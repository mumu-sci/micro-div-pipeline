# =============================================================================
# R/30_null_models.R -- 群落组装零模型引擎（betaNTI + RCbray）
# -----------------------------------------------------------------------------
# Stegen et al. (2013, ISME J; 2015, Front Microbiol) 联合零模型框架。
#
# betaNTI（系统发育 beta 周转的标准化效应量）
#   · 观察值：丰度加权 betaMNTD = picante::comdistnt(abundance.weighted = TRUE)
#     （comdistnt 的加权口径与 nifH 09_Supp_assembly.R 的 MNTD 手写实现同构：
#       对两样本各自取"到另一样本任一 OTU 的最小系统发育距离"，按各自相对丰度加权）
#   · 零模型：taxa.labels -- 每次打乱 ASV 标签（picante::taxaShuffle 对 cophenetic
#     距离矩阵行列施加同一置换）。本数据集每样本多为单例，"按 ASV 打乱标签"与
#     nifH 的 randomizeMatrix(null.model = "independentswap") 在该数据上数值等价
#     （差异见 params$bnti_null_model 记录）。
#   · betaNTI = (obs - null_mean) / null_sd，非有限值记 NA
#
# RCbray（Bray-Curtis 的标准化效应量）
#   · 观察值：vegan::vegdist(bray)
#   · 零模型：vegan::permatfull(fixedmar = CFG$assembly$fixedmar, times = runs)
#     -- 固定边际的列置换（保持各样本读长与丰富度）
#   · RCbray = (obs - null_mean) / null_sd，非有限值记 0，再截断到 [-1, 1]
#
# 三件套：
#   ① 分块：CFG$assembly$chunk 不为 NULL 时按"样本对"分块算完再合并
#      （每块的置换种子与块划分无关，故结果与不分块逐位一致；chunk 只影响内存
#        峰值与并行粒度）
#   ② 并行：md_parallel() / future.apply，每个任务一个由全包唯一种子派生的子种子，
#      故核数不影响结果；结束必 md_parallel_off()
#   ③ 缓存：整块结果走 md_cache(CFG, "30_null", ...)，落 cache/30_null.rds
#
# 大样本提速开关（CFG$assembly$subset_top）：只用丰度前 N 个 ASV -- **会改变结论**，
# 代价写进结果表 Note 与日志。
#
# 产物：cache/30_null.rds =
#   list(bnti = n×n 矩阵, rc = n×n 矩阵, pair_df = 长表, process = 因子向量,
#        params = list(...))            -- 见 MODULE_SPEC §8（31 / 33 依赖）
#       pair_df 列名固定为：
#   sample1, sample2, bray, bnti, rcbray, process, Zone1, Zone2, Vegetation1,
#   Vegetation2, Same_Zone, Same_Vegetation, Comparison_Type,
#   f1_1, f1_2, f2_1, f2_2, Same_f1, Same_f2
# =============================================================================

md_null_models <- function(CFG) {
  module <- "30"
  t_start <- Sys.time()
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少依赖包"))
  if (!md_precheck(CFG, module, need_response = FALSE, need_phylo = FALSE))
    return(md_result(module, "SKIP", "前置条件不满足"))

  md_load(c("ape", "vegan", "picante"))

  # ---- 参数（全部来自 CFG，缺则 %||% 默认）----
  A <- CFG$assembly %||% list()
  runs       <- as.integer(A$runs %||% 999L)
  if (!is.finite(runs) || runs < 1L) runs <- 999L
  null_model <- A$null_model %||% "independentswap"   # 仅记录进 params
  fixedmar   <- A$fixedmar %||% "columns"             # RCbray 的固定边际
  bnti_thr   <- A$bnti_thr %||% 2
  rc_thr     <- A$rcbray_thr %||% 0.95
  ncores     <- A$ncores %||% CFG$opts$ncores %||% 4L
  chunk      <- A$chunk
  subset_top <- A$subset_top
  factors    <- CFG$design$factors %||% character(0)
  smp_col    <- CFG$data$sample_col %||% "sample_id"
  detail     <- isTRUE(CFG$table$detail)
  th         <- local30_thresholds(CFG)

  if (!identical(null_model, "independentswap"))
    md_warn(CFG, "assembly$null_model = ", null_model, "：已按 CFG 记录，但 betaNTI 的",
            "零模型口径固定为 taxa.labels（cophenetic 标签置换）")

  # ---- 数据 ----
  m    <- md_core_asv(CFG, "assembly")
  meta <- md_core_meta(CFG)
  m    <- m[, colSums(m) > 0, drop = FALSE]
  ids  <- colnames(m)
  if (ncol(m) < 4L) {
    md_skip(CFG, module, "非空样本少于 4 个，成对比较无意义")
    return(md_result(module, "SKIP", "样本过少"))
  }
  mt <- meta[match(ids, meta[[smp_col]]), , drop = FALSE]

  # ---- subset_top 代价开关 ----
  subset_note <- character(0)
  if (!is.null(subset_top) && length(subset_top) == 1 && is.finite(subset_top) &&
      subset_top > 0 && subset_top < nrow(m)) {
    top <- order(rowSums(m), decreasing = TRUE)[seq_len(as.integer(subset_top))]
    m <- m[top, , drop = FALSE]
    m <- m[, colSums(m) > 0, drop = FALSE]
    ids <- colnames(m); mt <- mt[match(ids, meta[[smp_col]]), , drop = FALSE]
    subset_note <- sprintf("仅用丰度前 %d 个 ASV（**改变原始口径，会改变结论**）",
                           as.integer(subset_top))
    md_warn(CFG, "assembly$subset_top = ", subset_top, " -- ", subset_note)
  }
  m <- m[rowSums(m) > 0, , drop = FALSE]
  n <- length(ids)

  # ---- 内存 / 耗时预估（打印出来）----
  tri_cnt <- n * (n - 1) / 2
  bytes_null <- as.numeric(tri_cnt) * runs * 8
  n_block <- length(local30_split_pairs(tri_cnt, chunk))
  est_rc <- runs * 0.016          # 经验值：vegdist(36x321) 约 16 ms/次
  md_log(CFG, "INFO", sprintf(
    "零模型规模：%d ASV x %d 样本 | 成对 %d | runs %d | 核 %d | chunk %s | subset_top %s",
    nrow(m), n, tri_cnt, runs, ncores, as.character(chunk %||% "NULL"),
    as.character(subset_top %||% "NULL")))
  md_log(CFG, "INFO", sprintf(
    "内存预估：null 距离数组 %.1f MB（每样本对 %.3f MB）；RCbray 分块数 %d；顺序耗时约 %.1f min",
    bytes_null / 1024^2, runs * 8 / 1024^2, n_block, est_rc / 60))

  # ---- 树 ----
  tree <- md_core_tree(CFG)
  tr <- if (is.null(tree)) NULL else local30_prune_tree(tree, rownames(m), CFG)
  bnti_ok <- !is.null(tr) && length(tr$tip.label) >= 3L
  if (!bnti_ok) {
    md_warn(CFG, "无可用树 / 匹配 tip 不足，betaNTI 记 NA（只算 RCbray）")
    tr <- NULL
  }

  md_seed(CFG, th$seed30)
  n_used <- md_parallel(CFG, ncores, seed = CFG$opts$seed + th$seed30)
  md_seed(CFG, th$seed30)
  obs_bc <- as.matrix(vegan::vegdist(t(m), method = "bray"))

  res <- md_cache(CFG, "30_null", {
    # ---- 1. betaNTI ----
    if (bnti_ok) {
      t0 <- Sys.time()
      bnti_mat <- local30_bnti(CFG, t(m), tr, runs, n_used, th)
      md_log(CFG, "INFO", "betaNTI 完成，用时 ",
             md_fmt_ms(0, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    } else {
      bnti_mat <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
    }

    # ---- 2. RCbray ----
    t0 <- Sys.time()
    rc_mat <- local30_rcbray(CFG, t(m), obs_bc, runs, n_used, fixedmar, chunk, th)
    md_log(CFG, "INFO", "RCbray 完成，用时 ",
           md_fmt_ms(0, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

    # ---- 3. 长表 + 过程判定 ----
    loc <- local30_pair_indices(n)
    idx <- loc$row + (loc$col - 1) * n
    s1 <- ids[loc$row]; s2 <- ids[loc$col]
    bnti_v <- bnti_mat[idx]
    rc_v   <- rc_mat[idx]
    proc <- local30_process(bnti_v, rc_v, bnti_thr, rc_thr)
    proc <- factor(proc, levels = c("Variable selection", "Homogeneous selection",
                                    "Dispersal limitation", "Homogenizing dispersal",
                                    "Drift"))

    pdf <- data.frame(sample1 = s1, sample2 = s2,
                      bray = obs_bc[idx],
                      bnti = bnti_v, rcbray = rc_v, process = as.character(proc),
                      stringsAsFactors = FALSE)
    pdf$bray[!is.finite(pdf$bray)] <- NA_real_
    pdf <- cbind(pdf, local30_pair_meta(mt, ids, smp_col, factors, CFG$design$site_var))
    for (f in factors) {
      if (!f %in% names(mt)) next
      v <- as.character(mt[[f]][match(ids, mt[[smp_col]])])
      pdf[[paste0(f, "1")]] <- v[loc$row]
      pdf[[paste0(f, "2")]] <- v[loc$col]
    }
    pdf$Same_Zone <- pdf$Same_f1
    pdf$Same_Vegetation <- pdf$Same_f2
    pdf$Comparison_Type <- local30_comparison_type(pdf$Same_f1, pdf$Same_f2)
    pair_df <- pdf
    names(proc) <- as.character(pair_df$sample1)

    notes <- subset_note
    if (!bnti_ok) notes <- c(notes, "无树/树 tip 不足：betaNTI 全部为 NA")

    list(bnti = bnti_mat, rc = rc_mat, pair_df = pair_df, process = proc,
         params = list(
           runs = runs, null_model = null_model, bnti_null_model = "taxa.labels",
           fixedmar = fixedmar, bnti_thr = bnti_thr, rcbray_thr = rc_thr,
           ncores = n_used, chunk = chunk %||% NA_integer_,
           subset_top = subset_top %||% NA_integer_,
           n_asv = nrow(m), n_sample = n, n_pair = nrow(pair_df),
           has_tree = bnti_ok, seed = CFG$opts$seed,
           mem_est_mb = round(bytes_null / 1024^2, 1),
           est_rc_min = round(est_rc / 60, 1), n_block = n_block,
           notes = notes,
           generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  }, deps = md_dep(list(runs = runs, null_model = null_model, fixedmar = fixedmar,
                        chunk = chunk, subset_top = subset_top,
                        n_asv = nrow(m), ids = ids, seed = CFG$opts$seed)))

  md_parallel_off()

  bnti_mat <- res$bnti; rc_mat <- res$rc; pair_df <- res$pair_df; proc <- res$process
  params <- res$params
  md_log(CFG, "INFO", "betaNTI 非NA ", sum(is.finite(pair_df$bnti)), "/", nrow(pair_df),
         " | RCbray 非NA ", sum(is.finite(pair_df$rcbray)), "/", nrow(pair_df))
  print(table(proc))

  # ---- 4. 图：betaNTI + RCbray 直方图（标注锚在面板内）----
  # 画法逐项对照：
  #   面板 A ← nifH 11_shared.R §4 :: plot_bnti_hist()（Fig3C 与 FigS5A 共用同一函数）
  #              bins = 60、填充 steelblue、阈值线 ±2 红色虚线、两侧过程名标注锚在面板内部
  #              （plot_bnti_hist 的标注色 #2166AC / #E31A1C 正是 md_hist() 的默认值）
  #   面板 B ← nifH 09_Supp_assembly.R §5 的 RCbray 直方图（FigS5B）
  #              bins = 60、填充 steelblue、阈值线 ±0.95、标注色 #4DAF4A / #FF7F00
  # 用基元 md_hist()（13_plot_types.R）承载直方图/参考线/两侧标注；
  # md_hist() 没有暴露「阈值线颜色与线宽」「直方图透明度」这几个参数，
  # 故按原脚本那几行补一层 geom_vline()（样式常量登记在 local30_hist_style()）。
  hsty  <- local30_hist_style(CFG)
  hpt   <- md_pt(CFG)
  d_bnti <- data.frame(bnti = pair_df$bnti[is.finite(pair_df$bnti)])
  d_rc   <- data.frame(rc = pair_df$rcbray[is.finite(pair_df$rcbray)])

  if (nrow(d_bnti) >= 3) {
    pA <- md_hist(d_bnti, x = "bnti", bins = hsty$bins, ref = NULL, fill = hsty$fill,
                  left_label = "Homogeneous\nselection",
                  right_label = "Variable\nselection",
                  xlab = expression(beta * "NTI"), ylab = "Frequency",
                  tag = "A", pt = hpt) +
      local30_thr_line(c(-bnti_thr, bnti_thr), hsty) +
      ggplot2::labs(title = "Overall")
  } else {
    pA <- md_empty_plot(CFG, "beta*NTI unavailable")
  }

  if (nrow(d_rc) >= 3) {
    pB <- md_hist(d_rc, x = "rc", bins = hsty$bins, ref = NULL, fill = hsty$fill,
                  left_label = "Homogenizing\ndispersal",
                  right_label = "Dispersal\nlimitation",
                  left_col = hsty$rc_left, right_col = hsty$rc_right,
                  xlab = "RCbray", ylab = "Frequency", tag = "B", pt = hpt) +
      local30_thr_line(c(-rc_thr, rc_thr), hsty) +
      ggplot2::labs(title = "Overall")
  } else {
    pB <- md_empty_plot(CFG, "RCbray unavailable")
  }
  p <- md_compose(CFG, list(A = pA, B = pB), design = "AB", guides = "collect",
                  legend = "none")
  md_save_figure(CFG, p, "30_null_models", height = 90)

  # ---- 5. 表 ----
  bnti_stat <- local30_stat_rows(pair_df$bnti, "betaNTI", bnti_thr)
  rc_stat   <- local30_stat_rows(pair_df$rcbray, "RCbray", rc_thr)
  proc_tab <- as.data.frame(table(process = proc), stringsAsFactors = FALSE)
  names(proc_tab)[2] <- "n"
  proc_tab$prop <- proc_tab$n / sum(proc_tab$n)
  proc_tab <- proc_tab[proc_tab$n > 0, , drop = FALSE]
  proc_tab <- proc_tab[order(-proc_tab$n), , drop = FALSE]
  rownames(proc_tab) <- NULL

  f1 <- factors[1] %||% NA_character_
  f2 <- if (length(factors) > 1) factors[2] else NA_character_
  note <- data.frame(item = c(
    "betaNTI 口径", "RCbray 口径", "过程判定阈值", "成对比较", "已知代价 / 警示",
    "数据特征"),
    value = c(
      paste0("picante::comdistnt(abundance.weighted = TRUE) + 零模型 taxa.labels",
             "（cophenetic 标签置换），runs = ", runs),
      paste0("vegan::vegdist(bray) + vegan::permatfull(fixedmar = '", fixedmar,
             "', times = ", runs, ")"),
      sprintf("|betaNTI| > %g 判定确定性（< -%g 同质选择，> %g 异质选择）；|RCbray| > %g 判定扩散（> %g 扩散限制，< -%g 均质化扩散）；其余为漂变",
              bnti_thr, bnti_thr, bnti_thr, rc_thr, rc_thr, rc_thr),
      sprintf("upper.tri 长表，n = %d 对（%d 样本）", nrow(pair_df), n),
      paste(c(params$notes %||% character(0),
              "rc / bnti 矩阵数值已截断/置 NA：rc 截断到 [-1, 1]，SD = 0 的样本对 bnti 记 NA",
              "RCbray 分块只切样本对、不切 runs（切 runs 会打断 RNG 流破坏可复现性），",
              "故并发度上限 = 分块数；样本数很大时请显式设 chunk"),
            collapse = "；"),
      paste0("本数据集极端单例化，相邻样本对稀少 -> null_sd 可能为 0，",
             "此时 SES 记 NA / RCbray 记 0（不是计算失败）")),
    stringsAsFactors = FALSE)

  sheets <- list(
    Process_summary = proc_tab,
    BetaNTI_stats   = bnti_stat,
    RCbray_stats    = rc_stat,
    By_factor1      = local30_proc_by(pair_df, if (!is.na(f1)) paste0(f1, "1") else NULL,
                                      paste0("Same_f", 1)),
    By_factor2      = local30_proc_by(pair_df, if (!is.na(f2)) paste0(f2, "1") else NULL,
                                      paste0("Same_f", 2)),
    Note            = note)
  if (detail) sheets$Pairwise <- pair_df
  tp <- md_write_xlsx(CFG, "30_null_models", sheets,
                      meta_extra = list(runs = runs, bnti_null_model = "taxa.labels",
                                        fixedmar = fixedmar, n_pair = nrow(pair_df)))

  # 距离矩阵单独存 rds（绝不并进 data.frame）
  saveRDS(list(bnti = bnti_mat, rc = rc_mat, ids = ids),
          md_path(CFG, "cache", "30_null_matrices.rds"))

  md_log(CFG, "OK", "模块 30 用时 ",
         md_fmt_ms(0, as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  md_result(module, "OK", tables = basename(tp), figures = "figures/30_null_models.pdf",
            notes = sprintf("runs = %d | 成对 %d | betaNTI 非NA %d",
                            runs, nrow(pair_df), sum(is.finite(pair_df$bnti))))
}

# =============================================================================
# local30_ 辅助函数（共享库缺的，写在本文件内，不动 R/shared/）
# =============================================================================

#' 本模块用的 seed offset（可被 CFG$opts$md_seed_offsets 覆盖）
local30_thresholds <- function(CFG) {
  th <- CFG$opts$md_seed_offsets %||% list()
  list(seed30 = th$seed30 %||% 301L, seed30b = th$seed30b %||% 302L)
}

#' 直方图样式常量 —— **逐字照抄 nifH 原脚本**：
#'   bins = 60、fill = "steelblue"、阈值线 "red" / linewidth 0.8、RCbray 标注色
#'   #4DAF4A / #FF7F00 ← 11_shared.R::plot_bnti_hist() 与 09_Supp_assembly.R §5。
#' 这几个常量在原脚本里是内联字面量，12_plot_theme.R 没有对应键，
#' 故在此统一登记（并允许 CFG$figure$style 覆盖），模块内别再散落字面量。
local30_hist_style <- function(CFG) {
  s <- CFG$figure$style %||% list()
  list(bins     = as.integer(s$hist_bins %||% 60L),
       fill     = s$hist_fill        %||% "steelblue",
       thr_col  = s$threshold_colour %||% "red",
       thr_lw   = s$threshold_lw     %||% 0.8,
       rc_left  = s$rcbray_left_col  %||% "#4DAF4A",
       rc_right = s$rcbray_right_col %||% "#FF7F00")
}

#' 阈值线图层（照抄 09_Supp_assembly.R §5 的 geom_vline 那两行）
local30_thr_line <- function(xint, hsty) {
  ggplot2::geom_vline(xintercept = xint, linetype = "dashed",
                      colour = hsty$thr_col, linewidth = hsty$thr_lw)
}

#' 树修剪：保留与 ASV 的交集，确保 rooted + 二叉化
local30_prune_tree <- function(tree, asv_ids, CFG = NULL) {
  if (is.null(tree) || length(tree$tip.label) < 3L) return(NULL)
  tr <- tree
  tr$tip.label <- gsub('^["\']|["\']$', "", trimws(tr$tip.label))
  common <- intersect(asv_ids, tr$tip.label)
  if (length(common) < 3L) return(NULL)
  if (length(common) < length(tr$tip.label)) tr <- ape::keep.tip(tr, common)
  if (!is.null(tr$edge.length) && any(!is.finite(tr$edge.length)))
    tr$edge.length[!is.finite(tr$edge.length)] <- 0
  tr <- tryCatch(ape::multi2di(tr), error = function(e) tr)
  if (!ape::is.rooted(tr)) {
    rm <- CFG$phylo$root_method %||% "midpoint"
    tr <- tryCatch({
      if (identical(rm, "none")) tr
      else if (identical(rm, "outgroup")) local30_root_outgroup(tr)
      else if (requireNamespace("phangorn", quietly = TRUE)) phangorn::midpoint(tr)
      else local30_root_outgroup(tr)
    }, error = function(e) local30_root_outgroup(tr))
  }
  tr
}

local30_root_outgroup <- function(tr) {
  tryCatch({
    te <- which(tr$edge[, 2] <= length(tr$tip.label))
    if (!length(te)) return(tr)
    longest <- tr$edge[te, 2][which.max(tr$edge.length[te])]
    ape::root(tr, outgroup = longest, resolve.root = TRUE)
  }, error = function(e) tr)
}

#' upper.tri 的线性下标（row / col），与 nifH 09_Supp_assembly.R 的
#' which(upper.tri(m), arr.ind = TRUE) 完全同序
local30_pair_indices <- function(n) {
  if (n < 2L) return(data.frame(row = integer(0), col = integer(0)))
  rows <- integer(0); cols <- integer(0)
  for (i in seq_len(n - 1L)) {
    k <- n - i
    rows <- c(rows, rep.int(i, k))
    cols <- c(cols, (i + 1L):n)
  }
  data.frame(row = rows, col = cols)
}

#' 由全包唯一种子派生子种子序列
local30_perm_seeds <- function(seed, offset, n) {
  as.integer(seed) + as.integer(offset) + seq_len(as.integer(n))
}

#' 过程判定（nifH determine_process 的同序逻辑，向量化）
local30_process <- function(bnti, rc, bnti_thr, rc_thr) {
  out <- rep("Drift", length(bnti))
  okr <- is.finite(rc)
  out[okr & rc >  rc_thr] <- "Dispersal limitation"
  out[okr & rc < -rc_thr] <- "Homogenizing dispersal"
  okb <- is.finite(bnti)
  out[okb & bnti < -bnti_thr] <- "Homogeneous selection"
  out[okb & bnti >  bnti_thr] <- "Variable selection"
  out
}

#' 成对分组元数据（f1_1/Same_f1/... + Site1/Site2）。
#' 列顺序/语义与 22 的实现一致；这里用 local30_pair_indices()（upper.tri 同序），
#' 故分块/分文件都不会错位。
local30_pair_meta <- function(meta, ids, smp_col, factors, site_var = NULL) {
  n <- length(ids)
  loc <- local30_pair_indices(n)
  # 坑：0 列 data.frame 的 nrow() 是 0，必须先固定行数，否则第一次列赋值报
  # "replacement has n rows, data has 0"
  out <- data.frame(row.names = seq_len(nrow(loc)))
  for (k in seq_along(factors[seq_len(min(2, length(factors)))])) {
    f <- factors[k]
    v <- if (f %in% names(meta)) as.character(meta[[f]][match(ids, meta[[smp_col]])])
         else rep(NA_character_, n)
    out[[paste0("f", k, "_1")]] <- v[loc$row]
    out[[paste0("f", k, "_2")]] <- v[loc$col]
    out[[paste0("Same_f", k)]]  <- out[[paste0("f", k, "_1")]] == out[[paste0("f", k, "_2")]]
  }
  sv <- if (!is.null(site_var) && site_var %in% names(meta))
    as.character(meta[[site_var]][match(ids, meta[[smp_col]])]) else rep(NA_character_, n)
  out$Site1 <- sv[loc$row]; out$Site2 <- sv[loc$col]
  out
}

#' Comparison_Type（nifH 口径：Within_Group / Across_Vegetation / Across_Zone / Across_Both）
local30_comparison_type <- function(same1, same2) {
  out <- rep("Across_Both", length(same1))
  out[same1 & !same2] <- "Across_Vegetation"
  out[!same1 & same2] <- "Across_Zone"
  out[same1 & same2]  <- "Within_Group"
  out[is.na(same1) | is.na(same2)] <- NA_character_
  out
}

#' 统计汇总行（阈值判定用 CFG 的阈值，不写死）
local30_stat_rows <- function(x, label, thr = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(data.frame(metric = label, Mean = NA_real_, SD = NA_real_, Median = NA_real_,
                      Min = NA_real_, Max = NA_real_, N_total = 0L,
                      N_low = 0L, N_high = 0L, N_within = 0L, stringsAsFactors = FALSE))
  }
  data.frame(metric = label, Mean = mean(x), SD = stats::sd(x), Median = stats::median(x),
             Min = min(x), Max = max(x), N_total = length(x),
             N_low = sum(x < -thr), N_high = sum(x > thr),
             N_within = sum(abs(x) <= thr), stringsAsFactors = FALSE)
}

#' 按因子（同因子对）汇总过程比例
local30_proc_by <- function(pair_df, group_col, same_col) {
  if (is.null(group_col) || !group_col %in% names(pair_df)) return(NULL)
  d <- pair_df[!is.na(pair_df[[group_col]]), , drop = FALSE]
  if (same_col %in% names(pair_df)) d <- d[which(d[[same_col]] %in% TRUE), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  keys <- sort(unique(as.character(d[[group_col]])))
  out <- do.call(rbind, lapply(keys, function(k) {
    x <- d[as.character(d[[group_col]]) == k, , drop = FALSE]
    t <- table(factor(x$process, levels = c("Variable selection", "Homogeneous selection",
                                            "Dispersal limitation", "Homogenizing dispersal",
                                            "Drift")))
    data.frame(group = k, process = names(t), n = as.integer(t),
               prop = as.numeric(t) / sum(t), n_pairs = nrow(x), stringsAsFactors = FALSE)
  }))
  out[out$n > 0, , drop = FALSE]
}

#' betaNTI 零模型：并行 + 逐任务子种子。
#' 每个任务打乱 cophenetic 矩阵的标签（taxa.labels），对全部样本算丰度加权 betaMNTD。
#' @return n×n 矩阵；无树时全 NA
local30_bnti <- function(CFG, samp, tr, runs, ncores, th) {
  n <- nrow(samp)
  ids <- rownames(samp)
  coph <- as.matrix(stats::cophenetic(tr))
  coph[!is.finite(coph)] <- NA_real_
  # comdistnt 返回 dist；转成完整对称矩阵并把对角（自比较）置 0。
  # 坑 1：matrix(as.numeric(dist), n, n) 会把 n(n-1)/2 个值循环填充，对角变成别的
  #        样本对的距离、矩阵也不再对称；必须显式填上三角再镜像。
  # 坑 2：stats::as.dist(完整矩阵) 在矩阵含 NA 时会**丢掉全 NA 的行/列**，
  #        使取值向量长度 < n(n-1)/2，与 upper.tri 数量不匹配（表现为整列 NA）。
  #        故这里一律用 local30_upper_vals() 按 upper.tri 全量取值。
  obs <- as.matrix(picante::comdistnt(samp, coph, abundance.weighted = TRUE))
  obs[!is.finite(obs)] <- NA_real_
  diag(obs) <- NA_real_
  obs_v <- local30_upper_vals(obs)
  seeds <- local30_perm_seeds(CFG$opts$seed, th$seed30, runs)
  one <- function(r) {
    set.seed(seeds[r])
    m <- as.matrix(picante::comdistnt(samp, picante::taxaShuffle(coph),
                                      abundance.weighted = TRUE))
    diag(m) <- NA_real_
    local30_upper_vals(m)
  }
  mat <- if (ncores > 1L && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(seq_len(runs), one, future.seed = NULL)
  } else {
    lapply(seq_len(runs), one)
  }
  mat <- do.call(rbind, mat)
  mu  <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, stats::sd, na.rm = TRUE)
  zvals <- (obs_v - mu) / sdv
  out <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  out[lower.tri(out)] <- zvals
  out <- out + t(out)
  out[!is.finite(out)] <- NA_real_
  out
}

#' 取对称矩阵的上三角（行优先，与 which(upper.tri(.), arr.ind=TRUE) 同序）。
#' 不用 stats::as.dist()：它在含 NA 时会丢掉全 NA 的行，长度会缩水。
local30_upper_vals <- function(m) {
  n <- nrow(m)
  if (n < 2L) return(numeric(0))
  idx <- local30_pair_indices(n)
  m[cbind(idx$row, idx$col)]
}

#' RCbray 零模型：vegan::permatfull 固定边际 + 分块并行。
#' 关键点：每次 permatfull(times = runs) 对应一串固定种子，**与块的划分无关**，
#' 因此分块只是把样本对分批，不改变任何数值；块内返回 (sum, sumsq, n) 充分统计量，
#' 块间合并得到精确 mean / sd（结果与不分块逐位一致）。
#' 说明：不进一步切分 runs -- 切分 runs 会把 RNG 流打断，破坏可复现性；
#' 因此数据规模很大时并发度受"每样本对一块"限制（见 params 的 Note）。
local30_rcbray <- function(CFG, samp, obs_bc, runs, ncores, fixedmar, chunk, th) {
  n <- nrow(samp)
  loc <- local30_pair_indices(n)
  np <- nrow(loc)
  obs_v <- obs_bc[loc$row + (loc$col - 1) * n]

  blocks <- local30_split_pairs(np, chunk)
  nb <- length(blocks)
  seed0 <- as.integer(CFG$opts$seed) + as.integer(th$seed30b)
  md_log(CFG, "INFO", "RCbray：", np, " 个样本对 -> ", nb, " 块（每块约 ",
         round(np / nb), " 对），每块 ", runs, " 次置换（种子与块无关，可复现）")

  one_block <- function(b) {
    set.seed(seed0)                      # 同一个种子 -> 每一块看到同一批置换
    ix <- blocks[[b]]
    summ <- matrix(0, length(ix), runs)
    pf <- vegan::permatfull(samp, fixedmar = fixedmar, times = runs)
    for (k in seq_len(runs)) {
      d <- vegan::vegdist(pf$perm[[k]], method = "bray")
      summ[, k] <- d[loc$row[ix] + (loc$col[ix] - 1) * n]
    }
    summ[!is.finite(summ)] <- NA_real_
    list(ix = ix, sum = rowSums(summ, na.rm = TRUE),
         sq = rowSums(summ^2, na.rm = TRUE), n = rowSums(!is.na(summ)))
  }

  parts <- if (ncores > 1L && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(seq_len(nb), one_block, future.seed = NULL)
  } else {
    lapply(seq_len(nb), one_block)
  }

  # 块间合并充分统计量 -> 精确 mean / sd
  s1 <- s2 <- numeric(np); cnt <- integer(np)
  for (p in parts) {
    s1[p$ix]  <- p$sum
    s2[p$ix]  <- p$sq
    cnt[p$ix] <- p$n
  }
  means <- ifelse(cnt > 0, s1 / pmax(cnt, 1), 0)
  sds   <- sqrt(pmax((s2 - cnt * means^2) / pmax(cnt - 1, 1), 0))
  sds[!is.finite(sds)] <- 0

  rc <- (obs_v - means) / sds
  rc[!is.finite(rc)] <- 0
  rc <- pmax(pmin(rc, 1), -1)
  m <- matrix(0, n, n, dimnames = list(rownames(samp), rownames(samp)))
  m[loc$row + (loc$col - 1) * n] <- rc
  m[loc$col + (loc$row - 1) * n] <- rc
  diag(m) <- 0
  m
}

#' 把 np 个样本对切成块（chunk = 每块的对数；NULL = 全部一块）
local30_split_pairs <- function(np, chunk) {
  if (np < 1L) return(list())
  if (is.null(chunk) || !is.finite(chunk) || chunk < 1L) return(list(seq_len(np)))
  k <- as.integer(chunk)
  unname(split(seq_len(np), ceiling(seq_len(np) / k)))
}
