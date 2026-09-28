# =============================================================================
# R/22_phylo_beta.R -- 系统发育 beta 多样性分解（turnover / nestedness）
# -----------------------------------------------------------------------------
# 方法（与 nifH 论文 11_shared.R §3 compute_phylo_sorensen + 14_Supp_phylo_beta.R 对齐）：
#   · 先在 prev >= CFG$phylo$prev_min 的**非单例子集**上做计算 -- 全树上直接跑
#     phylo.beta.pair 是 O(n²) 内存爆炸（nifH 实测卡死）；子集内为空的样本剔除。
#   · 树：确保 rooted -> multi2di（多分叉按零长枝拆开），再 keep.tip 到子集。
#   · betapart::phylo.beta.pair(index.family = "sorensen")
#       -> phylo.beta.sor（总）/ phylo.beta.sim（周转）/ phylo.beta.sne（嵌套）
#   · 分类学对照：betapart::beta.pair(index.family = "sorensen")$beta.sim（同子集）
#       -- 分类学 vs 系统发育 的 Mantel（Spearman）与成对 Spearman +
#          **站点对整簇 bootstrap 置信区间**（md_boot_rho，控制成对伪重复）。
#
# 产物：results/22_phylo_beta.xlsx + figures/22_phylo_beta.pdf
#      （本模块不写 cache/22_*，跨模块契约里没有它；数值表即交付物）
#
# 无树 / phylo$enabled = FALSE 时 SKIP。
# =============================================================================

md_phylo_beta <- function(CFG) {
  module <- "22"
  t_start <- Sys.time()
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少依赖包"))
  if (!md_precheck(CFG, module, need_response = FALSE, need_phylo = TRUE))
    return(md_result(module, "SKIP", "前置条件不满足"))

  md_load(c("ape", "betapart", "vegan"))

  prev_min <- CFG$phylo$prev_min %||% 2
  factors  <- CFG$design$factors %||% character(0)
  smp_col  <- CFG$data$sample_col %||% "sample_id"
  perms    <- CFG$beta_tests$permutations %||% 999
  B        <- CFG$response_models$bootstrap_B %||% 1000

  m    <- md_core_asv(CFG, "phylo")
  meta <- md_core_meta(CFG)
  tree <- md_core_tree(CFG)
  if (is.null(tree)) {
    md_skip(CFG, module, "无可用系统发育树")
    return(md_result(module, "SKIP", "无系统发育树"))
  }

  # ---- 1. 树修剪 ----
  tr <- local22_phylo_prune(tree, rownames(m), CFG)
  if (is.null(tr) || length(tr$tip.label) < 3L) {
    md_skip(CFG, module, "ASV 与树 tip 匹配后类群不足 3 个")
    return(md_result(module, "SKIP", "树 tip 匹配不足"))
  }
  m <- m[tr$tip.label, , drop = FALSE]

  # ---- 2. prev >= prev_min 的非单例子集 ----
  pa   <- t((m > 0) * 1L)                 # samples x taxa
  prev <- colSums(pa)
  keep_tax <- names(prev)[prev >= prev_min]
  if (length(keep_tax) < 3L) {
    md_skip(CFG, module,
            paste0("出现 >= ", prev_min, " 个样本的 ASV 不足 3 个（本数据集极端单例化）"))
    return(md_result(module, "SKIP", "prev 过滤后类群不足"))
  }
  pa_sub <- pa[, keep_tax, drop = FALSE]
  rs <- rowSums(pa_sub)
  dropped <- names(rs)[rs == 0]
  pa_sub <- pa_sub[rs > 0, , drop = FALSE]
  storage.mode(pa_sub) <- "numeric"
  ids <- rownames(pa_sub)
  if (length(ids) < 3L) {
    md_skip(CFG, module, "prev 子集内非空样本少于 3 个")
    return(md_result(module, "SKIP", "子集内样本过少"))
  }
  tr_sub <- ape::keep.tip(tr, intersect(tr$tip.label, colnames(pa_sub)))
  tr_sub <- tryCatch(ape::multi2di(tr_sub), error = function(e) tr_sub)
  if (!ape::is.rooted(tr_sub)) tr_sub <- local22_root_outgroup(tr_sub)
  common <- intersect(tr_sub$tip.label, colnames(pa_sub))
  pa_sub <- pa_sub[, common, drop = FALSE]
  md_log(CFG, "INFO", "系统发育 beta：prev >= ", prev_min, " 子集 = ", length(common),
         " 类群 x ", length(ids), " 样本（原 ", nrow(m), " 类群 x ", ncol(m), " 样本；",
         "剔除空样本 ", length(dropped), " 个）")

  # ---- 3. 系统发育 Sørensen 分解 ----
  pbp <- tryCatch(betapart::phylo.beta.pair(pa_sub, tr_sub, index.family = "sorensen"),
                  error = function(e) {
                    md_warn(CFG, "phylo.beta.pair 失败：", conditionMessage(e)); NULL })
  if (is.null(pbp)) {
    md_skip(CFG, module, "betapart::phylo.beta.pair 计算失败")
    return(md_result(module, "SKIP", "phylo.beta.pair 失败"))
  }
  PO <- as.matrix(pbp$phylo.beta.sor)
  PT <- as.matrix(pbp$phylo.beta.sim)
  PN <- as.matrix(pbp$phylo.beta.sne)
  PO[!is.finite(PO)] <- NA_real_; PT[!is.finite(PT)] <- NA_real_; PN[!is.finite(PN)] <- NA_real_

  # ---- 4. 分类学对照（同子集，Sørensen 周转）----
  TT <- tryCatch({
    tb <- betapart::beta.pair(pa_sub, index.family = "sorensen")
    as.matrix(tb$beta.sim)
  }, error = function(e) NULL)
  if (!is.null(TT)) TT[!is.finite(TT)] <- NA_real_

  # ---- 5. 区域对聚合（upper.tri 与跨模块 pair_df 同序）----
  f1 <- factors[1] %||% NA_character_
  tri <- which(upper.tri(PO), arr.ind = TRUE)
  long <- data.frame(
    sample1 = ids[tri[, 1]], sample2 = ids[tri[, 2]],
    phylo_total = PO[tri], phylo_turnover = PT[tri], phylo_nested = PN[tri],
    stringsAsFactors = FALSE)
  if (!is.null(TT)) long$tax_turnover <- TT[tri]

  # 因子对列（CFG$design$factors 前两个）+ 站点对标签
  pair_meta <- local22_pair_meta(meta, ids, smp_col, factors,
                                 site_var = CFG$design$site_var)
  if (nrow(pair_meta) != nrow(long))
    stop("成对元数据行数不匹配：pair_meta = ", nrow(pair_meta), " vs long = ", nrow(long),
         call. = FALSE)
  long <- cbind(long, pair_meta)
  long$SitePair <- local22_site_pair(long$Site1, long$Site2)
  long$RegionPair <- if (!is.na(f1))
    mapply(function(a, b) paste(sort(c(a, b)), collapse = " - "),
           long$f1_1, long$f1_2, USE.NAMES = FALSE) else "All"
  long$RegionPair[is.na(long$f1_1) | is.na(long$f1_2)] <- NA_character_

  agg <- local22_pair_aggregate(long)
  if (is.null(agg) || nrow(agg) == 0) {
    md_skip(CFG, module, "无有效区域对可聚合")
    return(md_result(module, "SKIP", "无有效区域对"))
  }

  # ---- 6. 分类 vs 系统发育：Mantel + 成对 Spearman（整簇 bootstrap）----
  mantel_df <- NULL
  boot_df <- NULL
  if (!is.null(TT)) {
    ok <- is.finite(TT) & is.finite(PT)
    if (sum(ok) >= 10L) {
      mt <- md_mantel(CFG, stats::as.dist(TT), stats::as.dist(PT),
                      method = CFG$env$cor_method %||% "spearman",
                      permutations = perms,
                      label = sprintf("Tax turnover vs Phylo turnover (prev>=%d, n=%d)",
                                      prev_min, length(ids)))
      if (!is.null(mt)) mantel_df <- mt
      bd <- long[is.finite(long$tax_turnover) & is.finite(long$phylo_turnover), , drop = FALSE]
      if (nrow(bd) >= 5L) {
        boot_df <- tryCatch(
          md_boot_rho(CFG, bd, x = "tax_turnover", y = "phylo_turnover",
                      cluster = "SitePair", B = B,
                      label = "Tax turnover vs Phylo turnover (pairwise)"),
          error = function(e) NULL)
      }
    }
  }

  means_df <- data.frame(
    metric = c("phylo.total", "phylo.turnover", "phylo.nested",
               "tax.turnover(subset)", "n_samples", "n_taxa", "prev_min"),
    value = c(mean(long$phylo_total, na.rm = TRUE), mean(long$phylo_turnover, na.rm = TRUE),
              mean(long$phylo_nested, na.rm = TRUE),
              if (is.null(TT)) NA_real_ else mean(long$tax_turnover, na.rm = TRUE),
              length(ids), length(common), prev_min),
    stringsAsFactors = FALSE)

  # ---- 7. 图 ----
  pal <- md_pal(CFG, unique(stats::na.omit(agg$RegionPair)), "group")
  ba <- rbind(
    data.frame(RegionPair = agg$RegionPair, component = "Turnover",
               value = agg$phylo_turnover, stringsAsFactors = FALSE),
    data.frame(RegionPair = agg$RegionPair, component = "Nestedness",
               value = agg$phylo_nested, stringsAsFactors = FALSE))
  ba$RegionPair <- factor(ba$RegionPair, levels = agg$RegionPair)
  ba$component  <- factor(ba$component, levels = c("Turnover", "Nestedness"))
  pA <- md_bar_stack(ba, x = "RegionPair", y = "value", fill = "component",
                     colors = c(Turnover = "#E15759", Nestedness = "#4E79A7"),
                     xlab = NULL, ylab = "Phylogenetic beta (Sorensen)", tag = "A",
                     x_angle = 30, legend_pos = "bottom", ncol = 1)

  plots <- list(A = pA)
  if (!is.null(TT)) {
    dd <- long[is.finite(long$tax_turnover) & is.finite(long$phylo_turnover), , drop = FALSE]
    dd$RegionPair <- factor(dd$RegionPair, levels = agg$RegionPair)
    pB <- md_scatter_fit(dd, x = "tax_turnover", y = "phylo_turnover",
                         color = "RegionPair", colors = pal, line = "lm",
                         xlab = "Taxonomic turnover (Sorensen)",
                         ylab = "Phylogenetic turnover (Sorensen)", tag = "B")
    lab <- c(if (!is.null(mantel_df))
      sprintf("Mantel r = %.3f, P = %s", mantel_df$statistic[1], md_fmt_p(mantel_df$p[1])),
      if (!is.null(boot_df))
        sprintf("pairwise rho = %.3f [%.3f, %.3f]", boot_df$rho[1], boot_df$ci_lo[1], boot_df$ci_hi[1]),
      sprintf("n = %d pairs", nrow(dd)))
    pB <- local22_add_annot(pB, lab, x = -Inf, y = Inf, hjust = -0.06, vjust = 1.5)
    plots$B <- pB
  }
  design <- if (length(plots) > 1) "AB" else "A"
  p <- md_compose(CFG, plots, design = design, guides = "collect", legend = "bottom")
  md_save_figure(CFG, p, "22_phylo_beta", height = 105)

  # ---- 8. 表 ----
  note <- data.frame(item = c(
    "子集口径", "树处理", "分解指数", "分类学对照",
    "Mantel", "bootstrap", "已知代价"),
    value = c(
      sprintf("prev >= %d 的非单例子集（原 %d 类群 -> %d 类群；子集内空样本剔除 %d 个；n = %d）",
              prev_min, nrow(m), length(common), length(dropped), length(ids)),
      paste0("确保 rooted（root_method = ", CFG$phylo$root_method %||% "midpoint",
             "）-> multi2di"),
      "betapart::phylo.beta.pair(index.family = 'sorensen')：sor = sim + sne",
      "betapart::beta.pair(index.family = 'sorensen')$beta.sim（同子集）",
      sprintf("vegan::mantel(method = '%s', permutations = %d)",
              CFG$env$cor_method %||% "spearman", perms),
      sprintf("md_boot_rho：按站点对整簇重抽 B = %d（控制成对数据伪重复）", B),
      paste0("prev 过滤会剔除大量单例 ASV（本数据集极端单例化），",
             "系统发育 beta 数值来自该子集，不能外推到全部 ASV")),
    stringsAsFactors = FALSE)

  sheets <- list(
    Region_pair   = agg,
    Phylo_means   = means_df,
    Pairwise      = if (isTRUE(CFG$table$detail)) long else NULL,
    Mantel        = mantel_df,
    Pairwise_rho  = boot_df,
    Note          = note)
  tp <- md_write_xlsx(CFG, "22_phylo_beta", sheets)

  md_log(CFG, "OK", "模块 22 用时 ",
         md_fmt_ms(0, as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  md_result(module, "OK", tables = basename(tp), figures = "figures/22_phylo_beta.pdf",
            notes = sprintf("prev>=%d 子集 %d 类群 x %d 样本", prev_min, length(common), length(ids)))
}

# -----------------------------------------------------------------------------
# local22_ 辅助函数
# -----------------------------------------------------------------------------

#' 修剪树：保留与给定 ASV 的交集，确保 rooted + 二叉化（与 21 同口径，各模块自持）
local22_phylo_prune <- function(tree, asv_ids, CFG = NULL) {
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
      else if (identical(rm, "outgroup")) local22_root_outgroup(tr)
      else if (requireNamespace("phangorn", quietly = TRUE)) phangorn::midpoint(tr)
      else local22_root_outgroup(tr)
    }, error = function(e) local22_root_outgroup(tr))
  }
  tr
}

local22_root_outgroup <- function(tr) {
  tryCatch({
    te <- which(tr$edge[, 2] <= length(tr$tip.label))
    if (!length(te)) return(tr)
    longest <- tr$edge[te, 2][which.max(tr$edge.length[te])]
    ape::root(tr, outgroup = longest, resolve.root = TRUE)
  }, error = function(e) tr)
}

#' 生成成对比较的分组元数据列（f1_1/f1_2/Same_f1 ...），与 30 的 pair_df 同列名；
#' 另外给站点对（Site1/Site2）便于整簇 bootstrap。upper.tri 顺序与 30 一致。
local22_pair_meta <- function(meta, ids, smp_col, factors, site_var = NULL) {
  n <- length(ids)
  tri <- which(upper.tri(matrix(0, n, n), diag = FALSE), arr.ind = TRUE)
  # 坑：0 列 data.frame 的 nrow() 是 0（不是 0 行也不是 n 行），必须先固定行数，
  # 否则第一次 out[[col]] <- <n 元向量> 会报 "replacement has n rows, data has 0"
  out <- data.frame(row.names = seq_len(nrow(tri)))
  for (k in seq_along(factors[seq_len(min(2, length(factors)))])) {
    f <- factors[k]
    v <- if (f %in% names(meta)) as.character(meta[[f]][match(ids, meta[[smp_col]])])
         else rep(NA_character_, n)
    out[[paste0("f", k, "_1")]] <- v[tri[, 1]]
    out[[paste0("f", k, "_2")]] <- v[tri[, 2]]
    out[[paste0("Same_f", k)]]  <- out[[paste0("f", k, "_1")]] == out[[paste0("f", k, "_2")]]
  }
  sv <- if (!is.null(site_var) && site_var %in% names(meta))
    as.character(meta[[site_var]][match(ids, meta[[smp_col]])]) else rep(NA_character_, n)
  out$Site1 <- sv[tri[, 1]]
  out$Site2 <- sv[tri[, 2]]
  out
}

#' 站点对标签（整簇 bootstrap 的聚类单位）
local22_site_pair <- function(s1, s2) {
  ifelse(is.na(s1) | is.na(s2), NA_character_,
         ifelse(s1 == s2, s1,
                paste(pmin(s1, s2), pmax(s1, s2), sep = "_")))
}

#' 把文字标注加到**图对象**上。
#' 坑：共享库的 md_add_p_annot() 只返回 ggplot2::annotate(...) 的结果（一个 ggproto
#' layer），忘了 `p +`，调用方拿到的是 LayerInstance，再拼图会报
#' "Cannot add <ggproto> objects together"。模块内自持正确版本，不改 shared/。
local22_add_annot <- function(p, labels, x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
                              pt = MD_PT, sep = "\n") {
  p + ggplot2::annotate("text", x = x, y = y, label = paste(labels, collapse = sep),
                        hjust = hjust, vjust = vjust, size = pt2gg(pt$ann))
}

#' 按 RegionPair 聚合系统发育 beta 分解（周转 / 嵌套 / 分类学周转）
local22_pair_aggregate <- function(long) {
  d <- long[!is.na(long$RegionPair), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  keys <- sort(unique(d$RegionPair))
  out <- do.call(rbind, lapply(keys, function(k) {
    x <- d[d$RegionPair == k, , drop = FALSE]
    data.frame(
      RegionPair     = k,
      phylo_total    = mean(x$phylo_total,    na.rm = TRUE),
      phylo_turnover = mean(x$phylo_turnover, na.rm = TRUE),
      phylo_nested   = mean(x$phylo_nested,   na.rm = TRUE),
      tax_turnover   = if ("tax_turnover" %in% names(x)) mean(x$tax_turnover, na.rm = TRUE) else NA_real_,
      n_pairs        = nrow(x),
      stringsAsFactors = FALSE)
  }))
  out$turnover_ratio <- out$phylo_turnover / pmax(out$phylo_total, .Machine$double.eps)
  out
}
