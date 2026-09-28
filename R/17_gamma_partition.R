# =============================================================================
# R/17_gamma_partition.R ---- α / β / γ 多样性分区（Hill 数刻度）
# -----------------------------------------------------------------------------
# 引擎（CFG$partition$engine）：iNEXT（分类学 Hill 数）| iNEXT.3D（系统发育 PD）| both
#
# 口径要点：
#   1. α 与 γ 取自**同一条 incidence 稀有化曲线**：
#        α = t = CFG$partition$t_alpha 个采样单元的稀有化 Hill 数（默认 1）
#        γ = t = t_gamma 的稀有化 Hill 数（主口径 A），
#        另报 AsyEst 渐近估计量（口径 B）与实际观测采样量下的值（口径 C）。
#      Hill 刻度上 γ = α x β 严格成立，故 α 占比 = 100 / β_mult（见 nifH TableS12 口径说明）。
#   2. t_gamma = NULL 时**自动判定** = 各组采样单元数的中位数；对达不到 t_gamma 的组
#      **显式标注**（NeverExcluded / Status 列），并用该组可支持的最大 t（观测到的 t）
#      代替 γ 口径 A，绝不静默剔除组（nifH 原码手工剔除 Coastal_Shrub / Nearshore_Shrub）。
#   3. iNEXT / iNEXT.3D 都对数据量有硬性下限：iNEXT.3D 明确要求每组观测物种数 >= 5，
#      iNEXT 在"没有 f 落在 1..nT-1 区间"的数据上会抛 NA/NaN。这些组**记 NA 并写入
#      NA_reason**，模块不报错（示范数据高度单例化，交互组普遍不可估计，是数据真实特征）。
#   4. 分组 = 每个因子水平 + 因子交互组合（交互层为主层），另加全池化 All 一行作参照。
#   5. 重计算模块（heavy=TRUE）：nboot 用 CFG$partition$nboot；缺值时退化为较小值并记警告。
#
# 输出：figures/17_gamma_partition.pdf
#       results/17_gamma_partition.xlsx（Partition / Coverage / Note）
# =============================================================================

md_gamma_partition <- function(CFG) {
  module <- "17"
  part <- CFG$partition %||% list()
  engine <- part$engine %||% "iNEXT"
  if (!md_need(CFG, "core", module)) return(md_result(module, "SKIP", "缺少核心包"))
  want_inext <- engine %in% c("iNEXT", "both")
  want_3d    <- engine %in% c("iNEXT.3D", "both")
  if (!md_precheck(CFG, module)) return(md_result(module, "SKIP", "前置条件不满足"))

  have_inext <- requireNamespace("iNEXT", quietly = TRUE)
  have_3d    <- requireNamespace("iNEXT.3D", quietly = TRUE)
  if (want_inext && !have_inext) {
    md_warn(CFG, "缺 iNEXT 包，退化为只用 iNEXT.3D")
    want_inext <- FALSE
  }
  if (want_3d && !have_3d) {
    md_warn(CFG, "缺 iNEXT.3D 包，退化为只用 iNEXT")
    want_3d <- FALSE
  }
  if (!want_inext && !want_3d) {
    md_skip(CFG, module, "partition$engine 指定的引擎都不可用（缺 iNEXT / iNEXT.3D）")
    return(md_result(module, "SKIP", "无可用引擎"))
  }
  if (!requireNamespace("iNEXT", quietly = TRUE) && !have_3d) {
    md_skip(CFG, module, "iNEXT 系列包都不可用")
    return(md_result(module, "SKIP", "缺包"))
  }

  m    <- md_core_asv(CFG, "partition")
  meta <- md_core_meta(CFG)
  if (is.null(m) || !nrow(m) || ncol(m) < 2) {
    md_skip(CFG, module, "分区口径矩阵为空")
    return(md_result(module, "SKIP", "矩阵不可用"))
  }
  pa <- t((m > 0) * 1L)                     # 样本 x ASV（iNEXT incidence_freq 用）
  sp <- t(pa)                               # ASV x 样本（iNEXT.3D incidence_raw 用）
  nboot  <- as.integer(part$nboot %||% 50)
  conf   <- part$conf %||% 0.95
  qs     <- part$q %||% c(0, 1, 2)
  qs     <- sort(unique(as.numeric(qs)))
  qs     <- qs[is.finite(qs) & qs >= 0]
  if (!length(qs)) qs <- c(0, 1, 2)
  t_alpha <- as.numeric(part$t_alpha %||% 1)
  n_singleton <- sum(colSums(pa) == 1)
  notes <- c(paste0("P/A 矩阵：", nrow(pa), " 样本 x ", ncol(pa), " ASV；",
                    n_singleton, " 个 ASV 为单例（",
                    round(100 * n_singleton / max(1, ncol(pa)), 1), "%）"))
  if (n_singleton / max(1, ncol(pa)) > 0.5)
    notes <- c(notes, "数据高度单例化：样本覆盖度低、AsyEst 外推倍数大是真实特征，需与 SC 同表阅读")

  # ---- 分组定义：因子水平 + 因子交互 + 全池化 ----
  facs <- CFG$design$factors %||% character(0)
  facs <- facs[facs %in% names(meta)]
  grp_list <- list()   # list(label, ids, Level, Group)
  for (f in facs) {
    lv <- CFG$design$ordered_levels[[f]] %||% sort(unique(as.character(meta[[f]])))
    for (x in lv) {
      ids <- which(as.character(meta[[f]]) == x)
      if (length(ids)) grp_list[[length(grp_list) + 1L]] <-
        list(Group = x, Level = f, ids = ids)
    }
  }
  if (length(facs) >= 2) {
    inter <- as.character(interaction(meta[, facs, drop = FALSE], sep = "_",
                                      drop = TRUE, lex.order = TRUE))
    for (x in unique(inter)) {
      ids <- which(inter == x)
      if (length(ids)) grp_list[[length(grp_list) + 1L]] <-
        list(Group = x, Level = paste(facs[1:2], collapse = " x "), ids = ids)
    }
    # 交互层排到前面（主层）
    ord <- order(vapply(grp_list, function(z) z$Level == paste(facs[1:2], collapse = " x "),
                        logical(1)), decreasing = TRUE)
    grp_list <- grp_list[ord]
  }
  grp_list[[length(grp_list) + 1L]] <- list(Group = "All", Level = "Pooled",
                                            ids = seq_len(nrow(pa)))
  n_units <- vapply(grp_list, function(z) length(z$ids), integer(1))

  # ---- t_gamma 自动判定 ----
  grp_n   <- n_units[vapply(grp_list, function(z) z$Level != "Pooled", logical(1))]
  t_gamma_auto <- stats::median(grp_n, na.rm = TRUE)
  t_gamma <- part$t_gamma %||% NULL
  t_auto  <- is.null(t_gamma)
  if (t_auto) t_gamma <- as.numeric(t_gamma_auto)
  t_gamma <- max(1, as.numeric(t_gamma))
  notes <- c(notes, paste0("t_gamma = ", t_gamma,
                           if (t_auto) paste0("（自动判定 = 各组采样单元数中位数，中位数值 ",
                                              t_gamma_auto, "）") else "（来自 CFG$partition$t_gamma）"))
  short <- sum(grp_n < t_gamma)
  if (short > 0)
    notes <- c(notes, paste0(short, " 个组的采样单元数 < t_gamma=", t_gamma,
                             "：这些组的 γ 口径 A 用该组可支持的最大 t 代替，并在 Status 列标注（不剔除）"))

  # ---- 分类学引擎（iNEXT 最慢，整块走缓存：改样式重出图时不再重算）----
  part_rows <- list()
  cov_rows  <- list()
  if (want_inext && have_inext) {
    res <- md_cache(CFG, "17_inext",
                    local_partition_inext(CFG, sp, grp_list, qs, t_alpha, t_gamma,
                                          nboot, conf, notes))
    part_rows <- c(part_rows, res$rows)
    cov_rows  <- c(cov_rows, res$cov)
    notes <- c(notes, res$notes)
  }

  # ---- 系统发育引擎（iNEXT.3D 同级别昂贵，同样缓存）----
  if (want_3d && have_3d) {
    tr <- md_core_tree(CFG)
    if (is.null(tr)) {
      notes <- c(notes, "phylo$enabled=FALSE 或无树：iNEXT.3D 系统发育分区跳过")
    } else {
      res3 <- md_cache(CFG, "17_inext3d",
                       local_partition_inext3d(CFG, sp, tr, grp_list, qs, t_alpha, t_gamma,
                                               nboot, conf, notes))
      part_rows <- c(part_rows, res3$rows)
      cov_rows  <- c(cov_rows, res3$cov)
      notes <- c(notes, res3$notes)
    }
  }

  part_df <- if (length(part_rows)) {
    nc <- vapply(part_rows, length, integer(1))
    if (length(unique(nc)) > 1) {
      md_warn(CFG, "17：分区结果行宽不一致（", paste(unique(nc), collapse = "/"),
              " 列），按并集对齐后合并")
      allc <- unique(unlist(lapply(part_rows, names)))
      part_rows <- lapply(part_rows, function(r) {
        for (mm in setdiff(allc, names(r))) r[[mm]] <- NA_real_
        r[allc]
      })
    }
    do.call(rbind, lapply(part_rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  } else data.frame(Group = character(0), stringsAsFactors = FALSE)
  cov_df  <- if (length(cov_rows)) do.call(rbind, cov_rows) else
    data.frame(Group = character(0), stringsAsFactors = FALSE)
  if (nrow(part_df)) {
    part_df$Status <- ifelse(is.na(part_df$NA_reason), "OK", "NA")
  }

  # ---- 图 ----
  # 只出 γ 分解（= 原脚本 02_Fig02.R §E 的 fig2E，**单面板**）。
  # 原脚本**没有覆盖度面板**，所以这里不再画 B（iNEXT 的 SC）——
  # 覆盖率数字仍写在 Coverage 表里；画图函数 local_fig_coverage() 保留备用。
  plots <- list(A = local_fig_gamma(part_df, CFG))
  fig <- tryCatch(md_compose(CFG, plots, design = "A"), error = function(e) NULL)
  if (is.null(fig))
    fig <- md_compose(CFG, list(A = md_empty_plot(CFG, "partition not estimable")), design = "A")
  md_save_figure(CFG, fig, "17_gamma_partition", height = 110)

  # ---- 表 ----
  note_df <- data.frame(
    Item = c("engine", "q", "nboot", "conf", "t_alpha", "t_gamma", "t_gamma_source",
             "coverage_aware", "alpha_definition", "gamma_calibers", "hill_scale",
             "groups_defined", "Notes"),
    Value = c(paste(c(if (want_inext) "iNEXT", if (want_3d) "iNEXT.3D"), collapse = " + "),
              paste(qs, collapse = ", "), as.character(nboot), as.character(conf),
              as.character(t_alpha), as.character(t_gamma),
              if (t_auto) "自动（各组采样单元数中位数）" else "CFG$partition$t_gamma",
              as.character(isTRUE(part$coverage_aware %||% TRUE)),
              "α = 同一条 incidence 曲线在 t = t_alpha 的稀有化 Hill 数（默认 t=1 个采样单元）",
              "A = t_gamma 稀有化；B = iNEXT/iNEXT.3D AsyEst 渐近估计（附 SE）；C = 实际观测采样量",
              "Hill 数刻度：γ = α x β 严格成立，α 占比 = 100 / β_mult",
              paste0(length(grp_list), " 组（", paste(unique(vapply(grp_list, function(z) z$Level,
                                                                   character(1))), collapse = " / "), "）"),
              if (length(notes)) paste(unique(notes), collapse = " | ") else "无"),
    stringsAsFactors = FALSE)
  tp <- md_write_xlsx(CFG, "17_gamma_partition", list(
    Partition = part_df,
    Coverage = cov_df,
    Note = note_df),
    meta_extra = list(engine = engine, nboot = nboot, t_gamma = t_gamma,
                      t_gamma_source = if (t_auto) "auto" else "config"))
  if (is.null(tp)) return(md_result(module, "PARTIAL", "写表失败（缺 openxlsx？）"))

  # 列名是 Gamma_A_q<q>（q 为 Hill 数刻度），不是 Gamma_A —— 原先写死 Gamma_A 导致这里恒打印 0
  gamma_col <- grep("^Gamma_A_q[0-9.]+$", names(part_df), value = TRUE)[1]
  md_log(CFG, "INFO", "17 完成：", nrow(part_df), " 行分区结果（",
         if (is.na(gamma_col)) NA_integer_ else sum(!is.na(part_df[[gamma_col]])),
         " 行 γ 可估计，列 = ", gamma_col %||% "无", "），t_gamma = ", t_gamma)
  md_result(module, "OK", tables = basename(tp), figures = "figures/17_gamma_partition.pdf",
            notes = notes)
}

# =============================================================================
# 内部工具（local_ 前缀）
# =============================================================================

#' 构造 incidence_freq 向量：c(n_units, 每个物种被多少个采样单元检出)
#' @param sp 物种(行) x 采样单元(列) 的 0/1 矩阵
local_incidence <- function(sp, ids) {
  sub <- sp[, ids, drop = FALSE]
  f <- rowSums(sub)
  n_inc <- sum(f)                       # 物种-采样单元 出现次数合计（iNEXT 的 S.obs 口径）
  f <- f[f > 0]
  list(n = length(ids), freq = f, n_species = length(f), n_incidence = n_inc)
}

#' 该组的数据量是否满足 iNEXT incidence_freq 的最低要求。
#' iNEXT 在"没有任何物种的频率落在 2..(n-1)"时会抛 NA/NaN 参数错误，
#' 因此在调用前直接判定并给出可读原因。
local_inext_ok <- function(inc) {
  if (inc$n < 2) return("采样单元 < 2")
  if (inc$n_species < 2) return(paste0("只检出 ", inc$n_species, " 个物种"))
  f <- inc$freq
  if (inc$n < 3) return("采样单元 < 3")
  if (!any(f > 1)) return("所有物种都是单例（f = 1），iNEXT 无法估计")
  if (!any(f < inc$n)) return(paste0("所有物种都在全部 ", inc$n, " 个采样单元中出现，无稀有化空间"))
  NA_character_
}

#' 从 iNEXT 结果的 size_based 曲线里取 t 处的 Hill qD（含 LCL/UCL/SC）
local_curve_at <- function(sb, assemblage, q, t) {
  d <- sb[sb$Assemblage == assemblage & as.numeric(as.character(sb$Order.q)) == q, , drop = FALSE]
  if (!nrow(d)) return(list(qD = NA_real_, lo = NA_real_, hi = NA_real_, SC = NA_real_, t = NA_real_))
  tt <- as.numeric(d$t)
  i <- which.min(abs(tt - t))
  list(qD = as.numeric(d$qD[i]), lo = as.numeric(d$qD.LCL[i]), hi = as.numeric(d$qD.UCL[i]),
       SC = as.numeric(d$SC[i]), t = tt[i])
}

#' 选一个用于绘图/汇总的 Hill 阶：优先 q = 1，否则最小的可用 q
local_pick_q <- function(part_df) {
  cand <- as.numeric(sub("^Alpha_q", "", grep("^Alpha_q[0-9.]+$", names(part_df), value = TRUE)))
  cand <- sort(unique(cand[is.finite(cand)]))
  if (!length(cand)) return(NA_real_)
  if (1 %in% cand) 1 else cand[1]
}

#' 分类学（iNEXT）分区
local_partition_inext <- function(CFG, pa, grp_list, qs, t_alpha, t_gamma, nboot, conf, notes) {
  rows <- list(); cov <- list()
  allc <- local_partition_cols(qs)
  new_row <- function(Group, Level, Type, n_units, n_distinct, n_incidence,
                      NA_reason = NA_character_) {
    r <- stats::setNames(as.list(rep(NA_real_, length(allc))), allc)
    r$Group <- Group; r$Level <- Level; r$Type <- Type
    r$n_units <- n_units; r$n_species_distinct <- n_distinct
    r$n_species_incidence <- n_incidence; r$NA_reason <- NA_reason
    r
  }
  for (g in grp_list) {
    inc <- local_incidence(pa, g$ids)
    reason <- local_inext_ok(inc)
    if (!is.na(reason)) {
      rows[[length(rows) + 1L]] <- new_row(g$Group, g$Level, "Taxonomic",
                                           inc$n, inc$n_species, inc$n_incidence, reason)
      notes <- c(notes, paste0("iNEXT 跳过组 ", g$Group, "：", reason))
      next
    }
    md_seed(CFG, 1100 + nchar(g$Group) + inc$n)
    obj <- tryCatch(iNEXT::iNEXT(list(A = c(inc$n, inc$freq)), q = qs,
                                 datatype = "incidence_freq", nboot = nboot, conf = conf),
                    error = function(e) e)
    if (inherits(obj, "error")) {
      rows[[length(rows) + 1L]] <- new_row(g$Group, g$Level, "Taxonomic",
                                           inc$n, inc$n_species, inc$n_incidence,
                                           conditionMessage(obj))
      notes <- c(notes, paste0("iNEXT 失败组 ", g$Group, "：", conditionMessage(obj)))
      next
    }
    di <- as.data.frame(obj$DataInfo)
    sc <- if ("SC" %in% names(di)) di$SC[1] else NA_real_
    sb <- as.data.frame(obj$iNextEst$size_based)
    asy <- as.data.frame(obj$AsyEst)
    max_t <- max(as.numeric(sb$t), na.rm = TRUE)
    t_use <- if (inc$n >= t_gamma) t_gamma else max_t
    cov[[length(cov) + 1L]] <- data.frame(
      Group = g$Group, Level = g$Level, Type = "Taxonomic", n_units = inc$n,
      n_species_distinct = inc$n_species, n_species_incidence = inc$n_incidence, SC = sc,
      gamma_units_used = if (inc$n >= t_gamma) t_gamma else NA_real_,
      gamma_source = if (inc$n >= t_gamma) "t_gamma" else paste0("observed_max_t=", max_t),
      Extrapolated = inc$n < t_gamma, stringsAsFactors = FALSE)
    one <- new_row(g$Group, g$Level, "Taxonomic", inc$n, inc$n_species, inc$n_incidence)
    one$SC <- sc
    one$t_gamma_used <- t_use
    for (q in qs) {
      a <- local_curve_at(sb, "A", q, t_alpha)
      b <- local_curve_at(sb, "A", q, t_use)
      cc <- local_curve_at(sb, "A", q, max_t)
      ast <- local_asy_at(asy, q)
      one[[paste0("Alpha_q", q)]] <- a$qD
      one[[paste0("Alpha_q", q, "_LCL")]] <- a$lo
      one[[paste0("Alpha_q", q, "_UCL")]] <- a$hi
      one[[paste0("Gamma_A_q", q)]] <- b$qD
      one[[paste0("Gamma_A_q", q, "_LCL")]] <- b$lo
      one[[paste0("Gamma_A_q", q, "_UCL")]] <- b$hi
      one[[paste0("Gamma_B_q", q)]] <- ast$est
      one[[paste0("Gamma_B_q", q, "_SE")]] <- ast$se
      one[[paste0("Gamma_C_q", q)]] <- cc$qD
      one[[paste0("Beta_mult_q", q)]] <- ifelse(is.finite(b$qD) && is.finite(a$qD) && a$qD > 0,
                                                b$qD / a$qD, NA_real_)
      one[[paste0("Alpha_pct_q", q)]] <- ifelse(is.finite(b$qD) && is.finite(a$qD) && b$qD > 0,
                                                100 * a$qD / b$qD, NA_real_)
    }
    rows[[length(rows) + 1L]] <- one
  }
  list(rows = rows, cov = cov, notes = notes)
}

#' 从 AsyEst 里取某 q 的渐近估计（Diversity 列名随 q 变化）
local_asy_at <- function(asy, q) {
  dv <- as.character(asy$Diversity)
  pat <- if (q == 0) "Richness" else if (abs(q - 1) < 1e-8) "Shannon" else if (abs(q - 2) < 1e-8) "Simpson"
         else paste0("q = ", q)
  i <- grep(pat, dv, ignore.case = TRUE)[1]
  if (is.na(i)) return(list(est = NA_real_, se = NA_real_))
  list(est = as.numeric(asy$Estimator[i]), se = as.numeric(asy$s.e.[i]))
}

#' 分区表每一行的规范列（success / NA 两条路径共用，避免列数漂移）
local_partition_cols <- function(qs) {
  c("Group", "Level", "Type", "n_units", "n_species_distinct", "n_species_incidence",
    "SC", "t_gamma_used", "NA_reason",
    unlist(lapply(qs, function(q) c(
      paste0("Alpha_q", q), paste0("Alpha_q", q, "_LCL"), paste0("Alpha_q", q, "_UCL"),
      paste0("Gamma_A_q", q), paste0("Gamma_A_q", q, "_LCL"), paste0("Gamma_A_q", q, "_UCL"),
      paste0("Gamma_B_q", q), paste0("Gamma_B_q", q, "_SE"), paste0("Gamma_C_q", q),
      paste0("Beta_mult_q", q), paste0("Alpha_pct_q", q))), use.names = FALSE))
}

#' 系统发育（iNEXT.3D）分区。失败/数据不足一律记 NA_reason，不报错。
local_partition_inext3d <- function(CFG, pa, tr, grp_list, qs, t_alpha, t_gamma,
                                    nboot, conf, notes) {
  rows <- list(); cov <- list()
  # 树准备：去引号 -> 剪到当前 ASV -> 二分支化 -> midpoint 生根 -> 去 node.label
  prep <- local_phylo_ready(tr, rownames(pa), CFG)
  if (is.null(prep$tree)) {
    notes <- c(notes, paste0("iNEXT.3D 跳过：", prep$msg))
    return(list(rows = rows, cov = cov, notes = notes))
  }
  tree <- prep$tree
  common <- prep$common
  notes <- c(notes, paste0("iNEXT.3D 用树：", length(tree$tip.label), " tips（",
                           length(common), "/", ncol(pa), " ASV 匹配），生根方式 ",
                           prep$rooting))
  allc <- local_partition_cols(qs)
  new_row <- function(Group, Level, Type, n_units, n_distinct, n_incidence,
                      NA_reason = NA_character_) {
    r <- stats::setNames(as.list(rep(NA_real_, length(allc))), allc)
    r$Group <- Group; r$Level <- Level; r$Type <- Type
    r$n_units <- n_units; r$n_species_distinct <- n_distinct
    r$n_species_incidence <- n_incidence; r$NA_reason <- NA_reason
    r
  }
  for (g in grp_list) {
    inc <- local_incidence(pa, g$ids)
    # iNEXT.3D 的 incidence_raw 要求 data 是 **物种(行) x 采样单元(列)** 的矩阵，
    # 行名 = 物种、列名 = 采样单元（iNEXT.3D::check.tree 用 rownames(data) 匹配 tip）。
    mtx <- pa[tree$tip.label, g$ids, drop = FALSE]
    mtx <- mtx[rowSums(mtx) > 0, , drop = FALSE]
    n_species <- nrow(mtx)
    if (n_species < 5) {
      reason <- paste0("观测谱系数 < 5（", n_species, "），iNEXT.3D 硬性下限")
      rows[[length(rows) + 1L]] <- new_row(g$Group, g$Level, "Phylogenetic",
                                           inc$n, n_species, inc$n_incidence, reason)
      notes <- c(notes, paste0("iNEXT.3D 跳过组 ", g$Group, "：", reason))
      next
    }
    md_seed(CFG, 1300 + nchar(g$Group) + inc$n)
    obj <- tryCatch(iNEXT.3D::iNEXT3D(list(A = mtx), diversity = "PD", q = qs,
                                      datatype = "incidence_raw", PDtree = tree,
                                      PDtype = "meanPD", nboot = nboot,
                                      endpoint = max(t_gamma, inc$n)),
                    error = function(e) e)
    if (inherits(obj, "error")) {
      rows[[length(rows) + 1L]] <- new_row(g$Group, g$Level, "Phylogenetic",
                                           inc$n, n_species, inc$n_incidence,
                                           conditionMessage(obj))
      notes <- c(notes, paste0("iNEXT.3D 失败组 ", g$Group, "：", conditionMessage(obj)))
      next
    }
    sb  <- as.data.frame(obj$PDiNextEst$size_based)
    asy <- as.data.frame(obj$PDAsyEst)
    max_t <- max(as.numeric(sb$mT), na.rm = TRUE)
    t_use <- if (inc$n >= t_gamma) t_gamma else max_t
    sc <- if ("SC" %in% names(sb)) sb$SC[which.min(abs(as.numeric(sb$mT) - t_use))] else NA_real_
    cov[[length(cov) + 1L]] <- data.frame(
      Group = g$Group, Level = g$Level, Type = "Phylogenetic", n_units = inc$n,
      n_species_distinct = n_species, n_species_incidence = inc$n_incidence, SC = sc,
      gamma_units_used = if (inc$n >= t_gamma) t_gamma else NA_real_,
      gamma_source = if (inc$n >= t_gamma) "t_gamma" else paste0("observed_max_t=", max_t),
      Extrapolated = inc$n < t_gamma, stringsAsFactors = FALSE)
    one <- new_row(g$Group, g$Level, "Phylogenetic", inc$n, n_species, inc$n_incidence)
    one$SC <- sc
    one$t_gamma_used <- t_use
    for (q in qs) {
      a  <- local_pd_curve_at(sb, q, t_alpha)
      b  <- local_pd_curve_at(sb, q, t_use)
      cc <- local_pd_curve_at(sb, q, max_t)
      ast <- local_pd_asy(asy, q)
      one[[paste0("Alpha_q", q)]] <- a$qD
      one[[paste0("Alpha_q", q, "_LCL")]] <- a$lo
      one[[paste0("Alpha_q", q, "_UCL")]] <- a$hi
      one[[paste0("Gamma_A_q", q)]] <- b$qD
      one[[paste0("Gamma_A_q", q, "_LCL")]] <- b$lo
      one[[paste0("Gamma_A_q", q, "_UCL")]] <- b$hi
      one[[paste0("Gamma_B_q", q)]] <- ast$est
      one[[paste0("Gamma_B_q", q, "_SE")]] <- ast$se
      one[[paste0("Gamma_C_q", q)]] <- cc$qD
      one[[paste0("Beta_mult_q", q)]] <- ifelse(is.finite(b$qD) && is.finite(a$qD) && a$qD > 0,
                                                b$qD / a$qD, NA_real_)
      one[[paste0("Alpha_pct_q", q)]] <- ifelse(is.finite(b$qD) && is.finite(a$qD) && b$qD > 0,
                                                100 * a$qD / b$qD, NA_real_)
    }
    rows[[length(rows) + 1L]] <- one
  }
  list(rows = rows, cov = cov, notes = notes)
}

#' 树准备（与 20_phylo_prep 同一策略，但不依赖其缓存，模块可独立跑）
local_phylo_ready <- function(tr, asvs, CFG) {
  if (!requireNamespace("ape", quietly = TRUE)) return(list(tree = NULL, msg = "缺 ape 包"))
  tr$tip.label <- gsub('^["\']|["\']$', "", trimws(tr$tip.label))
  common <- intersect(tr$tip.label, asvs)
  if (length(common) < 5)
    return(list(tree = NULL, msg = paste0("树与 ASV 交集不足 5（", length(common), "）")))
  t2 <- tr
  t2$node.label <- NULL
  t2 <- tryCatch(ape::keep.tip(t2, common), error = function(e) NULL)
  if (is.null(t2)) return(list(tree = NULL, msg = "keep.tip 失败"))
  rooting <- CFG$phylo$root_method %||% "midpoint"
  t2 <- tryCatch(ape::multi2di(t2), error = function(e) t2)
  t2 <- tryCatch(local_root_tree(t2, rooting, CFG), error = function(e) t2)
  if (any(is.na(t2$edge.length))) t2$edge.length[is.na(t2$edge.length)] <- 0
  if (any(t2$edge.length < 0)) t2$edge.length[t2$edge.length < 0] <- 0
  if (!ape::is.rooted(t2)) {
    t2 <- tryCatch(local_root_tree(t2, "midpoint", CFG), error = function(e) t2)
  }
  list(tree = t2, common = common, rooting = rooting,
       msg = "")
}

#' 按 root_method 生根：midpoint（phangorn）| outgroup（最长 tip 支）| none
local_root_tree <- function(tr, method = "midpoint", CFG = NULL) {
  method <- match.arg(method, c("midpoint", "outgroup", "none"))
  if (method == "none") return(tr)
  if (method == "midpoint") {
    if (requireNamespace("phangorn", quietly = TRUE))
      return(phangorn::midpoint(tr))
    method <- "outgroup"
  }
  if (ape::is.rooted(tr) && method == "outgroup") {
    # 已生根则不重复操作
    return(tr)
  }
  tip_edges <- which(tr$edge[, 2] <= length(tr$tip.label))
  if (!length(tip_edges)) return(tr)
  out_group <- tr$edge[tip_edges, 2][which.max(tr$edge.length[tip_edges])]
  ape::root(tr, outgroup = out_group, resolve.root = TRUE)
}

#' iNEXT.3D size_based 曲线取值（列名 qPD / mT）
local_pd_curve_at <- function(sb, q, t) {
  d <- sb[as.numeric(sb$Order.q) == q, , drop = FALSE]
  if (!nrow(d)) return(list(qD = NA_real_, lo = NA_real_, hi = NA_real_, SC = NA_real_, t = NA_real_))
  tt <- as.numeric(d$mT)
  i <- which.min(abs(tt - t))
  list(qD = as.numeric(d$qPD[i]), lo = as.numeric(d$qPD.LCL[i]),
       hi = as.numeric(d$qPD.UCL[i]), SC = as.numeric(d$SC[i]), t = tt[i])
}

local_pd_asy <- function(asy, q) {
  if (!nrow(asy)) return(list(est = NA_real_, se = NA_real_))
  lab <- as.character(asy[[1]])
  i <- grep(paste0("q = ", q, " PD"), lab, fixed = TRUE)[1]
  if (is.na(i)) i <- which.min(abs(seq_len(nrow(asy)) - q - 1))
  est_col <- grep("PD_asy", names(asy), value = TRUE)[1]
  se_col  <- grep("s\\.e\\.", names(asy), value = TRUE)[1]
  list(est = if (!is.na(est_col)) as.numeric(asy[[est_col]][i]) else NA_real_,
       se  = if (!is.na(se_col)) as.numeric(asy[[se_col]][i]) else NA_real_)
}

# γ 分解的组分配色与图例文字 —— 照抄 02_Fig02.R §E 的 cols_ep4 / labs_ep4：
#   分类学整族浅色、系统发育整族深色；组分色相不变（turnover = 蓝系、alpha = 鲑鱼系）。
# 说明：R/shared/12_plot_theme.R 的 MD_BUILTIN_PAL 没有这一族颜色，
#       md_pal() 会落到 qualitative 兜底色、与原图不一致，故按原脚本照抄这 4 个 hex
#       （见最终报告里的 shared 缺口）。
local_gamma_cols <- c(Taxonomic_turnover    = "#AECFE6",
                      Taxonomic_alpha       = "#F7CCC5",
                      Phylogenetic_turnover = "#5B93C4",
                      Phylogenetic_alpha    = "#E0877C")
local_gamma_labels <- c(Taxonomic_turnover    = "Taxonomic turnover (gamma - alpha)",
                        Taxonomic_alpha       = "Taxonomic alpha",
                        Phylogenetic_turnover = "Phylogenetic turnover (gamma - alpha)",
                        Phylogenetic_alpha    = "Phylogenetic alpha")

#' 面板 A：γ 分解（Hill 刻度）—— 画法照抄 02_Fig02.R §E（fig2E）：
#'   · 每个分组两根柱：左 = 分类学、右 = 系统发育（本包 engine=iNEXT 时只有分类学，
#'     按原脚本同一套水平退化为一根柱）
#'   · 柱内堆叠 α 段（0..α）+ 周转段（α..γ），用 geom_rect 手工累计 xmin/xmax
#'     （原脚本注释：geom_col 的 position_dodge 无法同时堆叠 + 错位）
#'   · 组内两根柱中心错开 ±BAR_OFF、半宽 BAR_HW；白色描边 linewidth 0.15
#'   · 双 x 轴：主轴 = 有效物种数，副轴 = meanPD（系统发育柱按 K 放大，
#'     使系统发育 γmax 与分类学 γmax 物理等长）；分组走 scale_y_reverse
local_fig_gamma <- function(part_df, CFG) {
  if (is.null(part_df) || !nrow(part_df)) return(md_empty_plot(CFG, "no partition result"))
  qsel <- local_pick_q(part_df)
  if (is.na(qsel)) return(md_empty_plot(CFG, "no Hill q column"))
  a_col <- paste0("Alpha_q", qsel); g_col <- paste0("Gamma_A_q", qsel)
  if (!a_col %in% names(part_df) || !g_col %in% names(part_df))
    return(md_empty_plot(CFG, "no Hill q column"))
  d <- part_df[is.finite(part_df[[a_col]]) & is.finite(part_df[[g_col]]), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(CFG, "gamma not estimable in any group"))
  d$Alpha <- as.numeric(d[[a_col]])
  d$Gamma <- as.numeric(d[[g_col]])
  d$Type  <- factor(as.character(d$Type), levels = c("Taxonomic", "Phylogenetic"))
  d <- d[!is.na(d$Type), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(CFG, "gamma not estimable in any group"))
  have_phy <- all(c("Taxonomic", "Phylogenetic") %in% as.character(d$Type))

  # 双轴缩放 K：分类学 γmax / 系统发育 γmax（原脚本 K_EP）
  K <- 1
  if (have_phy) {
    gt <- max(d$Gamma[d$Type == "Taxonomic"], na.rm = TRUE)
    gp <- max(d$Gamma[d$Type == "Phylogenetic"], na.rm = TRUE)
    if (is.finite(gt) && is.finite(gp) && gp > 0) K <- gt / gp
  }
  d$Alpha_x <- ifelse(d$Type == "Phylogenetic", d$Alpha * K, d$Alpha)
  d$Gamma_x <- ifelse(d$Type == "Phylogenetic", d$Gamma * K, d$Gamma)

  # 分组顺序 = part_df 的行序（= CFG 定义的分组顺序，原脚本用 both_levels 定序）
  d$Lab <- factor(as.character(d$Group), levels = unique(as.character(d$Group)))
  grp_idx <- data.frame(Group = levels(d$Lab), idx = seq_along(levels(d$Lab)),
                        stringsAsFactors = FALSE)
  d <- merge(d, grp_idx, by.x = "Lab", by.y = "Group", all.x = TRUE, sort = FALSE)
  BAR_OFF <- 0.19   # 照抄 02_Fig02.R 顶部常量
  BAR_HW  <- 0.17
  d$xc <- d$idx + ifelse(d$Type == "Phylogenetic", BAR_OFF, -BAR_OFF)

  segs <- rbind(
    data.frame(xc = d$xc, Type = d$Type, xmin = 0, xmax = d$Alpha_x,
               Comp = "alpha", stringsAsFactors = FALSE),
    data.frame(xc = d$xc, Type = d$Type, xmin = d$Alpha_x,
               xmax = d$Alpha_x + pmax(d$Gamma_x - d$Alpha_x, 0),
               Comp = "turnover", stringsAsFactors = FALSE))
  segs$FC <- factor(paste(segs$Type, segs$Comp, sep = "_"),
                    levels = c("Taxonomic_turnover", "Taxonomic_alpha",
                               "Phylogenetic_turnover", "Phylogenetic_alpha"))

  brks <- pretty(c(0, max(segs$xmax, na.rm = TRUE) * 1.02), n = 5)
  XMAX <- max(brks)
  # 注意：ggplot2 4.x 里 sec.axis = NULL 会报错，必须只在有系统发育柱时才带上副轴
  x_scale <- if (have_phy && is.finite(K) && K > 0) {
    ggplot2::scale_x_continuous(
      limits = c(0, XMAX), expand = ggplot2::expansion(mult = c(0, 0)), breaks = brks,
      sec.axis = ggplot2::sec_axis(~ . / K, name = "phylogenetic diversity (meanPD)",
                                   breaks = pretty(c(0, XMAX / K), n = 4)))
  } else {
    ggplot2::scale_x_continuous(limits = c(0, XMAX),
                                expand = ggplot2::expansion(mult = c(0, 0)), breaks = brks)
  }
  n_leg <- length(intersect(names(local_gamma_labels), as.character(unique(segs$FC))))
  g <- ggplot2::ggplot(segs) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$xc - BAR_HW, ymax = .data$xc + BAR_HW,
                                    fill = .data$FC),
                       colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = local_gamma_cols,
                               breaks = names(local_gamma_labels),
                               labels = local_gamma_labels, name = NULL,
                               guide = ggplot2::guide_legend(nrow = if (n_leg > 2) 2 else 1,
                                                             byrow = TRUE)) +
    x_scale +
    ggplot2::scale_y_reverse(breaks = grp_idx$idx, labels = as.character(grp_idx$Group),
                             expand = ggplot2::expansion(add = 0.85)) +
    ggplot2::labs(x = "effective number of species (taxonomic)", y = NULL, tag = "A") +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical",
                   panel.grid.major.y = ggplot2::element_blank())
  g
}

#' 面板 B：样本覆盖度 SC（每个组一个点/柱），低 SC 明确提示外推风险。
#' 原脚本没有这个面板（覆盖度只在 TableS12 里报），故配色走 md_pal()、几何走 MD_* 档位。
local_fig_coverage <- function(part_df, cov_df, CFG) {
  d <- cov_df
  if (is.null(d) || !nrow(d) || !"SC" %in% names(d)) d <- part_df
  if (is.null(d) || !nrow(d) || !"SC" %in% names(d))
    return(md_empty_plot(CFG, "coverage not available"))
  d <- d[!is.na(d$SC), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(CFG, "coverage not available"))
  d$Lab <- paste0(d$Type, " | ", d$Level, " | ", d$Group)
  d$Lab <- factor(d$Lab, levels = d$Lab[order(d$SC)])
  gp <- md_pal(CFG, sort(unique(as.character(d$Type))))
  g <- ggplot2::ggplot(d, ggplot2::aes(y = .data$Lab, x = .data$SC)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$SC, y = .data$Lab, yend = .data$Lab),
                          colour = "grey80", linewidth = MD_LW$line) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Type), size = MD_PT_SIZE$jitter) +
    ggplot2::scale_colour_manual(values = gp, name = NULL) +
    ggplot2::geom_vline(xintercept = c(0.9, 0.95), linetype = "dashed",
                        linewidth = MD_LW$ref, colour = "grey55") +
    ggplot2::scale_x_continuous(limits = c(0, 1.02),
                                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::labs(x = "Sample coverage (SC)", y = NULL, tag = "B", colour = NULL,
                  subtitle = "AsyEst extrapolation grows sharply below SC ~ 0.95") +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical",
                   plot.subtitle = ggplot2::element_text(size = MD_PT$ann * 0.85, colour = "grey30"),
                   axis.text.y = ggplot2::element_text(size = MD_PT$tick * 0.8))
  g
}
