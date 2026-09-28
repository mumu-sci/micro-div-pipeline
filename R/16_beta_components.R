# =============================================================================
# R/16_beta_components.R ---- β 多样性的周转 / 嵌套分解（betapart）
# -----------------------------------------------------------------------------
# 依赖：13_distances（cache/13_dist.rds，用于 β 分类 vs 系统发育的 Mantel 对比）
#
# 口径（与 nifH 02_Fig02.R §D 同源）：
#   分类学 β：betapart::beta.pair(index.family = "sorensen") 在**抽平后 P/A 矩阵**上算，
#             β_total = β_turnover(sim) + β_nestedness(sne) 在配对级严格可加。
#   丰度型敏感性（可选）：betapart::beta.pair.abund(index.family = "bray")
#             的 balanced / gradient 分量；版本差异（balanced/gradient vs bal/gra）
#             用 local_ 取列函数兜住。
#   系统发育 β：本模块**不自己算** phylo Sorensen，而是用 13 的 wunifrac/unifrac
#             距离与分类学距离做 Mantel 对比（避免与 22_phylo_beta 口径漂移）。
#   分组聚合：先取配对级（下三角），再按因子水平对（无序水平对）求均值；
#             同时保留配对级长表供 detail=TRUE 时导出。
#
# 输出：figures/16_beta_components.pdf（一张 patchwork）
#       results/16_beta_components.xlsx（Turnover_nestedness / By_levelpair / Abundance /
#                                        Mantel / Note）
# =============================================================================

md_beta_components <- function(CFG) {
  module <- "16"
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少 vegan/ape/picante/betapart"))
  if (!md_precheck(CFG, module)) return(md_result(module, "SKIP", "前置条件不满足"))

  m    <- md_core_asv(CFG, "distances")
  meta <- md_core_meta(CFG)
  if (is.null(m) || !nrow(m) || ncol(m) < 2) {
    md_skip(CFG, module, "距离口径矩阵为空或样本不足")
    return(md_result(module, "SKIP", "矩阵不可用"))
  }
  # ASV x 样本 -> 样本 x ASV 的 0/1 矩阵（betapart 要 sites x species）
  pa <- t((m > 0) * 1L)
  n_singleton <- sum(colSums(pa) == 1)
  notes <- c(paste0("P/A 矩阵：", nrow(pa), " 样本 x ", ncol(pa), " ASV；其中 ",
                    n_singleton, " 个 ASV 只出现在 1 个样本（占 ",
                    round(100 * n_singleton / max(1, ncol(pa)), 1), "%）"))

  # ---- 1. 分类学 Sorensen 分解 ----
  bp <- tryCatch(betapart::beta.pair(pa, index.family = "sorensen"), error = function(e) e)
  if (inherits(bp, "error")) {
    md_skip(CFG, module, paste0("betapart::beta.pair 失败：", conditionMessage(bp)))
    return(md_result(module, "SKIP", "beta.pair 失败"))
  }
  b_tot <- as.matrix(bp$beta.sor)
  b_tur <- as.matrix(bp$beta.sim)
  b_nes <- as.matrix(bp$beta.sne)
  ids <- rownames(b_tot)

  pair_long <- local_mat_pairs(b_tot, b_tur, b_nes, ids)
  pair_long <- local_attach_pair_factors(pair_long, meta, CFG)

  # ---- 2. 丰度型敏感性（可选，失败即降级）----
  ab_long  <- NULL
  ab_note  <- "未计算"
  if (isTRUE(CFG$beta_components$abundance %||% TRUE)) {
    ab <- tryCatch(betapart::beta.pair.abund(pa, index.family = "bray"), error = function(e) e)
    if (inherits(ab, "error")) {
      ab_note <- paste0("beta.pair.abund 失败：", conditionMessage(ab))
      notes <- c(notes, ab_note)
    } else {
      b_bal <- local_first_matrix(ab, c("beta.bray.bal", "beta.bal"))
      b_gra <- local_first_matrix(ab, c("beta.bray.gra", "beta.gra"))
      b_abt <- local_first_matrix(ab, c("beta.bray", "beta.bray.total"))
      if (!is.null(b_bal) && !is.null(b_gra)) {
        ab_long <- local_mat_pairs_multi(ids, list(Balanced = b_bal, Gradient = b_gra,
                                                   Total = b_abt))
        ab_long <- local_attach_pair_factors(ab_long, meta, CFG)
        ab_note <- "betapart::beta.pair.abund(index.family = 'bray')"
      } else {
        ab_note <- "beta.pair.abund 返回结构与预期不符，已跳过"
        notes <- c(notes, ab_note)
      }
    }
  }

  # ---- 3. 按因子水平对聚合 ----
  by_lp <- local_aggregate_pairwise(pair_long, CFG)
  by_ab <- if (!is.null(ab_long)) local_aggregate_pairwise(ab_long, CFG, value_cols = c("Balanced", "Gradient", "Total")) else NULL

  # ---- 4. 分类 vs 系统发育 β 的 Mantel 对比 ----
  mantel_df <- local_beta_mantel(CFG, b_tot, b_tur)

  # ---- 5. 图 ----
  # 只画 **nifH Fig.2D** 对应的内容：A 全配对池化分解、B 按因子水平对的分解。
  # 原实现还有 C（丰度型敏感性）与 D（分类 vs 系统发育 Mantel 散点）两格，
  # 但 **Fig.2D 里没有这两格** —— 按用户「不得画其他图」的要求从图面移除；
  # 两者的**分析照常计算**并写进结果表（abundance / Mantel sheet），只是不再出图。
  # 相应的画图函数 local_fig_abund() / local_fig_mantel() 保留在文件末尾备用（默认不调用）。
  # 只保留**逐因子水平对**的分解（原 Fig.2D 的呈现口径）。
  # 池化单柱（A）与丰度敏感性（C）、Mantel 散点（D）都不是 Fig.2D 的内容：
  # A 是一根柱子、C/D 原图没有 —— 全部从图面移除，**分析照常计算并写进结果表**。
  # 备用画图函数 local_fig_pooled() / local_fig_abund() / local_fig_mantel() 保留在文件末尾。
  plots <- list(A = local_fig_levelpair(by_lp, CFG))
  fig <- tryCatch(md_compose(CFG, plots, design = "A"), error = function(e) NULL)
  if (is.null(fig))
    fig <- md_compose(CFG, list(A = md_empty_plot(CFG, "beta partition not estimable")), design = "A")
  md_save_figure(CFG, fig, "16_beta_components", height = 105)

  # ---- 6. 表 ----
  summ <- local_pooled_summary(pair_long)
  note_df <- data.frame(
    Item = c("taxonomic_partition", "abundance_partition", "phylo_comparison",
             "aggregation", "additivity", "singletons", "Notes"),
    Value = c(
      "betapart::beta.pair(index.family = 'sorensen')，在抽平后 P/A 矩阵上（与 nifH 02_Fig02.R 同源）",
      ab_note,
      "本模块只用 13 的 wunifrac/unifrac 与分类学距离做 Mantel 对比；系统发育 Sorensen 的分解见模块 22",
      "配对级（下三角）先算，再按因子水平对（无序）求均值；水平对标签按因子水平顺序拼接",
      "β_total = β_turnover + β_nestedness 在配对级严格可加；均值层面亦成立",
      paste0(n_singleton, " / ", ncol(pa), " 个 ASV 为单例（示范数据极端单例化，嵌套分量接近 0 是真实特征）"),
      if (length(notes)) paste(unique(notes), collapse = " | ") else "无"),
    stringsAsFactors = FALSE)

  sheets <- list(
    Pair_summary = summ,
    By_levelpair = by_lp,
    Abundance = by_ab,
    Mantel = mantel_df,
    Note = note_df)
  if (isTRUE(CFG$table$detail)) sheets$Pairs <- pair_long
  tp <- md_write_xlsx(CFG, "16_beta_components", sheets,
                      meta_extra = list(partition = "sorensen",
                                        n_pairs = nrow(pair_long)))
  if (is.null(tp)) return(md_result(module, "PARTIAL", "写表失败（缺 openxlsx？）"))

  md_log(CFG, "INFO", "16 完成：配对 ", nrow(pair_long), " 对（", nrow(pa), " 样本）",
         " | 周转均值 ", md_fmt_num(mean(pair_long$Turnover, na.rm = TRUE), 4),
         " / 嵌套均值 ", md_fmt_num(mean(pair_long$Nestedness, na.rm = TRUE), 4))
  md_result(module, "OK", tables = basename(tp), figures = "figures/16_beta_components.pdf",
            notes = notes)
}

# =============================================================================
# 内部工具（local_ 前缀）
# =============================================================================

# 周转 / 嵌套的组分配色与图例文字 —— 照抄 02_Fig02.R §D：
#   退化支路（无系统发育版时）写死 c(Turnover = "#5B93C4", Nestedness = "#E0877C")，
#   图例文字同 §D 的 lab_beta4（"Taxonomic turnover" / "Taxonomic nestedness"）。
# 说明：R/shared/12_plot_theme.R 的 MD_BUILTIN_PAL（zone/veg/assembly/quadrant）
#       没有这一对颜色，md_pal() 会落到 qualitative 兜底色、与原图不一致，
#       故按原脚本照抄这两个 hex（见最终报告里的 shared 缺口）。
local_beta_cols <- c(Turnover = "#5B93C4", Nestedness = "#E0877C")
local_beta_labels <- c(Turnover = "Taxonomic turnover", Nestedness = "Taxonomic nestedness")

#' 三个对称矩阵 -> 下三角配对长表（保持 sample1 < sample2 的行序）
local_mat_pairs <- function(m_tot, m_tur, m_nes, ids) {
  out <- local_mat_pairs_multi(ids, list(Total = m_tot, Turnover = m_tur, Nestedness = m_nes))
  out
}

#' 任意个对称矩阵 -> 配对长表
local_mat_pairs_multi <- function(ids, mats) {
  k <- which(lower.tri(matrix(0, length(ids), length(ids))), arr.ind = TRUE)
  out <- data.frame(sample1 = ids[k[, 1]], sample2 = ids[k[, 2]], stringsAsFactors = FALSE)
  for (nm in names(mats)) {
    mm <- mats[[nm]]
    if (is.null(mm)) { out[[nm]] <- NA_real_; next }
    out[[nm]] <- mm[k]
  }
  out
}

#' betapart 版本差异兜底：按候选列名取第一个命中的矩阵
local_first_matrix <- function(x, cand) {
  for (nm in cand) if (!is.null(x[[nm]])) return(as.matrix(x[[nm]]))
  NULL
}

#' 给配对长表补因子列：f1_1/f1_2/... + 每个因子的"水平对"列
local_attach_pair_factors <- function(df, meta, CFG) {
  smp <- CFG$data$sample_col
  facs <- CFG$design$factors %||% character(0)
  facs <- facs[facs %in% names(meta)]
  idx1 <- match(df$sample1, meta[[smp]])
  idx2 <- match(df$sample2, meta[[smp]])
  for (f in facs) {
    a <- as.character(meta[[f]][idx1]); b <- as.character(meta[[f]][idx2])
    lv <- CFG$design$ordered_levels[[f]] %||% sort(unique(c(a, b)))
    a <- factor(a, levels = lv); b <- factor(b, levels = lv)
    df[[paste0(f, "_1")]] <- as.character(a)
    df[[paste0(f, "_2")]] <- as.character(b)
    df[[paste0("Same_", f)]] <- as.character(a) == as.character(b)
    df[[paste0("Pair_", f)]] <- local_level_pair(a, b, lv)
  }
  if (length(facs) >= 2) {
    df$Pair_Interaction <- local_level_pair2(df[[paste0("Pair_", facs[1])]],
                                             df[[paste0("Pair_", facs[2])]])
  }
  df
}

#' 无序水平对标签（"A - B"，按 level 顺序）
local_level_pair <- function(a, b, lv) {
  a <- as.character(a); b <- as.character(b)
  out <- mapply(function(x, y) {
    if (is.na(x) || is.na(y)) return(NA_character_)
    if (x == y) return(x)
    v <- c(x, y)
    if (!is.null(lv)) v <- lv[lv %in% v]
    paste(sort(v), collapse = " - ")
  }, a, b, USE.NAMES = FALSE)
  out
}

local_level_pair2 <- function(pa, pb) {
  mapply(function(x, y) paste(sort(c(x, y)), collapse = " | "), pa, pb, USE.NAMES = FALSE)
}

#' 按因子（与因子交互）聚合配对级均值。
#' 用 tapply 逐列算，避免 stats::aggregate 在分组列与数据列同名/空水平时的静默失败。
local_aggregate_pairwise <- function(df, CFG,
                                     value_cols = c("Turnover", "Nestedness", "Total")) {
  facs <- CFG$design$factors %||% character(0)
  facs <- facs[paste0("Pair_", facs) %in% names(df)]
  cols <- intersect(value_cols, names(df))
  if (!length(facs) || !length(cols)) return(NULL)

  agg_by <- function(d, key, factor_label) {
    keep <- !is.na(key)
    d <- d[keep, , drop = FALSE]; key <- key[keep]
    if (!nrow(d) || length(unique(key)) == 0) return(NULL)
    key <- as.character(key)
    out <- data.frame(Factor = factor_label, Level_pair = unique(key), stringsAsFactors = FALSE)
    for (cc in cols) {
      mu <- tapply(d[[cc]], key, function(x) mean(x, na.rm = TRUE))
      out[[cc]] <- unname(mu[out$Level_pair])
    }
    out$n_pairs <- as.integer(table(key)[out$Level_pair])
    out[, c("Factor", "Level_pair", "n_pairs", cols), drop = FALSE]
  }

  rows <- list()
  for (f in facs) {
    col <- paste0("Pair_", f)
    if (!col %in% names(df)) next
    a <- agg_by(df, df[[col]], f)
    if (!is.null(a)) rows[[length(rows) + 1L]] <- a
  }
  if (length(facs) >= 2 && "Pair_Interaction" %in% names(df)) {
    a <- agg_by(df, df$Pair_Interaction, paste(facs[1:2], collapse = " x "))
    if (!is.null(a)) rows[[length(rows) + 1L]] <- a
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' 全配对池化后的三行摘要（Total / Turnover / Nestedness）
local_pooled_summary <- function(pair_long) {
  cols <- intersect(c("Total", "Turnover", "Nestedness"), names(pair_long))
  data.frame(
    Component = cols,
    mean = vapply(cols, function(cc) mean(pair_long[[cc]], na.rm = TRUE), numeric(1)),
    sd   = vapply(cols, function(cc) stats::sd(pair_long[[cc]], na.rm = TRUE), numeric(1)),
    median = vapply(cols, function(cc) stats::median(pair_long[[cc]], na.rm = TRUE), numeric(1)),
    n_pairs = nrow(pair_long),
    stringsAsFactors = FALSE)
}

#' 分类学 β vs 13 缓存里的系统发育距离做 Mantel
local_beta_mantel <- function(CFG, b_tot, b_tur) {
  rows <- list()
  cand <- intersect(CFG$distances$metrics %||% character(0),
                    c("wunifrac", "unifrac"))
  for (mt in cand) {
    d <- md_dist_get(CFG, mt)
    if (is.null(d)) next
    ids <- intersect(rownames(b_tot), attr(d, "Labels") %||% labels(d))
    if (length(ids) < 4) next
    tm <- stats::as.dist(b_tot[ids, ids, drop = FALSE])
    dm <- md_subset_dist(d, ids)
    for (lab in c("Sorensen total", "Sorensen turnover")) {
      tax_d <- if (lab == "Sorensen total") tm else stats::as.dist(b_tur[ids, ids, drop = FALSE])
      # **必须设种子**：vegan::mantel 的置换 p 是蒙特卡洛量。
      # 不设种子时它取决于进入本模块时的会话 RNG（是否先跑过 15/30/60），
      # 同一份数据会在两次运行间给出不同的 P（实测 0.003 → 0.001）—— 不可复现。
      # 种子按「距离度量 × 分解类型」派生，保证可复现且各次检验互不相关。
      md_seed(CFG, 1600L + match(mt, cand) * 10L +
                match(lab, c("Sorensen total", "Sorensen turnover")))
      r <- tryCatch(vegan::mantel(tax_d, dm, method = "spearman", permutations = 999),
                    error = function(e) NULL)
      if (is.null(r)) next
      rows[[length(rows) + 1L]] <- data.frame(
        tax_metric = lab, phylo_metric = mt, method = "Spearman",
        statistic = unname(r$statistic), p = r$signif, permutations = 999,
        n_samples = length(ids), stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(Note = "13 缓存中没有系统发育距离（unifrac/wunifrac），跳过", stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

#' 面板 A：全配对池化分解（堆叠柱，标出各分量占比）
#' 画法对齐 02_Fig02.R §D：geom_col(position = "stack", width = 0.54) + 白色描边 0.15，
#' 分量配色与图例文字同 §D。
local_fig_pooled <- function(pair_long) {
  cols <- intersect(c("Total", "Turnover", "Nestedness"), names(pair_long))
  if (!length(cols)) return(md_empty_plot(NULL, "no partition"))
  mu <- vapply(cols, function(cc) mean(pair_long[[cc]], na.rm = TRUE), numeric(1))
  tot <- mu["Total"]
  d <- data.frame(
    Component = factor(c("Turnover", "Nestedness"), levels = c("Nestedness", "Turnover")),
    Value = c(mu["Turnover"], mu["Nestedness"]),
    stringsAsFactors = FALSE)
  d$pct <- d$Value / sum(d$Value)
  g <- ggplot2::ggplot(d, ggplot2::aes(x = "All pairs", y = .data$Value, fill = .data$Component)) +
    ggplot2::geom_col(width = 0.54, colour = "white", linewidth = 0.15) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(md_fmt_num(.data$Value, 3), "\n(",
                                                   md_fmt_pct(.data$pct), ")")),
                       position = ggplot2::position_stack(vjust = 0.5),
                       size = pt2gg(MD_PT$ann * 0.9), colour = "grey15") +
    ggplot2::scale_fill_manual(values = local_beta_cols, labels = local_beta_labels, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, max(1, sum(d$Value)) * 1.06),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "Sorensen dissimilarity", tag = "A",
                  subtitle = paste0("n = ", nrow(pair_long), " pairs; total = ",
                                    md_fmt_num(tot, 3))) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "none",
                   plot.subtitle = ggplot2::element_text(size = MD_PT$ann * 0.9, colour = "grey30"))
  g
}

#' 面板 B：按因子水平对的周转/嵌套均值（堆叠柱；轴太长时只保留最高/最低若干）
#' 画法对齐 02_Fig02.R §D 的退化支路：堆叠柱 + 白色描边 0.15 + 图例在下方，
#' 轴标题文字同原图（"Sorensen dissimilarity"，原图用非 ASCII 的 ø，本包按铁律写 ASCII）。
local_fig_levelpair <- function(by_lp, CFG) {
  if (is.null(by_lp) || !nrow(by_lp) || !"Turnover" %in% names(by_lp))
    return(md_empty_plot(CFG, "no level-pair summary"))
  d <- by_lp[!is.na(by_lp$Turnover), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(CFG, "no level-pair summary"))
  # 只画**第一个因子**（RegionPair = Zone 的区域对）—— 原脚本 02_Fig02.R §D 就是
  # 一条 Zone RegionPair 的 x 轴，没有按因子分面。植被型的水平对不再出图（数字留在表里）。
  fac_keep <- unique(as.character(d$Factor))[1]
  d <- d[as.character(d$Factor) == fac_keep, , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(CFG, "no level-pair summary"))
  # 标签用**短的** Level_pair（如 "Coastal - Nearshore"）；
  # 竖柱 + 45° 轴标签 + 底部两行图例 —— 对齐原脚本 02_Fig02.R §D（L558-572）。
  d$Lab <- as.character(d$Level_pair)
  d$Tot <- d$Turnover + ifelse(is.na(d$Nestedness), 0, d$Nestedness)
  # 只画信息量最大的前 14 个（避免 x 轴标签爆掉）
  d <- d[order(-d$Tot), , drop = FALSE]
  if (nrow(d) > 14) d <- d[seq_len(14), , drop = FALSE]
  d$Lab <- factor(d$Lab, levels = unique(d$Lab))
  d$Facet <- factor(as.character(d$Factor), levels = unique(as.character(d$Factor)))
  d2 <- rbind(
    data.frame(Lab = d$Lab, Facet = d$Facet, Component = "Turnover",
               Value = d$Turnover, stringsAsFactors = FALSE),
    data.frame(Lab = d$Lab, Facet = d$Facet, Component = "Nestedness",
               Value = d$Nestedness, stringsAsFactors = FALSE))
  d2$Component <- factor(d2$Component, levels = c("Nestedness", "Turnover"))
  g <- ggplot2::ggplot(d2, ggplot2::aes(x = .data$Lab, y = .data$Value, fill = .data$Component)) +
    ggplot2::geom_col(width = 0.54, colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = local_beta_cols, labels = local_beta_labels, name = NULL,
                               guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "Sorensen dissimilarity", tag = "A") +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = MD_PT_COMP$tick))
  # **不分面**：原脚本 02_Fig02.R §D 只有 Zone 的 RegionPair 一条 x 轴，
  # 没有按因子分面的 strip。多因子时只画第一个因子（RegionPair）的水平对。
  g
}

#' 面板 C：丰度型分量（balanced / gradient）按水平对。
#' 本面板是原脚本没有的敏感性面板（Fig. 2D 只画分类学 Sørensen），
#' 故配色不照抄任何原图，走 md_pal()；几何与 B 面板同款。
local_fig_abund <- function(by_ab, CFG) {
  if (is.null(by_ab) || !nrow(by_ab) || !all(c("Balanced", "Gradient") %in% names(by_ab)))
    return(md_empty_plot(CFG, "abundance partition not computed"))
  d <- by_ab[!is.na(by_ab$Balanced), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(CFG, "abundance partition not computed"))
  d$Lab <- paste0(d$Factor, ": ", d$Level_pair)
  d$Tot <- d$Balanced + ifelse(is.na(d$Gradient), 0, d$Gradient)
  d <- d[order(-d$Tot), , drop = FALSE]
  if (nrow(d) > 14) d <- d[seq_len(14), , drop = FALSE]
  d$Lab <- factor(d$Lab, levels = rev(d$Lab))
  d2 <- rbind(
    data.frame(Lab = d$Lab, Component = "Balanced", Value = d$Balanced, stringsAsFactors = FALSE),
    data.frame(Lab = d$Lab, Component = "Gradient", Value = d$Gradient, stringsAsFactors = FALSE))
  d2$Component <- factor(d2$Component, levels = c("Gradient", "Balanced"))
  gp <- md_pal(CFG, c("Balanced", "Gradient"))
  ggplot2::ggplot(d2, ggplot2::aes(y = .data$Lab, x = .data$Value, fill = .data$Component)) +
    ggplot2::geom_col(width = 0.54, colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = gp, name = NULL) +
    ggplot2::labs(x = "Bray dissimilarity", y = NULL, tag = "C",
                  subtitle = "beta.pair.abund (abundance-based sensitivity)") +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical",
                   plot.subtitle = ggplot2::element_text(size = MD_PT$ann * 0.85, colour = "grey30"),
                   axis.text.y = ggplot2::element_text(size = MD_PT$tick * 0.85))
}

#' 面板 D：分类学 vs 系统发育 β 的 Mantel 散点（抽样 3000 对以内）
#' 本面板也是原脚本没有的（原图 §D 不画散点），故散点走 MD_PT_SIZE / MD_REG 的通用档位。
local_fig_mantel <- function(CFG, b_tot, b_tur) {
  dc <- tryCatch(md_dist_cache(CFG), error = function(e) NULL)
  if (is.null(dc)) return(NULL)
  cand <- intersect(names(dc$dist), c("wunifrac", "unifrac"))
  if (!length(cand)) return(NULL)
  mt <- cand[1]
  d <- dc$dist[[mt]]
  ids <- intersect(rownames(b_tot), attr(d, "Labels") %||% labels(d))
  if (length(ids) < 4) return(NULL)
  pl <- local_mat_pairs(base::as.matrix(md_subset_dist(d, ids)), NULL, NULL, ids)
  pl$Phylo <- pl$Total; pl$Total <- NULL
  pl$Tax <- b_tot[cbind(match(pl$sample1, ids), match(pl$sample2, ids))]
  pl <- pl[stats::complete.cases(pl), , drop = FALSE]
  if (!nrow(pl)) return(NULL)
  set.seed(1)
  if (nrow(pl) > 3000) pl <- pl[sample.int(nrow(pl), 3000), , drop = FALSE]
  r <- tryCatch(stats::cor.test(pl$Tax, pl$Phylo, method = "spearman"),
                error = function(e) NULL)
  g <- ggplot2::ggplot(pl, ggplot2::aes(x = .data$Tax, y = .data$Phylo)) +
    ggplot2::geom_point(size = MD_PT_SIZE$scatter, alpha = MD_PT_SIZE$point_alpha,
                        colour = MD_REG$col) +
    ggplot2::labs(x = "Taxonomic Sorensen", y = paste0("Phylogenetic (", mt, ")"),
                  tag = "D") +
    theme_md(CFG = NULL)
  if (!is.null(r)) {
    g <- g + ggplot2::annotate("text", x = -Inf, y = Inf,
                               label = paste0("Spearman rho = ", md_fmt_num(unname(r$estimate), 3),
                                              ", P = ", md_fmt_p(r$p.value)),
                               hjust = -0.06, vjust = 1.6, size = pt2gg(MD_PT$ann * 0.85))
  }
  g
}
