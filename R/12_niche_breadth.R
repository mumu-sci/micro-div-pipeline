# =============================================================================
# R/12_niche_breadth.R ---- 模块 12：生态位宽度与生态策略组成（原图 Fig1c）
# -----------------------------------------------------------------------------
# 源：毕业论文 1_community_composition.R 第 415-900 行
#       （原 6_niche_rarity.R 的样本级四象限组成 + §7 稀有物种清单 / Levins' B）
#     关键原文行：679-708（fig07 双因素堆叠柱）、850-867（p_07sa 四象限散点）、
#                 869-879（p_07sb Levins B + gam）、887-891（两行拼图 200x150 mm）
#
# 做什么（**版式完全对齐原图 Fig1c_niche_rarity**）
#   A  样本生态策略组成堆叠柱，x = **前两个设计因子的交互组**（原 fig07）
#   B  丰度 x 出现率四象限散点（原 p_07sa；2026-09-26 从模块 10_composition 移入）
#   C  Levins B vs 出现率 + **gam 平滑 + 95% CI**（原 p_07sb）
#   表 = Niche_Breadth / Cell_summary / Sample_strategy / Threshold_sensitivity /
#        Sample_stats / Gradient_test / Note
#
# 口径
#   · 计算走共享库：md_niche_table() + md_sample_strategy()（不重写公式）
#   · 四象限 = 平均相对丰度(%) x 出现率(%) 双阈值（默认 0.1 x 10，原 classify_2x2 口径）
#   · 输入矩阵：md_core(CFG)$asv_raw（未抽平；原文用 ps_rare，CFG$niche$input = "rarefied"
#     时改用 $asv_rare）
#   · index = "both" 时表里同时给 Levins_B 与 Shannon_B，**图仍以 Levins_B 为纵轴**
#   · 图上文字**纯 ASCII**：角标 Core / CondRare / TransRare / PermRare
#
# 依赖共享库：md_core / md_niche_table / md_sample_strategy / md_group_test /
#             md_pal / md_compose / md_save_figure / md_write_xlsx /
#             md_maybe_single / md_empty_plot / pt2gg / md_result
# =============================================================================

#' 模块 12 入口
#' @param CFG config 列表（见 config/config_template.R 的 niche 块）
md_niche_breadth <- function(CFG) {
  module <- "12"

  if (!isTRUE(CFG$niche$enabled %||% TRUE))
    return(md_result(module, "SKIP", "niche$enabled = FALSE"))
  if (!md_need(CFG, c("phyloseq"), module))
    return(md_result(module, "SKIP", "缺少依赖包（phyloseq）"))

  core <- tryCatch(md_core(CFG), error = function(e) NULL)
  if (is.null(core) || is.null(core$asv_raw)) {
    md_skip(CFG, module, "取不到 cache/05_core.rds")
    return(md_result(module, "SKIP", "缺 05 产物"))
  }

  # ---- 参数（全部来自 CFG，%||% 给默认）----
  idx      <- CFG$niche$index       %||% "levins"
  lvl      <- CFG$niche$level       %||% "both"
  inp      <- CFG$niche$input       %||% "raw"
  thr_a    <- CFG$niche$thr_abund   %||% 0.1
  thr_p    <- CFG$niche$thr_prev    %||% 10
  sens_occ <- CFG$niche$sensitivity_occ   %||% c(0.1, 0.3, 0.5, 0.7, 0.9)
  sens_ab  <- CFG$niche$sensitivity_abund %||% c(0.0001, 0.001, 0.01)
  st       <- CFG$niche$style %||% list()
  pt       <- md_pt(CFG, comp = FALSE)
  digits   <- CFG$table$digits %||% 3

  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  # ---- 输入矩阵口径（进缓存键：input 是计算参数，不是样式）----
  m <- if (identical(inp, "rarefied")) (core$asv_rare %||% core$asv_raw) else core$asv_raw
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
    md_skip(CFG, module, "生态位输入矩阵为空")
    return(md_result(module, "SKIP", "输入矩阵为空"))
  }
  m <- as.matrix(m)
  storage.mode(m) <- "double"

  meta <- core$meta
  smp  <- CFG$data$sample_col %||% "sample_id"
  if (!is.null(meta)) meta <- meta[match(colnames(m), meta[[smp]]), , drop = FALSE]

  # ---- 分组变量（样本级堆叠柱的横轴）----
  #   默认 = **前两个设计因子的交互组**（与原图 Fig1c 的 x 轴一致：
  #   1_community_composition.R:691-698 的 interaction_order 就是 Zone×Vegetation，
  #   顺序为「第一因子为主序、第二因子为次序」）；
  #   CFG$niche$groups 给字符串时强制单因子。
  facs <- CFG$design$factors %||% character(0)
  facs <- if (is.null(meta)) character(0) else facs[facs %in% names(meta)]
  grp_cfg <- CFG$niche$groups
  grp_facs <- if (!is.null(grp_cfg) && length(grp_cfg) == 1 && nzchar(as.character(grp_cfg))) {
    intersect(as.character(grp_cfg), names(meta))
  } else {
    utils::head(facs, 2)
  }
  fac_lv <- lapply(grp_facs, function(f) {
    ol  <- CFG$design$ordered_levels[[f]]
    obs <- unique(as.character(meta[[f]])); obs <- obs[!is.na(obs)]
    if (is.null(ol)) return(obs)
    hit <- intersect(as.character(ol), obs)
    if (length(hit)) hit else obs
  })
  # 交互水平顺序：第一因子为主序、第二因子为次序（t(outer(...)) 才是行优先）
  grp_lv <- if (length(grp_facs) >= 2) {
    as.vector(t(outer(fac_lv[[1]], fac_lv[[2]], paste, sep = "_")))
  } else if (length(grp_facs) == 1) {
    fac_lv[[1]]
  } else character(0)
  grp_var <- if (length(grp_facs)) paste(grp_facs, collapse = "_") else NULL
  grp_of  <- if (length(grp_facs) >= 2) {
    paste(as.character(meta[[grp_facs[1]]]), as.character(meta[[grp_facs[2]]]), sep = "_")
  } else if (length(grp_facs) == 1) {
    as.character(meta[[grp_facs[1]]])
  } else {
    rep("All samples", ncol(m))
  }

  # ===========================================================================
  # 分析段 —— 整体包进 md_stage()，图与表在段外
  # ===========================================================================
  md_stage(CFG, "12_analysis", {

    niche_df <- md_niche_table(m, thr_a, thr_p, index = idx)

    # 原文 classify_2x2() 的 Cell 顺序（00c_setup_and_data_preparation.R:1103-1112）
    cell_levels <- c("Core", "Permanently_rare", "Conditionally_rare", "Transient_rare")
    niche_df$Cell <- factor(niche_df$Cell, levels = cell_levels)

    cell_summary <- as.data.frame(table(cell = niche_df$Cell), stringsAsFactors = FALSE)
    names(cell_summary) <- c("Cell", "n_ASV")
    cell_summary$pct_ASV <- 100 * cell_summary$n_ASV / max(1L, sum(cell_summary$n_ASV))
    cell_summary$mean_Levins_B <- vapply(as.character(cell_summary$Cell), function(cl) {
      v <- niche_df$Levins_B[!is.na(niche_df$Cell) & as.character(niche_df$Cell) == cl]
      if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
    }, numeric(1))
    cell_summary$sum_Levins_B <- vapply(as.character(cell_summary$Cell), function(cl) {
      v <- niche_df$Levins_B[!is.na(niche_df$Cell) & as.character(niche_df$Cell) == cl]
      if (length(v) == 0) NA_real_ else sum(v, na.rm = TRUE)
    }, numeric(1))

    # ---- 样本级生态策略组成（原文 513-568；共享库逐字复刻该口径）----
    sample_strategy <- md_sample_strategy(m, niche_df, a_thr = thr_a, p_thr = thr_p)
    sample_strategy <- sample_strategy[match(colnames(m), sample_strategy$SampleID), , drop = FALSE]
    sample_strategy$SampleID <- colnames(m)

    # ---- 样本统计（原文 533-538 的 Richness / Shannon / Pielou）----
    rel_s <- sweep(m, 2, pmax(colSums(m), 1), "/")
    shan <- apply(rel_s, 2, function(p) {
      p <- p[p > 0]
      if (length(p) < 2) NA_real_ else -sum(p * log(p))
    })
    rich <- colSums(m > 0)
    sample_stats <- data.frame(
      SampleID = colnames(m),
      Richness = as.integer(rich),
      Shannon  = as.numeric(shan),
      Pielou_Evenness = ifelse(rich > 1, as.numeric(shan) / log(pmax(rich, 2)), NA_real_),
      stringsAsFactors = FALSE
    )
    sample_stats <- cbind(sample_stats,
                          sample_strategy[, setdiff(names(sample_strategy), "SampleID"),
                                          drop = FALSE])
    sample_stats$Group <- grp_of
    # 参与交互的两个因子各存一列（便于回溯），单因子时只存一列
    for (f in grp_facs) sample_stats[[f]] <- as.character(meta[[f]])
    if (!is.null(meta)) {
      cov_cols <- setdiff(intersect(CFG$design$factors %||% character(0), names(meta)), grp_facs)
      for (cc in cov_cols) sample_stats[[cc]] <- as.character(meta[[cc]])
    }

    # ---- 阈值敏感性扫描（原文 SCI-8；这里单位统一为百分比）----
    sens <- expand.grid(Occ_Threshold_pct = as.numeric(sens_occ) * 100,
                        Abund_Threshold_pct = as.numeric(sens_ab) * 100,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    sens$N_Core <- vapply(seq_len(nrow(sens)), function(i) {
      sum(niche_df$Occurrence > sens$Occ_Threshold_pct[i] &
            niche_df$Mean_Abund > sens$Abund_Threshold_pct[i], na.rm = TRUE)
    }, integer(1))
    sens$N_ASV_total <- nrow(niche_df)
    sens$Pct_Core <- 100 * sens$N_Core / max(1L, nrow(niche_df))

    # ---- 样本级指标的因子检验（走共享 md_group_test，不重写 KW/Dunn）----
    strat_pct <- grep("_pct$", names(sample_strategy), value = TRUE)
    sample_metrics <- c("Richness", "Shannon", "Pielou_Evenness", strat_pct)
    factor_vars <- intersect(CFG$design$factors %||% character(0),
                             names(sample_stats))
    gradient_test <- list()
    for (fv in factor_vars) {
      for (mv in sample_metrics) {
        if (!mv %in% names(sample_stats)) next
        gt <- tryCatch(md_group_test(CFG, sample_stats, mv, fv,
                                     test = CFG$alpha$test %||% "kruskal_dunn",
                                     p_adjust = CFG$alpha$p_adjust %||% "BH",
                                     min_n = 3),
                       error = function(e) NULL)
        if (is.null(gt) || is.null(gt$test) || nrow(gt$test) == 0) next
        tt <- gt$test
        tt$metric <- mv
        tt$factor <- fv
        tt$n_group <- nrow(gt$desc)
        gradient_test[[paste(fv, mv, sep = "::")]] <- tt
      }
    }
    gradient_test <- if (length(gradient_test)) do.call(rbind, gradient_test) else NULL

    # ---- 四象限散点数据（原文 824-841：只保留 Levins_B 非 NA 且平均丰度 > 0 的 ASV）----
    scatter_df <- niche_df[!is.na(niche_df$Levins_B) & !is.na(niche_df$Mean_Abund) &
                             niche_df$Mean_Abund > 0 & niche_df$Levins_B > 0, , drop = FALSE]
    if (nrow(scatter_df) > 0) scatter_df$Cell <- droplevels(scatter_df$Cell)

    # ---- 组平均样本策略组成（堆叠柱数据）----
    strat_long <- do.call(rbind, lapply(strat_pct, function(cc) {
      data.frame(Group = as.character(sample_stats$Group),
                 Cell = sub("_pct$", "", cc),
                 pct = as.numeric(sample_stats[[cc]]),
                 stringsAsFactors = FALSE)
    }))
    strat_long <- strat_long[!is.na(strat_long$Group), , drop = FALSE]
    stack_df <- if (nrow(strat_long)) {
      agg <- stats::aggregate(pct ~ Group + Cell, data = strat_long, FUN = mean)
      n_by <- table(strat_long$Group)
      agg$n_sample <- as.integer(n_by[agg$Group])
      agg$Cell <- factor(agg$Cell, levels = cell_levels)
      agg$Group <- factor(agg$Group, levels = if (length(grp_lv)) grp_lv else unique(agg$Group))
      agg[order(agg$Group, agg$Cell), , drop = FALSE]
    } else NULL
  })

  # ---- 配色：一律走 12_plot_theme.R 的 MD_BUILTIN_PAL$quadrant ----
  cell_cols <- CFG$niche$style$cell_colors %||%
    md_pal(CFG, c("Core", "Permanently_rare", "Conditionally_rare", "Transient_rare"),
           name = "quadrant")
  cell_cols <- cell_cols[c("Core", "Permanently_rare", "Conditionally_rare",
                           "Transient_rare")]

  # ---------------------------------------------------------------------------
  # 出图段（每次运行都执行；样式不进缓存键）
  #
  # 版式 = 毕业论文 Fig1c_niche_rarity（1_community_composition.R:679-891）：
  #   上排 A  样本生态策略组成堆叠柱（x = 前两个设计因子的交互组）  <- fig07
  #   下排 B  丰度 x 出现率四象限散点                               <- p_07sa（原在 10 号）
  #   下排 C  Levins B vs 出现率 + gam 平滑 + 95% CI                 <- p_07sb
  #   原图 200 x 150 mm；本模块同尺寸。
  # ---------------------------------------------------------------------------
  ps <- as.numeric(st$point_size %||% MD_PT_SIZE$scatter)
  pa <- as.numeric(st$point_alpha %||% MD_PT_SIZE$point_alpha)
  gam_cfg <- st$gam %||% list()

  # ---- A. 样本生态策略组成堆叠柱（原 fig07：alpha 0.85 + 黑边 0.3 + 45 度两行标签）----
  if (is.null(stack_df) || nrow(stack_df) == 0) {
    pA <- md_empty_plot(CFG, "No sample strategy composition")
  } else {
    pA <- ggplot2::ggplot(stack_df,
                          ggplot2::aes(x = .data$Group, y = .data$pct, fill = .data$Cell)) +
      ggplot2::geom_bar(stat = "identity", position = "stack",
                        alpha = as.numeric(st$inter_alpha %||% 0.85),
                        colour = "black",
                        linewidth = as.numeric(st$inter_lw %||% 0.3)) +
      ggplot2::scale_fill_manual(values = cell_cols, name = "Cell", drop = FALSE,
                                 guide = md_guide_flex(4, "bottom")) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 101),
                                  labels = function(x) paste0(x, "%")) +
      ggplot2::scale_x_discrete(labels = function(v) gsub("_", "\n", v, fixed = TRUE)) +
      ggplot2::labs(x = NULL, y = "Mean relative abundance (%)", tag = "A") +
      theme_md(CFG = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                         size = pt$tick),
                     legend.position = "bottom",
                     legend.box = "vertical")
  }

  # ---- B. 丰度 x 出现率四象限散点（原 p_07sa，逐项照抄 1_community_composition.R:850-867）----
  #   注意：原脚本这一格**没有**任何拟合曲线；阈值线固定 0.1% / 10%（不读 CFG 阈值 ——
  #   原脚本就是写死的 0.1 / 10，与四象限分类的默认阈值一致）。
  quad_ok <- niche_df[is.finite(niche_df$Mean_Abund) & niche_df$Mean_Abund > 0, , drop = FALSE]
  if (nrow(quad_ok) == 0) {
    pB <- md_empty_plot(CFG, "No ASV with positive mean abundance")
  } else {
    pB <- ggplot2::ggplot(quad_ok, ggplot2::aes(x = .data$Occurrence,
                                                y = .data$Mean_Abund,
                                                colour = .data$Cell)) +
      ggplot2::geom_vline(xintercept = 10, linetype = "dashed",
                          colour = "grey40", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed",
                          colour = "grey40", linewidth = 0.3) +
      ggplot2::geom_point(size = ps, alpha = pa) +
      ggplot2::scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1),
                             labels = c("0.001", "0.01", "0.1", "1")) +
      ggplot2::annotation_logticks(sides = "l", linewidth = 0.2) +
      ggplot2::scale_colour_manual(values = cell_cols, name = "Cell", drop = FALSE,
                                   guide = "none") +
      ggplot2::coord_cartesian(xlim = c(0, 40)) +
      ggplot2::annotate("text", x = 30, y = 2, label = "Core",
                        size = pt2gg(pt$ann), colour = unname(cell_cols["Core"])) +
      ggplot2::annotate("text", x = 4, y = 2, label = "CondRare",
                        size = pt2gg(pt$ann), colour = unname(cell_cols["Conditionally_rare"])) +
      ggplot2::annotate("text", x = 4, y = 0.004, label = "TransRare",
                        size = pt2gg(pt$ann), colour = unname(cell_cols["Transient_rare"])) +
      ggplot2::annotate("text", x = 30, y = 0.004, label = "PermRare (N~0)",
                        size = pt2gg(pt$ann), colour = unname(cell_cols["Permanently_rare"])) +
      ggplot2::labs(x = "Occurrence (%)", y = "Mean relative abundance (%)", tag = "B") +
      theme_md(CFG = NULL) +
      ggplot2::theme(legend.position = "none")
  }

  # ---- C. Levins B vs 出现率 + gam 平滑 + 95% CI（原 p_07sb，逐项照抄 1_community_composition.R:869-879）----
  #   ⚠ 原式 `s(x, bs = "cs")` 用的是 mgcv 默认 k = 10。单例化严重的数据里出现率只有
  #     少数几个离散取值时会报 "x has insufficient unique values to support 10 knots"。
  #     这里把 k 自适应到 `min(10, 唯一值个数 - 1)`；`CFG$niche$style$gam$k` 可显式覆盖。
  #     唯一取值 < 4 时**不画平滑**（画出来必是伪曲线），只留散点并在 notes 说明。
  n_ux <- if (nrow(scatter_df)) length(unique(scatter_df$Occurrence)) else 0L
  k_sm <- as.integer(gam_cfg$k %||% max(3L, min(10L, n_ux - 1L)))
  if (nrow(scatter_df) == 0) {
    pC <- md_empty_plot(CFG, "No ASV with measurable niche breadth")
  } else {
    pC <- ggplot2::ggplot(scatter_df,
                          ggplot2::aes(x = .data$Occurrence, y = .data$Levins_B,
                                       colour = .data$Cell)) +
      ggplot2::geom_vline(xintercept = 10, linetype = "dashed",
                          colour = "grey40", linewidth = 0.3) +
      ggplot2::geom_point(size = ps * 1.15, alpha = pa)
    if (n_ux >= 4L) {
      # 原脚本：geom_smooth(aes(group = 1), method = "gam", formula = y ~ s(x, bs = "cs"),
      #                     se = TRUE, color = "#D55E00", linewidth = 0.6, alpha = 0.15)
      pC <- pC + ggplot2::geom_smooth(
        ggplot2::aes(group = 1),
        method  = as.character(gam_cfg$method %||% "gam"),
        formula = stats::as.formula(sprintf("y ~ s(x, bs = '%s', k = %d)",
                                            as.character(gam_cfg$bs %||% "cs"), k_sm)),
        se = TRUE,
        colour = as.character(gam_cfg$colour %||% "#D55E00"),
        linewidth = as.numeric(gam_cfg$lw %||% 0.6),
        fill = as.character(gam_cfg$colour %||% "#D55E00"),
        alpha = as.numeric(gam_cfg$alpha %||% 0.15))
    }
    pC <- pC +
      ggplot2::scale_colour_manual(values = cell_cols, drop = FALSE, guide = "none") +
      ggplot2::scale_y_log10() +
      ggplot2::coord_cartesian(xlim = c(0, 40)) +
      ggplot2::labs(x = "Occurrence (%)", y = "Levins B (niche breadth)", tag = "C") +
      theme_md(CFG = NULL) +
      ggplot2::theme(legend.position = "none")
  }

  fig <- md_compose(CFG, list(A = pA, B = pB, C = pC), design = "AAAA\nBBCC",
                    guides = "keep")
  md_save_figure(CFG, fig, "12_niche_breadth", width = 200, height = 150)
  md_maybe_single(CFG, pA, "12_niche_breadth_A")
  md_maybe_single(CFG, pB, "12_niche_breadth_B")
  md_maybe_single(CFG, pC, "12_niche_breadth_C")

  # ---------------------------------------------------------------------------
  # 写表段
  # ---------------------------------------------------------------------------
  nb_out <- data.frame(
    ASV_ID     = niche_df$ASV_ID,
    Mean_Abund = round(niche_df$Mean_Abund, 6),
    Occurrence = round(niche_df$Occurrence, digits),
    Cell       = as.character(niche_df$Cell),
    Levins_B   = round(niche_df$Levins_B, digits),
    Shannon_B  = round(niche_df$Shannon_B, digits),
    N_samples  = niche_df$N_samples,
    stringsAsFactors = FALSE)
  if (!is.null(core$tax)) {
    nb_out$Phylum <- as.character(core$tax$Phylum[match(nb_out$ASV_ID, core$tax[[CFG$data$id_col %||% "asv_id"]])])
    if ("Genus" %in% names(core$tax))
      nb_out$Genus <- as.character(core$tax$Genus[match(nb_out$ASV_ID, core$tax[[CFG$data$id_col %||% "asv_id"]])])
  }

  ss_out <- sample_stats
  for (cc in grep("_pct$|_Abundance$", names(ss_out), value = TRUE))
    ss_out[[cc]] <- round(as.numeric(ss_out[[cc]]), digits)
  if ("Shannon" %in% names(ss_out)) ss_out$Shannon <- round(ss_out$Shannon, digits)
  if ("Pielou_Evenness" %in% names(ss_out)) ss_out$Pielou_Evenness <- round(ss_out$Pielou_Evenness, digits)

  note <- data.frame(
    item = c("Input matrix", "Normalization", "Index", "Level",
             "Abundance threshold (%)", "Occurrence threshold (%)",
             "Cell definitions", "Levins B formula", "Shannon B formula",
             "Panel A", "Panel B", "Grouping factor",
             "Sample strategy unit", "n_ASV", "n_sample",
             "Corner labels"),
    value = c(
      paste0(if (identical(inp, "rarefied")) "asv_rare (rarefied)" else "asv_raw (not rarefied)",
             " x ", nrow(m), " ASV x ", ncol(m), " samples"),
      as.character(CFG$normalize$by_module$composition %||% CFG$normalize$method),
      idx, lvl, as.character(thr_a), as.character(thr_p),
      paste0("mean rel. abund > ", thr_a, "% and occurrence > ", thr_p,
             "% -> Core; abund > thr & occ <= thr -> Conditionally_rare; ",
             "abund <= thr & occ > thr -> Permanently_rare; else Transient_rare"),
      "B = 1 / sum(p^2), p = ASV relative abundance across samples (NA if present in < 2 samples)",
      "B = -sum(p * log(p)) / log(n_samples) (Shannon-type evenness of the abundance distribution)",
      paste0("A = sample ecological-strategy composition stacked bars, x = interaction group (",
             grp_var %||% "All samples", "); ",
             "B = abundance x occurrence 4-quadrant scatter (moved here from module 10, ",
             "matches Fig1c p_07sa); ",
             "C = Levins B vs occurrence with GAM smooth + 95% CI (matches Fig1c p_07sb)"),
      "mean per-sample percentage of each Cell within each group (stack sums to ~100%)",
      grp_var %||% "All samples",
      paste0("per-sample relative abundance (%) of the ASVs assigned to each Cell; ",
             "Cell assignment is done at ASV level over all samples (原始脚本口径)"),
      as.character(nrow(niche_df)), as.character(ncol(m)),
      "Core / CondRare / TransRare / PermRare (N~0)"),
    stringsAsFactors = FALSE)
  note <- rbind(note, data.frame(item = c("Panel B thresholds", "Panel C smooth"),
    value = c("fixed at 0.1% abundance / 10% occurrence (hard-coded in the original p_07sa)",
              "geom_smooth(method='gam', formula=y~s(x,bs='cs'), se=TRUE, colour='#D55E00', alpha=0.15)"),
    stringsAsFactors = FALSE))
  if (identical(lvl, "sample") && nrow(scatter_df) > 0)
    note <- rbind(note, data.frame(item = "Note",
      value = paste0("niche$level = sample: panel A is still produced because it is the ",
                     "niche-breadth chart requested by the spec."),
      stringsAsFactors = FALSE))

  tp <- md_write_xlsx(CFG, "12_niche_breadth", list(
    Niche_Breadth         = nb_out,
    Cell_summary          = cell_summary,
    Sample_strategy       = sample_strategy,
    Threshold_sensitivity = sens,
    Sample_stats          = ss_out,
    Gradient_test         = gradient_test,
    Note                  = note))

  md_result(module, "OK", NULL,
            tables  = "12_niche_breadth.xlsx",
            figures = "figures/12_niche_breadth.pdf",
            notes = c(notes,
                      paste0("index=", idx, "; level=", lvl, "; input=", inp),
                      paste0("Core ASV = ", sum(niche_df$Cell == "Core", na.rm = TRUE),
                             " / ", nrow(niche_df),
                             " (nifH data is singleton-dominated; a near-empty Core ",
                             "quadrant is a real feature of the data, not a bug)")))
}
