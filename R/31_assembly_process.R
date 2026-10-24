# =============================================================================
# R/31_assembly_process.R -- 群落组装过程判定与比例汇总（Stegen 框架）
# -----------------------------------------------------------------------------
# 本模块**只读** cache/30_null.rds（经 md_null_cache(CFG)），不重算零模型。
# 内容：
#   1. Stegen 五过程判定复核 + 比例汇总，四种口径：
#        总体 / 按因子 1（同因子 1 的样本对）/ 按因子 2 / 按因子交互
#        -- 与 nifH 原实现 Overall / By_Zone / By_Vegetation / By_Both 对齐
#   2. betaNTI × RCbray 四象限散点（阈值来自 CFG$assembly$bnti_thr / rcbray_thr）
#   3. 图（画法逐项对照 nifH 原脚本）：
#        A = 过程比例堆叠柱      ← 原脚本 §B
#        B = betaNTI×RCbray 四象限 ← 原脚本（附图C）
#        C = **可选备用面板**（默认关）「样本级随机性 vs 平均群落周转」
#                                  ← 原脚本 §D的画法
#
# 可选备用面板（默认关，不开时本模块仍只出 A/B 两格，默认出图张数不变）：
#   CFG$spare$process_vs_turnover <- list(
#     enabled = FALSE,      # TRUE 才出这一格
#     x       = NULL,       # 必填：横轴列。样本级过程名（如 "Drift"，画成 %）或分析总表里的数值列
#     y       = "MeanBeta", # 纵轴：该样本与**其余全部样本**的平均群落距离（也可写分析总表里的列）
#     group   = NULL,       # 可选：分组着色（走 md_pal()）
#     covar   = NULL,       # 可选：加进 lm/gam 的协变量
#     fit     = "lm"        # lm | gam
#   )
#   x 的过程占比口径 = pair_df 展平到样本级后样本内归一化
#     （同 nifH 原脚本）；y 的平均周转取自模块 13 的距离矩阵
#     （md_dist_get，逐行对非对角元素取均值），**本模块不自己算距离、不重算零模型**。
#
# **距离衰减 DDR 不在本模块** -- 那是 33_distance_decay（用户明确要求分开）。
#
# 产物：results/31_assembly_process.xlsx + figures/31_assembly_process.pdf
# =============================================================================

md_assembly_process <- function(CFG) {
  module <- "31"
  t_start <- Sys.time()
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少依赖包"))
  if (!md_precheck(CFG, module, need_response = FALSE, need_phylo = FALSE))
    return(md_result(module, "SKIP", "前置条件不满足"))

  nc <- tryCatch(md_null_cache(CFG), error = function(e) NULL)
  if (is.null(nc) || is.null(nc$pair_df) || nrow(nc$pair_df) == 0) {
    md_skip(CFG, module, "缺少上游产物 cache/30_null.rds（请先跑模块 30）")
    return(md_result(module, "SKIP", "缺少 30 的产物"))
  }
  pair_df <- nc$pair_df
  params  <- nc$params %||% list()
  bnti_thr <- params$bnti_thr %||% CFG$assembly$bnti_thr %||% 2
  rc_thr   <- params$rcbray_thr %||% CFG$assembly$rcbray_thr %||% 0.95

  factors  <- CFG$design$factors %||% character(0)
  f1 <- factors[1] %||% NA_character_
  f2 <- if (length(factors) > 1) factors[2] else NA_character_
  proc_levels <- c("Variable selection", "Homogeneous selection",
                   "Dispersal limitation", "Homogenizing dispersal", "Drift")
  # 配色走主题：md_pal() 会自动命中 MD_BUILTIN_PAL$assembly（= nifH cols_assembly）
  proc_cols <- local31_process_cols(CFG, proc_levels)
  pt <- md_pt(CFG)

  # ---- 1. 过程判定复核（不改动 30 的判定，只做一致性检查，异常时按本模块阈值重判）----
  if (!"process" %in% names(pair_df) || all(is.na(pair_df$process))) {
    md_warn(CFG, "pair_df 缺 process 列，按本模块阈值重新判定")
  } else {
    expect <- local31_process(pair_df$bnti, pair_df$rcbray, bnti_thr, rc_thr)
    mism <- sum(as.character(pair_df$process) != expect, na.rm = TRUE)
    if (mism > 0)
      md_warn(CFG, "过程判定与 30 的记录不一致 ", mism, " 处，已按 CFG 阈值重判")
  }
  pair_df$process <- local31_process(pair_df$bnti, pair_df$rcbray, bnti_thr, rc_thr)
  n_pair <- nrow(pair_df)

  # ---- 2. 比例汇总：总体 / f1 / f2 / 交互 ----
  overall <- local31_prop_table(pair_df, "Overall", proc_levels)
  by_f1 <- NULL; by_f2 <- NULL; by_both <- NULL
  if (!is.na(f1)) {
    by_f1 <- local31_prop_by(pair_df, paste0(f1, "1"), paste0("Same_f", 1), proc_levels)
  }
  if (!is.na(f2)) {
    by_f2 <- local31_prop_by(pair_df, paste0(f2, "1"), paste0("Same_f", 2), proc_levels)
  }
  if (!is.na(f1) && !is.na(f2)) {
    both_col <- paste0(f1, "1"); veg_col <- paste0(f2, "1")
    d <- pair_df[which(pair_df$Same_f1 %in% TRUE & pair_df$Same_f2 %in% TRUE), , drop = FALSE]
    if (nrow(d) > 0 && both_col %in% names(d) && veg_col %in% names(d)) {
      d$.both <- paste(as.character(d[[both_col]]), as.character(d[[veg_col]]), sep = " | ")
      by_both <- local31_prop_by(d, ".both", NULL, proc_levels)
    }
  }

  # ---- 3. 统计：betaNTI / RCbray 描述 + 四象限计数 ----
  desc <- rbind(local31_desc(pair_df$bnti, "betaNTI"),
                local31_desc(pair_df$rcbray, "RCbray"))
  quad <- local31_quadrant_table(pair_df, bnti_thr, rc_thr)
  data_feat <- local31_feature_table(pair_df, bnti_thr)

  # ---- 4. 图 A：过程比例堆叠柱（总体 + 各口径）----
  # 画法沿用同一套样式 §B：geom_bar(stat="identity", width=0.75,
  #   color="white", linewidth=0.2) + scale_fill_manual(cols_assembly, drop=FALSE)
  #   + labs(x=..., y="Proportion (%)") + 45 度 x 轴标签。
  # （md_bar_stack() 把柱宽固定 0.65、描边固定 NA，表达不了原图的白描边，
  #   故这里按原脚本那几行写，其余常量仍全部走 12_plot_theme.R。）
  ba <- rbind(
    data.frame(group = "Overall", overall[, c("process", "prop")], stringsAsFactors = FALSE),
    if (!is.null(by_f1)) data.frame(group = paste0(f1, ": ", by_f1$group),
                                    by_f1[, c("process", "prop")], stringsAsFactors = FALSE),
    if (!is.null(by_f2)) data.frame(group = paste0(f2, ": ", by_f2$group),
                                    by_f2[, c("process", "prop")], stringsAsFactors = FALSE))
  ba$prop <- 100 * ba$prop
  ba$process <- factor(ba$process, levels = proc_levels)
  ba$group <- factor(ba$group, levels = unique(ba$group))
  xlab_bar <- paste(c(f1, f2)[!is.na(c(f1, f2))], collapse = " x ")
  pA <- ggplot2::ggplot(ba, ggplot2::aes(x = .data$group, y = .data$prop,
                                         fill = .data$process)) +
    ggplot2::geom_bar(stat = "identity", width = 0.75, colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = proc_cols, drop = FALSE,
                               # 图例在下方 → 横向排成一行（bottom_max = 项目数），别竖着一列浪费版面
                               guide = md_guide_flex(length(proc_cols), "bottom",
                                                     length(proc_cols))) +
    ggplot2::labs(x = xlab_bar, y = "Proportion (%)", tag = "A", fill = NULL) +
    theme_md(CFG = CFG) +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical",
                   # 沿用同一套样式 §B：x 轴标签一律 45 度（标签保持单行）
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = pt$tick))

  # ---- 图 B：betaNTI × RCbray 四象限散点 ----
  # 画法沿用同一套样式 的四象限散点（附图C）：
  #   只画 ±阈值 参考线（不画 0 线）、按过程着色（cols_assembly, drop=FALSE）、
  #   图例名 "Ecological Process"、右下角副标题给成对比较数。
  dd <- pair_df[is.finite(pair_df$bnti) & is.finite(pair_df$rcbray), , drop = FALSE]
  if (nrow(dd) >= 5) {
    dd$process <- factor(dd$process, levels = proc_levels)
    pB <- md_quadrant(dd, x = "bnti", y = "rcbray",
                      hline = c(-rc_thr, rc_thr), vline = c(-bnti_thr, bnti_thr),
                      color = "process", colors = NULL,
                      xlab = expression(beta * "NTI"), ylab = "RCbray", tag = "B",
                      pt = pt) +
      ggplot2::scale_colour_manual(values = proc_cols, drop = FALSE,
                                   name = "Ecological Process") +
      ggplot2::labs(subtitle = sprintf("n = %d pairwise comparisons", nrow(dd))) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = pt$title,
                                                           colour = "grey40"))
  } else {
    pB <- md_empty_plot(CFG, "beta*NTI x RCbray unavailable (n < 5)")
  }

  # ---- 图 C（可选备用面板）：任意样本级指标 vs 平均群落周转 ----
  # 对照 原脚本 §D（散点 + geom_smooth + 报 adj R2）。
  # **默认关**：CFG$spare$process_vs_turnover$enabled = FALSE 时不画这一格，
  # 本模块仍按「过程比例堆叠柱 + 四象限散点」两个面板出图。
  # 打开方式见本模块头部注释与 Note sheet。
  spare <- local31_spare_panel(CFG, pair_df, proc_levels, proc_cols, pt)
  pC <- if (is.null(spare)) NULL else spare$plot

  panels <- list(A = pA, B = pB)
  design <- "AB"
  if (!is.null(pC)) { panels$C <- pC; design <- "AB\nCC" }
  p <- md_compose(CFG, panels, design = design, guides = "collect", legend = "bottom")
  md_save_figure(CFG, p, "31_assembly_process",
                 height = if (is.null(pC)) 100 else 135)

  # ---- 5. 表 ----
  note <- data.frame(item = c(
    "数据来源", "过程判定", "比例口径", "成分数据警示", "样本量口径", "DDR",
    "可选备用面板"),
    value = c(
      "cache/30_null.rds（模块 30 的 betaNTI / RCbray），本模块不重算零模型",
      sprintf("|betaNTI| > %g：< -%g 同质选择 / > %g 异质选择；|RCbray| > %g：> %g 扩散限制 / < -%g 均质化扩散；其余漂变",
              bnti_thr, bnti_thr, bnti_thr, rc_thr, rc_thr, rc_thr),
      "总体 / 按因子1（仅同因子1的样本对）/ 按因子2 / 按因子交互（同因子1且同因子2）",
      "各口径内过程比例为成分数据（sum = 1），同一模型放入多个过程会完全共线，回归须用单焦点 + 控制变量",
      sprintf("成对口径 n = %d /（%d 样本）；betaNTI 有效 %d 对、RCbray 有效 %d 对",
              n_pair, params$n_sample %||% NA_integer_,
              sum(is.finite(pair_df$bnti)), sum(is.finite(pair_df$rcbray))),
      "距离衰减关系（DDR）**不在本模块**，见 33_distance_decay",
      local31_spare_note(CFG, spare)),
    stringsAsFactors = FALSE)

  sheets <- list(
    Process_overall = overall,
    By_factor1      = by_f1,
    By_factor2      = by_f2,
    By_interaction  = by_both,
    BetaNTI_RCbray  = desc,
    Quadrant        = quad,
    Data_features   = data_feat,
    Note            = note)
  # 备用面板开启时，把它用的样本级两列一并交出去（关闭时不出这个 sheet）
  if (!is.null(spare) && !is.null(spare$data)) sheets$Spare_process_vs_turnover <- spare$data
  tp <- md_write_xlsx(CFG, "31_assembly_process", sheets,
                      meta_extra = list(bnti_thr = bnti_thr, rcbray_thr = rc_thr,
                                        n_pair = n_pair))

  md_log(CFG, "OK", "模块 31 用时 ",
         md_fmt_ms(0, as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  md_result(module, "OK", tables = basename(tp), figures = "figures/31_assembly_process.pdf",
            notes = sprintf("成对 %d | 主导过程 %s", n_pair, overall$process[1]))
}

# -----------------------------------------------------------------------------
# local31_ 辅助函数
# -----------------------------------------------------------------------------

#' 组装过程配色 —— 走主题：md_pal(CFG, levels, "assembly") 会自动命中
#' MD_BUILTIN_PAL$assembly（**逐字照搬 原实现 cols_assembly**，见 12_plot_theme.R）。
#' 只有用户显式在 CFG$figure$palette$assembly / $process 里覆盖时才用自定义色。
local31_process_cols <- function(CFG, levels) {
  p <- CFG$figure$palette$assembly %||% CFG$figure$palette$process
  if (!is.null(p) && !is.null(names(p)) && all(levels %in% names(p))) return(p[levels])
  md_pal(CFG, levels, "assembly")
}

# -----------------------------------------------------------------------------
# 可选备用面板：任意**样本级**指标 vs 平均群落周转（默认关，不影响默认出图）
#   CFG$spare$process_vs_turnover = list(enabled = FALSE, x = NULL, y = "MeanBeta",
#                                        group = NULL, covar = NULL, fit = "lm")
# -----------------------------------------------------------------------------

#' 读备用面板配置（键缺省时用 NULL / FALSE；不改 config/*，只在模块内解释）
local31_spare_cfg <- function(CFG) {
  d <- CFG$spare$process_vs_turnover %||% list()
  list(enabled = isTRUE(d$enabled),
       x       = if (is.null(d$x)) NULL else as.character(d$x)[1],
       y       = if (is.null(d$y)) "MeanBeta" else as.character(d$y)[1],
       group   = if (is.null(d$group)) NULL else as.character(d$group)[1],
       covar   = if (is.null(d$covar)) NULL else as.character(d$covar),
       fit     = if (identical(d$fit, "gam")) "gam" else "lm")
}

#' 每个样本与**其余全部样本**的平均群落距离（面板 y 的默认口径）。
#' 距离矩阵来自模块 13（md_dist_get），本模块不自己算距离、不重算零模型。
local31_mean_dist <- function(CFG, metric) {
  d <- tryCatch(md_dist_get(CFG, metric), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  m <- as.matrix(d)
  if (nrow(m) < 3L) return(NULL)
  diag(m) <- NA_real_
  v <- rowMeans(m, na.rm = TRUE)
  v[!is.finite(v)] <- NA_real_
  data.frame(sample = rownames(m), .value = as.numeric(v), stringsAsFactors = FALSE)
}

#' 样本级过程占比（口径同 nifH 原脚本：
#' 把 pair_df 的 Sample1 / Sample2 展平后按样本计数、样本内归一化）
local31_sample_prop <- function(pair_df, proc_levels) {
  s <- c(as.character(pair_df$sample1), as.character(pair_df$sample2))
  p <- c(as.character(pair_df$process),  as.character(pair_df$process))
  keep <- !is.na(s) & !is.na(p)
  s <- s[keep]; p <- p[keep]
  if (!length(s)) return(NULL)
  tab <- table(factor(s, levels = unique(s)), factor(p, levels = proc_levels))
  out <- as.data.frame.matrix(tab, stringsAsFactors = FALSE)
  out$sample <- rownames(out)
  out$n_pairs <- rowSums(out[, proc_levels, drop = FALSE])
  for (k in proc_levels) out[[k]] <- out[[k]] / pmax(out$n_pairs, 1)
  out[, c("sample", "n_pairs", proc_levels), drop = FALSE]
}

#' 画备用面板。返回 NULL = 不出这一格（未启用 / 没给 x / 数据不足 / y 取不到）。
#' 画法照抄 nifH 原脚本 §D：散点 + geom_smooth + 面板内报 adj R2；
#' 分组着色走 md_pal()，线条/散点/字号全部走 12_plot_theme.R。
local31_spare_panel <- function(CFG, pair_df, proc_levels, proc_cols, pt) {
  sp <- local31_spare_cfg(CFG)
  if (!sp$enabled || is.null(sp$x)) return(NULL)

  metric <- CFG$response_models$pair_metric %||% CFG$distances$metrics[1]
  smp    <- CFG$data$sample_col
  prop   <- local31_sample_prop(pair_df, proc_levels)
  ana    <- tryCatch(md_analysis_with_derived(CFG), error = function(e) NULL)

  base <- if (!is.null(prop)) prop else
    data.frame(sample = sort(unique(c(as.character(pair_df$sample1),
                                      as.character(pair_df$sample2)))),
               stringsAsFactors = FALSE)
  if (!is.null(ana) && smp %in% names(ana)) {
    keep <- setdiff(names(ana), intersect(names(ana), names(base)))
    base <- merge(base, ana[, c(smp, keep), drop = FALSE],
                  by.x = "sample", by.y = smp, all.x = TRUE)
  }

  # ---- 解析 x：过程名（样本级占比，%）或样本级表里的任意数值列 ----
  x_lab <- sp$x
  if (sp$x %in% proc_levels) {
    base$.x <- 100 * suppressWarnings(as.numeric(base[[sp$x]]))
    x_lab   <- paste0(sp$x, " (%)")
  } else if (sp$x %in% names(base)) {
    base$.x <- suppressWarnings(as.numeric(base[[sp$x]]))
  } else {
    md_warn(CFG, "31 备用面板：x = '", sp$x,
            "' 既不是过程名，也不在样本级表（过程占比 + 分析总表）里，跳过该面板")
    return(NULL)
  }

  # ---- 解析 y：默认 MeanBeta = 该样本与其余样本的平均群落距离 ----
  y_lab <- sp$y
  if (identical(sp$y, "MeanBeta") || !sp$y %in% names(base)) {
    md_ <- local31_mean_dist(CFG, metric)
    if (is.null(md_)) return(NULL)
    names(md_) <- c("sample", ".y")
    base <- merge(base, md_, by = "sample", all.x = TRUE)
    y_lab <- sprintf("Mean community distance (%s)", metric)
  } else {
    base$.y <- suppressWarnings(as.numeric(base[[sp$y]]))
  }

  if (!is.null(sp$group) && sp$group %in% names(base))
    base$.g <- as.character(base[[sp$group]])
  cov_use <- intersect(sp$covar %||% character(0), names(base))
  for (v in cov_use) base[[v]] <- suppressWarnings(as.numeric(base[[v]]))

  need <- c(".x", ".y", cov_use, if (!is.null(base$.g)) ".g")
  d <- base[stats::complete.cases(base[, need, drop = FALSE]), , drop = FALSE]
  d <- d[is.finite(d$.x) & is.finite(d$.y), , drop = FALSE]
  if (nrow(d) < 5L || length(unique(d$.x)) < 2L) {
    md_warn(CFG, "31 备用面板：有效样本 ", nrow(d), " 个 / x 唯一值 ",
            length(unique(d$.x)), " 个，不足以画散点+拟合，跳过")
    return(NULL)
  }

  # ---- 拟合（只用于面板上的 adj R2 标注；不改变本模块任何已有统计量）----
  fd <- stats::setNames(d[, c(".x", ".y", cov_use), drop = FALSE],
                        c("x", "y", cov_use))
  rhs <- paste(c("x", cov_use), collapse = " + ")
  f <- tryCatch(
    if (identical(sp$fit, "gam") && requireNamespace("mgcv", quietly = TRUE))
      mgcv::gam(stats::as.formula(paste0("y ~ s(x)",
                                         if (length(cov_use)) paste0(" + ",
                                           paste(cov_use, collapse = " + ")) else "")),
                data = fd)
    else stats::lm(stats::as.formula(paste("y ~", rhs)), data = fd),
    error = function(e) NULL)
  r2 <- if (is.null(f)) NA_real_ else tryCatch(
    if (inherits(f, "gam")) summary(f)$r.sq else summary(f)$adj.r.squared,
    error = function(e) NA_real_)
  lab <- sprintf("R2 adj = %s", if (is.finite(r2)) sprintf("%.3f", r2) else "NA")
  if (length(cov_use))
    lab <- paste0(lab, " (adjusted for ", paste(cov_use, collapse = ", "), ")")

  # ---- 画：散点（可分组着色）+ 单条 geom_smooth + 面板内 adj R2 标注 ----
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$.x, y = .data$.y))
  g <- if (!is.null(d$.g))
    g + ggplot2::geom_point(ggplot2::aes(colour = .data$.g),
                            size = MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha)
  else
    g + ggplot2::geom_point(size = MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha,
                            colour = MD_REG$col)
  g <- g +
    ggplot2::geom_smooth(method = if (identical(sp$fit, "gam")) "gam" else "lm",
                         formula = y ~ x, se = TRUE,
                         colour = MD_REG$col, fill = MD_REG$col,
                         alpha = MD_REG$ribbon_alpha, linewidth = MD_REG$lw) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.4,
                      label = lab, size = pt2gg(pt$ann), colour = "grey25") +
    ggplot2::labs(x = x_lab, y = y_lab, tag = "C") +
    theme_md(CFG = CFG)
  g <- if (!is.null(d$.g))
    g + ggplot2::scale_colour_manual(values = md_pal(CFG, unique(d$.g), sp$group),
                                     name = sp$group) +
      ggplot2::theme(legend.position = "bottom")
  else
    g + ggplot2::theme(legend.position = "none")

  out <- data.frame(sample = d$sample, x = d$.x, y = d$.y,
                    group = if (!is.null(d$.g)) d$.g else NA_character_,
                    stringsAsFactors = FALSE)
  names(out)[2] <- sp$x
  names(out)[3] <- sp$y
  list(plot = g, data = out,
       note = sprintf("x = %s（%s）；y = %s（%s）%s；fit = %s",
                      sp$x, x_lab, sp$y, y_lab,
                      if (length(cov_use)) paste0("；协变量 = ", paste(cov_use, collapse = ", ")) else "",
                      sp$fit))
}

#' Note sheet 里那行说明：备用面板默认关，打开时才出图
local31_spare_note <- function(CFG, spare) {
  sp <- local31_spare_cfg(CFG)
  if (is.null(spare))
    return(paste0("【可选备用面板】未出图（CFG$spare$process_vs_turnover$enabled = ",
                  sp$enabled, "）：本模块图件只有「过程比例堆叠柱 + betaNTI x RCbray 四象限」",
                  "两个面板。想加「随机性 vs 群落周转」这一格，置 enabled = TRUE 并给 x",
                  "（样本级过程名或分析总表里的数值列），y 默认 MeanBeta = 该样本与其余样本的",
                  "平均群落距离（md_dist_get 的距离矩阵逐行平均）。画法照 原脚本 §D。"))
  paste0("【可选备用面板】已出图（enabled = TRUE）：", spare$note,
         "；散点 + 单条 geom_smooth + 面板内 adj R2 标注（照 原脚本 §D）；",
         "过程占比口径同 nifH 原脚本（pair_df 展平后样本内归一化）。")
}

local31_process <- function(bnti, rc, bnti_thr, rc_thr) {
  out <- rep("Drift", length(bnti))
  okr <- is.finite(rc)
  out[okr & rc >  rc_thr] <- "Dispersal limitation"
  out[okr & rc < -rc_thr] <- "Homogenizing dispersal"
  okb <- is.finite(bnti)
  out[okb & bnti < -bnti_thr] <- "Homogeneous selection"
  out[okb & bnti >  bnti_thr] <- "Variable selection"
  out
}

#' 单一分组水平下的过程比例表
local31_prop_table <- function(df, label, proc_levels) {
  d <- df[!is.na(df$process), , drop = FALSE]
  n <- nrow(d)
  out <- data.frame(group = label, process = proc_levels,
                    n = as.integer(table(factor(d$process, levels = proc_levels))),
                    stringsAsFactors = FALSE)
  out$prop <- if (n > 0) out$n / n else NA_real_
  out$n_pairs <- n
  out[order(-out$n), , drop = FALSE]
}

#' 按分组列汇总（same_col 不为 NULL 时只保留同该因子的样本对）
local31_prop_by <- function(df, group_col, same_col, proc_levels) {
  if (is.null(group_col) || !group_col %in% names(df)) return(NULL)
  d <- df[!is.na(df[[group_col]]), , drop = FALSE]
  if (!is.null(same_col) && same_col %in% names(df))
    d <- d[which(d[[same_col]] %in% TRUE), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  keys <- sort(unique(as.character(d[[group_col]])))
  out <- do.call(rbind, lapply(keys, function(k) {
    x <- d[as.character(d[[group_col]]) == k, , drop = FALSE]
    tab <- as.integer(table(factor(x$process, levels = proc_levels)))
    data.frame(group = k, process = proc_levels, n = tab,
               prop = tab / sum(tab), n_pairs = nrow(x), stringsAsFactors = FALSE)
  }))
  out$group <- as.character(out$group)
  out
}

local31_desc <- function(x, label) {
  x <- x[is.finite(x)]
  if (!length(x))
    return(data.frame(metric = label, n = 0L, mean = NA_real_, sd = NA_real_,
                      median = NA_real_, min = NA_real_, max = NA_real_,
                      stringsAsFactors = FALSE))
  data.frame(metric = label, n = length(x), mean = mean(x), sd = stats::sd(x),
             median = stats::median(x), min = min(x), max = max(x),
             stringsAsFactors = FALSE)
}

#' 四象限计数（按 CFG 阈值）
local31_quadrant_table <- function(pair_df, bnti_thr, rc_thr) {
  b <- pair_df$bnti; r <- pair_df$rcbray
  glev <- c("betaNTI NA", "|betaNTI| <= thr",
            "|betaNTI| > thr (negative)", "|betaNTI| > thr (positive)")
  rlev <- c("RCbray NA", "RCbray < -thr", "|RCbray| <= thr", "RCbray > +thr")
  cat_g <- ifelse(!is.finite(b), glev[1],
           ifelse(b < -bnti_thr, glev[3],
           ifelse(b >  bnti_thr, glev[4], glev[2])))
  cat_r <- ifelse(!is.finite(r), rlev[1],
           ifelse(r >  rc_thr, rlev[4],
           ifelse(r < -rc_thr, rlev[2], rlev[3])))
  tt <- table(factor(cat_g, levels = glev), factor(cat_r, levels = rlev))
  out <- as.data.frame(tt, stringsAsFactors = FALSE)
  names(out) <- c("BetaNTI_class", "RCbray_class", "n")
  out$prop <- out$n / sum(out$n)
  out[out$n > 0, , drop = FALSE]
}

#' 数据特征清单（本项目数据极端单例化，这些数字要随表一起交出去）
local31_feature_table <- function(pair_df, bnti_thr = 2) {
  data.frame(
    item = c("n_pairs", "betaNTI_valid", "betaNTI_NA", "RCbray_valid",
             "RCbray_at_bound", "n_x_Null", "n_x_nonNull",
             "n_both_available"),
    value = c(nrow(pair_df),
              sum(is.finite(pair_df$bnti)),
              sum(!is.finite(pair_df$bnti)),
              sum(is.finite(pair_df$rcbray)),
              sum(abs(pair_df$rcbray) >= 1 - 1e-12, na.rm = TRUE),
              sum(is.finite(pair_df$bnti) & abs(pair_df$bnti) <= bnti_thr),
              sum(is.finite(pair_df$bnti) & abs(pair_df$bnti) >  bnti_thr),
              sum(is.finite(pair_df$bnti) & is.finite(pair_df$rcbray))),
    stringsAsFactors = FALSE)
}
