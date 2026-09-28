# =============================================================================
# R/shared/13_plot_types.R —— 图形基元
# -----------------------------------------------------------------------------
# 库文件。每个函数只负责"把一张图画出来"，不做统计、不读文件。
# 所有 geom_text/annotate 的 size 一律走 pt2gg()，字号通过参数传入
# （nifH 论文踩过的坑：geom_text 的 size 在构造时固化，事后 + theme() 改不动）。
# =============================================================================

#' 面板为空时的占位图（模块降级时用，避免出 0 字节文件）
md_empty_plot <- function(CFG = NULL, msg = "Skipped") {
  # 注意：光靠 theme_void() 不够 —— md_compose() 末尾会把 axis.text 重新设回 element_text，
  # 占位面板就会重新长出 "-0.05 / 0.00 / 0.05" 这类坐标轴数字。
  # 把两个 scale 的 breaks 清空才是根治：没有 break 就没有刻度文字。
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, size = pt2gg(MD_PT$ann),
                      colour = "grey40") +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::theme_void()
}

# ---------------------------------------------------------------------------
# 箱线图
# ---------------------------------------------------------------------------
#' 单因子箱线 + 抖动点 + CLD 字母 + 可选的 KW 注释
#'
#' 默认样式**完全复刻 nifH**：空心箱体、按因子给箱体边框与抖动点上色、
#' 字母画在每组 max(y) * letters_at 处、y 轴从 0 起、刻度朝内、**不标样本量**。
#' 每一项都可用 `CFG$plot$box` 或本函数的具名参数覆盖（换论文时按需改）。
md_box_group <- function(df, y, x, fill = NULL, letters = NULL,
                         colors = NULL, ylab = NULL, xlab = NULL,
                         y_limits = NULL, y_expand = NULL,
                         x_angle = 0, show_n = NULL, tag = NULL,
                         point_size = NULL, pt = MD_PT, legend_pos = "none",
                         color_by = NULL, kw_text = NULL, kw_annot = NULL,
                         y_from_zero = NULL, inward_ticks = NULL,
                         letters_at = NULL, letters_face = NULL,
                         letters_colour = NULL) {
  b <- md_box_cfg(y_expand = y_expand, show_n = show_n, color_by = color_by,
                  kw_annot = kw_annot, y_from_zero = y_from_zero,
                  inward_ticks = inward_ticks, letters_at = letters_at,
                  letters_face = letters_face, letters_colour = letters_colour)
  point_size <- point_size %||% MD_PT_SIZE$jitter

  lv_want <- names(colors) %||% unique(as.character(df[[x]]))
  d <- df[!is.na(df[[y]]) & !is.na(df[[x]]), , drop = FALSE]
  d[[x]] <- droplevels(factor(as.character(d[[x]]), levels = intersect(lv_want, unique(as.character(d[[x]])))))
  lv <- levels(d[[x]])
  if (length(lv) == 0) return(md_empty_plot(NULL, "no data"))

  # ---- 着色因子：nifH 用固定因子给箱体边框与点上色 ----
  cb <- b$color_by
  if (identical(cb, "auto")) cb <- if (length(lv) >= 2) "f1" else "none"
  if (identical(cb, "none") || length(lv) < 2) {
    cb <- "none"
    if (is.null(colors)) colors <- stats::setNames(rep("grey25", length(lv)), lv)
  } else if (is.null(colors)) {
    colors <- md_pal(NULL, lv, name = x)
  }

  # ---- y 轴范围 ----
  yr <- range(d[[y]], na.rm = TRUE)
  base <- if (!is.null(y_limits)) y_limits else
    if (isTRUE(b$y_from_zero)) c(0, max(yr[2], 0) * b$y_expand)
    else c(yr[1] - diff(yr) * 0.05, yr[2] + diff(yr) * (b$y_expand - 1))

  aes_main <- if (identical(cb, "none"))
    ggplot2::aes(x = .data[[x]], y = .data[[y]]) else
    ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[x]])

  g <- ggplot2::ggplot(d, aes_main) +
    ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = MD_LW$box,
                          fill = b$fill) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = MD_PT_SIZE$jitter_width, seed = MD_JITTER_SEED), size = point_size,
                         alpha = MD_PT_SIZE$alpha, stroke = MD_PT_SIZE$stroke) +
    ggplot2::scale_y_continuous(limits = base,
                                expand = ggplot2::expansion(
                                  mult = c(b$y_mult %||% 0.05, b$y_mult %||% 0.05))) +
    ggplot2::labs(y = ylab %||% y, x = xlab, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = legend_pos,
                   axis.text.x = ggplot2::element_text(
                     angle = x_angle, hjust = if (x_angle > 0) 1 else 0.5, size = pt$tick),
                   axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
                   axis.title = ggplot2::element_text(size = pt$title))
  if (!identical(cb, "none")) {
    g <- g + ggplot2::scale_colour_manual(values = colors, name = if (identical(legend_pos, "none")) NULL else x)
  }
  if (isTRUE(b$inward_ticks)) {
    g <- g + ggplot2::theme(axis.ticks.length = grid::unit(-0.12, "cm"))
  }

  # ---- CLD 字母（nifH 口径：每组 max(y) * letters_at；CFG$plot$box$letters = FALSE 时不画）----
  if (isTRUE(b$letters) && !is.null(letters) && nrow(letters) > 0) {
    ld <- local_letters_df(letters, d, y, x, lv, b)
    if (!is.null(ld) && nrow(ld) > 0) {
      g <- g + ggplot2::geom_text(data = ld,
                                  ggplot2::aes(x = .data$g, y = .data$y, label = .data$lab),
                                  inherit.aes = FALSE, vjust = -0.5,
                                  size = pt2gg(pt$ann),
                                  fontface = b$letters_face, colour = b$letters_colour)
    }
  }

  # ---- 样本量（默认关；开了才画）----
  if (isTRUE(b$show_n)) {
    n_df <- as.data.frame(table(d[[x]]), stringsAsFactors = FALSE)
    names(n_df) <- c("g", "n")
    n_df$y <- base[2] - diff(base) * 0.04
    g <- g + ggplot2::geom_text(data = n_df,
                                ggplot2::aes(x = .data$g, y = .data$y, label = .data$n),
                                inherit.aes = FALSE, size = pt2gg(pt$ann * 0.9),
                                colour = "grey35")
  }

  # ---- KW 注释（右上角，锚在面板内）----
  if (isTRUE(b$kw_annot) && !is.null(kw_text) && nzchar(kw_text)) {
    g <- g + ggplot2::annotate("text", x = Inf, y = Inf, label = kw_text,
                               hjust = 1.1, vjust = 1.5,
                               size = pt2gg(pt$ann), colour = "grey30")
  }
  g
}

#' 把 CLD 字母表对齐到图上的 x 水平，并算好 y 位置。
#' **字母分组名与 x 水平对不上时会 md_warn 报出来并丢弃不匹配项**（不再静默丢字母）。
local_letters_df <- function(letters, d, y, x, lv, b) {
  if (!all(c("group", "letter") %in% names(letters))) {
    stop("letters 表必须含 `group` 与 `letter` 两列（现在有: ",
         paste(names(letters), collapse = ", "), "）", call. = FALSE)
  }
  g_all <- as.character(letters$group)
  miss <- setdiff(g_all, lv)
  if (length(miss)) {
    md_warn(NULL, "箱线图字母与 x 水平不匹配，已丢弃 ", length(miss), " 个: ",
            paste(utils::head(miss, 5), collapse = ", "),
            if (length(miss) > 5) " …" else "",
            " | 检查生成 letters 时用的 inter_sep 是否与画图时一致")
  }
  keep <- intersect(g_all, lv)
  if (length(keep) == 0) return(NULL)
  max_y <- tapply(d[[y]], as.character(d[[x]]), max, na.rm = TRUE)
  data.frame(g = keep,
             lab = as.character(letters$letter)[match(keep, g_all)],
             y = as.numeric(max_y[keep]) * b$letters_at,
             stringsAsFactors = FALSE)
}

#' 双因子箱线：x 轴 = 两个因子的组合（顺序由 levels 决定）
#'
#' 默认样式**完全复刻 nifH 的 plot_dual_factor_panel**：
#'   · 空心箱体（`fill = NA`）+ 按 `color_by` 给**箱体边框与抖动点**上色（nifH 用第二因子）
#'   · 字母画在每组 `max(y) * letters_at`（默认 1.08）处，加粗、grey20
#'   · y 轴从 0 起、上限 = max * 1.35；刻度朝内
#'   · **不标样本量**（nifH 没有这个标注）
#'   · 右上角 KW 注释 "Zone: p = 0.033 *; Vegetation: p = 0.486 ns"
#' 换论文时逐项可改：`CFG$plot$box` 或本函数的具名参数。
#'
#' ⚠ `inter_sep` 必须与生成 `letters` 时用的分隔符一致，否则字母会被丢弃
#'   （local_letters_df() 会 md_warn 报出来，不再静默）。
md_box_two_factor <- function(df, y, f1, f2, letters = NULL, colors = NULL,
                              f1_levels = NULL, f2_levels = NULL,
                              ylab = NULL, y_limits = NULL, y_expand = NULL,
                              x_angle = 45, show_n = NULL, tag = NULL,
                              pt = MD_PT, legend_pos = "none",
                              x_labels = NULL, x_labels_one_line = NULL,
                              point_size = NULL,
                              color_by = NULL, kw_text = NULL, kw_annot = NULL,
                              y_from_zero = NULL, inward_ticks = NULL,
                              inter_sep = NULL, letters_at = NULL,
                              letters_face = NULL, letters_colour = NULL) {
  if (!is.null(x_labels_one_line))
    x_labels <- if (isTRUE(x_labels_one_line)) "one_line" else "two_line"
  b <- md_box_cfg(y_expand = y_expand, show_n = show_n, color_by = color_by,
                  kw_annot = kw_annot, y_from_zero = y_from_zero,
                  inward_ticks = inward_ticks, inter_sep = inter_sep,
                  letters_at = letters_at, letters_face = letters_face,
                  letters_colour = letters_colour, x_labels = x_labels)
  point_size <- point_size %||% MD_PT_SIZE$jitter
  sep <- b$inter_sep %||% "_"

  d <- df
  d[[f1]] <- factor(as.character(d[[f1]]), levels = f1_levels %||% unique(as.character(d[[f1]])))
  d[[f2]] <- factor(as.character(d[[f2]]), levels = f2_levels %||% unique(as.character(d[[f2]])))
  join <- function(a, bb) paste(a, bb, sep = sep)
  d$.inter <- factor(join(d[[f1]], d[[f2]]),
                     levels = as.vector(t(outer(levels(d[[f1]]), levels(d[[f2]]), join))))
  d <- d[!is.na(d$.inter) & !is.na(d[[y]]), , drop = FALSE]
  if (nrow(d) == 0) return(md_empty_plot(NULL, "no data"))

  # ---- 着色因子（默认 = 第二因子，与 nifH 一致）----
  cb <- b$color_by
  n1 <- nlevels(d[[f1]]); n2 <- nlevels(d[[f2]])
  if (identical(cb, "auto")) cb <- if (n2 >= 2) "f2" else if (n1 >= 2) "f1" else "none"
  if (identical(cb, "f2") && n2 < 2) {
    md_warn(NULL, "md_box_two_factor(): f2 只有 ", n2, " 个水平，着色退到 '",
            if (n1 >= 2) "f1" else "none", "'"); cb <- if (n1 >= 2) "f1" else "none"
  }
  if (identical(cb, "f1") && n1 < 2) {
    md_warn(NULL, "md_box_two_factor(): f1 只有 ", n1, " 个水平，着色退到 'none'"); cb <- "none"
  }
  if (!cb %in% c("f1", "f2")) cb <- "none"
  colvar <- switch(cb, f2 = f2, f1 = f1, NULL)
  if (is.null(colvar)) {
    if (is.null(colors)) colors <- "grey25"
  } else if (is.null(colors)) {
    colors <- md_pal(NULL, levels(d[[colvar]]), name = colvar)
  }

  # ---- y 轴范围 ----
  yr <- range(d[[y]], na.rm = TRUE)
  base <- if (!is.null(y_limits)) y_limits else
    if (isTRUE(b$y_from_zero)) c(0, max(yr[2], 0) * b$y_expand)
    else c(yr[1] - diff(yr) * 0.05, yr[2] + diff(yr) * (b$y_expand - 1))

  aes_main <- if (is.null(colvar))
    ggplot2::aes(x = .data$.inter, y = .data[[y]]) else
    ggplot2::aes(x = .data$.inter, y = .data[[y]], colour = .data[[colvar]])

  g <- ggplot2::ggplot(d, aes_main) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = MD_LW$box,
                          fill = b$fill) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = MD_PT_SIZE$jitter_width, seed = MD_JITTER_SEED), size = point_size,
                         alpha = MD_PT_SIZE$alpha, stroke = MD_PT_SIZE$stroke) +
    ggplot2::scale_y_continuous(limits = base,
                                expand = ggplot2::expansion(
                                  mult = c(b$y_mult %||% 0.05, b$y_mult %||% 0.05))) +
    ggplot2::labs(y = ylab %||% y, x = NULL, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(
      legend.position = legend_pos,
      axis.text.x = ggplot2::element_text(
        angle = x_angle, hjust = 1, size = pt$tick,
        lineheight = if (identical(b$x_labels, "one_line")) 1 else 0.85),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
      axis.title = ggplot2::element_text(size = pt$title))
  if (!is.null(colvar)) {
    g <- g + ggplot2::scale_colour_manual(
      values = colors, name = if (identical(legend_pos, "none")) NULL else colvar)
  }
  lab_fun <- switch(b$x_labels,
                    one_line = function(v) gsub(sep, " ", v, fixed = TRUE),
                    two_line = function(v) gsub(sep, "\n", v, fixed = TRUE),
                    function(v) v)
  g <- g + ggplot2::scale_x_discrete(labels = lab_fun)
  if (isTRUE(b$inward_ticks)) {
    g <- g + ggplot2::theme(axis.ticks.length = grid::unit(-0.12, "cm"),
                            axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)))
  }

  # ---- CLD 字母（CFG$plot$box$letters = FALSE 时不画）----
  if (isTRUE(b$letters) && !is.null(letters) && nrow(letters) > 0) {
    ld <- local_letters_df(letters, d, y, ".inter", levels(d$.inter), b)
    if (!is.null(ld) && nrow(ld) > 0) {
      g <- g + ggplot2::geom_text(data = ld,
                                  ggplot2::aes(x = .data$g, y = .data$y, label = .data$lab),
                                  inherit.aes = FALSE, vjust = -0.5,
                                  size = pt2gg(pt$ann),
                                  fontface = b$letters_face, colour = b$letters_colour)
    }
  }

  # ---- 样本量（默认关）----
  if (isTRUE(b$show_n)) {
    n_df <- as.data.frame(table(d$.inter), stringsAsFactors = FALSE)
    names(n_df) <- c("g", "n")
    n_df$y <- base[2] - diff(base) * 0.04
    g <- g + ggplot2::geom_text(data = n_df,
                                ggplot2::aes(x = .data$g, y = .data$y, label = .data$n),
                                inherit.aes = FALSE, size = pt2gg(pt$ann * 0.9),
                                colour = "grey35")
  }

  # ---- KW 注释 ----
  if (isTRUE(b$kw_annot) && !is.null(kw_text) && nzchar(kw_text)) {
    g <- g + ggplot2::annotate("text", x = Inf, y = Inf, label = kw_text,
                               hjust = 1.1, vjust = 1.5,
                               size = pt2gg(pt$ann), colour = "grey30")
  }
  g
}

#' 在图上加一段 P 值注解
#'
#' **双模式**（兼容两种写法，避免踩坑）：
#'   · `md_add_p_annot(p, labels = "…")`  -> 返回 `p + annotate(...)`（传图对象）
#'   · `p + md_add_p_annot(labels = "…")` -> 返回 annotation layer（当图层工厂用）
md_add_p_annot <- function(g = NULL, labels, x = -Inf, y = Inf, hjust = -0.05,
                           vjust = 1.4, pt = MD_PT, sep = "\n",
                           colour = "grey30", size = NULL, lineheight = NULL) {
  layer <- ggplot2::annotate("text", x = x, y = y, label = paste(labels, collapse = sep),
                             hjust = hjust, vjust = vjust,
                             size = size %||% pt2gg(pt$ann), colour = colour,
                             lineheight = lineheight %||% 0.9)
  if (inherits(g, "gg") || inherits(g, "patchwork")) return(g + layer)
  layer
}

# ---------------------------------------------------------------------------
# 散点 + 回归/平滑带
# ---------------------------------------------------------------------------
md_scatter_fit <- function(df, x, y, color = NULL, colors = NULL,
                           line = c("lm", "gam", "none")[1], ribbon = TRUE,
                           point_size = MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha,
                           xlab = NULL, ylab = NULL, tag = NULL,
                           pt = MD_PT, se = TRUE, method_args = list()) {
  d <- df[stats::complete.cases(df[, c(x, y)]), , drop = FALSE]
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(color)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = .data[[color]]),
                                 size = point_size, alpha = alpha)
    if (!is.null(colors)) g <- g + ggplot2::scale_colour_manual(values = colors)
  } else {
    g <- g + ggplot2::geom_point(size = point_size, alpha = alpha, colour = "grey25")
  }
  if (line == "lm") {
    g <- g + ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = ribbon,
                                  linewidth = MD_REG$lw, colour = MD_REG$col,
                                  fill = MD_REG$col, alpha = MD_REG$ribbon_alpha)
  } else if (line == "gam") {
    g <- g + ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = ribbon,
                                  linewidth = MD_REG$lw, colour = MD_REG$col,
                                  fill = MD_REG$col, alpha = MD_REG$ribbon_alpha)
  }
  g + ggplot2::labs(x = xlab %||% x, y = ylab %||% y, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = if (is.null(color)) "none" else "bottom",
                   axis.title = ggplot2::element_text(size = pt$title))
}

#' 条件效应曲线（已算好预测值与置信带，模块只负责把数据整理好）
md_effect_curve <- function(df, x, fit, lo = NULL, hi = NULL, points = NULL,
                            y = NULL, ribbon = TRUE, xlab = NULL, ylab = NULL, tag = NULL,
                            pt = MD_PT, x_breaks = NULL) {
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[fit]]))
  if (ribbon && !is.null(lo) && !is.null(hi)) {
    g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]),
                                  fill = "#1E72B0", alpha = 0.12, colour = NA)
  }
  g <- g + ggplot2::geom_line(linewidth = MD_LW$line, colour = "#1E72B0")
  if (!is.null(points) && !is.null(y)) {
    g <- g + ggplot2::geom_point(data = points, ggplot2::aes(x = .data[[x]], y = .data[[y]]),
                                 inherit.aes = FALSE, size = MD_PT_SIZE$point * 0.8,
                                 alpha = MD_PT_SIZE$alpha, colour = "grey25")
  }
  if (!is.null(x_breaks)) g <- g + ggplot2::scale_x_continuous(breaks = x_breaks)
  g + ggplot2::labs(x = xlab %||% x, y = ylab %||% fit, tag = tag) +
    theme_md(CFG = NULL)
}

# ---------------------------------------------------------------------------
# 热图 / 柱状 / 森林 / 直方 / 四象限
# ---------------------------------------------------------------------------
md_heatmap <- function(df, x, y, fill, values = NULL, label = FALSE,
                       xlab = NULL, ylab = NULL, tag = NULL,
                       mid = NULL, pt = MD_PT, label_digits = 2,
                       x_angle = 45, border = "white") {
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    ggplot2::geom_tile(colour = border, linewidth = 0.25)
  if (!is.null(values)) {
    g <- g + ggplot2::scale_fill_gradient2(low = values[1], mid = values[2], high = values[3],
                                           midpoint = mid %||% 0)
  }
  if (label) {
    g <- g + ggplot2::geom_text(ggplot2::aes(label = formatC(.data[[fill]], format = "f",
                                                             digits = label_digits)),
                                size = pt2gg(pt$ann * 0.85))
  }
  g + ggplot2::labs(x = xlab, y = ylab, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1, size = pt$tick),
                   axis.ticks = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
}

#' 堆叠柱（组成图 / 组装过程比例等）。x 可为因子或连续。
md_bar_stack <- function(df, x, y, fill, colors = NULL, xlab = NULL, ylab = NULL,
                         tag = NULL, width = 0.65, pt = MD_PT, x_angle = 45,
                         legend_pos = "bottom", ncol = NULL,
                         outline_colour = NA, outline_lw = 0.15, drop = TRUE) {
  # outline_colour/lw：nifH Fig3B 用 `colour="white", linewidth=0.2` 做柱间白描边
  # drop：FALSE 时保留未出现的因子水平（原脚本的 scale_fill_manual(drop=FALSE)）
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    ggplot2::geom_col(width = width, colour = outline_colour, linewidth = outline_lw)
  if (!is.null(colors)) {
    g <- g + ggplot2::scale_fill_manual(values = colors, drop = drop,
                                        guide = md_guide_flex(length(colors), legend_pos, ncol %||% 4))
  }
  g + ggplot2::labs(x = xlab, y = ylab, tag = tag, fill = NULL) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = legend_pos,
                   legend.box = "vertical",
                   axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1, size = pt$tick))
}

#' 森林图：estimate + CI + P
#'
#' 点样式支持 nifH `make_coef_panel` 的「方向 × 显著性」三色：
#'  - `point_col`：df 里**逐行字面颜色**的列名（配合 `scale_colour_identity()`）
#'  - 或 `color_by` + `colors`：按某列映射（整列一套色板）
#'  - `point_shape` / `point_stroke` / `fill_by` + `fill_vals`：画 shape 21 空心点时按列给填充色
md_forest <- function(df, term, estimate, lo = NULL, hi = NULL, p = NULL,
                      ref = 0, xlab = NULL, tag = NULL, pt = MD_PT,
                      colors = NULL, color_by = NULL, label_p = TRUE,
                      point_col = NULL, point_shape = 19, point_stroke = 0.3,
                      point_size = 1.8, fill_by = NULL, fill_vals = NULL,
                      label_p_ascii = FALSE) {
  d <- df
  d$.term <- factor(as.character(d[[term]]), levels = rev(unique(as.character(d[[term]]))))
  g <- ggplot2::ggplot(d, ggplot2::aes(y = .data$.term, x = .data[[estimate]]))
  if (!is.null(lo) && !is.null(hi)) {
    g <- g + ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data[[lo]], xmax = .data[[hi]]),
      orientation = "y", width = 0.18, linewidth = MD_LW$line * 0.8,
      colour = "grey30")
  }
  if (!is.null(point_col)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = .data[[point_col]]),
                                 shape = point_shape, stroke = point_stroke,
                                 size = point_size) +
      ggplot2::scale_colour_identity()
  } else if (!is.null(color_by) && !is.null(colors)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = .data[[color_by]]),
                                 shape = point_shape, stroke = point_stroke,
                                 size = point_size)
    g <- g + ggplot2::scale_colour_manual(values = colors)
  } else {
    g <- g + ggplot2::geom_point(shape = point_shape, stroke = point_stroke,
                                 size = point_size, colour = "grey15")
  }
  if (!is.null(fill_by) && !is.null(fill_vals)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(fill = .data[[fill_by]]),
                                 shape = 21, stroke = point_stroke, size = point_size) +
      ggplot2::scale_fill_manual(values = fill_vals)
  }
  g <- g + ggplot2::geom_vline(xintercept = ref, linetype = "dashed",
                               linewidth = MD_LW$ref, colour = "grey55")
  if (label_p && !is.null(p)) {
    lab <- if (isTRUE(label_p_ascii)) md_fmt_p_fig(d[[p]]) else md_fmt_p(d[[p]])
    g <- g + ggplot2::geom_text(ggplot2::aes(label = lab),
                                hjust = -0.15, size = pt2gg(pt$ann * 0.85), colour = "grey30")
  }
  g + ggplot2::labs(x = xlab %||% "Estimate", y = NULL, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = pt$tick),
                   panel.grid.major.y = ggplot2::element_line(linewidth = 0.2, colour = "grey93"))
}

#' 直方图 + 双侧参考线标注（标注锚在面板内部，避免被裁）
#'
#' 原图（nifH 09_Supp_assembly.R §5 / 11_shared.R plot_bnti_hist）需要的样式都已暴露为参数：
#' `alpha`（原 0.7）、`outline_colour/lw`、`ref_col/ref_lw/ref_type`（原 red / 0.8）、
#' `label_vjust`（原 1.15）、`x_expand`（RCbray 原 expansion(mult = 0.22)）。
md_hist <- function(df, x, bins = 40, ref = NULL, fill = "steelblue",
                    left_label = NULL, right_label = NULL,
                    left_col = "#2166AC", right_col = "#E31A1C",
                    xlab = NULL, ylab = "Frequency", tag = NULL, pt = MD_PT,
                    alpha = 0.85, outline_colour = "white", outline_lw = 0.15,
                    ref_col = "grey45", ref_lw = MD_LW$ref, ref_type = "dashed",
                    label_vjust = 1.3, x_expand = NULL) {
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]])) +
    ggplot2::geom_histogram(bins = bins, fill = fill, colour = outline_colour,
                            linewidth = outline_lw, alpha = alpha)
  if (!is.null(ref)) {
    g <- g + ggplot2::geom_vline(xintercept = ref, linetype = ref_type,
                                 linewidth = ref_lw, colour = ref_col)
  }
  if (!is.null(x_expand)) {
    g <- g + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = x_expand))
  }
  if (!is.null(left_label)) {
    g <- g + ggplot2::annotate("text", x = -Inf, y = Inf, label = left_label,
                               hjust = -0.08, vjust = label_vjust,
                               size = pt2gg(pt$ann), colour = left_col)
  }
  if (!is.null(right_label)) {
    g <- g + ggplot2::annotate("text", x = Inf, y = Inf, label = right_label,
                               hjust = 1.08, vjust = label_vjust,
                               size = pt2gg(pt$ann), colour = right_col)
  }
  g + ggplot2::labs(x = xlab %||% x, y = ylab, tag = tag) + theme_md(CFG = NULL)
}

#' 四象限散点
md_quadrant <- function(df, x, y, hline = 0, vline = 0, color = NULL, colors = NULL,
                        xlab = NULL, ylab = NULL, tag = NULL,
                        labels = NULL, pt = MD_PT, point_size = MD_PT_SIZE$point) {
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  g <- g + ggplot2::geom_hline(yintercept = hline, linetype = "dashed",
                               linewidth = MD_LW$ref, colour = "grey55") +
    ggplot2::geom_vline(xintercept = vline, linetype = "dashed",
                        linewidth = MD_LW$ref, colour = "grey55")
  if (!is.null(color)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = .data[[color]]), size = point_size, alpha = 0.7)
    if (!is.null(colors)) g <- g + ggplot2::scale_colour_manual(values = colors)
  } else {
    g <- g + ggplot2::geom_point(size = point_size, alpha = 0.6, colour = "grey25")
  }
  if (!is.null(labels)) {
    g <- g + ggplot2::annotate("text", x = -Inf, y = Inf, label = labels[1],
                               hjust = -0.1, vjust = 1.4, size = pt2gg(pt$ann), colour = labels[3]) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = labels[2],
                        hjust = 1.1, vjust = 1.4, size = pt2gg(pt$ann), colour = labels[4])
  }
  g + ggplot2::labs(x = xlab %||% x, y = ylab %||% y, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = if (is.null(color)) "none" else "bottom")
}

#' 三圆 Venn（变差分解用）—— **逐行照抄 nifH `03_Fig03_assembly.R:288-329` 的 `draw_venn3()`**
#'
#' 几何：圆半径 1.15、圆心落在半径 0.85 的三个方向（90°/210°/330°）、标签 offset 1.55；
#' 填充 alpha 0.3 + 灰边 0.7；7 个百分比标签（3 个独立 + 3 个两两 + 1 个三交）；
#' 块名用该块的颜色加粗写在圈外；`coord_fixed(clip="off")` + `theme_void()`。
#' **负交叠项按原脚本 `clamp()` 截为 0**（图上显示 0.0%）；负值明细请放在结果表的 Note 里。
md_venn3 <- function(lab1, lab2, lab3,
                     f1, f2, f3, f12, f13, f23, f123,
                     col1, col2, col3, title_txt = NULL, subtitle = NULL,
                     pt = MD_PT, tag = NULL) {
  clamp <- function(x) pmax(x, 0)
  r  <- 1.15
  cx1 <- cos(pi / 2) * 0.85;               cy1 <- sin(pi / 2) * 0.85
  cx2 <- cos(pi / 2 + 2 * pi / 3) * 0.85;  cy2 <- sin(pi / 2 + 2 * pi / 3) * 0.85
  cx3 <- cos(pi / 2 + 4 * pi / 3) * 0.85;  cy3 <- sin(pi / 2 + 4 * pi / 3) * 0.85
  offset <- 1.55
  labels <- data.frame(
    x   = c(cx1 * offset, cx2 * offset, cx3 * offset,
            (cx1 + cx2) * 0.55, (cx1 + cx3) * 0.55, (cx2 + cx3) * 0.55, 0),
    y   = c(cy1 * offset - 0.35, cy2 * offset - 0.35, cy3 * offset - 0.35,
            (cy1 + cy2) * 0.55, (cy1 + cy3) * 0.55, (cy2 + cy3) * 0.55, 0),
    txt = sprintf("%.1f%%", 100 * clamp(c(f1, f2, f3, f12, f13, f23, f123))),
    size = c(rep(pt2gg(pt$ann), 3), rep(pt2gg(pt$ann) * 0.9, 3), pt2gg(pt$ann) * 0.9)
  )
  ggplot2::ggplot() +
    ggforce::geom_circle(ggplot2::aes(x0 = cx1, y0 = cy1, r = r),
                         fill = col1, alpha = 0.3, colour = "grey40", linewidth = 0.7) +
    ggforce::geom_circle(ggplot2::aes(x0 = cx2, y0 = cy2, r = r),
                         fill = col2, alpha = 0.3, colour = "grey40", linewidth = 0.7) +
    ggforce::geom_circle(ggplot2::aes(x0 = cx3, y0 = cy3, r = r),
                         fill = col3, alpha = 0.3, colour = "grey40", linewidth = 0.7) +
    ggplot2::annotate("text", x = cx1 * 1.55, y = cy1 * 1.55, label = lab1,
                      size = pt2gg(pt$legend), fontface = "bold", colour = col1) +
    ggplot2::annotate("text", x = cx2 * 1.55, y = cy2 * 1.55, label = lab2,
                      size = pt2gg(pt$legend), fontface = "bold", colour = col2) +
    ggplot2::annotate("text", x = cx3 * 1.55, y = cy3 * 1.55, label = lab3,
                      size = pt2gg(pt$legend), fontface = "bold", colour = col3) +
    ggplot2::geom_text(data = labels,
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt,
                                    size = I(.data$size)),
                       colour = "grey20") +
    # ⚠ theme_void() 挡不住 md_compose() 之后套上的 theme_md_comp()：
    #   后者会把 axis.text / axis.title 重新设成可见，于是 Venn 旁边会长出 -2..2 的坐标轴。
    #   清空两个 scale 的 breaks 才是根治（没有 break 就没有刻度文字）。
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(title = title_txt, subtitle = subtitle, tag = tag, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = pt$title, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = pt$ann, colour = "grey30", hjust = 0.5),
      plot.tag      = ggplot2::element_text(size = pt$tag, face = "bold"),
      plot.margin   = ggplot2::margin(4, 4, 4, 4, unit = "pt")
    )
}

#' NRI / NTI 面板（±2 参考带 + 分组箱线 + CLD 字母 + KW 注释）
#' 样式与 md_box_group 同源，便于 Fig3A / FigS4 口径一致（nifH 的 plot_nri_nti）。
md_nri_nti_panel <- function(df, metric, group, order_levels = NULL,
                             colors = NULL, ref = c(-2, 2), ylab = NULL, tag = NULL,
                             x_angle = 45, x_one_line = NULL, pt = MD_PT,
                             letters = NULL, kw_text = NULL,
                             color_by = NULL, inward_ticks = NULL,
                             letters_at = NULL,
                             # 以下默认值逐项 = nifH 11_shared.R §4.2 plot_nri_nti()
                             ref_col = "red", ref_lw = MD_LW$ref, ref_alpha = 0.5,
                             y_expand = c(0.05, 0.15),
                             letters_at_ori = 1.18,
                             show_kw_test = FALSE) {
  b <- md_box_cfg(color_by = color_by, inward_ticks = inward_ticks,
                  letters_at = letters_at)
  d <- df[!is.na(df[[metric]]), , drop = FALSE]
  d[[group]] <- droplevels(factor(as.character(d[[group]]),
                                  levels = order_levels %||% unique(as.character(d[[group]]))))
  lv <- levels(d[[group]])
  if (length(lv) == 0) return(md_empty_plot(NULL, "no data"))

  cb <- b$color_by
  if (identical(cb, "auto")) cb <- if (length(lv) >= 2) "f1" else "none"
  if (identical(cb, "none") || length(lv) < 2) {
    cb <- "none"
    if (is.null(colors)) colors <- stats::setNames(rep("grey25", length(lv)), lv)
  } else if (is.null(colors)) {
    colors <- md_pal(NULL, lv, name = group)
  }

  aes_main <- if (identical(cb, "none"))
    ggplot2::aes(x = .data[[group]], y = .data[[metric]]) else
    ggplot2::aes(x = .data[[group]], y = .data[[metric]], colour = .data[[group]])

  kmax <- tryCatch(stats::kruskal.test(d[[metric]] ~ d[[group]])$p.value,
                   error = function(e) NA_real_)
  # 默认值逐项对齐 nifH 11_shared.R §4.2 plot_nri_nti()：
  #   ±2 参考线 = red / alpha 0.5 的虚线；jitter 用整份 JIT_SIZE_STD（不缩放）；
  #   字母在 max_y*1.18、vjust 0、黑色；y 轴顶部留 15% 余量（原文注释说明必需）；
  #   **不在左下角画 Kruskal-Wallis**（原图没有，由 show_kw_test 显式打开）
  g <- ggplot2::ggplot(d, aes_main) +
    ggplot2::geom_hline(yintercept = ref, linetype = "dashed",
                        linewidth = ref_lw, colour = ref_col, alpha = ref_alpha) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, fill = b$fill,
                          linewidth = MD_LW$box) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = MD_PT_SIZE$jitter_width, seed = MD_JITTER_SEED), size = MD_PT_SIZE$jitter,
                         alpha = MD_PT_SIZE$alpha, stroke = MD_PT_SIZE$stroke) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = y_expand)) +
    ggplot2::labs(y = ylab %||% metric, x = NULL, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1, size = pt$tick))
  if (!identical(cb, "none")) {
    g <- g + ggplot2::scale_colour_manual(values = colors)
  }
  if (isTRUE(x_one_line)) {
    g <- g + ggplot2::scale_x_discrete(labels = function(x) gsub("_", " ", x, fixed = TRUE))
  }
  if (isTRUE(b$inward_ticks)) {
    g <- g + ggplot2::theme(axis.ticks.length = grid::unit(-0.12, "cm"))
  }

  # ---- CLD 字母（原图口径：max_y * 1.18、vjust 0、黑色；letters 开关见 CFG$plot$box$letters）----
  if (isTRUE(b$letters) && !is.null(letters) && nrow(letters) > 0) {
    ld <- local_letters_df(letters, d, metric, group, lv,
                           modifyList(b, list(letters_at = letters_at_ori)))
    if (!is.null(ld) && nrow(ld) > 0) {
      g <- g + ggplot2::geom_text(data = ld,
                                  ggplot2::aes(x = .data$g, y = .data$y, label = .data$lab),
                                  inherit.aes = FALSE, vjust = 0,
                                  size = pt2gg(pt$ann),
                                  fontface = "bold", colour = "black")
    }
  }

  # ---- 左下角 Kruskal-Wallis：**默认不画**（原图没有），需要时 show_kw_test = TRUE ----
  if (isTRUE(show_kw_test) && !is.na(kmax)) {
    g <- g + ggplot2::annotate("text", x = -Inf, y = Inf,
                               label = paste0("Kruskal-Wallis ", md_fmt_p(kmax)),
                               hjust = -0.06, vjust = 1.5, size = pt2gg(pt$ann * 0.9))
  }
  # 右上角 KW 注释（" ; " 连接，原图 plot_nri_nti 的写法）
  if (isTRUE(b$kw_annot) && !is.null(kw_text) && nzchar(kw_text)) {
    g <- g + ggplot2::annotate("text", x = Inf, y = Inf, label = kw_text,
                               hjust = 1.05, vjust = 1.5,
                               size = pt2gg(pt$ann * 0.9), colour = "grey30")
  }
  g
}

# =============================================================================
# 2026 毕业论文补充轮新增图型
#   · md_lollipop()           指示种 IndVal 点棒图（10_composition 的 indicator 面板 / 81_taxon_group）
#   · md_rf_importance_plot() 随机森林重要性横向条形（40_env_summary 的 rf 面板）
#   · md_gam_curve_plot()     GAM/GAMM 响应曲线（45_response_curves）
#   · md_venn2()              两圆 Venn（82_taxon_share 在 2 水平时的退化情形）
# 每一条的版式都逐项照抄毕业论文原脚本，参数只暴露"样式"不暴露"统计"。
# =============================================================================

#' IndVal 指示种点棒图（横向 lollipop）
#'
#' 源：`S4_supplementary_taxa.R:222-269`（Fig08a / Fig08b）
#' 原文版式：`geom_segment`(到 0) + `geom_point`(1.8) + 显著性星号文本
#' （y 偏移 = `ifelse(IndVal > 0.85, -0.06, 0.04)`）+ `coord_flip()` +
#' 组内按 IndVal 降序 + 图例关掉（共用图例在拼图时收集）+ 标题左对齐加粗。
#'
#' @param df     需含 `label` / `indval` / `level` 三列（`sig` 可选）
#' @param colors 按 `level` 命名的颜色向量（模块用 `md_pal()` 取）
#' @param title  面板标题（ASCII）
md_lollipop <- function(df, label = "label", indval = "indval", level = "level",
                        sig = "sig", colors = NULL, title = NULL,
                        xlab = "Indicator Value (IndVal)", ylab = NULL,
                        pt = MD_PT, tag = NULL, text_size = 2,
                        label_size_pt = 4.3, point_size = 1.8) {
  d <- as.data.frame(df)
  d <- d[is.finite(d[[indval]]), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(NULL, "no indicator taxa"))
  lv <- unique(as.character(d[[level]]))
  d$.lv <- factor(as.character(d[[level]]), levels = lv)
  d <- d[order(d$.lv, -d[[indval]]), , drop = FALSE]
  d$.lab <- factor(as.character(d[[label]]), levels = rev(unique(as.character(d[[label]]))))
  # 原文：IndVal > 0.85 时星号画在点下方，否则上方
  d$.sig <- if (sig %in% names(d)) as.character(d[[sig]]) else ""
  d$.sig[is.na(d$.sig)] <- ""
  d$.ty <- d[[indval]] + ifelse(d[[indval]] > 0.85, -0.06, 0.04)

  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$.lab, y = .data[[indval]])) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$.lab, yend = 0,
                                       colour = .data$.lv),
                          linewidth = 0.4, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.lv), size = point_size) +
    ggplot2::geom_text(ggplot2::aes(y = .data$.ty, label = .data$.sig),
                       size = text_size, colour = "black", hjust = 0.5) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = ylab, y = xlab, colour = NULL, title = title, tag = tag)
  if (!is.null(colors)) {
    g <- g + ggplot2::scale_colour_manual(values = colors, drop = FALSE)
  }
  g + theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = label_size_pt),
                   plot.title = ggplot2::element_text(size = pt$title, face = "bold", hjust = 0))
}

#' 随机森林变量重要性（横向条形，相对重要性缩放到 max = 1）
#'
#' 源：`3_environmental_drivers.R:1232-1250`（`p_rf`）—— 单色单柱、无 dodged、
#' 无 model 图例；x 轴题 "Relative permutation importance"；y 轴题关掉。
md_rf_importance_plot <- function(df, variable = "Variable",
                                  importance = "Importance_scaled",
                                  top_n = 10L, fill = "#3C5488",
                                  xlab = "Relative permutation importance",
                                  ylab = NULL, pt = MD_PT, tag = NULL,
                                  r2 = NULL, label_r2 = TRUE) {
  d <- as.data.frame(df)
  if (!nrow(d)) return(md_empty_plot(NULL, "no importance"))
  d <- d[is.finite(d[[importance]]), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(NULL, "no importance"))
  d <- utils::head(d[order(-d[[importance]]), , drop = FALSE], as.integer(top_n))
  d$.v <- factor(as.character(d[[variable]]), levels = rev(as.character(d[[variable]])))
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[importance]], y = .data$.v)) +
    ggplot2::geom_col(fill = fill, width = 0.68) +
    ggplot2::labs(x = xlab, y = ylab, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = pt$tick),
                   panel.grid.major.y = ggplot2::element_blank())
  if (isTRUE(label_r2) && !is.null(r2) && is.finite(r2)) {
    g <- g + ggplot2::annotate("text", x = Inf, y = -Inf,
                               label = sprintf("R2 = %.3f", r2),
                               hjust = 1.1, vjust = -0.8,
                               size = pt2gg(pt$ann * 0.9), colour = "grey30")
  }
  g
}

#' GAM / GAMM 响应曲线（点 + 拟合线 + 95% 置信带）
#'
#' 源：`3_environmental_drivers.R:1329-1343` —— 灰色 95% ribbon（fill grey80,
#' alpha 0.5）、`#0072B2` 拟合线、灰点 alpha 0.3；R²/P 以**面板内注释**给出
#' （原注：标题写法在 183 mm 三列排布下会互相压字）。
md_gam_curve_plot <- function(curve, points = NULL, title = NULL, subtitle = NULL,
                              xlab = NULL, ylab = NULL, tag = NULL,
                              line_col = "#0072B2", ribbon_fill = "grey80",
                              ribbon_alpha = 0.5, point_col = "grey40",
                              point_size = 1, point_alpha = 0.3,
                              r2 = NULL, p = NULL, pt = MD_PT) {
  if (is.null(curve) || !nrow(curve)) return(md_empty_plot(NULL, "not fitted"))
  g <- ggplot2::ggplot()
  if (!is.null(points) && nrow(points)) {
    g <- g + ggplot2::geom_point(data = points,
                                 ggplot2::aes(x = .data$x, y = .data$y),
                                 colour = point_col, size = point_size, alpha = point_alpha)
  }
  g <- g +
    ggplot2::geom_ribbon(data = curve,
                         ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
                         fill = ribbon_fill, alpha = ribbon_alpha) +
    ggplot2::geom_line(data = curve, ggplot2::aes(x = .data$x, y = .data$fit),
                       colour = line_col, linewidth = 1)
  if (!is.null(r2) || !is.null(p)) {
    lab <- paste(c(if (!is.null(r2)) sprintf("R2=%.3f", r2),
                   if (!is.null(p)) sprintf("P=%.4f", p)), collapse = "\n")
    g <- g + ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.12, vjust = 1.6,
                               label = lab, size = pt2gg(pt$ann * 0.85),
                               colour = "grey25", lineheight = 0.95)
  }
  g + ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = pt$title))
}

#' 两圆 Venn（82_taxon_share 在因子只有 2 水平时用）
#'
#' 几何口径沿用 `md_venn3()`（圆半径 1.15、alpha 0.3、grey40 边 0.7、
#' `coord_fixed(clip="off")` + `theme_void()`），只是圆数由 3 减到 2。
#' 覆盖度百分比与 `md_venn3` 同样的 5 位（A、B、A∩B、A\B、B\A）。
md_venn2 <- function(lab1, lab2, f1, f2, f12,
                     col1 = "#5B93C4", col2 = "#E0877C",
                     title_txt = NULL, subtitle = NULL, pt = MD_PT, tag = NULL) {
  r    <- 1.15
  dsep <- 0.72
  cx1  <- -dsep / 2; cx2 <- dsep / 2
  pct  <- function(x) sprintf("%.1f%%", 100 * pmax(x, 0))
  only1 <- max(f1 - f12, 0); only2 <- max(f2 - f12, 0)
  labels <- data.frame(
    x = c(cx1 - 0.62, cx2 + 0.62, 0, cx1 - 0.62, cx2 + 0.62),
    y = c(0, 0, 0, -0.62, -0.62),
    txt = c(pct(only1), pct(only2), pct(f12), "", ""),
    stringsAsFactors = FALSE)
  names_lab <- data.frame(
    x = c(cx1, cx2), y = c(r + 0.35, r + 0.35),
    txt = c(lab1, lab2), col = c(col1, col2), stringsAsFactors = FALSE)
  ggplot2::ggplot() +
    ggforce::geom_circle(ggplot2::aes(x0 = cx1, y0 = 0, r = r),
                         fill = col1, alpha = 0.3, colour = "grey40", linewidth = 0.7) +
    ggforce::geom_circle(ggplot2::aes(x0 = cx2, y0 = 0, r = r),
                         fill = col2, alpha = 0.3, colour = "grey40", linewidth = 0.7) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt),
                       size = pt2gg(pt$ann), colour = "grey20") +
    ggplot2::geom_text(data = names_lab,
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt),
                       size = pt2gg(pt$title), colour = names_lab$col, fontface = "bold") +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(title = title_txt, subtitle = subtitle, tag = tag) +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = pt$title),
                   plot.subtitle = ggplot2::element_text(size = pt$label, colour = "grey35"),
                   plot.margin = ggplot2::margin(6, 6, 6, 6))
}
