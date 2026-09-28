# =============================================================================
# R/14_ordination.R —— 模块 14：排序（PCoA / NMDS × 多距离度量，含 95% 椭圆）
# -----------------------------------------------------------------------------
# 导出：md_ordination(CFG)
#
# 输入：**必须通过** cache/13_dist.rds（md_dist_cache() / md_dist_get()），
#       即本模块依赖模块 13；单独跑 14 会因缺缓存报错，请与 13 一起跑。
#
# 方法（CFG$ordination$method）：pcoa | nmds | both
#   PCoA = ape::pcoa(d, correction = "cailliez")，轴解释率 = 特征值 / 正特征值之和
#   NMDS = vegan::metaMDS(d, trymax = CFG$ordination$nmds_trymax)
# 椭圆（CFG$ordination$ellipse）：ggplot2::stat_ellipse(type = "norm", level = CFG$ordination$ellipse_level)
#   每组样本数 < CFG$ordination$ellipse_min_n 或组内轴方差 ~ 0 时不画，并在表里记 NA + 原因
# 分组着色 = md_pal(CFG, levels, name)，name 用因子的列名（如 "Zone" / "Vegetation"）
#
# 出图：figures/14_ordination.pdf（每个距离度量一个面板，拼成一张图）
# 出表：results/14_ordination.xlsx（轴解释率 / 椭圆诊断 / 站点得分）
# 注：PERMDISP 不在此模块（属于模块 15）。
# =============================================================================

# ---------------------------------------------------------------------------
# 模块内小工具（local_ 前缀）
# ---------------------------------------------------------------------------
local_factor_levels <- function(CFG, meta, var) {
  ol <- CFG$design$ordered_levels[[var]]
  obs <- unique(as.character(meta[[var]]))
  obs <- obs[!is.na(obs)]
  if (!is.null(ol)) return(intersect(as.character(ol), obs))
  obs
}

#' 形状值：照抄 02_Fig02.R §B/C 的
#'   scale_shape_manual(values = c(Bare = 16, Herb = 17, Shrub = 15, Tree = 18))
#' 原脚本没有的水平按同一族的填充形状补位（本包示范数据的 Vegetation 只有 Bare/Herb/Tree）。
local_shape_values <- function(lv) {
  builtin <- c(Bare = 16, Herb = 17, Shrub = 15, Tree = 18)
  out <- unname(builtin[lv])
  if (anyNA(out)) {
    pool <- setdiff(c(16, 17, 15, 18, 8, 3, 7, 4, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14), out)
    out[is.na(out)] <- pool[seq_len(sum(is.na(out)))]
  }
  out
}

#' 取轴解释率（%）：优先用 Cailliez 校正特征值，退化时回退原始特征值
local_eig_pct <- function(pc, n_axis = 2) {
  v <- pc$values
  ev <- v$Corr_eig
  if (is.null(ev) || all(!is.finite(ev))) ev <- v$Eigenvalues
  ev <- as.numeric(ev)
  ev[!is.finite(ev)] <- 0
  denom <- sum(ev[ev > 0])
  if (!is.finite(denom) || denom <= 0) return(rep(NA_real_, n_axis))
  pct <- 100 * ev[seq_len(n_axis)] / denom
  as.numeric(pct)
}

#' 可以画椭圆的组（n >= min_n 且两轴方差都非 0）；返回 list(keep, dropped)
local_ellipse_groups <- function(df, group, min_n, sd_tol = 1e-6) {
  if (is.null(group) || !group %in% names(df)) return(list(keep = character(0), dropped = NULL))
  g <- as.character(df[[group]])
  lv <- unique(g[!is.na(g)])
  keep <- character(0); drop <- list()
  for (x in lv) {
    d <- df[!is.na(g) & g == x, , drop = FALSE]
    n <- nrow(d)
    s1 <- if (n > 1) stats::sd(d$Axis1) else 0
    s2 <- if (n > 1) stats::sd(d$Axis2) else 0
    if (n < min_n) {
      drop[[length(drop) + 1L]] <- data.frame(group = x, n = n, sd_axis1 = s1, sd_axis2 = s2,
                                              drawn = FALSE, reason = paste0("n < ellipse_min_n (", min_n, ")"),
                                              stringsAsFactors = FALSE)
    } else if (!is.finite(s1) || !is.finite(s2) || s1 < sd_tol || s2 < sd_tol) {
      drop[[length(drop) + 1L]] <- data.frame(group = x, n = n, sd_axis1 = s1, sd_axis2 = s2,
                                              drawn = FALSE, reason = "zero variance on an axis",
                                              stringsAsFactors = FALSE)
    } else {
      keep <- c(keep, x)
      drop[[length(drop) + 1L]] <- data.frame(group = x, n = n, sd_axis1 = s1, sd_axis2 = s2,
                                              drawn = TRUE, reason = "", stringsAsFactors = FALSE)
    }
  }
  list(keep = keep, dropped = if (length(drop)) do.call(rbind, drop) else NULL)
}

#' 单个排序面板 —— 画法逐行照抄 02_Fig02.R §B/C（fig2A / fig2B）：
#'   · 椭圆：stat_ellipse(geom = "polygon", alpha = 0.10 * 4/3, level = CFG 值,
#'            type = "norm", colour = NA, show.legend = FALSE)  —— 无边框、只填半透明多边形
#'   · 散点：geom_point(size = 2.2 * 2/3, alpha = 0.75, stroke = 0.5)，
#'            第一因子着色（原图 Zone）、第二因子给形状（原图 Vegetation）
#'   · 零轴参考线：geom_hline/vline(linetype = "dashed", colour = "grey70")
#'   · 轴标题 "PCoA1 (x%)" / "PCoA2 (y%)"，x/y 都用 round(特征值占比, 2)
#' 颜色走 md_pal()，字号/线宽走 MD_PT / MD_LW。
local_ord_panel <- function(df, group, pal, xlab, ylab, tag, ellipse_level, min_n,
                            method = "pcoa", stress = NA_real_,
                            shape_var = NULL, subtitle = NULL) {
  eg <- local_ellipse_groups(df, group, min_n)
  d <- df
  d$.g <- factor(as.character(d[[group]]), levels = names(pal))
  d <- d[!is.na(d$.g), , drop = FALSE]
  d_ell <- d[as.character(d$.g) %in% eg$keep, , drop = FALSE]

  use_shape <- !is.null(shape_var) && shape_var %in% names(d) &&
    !identical(shape_var, group) &&
    length(unique(stats::na.omit(as.character(d[[shape_var]])))) >= 2
  if (use_shape) {
    d$.s <- droplevels(factor(as.character(d[[shape_var]])))
    shape_vals <- local_shape_values(levels(d$.s))
  }

  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Axis1, y = .data$Axis2, colour = .data$.g))
  if (nrow(d_ell) > 0) {
    g <- g + ggplot2::stat_ellipse(
      data = d_ell, ggplot2::aes(group = .data$.g, fill = .data$.g),
      geom = "polygon", alpha = 0.10 * 4 / 3, level = ellipse_level,
      type = "norm", colour = NA, show.legend = FALSE)
  }
  g <- g +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70",
                        linewidth = MD_LW$ref) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70",
                        linewidth = MD_LW$ref)
  if (use_shape) {
    g <- g + ggplot2::geom_point(ggplot2::aes(shape = .data$.s),
                                 size = 2.2 * 2 / 3, alpha = 0.75, stroke = 0.5)
  } else {
    g <- g + ggplot2::geom_point(size = 2.2 * 2 / 3, alpha = 0.75, stroke = 0.5)
  }
  g <- g +
    ggplot2::scale_colour_manual(values = pal, name = group) +
    ggplot2::scale_fill_manual(values = pal, name = group) +
    ggplot2::labs(x = xlab, y = ylab, tag = tag, colour = group) +
    ggplot2::guides(colour = md_guide_flex(length(pal), "bottom", length(pal))) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  if (use_shape) {
    g <- g + ggplot2::scale_shape_manual(values = shape_vals, name = shape_var) +
      ggplot2::guides(shape = md_guide_flex(length(shape_vals), "bottom", length(shape_vals)))
  }
  if (!is.null(subtitle)) g <- g + ggplot2::labs(subtitle = subtitle)
  g
}

# ---------------------------------------------------------------------------
# 模块主体
# ---------------------------------------------------------------------------
md_ordination <- function(CFG) {
  module <- "14"
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少 ecology 包组"))

  dc <- tryCatch(md_dist_cache(CFG), error = function(e) NULL)
  if (is.null(dc) || is.null(dc$dist) || length(dc$dist) == 0) {
    msg <- "cache/13_dist.rds 不存在或为空 —— 请先运行模块 13（--only=13,14）"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }
  dist_list <- dc$dist
  meta <- dc$meta

  method <- tolower(CFG$ordination$method %||% "pcoa")
  do_pcoa <- method %in% c("pcoa", "both", "pcoa+nmds")
  do_nmds <- method %in% c("nmds", "both", "pcoa+nmds")
  if (!do_pcoa && !do_nmds) {
    msg <- paste0("ordination$method = '", method, "' 不可识别（pcoa|nmds|both）")
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }
  ell_on <- isTRUE(CFG$ordination$ellipse %||% TRUE)
  ell_level <- CFG$ordination$ellipse_level %||% 0.95
  ell_min_n <- as.integer(CFG$ordination$ellipse_min_n %||% 4)
  trymax <- as.integer(CFG$ordination$nmds_trymax %||% 200)
  k_axes <- as.integer(CFG$ordination$nmds_k %||% 2)

  factors <- intersect(CFG$design$factors, names(meta))
  factors <- factors[vapply(factors, function(f) length(unique(stats::na.omit(meta[[f]]))) >= 2, logical(1))]
  if (length(factors) == 0) {
    msg <- "元数据中没有可用的分组因子"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }

  axis_rows <- list(); ell_rows <- list(); score_rows <- list(); stress_rows <- list()
  note_extra <- character(0)

  md_seed(CFG, 140)
  for (mname in names(dist_list)) {
    d <- dist_list[[mname]]
    d <- stats::as.dist(as.matrix(d))
    if (length(d) < 2) { note_extra <- c(note_extra, paste0(mname, ": 样本不足")); next }

    # ---- PCoA ----
    if (do_pcoa) {
      pc <- tryCatch(ape::pcoa(d, correction = "cailliez"), error = function(e) NULL)
      if (is.null(pc)) {
        md_warn(CFG, "ape::pcoa 失败: ", mname)
        note_extra <- c(note_extra, paste0(mname, ": PCoA 失败"))
      } else {
        vec <- pc$vectors
        if (is.null(vec) || ncol(vec) < 1) {
          note_extra <- c(note_extra, paste0(mname, ": PCoA 特征向量为空"))
        } else {
          n_ax <- min(2L, ncol(vec))
          pct <- local_eig_pct(pc, n_axis = 2L)
          df <- data.frame(sample_id = rownames(vec), Axis1 = as.numeric(vec[, 1]),
                           Axis2 = if (n_ax >= 2) as.numeric(vec[, 2]) else 0,
                           stringsAsFactors = FALSE)
          df <- merge(df, meta, by.x = "sample_id", by.y = CFG$data$sample_col,
                      all.x = TRUE, sort = FALSE)
          df <- df[match(rownames(vec), df$sample_id), , drop = FALSE]
          df$metric <- mname; df$method <- "PCoA"
          axis_rows[[length(axis_rows) + 1L]] <- data.frame(
            metric = mname, method = "PCoA",
            axis = c("Axis1", "Axis2"),
            eigenvalue = c(pc$values$Corr_eig[1] %||% NA_real_,
                           if (!is.null(pc$values$Corr_eig) && length(pc$values$Corr_eig) >= 2)
                             pc$values$Corr_eig[2] else NA_real_),
            explained_pct = pct,
            stringsAsFactors = FALSE)
          score_rows[[length(score_rows) + 1L]] <- df
          for (f in factors) {
            eg <- local_ellipse_groups(df, f, ell_min_n)
            if (!is.null(eg$dropped)) {
              e <- eg$dropped; e$metric <- mname; e$method <- "PCoA"; e$factor <- f
              if (!ell_on) { e$drawn <- FALSE; e$reason <- "ellipse = FALSE" }
              ell_rows[[length(ell_rows) + 1L]] <- e
            }
          }
        }
      }
    }

    # ---- NMDS ----
    if (do_nmds) {
      md_seed(CFG, 141)
      nm <- tryCatch(vegan::metaMDS(d, k = k_axes, trymax = trymax, trace = FALSE),
                     error = function(e) { md_warn(CFG, "metaMDS 失败: ", mname, " ", conditionMessage(e)); NULL })
      if (is.null(nm)) {
        note_extra <- c(note_extra, paste0(mname, ": NMDS 失败"))
      } else {
        pts <- as.data.frame(nm$points, stringsAsFactors = FALSE)
        names(pts)[1:2] <- c("Axis1", "Axis2")
        df <- data.frame(sample_id = rownames(pts), Axis1 = pts$Axis1, Axis2 = pts$Axis2,
                         stringsAsFactors = FALSE)
        df <- merge(df, meta, by.x = "sample_id", by.y = CFG$data$sample_col,
                    all.x = TRUE, sort = FALSE)
        df <- df[match(rownames(pts), df$sample_id), , drop = FALSE]
        df$metric <- mname; df$method <- "NMDS"
        st <- nm$stress %||% NA_real_
        stress_rows[[length(stress_rows) + 1L]] <- data.frame(
          metric = mname, method = "NMDS", stress = as.numeric(st),
          stress_type = if (!is.null(st)) "vegan metaMDS stress" else NA_character_,
          stringsAsFactors = FALSE)
        axis_rows[[length(axis_rows) + 1L]] <- data.frame(
          metric = mname, method = "NMDS", axis = c("NMDS1", "NMDS2"),
          eigenvalue = c(NA_real_, NA_real_),
          explained_pct = c(NA_real_, NA_real_), stringsAsFactors = FALSE)
        score_rows[[length(score_rows) + 1L]] <- df
        for (f in factors) {
          eg <- local_ellipse_groups(df, f, ell_min_n)
          if (!is.null(eg$dropped)) {
            e <- eg$dropped; e$metric <- mname; e$method <- "NMDS"; e$factor <- f
            e$stress <- st
            if (!ell_on) { e$drawn <- FALSE; e$reason <- "ellipse = FALSE" }
            ell_rows[[length(ell_rows) + 1L]] <- e
          }
        }
        note_extra <- c(note_extra, paste0(mname, ": NMDS stress = ", round(st, 4)))
      }
    }
  }

  tab_axis <- if (length(axis_rows)) do.call(rbind, axis_rows) else NULL
  tab_ell <- if (length(ell_rows)) do.call(rbind, ell_rows) else NULL
  tab_stress <- if (length(stress_rows)) do.call(rbind, stress_rows) else NULL
  scores <- if (length(score_rows)) do.call(rbind, score_rows) else NULL
  if (!is.null(tab_axis)) tab_axis$explained_pct <- round(tab_axis$explained_pct, 2)
  if (!is.null(tab_ell)) {
    tab_ell$sd_axis1 <- round(tab_ell$sd_axis1, 6)
    tab_ell$sd_axis2 <- round(tab_ell$sd_axis2, 6)
    tab_ell <- tab_ell[, c("metric", "method", "factor", "group", "n",
                           "sd_axis1", "sd_axis2", "drawn", "reason")]
  }

  # ---- 出图：每个距离度量一个面板（照 02_Fig02.R §B/C 的构图）----
  #   原图是「一个距离度量 = 一个 PCoA 面板」，面板内按 Zone 着色、按 Vegetation 给形状；
  #   故这里不再按因子拆面板，颜色固定用第一个设计因子、形状用第二个（若存在）。
  panels <- list()
  colour_var <- factors[1]
  shape_var  <- if (length(factors) >= 2) factors[2] else NULL
  if (!is.null(scores)) {
    for (mname in unique(scores$metric)) {
      for (meth in unique(scores$method[scores$metric == mname])) {
        sub0 <- scores[scores$metric == mname & scores$method == meth, , drop = FALSE]
        pc1 <- if (meth == "PCoA" && !is.null(tab_axis)) {
          z <- tab_axis[tab_axis$metric == mname & tab_axis$method == "PCoA" & tab_axis$axis == "Axis1", ]
          if (nrow(z)) z$explained_pct[1] else NA_real_
        } else NA_real_
        pc2 <- if (meth == "PCoA" && !is.null(tab_axis)) {
          z <- tab_axis[tab_axis$metric == mname & tab_axis$method == "PCoA" & tab_axis$axis == "Axis2", ]
          if (nrow(z)) z$explained_pct[1] else NA_real_
        } else NA_real_
        # 轴标题文字照抄原脚本：round(特征值占比, 2) 直接进括号
        xlab <- if (is.na(pc1)) paste0(meth, "1") else
          paste0("PCoA1 (", as.character(round(pc1, 2)), "%)")
        ylab <- if (is.na(pc2)) paste0(meth, "2") else
          paste0("PCoA2 (", as.character(round(pc2, 2)), "%)")
        pal <- md_pal(CFG, local_factor_levels(CFG, meta, colour_var), colour_var)
        subt <- NULL
        if (meth == "NMDS" && !is.null(tab_stress)) {
          z <- tab_stress[tab_stress$metric == mname & tab_stress$method == "NMDS", , drop = FALSE]
          if (nrow(z) && is.finite(z$stress[1]))
            subt <- paste0("stress = ", formatC(z$stress[1], format = "f", digits = 3))
        }
        pl <- tryCatch(
          local_ord_panel(sub0, group = colour_var, pal = pal, xlab = xlab, ylab = ylab,
                          tag = NULL, ellipse_level = ell_level, min_n = ell_min_n,
                          method = meth, shape_var = shape_var, subtitle = subt),
          error = function(e) { md_warn(CFG, "排序面板失败 ", mname, "/", meth, ": ",
                                        conditionMessage(e)); md_empty_plot(CFG, mname) })
        panels[[length(panels) + 1L]] <- pl
      }
    }
  }
  if (length(panels) == 0) panels <- list(md_empty_plot(CFG, "Ordination not available"))

  fig <- md_compose(CFG, panels, ncol = min(2L, length(panels)))
  h <- if (length(panels) <= 2) 105 else if (length(panels) <= 4) 165 else 225
  md_save_figure(CFG, fig, "14_ordination", height = h)
  for (i in seq_along(panels)) md_maybe_single(CFG, panels[[i]], paste0("14_ordination_P", i))

  # ---- 表 ----
  note <- data.frame(
    item = c("method", "distances used", "ellipse", "ellipse level", "ellipse_min_n",
             "ellipse type", "PCoA correction", "axis % definition",
             "NMDS trymax", "PERMDISP", "Notes"),
    value = c(method, paste(names(dist_list), collapse = ", "), as.character(ell_on),
              as.character(ell_level), as.character(ell_min_n), "stat_ellipse(type = 'norm')",
              "cailliez", "eigenvalue / sum(positive eigenvalues)",
              as.character(trymax), "not computed here (module 15)",
              paste(note_extra, collapse = " | ")),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "14_ordination", list(
    Axis_variance = tab_axis,
    Ellipse_diagnostics = tab_ell,
    NMDS_stress = tab_stress,
    Site_scores = if (isTRUE(CFG$table$detail)) scores else
      if (!is.null(scores)) scores[, intersect(c("sample_id", "metric", "method", "Axis1", "Axis2",
                                                 factors), names(scores)), drop = FALSE] else NULL,
    Notes = note))

  md_result(module, "OK", tables = basename(tp), figures = "figures/14_ordination.pdf",
            notes = paste0("method=", method, "; metrics=", paste(names(dist_list), collapse = "+"),
                           if (length(note_extra)) paste0("; ", paste(note_extra, collapse = " | ")) else ""))
}
