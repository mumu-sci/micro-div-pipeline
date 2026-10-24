# =============================================================================
# R/70_response_pattern.R ---- 响应变量的格局（因子格局 + 连续环境梯度）
# -----------------------------------------------------------------------------
# 契约见 _build-notes/micro-div-pipeline/MODULE_SPEC.md
#
# 做什么（只做"格局"，不做条件效应、不做关联分析 ---- 那两件事属于 71）
#   A  双因子交互箱线图（CFG$response_models$pattern_test，KW + Dunn + CLD 字母）
#   B  单因子各自的箱线图（每个 design$factors 一个）
#   C  响应变量 vs 连续环境变量的散点 + 回归拟合
#   图 = 多面板拼一张；表 = 描述统计 / 检验 / 事后 / 交互模型 / 关联汇总 / 注释
#
# 参数来源（全部 CFG，未给的走 %||% 默认）
#   CFG$response$enabled / var / label / unit / trans
#   CFG$response_models$pattern_test
#   CFG$design$factors / ordered_levels
#   CFG$env$vars, CFG$model$covariates / focus / derived, CFG$model$standardize
#   CFG$alpha$p_adjust
#
# 依赖：05（cache/05_core.rds）。
# =============================================================================

md_response_pattern <- function(CFG) {
  module <- "70"

  if (!md_precheck(CFG, module, need_response = TRUE))
    return(md_result(module, "SKIP", "response$enabled = FALSE 或前置条件不满足"))
  if (!md_need(CFG, "ecology", module))
    return(md_result(module, "SKIP", "缺少统计依赖包"))

  d <- tryCatch(md_analysis_with_derived(CFG), error = function(e) NULL)
  if (is.null(d)) {
    md_skip(CFG, module, "取不到分析总表（cache/05_core.rds 缺失？）")
    return(md_result(module, "SKIP", "取不到分析总表"))
  }

  yv    <- CFG$response$var
  ylab  <- CFG$response$label %||% yv
  yunit <- CFG$response$unit %||% ""
  if (!yv %in% names(d)) {
    md_skip(CFG, module, paste0("响应变量列 ", yv, " 不在数据中"))
    return(md_result(module, "SKIP", paste0("缺响应变量列 ", yv)))
  }

  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  # ---- 响应变量变换 ----
  y <- suppressWarnings(as.numeric(d[[yv]]))
  tr <- CFG$response$trans %||% "none"
  y_t <- switch(tr,
                log  = if (all(y[is.finite(y)] > 0)) log(y) else { add_note(
                  "response$trans = log 但存在 <= 0 的值，本模块未做变换。"); y },
                sqrt = if (all(y[is.finite(y)] >= 0)) sqrt(y) else { add_note(
                  "response$trans = sqrt 但存在 < 0 的值，本模块未做变换。"); y },
                z    = as.numeric(scale(y)),
                y)
  d$.y <- y_t
  n_ok <- sum(is.finite(d$.y))
  if (n_ok < 4) {
    md_skip(CFG, module, paste0("响应变量有效观测只有 ", n_ok, " 个，无法做格局分析"))
    return(md_result(module, "SKIP", "响应变量有效观测不足"))
  }
  if (stats::sd(d$.y, na.rm = TRUE) <= 0) {
    md_skip(CFG, module, "响应变量零方差，无法做格局分析")
    return(md_result(module, "SKIP", "响应变量零方差"))
  }
  y_axis <- if (identical(tr, "none")) {
    if (nzchar(yunit)) paste0(ylab, " (", yunit, ")") else ylab
  } else paste0(ylab, " (", tr, "-transformed)")

  fac <- CFG$design$factors %||% character(0)
  fac <- fac[fac %in% names(d)]
  ord_lv <- CFG$design$ordered_levels %||% list()
  test <- CFG$response_models$pattern_test %||% "kruskal_dunn"
  padj <- CFG$model$p_adjust %||% CFG$alpha$p_adjust %||% "BH"

  # ---------------------------------------------------------------------------
  # A0. 备选口径：Scheirer-Ray-Hare 双因子非参数检验
  #   参考实现（rcompanion::scheirerRayHare；
  #       缺包退回 ART + 双因素 ANOVA，见 shared/11_stats.R 的 md_srh()）
  #   ⚠ 只在 pattern_test = "scheirer_ray_hare" 时执行；默认口径下图件与表**完全不变**。
  #     SRH 不做事后两两比较，所以单因子检验/CLD 字母仍走 KW + Dunn。
  # ---------------------------------------------------------------------------
  srh_tab <- NULL
  if (identical(as.character(test), "scheirer_ray_hare")) {
    if (length(fac) >= 2) {
      srh_tab <- tryCatch(md_srh(CFG, d, ".y", fac[1], fac[2]), error = function(e) NULL)
      if (is.null(srh_tab)) {
        add_note("Scheirer-Ray-Hare 未算出（rcompanion 缺失或水平/样本量不足），" ,
                 "已退回 KW + Dunn 口径。")
      } else {
        add_note("pattern_test = scheirer_ray_hare：双因子非参数检验见 SRH_two_way sheet；",
                 "面板 A 的右上角注释改用 SRH 的 P。SRH 不做事后两两比较，",
                 "字母仍来自 KW + Dunn。")
      }
    } else {
      add_note("scheirer_ray_hare 需要 2 个设计因子，本数据只有 ", length(fac),
               " 个，已退回 KW + Dunn 口径。")
    }
    test <- "kruskal_dunn"
  }

  # ---------------------------------------------------------------------------
  # A. 双因子交互格局
  # ---------------------------------------------------------------------------
  desc <- NULL; test_tab <- NULL; posthoc <- NULL; letters_df <- NULL
  letters_inter <- NULL; n_dropped <- character(0)
  if (length(fac) >= 2) {
    f1 <- fac[1]; f2 <- fac[2]
    d$.inter <- factor(paste(as.character(d[[f1]]), as.character(d[[f2]]), sep = "_"))
    gt <- md_group_test(CFG, d, ".y", ".inter", test = test, p_adjust = padj,
                        min_n = 3, order_levels = NULL)
    desc <- gt$desc; test_tab <- gt$test; posthoc <- gt$posthoc
    letters_df <- gt$letters; n_dropped <- gt$dropped
    if (length(n_dropped))
      add_note("交互组样本量 < 3 的组已从检验中剔除：",
               paste(n_dropped, collapse = ", "), "（仍出现在箱线图上）。")
    if (is.null(posthoc) || nrow(posthoc) == 0)
      add_note("交互组事后两两比较不可用（组数 < 3 或 rstatix 缺失），CLD 字母退化为单一字母。")
    if (!is.null(letters_df)) letters_inter <- letters_df
  } else {
    add_note("design$factors 少于 2 个，跳过双因子交互面板。")
  }

  # 单因子检验
  single <- list()
  for (f in fac) {
    g <- md_group_test(CFG, d, ".y", f, test = test, p_adjust = padj,
                       min_n = 3, order_levels = ord_lv[[f]])
    single[[f]] <- g
  }

  # ---------------------------------------------------------------------------
  # B. 响应变量 vs 连续环境变量
  # ---------------------------------------------------------------------------
  cont_vars <- unique(c(CFG$env$vars %||% character(0),
                        CFG$model$covariates %||% character(0),
                        CFG$model$focus %||% character(0)))
  cont_vars <- cont_vars[cont_vars %in% names(d)]
  cont_vars <- cont_vars[vapply(d[cont_vars], function(x) is.numeric(x) &&
                                  sum(is.finite(x)) >= 4 &&
                                  stats::sd(x, na.rm = TRUE) > 0, logical(1))]
  cont_vars <- setdiff(cont_vars, yv)
  if (length(cont_vars) == 0)
    add_note("没有可用的连续环境变量，散点面板记为空白。")

  cont_tab <- NULL
  if (length(cont_vars) > 0) {
    cont_tab <- do.call(rbind, lapply(cont_vars, function(v) {
      dd <- data.frame(x = as.numeric(d[[v]]), y = d$.y)
      dd <- dd[stats::complete.cases(dd), , drop = FALSE]
      if (nrow(dd) < 4) return(NULL)
      lmf <- tryCatch(stats::lm(y ~ x, data = dd), error = function(e) NULL)
      sp  <- tryCatch(suppressWarnings(stats::cor.test(dd$x, dd$y, method = "spearman")),
                      error = function(e) NULL)
      data.frame(variable = v, n = nrow(dd),
                 slope = if (is.null(lmf)) NA_real_ else unname(stats::coef(lmf)[2]),
                 r2 = if (is.null(lmf)) NA_real_ else summary(lmf)$r.squared,
                 p_lm = if (is.null(lmf)) NA_real_ else
                   unname(summary(lmf)$coefficients[2, "Pr(>|t|)"]),
                 rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
                 p_spearman = if (is.null(sp)) NA_real_ else sp$p.value,
                 stringsAsFactors = FALSE)
    }))
    if (!is.null(cont_tab) && nrow(cont_tab)) {
      cont_tab$p_lm_BH <- md_p_adjust(cont_tab$p_lm, padj)
      cont_tab$p_spearman_BH <- md_p_adjust(cont_tab$p_spearman, padj)
    }
  }

  # ---------------------------------------------------------------------------
  # C. 交互模型（含交互项的 lm；响应与预测变量标准化的双份口径）
  # ---------------------------------------------------------------------------
  inter_mod <- NULL
  if (length(fac) >= 2) {
    f1 <- fac[1]; f2 <- fac[2]
    dm <- d[stats::complete.cases(d[, c(".y", f1, f2)]), , drop = FALSE]
    dm[[f1]] <- factor(as.character(dm[[f1]]), levels = ord_lv[[f1]] %||% unique(as.character(dm[[f1]])))
    dm[[f2]] <- factor(as.character(dm[[f2]]), levels = ord_lv[[f2]] %||% unique(as.character(dm[[f2]])))
    if (nlevels(droplevels(dm[[f1]])) >= 2 && nlevels(droplevels(dm[[f2]])) >= 2 &&
        nrow(dm) >= 6) {
      fit <- tryCatch(stats::lm(stats::as.formula(paste0(".y ~ `", f1, "` * `", f2, "`")),
                                data = dm), error = function(e) NULL)
      if (!is.null(fit)) {
        sm <- summary(fit)
        co <- as.data.frame(sm$coefficients, stringsAsFactors = FALSE)
        names(co) <- c("Estimate", "Std_Error", "t_value", "P_value")
        co$Term <- rownames(co); rownames(co) <- NULL
        co$R2 <- sm$r.squared
        co$Adj_R2 <- sm$adj.r.squared
        co$F_interaction <- tryCatch({
          an <- stats::anova(fit)
          r <- grep(":", rownames(an))
          if (length(r)) an$F[r[1]] else NA_real_
        }, error = function(e) NA_real_)
        co$P_interaction <- tryCatch({
          an <- stats::anova(fit)
          r <- grep(":", rownames(an))
          if (length(r)) an$`Pr(>F)`[r[1]] else NA_real_
        }, error = function(e) NA_real_)
        co$Term_type <- ifelse(co$Term == "(Intercept)", "Intercept",
                               ifelse(grepl(":", co$Term), "Interaction", "Main effect"))
        inter_mod <- co[, c("Term", "Term_type", "Estimate", "Std_Error",
                            "t_value", "P_value", "R2", "Adj_R2",
                            "F_interaction", "P_interaction")]
      } else {
        add_note("交互项 lm 拟合失败，交互模型表记 NA。")
      }
    } else {
      add_note("某个因子不足 2 个水平或样本量 < 6，交互项 lm 未拟合。")
    }
  }

  # ---------------------------------------------------------------------------
  # D. 出图
  # ---------------------------------------------------------------------------
  pt <- md_pt(CFG)
  f1 <- if (length(fac) >= 1) fac[1] else NULL
  f2 <- if (length(fac) >= 2) fac[2] else NULL

  if (!is.null(f2)) {
    # inter_sep 必须与上面生成 .inter 用的 "_" 一致，否则字母与 x 水平对不上被丢弃
    d$.inter <- factor(paste(as.character(d[[f1]]), as.character(d[[f2]]), sep = "_"))
    kw <- md_kw_text(stats::setNames(
      lapply(c(f1, f2), function(ff)
        if (!is.null(single[[ff]]$test) && nrow(single[[ff]]$test) >= 1)
          single[[ff]]$test$p[1] else NA_real_),
      c(f1, f2)))
    # SRH 口径：右上角改用双因子非参数检验的 P（把 2 个主效应压成一行，纯 ASCII）
    if (!is.null(srh_tab) && nrow(srh_tab) >= 2) {
      fx <- c(f1, f2); fx <- gsub("[^ -~]", "", fx)
      kw <- paste(sprintf("%s %s", fx, md_fmt_p_fig(srh_tab$p.value[1:2])), collapse = "; ")
    }
    pA <- md_box_two_factor(d, y = ".y", f1 = f1, f2 = f2,
                            letters = letters_inter,
                            colors = md_pal(CFG, ord_lv[[f2]] %||% unique(as.character(d[[f2]])), f2),
                            f1_levels = ord_lv[[f1]], f2_levels = ord_lv[[f2]],
                            ylab = y_axis, tag = "A", x_angle = 45,
                            legend_pos = "none", inter_sep = "_",
                            x_labels = "two_line",
                            kw_text = kw, pt = pt)
  } else if (!is.null(f1)) {
    pA <- md_box_group(d, y = ".y", x = f1, letters = single[[f1]]$letters,
                       colors = md_pal(CFG, ord_lv[[f1]] %||% unique(as.character(d[[f1]])), f1),
                       ylab = y_axis, xlab = f1, tag = "A", pt = pt,
                       kw_text = md_kw_text(stats::setNames(
                         list(if (!is.null(single[[f1]]$test) && nrow(single[[f1]]$test) >= 1)
                           single[[f1]]$test$p[1] else NA_real_), f1)))
  } else {
    pA <- md_empty_plot(CFG, "no design factor") + ggplot2::labs(tag = "A")
  }
  pA <- pA + theme_md(CFG = CFG)

  # 单因子箱线（B/C）**不再出图** —— 用户口径：A 已经是两因子交互格局，
  # 单因子只是它的**子集**，重复呈现浪费版面。
  # （单因子检验的 p 值仍照常计算并写在 Tests / Posthoc 表里，只是不画。）
  single_plots <- list()

  # 连续环境变量散点（最多 2 个面板，按 |rho| 取最强）
  scat_plots <- list()
  if (length(cont_vars) > 0) {
    pick <- cont_vars
    if (!is.null(cont_tab) && nrow(cont_tab)) {
      ordc <- cont_tab$variable[order(-abs(cont_tab$rho))]
      pick <- head(unique(c(ordc, cont_vars)), 2)
    } else pick <- head(cont_vars, 2)
    # 着色因子沿用同一套样式 §B/D/E：按第二因子（Vegetation）上色 + cols_veg
    sc_col <- if (!is.null(f2)) f2 else f1
    sc_lv <- if (!is.null(sc_col))
      (ord_lv[[sc_col]] %||% unique(as.character(d[[sc_col]]))) else NULL
    for (i in seq_along(pick)) {
      v <- pick[i]
      g <- md_scatter_fit(d, x = v, y = ".y", color = sc_col,
                          colors = if (!is.null(sc_col)) md_pal(CFG, sc_lv, sc_col),
                          line = "lm", ribbon = TRUE,
                          xlab = v, ylab = y_axis,
                          tag = LETTERS[1 + i], pt = pt)
      if (!is.null(cont_tab)) {
        r <- cont_tab[cont_tab$variable == v, , drop = FALSE]
        if (nrow(r))
          # 标注锚在面板右上角（同 原脚本 §B/D/E 的 ann：
          # x=Inf / hjust=1 / vjust=1.3 / color="grey30"）
          g <- md_add_p_annot(g,
                              labels = sprintf("R2 = %.3f, P = %s",
                                               r$r2[1], md_fmt_p(r$p_lm[1])),
                              x = Inf, hjust = 1, vjust = 1.3,
                              colour = "grey30", pt = pt)
      }
      scat_plots[[v]] <- g + theme_md(CFG = CFG)
    }
  }

  all_panels <- c(list(main = pA), single_plots, scat_plots)
  n_bottom <- length(all_panels) - 1
  # 图例：原 Fig.5 全图**没有图例**（A 面板 legend.position="none"，其余面板同样 none，
  # 组合层只加 theme_comp()）—— 故组合层也不放图例，否则会把各面板自己的 "none"
  # 覆盖成底部共用图例。
  if (n_bottom == 0) {
    fig <- pA
  } else if (n_bottom <= 3) {
    # 顶行 A 整宽；底行按 **连续等分** 给 B..。
    # patchwork 要求 design 每行等宽，且每个面板占的格子必须构成**矩形**，
    # 所以不能用 "BCBCBC" 这种交错写法（会报 "Patch areas must be rectangular"）。
    len <- 12L
    per <- max(1L, len %/% n_bottom)
    # ⚠ 底行的字母必须从 **B** 开始（A 已被主面板占用），否则同一个字母会占两块
    #   不连续的区域，patchwork 报 "Patch areas must be rectangular"。
    bot_letters <- LETTERS[1L + seq_len(n_bottom)]
    design <- paste0(paste(rep("A", len), collapse = ""), "\n",
                     paste(rep(bot_letters, each = per), collapse = ""))
    fig <- md_compose(CFG, all_panels, design = design,
                      guides = "collect", legend = "none")
  } else {
    fig <- md_compose(CFG, all_panels, ncol = 3, guides = "collect", legend = "none")
  }
  md_save_figure(CFG, fig, "70_response_pattern",
                 height = if (n_bottom > 3) 150 else 110)

  # ---------------------------------------------------------------------------
  # E. 出表
  # ---------------------------------------------------------------------------
  desc_all <- NULL
  for (f in fac) {
    dd0 <- d[stats::complete.cases(d[, c(".y", f)]), , drop = FALSE]
    if (nrow(dd0) == 0) next
    tb <- do.call(rbind, lapply(split(dd0$.y, droplevels(factor(dd0[[f]]))), function(x) {
      data.frame(n = length(x), mean = mean(x), sd = stats::sd(x),
                 median = stats::median(x),
                 q1 = unname(stats::quantile(x, .25)),
                 q3 = unname(stats::quantile(x, .75)),
                 min = min(x), max = max(x), stringsAsFactors = FALSE)
    }))
    desc_all <- rbind(desc_all, data.frame(Factor = rep(f, nrow(tb)),
                                           Level = rownames(tb), tb,
                                           row.names = NULL, stringsAsFactors = FALSE))
  }
  if (!is.null(desc_all)) rownames(desc_all) <- NULL

  test_all <- NULL
  for (f in fac) {
    t0 <- single[[f]]$test
    if (!is.null(t0)) test_all <- rbind(test_all, data.frame(Factor = f, t0,
                                                             stringsAsFactors = FALSE))
  }
  if (!is.null(test_tab)) {
    test_all <- rbind(test_all, data.frame(Factor = paste(fac[1:2], collapse = " x "),
                                           test_tab, stringsAsFactors = FALSE))
  }
  if (is.null(test_all)) test_all <- data.frame(Note = "no test available")

  post_all <- NULL
  for (f in fac) {
    p0 <- single[[f]]$posthoc
    if (!is.null(p0)) post_all <- rbind(post_all, data.frame(Factor = f, p0,
                                                             stringsAsFactors = FALSE))
  }
  if (!is.null(posthoc)) {
    post_all <- rbind(post_all, data.frame(Factor = paste(fac[1:2], collapse = " x "),
                                           posthoc, stringsAsFactors = FALSE))
  }
  if (is.null(post_all)) post_all <- data.frame(Note = "no post-hoc comparison available")

  note_tab <- data.frame(
    Item = c("Response variable", "Response label", "Response unit",
             "Response transform", "Pattern test", "p adjust",
             "Design factors", "n usable", "CLD letters",
             "Scope of this module", "Notes"),
    Value = c(yv, ylab, yunit, tr, test, padj,
              paste(fac, collapse = " x "), n_ok,
              if (is.null(letters_inter)) "NA" else
                paste(sprintf("%s=%s", letters_inter$.inter, letters_inter$letter), collapse = "; "),
              "descriptive pattern only; conditional effects and community linkage belong to module 71",
              paste(notes, collapse = " ")),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "70_response_pattern", list(
    Descriptive = desc_all,
    Tests = test_all,
    Posthoc = post_all,
    Interaction_model = inter_mod,
    Continuity = cont_tab,
    SRH_two_way = srh_tab,
    Notes = note_tab
  ))

  md_result(module, "OK",
            tables = basename(tp), figures = "figures/70_response_pattern.pdf",
            notes = notes)
}
