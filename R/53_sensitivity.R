# =============================================================================
# R/53_sensitivity.R —— 敏感性分析：剔除特定组/站点/样本后重跑关键模型并并列对比
# -----------------------------------------------------------------------------
# 把 原始研究里写死的「剔除人工岛（南沙）敏感性分析」做成**通用机制**：
#   配置里声明"剔掉哪些因子水平 / 哪些站点 / 哪些样本"，
#   模块对每个子集重跑同一个模型，把关键系数与 P 值并列成表 + 森林图。
#
# 为什么需要：审稿人一定会问「结论是不是被某几个样本/某个特殊组驱动的」。
#
# 配置键（CFG$sensitivity，缺省值见下）：
#   enabled       TRUE/FALSE
#   drop_levels   list(因子名 = c("要剔除的水平"))   例如 list(Zone = "Xisha")
#   drop_sites    c("站点名")                        按 design$site_var 匹配
#   drop_samples  c("sample_id")
#   target        c("lmm") 或 c("lmm","sem")         要重跑哪类模型
#   max_subsets   最多展示几个子集（防止爆炸）
#   min_n         子集样本量低于此值就不跑该子集
# =============================================================================

md_sensitivity <- function(CFG) {
  module <- "53"

  sens <- CFG$sensitivity %||% list()
  enabled <- sens$enabled %||% TRUE
  drop_levels <- sens$drop_levels %||% list()
  drop_sites  <- sens$drop_sites %||% character(0)
  drop_samples <- sens$drop_samples %||% character(0)
  target      <- sens$target %||% c("lmm")
  max_subsets <- sens$max_subsets %||% 6L
  min_n       <- sens$min_n %||% 12L

  if (!isTRUE(enabled)) {
    md_skip(CFG, module, "sensitivity$enabled = FALSE")
    return(md_result(module, "SKIP", "sensitivity$enabled = FALSE"))
  }
  if (length(drop_levels) == 0 && length(drop_sites) == 0 && length(drop_samples) == 0) {
    md_skip(CFG, module, "未声明任何剔除项（sensitivity$drop_levels / drop_sites / drop_samples 都为空）")
    return(md_result(module, "SKIP", "未声明剔除项"))
  }
  if (!md_need(CFG, "model", module)) {
    return(md_result(module, "SKIP", "缺少建模依赖包"))
  }

  d_all <- tryCatch(md_analysis_with_derived(CFG), error = function(e) NULL)
  if (is.null(d_all)) {
    md_skip(CFG, module, "取不到分析总表（需要 05 与 40 的产物）")
    return(md_result(module, "SKIP", "上游产物缺失"))
  }

  smp     <- CFG$data$sample_col
  y       <- CFG$response$var
  focus   <- CFG$model$focus %||% character(0)
  covs    <- CFG$model$covariates %||% character(0)
  random  <- CFG$model$random
  sitevar <- CFG$design$site_var

  if (is.null(d_all[[smp]])) d_all[[smp]] <- rownames(d_all)
  if (isTRUE(CFG$response$enabled) && !y %in% names(d_all)) {
    md_skip(CFG, module, paste0("响应变量 ", y, " 不在分析总表里"))
    return(md_result(module, "SKIP", "响应变量缺失"))
  }
  keep_vars <- intersect(c(smp, y, focus, covs, sitevar, CFG$design$factors),
                         names(d_all))
  d_all <- d_all[, keep_vars, drop = FALSE]

  # ---- 1. 构造子集 ----
  subsets <- list(Full = d_all)
  for (f in names(drop_levels)) {
    if (!f %in% names(d_all)) {
      md_warn(CFG, "sensitivity$drop_levels 里的因子 ", f, " 不在数据中，跳过")
      next
    }
    lv <- as.character(drop_levels[[f]])
    d <- d_all[!as.character(d_all[[f]]) %in% lv, , drop = FALSE]
    nm <- paste0("drop_", f, "_", paste(lv, collapse = "+"))
    subsets[[nm]] <- d
  }
  if (length(drop_sites) && !is.null(sitevar) && sitevar %in% names(d_all)) {
    d <- d_all[!as.character(d_all[[sitevar]]) %in% as.character(drop_sites), , drop = FALSE]
    subsets[[paste0("drop_sites_", length(drop_sites))]] <- d
  }
  if (length(drop_samples)) {
    d <- d_all[!as.character(d_all[[smp]]) %in% as.character(drop_samples), , drop = FALSE]
    subsets[[paste0("drop_samples_", length(drop_samples))]] <- d
  }

  # 超上限时只保留 Full + 前几个
  if (length(subsets) > max_subsets + 1L) {
    md_warn(CFG, "子集数 ", length(subsets) - 1L, " 超过 max_subsets = ", max_subsets,
            "，只展示前 ", max_subsets, " 个")
    subsets <- subsets[seq_len(max_subsets + 1L)]
  }

  # ---- 2. 逐个子集重跑模型 ----
  coef_rows <- list(); fit_rows <- list(); skip_rows <- list()

  for (nm in names(subsets)) {
    d <- subsets[[nm]]
    d <- d[stats::complete.cases(d[, intersect(c(y, focus, covs), names(d)), drop = FALSE]), ,
           drop = FALSE]
    if (nrow(d) < min_n) {
      skip_rows[[length(skip_rows) + 1L]] <- data.frame(
        subset = nm, n = nrow(d), reason = paste0("样本量 < min_n(", min_n, ")"),
        stringsAsFactors = FALSE)
      next
    }
    if (length(unique(na.omit(d[[y]]))) < 3) {
      skip_rows[[length(skip_rows) + 1L]] <- data.frame(
        subset = nm, n = nrow(d), reason = "响应变量取值过少", stringsAsFactors = FALSE)
      next
    }

    # --- LMM / 普通回归 ---
    if ("lmm" %in% target) {
      rhs_terms <- unique(c(focus, covs))
      rhs_terms <- rhs_terms[rhs_terms %in% names(d)]
      if (length(rhs_terms) == 0) {
        skip_rows[[length(skip_rows) + 1L]] <- data.frame(
          subset = nm, n = nrow(d), reason = "焦点变量与协变量都不在数据中",
          stringsAsFactors = FALSE)
      } else {
        # 站点级共线守卫：站内唯一值 = 1 的变量不能与站点随机效应同时进模型
        use_re <- !is.null(random) && nzchar(random)
        if (use_re && !is.null(sitevar) && sitevar %in% names(d)) {
          bad <- character(0)
          for (v in rhs_terms) {
            u <- tryCatch(md_within_unique(d, v, sitevar), error = function(e) NA_integer_)
            if (!is.na(u) && u <= 1) bad <- c(bad, v)
          }
          if (length(bad)) {
            md_log(CFG, "INFO", "  [", nm, "] 站点级共线守卫剔除: ", paste(bad, collapse = ", "))
            rhs_terms <- setdiff(rhs_terms, bad)
            if (length(rhs_terms) == 0) {
              skip_rows[[length(skip_rows) + 1L]] <- data.frame(
                subset = nm, n = nrow(d), reason = "共线守卫后没有剩余预测变量",
                stringsAsFactors = FALSE)
            }
          }
        }
        if (length(rhs_terms) > 0) {
          fml <- stats::as.formula(paste(y, "~", paste(rhs_terms, collapse = " + "),
                                         if (use_re) paste("+", random) else ""))
          md_seed(CFG, 2001)
          fit <- tryCatch(
            lmerTest::lmer(fml, data = d, REML = CFG$model$reml %||% FALSE),
            error = function(e) { md_warn(CFG, "  [", nm, "] lmer 拟合失败: ",
                                          conditionMessage(e)); NULL })
          if (!is.null(fit)) {
            co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
            if (!is.null(co)) {
              for (tm in rownames(co)) {
                coef_rows[[length(coef_rows) + 1L]] <- data.frame(
                  subset = nm, n = nrow(d), term = tm,
                  estimate = co[tm, "Estimate"],
                  se = co[tm, "Std. Error"],
                  statistic = co[tm, "t value"],
                  df = if ("df" %in% colnames(co)) co[tm, "df"] else NA_real_,
                  p = co[tm, if ("Pr(>|t|)" %in% colnames(co)) "Pr(>|t|)" else "Pr(>|z|)"],
                  stringsAsFactors = FALSE)
              }
            }
            r2 <- tryCatch(performance::r2(fit), error = function(e) NULL)
            fit_rows[[length(fit_rows) + 1L]] <- data.frame(
              subset = nm, n = nrow(d), model = "lmm",
              singular = isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) FALSE)),
              R2m = if (!is.null(r2)) r2$R2_marginal else NA_real_,
              R2c = if (!is.null(r2)) r2$R2_conditional else NA_real_,
              AIC = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
              guarded = paste(setdiff(unique(c(focus, covs)), rhs_terms), collapse = ","),
              stringsAsFactors = FALSE)
          }
        }
      }
    }
  }

  coef_df <- if (length(coef_rows)) do.call(rbind, coef_rows) else
    data.frame(subset = character(0), n = integer(0), term = character(0),
               estimate = numeric(0), se = numeric(0), statistic = numeric(0),
               df = numeric(0), p = numeric(0), stringsAsFactors = FALSE)
  fit_df <- if (length(fit_rows)) do.call(rbind, fit_rows) else
    data.frame(subset = character(0), n = integer(0), model = character(0),
               singular = logical(0), R2m = numeric(0), R2c = numeric(0),
               AIC = numeric(0), guarded = character(0), stringsAsFactors = FALSE)
  skip_df <- if (length(skip_rows)) do.call(rbind, skip_rows) else
    data.frame(subset = character(0), n = integer(0), reason = character(0),
               stringsAsFactors = FALSE)

  if (nrow(coef_df) == 0) {
    md_warn(CFG, "所有子集都没能拟合出模型")
    md_write_xlsx(CFG, "53_sensitivity", list(Coefficients = coef_df, Fit = fit_df,
                                              Skipped = skip_df,
                                              Note = data.frame(
                                                note = "所有子集都未拟合成功，见 Skipped sheet",
                                                stringsAsFactors = FALSE)))
    return(md_result(module, "PARTIAL", reason = "所有子集拟合失败"))
  }

  # P 校正（同一模型内的多个路径 / 跨子集同一条路径都不该混着校正）
  coef_df$p_adj <- ave(coef_df$p, coef_df$term, FUN = function(x) md_p_adjust(x, "BH"))

  # ---- 3. 并列对比表（每条路径 × 每个子集 一行）----
  focus_terms <- unique(coef_df$term[coef_df$term != "(Intercept)"])
  wide <- do.call(rbind, lapply(focus_terms, function(tm) {
    sub <- coef_df[coef_df$term == tm, , drop = FALSE]
    full <- sub[sub$subset == "Full", , drop = FALSE]
    data.frame(
      term = tm,
      full_estimate = if (nrow(full)) full$estimate[1] else NA_real_,
      full_p = if (nrow(full)) full$p[1] else NA_real_,
      n_subsets = nrow(sub) - nrow(full),
      estimate_range = if (nrow(sub) > 1)
        sprintf("%.3f ~ %.3f", min(sub$estimate), max(sub$estimate)) else "—",
      p_range = if (nrow(sub) > 1)
        paste(md_fmt_p(min(sub$p)), " ~ ", md_fmt_p(max(sub$p)), sep = "") else "—",
      sign_consistent = if (nrow(sub) > 1)
        all(sign(sub$estimate) == sign(sub$estimate[1])) else TRUE,
      stringsAsFactors = FALSE)
  }))
  wide$verdict <- ifelse(wide$sign_consistent & wide$p_range != "—",
                         "方向一致（需看 P 是否跨过 0.05）", "方向不一致 —— 结论脆弱")

  # ---- 4. 图：剔组系数并列森林图（对 原脚本「剔除南沙敏感性」段的
  #         make_coef_panel 并列呈现方式：同一批系数、按数据集并列成行）
  #   画法口径与 原脚本 §C/D make_coef_panel 一致：
  #     误差棒（height 0.18）+ 实心白描边圆点 shape 21；参考线 0 灰色虚线；
  #     方向×显著性决定点色（正 / 负 / 不显著 grey60）；y 轴标签加粗。
  #   §A 的站点级共线守卫保留，故同一 term 在不同子集里可能只剩部分行。
  pt <- md_pt(CFG)
  plot_df <- coef_df[coef_df$term != "(Intercept)", , drop = FALSE]
  if (!nrow(plot_df)) plot_df <- coef_df
  # 子集顺序：Full 在最前，其余按名称（原脚本「全量 -> 剔除后」的并列次序）
  sub_lv <- unique(as.character(plot_df$subset))
  sub_lv <- c(intersect("Full", sub_lv), sort(setdiff(sub_lv, "Full")))
  plot_df$subset <- factor(as.character(plot_df$subset), levels = sub_lv)
  # y 轴 = term × subset（同一 term 的各子集并列成相邻的行）
  term_lv <- unique(as.character(plot_df$term))
  plot_df$row <- factor(paste(as.character(plot_df$term), as.character(plot_df$subset), sep = "\n"),
                        levels = as.vector(t(outer(term_lv, sub_lv, function(a, b)
                          paste(a, b, sep = "\n")))))
  plot_df <- plot_df[!is.na(plot_df$row), , drop = FALSE]
  plot_df$lo <- plot_df$estimate - 1.96 * plot_df$se
  plot_df$hi <- plot_df$estimate + 1.96 * plot_df$se
  plot_df$Col <- ifelse(is.na(plot_df$p), "grey60",
                        ifelse(plot_df$p < 0.05,
                               ifelse(plot_df$estimate > 0, md_pal_cont(CFG, "diverging")[3],
                                      md_pal_cont(CFG, "diverging")[1]),
                               "grey60"))
  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$estimate, y = .data$row)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey35",
                        linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                           orientation = "y", height = 0.18, linewidth = MD_LW$line,
                           colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$Col), size = MD_PT_SIZE$scatter * 2.0,
                        shape = 21, colour = "white", stroke = 0.3) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(clip = "off") +
    # 口径说明：本模块的估计量**不是**标准化系数（z 化在 50 里做），如实标注
    ggplot2::labs(x = "Coefficient (raw scale)", y = NULL, tag = "A") +
    theme_md(CFG = CFG) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold", size = pt$tick))

  p <- md_compose(CFG, list(A = p1), design = "AA\nAA", legend = "none")
  md_save_figure(CFG, p, "53_sensitivity", height = 120)

  tp <- md_write_xlsx(CFG, "53_sensitivity", list(
    Robustness_summary = wide,
    Coefficients = coef_df,
    Model_fit = fit_df,
    Skipped = skip_df,
    Design = data.frame(
      item = c("drop_levels", "drop_sites", "drop_samples", "target", "max_subsets", "min_n"),
      value = c(paste(names(drop_levels), vapply(drop_levels, paste, character(1), collapse = "+"),
                      sep = "=", collapse = "; "),
                paste(drop_sites, collapse = ",") %||% "—",
                paste(drop_samples, collapse = ",") %||% "—",
                paste(target, collapse = "+"), max_subsets, min_n),
      stringsAsFactors = FALSE),
    Note = data.frame(note = c(
      "每个子集重跑同一个模型；Full = 全部样本。",
      "sign_consistent = 各子集系数方向是否与 Full 一致；方向翻转说明该结论由被剔除的组驱动。",
      "P 值跨子集需各自解读，不要跨子集做多重校正。",
      "被共线守卫剔除的变量见 Model_fit 的 guarded 列。"),
      stringsAsFactors = FALSE)),
    meta_extra = list(n_subsets = nrow(fit_df), n_terms = length(focus_terms)))

  md_result(module, "OK",
            tables = basename(tp), figures = "figures/53_sensitivity.pdf",
            notes = paste0("比较了 ", nrow(fit_df), " 个子集 / ", length(focus_terms), " 条路径"))
}
