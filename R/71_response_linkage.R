# =============================================================================
# R/71_response_linkage.R ---- 响应变量与群落指标 / 群落距离的关联
# -----------------------------------------------------------------------------
# 契约见 _build-notes/micro-div-pipeline/MODULE_SPEC.md
#
# 按 CFG$response_models$linkage 的子集（alpha / ordination / pairwise）分别做：
#   alpha       响应 vs α 多样性（CFG$response_models$alpha_focus）
#               单回归 + 偏回归 + LMM / GAMM（协变量 CFG$model$covariates，随机效应 CFG$model$random）
#   ordination  响应 vs 排序轴（CFG$response_models$ordination_axes）的条件效应
#               （全变量 z 化网格 + 95% 置信带 + 反变换回原始单位）
#   pairwise    成对级 |delta响应| vs 群落距离（CFG$response_models$pair_metric）
#               MRM（对象级置换）+ Mantel（秩）+ md_boot_rho() 站点对整簇 bootstrap
#               + 饱和诊断（mean / sd / frac == 1）—— **只出表，不出图**
#               （图件版式属 43_mantel_mrm，对应论文 Fig.S12B）
#   图 = 两类各一个面板拼一张（A 响应 vs α、B PCoA1 条件效应）；表 = 每类一张 + 注释
#
# 数据来源（**只走访问器**）
#   md_analysis_with_derived() / md_derived()            ---- 40 号产出（缺则本地兜底）
#   md_phylo_cache()                                     ---- 21 号产出（PD / NRI / NTI）
#   md_dist_get(CFG, metric)                             ---- 13 号产出（缺则本地重算）
#
# 所有随机过程前调用 md_seed(CFG, offset)。
# =============================================================================

md_response_linkage <- function(CFG) {
  module <- "71"

  if (!md_precheck(CFG, module, need_response = TRUE))
    return(md_result(module, "SKIP", "response$enabled = FALSE 或前置条件不满足"))
  if (!md_need(CFG, c("model", "mantel"), module))
    return(md_result(module, "SKIP", "缺少依赖包（lme4 / mgcv / ecodist / geosphere）"))

  d <- tryCatch(md_analysis_with_derived(CFG), error = function(e) NULL)
  if (is.null(d)) {
    md_skip(CFG, module, "取不到分析总表（cache/05_core.rds 缺失？）；71 需要 13 与 40 的产物")
    return(md_result(module, "SKIP", "取不到分析总表"))
  }

  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  smp  <- CFG$data$sample_col
  yv   <- CFG$response$var
  ylab <- CFG$response$label %||% yv
  yunit <- CFG$response$unit %||% ""
  if (!yv %in% names(d)) {
    md_skip(CFG, module, paste0("响应变量列 ", yv, " 不在数据中"))
    return(md_result(module, "SKIP", paste0("缺响应变量列 ", yv)))
  }

  y <- suppressWarnings(as.numeric(d[[yv]]))
  tr <- CFG$response$trans %||% "none"
  y <- switch(tr,
              log  = if (all(y[is.finite(y)] > 0)) log(y) else y,
              sqrt = if (all(y[is.finite(y)] >= 0)) sqrt(y) else y,
              z    = as.numeric(scale(y)),
              y)
  d$.y <- y
  if (sum(is.finite(d$.y)) < 5 || stats::sd(d$.y, na.rm = TRUE) <= 0) {
    md_skip(CFG, module, "响应变量有效观测不足或零方差")
    return(md_result(module, "SKIP", "响应变量不可用"))
  }
  y_axis <- if (identical(tr, "none")) {
    if (nzchar(yunit)) paste0(ylab, " (", yunit, ")") else ylab
  } else paste0(ylab, " (", tr, "-transformed)")

  ids_all <- as.character(d[[smp]])
  linkage <- CFG$response_models$linkage %||% c("alpha", "ordination", "pairwise")

  # ---- 社区指标补齐 ----
  # 总表（md_analysis_with_derived）只含 meta + env + 派生变量，不含 Shannon / PD /
  # 排序轴。此处按两类关联的需要补齐：α 指数本地现算，PD 取 21 号缓存；
  # 排序轴在 ordination 分支里由 13 号距离（或本地重算）现做 PCoA。
  local_alpha <- function(cfg, ids, idx) {
    m <- tryCatch(md_core_asv(cfg, "alpha"), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    mm <- t(m)[ids, , drop = FALSE]
    out <- list()
    for (k in idx) {
      v <- switch(k,
                  Shannon  = tryCatch(vegan::diversity(mm, index = "shannon"), error = function(e) NULL),
                  Simpson  = tryCatch(vegan::diversity(mm, index = "simpson"), error = function(e) NULL),
                  Observed = rowSums(mm > 0),
                  Richness = rowSums(mm > 0),
                  Pielou   = {
                    ob <- rowSums(mm > 0)
                    sh <- tryCatch(vegan::diversity(mm, index = "shannon"), error = function(e) NULL)
                    if (is.null(sh)) NULL else sh / log(pmax(ob, 2))
                  },
                  NULL)
      if (!is.null(v)) out[[k]] <- as.numeric(v)
    }
    if (!length(out)) return(NULL)
    data.frame(sample_id = ids, out, check.names = FALSE, stringsAsFactors = FALSE)
  }
  local_fill <- function(df, tab, cols, why) {
    tab <- tab[match(as.character(df[[smp]]), as.character(tab$sample_id)), , drop = FALSE]
    rownames(tab) <- NULL; rownames(df) <- NULL
    df <- cbind(df, tab[, cols, drop = FALSE])
    add_note("已补齐社区指标 ", paste(cols, collapse = ", "), "（", why, "）。")
    df
  }
  a_need <- intersect(CFG$response_models$alpha_focus %||% character(0),
                      c("Shannon", "Simpson", "Observed", "Richness", "Pielou"))
  a_need <- setdiff(a_need, names(d))
  if (length(a_need)) {
    a <- local_alpha(CFG, ids_all, a_need)
    if (!is.null(a)) d <- local_fill(d, a, intersect(names(a), a_need),
                                     "由 md_core_asv(\"alpha\") 现算")
    else add_note("α 指标 ", paste(a_need, collapse = ", "),
                  " 无法现算（md_core_asv(\"alpha\") 不可用）。")
  }
  ph_need <- intersect(CFG$response_models$alpha_focus %||% character(0),
                       c("PD", "MPD", "MNTD", "NRI", "NTI"))
  ph_need <- setdiff(ph_need, names(d))
  if (length(ph_need)) {
    ph <- tryCatch(md_phylo_cache(CFG), error = function(e) NULL)
    if (is.data.frame(ph) && smp %in% names(ph)) {
      ex <- intersect(names(ph), ph_need)
      if (length(ex)) {
        d <- local_fill(d, ph[, c(smp, ex), drop = FALSE], ex,
                        "来源 21_phylo_diversity 缓存")
      } else add_note("21 号缓存里没有 ", paste(ph_need, collapse = ", "), "。")
    } else {
      add_note("取不到 21 号产物（md_phylo_cache），",
               paste(ph_need, collapse = ", "), " 不可用。")
    }
  }
  ids_all <- as.character(d[[smp]])

  # ---- 分组因子（图例用；水平过多则不上色）----
  # 上色因子照抄原脚本 §附图 与 原脚本 §附图
  # 的散点一律 `colour = Vegetation` + `scale_colour_manual(values = cols_veg)`，
  # 即**第二个设计因子**（本包 design$factors = c("Zone","Vegetation")）。
  .fac_all <- CFG$design$factors %||% character(0)
  gvar <- if (length(.fac_all) >= 2) .fac_all[2] else
    if (length(.fac_all) == 1) .fac_all[1] else NULL
  if (is.null(gvar) || !gvar %in% names(d)) gvar <- NULL
  if (!is.null(gvar) && length(unique(d[[gvar]])) > 6) gvar <- NULL
  gcols <- if (!is.null(gvar))
    md_pal(CFG, (CFG$design$ordered_levels[[gvar]] %||% unique(as.character(d[[gvar]]))), gvar) else NULL
  if (!is.null(gvar)) d[[gvar]] <- factor(as.character(d[[gvar]]),
                                          levels = CFG$design$ordered_levels[[gvar]] %||% NULL)

  pt <- md_pt(CFG)
  panels <- list(); tabs <- list()

  # ===========================================================================
  # 公共小工具（local_ 前缀，不改共享库）
  # ===========================================================================
  local_resid_ok <- function(df, vars) {
    if (length(vars) == 0) return(logical(0))
    vapply(vars, function(v) v %in% names(df) && is.numeric(df[[v]]) &&
             sum(is.finite(df[[v]])) >= 4, logical(1))
  }

  # 站点内无变异 => 与站点随机截距共线
  local_site_level <- function(df, vars) {
    sv <- CFG$design$site_var
    if (is.null(sv) || !sv %in% names(df) || length(vars) == 0) return(setNames(rep(FALSE, length(vars)), vars))
    out <- vapply(vars, function(v) {
      x <- df[[v]]
      if (!is.numeric(x)) x <- as.numeric(factor(x))
      u <- tapply(x, df[[sv]], function(z) length(unique(z[is.finite(z)])))
      isTRUE(max(u, na.rm = TRUE) == 1)
    }, logical(1))
    out
  }

  local_re_term <- function(df) {
    rs <- CFG$model$random
    if (is.null(rs) || !nzchar(rs)) return(NULL)
    m <- regmatches(rs, regexec("\\|\\s*([A-Za-z0-9_.]+)", rs))[[1]]
    if (length(m) < 2) return(NULL)
    v <- m[2]
    if (!v %in% names(df)) return(NULL)
    if (length(unique(stats::na.omit(df[[v]]))) < 3) return(NULL)
    v
  }

  local_stats <- function(fit, term) {
    na <- c(est = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_)
    if (is.null(fit)) return(na)
    cf <- tryCatch(as.data.frame(summary(fit)$coefficients), error = function(e) NULL)
    if (is.null(cf) || !term %in% rownames(cf)) return(na)
    pcol <- intersect(c("Pr(>|t|)", "p.value", "Pr(>|z|)", "P"), names(cf))[1]
    scol <- intersect(c("t value", "z value"), names(cf))[1]
    est  <- as.numeric(cf[term, "Estimate"])
    se   <- as.numeric(cf[term, intersect(c("Std. Error", "Std.Error"), names(cf))[1]])
    st   <- if (is.na(scol)) NA_real_ else as.numeric(cf[term, scol])
    if (is.na(pcol)) {
      pcol <- "Pr(>|t|)"
      cf[[pcol]] <- 2 * stats::pnorm(-abs(st))
    }
    p <- as.numeric(cf[term, pcol])
    c(est = est, se = se, stat = st, p = p)
  }

  local_coef <- function(fit) {
    if (is.null(fit)) return(c("(Intercept)" = NA_real_))
    if (inherits(fit, "merMod")) return(lme4::fixef(fit))
    stats::coef(fit)
  }
  local_vcov <- function(fit) {
    if (is.null(fit)) return(matrix(NA_real_))
    if (inherits(fit, "merMod")) return(as.matrix(stats::vcov(fit)))
    stats::vcov(fit)
  }
  local_back <- function(z, cen, sc) as.numeric(z) * sc + cen

  #' 条件效应：z 化网格 -> 点估计 + 95% 置信带 -> 反变换回原始单位
  #' 返回 list(curve, d, n, beta, SE, stat, p_lmm, edf, p_smooth, src, covs)
  local_cond_effect <- function(df, yvar, xvar, covars, allow_smooth = TRUE,
                                n_grid = NULL, edf_thr = NULL, re = NULL) {
    n_grid <- n_grid %||% CFG$model$grid_n %||% 200
    edf_thr <- edf_thr %||% CFG$model$edf_nonlin %||% 1.1
    vars <- unique(c(yvar, xvar, covars))
    need <- unique(c(vars, re))
    dd <- df[stats::complete.cases(df[, need, drop = FALSE]), , drop = FALSE]
    if (nrow(dd) < max(6, length(covars) + 3)) return(NULL)
    if (stats::sd(dd[[xvar]], na.rm = TRUE) <= 0 || stats::sd(dd[[yvar]], na.rm = TRUE) <= 0)
      return(NULL)
    if (nrow(dd) < 3 * length(vars) + 2) {
      # 样本量不足以支撑这么多参数 -> 只留最强的协变量（按 |Spearman| 排序）
      if (length(covars) > 0) {
        rs <- vapply(covars, function(v) {
          s <- suppressWarnings(stats::cor(dd[[v]], dd[[yvar]], method = "spearman",
                                           use = "complete.obs"))
          if (is.na(s)) 0 else abs(s)
        }, numeric(1))
        keep <- names(sort(rs, decreasing = TRUE))[seq_len(min(1, length(covars)))]
        covars <- keep
        vars <- unique(c(yvar, xvar, covars))
      }
    }
    z <- dd; cen <- list(); sc <- list()
    for (v in vars) {
      cen[[v]] <- mean(dd[[v]], na.rm = TRUE)
      sc[[v]]  <- stats::sd(dd[[v]], na.rm = TRUE)
      if (!is.finite(sc[[v]]) || sc[[v]] == 0) return(NULL)
      z[[v]] <- (dd[[v]] - cen[[v]]) / sc[[v]]
    }
    fx <- c(xvar, covars)
    re_term <- if (is.null(re)) "" else paste0(" + (1|", re, ")")
    f_lmm <- stats::as.formula(paste0(yvar, " ~ ", paste0("`", fx, "`", collapse = " + "), re_term))
    fit <- tryCatch({
      if (is.null(re)) stats::lm(f_lmm, data = z)
      else lmerTest::lmer(f_lmm, data = z, REML = FALSE)
    }, error = function(e) NULL)
    if (is.null(fit)) {
      fit <- tryCatch(stats::lm(stats::as.formula(paste0(yvar, " ~ ",
                                                         paste0("`", fx, "`", collapse = " + "))),
                                data = z), error = function(e) NULL)
      if (!is.null(fit)) re <- NULL
    }
    if (is.null(fit)) return(NULL)
    st <- local_stats(fit, xvar)

    # ---- 平滑项诊断 ----
    edf_v <- NA_real_; p_sm <- NA_real_; g <- NULL
    if (allow_smooth && length(unique(z[[xvar]])) >= 5 && nrow(z) >= 10) {
      k_s <- max(3, min(5, floor(length(unique(z[[xvar]])) / 2)))
      rhs <- c(paste0("s(`", xvar, "`, k = ", k_s, ")"),
               paste0("`", setdiff(fx, xvar), "`"))
      if (is.null(re)) {
        g <- tryCatch(mgcv::gam(stats::as.formula(paste0(yvar, " ~ ", paste(rhs, collapse = " + "))),
                                data = z, method = "REML"), error = function(e) NULL)
        if (!is.null(g)) {
          s <- summary(g)$s.table
          rn <- rownames(s)
          hit <- grep(paste0("s\\(", xvar), rn)
          if (length(hit)) { edf_v <- s[hit[1], "edf"]; p_sm <- s[hit[1], "p-value"] }
        }
      } else {
        gm <- tryCatch(mgcv::gamm(
          stats::as.formula(paste0(yvar, " ~ ", paste(rhs, collapse = " + "))),
          random = stats::as.formula(paste0("~ 1 | ", re)), data = z, method = "REML"),
          error = function(e) NULL)
        if (!is.null(gm)) {
          g <- gm$gam
          s <- summary(g)$s.table
          rn <- rownames(s)
          hit <- grep(paste0("s\\(", xvar), rn)
          if (length(hit)) { edf_v <- s[hit[1], "edf"]; p_sm <- s[hit[1], "p-value"] }
        }
      }
    }
    smooth <- is.finite(edf_v) && edf_v > edf_thr

    # ---- 预测网格（z 空间，其余变量取 0 = 均值）----
    xs <- seq(min(z[[xvar]]), max(z[[xvar]]), length.out = n_grid)
    nd <- data.frame(xs); names(nd) <- xvar
    for (v in setdiff(fx, xvar)) nd[[v]] <- 0
    nd <- nd[, fx, drop = FALSE]

    if (smooth && !is.null(g)) {
      pr <- tryCatch(stats::predict(g, newdata = nd, se.fit = TRUE, type = "response"),
                     error = function(e) NULL)
      if (is.null(pr)) { smooth <- FALSE } else {
        cv <- data.frame(x_raw = local_back(xs, cen[[xvar]], sc[[xvar]]),
                         fit = local_back(pr$fit, cen[[yvar]], sc[[yvar]]),
                         lower = local_back(pr$fit - 1.96 * pr$se.fit, cen[[yvar]], sc[[yvar]]),
                         upper = local_back(pr$fit + 1.96 * pr$se.fit, cen[[yvar]], sc[[yvar]]))
      }
    }
    if (!smooth) {
      X <- tryCatch(stats::model.matrix(
        stats::as.formula(paste0("~ ", paste0("`", fx, "`", collapse = " + "))), data = nd),
        error = function(e) NULL)
      if (is.null(X)) return(NULL)
      b <- local_coef(fit); V <- local_vcov(fit)
      cn <- intersect(colnames(X), names(b))
      X <- X[, cn, drop = FALSE]
      fv <- as.vector(X %*% b[cn])
      se <- sqrt(pmax(0, rowSums((X %*% V[cn, cn, drop = FALSE]) * X)))
      cv <- data.frame(x_raw = local_back(xs, cen[[xvar]], sc[[xvar]]),
                       fit = local_back(fv, cen[[yvar]], sc[[yvar]]),
                       lower = local_back(fv - 1.96 * se, cen[[yvar]], sc[[yvar]]),
                       upper = local_back(fv + 1.96 * se, cen[[yvar]], sc[[yvar]]))
    }
    list(curve = cv, d = dd, n = nrow(dd),
         beta = unname(st["est"]), SE = unname(st["se"]), stat = unname(st["stat"]),
         p_lmm = unname(st["p"]), edf = edf_v, p_smooth = p_sm,
         src = if (smooth) "GAMM" else if (is.null(re)) "OLS" else "LMM",
         covs = covars, re = re %||% "none")
  }

  local_marginal <- function(df, xvar) {
    dd <- df[stats::complete.cases(df[, c(yv, xvar)]), , drop = FALSE]
    if (nrow(dd) < 4) return(NULL)
    fit <- tryCatch(stats::lm(stats::as.formula(paste0("`", yv, "` ~ `", xvar, "`")), data = dd),
                    error = function(e) NULL)
    sp  <- tryCatch(suppressWarnings(stats::cor.test(dd[[xvar]], dd[[yv]], method = "spearman")),
                    error = function(e) NULL)
    st <- local_stats(fit, xvar)
    data.frame(predictor = xvar, n = nrow(dd),
               slope = unname(st["est"]), SE = unname(st["se"]), t = unname(st["stat"]),
               p = unname(st["p"]),
               r2 = if (is.null(fit)) NA_real_ else summary(fit)$r.squared,
               rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
               p_spearman = if (is.null(sp)) NA_real_ else sp$p.value,
               stringsAsFactors = FALSE)
  }

  local_annot <- function(label, edf, p_lmm, p_sm, src) {
    l <- c(label, paste0("edf = ", if (is.na(edf)) "n.i." else sprintf("%.2f", edf)),
           paste0(src, " P = ", md_fmt_p(p_lmm)))
    if (!is.na(edf) && !is.na(p_sm) && edf > (CFG$model$edf_nonlin %||% 1.1))
      l <- c(l, paste0("GAMM P = ", md_fmt_p(p_sm)))
    paste(l, collapse = "\n")
  }

  # ===========================================================================
  # 1. alpha ---- 响应 vs α 多样性
  # ===========================================================================
  alpha_tab <- NULL; partial_tab <- NULL
  if ("alpha" %in% linkage) {
    a_vars <- CFG$response_models$alpha_focus %||% c("Shannon")
    a_vars <- a_vars[a_vars %in% names(d)]
    a_vars <- a_vars[vapply(a_vars, function(v) v %in% names(d) && is.numeric(d[[v]]) &&
                              sum(is.finite(d[[v]])) >= 4 &&
                              stats::sd(d[[v]], na.rm = TRUE) > 0, logical(1))]
    if (length(a_vars) == 0) {
      add_note("linkage 含 alpha，但没有任何可用的 α 指标（alpha_focus 均不在数据中）。")
    } else {
      covs <- CFG$model$covariates %||% character(0)
      drv  <- CFG$model$derived %||% character(0)
      covs_all <- unique(c(covs, drv))
      covs_all <- covs_all[local_resid_ok(d, covs_all)]
      covs_all <- setdiff(covs_all, a_vars)
      site_all <- local_site_level(d, unique(c(a_vars, covs_all)))
      site_of  <- function(x) {
        i <- match(x, names(site_all))
        if (length(i) == 0 || is.na(i)) return(FALSE)
        isTRUE(unname(site_all[i]))
      }
      if (length(covs_all)) {
        bad <- covs_all[vapply(covs_all, site_of, logical(1))]
        if (length(bad))
          add_note("以下协变量在站点内无变异（与站点随机截距共线），未纳入 LMM：",
                   paste(bad, collapse = ", "), "。")
      }
      re <- local_re_term(d)

      alpha_tab <- do.call(rbind, lapply(a_vars, function(v) local_marginal(d, v)))
      alpha_tab$p_BH <- md_p_adjust(alpha_tab$p, CFG$model$p_adjust %||% "BH")

      rows <- list()
      for (v in a_vars) {
        cv <- setdiff(covs_all, v)
        cv <- cv[!vapply(cv, site_of, logical(1))]
        f <- local_cond_effect(d, ".y", v, cv,
                               allow_smooth = !site_of(v), re = re)
        if (is.null(f)) {
          rows[[length(rows) + 1L]] <- data.frame(
            metric = v, n = NA_integer_, beta = NA_real_, SE = NA_real_,
            statistic = NA_real_, p = NA_real_, edf = NA_real_, p_smooth = NA_real_,
            model = "not estimable", random = re %||% "none",
            covariates = paste(cv, collapse = " + "), stringsAsFactors = FALSE)
          next
        }
        rows[[length(rows) + 1L]] <- data.frame(
          metric = v, n = f$n, beta = f$beta, SE = f$SE, statistic = f$stat,
          p = f$p_lmm, edf = f$edf, p_smooth = f$p_smooth,
          model = f$src, random = f$re, covariates = paste(f$covs, collapse = " + "),
          stringsAsFactors = FALSE)
      }
      partial_tab <- do.call(rbind, rows)
      if (!is.null(partial_tab)) partial_tab$p_BH <- md_p_adjust(partial_tab$p,
                                                                 CFG$model$p_adjust %||% "BH")

      # 面板 A：取 alpha_focus 第一项画条件效应
      v1 <- a_vars[1]
      cv1 <- setdiff(covs_all, v1); cv1 <- cv1[!vapply(cv1, site_of, logical(1))]
      f1 <- local_cond_effect(d, ".y", v1, cv1, allow_smooth = TRUE, re = re)
      if (!is.null(f1)) {
        gA <- ggplot2::ggplot()
        if (!is.null(gvar))
          gA <- gA + ggplot2::geom_point(data = f1$d,
                                         ggplot2::aes(x = .data[[v1]], y = .data$.y,
                                                      colour = .data[[gvar]]),
                                         size = MD_PT_SIZE$point,
                                         alpha = MD_PT_SIZE$point_alpha)
        else
          gA <- gA + ggplot2::geom_point(data = f1$d,
                                         ggplot2::aes(x = .data[[v1]], y = .data$.y),
                                         size = MD_PT_SIZE$point,
                                         alpha = MD_PT_SIZE$point_alpha, colour = "grey30")
        # 曲线/置信带沿用同一套样式 §附图：fill = alpha(REG_LINE_COL, REG_RIBBON_A)
        gA <- gA +
          ggplot2::geom_ribbon(data = f1$curve,
                               ggplot2::aes(x = .data$x_raw, ymin = .data$lower,
                                            ymax = .data$upper),
                               fill = MD_REG$col, alpha = MD_REG$ribbon_alpha,
                               colour = NA) +
          ggplot2::geom_line(data = f1$curve,
                             ggplot2::aes(x = .data$x_raw, y = .data$fit),
                             colour = MD_REG$col, linewidth = MD_REG$lw)
        if (!is.null(gvar)) gA <- gA + ggplot2::scale_colour_manual(values = gcols, name = gvar)
        gA <- gA +
          # 标注锚在面板**左上**内侧（同 附图 的 -Inf / hjust = -0.08 / vjust = 1.3）
          ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.3,
                            size = pt2gg(pt$ann), colour = "grey25", lineheight = 0.9,
                            label = local_annot(v1, f1$edf, f1$p_lmm, f1$p_smooth, f1$src)) +
          ggplot2::labs(x = v1, y = y_axis, tag = "A") +
          theme_md(CFG = CFG) +
          ggplot2::theme(legend.position = if (is.null(gvar)) "none" else "bottom",
                         legend.key.size = grid::unit(0.28, "cm"))
        panels$A <- gA
      } else {
        add_note("α 面板的条件效应不可估计（样本量不足或零方差）。")
      }
    }
  }

  # ===========================================================================
  # 2. ordination ---- 响应 vs 排序轴的条件效应
  # ===========================================================================
  pcoa_axis_tab <- NULL; pcoa_curve <- NULL
  if ("ordination" %in% linkage) {
    axes <- CFG$response_models$ordination_axes %||% c("PCoA1")
    # 排序轴来源：优先 13 号缓存的距离矩阵；缺则本地用 md_core_asv 重算
    need_compute <- TRUE
    metric <- CFG$response_models$pair_metric %||% "bray"
    dmat <- tryCatch(md_dist_get(CFG, metric), error = function(e) NULL)
    dsrc <- paste0("13_distances cache (", metric, ")")
    if (is.null(dmat)) {
      md_warn(CFG, "[71] 取不到 13 号距离产物（cache/13_dist.rds），本地重算 ", metric, " 距离")
      add_note("取不到 13 号产物（md_dist_get），排序轴由本地 vegdist(", metric,
               ") + ape::pcoa 重算；若与 13/14 的轴不一致，请先运行模块 13。")
      m <- tryCatch(md_core_asv(CFG, "distances"), error = function(e) NULL)
      if (!is.null(m)) {
        mm <- t(m)                                   # 样本 x ASV
        mm <- mm[ids_all, , drop = FALSE]
        dmat <- tryCatch(vegan::vegdist(mm, method = metric), error = function(e) NULL)
        dsrc <- paste0("local recompute (", metric, ")")
      }
    }
    if (is.null(dmat)) {
      add_note("排序轴不可用：既无 13 号距离产物，本地重算也失败。ordination 面板跳过。")
    } else {
      pco <- tryCatch(ape::pcoa(dmat, correction = "cailliez"), error = function(e) NULL)
      if (is.null(pco)) {
        add_note("ape::pcoa 失败，ordination 面板跳过。")
      } else {
        eig <- pco$values$Corr_eig
        if (is.null(eig) || all(is.na(eig))) eig <- pco$values$Eigenvalues
        pct <- 100 * eig / sum(eig[eig > 0])
        ax_all <- as.data.frame(pco$vectors, stringsAsFactors = FALSE)
        ax_all[[smp]] <- rownames(ax_all)
        nm <- paste0("PCoA", seq_along(pct))
        names(ax_all) <- c(nm, smp)
        ax_sel <- intersect(axes, nm)
        if (length(ax_sel) == 0) ax_sel <- nm[1]
        ax_use <- ax_all[, c(smp, ax_sel), drop = FALSE]
        d <- merge(d, ax_use, by = smp, all.x = TRUE, sort = FALSE)
        d <- d[match(ids_all, as.character(d[[smp]])), , drop = FALSE]
        add_note("排序轴来源 = ", dsrc, "；轴 1 解释 ", sprintf("%.1f%%", pct[1]),
                 "。")

        covs_all <- unique(c(CFG$model$covariates %||% character(0),
                             CFG$model$derived %||% character(0)))
        covs_all <- covs_all[local_resid_ok(d, covs_all)]
        covs_all <- setdiff(covs_all, ax_sel)
        site_all <- local_site_level(d, unique(c(ax_sel, covs_all)))
        site_of  <- function(x) {
          i <- match(x, names(site_all))
          if (length(i) == 0 || is.na(i)) return(FALSE)
          isTRUE(unname(site_all[i]))
        }
        if (length(covs_all)) {
          bad <- covs_all[vapply(covs_all, site_of, logical(1))]
          if (length(bad))
            add_note("以下协变量在站点内无变异，未纳入排序轴模型：",
                     paste(bad, collapse = ", "), "。")
        }
        re <- local_re_term(d)

        rows <- list()
        for (ax in ax_sel) {
          cv <- setdiff(covs_all, ax); cv <- cv[!vapply(cv, site_of, logical(1))]
          f <- local_cond_effect(d, ".y", ax, cv, allow_smooth = !site_of(ax), re = re)
          if (is.null(f)) {
            rows[[length(rows) + 1L]] <- data.frame(
              axis = ax, axis_var_explained_pct = pct[match(ax, nm)],
              n = NA_integer_, beta = NA_real_, SE = NA_real_, statistic = NA_real_,
              p = NA_real_, edf = NA_real_, p_smooth = NA_real_, model = "not estimable",
              random = re %||% "none", covariates = paste(cv, collapse = " + "),
              stringsAsFactors = FALSE)
            next
          }
          rows[[length(rows) + 1L]] <- data.frame(
            axis = ax, axis_var_explained_pct = pct[match(ax, nm)],
            n = f$n, beta = f$beta, SE = f$SE, statistic = f$stat, p = f$p_lmm,
            edf = f$edf, p_smooth = f$p_smooth, model = f$src, random = f$re,
            covariates = paste(f$covs, collapse = " + "), stringsAsFactors = FALSE)
        }
        pcoa_axis_tab <- do.call(rbind, rows)
        if (!is.null(pcoa_axis_tab))
          pcoa_axis_tab$p_BH <- md_p_adjust(pcoa_axis_tab$p, CFG$model$p_adjust %||% "BH")

        ax1 <- ax_sel[1]
        cv1 <- setdiff(covs_all, ax1); cv1 <- cv1[!vapply(cv1, site_of, logical(1))]
        f1 <- local_cond_effect(d, ".y", ax1, cv1, allow_smooth = TRUE, re = re)
        if (!is.null(f1)) {
          pcoa_curve <- f1
          gB <- ggplot2::ggplot()
          if (!is.null(gvar))
            gB <- gB + ggplot2::geom_point(data = f1$d,
                                           ggplot2::aes(x = .data[[ax1]], y = .data$.y,
                                                        colour = .data[[gvar]]),
                                           size = MD_PT_SIZE$point,
                                           alpha = MD_PT_SIZE$point_alpha)
          else
            gB <- gB + ggplot2::geom_point(data = f1$d,
                                           ggplot2::aes(x = .data[[ax1]], y = .data$.y),
                                           size = MD_PT_SIZE$point,
                                           alpha = MD_PT_SIZE$point_alpha, colour = "grey30")
          # 曲线/置信带沿用同一套样式 §附图A
          gB <- gB +
            ggplot2::geom_ribbon(data = f1$curve,
                                 ggplot2::aes(x = .data$x_raw, ymin = .data$lower,
                                              ymax = .data$upper),
                                 fill = MD_REG$col, alpha = MD_REG$ribbon_alpha,
                                 colour = NA) +
            ggplot2::geom_line(data = f1$curve,
                               ggplot2::aes(x = .data$x_raw, y = .data$fit),
                               colour = MD_REG$col, linewidth = MD_REG$lw)
          if (!is.null(gvar)) gB <- gB + ggplot2::scale_colour_manual(values = gcols, name = gvar)
          gB <- gB +
            # 标注锚在面板**右上**内侧（同 附图A 的 Inf / hjust = 1 / vjust = 1.3）
            ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1.3,
                              size = pt2gg(pt$ann), colour = "grey25", lineheight = 0.95,
                              label = local_annot(ax1, f1$edf, f1$p_lmm, f1$p_smooth, f1$src)) +
            ggplot2::labs(x = sprintf("%s (%.2f%% of community variation)",
                                      ax1, pct[match(ax1, nm)]),
                          y = NULL, tag = "B") +
            theme_md(CFG = CFG) +
            ggplot2::theme(legend.position = if (is.null(gvar)) "none" else "bottom",
                           legend.key.size = grid::unit(0.28, "cm"))
          panels$B <- gB
        } else {
          add_note("排序轴面板的条件效应不可估计。")
        }
      }
    }
  }

  # ===========================================================================
  # 3. pairwise ---- 成对级 |delta响应| vs 群落距离（**只出表，不出图**）
  #    图件（|delta响应| vs 群落距离，Fig.S12B）归 43_mantel_mrm；
  #    本模块出统计量 + 饱和诊断表，供表件与 43 引用。
  # ===========================================================================
  pair_tab <- NULL; boot_tab <- NULL; sat_tab <- NULL; mrm_tab <- NULL; mant_tab <- NULL
  if ("pairwise" %in% linkage) {
    metric <- CFG$response_models$pair_metric %||% "bray"
    dmat <- tryCatch(md_dist_get(CFG, metric), error = function(e) NULL)
    if (is.null(dmat)) {
      md_warn(CFG, "[71] 成对分析：取不到 13 号产物，本地重算 ", metric, " 距离")
      add_note("成对分析的距离矩阵由本地重算（", metric, "），未复用 13 号缓存。")
      m <- tryCatch(md_core_asv(CFG, "distances"), error = function(e) NULL)
      if (!is.null(m)) {
        mm <- t(m)[ids_all, , drop = FALSE]
        dmat <- tryCatch(vegan::vegdist(mm, method = metric), error = function(e) NULL)
      }
    }
    if (is.null(dmat)) {
      add_note("成对分析跳过：没有可用的群落距离矩阵。")
    } else {
      M <- as.matrix(dmat)
      ids2 <- intersect(ids_all, rownames(M))
      M <- M[ids2, ids2, drop = FALSE]
      if (length(ids2) < 5) {
        add_note("成对分析跳过：可对齐的样本不足。")
      } else {
        yv2 <- d$.y[match(ids2, ids_all)]
        yy  <- stats::setNames(yv2, ids2)
        lt <- which(lower.tri(M), arr.ind = TRUE)
        P <- data.frame(sample1 = rownames(M)[lt[, 1]],
                        sample2 = colnames(M)[lt[, 2]],
                        comm = M[lt], stringsAsFactors = FALSE)
        P$abs_dy <- abs(yy[P$sample1] - yy[P$sample2])
        P <- P[is.finite(P$abs_dy) & is.finite(P$comm), , drop = FALSE]
        sv <- CFG$design$site_var
        if (!is.null(sv) && sv %in% names(d)) {
          s_of <- stats::setNames(as.character(d[[sv]])[match(ids2, ids_all)], ids2)
          s1 <- s_of[P$sample1]; s2 <- s_of[P$sample2]
          P$site_pair <- ifelse(s1 <= s2, paste(s1, s2, sep = "|"), paste(s2, s1, sep = "|"))
        } else {
          P$site_pair <- paste(pmin(P$sample1, P$sample2), pmax(P$sample1, P$sample2), sep = "|")
        }
        geo <- tryCatch(md_geo_dist(md_core_env(CFG), ids2), error = function(e) NULL)
        if (!is.null(geo)) {
          G <- as.matrix(geo)[ids2, ids2, drop = FALSE]
          P$geo <- G[cbind(match(P$sample1, ids2), match(P$sample2, ids2))]
        }

        if (nrow(P) >= 5 && stats::sd(P$comm, na.rm = TRUE) > 0 &&
            stats::sd(P$abs_dy, na.rm = TRUE) > 0) {
          P$comm_z  <- as.numeric(scale(P$comm))
          P$abs_dy_z <- as.numeric(scale(P$abs_dy))
          if ("geo" %in% names(P)) P$geo_z <- as.numeric(scale(P$geo))

          nperm <- CFG$spatial$permutations %||% CFG$beta_tests$permutations %||% 999
          mrm_tab <- tryCatch(md_mrm(CFG, P, "abs_dy_z", "comm_z",
                                     nperm = nperm, label = paste0("|dResponse| ~ ", metric)),
                              error = function(e) NULL)
          mant_tab <- tryCatch(md_mantel(CFG, stats::as.dist(M), md_dist_mat(
            data.frame(sample1 = P$sample1, sample2 = P$sample2, value = P$abs_dy),
            ids2), method = "spearman", permutations = nperm,
            label = paste0(metric, " vs |dResponse| (rank)")), error = function(e) NULL)

          B <- as.integer(CFG$response_models$bootstrap_B %||% 1000)
          boot_tab <- tryCatch(md_boot_rho(CFG, P, x = "comm", y = "abs_dy",
                                           cluster = "site_pair", B = B,
                                           label = metric), error = function(e) NULL)

          sat_tab <- data.frame(
            metric = metric, n_pairs = nrow(P),
            mean = mean(P$comm), sd = stats::sd(P$comm),
            frac_eq1 = mean(P$comm >= 0.9999),
            rho = suppressWarnings(stats::cor(P$comm, P$abs_dy, method = "spearman")),
            r_pearson = suppressWarnings(stats::cor(P$comm, P$abs_dy)),
            R2_MRM = if (!is.null(mrm_tab) && "R2" %in% names(mrm_tab))
              mrm_tab$R2[1] else NA_real_,
            P_MRM = if (!is.null(mrm_tab) && "p" %in% names(mrm_tab)) {
              cr <- mrm_tab[!(as.character(mrm_tab$term) %in% c("Int", "Intercept")), , drop = FALSE]
              if (nrow(cr)) cr$p[1] else NA_real_
            } else NA_real_,
            P_Mantel = if (!is.null(mant_tab)) mant_tab$p[1] else NA_real_,
            rho_ci_lo = if (!is.null(boot_tab)) boot_tab$ci_lo[1] else NA_real_,
            rho_ci_hi = if (!is.null(boot_tab)) boot_tab$ci_hi[1] else NA_real_,
            stringsAsFactors = FALSE)
          sat_tab$saturated <- sat_tab$frac_eq1 >= 0.5

          notes <- c(notes, if (sat_tab$saturated)
            paste0("成对距离度量 ", metric, " 高度饱和（均值 ", sprintf("%.3f", sat_tab$mean),
                   "，", sprintf("%.0f%%", 100 * sat_tab$frac_eq1),
                   " 的样本对 = 1.0）-> 其 R2 上限被结构性压低，",
                   "低 R2 应读作'距离方差近乎耗尽'，不能读作'没有关系'。")
            else paste0("成对距离度量 ", metric, " 未饱和（均值 ",
                        sprintf("%.3f", sat_tab$mean), "，仅 ",
                        sprintf("%.0f%%", 100 * sat_tab$frac_eq1),
                        " 的样本对 = 1.0），低 R2 无'距离方差耗尽'这一免责理由。"))

          # 成对级 **只出表，不出图**：图件清单里「|delta响应| vs 群落距离」属
          # 43_mantel_mrm（论文 Fig.S12B），本模块只负责 A 响应 vs α、
          # B PCoA1 条件效应两块面板（见 references/06-figure-map.md）。
          # 下面只登记表格与饱和诊断，不再构造面板对象。
          pair_tab <- P
        } else {
          add_note("成对分析跳过：有效样本对不足或群落距离零方差。")
        }
      }
    }
  }

  # ===========================================================================
  # 2b. procrustes ---- 群落排序 vs 响应+环境排序的 Procrustes 一致性
  #     参考实现
  #     ⚠ 只在 linkage 含 "procrustes" 时执行；默认配置下不产生任何变化。
  # ===========================================================================
  proc_res <- NULL; proc_sum <- NULL
  if ("procrustes" %in% linkage) {
    pro <- CFG$response_models$procrustes %||% list()
    nax <- as.integer(pro$axes %||% 5L)
    nax <- max(2L, nax)
    pm  <- CFG$response_models$pair_metric %||% "bray"
    dmat <- tryCatch(md_dist_get(CFG, pm), error = function(e) NULL)
    if (is.null(dmat)) {
      m <- tryCatch(md_core_asv(CFG, "distances"), error = function(e) NULL)
      if (!is.null(m)) {
        mm <- t(m)[ids_all, , drop = FALSE]
        dmat <- tryCatch(vegan::vegdist(mm, method = pm), error = function(e) NULL)
        if (!is.null(dmat))
          add_note("Procrustes 的群落距离由本地重算（", pm, "），未复用 13 号缓存。")
      }
    }
    if (is.null(dmat)) {
      add_note("Procrustes 跳过：没有可用的群落距离矩阵。")
    } else {
      M <- as.matrix(dmat)
      ids2 <- intersect(ids_all, rownames(M))
      M <- M[ids2, ids2, drop = FALSE]
      if (length(ids2) < 5) {
        add_note("Procrustes 跳过：可对齐的样本不足。")
      } else {
        pc <- tryCatch(ape::pcoa(stats::as.dist(M), correction = "cailliez"),
                       error = function(e) NULL)
        Y <- if (!is.null(pc)) pc$vectors[, seq_len(min(nax, ncol(pc$vectors))), drop = FALSE] else NULL
        envv <- unique(c(yv, CFG$model$covariates %||% character(0),
                         CFG$model$derived %||% character(0)))
        envv <- envv[envv %in% names(d)]
        envv <- envv[vapply(envv, function(v) is.numeric(d[[v]]), logical(1))]
        E <- if (length(envv)) {
          em <- apply(as.matrix(d[match(ids2, ids_all), envv, drop = FALSE]), 2, as.numeric)
          rownames(em) <- ids2
          em[stats::complete.cases(em), , drop = FALSE]
        } else NULL
        if (is.null(Y) || is.null(E) || nrow(E) < 5 || ncol(E) < 1) {
          add_note("Procrustes 跳过：排序轴或响应/协变量矩阵不可用。")
        } else {
          keep <- intersect(rownames(E), ids2)
          X <- tryCatch(stats::prcomp(E[keep, , drop = FALSE],
                                      center = TRUE, scale. = TRUE)$x, error = function(e) NULL)
          if (is.null(X)) {
            add_note("Procrustes 跳过：响应+协变量的 PCA 拟合失败。")
          } else {
            k <- min(nax, ncol(X), ncol(Y))
            pr <- tryCatch(md_procrustes_res(CFG, X[keep, seq_len(k), drop = FALSE],
                                             Y[keep, seq_len(k), drop = FALSE],
                                             permutations = as.integer(pro$permutations %||% 999)),
                           error = function(e) NULL)
            if (is.null(pr)) {
              add_note("Procrustes 未算出（vegan::procrustes / protest 失败）。")
            } else {
              proc_res <- pr$residuals
              proc_res$sample_id <- as.character(proc_res$sample_id)
              for (f in utils::head(CFG$design$factors %||% character(0), 2)) {
                if (f %in% names(d))
                  proc_res[[f]] <- d[[f]][match(proc_res$sample_id, ids_all)]
              }
              proc_sum <- pr$summary
              proc_sum$metric <- pm
              proc_sum$axes <- k
              add_note("Procrustes：m2 = ", md_fmt_num(pr$summary$m2, 4),
                       "，PROTEST r = ", md_fmt_num(pr$summary$correlation, 3),
                       "，", md_fmt_p_fig(pr$summary$p), "（", pr$summary$permutations,
                       " 次置换）。残差只反映一致性，不代表因果方向。")
            }
          }
        }
      }
    }
  }

  # ===========================================================================
  # 2c. rf ---- 群落组成（PCoA 轴 + 门水平相对丰度）-> 响应变量的随机森林
  #     参考实现
  #     预测变量 = PCoA1..k + 门（或 CFG$response_models$rf$rank）相对丰度；
  #     响应 = CFG$response$var。⚠ 只在 linkage 含 "rf" 时执行。
  # ===========================================================================
  rf_tab <- NULL; rf_fit <- NULL
  if ("rf" %in% linkage) {
    rc <- CFG$response_models$rf %||% list()
    nax <- as.integer(rc$pcoa_axes %||% 5L)
    rr  <- as.character(rc$rank %||% "Phylum")
    pm  <- CFG$response_models$pair_metric %||% "bray"
    dmat <- tryCatch(md_dist_get(CFG, pm), error = function(e) NULL)
    feat <- NULL
    if (!is.null(dmat)) {
      M <- as.matrix(dmat)
      ids2 <- intersect(ids_all, rownames(M))
      if (length(ids2) >= 5) {
        pc <- tryCatch(ape::pcoa(stats::as.dist(M[ids2, ids2, drop = FALSE]),
                                 correction = "cailliez"), error = function(e) NULL)
        if (!is.null(pc)) {
          k <- min(nax, ncol(pc$vectors))
          feat <- as.data.frame(pc$vectors[, seq_len(k), drop = FALSE])
          names(feat) <- paste0("PCoA", seq_len(k))
          feat[[smp]] <- ids2
        }
      }
    }
    m <- tryCatch(md_core_asv(CFG, "alpha"), error = function(e) NULL)
    tax <- tryCatch(md_core(CFG)$tax, error = function(e) NULL)
    if (!is.null(m) && !is.null(tax) && rr %in% names(tax)) {
      # ⚠ md_glom() 返回 list(mat=, group=)，取 $mat 才是"类群 × 样本"矩阵
      gl <- tryCatch(md_glom(m, tax, rr, CFG$data$id_col %||% "asv_id"),
                     error = function(e) NULL)
      rel <- if (is.null(gl)) NULL else md_rel(gl$mat) * 100
    }
    if (!is.null(rel)) {
      rk <- as.data.frame(t(rel))
      rk[[smp]] <- rownames(rk)
      nm <- names(rk); nm[nm != smp] <- paste0(rr, "_", gsub("[^A-Za-z0-9_.]", "_", nm[nm != smp]))
      names(rk) <- nm
      feat <- if (is.null(feat)) rk else merge(feat, rk, by = smp, all = FALSE)
    }
    if (is.null(feat) || ncol(feat) < 3) {
      add_note("RF 跳过：构建不出特征矩阵（PCoA 轴或 ", rr, " 水平丰度不可用）。")
    } else {
      dd <- merge(feat, d[, c(smp, yv), drop = FALSE], by = smp, all = FALSE)
      xs <- setdiff(names(dd), c(smp, yv))
      fit_r <- tryCatch(md_rf_fit(CFG, dd, yv, xs,
                                  engine = as.character(rc$engine %||% "randomForest"),
                                  num_trees = as.integer(rc$num_trees %||% 3000L),
                                  importance = as.character(rc$importance %||% "permutation"),
                                  seed_offset = 8150L),
                        error = function(e) NULL)
      if (is.null(fit_r)) {
        add_note("RF 未算出（ranger / randomForest 不可用或样本不足）。")
      } else {
        rf_fit <- fit_r
        rf_tab <- fit_r$importance
        rf_tab$response <- yv; rf_tab$engine <- fit_r$engine
        rf_tab$R2 <- fit_r$r2; rf_tab$n <- fit_r$n
        add_note("RF（群落组成 -> ", yv, "）：engine=", fit_r$engine, "，n=", fit_r$n,
                 "，R2=", md_fmt_num(fit_r$r2, 3), "，特征数=", length(fit_r$predictors), "。")
      }
    }
  }

  # ===========================================================================
  # 4. 拼图 + 出表
  # ===========================================================================
  # ---- Procrustes 面板：残差按第一个设计因子的箱线（原图分组箱线口径）----
  if (!is.null(proc_res) && nrow(proc_res) > 0) {
    gf <- (CFG$design$factors %||% character(0))
    gf <- gf[gf %in% names(proc_res)]
    if (length(gf)) {
      g <- gf[1]
      pr <- proc_res
      pr$.grp <- as.character(pr[[g]])
      pr <- pr[!is.na(pr$.grp) & is.finite(pr$Residual), , drop = FALSE]
      if (nrow(pr) >= 4 && length(unique(pr$.grp)) >= 2) {
        kt <- tryCatch(stats::kruskal.test(Residual ~ .grp, data = pr),
                       error = function(e) NULL)
        labs_df <- tryCatch({
          gt <- md_group_test(CFG, pr, "Residual", ".grp", test = "kruskal_dunn",
                              p_adjust = CFG$model$p_adjust %||% "BH", min_n = 3,
                              order_levels = CFG$design$ordered_levels[[g]])
          gt$letters
        }, error = function(e) NULL)
        panels$C <- md_box_group(pr, y = "Residual", x = ".grp", letters = labs_df,
                                 colors = md_pal(CFG, unique(pr$.grp), g),
                                 xlab = g, ylab = "Procrustes residual",
                                 tag = "C", pt = pt,
                                 kw_text = if (is.null(kt)) NULL else md_kw_text(
                                   stats::setNames(list(kt$p.value), g)))
      }
    }
  }

  # ---- RF 面板：变量重要性（原图 RF 口径）----
  if (!is.null(rf_fit) && !is.null(rf_fit$importance) && nrow(rf_fit$importance) > 0) {
    panels$D <- md_rf_importance_plot(
      rf_fit$importance,
      top_n = as.integer((CFG$response_models$rf %||% list())$top_n %||% 15L),
      fill = "#3C5488", r2 = rf_fit$r2, label_r2 = TRUE,
      xlab = "Relative permutation importance", pt = pt, tag = "D")
  }

  if (length(panels) == 0) {
    md_skip(CFG, module, "三类关联都没有产出可用面板（详见日志）")
    return(md_result(module, "SKIP", "无可出图内容", notes = notes))
  }
  # 重排 tag，使面板字母连续
  panels <- panels[order(match(names(panels), LETTERS))]
  for (i in seq_along(panels)) {
    panels[[i]] <- panels[[i]] + ggplot2::labs(tag = LETTERS[i])
  }
  design <- switch(as.character(length(panels)),
                   "1" = "A",
                   "2" = "AB",
                   "3" = "AB\nCC",
                   "4" = "AB\nCD",
                   "AB\nCD")
  n_row <- if (length(panels) >= 3) 2 else 1
  fig <- md_compose(CFG, panels, design = design, guides = "collect", legend = "bottom")
  md_save_figure(CFG, fig, "71_response_linkage", height = if (n_row >= 2) 115 else 95)

  note_tab <- data.frame(
    Item = c("Response variable", "Response label", "Response transform",
             "Linkage modalities", "alpha focus", "ordination axes",
             "pairwise metric", "bootstrap B", "Distance source",
             "Scope", "Notes"),
    Value = c(yv, ylab, tr, paste(linkage, collapse = ", "),
              paste(CFG$response_models$alpha_focus %||% character(0), collapse = ", "),
              paste(CFG$response_models$ordination_axes %||% character(0), collapse = ", "),
              CFG$response_models$pair_metric %||% "bray",
              CFG$response_models$bootstrap_B %||% 1000,
              "13_distances cache if present, else local vegdist + ape::pcoa",
              paste0("figure = A alpha (single + partial + LMM/GAMM) + B ordination ",
                     "(conditional effect); pairwise = MRM/Mantel + cluster bootstrap, ",
                     "tables only (its figure belongs to 43_mantel_mrm, Fig.S12B)"),
              paste(notes, collapse = " ")),
    stringsAsFactors = FALSE)

  tabs <- list(
    Alpha_marginal = alpha_tab,
    Alpha_partial = partial_tab,
    Ordination_axes = pcoa_axis_tab,
    Pairwise_MRM = mrm_tab,
    Pairwise_Mantel = mant_tab,
    Pairwise_effect = boot_tab,
    Pairwise_saturation = sat_tab,
    Pairwise_pairs = if (isTRUE(CFG$table$detail)) pair_tab else NULL,
    Procrustes = proc_res,
    Procrustes_summary = proc_sum,
    RF_importance = rf_tab,
    Notes = note_tab
  )
  tp <- md_write_xlsx(CFG, "71_response_linkage", tabs)

  md_result(module, "OK",
            tables = basename(tp), figures = "figures/71_response_linkage.pdf",
            notes = notes)
}
