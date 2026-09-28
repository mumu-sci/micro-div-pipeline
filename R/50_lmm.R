# =============================================================================
# R/50_lmm.R ---- LMM / GLMM / GAMM 与**条件效应**（合并）
# -----------------------------------------------------------------------------
# 一次做完四件事：
#   1. 拟合 lmerTest::lmer(REML = CFG$model$reml) -> Satterthwaite P；
#      CFG$model$family 非 gaussian 时改用 lme4::glmer（GLMM，Wald z）；
#      CFG$model$random = NULL 时退化为 lm/glm。
#      **优雅降级**：随机效应方差为 0 / singular fit / 设计矩阵不满秩 / 组内样本
#      太少 / 全 NA ---- 一律不崩，退化后把原因写进 Fit$note 与 Guard sheet。
#   2. mgcv::gamm（失败退 gam）按 edf > CFG$model$edf_nonlin 判非线性：
#      越阈值 -> 画平滑曲线，否则退回直线；图面标注随实际曲线来源联动
#      （nifH 踩过的坑：5D 画的是直线、标注却印 GAMM 平滑项 P）。
#   3. performance::r2() 报 R2m / R2c（lm/glm 只报 R2）；
#      CFG$model$robust_se = TRUE 时额外报**站点级整簇**稳健 SE（CR1 三明治，
#      df = 簇数 - 1），与 Satterthwaite 口径并列而不是替换。
#   4. 条件效应曲线：CFG$model$focus x CFG$model$covariates，全变量 z 化 ->
#      CFG$model$grid_n 网格 -> 其余协变量置 z = 0 -> 预测 +/- 1.96SE -> 反变换回原始单位。
#
# **站点级共线守卫**（nifH Fig. 5B / FigS9 的教训）：
#   焦点变量若在站点内唯一值 = 1（md_within_unique()），它就是站点级变量，与
#   (1|Site) 随机截距共线 -> 自动剔除同为站点级的协变量、禁用平滑（此时 s(x) 的
#   有效样本量 = 站点数，edf 不可识别），并在 Guard sheet 与图面注记里说明。
#
# 依赖：40（cache/40_env.rds）。Soil_PC1 只经 md_derived() 取，**不重算**；
#       40 缺失时自动剔除该变量并记录，不报错。
# =============================================================================

md_lmm <- function(CFG) {
  module <- "50"
  if (!md_need(CFG, "model", module))
    return(md_result(module, "SKIP", "缺少 lme4/lmerTest/performance"))
  if (!md_need(CFG, "gam", module))
    return(md_result(module, "SKIP", "缺少 mgcv"))
  if (!md_precheck(CFG, module, need_response = TRUE))
    return(md_result(module, "SKIP", "前置条件不满足"))
  md_load(c("lme4", "lmerTest", "performance"))

  # ---------------------------------------------------------------------------
  # 0. 参数
  # ---------------------------------------------------------------------------
  smp       <- CFG$data$sample_col
  y_var     <- CFG$model$response %||% CFG$response$var
  focus     <- as.character(CFG$model$focus %||% character(0))
  focus     <- unique(focus[nzchar(focus)])
  covs_spec <- as.character(CFG$model$covariates %||% character(0))
  covs_spec <- unique(covs_spec[nzchar(covs_spec)])
  site_var  <- CFG$design$site_var
  random    <- CFG$model$random
  if (!length(focus)) {
    md_skip(CFG, module, "CFG$model$focus 未设置（不许自动探索）")
    return(md_result(module, "SKIP", "CFG$model$focus 未设置"))
  }
  reml     <- isTRUE(CFG$model$reml)
  std      <- isTRUE(CFG$model$standardize %||% TRUE)
  edf_thr  <- as.numeric(CFG$model$edf_nonlin %||% 1.1)
  grid_n   <- as.integer(CFG$model$grid_n %||% 200)
  robust   <- isTRUE(CFG$model$robust_se)
  fam_name <- tolower(CFG$model$family %||% "gaussian")
  fam <- switch(fam_name,
                gaussian = stats::gaussian(), binomial = stats::binomial(),
                poisson  = stats::poisson(),  gamma = stats::Gamma(),
                stats::gaussian())

  d <- md_analysis_with_derived(CFG)
  if (is.null(d) || !nrow(d)) {
    md_skip(CFG, module, "analysis 表为空（先跑 05 / 40）")
    return(md_result(module, "SKIP", "analysis 表为空"))
  }
  ids      <- as.character(d[[smp]])
  has_site <- !is.null(site_var) && site_var %in% names(d)
  use_re   <- has_site && !is.null(random) && nzchar(as.character(random))

  # ---------------------------------------------------------------------------
  # 1. 变量解析：analysis 表 -> md_derived(40) -> 本地由 05/13/21 派生
  # ---------------------------------------------------------------------------
  cv <- new.env()
  local_z <- function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
  }
  local_alpha <- function() {
    if (!is.null(cv$alpha)) return(cv$alpha)
    m <- md_core_asv(CFG, "alpha")          # ASV x 样本
    rel <- md_rel(m)
    sh  <- -colSums(rel * log(ifelse(rel > 0, rel, 1)))
    ob  <- colSums(m > 0)
    cv$alpha <- data.frame(sample_id = colnames(m), Shannon = as.numeric(sh),
                           Observed = as.numeric(ob),
                           Simpson = as.numeric(colSums(rel^2)),
                           Pielou = as.numeric(sh / log(pmax(ob, 2))),
                           stringsAsFactors = FALSE)
    cv$alpha
  }
  # 备用面板与主分支共用解析器；用调用者环境接收（避免误写 cv$res 之外的键）
  md_lmm_env <- environment()
  local_phylo <- function() {
    if (!is.null(cv$phylo)) return(cv$phylo)
    cv$phylo <- tryCatch(md_phylo_cache(CFG), error = function(e) NULL)
    cv$phylo
  }
  local_pcoa1 <- function() {
    if (!is.null(cv$pcoa1)) return(cv$pcoa1)
    dd <- md_dist_get(CFG, CFG$response_models$pair_metric %||% CFG$distances$metrics[1])
    if (is.null(dd)) { cv$pcoa1 <- NULL; cv$pcoa1_ve <- NA_real_; return(NULL) }
    M <- as.matrix(dd); v <- NULL; ve <- NA_real_
    if (requireNamespace("ape", quietly = TRUE)) {
      p <- tryCatch(ape::pcoa(M, correction = "cailliez"), error = function(e) NULL)
      if (!is.null(p)) {
        eig <- p$values$Corr_eig
        if (is.null(eig) || all(is.na(eig))) eig <- p$values$Eigenvalues
        v <- stats::setNames(p$vectors[, 1], rownames(p$vectors))
        if (!all(is.na(eig))) ve <- eig[1] / sum(eig[eig > 0])
      }
    }
    if (is.null(v)) {
      cm <- tryCatch(stats::cmdscale(M, k = 1), error = function(e) NULL)
      if (!is.null(cm)) v <- stats::setNames(cm[, 1], rownames(M))
    }
    cv$pcoa1_ve <- ve; cv$pcoa1 <- v; v
  }

  prov <- list(); prov_seen <- character(0)
  add_prov <- function(v, source, n_used) {
    if (v %in% prov_seen) return(invisible(NULL))
    prov_seen <<- c(prov_seen, v)
    prov[[length(prov) + 1L]] <<- data.frame(
      variable = v, source = source, n_used = n_used, stringsAsFactors = FALSE)
    invisible(NULL)
  }
  local_resolve <- function(v) {
    if (!is.null(v) && v %in% names(d)) {
      x <- suppressWarnings(as.numeric(d[[v]]))
      if (!all(is.na(x)))
        return(list(value = x, source = "analysis 表（meta / env / 40 派生变量）"))
    }
    dv <- tryCatch(md_derived(CFG, v, ids), error = function(e) NULL)
    if (!is.null(dv) && any(!is.na(dv)))
      return(list(value = suppressWarnings(as.numeric(dv)),
                  source = "md_derived()（40 号派生变量）"))
    a <- local_alpha()
    if (v %in% names(a))
      return(list(value = a[[v]][match(ids, a$sample_id)],
                  source = "本地由 05 核心矩阵（alpha 口径）计算"))
    pc <- local_phylo()
    if (!is.null(pc) && v %in% names(pc))
      return(list(value = suppressWarnings(as.numeric(pc[[v]][match(ids, pc[[1]])])),
                  source = "md_phylo_cache()（21 号）"))
    po <- local_pcoa1()
    if (!is.null(po) && v %in% c("PCoA1", "PCoA_1", "PC1"))
      return(list(value = as.numeric(po[match(ids, names(po))]),
                  source = sprintf("本地由 13 距离矩阵 PCoA 轴 1 计算（解释率 %s）",
                                   if (is.na(cv$pcoa1_ve %||% NA_real_)) "NA"
                                   else sprintf("%.1f%%", 100 * cv$pcoa1_ve))))
    NULL
  }

  res <- list()
  ry <- local_resolve(y_var)
  if (is.null(ry)) {
    md_skip(CFG, module, paste0("响应变量 ", y_var, " 无法解析"))
    return(md_result(module, "SKIP", paste0("响应变量 ", y_var, " 不可用")))
  }
  res[[y_var]] <- ry; add_prov(y_var, ry$source, sum(!is.na(ry$value)))
  # 备用派生环境：把主分支已解析的变量（含响应变量）暴露给 §3b 的备用面板
  cv$res <- res
  dropped_var <- character(0)
  for (v in unique(c(focus, covs_spec))) {
    r <- local_resolve(v)
    if (is.null(r)) { dropped_var <- c(dropped_var, v); next }
    res[[v]] <- r; add_prov(v, r$source, sum(!is.na(r$value)))
  }
  focus_ok <- intersect(focus, names(res))
  if (length(dropped_var))
    md_warn(CFG, "50: 以下变量无法解析，已剔除：", paste(dropped_var, collapse = ", "))
  if (!length(focus_ok)) {
    md_skip(CFG, module, paste0("焦点变量全部不可用（", paste(focus, collapse = ", "),
                                "）；Soil_PC1 需先跑 40"))
    return(md_result(module, "SKIP", "焦点变量全部不可用"))
  }

  # ---------------------------------------------------------------------------
  # 2. 逐焦点拟合
  # ---------------------------------------------------------------------------
  local_cr_se <- function(fit, cluster) {
    X <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
    if (is.null(X) || !nrow(X)) return(NULL)
    b <- tryCatch(lme4::fixef(fit), error = function(e) NULL)
    if (is.null(b)) b <- tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    keep <- intersect(names(b), colnames(X))
    keep <- keep[!is.na(b[keep])]
    if (length(keep) < 2) return(NULL)
    X <- X[, keep, drop = FALSE]
    u <- tryCatch(stats::residuals(fit), error = function(e) NULL)
    if (is.null(u) || length(u) != nrow(X)) return(NULL)
    bread <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
    if (is.null(bread)) return(NULL)
    g  <- as.character(cluster)
    meat <- matrix(0, ncol(X), ncol(X), dimnames = list(keep, keep))
    for (k in unique(g[!is.na(g)])) {
      i <- which(g == k)
      s <- colSums(X[i, , drop = FALSE] * u[i])
      meat <- meat + outer(s, s)
    }
    V  <- bread %*% meat %*% bread
    se <- sqrt(pmax(diag(V), 0))
    G <- length(unique(g[!is.na(g)])); n <- nrow(X); k <- ncol(X)
    se * sqrt((G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1)))
  }

  coef_rows <- list(); fit_rows <- list(); nl_rows <- list(); guard_rows <- list()
  curve_rows <- list()
  # 图件构造需要的中间量（只影响画法，不参与任何统计）
  curve_df   <- list()   # focus -> curve data.frame（画曲线用，与 curve_rows 同源）
  pt_df      <- list()   # focus -> 原始单位散点（含 .row 以便挂分组因子）
  ux_list    <- list()   # focus -> 焦点变量在 z 空间的取值（决定 x 轴刻度）

  # 站内唯一值：派生变量不在 analysis 表里，必须先拼成 data.frame 再算
  local_wu <- function(v) {
    if (!has_site) return(NA_integer_)
    md_within_unique(data.frame(.v = res[[v]]$value,
                                .s = as.character(d[[site_var]])), ".v", ".s")
  }

  for (f in focus_ok) {
    cv_use <- setdiff(intersect(covs_spec, names(res)), c(f, y_var))
    wu_focus <- local_wu(f)
    site_level <- isTRUE(!is.na(wu_focus) && wu_focus == 1)
    dropped_cv <- character(0); guard_note <- character(0)

    if (site_level) {
      wu_cv <- vapply(cv_use, local_wu, integer(1))
      drop_now <- cv_use[!is.na(wu_cv) & wu_cv == 1]
      cv_use <- setdiff(cv_use, drop_now); dropped_cv <- c(dropped_cv, drop_now)
      guard_note <- c(guard_note,
                      paste0("焦点是站点级变量（站内唯一值 = 1）-> 剔除同为站点级的协变量：",
                             if (length(drop_now)) paste(drop_now, collapse = ", ") else "无"),
                      "禁用平滑：站点级焦点与 (1|Site) 共线，edf 不可识别")
    }
    const_cv <- cv_use[vapply(cv_use, function(v) {
      x <- res[[v]]$value; length(unique(x[!is.na(x)])) < 2 }, logical(1))]
    if (length(const_cv)) {
      cv_use <- setdiff(cv_use, const_cv); dropped_cv <- c(dropped_cv, const_cv)
      guard_note <- c(guard_note, paste0("常量协变量剔除：", paste(const_cv, collapse = ", ")))
    }
    if (length(cv_use) >= 2) {
      vd <- as.data.frame(stats::setNames(lapply(cv_use, function(v) res[[v]]$value), cv_use))
      vif <- md_vif(vd, cv_use); thr <- CFG$env$vif_thr %||% 10
      if (any(!is.na(vif$vif) & vif$vif > thr)) {
        bad <- vif$variable[which.max(vif$vif)]
        cv_use <- setdiff(cv_use, bad); dropped_cv <- c(dropped_cv, bad)
        guard_note <- c(guard_note, paste0("协变量 VIF > ", thr, "，剔除：", bad))
      }
    }

    # ---- 建模数据 ----
    dd_raw <- data.frame(.row = seq_len(nrow(d)))
    for (v in unique(c(y_var, f, cv_use))) dd_raw[[v]] <- res[[v]]$value
    if (use_re) dd_raw[[site_var]] <- as.character(d[[site_var]])
    dd_raw <- dd_raw[stats::complete.cases(dd_raw[, setdiff(names(dd_raw), ".row"), drop = FALSE]), ,
                     drop = FALSE]
    if (use_re) dd_raw[[site_var]] <- factor(dd_raw[[site_var]])
    if (nrow(dd_raw) < 6) {
      guard_rows[[length(guard_rows) + 1L]] <- data.frame(
        focus = f, within_unique_max = wu_focus, site_level = site_level,
        dropped_covariates = if (length(cv_use)) paste(cv_use, collapse = ", ") else "无",
        smooth_disabled = TRUE, note = "完整样本 < 6，跳过该焦点", stringsAsFactors = FALSE)
      next
    }
    dd <- dd_raw; cs <- list()
    for (v in unique(c(y_var, f, cv_use))) {
      cen <- mean(dd_raw[[v]], na.rm = TRUE); scl <- stats::sd(dd_raw[[v]], na.rm = TRUE)
      cs[[v]] <- c(center = cen, scale = scl)
      dd[[v]] <- if (std) local_z(dd_raw[[v]]) else dd_raw[[v]]
    }

    # ---- 拟合 ----
    rhs <- c(f, cv_use, if (use_re) random else NULL)
    fml <- stats::as.formula(paste(y_var, "~", paste(rhs, collapse = " + ")))
    note_fit <- character(0); fit <- NULL; est <- NA_character_
    muffle <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
    if (fam_name == "gaussian") {
      if (use_re)
        fit <- tryCatch(muffle(lmerTest::lmer(fml, data = dd, REML = reml)),
                        error = function(e) NULL)
      est <- if (use_re) "LMM (lmerTest)" else "LM"
      if (is.null(fit)) {
        fit <- tryCatch(stats::lm(fml, data = dd), error = function(e) NULL)
        if (use_re) note_fit <- c(note_fit, "lmer 拟合失败 -> 退化为 lm（随机效应被忽略）")
        est <- "LM (LMM 失败后降级)"
      }
    } else {
      if (use_re)
        fit <- tryCatch(muffle(lme4::glmer(fml, data = dd, family = fam)),
                        error = function(e) NULL)
      est <- if (use_re) paste0("GLMM (glmer, ", fam_name, ")")
             else paste0("GLM (", fam_name, ")")
      if (is.null(fit)) {
        fit <- tryCatch(stats::glm(fml, data = dd, family = fam), error = function(e) NULL)
        note_fit <- c(note_fit, "glmer 拟合失败 -> 退化为 glm")
        est <- paste0("GLM (", fam_name, "，GLMM 失败后降级)")
      }
    }
    if (is.null(fit)) {
      guard_rows[[length(guard_rows) + 1L]] <- data.frame(
        focus = f, within_unique_max = wu_focus, site_level = site_level,
        dropped_covariates = if (length(cv_use)) paste(cv_use, collapse = ", ") else "无",
        smooth_disabled = TRUE,
        note = "所有候选模型均拟合失败（设计矩阵退化 / 不满秩）", stringsAsFactors = FALSE)
      next
    }
    is_mer <- inherits(fit, "merMod")
    singular <- FALSE
    if (is_mer) {
      singular <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) FALSE))
      v0 <- tryCatch({
        vcd <- as.data.frame(lme4::VarCorr(fit))
        if (nrow(vcd) && "vcov" %in% names(vcd)) as.numeric(vcd$vcov[1]) else NA_real_
      }, error = function(e) NA_real_)
      if (singular || (is.finite(v0) && v0 < 1e-6))
        note_fit <- c(note_fit, sprintf("singular fit：随机效应方差 ~ %s（与 0 无法区分）",
                                        format(signif(as.numeric(v0), 3))))
    }

    # ---- 系数 ----
    co <- as.data.frame(summary(fit)$coefficients)
    cn <- names(co)
    pcol <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "p.value", "P"), cn)[1]
    tcol <- intersect(c("t value", "z value"), cn)[1]
    if (is.na(tcol)) tcol <- NULL
    if (is.na(pcol)) {
      co$.P <- 2 * stats::pt(-abs(co[[tcol]]), df = if ("df" %in% cn) co[["df"]] else Inf)
      pcol <- ".P"
    }
    dfv <- if ("df" %in% cn) as.numeric(co[["df"]])
           else if (inherits(fit, "lm")) rep(stats::df.residual(fit), nrow(co)) else Inf
    nco <- nrow(co)
    cdf <- data.frame(
      focus = rep(f, nco), term = rownames(co), is_focus = rownames(co) == f,
      estimate = as.numeric(co[["Estimate"]]),
      se = as.numeric(co[["Std. Error"]]),
      df = dfv,
      statistic = if (!is.null(tcol)) as.numeric(co[[tcol]]) else rep(NA_real_, nco),
      p = as.numeric(co[[pcol]]),
      estimator = rep(est, nco), n = rep(stats::nobs(fit), nco), stringsAsFactors = FALSE)
    crit <- ifelse(is.finite(cdf$df), stats::qt(0.975, pmax(cdf$df, 1)), 1.96)
    cdf$lo <- cdf$estimate - crit * cdf$se
    cdf$hi <- cdf$estimate + crit * cdf$se
    cdf$se_robust <- NA_real_; cdf$p_robust <- NA_real_
    if (robust && has_site) {
      se_r <- local_cr_se(fit, dd_raw[[site_var]])
      if (!is.null(se_r)) {
        m <- match(names(se_r), cdf$term)
        cdf$se_robust[m] <- as.numeric(se_r)
        G <- length(unique(dd_raw[[site_var]]))
        cdf$p_robust <- 2 * stats::pt(-abs(cdf$estimate / cdf$se_robust), df = max(G - 1, 1))
        note_fit <- c(note_fit, sprintf("稳健 SE：站点级整簇 CR1 三明治，df = %d", G - 1))
      } else note_fit <- c(note_fit, "稳健 SE 计算失败（设计矩阵不可逆）")
    }
    cdf$singular <- singular
    cdf$note <- paste(unique(note_fit), collapse = "；")
    coef_rows[[length(coef_rows) + 1L]] <- cdf

    # ---- 拟合优度 ----
    r2m <- NA_real_; r2c <- NA_real_; icc <- NA_real_; sig2 <- NA_real_
    rr <- tryCatch(performance::r2(fit), error = function(e) NULL)
    if (!is.null(rr)) {
      r2m <- as.numeric(rr$R2_marginal %||% NA_real_)
      r2c <- as.numeric(rr$R2_conditional %||% rr$R2 %||% NA_real_)
    }
    if (inherits(fit, "lm")) sig2 <- as.numeric(summary(fit)$sigma %||% NA_real_)
    if (is_mer) icc <- tryCatch(as.numeric(performance::icc(fit)$ICC_adjusted),
                                error = function(e) NA_real_)
    # 随机效应方差 ~ 0 时 performance::r2 不报 R2c：照实取边际 R2 并记明，
    # 不留下无法解释的 NA
    if (is_mer && is.na(r2c) && is.finite(r2m)) {
      r2c <- r2m
      note_fit <- c(note_fit, "随机效应方差 ~ 0，无法分解 -> R2c 取边际 R2")
    }
    fit_rows[[length(fit_rows) + 1L]] <- data.frame(
      focus = f, response = y_var, estimator = est,
      formula = paste(deparse(stats::formula(fit)), collapse = " "),
      n = stats::nobs(fit),
      n_site = if (has_site) length(unique(dd_raw[[site_var]])) else NA_integer_,
      R2_marginal = r2m, R2_conditional = r2c, ICC = icc,
      AIC = tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_),
      sigma = sig2, family = fam_name, singular = singular, standardize = std,
      note = paste(unique(note_fit), collapse = "；"), stringsAsFactors = FALSE)

    # ---- GAMM 非线性判定 ----
    allow_smooth <- !site_level && nrow(dd_raw) >= 8 &&
      length(unique(dd_raw[[f]])) >= 4
    edf <- NA_real_; p_smooth <- NA_real_; curve_src <- "linear"
    g <- NULL
    if (allow_smooth) {
      gam_rhs <- paste0("s(", f, ")",
                        if (length(cv_use)) paste0(" + ", paste(cv_use, collapse = " + ")) else "")
      gfml <- stats::as.formula(paste(y_var, "~", gam_rhs))
      if (use_re) {
        dd2 <- dd; dd2$md_site <- dd[[site_var]]
        g <- tryCatch(muffle(mgcv::gamm(gfml, random = list(md_site = ~1),
                                        data = dd2, method = "REML")),
                      error = function(e) NULL)
      }
      if (is.null(g))
        g <- tryCatch(muffle(mgcv::gam(gfml, data = dd, method = "REML")),
                      error = function(e) NULL)
      if (!is.null(g)) {
        st <- tryCatch(summary(g$gam)$s.table, error = function(e) NULL)
        key <- paste0("s(", f, ")")
        if (!is.null(st) && key %in% rownames(st)) {
          edf <- as.numeric(st[key, "edf"]); p_smooth <- as.numeric(st[key, "p-value"])
        }
        if (is.finite(edf) && edf > edf_thr) curve_src <- "GAMM"
      }
      if (is.null(g)) note_fit <- c(note_fit, "GAMM 拟合失败，非线性未检验")
    } else if (site_level) {
      note_fit <- c(note_fit, "站点级焦点：禁用平滑")
    } else {
      note_fit <- c(note_fit, "唯一值/样本量不足，禁用平滑")
    }
    guard_rows[[length(guard_rows) + 1L]] <- data.frame(
      focus = f, within_unique_max = wu_focus, site_level = site_level,
      dropped_covariates = if (length(dropped_cv)) paste(dropped_cv, collapse = ", ") else "无",
      covariates_used = if (length(cv_use)) paste(cv_use, collapse = ", ") else "无",
      smooth_disabled = !allow_smooth,
      note = if (length(guard_note)) paste(guard_note, collapse = "；") else "未触发守卫",
      stringsAsFactors = FALSE)

    # ---- 条件效应曲线（z 空间预测 -> 反变换回原始单位） ----
    back_y <- function(z) as.numeric(z) * cs[[y_var]]["scale"] + cs[[y_var]]["center"]
    ux <- sort(unique(dd[[f]]))
    xs <- if (length(ux) <= 10) ux else seq(min(dd[[f]]), max(dd[[f]]), length.out = grid_n)
    nd <- data.frame(x = xs); names(nd) <- f
    for (v in cv_use) nd[[v]] <- if (std) 0 else cs[[v]]["center"]
    nd <- nd[, c(f, cv_use), drop = FALSE]
    curve <- NULL
    if (identical(curve_src, "GAMM") && !is.null(g)) {
      pr <- tryCatch(stats::predict(g$gam, newdata = nd, se.fit = TRUE), error = function(e) NULL)
      if (!is.null(pr)) {
        fz <- as.numeric(pr$fit); sz <- as.numeric(pr$se.fit)
        if (fam_name == "gaussian")
          curve <- data.frame(x_z = xs, x_raw = xs * cs[[f]]["scale"] + cs[[f]]["center"],
                              fit = back_y(fz), lower = back_y(fz - 1.96 * sz),
                              upper = back_y(fz + 1.96 * sz))
        else
          curve <- data.frame(x_z = xs, x_raw = xs * cs[[f]]["scale"] + cs[[f]]["center"],
                              fit = fz, lower = fz - 1.96 * sz, upper = fz + 1.96 * sz)
      }
    }
    if (is.null(curve)) {
      curve_src <- if (site_level) "LMM (site-level focus)" else "LMM"
      X <- tryCatch(stats::model.matrix(
        stats::as.formula(paste0("~ ", paste(c(f, cv_use), collapse = " + "))), data = nd),
        error = function(e) NULL)
      if (!is.null(X)) {
        b <- if (is_mer) lme4::fixef(fit) else stats::coef(fit)
        V <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
        keep <- intersect(names(b), colnames(X)); keep <- keep[!is.na(b[keep])]
        if (!is.null(V) && length(keep) >= 2) {
          X <- X[, keep, drop = FALSE]; b <- b[keep]; V <- V[keep, keep, drop = FALSE]
          fz <- as.vector(X %*% b)
          sz <- sqrt(pmax(0, rowSums((X %*% V) * X)))
          if (fam_name == "gaussian")
            curve <- data.frame(x_z = xs, x_raw = xs * cs[[f]]["scale"] + cs[[f]]["center"],
                                fit = back_y(fz), lower = back_y(fz - 1.96 * sz),
                                upper = back_y(fz + 1.96 * sz))
          else
            curve <- data.frame(x_z = xs, x_raw = xs * cs[[f]]["scale"] + cs[[f]]["center"],
                                fit = fz, lower = fz - 1.96 * sz, upper = fz + 1.96 * sz)
        }
      }
    }
    p_focus <- cdf$p[cdf$term == f][1]
    # 非线性判定表在曲线**落地之后**才写：curve_drawn 一定等于实际画出来的曲线，
    # 不允许出现"画直线、标注 GAMM P"的口径错配（nifH Fig.5D 踩过）。
    nl_rows[[length(nl_rows) + 1L]] <- data.frame(
      focus = f, edf = edf, P_smooth = p_smooth, threshold = edf_thr,
      curve_drawn = if (identical(curve_src, "GAMM")) "GAMM smooth"
                    else if (is.null(curve)) "none" else "linear",
      smooth_disabled = !allow_smooth,
      note = if (!allow_smooth) "禁用平滑（站点级焦点或取值/样本量不足，edf 不可识别）"
             else if (is.na(edf) || is.null(g)) "edf 不可用（GAMM 拟合失败）-> 退回直线"
             else if (edf > edf_thr && !identical(curve_src, "GAMM"))
               "edf 越阈值但 GAMM 曲线预测失败 -> 退回直线"
             else if (identical(curve_src, "GAMM")) "edf 越阈值 -> 画 GAMM 平滑曲线"
             else "edf 未越阈值 -> 退回直线",
      stringsAsFactors = FALSE)
    if (!is.null(curve)) {
      curve$focus <- f
      curve_rows[[length(curve_rows) + 1L]] <- curve
      curve_df[[f]] <- curve
      # 原始单位散点（原脚本 Fig5 B-E：geom_point(x = 原始值, y = 原始响应)）
      pt_df[[f]] <- data.frame(.row = dd_raw$.row, x = dd_raw[[f]], y = dd_raw[[y_var]])
      names(pt_df[[f]])[2] <- f
      names(pt_df[[f]])[3] <- y_var
      ux_list[[f]] <- ux
    }
    # 把后补的 note 回写（GAMM/曲线阶段新增的 note）
    fit_rows[[length(fit_rows)]]$note <- paste(unique(note_fit), collapse = "；")
  }

  if (!length(coef_rows)) {
    md_skip(CFG, module, "所有焦点模型都无法拟合（设计矩阵退化 / 样本不足）")
    return(md_result(module, "SKIP", "所有焦点模型都无法拟合"))
  }
  coef_tab  <- do.call(rbind, coef_rows)
  coef_tab$P_BH <- md_p_adjust(coef_tab$p, CFG$model$p_adjust %||% "BH")
  coef_tab$sig  <- ifelse(is.na(coef_tab$p), "",
                          ifelse(coef_tab$p < 0.001, "***",
                                 ifelse(coef_tab$p < 0.01, "**",
                                        ifelse(coef_tab$p < 0.05, "*", "ns"))))
  fit_tab   <- do.call(rbind, fit_rows)
  nl_tab    <- do.call(rbind, nl_rows)
  guard_tab <- do.call(rbind, guard_rows)
  prov_tab  <- do.call(rbind, prov)

  # ---------------------------------------------------------------------------
  # 3. 图：系数森林图（A，对 04_Fig04_SEM.R §C/D make_coef_panel）
  #       + 条件效应曲线（B-E，对 05_Fig05_decoupling.R §B-E）
  # ---------------------------------------------------------------------------
  # 画法要点（逐条照抄原脚本）：
  #   A：误差棒（height 0.18）+ 实心白描边圆点 shape 21；颜色 = 方向×显著性
  #      （正 + md_pal_cont(diverging)[3] / 负 - [1] / 不显著 grey60）；参考线 0
  #      灰色虚线；右上角标注 R2m / R2c / n；面板 3:2（宽:高）；y 轴标签加粗；
  #      节点顺序按 |estimate| 降序（原脚本 order(-abs(Estimate))）。
  #   B-E：原始单位散点（按模型里另一个因子着色）+ 置信带 + 曲线（MD_REG 口径）；
  #        右上角标注效应 P 与曲线来源；x = 焦点变量名，y = NULL；一列面板。
  pt <- md_pt(CFG)
  # 面板字母：A = 系数森林图；B 起按 focus_ok 顺序给条件效应曲线，备用面板继续顺延
  tag_of <- function(i) if (i <= length(LETTERS)) LETTERS[i] else as.character(i)
  fr <- coef_tab[coef_tab$term != "(Intercept)", , drop = FALSE]
  if (!nrow(fr)) fr <- coef_tab
  focus_lab_of <- function(f) {
    hit <- coef_tab$term[coef_tab$focus == f & coef_tab$term != "(Intercept)"]
    if (length(hit)) hit[1] else f
  }
  fr$Label <- ifelse(fr$is_focus, vapply(fr$focus, focus_lab_of, character(1)),
                     as.character(fr$term))
  fr$Col <- ifelse(is.na(fr$p), "grey60",
                   ifelse(fr$p < 0.05,
                          ifelse(fr$estimate > 0, md_pal_cont(CFG, "diverging")[3],
                                 md_pal_cont(CFG, "diverging")[1]),
                          "grey60"))
  fr <- fr[order(-abs(fr$estimate)), , drop = FALSE]
  fr$Lbl <- factor(fr$Label, levels = rev(unique(fr$Label)))
  fit_main <- fit_tab[fit_tab$focus %in% fr$focus, , drop = FALSE]
  if (!nrow(fit_main)) fit_main <- fit_tab
  n_first <- suppressWarnings(as.numeric(fit_main$n[1]))
  ann_coef <- paste0("R2m = ", md_fmt_num(suppressWarnings(as.numeric(fit_main$R2_marginal[1]))),
                     "\nR2c = ", md_fmt_num(suppressWarnings(as.numeric(fit_main$R2_conditional[1]))),
                     "\nn = ", if (is.na(n_first)) "NA" else as.integer(n_first))

  p_forest <- ggplot2::ggplot(fr, ggplot2::aes(x = .data$estimate, y = .data$Lbl)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey35",
                        linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                           orientation = "y", height = 0.18, linewidth = MD_LW$line,
                           colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$Col), size = MD_PT_SIZE$scatter * 2.0,
                        shape = 21, colour = "white", stroke = 0.3) +
    ggplot2::scale_fill_identity() +
    ggplot2::annotate("text", x = Inf, y = Inf, label = ann_coef,
                      hjust = 1.05, vjust = 1.35, lineheight = 0.9,
                      size = pt2gg(pt$ann), colour = "grey30") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = if (std) "Standardized coefficient" else "Coefficient",
                  y = NULL, tag = "A") +
    theme_md(CFG = CFG) +
    ggplot2::theme(aspect.ratio = 2 / 3,
                   axis.text.y = ggplot2::element_text(face = "bold", size = pt$tick))

  # ---- 逐焦点条件效应曲线（B 起，按 focus_ok 顺序） ----
  p_list <- list()
  for (f in focus_ok) {
    cvd <- curve_df[[f]]
    if (is.null(cvd) || !nrow(cvd)) next
    nl <- nl_tab[nl_tab$focus == f, , drop = FALSE]
    src_lab <- if (nrow(nl)) as.character(nl$curve_drawn[1]) else "linear"
    est_lab <- if (identical(src_lab, "GAMM smooth")) "GAMM curve" else "linear"
    p_focus_p <- coef_tab$p[coef_tab$focus == f & coef_tab$is_focus][1]
    lbl <- sprintf("%s: P = %s\n(LMM, %s)", f, md_fmt_p(p_focus_p), est_lab)

    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = cvd, ggplot2::aes(x = .data$x_raw, ymin = .data$lower,
                                                    ymax = .data$upper),
                           fill = MD_REG$col, alpha = MD_REG$ribbon_alpha, colour = NA) +
      ggplot2::geom_line(data = cvd, ggplot2::aes(x = .data$x_raw, y = .data$fit),
                         colour = MD_REG$col, linewidth = MD_REG$lw)

    # 散点：x = 焦点变量原始值，y = 响应变量原始值；按模型里另一个因子着色
    dpt <- pt_df[[f]]
    if (!is.null(dpt) && nrow(dpt)) {
      # 着色因子**优先用设计因子**（Vegetation / Zone）——图例显示的是**分组名**；
      # 原稿 05_Fig05_decoupling.R §B–E 就是按 Vegetation 着色。
      # ⚠ res 里只放了 y / focus / covariates，设计因子要**按需再解析一次**（get_val），
      #   否则候选会被下面的过滤器全部滤掉，导致全图根本没有图例。
      # ⚠ res/local_resolve 只走 as.numeric()，**分类的设计因子会被解析成全 NA**，
      #   所以这里直接从模型数据 d 里取原值（分类列原样返回字符）。
      get_val <- function(v) {
        if (!is.null(v) && v %in% names(d)) {
          raw <- d[[v]]
          if (!is.numeric(raw)) return(as.character(raw))
          x <- suppressWarnings(as.numeric(raw))
          if (!all(is.na(x))) return(x)
        }
        x <- res[[v]]$value
        if (is.null(x)) x <- tryCatch(local_resolve(v)$value, error = function(e) NULL)
        x
      }
      # rev() 是为了让第二个设计因子（Vegetation）优先于第一个（Zone）
      gvars <- setdiff(unique(c(rev(CFG$design$factors %||% character(0)), covs_spec, focus_ok)),
                       c(f, y_var))
      gvars <- gvars[vapply(gvars, function(v) {
        x <- get_val(v)
        !is.null(x) && length(unique(x[!is.na(x)])) >= 2
      }, logical(1))]
      if (length(gvars)) {
        gv <- gvars[1]
        dpt$.g <- factor(as.character(get_val(gv))[match(dpt$.row, seq_len(nrow(d)))])
        dpt <- dpt[!is.na(dpt$.g), , drop = FALSE]
        p <- p + ggplot2::geom_point(data = dpt,
                                     ggplot2::aes(x = .data[[f]], y = .data[[y_var]],
                                                  colour = .data$.g),
                                     alpha = MD_PT_SIZE$point_alpha,
                                     size = MD_PT_SIZE$point) +
          ggplot2::scale_colour_manual(values = md_pal(CFG, levels(dpt$.g), gv),
                                       name = gv)
      } else {
        p <- p + ggplot2::geom_point(data = dpt,
                                     ggplot2::aes(x = .data[[f]], y = .data[[y_var]]),
                                     colour = "grey30", alpha = MD_PT_SIZE$point_alpha,
                                     size = MD_PT_SIZE$point)
      }
    }

    p <- p +
      ggplot2::annotate("text", x = Inf, y = Inf, label = lbl, hjust = 1, vjust = 1.3,
                        lineheight = 0.9, size = pt2gg(pt$ann), colour = "grey30") +
      ggplot2::labs(x = f, y = NULL, tag = tag_of(1L + match(f, focus_ok))) +
      theme_md(CFG = CFG) +
      # 图例**必须留着**：组合层用 guides="collect" 把 B–E 的同一个分组图例合并成
      # 底部一份；若这里写 "none"，patchwork 收不到该 guide，全图就完全没有图例。
      ggplot2::theme(legend.position = "bottom")
    ux <- ux_list[[f]]
    if (!is.null(ux) && length(ux) <= 10)
      p <- p + ggplot2::scale_x_continuous(breaks = sort(unique(ux)))
    p_list[[f]] <- p
  }

  # ---------------------------------------------------------------------------
  # 3b. 可选备用面板：标准化 β 森林图（CFG$spare$std_beta；**默认关**）
  # ---------------------------------------------------------------------------
  # 通用化：对 focal 里每个变量拟合 response ~ focal_i + covariates [+ random]，
  # 取【z 化后】的标准化系数 β / SE / P，画一张森林图。by 非空时按该因子分组
  # 各拟合一套，做并列森林图。该分支与上面的主图完全独立，enabled = FALSE 时
  # 一行图代码都不执行，对 50 的现有数值与图件零影响。
  spare_spec <- CFG$spare$std_beta %||% list()
  spare_on   <- isTRUE(spare_spec$enabled)
  spare_beta <- NULL
  p_spare    <- NULL
  spare_note <- if (spare_on) character(0) else
    "可选备用面板 CFG$spare$std_beta 未开启（enabled = FALSE），本次未出图。"
  if (spare_on) {
    sp_y    <- spare_spec$response %||% y_var
    sp_foc  <- as.character(spare_spec$focal %||% character(0))
    sp_foc  <- unique(sp_foc[nzchar(sp_foc)])
    sp_cov  <- as.character(spare_spec$covariates %||% character(0))
    sp_cov  <- unique(sp_cov[nzchar(sp_cov)])
    sp_rand <- spare_spec$random %||% random
    sp_by   <- spare_spec$by %||% NULL
    if (!length(sp_foc)) {
      spare_note <- c(spare_note,
        "CFG$spare$std_beta$enabled = TRUE 但 focal 为空，未出图。")
      spare_on <- FALSE
    }
    # 逐变量解析（复用主分支的 local_resolve；响应变量必须在 analysis 表里）
    sp_need <- unique(c(sp_y, sp_foc, sp_cov, if (!is.null(sp_by)) sp_by else NULL))
    sp_val  <- list(); sp_miss <- character(0)
    # 主分支的解析结果（cv$res，含响应变量）优先；其余走同一个 local_resolve()
    sp_base <- cv$res
    if (is.null(sp_base)) sp_base <- list()
    for (v in sp_need) {
      r <- if (v %in% names(sp_base)) sp_base[[v]] else local_resolve(v)
      if (is.null(r)) sp_miss <- c(sp_miss, v) else sp_val[[v]] <- r$value
    }
    if (length(sp_miss)) {
      spare_note <- c(spare_note, paste0(
        "CFG$spare$std_beta 引用的变量无法解析，已跳过备用面板：",
        paste(sp_miss, collapse = ", ")))
      spare_on <- FALSE
    }
    if (spare_on) {
      sp_use_re <- !is.null(sp_rand) && nzchar(as.character(sp_rand)) && has_site
      sp_by_lv  <- if (!is.null(sp_by) && sp_by %in% names(sp_val))
        sort(unique(as.character(sp_val[[sp_by]]))) else "All"
      sp_rows <- list()
      for (bl in sp_by_lv) {
        sel <- if (identical(sp_by_lv, "All")) seq_along(ids)
               else which(as.character(sp_val[[sp_by]]) == bl)
        for (f in sp_foc) {
          keep_v <- unique(c(sp_y, f, sp_cov))
          spd <- data.frame(.row = sel)
          for (v in keep_v) spd[[v]] <- sp_val[[v]][sel]
          if (sp_use_re) spd[[site_var]] <- as.character(d[[site_var]])[sel]
          spd <- spd[stats::complete.cases(spd[, setdiff(names(spd), ".row"), drop = FALSE]), ,
                     drop = FALSE]
          if (nrow(spd) < 6 || length(unique(stats::na.omit(spd[[f]]))) < 2) {
            spare_note <- c(spare_note, sprintf(
              "备用面板 [%s | %s]：完整样本不足或焦点变量无变化，跳过。", bl, f))
            next
          }
          # 站点级共线守卫（复用 md_within_unique）：焦点与站点随机效应共线时去随机效应
          sp_re <- sp_use_re
          if (sp_use_re && !is.null(site_var) && site_var %in% names(spd)) {
            u <- tryCatch(md_within_unique(spd, f, site_var), error = function(e) NA_integer_)
            if (!is.na(u) && u <= 1) {
              sp_re <- FALSE
              spare_note <- c(spare_note, sprintf(
                "备用面板 [%s | %s]：焦点站内唯一值 = 1（站点级），已去掉随机效应避免共线。",
                bl, f))
            }
          }
          # z 化（CFG$model$standardize；响应变量一并 z 化 → 系数为标准化 β）
          spd_z <- spd
          if (!sp_y %in% names(spd_z)) {
            spare_note <- c(spare_note, sprintf(
              "备用面板 [%s | %s]：响应变量 %s 未进入建模数据，跳过。", bl, f, sp_y))
            next
          }
          for (v in keep_v) if (isTRUE(std)) spd_z[[v]] <- local_z(spd_z[[v]])
          if (isTRUE(std)) spd_z[[sp_y]] <- local_z(spd_z[[sp_y]])
          rhs <- c(f, sp_cov, if (sp_re) sp_rand else NULL)
          fml <- stats::as.formula(paste(sp_y, "~", paste(rhs, collapse = " + ")))
          muffle <- function(e) withCallingHandlers(
            e, warning = function(w) invokeRestart("muffleWarning"))
          fit <- if (sp_re)
            tryCatch(muffle(lmerTest::lmer(fml, data = spd_z, REML = reml)),
                     error = function(e) NULL)
          else NULL
          if (is.null(fit))
            fit <- tryCatch(stats::lm(fml, data = spd_z), error = function(e) NULL)
          if (is.null(fit)) {
            spare_note <- c(spare_note, sprintf("备用面板 [%s | %s]：模型拟合失败。", bl, f))
            next
          }
          co <- as.data.frame(summary(fit)$coefficients)
          cn <- names(co)
          pcol <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "p.value", "P"), cn)[1]
          if (is.na(pcol)) next
          sp_rows[[length(sp_rows) + 1L]] <- data.frame(
            group = bl, focus = f, term = rownames(co),
            estimate = as.numeric(co[["Estimate"]]),
            se = as.numeric(co[["Std. Error"]]),
            statistic = if ("t value" %in% cn) as.numeric(co[["t value"]])
                        else if ("z value" %in% cn) as.numeric(co[["z value"]]) else NA_real_,
            p = as.numeric(co[[pcol]]),
            n = stats::nobs(fit),
            estimator = if (inherits(fit, "merMod")) "LMM (standardized)" else "LM (standardized)",
            stringsAsFactors = FALSE)
        }
      }
      if (length(sp_rows)) {
        spare_beta <- do.call(rbind, sp_rows)
        spare_beta <- spare_beta[spare_beta$term != "(Intercept)", , drop = FALSE]
      }
      if (!is.null(spare_beta) && nrow(spare_beta)) {
        spare_beta$sig <- ifelse(is.na(spare_beta$p), "",
          ifelse(spare_beta$p < 0.001, "***", ifelse(spare_beta$p < 0.01, "**",
                 ifelse(spare_beta$p < 0.05, "*", "ns"))))
        spare_beta$lo <- spare_beta$estimate - 1.96 * spare_beta$se
        spare_beta$hi <- spare_beta$estimate + 1.96 * spare_beta$se
        spare_beta$Label <- paste0(spare_beta$focus, ": ", spare_beta$term,
                                   if (length(sp_by_lv) > 1) paste0("  [", spare_beta$group, "]") else "")
        spare_beta$Col <- ifelse(spare_beta$sig %in% c("*", "**", "***"),
                                 md_pal_cont(CFG, "diverging")[3], "grey60")
        p_spare <- md_forest(spare_beta, term = "Label", estimate = "estimate",
                             lo = "lo", hi = "hi", p = "p", ref = 0,
                             xlab = "Standardized beta",
                             tag = tag_of(2L + length(p_list)), pt = pt)
        # 点色 = 方向×显著性（md_pal_cont(diverging)：正 [3] / 不显著 grey60），
        # 与 A 面板同口径；md_forest() 的默认灰点用白色小点遮盖后再叠彩点。
        p_spare <- p_spare +
          ggplot2::geom_point(ggplot2::aes(fill = .data$Col), size = MD_PT_SIZE$scatter * 2.0,
                              shape = 21, colour = "white", stroke = 0.3) +
          ggplot2::scale_fill_identity() +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = pt$tick))
        spare_note <- c(spare_note, sprintf(
          "备用面板已出图：response = %s；focal = %s；covariates = %s；random = %s；by = %s。",
          sp_y, paste(sp_foc, collapse = ", "),
          if (length(sp_cov)) paste(sp_cov, collapse = ", ") else "无",
          if (sp_use_re) as.character(sp_rand) else "无（普通回归）",
          if (is.null(sp_by)) "无（不分亚组）" else sp_by))
      } else {
        spare_note <- c(spare_note, "CFG$spare$std_beta 开启但没有任何模型拟合成功，未出图。")
      }
    }
  }
  if (length(p_list)) {
    # A 全宽；B 起一行并列（对原脚本 Fig5 的 design "AAAA\nBCDE" 口径）。
    # patchwork 要求同一字母占的格子必须是矩形，故这行按面板数排布：
    #   n < 4  -> 该行就是 n 格宽的网格；n >= 4 -> 按 4 列网格换行（每格 1 个面板）。
    n_curve <- length(p_list)
    ncol_g  <- if (n_curve < 4) n_curve else 4L
    row2 <- paste(rep(LETTERS[2:(n_curve + 1)], length.out = n_curve), collapse = "")
    if (nchar(row2) < ncol_g)
      row2 <- paste0(row2, paste(rep(LETTERS[n_curve + 1], ncol_g - nchar(row2)), collapse = ""))
    a_row <- paste(rep("A", ncol_g), collapse = "")
    plot_list <- c(list(A = p_forest), p_list)
    design_str <- paste0(a_row, "\n", row2)
    if (!is.null(p_spare)) {
      # 备用面板另起一行，行宽与上一行一致（patchwork design 要求每行等宽）
      plot_list[["SPARE"]] <- p_spare
      design_str <- paste0(design_str, "\n", paste(rep("S", ncol_g), collapse = ""))
    }
    fig <- md_compose(CFG, plot_list, design = design_str,
                      # collect：B–E 四个面板共用同一个 Vegetation 图例，只出底部一份
                      guides = "collect", legend = "bottom")
  } else {
    fig <- md_compose(CFG, list(A = p_forest), design = "AA\nAA",
                      guides = "keep", legend = NULL)
  }
  md_save_figure(CFG, fig, "50_lmm", height = 130)

  # ---------------------------------------------------------------------------
  # 4. 表
  # ---------------------------------------------------------------------------
  coef_out <- coef_tab[, c("focus", "term", "is_focus", "estimate", "se", "se_robust",
                           "df", "statistic", "p", "p_robust", "P_BH", "sig",
                           "lo", "hi", "estimator", "n", "singular", "note")]
  sheets <- list(
    Coefficients = coef_out,
    Fit = fit_tab,
    Nonlinearity = nl_tab,
    Collinearity_guard = guard_tab,
    Variables = prov_tab,
    Note = data.frame(Item = c(
      "响应变量", "焦点变量", "协变量", "随机效应", "估计方法", "P 口径", "P 校正",
      "稳健 SE", "非线性判据", "条件效应口径", "站点级共线守卫", "降级变量", "退化处理"),
      Value = c(y_var, paste(focus_ok, collapse = ", "), paste(covs_spec, collapse = ", "),
                if (use_re) as.character(random) else "无（退化为普通回归）",
                paste(unique(fit_tab$estimator), collapse = "; "),
                if (fam_name == "gaussian") "lmerTest Satterthwaite t 检验（GLMM 为 Wald z）"
                else paste0("Wald z（family = ", fam_name, "）"),
                paste0(CFG$model$p_adjust %||% "BH", "（对同一模型内全部项一起校正）"),
                if (robust) "已开启：站点级整簇 CR1 三明治（df = 簇数 - 1），与 Satterthwaite 并列"
                else "未开启（CFG$model$robust_se = FALSE）",
                paste0("edf > ", edf_thr, " 判为非线性并画 GAMM 平滑曲线，否则退回直线"),
                paste0("全变量 z 化（standardize = ", std, "）-> ", grid_n,
                       " 点网格 -> 其余协变量置 z = 0 -> 预测 +/- 1.96SE -> 反变换回原始单位"),
                "焦点站内唯一值 = 1 -> 剔站点级协变量 + 禁平滑（见 Collinearity_guard）",
                if (length(dropped_var)) paste(dropped_var, collapse = ", ") else "无",
                "随机效应方差 ~ 0 / singular fit / 不满秩：照常报系数并在 Fit$note 标注，不中断"),
      stringsAsFactors = FALSE))
  if (length(curve_rows) && isTRUE(CFG$table$detail))
    sheets$Curve_data <- do.call(rbind, curve_rows)
  sheets$Note <- rbind(
    sheets$Note,
    data.frame(Item = "可选备用面板（CFG$spare$std_beta）",
               Value = paste(spare_note, collapse = " | "),
               stringsAsFactors = FALSE))
  if (!is.null(spare_beta) && nrow(spare_beta))
    sheets$Std_beta_optional <- spare_beta[, c("group", "focus", "term", "estimate", "se",
                                               "statistic", "p", "sig", "n", "estimator"),
                                           drop = FALSE]
  sheets <- Filter(function(x) !is.null(x) && nrow(x) > 0, sheets)
  tp <- md_write_xlsx(CFG, "50_lmm", sheets)
  md_log(CFG, "INFO", "50: ", length(focus_ok), " 个焦点模型 | 系数 ", nrow(coef_out),
         " 行 | singular ", sum(fit_tab$singular, na.rm = TRUE), " 个")
  md_result(module, "OK", tables = basename(tp), figures = "figures/50_lmm.pdf")
}
