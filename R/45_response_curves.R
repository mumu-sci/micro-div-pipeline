# =============================================================================
# R/45_response_curves.R ---- 模块 45：GAM / GAMM 响应曲线
# -----------------------------------------------------------------------------
# 参考实现（§6 GAM 非线性响应，原图）
#     契约见 _build-notes/micro-div-pipeline/ROUND5_MODULE_SPEC.md
#
# 导出：md_response_curves(CFG)
#
# 通用化（相对原文的关键区别）：
#   · 原文写死 gam_data 里的 PCoA1 / Shannon 两个响应 + RF top 3 预测变量；
#     这里 y 可以是**任意连续列**（排序轴 / α 指数 / 环境量 / 响应变量），
#     x 是任意列名向量（NULL -> CFG$env$vars）；x_from = "rf_top" 时改用
#     模块 40 / 71 的随机森林重要性前 max_panels 个变量（读不到就回退 x）。
#   · 每个变量调一次共享 md_gam_curve()（REML，k = CFG$response_curves$k）；
#     random 非空且 engine = "gam" 时共享函数内部自动走 mgcv::gamm。
#   · 只保留 p_smooth <= p_max 的变量（原文 1320 行口径）；
#     **全部都超阈值时保留 P 最小的那个**（避免出空图），并在 notes 说明。
#
# 图：figures/45_response_curves.pdf（183 × 100 mm）
#     A = 各变量的响应曲线（每变量一格，ncol 自适应）
#     B = EDF / P 值条形（横向，按 edf 排序）
# 表：results/45_response_curves.xlsx（Curve_fits / Curve_grid / Note）
# =============================================================================

# ---------------------------------------------------------------------------
# local_* 辅助（只在本模块内使用）
# ---------------------------------------------------------------------------

#' 图上文字一律纯 ASCII
local_rc_ascii <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- suppressWarnings(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"))
  x[is.na(x)] <- ""
  gsub("[^ -~]", "_", x)
}

#' α 指数（Shannon / Observed / Simpson / Pielou）——与 50_lmm 的本地口径一致，
#' 直接由 05 核心矩阵算，不依赖模块 11 的产物。
local_rc_alpha <- function(CFG) {
  m <- tryCatch(md_core_asv(CFG, "alpha"), error = function(e) NULL)
  if (is.null(m) || !ncol(m)) return(NULL)
  rel <- md_rel(m)
  sh <- -colSums(rel * log(ifelse(rel > 0, rel, 1)))
  ob <- colSums(m > 0)
  data.frame(sample_id = colnames(m),
             Shannon = as.numeric(sh), Observed = as.numeric(ob),
             Simpson = as.numeric(colSums(rel^2)),
             Pielou  = as.numeric(sh / log(pmax(ob, 2))),
             stringsAsFactors = FALSE)
}

#' PCoA 轴（PCoA1..PCoA3）——从 13 号距离缓存现算，与 50_lmm 的 local_pcoa1() 同源。
#' 取不到距离 / 没有 ape 时退回 cmdscale；再不行返回 NULL。
local_rc_pcoa <- function(CFG) {
  metric <- CFG$response_models$pair_metric %||%
    (CFG$distances$metrics %||% "bray")[1] %||% "bray"
  d <- tryCatch(md_dist_get(CFG, metric), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  M <- as.matrix(d)
  if (nrow(M) < 3) return(NULL)
  if (requireNamespace("ape", quietly = TRUE)) {
    p <- tryCatch(ape::pcoa(M, correction = "cailliez"), error = function(e) NULL)
    if (!is.null(p) && !is.null(p$vectors) && ncol(p$vectors) >= 1) {
      v <- p$vectors
      out <- list()
      for (k in seq_len(min(3L, ncol(v)))) {
        out[[paste0("PCoA", k)]] <- stats::setNames(as.numeric(v[, k]), rownames(v))
      }
      return(out)
    }
  }
  cm <- tryCatch(stats::cmdscale(M, k = min(2L, nrow(M) - 1L)), error = function(e) NULL)
  if (is.null(cm) || !ncol(cm)) return(NULL)
  out <- list()
  for (k in seq_len(ncol(cm)))
    out[[paste0("PCoA", k)]] <- stats::setNames(as.numeric(cm[, k]), rownames(M))
  out
}

#' 变量解析：analysis 表 -> 40 号派生变量 -> 本地 α 指数 -> 13 号距离的 PCoA 轴
local_rc_resolve <- function(CFG, ana, ids, v, ai = NULL, po = NULL) {
  if (is.null(v) || !nzchar(v)) return(NULL)
  if (v %in% names(ana)) {
    x <- suppressWarnings(as.numeric(ana[[v]]))
    if (any(is.finite(x)))
      return(list(value = x, source = "analysis table (meta / env / derived)"))
  }
  dv <- tryCatch(md_derived(CFG, v, ids), error = function(e) NULL)
  if (!is.null(dv) && any(is.finite(dv)))
    return(list(value = as.numeric(dv), source = "md_derived() (module 40)"))
  if (!is.null(ai) && v %in% names(ai))
    return(list(value = ai[[v]][match(ids, ai$sample_id)],
                source = "alpha index computed from 05 core (alpha normalisation)"))
  if (!is.null(po) && v %in% names(po))
    return(list(value = as.numeric(po[[v]][match(ids, names(po[[v]]))]),
                source = "PCoA axis computed from 13 distance cache"))
  NULL
}

#' 取"随机森林重要性前 n 个变量"，对应参考实现
#' `gam_preds <- names(sort(rf_A_pcoa1$variable.importance, decreasing = TRUE))[1:3]`
#' （原脚本）。
#'
#' 三级来源，先到先用：
#'   1. 模块 40 落盘的 `cache/40_env.rds`（`env$rf$enabled = TRUE` 时才有）；
#'   2. 模块 71 落盘的 `cache/71_response.rds`；
#'   3. **都没有就自己跑一次** `md_rf_fit()`（与原脚本同一引擎、同一条流水线）——
#'      这样 45 号是自足的，不依赖用户是否顺手开了 40 号的 RF 面板。
local_rc_rf_top <- function(CFG, n, y = NULL, ana = NULL) {
  pick <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (is.data.frame(obj) && "Variable" %in% names(obj) &&
        any(c("Importance", "Importance_scaled") %in% names(obj))) return(obj)
    if (is.list(obj)) {
      for (k in c("rf", "rf_importance", "importance", "imp")) {
        z <- obj[[k]]
        if (is.null(z)) next
        if (is.list(z) && !is.null(z$importance)) z <- z$importance
        if (is.data.frame(z) && "Variable" %in% names(z) &&
            any(c("Importance", "Importance_scaled") %in% names(z))) return(z)
      }
    }
    NULL
  }
  imp <- pick(tryCatch(md_env_cache(CFG), error = function(e) NULL))
  if (is.null(imp)) {
    f <- md_path(CFG, "cache", "71_response.rds")
    if (file.exists(f)) imp <- pick(tryCatch(readRDS(f), error = function(e) NULL))
  }
  if (is.null(imp) && !is.null(y) && !is.null(ana) && y %in% names(ana)) {
    rfc <- CFG$response_curves$rf %||% list()
    fit <- tryCatch(
      md_rf_fit(CFG, ana, y, as.character(CFG$env$vars %||% character(0)),
                engine = as.character(rfc$engine %||% "ranger"),
                num_trees = as.integer(rfc$num_trees %||% 2000L),
                importance = "permutation", seed_offset = 8420L),
      error = function(e) NULL)
    if (!is.null(fit)) {
      imp <- fit$importance
      attr(imp, "source") <- paste0("md_rf_fit(", fit$engine, ") computed inside module 45")
    }
  }
  if (is.null(imp) || !nrow(imp)) return(NULL)
  col <- if ("Importance_scaled" %in% names(imp)) "Importance_scaled" else "Importance"
  imp <- imp[order(-imp[[col]]), , drop = FALSE]
  out <- utils::head(as.character(imp$Variable), as.integer(n))
  attr(out, "source") <- attr(imp, "source") %||%
    "random-forest importance from cache (module 40 / 71)"
  out
}

# ---------------------------------------------------------------------------
# 模块主体
# ---------------------------------------------------------------------------
md_response_curves <- function(CFG) {
  module <- "45"
  if (!isTRUE(CFG$response_curves$enabled %||% TRUE))
    return(md_result(module, "SKIP", "response_curves$enabled = FALSE"))
  if (!md_need(CFG, "gam", module))
    return(md_result(module, "SKIP", "缺少 mgcv"))

  smp <- CFG$data$sample_col
  ana <- tryCatch(md_core_analysis(CFG), error = function(e) NULL)
  if (is.null(ana) || nrow(ana) < 10 || !smp %in% names(ana)) {
    msg <- "分析总表不可用（先跑模块 05；或样本数 < 10）"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }
  ids <- as.character(ana[[smp]])

  notes <- character(0)

  # ---- 分析段 ----
  md_stage(CFG, "45_analysis", {
    md_seed(CFG, md_offset(CFG, "response_curves", 8410L))
    rc_notes <- character(0)

    rc_y          <- as.character(CFG$response_curves$y %||% "PCoA1")
    rc_engine     <- tolower(as.character(CFG$response_curves$engine %||% "gam"))
    rc_k          <- as.integer(CFG$response_curves$k %||% 4L)
    rc_random     <- CFG$response_curves$random
    rc_grid       <- as.integer(CFG$response_curves$grid_n %||% 100L)
    rc_edf_nonlin <- CFG$response_curves$edf_nonlin %||% 1.1
    rc_pmax       <- CFG$response_curves$p_max %||% 0.1
    rc_maxp       <- max(1L, as.integer(CFG$response_curves$max_panels %||% 6L))
    rc_xfrom      <- as.character(CFG$response_curves$x_from %||% "")
    rc_xs         <- as.character(CFG$response_curves$x %||% CFG$env$vars %||% character(0))
    if (!rc_engine %in% c("gam", "gamm", "lm")) {
      rc_notes <- c(rc_notes, paste0("engine '", rc_engine, "' 不可识别 —— 已退回 gam"))
      rc_engine <- "gam"
    }

    # ---- 自变量来源：x_from = "rf_top" 时**在因变量解析之后**再算 ----
    #   原脚本流水线：随机森林给环境变量排序 -> 取前 3 -> 逐个拟合 GAM
    #   （原脚本）。RF 需要因变量那一列，所以放到下面 rc_d 建好之后。
    rc_use_rf <- identical(rc_xfrom, "rf_top")


    # ---- 变量解析 ----
    rc_ai <- tryCatch(local_rc_alpha(CFG), error = function(e) NULL)
    rc_po <- tryCatch(local_rc_pcoa(CFG), error = function(e) NULL)
    rc_y_src <- NA_character_
    ry <- local_rc_resolve(CFG, ana, ids, rc_y, rc_ai, rc_po)
    if (!is.null(ry)) rc_y_src <- ry$source

    rc_d <- data.frame(.sid = ids, stringsAsFactors = FALSE)
    if (!is.null(ry)) rc_d[[rc_y]] <- as.numeric(ry$value)

    # ---- x_from = "rf_top"：现在才有多一列因变量的表可供随机森林使用 ----
    if (rc_use_rf) {
      rf_df <- ana
      if (!is.null(ry)) rf_df[[rc_y]] <- as.numeric(ry$value)
      rf <- local_rc_rf_top(CFG, rc_maxp, y = if (is.null(ry)) NULL else rc_y, ana = rf_df)
      if (length(rf)) {
        rc_xs <- rf
        rc_notes <- c(rc_notes, sprintf(
          "x_from = 'rf_top'（对齐原脚本）：取随机森林重要性前 %d 个变量（%s）；来源 = %s。",
          length(rf), paste(rf, collapse = ", "),
          attr(rf, "source") %||% "unknown"))
      } else {
        rc_notes <- c(rc_notes, paste0(
          "x_from = 'rf_top' 但随机森林排不出序（因变量无法解析 / 样本不足 / ",
          "ranger 与 randomForest 都不可用）—— 已回退用 response_curves$x。"))
      }
    }
    if (!length(rc_xs))
      rc_notes <- c(rc_notes, "没有可用的自变量（x / env$vars / rf_top 均为空）")
    rc_xs <- unique(rc_xs[nzchar(rc_xs)])
    rc_used_xs <- character(0)
    rc_miss_xs <- character(0)
    for (xv in rc_xs) {
      if (identical(xv, rc_y)) {
        rc_notes <- c(rc_notes, paste0("自变量 '", xv, "' 与因变量同名 —— 已跳过"))
        next
      }
      r <- local_rc_resolve(CFG, ana, ids, xv, rc_ai, rc_po)
      if (is.null(r)) { rc_miss_xs <- c(rc_miss_xs, xv); next }
      rc_d[[xv]] <- as.numeric(r$value)
      rc_used_xs <- c(rc_used_xs, xv)
    }
    if (length(rc_miss_xs))
      rc_notes <- c(rc_notes, paste0("自变量无法解析，已剔除：", paste(rc_miss_xs, collapse = ", ")))

    # ---- 逐变量拟合 ----
    rc_fits <- list(); rc_curves <- list(); rc_points <- list()
    rc_nfail <- 0L
    if (!is.null(ry) && length(rc_used_xs)) {
      for (xv in rc_used_xs) {
        g <- tryCatch(
          md_gam_curve(CFG, rc_d, rc_y, xv, engine = rc_engine, k = rc_k,
                       random = rc_random, grid_n = rc_grid, edf_nonlin = rc_edf_nonlin),
          error = function(e) {
            rc_notes <<- c(rc_notes, paste0("'", xv, "' 拟合报错：", conditionMessage(e)))
            NULL
          })
        if (is.null(g)) { rc_nfail <- rc_nfail + 1L; next }
        rc_fits[[xv]]   <- g$fit
        rc_curves[[xv]] <- g$curve
        rc_points[[xv]] <- g$points
      }
    }
    if (rc_nfail)
      rc_notes <- c(rc_notes, sprintf("%d 个自变量拟合失败（样本不足 / 常数列）已跳过。", rc_nfail))

    # ---- 只留 p_smooth <= p_max；全超阈值时保 P 最小的一个 ----
    rc_fit_tab <- NULL; rc_curve_tab <- NULL; rc_use <- character(0); rc_ncol <- 1L
    if (length(rc_fits)) {
      ft <- do.call(rbind, rc_fits)
      rownames(ft) <- NULL
      # 全部候选（含没通过 p_max 的）留一份，写进 Curve_candidates sheet ——
      # 被跳过的变量与它们的 P 必须能查到（只画一个时尤其重要）。
      rc_cand_tab <- ft
      keep <- !is.na(ft$p_smooth) & ft$p_smooth <= rc_pmax
      rc_fallback <- FALSE
      if (!any(keep)) {
        # 原脚本口径（原脚本 `if (gr$p_smooth > 0.1) next`）是**不画**。
        # 但一个都不画会得到一张空图，读者看不到"跑过了、P 是多少"。
        # 与模块 83 的 `min_sig` 回退同一口径：**显示 P 最小的一个，并如实标注它没通过阈值**
        # （每格的 R2/P 注释本来就在图面上），同时在 Note 里写明这是回退。
        w <- if (any(!is.na(ft$p_smooth))) which.min(ifelse(is.na(ft$p_smooth), Inf, ft$p_smooth)) else 1L
        keep[w] <- TRUE
        rc_fallback <- TRUE
        ptxt <- paste(sprintf("%s = %.4f", ft$term,
                              ifelse(is.na(ft$p_smooth), NA_real_, ft$p_smooth)),
                      collapse = "; ")
        rc_notes <- c(rc_notes, paste0(
          "所有变量的 p_smooth 都 > ", format(rc_pmax), "（原脚本此处不画）—— ",
          "按与模块 83 相同的回退口径，只显示 P 最小的 '", ft$term[w],
          "'（P = ", format(ft$p_smooth[w]), "），图上注释如实给出该 P；各变量 P 值：", ptxt))
      }
      sel <- ft$term[keep]
      sel <- sel[order(match(sel, rc_used_xs))]          # 保持用户 / RF 给定的顺序
      drop_n <- max(0L, length(sel) - rc_maxp)
      if (drop_n > 0)
        rc_notes <- c(rc_notes, sprintf("变量数超过 max_panels = %d，已按给定顺序截断 %d 个。",
                                        rc_maxp, drop_n))
      rc_use <- utils::head(sel, rc_maxp)
      # 一个都没通过 p_max 时，把**全部**拟合结果留在表里（下游出说明面板 + PARTIAL），
      # 而不是把它清空成 SKIP —— 用户要能看到"跑过了、P 值是多少"。
      rc_fit_tab <- if (length(rc_use)) ft[match(rc_use, ft$term), , drop = FALSE] else ft
      rownames(rc_fit_tab) <- NULL
      rc_curve_tab <- if (length(rc_use)) {
        do.call(rbind, lapply(rc_use, function(t)
          cbind(term = t, rc_curves[[t]], stringsAsFactors = FALSE)))
      } else NULL
      if (!is.null(rc_curve_tab)) rownames(rc_curve_tab) <- NULL
      rc_ncol <- if (length(rc_use) <= 1) 1L else if (length(rc_use) == 2) 2L else 3L
      rc_notes <- c(rc_notes, sprintf(
        "保留 %d / %d 个变量出图（p_smooth <= %s）；engine = %s，k = %d，edf_nonlin = %s，grid_n = %d。",
        length(rc_use), length(ft$term), format(rc_pmax), rc_engine, rc_k,
        format(rc_edf_nonlin), rc_grid))
    } else {
      rc_notes <- c(rc_notes, "没有任何变量拟合成功 —— 本模块将 SKIP。")
    }
    rc_y_used <- rc_y
  })
  if (length(rc_notes)) notes <- c(notes, rc_notes)

  if (is.null(rc_fit_tab) || !nrow(rc_fit_tab)) {
    msg <- if (length(rc_notes)) paste(rc_notes, collapse = " | ") else
      "没有任何响应曲线可以拟合"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", "response curves unavailable", notes = rc_notes))
  }

  # ---------------------------------------------------------------------------
  # 出图：只出原脚本那一组 GAM 响应曲线（原脚本）
  #   ⚠ 2026-09-26（第六轮）：原实现里的 **B 面板（EDF / -log10(P) 横向条形）**
  #     是自造的，原代码没有这张图 —— 已删除。本模块现在只有这一组曲线。
  # ---------------------------------------------------------------------------
  pt    <- md_pt(CFG)
  style <- CFG$response_curves$style %||% list()

  # 没有任何变量的 p_smooth 通过 p_max：按原脚本就是"不画"。
  # 这里出一个说明用的空面板并把 P 值写进 notes，模块状态记 PARTIAL（不是 OK）。
  if (!length(rc_use)) {
    fig <- md_compose(CFG, list(A = md_empty_plot(
      CFG, paste0("No smooth with p <= ", format(rc_pmax)))))
    md_save_figure(CFG, fig, "45_response_curves", width = 90, height = 60)
    fit_out0 <- rc_fit_tab %||% data.frame()
    tp0 <- md_write_xlsx(CFG, "45_response_curves", list(
      Curve_fits = if (nrow(fit_out0)) fit_out0 else NULL,
      Note = data.frame(Note = c(
        "口径参考实现（原图 GAM 响应曲线）。",
        sprintf("y = %s；engine = %s；k = %d。", rc_y_used, rc_engine, rc_k),
        sprintf("筛选口径照原脚本：只画 p_smooth <= %s 的平滑项；本次一个都没有通过。",
                format(rc_pmax)),
        notes)),
      meta_extra = list(y = rc_y_used, engine = rc_engine, k = rc_k, n_terms = 0)))
    return(md_result(module, "PARTIAL",
                     paste0("无平滑项通过 p <= ", format(rc_pmax), "，按原脚本口径不出曲线"),
                     tables = basename(tp0), figures = "figures/45_response_curves.pdf",
                     notes = notes))
  }

  panelsA <- lapply(rc_use, function(t) {
    f <- rc_fit_tab[rc_fit_tab$term == t, , drop = FALSE]
    g <- md_gam_curve_plot(
      curve = rc_curves[[t]], points = rc_points[[t]],
      title = local_rc_ascii(sprintf("%s ~ s(%s)", rc_y_used, t)),
      xlab = local_rc_ascii(t), ylab = local_rc_ascii(rc_y_used),
      r2 = f$r2[1], p = f$p_smooth[1],
      line_col = style$line_col %||% "#0072B2",
      ribbon_fill = style$ribbon_fill %||% "grey80",
      ribbon_alpha = style$ribbon_alpha %||% 0.5,
      point_alpha = style$point_alpha %||% 0.3,
      pt = pt)
    # 回退情形：在面板内补一行醒目说明（纯 ASCII），读者不会把不显著的曲线当结论
    if (isTRUE(rc_fallback)) {
      g <- g + ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8,
                                 label = sprintf("p > %.2f (not significant)", rc_pmax),
                                 size = pt2gg(pt$ann * 0.9), colour = "grey35")
    }
    g
  })
  n_panel <- length(panelsA)
  fig <- if (n_panel == 1) {
    md_compose(CFG, list(A = panelsA[[1]] + ggplot2::labs(tag = NULL)))
  } else {
    md_compose(CFG, panelsA, ncol = rc_ncol, tag_levels = NULL)
  }
  fig_h <- if (n_panel <= 3) 70 else if (n_panel <= 6) 130 else 190
  md_save_figure(CFG, fig, "45_response_curves", width = 183, height = fig_h)
  for (i in seq_along(panelsA))
    md_maybe_single(CFG, panelsA[[i]], paste0("45_response_curves_", LETTERS[i]))

  # ---------------------------------------------------------------------------
  # 结果表
  # ---------------------------------------------------------------------------
  fit_out <- rc_fit_tab[, c("term", "engine", "edf", "p_smooth", "r2", "n", "nonlinear"),
                        drop = FALSE]
  names(fit_out) <- c("term", "engine", "edf", "p_smooth", "r2", "n", "nonlinear")

  note_df <- data.frame(Note = c(
    sprintf("口径参考实现（原图 GAM 响应曲线）。"),
    sprintf("因变量 y = %s（来源：%s）；engine = %s；k = %d；REML。",
            rc_y_used, rc_y_src %||% "NA", rc_engine, rc_k),
    sprintf("random = %s（非空且 engine = gam 时共享 md_gam_curve() 内部走 mgcv::gamm；GAMM 的 summary() 不提供 s.pv，故 p_smooth 可能为 NA）。",
            if (is.null(rc_random) || !nzchar(rc_random)) "NULL" else rc_random),
    sprintf("筛选：只画 p_smooth <= %s 的变量（原脚本 口径）；一个都没通过时按原脚本不出曲线，只在图里给说明。",
            format(rc_pmax)),
    sprintf("每格标题 %s ~ s(x)；R2/P 以面板内注释给出（原文 1334-1339 行的版式）。",
            rc_y_used),
    "edf > edf_nonlin 记为非线性（nonlinear 列）。",
    "本模块只有这一组 GAM 曲线面板；原实现里的 EDF / -log10(P) 条形图是自造的，已于 2026-09-26 删除。",
    notes))

  tp <- md_write_xlsx(CFG, "45_response_curves", list(
    Curve_fits = fit_out,
    Curve_candidates = rc_cand_tab,
    Curve_grid = rc_curve_tab,
    Note       = note_df),
    meta_extra = list(y = rc_y_used, engine = rc_engine, k = rc_k,
                      n_terms = nrow(fit_out)))

  md_result(module, "OK", tables = basename(tp),
            figures = "figures/45_response_curves.pdf",
            notes = paste(notes, collapse = " | "))
}
