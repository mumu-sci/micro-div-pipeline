# =============================================================================
# R/52_sem.R ---- 声明式边表 -> piecewiseSEM -> d-separation -> 路径系数 -> 分层路径图
# -----------------------------------------------------------------------------
# 结构**完全来自用户声明**（CFG$sem$edges：from / to / layer / model）：
#   边表为空 -> md_skip("用户未声明 DAG")，**绝不自动探索结构**。
#   layer 决定路径图纵向分层；CFG$sem$layers 给层的上下顺序；
#   CFG$sem$labels 给显示名；CFG$sem$dsep 控制 Fisher's C；CFG$sem$p_adjust
#   控制路径系数多重校正（BH，另列 P_BH）。
#
# 拟合：按 `to` 分组，每个内生变量一个模型（lm 或 lmer），
#       do.call(piecewiseSEM::psem, ...) -> fisherC() d-separation ->
#       summary(sem) 路径系数（含 Std.Estimate）-> 各内生变量 R2。
#       标准化（CFG$sem$standardize，缺省跟随 CFG$model$standardize）后
#       Std.Estimate ~ Estimate，系数可横向比较。
#
# 照实报告而不是硬撑：
#   · 抽样量不足（完整样本 < 10）-> SKIP 并写明；
#   · 边表引用的变量解析不出来（例如 40 号未跑导致 Soil_PC1 缺失）->
#     默认 SKIP 并列出变量名（CFG$sem$on_missing = "drop" 时才降集边表并记录）；
#   · 模型不可识别（共线 / 不满秩 / lmer 失败退化 lm）-> 记录在 Fit 与 Note；
#   · Fisher's C 拒绝（P < 0.05）-> 在 Fit_Interpretation 里写
#     "REJECTED by d-separation"，不得只留一个数字。
#
# 路径图：复用 nifH make_sem_path_nature 的思路 ---- 分层布局、跨层贝塞尔绕行、
#   系数标签避让；正系数 #2166AC、负系数 #B2182B，边上标标准化系数与显著性星号。
#
# 依赖：40（Soil_PC1 等派生变量，只经 md_derived() 取）。
# =============================================================================

md_sem <- function(CFG) {
  module <- "52"
  if (!md_need(CFG, "sem", module))
    return(md_result(module, "SKIP", "缺少 piecewiseSEM"))
  if (!md_need(CFG, "model", module))
    return(md_result(module, "SKIP", "缺少 lme4/lmerTest/performance"))
  if (!md_precheck(CFG, module, need_response = TRUE))
    return(md_result(module, "SKIP", "前置条件不满足"))
  md_load(c("lme4", "lmerTest", "performance"))

  # ---------------------------------------------------------------------------
  # 0. 边表
  # ---------------------------------------------------------------------------
  edges <- CFG$sem$edges
  if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0) {
    md_skip(CFG, module, "用户未声明 DAG（CFG$sem$edges 为空）")
    return(md_result(module, "SKIP", "用户未声明 DAG"))
  }
  need_cols <- c("from", "to")
  if (!all(need_cols %in% names(edges))) {
    md_skip(CFG, module, "CFG$sem$edges 缺少 from / to 列")
    return(md_result(module, "SKIP", "边表列不完整"))
  }
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  if (!"layer" %in% names(edges)) edges$layer <- ""
  if (!"model" %in% names(edges)) edges$model <- "lm"
  edges$from <- as.character(edges$from); edges$to <- as.character(edges$to)
  edges$layer <- as.character(edges$layer); edges$model <- tolower(as.character(edges$model))
  edges <- edges[nzchar(edges$from) & nzchar(edges$to), , drop = FALSE]
  if (nrow(edges) == 0) {
    md_skip(CFG, module, "边表里没有有效边（from/to 为空）")
    return(md_result(module, "SKIP", "边表为空"))
  }
  layers   <- as.character(CFG$sem$layers %||% unique(edges$layer))
  dsep     <- isTRUE(CFG$sem$dsep %||% TRUE)
  p_adj    <- CFG$sem$p_adjust
  std      <- isTRUE(CFG$sem$standardize %||% CFG$model$standardize %||% TRUE)
  on_miss  <- CFG$sem$on_missing %||% "skip"
  smp      <- CFG$data$sample_col
  site_var <- CFG$design$site_var
  random   <- CFG$model$random

  d <- md_analysis_with_derived(CFG)
  if (is.null(d) || !nrow(d)) {
    md_skip(CFG, module, "analysis 表为空（先跑 05 / 40）")
    return(md_result(module, "SKIP", "analysis 表为空"))
  }
  ids <- as.character(d[[smp]])

  # ---------------------------------------------------------------------------
  # 1. 变量解析
  # ---------------------------------------------------------------------------
  cv <- new.env()
  local_z <- function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
  }
  local_alpha <- function() {
    if (!is.null(cv$alpha)) return(cv$alpha)
    m <- md_core_asv(CFG, "alpha"); rel <- md_rel(m)
    sh <- -colSums(rel * log(ifelse(rel > 0, rel, 1))); ob <- colSums(m > 0)
    cv$alpha <- data.frame(sample_id = colnames(m), Shannon = as.numeric(sh),
                           Observed = as.numeric(ob), Simpson = as.numeric(colSums(rel^2)),
                           Pielou = as.numeric(sh / log(pmax(ob, 2))), stringsAsFactors = FALSE)
    cv$alpha
  }
  local_pcoa1 <- function() {
    if (!is.null(cv$pcoa1)) return(cv$pcoa1)
    dd <- md_dist_get(CFG, CFG$response_models$pair_metric %||% CFG$distances$metrics[1])
    if (is.null(dd)) { cv$pcoa1 <- NULL; return(NULL) }
    M <- as.matrix(dd); v <- NULL
    if (requireNamespace("ape", quietly = TRUE)) {
      p <- tryCatch(ape::pcoa(M, correction = "cailliez"), error = function(e) NULL)
      if (!is.null(p)) v <- stats::setNames(p$vectors[, 1], rownames(p$vectors))
    }
    if (is.null(v)) {
      cm <- tryCatch(stats::cmdscale(M, k = 1), error = function(e) NULL)
      if (!is.null(cm)) v <- stats::setNames(cm[, 1], rownames(M))
    }
    cv$pcoa1 <- v; v
  }
  prov <- list(); prov_seen <- character(0)
  local_resolve <- function(v) {
    src <- NULL; val <- NULL
    if (v %in% names(d)) {
      x <- suppressWarnings(as.numeric(d[[v]]))
      if (!all(is.na(x))) { val <- x; src <- "analysis 表（meta / env / 40 派生变量）" }
    }
    if (is.null(val)) {
      dv <- tryCatch(md_derived(CFG, v, ids), error = function(e) NULL)
      if (!is.null(dv) && any(!is.na(dv))) {
        val <- suppressWarnings(as.numeric(dv)); src <- "md_derived()（40 号派生变量）"
      }
    }
    if (is.null(val)) {
      a <- local_alpha()
      if (v %in% names(a)) {
        val <- a[[v]][match(ids, a$sample_id)]
        src <- "本地由 05 核心矩阵（alpha 口径）计算"
      }
    }
    if (is.null(val)) {
      pc <- tryCatch(md_phylo_cache(CFG), error = function(e) NULL)
      if (!is.null(pc) && v %in% names(pc)) {
        val <- suppressWarnings(as.numeric(pc[[v]][match(ids, pc[[1]])]))
        src <- "md_phylo_cache()（21 号）"
      }
    }
    if (is.null(val)) {
      po <- local_pcoa1()
      if (!is.null(po) && v %in% c("PCoA1", "PCoA_1", "PC1")) {
        val <- as.numeric(po[match(ids, names(po))]); src <- "本地由 13 距离矩阵 PCoA 轴 1 计算"
      }
    }
    if (is.null(val)) return(NULL)
    if (!v %in% prov_seen) {
      prov_seen <<- c(prov_seen, v)
      prov[[length(prov) + 1L]] <<- data.frame(variable = v, source = src,
                                               n_used = sum(!is.na(val)), stringsAsFactors = FALSE)
    }
    list(value = val, source = src)
  }

  vars <- unique(c(edges$from, edges$to))
  res  <- list(); missing_vars <- character(0)
  for (v in vars) {
    r <- local_resolve(v)
    if (is.null(r)) missing_vars <- c(missing_vars, v) else res[[v]] <- r
  }
  if (length(missing_vars)) {
    msg <- paste0("边表引用的变量无法解析：", paste(missing_vars, collapse = ", "))
    if (identical(on_miss, "drop")) {
      edges <- edges[!(edges$from %in% missing_vars | edges$to %in% missing_vars), , drop = FALSE]
      md_warn(CFG, "52: ", msg, " -> 已按 CFG$sem$on_missing='drop' 降集边表（剩余 ",
              nrow(edges), " 条边）")
      if (nrow(edges) == 0) {
        md_skip(CFG, module, paste0(msg, "；降集后无剩余边"))
        return(md_result(module, "SKIP", msg))
      }
    } else {
      md_skip(CFG, module, paste0(msg, "（若确有缺失，可显式设 CFG$sem$on_missing='drop' 降集）"))
      return(md_result(module, "SKIP", msg))
    }
  }

  # ---------------------------------------------------------------------------
  # 2. 建模数据（完整样本 + 可选 z 化）
  # ---------------------------------------------------------------------------
  use_re <- !is.null(random) && nzchar(as.character(random)) &&
    !is.null(site_var) && site_var %in% names(d)
  d_use <- data.frame(.row = seq_len(nrow(d)))
  for (v in unique(c(edges$from, edges$to))) d_use[[v]] <- res[[v]]$value
  if (use_re) d_use[[site_var]] <- as.character(d[[site_var]])
  d_use <- d_use[stats::complete.cases(d_use[, setdiff(names(d_use), ".row"), drop = FALSE]), ,
                 drop = FALSE]
  if (use_re) d_use[[site_var]] <- factor(d_use[[site_var]])
  n_use <- nrow(d_use)
  if (n_use < 10) {
    md_skip(CFG, module, paste0("完整样本仅 ", n_use, " 个，SEM 抽样量不足"))
    return(md_result(module, "SKIP", paste0("抽样量不足 (n = ", n_use, ")")))
  }
  centers <- data.frame(variable = unique(c(edges$from, edges$to)), center = NA_real_,
                        scale = NA_real_, stringsAsFactors = FALSE)
  if (std) {
    for (v in centers$variable) {
      centers$center[centers$variable == v] <- mean(d_use[[v]], na.rm = TRUE)
      centers$scale[centers$variable == v]  <- stats::sd(d_use[[v]], na.rm = TRUE)
      d_use[[v]] <- local_z(d_use[[v]])
    }
  }

  # ---------------------------------------------------------------------------
  # 3. 逐内生变量拟合 -> psem
  # ---------------------------------------------------------------------------
  resp_list <- unique(edges$to)
  models <- list(); fit_note <- list(); aliased <- character(0)
  muffle <- function(e) withCallingHandlers(e, warning = function(w) invokeRestart("muffleWarning"))
  for (r in resp_list) {
    pr  <- unique(edges$from[edges$to == r])
    mdls<- unique(edges$model[edges$to == r])
    kind<- if (any(mdls == "lmer")) "lmer" else "lm"
    note <- character(0)
    if (kind == "lmer" && !use_re) {
      kind <- "lm"; note <- c(note, "声明为 lmer 但无可用随机效应（CFG$model$random / site_var）-> 降为 lm")
    }
    rhs <- c(pr, if (kind == "lmer") random else NULL)
    fml <- stats::as.formula(paste(r, "~", paste(rhs, collapse = " + ")))
    fit <- NULL
    if (kind == "lmer") {
      fit <- tryCatch(suppressMessages(muffle(lmerTest::lmer(fml, data = d_use, REML = FALSE))),
                      error = function(e) NULL)
      if (is.null(fit)) { kind <- "lm"; note <- c(note, "lmer 拟合失败 -> 降为 lm") }
    }
    if (is.null(fit)) fit <- tryCatch(stats::lm(fml, data = d_use), error = function(e) NULL)
    if (is.null(fit)) {
      md_skip(CFG, module, paste0("内生变量 ", r, " 的模型无法拟合"))
      return(md_result(module, "SKIP", paste0("模型无法拟合：", r)))
    }
    co <- as.data.frame(summary(fit)$coefficients)
    if (any(is.na(co[["Estimate"]]))) aliased <- c(aliased, rownames(co)[is.na(co[["Estimate"]])])
    if (inherits(fit, "merMod") && isTRUE(tryCatch(lme4::isSingular(fit),
                                                   error = function(e) FALSE)))
      note <- c(note, "singular fit（随机效应方差 ~ 0）")
    models[[r]] <- fit; fit_note[[r]] <- paste(unique(note), collapse = "；")
  }
  if (length(aliased))
    md_warn(CFG, "52: 模型不可识别（存在共线/别名项）：", paste(unique(aliased), collapse = ", "))

  sem <- tryCatch(suppressMessages(do.call(piecewiseSEM::psem, c(models, list(data = d_use)))),
                  error = function(e) NULL)
  if (is.null(sem)) {
    md_skip(CFG, module, "piecewiseSEM::psem 拟合失败（结构或数据不满足要求）")
    return(md_result(module, "SKIP", "psem 拟合失败"))
  }

  # ---- d-separation（Fisher's C） ----
  fc_tab <- NULL; fc <- NA_real_; fc_p <- NA_real_
  if (dsep) {
    fc_tab <- tryCatch(piecewiseSEM::fisherC(sem), error = function(e) NULL)
    if (!is.null(fc_tab) && nrow(fc_tab)) {
      fc   <- suppressWarnings(as.numeric(fc_tab$Fisher.C[1]))
      fc_p <- suppressWarnings(as.numeric(fc_tab$P.Value[1]))
    }
  }

  # ---- 路径系数 ----
  sm <- tryCatch({
    junk <- utils::capture.output(sm0 <- summary(sem, .progressBar = FALSE)); sm0
  }, error = function(e) NULL)
  if (is.null(sm) || is.null(sm$coefficients)) {
    md_skip(CFG, module, "summary(psem) 未返回路径系数")
    return(md_result(module, "SKIP", "summary(psem) 无系数"))
  }
  co <- as.data.frame(sm$coefficients)
  names(co) <- make.names(names(co), unique = TRUE)
  pcol <- intersect(c("P.Value", "P.value", "P"), names(co))[1]
  if (!"Std.Estimate" %in% names(co)) co$Std.Estimate <- co$Estimate
  co$P_BH <- md_p_adjust(co[[pcol]], p_adj)
  co$sig  <- ifelse(is.na(co[[pcol]]), "",
                    ifelse(co[[pcol]] < 0.001, "***", ifelse(co[[pcol]] < 0.01, "**",
                           ifelse(co[[pcol]] < 0.05, "*", "ns"))))
  is_path <- !grepl("~~", co$Predictor) & !grepl("\\(Intercept\\)", co$Predictor) &
    !grepl("~~", co$Response)
  path_tab <- co[is_path, , drop = FALSE]

  # ---- 各内生变量 R2 ----
  r2_tab <- do.call(rbind, lapply(names(models), function(r) {
    m <- models[[r]]
    if (inherits(m, "merMod")) {
      rr <- tryCatch(performance::r2(m), error = function(e) NULL)
      r2c <- as.numeric(rr$R2_conditional %||% NA_real_)
      r2m <- as.numeric(rr$R2_marginal %||% NA_real_)
      # 随机效应方差 ~ 0 时无法分解 R2 -> 取边际 R2 并记明，不留无法解释的 NA
      nt <- fit_note[[r]] %||% ""
      if (is.na(r2c) && is.finite(r2m)) {
        r2c <- r2m
        nt <- paste(c(nt[nzchar(nt)], "随机效应方差 ~ 0，无法分解 -> R2 取边际 R2"), collapse = "；")
      }
      data.frame(response = r, estimator = "lmer", n = stats::nobs(m),
                 R2 = if (is.finite(r2c)) r2c else r2m,
                 R2_marginal = r2m, R2_conditional = r2c,
                 note = nt, stringsAsFactors = FALSE)
    } else {
      sm0 <- tryCatch(summary(m), error = function(e) NULL)
      data.frame(response = r, estimator = "lm", n = stats::nobs(m),
                 R2 = as.numeric(sm0$r.squared %||% NA_real_),
                 R2_marginal = as.numeric(sm0$r.squared %||% NA_real_),
                 R2_conditional = NA_real_,
                 note = fit_note[[r]] %||% "", stringsAsFactors = FALSE)
    }
  }))

  # ---------------------------------------------------------------------------
  # 4. 分层路径图（逐项对 04_Fig04_SEM.R §A/B make_sem_path_nature）
  # ---------------------------------------------------------------------------
  # 布局：层号决定纵向位置；层内按 CFG$sem$layers 顺序左右均分。
  # 节点：**两行标签**（显示名 + "R2 = x.xx" 斜体），框宽由两行里较长的一行决定，
  #       避免 R2 被框压住（原脚本 node_box 口径）。
  # 边：正 #2166AC / 负 #B2182B（只给显著路径上色），不显著 grey40 点线且**不标系数**；
  #     线宽 ∝ |beta|；跨层直线会穿过别的节点框时自动贝塞尔绕行。
  xs_lim <- c(0.35, 9.65); x_lo <- 2.4; x_hi <- 7.6; x_mid <- 5.0
  y_step <- 2.75; y_base <- 0.3
  # 颜色（原脚本字面值，与 12_plot_theme.R 的 MD_PAL_DEFAULT$sem_pos/sem_neg 同值）
  col_pos <- md_pal_cont(CFG, "diverging")[3]   # 正系数
  col_neg <- md_pal_cont(CFG, "diverging")[1]   # 负系数
  col_ns  <- "grey40"                            # 不显著
  lay_of <- function(nm) match(nm, layers)
  nodes <- unique(c(as.character(path_tab$Predictor), as.character(path_tab$Response)))
  if (!length(nodes)) nodes <- unique(c(edges$from, edges$to))
  node_lay <- stats::setNames(rep(NA_integer_, length(nodes)), nodes)
  for (i in seq_len(nrow(edges))) {
    li <- lay_of(edges$layer[i]); if (is.na(li)) li <- 1L
    node_lay[edges$to[i]] <- li
  }
  for (v in nodes) {
    if (is.na(node_lay[[v]])) {
      cons <- lay_of(edges$layer[edges$from == v])
      cons <- cons[!is.na(cons)]
      node_lay[[v]] <- if (length(cons)) max(1L, min(cons) - 1L) else 1L
    }
  }
  n_lay  <- max(node_lay, na.rm = TRUE)
  y_of   <- y_base + (n_lay - node_lay) * y_step
  x_of   <- rep(x_mid, length(nodes)); names(x_of) <- nodes
  for (li in sort(unique(node_lay))) {
    idx <- which(node_lay == li)
    if (!length(idx)) next
    # 层内次序与 CFG$sem$layers / 边表一致（原脚本 SEM_ORDER 的口径）
    ord <- order(match(nodes[idx], unique(c(nodes, edges$from, edges$to))))
    idx <- idx[ord]
    x_of[idx] <- if (length(idx) == 1) x_mid else seq(x_lo, x_hi, length.out = length(idx))
  }

  lab_map <- CFG$sem$labels %||% character(0)
  nd <- data.frame(name = nodes, x = as.numeric(x_of[nodes]), y = as.numeric(y_of[nodes]),
                   stringsAsFactors = FALSE)
  nd$first <- ifelse(nd$name %in% names(lab_map), as.character(lab_map[nd$name]), nd$name)
  nd$r2 <- vapply(nd$name, function(v) {
    i <- match(v, r2_tab$response)
    if (is.na(i) || is.na(r2_tab$R2[i])) "" else sprintf("R2 = %.2f", r2_tab$R2[i])
  }, character(1))
  nd$label <- ifelse(nzchar(nd$r2), paste0(nd$first, "\n", nd$r2), nd$first)
  # 框宽由两行里较长的一行决定（原脚本 node_box / lab_nchar）
  nd$lab_nchar <- pmax(nchar(nd$first), ifelse(nzchar(nd$r2), nchar(nd$r2), 0L))
  nd$hw <- 0.45 + 0.050 * nd$lab_nchar
  nd$hh <- ifelse(nzchar(nd$r2), 0.60, 0.42)
  # 模型拟合注（Fisher's C）当成一个「区块」参与标签避让
  fc_str <- if (is.na(fc)) sprintf("Saturated model   n = %d", n_use)
            else sprintf("Fisher's C = %.1f,  P = %.3f,  n = %d", fc, fc_p, n_use)
  fc_hw <- 2.2

  y_lim <- c(min(nd$y) - 1.55, max(nd$y) + 1.05)
  fc_pt <- c(x_mid, y_lim[1] + 0.45)
  panel_wh <- 3 / 2
  ASPECT <- (1 / panel_wh) * (diff(xs_lim) / diff(y_lim))

  bezier_quad <- function(p0, p1, p2, n = 90L) {
    t <- seq(0, 1, length.out = n)
    cbind(x = (1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1],
          y = (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2])
  }
  resample_path <- function(pts, n = 120L) {
    if (nrow(pts) < 2) return(pts)
    d <- sqrt(rowSums(diff(pts)^2)); s <- c(0, cumsum(d))
    if (s[length(s)] < 1e-9) return(pts[c(1, nrow(pts)), , drop = FALSE])
    s <- s / s[length(s)]; tq <- seq(0, 1, length.out = n)
    cbind(x = stats::approx(s, pts[, 1], tq)$y, y = stats::approx(s, pts[, 2], tq)$y)
  }
  pt_box_dist <- function(px, py, blk) {
    if (is.null(blk) || !nrow(as.matrix(blk))) return(Inf)
    dx <- pmax(0, abs(px - blk[, 1]) - blk[, 3])
    dy <- pmax(0, abs(py - blk[, 2]) - blk[, 4])
    min(sqrt(dx^2 + (dy * ASPECT)^2))
  }
  path_clearance <- function(pts, blk, trim = 0.08) {
    n <- nrow(pts); i0 <- max(1, floor(trim * n)); i1 <- min(n, ceiling((1 - trim) * n))
    if (i1 <= i0) { i0 <- 1L; i1 <- n }
    min(vapply(i0:i1, function(i) pt_box_dist(pts[i, 1], pts[i, 2], blk), numeric(1)))
  }
  path_in_bounds <- function(pts, m = 0.10)
    all(pts[, 1] >= xs_lim[1] + m, pts[, 1] <= xs_lim[2] - m,
        pts[, 2] >= y_lim[1] + m, pts[, 2] <= y_lim[2] - m)
  trim_path_tail <- function(pts, cx, cy, hw, hh, gap = 0.02) {
    inside <- abs(pts[, 1] - cx) <= (hw + gap) & abs(pts[, 2] - cy) <= (hh + gap)
    if (!any(inside)) return(pts)
    keep <- max(2L, min(which(inside)) - 1L)
    pts[seq_len(min(keep, nrow(pts))), , drop = FALSE]
  }
  choose_curvature <- function(p0, p2, blk, target = 0.22,
                               cand_k = seq(-0.40, 0.40, by = 0.02)) {
    L <- sqrt(sum((p2 - p0)^2)); if (L < 1e-9) return(list(k = 0, pts = rbind(p0, p2), cl = Inf))
    mid <- (p0 + p2) / 2; u <- (p2 - p0) / L; nrm <- c(-u[2], u[1])
    cl <- rep(NA_real_, length(cand_k)); ok <- rep(FALSE, length(cand_k))
    paths <- vector("list", length(cand_k))
    for (i in seq_along(cand_k)) {
      pts <- bezier_quad(p0, mid + nrm * (cand_k[i] * L), p2)
      paths[[i]] <- pts; cl[i] <- path_clearance(pts, blk); ok[i] <- path_in_bounds(pts)
    }
    sel <- which(ok & cl >= target)
    i <- if (length(sel)) sel[which.min(abs(cand_k[sel]))] else {
      ok2 <- which(ok)
      if (length(ok2)) ok2[which.max(cl[ok2])] else which.min(abs(cand_k))
    }
    list(k = cand_k[i], pts = paths[[i]], cl = cl[i])
  }

  blocks <- rbind(cbind(nd$x, nd$y, nd$hw, nd$hh),
                  c(fc_pt[1], fc_pt[2], fc_hw, 0.24))
  pe <- data.frame(
    from = as.character(path_tab$Predictor), to = as.character(path_tab$Response),
    est = suppressWarnings(as.numeric(path_tab$Std.Estimate)),
    p = suppressWarnings(as.numeric(path_tab[[pcol]])),
    stringsAsFactors = FALSE)
  pe$est[is.na(pe$est)] <- suppressWarnings(as.numeric(path_tab$Estimate))[is.na(pe$est)]
  pe$sig <- !is.na(pe$p) & pe$p < 0.05
  pe$star <- ifelse(is.na(pe$p), "", ifelse(pe$p < 0.001, "***", ifelse(pe$p < 0.01, "**",
                     ifelse(pe$p < 0.05, "*", ""))))
  # 只有显著路径标系数（原脚本口径），不显著的只留点线
  pe$lab <- ifelse(pe$sig, sprintf("%.2f%s", pe$est, pe$star), NA_character_)
  pe$col <- ifelse(pe$sig, ifelse(pe$est > 0, col_pos, col_neg), col_ns)
  pe$lty <- ifelse(pe$sig, "Significant paths", "Non-significant paths")
  pe$lty <- factor(pe$lty, levels = c("Significant paths", "Non-significant paths"))
  # |beta| -> 线宽（原脚本：显著 0.35 + 0.85 * |b|/max|b|；不显著 0.40）
  mx_b <- suppressWarnings(max(abs(pe$est), na.rm = TRUE))
  if (!is.finite(mx_b) || mx_b <= 0) mx_b <- 1
  pe$lw <- ifelse(pe$sig, 0.35 + 0.85 * (abs(pe$est) / mx_b), 0.40)

  paths <- vector("list", nrow(pe))
  for (i in seq_len(nrow(pe))) {
    si <- match(pe$from[i], nd$name); ti <- match(pe$to[i], nd$name)
    others <- blocks[-c(si, ti), , drop = FALSE]
    cc <- choose_curvature(c(nd$x[si], nd$y[si]), c(nd$x[ti], nd$y[ti]), others)
    paths[[i]] <- trim_path_tail(cc$pts, nd$x[ti], nd$y[ti], nd$hw[ti], nd$hh[ti])
  }
  # ---- 系数标签避让（只给显著路径摆标签） ----
  pe$lx <- NA_real_; pe$ly <- NA_real_
  placed <- matrix(numeric(0), ncol = 4)
  gap_to <- function(cx, cy, hw, hh, blk) {
    if (is.null(blk) || !nrow(as.matrix(blk))) return(Inf)
    dx <- pmax(0, abs(cx - blk[, 1]) - (hw + blk[, 3]) - 0.06)
    dy <- pmax(0, abs(cy - blk[, 2]) - (hh + blk[, 4]) - 0.06)
    min(sqrt(dx^2 + (dy * ASPECT)^2))
  }
  for (i in order(-abs(pe$est))) {
    if (!isTRUE(pe$sig[i]) || is.na(pe$lab[i])) next
    lab <- pe$lab[i]
    hw <- 0.15 + 0.038 * nchar(lab); hh <- 0.22
    Pr <- resample_path(paths[[i]], 140L)
    d <- sqrt(rowSums(diff(Pr)^2)); s <- c(0, cumsum(d))
    s <- if (s[length(s)] > 1e-9) s / s[length(s)] else seq(0, 1, length.out = nrow(Pr))
    at <- function(u) {
      xi <- stats::approx(s, Pr[, 1], u)$y; yi <- stats::approx(s, Pr[, 2], u)$y
      du <- 0.02
      x1 <- stats::approx(s, Pr[, 1], max(0, u - du))$y; y1 <- stats::approx(s, Pr[, 2], max(0, u - du))$y
      x2 <- stats::approx(s, Pr[, 1], min(1, u + du))$y; y2 <- stats::approx(s, Pr[, 2], min(1, u + du))$y
      tx <- x2 - x1; ty <- y2 - y1; Lt <- sqrt(tx^2 + ty^2)
      c(xi, yi, if (Lt > 1e-9) -ty / Lt else 0, if (Lt > 1e-9) tx / Lt else 0)
    }
    best <- NULL; best_s <- -Inf
    for (u in seq(0.16, 0.84, by = 0.035)) {
      q <- at(u)
      for (o in c(-0.55, -0.30, 0, 0.30, 0.55)) {
        cx <- q[1] + q[3] * o; cy <- q[2] + q[4] * o
        g <- gap_to(cx, cy, hw, hh, blocks)
        if (nrow(placed)) g <- min(g, gap_to(cx, cy, hw, hh, placed))
        pen <- min(cx - xs_lim[1] - hw, xs_lim[2] - cx - hw,
                   cy - y_lim[1] - hh, y_lim[2] - cy - hh)
        if (pen < 0.20) g <- min(g, pen)
        sc <- g - 0.55 * abs(o) - 0.45 * abs(u - 0.5)
        if (sc > best_s + 1e-9) { best_s <- sc; best <- c(cx, cy) }
      }
    }
    pe$lx[i] <- best[1]; pe$ly[i] <- best[2]
    placed <- rbind(placed, c(best[1], best[2], hw, hh))
  }

  pd <- do.call(rbind, lapply(seq_len(nrow(pe)), function(i) {
    P <- paths[[i]]
    data.frame(id = i, x = P[, 1], y = P[, 2], col = pe$col[i], lw = pe$lw[i],
               lty = pe$lty[i], sig = pe$sig[i],
               stringsAsFactors = FALSE)
  }))
  pd_sig <- pd[pd$sig, , drop = FALSE]
  pd_ns  <- pd[!pd$sig, , drop = FALSE]

  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(limits = xs_lim, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_family = CFG$figure$font %||% "Arial") +
    ggplot2::theme(aspect.ratio = 1 / panel_wh)
  # 不显著路径：灰色点线、无箭头（不主张方向性）
  if (nrow(pd_ns))
    p <- p + ggplot2::geom_path(data = pd_ns,
                                ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                             linetype = .data$lty),
                                colour = col_ns, linewidth = 0.45)
  # 显著路径：实线、蓝=正 / 红=负、线宽 ∝ |beta|、带箭头
  if (nrow(pd_sig))
    p <- p + ggplot2::geom_path(data = pd_sig,
                                ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                             colour = .data$col, linewidth = .data$lw,
                                             linetype = .data$lty),
                                arrow = ggplot2::arrow(length = grid::unit(1.4, "mm"),
                                                       type = "closed")) +
      ggplot2::scale_colour_manual(values = c(col_pos, col_neg), guide = "none") +
      ggplot2::scale_linewidth_identity()
  # 线型图例（与 04_Fig04_SEM.R 的图例口径一致；放在所有层之后以免被 scale 重置）
  p <- p + ggplot2::scale_linetype_manual(
    name = NULL,
    values = c("Significant paths" = "solid", "Non-significant paths" = "dotted"),
    guide = ggplot2::guide_legend(override.aes = list(linewidth = 0.6, colour = "grey20")))
  e_sig <- pe[pe$sig & !is.na(pe$lab), , drop = FALSE]
  if (nrow(e_sig))
    p <- p + ggplot2::geom_label(data = e_sig, ggplot2::aes(x = .data$lx, y = .data$ly,
                                                            label = .data$lab),
                                 size = pt2gg(MD_PT$ann * 0.9), colour = e_sig$col,
                                 fill = "white", linewidth = 0,
                                 label.padding = grid::unit(0.05, "lines"),
                                 fontface = "plain")
  p <- p + ggplot2::geom_label(data = nd, ggplot2::aes(x = .data$x, y = .data$y,
                                                       label = .data$label),
                               size = pt2gg(MD_PT$ann), lineheight = 0.9, fill = "white",
                               colour = "black", label.padding = grid::unit(0.22, "lines"),
                               label.r = grid::unit(0.16, "lines"), linewidth = 0.35) +
    ggplot2::annotate("text", x = fc_pt[1], y = fc_pt[2], label = fc_str,
                      size = pt2gg(MD_PT$ann * 0.9), colour = "grey30") +
    ggplot2::theme(legend.position = "bottom")

  md_save_figure(CFG, p, "52_sem", height = 150)

  # ---------------------------------------------------------------------------
  # 5. 表
  # ---------------------------------------------------------------------------
  co_out <- data.frame(
    Response = path_tab$Response, Predictor = path_tab$Predictor,
    Estimate = suppressWarnings(as.numeric(path_tab$Estimate)),
    Std_Estimate = suppressWarnings(as.numeric(path_tab$Std.Estimate)),
    Std_Error = suppressWarnings(as.numeric(path_tab$Std.Error %||% NA_real_)),
    DF = suppressWarnings(as.numeric(path_tab$DF %||% NA_real_)),
    Crit_Value = suppressWarnings(as.numeric(path_tab$Crit.Value %||% NA_real_)),
    P = suppressWarnings(as.numeric(path_tab[[pcol]])),
    P_BH = path_tab$P_BH, sig = path_tab$sig,
    layer = edges$layer[match(paste(path_tab$Predictor, path_tab$Response),
                              paste(edges$from, edges$to))],
    stringsAsFactors = FALSE)
  if (!"Std_Error" %in% names(co_out) || all(is.na(co_out$Std_Error)))
    co_out$Std_Error <- suppressWarnings(as.numeric(path_tab$Std.Error))

  fit_out <- data.frame(
    Item = c("n_samples", "n_paths", "n_endogenous", "estimators", "standardized",
             "dsep_tested", "Fisher_C", "Fisher_C_P", "dsep_df", "Verdict",
             "AIC_total", "aliased_terms", "on_missing_policy"),
    Value = c(n_use, nrow(path_tab), length(models),
              paste(unique(r2_tab$estimator), collapse = "; "), std,
              dsep, ifelse(is.na(fc), NA, fc), ifelse(is.na(fc_p), NA, fc_p),
              if (!is.null(fc_tab) && "df" %in% names(fc_tab)) fc_tab$df[1] else NA,
              if (is.na(fc_p)) "no d-separation test (saturated / disabled)"
              else if (fc_p > 0.05) "fits (not rejected)" else "REJECTED by d-separation",
              tryCatch(as.numeric(sum(vapply(models, stats::AIC, numeric(1)))),
                       error = function(e) NA_real_),
              if (length(aliased)) paste(unique(aliased), collapse = ", ") else "none",
              on_miss),
    stringsAsFactors = FALSE)
  dsep_out <- if (!is.null(fc_tab) && nrow(fc_tab)) fc_tab else
    data.frame(Note = "未做 d-separation（CFG$sem$dsep = FALSE 或模型饱和）")
  if (!is.null(sm$dTable) && nrow(sm$dTable)) dsep_out <- as.data.frame(sm$dTable)

  note_txt <- c(
    paste0("结构由用户声明：", nrow(edges), " 条边 / ", length(layers), " 层（",
           paste(layers, collapse = " -> "), "）"),
    paste0("全变量 z 化 = ", std, "；路径系数取 Std.Estimate，可作为标准化路径系数"),
    paste0("P 校正：", p_adj %||% "none", "（P_BH 列）"),
    if (is.na(fc_p)) "Fisher's C 不可用（模型饱和或未做 d-separation）"
    else if (fc_p > 0.05) "Fisher's C 不拒绝该 DAG（P > 0.05）"
    else "Fisher's C 拒绝该 DAG（P < 0.05）：声明的因果结构不足以重现数据的条件独立性，结论须谨慎",
    if (length(aliased)) paste0("模型不可识别：存在共线/别名项 ", paste(unique(aliased), collapse = ", "))
    else "模型可识别（无别名项）",
    if (any(grepl("singular", unlist(fit_note)))) "存在 singular fit：随机效应方差 ~ 0，随机结构信息量有限" else "",
    "路径图按 CFG$sem$layers 分层；跨层边若不穿过其它节点框走直线，否则自动贝塞尔绕行",
    "正系数 #2166AC / 负系数 #B2182B；仅显著路径标系数，不显著路径为灰色点线")
  note_txt <- note_txt[nzchar(note_txt)]

  sheets <- list(
    Path_coefficients = co_out,
    Fit = r2_tab,
    Fit_Interpretation = fit_out,
    dSeparation = dsep_out,
    Variables = do.call(rbind, prov),
    Note = data.frame(Note = note_txt, stringsAsFactors = FALSE))
  sheets <- Filter(function(x) !is.null(x) && nrow(x) > 0, sheets)
  tp <- md_write_xlsx(CFG, "52_sem", sheets)
  md_log(CFG, "OK", "52: ", nrow(path_tab), " 条路径 | Fisher's C = ",
         ifelse(is.na(fc), "NA", sprintf("%.2f (P = %.3f)", fc, fc_p)),
         " | n = ", n_use)
  md_result(module, "OK", tables = basename(tp), figures = "figures/52_sem.pdf",
            notes = if (is.na(fc_p)) NULL else if (fc_p > 0.05) "d-sep 不拒绝" else
              "d-sep 拒绝该 DAG（见表注）")
}
