# =============================================================================
# R/60_spatial_autocorr.R ---- 空间自相关：Moran's I + correlogram（合并为一张图）
# -----------------------------------------------------------------------------
# 契约见 _build-notes/micro-div-pipeline/MODULE_SPEC.md
#
# 做什么
#   面板 A  Moran's I 条形图（k-NN 反距离权重，置换检验，显著性星号）
#   面板 B  correlogram（距离分段二值权重；x 轴 = 等距【类别】轴，不是数值轴）
#   两者都在 CFG$spatial$level 指定的分析层级上做，并对原始值 + 环境残差各做一遍。
#
# 分析层级（CFG$spatial$level）
#   auto   ---- 坐标有重复（站内重复）则**强制站点级聚合**；坐标全唯一才做样本级
#   sample ---- 强制样本级；若坐标有重复，如实登记「站内距离恒为 0 -> 不可检验」并降级
#   site   ---- 强制站点级
#
# 统计口径
#   - Moran's I：k-NN 反距离权重，k = CFG$spatial$moran_k 或 max(3, min(n/4, 10))
#   - 置换检验：CFG$spatial$permutations（双侧，以 E[I] 为中心）
#   - correlogram：分段【二值】权重（段内为 1），置换同上
#   - 残差 = 变量对 CFG$model$covariates（+ CFG$model$derived）的 lm 残差
#
# 参数来源（全部 CFG，未给的走 %||% 默认）
#   CFG$spatial$level / moran_k / permutations / dist_classes / min_links
#   CFG$design$site_var, CFG$response$var, CFG$model$covariates, CFG$model$derived,
#   CFG$response_models$alpha_focus, CFG$response_models$ordination_axes
#
# 依赖：05（cache/05_core.rds）；残差协变量若含派生变量，则优先 40，缺则本地兜底。
# =============================================================================

md_spatial_autocorr <- function(CFG) {
  module <- "60"

  if (!md_need(CFG, c("vegan", "geosphere"), module))
    return(md_result(module, "SKIP", "缺少包 vegan / geosphere"))
  if (!md_precheck(CFG, module))
    return(md_result(module, "SKIP", "前置条件不满足"))

  d_raw <- tryCatch(md_analysis_with_derived(CFG), error = function(e) NULL)
  if (is.null(d_raw)) {
    md_skip(CFG, module, "取不到分析总表（cache/05_core.rds 缺失？）")
    return(md_result(module, "SKIP", "取不到分析总表"))
  }
  core <- md_core(CFG)
  smp  <- CFG$data$sample_col
  site_var <- CFG$design$site_var
  lat <- "Latitude"; lon <- "Longitude"

  if (!all(c(lat, lon) %in% names(d_raw))) {
    md_skip(CFG, module, paste0("数据里没有经纬度列（", lat, " / ", lon, "）"))
    return(md_result(module, "SKIP", "缺经纬度列"))
  }

  # ---------------------------------------------------------------------------
  # 1. 决定分析层级
  # ---------------------------------------------------------------------------
  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  ok_ll <- !is.na(d_raw[[lat]]) & !is.na(d_raw[[lon]])
  n_coord_unique <- length(unique(paste(round(d_raw[[lat]][ok_ll], 6),
                                        round(d_raw[[lon]][ok_ll], 6))))
  dup_coord <- n_coord_unique < sum(ok_ll)
  lvl_req <- (CFG$spatial$level %||% "auto")
  level <- switch(lvl_req,
                  site   = "site",
                  sample = "sample",
                  if (dup_coord) "site" else "sample")

  if (identical(lvl_req, "auto") && dup_coord) {
    add_note("spatial$level = auto：坐标存在重复（", sum(ok_ll), " 个样本仅 ",
             n_coord_unique, " 个唯一坐标），按站点级分析。")
    md_log(CFG, "INFO", "[60] auto + 坐标重复 -> 站点级聚合")
  }
  if (identical(level, "sample")) {
    # 样本级：站内距离恒为 0 时不可检验，如实登记并降级为站点级
    within0 <- FALSE
    if (dup_coord) within0 <- TRUE
    if (within0) {
      add_note("样本级空间自相关因站内距离恒为 0 而不可检验（坐标是站点级，",
               sum(ok_ll), " 个样本仅 ", n_coord_unique, " 个唯一坐标）；",
               "已自动改用站点级聚合，样本级 Moran's I 记为 NA。")
      md_warn(CFG, "[60] 样本级空间自相关不可检验（站内距离恒为 0），降级为站点级")
      level <- "site"
    }
  }
  if (identical(level, "site")) {
    if (is.null(site_var) || !site_var %in% names(d_raw)) {
      md_skip(CFG, module, paste0("需要站点级分析，但 design$site_var（",
                                  site_var %||% "NULL", "）不在数据中"))
      return(md_result(module, "SKIP", "站点级分析缺 site_var"))
    }
  }

  # ---------------------------------------------------------------------------
  # 2. 变量清单（响应 / α / 排序轴）
  # ---------------------------------------------------------------------------
  # 社区指标（Shannon / PD / 排序轴）不在 md_analysis_with_derived() 的总表里：
  # 本模块按下述顺序自行补齐，缺什么就补什么 ----
  #   Shannon / Observed / Simpson / Pielou -> 本地由 md_core_asv("alpha") 现算
  #   PD                                     -> 21_phylo_diversity 的缓存
  #   PCoA*                                  -> 13_distances 的 bray 距离现算
  local_alpha <- function(cfg, ids, idx) {
    m <- tryCatch(md_core_asv(cfg, "alpha"), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    mm <- t(m)[ids, , drop = FALSE]
    out <- data.frame(sample_id = as.character(ids), stringsAsFactors = FALSE)
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
      if (!is.null(v) && length(v) == nrow(out)) out[[k]] <- as.numeric(v)
    }
    if (ncol(out) < 2) return(NULL)
    out
  }
  local_pcoa <- function(cfg, ids, metric = "bray", axes = 2L) {
    dmat <- tryCatch(md_dist_get(cfg, metric), error = function(e) NULL)
    if (is.null(dmat)) {
      m <- tryCatch(md_core_asv(cfg, "distances"), error = function(e) NULL)
      if (is.null(m)) return(NULL)
      mm <- t(m)[ids, , drop = FALSE]
      dmat <- tryCatch(vegan::vegdist(mm, method = metric), error = function(e) NULL)
    }
    if (is.null(dmat)) return(NULL)
    M <- as.matrix(dmat)
    keep <- intersect(ids, rownames(M))
    if (length(keep) < 3) return(NULL)
    p <- tryCatch(ape::pcoa(stats::as.dist(M[keep, keep, drop = FALSE]),
                            correction = "cailliez"), error = function(e) NULL)
    if (is.null(p)) return(NULL)
    eig <- p$values$Corr_eig
    if (is.null(eig) || all(is.na(eig))) eig <- p$values$Eigenvalues
    pct <- 100 * eig[seq_len(axes)] / sum(eig[eig > 0])
    out <- as.data.frame(p$vectors[, seq_len(axes), drop = FALSE])
    names(out) <- paste0("PCoA", seq_len(axes))
    out$sample_id <- rownames(p$vectors)
    attr(out, "pct") <- pct
    out
  }

  yv <- CFG$response$var
  want_all <- unique(c(yv,
                       CFG$response_models$alpha_focus %||% character(0),
                       CFG$response_models$ordination_axes %||% character(0),
                       intersect(c("Shannon", "PD", "Observed", "Pielou", "PCoA1"), names(d_raw))))
  cand <- want_all[!is.na(want_all) & nzchar(want_all) & want_all %in% names(d_raw)]

  miss <- setdiff(want_all[!is.na(want_all) & nzchar(want_all)], names(d_raw))
  miss_alpha <- intersect(miss, c("Shannon", "Simpson", "Observed", "Richness", "Pielou"))
  if (length(miss_alpha)) {
    a <- local_alpha(CFG, as.character(d_raw[[smp]]), miss_alpha)
    if (!is.null(a)) {
      newc <- intersect(setdiff(names(a), "sample_id"), miss_alpha)
      ix <- match(as.character(d_raw[[smp]]), a$sample_id)
      for (cl in newc) d_raw[[cl]] <- a[[cl]][ix]
      add_note("α 指标 ", paste(newc, collapse = ", "),
               " 不在总表里，已由 md_core_asv(\"alpha\") 本地现算（",
               sum(is.finite(d_raw[[newc[1]]])), "/", nrow(d_raw), " 个样本有值）。")
      cand <- unique(c(cand, newc))
    }
  }
  miss_pd <- intersect(miss, c("PD", "MPD", "MNTD", "NRI", "NTI"))
  if (length(miss_pd)) {
    ph <- tryCatch(md_phylo_cache(CFG), error = function(e) NULL)
    if (is.data.frame(ph) && smp %in% names(ph)) {
      ex <- intersect(names(ph), miss_pd)
      if (length(ex)) {
        ix <- match(as.character(d_raw[[smp]]), as.character(ph[[smp]]))
        for (cl in ex) d_raw[[cl]] <- ph[[cl]][ix]
        add_note("系统发育指标 ", paste(ex, collapse = ", "),
                 " 已并入（来源：21_phylo_diversity 缓存；",
                 sum(is.finite(d_raw[[ex[1]]])), "/", nrow(d_raw), " 个样本有值）。")
        cand <- unique(c(cand, ex))
      }
    } else {
      add_note("取不到 21 号产物（md_phylo_cache），",
               paste(miss_pd, collapse = ", "), " 不可用。")
    }
  }
  miss_pc <- intersect(miss, c("PCoA1", "PCoA2", "PCoA3"))
  if (length(miss_pc)) {
    pc <- local_pcoa(CFG, as.character(d_raw[[smp]]),
                     CFG$response_models$pair_metric %||% "bray", 3L)
    if (!is.null(pc)) {
      ex <- intersect(names(pc)[names(pc) != smp], miss_pc)
      if (length(ex)) {
        ix <- match(as.character(d_raw[[smp]]), as.character(pc[[smp]]))
        for (cl in ex) d_raw[[cl]] <- pc[[cl]][ix]
        note_pc <- if (!is.null(attr(pc, "pct"))) sprintf("%.1f%%", attr(pc, "pct")[1]) else "NA"
        add_note("排序轴 ", paste(ex, collapse = ", "), " 本地重算（轴 1 解释 ",
                 note_pc, "）。")
        cand <- unique(c(cand, ex))
      }
    }
  }

  vars <- cand[vapply(cand, function(v) v %in% names(d_raw) && is.numeric(d_raw[[v]]) &&
                        sum(is.finite(d_raw[[v]])) >= 5 &&
                        stats::sd(d_raw[[v]], na.rm = TRUE) > 0, logical(1))]
  if (length(vars) < length(cand)) {
    add_note("以下变量在数据中不可用或零方差，已剔除：",
             paste(setdiff(cand, vars), collapse = ", "), "。")
  }
  if (length(vars) == 0) {
    md_skip(CFG, module, "没有任何可用于空间自相关的数值变量")
    return(md_result(module, "SKIP", "无可用变量"))
  }

  # ---------------------------------------------------------------------------
  # 3. 组装分析表 + 坐标 + 站点聚合
  # ---------------------------------------------------------------------------
  # 残差协变量（CFG$model$covariates + CFG$model$derived）必须一起进分析表，
  # 否则站点聚合后无处可取 ---- 曾导致"协变量集为空"而残差全部记 NA。
  cf_covs <- unique(c(CFG$model$covariates %||% character(0),
                      CFG$model$derived %||% character(0)))
  cf_covs <- cf_covs[cf_covs %in% names(d_raw)]
  keep_cols <- unique(c(smp, site_var, lat, lon, vars, cf_covs))
  d <- d_raw[, keep_cols, drop = FALSE]

  if (identical(level, "site")) {
    # md_site_aggregate() 只保留 site_var + vars，故协变量必须一并列进 vars
    agg <- md_site_aggregate(d, site_var, vars = unique(c(vars, cf_covs)), id_col = smp)
    # 站点坐标取该站第一行（坐标本身是站点级常量）
    cxy <- do.call(rbind, lapply(split(seq_len(nrow(d)), d[[site_var]]), function(ix) {
      data.frame(site = as.character(d[[site_var]][ix[1]]),
                 lat = d[[lat]][ix[1]], lon = d[[lon]][ix[1]],
                 stringsAsFactors = FALSE)
    }))
    names(cxy) <- c(site_var, lat, lon)
    agg <- merge(agg, cxy, by = site_var, all.x = TRUE, sort = FALSE)
    agg <- agg[match(cxy[[site_var]], agg[[site_var]]), , drop = FALSE]
    d_use <- agg
    unit  <- "site"
    unit_id <- site_var
  } else {
    d_use <- d
    unit  <- "sample"
    unit_id <- smp
  }

  coords <- as.matrix(d_use[, c(lon, lat), drop = FALSE])
  coords <- apply(coords, 2, as.numeric)
  keep <- stats::complete.cases(coords)
  if (sum(keep) < 5) {
    md_skip(CFG, module, paste0("有效坐标点仅 ", sum(keep), " 个（少于 5），无法做空间自相关"))
    return(md_result(module, "SKIP", "有效坐标点不足"))
  }
  d_use <- d_use[keep, , drop = FALSE]
  coords <- coords[keep, , drop = FALSE]
  n_unit <- nrow(d_use)

  dm_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
  dimnames(dm_km) <- list(d_use[[unit_id]], d_use[[unit_id]])
  upper_km <- dm_km[upper.tri(dm_km)]
  add_note("分析层级 = ", unit, "；n = ", n_unit, "；站点/样本间地理距离 ",
           sprintf("%.2f", min(upper_km)), " - ", sprintf("%.1f", max(upper_km)),
           " km（中位 ", sprintf("%.1f", stats::median(upper_km)), " km）。")

  # ---------------------------------------------------------------------------
  # 4. Moran's I（k-NN 反距离权重）
  # ---------------------------------------------------------------------------
  k_nn <- CFG$spatial$moran_k
  if (is.null(k_nn)) k_nn <- max(3, min(floor(n_unit / 4), 10))
  B <- as.integer(CFG$spatial$permutations %||% 999)

  local_moran_knn <- function(v, label) {
    okv <- is.finite(v)
    if (sum(okv) < 5) return(NULL)
    vv <- v[okv]
    n  <- length(vv)
    k  <- max(1L, min(as.integer(k_nn), n - 1L))
    dd <- dm_km[okv, okv, drop = FALSE]
    w  <- matrix(0, n, n)
    for (i in seq_len(n)) {
      ord <- order(dd[i, ])
      nn  <- ord[ord != i][seq_len(k)]
      dnn <- dd[i, nn]
      if (any(dnn == 0)) dnn[dnn == 0] <- min(dnn[dnn > 0], na.rm = TRUE) / 2
      w[i, nn] <- 1 / dnn
    }
    rs <- rowSums(w); rs[rs == 0] <- 1
    w  <- w / rs
    z  <- vv - mean(vv)
    S0 <- sum(w)
    if (!is.finite(S0) || S0 <= 0 || sum(z^2) <= 0) return(NULL)
    moran_stat <- function(zz) (n / S0) * sum(w * outer(zz, zz)) / sum(zz^2)
    I  <- moran_stat(z)
    EI <- -1 / (n - 1)
    md_seed(CFG, 601)
    Isim <- replicate(B, moran_stat(sample(z)))
    p <- (sum(abs(Isim - EI) >= abs(I - EI)) + 1) / (B + 1)
    data.frame(Variable = label, n = n, Moran_I = I, Expected_I = EI,
               P_perm = p, k_nn = k, level = unit,
               stringsAsFactors = FALSE)
  }

  raw_tab <- do.call(rbind, lapply(vars, function(v) local_moran_knn(d_use[[v]], v)))

  # ---- 环境残差 ----
  covs <- cf_covs[cf_covs %in% names(d_use)]
  covs <- covs[vapply(d_use[covs], is.numeric, logical(1))]
  covs <- covs[vapply(d_use[covs], function(x) sum(is.finite(x)) >= 5, logical(1))]
  if (length(covs) == 0) {
    add_note("协变量集为空（model$covariates / model$derived 均不在数据中或不可用），",
             "环境残差的 Moran's I 记 NA。")
  }

  resid_cols <- character(0)
  resid_tab <- NULL
  if (length(covs) > 0) {
    for (v in vars) {
      ddv <- data.frame(y = d_use[[v]], d_use[, covs, drop = FALSE])
      okv <- stats::complete.cases(ddv)
      n_ok <- sum(okv); p_cov <- length(covs)
      if (n_ok < max(6, p_cov + 2)) next
      if (stats::sd(ddv$y[okv]) <= 0) next
      r <- tryCatch(stats::residuals(stats::lm(
        stats::as.formula(paste0("y ~ ", paste0("`", covs, "`", collapse = " + "))),
        data = ddv)), error = function(e) NULL)
      if (is.null(r) || stats::sd(r, na.rm = TRUE) <= 0) next
      col <- paste0(v, "_resid")
      d_use[[col]] <- NA_real_
      d_use[[col]][okv] <- as.numeric(r)
      resid_cols <- c(resid_cols, col)
    }
    if (length(resid_cols) == 0)
      add_note("所有变量都无法拟合协变量 lm 残差（样本量不足 / 零方差），",
               "环境残差的 Moran's I 记 NA。")
    if (length(resid_cols) > 0) {
      resid_tab <- do.call(rbind, lapply(resid_cols, function(cl) {
        g <- local_moran_knn(d_use[[cl]], paste0(sub("_resid$", "", cl), " (env. resid.)"))
        if (is.null(g)) return(NULL)
        g
      }))
    }
  }

  moran_all <- rbind(raw_tab, resid_tab)
  if (is.null(moran_all) || nrow(moran_all) == 0) {
    md_skip(CFG, module, "所有变量的 Moran's I 都无法计算（样本量不足）")
    return(md_result(module, "SKIP", "Moran's I 无法计算"))
  }
  moran_all$Value_type <- ifelse(grepl("\\(env. resid\\.\\)$", moran_all$Variable),
                                 "Env. residual", "Raw")
  add_note("Moran's I 权重 = k 近邻反距离（k = ", k_nn, "）；置换次数 = ", B, "。")
  add_note("环境残差 = 变量对 ", paste(covs, collapse = " + "),
           " 的 lm 残差。")

  # ---------------------------------------------------------------------------
  # 5. correlogram（距离分段二值权重）
  # ---------------------------------------------------------------------------
  brks <- CFG$spatial$dist_classes %||% c(0, 10, 100, 600, 1600)
  brks <- sort(unique(as.numeric(brks)))
  if (length(brks) < 2) brks <- c(0, 10, 100, 600, 1600)
  classes <- lapply(seq_len(length(brks) - 1), function(i) c(brks[i], brks[i + 1]))
  min_links <- as.integer(CFG$spatial$min_links %||% 3)

  n_links_tot <- vapply(classes, function(cl) sum(upper_km > cl[1] & upper_km <= cl[2]),
                        integer(1))
  if (all(n_links_tot == 0)) {
    brks <- unique(stats::quantile(upper_km, probs = seq(0, 1, length.out = 5)))
    classes <- lapply(seq_len(length(brks) - 1), function(i) c(brks[i], brks[i + 1]))
    n_links_tot <- vapply(classes, function(cl) sum(upper_km > cl[1] & upper_km <= cl[2]),
                          integer(1))
    add_note("CFG$spatial$dist_classes 与实测距离范围不重叠，correlogram 已改用四分位分段。")
  }

  local_correlogram <- function(v, label) {
    okv <- is.finite(v)
    if (sum(okv) < 5) return(NULL)
    vv <- v[okv]; dd <- dm_km[okv, okv, drop = FALSE]
    n  <- length(vv); z <- vv - mean(vv)
    if (sum(z^2) <= 0) return(NULL)
    EI <- -1 / (n - 1)
    stat_of <- function(w, zz) {
      rs <- rowSums(w); rs[rs == 0] <- 1
      wn <- w / rs; S0 <- sum(wn)
      if (!is.finite(S0) || S0 <= 0) return(NA_real_)
      (n / S0) * sum(wn * outer(zz, zz)) / sum(zz^2)
    }
    out <- lapply(classes, function(cl) {
      w <- (dd > cl[1] & dd <= cl[2]) * 1
      diag(w) <- 0
      nl <- sum(w > 0)
      if (nl < min_links) return(NULL)
      I <- stat_of(w, z)
      if (!is.finite(I)) return(NULL)
      md_seed(CFG, 602)
      Isim <- replicate(B, stat_of(w, sample(z)))
      data.frame(Variable = label, lo = cl[1], hi = cl[2],
                 Dist_mid = mean(cl), Moran_I = I, Expected_I = EI,
                 P_perm = (sum(abs(Isim - EI) >= abs(I - EI)) + 1) / (B + 1),
                 n_links = nl, level = unit, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }

  cg <- do.call(rbind, lapply(vars, function(v) local_correlogram(d_use[[v]], v)))
  if (is.null(cg) || nrow(cg) == 0) {
    add_note("所有距离分段的连接数都少于 min_links = ", min_links,
             "，correlogram 无可用点（面板 B 记为空白）。")
  } else {
    cg$Class <- sprintf("%.0f-%.0f", cg$lo, cg$hi)
    lv <- unique(cg$Class[order(cg$Dist_mid)])
    cg$Class <- factor(cg$Class, levels = lv)
    add_note("correlogram 距离段（km）= ",
             paste(sprintf("%.0f-%.0f (n=%d)", cg$lo, cg$hi, cg$n_links),
                   collapse = ", "), "；分段内连接数 < ", min_links, " 的段不出点。")
  }

  # ---------------------------------------------------------------------------
  # 6. 出图
  # ---------------------------------------------------------------------------
  pt <- md_pt(CFG)

  # ---- 画法对齐 18_Supp_spatial_autocorr.R §FigS13 ----
  # A：单一色的 geom_col（原脚本用 cols_zone[["Xisha"]]），**没有** fill 图例；
  #    x 轴按 Moran's I 排序（原脚本 aes(x = reorder(Variable, Moran_I))）
  # B：correlogram，等距【类别】轴；颜色 = 内置 zone 色板按顺序取（原脚本把
  #    Shannon/Faith_PD/PCoA1/ARA 逐一手写成 cols_zone 的四个值，顺序完全一致）
  bar_col  <- unname(md_pal(CFG, "Xisha", "zone")["Xisha"])
  zone_cols <- unname(md_pal(CFG, names(MD_BUILTIN_PAL$zone), "zone"))

  ma <- moran_all
  ma$lab <- sprintf("%.3f%s", ma$Moran_I, ifelse(ma$P_perm < 0.05, "*", ""))
  ma$Var_f <- stats::reorder(ma$Variable, ma$Moran_I)

  pA <- ggplot2::ggplot(ma, ggplot2::aes(x = Var_f, y = Moran_I)) +
    ggplot2::geom_col(fill = bar_col, width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linewidth = MD_LW$ref) +
    ggplot2::geom_text(ggplot2::aes(label = lab), hjust = -0.15,
                       size = pt2gg(pt$ann)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.30))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = paste0("Moran's I (", unit, " level)"), tag = "A") +
    theme_md(CFG = CFG)

  if (!is.null(cg) && nrow(cg) > 0) {
    cg$Var_f <- factor(cg$Variable, levels = levels(ma$Var_f))
    cg_cols <- stats::setNames(rep_len(zone_cols, nlevels(cg$Var_f)), levels(cg$Var_f))
    pB <- ggplot2::ggplot(cg, ggplot2::aes(x = Class, y = Moran_I,
                                           colour = Var_f, group = Var_f)) +
      ggplot2::geom_hline(yintercept = unique(cg$Expected_I), linetype = "dashed",
                          linewidth = MD_LW$ref, colour = "grey50") +
      ggplot2::geom_line(linewidth = MD_LW$line) +
      ggplot2::geom_point(ggplot2::aes(shape = P_perm < 0.05),
                          size = MD_PT_SIZE$point * 2.4) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                                  name = "Permutation P", labels = c("n.s.", "< 0.05")) +
      ggplot2::scale_colour_manual(values = cg_cols, name = NULL) +
      ggplot2::labs(x = "Distance class (km)", y = "Moran's I", tag = "B") +
      theme_md(CFG = CFG) +
      ggplot2::theme(legend.position = "bottom",
                     legend.key.size = grid::unit(0.28, "cm"))
  } else {
    pB <- md_empty_plot(CFG, "correlogram: no distance class has enough links")
    pB <- pB + ggplot2::labs(tag = "B")
  }

  # 版式对齐 FIGURE_PANELS$FigS13：183 x 85 mm，1 行 2 列等宽，
  # legend = "none"（B 的图例在面板右侧，故组合层不覆盖 legend.position）
  fig <- md_compose(CFG, list(A = pA, B = pB), design = "AB",
                    guides = "keep", legend = NULL)
  md_save_figure(CFG, fig, "60_spatial_autocorr", height = 85)

  # ---------------------------------------------------------------------------
  # 7. 出表
  # ---------------------------------------------------------------------------
  coords_tab <- data.frame(unit = d_use[[unit_id]],
                           Longitude = coords[, 1], Latitude = coords[, 2],
                           stringsAsFactors = FALSE)
  names(coords_tab)[1] <- unit_id
  if (identical(level, "site")) {
    for (v in vars) coords_tab[[v]] <- d_use[[v]]
  }

  note_tab <- data.frame(Item = c("Analysis level (requested)", "Analysis level (used)",
                                  "n units", "Unique coordinates in data",
                                  "Samples in data", "k (k-NN inverse distance)",
                                  "Permutations", "Distance classes (km)",
                                  "min_links", "Residual covariates", "Notes"),
                         Value = c(lvl_req, unit, n_unit, n_coord_unique,
                                   sum(ok_ll), k_nn, B,
                                   paste(sprintf("%.0f-%.0f", brks[-length(brks)], brks[-1]),
                                         collapse = ", "),
                                   min_links,
                                   if (length(covs)) paste(covs, collapse = " + ") else "NA",
                                   paste(notes, collapse = " ")),
                         stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "60_spatial_autocorr", list(
    Moran_I = moran_all,
    Correlogram = if (!is.null(cg) && nrow(cg) > 0)
      cg[, setdiff(names(cg), "Var_f"), drop = FALSE] else
      data.frame(Note = "no distance class reached min_links"),
    Coordinates = coords_tab,
    Notes = note_tab
  ))

  md_result(module, "OK",
            tables = basename(tp), figures = "figures/60_spatial_autocorr.pdf",
            notes = notes)
}
