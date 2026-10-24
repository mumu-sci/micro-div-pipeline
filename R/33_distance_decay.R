# =============================================================================
# R/33_distance_decay.R —— 距离衰减关系（Distance-Decay Relationship, DDR）
# -----------------------------------------------------------------------------
# 口径来源（原始研究 原脚本，只读参考，不修改）：
#   · 地理距离   = Haversine（km）—— 本包统一走 md_geo_dist()
#   · 群落相似度 = 1 - Bray-Curtis（来自 md_dist_get(CFG, "bray")，即模块 13）
#   · 拟合模型   = lm(log(Similarity) ~ log(Dist_km))，报斜率 + P + R²
#   · 画图       = 距离衰减散点 + 拟合线，log-log 双轴
# 本包相对 原实现扩展：
#   · 「同组 vs 跨组」按 CFG$design$factors 的第一个因子 + 第二个因子各自判定
#     （nifH 只按 Vegetation 分）；两组曲线同图对比。
#   · 三口径都出：地理距离 / 环境距离 / 群落距离 vs 相似度。
#   · 样本对最少对数门限参数化（CFG$distance_decay$min_pairs，默认 10）。
#
# 层级与伪重复（示范数据关键事实）：
#   11 个站点的经纬度是**站点级**的（站内 4 个重复共享同一坐标），因此站内样本对
#   的地理距离恒为 0。这些对全部被 Dist_km > 0 过滤掉，剩下的是**跨站点对**，
#   故本模块同时给出「站点对」聚合口径（站对间 Bray 取均值），作为稳健性对照。
#
# 输出：figures/33_distance_decay.pdf（一张 patchwork，4 面板）
#       results/33_distance_decay.xlsx
#       cache/33_ddr.rds（距离矩阵单独存 rds，绝不并进 data.frame）
#
# 依赖：13（Bray 距离）。若 13 的缓存不存在，回退到本模块自算 Bray（见下方
#       local_bray_from_core），并在 Note sheet 与日志里写明回退口径。
# =============================================================================

# ---------------------------------------------------------------------------
# local_* 辅助（共享库缺的写在自己文件里）
# ---------------------------------------------------------------------------

#' log-log 距离衰减拟合：log(Similarity) ~ log(Dist)。
#' 优雅处理：对数 < min_pairs、自变量零方差、非有限值、拟合失败 → 返回 NA 行。
local_ddr_fit <- function(df, x = "Dist", y = "Similarity",
                          min_pairs = 10L, label = "", metric = NA_character_,
                          group_var = "", group_level = "", level = "sample_pair") {
  out <- data.frame(
    metric = metric, level = level, group_var = group_var, group_level = group_level,
    label = label, n_pairs = 0L,
    slope = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
    r2 = NA_real_, adj_r2 = NA_real_, p = NA_real_, intercept = NA_real_,
    mean_dist = NA_real_, median_dist = NA_real_,
    mean_sim = NA_real_, median_sim = NA_real_,
    note = NA_character_, stringsAsFactors = FALSE)

  if (is.null(df) || nrow(df) == 0 ||
      !all(c(x, y) %in% names(df))) {
    out$note <- "无输入数据"
    return(out)
  }
  d <- df[is.finite(df[[x]]) & is.finite(df[[y]]) & df[[x]] > 0 & df[[y]] > 0,
          c(x, y), drop = FALSE]
  names(d) <- c("x", "y")
  out$n_pairs     <- nrow(d)
  out$mean_dist   <- if (nrow(d)) mean(d$x) else NA_real_
  out$median_dist <- if (nrow(d)) stats::median(d$x) else NA_real_
  out$mean_sim    <- if (nrow(d)) mean(d$y) else NA_real_
  out$median_sim  <- if (nrow(d)) stats::median(d$y) else NA_real_

  if (nrow(d) < min_pairs) {
    out$note <- sprintf("对数 %d < min_pairs %d，不拟合", nrow(d), min_pairs)
    return(out)
  }
  d$lx <- log(d$x); d$ly <- log(d$y)
  d <- d[is.finite(d$lx) & is.finite(d$ly), , drop = FALSE]
  if (nrow(d) < min_pairs || length(unique(d$lx)) < 2) {
    out$note <- "自变量 log(Dist) 零方差或有效对数不足"
    return(out)
  }
  fit <- tryCatch(stats::lm(ly ~ lx, data = d), error = function(e) NULL)
  if (is.null(fit)) { out$note <- "lm 拟合失败"; return(out) }
  sm <- summary(fit)
  cf <- sm$coefficients
  if (nrow(cf) < 2 || !is.finite(cf[2, 1])) { out$note <- "斜率不可估"; return(out) }
  ci <- tryCatch(stats::confint(fit)[2, ], error = function(e) c(NA_real_, NA_real_))
  out$slope     <- unname(cf[2, 1])
  out$se        <- unname(cf[2, 2])
  out$ci_lo     <- unname(ci[1]); out$ci_hi <- unname(ci[2])
  out$intercept <- unname(cf[1, 1])
  out$r2        <- unname(sm$r.squared)
  out$adj_r2    <- unname(sm$adj.r.squared)
  out$p         <- unname(cf[2, 4])
  out$note      <- ""
  out
}

#' 地理距离矩阵（优先共享 md_geo_dist()；缺列/缺包返回 NULL）
local_geo <- function(env, ids, CFG = NULL) {
  d <- tryCatch(md_geo_dist(env, ids), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  m <- as.matrix(d)
  if (nrow(m) != length(ids) || any(!is.finite(m))) return(NULL)
  m
}

#' 环境距离矩阵：对标准化后的环境变量取欧氏距离（stats::dist）。
#' 常量列 / 非数值列自动剔除；有效变量 < 2 返回 NULL。
local_env_dist <- function(df, ids, vars) {
  vars <- intersect(vars, names(df))
  vars <- vars[vapply(df[vars], is.numeric, logical(1))]
  if (length(vars) < 2) return(NULL)
  X <- as.data.frame(df[match(ids, df[[1]]), vars, drop = FALSE])
  ok <- vapply(X, function(v) is.finite(stats::sd(v, na.rm = TRUE)) &&
                 stats::sd(v, na.rm = TRUE) > 0, logical(1))
  vars <- vars[ok]
  if (length(vars) < 2) return(NULL)
  X <- X[, vars, drop = FALSE]
  # 缺行整行剔除（宁可少几个样本，也不用中位数填补 —— 补出来的是假数据）
  cc <- stats::complete.cases(X)
  if (sum(cc) < 3) return(NULL)
  Xc <- X[cc, , drop = FALSE]
  Z  <- scale(as.matrix(Xc))
  if (!nrow(Z)) return(NULL)
  rownames(Z) <- ids[cc]
  m <- as.matrix(stats::dist(Z))
  dimnames(m) <- list(ids[cc], ids[cc])
  attr(m, "n_used") <- sum(cc)
  attr(m, "vars") <- vars
  m
}

#' 样本对长表（上三角），携带指定分组列
local_pairs <- function(mats, meta, ids, group_vars) {
  k <- length(ids)
  ut <- which(upper.tri(matrix(0, k, k)), arr.ind = TRUE)
  d <- data.frame(sample1 = ids[ut[, 1]], sample2 = ids[ut[, 2]],
                  stringsAsFactors = FALSE)
  for (nm in names(mats)) d[[nm]] <- mats[[nm]][ut]
  for (g in group_vars) {
    d[[paste0(g, "1")]] <- as.character(meta[[g]][ut[, 1]])
    d[[paste0(g, "2")]] <- as.character(meta[[g]][ut[, 2]])
    d[[paste0("Same_", g)]] <- d[[paste0(g, "1")]] == d[[paste0(g, "2")]]
  }
  d
}

#' 站点对聚合：把样本对按「站点-站点」分组取均值（伪重复对照口径）
local_site_pairs <- function(d, meta, site_col) {
  s1 <- as.character(meta[[site_col]])[match(d$sample1, meta[[1]])]
  s2 <- as.character(meta[[site_col]])[match(d$sample2, meta[[1]])]
  key <- ifelse(s1 < s2, paste(s1, s2, sep = "|"), paste(s2, s1, sep = "|"))
  num_cols <- names(d)[vapply(d, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, "Dist_km")   # 地理距离取组内均值即可，与其余数值列同样处理
  sp <- split(seq_len(nrow(d)), key)
  rows <- lapply(sp, function(i) {
    row <- as.data.frame(lapply(d[i, num_cols, drop = FALSE], function(x) mean(x, na.rm = TRUE)),
                         stringsAsFactors = FALSE)
    row$site1 <- s1[i[1]]; row$site2 <- s2[i[1]]
    row$n_sample_pairs <- length(i)
    row
  })
  do.call(rbind, rows)
}

#' 反查某个样本对属于哪个「站点对」（供聚合列回填，仅用于图注）
local_bray_from_core <- function(CFG) {
  m <- tryCatch(md_core_asv(CFG, "distances"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  m <- m[rowSums(m) > 0, , drop = FALSE]
  if (ncol(m) < 3 || nrow(m) < 2) return(NULL)
  tryCatch(vegan::vegdist(t(m), method = "bray"), error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------
md_distance_decay <- function(CFG) {
  module <- "33"
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少依赖包"))
  if (!md_precheck(CFG, module)) return(md_result(module, "SKIP", "前置条件不满足"))

  dd_cfg   <- CFG$distance_decay %||% list()
  min_pairs <- as.integer(dd_cfg$min_pairs %||% 10L)
  env_vars <- (dd_cfg$env_vars %||% CFG$env$pca_vars %||% CFG$env$vars)
  group_vars <- (dd_cfg$group_vars %||% CFG$design$factors) %||% character(0)
  group_vars <- group_vars[group_vars %in% names(md_core_meta(CFG))]

  # ---- 1. 数据 ----
  ana  <- tryCatch(md_core_analysis(CFG), error = function(e) NULL)
  if (is.null(ana) || nrow(ana) < 4) {
    md_skip(CFG, module, "分析总表不可用或样本过少")
    return(md_result(module, "SKIP", "分析总表不可用"))
  }
  smp  <- CFG$data$sample_col
  ids  <- as.character(ana[[smp]])

  # ---- 2. 群落距离（优先模块 13；缺则回退自算，并登记 Note）----
  notes <- character(0)
  bray <- tryCatch(md_dist_get(CFG, "bray"), error = function(e) NULL)
  bray_src <- "模块 13（cache/13_dist.rds, bray）"
  if (is.null(bray)) {
    bray <- local_bray_from_core(CFG)
    bray_src <- "回退：本模块用 vegan::vegdist 自算（模块 13 缓存缺失）"
    notes <- c(notes, "模块 13 的 Bray 距离不可用，本模块回退为对 md_core_asv(CFG,'distances') 自算 Bray-Curtis。")
    md_warn(CFG, module, " 独立：", bray_src)
  }
  if (is.null(bray)) {
    md_skip(CFG, module, "Bray 距离不可用（模块 13 未跑且无法自算）")
    return(md_result(module, "SKIP", "Bray 距离不可用"))
  }
  bm <- as.matrix(bray)
  if (is.null(rownames(bm))) { md_skip(CFG, module, "Bray 矩阵无样本名"); return(md_result(module, "SKIP", "Bray 无样本名")) }
  common <- intersect(ids, rownames(bm))
  if (length(common) < 4) {
    md_skip(CFG, module, "Bray 矩阵与分析总表交集过少")
    return(md_result(module, "SKIP", "样本交集过少"))
  }
  ana  <- ana[match(common, ana[[smp]]), , drop = FALSE]
  ids  <- common
  bm   <- bm[ids, ids, drop = FALSE]

  # ---- 3. 地理距离 ----
  geo <- local_geo(ana, ids)
  has_geo <- !is.null(geo)

  # ---- 4. 环境距离 ----
  edist <- local_env_dist(ana, ids, env_vars)
  if (is.null(edist)) notes <- c(notes, "环境距离不可用（有效环境变量 < 2）。")

  # ---- 5. 站点级判定 ----
  site_col <- CFG$design$site_var
  lat_lon  <- c(dd_cfg$lat_col %||% "Latitude", dd_cfg$lon_col %||% "Longitude")
  site_level <- FALSE
  if (!is.null(site_col) && site_col %in% names(ana)) {
    uq <- vapply(intersect(lat_lon, names(ana)), function(v)
      md_within_unique(ana, v, site_col), integer(1))
    site_level <- length(uq) > 0 && all(uq == 1)
  }
  if (site_level) {
    notes <- c(notes, sprintf(
      "坐标在 %s 内唯一（%s 为站点级变量）→ 站内样本对地理距离恒为 0；样本级口径下这些对已被 Dist>0 过滤，实际比较的是跨站点对。同时给出站点对聚合口径作稳健性对照。",
      site_col, paste(intersect(lat_lon, names(ana)), collapse = "/")))
  }

  # ---- 6. 样本对长表 ----
  # 环境距离可能因缺失值只覆盖部分样本 → 统一把 ids 收缩到「地理 × 群落 × 环境」的公共集合，
  # 保证所有矩阵同维（local_pairs 才能安全索引）。
  if (!is.null(edist)) {
    common2 <- intersect(ids, rownames(edist))
    if (has_geo) common2 <- intersect(common2, rownames(geo))
    if (length(common2) < length(ids)) {
      notes <- c(notes, sprintf(
        "环境变量存在缺失 → 统一样本集收缩到 %d/%d；环境距离整行剔除（不插补）。",
        length(common2), length(ids)))
    }
    ids <- common2
    ana <- ana[match(ids, ana[[smp]]), , drop = FALSE]
    bm  <- bm[ids, ids, drop = FALSE]
    if (has_geo) geo <- geo[ids, ids, drop = FALSE]
    edist <- edist[ids, ids, drop = FALSE]
  }
  if (length(ids) < 4) {
    md_skip(CFG, module, "公共样本集过少")
    return(md_result(module, "SKIP", "公共样本集过少"))
  }
  # ★ local_pairs 的 4 个位置参数：Geo / Env / Bray / Bray(未用)。
  #   列名必须与下面引用的一致：Geo_dist_km / Env_dist / Bray_dist。
  na_mat <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  pair <- local_pairs(list(
    Geo_dist_km = if (has_geo) geo[ids, ids, drop = FALSE] else na_mat,
    Env_dist    = if (!is.null(edist)) edist[ids, ids, drop = FALSE] else na_mat,
    Bray_dist   = bm[ids, ids, drop = FALSE],
    Bray_unused = bm[ids, ids, drop = FALSE]), ana, ids, group_vars)
  pair$Similarity <- 1 - pair$Bray_dist
  pair$Dist_km <- pair$Geo_dist_km
  # 环境距离可能因缺失值只覆盖部分样本（local_env_dist 整行剔除）→ 此处应为全覆盖
  if (!is.null(pair$Env_dist)) {
    n_env <- sum(is.finite(pair$Env_dist))
    if (n_env < nrow(pair)) {
      notes <- c(notes, "环境距离矩阵仍有 NA（不应出现，请检查 local_env_dist）。")
    }
  }

  n_pairs_all <- nrow(pair)

  # ---- 7. 三口径拟合 ----
  fits <- list()

  # 7a 地理 vs 相似度
  if (has_geo) {
    d0 <- pair[pair$Dist_km > 0, , drop = FALSE]
    fits[["Overall"]] <- local_ddr_fit(d0, "Dist_km", "Similarity", min_pairs,
                                       "Overall", metric = "Geographic")
    if (length(group_vars)) {
      gf <- do.call(rbind, lapply(group_vars, function(g) {
        lv <- unique(c(d0[[paste0(g, "1")]], d0[[paste0(g, "2")]]))
        lv <- lv[!is.na(lv)]
        same <- rbind(
          local_ddr_fit(d0[d0[[paste0("Same_", g)]] %in% TRUE, , drop = FALSE],
                        "Dist_km", "Similarity", min_pairs,
                        "Same", metric = "Geographic", group_var = g, group_level = "Same"),
          local_ddr_fit(d0[d0[[paste0("Same_", g)]] %in% FALSE, , drop = FALSE],
                        "Dist_km", "Similarity", min_pairs,
                        "Across", metric = "Geographic", group_var = g, group_level = "Across"))
        per <- do.call(rbind, lapply(lv, function(x) {
          ii <- d0[[paste0(g, "1")]] == x & d0[[paste0(g, "2")]] == x
          local_ddr_fit(d0[ii, , drop = FALSE], "Dist_km", "Similarity", min_pairs,
                        x, metric = "Geographic", group_var = g, group_level = x)
        }))
        rbind(same, per)
      }))
      fits[["By_group"]] <- gf
    }
  } else {
    notes <- c(notes, "分析总表缺经纬度列（或 geosphere 缺失）→ 地理距离口径跳过。")
  }

  # 7b 环境 vs 相似度
  if (!is.null(edist)) {
    fits[["Env"]] <- local_ddr_fit(pair[pair$Env_dist > 0, , drop = FALSE],
                                   "Env_dist", "Similarity", min_pairs,
                                   "Overall", metric = "Environmental")
  }
  # 7c 群落距离 vs 相似度（Bray 距离 vs 1-Bray，恒等于 1-x，作退化诊断口径）
  fits[["Community"]] <- local_ddr_fit(pair[pair$Bray_dist > 0, , drop = FALSE],
                                       "Bray_dist", "Similarity", min_pairs,
                                       "Overall", metric = "Community")

  # ---- 8. 站点对聚合口径（伪重复对照）----
  site_fits <- NULL
  if (has_geo && site_level && !is.null(site_col) && site_col %in% names(ana)) {
    sp <- tryCatch(local_site_pairs(pair, ana, site_col), error = function(e) NULL)
    if (!is.null(sp) && nrow(sp) >= 3) {
      sp$Similarity <- 1 - sp$Bray_dist
      site_fits <- local_ddr_fit(sp, "Geo_dist_km", "Similarity", min_pairs,
                                 "Overall (site-pair mean)", metric = "Geographic",
                                 level = "site_pair")
      site_fits$n_sites <- length(unique(c(sp$site1, sp$site2)))
    }
  }

  fit_tab <- do.call(rbind, fits)
  if (is.list(site_fits) && !is.data.frame(site_fits)) site_fits <- as.data.frame(site_fits)
  if (!is.null(site_fits)) {
    # 站点对口径带一个专列 n_sites；样本对口径补 NA 以保证 rbind 列数一致
    fit_tab$n_sites <- NA_integer_
    if (is.null(site_fits$n_sites)) site_fits$n_sites <- NA_integer_
    fit_tab <- rbind(fit_tab, site_fits)
  }
  rownames(fit_tab) <- NULL
  fit_tab$sig <- ifelse(is.na(fit_tab$p), "—",
                 ifelse(fit_tab$p < 0.001, "***",
                 ifelse(fit_tab$p < 0.01, "**",
                 ifelse(fit_tab$p < 0.05, "*", "ns"))))

  # ---- 9. 描述统计 ----
  desc_tab <- data.frame(
    metric = c("Geographic (km)", "Environmental (z-Euclidean)", "Community (Bray)")[c(has_geo, !is.null(edist), TRUE)],
    n_pairs = c(if (has_geo) sum(pair$Dist_km > 0) else NA_integer_,
                if (!is.null(edist)) sum(pair$Env_dist > 0) else NA_integer_,
                sum(pair$Bray_dist > 0)),
    min = c(if (has_geo) min(pair$Dist_km, na.rm = TRUE) else NA_real_,
            if (!is.null(edist)) min(pair$Env_dist, na.rm = TRUE) else NA_real_,
            min(pair$Bray_dist, na.rm = TRUE)),
    median = c(if (has_geo) stats::median(pair$Dist_km, na.rm = TRUE) else NA_real_,
               if (!is.null(edist)) stats::median(pair$Env_dist, na.rm = TRUE) else NA_real_,
               stats::median(pair$Bray_dist, na.rm = TRUE)),
    mean = c(if (has_geo) mean(pair$Dist_km, na.rm = TRUE) else NA_real_,
             if (!is.null(edist)) mean(pair$Env_dist, na.rm = TRUE) else NA_real_,
             mean(pair$Bray_dist, na.rm = TRUE)),
    max = c(if (has_geo) max(pair$Dist_km, na.rm = TRUE) else NA_real_,
            if (!is.null(edist)) max(pair$Env_dist, na.rm = TRUE) else NA_real_,
            max(pair$Bray_dist, na.rm = TRUE)),
    stringsAsFactors = FALSE)
  desc_tab$Sim_mean <- 1 - desc_tab$mean
  desc_tab$n_pairs_ge0_sim <- sum(pair$Similarity > 0)
  desc_tab$pct_pairs_sim0 <- 100 * mean(pair$Similarity <= 0)


  # ---- 10. 出图 ----
  # 画法逐项对照 原实现 附图（make_ddr_panel）：
  #   散点 size = PT_SIZE_STD(MD_PT_SIZE$point)、alpha = PT_ALPHA_STD(MD_PT_SIZE$point_alpha)、
  #   颜色 = 该分组的 cols_veg / cols_zone（走 md_pal()）；
  #   拟合线 geom_smooth(method="lm", se=TRUE) 用 REG_LINE_COL/REG_LINE_LW/REG_RIBBON_A(MD_REG)；
  #   面板内右上角标注 "R2 = ..., slope = ..." + P（hjust=1, vjust=1.3, lineheight=0.9, grey30）；
  #   双对数轴 scale_x_log10(labels=scales::comma) + scale_y_log10()；
  #   轴题照抄原图英文："Geographic Distance (km)" / "Bray-Curtis Similarity"。
  pt <- md_pt(CFG)
  col_grp <- md_pal(CFG, c("Same", "Across"))

  # in-panel 的 P 文案：照抄 nifH 原脚本（ASCII；"NA" 不写 "—"）
  fmt_p_fig <- function(p) ifelse(is.na(p), "NA",
                           ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))

  # 面板内右上角标注：照抄 附图 的 annotate 位置与字号
  ann_layer <- function(r) {
    lab <- if (is.null(r) || is.na(r$slope))
      sprintf("n = %d\nnot fitted", r$n_pairs %||% 0L)
    else
      sprintf("R2 = %.3f, slope = %.4f\n%s", r$r2, r$slope, fmt_p_fig(r$p))
    ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1.3,
                      label = lab, size = pt2gg(pt$ann), colour = "grey30",
                      lineheight = 0.9)
  }

  # A：总体地理距离衰减（log-log）
  pA <- md_empty_plot(CFG, "No geographic coordinates")
  if (has_geo) {
    d0 <- pair[pair$Dist_km > 0 & pair$Similarity > 0, , drop = FALSE]
    pA <- ggplot2::ggplot(d0, ggplot2::aes(x = .data$Dist_km, y = .data$Similarity)) +
      ggplot2::geom_point(size = MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha,
                           colour = MD_REG$col) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                           linewidth = MD_REG$lw, colour = MD_REG$col,
                           fill = MD_REG$col, alpha = MD_REG$ribbon_alpha) +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::scale_y_log10() +
      ann_layer(fits[["Overall"]]) +
      ggplot2::labs(x = "Geographic Distance (km)", y = "Bray-Curtis Similarity", tag = "A") +
      theme_md(CFG = CFG)
  }

  # B：同组 vs 跨组（按第一个设计因子）
  pB <- md_empty_plot(CFG, "No factor for grouping")
  if (has_geo && length(group_vars)) {
    g1 <- group_vars[1]
    d0 <- pair[pair$Dist_km > 0 & pair$Similarity > 0, , drop = FALSE]
    d0$SameGroup <- factor(ifelse(d0[[paste0("Same_", g1)]], "Same", "Across"),
                           levels = c("Same", "Across"))
    f_same <- fits[["By_group"]]
    f_same <- f_same[f_same$group_var == g1 & f_same$group_level %in% c("Same", "Across"), , drop = FALSE]
    pB <- ggplot2::ggplot(d0, ggplot2::aes(x = .data$Dist_km, y = .data$Similarity,
                                           colour = .data$SameGroup)) +
      ggplot2::geom_point(size = MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                           linewidth = MD_REG$lw, fill = MD_REG$col,
                           alpha = MD_REG$ribbon_alpha) +
      ggplot2::scale_colour_manual(values = col_grp, name = g1) +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Geographic Distance (km)", y = "Bray-Curtis Similarity", tag = "B") +
      theme_md(CFG = CFG) +
      ggplot2::theme(legend.position = "bottom")
    # 每条拟合线各自的 R2 / slope / P：锚在面板内部两侧（与 附图 同一标注格式）
    for (i in seq_len(nrow(f_same))) {
      lvl <- as.character(f_same$group_level[i])
      pB <- pB + ggplot2::annotate("text",
        x = if (identical(lvl, "Same")) -Inf else Inf, y = Inf,
        hjust = if (identical(lvl, "Same")) -0.08 else 1.08, vjust = 1.3,
        label = sprintf("R2 = %.3f, slope = %.4f\n%s",
                        f_same$r2[i], f_same$slope[i], fmt_p_fig(f_same$p[i])),
        size = pt2gg(pt$ann), colour = unname(col_grp[lvl]), lineheight = 0.9)
    }
  }

  # C：按第一个因子的各水平单独拟合（对应 附图 的「按植被型分组」2x2 版式）
  pC <- md_empty_plot(CFG, "No per-level fit")
  if (has_geo && length(group_vars)) {
    g1 <- group_vars[1]
    d0 <- pair[pair$Dist_km > 0 & pair$Similarity > 0, , drop = FALSE]
    d0[[g1]] <- d0[[paste0(g1, "1")]]
    f_lv <- fits[["By_group"]]
    lv <- unique(f_lv$group_level[f_lv$group_var == g1 &
                                  !f_lv$group_level %in% c("Same", "Across")])
    if (length(lv)) {
      d0 <- d0[d0[[g1]] %in% lv, , drop = FALSE]
      d0[[g1]] <- factor(as.character(d0[[g1]]), levels = lv)
      pC <- ggplot2::ggplot(d0, ggplot2::aes(x = .data$Dist_km, y = .data$Similarity,
                                             colour = .data[[g1]])) +
        ggplot2::geom_point(size = MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                             linewidth = MD_REG$lw, fill = MD_REG$col,
                             alpha = MD_REG$ribbon_alpha) +
        ggplot2::scale_colour_manual(values = md_pal(CFG, lv, g1)) +
        ggplot2::scale_x_log10(labels = scales::comma) +
        ggplot2::scale_y_log10() +
        ggplot2::facet_wrap(stats::as.formula(paste0("~ `", g1, "`")), nrow = 1) +
        ggplot2::labs(x = "Geographic Distance (km)", y = "Bray-Curtis Similarity", tag = "C") +
        theme_md(CFG = CFG) +
        ggplot2::theme(legend.position = "none")
    }
  }

  # D（斜率森林图）**原脚本没有** —— nifH 原实现 附图 是
  # (Bare|Herb)/(Shrub|Tree) 的 2×2 DDR 四格，没有任何森林图。故从图面移除；
  # 各口径的斜率/CI/P 仍完整写在 DDR_fits 表里。

  fig <- md_compose(CFG, list(A = pA, B = pB, C = pC),
                    design = "AB\nCC", guides = "keep", legend = "bottom")
  md_save_figure(CFG, fig, "33_distance_decay", height = 150)

  # ---- 11. 写表 + 缓存 ----
  note_df <- data.frame(Note = c(
    sprintf("地理距离：Haversine（km），来源 md_geo_dist()；群落距离来源：%s", bray_src),
    sprintf("拟合模型：lm(log(Similarity) ~ log(Geographic distance))，Similarity = 1 - Bray-Curtis。"),
    sprintf("样本对最少对数门限 min_pairs = %d（CFG$distance_decay$min_pairs）。", min_pairs),
    sprintf("样本对数（全部）= %d；站内样本对地理距离为 0 已被过滤。", n_pairs_all),
    if (site_level) sprintf("站点级变量：%s（坐标站内唯一）→ 站内样本对 Dist = 0；已额外给出站点对聚合口径。", site_col) else "未检出站点级坐标。",
    sprintf("环境距离：CFG$env$pca_vars 交集共 %d 个变量，z-score 后取欧氏距离（stats::dist）。", length(env_vars)),
    "群落距离口径（Bray vs 1-Bray）在数学上恒为 y = 1 - x，仅作退化诊断，不解释为生态学信号。",
    notes))

  tp <- md_write_xlsx(CFG, "33_distance_decay", list(
    DDR_fits = fit_tab,
    Distance_summary = desc_tab,
    Note = note_df),
    meta_extra = list(n_pairs = n_pairs_all, min_pairs = min_pairs,
                      bray_source = bray_src, site_level = site_level))

  saveRDS(list(pair = if (isTRUE(CFG$table$detail)) pair else NULL,
               fits = fit_tab,
               dist = list(geo = if (has_geo) stats::as.dist(geo) else NULL,
                           env = if (is.null(edist)) NULL else stats::as.dist(edist),
                           bray = as.dist(bm)),
               params = list(min_pairs = min_pairs, env_vars = env_vars,
                             group_vars = group_vars, site_level = site_level,
                             bray_source = bray_src, seed = CFG$opts$seed)),
          md_path(CFG, "cache", "33_ddr.rds"))

  md_result(module, "OK",
            tables = basename(tp),
            figures = "figures/33_distance_decay.pdf",
            notes = paste(notes, collapse = " | "))
}
