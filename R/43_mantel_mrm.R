# =============================================================================
# R/43_mantel_mrm.R ---- Mantel / partial Mantel / MRM（**双口径**）+ 饱和诊断
# -----------------------------------------------------------------------------
# 口径（**必须写进表里，口径错了整个 P 都错**）：
#   · Mantel / MRM 的置换发生在【对象级】：置换的是 n 个样本（口径 pair）或
#     n 个站点（口径 site）的标签，**不是**把每一对当独立观测。
#     原实现实测诊断：对级置换会把 P 系统压小约 18 倍。
#   · 双口径：
#       pair ---- 对象 = 样本（n = 样本数），观测 = 全部样本对；
#       site ---- 对象 = 站点（n = 站点数），观测 = 全部站点对；
#              站间群落距离 = 站内所有样本对的平均，响应距离 = |delta(站均响应)|。
#       site 口径是控制"同站点多对"伪重复的保守口径。
#   · 效应量统一用 md_boot_rho() 的【站点对整簇】bootstrap 95% CI（口径 pair）
#     或【站点整簇】（口径 site）。
#
# 饱和 / 膨胀诊断（nifH 原在 16 号脚本，本包并入本模块）：
#   群落距离类变量有数学上限 1.0。当均值 ~ 1 或大量对顶到 1.0 时，R2 的上限被
#   结构性压低 ---- 低 R2 应读作"β 方差近乎耗尽"，**不得**读作"无效应"。
#   诊断结果写进 Saturation sheet，并在 Note 里写明"该关系受数学上限约束，
#   不宜解读为'无效应'"。
#
# 参数（全部从 CFG 读，缺键用 %||% 默认）：
#   CFG$mantel$y            响应变量（样本级数值列），默认 CFG$response$var
#   CFG$mantel$xs           解释距离，默认 CFG$response_models$pair_metric
#   CFG$mantel$partial      偏相关控制变量
#   CFG$mantel$permutations 置换次数，默认 CFG$beta_tests$permutations
#   CFG$mantel$method       "spearman" / "pearson" / "kendall"
#   CFG$mantel$bootstrap_B  bootstrap 重复数
#   CFG$mantel$levels       c("pair","site")
#   xs / partial 的每一项可以是：13 号缓存里的距离度量名（bray/jaccard/…）、
#   "geo"（地理距离）或 analysis 表里的数值列（取 |delta| 作为距离）。
#
# 依赖：13（cache/13_dist.rds）---- 缺失时用 05 的矩阵本地重算并记进 Note。
# =============================================================================

md_mantel_mrm <- function(CFG) {
  module <- "43"
  if (!md_need(CFG, "ecology", module))
    return(md_result(module, "SKIP", "缺少 vegan/ape/picante/betapart"))
  if (!md_need(CFG, "mantel", module))
    return(md_result(module, "SKIP", "缺少 ecodist/geosphere"))
  if (!md_precheck(CFG, module, need_response = TRUE))
    return(md_result(module, "SKIP", "前置条件不满足"))

  # ---------------------------------------------------------------------------
  # 0. 参数与数据
  # ---------------------------------------------------------------------------
  smp      <- CFG$data$sample_col
  y_var    <- CFG$mantel$y %||% CFG$response$var
  y_label  <- CFG$response$label %||% y_var
  xs_spec  <- as.character(CFG$mantel$xs %||%
                             (CFG$response_models$pair_metric %||% CFG$distances$metrics[1]))
  xs_spec  <- unique(xs_spec[nzchar(xs_spec)])
  pt_spec  <- as.character(CFG$mantel$partial %||% character(0))
  pt_spec  <- unique(pt_spec[nzchar(pt_spec)])
  nperm    <- as.integer(CFG$mantel$permutations %||% CFG$beta_tests$permutations %||% 999)
  m_method <- CFG$mantel$method %||% "spearman"
  mrm_z    <- isTRUE(CFG$mantel$standardize %||% TRUE)
  B_boot   <- as.integer(CFG$mantel$bootstrap_B %||%
                           CFG$response_models$bootstrap_B %||% 1000)
  levels   <- as.character(CFG$mantel$levels %||% c("pair", "site"))
  detail   <- isTRUE(CFG$table$detail)

  d <- md_analysis_with_derived(CFG)
  if (is.null(d) || !nrow(d)) {
    md_skip(CFG, module, "analysis 表为空")
    return(md_result(module, "SKIP", "analysis 表为空"))
  }
  ids      <- as.character(d[[smp]])
  site_var <- CFG$design$site_var
  site_of  <- if (!is.null(site_var) && site_var %in% names(d))
    as.character(d[[site_var]]) else NULL
  if (is.null(site_of)) {
    levels <- "pair"
    md_warn(CFG, "design$site_var 缺失 ---- 只跑 pair 口径")
  }

  # ---- 局部工具（shared 里没有的，一律 local_ 前缀） ----
  local_z <- function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
  }
  local_num_col <- function(dat, v) {
    if (is.null(v) || !v %in% names(dat)) return(NULL)
    x <- suppressWarnings(as.numeric(dat[[v]]))
    if (all(is.na(x))) NULL else x
  }
  local_canon <- function(M) {
    tri <- which(lower.tri(M), arr.ind = TRUE)
    data.frame(i = tri[, 1], j = tri[, 2],
               row = rownames(M)[tri[, 1]], col = colnames(M)[tri[, 2]],
               stringsAsFactors = FALSE)
  }
  local_attach <- function(df, M, name, ids_use) {
    dimnames(M) <- list(ids_use, ids_use)
    df[[name]] <- M[cbind(df$i, df$j)]
    df
  }
  local_abs_mat <- function(x, ids_use) {
    M <- abs(outer(x, x, "-"))
    dimnames(M) <- list(ids_use, ids_use)
    M
  }
  local_geo_mat <- function(dat, ids_use) {
    if (!all(c("Latitude", "Longitude") %in% names(dat))) return(NULL)
    g <- md_geo_dist(data.frame(id = as.character(dat[[smp]]),
                                Latitude = dat$Latitude, Longitude = dat$Longitude),
                     ids_use)
    if (is.null(g)) return(NULL)
    as.matrix(g)
  }
  local_site_mat <- function(M, site_of_use, sites) {
    out <- matrix(0, length(sites), length(sites), dimnames = list(sites, sites))
    for (a in sites) for (b in sites) {
      if (a == b) next
      ia <- which(site_of_use == a); ib <- which(site_of_use == b)
      v <- M[ia, ib, drop = FALSE]
      out[a, b] <- mean(v, na.rm = TRUE)
    }
    out[!is.finite(out)] <- NA_real_
    out
  }
  local_pair_label <- function(row_ids, col_ids, loc_of) {
    l1 <- loc_of[row_ids]; l2 <- loc_of[col_ids]
    ifelse(l1 == l2, paste0("within:", l1),
           paste0(pmin(l1, l2), "__", pmax(l1, l2)))
  }

  # 距离度量名字表（用于判断某字符串是不是"距离度量"）
  VEGDIST_METHODS <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski",
                       "jaccard", "gower", "altgower", "morisita", "horn",
                       "mountford", "raup", "binomial", "chao", "cao", "mahalanobis")
  dc <- tryCatch(md_dist_cache(CFG), error = function(e) NULL)
  dc_names <- if (!is.null(dc) && !is.null(dc$dist)) names(dc$dist) else character(0)

  spec_kind <- function(v) {
    if (tolower(v) %in% c("geo", "geography", "geo_dist", "geo_distance")) return("geo")
    if (length(dc_names) && any(tolower(dc_names) == tolower(v))) return("metric")
    if (tolower(v) %in% VEGDIST_METHODS) return("metric")
    if (!is.null(local_num_col(d, v))) return("numeric")
    "none"
  }
  get_metric_mat <- function(v, ids_use) {
    M <- NULL; src <- NA_character_
    if (length(dc_names)) {
      hit <- dc_names[tolower(dc_names) == tolower(v)]
      if (!length(hit)) hit <- dc_names[grepl(v, dc_names, ignore.case = TRUE)]
      if (length(hit)) {
        M <- as.matrix(dc$dist[[hit[1]]])
        src <- paste0("13_dist cache (", hit[1], ")")
      }
    }
    if (is.null(M) && tolower(v) %in% VEGDIST_METHODS) {
      m <- md_core_asv(CFG, "distances")
      bin <- identical(tolower(v), "jaccard") && isTRUE(CFG$distances$binary_jaccard)
      M <- as.matrix(vegan::vegdist(t(m), method = tolower(v), binary = bin))
      src <- "本地重算（13 号缓存不可用 / 度量不在其中）"
    }
    if (is.null(M)) return(NULL)
    keep <- intersect(ids_use, rownames(M))
    if (length(keep) < 3) return(NULL)
    M <- M[keep, keep, drop = FALSE]
    attr(M, "md_src") <- src
    M
  }

  # ---------------------------------------------------------------------------
  # 1. 变量解析 + 公共对象集
  # ---------------------------------------------------------------------------
  ynum <- local_num_col(d, y_var)
  if (is.null(ynum)) {
    md_skip(CFG, module, paste0("响应变量 ", y_var, " 不是数值列"))
    return(md_result(module, "SKIP", paste0("响应变量 ", y_var, " 不是数值列")))
  }

  prov <- list()   # 变量来源（写进 Variables sheet）
  add_prov <- function(variable, role, source) {
    prov[[length(prov) + 1L]] <<- data.frame(
      variable = variable, role = role, source = source, stringsAsFactors = FALSE)
  }

  # 逐项判定；不认识的项登记后剔除（不报错）
  dropped_spec <- character(0)
  spec_info <- do.call(rbind, lapply(unique(c(xs_spec, pt_spec)), function(v) {
    k <- spec_kind(v)
    if (k == "none") dropped_spec <<- c(dropped_spec, v)
    data.frame(spec = v, kind = k,
               role = if (v %in% xs_spec) "x" else "partial",
               stringsAsFactors = FALSE)
  }))
  if (length(dropped_spec))
    md_warn(CFG, "43: 以下解释变量既不是距离度量也不是数值列，已剔除：",
            paste(dropped_spec, collapse = ", "))
  spec_info <- spec_info[spec_info$kind != "none", , drop = FALSE]

  # 公共对象集：各距离矩阵 rownames 的交集，再叠加数值列的非缺失
  ids_use <- ids
  mats_pair <- list()
  for (i in seq_len(nrow(spec_info))) {
    v <- spec_info$spec[i]
    if (spec_info$kind[i] == "metric") {
      M <- get_metric_mat(v, ids)
      if (is.null(M)) {
        md_warn(CFG, "43: 距离度量 ", v, " 不可用，已剔除")
        spec_info$kind[i] <- "none"
        next
      }
      mats_pair[[v]] <- M
      ids_use <- intersect(ids_use, rownames(M))
      add_prov(v, spec_info$role[i], attr(M, "md_src"))
    } else if (spec_info$kind[i] == "geo") {
      M <- local_geo_mat(d, ids)
      if (is.null(M)) {
        md_warn(CFG, "43: 地理距离不可用（缺经纬度 / geosphere），已剔除 ", v)
        spec_info$kind[i] <- "none"
        next
      }
      mats_pair[[v]] <- M
      ids_use <- intersect(ids_use, rownames(M))
      add_prov(v, spec_info$role[i], "geosphere distHaversine (km)")
    } else {
      add_prov(v, spec_info$role[i], "analysis 表数值列 -> 成对 |delta|")
    }
  }
  spec_info <- spec_info[spec_info$kind != "none", , drop = FALSE]
  add_prov(y_var, "y", "analysis 表数值列 -> 成对 |delta|")
  if (nrow(spec_info) == 0) {
    md_skip(CFG, module, "没有任何可用的解释距离变量")
    return(md_result(module, "SKIP", "没有任何可用的解释距离变量"))
  }

  # 完备样本（在对象层面剔除，保证成对表始终是完整三角）
  keep_cc <- rep(TRUE, length(ids_use))
  for (i in seq_len(nrow(spec_info))) {
    if (spec_info$kind[i] == "numeric") {
      v <- local_num_col(d, spec_info$spec[i])[match(ids_use, ids)]
      keep_cc <- keep_cc & !is.na(v)
    }
  }
  keep_cc <- keep_cc & !is.na(ynum[match(ids_use, ids)])
  n_drop_cc <- sum(!keep_cc)
  ids_use <- ids_use[keep_cc]
  if (n_drop_cc > 0)
    md_log(CFG, "INFO", "43: 因缺失剔除 ", n_drop_cc, " 个样本（对象），剩余 ", length(ids_use))
  if (length(ids_use) < 4) {
    md_skip(CFG, module, paste0("可用样本太少（n = ", length(ids_use), "）"))
    return(md_result(module, "SKIP", "可用样本太少"))
  }
  pos <- match(ids_use, ids)
  # 站点向量必须**按样本名命名**：成对表的 row/col 是样本名，靠名字取站点对标签
  site_of_use <- if (!is.null(site_of)) stats::setNames(site_of[pos], ids_use) else NULL

  # 把每个度量矩阵都裁到 ids_use
  for (v in names(mats_pair)) {
    M <- mats_pair[[v]]
    mats_pair[[v]] <- M[ids_use, ids_use, drop = FALSE]
  }

  # ---------------------------------------------------------------------------
  # 2. 逐口径构建成对表 + Mantel/MRM
  # ---------------------------------------------------------------------------
  build_level <- function(level) {
    if (level == "site") {
      sites <- unique(site_of_use)
      lid   <- sites
      lsite <- sites
    } else {
      lid   <- ids_use
      lsite <- site_of_use
    }
    # 响应距离矩阵
    yv <- ynum[pos]
    if (level == "site") yv <- as.numeric(tapply(yv, site_of_use, mean, na.rm = TRUE)[lid])
    Dy <- local_abs_mat(yv, lid)
    # 解释距离矩阵
    Dx <- list()
    for (v in names(mats_pair)) {
      M <- mats_pair[[v]]
      if (level == "site") M <- local_site_mat(M, site_of_use, lid)
      Dx[[v]] <- M
    }
    for (i in seq_len(nrow(spec_info))) {
      if (spec_info$kind[i] != "numeric") next
      v <- spec_info$spec[i]
      xv <- local_num_col(d, v)[pos]
      if (level == "site") xv <- as.numeric(tapply(xv, site_of_use, mean, na.rm = TRUE)[lid])
      Dx[[v]] <- local_abs_mat(xv, lid)
    }
    # 成对表（顺序 = which(lower.tri) 的列优先顺序，ecodist::MRM 依赖这个顺序）
    P <- local_canon(Dy)
    P$y <- Dy[cbind(P$i, P$j)]
    for (v in names(Dx)) P <- local_attach(P, Dx[[v]], make.names(v), lid)
    if (is.null(lsite)) {
      P$cluster <- rep("all", nrow(P))
    } else if (level == "site") {
      P$cluster <- as.character(P$row)
    } else {
      P$cluster <- local_pair_label(P$row, P$col, lsite)
      if (all(is.na(P$cluster))) P$cluster <- rep("all", nrow(P))
    }
    # 供 MRM 用的 z 化列（线性变换，不改变 R2/P，只让 Beta 可横向比较）
    P$y_z <- local_z(P$y)
    for (v in names(Dx)) {
      cz <- make.names(v)
      if (cz %in% names(P)) P[[paste0(cz, "_z")]] <- local_z(P[[cz]])
    }
    list(id = lid, dist = Dx, Dy = Dy, P = P, yv = yv)
  }

  lv_ok <- intersect(levels, c("pair", "site"))
  if (!length(lv_ok)) lv_ok <- "pair"
  L <- lapply(lv_ok, build_level)
  names(L) <- lv_ok

  x_cols_use <- intersect(make.names(spec_info$spec), names(L[[1]]$P))

  mant_rows <- list(); mrm_rows <- list(); es_rows <- list(); sat_rows <- list()
  add_row <- function(lst, df) { if (!is.null(df) && nrow(df)) lst[[length(lst) + 1L]] <- df; lst }

  for (lv in lv_ok) {
    Li <- L[[lv]]
    P  <- Li$P
    n_obj <- length(Li$id)
    n_pair <- nrow(P)
    perm_txt <- paste0("object-level (n = ", n_obj,
                       if (lv == "site") " sites)" else " samples)")
    for (v in names(Li$dist)) {
      cx <- make.names(v)
      if (!cx %in% names(P)) next
      lab <- paste0(v, " vs |delta ", y_var, "| [", lv, "]")
      m1 <- md_mantel(CFG, as.dist(Li$dist[[v]]), Li$Dy,
                      method = m_method, permutations = nperm, label = lab)
      if (!is.null(m1)) {
        m1$level <- lv; m1$x <- v; m1$y <- y_var; m1$partial <- "--"
        m1$n_objects <- n_obj; m1$n_pairs <- n_pair; m1$permutation_level <- perm_txt
        mant_rows <- add_row(mant_rows, m1)
      }
      # partial Mantel：每个偏变量一条
      for (z in setdiff(names(Li$dist), v)) {
        if (z %in% xs_spec && !(z %in% pt_spec)) next   # partial 只对 CFG$mantel$partial 里的变量做
        m2 <- md_mantel(CFG, as.dist(Li$dist[[v]]), Li$Dy,
                        method = m_method, permutations = nperm,
                        partial = as.dist(Li$dist[[z]]),
                        label = paste0(v, " vs |delta ", y_var, "| | ", z, " [", lv, "]"))
        if (!is.null(m2)) {
          m2$level <- lv; m2$x <- v; m2$y <- y_var; m2$partial <- z
          m2$n_objects <- n_obj; m2$n_pairs <- n_pair; m2$permutation_level <- perm_txt
          mant_rows <- add_row(mant_rows, m2)
        }
      }
      # MRM：单解释变量（z 化口径：Beta 可横向比较，R2/P 不变）
      ycol <- if (mrm_z) "y_z" else "y"
      xcol <- if (mrm_z) paste0(cx, "_z") else cx
      mm <- md_mrm(CFG, P, ycol, xcol, nperm = nperm,
                   label = paste0("|delta ", y_var, "|", if (mrm_z) "(z)" else "",
                                  " ~ ", v, if (mrm_z) "(z)" else "", " [", lv, "]"))
      if (!is.null(mm)) {
        mm$level <- lv; mm$x <- v
        mm$n_objects <- n_obj; mm$permutation_level <- perm_txt
        mrm_rows <- add_row(mrm_rows, mm)
      }
      # 效应量：站点（对）整簇 bootstrap
      es <- md_boot_rho(CFG, P, x = cx, y = "y", cluster = "cluster",
                        B = B_boot, label = paste0(v, " [", lv, "]"))
      if (!is.null(es)) {
        es$level <- lv; es$x <- v
        es$p_mantel <- if (!is.null(m1)) m1$p[1] else NA_real_
        es_rows <- add_row(es_rows, es)
      }
      # 饱和诊断
      M <- Li$dist[[v]]
      vv <- M[lower.tri(M)]
      bounded <- (max(vv, na.rm = TRUE) <= 1 + 1e-8 && min(vv, na.rm = TRUE) >= -1e-8)
      sat_rows <- add_row(sat_rows, data.frame(
        level = lv, variable = v, kind = if (bounded) "dissimilarity [0,1]" else "abs-difference",
        bounded = bounded, n_pairs = length(vv),
        mean = mean(vv, na.rm = TRUE), sd = stats::sd(vv, na.rm = TRUE),
        min = min(vv, na.rm = TRUE), max = max(vv, na.rm = TRUE),
        frac_at_upper = if (bounded) mean(vv >= 0.9999) else NA_real_,
        frac_near_zero = mean(vv <= 0.01, na.rm = TRUE),
        stringsAsFactors = FALSE))
    }
    # 响应距离本身的极值（供饱和/膨胀诊断对照）
    yvv <- Li$Dy[lower.tri(Li$Dy)]
    sat_rows <- add_row(sat_rows, data.frame(
      level = lv, variable = paste0("|d", y_var, "|"),
      kind = "abs-difference (response)", bounded = FALSE, n_pairs = length(yvv),
      mean = mean(yvv), sd = stats::sd(yvv), min = min(yvv), max = max(yvv),
      frac_at_upper = NA_real_, frac_near_zero = mean(yvv <= 0.01),
      stringsAsFactors = FALSE))
    # MRM：全部解释变量 + 偏变量一起入模
    all_x <- intersect(make.names(spec_info$spec), names(P))
    all_x <- if (mrm_z) paste0(all_x, "_z") else all_x
    all_x <- all_x[all_x %in% names(P)]
    if (length(all_x) >= 2) {
      mm <- md_mrm(CFG, P, if (mrm_z) "y_z" else "y", all_x, nperm = nperm,
                   label = paste0("|delta ", y_var, "|", if (mrm_z) "(z)" else "", " ~ ",
                                  paste(all_x, collapse = " + "), " [", lv, "]"))
      if (!is.null(mm)) {
        mm$level <- lv; mm$x <- paste(all_x, collapse = " + ")
        mm$n_objects <- n_obj; mm$permutation_level <- perm_txt
        mrm_rows <- add_row(mrm_rows, mm)
      }
    }
  }

  mant <- if (length(mant_rows)) do.call(rbind, mant_rows) else NULL
  mrm  <- if (length(mrm_rows)) do.call(rbind, mrm_rows) else NULL
  es   <- if (length(es_rows)) do.call(rbind, es_rows) else NULL
  sat  <- if (length(sat_rows)) do.call(rbind, sat_rows) else NULL
  if (!is.null(mrm)) mrm <- mrm[as.character(mrm$term) != "Int", , drop = FALSE]  # ecodist 截距名是 Int

  if (is.null(es) || !nrow(es)) {
    md_skip(CFG, module, "成对表不足，无法给出效应量（对数 < 5）")
    return(md_result(module, "SKIP", "成对表不足"))
  }

  # 饱和度结论
  sat$note <- ifelse(!is.na(sat$frac_at_upper) &
                       (sat$frac_at_upper >= 0.5 | sat$mean >= 0.9),
                     "受数学上限约束：多数对已顶到 1.0，R2 上限被结构性压低", "")
  sat_note <- any(nzchar(sat$note))
  sat_txt <- if (sat_note)
    "该关系受数学上限约束（群落距离有上界），不宜解读为'无效应'：低 R2 应读作距离方差近乎耗尽。"
  else "本次分析的距离变量未检出明显饱和。"

  # ---------------------------------------------------------------------------
  # 3. 距离矩阵单独存 rds（绝不并进 data.frame）
  # ---------------------------------------------------------------------------
  dist_path <- md_path(CFG, "cache", "43_mantel_mrm_dists.rds")
  saveRDS(list(levels = lapply(L, function(x) x$dist),
               response = lapply(L, function(x) x$Dy),
               ids = lapply(L, function(x) x$id),
               params = list(y = y_var, xs = xs_spec, partial = pt_spec,
                             nperm = nperm, method = m_method, seed = CFG$opts$seed)),
          dist_path)

  # ---------------------------------------------------------------------------
  # 4. 图：成对散点（hexbin）+ 效应量森林图
  #    面板 A = 论文 Fig.S12B「区域 |delta 响应| vs 群落距离」，
  #             画法逐项照抄 nifH 原脚本 §B/§C 的 mk_hex：
  #      geom_hex(bins = 55) + 共用色标（low "grey90" / high "#2B6AA0" / trans "log10"
  #      / name "pairs" / limits = c(1, CNT_MAX) / breaks c(1,10,100) / oob squish）
  #      + geom_smooth(method = "lm", se = FALSE, colour "#B2182B", linewidth REG_LINE_LW)
  #      + 面板内右上角多行标注（hjust = 1, vjust = 1.3, lineheight = 0.95, grey25）：
  #          MRM 的 beta 与 P、R2；Mantel 的 rho 与 P；饱和度诊断 mean / sd / frac >= 0.9999
  #      + 副标题给量纲/口径（S_ANN_COMP 字号、grey35、hjust = 0）。
  #      hexbin 缺失时 md_warn 并降级为普通散点（不报错、不静默）。
  #    面板 B = 效应量森林图，走基元 md_forest()
  #      （原实现标准化 beta 森林图版式；坐标轴方向由基元给出）。
  # ---------------------------------------------------------------------------
  prim_lv <- lv_ok[1]
  prim_x  <- names(L[[prim_lv]]$dist)[1]
  P_prim  <- L[[prim_lv]]$P
  cx_prim <- make.names(prim_x)
  pt <- md_pt(CFG)

  # P 文案照抄 nifH 原脚本（ASCII；不用 md_fmt_p 的 "—"）
  fmt_p_fig <- function(p) ifelse(is.na(p) | length(p) == 0, "NA",
                           ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))
  # 取「某个 x 在该口径下」的那一行
  pick_row <- function(df, key, val, lvl = NULL) {
    if (is.null(df) || !nrow(df) || !key %in% names(df)) return(NULL)
    d <- df[!is.na(df[[key]]) & as.character(df[[key]]) == val, , drop = FALSE]
    if (!is.null(lvl) && "level" %in% names(d))
      d <- d[!is.na(d$level) & as.character(d$level) == lvl, , drop = FALSE]
    if (!nrow(d)) NULL else d[1, , drop = FALSE]
  }

  # 面板 A 的图内标注：MRM / Mantel / 饱和诊断三行（照抄 mk_lab_pairs() 的写法）
  a_mrm  <- pick_row(mrm,  "x",        prim_x, prim_lv)
  a_mant <- pick_row(mant, "x",        prim_x, prim_lv)
  a_sat  <- pick_row(sat,  "variable", prim_x, prim_lv)
  lab_lines <- character(0)
  if (!is.null(a_mrm))
    lab_lines <- c(lab_lines,
                   sprintf("MRM (linear): beta = %.2f, %s", a_mrm$estimate, fmt_p_fig(a_mrm$p)),
                   sprintf("R2 = %.4f", a_mrm$R2))
  if (!is.null(a_mant))
    lab_lines <- c(lab_lines,
                   sprintf("Mantel (rank): rho = %.3f, %s", a_mant$statistic, fmt_p_fig(a_mant$p)))
  if (!is.null(a_sat))
    lab_lines <- c(lab_lines,
                   sprintf("%s: mean %.3f, sd %.3f, %.0f%% >= 0.9999",
                           if (isTRUE(a_sat$bounded)) "beta saturated" else "beta not saturated",
                           a_sat$mean, a_sat$sd, 100 * a_sat$frac_at_upper))
  lab_a <- if (length(lab_lines)) paste(lab_lines, collapse = "\n") else
    sprintf("n = %d pairs", nrow(P_prim))
  sub_a <- sprintf("%s vs |delta %s|, %s level (%d objects, %d pairs)",
                   prim_x, y_var, prim_lv, length(L[[prim_lv]]$id), nrow(P_prim))

  if (requireNamespace("hexbin", quietly = TRUE)) {
    hx <- suppressWarnings(as.numeric(P_prim[[cx_prim]]))
    hy <- suppressWarnings(as.numeric(P_prim$y))
    ok <- is.finite(hx) & is.finite(hy)
    cnt_max <- if (any(ok)) max(hexbin::hexbin(hx[ok], hy[ok], xbins = 55)@count) else NA_real_
    pA <- ggplot2::ggplot(P_prim, ggplot2::aes(x = .data[[cx_prim]], y = .data$y)) +
      ggplot2::geom_hex(bins = 55) +
      ggplot2::scale_fill_gradient(low = "grey90", high = "#2B6AA0", trans = "log10",
                                   name = "pairs",
                                   limits = c(1, max(cnt_max, 2, na.rm = TRUE)),
                                   breaks = c(1, 10, 100), oob = scales::squish) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                           colour = "#B2182B", linewidth = MD_REG$lw)
  } else {
    # hexbin 缺失：显式告警后降级为普通散点（不报错、不静默）
    md_warn(CFG, "43: 缺少 hexbin 包，Fig.S12B 的六边形密度图降级为普通散点")
    pA <- ggplot2::ggplot(P_prim, ggplot2::aes(x = .data[[cx_prim]], y = .data$y)) +
      ggplot2::geom_point(size = MD_PT_SIZE$point * 0.6, alpha = MD_PT_SIZE$point_alpha,
                          colour = "grey30") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                           colour = "#B2182B", linewidth = MD_REG$lw)
  }
  pA <- pA +
    ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1.3,
                      label = lab_a, size = pt2gg(pt$ann), colour = "grey25",
                      lineheight = 0.95) +
    ggplot2::labs(x = prim_x, y = paste0("|delta ", y_var, "|"),
                  subtitle = sub_a, tag = "A") +
    theme_md(CFG = CFG) +
    ggplot2::theme(legend.position = "bottom",
                   plot.subtitle = ggplot2::element_text(size = MD_PT_COMP$ann,
                                                         colour = "grey35", hjust = 0))

  # B（效应量森林图）**原脚本已明确移除** —— 原脚本 写着
  # 效应量面板"→ 从 附图 移除，数字完整保留在 TableS15 的 B_effect_sizes 工作表"。
  # 故这里只出 A（hexbin |Δ响应| vs 群落距离 + MRM/Mantel/饱和度标注）；
  # 效应量数字仍完整写在 Effect_sizes 表里。

  md_save_figure(CFG, md_compose(CFG, list(A = pA), design = "A"),
                 "43_mantel_mrm", height = 105)

  # ---------------------------------------------------------------------------
  # 5. 表
  # ---------------------------------------------------------------------------
  mant_out <- if (!is.null(mant)) mant[, c("label", "level", "x", "y", "partial", "method",
                                           "statistic", "p", "permutations",
                                           "n_objects", "n_pairs", "permutation_level")]
              else data.frame(label = "--", level = "--", x = "--", y = y_var, partial = "--",
                              method = m_method, statistic = NA_real_, p = NA_real_,
                              permutations = nperm, n_objects = NA_integer_,
                              n_pairs = NA_integer_, permutation_level = "--")
  names(mant_out)[names(mant_out) == "statistic"] <- "r"
  names(mant_out)[names(mant_out) == "p"] <- "P"

  mrm_out <- if (!is.null(mrm) && nrow(mrm))
    mrm[, c("label", "level", "term", "estimate", "p", "R2", "n_perm", "n_obs",
            "n_objects", "permutation_level")]
  else data.frame(label = "--", level = "--", term = "--", estimate = NA_real_,
                  p = NA_real_, R2 = NA_real_, n_perm = nperm, n_obs = NA_integer_,
                  n_objects = NA_integer_, permutation_level = "--")
  names(mrm_out)[names(mrm_out) == "estimate"] <- "Beta"
  names(mrm_out)[names(mrm_out) == "p"] <- "P"

  es_out <- es[, c("label", "level", "x", "rho", "ci_lo", "ci_hi", "n_pairs",
                   "n_cluster", "B", "p_mantel")]
  names(es_out)[names(es_out) == "p_mantel"] <- "P_Mantel"

  prov_out <- do.call(rbind, prov)
  prov_out <- prov_out[!duplicated(prov_out[, c("variable", "role")]), , drop = FALSE]

  sheets <- list(
    Mantel = mant_out,
    MRM = mrm_out,
    Effect_sizes = es_out,
    Saturation = sat[, c("level", "variable", "kind", "bounded", "n_pairs", "mean", "sd",
                         "min", "max", "frac_at_upper", "frac_near_zero", "note")],
    Pair_levels = data.frame(
      level = lv_ok,
      unit = vapply(lv_ok, function(k) if (k == "site") "site (样本级聚合)" else "sample",
                    character(1)),
      n_objects = vapply(lv_ok, function(k) length(L[[k]]$id), integer(1)),
      n_pairs = vapply(lv_ok, function(k) nrow(L[[k]]$P), integer(1)),
      permutation_level = "object-level (置换对象标签，不是置换对)",
      stringsAsFactors = FALSE),
    Variables = prov_out,
    Note = data.frame(Item = c(
      "响应距离", "解释距离", "偏变量", "置换层级", "效应量 CI", "饱和诊断",
      "ecodist 截距", "距离矩阵", "未识别变量", "口径 pair", "口径 site", "饱和结论"),
      Value = c(
        paste0("|delta", y_var, "|（样本级", y_label, "的成对绝对差）"),
        paste(xs_spec, collapse = ", "),
        if (length(pt_spec)) paste(pt_spec, collapse = ", ") else "未指定（CFG$mantel$partial = NULL）",
        paste0("Mantel/MRM 均在【对象级】置换：", nperm, " 次，对象 = 样本(pair 口径) 或站点(site 口径)；报告 n = 对象数，不报对数"),
        paste0("md_boot_rho() 整簇 bootstrap，B = ", B_boot,
               "；pair 口径按【站点对】整簇，site 口径按【站点】整簇"),
        "群落距离类变量有数学上限 1.0；报告均值/极值/顶到上限的对数占比",
        "ecodist::MRM 把截距命名为 Int（不是 Intercept），已在 MRM 表中剔除",
        paste0("距离矩阵单独存 ", basename(dist_path), "（未并进 data.frame）"),
        if (length(dropped_spec)) paste(dropped_spec, collapse = ", ") else "无",
        "全部样本对，统计单元 = 样本对，置换对象 = 样本",
        "站点间平均距离，统计单元 = 站点对，置换对象 = 站点（控制站内伪重复的保守口径）",
        sat_txt),
      stringsAsFactors = FALSE))

  if (detail) {
    sheets$Pair_detail <- do.call(rbind, lapply(lv_ok, function(k) {
      P <- L[[k]]$P; P$level <- k; P
    }))
  }
  sheets <- Filter(function(x) !is.null(x) && nrow(x) > 0, sheets)
  tp <- md_write_xlsx(CFG, "43_mantel_mrm", sheets)

  md_log(CFG, "INFO", "43: ", nrow(es_out), " 个效应量 | ", nrow(mant_out), " 个 Mantel | ",
         nrow(mrm_out), " 个 MRM；饱和诊断：", sat_txt)
  md_result(module, "OK",
            tables = basename(tp), figures = "figures/43_mantel_mrm.pdf",
            notes = sat_txt)
}
