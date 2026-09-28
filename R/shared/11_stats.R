# =============================================================================
# R/shared/11_stats.R —— 统计与距离工具
# -----------------------------------------------------------------------------
# 库文件，无分析逻辑。所有随机过程必须先 md_seed(CFG, offset)。
# =============================================================================

# ---------------------------------------------------------------------------
# 多重校正
# ---------------------------------------------------------------------------
md_p_adjust <- function(p, method = "BH") {
  if (is.null(method) || identical(method, "none") || length(p) == 0) return(p)
  out <- p.adjust(p, method = method)
  out[is.na(p)] <- NA
  out
}

# ---------------------------------------------------------------------------
# 组间比较：描述统计 + 主检验 + 事后 + 显著性字母
# ---------------------------------------------------------------------------
#' @param df data.frame
#' @param value 因变量列名
#' @param group 分组列名（单因子）
#' @param test "anova_tukey" | "kruskal_dunn" | "lmm" | "none"
#' @param p_adjust 事后检验的校正方法
md_group_test <- function(CFG, df, value, group, test = "kruskal_dunn",
                          p_adjust = "BH", min_n = 3, order_levels = NULL) {
  d <- df[!is.na(df[[value]]) & !is.na(df[[group]]), c(value, group), drop = FALSE]
  names(d) <- c("y", "g")
  d$g <- as.character(d$g)
  if (!is.null(order_levels)) d$g <- factor(d$g, levels = order_levels)
  d$g <- droplevels(factor(d$g))
  # 样本量不足的组剔除
  n_tab <- table(d$g)
  drop_g <- names(n_tab)[n_tab < max(2, min_n)]
  if (length(drop_g)) {
    d <- d[!d$g %in% drop_g, , drop = FALSE]
    d$g <- droplevels(d$g)
  }
  lv <- levels(d$g)

  desc <- do.call(rbind, lapply(lv, function(g) {
    x <- d$y[d$g == g]
    data.frame(group = g, n = length(x), mean = mean(x), sd = stats::sd(x),
               se = stats::sd(x) / sqrt(length(x)),
               median = stats::median(x),
               q1 = unname(stats::quantile(x, .25)), q3 = unname(stats::quantile(x, .75)),
               min = min(x), max = max(x), stringsAsFactors = FALSE)
  }))

  trend <- NULL
  if (length(lv) == 2) {
    tr <- tryCatch(stats::cor.test(as.numeric(d$g), d$y, method = "spearman"),
                   error = function(e) NULL)
    if (!is.null(tr))
      trend <- data.frame(test = "Spearman trend (2 groups)", statistic = unname(tr$estimate),
                          df = NA_real_, p = tr$p.value, method = "spearman",
                          effect = unname(tr$estimate), stringsAsFactors = FALSE)
  }

  res <- list(desc = desc, test = NULL, posthoc = NULL, letters = NULL, dropped = drop_g)
  if (test == "none" || length(lv) < 2) return(res)

  if (test == "anova_tukey") {
    fit <- stats::aov(y ~ g, data = d)
    sm  <- summary(fit)[[1]]
    res$test <- data.frame(test = "One-way ANOVA", statistic = sm[["F value"]][1],
                           df = sm[["Df"]][1], p = sm[["Pr(>F)"]][1],
                           method = "aov", effect = NA_real_, stringsAsFactors = FALSE)
    if (requireNamespace("multcomp", quietly = TRUE)) {
      th <- multcomp::TukeyHSD(fit)$g
      res$posthoc <- data.frame(group1 = sub("-.*$", "", rownames(th)),
                                group2 = sub("^.*-", "", rownames(th)),
                                diff = th[, "diff"], p = th[, "p adj"],
                                p_adj = th[, "p adj"], stringsAsFactors = FALSE)
    }
  } else if (test == "kruskal_dunn") {
    kt <- stats::kruskal.test(y ~ g, data = d)
    n  <- nrow(d); k <- length(lv)
    eps2 <- unname(kt$statistic) / ((n^2 - 1) / (n + 1))   # 与 nifH 论文口径一致
    res$test <- data.frame(test = "Kruskal-Wallis", statistic = unname(kt$statistic),
                           df = k - 1, p = kt$p.value, method = "kruskal.test",
                           effect = eps2, stringsAsFactors = FALSE)
    if (requireNamespace("rstatix", quietly = TRUE) && k >= 3) {
      ph <- tryCatch(rstatix::dunn_test(d, y ~ g, p.adjust.method = p_adjust),
                     error = function(e) NULL)
      if (!is.null(ph)) {
        res$posthoc <- data.frame(group1 = ph$group1, group2 = ph$group2,
                                  statistic = ph$statistic, p = ph$p,
                                  p_adj = ph$p.adj, stringsAsFactors = FALSE)
      }
    }
  } else if (test == "lmm") {
    if (!md_need(CFG, "model", "11_alpha/lmm")) return(res)
    md_load(c("lme4", "lmerTest"))
    re <- CFG$design$site_var
    fml <- if (is.null(re)) stats::as.formula("y ~ g")
           else stats::as.formula(paste0("y ~ g + (1|", re, ")"))
    fit <- tryCatch(lmerTest::lmer(fml, data = d, REML = FALSE), error = function(e) NULL)
    if (!is.null(fit)) {
      co <- summary(fit)$coefficients
      i <- grep("^g", rownames(co))
      res$test <- data.frame(test = "LMM (group fixed effect)",
                             statistic = co[i, "t value"], df = co[i, "df"],
                             p = co[i, "Pr(>|t|)"], method = "lmerTest",
                             effect = co[i, "Estimate"], stringsAsFactors = FALSE)
    }
  }
  if (!is.null(trend)) res$test <- rbind(res$test, trend)

  # ---- 显著性字母 ----
  res$letters <- md_sig_letters(res$posthoc, lv, p_col = "p_adj")
  res
}

#' 由两两比较表生成 CLD 字母。无两两结果时返回单字母。
md_sig_letters <- function(posthoc, levels, p_col = "p_adj", alpha = 0.05) {
  lv <- as.character(levels)
  if (length(lv) < 2) return(data.frame(group = lv, letter = "a", stringsAsFactors = FALSE))
  if (is.null(posthoc) || nrow(posthoc) == 0 ||
      !requireNamespace("multcompView", quietly = TRUE)) {
    return(data.frame(group = lv, letter = "a", stringsAsFactors = FALSE))
  }
  p <- posthoc[[p_col]]
  M <- matrix(1, length(lv), length(lv), dimnames = list(lv, lv))
  for (i in seq_len(nrow(posthoc))) {
    g1 <- as.character(posthoc$group1[i]); g2 <- as.character(posthoc$group2[i])
    if (!g1 %in% lv || !g2 %in% lv) next
    v <- if (is.na(p[i])) 1 else if (p[i] < alpha) 0 else 1
    M[g1, g2] <- v; M[g2, g1] <- v
  }
  diag(M) <- 1
  out <- tryCatch(multcompView::multcompLetters(M, threshold = alpha),
                  error = function(e) NULL)
  if (is.null(out)) {
    return(data.frame(group = lv, letter = "a", stringsAsFactors = FALSE))
  }
  data.frame(group = names(out$Letters), letter = as.character(out$Letters),
             stringsAsFactors = FALSE)
}

#' 把字母并回描述统计表（保持因子顺序）
md_attach_letters <- function(desc, letters, group_col = "group") {
  if (is.null(letters)) { desc$letter <- ""; return(desc) }
  desc$letter <- letters$letter[match(desc[[group_col]], letters$group)]
  desc$letter[is.na(desc$letter)] <- ""
  desc
}

#' 每格的样本量（供样本量标注）
md_n_label <- function(df, group) {
  t <- table(df[[group]])
  data.frame(group = names(t), n = as.integer(t), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 距离矩阵工具
# ---------------------------------------------------------------------------
#' dist 对象 -> 下三角长表
md_dist_long <- function(d) {
  m <- as.matrix(d)
  ids <- rownames(m)
  ut <- which(lower.tri(m), arr.ind = TRUE)
  data.frame(sample1 = ids[ut[, 1]], sample2 = ids[ut[, 2]],
             value = m[ut], stringsAsFactors = FALSE)
}

#' 长表 -> dist 对象
md_dist_mat <- function(long, ids, value = "value") {
  m <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  i <- match(long$sample1, ids); j <- match(long$sample2, ids)
  m[cbind(i, j)] <- long[[value]]; m[cbind(j, i)] <- long[[value]]
  stats::as.dist(m)
}

#' 在保持样本顺序的前提下取距离子集
md_subset_dist <- function(d, ids) {
  m <- as.matrix(d)
  keep <- intersect(ids, rownames(m))
  stats::as.dist(m[keep, keep, drop = FALSE])
}

#' 地理距离矩阵（km）。需 env 里有经纬度列。
md_geo_dist <- function(env, ids, lat = "Latitude", lon = "Longitude") {
  if (!requireNamespace("geosphere", quietly = TRUE)) return(NULL)
  if (!all(c(lat, lon) %in% names(env))) return(NULL)
  co <- as.matrix(env[match(ids, env[[1]]), c(lon, lat), drop = FALSE])
  co <- apply(co, 2, as.numeric)
  if (anyNA(co)) return(NULL)
  d <- geosphere::distm(co, fun = geosphere::distHaversine) / 1000
  dimnames(d) <- list(ids, ids)
  stats::as.dist(d)
}

#' 站点级聚合（坐标重复时用；避免"站内距离恒为 0"的伪重复）
#' ⚠ `vars` 语义：**不给** = 自动取所有数值列；**显式给了** = 只保留这些列，
#'   其余列（包括残差协变量）一律丢弃。需要保留别的列请自己并进 `vars`。
md_site_aggregate <- function(df, site_var, vars = NULL, id_col = "sample_id") {
  vars <- vars %||% setdiff(names(df), c(id_col, site_var))
  vars <- vars[vapply(df[vars], is.numeric, logical(1))]
  if (is.null(site_var) || !site_var %in% names(df)) return(df)
  g <- df[[site_var]]
  out <- data.frame(site = unique(as.character(g)), stringsAsFactors = FALSE)
  names(out) <- site_var
  for (v in vars) out[[v]] <- as.numeric(tapply(df[[v]], g, mean, na.rm = TRUE)[out[[site_var]]])
  out
}

#' 某变量在某分组内是否恒为常数（站点级变量的共线守卫）
md_within_unique <- function(df, var, site_var) {
  if (is.null(site_var) || !site_var %in% names(df)) return(NA_integer_)
  x <- df[[var]]
  if (!is.numeric(x)) x <- as.numeric(factor(x))
  u <- tapply(x, df[[site_var]], function(v) length(unique(v[!is.na(v)])))
  max(u, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# 距离-矩阵关联：Mantel / partial Mantel / MRM
# ---------------------------------------------------------------------------
md_mantel <- function(CFG, d1, d2, method = "spearman", permutations = 999,
                      partial = NULL, label = "") {
  if (!md_need(CFG, "ecology", "mantel")) return(NULL)
  md_seed(CFG, 101)
  z <- if (is.null(partial)) NULL else as.dist(as.matrix(partial))
  r <- tryCatch(
    if (is.null(z)) vegan::mantel(d1, d2, method = method, permutations = permutations)
    else vegan::mantel.partial(d1, d2, z, method = method, permutations = permutations),
    error = function(e) NULL)
  if (is.null(r)) return(NULL)
  data.frame(label = label, method = paste0("Mantel", method, if (!is.null(z)) " (partial)"),
             statistic = unname(r$statistic), p = r$signif,
             permutations = permutations, stringsAsFactors = FALSE)
}

#' MRM：成对距离 ~ 解释变量。置换发生在**对象级**（n = 样本数），不是对数级。
md_mrm <- function(CFG, df, y, xs, nperm = 999, label = "") {
  if (!md_need(CFG, "mantel", "mrm")) return(NULL)
  md_seed(CFG, 102)
  fml <- stats::as.formula(paste0(y, " ~ ", paste(xs, collapse = " + ")))
  r <- tryCatch(ecodist::MRM(fml, data = df, nperm = nperm), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  co <- as.data.frame(r$coef)
  co$term <- rownames(co)
  names(co)[1:2] <- c("estimate", "p")
  co$label <- label; co$n_perm <- nperm; co$n_obs <- nrow(df)
  co$R2 <- r$r.squared[1] %||% NA_real_
  co
}

#' 成对级 Spearman rho + **站点对整簇 bootstrap** 置信区间
#' （控制成对数据点的伪重复：以"站点对"为聚类单位整簇重抽）
md_boot_rho <- function(CFG, df, x, y, cluster = NULL, B = 1000, label = "") {
  d <- df[stats::complete.cases(df[, c(x, y)]), , drop = FALSE]
  if (nrow(d) < 5) return(NULL)
  rho0 <- suppressWarnings(stats::cor(d[[x]], d[[y]], method = "spearman"))
  md_seed(CFG, 103)
  if (is.null(cluster) || !cluster %in% names(d)) {
    idx <- replicate(B, sample.int(nrow(d), nrow(d), replace = TRUE), simplify = FALSE)
  } else {
    cl <- split(seq_len(nrow(d)), d[[cluster]])
    idx <- replicate(B, {
      pick <- sample(names(cl), length(cl), replace = TRUE)
      unlist(cl[pick], use.names = FALSE)
    }, simplify = FALSE)
  }
  rs <- vapply(idx, function(i) {
    suppressWarnings(stats::cor(d[[x]][i], d[[y]][i], method = "spearman"))
  }, numeric(1))
  rs <- rs[is.finite(rs)]
  data.frame(label = label, rho = rho0, n_pairs = nrow(d), n_cluster = if (is.null(cluster)) NA_integer_ else length(unique(d[[cluster]])),
             ci_lo = unname(stats::quantile(rs, .025)), ci_hi = unname(stats::quantile(rs, .975)),
             B = length(rs), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 共线诊断
# ---------------------------------------------------------------------------
#' VIF = 1/(1-R^2)（手写版，与 nifH 论文口径一致）
md_vif <- function(df, vars) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) < 2) return(data.frame(variable = vars, vif = NA_real_))
  do.call(rbind, lapply(vars, function(v) {
    others <- setdiff(vars, v)
    f <- stats::as.formula(paste0("`", v, "` ~ ", paste0("`", others, "`", collapse = " + ")))
    r2 <- tryCatch(summary(stats::lm(f, data = df))$r.squared, error = function(e) NA_real_)
    data.frame(variable = v, vif = 1 / (1 - r2), stringsAsFactors = FALSE)
  }))
}

#' 方差膨胀（供 model 模块挑选协变量）
md_drop_collinear <- function(df, vars, thr = 10) {
  keep <- vars
  repeat {
    v <- md_vif(df, keep)
    if (all(is.na(v$vif)) || max(v$vif, na.rm = TRUE) <= thr) break
    drop <- v$variable[which.max(v$vif)]
    keep <- setdiff(keep, drop)
    if (length(keep) < 2) break
  }
  keep
}

# ---------------------------------------------------------------------------
# 抽平 / 亚抽样
# ---------------------------------------------------------------------------
#' 抽平到统一深度（无放回超几何抽样，等价于 phyloseq::rarefy_even_depth 的结果口径）
#' 自实现以便：不依赖 phyloseq、可由 CFG$opts$seed 完全复现、速度更快。
#' @param m ASV x 样本 的整数矩阵
#' @param depth 目标深度
#' @param seed 随机种子
md_rarefy <- function(m, depth, seed = 123L) {
  depth <- as.integer(depth)
  if (depth < 1) stop("抽平深度必须 >= 1", call. = FALSE)
  cs <- colSums(m)
  if (any(cs < depth)) stop(sprintf("有 %d 个样本读长低于抽平深度 %d",
                                    sum(cs < depth), depth), call. = FALSE)
  set.seed(as.integer(seed))
  out <- vapply(seq_len(ncol(m)), function(j) {
    x <- m[, j]
    pool <- rep.int(seq_along(x), as.integer(x))
    sel <- sample.int(length(pool), depth, replace = FALSE)
    tabulate(pool[sel], nbins = length(x))
  }, numeric(nrow(m)))
  dimnames(out) <- dimnames(m)
  out
}

#' 多次亚抽样取均值（降低单次抽平的随机性；结果为非整数）
md_subsample_avg <- function(m, depth, reps = 20L, seed = 123L) {
  reps <- as.integer(reps)
  acc <- matrix(0, nrow(m), ncol(m), dimnames = dimnames(m))
  for (r in seq_len(reps)) acc <- acc + md_rarefy(m, depth, seed = as.integer(seed) + r)
  acc / reps
}

# =============================================================================
# 以下为 2026 毕业论文补充轮新增（模块 10 指示种 / 12 生态位 / 40 RF /
# 44 约束排序 / 45 响应曲线 / 70 SRH / 71 Procrustes / 83 差异类群）
#
# 纪律：
#   · 全部纯计算，不绘图、不写文件；
#   · 所有随机过程先 md_seed(CFG, md_offset(...))，模块里不出现字面量种子；
#   · 统计口径逐项照抄毕业论文原脚本（源文件与行号写在每个函数的注释里）；
#   · 每个函数都对"输入不足"给出 NA / NULL 而不是报错，由模块决定 SKIP 还是降级。
# =============================================================================

#' 取某个随机过程的种子偏移量。可被 `CFG$opts$md_seed_offsets$<name>` 覆盖。
md_offset <- function(CFG, name, default) {
  th <- CFG$opts$md_seed_offsets %||% list()
  as.integer(th[[name]] %||% default)
}

# ---------------------------------------------------------------------------
# 生态位宽度与稀有性（12_niche_breadth；源：1_community_composition.R §7）
# ---------------------------------------------------------------------------
#' Levins' 生态位宽度 B（源：`1_community_composition.R:815-817`）
#'   B = 1 / sum(p^2)，p 为该 ASV 在各样本中的相对丰度构成
#'   `sum(p > 0) < 2` 时返回 NA（原脚本口径：只在 1 个样本出现无法定义宽度）
md_levins_b <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x)) || sum(x, na.rm = TRUE) <= 0) return(NA_real_)
  p <- x / sum(x, na.rm = TRUE)
  if (sum(p > 0, na.rm = TRUE) < 2) return(NA_real_)
  1 / sum(p^2, na.rm = TRUE)
}

#' Shannon 型生态位宽度（备选口径；`niche$index = "shannon"` 时用）
md_shannon_breadth <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x)) || sum(x, na.rm = TRUE) <= 0) return(NA_real_)
  p <- x / sum(x, na.rm = TRUE); p <- p[p > 0]
  if (length(p) < 2) return(NA_real_)
  -sum(p * log(p)) / log(length(p))
}

#' ASV × 样本 相对丰度矩阵 -> ASV 级生态位/稀有性表
#' 源：`1_community_composition.R:773-836`（`classify_2x2` 的 4 象限 + Levins' B）
#' @param m  ASV × 样本 **原始 reads**（未抽平；原脚本用 `ps_rare`）
#' @param a_thr 平均相对丰度阈值（%，原脚本 0.1）
#' @param p_thr 出现率阈值（%，原脚本 10）
#' @param index "levins" | "shannon" | "both"
md_niche_table <- function(m, a_thr = 0.1, p_thr = 10, index = "levins") {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  n_s <- ncol(m)
  rel <- sweep(m, 2, pmax(colSums(m), 1), "/") * 100
  mean_rel <- rowMeans(rel)
  occ <- rowSums(m > 0) / n_s * 100
  cell <- ifelse(mean_rel > a_thr,
                 ifelse(occ > p_thr, "Core", "Conditionally_rare"),
                 ifelse(occ > p_thr, "Permanently_rare", "Transient_rare"))
  out <- data.frame(ASV_ID = rownames(m), Mean_Abund = as.numeric(mean_rel),
                    Occurrence = as.numeric(occ), Cell = as.character(cell),
                    stringsAsFactors = FALSE)
  out$Levins_B <- if (index %in% c("levins", "both")) {
    apply(rel, 1, md_levins_b)
  } else NA_real_
  out$Shannon_B <- if (index %in% c("shannon", "both")) {
    apply(rel, 1, md_shannon_breadth)
  } else NA_real_
  out$N_samples <- as.integer(rowSums(m > 0))
  out
}

#' 样本级"生态策略组成"（源：`1_community_composition.R:513-568`）
#' 每个样本内四类细胞的**相对丰度占比**（%）
md_sample_strategy <- function(m, niche_df, a_thr = 0.1, p_thr = 10) {
  m <- as.matrix(m)
  cells <- stats::setNames(niche_df$Cell, niche_df$ASV_ID)
  lvs <- c("Core", "Permanently_rare", "Conditionally_rare", "Transient_rare")
  tot <- pmax(colSums(m), 1)
  out <- data.frame(SampleID = colnames(m), stringsAsFactors = FALSE)
  for (cl in lvs) {
    idx <- names(cells)[!is.na(cells) & cells == cl]
    idx <- intersect(idx, rownames(m))
    sub <- if (length(idx) == 0) rep(0, ncol(m)) else colSums(m[idx, , drop = FALSE])
    out[[paste0(cl, "_pct")]] <- as.numeric(sub) / tot * 100
    out[[paste0(cl, "_Abundance")]] <- as.numeric(sub)
  }
  out$Richness <- as.integer(colSums(m > 0))
  out
}

# ---------------------------------------------------------------------------
# 差异类群（83_taxon_diff；源：S4_supplementary_taxa.R §3 / 原 9_differential_analysis.R）
# ---------------------------------------------------------------------------
#' 中心对数比变换（源：`S4_supplementary_taxa.R:781-791`）
#' 伪计数缺省 = 最小正值的 1/2（与原脚本一致）
md_clr <- function(mat, pseudocount = NULL) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  if (is.null(pseudocount)) {
    nz <- mat[mat > 0]
    pseudocount <- if (length(nz)) min(nz) / 2 else 0.5
  }
  p <- mat + pseudocount
  rel <- p / rowSums(p)
  lg <- log(rel)
  lg - rowMeans(lg)
}

#' Cliff's delta（源：`S4_supplementary_taxa.R:809-824`）
#' 返回"绝对效应量最大"的那一对的 delta（原脚本口径，不是加权平均）
md_cliff_delta <- function(x, group) {
  x <- as.numeric(x); group <- as.character(group)
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]; group <- group[ok]
  gs <- unique(group)
  if (length(gs) < 2) return(NA_real_)
  deltas <- numeric(0)
  for (i in seq_len(length(gs) - 1)) {
    for (j in seq_len(length(gs) - i) + i) {
      g1 <- x[group == gs[i]]; g2 <- x[group == gs[j]]
      if (length(g1) < 2 || length(g2) < 2) next
      comp <- outer(g1, g2, "-")
      deltas <- c(deltas, (sum(comp > 0) - sum(comp < 0)) / (length(g1) * length(g2)))
    }
  }
  if (!length(deltas)) return(0)
  as.numeric(deltas[which.max(abs(deltas))])
}

#' 类群级差异分析：CLR + Kruskal-Wallis + BH + Cliff's delta + 富集组
#' 源：`S4_supplementary_taxa.R:828-872`（`run_diff_analysis`）
#' （函数名带 `_kw` 后缀，避免与模块 83 的入口 `md_taxon_diff()` 撞名）
#' @param clr_mat 类群 × 样本 的 CLR 矩阵
#' @param rel_mat 类群 × 样本 的相对丰度（%）矩阵，同维度
#' @param meta_df 含 `sample_id` 与分组列的 data.frame
md_diff_taxa_kw <- function(clr_mat, rel_mat, meta_df, group_var,
                            sample_col = "sample_id", min_n = 3, adjust = "BH") {
  ids <- as.character(meta_df[[sample_col]])
  g   <- as.character(meta_df[[group_var]])
  taxa <- intersect(rownames(clr_mat), rownames(rel_mat))
  if (length(taxa) == 0) return(data.frame())
  rows <- lapply(taxa, function(tx) {
    v <- as.numeric(clr_mat[tx, ids])
    d <- data.frame(y = v, g = g)
    d <- d[is.finite(d$y) & !is.na(d$g), , drop = FALSE]
    if (nrow(d) < 2) return(NULL)
    if (any(table(d$g) < min_n)) return(NULL)
    if (length(unique(d$g)) < 2) return(NULL)
    kw <- tryCatch(stats::kruskal.test(y ~ g, data = d), error = function(e) NULL)
    if (is.null(kw)) return(NULL)
    data.frame(Taxon = tx, KW_statistic = unname(kw$statistic), KW_p = kw$p.value,
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  res <- do.call(rbind, rows)
  res$KW_padj <- md_p_adjust(res$KW_p, adjust)
  res$Cliff_delta <- vapply(res$Taxon, function(tx)
    md_cliff_delta(as.numeric(clr_mat[tx, ids]), g), numeric(1))
  grp_lv <- unique(g[!is.na(g)])
  for (gv in grp_lv) {
    res[[paste0("Rel_", gv)]] <- vapply(res$Taxon, function(tx) {
      i <- which(g == gv & is.finite(as.numeric(rel_mat[tx, ids])))
      if (!length(i)) return(NA_real_)
      mean(as.numeric(rel_mat[tx, ids])[i], na.rm = TRUE)
    }, numeric(1))
  }
  rel_cols <- paste0("Rel_", grp_lv)
  if (length(rel_cols) == 1) {
    res$Enrich_group <- grp_lv[1]
    res$Mean_Rel_max  <- res[[rel_cols[1]]]
  } else {
    m <- as.matrix(res[, rel_cols, drop = FALSE])
    mi <- apply(m, 1, function(z) if (all(is.na(z))) NA_integer_ else which.max(z))
    res$Enrich_group <- ifelse(is.na(mi), NA_character_, grp_lv[mi])
    res$Mean_Rel_max <- ifelse(is.na(mi), NA_real_, m[cbind(seq_len(nrow(m)), pmax(mi, 1))])
  }
  res[order(res$KW_padj), , drop = FALSE]
}

#' 差异类群的 Dunn 事后检验（源：`S4_supplementary_taxa.R:874-889`）
md_taxon_dunn <- function(clr_mat, meta_df, group_var, taxa,
                          sample_col = "sample_id", adjust = "BH") {
  if (!length(taxa) || !requireNamespace("rstatix", quietly = TRUE)) return(data.frame())
  ids <- as.character(meta_df[[sample_col]]); g <- as.character(meta_df[[group_var]])
  out <- lapply(taxa, function(tx) {
    v <- as.numeric(clr_mat[tx, ids])
    d <- data.frame(y = v, g = g)
    d <- d[is.finite(d$y) & !is.na(d$g), , drop = FALSE]
    if (length(unique(d$g)) < 2) return(NULL)
    r <- tryCatch(rstatix::dunn_test(d, y ~ g, p.adjust.method = adjust),
                  error = function(e) NULL)
    if (is.null(r) || !nrow(r)) return(NULL)
    r$Taxon <- tx
    r
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# 指示种 IndVal（10_composition 的 indicator 面板 / 81_taxon_group）
# 源：S4_supplementary_taxa.R:117-171（indicspecies::multipatt, IndVal.g, 999 置换）
# ---------------------------------------------------------------------------
#' @param comm  ASV × 样本 的**相对丰度**矩阵（行 = ASV，列 = 样本）
#' @param group 长度 = ncol(comm) 的分组向量
#' @param func   "IndVal.g"（推荐，原脚本）/ "IndVal" / "r.g" / "r"
#' @param nperm  置换次数（原脚本 999）
#' @param adjust p 值校正
md_indval <- function(CFG, comm, group, func = "IndVal.g", nperm = 999,
                      adjust = "BH", seed_offset = 8100L) {
  if (!md_need(CFG, "indicator", "indval")) return(NULL)
  grp <- as.character(group); grp[is.na(grp)] <- NA_character_
  keep <- !is.na(grp)
  if (sum(keep) < 4 || length(unique(grp[keep])) < 2) return(NULL)
  m <- as.matrix(comm[, keep, drop = FALSE])
  storage.mode(m) <- "double"
  # multipatt 要求行 = 样本
  X <- t(m)
  cl <- grp[keep]
  md_seed(CFG, md_offset(CFG, "indval", seed_offset))
  mp <- tryCatch(indicspecies::multipatt(X, cluster = cl, func = func,
                                         control = permute::how(nperm = nperm)),
                 error = function(e) NULL)
  if (is.null(mp) || is.null(mp$sign)) return(NULL)
  s <- mp$sign
  out <- data.frame(ASV_ID = rownames(s), Level = NA_character_,
                    IndVal = as.numeric(s$stat), p.value = as.numeric(s$p.value),
                    stringsAsFactors = FALSE)
  # multipatt 的 s.index 列名形如 "1"、"1+2"（0/1 指示向量编码）
  dec <- attr(mp, "group.levels")
  if (is.null(dec)) dec <- levels(factor(cl))
  idx_cols <- setdiff(names(s), c("stat", "p.value"))
  out$Level <- vapply(seq_len(nrow(s)), function(i) {
    v <- as.numeric(s[i, idx_cols])
    paste(dec[which(v == 1)], collapse = "+")
  }, character(1))
  # 原脚本只保留单水平 + 多水平组合都存在的表；这里一律保留并标注
  out$fdr <- md_p_adjust(out$p.value, adjust)
  out$sig <- ifelse(out$fdr < 0.05, "**", "")
  out[order(out$Level, -out$IndVal), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Scheirer-Ray-Hare 双因子非参数检验（70_response_pattern 的备选口径）
# 源：8_ARA_function.R:200-236 —— rcompanion 优先，缺包退回 ART + 双因素 ANOVA
# ---------------------------------------------------------------------------
#' @param df data.frame
#' @param y  因变量列名（连续）
#' @param f1,f2 两个分组因子列名
md_srh <- function(CFG, df, y, f1, f2) {
  d <- df[is.finite(suppressWarnings(as.numeric(df[[y]]))) &
            !is.na(df[[f1]]) & !is.na(df[[f2]]), , drop = FALSE]
  if (nrow(d) < 6) return(NULL)
  d$.y <- as.numeric(d[[y]])
  d$.f1 <- factor(as.character(d[[f1]]))
  d$.f2 <- factor(as.character(d[[f2]]))
  if (nlevels(d$.f1) < 2 || nlevels(d$.f2) < 2) return(NULL)
  fml <- stats::as.formula(".y ~ .f1 + .f2")
  if (requireNamespace("rcompanion", quietly = TRUE)) {
    s <- tryCatch(rcompanion::scheirerRayHare(fml, data = d), error = function(e) NULL)
    if (!is.null(s)) {
      o <- as.data.frame(s)
      o$term <- rownames(o)
      o$epsilon_sq <- o$H / (o$H + nrow(d) - o$Df - 1)   # 原脚本口径
      o$method <- "rcompanion::scheirerRayHare"
      o$n <- nrow(d)
      return(o[, c("term", "Df", "H", "p.value", "epsilon_sq", "method", "n")])
    }
  }
  # 备选：对齐秩转换（ART）+ 双因素 ANOVA（原脚本注释里的 fallback）
  d$.ar <- rank(d$.y - stats::ave(d$.y, d$.f1, FUN = stats::median) -
                  stats::ave(d$.y, d$.f2, FUN = stats::median) + stats::median(d$.y))
  a <- tryCatch(stats::anova(stats::lm(.ar ~ .f1 * .f2, data = d)), error = function(e) NULL)
  if (is.null(a)) return(NULL)
  data.frame(term = c(paste0(f1), paste0(f2), paste0(f1, ":", f2)),
             Df = a$Df[1:3], H = a$`Sum Sq`[1:3], p.value = a$`Pr(>F)`[1:3],
             epsilon_sq = a$`Sum Sq`[1:3] / sum(a$`Sum Sq`),
             method = "ART + two-way ANOVA (fallback)", n = nrow(d),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Procrustes（71_response_linkage 的 "procrustes" 口径）
# 源：8_ARA_function.R:1173-1229
# ---------------------------------------------------------------------------
#' @param X,Y 两个样本 × 维度的排序得分矩阵（行名 = sample_id，顺序需一致或可对齐）
#' @param permutations PROTEST 置换次数
md_procrustes_res <- function(CFG, X, Y, permutations = 999, seed_offset = 8120L) {
  ids <- intersect(rownames(X), rownames(Y))
  if (length(ids) < 5) return(NULL)
  A <- as.matrix(X[ids, , drop = FALSE]); B <- as.matrix(Y[ids, , drop = FALSE])
  k <- min(ncol(A), ncol(B))
  A <- A[, seq_len(k), drop = FALSE]; B <- B[, seq_len(k), drop = FALSE]
  A <- A[, colSums(is.na(A)) == 0, drop = FALSE]
  if (!ncol(A)) return(NULL)
  md_seed(CFG, md_offset(CFG, "procrustes", seed_offset))
  pr <- tryCatch(vegan::procrustes(A, B, symmetric = TRUE), error = function(e) NULL)
  if (is.null(pr)) return(NULL)
  pt <- tryCatch(vegan::protest(A, B, permutations = permutations, symmetric = TRUE),
                 error = function(e) NULL)
  res <- sqrt(rowSums((pr$Yrot - B)^2))
  list(residuals = data.frame(sample_id = ids, Residual = as.numeric(res),
                              stringsAsFactors = FALSE),
       summary = data.frame(m2 = unname(pr$ss),
                            correlation = if (is.null(pt)) NA_real_ else unname(pt$t0),
                            p = if (is.null(pt)) NA_real_ else unname(pt$signif),
                            permutations = permutations, n = length(ids),
                            stringsAsFactors = FALSE))
}

# ---------------------------------------------------------------------------
# 随机森林重要性（40_env_summary 的 rf 面板 / 71 的 "rf" 口径）
# 源：3_environmental_drivers.R:1199-1258（ranger，permutation importance，缩放到 max=1）
#     8_ARA_function.R:1130+（randomForest 口径）
# ---------------------------------------------------------------------------
#' @param engine "ranger"（推荐，与毕业论文 Fig3b-B 一致）| "randomForest" | "cforest"
#' 本函数**永不抛错**：任何异常都降级为 `md_warn()` + 返回 NULL，由调用模块决定跳过该面板。
md_rf_fit <- function(CFG, df, y, xs, engine = "ranger", num_trees = 2000L,
                      importance = "permutation", seed_offset = 8130L) {
  tryCatch(
    md_rf_fit_impl(CFG, df, y, xs, engine = engine, num_trees = num_trees,
                   importance = importance, seed_offset = seed_offset),
    error = function(e) {
      md_warn(CFG, "md_rf_fit() 失败，该面板跳过：", conditionMessage(e))
      NULL
    })
}

#' `md_rf_fit()` 的实现体（不直接调用；见上面的包装）
md_rf_fit_impl <- function(CFG, df, y, xs, engine = "ranger", num_trees = 2000L,
                           importance = "permutation", seed_offset = 8130L) {
  xs <- intersect(xs, names(df))
  d <- df[stats::complete.cases(df[, c(y, xs), drop = FALSE]), c(y, xs), drop = FALSE]
  if (nrow(d) < 10 || length(xs) < 1) return(NULL)
  d <- d[is.finite(suppressWarnings(as.numeric(d[[y]]))), , drop = FALSE]
  if (nrow(d) < 10) return(NULL)
  fml <- stats::as.formula(paste0("`", y, "` ~ ."))
  md_seed(CFG, md_offset(CFG, "rf", seed_offset))
  imp <- NULL; r2 <- NA_real_; fit <- NULL; used <- engine

  if (engine == "ranger" && md_need(CFG, "rf", "rf/ranger")) {
    fit <- tryCatch(ranger::ranger(fml, data = d, num.trees = as.integer(num_trees),
                                   importance = importance, seed = md_offset(CFG, "rf", seed_offset)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      imp <- data.frame(Variable = names(fit$variable.importance),
                        Importance = as.numeric(fit$variable.importance),
                        stringsAsFactors = FALSE)
      r2 <- unname(fit$r.squared) %||% NA_real_
    }
  }
  if (is.null(imp) && engine != "cforest" && md_need(CFG, "rf", "rf/randomForest")) {
    used <- "randomForest"
    fit <- tryCatch(randomForest::randomForest(fml, data = d,
                                               ntree = as.integer(num_trees),
                                               importance = TRUE), error = function(e) NULL)
    if (!is.null(fit)) {
      ii <- tryCatch(randomForest::importance(fit), error = function(e) NULL)
      if (is.null(ii)) { imp <- NULL; fit <- NULL }
      else {
        col <- if ("%IncMSE" %in% colnames(ii)) "%IncMSE" else colnames(ii)[1]
        imp <- data.frame(Variable = rownames(ii), Importance = as.numeric(ii[, col]),
                          stringsAsFactors = FALSE)
        r2 <- unname(fit$rsq[length(fit$rsq)]) %||% NA_real_
      }
    }
  }
  if (is.null(imp) && md_need(CFG, "rf", "rf/cforest")) {
    used <- "cforest"
    fit <- tryCatch(party::cforest(fml, data = d,
                                   controls = party::cforest_unbiased(ntree = as.integer(num_trees))),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      vi <- tryCatch(party::varimp(fit), error = function(e) NULL)
      if (!is.null(vi)) {
        imp <- data.frame(Variable = names(vi), Importance = as.numeric(vi),
                          stringsAsFactors = FALSE)
        r2 <- NA_real_
      }
    }
  }
  if (is.null(imp)) return(NULL)
  imp$Importance_scaled <- imp$Importance / max(imp$Importance, na.rm = TRUE)
  imp <- imp[order(-imp$Importance), , drop = FALSE]
  rownames(imp) <- NULL
  # 部分依赖（可选，只有装了 pdp 才算；原论文 Fig3b 没有 PDP，属本包可选面板）
  pdp_df <- NULL
  if (isTRUE(CFG$env$rf$pdp %||% FALSE) && !is.null(fit) &&
      used %in% c("ranger", "randomForest") && md_need(CFG, "rf", "rf/pdp")) {
    top <- head(imp$Variable, as.integer(CFG$env$rf$pdp_top %||% 3L))
    pdp_df <- lapply(top, function(v) {
      g <- tryCatch(pdp::partial(fit, pred.var = v, train = d,
                                 grid.resolution = 30, plot = FALSE),
                    error = function(e) NULL)
      if (is.null(g) || !nrow(g)) return(NULL)
      # ⚠ 每个变量的 partial 输出列名不同（列名 = 该预测变量名 + yhat），
      #   直接 rbind 会报 "names do not match previous names"。
      #   统一成 (Variable, x, yhat) 三列再合并。
      pv <- setdiff(names(g), "yhat")[1]
      if (is.na(pv)) return(NULL)
      data.frame(Variable = v, x = as.numeric(g[[pv]]), yhat = as.numeric(g$yhat),
                 stringsAsFactors = FALSE)
    })
    pdp_df <- pdp_df[!vapply(pdp_df, is.null, logical(1))]
    pdp_df <- if (length(pdp_df)) do.call(rbind, pdp_df) else NULL
  }
  list(importance = imp, r2 = r2, n = nrow(d), engine = used,
       predictors = xs, pdp = pdp_df)
}

# ---------------------------------------------------------------------------
# GAM / GAMM 响应曲线（45_response_curves；源：3_environmental_drivers.R:1261-1345）
# ---------------------------------------------------------------------------
#' @param engine  "gam"（推荐，论文 Fig3b-C 口径）| "gamm" | "lm"
#' @param k       平滑基维数（原脚本 k = 4）
#' @param random  随机效应项字符串，例如 "(1|Site)"；给出则走 gamm
#' @param grid_n  预测网格点数（原脚本 100）
#' @param edf_nonlin 线性判据：edf > 该值才认为非线性
md_gam_curve <- function(CFG, df, y, x, engine = "gam", k = 4L, random = NULL,
                         grid_n = 100L, edf_nonlin = 1.1) {
  d <- df[stats::complete.cases(df[, c(y, x), drop = FALSE]), c(y, x), drop = FALSE]
  d <- d[is.finite(suppressWarnings(as.numeric(d[[y]]))), , drop = FALSE]
  if (nrow(d) < 10) return(NULL)
  d$.y <- as.numeric(d[[y]]); d$.x <- as.numeric(d[[x]])
  eng <- engine
  if (!is.null(random) && nzchar(random) && engine == "gam") eng <- "gamm"
  # 原脚本：gam(s(pred, k = k), method = "REML")
  f <- NULL
  if (eng %in% c("gam", "gamm") && md_need(CFG, "gam", "gam_curve")) {
    rhs <- sprintf("s(.x, k = %d)", as.integer(k))
    if (eng == "gamm") {
      if (!is.null(random) && nzchar(random))
        f <- tryCatch(mgcv::gamm(stats::as.formula(paste0(".y ~ ", rhs, " + ", random)),
                                 data = d, method = "REML"), error = function(e) NULL)
      if (is.null(f))
        f <- tryCatch(mgcv::gam(stats::as.formula(paste0(".y ~ ", rhs)),
                                data = d, method = "REML"), error = function(e) NULL)
    } else {
      f <- tryCatch(mgcv::gam(stats::as.formula(paste0(".y ~ ", rhs)),
                              data = d, method = "REML"), error = function(e) NULL)
    }
  }
  if (is.null(f)) {
    f <- stats::lm(.y ~ .x, data = d)
    eng <- "lm"
  }
  sm <- summary(f)
  edf  <- if (!is.null(sm$edf)) as.numeric(sm$edf[1]) else 1
  psmo <- if (!is.null(sm$s.pv)) as.numeric(sm$s.pv[1]) else NA_real_
  r2   <- as.numeric(sm$r.sq %||% sm$adj.r.squared %||% NA_real_)
  gx <- seq(min(d$.x, na.rm = TRUE), max(d$.x, na.rm = TRUE), length.out = as.integer(grid_n))
  nd <- data.frame(.x = gx)
  if (!is.null(random) && nzchar(random) && inherits(f, "gamm")) {
    pv <- tryCatch(stats::predict(f$gam, newdata = nd, se.fit = TRUE), error = function(e) NULL)
  } else {
    pv <- tryCatch(stats::predict(f, newdata = nd, se.fit = TRUE), error = function(e) NULL)
  }
  if (is.null(pv)) return(NULL)
  curve <- data.frame(x = gx, fit = as.numeric(pv$fit),
                      lo = as.numeric(pv$fit - 1.96 * pv$se.fit),
                      hi = as.numeric(pv$fit + 1.96 * pv$se.fit))
  list(curve = curve, points = data.frame(x = d$.x, y = d$.y),
       fit = data.frame(term = x, engine = eng, edf = edf, p_smooth = psmo,
                        r2 = r2, n = nrow(d), nonlinear = isTRUE(edf > edf_nonlin),
                        stringsAsFactors = FALSE),
       model = f)
}

# ---------------------------------------------------------------------------
# 约束排序 db-RDA / RDA / CCA / capscale（44_constrained_ord）
# 源：3_environmental_drivers.R:855-928
# ---------------------------------------------------------------------------
#' @param Y  样本 × 变量 的响应矩阵（`response = "community"` 时），
#'           或已是 `dist` 对象（`response = "distance"` 时用 `Y` 作距离）
#' @param X  样本 × 解释变量 的 data.frame（行名 = 样本 ID）
#' @param method  "dbrda"（推荐，论文口径）| "rda" | "cca" | "capscale"
#' @param transform 群落矩阵预处理："hellinger"（推荐）| "none"
#' @param vif_thr VIF 逐步剔除阈值（原脚本 15）
#' （函数名带 `_fit` 后缀，避免与模块 44 的入口 `md_constrained_ord()` 撞名）
md_cord_fit <- function(CFG, Y, X, method = "dbrda", transform = "hellinger",
                        vif_thr = 15, permutations = 999,
                        seed_offset = 8140L) {
  if (!md_need(CFG, "ecology", "constrained_ord")) return(NULL)
  ids <- intersect(rownames(X), if (inherits(Y, "dist")) labels(Y) else rownames(Y))
  if (length(ids) < 5) return(NULL)
  Xd <- as.data.frame(X[ids, , drop = FALSE])
  Xd <- Xd[, vapply(Xd, function(z) length(unique(z[!is.na(z)])) > 1, logical(1)), drop = FALSE]
  Xd <- Xd[, vapply(Xd, is.numeric, logical(1)), drop = FALSE]
  Xd <- Xd[stats::complete.cases(Xd), , drop = FALSE]
  if (ncol(Xd) < 1 || nrow(Xd) < 5) return(NULL)
  ids2 <- rownames(Xd)

  if (inherits(Y, "dist")) {
    d <- stats::as.dist(as.matrix(Y)[ids2, ids2])
    resp_label <- "distance"
  } else {
    M <- as.matrix(Y[ids2, , drop = FALSE])
    storage.mode(M) <- "double"
    M <- M[, colSums(M) > 0, drop = FALSE]
    if (transform == "hellinger") M <- vegan::decostand(M, method = "hellinger")
    resp_label <- "community"
  }
  fml <- if (inherits(Y, "dist")) {
    stats::as.formula(paste("d ~", paste(colnames(Xd), collapse = " + ")))
  } else {
    stats::as.formula(paste("M ~", paste(colnames(Xd), collapse = " + ")))
  }
  fit_once <- function(xdf, M = NULL) {
    dat <- as.data.frame(xdf)
    if (!is.null(M)) dat$M <- M
    else dat$d <- d
    md_seed(CFG, md_offset(CFG, "constrained_ord", seed_offset))
    switch(method,
           dbrda    = vegan::dbrda(fml, data = dat),
           capscale = vegan::capscale(fml, data = dat),
           rda      = vegan::rda(fml, data = dat),
           cca      = vegan::cca(fml, data = dat),
           vegan::dbrda(fml, data = dat))
  }
  env_vars <- colnames(Xd)
  fit <- tryCatch(if (inherits(Y, "dist")) fit_once(Xd) else fit_once(Xd, M),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  vif_vals <- tryCatch(vegan::vif.cca(fit), error = function(e) NULL)
  guard <- 0L
  while (!is.null(vif_vals) && any(vif_vals > vif_thr, na.rm = TRUE) && length(env_vars) > 2 && guard < 50L) {
    guard <- guard + 1L
    worst <- names(which.max(vif_vals))
    env_vars <- setdiff(env_vars, worst)
    Xd <- Xd[, env_vars, drop = FALSE]
    fml <- if (inherits(Y, "dist")) {
      stats::as.formula(paste("d ~", paste(env_vars, collapse = " + ")))
    } else {
      stats::as.formula(paste("M ~", paste(env_vars, collapse = " + ")))
    }
    fit <- tryCatch(if (inherits(Y, "dist")) fit_once(Xd) else fit_once(Xd, M),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    vif_vals <- tryCatch(vegan::vif.cca(fit), error = function(e) NULL)
  }
  md_seed(CFG, md_offset(CFG, "constrained_ord", seed_offset) + 1L)
  an <- tryCatch(stats::anova(fit, permutations = permutations), error = function(e) NULL)
  tr <- tryCatch(stats::anova(fit, by = "terms", permutations = permutations),
                 error = function(e) NULL)
  r2 <- tryCatch(vegan::RsquareAdj(fit), error = function(e) NULL)
  sites <- tryCatch(as.data.frame(vegan::scores(fit, display = "sites", choices = 1:2)),
                    error = function(e) NULL)
  if (is.null(sites)) return(NULL)
  sites$sample_id <- rownames(sites)
  envsc <- tryCatch(as.data.frame(vegan::scores(fit, display = "bp", choices = 1:2)),
                    error = function(e) NULL)
  if (!is.null(envsc)) { envsc$env_var <- rownames(envsc); rownames(envsc) <- NULL }
  eig <- tryCatch(as.numeric(fit$CCA$eig[1:2] / fit$tot.chi * 100), error = function(e) c(NA, NA))
  list(sites = sites, env_scores = envsc, vif = data.frame(variable = names(vif_vals),
                                                           vif = as.numeric(vif_vals),
                                                           stringsAsFactors = FALSE),
       anova = an, terms = tr,
       R2adj = as.numeric(r2$adj.r.squared %||% r2$adj.r %||% NA_real_),
       R2 = as.numeric(r2$r.squared %||% r2$r %||% NA_real_),
       var_exp = eig, method = method, response = resp_label,
       variables = env_vars, n = nrow(Xd))
}
