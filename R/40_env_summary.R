# =============================================================================
# R/40_env_summary.R —— 环境变量表征（描述统计 / VIF / 相关矩阵 / 组均值热图 / PCA 双标图）
# -----------------------------------------------------------------------------
# 口径来源（nifH 论文 10_Supp_env_SEM.R 与 11_shared.R §1，只读参考，不修改）：
#   · 组均值热图：行 = 分组（Zone×Vegetation），列 = 环境变量，值 = 各组均值的列内 z-score
#   · 相关矩阵：Spearman，**只对上三角独立配对做 BH 校正，再镜像**回下三角（nifH 口径）
#   · VIF：手写 1/(1-R²)，与 nifH 口径一致（本包走共享 md_vif()）
#
# ★ 跨模块产物契约（MODULE_SPEC §8）：
#   必须写 cache/40_env.rds = list(env, vif, cor, pca, soil_pca, derived, ...)
#   其中 derived 至少含 sample_id 与 Soil_PC1（50/52/71 通过 md_derived() 取）。
#   Soil_PC1 的定义全项目唯一：prcomp(scale(X)) 的第一主轴得分，
#   X = CFG$env$pca_vars，是否 log 由 CFG$env$pca_log 决定，pH 在 log 时豁免。
#   —— 本模块是全项目**唯一**计算它的地方；其它模块一律用 md_derived()。
#
# 输出：figures/40_env_summary.pdf（A PCA 双标图 / B 相关矩阵 / C 组均值热图）
#       results/40_env_summary.xlsx（Descriptive / VIF / Correlations / Soil_PCA / Note）
#       cache/40_env.rds
# =============================================================================

# ---------------------------------------------------------------------------
# local_* 辅助
# ---------------------------------------------------------------------------

#' 与 nifH make_soil_pca() 等价的唯一定义：prcomp(scale(X))，log 时 pH 豁免
local_pca <- function(df, vars, log_transform = FALSE, eps = 1e-4) {
  vars <- intersect(vars, names(df))
  vars <- vars[vapply(df[vars], is.numeric, logical(1))]
  if (length(vars) < 2) return(NULL)
  X <- as.data.frame(df[, vars, drop = FALSE])
  if (isTRUE(log_transform)) {
    lg <- setdiff(vars, "pH")
    if (length(lg)) X[lg] <- log10(as.matrix(X[lg]) + eps)
  }
  cc <- stats::complete.cases(X)
  if (sum(cc) < 3) return(NULL)
  Xc <- X[cc, , drop = FALSE]
  # 常量列会让 prcomp 失效，先剔除（并记录，写在 Note 里）
  sds <- vapply(Xc, stats::sd, numeric(1))
  dropped <- names(sds)[!is.finite(sds) | sds <= 0]
  keep <- setdiff(names(Xc), dropped)
  if (length(keep) < 2) return(NULL)
  p <- stats::prcomp(scale(as.matrix(Xc[, keep, drop = FALSE])))
  ve <- summary(p)$importance[2, ]
  scores <- matrix(NA_real_, nrow(df), ncol(p$x),
                   dimnames = list(NULL, paste0("PC", seq_len(ncol(p$x)))))
  scores[cc, ] <- p$x
  list(
    scores  = scores,
    ve      = ve,
    rotation = p$rotation,
    sdev    = p$sdev,
    vars    = keep,
    dropped = dropped,
    log     = log_transform,
    n_used  = sum(cc),
    complete = cc,
    pca_obj = p)
}

#' 上三角独立配对 BH 校正后镜像的相关矩阵（nifH 口径）
local_cor_matrix <- function(X, method = "spearman", padj = "BH", digits = 6) {
  X <- as.matrix(X)
  k <- ncol(X)
  if (k < 2) return(NULL)
  res <- suppressWarnings(Hmisc::rcorr(X, type = method))
  R <- res$r; P <- res$P
  dimnames(R) <- list(colnames(X), colnames(X))
  dimnames(P) <- dimnames(R)
  ut <- upper.tri(P)
  Padj <- matrix(NA_real_, k, k, dimnames = dimnames(P))
  if (!is.null(padj) && !identical(padj, "none")) {
    Padj[ut] <- stats::p.adjust(P[ut], method = padj)
  } else {
    Padj[ut] <- P[ut]
  }
  Padj[lower.tri(Padj)] <- t(Padj)[lower.tri(Padj)]
  diag(Padj) <- NA_real_
  stars <- matrix(ifelse(is.na(Padj), "",
                  ifelse(Padj < 0.001, "***",
                  ifelse(Padj < 0.01, "**",
                  ifelse(Padj < 0.05, "*", "")))), k, k, dimnames = dimnames(P))
  diag(stars) <- ""
  # 长表：只保留下三角（每对一次），P_adj 用镜像后的值
  lt <- which(upper.tri(R), arr.ind = TRUE)
  long <- data.frame(
    Var1 = colnames(X)[lt[, 1]], Var2 = colnames(X)[lt[, 2]],
    r = R[lt], P = P[lt], P_adj = Padj[lt],
    sig = stars[lt], n = res$n[1], method = method,
    stringsAsFactors = FALSE)
  long <- long[order(long$P), , drop = FALSE]
  list(r = R, p = P, p_adj = Padj, stars = stars, n = res$n[1],
       method = method, p_adjust = padj, long = long)
}

#' 分组均值矩阵 + 列内 z-score
local_group_means <- function(df, group_vars, vars, order_levels = NULL) {
  g <- interaction(df[, group_vars, drop = FALSE], sep = "_", drop = TRUE, lex.order = TRUE)
  lv <- levels(g)
  if (!is.null(order_levels)) {
    lv <- as.vector(outer(order_levels[[group_vars[1]]] %||% levels(factor(df[[group_vars[1]]])),
                          order_levels[[group_vars[2]]] %||% levels(factor(df[[group_vars[2]]])),
                          function(a, b) paste(a, b, sep = "_")))
    lv <- lv[lv %in% levels(g)]
  }
  vars <- intersect(vars, names(df))
  vars <- vars[vapply(df[vars], is.numeric, logical(1))]
  if (!length(vars)) return(NULL)
  gm <- matrix(NA_real_, length(lv), length(vars), dimnames = list(lv, vars))
  for (i in seq_along(lv)) {
    idx <- as.character(g) == lv[i]
    for (v in vars) gm[i, v] <- mean(df[[v]][idx], na.rm = TRUE)
  }
  gm[!is.finite(gm)] <- NA_real_
  z <- gm
  for (j in seq_along(vars)) {
    x <- gm[, j]; s <- stats::sd(x, na.rm = TRUE)
    z[, j] <- if (is.finite(s) && s > 0) (x - mean(x, na.rm = TRUE)) / s else 0
  }
  n <- as.integer(table(g)[lv])
  list(means = gm, z = z, n = n, levels = lv)
}

#' 现算若干 α 指数（供 `env$rf` 的因变量用），返回 data.frame(sample_id, <指数>)
#'
#' 为什么需要它：`md_core_analysis()` 只含 meta + env + 派生变量，**不含** α 指数
#' （Shannon / Simpson / Observed / Pielou 由模块 11 算过但不并表）。
#' 于是"环境 -> α 多样性"的随机森林如果直接写 `response = "Shannon"` 会取不到列。
#' 这里按需现算，口径与模块 11 一致（同一 `md_core_asv(CFG, "alpha")`）。
local_rf_alpha <- function(CFG, want) {
  want <- intersect(as.character(want), c("Observed", "Shannon", "Simpson", "Pielou"))
  if (!length(want)) return(NULL)
  m <- tryCatch(md_core_asv(CFG, "alpha"), error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0) return(NULL)
  X <- t(as.matrix(m))                     # 样本 x ASV
  idx <- list(
    Observed = rowSums(X > 0),
    Shannon  = vegan::diversity(X, index = "shannon"),
    Simpson  = vegan::diversity(X, index = "simpson"),
    Pielou   = vegan::diversity(X, index = "shannon") / log(pmax(rowSums(X > 0), 1))
  )
  smp <- CFG$data$sample_col %||% "sample_id"
  out <- data.frame(sample_id = rownames(X), stringsAsFactors = FALSE)
  names(out)[1] <- smp
  for (nm in want) out[[nm]] <- as.numeric(idx[[nm]])
  out
}

# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------
md_env_summary <- function(CFG) {
  module <- "40"
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少依赖包"))
  if (!md_precheck(CFG, module)) return(md_result(module, "SKIP", "前置条件不满足"))

  vars     <- CFG$env$vars %||% character(0)
  pca_vars <- CFG$env$pca_vars %||% vars
  pca_log  <- isTRUE(CFG$env$pca_log)
  pca_show <- CFG$env$pca_show %||% "sites"
  cor_method <- CFG$env$cor_method %||% "spearman"
  cor_padj   <- CFG$env$cor_p_adjust %||% "BH"
  vif_thr    <- CFG$env$vif_thr %||% 10
  eps        <- CFG$env$pca_eps %||% 1e-4
  group_vars <- CFG$design$factors %||% character(0)
  smp        <- CFG$data$sample_col

  ana <- tryCatch(md_core_analysis(CFG), error = function(e) NULL)
  if (is.null(ana) || nrow(ana) < 3) {
    md_skip(CFG, module, "分析总表不可用或样本过少")
    return(md_result(module, "SKIP", "分析总表不可用"))
  }
  ids <- as.character(ana[[smp]])
  notes <- character(0)

  in_ana   <- vars %in% names(ana)
  num_vars <- vars[in_ana & vapply(ana[vars[in_ana]], is.numeric, logical(1))]
  miss_vars <- vars[!in_ana]
  if (length(miss_vars)) notes <- c(notes, paste0("缺失/非数值环境列已跳过：", paste(miss_vars, collapse = ", ")))
  if (!length(num_vars)) {
    md_skip(CFG, module, "CFG$env$vars 在数据里一个都没有")
    return(md_result(module, "SKIP", "环境变量缺失"))
  }
  # 强制数值型（防 read.csv 把某列读成 character 导致 scale() 报错）
  for (v in num_vars) ana[[v]] <- suppressWarnings(as.numeric(ana[[v]]))
  env_df <- ana[, c(smp, num_vars), drop = FALSE]

  # ---- 1. 描述统计（总体 + 按第一个因子）----
  desc_overall <- do.call(rbind, lapply(num_vars, function(v) {
    x <- ana[[v]]; x <- x[is.finite(x)]
    data.frame(variable = v, group = "Overall", n = length(x),
               mean = mean(x), sd = stats::sd(x),
               se = stats::sd(x) / sqrt(length(x)),
               median = stats::median(x),
               q1 = unname(stats::quantile(x, .25)), q3 = unname(stats::quantile(x, .75)),
               min = min(x), max = max(x),
               cv = stats::sd(x) / mean(x),
               n_missing = sum(!is.finite(ana[[v]])),
               stringsAsFactors = FALSE)
  }))
  desc_by <- NULL
  g1 <- group_vars[1]
  if (!is.null(g1) && g1 %in% names(ana)) {
    desc_by <- do.call(rbind, lapply(num_vars, function(v) {
      do.call(rbind, lapply(sort(unique(as.character(ana[[g1]]))), function(lv) {
        x <- ana[[v]][as.character(ana[[g1]]) == lv]; x <- x[is.finite(x)]
        data.frame(variable = v, group = lv, n = length(x), mean = mean(x),
                           sd = stats::sd(x), se = stats::sd(x) / sqrt(max(length(x), 1)),
                           median = stats::median(x),
                           q1 = unname(stats::quantile(x, .25)),
                           q3 = unname(stats::quantile(x, .75)),
                           min = min(x), max = max(x), cv = stats::sd(x) / mean(x),
                           n_missing = sum(!is.finite(ana[[v]][as.character(ana[[g1]]) == lv])),
                           stringsAsFactors = FALSE)
      }))
    }))
  }
  desc <- rbind(desc_overall, desc_by)

  # ---- 2. VIF（手写 1/(1-R²)）----
  vif <- tryCatch(md_vif(ana, num_vars), error = function(e) NULL)
  if (is.null(vif)) vif <- data.frame(variable = num_vars, vif = NA_real_)
  vif$vif[!is.finite(vif$vif)] <- NA_real_
  vif$flag <- ifelse(is.na(vif$vif), "—", ifelse(vif$vif > vif_thr, "collinear", "ok"))
  vif$threshold <- vif_thr
  keep_collinear <- tryCatch(md_drop_collinear(ana, num_vars, thr = vif_thr),
                             error = function(e) num_vars)
  vif$retained_after_pruning <- vif$variable %in% keep_collinear

  # ---- 3. 相关矩阵（上三角独立配对 BH 校正后镜像）----
  cor_in_vars <- num_vars[!vapply(ana[num_vars], function(x) {
    s <- stats::sd(x, na.rm = TRUE); !is.finite(s) || s <= 0
  }, logical(1))]
  if (length(cor_in_vars) < length(num_vars))
    notes <- c(notes, paste0("零方差列已从相关矩阵剔除：",
                             paste(setdiff(num_vars, cor_in_vars), collapse = ", ")))
  corm <- tryCatch(local_cor_matrix(as.matrix(ana[, cor_in_vars, drop = FALSE]),
                                    method = cor_method, padj = cor_padj),
                   error = function(e) NULL)
  if (is.null(corm)) {
    md_warn(CFG, module, "：相关矩阵计算失败（Hmisc/rcorr）")
    notes <- c(notes, "相关矩阵计算失败。")
  } else {
    notes <- c(notes, sprintf(
      "相关矩阵：%s，仅对上三角 %d 个独立配对做 %s 校正后再镜像（nifH 口径）。",
      cor_method, length(cor_in_vars) * (length(cor_in_vars) - 1) / 2, cor_padj))
  }

  # ---- 4. 环境 PCA（全变量）+ 土壤 PCA（Soil_PC1 唯一定义）----
  pca <- local_pca(ana, pca_vars, log_transform = pca_log, eps = eps)
  if (is.null(pca)) {
    md_warn(CFG, module, "：PCA 不可用（pca_vars 不足或全缺失）")
    notes <- c(notes, "环境 PCA 不可用。")
  } else {
    # 常数列剔除记录
    if (length(pca$dropped)) notes <- c(notes, paste0("PCA 输入剔除零方差列：",
                                                      paste(pca$dropped, collapse = ", ")))
    notes <- c(notes, sprintf(
      "环境 PCA：vars = %s，log = %s（pH %s），n = %d，PC1 = %.2f%%，PC2 = %.2f%%。",
      paste(pca$vars, collapse = "+"), pca_log,
      if (pca_log) "豁免" else "未变换（无需豁免）",
      pca$n_used, 100 * pca$ve[1], 100 * pca$ve[2]))
  }

  # Soil_PC1：CFG$env$pca_vars 与 pca_vars 是同一个键 → 与上面的 pca 是同一组数，
  # 直接复用（不再算第二遍，杜绝"两份 Soil_PC1"的漂移风险）。
  soil <- pca
  if (!is.null(soil)) {
    soil$var_names <- soil$vars
    notes <- c(notes, sprintf(
      "Soil_PC1 定义（全项目唯一）：prcomp(scale(X)) 第一主轴得分，X = CFG$env$pca_vars（%s），log = %s，pH %s。PC1 = %.2f%%。",
      paste(soil$vars, collapse = "+"), pca_log,
      if (pca_log) "豁免 log" else "为原始刻度", 100 * soil$ve[1]))
  } else {
    notes <- c(notes, "土壤 PCA 不可用 → Soil_PC1 全为 NA（50/52/71 会降级）。")
    md_warn(CFG, module, "：土壤 PCA 不可用，Soil_PC1 = NA")
  }

  # ---- 5. derived：sample_id + PC 列（50/52/71 靠 md_derived() 取）----
  derived <- data.frame(sample_id = ids, stringsAsFactors = FALSE)
  if (!is.null(soil)) {
    ncomp <- min(2L, ncol(soil$scores))
    for (k in seq_len(ncomp)) derived[[paste0("Soil_PC", k)]] <- soil$scores[, k]
  } else {
    derived$Soil_PC1 <- NA_real_
  }
  if (!is.null(pca)) {
    ncomp <- min(2L, ncol(pca$scores))
    for (k in seq_len(ncomp)) derived[[paste0("Env_PC", k)]] <- pca$scores[, k]
  }
  # CFG$model$derived 里额外点名的派生变量（目前只支持 PCA 分数）
  want_extra <- setdiff(CFG$model$derived %||% character(0), names(derived))
  if (length(want_extra)) {
    notes <- c(notes, paste0("CFG$model$derived 中点名但本模块未提供的派生变量：",
                             paste(want_extra, collapse = ", ")))
  }

  # ---- 6. 图 ----
  # 画法逐项对照原脚本：
  #   A 土壤 PCA 双标图  <- 10_Supp_env_SEM.R 「2a FigS7a 土壤 PCA 双标图」:
  #        Zone 着色 + Vegetation 定形（16/17/15/18）、点 2.2*2/3 / alpha 0.75 /
  #        stroke 0.5、变量箭头统一缩放（最长箭头 = 较小轴范围的 45%）、
  #        ggrepel 标签、过原点虚线、x/y expand 0.18、图例底部。
  #   B 相关矩阵         <- 10_Supp_env_SEM.R 「2d.3 FigS8B」:
  #        同一批指标顺序、geom_tile 白边 0.25、格内星号、蓝-白-绿同色带（±1）。
  #   C 组均值热图       <- 10_Supp_env_SEM.R 「2d.2 FigS8A」:
  #        列内 z-score、对称 limits ±max|z|、同色带、最上一行 Zone 注释条。
  pt <- md_pt(CFG)
  group_cols <- NULL
  if (length(group_vars) && all(group_vars %in% names(ana))) {
    group_cols <- lapply(group_vars, function(g) as.character(ana[[g]]))
    names(group_cols) <- group_vars
  }

  # A/B 共用的蓝-白-绿色带（原脚本 my_colors；端点各向白混 30%）。
  # 原脚本写字面色值，模块里不做常量表，故这两行按原脚本照抄（属于「画法」）。
  local_blend_white <- function(hex, f = 0.3) {
    col <- grDevices::col2rgb(hex)[, 1] / 255
    grDevices::rgb(col[1] + (1 - col[1]) * f,
                   col[2] + (1 - col[2]) * f,
                   col[3] + (1 - col[3]) * f)
  }
  local_my_colors <- function() {
    grDevices::colorRampPalette(c(local_blend_white("#3A558C"), "white",
                                  local_blend_white("#75D155")))(100)
  }
  # 变量显示顺序（原脚本「气候在前，其后环境变量」；这里是图面画法，不参与任何统计）
  var_order <- unique(c(intersect(c("MAT", "MAP"), num_vars),
                        intersect(c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP",
                                    "EC", "TK", "Ca", "Fe"), num_vars),
                        setdiff(num_vars, c("MAT", "MAP", "SWC", "NH4", "NO3", "pH",
                                            "OM", "TN", "TP", "EC", "TK", "Ca", "Fe"))))
  var_order <- var_order[var_order %in% num_vars]

  # A：土壤 PCA 双标图
  pA <- md_empty_plot(CFG, "PCA unavailable")
  if (!is.null(pca)) {
    # sc$sample_id 不存在（scores 是纯数值矩阵）→ 用行号对齐：
    # 先取完整行，再把 ids 收缩到同一子集，保证逐行对应。
    sc <- as.data.frame(pca$scores)
    cc <- stats::complete.cases(sc[, c("PC1", "PC2")])
    sc <- sc[cc, , drop = FALSE]
    sc_ids <- ids[cc]
    if (nrow(sc) >= 3) {
      sc <- cbind(sc, as.data.frame(lapply(group_cols, function(x) factor(x[cc]))))
      aes_pca <- if (!is.null(group_cols))
        ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data[[names(group_cols)[1]]],
                     shape = if (length(group_cols) > 1) .data[[names(group_cols)[2]]] else NULL)
      else ggplot2::aes(x = .data$PC1, y = .data$PC2)
      pA <- ggplot2::ggplot(sc, aes_pca) +
        ggplot2::geom_point(size = 2.2 * 2 / 3, alpha = 0.75, stroke = 0.5)
      if (!is.null(group_cols)) {
        pA <- pA + ggplot2::scale_colour_manual(
          values = md_pal(CFG, levels(sc[[names(group_cols)[1]]]), names(group_cols)[1]),
          name = names(group_cols)[1])
        if (length(group_cols) > 1)
          pA <- pA + ggplot2::scale_shape_manual(
            values = stats::setNames(c(16, 17, 15, 18, 3, 4)[seq_along(levels(sc[[names(group_cols)[2]]]))],
                                     levels(sc[[names(group_cols)[2]]])),
            name = names(group_cols)[2])
      }
      # 变量箭头（pca_show = species / both）
      if (pca_show %in% c("species", "both") && !is.null(pca$rotation)) {
        ld <- as.matrix(pca$rotation[, seq_len(min(2, ncol(pca$rotation))), drop = FALSE])
        if (ncol(ld) == 1) ld <- cbind(ld, 0)
        rng <- c(diff(range(sc$PC1)), diff(range(sc$PC2)))
        len <- sqrt(ld[, 1]^2 + ld[, 2]^2)
        mx  <- max(len, na.rm = TRUE)
        # 原脚本口径：最长箭头 = 较小轴范围的 45%
        scl <- if (is.finite(mx) && mx > 0 && all(is.finite(rng))) 0.45 * min(rng) / mx else 0
        arr <- data.frame(x = ld[, 1] * scl, y = ld[, 2] * scl,
                          lab = rownames(ld), stringsAsFactors = FALSE)
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          pA <- pA +
            ggplot2::geom_text_repel(data = arr, ggplot2::aes(x = .data$x, y = .data$y,
                                                             label = .data$lab),
                                     inherit.aes = FALSE, size = pt2gg(pt$ann),
                                     colour = "grey15", fontface = "bold", seed = CFG$opts$seed,
                                     min.segment.length = 0, segment.color = "grey55",
                                     segment.size = 0.25, box.padding = 0.35, force = 10,
                                     max.overlaps = Inf)
        } else {
          pA <- pA +
            ggplot2::geom_text(data = arr, ggplot2::aes(x = .data$x, y = .data$y,
                                                        label = .data$lab),
                               inherit.aes = FALSE, size = pt2gg(pt$ann),
                               colour = "grey15", fontface = "bold", vjust = -0.2)
        }
        pA <- pA +
          ggplot2::geom_segment(data = arr, ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
                                inherit.aes = FALSE, colour = "grey25",
                                linewidth = MD_LW$axis,
                                arrow = ggplot2::arrow(length = grid::unit(0.16, "cm"), type = "closed"))
        notes <- c(notes, sprintf("PCA 双标图显示模式 pca_show = %s（含变量箭头）。", pca_show))
      } else {
        notes <- c(notes, sprintf("PCA 双标图显示模式 pca_show = %s（只出样方点）。", pca_show))
      }
      # 过原点的参考线（原脚本画在散点之后，故这里同样置于边界线图层之后）
      pA <- pA +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70",
                            linewidth = MD_LW$ref) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70",
                            linewidth = MD_LW$ref) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.18)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.18)) +
        ggplot2::labs(x = sprintf("PC1 (%.1f%%)", 100 * pca$ve[1]),
                      y = sprintf("PC2 (%.1f%%)", 100 * pca$ve[2]), tag = "A") +
        theme_md(CFG = CFG) +
        ggplot2::theme(legend.position = "bottom",
                       legend.box = "vertical",
                       legend.box.just = "center",
                       legend.margin = ggplot2::margin(0, 0, 0, 0))
    }
  }

  # B：相关矩阵热图（上三角校正后镜像）
  pB <- md_empty_plot(CFG, "Correlation unavailable")
  if (!is.null(corm)) {
    k <- ncol(corm$r)
    lv <- var_order[var_order %in% colnames(corm$r)]
    if (length(lv) < 2) lv <- colnames(corm$r)
    df_cor <- expand.grid(Var = lv, Indicator = lv, stringsAsFactors = FALSE)
    df_cor$r <- corm$r[cbind(df_cor$Indicator, df_cor$Var)]
    df_cor$s <- corm$stars[cbind(df_cor$Indicator, df_cor$Var)]
    df_cor$Var <- factor(df_cor$Var, levels = lv)
    df_cor$Indicator <- factor(df_cor$Indicator, levels = rev(lv))
    pB <- ggplot2::ggplot(df_cor, ggplot2::aes(x = .data$Var, y = .data$Indicator)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$r), colour = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = .data$s), size = pt2gg(4.2), colour = "grey10") +
      ggplot2::scale_fill_gradientn(colours = local_my_colors(), name = "Spearman r",
                                    limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                                    guide = ggplot2::guide_colourbar(
                                      barwidth = grid::unit(30, "mm"),
                                      barheight = grid::unit(2.4, "mm"))) +
      ggplot2::scale_x_discrete(limits = lv) +
      ggplot2::scale_y_discrete(limits = rev(lv)) +
      ggplot2::labs(x = NULL, y = NULL, tag = "B") +
      theme_md(CFG = CFG) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = pt$tick),
                     axis.text.y = ggplot2::element_text(size = pt$tick),
                     axis.ticks = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     legend.position = "bottom",
                     legend.title = ggplot2::element_text(size = MD_PT_COMP$legend),
                     legend.text = ggplot2::element_text(size = MD_PT_COMP$legend))
  }

  # C：组均值热图（各组均值的列内 z-score + 最上一行 Zone 注释条）
  pC <- md_empty_plot(CFG, "Group means unavailable")
  gm <- NULL
  if (length(group_vars) >= 2) {
    gm <- local_group_means(ana, group_vars[1:2], num_vars, CFG$design$ordered_levels)
  } else if (length(group_vars) == 1) {
    gm <- local_group_means(ana, group_vars[1], num_vars, CFG$design$ordered_levels)
  }
  if (!is.null(gm)) {
    ind_levels <- colnames(gm$z)
    ind_use <- var_order[var_order %in% ind_levels]
    ind_use <- if (length(ind_use) >= 2) ind_use else ind_levels
    # y 轴自上而下 = Zone 注释条 + 全部指标（与 B 面板同序，便于上下对齐）
    ylev <- c("Zone", ind_use)
    df_gm <- expand.grid(group = gm$levels, variable = ind_use, stringsAsFactors = FALSE)
    df_gm$z <- as.vector(gm$z[cbind(match(df_gm$group, gm$levels),
                                    match(df_gm$variable, colnames(gm$z)))])
    df_gm$group <- factor(df_gm$group, levels = gm$levels)
    df_gm$variable <- factor(df_gm$variable, levels = ylev)
    zlim <- max(abs(df_gm$z), na.rm = TRUE)
    if (!is.finite(zlim) || zlim == 0) zlim <- 1
    pC <- ggplot2::ggplot(df_gm, ggplot2::aes(x = .data$group, y = .data$variable)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$z), colour = "white", linewidth = 0.25)
    # Zone 注释条：每个组水平取 Zone 色（md_pal() 与 cols_zone 同源）
    if (length(group_vars) >= 2 && group_vars[1] %in% names(ana)) {
      zone_of <- vapply(gm$levels, function(g) {
        i <- which(as.character(ana[[group_vars[1]]]) == sub("_.*$", "", g))
        if (length(i)) as.character(ana[[group_vars[1]]])[i[1]] else NA_character_
      }, character(1))
      ann_df <- data.frame(Group = factor(gm$levels, levels = gm$levels),
                           Zone = factor(zone_of, levels = levels(factor(ana[[group_vars[1]]]))),
                           stringsAsFactors = FALSE)
      ann_df <- ann_df[!is.na(ann_df$Zone), , drop = FALSE]
      if (nrow(ann_df)) {
        # Zone 注释条用离散色标叠在连续色标之上 —— 必须开一个独立的 fill scale
        # （原脚本 04/10 用 cols_zone 直接顶掉主 fill，本包用 ggnewscale 保持两个色标并存；
        #  缺 ggnewscale 时退化为只画 Z-score 热图，不静默画错）
        if (requireNamespace("ggnewscale", quietly = TRUE)) {
          zc <- md_pal(CFG, levels(ann_df$Zone), group_vars[1])
          pC <- pC +
            ggplot2::scale_fill_gradientn(colours = local_my_colors(), name = "Z-score",
                                          limits = c(-zlim, zlim),
                                          breaks = scales::pretty_breaks(n = 5),
                                          guide = ggplot2::guide_colourbar(
                                            barwidth = grid::unit(30, "mm"),
                                            barheight = grid::unit(2.4, "mm"))) +
            ggnewscale::new_scale_fill() +
            ggplot2::geom_tile(data = ann_df, ggplot2::aes(x = .data$Group, y = "Zone",
                                                          fill = .data$Zone),
                               inherit.aes = FALSE, colour = "white", linewidth = 0.25) +
            ggplot2::scale_fill_manual(values = zc, name = group_vars[1],
                                       guide = ggplot2::guide_legend(
                                         override.aes = list(colour = NA)))
        } else {
          notes <- c(notes, "未装 ggnewscale -> 组均值热图未画 Zone 注释条。")
          pC <- pC +
            ggplot2::scale_fill_gradientn(colours = local_my_colors(), name = "Z-score",
                                          limits = c(-zlim, zlim),
                                          breaks = scales::pretty_breaks(n = 5),
                                          guide = ggplot2::guide_colourbar(
                                            barwidth = grid::unit(30, "mm"),
                                            barheight = grid::unit(2.4, "mm")))
        }
      } else {
        pC <- pC +
          ggplot2::scale_fill_gradientn(colours = local_my_colors(), name = "Z-score",
                                        limits = c(-zlim, zlim),
                                        breaks = scales::pretty_breaks(n = 5),
                                        guide = ggplot2::guide_colourbar(
                                          barwidth = grid::unit(30, "mm"),
                                          barheight = grid::unit(2.4, "mm")))
      }
    } else {
      pC <- pC +
        ggplot2::scale_fill_gradientn(colours = local_my_colors(), name = "Z-score",
                                      limits = c(-zlim, zlim),
                                      breaks = scales::pretty_breaks(n = 5),
                                      guide = ggplot2::guide_colourbar(
                                        barwidth = grid::unit(30, "mm"),
                                        barheight = grid::unit(2.4, "mm")))
    }
    pC <- pC +
      ggplot2::scale_y_discrete(limits = ylev,
                                labels = c(Zone = "",
                                           stats::setNames(ind_use, ind_use))) +
      ggplot2::labs(x = NULL, y = NULL, tag = "C") +
      theme_md(CFG = CFG) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                         size = MD_PT_COMP$tick),
                     axis.text.y = ggplot2::element_text(size = pt$tick),
                     axis.ticks = ggplot2::element_blank(),
                     legend.position = "bottom",
                     legend.title = ggplot2::element_text(size = MD_PT_COMP$legend),
                     legend.text = ggplot2::element_text(size = MD_PT_COMP$legend),
                     panel.border = ggplot2::element_blank())
  }

  # ---- 6b. 随机森林变量重要性（**可选面板，默认关**；env$rf$enabled = TRUE 才出）----
  #   源：毕业论文 3_environmental_drivers.R:1199-1258（ranger + permutation importance，
  #       相对重要性缩放到 max = 1，单色单柱、无 dodged、无 model 图例）
  #   通用化：因变量与预测变量都来自 config，**不限于"环境 -> 群落"**。
  #   ⚠ 设计串 "ADD\nBCC" 里的 **D 区本来就是预留的空位** —— 打开本开关正好把它填满，
  #     不会挤动 A / B / C 的位置（默认关闭时图件字节与之前完全一致）。
  rf_cfg <- CFG$env$rf %||% list()
  p_rf <- NULL; tab_rf <- NULL; pdp_plot <- NULL; rf_notes <- character(0)
  if (isTRUE(rf_cfg$enabled)) {
    yv <- as.character(rf_cfg$response %||% "")
    # 因变量不在分析总表里时，先试"是不是 α 指数"（Shannon / Simpson / Observed / Pielou）——
    # 这些由模块 11 算过、但没有并进 md_core_analysis()。能现算就补上，算不出才跳过。
    if (nzchar(yv) && !yv %in% names(ana)) {
      av <- tryCatch(local_rf_alpha(CFG, yv), error = function(e) NULL)
      if (!is.null(av) && yv %in% names(av)) {
        ana <- merge(ana, av, by = CFG$data$sample_col, all.x = TRUE)
        rf_notes <- c(rf_notes, paste0("rf: 因变量 '", yv,
                                       "' 不在分析总表里，已由 md_core_asv(\"alpha\") 现算补齐。"))
      }
    }
    xs <- as.character(rf_cfg$predictors %||% num_vars)
    xs <- setdiff(intersect(xs, names(ana)), yv)
    if (!nzchar(yv) || !yv %in% names(ana)) {
      rf_notes <- c(rf_notes, paste0("rf: 因变量 '", yv, "' 不在分析总表里，跳过该面板"))
    } else if (length(xs) < 2) {
      rf_notes <- c(rf_notes, "rf: 可用预测变量少于 2 个，跳过该面板")
    } else {
      fit_rf <- tryCatch(
        md_rf_fit(CFG, ana, yv, xs,
                  engine = as.character(rf_cfg$engine %||% "ranger"),
                  num_trees = as.integer(rf_cfg$num_trees %||% 2000L),
                  importance = as.character(rf_cfg$importance %||% "permutation")),
        error = function(e) { rf_notes <<- c(rf_notes, paste0("rf: ", conditionMessage(e))); NULL })
      if (is.null(fit_rf)) {
        rf_notes <- c(rf_notes, "rf: 拟合失败或所依赖的包不可用，跳过该面板")
      } else {
        tab_rf <- fit_rf$importance
        tab_rf$response <- yv
        tab_rf$engine   <- fit_rf$engine
        tab_rf$R2       <- fit_rf$r2
        tab_rf$n        <- fit_rf$n
        st_rf <- rf_cfg$style %||% list()
        p_rf <- md_rf_importance_plot(
          fit_rf$importance, top_n = as.integer(rf_cfg$top_n %||% 10L),
          fill = st_rf$fill %||% "#3C5488", r2 = fit_rf$r2,
          label_r2 = isTRUE(st_rf$label_r2 %||% TRUE),
          xlab = "Relative permutation importance",
          pt = md_pt(CFG, comp = TRUE), tag = "D")
        rf_notes <- c(rf_notes, paste0("rf: engine=", fit_rf$engine, ", n=", fit_rf$n,
                                       ", R2=", md_fmt_num(fit_rf$r2, 3),
                                       "（面板 D）"))
        if (!is.null(fit_rf$pdp) && nrow(fit_rf$pdp) > 0) {
          pd <- fit_rf$pdp
          xcol <- setdiff(names(pd), c("yhat", "Variable"))
          if (length(xcol)) {
            pd$.x <- pd[[xcol[1]]]
            pdp_plot <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$.x, y = .data$yhat)) +
              ggplot2::geom_line(colour = MD_REG$col, linewidth = MD_REG$lw * 2) +
              ggplot2::facet_wrap(~ Variable, scales = "free_x", nrow = 1) +
              ggplot2::labs(x = NULL, y = paste0("Partial dependence of ", yv), tag = "E") +
              theme_md(CFG = NULL)
          }
        }
      }
    }
  }

  all_panels <- list(A = pA, B = pB, C = pC)
  des <- "ADD\nBCC"
  if (!is.null(p_rf)) all_panels$D <- p_rf
  fig_h <- 190
  if (!is.null(pdp_plot)) {
    all_panels$E <- pdp_plot
    des <- "ADD\nBCC\nEEE"
    fig_h <- 245
  }
  fig <- md_compose(CFG, all_panels, design = des, guides = "keep", legend = NULL)
  md_save_figure(CFG, fig, "40_env_summary", height = fig_h)

  # ---- 7. 写 cache（跨模块契约，先写缓存再写表）----
  cache_obj <- list(
    env      = env_df,
    desc     = desc,
    vif      = vif,
    cor      = corm,
    pca      = pca,
    soil_pca = soil,
    derived  = derived,
    group_means = gm,
    params   = list(vars = num_vars, pca_vars = pca_vars, pca_log = pca_log,
                    pca_show = pca_show, cor_method = cor_method,
                    cor_p_adjust = cor_padj, vif_thr = vif_thr,
                    n_samples = nrow(ana), seed = CFG$opts$seed),
    notes    = notes)
  saveRDS(cache_obj, md_path(CFG, "cache", "40_env.rds"))
  md_log(CFG, "OK", "写出 40_env.rds（derived: ", paste(names(derived), collapse = ", "), "）")

  # ---- 8. 写表 ----
  soil_tab <- if (is.null(soil)) NULL else data.frame(
    Variable = colnames(soil$rotation),
    Loading_PC1 = soil$rotation[, 1],
    Loading_PC2 = if (ncol(soil$rotation) >= 2) soil$rotation[, 2] else NA_real_,
    stringsAsFactors = FALSE)
  env_pca_tab <- if (is.null(pca)) NULL else data.frame(
    Variable = colnames(pca$rotation),
    Loading_PC1 = pca$rotation[, 1],
    Loading_PC2 = if (ncol(pca$rotation) >= 2) pca$rotation[, 2] else NA_real_,
    stringsAsFactors = FALSE)

  note_df <- data.frame(Note = c(
    sprintf("环境变量（CFG$env$vars）：%s。", paste(num_vars, collapse = ", ")),
    sprintf("VIF = 1/(1-R²)（共享 md_vif()，与 nifH 手写口径一致）；阈值 %s，超阈值标记 collinear，retained_after_pruning 为逐步剔除后保留的集合。", vif_thr),
    notes,
    "Soil_PC1 是本包 50/52/71 唯一的来源（md_derived(CFG, \"Soil_PC1\")），其它模块不得自行重算。",
    "组均值热图的值 = 各组均值的**列内 z-score**（nifH FigS8A 口径）。"))

  tp <- md_write_xlsx(CFG, "40_env_summary", list(
    Descriptive = desc,
    VIF = vif,
    Correlations = if (is.null(corm)) NULL else corm$long,
    Cor_matrix = if (is.null(corm)) NULL else data.frame(
      variable = rownames(corm$r), as.data.frame(corm$r), check.names = FALSE,
      stringsAsFactors = FALSE),
    P_adjusted = if (is.null(corm)) NULL else data.frame(
      variable = rownames(corm$p_adj), as.data.frame(corm$p_adj), check.names = FALSE,
      stringsAsFactors = FALSE),
    Soil_PCA = soil_tab,
    Env_PCA = env_pca_tab,
    Group_means = if (is.null(gm)) NULL else data.frame(
      group = gm$levels, n = gm$n, as.data.frame(gm$means), check.names = FALSE,
      stringsAsFactors = FALSE),
    Derived = derived,
    RF_importance = tab_rf,
    Note = rbind(note_df, if (length(rf_notes)) data.frame(Note = rf_notes) else NULL)),
    meta_extra = list(n_samples = nrow(ana), pca_log = pca_log,
                      soil_pc1_ve = if (is.null(soil)) NA_real_ else soil$ve[1]))

  md_result(module, "OK", tables = basename(tp),
            figures = "figures/40_env_summary.pdf",
            notes = paste(c(notes, rf_notes), collapse = " | "))
}
