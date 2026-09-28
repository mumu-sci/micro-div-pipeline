# =============================================================================
# R/15_beta_tests.R ---- 多元组间差异检验（PERMANOVA / ANOSIM / PERMDISP / MRPP）
# -----------------------------------------------------------------------------
# 依赖：13_distances（cache/13_dist.rds），距离一律经 md_dist_get(CFG, metric) 取。
#
# 口径要点（nifH 原码踩过的坑，逐条保留）：
#   1. PERMANOVA 必须 by = "terms"：本机 vegan 的 by = "margin" 只返回末项 + Residual + Total
#      （02_Fig02.R §D 注释实测）。交互模型因此用 by = "terms" 拿全部项。
#   2. pairwise adonis2 = combn(sort(levels), 2) 逐对子集重跑；任一组 n < min_group_n
#      记 NA 并在 Note 列写 "n<min_group_n"（08_Supp_diversity.R::pairwise_adonis2 口径）。
#      p 校正在**全部对**上做（p.adjust 对 NA 返回 NA，不影响其余值）。
#   3. PERMDISP：betadisper + permutest，P = (#{perm >= stat} + 1) / (n + 1)
#      （02_Fig02.R 手写口径，而非 permutest 的打印值）。
#   4. 置换前 md_seed(CFG, offset)；不同检验/不同对用不同 offset，保证同轮可复现。
#   5. 每个因子的 PERMDISP 单独拟合（各自取该因子可用水平），不跨因子混合。
#      ⚠ 因此 A 面板里并列的 PERMDISP 点来自**不同的 betadisper 拟合**，不能横向当作同一尺度比较
#        （已在 Note sheet 写明）。
#
# 输出：figures/15_beta_tests.pdf（一张 patchwork；PERMDISP 不可估计时自动只出 A/B 两面板）
#       results/15_beta_tests.xlsx（Effect_sizes / Pairwise / PERMDISP / Note）
# =============================================================================

md_beta_tests <- function(CFG) {
  module <- "15"
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少 vegan/ape/picante/betapart"))
  if (!md_precheck(CFG, module)) return(md_result(module, "SKIP", "前置条件不满足"))

  bt      <- CFG$beta_tests %||% list()
  nperm   <- as.integer(bt$permutations %||% 999)
  padj_m  <- bt$p_adjust %||% "fdr"
  min_gn  <- as.integer(bt$min_group_n %||% 3)
  do_anosim <- isTRUE(bt$anosim %||% TRUE)
  do_permd  <- isTRUE(bt$permdisp %||% TRUE)
  do_mrpp   <- isTRUE(bt$mrpp %||% FALSE)
  do_pw     <- isTRUE(bt$pairwise %||% TRUE)
  do_terms  <- isTRUE(bt$permanova_terms %||% TRUE)

  meta_all <- md_core_meta(CFG)
  smp_col  <- CFG$data$sample_col
  factors  <- CFG$design$factors %||% character(0)
  factors  <- factors[factors %in% names(meta_all)]
  if (length(factors) == 0) {
    md_skip(CFG, module, "design$factors 为空或不在元数据中")
    return(md_result(module, "SKIP", "无可用分组因子"))
  }

  # ---- 距离：只用 13 的缓存 ----
  metrics <- CFG$distances$metrics %||% "bray"
  d_cache <- tryCatch(md_dist_cache(CFG), error = function(e) NULL)
  if (is.null(d_cache)) {
    md_skip(CFG, module, "cache/13_dist.rds 不存在（请先跑模块 13）")
    return(md_result(module, "SKIP", "缺少 13 的距离缓存"))
  }
  dlist <- list()
  for (mt in metrics) {
    d <- md_dist_get(CFG, mt)
    if (!is.null(d) && inherits(d, "dist") && length(d) > 0) dlist[[mt]] <- d
  }
  if (length(dlist) == 0) {
    md_skip(CFG, module, "13 的缓存里没有可用 dist（metrics: ", paste(metrics, collapse = ", "), "）")
    return(md_result(module, "SKIP", "无可用的距离矩阵"))
  }
  metrics_ok <- names(dlist)

  d_main <- CFG$response_models$pair_metric %||% "bray"
  if (!d_main %in% metrics_ok) d_main <- metrics_ok[1]

  levels_of <- function(f) {
    lv  <- CFG$design$ordered_levels[[f]] %||% NULL
    obs <- unique(as.character(meta_all[[f]]))
    obs <- obs[!is.na(obs) & nzchar(obs)]
    if (is.null(lv)) return(sort(obs))
    c(intersect(lv, obs), setdiff(sort(obs), lv))
  }
  levs <- stats::setNames(lapply(factors, levels_of), factors)

  # meta 行序对齐到距离 labels（13 的契约是列序一致，这里再显式对齐一次更稳）
  align_meta <- function(ids, meta) {
    m <- meta[match(ids, meta[[smp_col]]), , drop = FALSE]
    rownames(m) <- ids
    m
  }

  eff_rows <- list(); pw_rows <- list(); disp_rows <- list()
  notes <- character(0)
  .off <- 0L

  permanova_one <- function(d, m, rhs, metric, model) {
    .off <<- .off + 1L
    md_seed(CFG, 500 + .off)
    form <- stats::as.formula(paste0("d ~ ", rhs))
    r <- tryCatch(vegan::adonis2(form, data = m, permutations = nperm, by = "terms"),
                  error = function(e) e)
    if (inherits(r, "error")) {
      notes <<- c(notes, paste0("PERMANOVA [", metric, "] ", rhs, " 失败：", conditionMessage(r)))
      return(NULL)
    }
    df <- as.data.frame(r, stringsAsFactors = FALSE)
    df$Term  <- rownames(df)
    df$Metric <- metric
    df$Model  <- model
    df$n      <- nrow(m)
    df$Permutations <- nperm
    df
  }

  for (mt in metrics_ok) {
    d   <- dlist[[mt]]
    ids <- attr(d, "Labels") %||% labels(d)
    m   <- align_meta(ids, meta_all)
    ok  <- stats::complete.cases(m[, factors, drop = FALSE])
    if (any(!ok)) {
      d   <- md_subset_dist(d, ids[ok])
      ids <- attr(d, "Labels") %||% labels(d)
      m   <- m[ok, , drop = FALSE]
      notes <- c(notes, paste0("距离 ", mt, "：", sum(!ok),
                               " 个样本的因子含缺失，已从该距离的全部检验中剔除"))
    }
    fac_use <- factors[vapply(factors, function(f) {
      t <- table(droplevels(factor(m[[f]])))
      length(t) >= 2 && min(t) >= 2
    }, logical(1))]
    dropped <- setdiff(factors, fac_use)
    if (length(dropped))
      notes <- c(notes, paste0("距离 ", mt, "：因子 ", paste(dropped, collapse = ", "),
                               " 水平数 < 2 或存在 n < 2 的水平，未做 PERMANOVA"))
    if (length(fac_use) == 0) {
      notes <- c(notes, paste0("距离 ", mt, "：没有任何因子满足最小样本量，跳过"))
      next
    }

    # (a) 各因子单独模型
    for (f in fac_use) {
      r <- permanova_one(d, m, f, mt, "single")
      if (!is.null(r) && nrow(r)) eff_rows[[length(eff_rows) + 1L]] <- r
    }
    # (b) 全域模型 ~ f1 + f2（各项 = 边际贡献）
    if (length(fac_use) >= 2) {
      r <- permanova_one(d, m, paste(fac_use, collapse = " + "), mt, "single_all")
      if (!is.null(r) && nrow(r)) eff_rows[[length(eff_rows) + 1L]] <- r
    }
    # (c) 交互模型 ~ f1 * f2，by = "terms"（**绝不用 by = "margin"**）
    if (do_terms) {
      rhs <- paste(fac_use, collapse = " * ")
      r <- permanova_one(d, m, rhs, mt, "interaction")
      if (!is.null(r) && nrow(r)) eff_rows[[length(eff_rows) + 1L]] <- r
    }

    # ---- 全局 ANOSIM ----
    if (do_anosim) {
      for (f in fac_use) {
        .off <<- .off + 1L
        md_seed(CFG, 600 + .off)
        a <- tryCatch(vegan::anosim(d, grouping = m[[f]], permutations = nperm),
                      error = function(e) NULL)
        if (is.null(a)) next
        eff_rows[[length(eff_rows) + 1L]] <- data.frame(
          Term = paste0("ANOSIM(", f, ")"), Df = NA_integer_, SumOfSqs = NA_real_,
          R2 = unname(a$statistic), F = NA_real_, `Pr(>F)` = a$signif,
          Metric = mt, Model = "anosim", n = nrow(m), Permutations = nperm,
          check.names = FALSE, stringsAsFactors = FALSE)
      }
    }

    # ---- 全局 MRPP（可选）----
    if (do_mrpp) {
      for (f in fac_use) {
        .off <<- .off + 1L
        md_seed(CFG, 700 + .off)
        mp <- tryCatch(vegan::mrpp(d, grouping = m[[f]], permutations = nperm),
                       error = function(e) NULL)
        if (is.null(mp)) next
        eff_rows[[length(eff_rows) + 1L]] <- data.frame(
          Term = paste0("MRPP(", f, ")"), Df = NA_integer_, SumOfSqs = unname(mp$delta),
          R2 = unname(mp$A), F = NA_real_, `Pr(>F)` = mp$Pvalue,
          Metric = mt, Model = "mrpp", n = nrow(m), Permutations = nperm,
          check.names = FALSE, stringsAsFactors = FALSE)
      }
    }

    # ---- PERMDISP（每因子单独拟合）----
    if (do_permd) {
      for (f in fac_use) {
        bd <- local_betadisper(CFG, d, m[[f]], min_gn)
        if (is.null(bd)) {
          notes <- c(notes, paste0("PERMDISP [", mt, "] ", f,
                                   "：满足 n >= ", min_gn, " 的水平不足 2 个，跳过"))
          next
        }
        disp_rows[[length(disp_rows) + 1L]] <- cbind(
          bd$test,
          data.frame(Metric = mt, Factor = f, stringsAsFactors = FALSE))
      }
    }

    # ---- 两两 PERMANOVA + 两两 ANOSIM ----
    if (do_pw) {
      for (f in fac_use) {
        x <- local_pairwise_beta(CFG, d, m, f, mt, nperm, padj_m, min_gn, do_anosim)
        if (!is.null(x) && nrow(x)) pw_rows[[length(pw_rows) + 1L]] <- x
      }
    }
  }

  eff_df  <- if (length(eff_rows)) do.call(rbind, eff_rows) else
    data.frame(Term = character(0), Metric = character(0), stringsAsFactors = FALSE)
  pw_df   <- if (length(pw_rows))  do.call(rbind, pw_rows)  else
    data.frame(Metric = character(0), Factor = character(0), stringsAsFactors = FALSE)
  disp_df <- if (length(disp_rows)) do.call(rbind, disp_rows) else
    data.frame(Factor = character(0), Metric = character(0), stringsAsFactors = FALSE)

  # ---- p 校正（按 距离 x 模型 分组，避免把不同模型混在一起校正）----
  if (nrow(eff_df) && "Pr(>F)" %in% names(eff_df)) {
    eff_df$p_adjusted <- NA_real_
    key <- paste(eff_df$Metric, eff_df$Model, sep = "|")
    for (k in unique(key)) {
      i <- which(key == k & !is.na(eff_df$`Pr(>F)`))
      if (length(i)) eff_df$p_adjusted[i] <- md_p_adjust(eff_df$`Pr(>F)`[i], padj_m)
    }
    eff_df$p_signif <- local_sig(eff_df$p_adjusted)
  }
  if (nrow(disp_df) && "P" %in% names(disp_df)) {
    disp_df$p_adjusted <- md_p_adjust(disp_df$P, padj_m)
    disp_df$p_signif   <- local_sig(disp_df$p_adjusted)
  }

  # ---- 图 ----
  d_main_group <- align_meta(attr(dlist[[d_main]], "Labels") %||% labels(dlist[[d_main]]), meta_all)
  gA <- local_fig_effect(eff_df)
  gB <- local_fig_disp(disp_df, min_gn)
  plots <- list(A = gA)
  if (!inherits(gB, "gg")) gB <- md_empty_plot(CFG, "PERMDISP not estimable")
  plots$B <- gB
  design <- "AA\nBB"
  fig <- tryCatch(md_compose(CFG, plots, design = design), error = function(e) NULL)
  if (is.null(fig)) fig <- md_compose(CFG, list(A = md_empty_plot(CFG, "no estimable result")), design = "A")
  md_save_figure(CFG, fig, "15_beta_tests", height = 150)

  # ---- 表 ----
  note_df <- data.frame(
    Item = c("distance_metrics", "primary_distance", "permutations", "p_adjust",
             "min_group_n", "permanova_by", "permanova_models",
             "permdisp_P", "permdisp_independence", "anosim", "mrpp", "empty_cells",
             "Notes"),
    Value = c(
      paste(metrics_ok, collapse = ", "),
      d_main,
      as.character(nperm),
      as.character(padj_m),
      as.character(min_gn),
      "terms",
      paste(c("single", "single_all (~ f1 + f2)", "interaction (~ f1 * f2)")[
        c(TRUE, length(factors) >= 2, do_terms)], collapse = "; "),
      "P = (#{perm >= statistic} + 1) / (n_perm + 1)（与 nifH 02_Fig02.R 手写口径一致）",
      "每个因子的 PERMDISP 各自拟合 betadisper（各自取该因子 n >= min_group_n 的水平）；并列展示不等于同一尺度可比",
      if (do_anosim) "vegan::anosim (R statistic)" else "关闭",
      if (do_mrpp) "vegan::mrpp (delta, A)" else "关闭",
      "组内样本不足 min_group_n 的两两比较一律填 NA，并在 Note 列注明（不删行）",
      if (length(notes)) paste(unique(notes), collapse = " | ") else "无"),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "15_beta_tests", list(
    Effect_sizes = eff_df,
    Pairwise = pw_df,
    PERMDISP = disp_df,
    Note = note_df),
    meta_extra = list(distance_primary = d_main, permanova_by = "terms",
                      permutations = nperm, min_group_n = min_gn))
  if (is.null(tp)) return(md_result(module, "PARTIAL", "写表失败（缺 openxlsx？）"))

  md_log(CFG, "INFO", "15 完成：效应量 ", nrow(eff_df), " 行 / 两两 ", nrow(pw_df),
         " 行 / PERMDISP ", nrow(disp_df), " 行")
  md_result(module, "OK", tables = basename(tp), figures = "figures/15_beta_tests.pdf",
            notes = notes)
}

# =============================================================================
# 内部工具（local_ 前缀，不进共享库）
# =============================================================================

local_sig <- function(p) {
  ifelse(is.na(p), "--",
         ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                         ifelse(p < 0.05, "*", "ns"))))
}

#' betadisper + permutest，P 用手写口径，并返回每组离散度摘要。
#' 返回 list(test = 1 行 data.frame(F, P, n_perm, n_groups, n, excl_groups),
#'           group = data.frame(Factor, Group, n, mean, sd, median))
local_betadisper <- function(CFG, d, g, min_group_n = 3) {
  g <- as.character(g)
  tab <- table(g[!is.na(g)])
  if (length(tab) < 2) return(NULL)
  ok_g <- names(tab)[tab >= min_group_n]
  if (length(ok_g) < 2) return(NULL)
  ids <- attr(d, "Labels") %||% labels(d)
  idx <- which(g %in% ok_g)
  d2 <- md_subset_dist(d, ids[idx])
  g2 <- droplevels(factor(g[idx], levels = intersect(names(tab), ok_g)))
  md_seed(CFG, 810 + length(unique(g2)) * 7L)
  bd <- tryCatch(vegan::betadisper(d2, group = g2, type = "median"), error = function(e) NULL)
  if (is.null(bd)) return(NULL)
  pt <- tryCatch(vegan::permutest(bd, permutations = 999), error = function(e) NULL)
  if (is.null(pt)) return(NULL)
  perm <- pt$perm %||% numeric(0)
  stat <- unname(pt$statistic)
  pval <- if (length(perm)) (sum(perm >= stat, na.rm = TRUE) + 1) / (length(perm) + 1) else NA_real_
  dists <- bd$distances
  gg <- as.character(bd$group)
  lv <- levels(bd$group)
  grp <- do.call(rbind, lapply(lv, function(l) {
    x <- dists[gg == l]
    data.frame(Group = l, n = length(x), mean = mean(x), sd = stats::sd(x),
               median = stats::median(x), stringsAsFactors = FALSE)
  }))
  list(test = data.frame(F = stat, P = pval, n_perm = length(perm),
                         n_groups = length(lv), n = length(g2),
                         excl_groups = paste(setdiff(names(tab), ok_g), collapse = ", "),
                         stringsAsFactors = FALSE),
       group = grp)
}

#' 两两 PERMANOVA / ANOSIM（combn + FDR；组内 n < min_group_n 记 NA，行保留）
local_pairwise_beta <- function(CFG, d, m, f, metric, nperm, padj_m, min_group_n,
                                do_anosim = TRUE) {
  mtx <- as.matrix(d)
  ids <- rownames(mtx)
  g   <- as.character(m[[f]][match(ids, rownames(m))])
  lv_ok <- sort(unique(g[!is.na(g) & nzchar(g)]))
  if (length(lv_ok) < 2) return(NULL)
  ord <- CFG$design$ordered_levels[[f]] %||% character(0)
  lv_ord <- c(intersect(ord, lv_ok), setdiff(lv_ok, ord))
  cmb <- utils::combn(lv_ord, 2, simplify = FALSE)
  out <- lapply(cmb, function(pr) {
    idx <- which(g %in% pr)
    n1 <- sum(g[idx] == pr[1]); n2 <- sum(g[idx] == pr[2])
    base <- data.frame(Metric = metric, Factor = f, Group1 = pr[1], Group2 = pr[2],
                       n1 = n1, n2 = n2, stringsAsFactors = FALSE)
    if (n1 < min_group_n || n2 < min_group_n || (n1 + n2) < 4) {
      return(cbind(base, data.frame(R2 = NA_real_, F = NA_real_, p_permanova = NA_real_,
                                    R_anosim = NA_real_, p_anosim = NA_real_,
                                    Note = paste0("n<", min_group_n), stringsAsFactors = FALSE)))
    }
    sid <- ids[idx]
    sub_d <- stats::as.dist(mtx[sid, sid, drop = FALSE])
    sub_g <- droplevels(factor(g[idx], levels = pr))
    .o <- (n1 * 7919L + n2 * 104729L + nchar(metric) * 31L) %% 100000L
    md_seed(CFG, 900 + .o)
    a <- tryCatch(vegan::adonis2(sub_d ~ sub_g, permutations = nperm, by = "terms"),
                  error = function(e) NULL)
    R2 <- Fs <- pp <- NA_real_
    if (!is.null(a)) {
      ad <- as.data.frame(a)
      i <- which(rownames(ad) == "sub_g")[1]
      if (length(i) && !is.na(i)) { R2 <- ad$R2[i]; Fs <- ad$F[i]; pp <- ad$`Pr(>F)`[i] }
    }
    Ra <- pa <- NA_real_
    if (do_anosim) {
      md_seed(CFG, 950 + .o)
      an <- tryCatch(vegan::anosim(sub_d, grouping = sub_g, permutations = nperm),
                     error = function(e) NULL)
      if (!is.null(an)) { Ra <- unname(an$statistic); pa <- an$signif }
    }
    cbind(base, data.frame(R2 = R2, F = Fs, p_permanova = pp,
                           R_anosim = Ra, p_anosim = pa, Note = "",
                           stringsAsFactors = FALSE))
  })
  df <- do.call(rbind, out)
  df$p_permanova_adjusted <- md_p_adjust(df$p_permanova, padj_m)
  df$p_anosim_adjusted    <- md_p_adjust(df$p_anosim, padj_m)
  df
}

#' 面板 A：效应量条形图（PERMANOVA R2 / ANOSIM R / MRPP A），按距离 x 模型 x 项
local_fig_effect <- function(eff_df) {
  if (is.null(eff_df) || !nrow(eff_df) || !"R2" %in% names(eff_df))
    return(md_empty_plot(NULL, "no PERMANOVA / ANOSIM result"))
  d <- eff_df[!is.na(eff_df$R2) & !eff_df$Term %in% c("Total", "Residual"), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(NULL, "no estimable term"))
  d$Lab <- paste0(d$Metric, " | ", d$Model, " | ", d$Term)
  d$Lab <- factor(d$Lab, levels = unique(d$Lab[order(d$R2)]))
  d$Sig <- local_sig(d$p_adjusted)
  d$Sig[d$Sig == "--"] <- ""
  g <- ggplot2::ggplot(d, ggplot2::aes(y = .data$Lab, x = .data$R2, fill = .data$Model)) +
    ggplot2::geom_col(width = 0.62, colour = "grey25", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Sig), hjust = -0.15,
                       size = pt2gg(MD_PT$ann * 0.85), colour = "grey25") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(x = "Effect size (PERMANOVA R2 / ANOSIM R / MRPP A)", y = NULL,
                  tag = "A", fill = NULL) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = MD_PT$tick * 0.85))
  g
}

#' 面板 B：组内离散度（betadisper 到组中位数的距离）箱线图，按因子分面。
#' 说明：箱体来自 PERMDISP 的那次 betadisper 拟合（各组 n >= min_group_n）
local_fig_disp <- function(disp_df, min_group_n) {
  if (is.null(disp_df) || !nrow(disp_df) ||
      !all(c("Group", "Factor", "Metric") %in% names(disp_df)))
    return(md_empty_plot(NULL, "PERMDISP not estimable"))
  d <- disp_df[!is.na(disp_df$median), , drop = FALSE]
  if (!nrow(d)) return(md_empty_plot(NULL, "PERMDISP not estimable"))
  # 箱线：用 mean +/- sd（betadisper 摘要只有汇总量，这里以均值为中心、sd 为须）
  d$Factor <- factor(as.character(d$Factor), levels = unique(as.character(d$Factor)))
  d$Group  <- factor(as.character(d$Group), levels = unique(as.character(d$Group)))
  d$lab <- paste0(d$Metric, " | ", d$Factor)
  d$lab <- factor(d$lab, levels = unique(d$lab))
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Group, y = .data$mean)) +
    ggplot2::geom_col(fill = "grey80", colour = "grey25", linewidth = 0.2, width = 0.62) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(.data$mean - .data$sd, 0),
                                        ymax = .data$mean + .data$sd),
                           width = 0.22, linewidth = MD_LW$errorbar, colour = "grey20") +
    ggplot2::geom_text(ggplot2::aes(label = .data$n), vjust = -0.35,
                       size = pt2gg(MD_PT$ann * 0.8), colour = "grey35") +
    ggplot2::facet_wrap(~ lab, scales = "free_x", nrow = 1) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "Distance to group median", tag = "B") +
    theme_md(CFG = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1,
                                                       size = MD_PT$tick * 0.9))
  ann <- do.call(rbind, lapply(split(d, d$lab), function(x)
    data.frame(lab = x$lab[1], txt = paste0("PERMDISP F = ", md_fmt_num(x$F[1], 2),
                                            ", P = ", md_fmt_p(x$P[1])),
               stringsAsFactors = FALSE)))
  ann$lab <- factor(ann$lab, levels = levels(d$lab))
  g <- g + ggplot2::geom_text(data = ann, ggplot2::aes(x = -Inf, y = Inf, label = .data$txt),
                              inherit.aes = FALSE, hjust = -0.06, vjust = 1.6,
                              size = pt2gg(MD_PT$ann * 0.85))
  g
}
