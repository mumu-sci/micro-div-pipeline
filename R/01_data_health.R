# =============================================================================
# R/01_data_health.R —— 数据体检（**硬门**：不通过就停）
# -----------------------------------------------------------------------------
# 检查：五文件完整性 / 维度 / 列名 / 缺失 / sample_id 与树 tip 对齐 /
#       taxonomy 层级完整度 / 每因子水平与交叉格样本量 / 读长分布 /
#       响应变量诊断（Shapiro / Bartlett / Levene / 离群点）
#
# 输出：results/01_data_health.xlsx（多 sheet）、figures/01_data_health.pdf
# 缓存：cache/01_health.rds（含 issues，供下游判断）
# =============================================================================

md_data_health <- function(CFG) {
  module <- "01"
  if (!md_need(CFG, "io", module)) return(md_result(module, "SKIP", "缺少 openxlsx"))

  inp <- md_read_inputs(CFG)
  asv <- inp$asv; tax <- inp$tax; meta <- inp$meta; env <- inp$env; tree <- inp$tree
  facs <- CFG$design$factors
  facs <- facs[facs %in% names(meta)]
  smp  <- inp$sample_col
  issues <- list()
  add_issue <- function(level, msg) issues[[length(issues) + 1L]] <<-
    data.frame(level = level, issue = msg, stringsAsFactors = FALSE)

  # ---- 1. 总览 ----
  readsum <- colSums(asv)
  overview <- data.frame(
    item = c("ASV 数", "样本数", "读长总量", "每样本读长 中位数", "每样本读长 最小",
             "每样本读长 最大", "每样本读长 最大值/最小值",
             "单样本检出 ASV 中位数", "零 reads 样本数", "全零 ASV 数",
             "ASV 出现率 中位数", "仅出现 1 个样本的 ASV 数"),
    value = c(nrow(asv), ncol(asv), format(sum(asv), big.mark = ","),
              format(round(stats::median(readsum)), big.mark = ","),
              format(min(readsum), big.mark = ","),
              format(max(readsum), big.mark = ","),
              sprintf("%.1f", max(readsum) / max(min(readsum), 1)),
              round(stats::median(colSums(asv > 0))),
              sum(readsum == 0), sum(rowSums(asv) == 0),
              sprintf("%.1f%%", 100 * stats::median(rowSums(asv > 0) / ncol(asv))),
              sum(rowSums(asv > 0) == 1)),
    stringsAsFactors = FALSE)

  if (max(readsum) / max(min(readsum), 1) > 20)
    add_issue("WARN", sprintf("读长跨度达 %.1f 倍，抽平会损失大量信息；建议考虑 CSS/相对丰度",
                              max(readsum) / max(min(readsum), 1)))
  if (any(readsum == 0)) add_issue("BLOCK", "存在零 reads 样本，必须剔除")

  # ---- 2. 文件指纹 ----
  fp <- md_data_md5(CFG)
  files_df <- data.frame(
    key = names(fp),
    path = vapply(names(fp), function(k) md_data_path(CFG, k) %||% "—", character(1)),
    md5 = unlist(fp), stringsAsFactors = FALSE)

  # ---- 3. taxonomy 完整度 ----
  tax_df <- NULL
  if (!is.null(tax)) {
    ranks <- intersect(CFG$data$tax_ranks, names(tax))
    tax_df <- do.call(rbind, lapply(ranks, function(r) {
      x <- as.character(tax[[r]][match(rownames(asv), tax[[inp$id_col]])])
      n_un <- sum(is.na(x) | x == "" | grepl("^(Unclassified|unclassified|unidentified|NA)$", x))
      data.frame(rank = r, n_unique = length(unique(x[!is.na(x) & x != ""])),
                 n_unclassified = n_un, pct_unclassified = n_un / length(x),
                 stringsAsFactors = FALSE)
    }))
    if (nrow(tax_df) && max(tax_df$pct_unclassified) > 0.5)
      add_issue("WARN", sprintf("层级 %s 的未分类比例达 %.0f%%，组成图结论需谨慎",
                                tax_df$rank[which.max(tax_df$pct_unclassified)],
                                100 * max(tax_df$pct_unclassified)))
  } else {
    add_issue("WARN", "未提供 taxonomy 表 —— 组成图模块将跳过")
  }

  # ---- 4. 因子水平与交叉格 ----
  fac_df <- NULL; cross_df <- NULL
  for (f in facs) {
    lv <- CFG$design$ordered_levels[[f]]
    x <- as.character(meta[[f]])
    if (!is.null(lv)) {
      not_in <- setdiff(unique(x), lv)
      if (length(not_in)) add_issue("WARN", paste0(
        "因子 ", f, " 里有未在 ordered_levels 中声明的水平: ", paste(not_in, collapse = ", ")))
    }
    t <- table(x)
    fac_df <- rbind(fac_df, data.frame(
      factor = f, level = names(t), n = as.integer(t),
      pct = as.numeric(t) / sum(t), stringsAsFactors = FALSE))
    small <- names(t)[t < 2]
    if (length(small)) add_issue("BLOCK", paste0(
      "因子 ", f, " 的水平样本量 < 2: ", paste(small, collapse = ", "),
      " —— 该水平的组间比较不可做"))
    else if (any(t < 4)) add_issue("WARN", paste0(
      "因子 ", f, " 有水平样本量 < 4: ", paste(names(t)[t < 4], collapse = ", "),
      " —— 组内统计量不稳定"))
  }
  if (length(facs) >= 2) {
    tb <- table(meta[[facs[1]]], meta[[facs[2]]])
    cross_df <- as.data.frame(tb, stringsAsFactors = FALSE)
    names(cross_df) <- c(facs[1], facs[2], "n")
    missing_cells <- cross_df[cross_df$n == 0, , drop = FALSE]
    if (nrow(missing_cells)) add_issue("WARN", sprintf(
      "%s x %s 有 %d 个空格（缺组合），交互作用分析的可估性受影响",
      facs[1], facs[2], nrow(missing_cells)))
    thin <- cross_df[cross_df$n > 0 & cross_df$n < 3, , drop = FALSE]
    if (nrow(thin)) add_issue("WARN", sprintf(
      "%s x %s 有 %d 个格样本量 < 3，两两比较将跳过这些格",
      facs[1], facs[2], nrow(thin)))
  }

  # ---- 5. 树对齐 ----
  tree_df <- data.frame(
    item = c("是否提供树", "tip 数", "是否 rooted", "被 ASV 表覆盖的 tip 比例",
             "被树覆盖的 ASV 比例"),
    value = c(if (is.null(tree)) "否" else "是",
              if (is.null(tree)) "—" else as.character(length(tree$tip.label)),
              if (is.null(tree)) "—" else as.character(ape::is.rooted(tree)),
              if (is.null(tree)) "—" else sprintf("%.1f%%", 100 * mean(rownames(asv) %in% tree$tip.label)),
              if (is.null(tree)) "—" else sprintf("%.1f%%", 100 * mean(tree$tip.label %in% rownames(asv)))),
    stringsAsFactors = FALSE)
  if (!is.null(tree)) {
    cov <- mean(rownames(asv) %in% tree$tip.label)
    if (cov < 0.5) add_issue("BLOCK", sprintf(
      "树只覆盖 %.1f%% 的 ASV —— 建议 phylo$enabled=FALSE，系统发育模块全部降级",
      100 * cov))
    else if (cov < 0.95) add_issue("WARN", sprintf(
      "树覆盖 %.1f%% 的 ASV（%.0f 个 ASV 无对应 tip）—— 系统发育指标用子集口径，须在表注披露",
      100 * cov, sum(!rownames(asv) %in% tree$tip.label)))
  } else if (isTRUE(CFG$phylo$enabled)) {
    add_issue("WARN", "phylo$enabled=TRUE 但没有树文件 —— 系统发育模块将降级")
  }

  # ---- 6. 样本读长 ----
  reads_df <- data.frame(sample_id = names(readsum), reads = as.integer(readsum),
                         n_asv = as.integer(colSums(asv > 0)),
                         stringsAsFactors = FALSE)
  reads_df <- merge(reads_df, meta, by = smp, all.x = TRUE, sort = FALSE)
  reads_df <- reads_df[match(names(readsum), reads_df[[smp]]), , drop = FALSE]
  q1 <- stats::quantile(readsum, .25); q3 <- stats::quantile(readsum, .75)
  iqr <- q3 - q1
  reads_df$outlier <- readsum < (q1 - 3 * iqr) | readsum > (q3 + 3 * iqr)
  if (any(reads_df$outlier)) add_issue("WARN", sprintf(
    "%d 个样本读长极端（3xIQR 外）：%s", sum(reads_df$outlier),
    paste(reads_df$sample_id[reads_df$outlier], collapse = ", ")))

  # ---- 7. 响应变量诊断 ----
  resp_df <- NULL; resp <- NULL
  if (isTRUE(CFG$response$enabled) && !is.null(env)) {
    rv <- CFG$response$var
    if (!rv %in% names(env)) {
      add_issue("BLOCK", paste0("response$var = '", rv, "' 不在环境/响应变量表里"))
    } else {
      y <- as.numeric(env[[rv]])
      d <- data.frame(y = y)
      for (f in facs) d[[f]] <- meta[[f]][match(env[[smp]], meta[[smp]])]
      sh <- tryCatch(stats::shapiro.test(d$y[!is.na(d$y)]), error = function(e) NULL)
      bar <- tryCatch(stats::bartlett.test(y ~ d[[facs[1]]], data = d), error = function(e) NULL)
      lev <- tryCatch(car::leveneTest(y ~ factor(d[[facs[1]]]), data = d), error = function(e) NULL)
      z <- scale(d$y)[, 1]
      resp_df <- data.frame(
        item = c("n", "mean", "sd", "median", "min", "max", "skewness",
                 "零值数", "Shapiro W", "Shapiro P", "Bartlett P (factor1)",
                 "Levene P (factor1)", "3xIQR 外样本数"),
        value = c(sum(!is.na(d$y)), sprintf("%.3f", mean(d$y, na.rm = TRUE)),
                  sprintf("%.3f", stats::sd(d$y, na.rm = TRUE)),
                  sprintf("%.3f", stats::median(d$y, na.rm = TRUE)),
                  sprintf("%.3f", min(d$y, na.rm = TRUE)),
                  sprintf("%.3f", max(d$y, na.rm = TRUE)),
                  sprintf("%.2f", mean(z^3, na.rm = TRUE)),
                  sum(d$y == 0, na.rm = TRUE),
                  if (is.null(sh)) "—" else sprintf("%.3f", sh$statistic),
                  if (is.null(sh)) "—" else md_fmt_p(sh$p.value),
                  if (is.null(bar)) "—" else md_fmt_p(bar$p.value),
                  if (is.null(lev)) "—" else md_fmt_p(lev[["Pr(>F)"]][1]),
                  sum(abs(z) > 3, na.rm = TRUE)),
        stringsAsFactors = FALSE)
      if (!is.null(sh) && sh$p.value < 0.05)
        add_issue("INFO", sprintf("响应变量 %s 非正态（Shapiro P = %s）—— 组间比较建议用非参数法，或先做变换",
                                  rv, md_fmt_p(sh$p.value)))
      resp <- d
    }
  }

  # ---- 8. 环境变量缺失 ----
  envmiss <- NULL
  if (!is.null(env)) {
    nm <- vapply(env, function(x) mean(is.na(x)), numeric(1))
    nm <- nm[setdiff(names(nm), smp)]
    envmiss <- data.frame(variable = names(nm), missing_rate = as.numeric(nm),
                          stringsAsFactors = FALSE)
    bad <- envmiss$variable[envmiss$missing_rate > 0.2]
    if (length(bad)) add_issue("WARN", paste0(
      "环境变量缺失率 > 20%: ", paste(bad, collapse = ", ")))
  }

  # ---- 9. 硬门 ----
  issues_df <- if (length(issues)) do.call(rbind, issues) else
    data.frame(level = "OK", issue = "未发现阻断性问题", stringsAsFactors = FALSE)
  n_block <- sum(issues_df$level == "BLOCK")
  if (n_block > 0) {
    md_log(CFG, "ERROR", "数据体检发现 ", n_block, " 个 BLOCK 级问题：")
    print(issues_df[issues_df$level == "BLOCK", , drop = FALSE], row.names = FALSE)
    if (isTRUE(CFG$opts$strict)) {
      stop("数据体检未通过（", n_block, " 个 BLOCK 级问题）。",
           "修正数据，或把 CFG$opts$strict 设为 FALSE 以只警告继续。", call. = FALSE)
    }
    md_warn(CFG, "体检有 BLOCK 级问题，但 strict=FALSE，继续执行")
  } else {
    md_log(CFG, "INFO", "数据体检通过（", sum(issues_df$level == "WARN"), " 个 WARN，无 BLOCK）")
  }

  # ---- 10. 出图 ----
  figs <- character(0)
  plots <- list()
  if (length(facs) >= 1) {
    plots[["A"]] <- md_box_group(reads_df, "reads", facs[1], ylab = "Reads per sample",
                                 xlab = NULL, x_angle = 0, show_n = TRUE)
  }
  st <- md_asv_stats(asv)
  plots[["B"]] <- ggplot2::ggplot(st, ggplot2::aes(x = mean_rel + 1e-6, y = n_present)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.5, colour = "grey30") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Mean relative abundance (log10)", y = "Occurrence (n samples)") +
    theme_md(CFG = NULL)
  if (!is.null(envmiss) && nrow(envmiss) > 0) {
    plots[["C"]] <- ggplot2::ggplot(envmiss, ggplot2::aes(x = stats::reorder(variable, missing_rate),
                                                          y = missing_rate)) +
      ggplot2::geom_col(fill = "#8FAE5B", width = 0.65) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(x = NULL, y = "Missing rate") + theme_md(CFG = NULL)
  }
  p <- md_compose(CFG, plots, design = if (length(plots) == 3) "AB\nCC" else NULL,
                  ncol = min(2, length(plots)), legend = NULL)
  md_save_figure(CFG, p, "01_data_health", width = CFG$figure$width_double, height = 110)
  figs <- "figures/01_data_health.pdf"

  # ---- 11. 写表 ----
  tp <- md_write_xlsx(CFG, "01_data_health", list(
    Overview = overview, Files = files_df, Issues = issues_df,
    Taxonomy = tax_df, Factors = fac_df, Factor_crosstab = cross_df,
    Tree_check = tree_df, Sample_reads = reads_df,
    Response_diagnostics = resp_df, Env_missing = envmiss),
    meta_extra = list(block = n_block, warn = sum(issues_df$level == "WARN")))

  saveRDS(list(overview = overview, issues = issues_df, n_block = n_block,
               facs = facs, samples = inp$samples, asvs = inp$asvs,
               fingerprint = fp, input_report = inp$report),
          md_path(CFG, "cache", "01_health.rds"))

  md_result(module, if (n_block > 0) "PARTIAL" else "OK",
            reason = if (n_block > 0) paste0(n_block, " 个 BLOCK 级问题（strict=FALSE 继续）") else NULL,
            tables = basename(tp), figures = figs)
}
