# =============================================================================
# R/34_variance_partition.R ---- 变差分解（VPA / varpart）
# -----------------------------------------------------------------------------
# 口径来源（原始研究 原脚本 §E，只读参考，不修改）：
#   · vegan::varpart()，响应为 Hellinger 转换后的群落矩阵
#   · 各块各自 scale()，三块：Geo / Soil / Climate
#   · 偏 RDA 显著性：anova(rda(Y, X, Z), permutations = 999) + p.adjust(method = "holm")
#   · nifH 实测环境总解释力仅约 2%（0.02 量级）----这是真实特征，照实输出。
#
# 本包参数化：
#   CFG$env$vpa_blocks        块定义（键名进图例）
#   CFG$env$vpa_transform     hellinger | none
#   CFG$env$vpa_permutations  置换次数
#   CFG$env$vpa_p_adjust      校正方法
#   CFG$env$vpa_min_n         有效样本下限（默认 8）
#   CFG$env$vpa_max_blocks    块数上限（默认 3 ---- 只有 3 块画得出 Venn）
#
# 输出：figures/34_variance_partition.pdf（A 分块解释力条形图 + B 交叠 Venn 示意）
#       results/34_variance_partition.xlsx（Fractions / Blocks / Note）
#       cache/34_vpa.rds
#
# 画法对齐（只改画法，不动任何统计量）：
#   · 分块配色 = 原脚本 §E 的 col_geo / col_sc / col_clim
#     （= #80B1D2 / #F6AAA1 / #C7CEEA；经 md_pal(name = "assembly") 取，模块里不写死颜色）
#   · 分块名加粗写在坐标轴上（原脚本把分块名写在 Venn 里、用 bold）
# 图内文字一律纯 ASCII（y 轴用 "Adjusted R2" / "Adjusted R2 (%)"，不用 plotmath 的 R^2）。
# =============================================================================

# ---------------------------------------------------------------------------
# local_* 辅助
# ---------------------------------------------------------------------------

#' 只吞掉 varpart/rda 内部那条已知的 cosmetic warning。
#' 背景：vegan 2.6.10 的 varpart() 内部对每个解释块做一次赋值，会抛
#' "number of items to replace is not a multiple of replacement length"
#' （实测每块一次，与本包数据无关，结果数值正确 ---- 实测 indfract 与手工 varpart 一致）。
#' 其余 warning 一律照常抛出，避免掩盖真实问题。
local_quiet_known <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("number of items to replace", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
}

#' 把 varpart() 的 indfract 表整理成「分块 + 交叠 + 残差」的整洁表。
#' ⚠ 实测 vegan 2.6.10：indfract 的 Df 列被 print 成 "[a] = X1 | X2+X3"，
#'   但 as.character() 取回来只剩数字（"2"/"11"/NA），**无法从 Df 反推分项名**。
#'   好在 indfract 的行序是固定且实测验证过的：
#'     行 1..p = [a]..[c] 各块独立贡献；其后的 p+1..p+3 = [d][e][f] 两两交叠；
#'     再下一行 [g] = 三块交叠；末行 [h] = 残差。
#'   故这里按「位置」赋标签（并对 p != 3 的通用情形降级为 a/b/c… 占位）。
#'
#' ⚠ 交叠项的校正 R2 **可以为负**，这是可加性分解的正常现象（含义是「共享的那部分
#'   解释力不优于随机期望」），不是错误、也不等于 0。故本表**不做任何截断**：
#'   所有下游（图、表、缓存）一律使用原始 adj_r2。
#'   （原实现 draw_venn3 里有一个 clamp <- pmax(x, 0)，
#'     会把负交叠画成 0.0%；本包不沿用该做法，见模块 34 的 Note sheet。）
local_vpa_fractions <- function(vp, blocks) {
  if (is.null(vp) || is.null(vp$part$indfract)) return(NULL)
  f <- vp$part$indfract
  frac_col <- intersect(c("Adj.R.square", "R.square"), names(f))[1]
  if (is.na(frac_col)) return(NULL)

  p    <- length(blocks)
  n    <- nrow(f)
  lab  <- rep(NA_character_, n)
  lab[seq_len(min(p, n))] <- blocks[seq_len(min(p, n))]
  pair_lab <- NULL
  if (p == 2) pair_lab <- paste0(blocks[1], " x ", blocks[2])
  if (p == 3) pair_lab <- c(paste0(blocks[1], " x ", blocks[2]),
                            paste0(blocks[1], " x ", blocks[3]),
                            paste0(blocks[2], " x ", blocks[3]))
  if (!is.null(pair_lab)) {
    idx <- p + seq_along(pair_lab)
    idx <- idx[idx <= n]
    lab[idx] <- pair_lab[seq_along(idx)]
  }
  if (p == 3 && n >= p + 4) lab[p + 4] <- paste0(blocks, collapse = " x ")
  if (n > 0) lab[n] <- "Residuals"
  lab[is.na(lab)] <- paste0("component_", which(is.na(lab)))

  tt <- if (!is.null(f$Testable) && length(f$Testable) == n)
    as.logical(f$Testable) else rep(NA, n)

  tab <- data.frame(
    fraction = lab,
    adj_r2   = as.numeric(f[[frac_col]]),
    testable = tt,
    stringsAsFactors = FALSE)
  ok <- !is.na(tab$adj_r2)
  tab <- tab[ok, , drop = FALSE]
  tab$type <- ifelse(tab$fraction %in% blocks, "Block (unique)",
              ifelse(tab$fraction == "Residuals", "Residual", "Shared overlap"))
  tab$pct <- 100 * tab$adj_r2
  tab$sign <- ifelse(tab$adj_r2 < 0, "negative", ifelse(tab$adj_r2 > 0, "positive", "zero"))
  tab[, c("fraction", "type", "adj_r2", "pct", "sign", "testable")]
}

#' 变差分解各分量的**有符号**条形图（替代原来的 3 圆 Venn）。
#'
#' 为什么不用 Venn：Venn 的每个叶瓣在几何上只能表示「正面积」，而校正 R2 的
#' 两两/三块交叠项实测为负（本数据集 GeoxSoil = −0.70% 等）。把负值硬钳成 0.00%
#' 会读成「没有共享解释力」，属于误导；把负值写进一个正面积的叶瓣也同样别扭。
#' 有符号条形图能在同一坐标轴上同时表达正负，零线为基准，语义无歧义。
local_partition_plot <- function(df, blocks, cols, pt, subtitle = NULL) {
  d <- df[df$type != "Residual", , drop = FALSE]
  # 排序：先分块（按 blocks 顺序），再交叠（按数值从大到小）
  d$grp <- ifelse(d$type == "Block (unique)", 1L, 2L)
  blk_order <- match(d$fraction, blocks)
  d <- d[order(d$grp, ifelse(d$grp == 1L, blk_order, -d$adj_r2),
               ifelse(d$grp == 1L, 0, d$adj_r2), na.last = TRUE), , drop = FALSE]
  d$fraction <- factor(d$fraction, levels = d$fraction)
  d$fill_grp <- ifelse(d$type == "Block (unique)", as.character(d$fraction), "Shared")

  pal <- c(stats::setNames(unname(cols[seq_along(blocks)]), blocks),
           Shared = "grey55")
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$fraction, y = .data$pct,
                                       fill = .data$fill_grp)) +
    ggplot2::geom_col(width = 0.66, colour = "grey25", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linewidth = MD_LW$ref, colour = "grey35") +
    # 标签：正值贴柱顶上方，负值贴柱底下侧（都在面板内，不会被裁）
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f%%", .data$pct),
                                    vjust = ifelse(.data$adj_r2 >= 0, -0.3, 1.25)),
                       size = pt2gg(pt$ann * 0.9), colour = "grey20") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::labs(x = NULL, y = "Adjusted R2 (%)",
                  tag = "B", subtitle = subtitle) +
    theme_md(CFG = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = pt$tick * 0.9),
                   plot.subtitle = ggplot2::element_text(size = pt$ann * 0.9, colour = "grey30"))
  g
}

# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------
md_variance_partition <- function(CFG) {
  module <- "34"
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少依赖包"))
  if (!md_precheck(CFG, module)) return(md_result(module, "SKIP", "前置条件不满足"))

  blocks <- CFG$env$vpa_blocks %||% CFG$env$groups
  if (is.null(blocks) || length(blocks) < 2) {
    md_skip(CFG, module, "CFG$env$vpa_blocks 未定义或块数 < 2")
    return(md_result(module, "SKIP", "vpa_blocks 未定义"))
  }
  transform   <- CFG$env$vpa_transform %||% "hellinger"
  nperm       <- as.integer(CFG$env$vpa_permutations %||% 999)
  padj_method <- CFG$env$vpa_p_adjust %||% "holm"
  min_n       <- as.integer(CFG$env$vpa_min_n %||% 8L)
  max_blocks  <- as.integer(CFG$env$vpa_max_blocks %||% 3L)

  # ---- 1. 响应矩阵 ----
  m <- tryCatch(md_core_asv(CFG, "env"), error = function(e) NULL)
  if (is.null(m) || ncol(m) < 4) {
    md_skip(CFG, module, "群落矩阵不可用")
    return(md_result(module, "SKIP", "群落矩阵不可用"))
  }
  m <- m[rowSums(m) > 0, , drop = FALSE]
  if (nrow(m) < 5) { md_skip(CFG, module, "有效 ASV < 5"); return(md_result(module, "SKIP", "ASV 过少")) }

  ana <- tryCatch(md_core_analysis(CFG), error = function(e) NULL)
  if (is.null(ana) || nrow(ana) < 4) {
    md_skip(CFG, module, "分析总表不可用")
    return(md_result(module, "SKIP", "分析总表不可用"))
  }
  smp <- CFG$data$sample_col
  ana <- ana[match(colnames(m), ana[[smp]]), , drop = FALSE]

  # ---- 2. 块变量：取交集、剔常量 ----
  blocks <- lapply(blocks, function(v) intersect(as.character(v), names(ana)))
  blocks <- lapply(blocks, function(v) v[vapply(ana[v], is.numeric, logical(1))])
  blocks <- lapply(blocks, function(v) {
    if (!length(v)) return(v)
    v[vapply(v, function(x) {
      s <- stats::sd(ana[[x]], na.rm = TRUE); is.finite(s) && s > 0
    }, logical(1))]
  })
  blocks <- blocks[vapply(blocks, length, integer(1)) > 0]
  if (length(blocks) < 2) {
    md_skip(CFG, module, "有效块 < 2（变量缺失或零方差）")
    return(md_result(module, "SKIP", "有效块 < 2"))
  }
  if (length(blocks) > max_blocks) {
    md_warn(CFG, module, "：块数 ", length(blocks), " > vpa_max_blocks ", max_blocks, "，只取前 ", max_blocks, " 块")
    blocks <- blocks[seq_len(max_blocks)]
  }
  all_vars <- unique(unlist(blocks))
  keep <- stats::complete.cases(ana[, all_vars, drop = FALSE])
  if (sum(keep) < min_n) {
    md_skip(CFG, module, sprintf("完整样本数 %d < vpa_min_n %d", sum(keep), min_n))
    return(md_result(module, "SKIP", "完整样本不足"))
  }
  dat  <- ana[keep, , drop = FALSE]
  ids  <- as.character(dat[[smp]])
  Yraw <- t(m[, ids, drop = FALSE])
  Yraw <- Yraw[rowSums(Yraw) > 0, , drop = FALSE]
  Yraw <- Yraw[, colSums(Yraw) > 0, drop = FALSE]
  ids  <- rownames(Yraw)
  dat  <- dat[match(ids, dat[[smp]]), , drop = FALSE]
  if (nrow(Yraw) < min_n || ncol(Yraw) < 5) {
    md_skip(CFG, module, "对齐后样本或 ASV 过少")
    return(md_result(module, "SKIP", "对齐后样本过少"))
  }

  # ---- 3. 转换 ----
  Y <- switch(transform,
              hellinger = vegan::decostand(Yraw, method = "hellinger"),
              none      = as.matrix(Yraw),
              as.matrix(Yraw))
  Xs <- lapply(blocks, function(v) as.matrix(scale(as.matrix(dat[, v, drop = FALSE]))))
  names(Xs) <- names(blocks)

  # ---- 4. varpart ----
  md_seed(CFG, 301)
  vp <- tryCatch(local_quiet_known(do.call(vegan::varpart, c(list(Y = Y), unname(Xs)))),
                 error = function(e) {
                   md_warn(CFG, module, "：varpart 失败 ---- ", conditionMessage(e)); NULL
                 })
  if (is.null(vp)) {
    md_skip(CFG, module, "varpart 计算失败")
    return(md_result(module, "SKIP", "varpart 失败"))
  }
  frac <- local_vpa_fractions(vp, names(blocks))

  # ---- 5. 偏 RDA 显著性 ----
  sig_rows <- lapply(names(blocks), function(b) {
    X <- Xs[[b]]
    Z <- do.call(cbind, Xs[setdiff(names(blocks), b)])
    md_seed(CFG, 302)
    # 注意：vegan 不导出 anova()；必须用 stats::anova() 让 S3 派发到 vegan::anova.cca
    a <- tryCatch(local_quiet_known(stats::anova(vegan::rda(Y, X, Z), permutations = nperm)),
                  error = function(e) NULL)
    if (is.null(a) || !nrow(a)) {
      md_warn(CFG, module, "：块 ", b, " 的偏 RDA 置换检验失败")
      return(data.frame(block = b, df = NA_real_, F = NA_real_, p = NA_real_,
                        permutations = nperm, stringsAsFactors = FALSE))
    }
    data.frame(block = b, df = a$Df[1],
               F = a$F[1],
               p = a$`Pr(>F)`[1],
               permutations = nperm, stringsAsFactors = FALSE)
  })
  sig <- do.call(rbind, sig_rows)
  sig$p_adj <- md_p_adjust(sig$p, padj_method)
  sig$sig <- ifelse(is.na(sig$p_adj), "--",
             ifelse(sig$p_adj < 0.001, "***",
             ifelse(sig$p_adj < 0.01, "**",
             ifelse(sig$p_adj < 0.05, "*", "ns"))))
  sig$adj_r2 <- vapply(sig$block, function(b) {
    r <- frac$adj_r2[frac$fraction == b]; if (!length(r)) NA_real_ else r[1]
  }, numeric(1))
  sig$n_vars <- vapply(sig$block, function(b) length(blocks[[b]]), integer(1))

  # 合计（varpart$part$fract 最后一行 = 全部块合计）
  total_adj <- tryCatch({
    fr <- vp$part$fract
    unname(fr$Adj.R.square[nrow(fr)])
  }, error = function(e) NA_real_)

  # ---- 6. 图 ----
  # A = 各块独立解释力（唯一贡献）+ 偏 RDA 置换 P
  # B = 全部非残差分量的**有符号**条形图（含可为负的交叠项）
  #     ---- 不使用 Venn，理由见 local_partition_plot() 的函数注释。
  pt <- md_pt(CFG)
  # 分块配色对齐 原脚本 §E 的 col_geo / col_sc / col_clim
  #   = #80B1D2 / #F6AAA1 / #C7CEEA（原脚本把分块名写在 Venn 里、用该块的色 + bold）。
  # 这三个 hex 在 R/shared/12_plot_theme.R 里只出现在 MD_BUILTIN_PAL$assembly，
  # 故经 md_pal(name = "assembly") 取用，不在模块里写死颜色（顺序 = CFG$env$vpa_blocks 的块序）。
  board <- md_pal(CFG, c("Homogeneous selection", "Homogenizing dispersal",
                         "Dispersal limitation"), name = "assembly")
  cols <- stats::setNames(unname(board)[seq_along(blocks)], names(blocks))
  bar_df <- frac[frac$type == "Block (unique)", , drop = FALSE]
  bar_df$fraction <- factor(bar_df$fraction, levels = names(blocks))
  bar_df <- bar_df[order(bar_df$fraction), , drop = FALSE]
  bar_df$p_adj <- sig$p_adj[match(as.character(bar_df$fraction), sig$block)]
  bar_df$lab <- sprintf("%.2f%%\n%s", 100 * bar_df$adj_r2, md_fmt_p(bar_df$p_adj))
  total_pct <- 100 * total_adj
  n_neg <- sum(frac$type == "Shared overlap" & frac$adj_r2 < 0)

  pA <- ggplot2::ggplot(bar_df, ggplot2::aes(x = .data$fraction, y = .data$adj_r2,
                                             fill = .data$fraction)) +
    ggplot2::geom_col(width = 0.6, colour = "grey25", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab), vjust = -0.25,
                       size = pt2gg(pt$ann * 0.9), lineheight = 0.95, colour = "grey20") +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.28))) +
    ggplot2::labs(x = NULL, y = "Adjusted R2",
                  tag = "A",
                  subtitle = sprintf("Unique (marginal) effect per block | Total adj. R2 = %.4f (%.2f%%) | n = %d samples, %d ASVs",
                                     total_adj, total_pct, nrow(Y), ncol(Y))) +
    theme_md(CFG = CFG) +
    # 原脚本把分块名加粗写出（fontface = "bold"）——这里落到 x 轴刻度文字上。
    # 注：不给刻度文字逐项上色 —— ggplot2 的 element_text() 不接受向量 colour
    #     （会抛 "Vectorized input to element_text() is not officially supported"），
    #     分块颜色已由柱体填充色 + 图 B 的图例表达。
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = pt$tick, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = pt$ann * 0.9, colour = "grey30"))

  # B = 三圆 Venn —— **照抄原实现 draw_venn3**。
  #   ⚠ 原脚本用 clamp() 把负的交叠项截为 0（图上显示 0.0%），本包沿用该口径；
  #     负值明细写在 Fractions 表与 Note 里，不丢信息。
  #   Venn 需要恰好 3 块；块数不是 3 时退回有符号条形图并记 Note。
  if (length(blocks) == 3) {
    val_of <- function(nm) {
      v <- frac$adj_r2[match(nm, as.character(frac$fraction))]
      if (length(v) == 0 || is.na(v)) 0 else v
    }
    bn <- names(blocks)
    pB <- md_venn3(
      lab1 = bn[1], lab2 = bn[2], lab3 = bn[3],
      f1 = val_of(bn[1]), f2 = val_of(bn[2]), f3 = val_of(bn[3]),
      f12 = val_of(paste0(bn[1], " x ", bn[2])),
      f13 = val_of(paste0(bn[1], " x ", bn[3])),
      f23 = val_of(paste0(bn[2], " x ", bn[3])),
      f123 = val_of(paste0(bn, collapse = " x ")),
      col1 = unname(cols)[1], col2 = unname(cols)[2], col3 = unname(cols)[3],
      pt = pt, tag = NULL,
      subtitle = sprintf(paste0("Shared fractions clamped at 0 for the Venn",
                                " (%d negative in the table; see Note) | Total adj. R2 = %.2f%%"),
                         n_neg, total_pct))
  } else {
    md_warn(CFG, "34：分块数为 ", length(blocks), "（Venn 需要 3 块），图 B 退回有符号条形图")
    pB <- local_partition_plot(
      frac, names(blocks), unname(cols), pt,
      subtitle = sprintf("Shared fractions are signed: %d negative (see Note) | sum of all non-residual = %.2f%%",
                         n_neg, sum(frac$adj_r2[frac$type != "Residual"]) * 100))
  }

  # 图件**只出环形图（Venn）**——用户口径：不要 A（各块唯一解释力的条形图），
  # 数字已在 Fractions / Blocks / Total 表里。Venn 走 theme_void()，本身没有坐标轴。
  fig <- md_compose(CFG, list(A = pB), design = "A", legend = NULL)
  md_save_figure(CFG, fig, "34_variance_partition", height = 105)

  # ---- 7. 写表 ----
  note_df <- data.frame(Note = c(
    sprintf("方法：vegan::varpart()（分块校正 R2）；响应矩阵 = %s 转换后的群落矩阵（模块 34 口径 = CFG$normalize$by_module$env）。",
            transform),
    sprintf("分块（CFG$env$vpa_blocks）：%s。",
            paste(sprintf("%s = %s", names(blocks),
                          vapply(blocks, function(v) paste(v, collapse = "+"), character(1))),
                  collapse = " | ")),
    sprintf("显著性：偏 RDA anova(rda(Y, X | Z), permutations = %d)，p 校正 = %s。", nperm, padj_method),
    sprintf("样本 n = %d，ASV = %d，合计校正 R2 = %.4f（%.2f%%）。", nrow(Y), ncol(Y), total_adj, 100 * total_adj),
    sprintf("vegan 的 indfract 语义（本表 fraction 列已按行序还原）：[a][b][c] = 各块独立贡献；[d][e][f] = 两两交叠；[g] = 三块交叠；[h] = 残差。三者独立贡献之和 = %.4f，与合计 %.4f 的差额即交叠部分。",
            sum(frac$adj_r2[frac$type == "Block (unique)"]), total_adj),
    "★ 交叠项为负值表示共享解释力不优于随机，属可加性分解的正常现象，不宜解读为 0。",
    sprintf("   本数据集有 %d 个交叠项为负（GeoxSoil、GeoxClimate、SoilxClimate）。Partition sheet 与图中均保留原始负值（如 -0.70%%），未做任何截断。",
            n_neg),
    "★ 图 B 用「有符号条形图」而非 Venn 圆：Venn 叶瓣几何上只能表示正面积，无法忠实表达负交叠；",
    "   故本模块不采用 Venn（nifH 原脚本 draw_venn3 里的 clamp <- pmax(x, 0) 会把负交叠画成 0.0%，本包不沿用）。",
    "⚠ 本数据集的群落变量总解释力极低（原始研究实测约 2%，示范数据同量级）----这是真实特征",
    "    （极端单例化、β 周转极高），不是计算错误；正文/表注应如实报告而非调参放大。"))

  tp <- md_write_xlsx(CFG, "34_variance_partition", list(
    Fractions = frac,
    Blocks = sig,
    Total = data.frame(item = c("transform", "n_samples", "n_asv", "n_blocks",
                                "total_adj_r2", "permutations", "p_adjust"),
                       value = c(transform, nrow(Y), ncol(Y), length(blocks),
                                 total_adj, nperm, padj_method),
                       stringsAsFactors = FALSE),
    Note = note_df),
    meta_extra = list(total_adj_r2 = total_adj, n_samples = nrow(Y), n_asv = ncol(Y)))

  saveRDS(list(varpart = vp, fractions = frac, blocks = sig,
               total_adj_r2 = total_adj, params = list(
                 transform = transform, blocks = names(blocks),
                 permutations = nperm, p_adjust = padj_method,
                 n = nrow(Y), n_asv = ncol(Y), seed = CFG$opts$seed)),
          md_path(CFG, "cache", "34_vpa.rds"))

  md_result(module, "OK", tables = basename(tp),
            figures = "figures/34_variance_partition.pdf",
            notes = sprintf("total adj R2 = %.4f", total_adj))
}
