# =============================================================================
# R/20_phylo_prep.R —— 系统发育树准备与可视化
# -----------------------------------------------------------------------------
# 做：tip 与 ASV 对齐检查 -> 按 CFG$phylo$root_method 生根 -> 二叉化 ->
#     剪到当前 ASV 集 -> 出 fan 树图（ggtreeExtra::geom_fruit 多环：
#     分类色环 -> 各设计因子色环 -> 最外丰度条形环），画法对照 01_Fig01.R §A
#
# 无树 / phylo$enabled = FALSE 时优雅 SKIP（ITS 用户的典型场景）。
#
# 产物：cache/20_phylo_prep.rds（供 21 / 22 与需要树的模块取用）
#       results/20_phylo_prep.xlsx、figures/20_phylo_prep.pdf
#
# 配置键：
#   phylo$enabled / root_method(midpoint|outgroup|none) / outgroup_tip
#   phylo$ring_rank   树外环按哪个分类层级着色（默认 Phylum）
#   phylo$ring_top_n  环上保留多少个分类（其余归 Other，默认 10）
# =============================================================================

md_phylo_prep <- function(CFG) {
  module <- "20"
  if (!md_precheck(CFG, module, need_phylo = TRUE)) {
    return(md_result(module, "SKIP", "未提供系统发育树或 phylo$enabled = FALSE"))
  }
  if (!md_need(CFG, c("ecology", "ggtree", "ggtreeExtra", "ggnewscale", "unifrac"), module)) {
    return(md_result(module, "SKIP", "缺少 ape / ggtree / ggtreeExtra / ggnewscale / phangorn"))
  }
  md_load(c("ecology", "ggtree", "ggtreeExtra", "ggnewscale", "unifrac"))

  core <- md_core(CFG)
  tr   <- core$tree
  m    <- md_core_asv(CFG, "phylo")
  tax  <- core$tax
  id_col <- core$id_col %||% CFG$data$id_col

  if (is.null(tr)) return(md_result(module, "SKIP", "核心对象里没有树"))

  raw_tip <- length(tr$tip.label)
  asvs    <- rownames(m)

  # ---- 1. tip 对齐 ----
  overlap   <- intersect(tr$tip.label, asvs)
  only_tree <- setdiff(tr$tip.label, asvs)
  only_asv  <- setdiff(asvs, tr$tip.label)
  cov_asv   <- length(overlap) / length(asvs)

  align_df <- data.frame(
    item = c("树的 tip 数", "分析的 ASV 数", "交集", "只在树中", "只在 ASV 表中",
             "ASV 被树覆盖的比例", "树 tip 被 ASV 覆盖的比例"),
    value = c(raw_tip, length(asvs), length(overlap), length(only_tree),
              length(only_asv),
              md_fmt_pct(cov_asv), md_fmt_pct(length(overlap) / raw_tip)),
    stringsAsFactors = FALSE)

  if (length(overlap) < 4) {
    md_skip(CFG, module, paste0("树与 ASV 的交集只有 ", length(overlap), " 个，无法做系统发育分析"))
    return(md_result(module, "SKIP", "树与 ASV 交集过少"))
  }
  if (cov_asv < 0.5) {
    md_warn(CFG, "树只覆盖 ", md_fmt_pct(cov_asv),
            " 的 ASV —— 系统发育指标属于子集口径，须在表注披露")
  }

  # ---- 2. 修剪到交集 ----
  tr <- ape::keep.tip(tr, overlap)
  tr$tip.label <- trimnames(tr$tip.label)

  # ---- 3. 生根 ----
  root_method <- CFG$phylo$root_method %||% "midpoint"
  root_used <- root_method
  if (ape::is.rooted(tr)) {
    root_used <- "already rooted"
  } else if (identical(root_method, "midpoint")) {
    tr <- tryCatch(phangorn::midpoint(tr), error = function(e) {
      md_warn(CFG, "midpoint 生根失败（", conditionMessage(e),
              "），保持无根并在此后按无根处理")
      root_used <<- "none (midpoint failed)"
      tr
    })
  } else if (identical(root_method, "outgroup")) {
    og <- CFG$phylo$outgroup_tip
    if (is.null(og) || !og %in% tr$tip.label) {
      og <- tr$tip.label[which.max(ape::cophenetic.phylo(tr)[, 1])]
      md_log(CFG, "INFO", "未指定外群，自动取 ", og)
    }
    tr <- tryCatch(ape::root(tr, outgroup = og, resolve.root = TRUE),
                   error = function(e) {
                     md_warn(CFG, "外群生根失败：", conditionMessage(e)); tr
                   })
  } else {
    root_used <- "none (保持无根)"
  }

  # ---- 4. 二叉化 ----
  n_poly <- 0L
  if (!ape::is.binary(tr)) {
    n_poly <- sum(ape::tabulate(tr$edge[, 1]) > 2)
    tr <- ape::multi2di(tr)
    md_log(CFG, "INFO", "多分支节点 ", n_poly, " 个，已 multi2di 二叉化")
  }
  tr <- ape::collapse.singles(tr)
  tr$node.label <- NULL

  # ---- 5. 缓存 ----
  prep <- list(
    tree = tr, root_method = root_used, n_tip = length(tr$tip.label),
    n_tip_raw = raw_tip, coverage_asv = cov_asv,
    n_only_tree = length(only_tree), n_only_asv = length(only_asv),
    n_polytomy = n_poly, is_binary = ape::is.binary(tr), is_rooted = ape::is.rooted(tr),
    total_bl = sum(tr$edge.length %||% NA_real_)
  )
  saveRDS(prep, md_path(CFG, "cache", "20_phylo_prep.rds"))

  # ---- 6. 分类环数据（照抄 01_Fig01.R §A 1.2/1.3/1.7，L120-147 + L217-227）----
  ring_rank <- CFG$phylo$ring_rank %||% "Phylum"
  top_n     <- CFG$phylo$ring_top_n %||% 10L
  ring_df <- NULL; ring_note <- "未提供 taxonomy，未画分类环"
  if (!is.null(tax) && ring_rank %in% names(tax)) {
    g <- as.character(tax[[ring_rank]][match(tr$tip.label, tax[[id_col]])])
    g[is.na(g) | g == ""] <- "Unclassified"
    tb <- sort(table(g), decreasing = TRUE)
    keep_g <- names(tb)[seq_len(min(top_n, length(tb)))]
    g2 <- ifelse(g %in% keep_g, g, "Other")
    lv <- c(keep_g, if (any(g2 == "Other")) "Other")
    ring_df <- data.frame(label = tr$tip.label, taxon = factor(g2, levels = lv),
                          RingLabel = ring_rank, stringsAsFactors = FALSE)
    ring_note <- paste0("分类环层级 = ", ring_rank, "；共 ", length(tb),
                        " 个分类，保留前 ", length(keep_g), " 个，其余归 Other")
  }

  # ---- 7. 分组环 / 丰度环的数据（照抄 01_Fig01.R §A 1.4-1.6，L149-211）----
  smp20  <- CFG$data$sample_col
  meta20 <- md_core_meta(CFG)
  meta20 <- meta20[match(colnames(m), meta20[[smp20]]), , drop = FALSE]
  m_ids  <- colnames(m)
  rel20  <- md_rel(m) * 100
  rel20[!is.finite(rel20)] <- 0

  # 外环：每 ASV 平均相对丰度 -> log1p(99% 分位缩尾)
  asv_mean <- rowMeans(rel20, na.rm = TRUE)
  cap      <- stats::quantile(asv_mean, 0.99, na.rm = TRUE)
  abund_df <- data.frame(label = rownames(rel20),
                         Abundance = log1p(pmin(asv_mean, cap)),
                         RingLabel = "Abundance", stringsAsFactors = FALSE)
  abund_df <- abund_df[match(tr$tip.label, abund_df$label), , drop = FALSE]

  facs20 <- intersect(CFG$design$factors %||% character(0), names(meta20))

  # 树环上的文字/轴线尺寸：照抄 01_Fig01.R L80 的 RING_TXT_MM <- pt2gg(S_LEG_COMP) * 0.5
  RING_TXT_MM  <- pt2gg(MD_PT_COMP$legend) * 0.5
  # 外环单色：照抄 01_Fig01.R L332 的 fill = colour = "#80B1D2"
  RING_BAR_COL <- "#80B1D2"
  # ⚠ axis.params 必须**就地写成字面 list**（照抄原脚本的写法）。实测 ggtreeExtra 1.20.1
  #   把 axis.params 存进变量再传进去时会静默丢掉整条环轴（图中的环名与刻度全部不出现）。

  # ---- 8. 出图：fan 树 + 多环（逐行照抄 01_Fig01.R §A L276-344）----
  #   基础树 ggtree(layout="fan", open.angle=20, size=0.04, colour="black")
  #   -> new_scale_fill() + geom_fruit(geom_tile) 逐环叠加（分类环 width 0.10 / 分组环 0.15）
  #   -> new_scale_fill() + geom_fruit(geom_bar) 丰度条形环（pwidth 0.22, offset 0.05）
  #   -> theme(legend.position="none", plot.margin=0)
  p <- ggtree::ggtree(tr, layout = "fan", open.angle = 20,
                      linewidth = 0.04, colour = "black")

  if (!is.null(ring_df)) {
    pal <- md_pal(CFG, levels(ring_df$taxon), name = ring_rank)
    if ("Other" %in% names(pal)) pal[["Other"]] <- "grey75"
    p <- p + ggtreeExtra::geom_fruit(
      data = ring_df, geom = geom_tile,
      mapping = ggplot2::aes(y = label, x = RingLabel, fill = taxon),
      width = 0.10, offset = 0.03,
      axis.params = list(axis = "x", text.size = RING_TXT_MM, text.angle = -90,
                         hjust = 0.5, vjust = 0.5)) +
      ggplot2::scale_fill_manual(values = pal, name = ring_rank,
                                 guide = md_guide_flex(length(pal), "bottom"))
  }

  # 中间各环：每个设计因子一环（每个 ASV 取该因子内平均相对丰度最高的水平）
  for (f in facs20) {
    lv <- CFG$design$ordered_levels[[f]] %||% unique(as.character(meta20[[f]]))
    lv <- intersect(as.character(lv), unique(as.character(meta20[[f]])))
    lv <- lv[!is.na(lv)]
    if (length(lv) < 2) next
    df_f <- data.frame(label = rownames(rel20),
                       g = factor(local20_dominant(rel20, meta20, m_ids, f, lv), levels = lv),
                       RingLabel = f, stringsAsFactors = FALSE)
    df_f <- df_f[match(tr$tip.label, df_f$label), , drop = FALSE]
    p <- p + ggnewscale::new_scale_fill() +
      ggtreeExtra::geom_fruit(
        data = df_f, geom = geom_tile,
        mapping = ggplot2::aes(y = label, x = RingLabel, fill = g),
        width = 0.15, offset = 0.03,
        axis.params = list(axis = "x", text.size = RING_TXT_MM, text.angle = -90,
                           hjust = 0.5, vjust = 0.5)) +
      ggplot2::scale_fill_manual(values = md_pal(CFG, lv, f), name = f,
                                 guide = md_guide_flex(length(lv), "bottom"))
  }

  # 最外环：丰度条形环（单色）
  p <- p + ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = abund_df, geom = geom_bar,
      mapping = ggplot2::aes(y = label, x = Abundance),
      fill = RING_BAR_COL, colour = RING_BAR_COL,
      stat = "identity", pwidth = 0.22, offset = 0.05, orientation = "y",
      axis.params = list(axis = "x", text.size = RING_TXT_MM, text.angle = -90,
                         hjust = 0.5, vjust = 0.5, nbreak = 3),
      grid.params = list())

  # 树面板收尾：**本模块是独立成图，必须自己出图例**（原 01_Fig01.R 的面板内写
  #   legend.position="none"，是因为整页 Fig1 底部另有一条共用图例行提供图例；
  #   本包按模块拆图后没有那条行，所以这里改为底部显示各环图例）。
  p <- p + ggplot2::theme(legend.position = "bottom",
                          legend.direction = "horizontal",
                          legend.box = "vertical",
                          legend.key.size = grid::unit(0.28, "cm"),
                          plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"))
  md_save_figure(CFG, p, "20_phylo_prep", height = 185)

  # ---- 8. 写表 ----
  summary_df <- data.frame(
    item = c("原始 tip 数", "修剪后 tip 数", "生根方式", "是否 rooted", "是否二叉",
             "多分支节点数（已处理）", "枝长总和", "ASV 覆盖率", "环层级", "环分类数"),
    value = c(raw_tip, length(tr$tip.label), root_used, as.character(ape::is.rooted(tr)),
              as.character(ape::is.binary(tr)), n_poly,
              if (is.null(tr$edge.length)) "无枝长" else md_fmt_num(sum(tr$edge.length), 3),
              md_fmt_pct(cov_asv), ring_rank %||% "—",
              if (is.null(ring_df)) "—" else length(levels(ring_df$taxon))),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "20_phylo_prep", list(
    Tree_summary = summary_df,
    Tip_alignment = align_df,
    Taxon_ring = if (!is.null(ring_df)) as.data.frame(table(ring_df$taxon)) else NULL,
    Note = data.frame(note = c(
      ring_note,
      paste0("系统发育指标的 ASV 覆盖率 = ", md_fmt_pct(cov_asv),
             if (cov_asv < 0.95) "（< 95%，属于子集口径，须在表注披露）" else ""),
      "生根方式影响 NRI/NTI 与系统发育 β 的数值，但本包不在此处做任何统计。"),
      stringsAsFactors = FALSE)),
    meta_extra = list(n_tip = length(tr$tip.label), root_method = root_used,
                      coverage = round(cov_asv, 4)))

  md_result(module, "OK", tables = basename(tp), figures = "figures/20_phylo_prep.pdf",
            notes = paste0(length(tr$tip.label), " tip | ", root_used,
                           " | 覆盖 ", md_fmt_pct(cov_asv)))
}

#' 去掉 Newick tip 名两端的引号/空白（ape 有时会保留）
trimnames <- function(x) trimws(gsub("^'|'$", "", gsub('^"|"$', "", x)))

#' 每个 ASV 在因子各水平上的平均相对丰度，取最高的一组（用于树的分组环）。
#' 口径照抄 01_Fig01.R §A 1.5/1.6（L173-211）的 veg_max / zone_max：
#'   rowMeans(rel[samples_of_level, ]) -> which.max -> 该水平名（并列取第一个）。
local20_dominant <- function(rel, meta, ids, var, levels) {
  mat <- vapply(levels, function(lv) {
    cols <- ids[!is.na(meta[[var]]) & as.character(meta[[var]]) == lv]
    cols <- intersect(cols, colnames(rel))
    if (!length(cols)) return(rep(0, nrow(rel)))
    rowMeans(rel[, cols, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(rel)))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
  mat[!is.finite(mat)] <- 0
  levels[max.col(mat, ties.method = "first")]
}
