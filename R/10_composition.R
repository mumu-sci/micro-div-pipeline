# =============================================================================
# R/10_composition.R —— 模块 10：类群组成（任意层级）
# -----------------------------------------------------------------------------
# 导出：md_composition(CFG)
#
# 口径（严格按 CFG$normalize$by_module，**不写死**）：
#   组成堆叠柱 = 未抽平口径 md_core_asv(CFG, "composition")（原 原始研究口径）
#   类群聚合 = 按 CFG$composition$rank（默认 Phylum）对 ASV 求和后再逐样本归一
#   Top-N    = 按"样本平均相对丰度"降序取 CFG$composition$top_n（默认 15），其余并入 Other
#   四象限   = 平均相对丰度(%) × 出现率(%) 双阈值，默认 0.1 × 10（同 nifH classify_2x2）
#              **只出表、不出图**（图已移交模块 12_niche_breadth，见下方说明）
#
# 出图：figures/10_composition.pdf —— 一排 = 各分类层级的组成堆叠柱（A = 门 / B = 属 / …）
#                                   其下 = 指示种（仅 composition$indicator$enabled = TRUE）
#       ⚠ 2026-09-26（第六轮）：原 B 面板「丰度×出现率四象限」已移入 12_niche_breadth
#         （原实现）。
# 出表：results/10_composition.xlsx（含 ASV_quadrant / Quadrant_counts / Note）
#
# 依赖共享库：md_core_asv / md_pal / md_bar_stack / md_quadrant / md_compose /
#             md_save_figure / md_write_xlsx / pt2gg / md_result / md_need
# =============================================================================

# ---------------------------------------------------------------------------
# 模块内小工具（local_ 前缀，不进共享库）
# ---------------------------------------------------------------------------
local_glue <- function(parts, sep = "_") {
  if (length(parts) == 0) return(character(0))
  parts <- lapply(parts, function(x) as.character(x))
  n <- max(vapply(parts, length, integer(1)))
  parts <- lapply(parts, function(x) rep(x, length.out = n))
  do.call(paste, c(parts, list(sep = sep)))
}

#' 取某个因子的有序水平（先按 CFG$design$ordered_levels，再按数据出现顺序）
local_factor_levels <- function(CFG, meta, var) {
  ol <- CFG$design$ordered_levels[[var]]
  obs <- unique(as.character(meta[[var]]))
  obs <- obs[!is.na(obs)]
  if (!is.null(ol)) return(intersect(as.character(ol), obs))
  obs
}

#' 用若干因子列拼出一个分组列（levels 与拼接串一一对应）
local_group_var <- function(CFG, meta, factors = NULL, name = NULL) {
  factors <- factors %||% CFG$design$factors
  factors <- intersect(factors, names(meta))
  if (length(factors) == 0) return(NULL)
  lv <- lapply(factors, function(f) local_factor_levels(CFG, meta, f))
  combos <- expand.grid(lv, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # ⚠ 不要对 combos 的列做 rev()：levels（水平名）与 vals（取值）必须用**同一个因子顺序**
  #   拼接，否则 "Bare_Coastal" 与 "Coastal_Bare" 对不上，factor() 全变 NA，图上分组消失。
  labs <- local_glue(combos, sep = "_")
  vals <- local_glue(lapply(factors, function(f) meta[[f]]), sep = "_")
  names(vals) <- as.character(meta[[CFG$data$sample_col]] %||% seq_along(vals))
  list(vals = vals, levels = labs, factors = factors,
       name = name %||% paste(factors, collapse = "_"))
}

#' 按分类层级聚合后逐样本归一为百分比（ASV × 样本 -> taxon × 样本）
local_taxon_rel <- function(m, tax, rank, id_col = "asv_id") {
  g <- as.character(tax[[rank]][match(rownames(m), tax[[id_col]])])
  g[is.na(g) | !nzchar(g)] <- "Unclassified"
  agg <- rowsum(m, group = g, reorder = TRUE)
  rel <- md_rel(agg) * 100
  rel[!is.finite(rel)] <- 0
  rel
}

#' ASV 级相对丰度（百分比，逐样本归一；rownames 保持 asv_id —— 四象限用）
local_asv_rel <- function(m) {
  rel <- md_rel(m) * 100
  rel[!is.finite(rel)] <- 0
  rel
}

#' 堆叠柱的长表：每个分组内先按样本取均值，再按 taxon 求和（= 组平均相对丰度）
local_stack_by_group <- function(rel, gv) {
  ids <- colnames(rel)
  if (is.null(gv)) {
    g <- rep("All samples", length(ids))
    lv <- "All samples"
  } else {
    g <- as.character(gv$vals[match(ids, names(gv$vals))])
    lv <- gv$levels
  }
  keep <- !is.na(g)
  out <- lapply(unique(g[keep]), function(gi) {
    cols <- ids[keep & g == gi]
    if (length(cols) == 0) return(NULL)
    sub <- rel[, cols, drop = FALSE]
    data.frame(group = gi, taxon = rownames(sub), value = rowMeans(sub, na.rm = TRUE),
               n_sample = length(cols), stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  if (is.null(d)) return(NULL)
  d$group <- factor(d$group, levels = intersect(lv, unique(d$group)))
  d
}

#' 丰度 × 出现率四象限分类（与 nifH classify_2x2 同口径）
local_classify_2x2 <- function(rel, abund_t = 0.1, freq_t = 10, asv_stats = NULL) {
  mean_ra <- rowMeans(rel, na.rm = TRUE)
  occ_pct <- rowSums(rel > 0.01) / ncol(rel) * 100
  hi_a <- mean_ra > abund_t
  hi_f <- occ_pct > freq_t
  cell <- ifelse(hi_a & hi_f, "Core",
          ifelse(hi_a & !hi_f, "Conditionally_rare",
          ifelse(!hi_a & hi_f, "Permanently_rare", "Transient_rare")))
  d <- data.frame(asv_id = rownames(rel), mean_rel_pct = mean_ra,
                  occurrence_pct = occ_pct, cell = cell, stringsAsFactors = FALSE)
  if (!is.null(asv_stats)) {
    i <- match(d$asv_id, asv_stats$asv_id)
    d$total_reads <- as.numeric(asv_stats$total_reads[i])
    d$n_present   <- as.integer(asv_stats$n_present[i])
  }
  d$cell <- factor(d$cell, levels = c("Core", "Conditionally_rare",
                                      "Permanently_rare", "Transient_rare"))
  d$interpretation <- if ("n_present" %in% names(d)) {
    ifelse(d$n_present <= 1, "singleton (1 sample)",
           ifelse(d$n_present <= 2, "doubleton (<=2 samples)", ""))
  } else ""
  d[order(-d$mean_rel_pct), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 模块主体
# ---------------------------------------------------------------------------
md_composition <- function(CFG) {
  module <- "10"
  if (!md_need(CFG, "core", module)) return(md_result(module, "SKIP", "缺少核心绘图包"))

  raw <- md_core_asv(CFG, "composition")            # ASV × 样本（按 by_module 口径）
  core <- md_core(CFG)
  tax <- core$tax
  meta <- md_core_meta(CFG)
  smp <- CFG$data$sample_col
  id_col <- CFG$data$id_col %||% "asv_id"

  if (is.null(raw) || nrow(raw) == 0 || ncol(raw) == 0)
    return(md_result(module, "SKIP", "组成矩阵为空"))

  # meta 行序严格对齐矩阵列序
  meta <- meta[match(colnames(raw), meta[[smp]]), , drop = FALSE]
  if (anyNA(meta[[smp]])) meta[[smp]] <- colnames(raw)[is.na(meta[[smp]])]

  # ---- 参数（全部来自 CFG，%||% 给默认）----
  rank_all <- CFG$data$tax_ranks %||% c("Phylum")
  rank <- CFG$composition$rank %||% "Phylum"
  top_n <- as.integer(CFG$composition$top_n %||% 15)
  abund_t <- CFG$composition$abund_threshold %||% 0.1
  freq_t  <- CFG$composition$freq_threshold %||% 10
  # 四象限配色一律走 12_plot_theme.R 的 md_pal()（= MD_BUILTIN_PAL$quadrant，
  # 逐字等于 nifH 原脚本 L823-826 的 cols_2x2）；模块内不写死颜色。
  cells_cols <- CFG$composition$cell_colors %||%
    md_pal(CFG, c("Core", "Permanently_rare", "Conditionally_rare", "Transient_rare"))

  if (is.null(tax) || !rank %in% names(tax)) {
    alt <- intersect(rank_all, names(tax))
    msg <- paste0("taxonomy 表里没有层级 '", rank, "'",
                  if (length(alt)) paste0("；回退到 '", alt[1], "'") else "；无可用层级")
    md_warn(CFG, msg)
    if (length(alt) == 0) return(md_result(module, "SKIP", msg))
    rank <- alt[1]
  }

  # ---- 1. 分类聚合：ASV -> rank 求和，再逐样本闭合到 100% ----
  #    （与 原实现 tax_glom(taxrank) -> transform_sample_counts(闭合) 同口径）
  rel_asv  <- local_asv_rel(raw)                        # ASV × 样本，%（四象限用）
  rel_tax  <- local_taxon_rel(raw, tax, rank, id_col)   # taxon × 样本，%（组成图用）
  tax_of_asv <- as.character(tax[[rank]][match(rownames(raw), tax[[id_col]])])
  tax_of_asv[is.na(tax_of_asv) | !nzchar(tax_of_asv)] <- "Unclassified"

  # ---- 2. Top-N + Other ----
  mean_ra <- rowMeans(rel_tax, na.rm = TRUE)
  ord <- order(mean_ra, decreasing = TRUE)
  # 每个层级各自的 Top-N：原脚本门用 top12（原脚本 slice_max n=12），
  # 属用 top20（原脚本 slice_max n=20）—— 数量不同直接决定
  # 属的配色分配（色板 names = rev(top20_genera)），所以必须分层给。
  top_n_rank <- as.integer(CFG$composition$top_n_by_rank[[rank]] %||% top_n)
  n_keep <- max(1L, min(top_n_rank, nrow(rel_tax)))
  top_taxa <- rownames(rel_tax)[ord][seq_len(n_keep)]
  plot_mat <- rel_tax[top_taxa, , drop = FALSE]
  rest <- setdiff(rownames(rel_tax), top_taxa)
  if (length(rest)) plot_mat <- rbind(plot_mat, Other = colSums(rel_tax[rest, , drop = FALSE]))
  # 表格/中间对象行序（本模块既有口径，不动数值）
  taxa_lv <- c(rev(top_taxa), if (length(rest)) "Other" else NULL)
  plot_mat <- plot_mat[taxa_lv, , drop = FALSE]
  # 堆叠顺序沿用同一套样式 §B L380 的
  #   factor(taxon, levels = c(top12_phyla, "Unclassified", "Other"))
  # ggplot position="stack" 把**第一个水平画在最底** -> 原实现口径是
  # 高丰度在底、Other 在顶（原模块写成 rev(top) 恰好反了）。
  stack_lv <- c(top_taxa, if (length(rest)) "Other" else NULL)

  # ---- 3. 分组变量 ----
  # 默认与箱线图一致：**前两个设计因子的交互组**（Zone_Vegetation，3×3 那种），
  # 这正是原脚本 / 原实现
  #   Group = factor(paste(Zone, Vegetation, sep = "_"), levels = both_levels)
  # 想按单因子分组就显式写 CFG$composition$group_by = "Zone"。
  g_by <- CFG$composition$group_by
  if (is.null(g_by)) {
    fd <- intersect(CFG$design$factors, names(meta))
    g_by <- if (length(fd) >= 2) fd[1:2] else if (length(fd) == 1) fd[1] else character(0)
  }
  gv <- if (length(g_by)) local_group_var(CFG, meta, g_by) else NULL
  if (is.null(gv)) {
    f1 <- intersect(CFG$design$factors, names(meta))[1]
    if (!is.na(f1) && !is.null(f1)) gv <- local_group_var(CFG, meta, f1)
  }
  if (is.null(gv)) {
    gv <- list(vals = stats::setNames(rep("All samples", nrow(meta)),
                                      as.character(meta[[smp]])),
               levels = "All samples", factors = character(0), name = "Group")
    md_log(CFG, "INFO", "没有可用分组因子，组成图按单一分组 \"All samples\" 出图")
  }
  group_var <- gv$name
  meta$.group <- gv$vals

  # ---- 4. 四象限 ----
  quad <- local_classify_2x2(rel_asv, abund_t, freq_t, md_asv_stats(raw, meta))
  cell_tab <- as.data.frame(table(cell = quad$cell), stringsAsFactors = FALSE)
  names(cell_tab) <- c("cell", "n_asv")
  cell_tab$pct_asv <- 100 * cell_tab$n_asv / sum(cell_tab$n_asv)

  # ---- 5. 堆叠柱数据 ----
  stack_df <- local_stack_by_group(plot_mat, gv)
  if (is.null(stack_df) || nrow(stack_df) == 0) {
    md_log(CFG, "WARN", "组成堆叠柱数据为空，出占位图")
    stack_df <- data.frame(group = factor("NA"), taxon = factor("NA"),
                           value = 0, n_sample = 0)
  }
  stack_df$taxon <- factor(stack_df$taxon, levels = stack_lv)
  stack_df$group <- factor(as.character(stack_df$group), levels = gv$levels)

  # 分类单元配色：**逐项照抄原脚本**（门 = Set3 插值 + Unclassified/Other，
  # 属 = 20 色硬编码 + Others），不再走 md_pal() 的兜底定性板
  other_col <- CFG$figure$palette$other %||%
    (if (identical(rank, "Genus")) MD_TAXON_OTHER else MD_TAXON_PHYLUM_OTHER)
  # 门配色复刻原脚本的关键细节：色板**名字**按 ASV 数量降序（top_phyla），
  # 堆叠**顺序**按平均丰度降序（top12_phyla）。scale_fill_manual 按名字取色，
  # 所以要传 color_order 才能把颜色对到正确类群上。
  cnt_all <- sort(table(tax_of_asv), decreasing = TRUE)
  color_order <- names(cnt_all)[seq_len(min(top_n_rank, max(1L, length(cnt_all))))]
  pal <- md_taxon_pal(CFG, rank, top_taxa, other_name = "Other", color_order = color_order)
  pal["Other"] <- other_col

  # ---- 6. 两个面板 ----
  # A. 组成堆叠柱：逐行沿用同一套样式 §B（L393-405）与 §C（L465-477）、
  #    属 Top-N 同款（原脚本 L285-296）——
  #      geom_bar(stat="identity", position="stack", width=0.65, color="white", linewidth=0.2)
  #      scale_y_continuous(expand=c(0,0), limits=c(0,101), labels = "<x>%")
  #      labs(y = "Relative Abundance (%)")  / 图例 guide_flex(n, "bottom") 置底
  #    基元 md_bar_stack() 无白描边、无百分比轴、无固定 y 上限，故按铁律 3 照抄原脚本那几行。
  #    x 轴分组名沿用原图箱线图（原脚本 L1672）的 "Zone\nVegetation" 两行写法
  #    （原 §B/C 置 element_blank() 是因为同列最下方的 E 面板已提供同一套标签）。
  pA <- ggplot2::ggplot(stack_df, ggplot2::aes(x = .data$group, y = .data$value,
                                               fill = .data$taxon)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.65,
                      colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = pal, name = rank,
                               guide = md_guide_flex(length(pal), "bottom")) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 101),
                                labels = function(x) paste0(x, "%")) +
    ggplot2::scale_x_discrete(labels = function(v) gsub("_", "\n", v, fixed = TRUE)) +
    ggplot2::labs(x = NULL, y = "Relative Abundance (%)", tag = "A") +
    theme_md(CFG = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = MD_PT$tick),
                   legend.position = "bottom", legend.box = "vertical")

  # A2..An. **额外分类层级**的组成堆叠柱（只出图、不进表 → 结果表零变化）
  #   为什么要这个：你论文的 Fig.1B（门组成）与 Fig.S3（属组成）是**两张图**；
  #   `composition$rank` 只能给一个层级，只出单个层级就会漏掉其中一张。
  #   层级由 CFG$composition$extra_ranks 指定（默认 "Genus"）；taxonomy 表里没有就跳过。
  #   画法与主面板 A 逐行一致（同 原脚本 §B / 原实现属 Top-N 口径）。
  extra_ranks <- setdiff(CFG$composition$extra_ranks %||% "Genus", rank)
  comp_panels <- list(A = pA)
  for (rk in extra_ranks) {
    if (!rk %in% names(tax)) {
      md_warn(CFG, "组成图：taxonomy 表里没有层级 '", rk, "'，跳过该面板")
      next
    }
    rel_r <- local_taxon_rel(raw, tax, rk, id_col)
    if (is.null(rel_r) || !nrow(rel_r)) next
    mr_r <- rowMeans(rel_r, na.rm = TRUE)
    ord_r <- order(mr_r, decreasing = TRUE)
    top_n_r <- as.integer(CFG$composition$top_n_by_rank[[rk]] %||% top_n)
    nk_r <- max(1L, min(top_n_r, nrow(rel_r)))
    top_r <- rownames(rel_r)[ord_r][seq_len(nk_r)]
    pm_r <- rel_r[top_r, , drop = FALSE]
    rest_r <- setdiff(rownames(rel_r), top_r)
    if (length(rest_r))
      pm_r <- rbind(pm_r, Other = colSums(rel_r[rest_r, , drop = FALSE]))
    lv_r <- c(top_r, if (length(rest_r)) "Other" else NULL)
    pm_r <- pm_r[lv_r, , drop = FALSE]
    sd_r <- local_stack_by_group(pm_r, gv)
    if (is.null(sd_r) || !nrow(sd_r)) next
    sd_r$taxon <- factor(sd_r$taxon, levels = lv_r)
    sd_r$group <- factor(as.character(sd_r$group), levels = gv$levels)
    tax_r <- as.character(tax[[rk]][match(rownames(raw), tax[[id_col]])])
    tax_r[is.na(tax_r) | !nzchar(tax_r)] <- "Unclassified"
    cnt_r <- sort(table(tax_r), decreasing = TRUE)
    ord_r2 <- names(cnt_r)[seq_len(min(top_n_r, max(1L, length(cnt_r))))]
    pal_r <- md_taxon_pal(CFG, rk, top_r, other_name = "Other", color_order = ord_r2)
    pal_r["Other"] <- other_col
    tg_r <- LETTERS[length(comp_panels) + 1L]
    comp_panels[[tg_r]] <- ggplot2::ggplot(sd_r, ggplot2::aes(x = .data$group, y = .data$value,
                                                              fill = .data$taxon)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.65,
                        colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = pal_r, name = rk,
                                 guide = md_guide_flex(length(pal_r), "bottom")) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 101),
                                  labels = function(x) paste0(x, "%")) +
      ggplot2::scale_x_discrete(labels = function(v) gsub("_", "\n", v, fixed = TRUE)) +
      ggplot2::labs(x = NULL, y = "Relative Abundance (%)", tag = tg_r) +
      theme_md(CFG = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                         size = MD_PT$tick),
                     legend.position = "bottom", legend.box = "vertical")
  }
  n_comp <- length(comp_panels)
  # ---- 2026-09-26（第六轮）：四象限面板**移交给模块 12_niche_breadth** ----
  #   原图下排左格就是「丰度 x 出现率四象限散点」
  #   （原实现 p_07sa），用户要求把它挪进 12 号，
  #   使 12 号与原图 完全对应。
  #   ⚠ 本模块只是**少一格图**：四象限的**结果表**（ASV_quadrant / Quadrant_counts）
  #     照旧输出，`quad` 在下面的表段仍被使用，统计口径与数值一个字没动。
  # ---- 6b. 指示种 IndVal（**可选面板，默认关**；composition$indicator$enabled = TRUE 才出）----
  #   参考实现（indicspecies::multipatt, IndVal.g, 999 置换）
  #   注意：本段在默认配置下完全不执行 -> 已有图的字节内容不变。
  ind_cfg <- CFG$composition$indicator %||% list()
  ind_panels <- list()
  tab_ind <- NULL
  ind_notes <- character(0)
  if (isTRUE(ind_cfg$enabled)) {
    facs <- as.character(ind_cfg$factors %||% CFG$design$factors %||% character(0))
    facs <- intersect(facs, names(meta))
    facs <- facs[vapply(facs, function(f) length(unique(stats::na.omit(meta[[f]]))) >= 2,
                        logical(1))]
    if (length(facs) > 2) {
      ind_notes <- c(ind_notes, paste0("indicator: 设计因子多于 2 个，只出前 2 个（",
                                       paste(facs[1:2], collapse = ", "), "）"))
      facs <- facs[1:2]
    }
    ind_rank  <- as.character(ind_cfg$rank %||% "ASV")
    ind_func  <- as.character(ind_cfg$func %||% "IndVal.g")
    ind_nperm <- as.integer(ind_cfg$nperm %||% 999)
    ind_adj   <- as.character(ind_cfg$adjust %||% "BH")
    ind_topn  <- as.integer(ind_cfg$top_n %||% 20)
    ind_thr   <- as.numeric(ind_cfg$fdr_thr %||% 0.05)
    ind_st    <- ind_cfg$style %||% list()

    # 在 ASV 级用 rel_asv（未抽平，与组成口径一致）；其它层级先用 local_taxon_rel 聚合
    if (identical(ind_rank, "ASV")) {
      comm_i <- rel_asv
      lab_gen <- as.character(tax[["Genus"]][match(rownames(rel_asv), tax[[id_col]])])
      lab_phy <- as.character(tax[["Phylum"]][match(rownames(rel_asv), tax[[id_col]])])
    } else if (ind_rank %in% names(tax)) {
      comm_i <- local_taxon_rel(raw, tax, ind_rank, id_col)
      lab_gen <- rownames(comm_i); lab_phy <- ""
    } else {
      comm_i <- NULL
      ind_notes <- c(ind_notes, paste0("indicator: taxonomy 无层级 '", ind_rank, "'，跳过"))
    }

    if (!is.null(comm_i) && nrow(comm_i) > 1) {
      for (f in facs) {
        iv <- tryCatch(md_indval(CFG, comm_i, meta[[f]], func = ind_func,
                                 nperm = ind_nperm, adjust = ind_adj),
                       error = function(e) NULL)
        if (is.null(iv) || !nrow(iv)) {
          ind_notes <- c(ind_notes, paste0("indicator(", f, "): 无可用结果（样本量或水平数不足）"))
          next
        }
        iv$Indicator_factor <- f
        # 标签：ASV 级 = "Genus\n(Phylum)"（原脚本写法），拿不到就用 ASV_ID
        if (identical(ind_rank, "ASV")) {
          gg <- lab_gen[match(iv$ASV_ID, rownames(comm_i))]
          pp <- lab_phy[match(iv$ASV_ID, rownames(comm_i))]
          lb <- ifelse(is.na(gg) | !nzchar(gg) | gg == "Unclassified",
                       iv$ASV_ID, paste0(gg, "\n(", ifelse(is.na(pp) | !nzchar(pp), "Unclassified", pp), ")"))
        } else {
          lb <- iv$ASV_ID
        }
        iv$Label <- lb
        tab_ind <- rbind(tab_ind, iv)
        keep <- iv[is.finite(iv$fdr) & iv$fdr < ind_thr, , drop = FALSE]
        if (!nrow(keep)) {
          ind_notes <- c(ind_notes, paste0("indicator(", f, "): 无 fdr < ", ind_thr, " 的指示种"))
          next
        }
        keep <- keep[order(-keep$IndVal), , drop = FALSE]
        keep <- utils::head(keep, ind_topn)
        keep <- keep[order(keep$Level, -keep$IndVal), , drop = FALSE]
        lv_i <- unique(as.character(keep$Level))
        pl <- md_lollipop(keep, label = "Label", indval = "IndVal", level = "Level",
                          sig = "sig", colors = md_pal(CFG, lv_i, name = f),
                          title = paste0(f, " indicator taxa"),
                          xlab = "Indicator Value (IndVal)",
                          pt = md_pt(CFG, comp = TRUE), tag = NULL,
                          text_size = ind_st$text_size %||% 2,
                          label_size_pt = ind_st$label_size_pt %||% 4.3,
                          point_size = ind_st$point_size %||% 1.8)
        ind_panels[[f]] <- pl
      }
    }
  }

  # 组图：一排 = 各分类层级的组成堆叠柱；其下 = 指示种（若有）。
  #   ⚠ 2026-09-26（第六轮）：四象限面板**已移交模块 12_niche_breadth**
  #     （原图下排左格就是它）。
  #     本模块只少一格图；四象限的结果表（ASV_quadrant / Quadrant_counts）照旧输出。
  all_plots <- comp_panels
  n_col <- n_comp
  used <- names(comp_panels)
  next_letters <- setdiff(LETTERS, used)
  for (i in seq_along(ind_panels)) {
    nm <- paste0("ind", i)
    all_plots[[nm]] <- ind_panels[[i]]
    names(all_plots)[names(all_plots) == nm] <- next_letters[i]
  }
  ind_tags <- next_letters[seq_along(ind_panels)]

  rows <- paste(LETTERS[seq_len(n_comp)], collapse = "")
  if (length(ind_tags)) {
    rows <- c(rows, vapply(ind_tags, function(t) paste(rep(t, n_col), collapse = ""),
                           character(1)))
  }
  design <- paste(rows, collapse = "\n")
  fig_h <- 150 + 70 * length(ind_tags)
  fig <- md_compose(CFG, all_plots, design = design)
  md_save_figure(CFG, fig, "10_composition", height = fig_h)
  for (nm in names(comp_panels)) md_maybe_single(CFG, comp_panels[[nm]], paste0("10_composition_", nm))

  # ---- 7. 表 ----
  tab_top <- data.frame(
    rank = rank,
    taxon = taxa_lv,
    mean_rel_pct_all = c(round(mean_ra[top_taxa], CFG$table$digits %||% 3),
                         if (length(rest)) round(sum(mean_ra[rest]), CFG$table$digits %||% 3) else NULL),
    n_asv = c(vapply(top_taxa, function(t) sum(tax_of_asv == t), integer(1)),
              if (length(rest)) sum(!tax_of_asv %in% top_taxa) else NULL),
    role = c(rep("Top", length(top_taxa)), if (length(rest)) "Other" else NULL),
    stringsAsFactors = FALSE)

  tab_group <- data.frame(
    rank = rank,
    factor = group_var,
    group = as.character(stack_df$group),
    taxon = as.character(stack_df$taxon),
    mean_rel_pct = round(stack_df$value, CFG$table$digits %||% 3),
    n_sample = stack_df$n_sample,
    stringsAsFactors = FALSE)

  tab_cell <- quad
  tab_cell$mean_rel_pct <- round(tab_cell$mean_rel_pct, 4)
  tab_cell$occurrence_pct <- round(tab_cell$occurrence_pct, 2)

    note <- data.frame(
    item = c("Taxonomic rank", "Top N", "Others pooled", "Abundance threshold (%)",
             "Occurrence threshold (%)", "Normalization (composition)",
             "Aggregation", "Grouping factor", "Bars show",
             "Core taxa note", "ASV classification unit", "n_taxon_summed"),
    value = c(rank, as.character(n_keep), if (length(rest)) "yes" else "no",
              as.character(abund_t), as.character(freq_t),
              as.character(CFG$normalize$by_module$composition %||% CFG$normalize$method),
              paste0("ASV -> ", rank, " rowsum, then per-sample closure to 100%"),
              group_var,
              "mean per-sample relative abundance within each group (stack sums to ~100%)",
              paste0("n = ", sum(quad$cell == "Core"), " ASV. This dataset is extremely ",
                     "singleton-dominated, so high-abundance ASVs are mostly ASVs that ",
                     "dominate the 1-3 samples they appear in (Conditionally_rare = ",
                     sum(quad$cell == "Conditionally_rare"), "), while Transient_rare = ",
                     sum(quad$cell == "Transient_rare"),
                     ". 'Core is nearly empty' is a real feature of the data, not a bug."),
              "ASV x sample (relative abundance in %, thresholds applied to ASV level)",
              as.character(nrow(rel_tax))),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "10_composition", list(
    Taxa_top = tab_top,
    Group_means = tab_group,
    ASV_quadrant = tab_cell,
    Quadrant_counts = cell_tab,
    Indicator_species = tab_ind,
    Notes = note,
    ASV_quadrant_detail = if (isTRUE(CFG$table$detail)) tab_cell else NULL))

  md_result(module, "OK",
            tables = basename(tp),
            figures = "figures/10_composition.pdf",
            notes = c(paste0("rank=", rank, "; top_n=", n_keep,
                             "; 核心类群 ", sum(quad$cell == "Core"), " 个 ASV（真实特征）"),
                      ind_notes))
}
