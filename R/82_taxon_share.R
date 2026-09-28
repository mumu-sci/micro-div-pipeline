# =============================================================================
# R/82_taxon_share.R ---- 模块 82：类群共享（UpSet）
# -----------------------------------------------------------------------------
# 源：毕业论文 S4_supplementary_taxa.R 第 296-690 行（原 8_asv_sharing.R）
#     关键原文行：360-428（get_group_asvs / build_shared_df 的共享-特有统计）、
#                 489-560（make_upset_hm：ComplexHeatmap::make_comb_mat + UpSet）
#
# 做什么
#   按 CFG$taxon_share$factor（NULL -> design$factors[1]）的水平把样本分组，
#   每组取"至少在 min_prev 个样本中出现"的 rank 层类群集合，做交集统计。
#   A  UpSet（ComplexHeatmap::UpSet，ggplotify::as.ggplot() 包装后拼 patchwork）
#   B  共享/独有数量条形；水平数 = 2 或 3 且 show_venn = TRUE 时
#      B 改为 md_venn2()（2 水平）/ md_venn3()（3 水平）
#   表 = Set_size / Intersection / Sharing_stats / Note
#
# 依赖共享库：md_core / md_core_asv / md_rel / md_pal / md_compose / md_venn2 /
#             md_venn3 / md_save_figure / md_write_xlsx / md_empty_plot / md_result
#
# 说明：本模块**只出 UpSet 一个面板**（原脚本 FigS3b 就是一张 UpSet）。
#   2026-09-26（第六轮）按用户裁定删掉原 B 面板（Venn / 共享条形）；
#   show_venn 配置键保留但不再生效。单面板不进 patchwork，因此也不需要
#   "定宽高比矢量包装" 那一层 —— 图是**全矢量**的，字号走 pt$tick / pt$legend。
# =============================================================================

#' 模块 82 入口
#' @param CFG config 列表（见 config/config_template.R 的 taxon_share 块）
md_taxon_share <- function(CFG) {
  module <- "82"

  if (!isTRUE(CFG$taxon_share$enabled %||% TRUE))
    return(md_result(module, "SKIP", "taxon_share$enabled = FALSE"))
  if (!md_need(CFG, c("upset"), module))
    return(md_result(module, "SKIP", "缺少依赖包（ComplexHeatmap / circlize / ggplotify）"))

  core <- tryCatch(md_core(CFG), error = function(e) NULL)
  if (is.null(core) || is.null(core$asv_raw)) {
    md_skip(CFG, module, "取不到 cache/05_core.rds")
    return(md_result(module, "SKIP", "缺 05 产物"))
  }

  ts       <- CFG$taxon_share
  rank     <- ts$rank %||% "ASV"
  min_prev <- as.integer(ts$min_prev %||% 1)
  max_comb <- as.integer(ts$max_comb %||% 20)
  top_n    <- as.integer(ts$upset_top_n %||% 20)
  show_venn <- isTRUE(ts$show_venn %||% TRUE)
  st       <- ts$style %||% list()
  pt       <- md_pt(CFG, comp = FALSE)
  id_col   <- CFG$data$id_col %||% "asv_id"
  smp      <- CFG$data$sample_col %||% "sample_id"

  # ---- 组合口径矩阵（默认 composition 口径：未抽平，与 10 / 81 一致）----
  m <- md_core_asv(CFG, ts$normalize_key %||% "composition")
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
    md_skip(CFG, module, "类群共享输入矩阵为空")
    return(md_result(module, "SKIP", "矩阵为空"))
  }
  m <- as.matrix(m)
  storage.mode(m) <- "double"

  meta <- core$meta
  if (is.null(meta)) {
    md_skip(CFG, module, "取不到样本元数据")
    return(md_result(module, "SKIP", "缺 meta"))
  }
  meta <- meta[match(colnames(m), meta[[smp]]), , drop = FALSE]

  # ---- 分组因子（NULL -> design$factors[1]）----
  fac <- ts$factor %||% (CFG$design$factors %||% character(0))[1]
  if (is.null(fac) || !fac %in% names(meta)) {
    msg <- paste0("分组因子 '", fac %||% "NULL", "' 不在样本元数据里")
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }
  ol  <- CFG$design$ordered_levels[[fac]]
  obs <- unique(as.character(meta[[fac]]))
  obs <- obs[!is.na(obs)]
  lv  <- if (is.null(ol)) obs else intersect(as.character(ol), obs)
  if (length(lv) == 0) lv <- obs
  if (length(lv) < 2) {
    msg <- paste0("因子 '", fac, "' 只有 ", length(lv), " 个有效水平，无法做共享分析")
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }

  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  # ===========================================================================
  # 分析段 —— 整体包进 md_stage()，图与表在段外
  # ===========================================================================
  md_stage(CFG, "82_analysis", {

    # ---- 1. 把 ASV 聚合到 rank 层（ASV 层直接用原表）----
    if (identical(toupper(rank), "ASV")) {
      mat_rank <- m
      tax_label <- rownames(m)
    } else {
      if (is.null(core$tax) || !rank %in% names(core$tax)) {
        mat_rank <- NULL
      } else {
        g <- as.character(core$tax[[rank]][match(rownames(m), core$tax[[id_col]])])
        g[is.na(g) | !nzchar(g)] <- "Unclassified"
        mat_rank <- rowsum(m, group = g, reorder = TRUE)
      }
      tax_label <- if (is.null(mat_rank)) character(0) else rownames(mat_rank)
    }

    if (is.null(mat_rank) || nrow(mat_rank) == 0) {
      sets <- list(); sets_bin <- list(); present_mat <- NULL
    } else {
      present_mat <- (mat_rank > 0)
      # 每组：至少在 min_prev 个样本中出现的类群（原文 360-369 的 get_group_asvs）
      sets <- lapply(lv, function(g) {
        idx <- which(as.character(meta[[fac]]) == g)
        if (length(idx) == 0) return(character(0))
        keep <- rowSums(present_mat[, idx, drop = FALSE]) >= min_prev
        rownames(mat_rank)[keep]
      })
      names(sets) <- lv
      # 稀疏组的"有无"集合（用于共享/独有表）
      sets_bin <- lapply(lv, function(g) {
        idx <- which(as.character(meta[[fac]]) == g)
        if (length(idx) == 0) return(character(0))
        rownames(mat_rank)[rowSums(present_mat[, idx, drop = FALSE]) > 0]
      })
      names(sets_bin) <- lv
    }
    sets <- sets[lengths(sets) > 0]
    sets_bin <- sets_bin[lengths(sets_bin) > 0]

    # ---- 2. 共享 / 独有统计（原文 403-454）----
    uniq_list <- lapply(names(sets_bin), function(nm) {
      setdiff(sets_bin[[nm]], unlist(sets_bin[names(sets_bin) != nm]))
    })
    names(uniq_list) <- names(sets_bin)
    shared_all <- if (length(sets_bin) >= 2) Reduce(intersect, sets_bin) else character(0)

    set_size <- data.frame(
      Factor = fac, Level = names(sets_bin), N_taxon = as.integer(lengths(sets_bin)),
      N_taxon_min_prev = as.integer(lengths(sets)[names(sets_bin)]),
      N_unique = as.integer(lengths(uniq_list)), stringsAsFactors = FALSE)
    set_size$Pct_of_pool <- 100 * set_size$N_taxon / max(1L, nrow(mat_rank))

    sharing_stats <- do.call(rbind, lapply(names(sets_bin), function(nm) {
      data.frame(Factor = fac, Level = nm, Type = "Unique", N_taxon = length(uniq_list[[nm]]),
                 stringsAsFactors = FALSE)
    }))
    if (length(shared_all) > 0) {
      sharing_stats <- rbind(sharing_stats, data.frame(
        Factor = fac, Level = paste0("Shared_", fac), Type = "Shared_All",
        N_taxon = length(shared_all), stringsAsFactors = FALSE))
    }
    sharing_stats <- rbind(sharing_stats, data.frame(
      Factor = fac, Level = names(sets_bin), Type = "Total",
      N_taxon = as.integer(lengths(sets_bin)), stringsAsFactors = FALSE))

    # ---- 3. UpSet 组合矩阵（原文 489-582 的 make_upset_hm 口径）----
    comb_mat <- NULL; int_tab <- NULL
    if (length(sets) >= 2 &&
        requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      comb_mat <- tryCatch(
        ComplexHeatmap::make_comb_mat(sets, mode = "distinct"),
        error = function(e) NULL)
      if (!is.null(comb_mat)) {
        cs <- ComplexHeatmap::comb_size(comb_mat)
        comb_mat <- comb_mat[, which(cs > 1), drop = FALSE]
        if (ncol(comb_mat) > 0) {
          ord <- order(ComplexHeatmap::comb_size(comb_mat), decreasing = TRUE)
          comb_mat <- comb_mat[, ord, drop = FALSE]
          if (ncol(comb_mat) > max_comb) comb_mat <- comb_mat[, seq_len(max_comb), drop = FALSE]
        }
      }
      if (!is.null(comb_mat) && ncol(comb_mat) > 0) {
        cs <- ComplexHeatmap::comb_size(comb_mat)
        int_tab <- data.frame(
          Intersection = names(cs),
          N_taxon = as.integer(cs),
          N_levels = vapply(ComplexHeatmap::comb_degree(comb_mat), as.integer, integer(1)),
          stringsAsFactors = FALSE)
        int_tab$Pct_of_pool <- 100 * int_tab$N_taxon / max(1L, nrow(mat_rank))
      } else {
        comb_mat <- NULL
      }
    }

    # ---- 4. 各类群出现率明细（原文 428-460 的 ASV_list 列在这里按类群给出现率）----
    prev_by_level <- NULL
    if (!is.null(present_mat)) {
      prev_by_level <- do.call(rbind, lapply(lv, function(g) {
        idx <- which(as.character(meta[[fac]]) == g)
        if (length(idx) == 0) return(NULL)
        n_in <- rowSums(present_mat[, idx, drop = FALSE])
        data.frame(Factor = fac, Level = g, Taxon = rownames(mat_rank),
                   N_sample_present = as.integer(n_in),
                   N_sample_level = length(idx),
                   In_set = n_in >= min_prev,
                   stringsAsFactors = FALSE)
      }))
    }
  })

  if (is.null(present_mat) || nrow(present_mat) == 0) {
    md_skip(CFG, module, paste0("层级 '", rank, "' 在 taxonomy 表里不存在或聚合后为空"))
    return(md_result(module, "SKIP", "聚合后为空"))
  }
  if (length(sets) < 2) {
    md_skip(CFG, module, "有效水平不足 2 个（或所有分组集合为空），无法做共享分析")
    return(md_result(module, "SKIP", "有效集合不足"))
  }

  # ---- 配色：走 md_pal()（CFG$figure$palette[[fac]] -> 内置板 -> 定性板兜底）----
  set_cols <- st$set_colors %||% md_pal(CFG, names(sets_bin), name = fac)
  set_cols <- set_cols[names(sets_bin)]
  if (anyNA(set_cols)) set_cols[is.na(set_cols)] <- "grey70"
  bar_col <- st$bar_color %||% "grey35"

  # ---------------------------------------------------------------------------
  # 出图段：**只有 UpSet 一个面板**（原脚本 FigS3b 就是一张 UpSet）
  #   2026-09-26（第六轮）：按用户裁定删掉 B 面板（Venn / 共享条形）。
  #   单面板不再进 patchwork，所以也不需要"定宽高比矢量包装"那一层 ——
  #   直接把 draw() 出来的 grob 用 ggplotify::as.ggplot() 包一次即可，全矢量。
  #   字号统一走 pt$tick / pt$legend（原来用 pt$ann - 3 = 4.5 pt，明显小于其余图）。
  # ---------------------------------------------------------------------------
  n_levels <- length(sets_bin)

  if (is.null(comb_mat) || ncol(comb_mat) == 0) {
    pA <- md_empty_plot(CFG, "No shared intersection (>= 2 levels)")
    add_note("没有 >= 2 个水平共享的组合，UpSet 面板退化为占位图。")
  } else {
    pA <- tryCatch({
      ta <- ComplexHeatmap::upset_top_annotation(
        comb_mat, add_numbers = TRUE, annotation_name_side = "left",
        annotation_name_gp = grid::gpar(fontsize = pt$legend),
        axis_param = list(gp = grid::gpar(fontsize = pt$tick)),
        numbers_gp = grid::gpar(fontsize = pt$tick),
        bar_width = 0.7, height = grid::unit(3.2, "cm"),
        gp = grid::gpar(fill = bar_col, col = bar_col))
      set_sizes <- as.integer(lengths(sets_bin)[names(sets)])
      ax_at <- pretty(c(0, max(set_sizes)), n = 3)
      ax_lab <- ifelse(ax_at >= 1000,
                       paste0(format(ax_at / 1000, trim = TRUE, scientific = FALSE), "k"),
                       format(ax_at, trim = TRUE, scientific = FALSE))
      la <- ComplexHeatmap::upset_left_annotation(
        comb_mat, width = grid::unit(1.6, "cm"),
        gp = grid::gpar(fill = unname(set_cols[names(sets)])),
        axis_param = list(at = ax_at, labels = ax_lab,
                          gp = grid::gpar(fontsize = pt$tick)),
        annotation_name_gp = grid::gpar(fontsize = pt$legend))
      ht <- ComplexHeatmap::UpSet(
        comb_mat, top_annotation = ta, left_annotation = la, right_annotation = NULL,
        set_order = names(sets),
        row_names_gp = grid::gpar(fontsize = pt$legend, fontface = "bold",
                                  col = unname(set_cols[names(sets)])),
        column_names_gp = grid::gpar(fontsize = pt$legend),
        column_names_rot = if (ncol(comb_mat) > 8) 45 else 0,
        border = TRUE, pt_size = grid::unit(as.numeric(st$point_size %||% 2.2), "mm"),
        bg_col = "#F2F2F2", lwd = as.numeric(st$line_width %||% 0.6),
        show_row_names = TRUE, row_title = NULL, column_title = NULL)
      gg <- ggplotify::as.ggplot(
        md_grid_grab(ComplexHeatmap::draw(ht, newpage = TRUE))) +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
      gg + ggplot2::labs(tag = "A")
    }, error = function(e) {
      md_warn(CFG, "UpSet 绘制失败：", conditionMessage(e))
      md_empty_plot(CFG, "UpSet failed")
    })
  }

  # 单面板图：**不经过 md_compose()** —— 既不套 theme_md_comp()（会覆盖 ComplexHeatmap
  # 的字号），也不受拼图缩放影响。
  md_save_figure(CFG, pA, "82_taxon_share", width = 183, height = 110)
  md_maybe_single(CFG, pA, "82_taxon_share_A")
  add_note("本模块只出 UpSet 一个面板（原脚本 FigS3b 就是一张 UpSet）；",
           "Venn 面板已按用户裁定删除（show_venn 键保留但不再生效）。")
  # ---------------------------------------------------------------------------
  # 写表段
  # ---------------------------------------------------------------------------
  note <- data.frame(
    item = c("Factor", "Levels", "Rank", "min_prev", "max_comb", "upset_top_n",
             "show_venn (deprecated)", "UpSet mode", "Input matrix",
             "Normalization (composition)",
             "Set definition", "Panel", "Original script difference"),
    value = c(fac, paste(lv, collapse = ", "), rank, as.character(min_prev),
              as.character(max_comb), as.character(top_n),
              paste0(as.character(show_venn), " —— 该键已弃用：本模块不再出 Venn"),
              "distinct (each taxon assigned to its most specific combination)",
              paste0(nrow(mat_rank), " ", rank, " x ", ncol(mat_rank), " samples"),
              as.character(CFG$normalize$by_module$composition %||% CFG$normalize$method),
              paste0("taxa present in >= ", min_prev, " samples of a level"),
              "single UpSet panel (ComplexHeatmap::UpSet, drawn grob wrapped by ggplotify::as.ggplot)",
              paste0("The original script drew three hard-wired UpSets (Zone, Vegetation, ",
                     "Zone x Vegetation at 508 x 406 mm). Here the grouping factor is ",
                     "configurable (CFG$taxon_share$factor) and one UpSet at 183 x 110 mm ",
                     "is produced. The former Venn / shared-unique bar panel was removed ",
                     "on 2026-09-26 by user decision.")),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "82_taxon_share", list(
    Set_size      = set_size,
    Intersection  = int_tab,
    Sharing_stats = sharing_stats,
    Prevalence    = if (isTRUE(CFG$table$detail)) prev_by_level else NULL,
    Note          = note))

  md_result(module, "OK", NULL,
            tables  = "82_taxon_share.xlsx",
            figures = "figures/82_taxon_share.pdf",
            notes = c(notes,
                      paste0("factor=", fac, "; rank=", rank, "; levels=", n_levels,
                             "; min_prev=", min_prev,
                             "; pooled ", rank, " = ", nrow(mat_rank),
                             "; shared by all = ", length(shared_all))))
}
