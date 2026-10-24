# =============================================================================
# R/83_taxon_diff.R ---- 模块 83：差异类群（CLR + Kruskal-Wallis + Dunn + Cliff's delta）
# -----------------------------------------------------------------------------
# 参考实现（原脚本）
#     契约见 _build-notes/micro-div-pipeline/ROUND5_MODULE_SPEC.md
#
# 导出：md_taxon_diff(CFG)
#
# 计算（逐项照抄原文口径）：
#   1) 取 composition 口径的丰度矩阵（= 未抽平原始 reads，原脚本用 ps_raw_used）
#   2) 出现率过滤：rowSums(m > 0) >= max(3, round(min_prev_frac * n_samples))
#      （原脚本 PREVALENCE_THRESHOLD = 0.05）
#   3) md_glom() 聚合到 taxon_diff$rank（默认 Genus），空/NA 归 "Unclassified"
#   4) rel = 相对丰度（%）；clr = md_clr()（伪计数缺省 = 最小正值的 1/2）
#      ⚠ md_clr() 的 "行 = 组成单元（样本）"，而 md_diff_taxa_kw() 要"类群 × 样本"，
#        所以这里 t() 进去、t() 回来。
#   5) md_diff_taxa_kw()：KW + BH + Cliff's delta + 富集组 + 各组平均相对丰度
#   6) md_taxon_dunn()：显著类群的 Dunn 事后检验
#
# 图：figures/83_taxon_diff.pdf（183 × 150 mm）
#     A = 差异类群热图（行 = Top top_n 显著类群，列 = 分组水平，
#         值 = heatmap_scale：row_z（推荐）/ rel / clr）
#     B = （仅 show_box = TRUE）Top 差异类群的 CLR 箱线图（md_box_group()）
# 表：results/83_taxon_diff.xlsx（Top_taxa / Kruskal / Dunn / Cliff_delta / Note）
# =============================================================================

# ---------------------------------------------------------------------------
# local_* 辅助（只在本模块内使用）
# ---------------------------------------------------------------------------

#' 图上文字一律纯 ASCII（类群名是拉丁名，但分组水平名不保证；一律过一遍）
local_td_ascii <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- suppressWarnings(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"))
  x[is.na(x)] <- ""
  gsub("[^ -~]", "_", x)
}

#' 把矩阵画成热图，返回 ggplot 对象（好让 md_save_figure() 统一落盘）。
#'   engine = "ComplexHeatmap" 时用 ComplexHeatmap::Heatmap + ggplotify::as.ggplot()
#'   engine = "pheatmap"       时用 pheatmap() + ggplotify::as.ggplot()
#'   两条路都走不通时退回 geom_tile() 画法（保证一定有一张 ggplot）
local_td_heatmap <- function(CFG, mat, engine, style, title = NULL,
                             subtitle = NULL, legend_title = "Z-score") {
  mat <- as.matrix(mat)
  if (!nrow(mat) || !ncol(mat)) return(NULL)
  row_pt <- style$row_label_pt %||% 5
  col_pt <- style$col_label_pt %||% 7
  ramp3  <- style$heatmap_col %||% c("#2166AC", "#F7F7F7", "#B2182B")
  zr <- suppressWarnings(max(abs(range(mat, na.rm = TRUE))))
  if (!is.finite(zr) || zr <= 0) zr <- 1
  mat[!is.finite(mat)] <- 0
  tt <- if (is.null(subtitle) || !nzchar(subtitle)) title else
    paste0(title %||% "", "\n", subtitle)

  # 包装后的有效性护栏：ggplotify 在失败时可能返回一个"没有图层数据"的空壳 ggplot，
  # 之前正是它让 83 的热图渲染成一张带 x/y 数值轴的空图。
  .usable <- function(g) {
    if (is.null(g)) return(FALSE)
    ok <- tryCatch(length(ggplot2::ggplot_build(g)$data) > 0, error = function(e) FALSE)
    isTRUE(ok)
  }

  if (identical(engine, "ComplexHeatmap") &&
      requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE) &&
      requireNamespace("ggplotify", quietly = TRUE)) {
    ht <- tryCatch({
      col_fun <- circlize::colorRamp2(c(-zr, 0, zr), ramp3)
      ComplexHeatmap::Heatmap(
        mat, name = legend_title, col = col_fun,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = row_pt),
        column_names_gp = grid::gpar(fontsize = col_pt),
        column_names_rot = 0,
        border = TRUE,
        column_title = tt,
        column_title_gp = grid::gpar(fontsize = col_pt, fontface = "bold"),
        heatmap_legend_param = list(
          title_gp = grid::gpar(fontsize = row_pt),
          labels_gp = grid::gpar(fontsize = row_pt),
          direction = "vertical"))
    }, error = function(e) { md_warn(CFG, "83: ComplexHeatmap 绘制失败 —— ", conditionMessage(e)); NULL })
    if (!is.null(ht)) {
      # ⚠ `ggplotify::as.ggplot(ht)` 对 Heatmap 对象**会返回空图** —— 必须先 draw()
      #   再抓成 grob（这正是之前"热图只有两个色块 + x/y 数值轴"的根因）。
      g <- tryCatch({
        gr <- md_grid_grab(ComplexHeatmap::draw(ht, newpage = TRUE))
        ggplotify::as.ggplot(gr) +
          ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
      }, error = function(e) NULL)
      if (.usable(g)) return(g)
      md_warn(CFG, "83: ComplexHeatmap 包装后仍为空图，退回下一档渲染。")
    }
  }

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ph <- tryCatch(pheatmap::pheatmap(
      mat, cluster_rows = FALSE, cluster_cols = FALSE, silent = TRUE,
      color = grDevices::colorRampPalette(ramp3)(100),
      breaks = seq(-zr, zr, length.out = 101),
      fontsize_row = row_pt, fontsize_col = col_pt,
      border_color = "white", main = tt), error = function(e) NULL)
    if (!is.null(ph) && requireNamespace("ggplotify", quietly = TRUE)) {
      g <- tryCatch(ggplotify::as.ggplot(ph$gtable), error = function(e) NULL)
      if (.usable(g)) return(g)
    }
  }

  # 兜底：geom_tile（不依赖任何热图包；离散轴，永远不会出现 0.00-1.00 的数值轴）
  long <- data.frame(
    Taxa  = rep(rownames(mat), times = ncol(mat)),
    Group = rep(colnames(mat), each = nrow(mat)),
    Value = as.numeric(mat), stringsAsFactors = FALSE)
  long$Taxa  <- factor(long$Taxa, levels = rev(rownames(mat)))
  long$Group <- factor(long$Group, levels = colnames(mat))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$Group, y = .data$Taxa)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$Value), colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(
      colours = ramp3, limits = c(-zr, zr), name = legend_title,
      guide = ggplot2::guide_colourbar(barwidth = grid::unit(4, "mm"),
                                       barheight = grid::unit(25, "mm"))) +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_md(CFG = CFG) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = col_pt),
                   axis.text.y = ggplot2::element_text(size = row_pt),
                   axis.ticks = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = col_pt, face = "bold", hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = row_pt, colour = "grey35", hjust = 0),
                   legend.position = "right")
}

# ---------------------------------------------------------------------------
# 模块主体
# ---------------------------------------------------------------------------
md_taxon_diff <- function(CFG) {
  module <- "83"
  if (!isTRUE(CFG$taxon_diff$enabled %||% TRUE))
    return(md_result(module, "SKIP", "taxon_diff$enabled = FALSE"))

  engine <- as.character(CFG$taxon_diff$heatmap_engine %||% "ComplexHeatmap")
  have_pheatmap <- md_need(CFG, "taxon_heat", module)
  have_ch <- FALSE
  if (identical(engine, "ComplexHeatmap")) have_ch <- md_need(CFG, "upset", module)
  if (!have_pheatmap && !have_ch) {
    msg <- "缺少热图依赖（pheatmap 与 ComplexHeatmap/circlize/ggplotify 都不可用）"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }

  core <- tryCatch(md_core(CFG), error = function(e) NULL)
  if (is.null(core)) {
    msg <- "取不到 cache/05_core.rds（先跑模块 05）"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }

  notes <- character(0)
  if (identical(engine, "ComplexHeatmap") && !have_ch) {
    engine <- "pheatmap"
    notes <- c(notes, "heatmap_engine = 'ComplexHeatmap' 但缺 ComplexHeatmap/circlize/ggplotify —— 已退回 pheatmap。")
  }

  # ---- 分析段 ----
  #   td_fallback_txt：显著类群不足 min_sig 时的副标题文案，在段内赋值、段外用在图面上。
  #   先初始化为 NULL，保证缓存命中/走别的分支时它也存在。
  td_fallback_txt <- NULL
  md_stage(CFG, "83_analysis", {
    md_seed(CFG, md_offset(CFG, "taxon_diff", 8420L))
    td_notes <- character(0)
    smp <- CFG$data$sample_col

    td_rank    <- as.character(CFG$taxon_diff$rank %||% "Genus")
    td_clr     <- isTRUE(CFG$taxon_diff$clr %||% TRUE)
    td_pc      <- CFG$taxon_diff$pseudocount
    td_prev    <- CFG$taxon_diff$min_prev_frac %||% 0.05
    td_min_n   <- as.integer(CFG$taxon_diff$min_n %||% 3L)
    td_test    <- tolower(as.character(CFG$taxon_diff$test %||% "kruskal"))
    td_posthoc <- tolower(as.character(CFG$taxon_diff$posthoc %||% "dunn"))
    td_adjust  <- as.character(CFG$taxon_diff$adjust %||% "BH")
    td_effect  <- tolower(as.character(CFG$taxon_diff$effect %||% "cliff"))
    td_pthr    <- CFG$taxon_diff$p_thr %||% 0.05
    td_minsig  <- max(1L, as.integer(CFG$taxon_diff$min_sig %||% 3L))
    td_topn    <- max(td_minsig, as.integer(CFG$taxon_diff$top_n %||% 25L))
    td_scale   <- as.character(CFG$taxon_diff$heatmap_scale %||% "row_z")
    if (!identical(td_test, "kruskal"))
      td_notes <- c(td_notes, paste0(
        "test = '", td_test, "' 在本包里没有独立实现（共享 md_diff_taxa_kw() 只做 Kruskal-Wallis）—— 已按 Kruskal-Wallis 执行。"))
    if (!td_scale %in% c("row_z", "rel", "clr")) {
      td_notes <- c(td_notes, paste0("heatmap_scale = '", td_scale, "' 不可识别 —— 已退回 row_z"))
      td_scale <- "row_z"
    }

    # ---- 1. 分组变量（原文做了 Zone 与 Vegetation 两套；本包取第一个设计因子）----
    td_meta <- core$meta
    facs <- intersect(CFG$design$factors %||% character(0), names(td_meta))
    td_group <- if (length(facs)) facs[1] else NA_character_
    if (is.na(td_group)) td_notes <- c(td_notes, "design$factors 里没有可用的分组列")

    # ---- 2. 丰度矩阵（composition 口径 = 未抽平）----
    td_m <- tryCatch(md_core_asv(CFG, "composition"), error = function(e) NULL)
    td_ok <- FALSE
    td_res <- NULL; td_top <- NULL; td_dunn <- NULL
    td_mat <- NULL; td_box <- NULL; td_grp_lv <- character(0)
    td_mat_label <- ""
    td_n_taxa_before <- NA_integer_; td_n_taxa_after <- NA_integer_

    if (!is.null(td_m) && nrow(td_m) >= 2 && ncol(td_m) >= 6 && !is.na(td_group)) {
      ids <- intersect(colnames(td_m), as.character(td_meta[[smp]]))
      if (length(ids) >= 6) {
        td_m <- td_m[, ids, drop = FALSE]
        td_meta <- td_meta[match(ids, as.character(td_meta[[smp]])), , drop = FALSE]
        # 组水平（尊重 design$ordered_levels）
        ol <- CFG$design$ordered_levels[[td_group]]
        obs <- unique(as.character(td_meta[[td_group]]))
        obs <- obs[!is.na(obs)]
        td_grp_lv <- if (!is.null(ol)) intersect(as.character(ol), obs) else obs
        if (length(td_grp_lv) < 2) {
          td_notes <- c(td_notes, paste0("分组变量 '", td_group, "' 少于 2 个水平"))
        } else {
          n_s <- length(ids)
          min_ns <- max(3L, as.integer(round(td_prev * n_s)))
          td_n_taxa_before <- nrow(td_m)
          keep <- rowSums(td_m > 0) >= min_ns
          td_m <- td_m[keep, , drop = FALSE]
          td_n_taxa_after <- nrow(td_m)
          td_notes <- c(td_notes, sprintf(
            "出现率过滤：>= %d 个样本（max(3, round(%s x %d))），%d -> %d 个单元。",
            min_ns, format(td_prev), n_s, td_n_taxa_before, td_n_taxa_after))

          if (nrow(td_m) >= 2) {
            # ---- 3. 聚合到目标层级 ----
            if (!identical(td_rank, "ASV")) {
              id_col <- CFG$data$id_col %||% "asv_id"
              if (!id_col %in% names(core$tax)) id_col <- names(core$tax)[1]
              if (!td_rank %in% names(core$tax)) {
                td_notes <- c(td_notes, paste0("taxonomy 表里没有层级 '", td_rank, "' —— 已在 ASV 层级继续"))
              } else {
                gl <- tryCatch(md_glom(td_m, core$tax, td_rank, id_col), error = function(e) NULL)
                if (is.null(gl)) {
                  td_notes <- c(td_notes, paste0("md_glom() 聚合到 '", td_rank, "' 失败 —— 已在 ASV 层级继续"))
                } else {
                  td_m <- gl$mat
                  td_notes <- c(td_notes, sprintf("聚合层级 %s：%d 个类群。", td_rank, nrow(td_m)))
                }
              }
            }

            rel <- md_rel(td_m) * 100
            clr_taxa_sample <- if (isTRUE(td_clr)) {
              t(md_clr(t(td_m), pseudocount = td_pc))
            } else {
              td_notes <- c(td_notes, "clr = FALSE —— 差异检验与热图的 clr 口径改用相对丰度（%）。")
              rel
            }
            rownames(clr_taxa_sample) <- rownames(td_m)
            colnames(clr_taxa_sample) <- colnames(td_m)

            # ---- 4. KW + BH + Cliff's delta ----
            td_res <- tryCatch(
              md_diff_taxa_kw(clr_taxa_sample, rel, td_meta, td_group,
                              sample_col = smp, min_n = td_min_n, adjust = td_adjust),
              error = function(e) {
                td_notes <<- c(td_notes, paste0("md_diff_taxa_kw() 报错：", conditionMessage(e)))
                NULL
              })
            if (!is.null(td_res) && nrow(td_res)) {
              # ---- 5. 显著类群（不足 min_sig 个时回退为按 KW_padj 排序的 Top-N，并如实标注）----
              n_sig <- sum(!is.na(td_res$KW_padj) & td_res$KW_padj < td_pthr)
              if (n_sig < td_minsig) {
                td_fallback_txt <- sprintf(
                  "None/only %d taxon(s) passed FDR < %s; showing top %d by KW p",
                  n_sig, format(td_pthr), min(td_topn, nrow(td_res)))
                td_notes <- c(td_notes, sprintf(
                  "只有 %d / %d 个类群通过 KW_padj < %s（< taxon_diff$min_sig = %d）—— ",
                  n_sig, nrow(td_res), format(td_pthr), td_minsig),
                  "热图回退为按 KW_padj 排序的 Top ", min(td_topn, nrow(td_res)),
                  "，并在面板副标题与 Note 里写明；没有把不显著的类群当显著结果。")
                sig <- td_res
              } else {
                sig <- td_res[!is.na(td_res$KW_padj) & td_res$KW_padj < td_pthr, , drop = FALSE]
                td_notes <- c(td_notes, sprintf("%d / %d 个类群通过 KW_padj < %s。",
                                                n_sig, nrow(td_res), format(td_pthr)))
              }
              td_top <- utils::head(sig[order(sig$KW_padj), , drop = FALSE], td_topn)
              if (nrow(td_top) < nrow(sig))
                td_notes <- c(td_notes, sprintf("热图只显示 Top %d 个类群（top_n）。", nrow(td_top)))

              # ---- 6. Dunn 事后 ----
              if (identical(td_posthoc, "dunn")) {
                td_dunn <- tryCatch(
                  md_taxon_dunn(clr_taxa_sample, td_meta, td_group,
                                taxa = as.character(td_top$Taxon),
                                sample_col = smp, adjust = td_adjust),
                  error = function(e) {
                    td_notes <<- c(td_notes, paste0("md_taxon_dunn() 报错：", conditionMessage(e)))
                    NULL
                  })
                if (is.null(td_dunn) || !nrow(td_dunn)) {
                  td_dunn <- NULL
                  td_notes <- c(td_notes, "Dunn 事后检验没有结果（rstatix 缺失或分组退化）。")
                }
              } else {
                td_notes <- c(td_notes, "posthoc != 'dunn' —— 跳过 Dunn 事后检验。")
              }

              # ---- 7. 热图矩阵：行 = Top 类群，列 = 分组水平 ----
              tax <- as.character(td_top$Taxon)
              M_rel <- matrix(NA_real_, length(tax), length(td_grp_lv),
                              dimnames = list(tax, td_grp_lv))
              M_clr <- M_rel
              for (j in seq_along(td_grp_lv)) {
                idx <- which(as.character(td_meta[[td_group]]) == td_grp_lv[j])
                if (!length(idx)) next
                M_rel[, j] <- rowMeans(rel[tax, ids[idx], drop = FALSE], na.rm = TRUE)
                M_clr[, j] <- rowMeans(clr_taxa_sample[tax, ids[idx], drop = FALSE], na.rm = TRUE)
              }
              if (identical(td_scale, "row_z")) {
                M <- t(scale(t(M_rel)))
                M[!is.finite(M)] <- 0
                rsd <- apply(M, 1, stats::sd, na.rm = TRUE)
                if (all(!is.finite(rsd) | rsd < 1e-6)) {
                  td_notes <- c(td_notes, "row_z 退化（各类群组间模式几乎相同）—— 热图改用组平均相对丰度（%）。")
                  M <- M_rel
                  td_mat_label <- "Mean relative abundance (%)"
                } else {
                  td_mat_label <- "Z-score"
                }
              } else if (identical(td_scale, "clr")) {
                M <- M_clr
                td_mat_label <- "Mean CLR"
              } else {
                M <- M_rel
                td_mat_label <- "Mean relative abundance (%)"
              }
              rownames(M) <- local_td_ascii(tax)
              colnames(M) <- local_td_ascii(td_grp_lv)
              M[!is.finite(M)] <- 0
              td_mat <- M

              # ---- 8. 箱线图数据（CLR，每样本一行）----
              # Taxon 保留**原始类群名**（与 td_top 对得上）；只在图面文字上用 ASCII
              box_rows <- lapply(tax, function(tx) {
                data.frame(Taxon = tx,
                           Group = as.character(td_meta[[td_group]]),
                           CLR   = as.numeric(clr_taxa_sample[tx, ids]),
                           stringsAsFactors = FALSE)
              })
              td_box <- do.call(rbind, box_rows)
              td_box <- td_box[is.finite(td_box$CLR) & !is.na(td_box$Group), , drop = FALSE]
              td_ok <- TRUE
            } else {
              td_notes <- c(td_notes, "md_diff_taxa_kw() 没有返回任何类群（每组样本数 < min_n 或全为 NA）。")
            }
          } else {
            td_notes <- c(td_notes, "出现率过滤后不足 2 个单元 —— 无法做差异分析。")
          }
        }
      } else {
        td_notes <- c(td_notes, "丰度矩阵与元数据的样本交集 < 6")
      }
    } else {
      td_notes <- c(td_notes, "composition 口径的丰度矩阵不可用（05 核心对象缺失或过小）")
    }
    td_rank_used <- td_rank
    td_group_used <- td_group
    td_scale_used <- td_scale
  })
  if (length(td_notes)) notes <- c(notes, td_notes)

  if (!isTRUE(td_ok) || is.null(td_res) || !nrow(td_res) || is.null(td_mat)) {
    msg <- if (length(td_notes)) paste(td_notes, collapse = " | ") else "差异分析不可用"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", "differential taxa unavailable", notes = td_notes))
  }

  # ---------------------------------------------------------------------------
  # 出图
  # ---------------------------------------------------------------------------
  pt    <- md_pt(CFG)
  style <- CFG$taxon_diff$style %||% list()
  show_box <- isTRUE(CFG$taxon_diff$show_box %||% FALSE)

  pA <- tryCatch(
    local_td_heatmap(CFG, td_mat, engine, style,
                     title = local_td_ascii(sprintf("Differential taxa by %s", td_group_used)),
                     subtitle = if (is.null(td_fallback_txt)) NULL else local_td_ascii(td_fallback_txt),
                     legend_title = td_mat_label),
    error = function(e) {
      md_warn(CFG, "83: 热图绘制失败 —— ", conditionMessage(e)); NULL
    })
  if (is.null(pA)) pA <- md_empty_plot(CFG, "Heatmap unavailable")

  panels <- list(A = pA)
  if (show_box && !is.null(td_box) && nrow(td_box)) {
    pal <- md_pal(CFG, td_grp_lv, td_group_used)
    box_taxa <- unique(td_box$Taxon)
    bp <- list()
    for (tx in box_taxa) {
      dd <- td_box[td_box$Taxon == tx, , drop = FALSE]
      if (nrow(dd) < 3 || length(unique(dd$Group)) < 2) next
      lt <- NULL; kw <- NULL
      if (!is.null(td_dunn) && nrow(td_dunn)) {
        dt <- td_dunn[as.character(td_dunn$Taxon) == tx, , drop = FALSE]
        if (nrow(dt) && all(c("group1", "group2") %in% names(dt))) {
          pcol <- if ("p.adj" %in% names(dt)) "p.adj" else "p"
          lt <- tryCatch(md_sig_letters(dt, unique(dd$Group), p_col = pcol),
                         error = function(e) NULL)
        }
      }
      prow <- td_top[td_top$Taxon == tx, , drop = FALSE]
      if (nrow(prow)) kw <- md_kw_text(stats::setNames(prow$KW_padj[1], td_group_used))
      bp[[tx]] <- tryCatch(
        md_box_group(dd, y = "CLR", x = "Group", colors = pal, letters = lt,
                     ylab = local_td_ascii(tx), xlab = NULL, y_from_zero = FALSE,
                     x_angle = 45, legend_pos = "none", kw_text = kw, kw_annot = TRUE,
                     tag = NULL, pt = pt),
        error = function(e) {
          md_warn(CFG, "83: 箱线图失败 ", tx, " —— ", conditionMessage(e))
          md_empty_plot(CFG, local_td_ascii(tx))
        })
    }
    if (length(bp)) {
      pB <- if (length(bp) == 1) bp[[1]] else
        patchwork::wrap_plots(bp, ncol = min(3L, length(bp)))
      panels$B <- pB
    } else {
      notes <- c(notes, "show_box = TRUE 但没有可画的类群（每组样本不足）—— 只出热图面板。")
    }
  }

  fig <- tryCatch({
    if (length(panels) == 1) {
      # ⚠ 单面板**不走 md_compose()**：它会把 theme_md_comp() 套上去，
      #   而 ComplexHeatmap 经 ggplotify 包装后是"0-1 坐标系 + 空面板"，
      #   一旦被重新主题化就会长出 0.00-1.00 的 x/y 数值轴（34 号踩过同一个坑）。
      pA
    } else {
      md_compose(CFG, panels, design = "AAAA\nAAAA\nBBBB")
    }
  }, error = function(e) pA)
  # 高度按行数自适应：2-3 行时 150 mm 太空，行多时再拉高
  n_row <- nrow(td_mat)
  fig_h <- if (n_row <= 2) 80 else if (n_row <= 6) 110 else 150
  md_save_figure(CFG, fig, "83_taxon_diff", height = fig_h)
  md_maybe_single(CFG, pA, "83_taxon_diff_A")
  if (!is.null(panels$B)) md_maybe_single(CFG, panels$B, "83_taxon_diff_B")

  # ---------------------------------------------------------------------------
  # 结果表
  # ---------------------------------------------------------------------------
  top_out <- td_top
  if (identical(td_effect, "none") && "Cliff_delta" %in% names(top_out))
    top_out$Cliff_delta <- NULL
  if (!identical(td_effect, "none"))
    notes <- c(notes, "Cliff_delta 是**绝对效应量最大**的那一对的 delta（原脚本 809-824 行口径，不是加权平均）。")

  cliff_tab <- NULL
  if (!identical(td_effect, "none")) {
    keep_c <- intersect(c("Taxon", "Cliff_delta", "Enrich_group", "Mean_Rel_max",
                          grep("^Rel_", names(td_top), value = TRUE)), names(td_top))
    cliff_tab <- td_top[, keep_c, drop = FALSE]
    cliff_tab$Abs_delta <- abs(cliff_tab$Cliff_delta)
    cliff_tab$Effect_size <- ifelse(is.na(cliff_tab$Abs_delta), NA_character_,
                                    ifelse(cliff_tab$Abs_delta > 0.474, "large",
                                           ifelse(cliff_tab$Abs_delta > 0.33, "medium",
                                                  ifelse(cliff_tab$Abs_delta > 0.147, "small", "negligible"))))
    cliff_tab <- cliff_tab[order(-cliff_tab$Abs_delta), , drop = FALSE]
  }

  dunn_out <- td_dunn
  if (!is.null(dunn_out) && ".y." %in% names(dunn_out))
    names(dunn_out)[names(dunn_out) == ".y."] <- "Response"

  note_df <- data.frame(Note = c(
    "口径参考实现（原脚本）。",
    sprintf("层级 = %s；CLR = %s（伪计数 %s，NULL = 最小正值的 1/2）；调整方法 = %s。",
            td_rank_used, as.character(td_clr),
            if (is.null(td_pc)) "auto (min positive / 2)" else as.character(td_pc),
            td_adjust),
    sprintf("场景：按 design$factors 的第一个因子 '%s' 分组（原文同时做了 Zone 与 Vegetation 两套图/表）。",
            td_group_used),
    sprintf("热图：行 = Top %d 差异类群，列 = %s 的各水平，值 = %s（heatmap_scale = %s）。",
            nrow(td_mat), td_group_used, td_mat_label, td_scale_used),
    sprintf("原文热图是**样本级**（列 = 样本，column_split 按分组 + 流行率/显著性行注释）；本模块按契约改为**组水平均值矩阵**。"),
    sprintf("热图引擎 = %s；box 面板 = %s。", engine, as.character(show_box)),
    "显著字母（箱线图）来自 Dunn 的 p.adj 经 md_sig_letters() 转成的 CLD 字母。",
    notes))

  tp <- md_write_xlsx(CFG, "83_taxon_diff", list(
    Top_taxa    = top_out,
    Kruskal     = td_res,
    Dunn        = dunn_out,
    Cliff_delta = cliff_tab,
    Note        = note_df),
    meta_extra = list(rank = td_rank_used, group = td_group_used,
                      n_taxa_tested = nrow(td_res), n_taxa_plotted = nrow(td_mat),
                      heatmap_engine = engine))

  md_result(module, "OK", tables = basename(tp),
            figures = "figures/83_taxon_diff.pdf",
            notes = paste(notes, collapse = " | "))
}
