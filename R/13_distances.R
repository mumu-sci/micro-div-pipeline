# =============================================================================
# R/13_distances.R —— 模块 13：β 距离矩阵（bray / jaccard / unifrac / wunifrac）
# -----------------------------------------------------------------------------
# 导出：md_distances(CFG)
#
# 口径：CFG$normalize$by_module$distances（默认 rarefy，抽平到 CFG$normalize$depth）
# 度量：CFG$distances$metrics（bray | jaccard | unifrac | wunifrac）
#       jaccard 是否二值化 = CFG$distances$binary_jaccard（默认 TRUE，与 nifH 一致）
#       unifrac / wunifrac 需要系统发育树；无树自动跳过并打印原因（不报错）
#
# **跨模块契约**：写出 cache/13_dist.rds，
#   structure = list(dist = <具名 dist 列表>, meta = data.frame, params = list(...))
#   距离矩阵只存 rds，绝不并入 data.frame（供 14 / 15 / 16 / 33 / 43 / 71 读取）。
#
# 出图：figures/13_distances.pdf —— A = 样本间距离热图（按分组排序）
#                                   B = 各度量成对距离分布箱线
# 出表：results/13_distances.xlsx
# =============================================================================

# ---------------------------------------------------------------------------
# 模块内小工具（local_ 前缀）
# ---------------------------------------------------------------------------
local_factor_levels <- function(CFG, meta, var) {
  ol <- CFG$design$ordered_levels[[var]]
  obs <- unique(as.character(meta[[var]]))
  obs <- obs[!is.na(obs)]
  if (!is.null(ol)) return(intersect(as.character(ol), obs))
  obs
}

local_order_by_group <- function(CFG, meta, f) {
  if (is.null(f) || !f %in% names(meta)) return(seq_len(nrow(meta)))
  x <- factor(as.character(meta[[f]]), levels = local_factor_levels(CFG, meta, f))
  order(x, seq_along(x))
}

#' 长表 -> 上三角 dist（供后续模块使用，不依赖共享库的实现细节）
local_dist_from_long <- function(long, ids) {
  m <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  i <- match(long$sample1, ids); j <- match(long$sample2, ids)
  ok <- !is.na(i) & !is.na(j)
  m[cbind(i[ok], j[ok])] <- long$value[ok]
  m[cbind(j[ok], i[ok])] <- long$value[ok]
  stats::as.dist(m)
}

#' 把距离矩阵的样本子集与分组顺序对齐（保持对称、对角线为 0）
local_dist_reorder <- function(d, ids) {
  m <- as.matrix(d)
  ids <- intersect(ids, rownames(m))
  if (length(ids) < 2) return(NULL)
  stats::as.dist(m[ids, ids, drop = FALSE])
}

# ---------------------------------------------------------------------------
# 模块主体
# ---------------------------------------------------------------------------
md_distances <- function(CFG) {
  module <- "13"
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少 ecology 包组"))

  m <- md_core_asv(CFG, "distances")
  if (is.null(m) || nrow(m) < 2 || ncol(m) < 2)
    return(md_result(module, "SKIP", "距离口径矩阵为空或样本不足"))
  # vegdist 针对样方 × 物种；本包矩阵为 ASV × 样本
  x <- t(m)
  meta <- md_core_meta(CFG)
  smp <- CFG$data$sample_col
  meta <- meta[match(rownames(x), meta[[smp]]), , drop = FALSE]
  if (anyNA(meta[[smp]])) meta[[smp]] <- rownames(x)[is.na(meta[[smp]])]
  rownames(meta) <- meta[[smp]]
  ids <- meta[[smp]]

  metrics <- tolower(CFG$distances$metrics %||% c("bray", "jaccard", "wunifrac"))
  metrics <- unique(metrics[nzchar(metrics)])
  binary_jac <- isTRUE(CFG$distances$binary_jaccard %||% TRUE)
  tree <- tryCatch(md_core_tree(CFG), error = function(e) NULL)
  has_tree <- !is.null(tree)

  dist_list <- list()
  skip_msg <- character(0)

  md_seed(CFG, 130)

  # ---- 1. 分类学距离 ----
  if ("bray" %in% metrics) {
    d <- tryCatch(vegan::vegdist(x, method = "bray"), error = function(e) NULL)
    if (is.null(d)) skip_msg <- c(skip_msg, "bray: vegdist 失败")
    else dist_list$bray <- d
  }
  if ("jaccard" %in% metrics) {
    d <- tryCatch(vegan::vegdist(x, method = "jaccard", binary = binary_jac),
                  error = function(e) NULL)
    if (is.null(d)) skip_msg <- c(skip_msg, "jaccard: vegdist 失败")
    else dist_list$jaccard <- d
  }
  # 兼容别名
  if ("bray-curtis" %in% metrics && is.null(dist_list$bray)) {
    d <- tryCatch(vegan::vegdist(x, method = "bray"), error = function(e) NULL)
    if (!is.null(d)) dist_list$bray <- d
  }

  # ---- 2. 系统发育距离 ----
  want_phy <- intersect(c("unifrac", "wunifrac"), metrics)
  if (length(want_phy)) {
    if (!has_tree) {
      skip_msg <- c(skip_msg, paste0(paste(want_phy, collapse = "/"),
                                     ": 无系统发育树（data$tree 为 NULL 或 phylo$enabled=FALSE）"))
      md_log(CFG, "INFO", "跳过 ", paste(want_phy, collapse = "/"),
             " —— 无系统发育树（data$tree 为空或 phylo$enabled=FALSE）")
    } else if (!md_need(CFG, "unifrac", module)) {
      skip_msg <- c(skip_msg, paste0(paste(want_phy, collapse = "/"), ": 缺少 GUniFrac/phangorn"))
    } else {
      tips <- intersect(colnames(x), tree$tip.label)
      if (length(tips) < 3) {
        skip_msg <- c(skip_msg, paste0(paste(want_phy, collapse = "/"),
                                       ": 树与 ASV 匹配过少（", length(tips), " 个 tip）"))
      } else {
        sub <- tryCatch(ape::drop.tip(tree, setdiff(tree$tip.label, tips)),
                        error = function(e) NULL)
        if (is.null(sub)) {
          skip_msg <- c(skip_msg, paste0(paste(want_phy, collapse = "/"), ": 树剪枝失败"))
        } else {
          if (!is.null(sub$edge.length) && any(sub$edge.length < 0)) sub$edge.length <- pmax(sub$edge.length, 0)
          if (is.null(sub$edge.length)) sub$edge.length <- rep(1, nrow(sub$edge))
          if (!ape::is.rooted(sub)) {
            rr <- tryCatch(if (requireNamespace("phangorn", quietly = TRUE))
              phangorn::midpoint(sub) else ape::root(sub, outgroup = sub$tip.label[1], resolve.root = TRUE),
              error = function(e) NULL)
            if (!is.null(rr)) sub <- rr
            else skip_msg <- c(skip_msg, paste0(paste(want_phy, collapse = "/"), ": 无法生根"))
          }
          md_seed(CFG, 131)
          guf <- tryCatch(GUniFrac::GUniFrac(x[, tips, drop = FALSE], sub,
                                             alpha = c(0, 0.5, 1)),
                          error = function(e) { md_warn(CFG, "GUniFrac 失败: ", conditionMessage(e)); NULL })
          if (is.null(guf)) {
            skip_msg <- c(skip_msg, paste0(paste(want_phy, collapse = "/"), ": GUniFrac 计算失败"))
          } else {
            uf <- guf$unifracs
            if ("unifrac" %in% metrics && "d_0" %in% dimnames(uf)[[3]])
              dist_list$unifrac <- stats::as.dist(uf[, , "d_0"])
            if ("wunifrac" %in% metrics && "d_1" %in% dimnames(uf)[[3]])
              dist_list$wunifrac <- stats::as.dist(uf[, , "d_1"])
          }
        }
      }
    }
  }

  if (length(dist_list) == 0) {
    msg <- paste0("没有任何距离度量可用：", paste(skip_msg, collapse = " | "))
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }

  # ---- 3. 统一样本顺序 + 抽样自检 ----
  for (k in names(dist_list)) {
    d <- dist_list[[k]]
    dm <- as.matrix(d)
    if (!identical(rownames(dm), ids)) {
      d <- local_dist_reorder(d, ids)
      if (is.null(d)) { dist_list[[k]] <- NULL; next }
    }
    v <- as.numeric(d)
    if (any(!is.finite(v))) {
      md_warn(CFG, "距离矩阵 ", k, " 存在非有限值，已置为 NA")
      d[!is.finite(v)] <- NA
    }
    if (any(v < 0, na.rm = TRUE)) {
      md_warn(CFG, "距离矩阵 ", k, " 存在负值（可能是 Cailliez/UniFrac 数值噪声），保留原值")
    }
    dist_list[[k]] <- d
  }
  if (length(dist_list) == 0) return(md_result(module, "SKIP", "所有距离矩阵重排后失效"))

  # ---- 4. 写 cache/13_dist.rds（**跨模块契约**）----
  params <- list(
    metrics = names(dist_list),
    metrics_requested = metrics,
    skipped = skip_msg,
    binary_jaccard = binary_jac,
    normalize = CFG$normalize$by_module$distances %||% CFG$normalize$method,
    depth = CFG$normalize$depth %||% (tryCatch(md_core(CFG)$depth, error = function(e) NA)),
    tree_used = has_tree && any(c("unifrac", "wunifrac") %in% names(dist_list)),
    n_tree_tips = if (has_tree) length(tree$tip.label) else 0L,
    seed = CFG$opts$seed,
    n_sample = length(ids),
    sample_col = smp,
    md_version = MD_VERSION)
  cache_obj <- list(dist = dist_list, meta = meta, params = params)
  saveRDS(cache_obj, md_path(CFG, "cache", "13_dist.rds"))
  md_log(CFG, "OK", "写出 cache/13_dist.rds（", paste(names(dist_list), collapse = ", "), "）")

  # ---- 5. 汇总表 ----
  tab_sum <- do.call(rbind, lapply(names(dist_list), function(k) {
    v <- as.numeric(dist_list[[k]]); v <- v[is.finite(v)]
    if (length(v) == 0) return(NULL)
    q <- unname(stats::quantile(v, c(0, .25, .5, .75, 1)))
    data.frame(metric = k, n_pairs = length(v),
               min = q[1], q1 = q[2], median = q[3], q3 = q[4], max = q[5],
               mean = mean(v), sd = stats::sd(v),
               pct_zero = 100 * mean(v == 0), pct_gt_0.9 = 100 * mean(v > 0.9),
               stringsAsFactors = FALSE)
  }))
  if (!is.null(tab_sum)) tab_sum[setdiff(names(tab_sum)[vapply(tab_sum, is.numeric, TRUE)], "n_pairs")] <-
    lapply(tab_sum[setdiff(names(tab_sum)[vapply(tab_sum, is.numeric, TRUE)], "n_pairs")],
           function(z) round(z, CFG$table$digits %||% 4))

  # ---- 6. 出图 ----
  f1 <- intersect(CFG$design$factors, names(meta))[1]
  ord <- local_order_by_group(CFG, meta, f1)
  key <- names(dist_list)[1]
  dm <- as.matrix(dist_list[[key]])[ord, ord, drop = FALSE]
  hm <- data.frame(
    sample1 = factor(rep(rownames(dm), times = ncol(dm)), levels = rownames(dm)),
    sample2 = factor(rep(colnames(dm), each = nrow(dm)), levels = colnames(dm)),
    value = as.numeric(dm), stringsAsFactors = FALSE)
  grp <- if (!is.null(f1)) as.character(meta[[f1]][ord]) else rep("All", nrow(dm))
  pA <- md_heatmap(hm, x = "sample2", y = "sample1", fill = "value", values = md_pal_cont(CFG, "sequential"),
                   xlab = NULL, ylab = NULL, tag = "A", label = FALSE,
                   mid = NULL, pt = MD_PT, x_angle = 90) +
    ggplot2::labs(fill = paste0("Distance (", key, ")")) +
    ggplot2::annotate("text", x = 1, y = nrow(dm) + 0.6, label = "samples ordered by group",
                      hjust = 0, vjust = 0, size = pt2gg(MD_PT$ann * 0.85), colour = "grey35") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "bottom")
  if (!is.null(f1) && length(unique(grp)) <= 8) {
    bounds <- cumsum(rle(grp)$lengths)
    if (length(bounds) > 1)
      pA <- pA + ggplot2::geom_hline(yintercept = bounds[-length(bounds)] + 0.5,
                                     linewidth = 0.25, colour = "grey25") +
        ggplot2::geom_vline(xintercept = bounds[-length(bounds)] + 0.5,
                            linewidth = 0.25, colour = "grey25")
  }

  ok_metrics <- vapply(dist_list, function(d) any(is.finite(as.numeric(d))), logical(1))
  pB <- if (any(ok_metrics)) {
    long <- do.call(rbind, lapply(names(dist_list)[ok_metrics], function(k) {
      v <- as.numeric(dist_list[[k]]); v <- v[is.finite(v)]
      if (length(v) == 0) return(NULL)
      data.frame(metric = k, value = v, stringsAsFactors = FALSE)
    }))
    md_seed(CFG, 132)
    med <- stats::aggregate(value ~ metric, long, stats::median)
    ggplot2::ggplot(long, ggplot2::aes(x = .data$metric, y = .data$value)) +
      ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA, fill = "white",
                            linewidth = MD_LW$box, colour = "grey25") +
      ggplot2::geom_jitter(position = ggplot2::position_jitter(width = 0.12, seed = MD_JITTER_SEED), size = MD_PT_SIZE$jitter * 0.5,
                           alpha = 0.18, colour = "#3B6E8F") +
      ggplot2::geom_text(data = med, ggplot2::aes(x = .data$metric, y = .data$value,
                                                  label = md_fmt_num(.data$value, 2)),
                         vjust = -1.2, size = pt2gg(MD_PT$ann * 0.85), colour = "grey30") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
      ggplot2::labs(x = NULL, y = "Pairwise distance", tag = "B") +
      theme_md(CFG = NULL)
  } else md_empty_plot(CFG, "No finite distance value")

  pC <- if (length(skip_msg)) {
    txt <- paste(strwrap(paste(skip_msg, collapse = "\n"), width = 70), collapse = "\n")
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = paste0("Skipped metrics:\n", txt),
                        size = pt2gg(MD_PT$ann * 0.9), colour = "grey35", hjust = 0.5, vjust = 0.5) +
      ggplot2::theme_void()
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
                        label = paste0("Metrics computed: ", paste(names(dist_list), collapse = ", "),
                                       "\n(", params$normalize, " normalization)"),
                        size = pt2gg(MD_PT$ann * 0.9), colour = "grey35") +
      ggplot2::theme_void()
  }

  fig <- md_compose(CFG, list(A = pA, B = pB, C = pC), design = "AABB\nAACC\nAACC")
  md_save_figure(CFG, fig, "13_distances", height = 165)
  md_maybe_single(CFG, pA, "13_distances_A")
  md_maybe_single(CFG, pB, "13_distances_B")

  # ---- 7. 表 ----
  tab_params <- data.frame(
    item = c("metrics_requested", "metrics_computed", "metrics_skipped", "binary_jaccard",
             "normalization (distances)", "rarefaction depth", "tree used", "n_tree_tips",
             "n samples", "seed", "cache file"),
    value = c(paste(metrics, collapse = ", "), paste(names(dist_list), collapse = ", "),
              if (length(skip_msg)) paste(skip_msg, collapse = " | ") else "none",
              as.character(binary_jac), as.character(params$normalize),
              as.character(params$depth), as.character(params$tree_used),
              as.character(params$n_tree_tips), as.character(length(ids)),
              as.character(CFG$opts$seed), "cache/13_dist.rds"),
    stringsAsFactors = FALSE)
  tab_long <- if (isTRUE(CFG$table$detail)) {
    do.call(rbind, lapply(names(dist_list), function(k) {
      d <- as.matrix(dist_list[[k]])
      ut <- which(upper.tri(d), arr.ind = TRUE)
      data.frame(metric = k, sample1 = rownames(d)[ut[, 1]], sample2 = rownames(d)[ut[, 2]],
                 distance = d[ut], stringsAsFactors = FALSE)
    }))
  } else NULL

  tp <- md_write_xlsx(CFG, "13_distances", list(
    Distance_summary = tab_sum,
    Parameters = tab_params,
    Pairwise_detail = tab_long))

  md_result(module, "OK", tables = basename(tp), figures = "figures/13_distances.pdf",
            notes = paste0("metrics=", paste(names(dist_list), collapse = "+"),
                           if (length(skip_msg)) paste0("; skipped: ", paste(skip_msg, collapse = " | ")) else ""))
}
