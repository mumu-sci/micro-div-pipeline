# =============================================================================
# R/21_phylo_diversity.R -- 系统发育 alpha 结构：Faith's PD + MPD/MNTD + NRI/NTI
# -----------------------------------------------------------------------------
# 本模块把「系统发育多样性（PD）」与「系统发育群落结构（NRI/NTI）」合并在一处：
# 两者都只依赖同一棵修剪后的树与同一个丰度矩阵，分开算会重复修剪、重复零模型。
#
# 方法（与 nifH 论文 13_Supp_NRI_NTI.R 逐条对齐）：
#   · 修剪：tree tip ∩ ASV 交集 -> 确保 rooted -> multi2di 处理多分叉
#   · Faith's PD：picante::pd(include.root = FALSE)（与 nifH 口径一致）
#   · MPD / MNTD：丰度加权（abundance.weighted = TRUE）
#   · 零模型：taxa.labels（等价于 nifH 的 sample(rownames) 逐次打乱 ASV 标签：
#             taxaShuffle 对行列同时施加同一置换，故与手写循环数值等价，
#             已用示范数据实测：|ΔMPD| < 5e-16、|ΔMNTD| < 3e-16）
#   · SES = (obs - null_mean) / null_sd；NRI = -SES_MPD，NTI = -SES_MNTD
#
# 计算走 picante::ses.mpd / ses.mntd（一次跑完全部 runs 与全部样本，比逐次
# 手写循环省内存、且自带 rank/p 值）；若该版本行为不符（缺函数 / 报错 / 结果
# 非法）则回退到本文件内的 local21_ 手写循环实现。
#
# 产物：cache/21_phylo.rds = data.frame(sample_id, PD, MPD, MNTD, NRI, NTI, <因子列>)
#       —— 见 MODULE_SPEC §8，71 等模块经 md_phylo_cache() 读取。
#
# 依赖：20（树准备）；参数全部来自 CFG$phylo / CFG$assembly / CFG$design / CFG$opts。
# =============================================================================

md_phylo_diversity <- function(CFG) {
  module <- "21"
  t_start <- Sys.time()
  if (!md_need(CFG, "ecology", module)) return(md_result(module, "SKIP", "缺少依赖包"))
  if (!md_precheck(CFG, module, need_response = FALSE, need_phylo = TRUE))
    return(md_result(module, "SKIP", "前置条件不满足"))

  md_load(c("ape", "picante"))

  runs <- as.integer(CFG$assembly$runs %||% 999L)
  if (!is.finite(runs) || runs < 1L) runs <- 999L
  prev_min <- CFG$phylo$prev_min %||% 2
  factors <- CFG$design$factors %||% character(0)
  smp_col <- CFG$data$sample_col %||% "sample_id"

  m    <- md_core_asv(CFG, "phylo")
  meta <- md_core_meta(CFG)
  tree <- md_core_tree(CFG)

  if (is.null(m) || ncol(m) < 3L || nrow(m) < 3L) {
    md_skip(CFG, module, "ASV 矩阵样本或类群过少")
    return(md_result(module, "SKIP", "样本或类群过少"))
  }
  if (is.null(tree)) {
    md_skip(CFG, module, "无可用系统发育树")
    return(md_result(module, "SKIP", "无系统发育树"))
  }

  # ---- 1. 树修剪 + 对齐 ----
  tr <- local21_phylo_prune(tree, rownames(m), CFG)
  if (is.null(tr) || length(tr$tip.label) < 3L) {
    md_skip(CFG, module, "ASV 与树 tip 匹配后类群不足 3 个")
    return(md_result(module, "SKIP", "树 tip 匹配不足"))
  }
  taxa <- tr$tip.label
  m    <- m[taxa, , drop = FALSE]
  m    <- m[, colSums(m) > 0, drop = FALSE]      # 全零样本剔除
  ids  <- colnames(m)
  if (length(ids) < 3L) {
    md_skip(CFG, module, "非空样本少于 3 个")
    return(md_result(module, "SKIP", "非空样本过少"))
  }
  md_log(CFG, "INFO", "系统发育修剪：", length(taxa), " 类群 x ", length(ids),
         " 样本 | 零模型 ", runs, " 次 | null.model = taxa.labels")

  samp <- t(m)                                    # samples x taxa
  coph <- stats::cophenetic(tr)
  coph <- as.matrix(coph)
  coph[!is.finite(coph)] <- NA_real_

  # ---- 2. Faith's PD ----
  pd_res <- tryCatch(picante::pd(samp, tr, include.root = FALSE), error = function(e) NULL)
  if (is.null(pd_res)) {
    md_warn(CFG, "picante::pd 失败，改用 local21_phylo_pd 手写实现")
    pd_res <- local21_phylo_pd(samp, tr)
  }
  PD <- as.numeric(pd_res$PD[match(ids, rownames(pd_res))])

  # ---- 3. MPD / MNTD 观察值 ----
  MPD  <- suppressWarnings(as.numeric(picante::mpd(samp, coph, abundance.weighted = TRUE)))
  MNTD <- suppressWarnings(as.numeric(picante::mntd(samp, coph, abundance.weighted = TRUE)))

  # ---- 4. 零模型 SES（并行；**整块缓存** —— 改样式重出图时不再重算）----
  ncores <- CFG$assembly$ncores %||% CFG$opts$ncores %||% 4L
  n_used <- max(1L, min(as.integer(ncores),
                        parallel::detectCores(logical = FALSE)))
  t_null <- Sys.time()
  ses_list <- md_cache(CFG, "21_ses", {
    md_parallel(CFG, n_used, seed = CFG$opts$seed + 211L)
    md_seed(CFG, 211)
    on.exit(md_parallel_off(), add = TRUE)
    local21_ses_mpd_mntd(CFG, samp, coph, runs = runs, ncores = n_used)
  })
  md_log(CFG, "INFO", "零模型就绪，用时 ",
         md_fmt_ms(0, as.numeric(difftime(Sys.time(), t_null, units = "secs"))),
         "（", runs, " 次 x 2 个指数，", ncores, " 核）")

  if (is.null(ses_list)) {
    md_warn(CFG, "零模型计算失败，NRI/NTI 记为 NA")
    ses_list <- list(NRI = rep(NA_real_, length(ids)), NTI = rep(NA_real_, length(ids)),
                     source = "failed", secs = NA_real_)
  }
  NRI <- ses_list$NRI
  NTI <- ses_list$NTI
  # 非有限 SES（null_sd = 0 / 单类群样本）一律记 NA，下游照 NA 处理
  NRI[!is.finite(NRI)] <- NA_real_
  NTI[!is.finite(NTI)] <- NA_real_

  md_log(CFG, "INFO", "PD 非NA ", sum(!is.na(PD)), "/", length(PD),
         " | MPD ", sum(!is.na(MPD)), " | MNTD ", sum(!is.na(MNTD)),
         " | NRI ", sum(!is.na(NRI)), " | NTI ", sum(!is.na(NTI)))

  # ---- 5. 组装跨模块契约表并落 cache/21_phylo.rds ----
  out <- data.frame(sample_id = ids, PD = PD, MPD = MPD, MNTD = MNTD,
                    NRI = NRI, NTI = NTI, stringsAsFactors = FALSE)
  mdl <- meta[match(ids, meta[[smp_col]]), , drop = FALSE]
  for (f in factors) if (f %in% names(mdl)) out[[f]] <- as.character(mdl[[f]])
  if (!is.null(CFG$design$site_var) && CFG$design$site_var %in% names(mdl))
    out[[CFG$design$site_var]] <- as.character(mdl[[CFG$design$site_var]])
  out$.n_taxa_pd <- as.integer(rowSums(samp > 0))
  attr(out, "params") <- list(runs = runs, null_model = "taxa.labels",
                              abundance_weighted = TRUE, include_root = FALSE,
                              n_taxa = length(taxa), n_sample = length(ids),
                              ncores = n_used, prev_min = prev_min,
                              ses_source = ses_list$source %||% "picante",
                              secs_null = as.numeric(difftime(Sys.time(), t_null, units = "secs")),
                              generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  saveRDS(out, md_path(CFG, "cache", "21_phylo.rds"))
  md_log(CFG, "OK", "写出 cache/21_phylo.rds（", nrow(out), " 样本 x ",
         ncol(out), " 列）")

  # ---- 6. 分组检验 ----
  test  <- CFG$alpha$test %||% "kruskal_dunn"
  padj  <- CFG$alpha$p_adjust %||% "BH"
  desc_list <- list(); test_list <- list(); ltr_list <- list()
  for (f in factors) {
    if (!f %in% names(out)) next
    ol <- CFG$design$ordered_levels[[f]]
    for (v in c("PD", "NRI", "NTI")) {
      g <- tryCatch(md_group_test(CFG, out, v, f, test = test, p_adjust = padj,
                                  min_n = 3, order_levels = ol),
                    error = function(e) NULL)
      if (is.null(g) || is.null(g$desc) || nrow(g$desc) == 0) next
      d <- g$desc; d$index <- v; d$factor <- f
      desc_list[[paste(v, f, sep = "_")]] <- d
      if (!is.null(g$letters) && nrow(g$letters) > 0)
        ltr_list[[paste(v, f, sep = "_")]] <- g$letters
      tt <- g$test
      if (!is.null(tt) && nrow(tt)) {
        tt$index <- v; tt$factor <- f
        test_list[[paste(v, f, sep = "_")]] <- tt
      }
    }
  }
  # 取某指标在某个因子下的 CLD 字母 / KW 的 P 值（供箱线图用）
  get_letters <- function(v, f) ltr_list[[paste(v, f, sep = "_")]]
  get_kw_p <- function(v, f) {
    tt <- test_list[[paste(v, f, sep = "_")]]
    if (is.null(tt) || !nrow(tt)) return(NULL)
    tt$p[1]
  }
  desc_df <- if (length(desc_list)) do.call(rbind, desc_list) else NULL
  test_df <- if (length(test_list)) do.call(rbind, test_list) else NULL
  if (!is.null(desc_df)) rownames(desc_df) <- NULL
  if (!is.null(test_df)) rownames(test_df) <- NULL

  # ---- 7. 图 ----
  # 面板构成按 references/06-figure-map.md：A = Faith's PD | B = NRI | C = NTI
  # 画法对照：01_Fig01.R §E 的 PD 箱线（走 11_shared.R §4.1 plot_dual_factor_panel()）
  #   与 13_Supp_NRI_NTI.R §6 的 NRI 箱线（走 11_shared.R §4.2 plot_nri_nti()）。
  # 字号走拼图轨道 MD_PT_COMP，与 md_compose() 末端套的 theme_md_comp() 同一条轨道。
  pt  <- md_pt(CFG, comp = TRUE)
  f1  <- factors[1] %||% NA_character_
  f2  <- if (length(factors) > 1) factors[2] else NA_character_
  lv_of <- function(v) {
    lv <- CFG$design$ordered_levels[[v]] %||% unique(as.character(out[[v]]))
    intersect(as.character(lv), unique(as.character(out[[v]])))
  }

  # 交互组（Zone_Vegetation）——NRI/NTI 与 PD 的 x 轴都用它（与 nifH 一致）。
  # ⚠ 只放在**画图副本** out_plot 里，绝不并进 out —— out 会整表写进
  #   results/21_phylo_diversity.xlsx 的 Phylo_indices sheet（不许改动结果表结构）。
  has2 <- !is.na(f1) && !is.na(f2)
  inter_vals <- if (has2) paste(out[[f1]], out[[f2]], sep = "_") else rep(NA_character_, nrow(out))
  inter_lv <- if (has2) {
    as.vector(t(outer(lv_of(f1), lv_of(f2), function(a, b) paste(a, b, sep = "_"))))
  } else character(0)
  out_plot <- out
  out_plot$.inter <- factor(inter_vals, levels = intersect(inter_lv, unique(inter_vals)))

  # 指标在某因子上的 KW P 值 -> 面板右上角文案
  kw_ps <- function(v) {
    ps <- list()
    for (f in c(f1, f2)) {
      if (is.na(f)) next
      p <- get_kw_p(v, f)
      if (!is.null(p)) ps[[f]] <- p
    }
    ps
  }
  # plot_dual_factor_panel()（PD 面板）用 "; " 连接（11_shared.R L1676）
  kw_pd <- function(v) md_kw_text(kw_ps(v))
  # plot_nri_nti()（NRI/NTI 面板）用 " ; " 连接（11_shared.R L1823，原脚本自身写法）
  kw_nri <- function(v) {
    ps <- kw_ps(v)
    if (!length(ps)) return(NULL)
    paste(vapply(names(ps), function(nm) paste0(nm, ": ", md_kw_fmt(ps[[nm]])),
                 character(1)), collapse = " ; ")
  }
  # 交互组 CLD 字母（原图：NRI/NTI 与 PD 都在 Zone×Vegetation 交互组上做 Dunn + multcompView）
  ltr_inter <- function(v) {
    g <- tryCatch(md_group_test(CFG, out_plot, v, ".inter", test = test, p_adjust = padj,
                                min_n = 3, order_levels = levels(out_plot$.inter)),
                  error = function(e) NULL)
    if (is.null(g)) NULL else g$letters
  }

  col_var <- if (has2) f2 else if (!is.na(f1)) f1 else NA_character_
  colours <- if (!is.na(col_var)) md_pal(CFG, lv_of(col_var), col_var) else NULL

  plots <- list()
  design <- NULL

  # A. Faith's PD（原 Fig.1E）——双因子箱线，按 Vegetation 上色
  if (!is.na(f1) && any(is.finite(out$PD))) {
    dd <- out
    dd[[f1]] <- as.character(dd[[f1]])
    if (!is.na(f2)) dd[[f2]] <- as.character(dd[[f2]])
    plots$A <- if (has2 && length(lv_of(f2)) >= 2) {
      tryCatch(
        md_box_two_factor(dd, y = "PD", f1 = f1, f2 = f2,
                          letters = ltr_inter("PD"), colors = colours,
                          f1_levels = lv_of(f1), f2_levels = lv_of(f2),
                          ylab = "PD", x_angle = 45, tag = "A",
                          pt = pt, legend_pos = "none",
                          inter_sep = "_", x_labels = "two_line",
                          kw_text = kw_pd("PD")),
        error = function(e) {
          md_warn(CFG, "PD 箱线图失败：", conditionMessage(e)); md_empty_plot(CFG, "PD")
        })
    } else {
      tryCatch(
        md_box_group(dd, y = "PD", x = f1, colors = colours, ylab = "PD",
                     xlab = NULL, letters = get_letters("PD", f1), tag = "A",
                     pt = pt, legend_pos = "none", kw_text = kw_pd("PD")),
        error = function(e) md_empty_plot(CFG, "PD"))
    }
  }

  # B / C. NRI / NTI（原 Fig.S4C / Fig.3A）——交互组箱线 + ±2 参考虚线 + KW 注释
  for (mv in c("NRI", "NTI")) {
    if (!mv %in% names(out) || all(is.na(out[[mv]]))) next
    tg <- c(NRI = "B", NTI = "C")[[mv]]
    plots[[tg]] <- tryCatch(
      local21_nri_nti_panel(out_plot, mv, ylab = mv, colour_var = col_var,
                            colours = colours, order_levels = levels(out_plot$.inter),
                            ref = c(-2, 2), letters = ltr_inter(mv),
                            kw_text = kw_nri(mv), tag = tg, x_one_line = TRUE,
                            pt = pt),
      error = function(e) {
        md_warn(CFG, mv, " 箱线图失败：", conditionMessage(e)); md_empty_plot(CFG, mv)
      })
  }

  if (length(plots) == 0) plots <- list(A = md_empty_plot(CFG, "No phylo index available"))
  ord <- intersect(c("A", "B", "C"), names(plots))
  plots <- plots[ord]
  design <- if (length(plots) >= 3) "AABB\nAACC" else if (length(plots) == 2) "AABB" else "AAAA"
  # 图例：照抄 13_Supp_NRI_NTI.R §6 L279-286（guides = "collect"，Vegetation 图例置底）
  p <- md_compose(CFG, plots, design = design, guides = "collect", legend = "bottom")
  md_save_figure(CFG, p, "21_phylo_diversity", height = 135)

  # ---- 8. 表 ----
  # ⚠ 这里**不再写"零模型耗时"** —— 那是每次运行都变的墙钟时间，写进结果表会让
  #   同一份数据两次运行的 xlsx 指纹对不上（可复现性校验会误报）。运行耗时统一进
  #   `_meta` sheet（`md_write_xlsx` 自动追加），结果表只放可复现的统计量。
  note <- data.frame(
    item = c("指标口径", "零模型", "SES 标准化", "PD", "非有限 SES"),
    value = c("MPD/MNTD 丰度加权（abundance.weighted = TRUE）",
              paste0("picante null.model = taxa.labels（等价于逐次打乱 ASV 标签），runs = ", runs),
              "SES = (obs - null_mean) / null_sd；NRI = -SES_MPD，NTI = -SES_MNTD",
              "picante::pd(include.root = FALSE)，与 nifH 论文口径一致",
              "null_sd = 0（样本内单类群 / 该样本 SES 退化）记为 NA，不参与检验"),
    stringsAsFactors = FALSE)
  sheets <- list(
    Phylo_indices = out,
    Group_desc    = desc_df,
    Group_test    = test_df,
    Note          = note)
  if (isTRUE(CFG$table$detail)) sheets$MPD_MNTD_obs <- out[, c("sample_id", "MPD", "MNTD", "NRI", "NTI")]
  tp <- md_write_xlsx(CFG, "21_phylo_diversity", sheets)

  md_log(CFG, "OK", "模块 21 用时 ",
         md_fmt_ms(0, as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  md_result(module, "OK", tables = basename(tp), figures = "figures/21_phylo_diversity.pdf",
            notes = paste("零模型", runs, "次 | SES 源:", ses_list$source %||% "picante"))
}

# -----------------------------------------------------------------------------
# local21_ 辅助函数（共享库缺的，写在本文件内，不动 R/shared/）
# -----------------------------------------------------------------------------

#' NRI / NTI 面板：逐行照抄 11_shared.R §4.2 的 plot_nri_nti()（L1770-1825）。
#'
#' 为什么在本文件内重写而不是用基元 md_nri_nti_panel()：
#' 该基元与原图有四处可见差异 ——
#'   (1) 多画一条左下角 "Kruskal-Wallis p = ..." 注释（原图没有）；
#'   (2) CLD 字母画在 max*1.08 / grey20 / vjust =-0.5（原图是 max*1.18 / black / vjust = 0）；
#'   (3) 参考线用 grey55 + MD_LW$ref（原图是红色虚线 alpha 0.5）；
#'   (4) jitter 缩到 0.8 倍（原图用整份 JIT_SIZE_STD）。
#' 按 FIGURE_ALIGN_SPEC 铁律 3「基元不够就照抄原脚本那几行」，在此照抄一份；
#' 共享库的问题只在报告里写，不改 R/shared/。
#'
#' @param df         含 metric / .inter / colour_var 的 data.frame
#' @param metric     指标列名（NRI / NTI）
#' @param colour_var 上色因子列名（原图 = Vegetation）
local21_nri_nti_panel <- function(df, metric, ylab, colour_var,
                                  order_levels = NULL, colours = NULL,
                                  ref = c(-2, 2), letters = NULL, kw_text = NULL,
                                  x_one_line = FALSE, tag = NULL, pt = MD_PT) {
  d <- df[!is.na(df[[metric]]), , drop = FALSE]
  d <- d[!is.na(d$.inter), , drop = FALSE]
  d$.inter <- droplevels(factor(as.character(d$.inter),
                                levels = intersect(as.character(order_levels),
                                                   unique(as.character(d$.inter)))))
  if (nrow(d) == 0 || nlevels(d$.inter) == 0) return(md_empty_plot(NULL, metric))

  # CLD 字母：原图画在每组 max(y) * 1.18、vjust = 0、加粗、黑色（L1786 / L1809）
  ld <- NULL
  if (!is.null(letters) && nrow(letters) > 0 && all(c("group", "letter") %in% names(letters))) {
    g_all <- as.character(letters$group)
    keep  <- intersect(g_all, levels(d$.inter))
    if (length(keep)) {
      max_y <- tapply(d[[metric]], as.character(d$.inter), max, na.rm = TRUE)
      ld <- data.frame(g = keep,
                       lab = as.character(letters$letter)[match(keep, g_all)],
                       y = as.numeric(max_y[keep]) * 1.18,
                       stringsAsFactors = FALSE)
    }
  }

  aes_main <- if (is.null(colour_var) || !colour_var %in% names(d))
    ggplot2::aes(x = .data$.inter, y = .data[[metric]]) else
    ggplot2::aes(x = .data$.inter, y = .data[[metric]], colour = .data[[colour_var]])

  lab_fun <- if (isTRUE(x_one_line)) function(v) v else function(v) gsub("_", "\n", v, fixed = TRUE)

  g <- ggplot2::ggplot(d, aes_main) +
    # 原图先 jitter 后 boxplot（点被箱体压住），顺序不要换
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = MD_PT_SIZE$jitter_width, seed = MD_JITTER_SEED), size = MD_PT_SIZE$jitter,
                         alpha = MD_PT_SIZE$alpha) +
    ggplot2::geom_boxplot(fill = NA, width = 0.55, outlier.shape = NA,
                          linewidth = MD_LW$box) +
    ggplot2::scale_x_discrete(limits = levels(d$.inter), drop = TRUE, labels = lab_fun) +
    ggplot2::labs(y = ylab, x = NULL, tag = tag) +
    theme_md(CFG = NULL) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = pt$tick))
  if (!is.null(colours) && !is.null(colour_var) && colour_var %in% names(d)) {
    g <- g + ggplot2::scale_colour_manual(values = colours, guide = "none")
  }
  if (!is.null(ld) && nrow(ld) > 0) {
    g <- g + ggplot2::geom_text(data = ld,
                                ggplot2::aes(x = .data$g, y = .data$y, label = .data$lab),
                                inherit.aes = FALSE, vjust = 0,
                                size = pt2gg(pt$ann), fontface = "bold", colour = "black")
  }
  if (!is.null(ref)) {
    g <- g + ggplot2::geom_hline(yintercept = ref, linetype = "dashed",
                                 colour = "red", alpha = 0.5)
  }
  # 顶部留白 15% 是必需的：字母在 max*1.18、KW 文案锚在面板顶部（原注释 L1817-1819）
  g <- g + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15)))
  if (!is.null(kw_text) && nzchar(kw_text)) {
    g <- g + ggplot2::annotate("text", x = Inf, y = Inf, label = kw_text,
                               hjust = 1.05, vjust = 1.5, size = pt2gg(pt$ann),
                               colour = "grey30")
  }
  g
}

#' 树修剪：保留与给定 ASV 的交集 -- 确保 rooted + 二叉化。
#' @return phylo 或 NULL（交集 < 3 个 tip）
local21_phylo_prune <- function(tree, asv_ids, CFG = NULL) {
  if (is.null(tree) || length(tree$tip.label) < 3L) return(NULL)
  tr <- tree
  tr$tip.label <- gsub('^["\']|["\']$', "", trimws(tr$tip.label))
  common <- intersect(asv_ids, tr$tip.label)
  if (length(common) < 3L) return(NULL)
  if (length(common) < length(tr$tip.label)) tr <- ape::keep.tip(tr, common)
  # 有非有限枝长时先修一下，否则 cophenetic 会出 NaN
  if (!is.null(tr$edge.length) && any(!is.finite(tr$edge.length))) {
    tr$edge.length[!is.finite(tr$edge.length)] <- 0
  }
  tr <- tryCatch(ape::multi2di(tr), error = function(e) tr)
  if (!ape::is.rooted(tr)) {
    rm <- CFG$phylo$root_method %||% "midpoint"
    tr <- tryCatch({
      if (identical(rm, "none")) tr
      else if (identical(rm, "outgroup")) local21_root_outgroup(tr)
      else if (requireNamespace("phangorn", quietly = TRUE)) phangorn::midpoint(tr)
      else local21_root_outgroup(tr)
    }, error = function(e) local21_root_outgroup(tr))
  }
  tr
}

#' 以最长枝 tip 为外群做根化（nifH 13_Supp_NRI_NTI.R 的口径）
local21_root_outgroup <- function(tr) {
  tryCatch({
    te <- which(tr$edge[, 2] <= length(tr$tip.label))
    if (!length(te)) return(tr)
    longest <- tr$edge[te, 2][which.max(tr$edge.length[te])]
    ape::root(tr, outgroup = longest, resolve.root = TRUE)
  }, error = function(e) {
    tryCatch(ape::root(tr, outgroup = tr$tip.label[1], resolve.root = TRUE),
             error = function(e2) tr)
  })
}

#' Faith's PD 手写回退实现（picante::pd 不可用时）。
#' 口径与 picante::pd(include.root = FALSE) 完全一致：存在类群诱导子树的枝长之和
#' （示范数据实测 |Δ| = 0）。
local21_phylo_pd <- function(samp, tree) {
  tip <- tree$tip.label
  pdv <- vapply(seq_len(nrow(samp)), function(i) {
    pres <- intersect(colnames(samp)[samp[i, ] > 0], tip)
    if (!length(pres)) return(0)
    sub <- tryCatch(ape::keep.tip(tree, pres), error = function(e) NULL)
    if (is.null(sub)) return(NA_real_)
    sum(sub$edge.length)
  }, numeric(1))
  data.frame(PD = pdv, SR = rowSums(samp > 0), row.names = rownames(samp))
}

#' 从 picante ses.* 结果里取 SES 列。列名跨版本可能不同，按候选名依次尝试。
#' 坑：ses.mpd 的 SES 列名是 `mpd.obs.z`，ses.mntd 的是 `mntd.obs.z`
#' （并非同名），用错会得到 numeric(0)，match 之后整列变 NA。
local21_z_col <- function(x, prefer) {
  cand <- unique(c(prefer, "obs.z", grep("\\.obs\\.z$", names(x), value = TRUE)))
  for (nm in cand) if (!is.null(nm) && !is.na(nm) && nm %in% names(x)) return(x[[nm]])
  NULL
}

#' SES 零模型：优先 picante::ses.mpd / ses.mntd，失败或结果非法时回退手写循环。
#' @return list(NRI, NTI, source, secs)
local21_ses_mpd_mntd <- function(CFG, samp, coph, runs, ncores = 1L) {
  t0 <- Sys.time()
  res <- tryCatch({
    a <- picante::ses.mpd(samp, coph, null.model = "taxa.labels",
                          abundance.weighted = TRUE, runs = runs)
    b <- picante::ses.mntd(samp, coph, null.model = "taxa.labels",
                           abundance.weighted = TRUE, runs = runs)
    za <- local21_z_col(a, "mpd.obs.z")
    zb <- local21_z_col(b, "mntd.obs.z")
    if (is.null(za) || is.null(zb))
      stop("picante 返回结构不含预期的 SES 列：", paste(names(b), collapse = ", "))
    list(NRI = -as.numeric(za)[match(rownames(samp), rownames(a))],
         NTI = -as.numeric(zb)[match(rownames(samp), rownames(b))],
         source = "picante::ses.mpd/ses.mntd")
  }, error = function(e) {
    md_warn(CFG, "picante ses.* 失败（", conditionMessage(e), "），回退手写零模型")
    NULL
  })
  if (!is.null(res) && length(res$NRI) == nrow(samp) &&
      (any(is.finite(res$NRI)) || any(is.finite(res$NTI)))) {
    res$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    return(res)
  }
  # ---- 手写回退：与 nifH 13_Supp_NRI_NTI.R 完全同构 ----
  md_log(CFG, "INFO", "零模型回退：手写循环（taxa.labels，", runs, " 次）")
  obs_mpd  <- as.numeric(picante::mpd(samp, coph, abundance.weighted = TRUE))
  obs_mntd <- as.numeric(picante::mntd(samp, coph, abundance.weighted = TRUE))
  seeds <- local21_perm_seeds(CFG$opts$seed, 211, runs)
  use_future <- requireNamespace("future.apply", quietly = TRUE) && ncores > 1L
  one <- function(r) {
    set.seed(seeds[r])
    ds <- picante::taxaShuffle(coph)
    c(picante::mpd(samp, ds, abundance.weighted = TRUE),
      picante::mntd(samp, ds, abundance.weighted = TRUE))
  }
  mat <- if (use_future) {
    future.apply::future_lapply(seq_len(runs), one, future.seed = NULL)
  } else {
    lapply(seq_len(runs), one)
  }
  mat <- do.call(rbind, mat)
  n <- nrow(samp)
  m_mpd  <- mat[, seq_len(n), drop = FALSE]
  m_mntd <- mat[, n + seq_len(n), drop = FALSE]
  ses_mpd  <- (obs_mpd  - colMeans(m_mpd,  na.rm = TRUE)) / apply(m_mpd,  2, stats::sd, na.rm = TRUE)
  ses_mntd <- (obs_mntd - colMeans(m_mntd, na.rm = TRUE)) / apply(m_mntd, 2, stats::sd, na.rm = TRUE)
  list(NRI = -ses_mpd, NTI = -ses_mntd, source = "manual (taxa.labels shuffle)",
       secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#' 由全包唯一种子派生一串子种子（同一个 offset 永远得到同一串）
local21_perm_seeds <- function(seed, offset, n) {
  as.integer(seed) + as.integer(offset) + seq_len(as.integer(n))
}
