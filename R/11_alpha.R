# =============================================================================
# R/11_alpha.R —— 模块 11：α 多样性（指数 + 组间检验 + CLD 字母 + 稀释曲线）
# -----------------------------------------------------------------------------
# 导出：md_alpha(CFG)
#
# 口径：
#   指数计算 = CFG$normalize$by_module$alpha（默认 rarefy，抽平到 CFG$normalize$depth）
#   稀释曲线 = CFG$normalize$by_module$rarefaction（默认 none，即原始 reads）
#
# 指数（CFG$alpha$indices）：Observed / Chao1 / ACE / Shannon / Simpson / Pielou / PD
#   Pielou = Shannon / log(Observed)；PD = 系统发育树 Faith's PD（无树则跳过并在表里说明）
# 检验（CFG$alpha$test）：anova_tukey | kruskal_dunn | lmm | none
#   显著性字母 = md_group_test() -> md_attach_letters()（CLD，multcompView）
# 稀释曲线（CFG$alpha$rarefaction$curve）：vegan | inext | both
#
# 出图：figures/11_alpha.pdf（A–G 合成 patchwork）
# 出表：results/11_alpha.xlsx
# =============================================================================

# ---------------------------------------------------------------------------
# 模块内小工具（local_ 前缀）
# ---------------------------------------------------------------------------
local_glue <- function(parts, sep = "_") {
  if (length(parts) == 0) return(character(0))
  parts <- lapply(parts, function(x) as.character(x))
  n <- max(vapply(parts, length, integer(1)))
  parts <- lapply(parts, function(x) rep(x, length.out = n))
  do.call(paste, c(parts, list(sep = sep)))
}

local_factor_levels <- function(CFG, meta, var) {
  ol <- CFG$design$ordered_levels[[var]]
  obs <- unique(as.character(meta[[var]]))
  obs <- obs[!is.na(obs)]
  if (!is.null(ol)) return(intersect(as.character(ol), obs))
  obs
}

local_group_var <- function(CFG, meta, factors = NULL, name = NULL) {
  factors <- factors %||% CFG$design$factors
  factors <- intersect(factors, names(meta))
  if (length(factors) == 0) return(NULL)
  lv <- lapply(factors, function(f) local_factor_levels(CFG, meta, f))
  combos <- expand.grid(lv, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # ⚠ 这里**不能**对 combos 的列做 rev()：levels（水平名）与 vals（取值）必须用
  #   **同一个因子顺序**拼接，否则 "Bare_Coastal" 与 "Coastal_Bare" 对不上，
  #   factor(vals, levels = levels) 会全部变成 NA —— 表现为「交互检验在表里缺失」
  #   且「箱线图上没有显著性字母」（字母分组名与图上级别不匹配被丢弃）。
  labs <- local_glue(combos, sep = "_")
  vals <- local_glue(lapply(factors, function(f) meta[[f]]), sep = "_")
  names(vals) <- as.character(meta[[CFG$data$sample_col]] %||% seq_along(vals))
  list(vals = vals, levels = labs, factors = factors, name = name %||% paste(factors, collapse = "_"))
}

#' 指标名 -> 纵轴题；沿用同一套样式（L1602-1609）：
#'   Observed 在图上写作 Richness，Faith's PD 写作 PD（nifH 原图口径，不要自己改写）
local_index_label <- function(ix) {
  m <- c(Pielou = "Pielou", Observed = "Richness", Shannon = "Shannon",
         PD = "PD", Faith_PD = "PD", Chao1 = "Chao1", ACE = "ACE", Simpson = "Simpson")
  v <- unname(m[ix])
  if (length(v) == 0 || is.na(v)) ix else v
}

#' 取某个指数在指定因子上的 KW 注释文案（供箱线图右上角用，纯 ASCII）
local_kw_text <- function(tab_test, facs, ix) {
  if (is.null(tab_test) || !nrow(tab_test)) return(NULL)
  sub <- tab_test[tab_test$index == ix & tab_test$factor %in% facs, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  md_kw_text(stats::setNames(as.list(sub$p_adj), sub$factor))
}

local_render <- function(x, th) {  if (inherits(x, "patchwork_list") || inherits(x, "patchwork")) {
    nm <- names(x)
    for (i in seq_along(x)) x[[i]] <- local_render(x[[i]], th)
    return(x)
  }
  x + th
}

#' 面板列表 -> 单张图（面板数 <= 3 用一行，否则两行）
local_grid <- function(plots, theme_fn = NULL) {
  plots <- Filter(function(p) inherits(p, "gg") || inherits(p, "patchwork"), plots)
  if (length(plots) == 0) return(NULL)
  p <- patchwork::wrap_plots(plots, ncol = if (length(plots) <= 3) length(plots) else 2)
  if (!is.null(theme_fn)) p <- local_render(p, theme_fn)
  p
}

#' Faith's PD；失败返回 NULL（调用方优雅降级）
#' 注：picante::pd(include.root=TRUE) 要求有根树，无根树时先 midpoint rooting
local_pd <- function(m, tree) {
  if (is.null(tree) || !requireNamespace("picante", quietly = TRUE)) return(NULL)
  tips <- intersect(rownames(m), tree$tip.label)
  if (length(tips) < 2) return(NULL)
  sub <- tryCatch(ape::drop.tip(tree, setdiff(tree$tip.label, tips)), error = function(e) NULL)
  if (is.null(sub) || is.null(sub$edge.length)) return(NULL)
  if (!ape::is.rooted(sub)) {
    rr <- tryCatch(if (requireNamespace("phangorn", quietly = TRUE))
      phangorn::midpoint(sub)
      else ape::root(sub, outgroup = sub$tip.label[1], resolve.root = TRUE),
      error = function(e) NULL)
    if (is.null(rr) || !ape::is.rooted(rr)) return(NULL)
    sub <- rr
  }
  comm <- t(m[tips, , drop = FALSE])
  out <- tryCatch(picante::pd(comm, sub, include.root = TRUE), error = function(e) NULL)
  if (is.null(out) || !"PD" %in% names(out)) return(NULL)
  stats::setNames(as.numeric(out$PD), rownames(out))
}

#' 稳健的行合并：自动跳过 NULL / 零行元素，并按并集补列
#' （base::rbind 无法合并 0x13 与 3x13，全 NA 指数如无树时的 PD 会踩这个坑）
local_rbind <- function(x) {
  x <- Filter(function(z) !is.null(z) && (is.data.frame(z) || is.matrix(z)) && NROW(z) > 0, x)
  if (length(x) == 0) return(NULL)
  if (length(x) == 1) return(as.data.frame(x[[1]], stringsAsFactors = FALSE))
  cols <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(z) {
    z <- as.data.frame(z, stringsAsFactors = FALSE)
    for (cc in setdiff(cols, names(z))) z[[cc]] <- NA
    z[, cols, drop = FALSE]
  })
  out <- x[[1]]
  for (i in seq_along(x)[-1]) out <- rbind(out, x[[i]])
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# 模块主体
# ---------------------------------------------------------------------------
md_alpha <- function(CFG) {
  module <- "11"
  if (!md_need(CFG, c("core", "ecology"), module))
    return(md_result(module, "SKIP", "缺少依赖包（core/ecology）"))
  m <- md_core_asv(CFG, "alpha")
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0)
    return(md_result(module, "SKIP", "α 口径矩阵为空"))
  meta <- md_core_meta(CFG)
  smp <- CFG$data$sample_col
  meta <- meta[match(colnames(m), meta[[smp]]), , drop = FALSE]

  want <- CFG$alpha$indices %||% c("Observed", "Shannon", "Pielou", "PD")
  want <- unique(as.character(want))
  if (length(want) == 0) want <- c("Observed", "Shannon")

  notes <- character(0)
  idx <- data.frame(sample_id = colnames(m), stringsAsFactors = FALSE)
  .mk <- local_marker(CFG)

  # ---- 1. 指数计算 ----
  pa <- m > 0
  n_present <- colSums(pa); tot <- colSums(m)
  if ("Observed" %in% want) idx$Observed <- as.numeric(n_present)
  if ("Shannon" %in% want) {
    idx$Shannon <- apply(m, 2, function(x) {
      p <- x / sum(x); p <- p[p > 0]; -sum(p * log(p))
    })
  }
  if ("Simpson" %in% want) {
    idx$Simpson <- apply(m, 2, function(x) {
      p <- x / sum(x); 1 - sum(p^2)
    })
  }
  if ("Pielou" %in% want) {
    sh <- if ("Shannon" %in% names(idx)) idx$Shannon else
      apply(m, 2, function(x) { p <- x / sum(x); p <- p[p > 0]; -sum(p * log(p)) })
    obs <- if ("Observed" %in% names(idx)) idx$Observed else as.numeric(n_present)
    idx$Pielou <- ifelse(obs > 1, sh / log(obs), NA_real_)
  }
  if (any(c("Chao1", "ACE") %in% want)) {
    est <- tryCatch({
      z <- suppressWarnings(phyloseq::estimate_richness(
        phyloseq::otu_table(m, taxa_are_rows = TRUE),
        measures = intersect(c("Chao1", "ACE"), want)))
      if (is.null(z) || nrow(z) == 0 || ncol(z) == 0) NULL else z
    }, error = function(e) NULL)
    # 样本名规范化：phyloseq 会把 '-' 转成 '.'（nifH 原代码同样处理）
    key <- function(x) gsub("-", ".", as.character(x))
    if (is.null(est)) {
      md_warn(CFG, "phyloseq::estimate_richness 失败，改用 vegan::estimateR")
      est <- tryCatch({
        z <- vegan::estimateR(t(m))
        d <- as.data.frame(t(z), stringsAsFactors = FALSE)
        if ("S.chao1" %in% names(d)) names(d)[names(d) == "S.chao1"] <- "Chao1"
        d
      }, error = function(e) NULL)
    }
    for (nm in intersect(c("Chao1", "ACE"), want)) {
      if (!is.null(est) && nm %in% names(est)) {
        v <- suppressWarnings(as.numeric(est[[nm]][match(key(colnames(m)), key(rownames(est)))]))
        v[is.nan(v) | !is.finite(v)] <- NA_real_
        idx[[nm]] <- v
        if (all(is.na(v))) notes <- c(notes, paste0(nm, " 全为 NA"))
      } else {
        idx[[nm]] <- NA_real_
        notes <- c(notes, paste0(nm, " 计算失败/不可用（列全为 NA）"))
      }
    }
  }
  tree <- md_core_tree(CFG)
  if ("PD" %in% want) {
    if (is.null(tree)) {
      idx$PD <- NA_real_
      notes <- c(notes, "无系统发育树（或 phylo$enabled=FALSE），Faith's PD 已跳过（全 NA）")
      md_log(CFG, "WARN", "无树 —— α 的 PD 跳过")
    } else {
      pv <- local_pd(m, tree)
      if (is.null(pv)) {
        idx$PD <- NA_real_
        notes <- c(notes, "Faith's PD 计算失败（树与 ASV 不匹配或退化），列全为 NA")
        md_warn(CFG, "Faith's PD 计算失败")
      } else {
        idx$PD <- as.numeric(pv[match(colnames(m), names(pv))])
      }
    }
  }

  # 未请求但用于下游的列不补；未在 want 里出现的指数直接不显示
  idx_idx <- intersect(want, names(idx))
  if (length(idx_idx) == 0) return(md_result(module, "SKIP", "没有任何可计算的 α 指数"))

  idx_m <- data.frame(idx, meta[, setdiff(names(meta), smp), drop = FALSE],
                      stringsAsFactors = FALSE)
  .mk("indices-done cols=", paste(names(idx_m), collapse = "/"))

  # ---- 2. 组间检验 + CLD ----
  tests <- CFG$alpha$test %||% "kruskal_dunn"
  p_adj <- CFG$alpha$p_adjust %||% "BH"
  min_n <- CFG$alpha$min_n %||% 3
  factors <- intersect(CFG$design$factors, names(meta))
  # 显式的空分组因子（design$factors = character(0)）表示"不分组"，不要回退到默认
  if (!is.null(CFG$design$factors) && length(CFG$design$factors) == 0) factors <- character(0)
  gvars <- list()
  for (f in factors) gvars[[f]] <- local_group_var(CFG, meta, f, name = f)
  if (length(factors) >= 2) {
    nm <- paste(factors, collapse = "_")
    gvars[[nm]] <- local_group_var(CFG, meta, factors[1:2], name = nm)
  }
  gvars <- Filter(function(g) !is.null(g) && length(unique(stats::na.omit(g$vals))) >= 2, gvars)

  test_rows <- list(); ph_rows <- list(); desc_rows <- list()
  for (gk in names(gvars)) {
    gv <- gvars[[gk]]
    idx_m$.g <- factor(gv$vals, levels = gv$levels)
    for (ix in idx_idx) {
      gt <- tryCatch(md_group_test(CFG, idx_m, ix, ".g", test = tests, p_adjust = p_adj,
                                   min_n = min_n, order_levels = gv$levels),
                     error = function(e) { md_warn(CFG, "md_group_test 失败 ", gk, "/", ix,
                                                   ": ", conditionMessage(e)); NULL })
      if (is.null(gt)) next
      if (!is.null(gt$test) && nrow(gt$test) > 0)
        test_rows[[length(test_rows) + 1L]] <- data.frame(
          factor = gk, index = ix, gt$test, stringsAsFactors = FALSE)
      if (!is.null(gt$posthoc) && nrow(gt$posthoc) > 0)
        ph_rows[[length(ph_rows) + 1L]] <- data.frame(
          factor = gk, index = ix, gt$posthoc, stringsAsFactors = FALSE)
      dd <- md_attach_letters(gt$desc, gt$letters)
      dd$factor <- gk; dd$index <- ix
      desc_rows[[length(desc_rows) + 1L]] <- dd
    }
  }
  idx_m$.g <- NULL
  .mk("tests-done")

  tab_test <- local_rbind(test_rows)
  if (!is.null(tab_test) && nrow(tab_test)) {
    tab_test$p_adj <- md_p_adjust(tab_test$p, p_adj)
    tab_test$sig <- ifelse(is.na(tab_test$p_adj), "", ifelse(tab_test$p_adj < 0.001, "***",
                       ifelse(tab_test$p_adj < 0.01, "**", ifelse(tab_test$p_adj < 0.05, "*", "ns"))))
  }
  tab_ph <- local_rbind(ph_rows)
  tab_desc <- local_rbind(desc_rows)
  .mk("tabs-done")
  if (!is.null(tab_ph)) tab_ph$sig <- ifelse(is.na(tab_ph$p_adj), "", ifelse(tab_ph$p_adj < 0.05, "*", "ns"))
  if (!is.null(tab_desc)) {
    tab_desc <- tab_desc[, c("factor", "index", "group", "n", "mean", "sd",
                             "se", "median", "q1", "q3", "min", "max", "letter"), drop = FALSE]
    num <- vapply(tab_desc, is.numeric, logical(1))
    num["n"] <- FALSE
    tab_desc[num] <- lapply(tab_desc[num], function(z) round(z, CFG$alpha$rounding %||% 3))
  }

  # ---- 3. 图：每个指数两因子箱线（A–F/G）----
  # 字号走**拼图轨道** MD_PT_COMP（12_plot_theme.R）—— 本模块把多个面板拼成一张图，
  # 与 local_grid()/theme_md_comp() 对 element_text 的处理保持同一条轨道
  # （原图右列小面板同样走 *_COMP：AX_TICK_PT / AX_TTL_PT / S_ANN_COMP）。
  pt <- md_pt(CFG, comp = TRUE)
  plots <- list()
  f1 <- factors[1]
  f2 <- if (length(factors) >= 2) factors[2] else NULL   # 注意 factors[2] 在单因子时是 NA 而非 NULL
  gv2 <- if (!is.null(f2)) local_group_var(CFG, meta, f2, name = f2) else NULL
  df_plot <- idx_m
  if (!is.null(f1)) df_plot[[f1]] <- factor(df_plot[[f1]],
                                            levels = local_factor_levels(CFG, meta, f1))
  if (!is.null(f2)) df_plot[[f2]] <- factor(df_plot[[f2]],
                                            levels = local_factor_levels(CFG, meta, f2))
  ix_plot <- utils::head(idx_idx, 6)
  if (length(ix_plot) == 0) ix_plot <- idx_idx[1]
  # 全 NA 的指数（如无树时的 PD）不出箱线面板，改在 Notes / 文字里说明
  ok_ix <- vapply(ix_plot, function(z) any(is.finite(idx_m[[z]])), logical(1))
  if (any(!ok_ix)) md_warn(CFG, "以下指数全为 NA，已跳过箱线面板：",
                           paste(ix_plot[!ok_ix], collapse = ", "))
  ix_plot <- ix_plot[ok_ix]
  if (length(ix_plot) == 0) {
    cand <- idx_idx[vapply(idx_idx, function(z) any(is.finite(idx_m[[z]])), logical(1))]
    ix_plot <- if (length(cand)) cand[1] else idx_idx[1]
  }

  for (ix in ix_plot) {
    .mk("loop ", ix)
    ylab <- local_index_label(ix)
    panel_tag <- LETTERS[length(plots) + 1L]
    if (!is.null(f1) && !is.null(f2) && length(unique(stats::na.omit(df_plot[[f1]]))) >= 2 &&
        length(unique(stats::na.omit(df_plot[[f2]]))) >= 2) {
      # 双因子 x 轴 = Zone|Vegetation 组合（组内 n >= min_n 才参与两两比较）
      gvc <- local_group_var(CFG, meta, c(f1, f2))
      dd <- df_plot
      dd$.inter <- factor(gvc$vals[match(dd[[smp]], names(gvc$vals))], levels = gvc$levels)
      gt <- tryCatch(md_group_test(CFG, dd, ix, ".inter", test = tests, p_adjust = p_adj,
                                   min_n = min_n, order_levels = gvc$levels),
                     error = function(e) NULL)
      gt2 <- tryCatch(md_group_test(CFG, df_plot, ix, f2, test = tests, p_adjust = p_adj,
                                    min_n = min_n,
                                    order_levels = local_factor_levels(CFG, meta, f2)),
                      error = function(e) NULL)
      cols <- md_pal(CFG, local_factor_levels(CFG, meta, f2), f2)
      dd[[f2]] <- as.character(dd[[f2]])
      kw <- local_kw_text(tab_test, c(f1, f2), ix)
      gp <- tryCatch(
        # inter_sep 必须与 local_group_var() 生成字母时用的分隔符一致（"_"），
        # 否则字母与 x 水平对不上而被丢弃（这是"箱线上只有样本量 4、没有字母"的根因）
        md_box_two_factor(dd, y = ix, f1 = f1, f2 = f2,
                          letters = if (!is.null(gt)) gt$letters else NULL,
                          colors = cols,
                          f1_levels = local_factor_levels(CFG, meta, f1),
                          f2_levels = local_factor_levels(CFG, meta, f2),
                          ylab = ylab, x_angle = 45, tag = panel_tag,
                          pt = pt, legend_pos = "none",
                          inter_sep = "_", x_labels = "two_line",
                          kw_text = kw),
        error = function(e) { md_warn(CFG, "箱线图失败 ", ix, ": ", conditionMessage(e));
                              md_empty_plot(CFG, ix) })
      plots[[length(plots) + 1L]] <- gp
    } else if (!is.null(f1)) {
      dd <- df_plot; dd[[f1]] <- as.character(dd[[f1]])
      cols <- md_pal(CFG, local_factor_levels(CFG, meta, f1), f1)
      gt1 <- tryCatch(md_group_test(CFG, df_plot, ix, f1, test = tests, p_adjust = p_adj,
                                    min_n = min_n,
                                    order_levels = local_factor_levels(CFG, meta, f1)),
                      error = function(e) NULL)
      gp <- tryCatch(
        md_box_group(dd, y = ix, x = f1, colors = cols, ylab = ylab,
                     xlab = NULL, letters = if (!is.null(gt1)) gt1$letters else NULL,
                     tag = panel_tag, pt = pt, legend_pos = "none",
                     kw_text = local_kw_text(tab_test, f1, ix)),
        error = function(e) md_empty_plot(CFG, ix))
      plots[[length(plots) + 1L]] <- gp
    }
  }

  # ---- 4. 稀释曲线 ----
  raw_curve <- tryCatch(md_core_asv(CFG, "rarefaction"), error = function(e) NULL)
  if (!is.null(raw_curve) && is.numeric(raw_curve) &&
      any(abs(raw_curve - round(raw_curve)) > 1e-9)) {
    md_log(CFG, "INFO", "rarefaction 口径非整数计数，已四舍五入后再算稀释曲线")
    raw_curve <- round(raw_curve)
  }
  kurve <- CFG$alpha$rarefaction$curve %||% "vegan"
  step <- as.integer(CFG$alpha$rarefaction$step %||% 200)
  nboot <- as.integer(CFG$alpha$rarefaction$nboot %||% 50)
  ci_on <- isTRUE(CFG$alpha$rarefaction$ci %||% TRUE)
  # 注：原图 Part C 的稀释曲线上**没有**抽平深度竖线，
  #     故本模块不再取 ref_depth、也不画 vline（只改画法，抽平深度的口径未动）。

  rare_df <- NULL; ci_df <- NULL; inext_ok <- FALSE
  if (kurve %in% c("inext", "both") && requireNamespace("iNEXT", quietly = TRUE)) {
    md_seed(CFG, 111)
    xi <- t(raw_curve)
    # 说明：iNEXT 对"极端单例化 + 样本多"的数据集很慢（实测 36 样本 x 321 ASV、
    # nboot=50 需数分钟），所以默认 curve="vegan"；here 只在显式选 inext/both 时跑。
    md_log(CFG, "INFO", "iNEXT 稀疏外推中（nboot = ", nboot, "）……")
    ine <- tryCatch(iNEXT::iNEXT(xi, q = 0, datatype = "abundance",
                                 nboot = nboot, conf = 0.95), error = function(e) NULL)
    if (!is.null(ine)) {
      rg <- tryCatch(as.data.frame(ine$iNextEst$size_based, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (!is.null(rg) && nrow(rg) > 0 && all(c("m", "qD") %in% names(rg))) {
        sid <- if ("Assemblage" %in% names(rg)) rg$Assemblage else if ("site" %in% names(rg)) rg$site else
          rep(NA_character_, nrow(rg))
        rare_df <- data.frame(sample_id = as.character(sid),
                              depth = as.numeric(rg$m), richness = as.numeric(rg$qD),
                              lower = if ("qD.LCL" %in% names(rg)) as.numeric(rg$qD.LCL) else NA_real_,
                              upper = if ("qD.UCL" %in% names(rg)) as.numeric(rg$qD.UCL) else NA_real_,
                              stringsAsFactors = FALSE)
        inext_ok <- TRUE
      }
    }
    if (!inext_ok) md_warn(CFG, "iNEXT 稀释曲线不可用，回退 vegan::rarecurve")
  }
  if (is.null(rare_df)) {
    rare_df <- tryCatch({
      # ⚠ vegan::rarecurve() **即使 tidy = TRUE 也会打开图形设备** —— Rscript 下没有
      #   设备时就会在包根掉一个 Rplots.pdf。所以必须套 md_offscreen()（见 pitfalls C39）。
      rc <- md_offscreen(vegan::rarecurve(t(raw_curve), step = step,
                             sample = min(colSums(raw_curve)), label = FALSE,
                             tidy = TRUE))
      data.frame(sample_id = as.character(rc$Site), depth = as.numeric(rc$Sample),
                 richness = as.numeric(rc$Species), lower = NA_real_, upper = NA_real_,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      md_warn(CFG, "vegan::rarecurve(tidy=TRUE) 失败：", conditionMessage(e))
      NULL
    })
    if (is.null(rare_df) || nrow(rare_df) == 0) {
      rare_df <- tryCatch({
        # 这条不是 tidy 版：rarecurve 会**真的画一遍**再把每条曲线当返回值给出。
        # 图我们不要，但副作用躲不掉 —— 同样套 md_offscreen()，让它画进不落盘的设备。
        rl <- md_offscreen(vegan::rarecurve(t(raw_curve), step = step,
                               sample = min(colSums(raw_curve)), label = FALSE))
        do.call(rbind, lapply(seq_along(rl), function(i) {
          data.frame(sample_id = rownames(raw_curve)[i],
                     depth = as.numeric(attr(rl[[i]], "Subsample")),
                     richness = as.numeric(rl[[i]]), lower = NA_real_, upper = NA_real_,
                     stringsAsFactors = FALSE)
        }))
      }, error = function(e) { md_warn(CFG, "vegan::rarecurve(list) 失败：", conditionMessage(e)); NULL })
    }
    if (is.null(rare_df) || nrow(rare_df) == 0) {
      md_warn(CFG, "稀释曲线引擎全部不可用，面板降级为占位图")
      rare_df <- NULL
    }
  }
  if (!is.null(rare_df) && length(gvars)) {
    gv1 <- gvars[[1]]
    rare_df$group <- as.character(gv1$vals[match(rare_df$sample_id, names(gv1$vals))])
    if (ci_on && any(!is.na(rare_df$group))) {
      brk <- sort(unique(rare_df$depth))
      ci_list <- lapply(split(rare_df, rare_df$group), function(z) {
        rows <- lapply(brk, function(d) {
          v <- z$richness[z$depth == d]
          v <- v[is.finite(v)]
          n <- length(v)
          if (n < 3) return(NULL)
          se <- stats::sd(v) / sqrt(n)
          data.frame(group = z$group[1], depth = d, mean = mean(v),
                     lower = mean(v) - stats::qt(0.975, n - 1) * se,
                     upper = mean(v) + stats::qt(0.975, n - 1) * se,
                     stringsAsFactors = FALSE)
        })
        rows <- Filter(Negate(is.null), rows)
        if (length(rows) == 0) return(NULL)
        do.call(rbind, rows)
      })
      ci_list <- Filter(Negate(is.null), ci_list)
      ci_df <- if (length(ci_list)) do.call(rbind, ci_list) else NULL
    }
  }
  md_log(CFG, "INFO", "稀释曲线：", if (is.null(rare_df)) "无" else nrow(rare_df), " 点 | CI 行数 ",
         if (is.null(ci_df)) 0 else nrow(ci_df))

  cols1 <- if (length(gvars)) md_pal(CFG, gvars[[1]]$levels, gvars[[1]]$name) else NULL
  # 稀释曲线画法沿用同一套样式 Part C（L596-643）：
  #   geom_line(linewidth = 0.3 (= MD_LW$line * 0.75), alpha = 0.7 (= MD_PT_SIZE$alpha))
  #   均值线 linewidth = REG_LINE_LW (= MD_LW$line)、置信带 alpha = REG_RIBBON_A (= MD_REG$ribbon_alpha)、colour = NA
  #   labs(x = "Reads sampled", y = "ASV richness")、图例置底
  #   **原图没有"抽平深度"竖线**，故不画（原模块自加的 geom_vline 已移除）。
  pR <- if (is.null(rare_df)) {
    md_empty_plot(CFG, "Rarefaction curve not available")
  } else {
    g <- ggplot2::ggplot(rare_df, ggplot2::aes(x = .data$depth, y = .data$richness,
                                               group = .data$sample_id,
                                               colour = .data$group)) +
      ggplot2::geom_line(linewidth = MD_LW$line * 0.75, alpha = MD_PT_SIZE$alpha)
    if (!is.null(ci_df) && nrow(ci_df) > 0) {
      g <- g + ggplot2::geom_ribbon(data = ci_df,
                                    ggplot2::aes(x = .data$depth, y = .data$mean,
                                                 ymin = .data$lower, ymax = .data$upper,
                                                 fill = .data$group, group = .data$group),
                                    inherit.aes = FALSE, alpha = MD_REG$ribbon_alpha,
                                    colour = NA) +
        ggplot2::geom_line(data = ci_df,
                           ggplot2::aes(x = .data$depth, y = .data$mean, colour = .data$group),
                           inherit.aes = FALSE, linewidth = MD_LW$line)
    }
    if (!is.null(cols1)) {
      # scale_fill 只在真的画了置信带时加，否则 ggplot 会报
      # "No shared levels found between names(values) ... and the data's fill values"
      g <- g + ggplot2::scale_colour_manual(values = cols1, name = gvars[[1]]$name)
      if (!is.null(ci_df) && nrow(ci_df) > 0)
        g <- g + ggplot2::scale_fill_manual(values = cols1, name = gvars[[1]]$name)
    }
    g + ggplot2::labs(x = "Reads sampled", y = "ASV richness",
                      tag = LETTERS[length(plots) + 1]) +
      theme_md(CFG = NULL) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(colour = md_guide_flex(length(cols1 %||% character(1)), "bottom"))
  }

  # ---- 5. 拼图 ----
  n_pan <- length(plots)
  if (n_pan == 0) plots <- list(md_empty_plot(CFG, "No index available"))
  body <- local_grid(plots, theme_md_comp(CFG))
  fig <- if (!is.null(body)) body / pR + patchwork::plot_layout(heights = c(3, 1.15)) else pR
  h <- if (length(ix_plot) >= 6) 190 else if (length(ix_plot) >= 4) 165 else 130
  md_save_figure(CFG, fig, "11_alpha", height = h)

  # ---- 6. 表 ----
  tab_rare <- if (!is.null(ci_df) && nrow(ci_df)) {
    ci_df$factor <- names(gvars)[1]
    ci_df[, c("factor", "group", "depth", "mean", "lower", "upper")]
  } else NULL
  note <- data.frame(
    item = c("Normalization (alpha)", "Normalization (rarefaction)", "Rarefaction depth (alpha)",
             "Test", "p adjust", "min_n", "Pielou definition", "PD", "Rarefaction engine",
             "Rarefaction engine used", "Curve step", "n samples"),
    value = c(as.character(CFG$normalize$by_module$alpha %||% CFG$normalize$method),
              as.character(CFG$normalize$by_module$rarefaction %||% CFG$normalize$method),
              as.character(CFG$normalize$depth %||% "auto (min positive sample reads)"),
              as.character(tests), as.character(p_adj), as.character(min_n),
              "Shannon / log(Observed)", if ("PD" %in% idx_idx && all(is.na(idx$PD))) "skipped (no tree)" else "Faith's PD (picante::pd, include.root=TRUE)",
              as.character(kurve), if (inext_ok) "iNEXT" else if (!is.null(rare_df)) "vegan::rarecurve" else "none",
              as.character(step), as.character(ncol(m))),
    stringsAsFactors = FALSE)
  if (length(notes)) note <- rbind(note, data.frame(item = "Notes", value = paste(notes, collapse = " | ")))

  tp <- md_write_xlsx(CFG, "11_alpha", list(
    Descriptive = tab_desc,
    Tests = tab_test,
    Posthoc = tab_ph,
    Alpha_values = if (isTRUE(CFG$table$detail)) idx_m else NULL,
    Rarefaction_curve = tab_rare,
    Notes = note))

  md_result(module, "OK", tables = basename(tp), figures = "figures/11_alpha.pdf",
            notes = paste0("indices=", paste(idx_idx, collapse = "+"),
                           "; test=", tests,
                           if (length(notes)) paste0("; ", paste(notes, collapse = "; ")) else ""))
}

# ---------------------------------------------------------------------------
# local_marker()：模块内轻量调试钩子
# 默认**完全静默**；把 CFG$opts$debug 设为 TRUE 才会把进度打到日志。
# 保留它是为了让排查长流程（指数 -> 检验 -> 事后 -> 稀释曲线）时能定位卡在哪一步。
# ---------------------------------------------------------------------------
local_marker <- function(CFG) {
  dbg <- isTRUE(CFG$opts$debug)
  function(...) if (dbg) md_log(CFG, "DEBUG", "11_alpha: ", ...)
}
