# =============================================================================
# R/05_assemble.R —— 数据装配与标准化（全下游共用的核心对象）
# -----------------------------------------------------------------------------
# 流程：读五文件 -> 过滤（ASV 读长/出现率、样本读长、类群子集）
#       -> 标准化（逐模块口径，自实现抽平）-> Good's coverage
#       -> 落 cache/05_core.rds
#
# **本模块只做"把原始数据变成规范丰度表"，不计算任何多样性指数或距离。**
# 下游读写一律走访问器：md_core() / md_core_asv(CFG, "alpha") / md_core_analysis()
#
# 标准化口径（CFG$normalize$by_module）：
#   composition / rarefaction -> none（未抽平）
#   alpha / distances / phylo / assembly / partition / env / model -> rarefy
# =============================================================================

md_assemble <- function(CFG) {
  module <- "05"
  if (!md_need(CFG, "io", module)) return(md_result(module, "SKIP", "缺少 openxlsx"))

  inp <- md_read_inputs(CFG)
  asv <- inp$asv; meta <- inp$meta; env <- inp$env; tree <- inp$tree; tax <- inp$tax
  smp <- inp$sample_col
  filt_log <- list()
  # step      = 中文步骤名（给人看结果表用）
  # step_label= **纯 ASCII** 步骤名（**只有它允许进图**；图中文字必须纯 ASCII，见 pitfalls B12）
  add_log <- function(step, step_label, n_before, n_after, unit) {
    filt_log[[length(filt_log) + 1L]] <<- data.frame(
      step = step, step_label = step_label,
      n_before = n_before, n_after = n_after,
      removed = n_before - n_after, unit = unit, stringsAsFactors = FALSE)
  }

  # ---- 1. ASV 读长过滤 ----
  n0 <- nrow(asv)
  asv <- asv[rowSums(asv) >= CFG$filter$min_asv_reads, , drop = FALSE]
  add_log(paste0("ASV 总读长 >= ", CFG$filter$min_asv_reads),
          paste0("ASV reads >= ", CFG$filter$min_asv_reads), n0, nrow(asv), "ASV")

  # ---- 2. ASV 出现率过滤 ----
  n0 <- nrow(asv)
  asv <- asv[rowSums(asv > 0) >= CFG$filter$min_prev_samples, , drop = FALSE]
  add_log(paste0("出现样本数 >= ", CFG$filter$min_prev_samples),
          paste0("Occurrence >= ", CFG$filter$min_prev_samples), n0, nrow(asv), "ASV")

  # ---- 3. 样本读长过滤 ----
  n0 <- ncol(asv)
  cs <- colSums(asv)
  keep_s <- names(cs)[cs >= CFG$filter$min_sample_reads]
  if (length(keep_s) == 0) stop("所有样本读长都低于 min_sample_reads", call. = FALSE)
  asv <- asv[, keep_s, drop = FALSE]
  meta <- meta[match(keep_s, meta[[smp]]), , drop = FALSE]
  if (!is.null(env)) env <- env[match(keep_s, env[[smp]]), , drop = FALSE]
  rownames(meta) <- keep_s
  if (!is.null(env)) rownames(env) <- keep_s
  add_log(paste0("样本读长 >= ", CFG$filter$min_sample_reads),
          paste0("Sample reads >= ", CFG$filter$min_sample_reads), n0, ncol(asv), "样本")

  # ---- 4. 类群子集过滤（通用化的类群筛选机制）----
  sub_log <- NULL
  if (!is.null(CFG$filter$taxon_subset) && !is.null(tax)) {
    ts <- CFG$filter$taxon_subset
    rank <- ts$rank %||% names(tax)[2]
    ids  <- rownames(asv)
    x <- as.character(tax[[rank]][match(ids, tax[[inp$id_col]])])
    x[is.na(x)] <- "Unclassified"
    drop_ids <- character(0)
    if (!is.null(ts$drop)) {
      hit <- Reduce(`|`, lapply(ts$drop, function(k) grepl(k, x, ignore.case = TRUE)))
      drop_ids <- c(drop_ids, ids[hit])
    }
    if (!is.null(ts$keep)) {
      hit <- Reduce(`|`, lapply(ts$keep, function(k) grepl(k, x, ignore.case = TRUE)))
      drop_ids <- c(drop_ids, ids[!hit])
    }
    drop_ids <- unique(drop_ids)
    n0 <- nrow(asv)
    asv <- asv[setdiff(ids, drop_ids), , drop = FALSE]
    add_log(paste0("类群子集过滤 @", rank),
            paste0("Taxon subset @", rank), n0, nrow(asv), "ASV")
    sub_log <- data.frame(rank = rank,
                          keep = paste(ts$keep %||% "—", collapse = ","),
                          drop = paste(ts$drop %||% "—", collapse = ","),
                          n_removed = n0 - nrow(asv), stringsAsFactors = FALSE)
    if (nrow(asv) < 10)
      md_warn(CFG, "类群子集过滤后只剩 ", nrow(asv), " 个 ASV —— 检查 filter$taxon_subset 是否过严")
  }
  if (nrow(asv) < 5) stop("过滤后 ASV 少于 5 个，无法继续", call. = FALSE)

  # ---- 5. 标准化 ----
  md_seed(CFG, 201)
  raw <- asv
  method <- CFG$normalize$method
  depth <- CFG$normalize$depth
  if (is.null(depth)) depth <- min(colSums(raw)[colSums(raw) > 0])
  depth <- as.integer(depth)
  if (depth < 1) stop("抽平深度无效", call. = FALSE)
  if (depth > min(colSums(raw)))
    stop(sprintf("指定的抽平深度 %d 大于最小样本读长 %d", depth, min(colSums(raw))), call. = FALSE)

  asv_rare <- raw; asv_subavg <- NULL
  if (method == "rarefy") {
    md_seed(CFG, 202)
    asv_rare <- md_rarefy(raw, depth, seed = CFG$opts$seed)
    asv_rare <- asv_rare[rowSums(asv_rare) > 0, , drop = FALSE]
  } else if (method == "subsample_avg") {
    reps <- CFG$normalize$subsample_reps %||% 20
    md_log(CFG, "INFO", "多次亚抽样取均值：", reps, " 次 x 深度 ", depth)
    asv_subavg <- md_subsample_avg(raw, depth, reps, seed = CFG$opts$seed)
    asv_rare <- asv_subavg
  } else if (method == "relative") {
    asv_rare <- md_rel(raw)
  } else if (method == "none") {
    asv_rare <- raw
  }

  asv_css <- NULL
  if (method == "css" || "css" %in% unlist(CFG$normalize$by_module)) {
    asv_css <- tryCatch(md_css(raw), error = function(e) {
      md_warn(CFG, "CSS 不可用（", conditionMessage(e), "），回退相对丰度")
      md_rel(raw)
    })
  }

  # ---- 6. Good's coverage ----
  cov_df <- data.frame(
    sample_id = colnames(raw),
    reads     = as.integer(colSums(raw)),
    n_asv     = as.integer(colSums(raw > 0)),
    singletons = as.integer(colSums(raw == 1)),
    goods_coverage = 1 - colSums(raw == 1) / pmax(colSums(raw), 1),
    stringsAsFactors = FALSE)
  cov_df <- merge(cov_df, meta, by = smp, all.x = TRUE, sort = FALSE)
  cov_df <- cov_df[match(colnames(raw), cov_df[[smp]]), , drop = FALSE]

  # ---- 7. 落盘 ----
  core <- list(
    asv_raw = raw, asv_rare = asv_rare, asv_css = asv_css, asv_subavg = asv_subavg,
    depth = depth, normalize_method = method, seed = CFG$opts$seed,
    meta = meta, env = env, tree = tree, tax = tax,
    sample_col = smp, id_col = inp$id_col,
    coverage = cov_df,
    filter_log = do.call(rbind, filt_log),
    taxon_subset_log = sub_log,
    input_report = inp$report
  )
  saveRDS(core, md_path(CFG, "cache", "05_core.rds"))
  md_log(CFG, "INFO", "核心对象：", nrow(raw), " ASV x ", ncol(raw), " 样本",
         " | 抽平深度 ", depth, " | 口径 ", method)

  # ---- 8. 出图（x 轴用纯 ASCII 的 step_label；中文 step 只留在结果表里）----
  ic <- data.frame(stage = c("Raw", core$filter_log$step_label),
                   n_asv = c(inp$report$n_asv, core$filter_log$n_after),
                   stringsAsFactors = FALSE)
  pA <- ggplot2::ggplot(ic, ggplot2::aes(x = stats::reorder(stage, -n_asv), y = n_asv)) +
    ggplot2::geom_col(fill = "#4E79A7", width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = n_asv), vjust = -0.4, size = pt2gg(MD_PT$ann * 0.9)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = NULL, y = "Number of ASVs", tag = "A") +
    theme_md(CFG = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  pB <- ggplot2::ggplot(cov_df, ggplot2::aes(x = reads, y = goods_coverage)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.55, colour = "#59A14F") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Reads per sample (log10)", y = "Good's coverage", tag = "B") +
    theme_md(CFG = NULL)
  pC <- ggplot2::ggplot(cov_df, ggplot2::aes(y = goods_coverage)) +
    ggplot2::geom_histogram(bins = 20, fill = "#8FAE5B", colour = "white", linewidth = 0.15) +
    ggplot2::labs(x = "Good's coverage", y = "Count", tag = "C") + theme_md(CFG = NULL)
  md_save_figure(CFG, md_compose(CFG, list(A = pA, B = pB, C = pC), design = "ABB\nCCC"),
                 "05_assemble", height = 120)

  f1 <- CFG$design$factors[1]
  tp <- md_write_xlsx(CFG, "05_assemble", list(
    Filter_log = core$filter_log,
    Taxon_subset = sub_log,
    Normalization = data.frame(
      item = c("全局默认口径",
               "composition", "rarefaction", "alpha", "distances", "phylo",
               "assembly", "partition", "env", "model",
               "抽平深度", "抽平前最小样本读长", "seed"),
      value = c(CFG$normalize$method,
                vapply(c("composition", "rarefaction", "alpha", "distances", "phylo",
                         "assembly", "partition", "env", "model"),
                       function(k) as.character(CFG$normalize$by_module[[k]] %||% "—"), character(1)),
                depth, min(colSums(raw)), CFG$opts$seed),
      stringsAsFactors = FALSE),
    Coverage = cov_df,
    Factor_sizes = if (!is.null(f1)) as.data.frame(table(meta[[f1]])) else NULL,
    Depth_sensitivity = md_depth_sensitivity(CFG, raw)))

  md_result(module, "OK", tables = basename(tp), figures = "figures/05_assemble.pdf")
}

#' 抽平深度敏感性：不同深度下保留的 ASV 数（供挑选深度时参考）
md_depth_sensitivity <- function(CFG, raw) {
  cs <- colSums(raw)
  cand <- unique(as.integer(stats::quantile(cs, probs = c(0, .05, .1, .25, .5))))
  cand <- sort(cand[cand >= 1 & cand <= min(cs)])
  if (length(cand) == 0) return(NULL)
  do.call(rbind, lapply(cand, function(d) {
    md_seed(CFG, 400 + d)
    r <- tryCatch(md_rarefy(raw, d, seed = CFG$opts$seed + d), error = function(e) NULL)
    data.frame(depth = d,
               reads_used_pct = d / stats::median(cs),
               n_asv_retained = if (is.null(r)) NA_integer_ else sum(rowSums(r) > 0),
               asv_retained_pct = if (is.null(r)) NA_real_ else mean(rowSums(r) > 0),
               stringsAsFactors = FALSE)
  }))
}

#' CSS（累积和标度）标准化。缺 metagenomeSeq 时由调用方回退。
md_css <- function(m) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE))
    stop("未安装 metagenomeSeq", call. = FALSE)
  x <- metagenomeSeq::newMRexperiment(m)
  x <- metagenomeSeq::cumNorm(x)
  as.matrix(metagenomeSeq::MRcounts(x, norm = TRUE))
}
