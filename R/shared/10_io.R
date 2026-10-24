# =============================================================================
# R/shared/10_io.R —— 输入读取 / 对齐 / 写表 / 格式化器 / 数据护栏
# -----------------------------------------------------------------------------
# 库文件，无分析逻辑。
# 输入契约见 references/00-data-contract.md。
# =============================================================================

# ---------------------------------------------------------------------------
# 数据护栏
# ---------------------------------------------------------------------------
#' 禁止向 data/ 目录写入任何东西（源数据永不改写）
md_assert_not_data <- function(CFG, path) {
  d <- normalizePath(CFG$paths$abs$data, mustWork = FALSE)
  p <- normalizePath(path, mustWork = FALSE)
  if (startsWith(p, d)) {
    stop("禁止写入 data/ 目录（源数据只读）：", p, call. = FALSE)
  }
  invisible(TRUE)
}

#' 输入五文件的 md5，用于缓存签名与"数据未被改写"校验
md_data_md5 <- function(CFG) {
  keys <- c("asv", "taxonomy", "meta", "env", "tree")
  out <- list()
  for (k in keys) {
    p <- md_data_path(CFG, k)
    if (is.null(p)) { out[[k]] <- "NULL"; next }
    out[[k]] <- if (file.exists(p)) unname(tools::md5sum(p)) else "MISSING"
  }
  out
}

# ---------------------------------------------------------------------------
# 读取
# ---------------------------------------------------------------------------
.md_read_csv <- function(p) {
  if (is.null(p) || !file.exists(p)) return(NULL)
  utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE,
                  fileEncoding = "UTF-8")
}

.md_read_nwk <- function(p, pkgs_ok) {
  if (is.null(p) || !file.exists(p) || !pkgs_ok) return(NULL)
  tr <- ape::read.tree(p)
  tr$tip.label <- gsub('^["\']|["\']$', "", trimws(tr$tip.label))
  tr
}

#' 读取并对齐五个输入文件。
#' 返回 list(asv, tax, meta, env, tree, samples, asvs, report)
md_read_inputs <- function(CFG, quiet = FALSE) {
  p_asv  <- md_data_path(CFG, "asv")
  p_tax  <- md_data_path(CFG, "taxonomy")
  p_meta <- md_data_path(CFG, "meta")
  p_env  <- md_data_path(CFG, "env")
  p_tree <- md_data_path(CFG, "tree")
  if (!file.exists(p_asv))  stop("丰度表不存在: ", p_asv, call. = FALSE)
  if (!file.exists(p_meta)) stop("样本元数据不存在: ", p_meta, call. = FALSE)
  if (!file.exists(p_tax))  stop("taxonomy 表不存在: ", p_tax, call. = FALSE)

  id_col  <- CFG$data$id_col
  smp_col <- CFG$data$sample_col

  d_asv <- .md_read_csv(p_asv)
  if (!id_col %in% names(d_asv)) {
    id_col <- names(d_asv)[1]
    if (!quiet) message("  [io] 丰度表首列名不是 '", CFG$data$id_col, "'，改用 '", id_col, "'")
  }
  asvs <- trimws(as.character(d_asv[[id_col]]))
  d_asv[[id_col]] <- NULL
  m <- as.matrix(d_asv)
  storage.mode(m) <- "double"
  rownames(m) <- asvs
  if (anyDuplicated(asvs)) stop("丰度表 asv_id 有重复", call. = FALSE)

  tax <- .md_read_csv(p_tax)
  if (!is.null(tax)) {
    if (!id_col %in% names(tax)) names(tax)[1] <- id_col
    tax[[id_col]] <- trimws(as.character(tax[[id_col]]))
  }

  meta <- .md_read_csv(p_meta)
  if (!smp_col %in% names(meta)) {
    smp_col <- names(meta)[1]
    if (!quiet) message("  [io] 元数据首列名不是 '", CFG$data$sample_col, "'，改用 '", smp_col, "'")
  }
  meta[[smp_col]] <- trimws(as.character(meta[[smp_col]]))
  if (!quiet && smp_col != CFG$data$sample_col)
    meta[[CFG$data$sample_col]] <- meta[[smp_col]]

  env <- .md_read_csv(p_env)
  if (!is.null(env)) {
    if (!smp_col %in% names(env)) names(env)[1] <- smp_col
    env[[smp_col]] <- trimws(as.character(env[[smp_col]]))
  }

  pkgs_ok <- requireNamespace("ape", quietly = TRUE)
  tree <- .md_read_nwk(p_tree, pkgs_ok)

  # ---- 对齐样本 ----
  in_asv  <- colnames(m)
  in_meta <- meta[[smp_col]]
  common  <- intersect(in_asv, in_meta)
  only_asv  <- setdiff(in_asv, in_meta)
  only_meta <- setdiff(in_meta, in_asv)
  if (length(common) == 0) stop("丰度表列名与元数据 sample_id 没有任何交集", call. = FALSE)

  meta <- meta[match(common, meta[[smp_col]]), , drop = FALSE]
  m    <- m[, common, drop = FALSE]
  if (!is.null(env)) env <- env[match(common, env[[smp_col]]), , drop = FALSE]
  rownames(meta) <- common
  if (!is.null(env)) rownames(env) <- common

  # ---- 树 tip 对齐 ----
  tree_cov <- NA_real_
  if (!is.null(tree)) {
    tree_cov <- mean(rownames(m) %in% tree$tip.label)
    if (tree_cov < 0.5) {
      md_warn(CFG, "树的 tip 只覆盖 ", sprintf("%.1f%%", 100 * tree_cov),
              " 的 ASV —— 系统发育模块建议关闭（phylo$enabled=FALSE）")
    }
  }

  report <- list(
    n_asv = nrow(m), n_sample = ncol(m),
    samples_only_in_asv = only_asv, samples_only_in_meta = only_meta,
    tree_tip_cov = tree_cov,
    tax_ranks_present = if (is.null(tax)) character(0) else
      intersect(CFG$data$tax_ranks, names(tax))
  )
  if (!quiet) {
    md_log(CFG, "INFO", "读入: ", report$n_asv, " ASV x ", report$n_sample, " 样本",
           if (!is.null(tree)) paste0(" | 树 tip ", length(tree$tip.label),
                                      "（覆盖 ", sprintf("%.1f%%", 100 * tree_cov), "）") else " | 无树")
    if (length(only_asv))  md_warn(CFG, "丰度表有 ", length(only_asv), " 个样本不在元数据中，已剔除")
    if (length(only_meta)) md_warn(CFG, "元数据有 ", length(only_meta), " 个样本不在丰度表中，已剔除")
  }
  list(asv = m, tax = tax, meta = meta, env = env, tree = tree,
       samples = common, asvs = rownames(m),
       id_col = id_col, sample_col = smp_col, report = report)
}

# ---------------------------------------------------------------------------
# 丰度表工具
# ---------------------------------------------------------------------------
#' 相对丰度（按列归一）
md_rel <- function(m) sweep(m, 2, pmax(colSums(m), 1), "/")

#' 按分类层级聚合（求和）
md_glom <- function(m, tax, rank, id_col = "asv_id") {
  if (is.null(tax) || !rank %in% names(tax)) {
    stop("taxonomy 表里没有层级: ", rank, call. = FALSE)
  }
  g <- as.character(tax[[rank]][match(rownames(m), tax[[id_col]])])
  g[is.na(g) | g == ""] <- "Unclassified"
  agg <- rowsum(m, group = g, reorder = TRUE)
  list(mat = agg, group = rownames(agg))
}

#' 矩阵 -> 长表
md_asv_long <- function(m, meta = NULL, sample_col = "sample_id", id_col = "asv_id") {
  d <- as.data.frame(m, check.names = FALSE)
  d[[id_col]] <- rownames(m)
  long <- tidyr::pivot_longer(d, cols = -dplyr::all_of(id_col),
                              names_to = sample_col, values_to = "abundance")
  if (!is.null(meta)) long <- dplyr::left_join(long, meta, by = sample_col)
  long
}

#' ASV 表：总丰度 / 出现样本数 / 出现率 / 平均相对丰度
md_asv_stats <- function(m, meta = NULL) {
  rel <- md_rel(m)
  data.frame(
    asv_id        = rownames(m),
    total_reads   = rowSums(m),
    mean_rel      = rowMeans(rel),
    n_present     = rowSums(m > 0),
    prev_rate     = rowSums(m > 0) / ncol(m),
    max_rel       = apply(rel, 1, max),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 写表
# ---------------------------------------------------------------------------
#' 写多 sheet xlsx 到 results/。sheets = 具名 list(data.frame)。
#' 自动追加一个 `_meta` sheet（运行信息），保证每个结果表可回溯。
md_write_xlsx <- function(CFG, name, sheets, meta_extra = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    md_log(CFG, "WARN", "缺 openxlsx，跳过写出 ", name)
    return(invisible(NULL))
  }
  sheets <- Filter(Negate(is.null), sheets)
  # 丢掉零行 sheet：合法的"本论文没有这类结果"不应该让 95 的产物校验报 FAIL。
  # 丢掉的名字会记进 `_meta`，可追溯。
  empty <- names(sheets)[vapply(sheets, function(x) nrow(as.data.frame(x)) == 0, logical(1))]
  if (length(empty)) {
    md_log(CFG, "INFO", "  跳过 ", length(empty), " 个空 sheet（无结果）: ",
           paste(empty, collapse = ", "))
    sheets <- sheets[!names(sheets) %in% empty]
  }
  if (length(sheets) == 0) return(invisible(NULL))
  sheets <- lapply(sheets, function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x[] <- lapply(x, function(col) {
      if (is.numeric(col)) col else as.character(col)
    })
    x
  })
  mi <- md_session_info(CFG)
  if (length(empty)) {
    mi <- rbind(mi, data.frame(item = "empty_sheets_dropped",
                               value = paste(empty, collapse = ", "),
                               stringsAsFactors = FALSE))
  }
  if (!is.null(meta_extra)) {
    mi <- rbind(mi, data.frame(item = names(meta_extra),
                               value = as.character(unlist(meta_extra)),
                               stringsAsFactors = FALSE))
  }
  sheets[["_meta"]] <- mi

  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, substr(nm, 1, 31))
    openxlsx::writeData(wb, substr(nm, 1, 31), sheets[[nm]], keepNA = TRUE, na.string = "NA")
    openxlsx::freezePane(wb, substr(nm, 1, 31), firstRow = TRUE)
  }
  p <- md_path(CFG, "results", paste0(name, ".xlsx"))
  md_assert_not_data(CFG, p)
  openxlsx::saveWorkbook(wb, p, overwrite = TRUE)
  md_log(CFG, "OK", "写出 ", basename(p), "（", length(sheets), " 个 sheet）")
  if (isTRUE(CFG$table$write_csv)) {
    for (nm in names(sheets)) {
      utils::write.csv(sheets[[nm]],
                       file.path(CFG$paths$abs$results, paste0(name, "__", nm, ".csv")),
                       row.names = FALSE, fileEncoding = "UTF-8")
    }
  }
  p
}

#' 按 **sheet 名** 读 xlsx（绝不按 sheetN 序号对齐）
md_read_xlsx <- function(path, sheet) {
  stopifnot(file.exists(path))
  openxlsx::read.xlsx(path, sheet = sheet)
}

# ---------------------------------------------------------------------------
# 格式化器（写论文表格用）
# ---------------------------------------------------------------------------
md_fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "—", formatC(round(as.numeric(x), digits),
                                format = "f", digits = digits))
}

md_fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "—",
         ifelse(p < 0.001, "< 0.001",
                ifelse(p < 0.01, formatC(p, format = "f", digits = 3),
                       formatC(p, format = "f", digits = digits))))
}

#' **图内用**的 P 值文案（纯 ASCII）—— 对应 nifH 原脚本
#'
#' 与 `md_fmt_p()` 的区别（`md_fmt_p` 是给结果表用的，含非 ASCII 破折号且无前缀）：
#'  - NA → `"NA"`（不是 U+2014 破折号，否则图会被 95 的 ASCII 审计判 FAIL）
#'  - 带前缀：`"P < 0.001"` / `"P = 0.020"` / `"P = 0.512"`
md_fmt_p_fig <- function(p, digits = 3, prefix = TRUE) {
  body <- ifelse(is.na(p), "NA",
                 ifelse(p < 0.001, "< 0.001",
                        formatC(p, format = "f", digits = digits)))
  if (!isTRUE(prefix)) return(body)
  ifelse(is.na(p), "P = NA",
         ifelse(p < 0.001, "P < 0.001", paste0("P = ", body)))
}

md_fmt_ms <- function(mean, sd, digits = 2) {
  ifelse(is.na(mean), "—",
         ifelse(is.na(sd), md_fmt_num(mean, digits),
                paste0(md_fmt_num(mean, digits), " ± ", md_fmt_num(sd, digits))))
}

md_fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "—", paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

# ---------------------------------------------------------------------------
# 核心对象访问器（05_assemble.R 的产物，全下游共用）
# ---------------------------------------------------------------------------
#' 读取 cache/05_core.rds。缺失时给出明确指引（而不是 NULL 引发下游怪错）。
md_core <- function(CFG) {
  f <- md_path(CFG, "cache", "05_core.rds")
  if (!file.exists(f)) {
    stop("核心对象不存在：", f,
         "\n请先运行模块 05（Rscript --vanilla R/90_run_all.R --only=01,05）",
         call. = FALSE)
  }
  readRDS(f)
}

#' 取某模块口径下的 ASV 矩阵。
#' @param module_key normalize$by_module 的键名，如 "alpha" / "composition"
#' @return 矩阵（ASV x 样本，列为原始 reads 或相对丰度 / 抽平值，取决于口径）
md_core_asv <- function(CFG, module_key = "alpha") {
  core <- md_core(CFG)
  method <- (CFG$normalize$by_module[[module_key]]) %||% CFG$normalize$method
  switch(method,
         none     = core$asv_raw,
         relative = md_rel(core$asv_raw),
         rarefy   = core$asv_rare,
         css      = core$asv_css %||% core$asv_raw,
         subsample_avg = core$asv_subavg %||% core$asv_rare,
         core$asv_rare)
}

#' 取某模块口径下的样本元数据（已按 ASV 矩阵列顺序对齐）
md_core_meta <- function(CFG) {
  core <- md_core(CFG)
  core$meta
}

#' 取某模块口径下的环境/响应变量表（已对齐）
md_core_env <- function(CFG) {
  core <- md_core(CFG)
  core$env
}

#' 合并 meta + env 的分析用总表（行 = 样本，列 = 因子 + 环境 + 响应）
md_core_analysis <- function(CFG) {
  core <- md_core(CFG)
  m <- core$meta
  e <- core$env
  if (is.null(e)) return(m)
  smp <- CFG$data$sample_col
  drop <- intersect(names(e), names(m))
  drop <- setdiff(drop, smp)
  e <- e[, setdiff(names(e), drop), drop = FALSE]
  out <- merge(m, e, by = smp, all.x = TRUE, sort = FALSE)
  out[match(m[[smp]], out[[smp]]), , drop = FALSE]
}

#' 取树（无树或 phylo$enabled=FALSE 时返回 NULL）
md_core_tree <- function(CFG) {
  if (!isTRUE(CFG$phylo$enabled)) return(NULL)
  core <- md_core(CFG)
  core$tree
}

#' 把核心对象（缓存）里的样本顺序对齐到给定 id 向量
md_align <- function(df_or_vec, ids, col = NULL) {
  if (is.null(col)) return(df_or_vec[ids])
  df_or_vec[match(ids, df_or_vec[[col]]), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 跨模块产物的访问器（**模块间唯一接口**，不要绕过去直接 readRDS）
# ---------------------------------------------------------------------------
.md_cache_get <- function(CFG, file, producer) {
  f <- md_path(CFG, "cache", file)
  if (!file.exists(f)) {
    stop("上游产物不存在：", f, "\n请先运行模块 ", producer,
         "（Rscript --vanilla R/90_run_all.R --only=", producer, "）", call. = FALSE)
  }
  readRDS(f)
}

#' 13_distances 的产物：list(dist = 具名 dist 列表, meta = data.frame, params = list)
#' **13 必须写 cache/13_dist.rds**，结构 = list(dist=<具名 dist 列表>, meta=, params=)
md_dist_cache <- function(CFG) .md_cache_get(CFG, "13_dist.rds", "13")
#' 取某个距离度量的 dist 对象。metric 见 CFG$distances$metrics（bray/jaccard/unifrac/wunifrac）。
#' 取不到返回 NULL（而不是报错），方便调用方优雅降级。
md_dist_get <- function(CFG, metric = NULL) {
  dc <- tryCatch(md_dist_cache(CFG), error = function(e) NULL)
  if (is.null(dc)) return(NULL)
  metric <- metric %||% names(dc$dist)[1]
  if (is.null(metric) || !metric %in% names(dc$dist)) return(NULL)
  dc$dist[[metric]]
}

#' 30_null_models 的产物：list(bnti = 矩阵, rc = 矩阵, pair_df = 长表, process = 因子向量,
#'                              params = list(runs, thr, ...))
md_null_cache <- function(CFG) .md_cache_get(CFG, "30_null.rds", "30")

#' 21_phylo_diversity 的产物：data.frame(sample_id, PD, MPD, MNTD, NRI, NTI, ...)
md_phylo_cache <- function(CFG) .md_cache_get(CFG, "21_phylo.rds", "21")

#' 40_env_summary 的产物：list(env, vif, cor, pca, soil_pca, derived, ...)
#' derived 里含 Soil_PC1 等派生变量，供 50_lmm / 52_sem / 71 引用（见 CFG$model$derived）
md_env_cache <- function(CFG) .md_cache_get(CFG, "40_env.rds", "40")

#' 取派生变量向量（如 Soil_PC1），按给定 sample_id 顺序对齐；取不到返回 NULL
md_derived <- function(CFG, var = "Soil_PC1", ids = NULL) {
  ec <- tryCatch(md_env_cache(CFG), error = function(e) NULL)
  if (is.null(ec) || is.null(ec$derived) || !var %in% names(ec$derived)) return(NULL)
  d <- ec$derived
  if (is.null(ids)) return(d[[var]])
  d[[var]][match(ids, d[[1]])]
}

#' 把派生变量并入分析总表（50 / 52 / 71 用）
md_analysis_with_derived <- function(CFG) {
  d <- md_core_analysis(CFG)
  ec <- tryCatch(md_env_cache(CFG), error = function(e) NULL)
  if (!is.null(ec$derived) && ncol(ec$derived) > 1) {
    smp <- CFG$data$sample_col
    extra <- setdiff(names(ec$derived), names(d))
    if (length(extra)) {
      d <- merge(d, ec$derived[, c(names(ec$derived)[1], extra), drop = FALSE],
                 by.x = smp, by.y = names(ec$derived)[1], all.x = TRUE, sort = FALSE)
      d <- d[match(md_core_meta(CFG)[[smp]], d[[smp]]), , drop = FALSE]
    }
  }
  d
}
