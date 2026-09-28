# =============================================================================
# scripts/snapshot_results.R —— 结果表逐 sheet 元素指纹（只读；改包回归用）
# -----------------------------------------------------------------------------
# 对 results/*.xlsx 的**每个 sheet** 算 (行数, 列数, 内容 digest)，写一张 CSV。
# 改包前后各跑一次，diff 两张 CSV 就知道"哪个模块的哪张表变了"。
#
# 为什么不用文件 MD5、也不用整本 xlsx 的 MD5：
#   openxlsx 写出的 zip 里带创建时间戳，同一个 data.frame 每次写出来的字节都不同。
#   所以必须按 **sheet 内容** 比对（见 references/08-maintenance.md）。
#
# 用法：
#   Rscript --vanilla scripts/snapshot_results.R --out=logs/results_fp_before.csv
#   # ……改包、重跑 --force ……
#   Rscript --vanilla scripts/snapshot_results.R --out=logs/results_fp_after.csv
#
# 比对（不需要额外脚本，一行就够）：
#   Rscript -e "a<-read.csv('logs/results_fp_before.csv');b<-read.csv('logs/results_fp_after.csv');\
#     m<-merge(a,b,by=c('file','sheet'),suffixes=c('.a','.b'));\
#     print(m[m$fp.a!=m$fp.b,c('file','sheet','nrow.a','nrow.b')])"
#
# ⚠ 有两类 sheet **每次跑都必然不同**，比对时请排除，否则会淹没真正的差异：
#     * 每个 xlsx 的 `_meta` sheet（版本 / seed / 运行时间）
#     * `00_index.xlsx` 的 `Index_figures`（文件字节数）/ `Run_summary`（各模块耗时）/ `Session`
#   剔除它们之后剩下的，才是"分析结果有没有被改动"。
#
# 用法补充：--dir= 指定 results 目录；--pkg= 指定包根。
# =============================================================================

.args <- commandArgs(TRUE)
.get_arg <- function(name, default = NULL) {
  hits <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hits)) return(sub(paste0("^--", name, "="), "", hits[1]))
  default
}

.self <- tryCatch({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE) else NA_character_
}, error = function(e) NA_character_)
PKG_DIR <- .get_arg("pkg", if (!is.na(.self)) dirname(dirname(.self)) else getwd())
DIR     <- .get_arg("dir", file.path(PKG_DIR, "results"))
OUT     <- .get_arg("out", file.path(PKG_DIR, "logs", "results_fp.csv"))

for (pk in c("openxlsx", "digest"))
  if (!requireNamespace(pk, quietly = TRUE))
    stop("需要 ", pk, "：install.packages('", pk, "')", call. = FALSE)
if (!dir.exists(DIR)) stop("结果目录不存在: ", DIR, call. = FALSE)

xlsx <- sort(list.files(DIR, "[.]xlsx$", full.names = TRUE))
if (!length(xlsx)) stop("目录里没有 xlsx: ", DIR, call. = FALSE)

rows <- list()
for (p in xlsx) {
  sheets <- tryCatch(openxlsx::getSheetNames(p), error = function(e) character(0))
  for (s in sheets) {
    d <- tryCatch(openxlsx::read.xlsx(p, sheet = s), error = function(e) NULL)
    rows[[length(rows) + 1L]] <- data.frame(
      file  = basename(p),
      sheet = s,
      nrow  = if (is.null(d)) -1L else nrow(d),
      ncol  = if (is.null(d)) -1L else ncol(d),
      fp    = digest::digest(d, algo = "xxhash64"),
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
# 排序保证两次快照可比（getSheetNames 的顺序是稳定的，但仍显式排一次）
res <- res[order(res$file, res$sheet), , drop = FALSE]
rownames(res) <- NULL

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, OUT, row.names = FALSE, fileEncoding = "UTF-8")

cat(sprintf("\n结果表 %d 个 | sheet %d 张 | 目录 %s\n", length(xlsx), nrow(res), DIR))
cat(sprintf("组合指纹（全部 sheet 的 digest 再 digest）: %s\n",
            digest::digest(res[, c("file", "sheet", "fp")], algo = "xxhash64")))
cat("报告: ", OUT, "\n\n", sep = "")
