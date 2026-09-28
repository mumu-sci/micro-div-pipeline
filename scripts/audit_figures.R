# =============================================================================
# scripts/audit_figures.R —— 图件体检（只读；改包回归与交付前自检用）
# -----------------------------------------------------------------------------
# 一次跑完四件事，输出一张 CSV + 控制台汇总：
#   1. **纯 ASCII 审计** —— 用 pdftools::pdf_text() 抽文字层，任何码点 > 0x7E 即计数
#      （硬约束 11：图中文字必须纯 ASCII；95_export_tables 会对交付 PDF 做同样的审计）
#   2. **位图 XObject 计数** —— 数 PDF 里的 /Subtype /Image。
#      正常情况下只有既有那两三张老图含位图；新图应当全是矢量。
#   3. **页尺寸** —— 毫米；对照 references/06-figure-map.md 里每张图的画幅
#   4. **文字层是否为空** —— 空的说明整幅被栅格化了（见 pitfalls C31/C32）
#
# 用法：
#   Rscript --vanilla scripts/audit_figures.R
#   Rscript --vanilla scripts/audit_figures.R --dir=figures --out=logs/figure_audit.csv
#   Rscript --vanilla scripts/audit_figures.R --pkg=D:/path/to/micro-div-pipeline
#
# 退出码：字符层出现非 ASCII，或某张图文字层为空 -> 1；否则 0。
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
DIR     <- .get_arg("dir", file.path(PKG_DIR, "figures"))
OUT     <- .get_arg("out", file.path(PKG_DIR, "logs", "figure_audit.csv"))

if (!requireNamespace("pdftools", quietly = TRUE))
  stop("需要 pdftools：install.packages('pdftools')", call. = FALSE)
if (!dir.exists(DIR)) stop("图件目录不存在: ", DIR, call. = FALSE)

pdfs <- sort(list.files(DIR, "[.]pdf$", full.names = TRUE))
if (!length(pdfs)) stop("目录里没有 pdf: ", DIR, call. = FALSE)

rows <- vector("list", length(pdfs))
for (i in seq_along(pdfs)) {
  p <- pdfs[i]
  txt <- tryCatch(pdftools::pdf_text(p), error = function(e) NULL)
  if (is.null(txt)) {
    rows[[i]] <- data.frame(file = basename(p), ok = FALSE, mm_w = NA_real_,
                            mm_h = NA_real_, non_ascii = NA_integer_,
                            chars = NA_integer_, images = NA_integer_,
                            note = "pdf_text() 读不了", stringsAsFactors = FALSE)
    next
  }
  joined <- paste(txt, collapse = "")
  ch <- unlist(strsplit(joined, ""))
  ch <- ch[nzchar(ch)]
  non_ascii <- if (length(ch)) sum(vapply(ch, function(z) utf8ToInt(z) > 0x7E, logical(1))) else 0L

  raw <- readBin(p, "raw", file.info(p)$size)
  n_img <- length(grepRaw(charToRaw("/Subtype /Image"), raw, all = TRUE, fixed = TRUE))

  ps <- tryCatch(pdftools::pdf_pagesize(p), error = function(e) NULL)
  mm_w <- if (!is.null(ps) && nrow(ps)) ps$width[1] / 72 * 25.4 else NA_real_
  mm_h <- if (!is.null(ps) && nrow(ps)) ps$height[1] / 72 * 25.4 else NA_real_

  note <- character(0)
  if (non_ascii > 0) note <- c(note, "含非 ASCII 字符")
  if (n_img > 0)     note <- c(note, sprintf("含 %d 个位图 XObject", n_img))
  if (!nchar(trimws(joined))) note <- c(note, "文字层为空（可能整幅被栅格化）")
  if (length(txt) > 1) note <- c(note, sprintf("%d 页", length(txt)))

  rows[[i]] <- data.frame(
    file = basename(p), ok = !length(note), mm_w = round(mm_w, 1), mm_h = round(mm_h, 1),
    non_ascii = non_ascii, chars = length(ch), images = n_img,
    note = paste(note, collapse = "; "), stringsAsFactors = FALSE)
}

res <- do.call(rbind, rows)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, OUT, row.names = FALSE, fileEncoding = "UTF-8")

cat(sprintf("\n图件 %d 张 | 目录 %s\n", nrow(res), DIR))
cat(sprintf("非 ASCII 字符总数 : %d\n", sum(res$non_ascii, na.rm = TRUE)))
cat(sprintf("文字层为空的图    : %d\n", sum(grepl("文字层为空", res$note))))
cat(sprintf("含位图的图        : %s\n",
            if (any(res$images > 0, na.rm = TRUE))
              paste(sprintf("%s(%d)", res$file[res$images > 0], res$images[res$images > 0]),
                    collapse = ", ") else "无（全部矢量）"))
bad <- res[!res$ok, , drop = FALSE]
if (nrow(bad)) {
  cat("\n需要注意的图：\n")
  print(bad[, c("file", "mm_w", "mm_h", "non_ascii", "images", "note")], row.names = FALSE)
} else {
  cat("\n全部图件通过（无非 ASCII、无空文字层）\n")
}
cat("\n报告: ", OUT, "\n", sep = "")
quit(status = if (any(res$non_ascii > 0, na.rm = TRUE) ||
                  any(grepl("文字层为空", res$note))) 1L else 0L, save = "no")
