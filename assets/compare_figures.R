# =============================================================================
# assets/compare_figures.R —— 两轮跑图后的逐像素比对（验收工具）
# -----------------------------------------------------------------------------
# 为什么需要它：PDF 内嵌 CreationDate，重跑后**未改动的图 MD5 也会变**，
# 所以不能用 MD5 判断有没有误伤（nifH 论文踩过这个坑）。
# 正解：渲染成位图逐像素比对。
#
# 用法：先跑一轮，把 figures/ 复制到别处；再跑一轮（例如无树配置），然后：
#   Rscript --vanilla assets/compare_figures.R <基准目录> <对比目录> [输出csv]
#
# 只比对**两边都存在**的图；只在一侧存在的图会被单独列出来。
# =============================================================================

.libPaths(c("D:/R-4.5.2/library", .libPaths()))
Sys.setlocale("LC_CTYPE", "")
suppressPackageStartupMessages({ library(pdftools); library(magick) })

args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("用法: Rscript --vanilla assets/compare_figures.R <基准figures目录> <对比figures目录> [输出csv]")
}
dir_a <- args[1]; dir_b <- args[2]
out_csv <- if (length(args) >= 3) args[3] else "figure_compare.csv"

fa <- list.files(dir_a, "\\.pdf$", full.names = TRUE)
fb <- list.files(dir_b, "\\.pdf$", full.names = TRUE)
na <- basename(fa); nb <- basename(fb)
common <- intersect(na, nb)

cat("基准目录:", dir_a, "->", length(fa), "张\n")
cat("对比目录:", dir_b, "->", length(fb), "张\n")
cat("共同:", length(common), "张 | 只在基准:", length(setdiff(na, nb)),
    "| 只在对比:", length(setdiff(nb, na)), "\n\n")

render <- function(p, dpi = 150) {
  img <- pdftools::pdf_render_page(p, page = 1, dpi = dpi)
  magick::image_read(img)
}

res <- do.call(rbind, lapply(common, function(n) {
  pa <- file.path(dir_a, n); pb <- file.path(dir_b, n)
  r <- tryCatch({
    ia <- render(pa); ib <- render(pb)
    # 对齐到较小尺寸（页面尺寸变化时也能比）
    w <- min(magick::image_info(ia)$width,  magick::image_info(ib)$width)
    h <- min(magick::image_info(ia)$height, magick::image_info(ib)$height)
    ia <- magick::image_crop(ia, paste0(w, "x", h, "+0+0"))
    ib <- magick::image_crop(ib, paste0(w, "x", h, "+0+0"))
    da <- as.numeric(magick::image_data(ia, "gray"))
    db <- as.numeric(magick::image_data(ib, "gray"))
    diff <- abs(da - db)
    data.frame(file = n, width = w, height = h,
               mean_abs_diff = mean(diff) / 255,
               max_abs_diff  = max(diff) / 255,
               frac_gt5pct   = mean(diff > 12.75),
               status = if (max(diff) == 0) "IDENTICAL" else "DIFF",
               stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(file = n, width = NA, height = NA, mean_abs_diff = NA,
               max_abs_diff = NA, frac_gt5pct = NA,
               status = "RENDER_FAIL", stringsAsFactors = FALSE)
  })
  r
}))

if (is.null(res)) res <- data.frame(file = character(0), status = character(0))
print(res, row.names = FALSE)
utils::write.csv(res, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat("\n明细已写:", out_csv, "\n")
cat("判定标准：两轮都**不该变**的图，status 应为 IDENTICAL（max_abs_diff = 0）。\n")
