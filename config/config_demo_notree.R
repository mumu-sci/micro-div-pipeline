# =============================================================================
# config/config_demo_notree.R —— 无系统发育树的降级验证配置
# -----------------------------------------------------------------------------
# 用途：只改「有没有树」这一个开关，验证全部系统发育相关模块能优雅降级，
#       而不是报错。ITS 用户的典型场景。
#
# 运行：
#   Rscript --vanilla R/90_run_all.R --config=config/config_demo_notree.R
#
# 预期：
#   13_distances  只算 bray / jaccard，跳过 unifrac / wunifrac 并打印原因
#   14_ordination 不画 wUniFrac 面板
#   20 / 21 / 22  SKIP，理由为「未提供系统发育树」
#   11_alpha      跳过 PD 指数
#   17_gamma_partition 不算系统发育分区
#   其余模块照常出图出表（42 / 50 / 52 / 60 / 70 / 71 等）
# =============================================================================

# 定位 config_demo.R：优先用 runner 传进来的 PKG_DIR，其次按工作目录/脚本位置推断
.cands <- c(
  if (exists("PKG_DIR", inherits = TRUE)) file.path(PKG_DIR, "config", "config_demo.R") else NULL,
  file.path(getwd(), "config", "config_demo.R"),
  "config/config_demo.R",
  file.path(tryCatch({
    fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)) else "."
  }, error = function(e) "."), "config_demo.R")
)
.p <- .cands[file.exists(.cands)][1]
if (is.na(.p) || is.null(.p))
  stop("找不到 config_demo.R，试过：\n  ", paste(.cands, collapse = "\n  "))
source(.p)

CFG$project$label <- "demo-36x300-notree"

# ---- 唯一的改动：没有树 ----
CFG$data$tree     <- NULL
CFG$phylo$enabled <- FALSE

# 距离口径里保留 unifrac/wunifrac（让 13 自己发现没树并跳过，才叫真正验证降级）
CFG$distances$metrics <- c("bray", "jaccard", "unifrac", "wunifrac")

# α 指数里保留 PD（让 11 自己发现没树并跳过）
CFG$alpha$indices <- c("Observed", "Chao1", "Shannon", "Simpson", "Pielou", "PD")

# 分区不做系统发育版
CFG$partition$engine <- "iNEXT"
