# =============================================================================
# config/config_demo_extras.R —— 示范数据的「全能力」配置
# -----------------------------------------------------------------------------
# 用途：在 config_demo.R 的基础上，把**往已有模块里加的增量**也全部打开，
#       用于验证这些可选路径确实能跑通，以及给用户看"全打开长什么样"。
#
# 与 config_demo.R 的唯一区别 = 下面这 5 个开关：
#   composition$indicator     10 号多一个 D 面板（指示种 IndVal）
#   env$rf                    40 号把预留的 D 区填上（随机森林重要性）
#   response_models$pattern   70 号改用 Scheirer-Ray-Hare 双因子非参数检验
#   response_models$linkage   71 号追加 Procrustes 与 RF 两个口径（C / D 面板）
#   taxon_group$indicator     81 号追加类群指示种
#
# ⚠ 这个配置**会改变已有模块的图与表**（多面板 / 多 sheet），
#   所以它只用于验证与演示，不作为默认交付口径。
# =============================================================================

# 复用 demo 配置（它会自己 source config_template.R）
.cfg_demo <- file.path(dirname(tryCatch({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE) else "."
}, error = function(e) ".")), "config_demo.R")
if (!file.exists(.cfg_demo)) .cfg_demo <- "config/config_demo.R"
source(.cfg_demo)

CFG$project$label <- "demo-36x300-extras"

# ---- 10 号：指示种 IndVal ----
CFG$composition$indicator$enabled <- TRUE
CFG$composition$indicator$rank    <- "ASV"
CFG$composition$indicator$func    <- "IndVal.g"
CFG$composition$indicator$nperm   <- 999
CFG$composition$indicator$top_n   <- 20

# ---- 40 号：随机森林重要性（填上设计串里预留的 D 区）----
CFG$env$rf$enabled    <- TRUE
# ⚠ 因变量必须是 `md_core_analysis()` 里有的列，或者是 α 指数
#   （Observed / Shannon / Simpson / Pielou —— 由模块 40 现算补齐）。
#   示范数据里 EF（生态系统功能）本来就在总表里，用它最直观。
CFG$env$rf$response   <- "EF"
CFG$env$rf$predictors <- NULL          # NULL = 用 env$vars
CFG$env$rf$engine     <- "ranger"
CFG$env$rf$num_trees  <- 2000
CFG$env$rf$top_n      <- 10
CFG$env$rf$pdp        <- TRUE          # 额外出 PDP 面板 E
CFG$env$rf$pdp_top    <- 3

# ---- 70 号：Scheirer-Ray-Hare 双因子非参数检验 ----
CFG$response_models$pattern_test <- "scheirer_ray_hare"

# ---- 71 号：追加 Procrustes 与 RF ----
CFG$response_models$linkage <- c("alpha", "ordination", "pairwise",
                                 "procrustes", "rf")
CFG$response_models$procrustes <- list(axes = 5, permutations = 999)
CFG$response_models$rf <- list(engine = "randomForest", num_trees = 3000,
                               importance = "permutation", pcoa_axes = 5,
                               rank = "Phylum", top_n = 15)

# ---- 81 号：类群指示种 ----
CFG$taxon_group$indicator$enabled <- TRUE
CFG$taxon_group$indicator$nperm   <- 999
CFG$taxon_group$indicator$top_n   <- 20
