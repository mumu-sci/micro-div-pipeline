# =============================================================================
# config/config_demo.R —— 示范数据配置（36 样本 x 300 ASV）
# -----------------------------------------------------------------------------
# 用法：Rscript --vanilla R/90_run_all.R --config=config/config_demo.R
# 这是「一份可直接跑的完整配置」，也是新论文配置的填空范例。
# =============================================================================

# 定位 config_template.R：优先用 runner 传进来的 PKG_DIR，其次按脚本/工作目录推断
.cfg_cands <- c(
  if (exists("PKG_DIR", inherits = TRUE)) file.path(PKG_DIR, "config", "config_template.R") else NULL,
  file.path(getwd(), "config", "config_template.R"),
  "config/config_template.R",
  file.path(tryCatch({
    fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)) else "."
  }, error = function(e) "."), "config_template.R")
)
.cfg_tpl <- .cfg_cands[file.exists(.cfg_cands)][1]
if (is.na(.cfg_tpl) || is.null(.cfg_tpl))
  stop("找不到 config_template.R，试过：\n  ", paste(.cfg_cands, collapse = "\n  "))
source(.cfg_tpl)

CFG$project$label <- "demo-36x300"
CFG$project$dir   <- md_pkg_root()   # 包根目录（config/ 的上一级）

# ---- 设计 ----
CFG$design$factors        <- c("Zone", "Vegetation")
CFG$design$ordered_levels <- list(
  Zone       = c("Coastal", "Nearshore", "Xisha"),
  Vegetation = c("Bare", "Herb", "Tree")
)
CFG$design$sampling    <- "independent"
CFG$design$site_var    <- "Site"
CFG$design$replicate_var <- NULL

# ---- 过滤 / 标准化 ----
CFG$filter$min_asv_reads    <- 10
# 【示范数据说明】这批 nifH 数据极端单例化（23,251/23,873 个 ASV 只出现 1 次），
# ASV 级出现率过滤会直接砍掉大半数据；论文原流程用的是 reads>=10 + 类群子集过滤，
# 没有做出现率过滤。因此示范配置里把出现率下限设为 1。
# 换成 16S/ITS 等常见数据时，这里通常应设为 2（可在 references/02-params.md 查推荐）。
CFG$filter$min_prev_samples <- 1
CFG$filter$min_sample_reads <- 1000
CFG$filter$taxon_subset     <- NULL

CFG$normalize$method <- "rarefy"
CFG$normalize$depth  <- NULL          # NULL = 最小正样本 reads

# ---- 响应变量（示范数据里是生态系统功能 EF，原 nifH 的 ARA）----
CFG$response <- list(
  enabled = TRUE, var = "EF", label = "Ecosystem function",
  unit = "pmol C2H4 g-1 h-1", trans = "none"
)

# ---- 系统发育 ----
CFG$phylo$enabled     <- TRUE
CFG$phylo$prev_min    <- 2
CFG$phylo$root_method <- "midpoint"

# ---- α ----
CFG$alpha$indices  <- c("Observed", "Chao1", "Shannon", "Simpson", "Pielou", "PD")
CFG$alpha$test     <- "kruskal_dunn"
CFG$alpha$p_adjust <- "BH"

# ---- β ----
CFG$distances$metrics <- c("bray", "jaccard", "wunifrac")
CFG$ordination$method <- "pcoa"
CFG$beta_tests$permutations <- 999

# ---- 分区 ----
CFG$partition$engine  <- "iNEXT"
CFG$partition$q       <- c(0, 1, 2)
CFG$partition$t_gamma <- NULL        # 自动判定

# ---- 组装零模型 ----
CFG$assembly$runs     <- 999
CFG$assembly$ncores   <- 4
CFG$assembly$chunk    <- NULL
CFG$assembly$subset_top <- NULL

# ---- 环境 ----
CFG$env$vars <- c("MAT", "MAP", "SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe")
CFG$env$groups <- list(
  Climate = c("MAT", "MAP"),
  Soil    = c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe")
)
CFG$env$pca_vars  <- c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe")
CFG$env$pca_log   <- FALSE
CFG$env$pca_show  <- "sites"
CFG$env$vpa_blocks <- list(
  Geo     = c("Latitude", "Longitude"),
  Soil    = c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe"),
  Climate = c("MAT", "MAP")
)

# ---- 建模 ----
CFG$model$random      <- "(1|Site)"
CFG$model$reml        <- FALSE
CFG$model$focus       <- c("Shannon", "Soil_PC1")
CFG$model$covariates  <- c("MAT", "Soil_PC1")
CFG$model$edf_nonlin  <- 1.1
CFG$model$grid_n      <- 200

# ---- SEM：声明式边表（4 层：气候 -> 土壤 -> 多样性 -> 功能）----
CFG$sem$edges <- data.frame(
  from  = c("MAT", "MAT", "Soil_PC1", "Shannon", "Soil_PC1", "Shannon"),
  to    = c("Soil_PC1", "Shannon", "Shannon", "EF", "EF", "EF"),
  layer = c("soil", "diversity", "diversity", "response", "response", "response"),
  model = c("lm", "lm", "lm", "lm", "lm", "lm"),
  stringsAsFactors = FALSE
)
CFG$sem$layers <- c("climate", "soil", "diversity", "response")
CFG$sem$labels <- c(MAT = "MAT", Soil_PC1 = "Soil PC1", Shannon = "Shannon",
                    EF = "Ecosystem function")
CFG$sem$dsep      <- TRUE
CFG$sem$p_adjust  <- "BH"

# ---- 敏感性分析（53）：示范数据里剔掉一个分区看结论稳不稳 ----
# 原 nifH 论文做过「剔南沙（人工岛）」；示范数据不含 Nansha，这里改成剔 Xisha（最远的岛群）。
CFG$sensitivity <- list(
  enabled      = TRUE,
  drop_levels  = list(Zone = "Xisha"),
  drop_sites   = character(0),
  drop_samples = character(0),
  target       = c("lmm"),
  max_subsets  = 6,
  min_n        = 12
)

# ---- 空间（示范数据 11 个站点，坐标站点级）----
CFG$spatial$level        <- "auto"
CFG$spatial$permutations <- 999
CFG$spatial$dist_classes <- c(0, 10, 100, 600, 1600)

# ---- 地图（示范数据有经纬度，但出图需联网取 DEM；默认关闭，见 DECISIONS.md）----
CFG$map$enabled <- FALSE
CFG$map$bbox    <- c(105, 4, 121, 24)

# ---- 响应变量关系 ----
CFG$response_models$linkage         <- c("alpha", "ordination", "pairwise")
CFG$response_models$alpha_focus     <- c("Shannon", "PD")
CFG$response_models$ordination_axes <- c("PCoA1")
CFG$response_models$pair_metric     <- "bray"

# =============================================================================
# 2026 毕业论文补充轮：7 个新模块（**示范配置全部打开**）
# =============================================================================
# 注意：往已有模块里"加面板"的四项增量（composition$indicator / env$rf /
# response_models$pattern_test = scheirer_ray_hare / linkage 加 procrustes+rf）
# 在这里**一律保持关闭**，以保证已有 19 张图与 135 张数据表零变化。
# 想看全部增量同时打开，用 config/config_demo_extras.R。

# ---- 12 生态位宽度与稀有性 ----
CFG$niche$enabled    <- TRUE
CFG$niche$index      <- "levins"
CFG$niche$level      <- "both"
CFG$niche$input      <- "raw"          # 未抽平（原脚本用 ps_rare）
CFG$niche$thr_abund  <- 0.1
CFG$niche$thr_prev   <- 10
# NULL = 前两个设计因子的交互组（Zone×Vegetation），与原图 Fig1c 的 x 轴一致
CFG$niche$groups     <- NULL

# ---- 44 约束排序（db-RDA）----
CFG$constrained_ord$enabled     <- TRUE
CFG$constrained_ord$method      <- "dbrda"
CFG$constrained_ord$response    <- "community"
CFG$constrained_ord$constraints <- c("MAT", "MAP", "SWC", "pH", "OM", "TN", "TP")
CFG$constrained_ord$transform   <- "hellinger"
CFG$constrained_ord$vif_thr     <- 15
CFG$constrained_ord$permutations <- 999

# ---- 45 GAM / GAMM 响应曲线（因变量不必是群落；预测变量 = 随机森林 top-N）----
CFG$response_curves$enabled    <- TRUE
CFG$response_curves$y          <- "PCoA1"
CFG$response_curves$x          <- c("MAT", "MAP", "SWC", "pH", "OM", "TN")
CFG$response_curves$x_from     <- "rf_top"   # 对齐原脚本：RF 排序后取前 3 个
CFG$response_curves$engine     <- "gam"
CFG$response_curves$k          <- 4
CFG$response_curves$random     <- NULL
CFG$response_curves$p_max      <- 0.1
CFG$response_curves$max_panels <- 3

# ---- 80 共现网络（会出 4 张图，是最重的一个模块）----
CFG$network$enabled        <- TRUE
CFG$network$method         <- "cor"     # 论文的"SparCC 近似"= 相关 + 阈值
CFG$network$cor_method     <- "spearman"
CFG$network$min_cor        <- 0.5
CFG$network$sub_min_cor    <- 0.5
CFG$network$min_prev       <- 0.02
CFG$network$max_asv        <- 1500      # 示范数据只有 327 个 ASV，不触发上限
CFG$network$module_algo    <- "louvain"
CFG$network$random_n       <- 100
CFG$network$robustness_nsim <- 30
CFG$network$n_boot         <- 1000
CFG$network$group_factors  <- c("Zone", "Vegetation")
CFG$network$sub_network_by <- "abundance_rank"
CFG$network$turnover_by    <- "Vegetation"
CFG$network$metrics_supp_figure <- TRUE

# ---- 81 通用类群分组（读 data/05_taxon_group.csv）----
CFG$taxon_group$enabled     <- TRUE
CFG$taxon_group$file        <- "data/05_taxon_group.csv"
CFG$taxon_group$match_by    <- "Genus"
CFG$taxon_group$group_label <- "Functional group"
CFG$taxon_group$drop_unmapped <- FALSE
CFG$taxon_group$indicator$enabled <- FALSE

# ---- 82 类群共享 UpSet ----
CFG$taxon_share$enabled    <- TRUE
CFG$taxon_share$factor     <- "Zone"    # 3 个水平 -> UpSet + 三圆 Venn
CFG$taxon_share$rank       <- "ASV"
CFG$taxon_share$show_venn  <- TRUE

# ---- 83 差异类群热图 ----
CFG$taxon_diff$enabled        <- TRUE
CFG$taxon_diff$rank           <- "Genus"
CFG$taxon_diff$clr            <- TRUE
CFG$taxon_diff$min_prev_frac  <- 0.05
CFG$taxon_diff$p_thr          <- 0.05
CFG$taxon_diff$top_n          <- 25
CFG$taxon_diff$heatmap_engine <- "ComplexHeatmap"

# ---- 往已有模块里加的增量：**全部保持关闭**（保住已有图与表）----
CFG$composition$indicator$enabled <- FALSE
CFG$env$rf$enabled                <- FALSE
CFG$response_models$pattern_test  <- "kruskal_dunn"

# ---- 图件 ----
CFG$figure$width_double  <- 183
CFG$figure$base_pt       <- 9
CFG$figure$cjk_font      <- NULL      # 示范包只出英文，避免无中文字体环境报错
CFG$figure$save_singles  <- FALSE
# 配色**不再自造** —— 用 12_plot_theme.R 里逐字照搬 nifH 的内置色板：
#   zone = #89C4D7 / #4DAFB5 / #2B6AA0 / #3B2C6B
#   veg  = #989898 / #959BB0 / #60859E / #B98079
# 需要覆盖时再显式写 CFG$figure$palette$Zone <- c(...)

# 默认不出图的模块 —— 只限「论文里没有对应图」的 5 个（见 references/06-figure-map.md）
#   01 数据体检 / 05 过滤与覆盖度 / 13 距离矩阵 / 15 多元检验（论文只有表）
#   22 系统发育 β（论文只有 Table S11）
# 其余模块各自出图。2026 补充轮之后，默认图件 = 19 + 10 = **29 张**：
#   12(1) + 44(1) + 45(1) + 80(4) + 81(1) + 82(1) + 83(1) = 10 张新图。
# 想把总数压到 26 张，把 CFG$network$metrics_supp_figure 设为 FALSE。
CFG$figure$skip_modules <- c("01", "05", "13", "15", "22")

# 箱线图样式：默认完全复刻 nifH 的 plot_dual_factor_panel（空心箱 + 第二因子上色 +
# 字母在 max*1.08 + KW 注释 + y 从 0 起 + 刻度朝内 + 不标样本量）
CFG$plot$box <- list(
  style          = "nifh",
  color_by       = "auto",   # 双因子 -> 第二因子着色（nifH 口径）
  fill           = NA,       # 空心箱体
  letters_at     = 1.08,
  letters_face   = "bold",
  letters_colour = "grey20",
  show_n         = FALSE,    # 不标样本量
  kw_annot       = TRUE,
  y_from_zero    = TRUE,
  y_expand       = 1.35,
  inward_ticks   = TRUE,
  inter_sep      = "_",
  x_labels       = "two_line"
)

# ---- 结果表 / 运行 ----
CFG$table$detail      <- FALSE
CFG$table$write_csv   <- FALSE
CFG$opts$seed         <- 123
CFG$opts$force_recalc <- FALSE
CFG$opts$ncores       <- 4
CFG$opts$strict       <- TRUE
