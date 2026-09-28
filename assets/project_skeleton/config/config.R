# =============================================================================
# <论文项目>/config/config.R
# -----------------------------------------------------------------------------
# 从包的 config/config_template.R 复制本文件，然后**只改下面标了 ★ 的地方**。
# 其余参数保持模板默认；要改先查 references/02-params.md。
# =============================================================================

.tpl <- "D:/01-Academic/00-学术文件/analysis-packages/micro-div-pipeline/config/config_template.R"
if (!file.exists(.tpl)) stop("找不到配置模板：", .tpl, "\n请把 .tpl 改成你机器上 micro-div-pipeline 的实际位置")
source(.tpl)

# ---- ★ 项目 ----
CFG$project$label <- "my-paper-2026"
CFG$project$dir   <- md_pkg_root()   # runner 会用 --project 覆盖

# ---- ★ 输入（相对项目目录）----
CFG$data$asv      <- "data/01_asv_table.csv"
CFG$data$taxonomy <- "data/02_taxonomy.csv"
CFG$data$meta     <- "data/03_sample_meta.csv"
CFG$data$env      <- "data/04_env_response.csv"
CFG$data$tree     <- "data/99_tree.nwk"      # 无树写 NULL
# 可选第 6 文件：类群 -> 分组映射（模块 81 用；不做这个分析就写 NULL）
CFG$data$taxon_group <- "data/05_taxon_group.csv"

# ---- ★ 硬门 1：分组因子与水平顺序 ----
CFG$design$factors <- c("FactorA", "FactorB")
CFG$design$ordered_levels <- list(
  FactorA = c("L1", "L2"),
  FactorB = c("G1", "G2")
)
CFG$design$sampling   <- "independent"       # independent | fixed_site | fixed_plot
CFG$design$site_var   <- "Site"

# ---- ★ 硬门 2：标准化 ----
CFG$normalize$method <- "rarefy"             # rarefy | relative | css | none | subsample_avg
CFG$normalize$depth  <- NULL                 # NULL = 最小正样本 reads

# ---- ★ 硬门 3：响应变量 ----
CFG$response <- list(enabled = TRUE, var = "EF", label = "Ecosystem function",
                     unit = "", trans = "none")
# 没有响应变量就写：CFG$response$enabled <- FALSE

# ---- ★ 硬门 4：系统发育树 ----
CFG$phylo$enabled <- TRUE                    # 无树改 FALSE

# ---- ★ 环境变量 ----
CFG$env$vars   <- c("MAT", "MAP", "pH", "SWC")
CFG$env$groups <- list(Climate = c("MAT", "MAP"), Soil = c("pH", "SWC"))
CFG$env$pca_vars <- c("pH", "SWC")
CFG$env$vpa_blocks <- list(
  Geo     = c("Latitude", "Longitude"),
  Soil    = c("pH", "SWC"),
  Climate = c("MAT", "MAP")
)

# ---- ★ 硬门 5：过滤 ----
CFG$filter$min_asv_reads    <- 10
CFG$filter$min_prev_samples <- 2             # 高周转数据可能要设 1（先看 01_data_health.xlsx）
CFG$filter$min_sample_reads <- 1000
# 16S 建议加：剔除叶绿体与线粒体
# CFG$filter$taxon_subset <- list(rank = "Order", keep = NULL,
#                                 drop = c("Chloroplast", "Mitochondria"))

# ---- ★ 硬门 6：零模型 ----
CFG$assembly$runs <- 999                     # 试跑可设 99；写论文前必须确认 999
CFG$assembly$ncores <- 4

# ---- ★ 硬门 7：SEM 边表（要跑 52 才填）----
CFG$sem$edges <- data.frame(
  from  = c("MAT", "Shannon"),
  to    = c("Shannon", "EF"),
  layer = c("diversity", "response"),
  model = c("lm", "lm"),
  stringsAsFactors = FALSE
)
CFG$sem$layers <- c("climate", "soil", "diversity", "response")
CFG$sem$labels <- c(MAT = "MAT", Shannon = "Shannon", EF = "Ecosystem function")

# ---- ★ 建模焦点与协变量 ----
CFG$model$random     <- "(1|Site)"
CFG$model$focus      <- c("Shannon")
CFG$model$covariates <- c("MAT")

# ---- ★ 空间（有经纬度才用）----
CFG$spatial$level        <- "auto"
CFG$spatial$dist_classes <- c(0, 1, 10, 100, 1000)

# ---- ★ 地图（要出研究区图才开；需联网取 DEM）----
CFG$map$enabled <- FALSE
CFG$map$bbox    <- c(100, 20, 120, 40)

# ---- ★ 2026 补充轮新增的 7 个模块（不用就把 enabled 设 FALSE）----
# 12 生态位宽度（Levins' B）
CFG$niche$enabled <- TRUE
CFG$niche$groups  <- "FactorA"               # 样本级策略组成的横轴分组
# 44 约束排序（db-RDA / RDA / CCA / capscale）
CFG$constrained_ord$enabled     <- TRUE
CFG$constrained_ord$constraints <- c("MAT", "pH", "SWC")
# 45 GAM / GAMM 响应曲线（因变量不必是群落）
CFG$response_curves$enabled <- TRUE
CFG$response_curves$y       <- "Shannon"
CFG$response_curves$x       <- c("MAT", "pH", "SWC")
# 80 共现网络（会出 4 张图）
CFG$network$enabled     <- TRUE
CFG$network$min_prev    <- 2 / 36            # 出现率下限 = 最少在几个样本里出现
CFG$network$group_factors <- c("FactorA", "FactorB")
CFG$network$turnover_by <- "FactorB"
# 81 通用类群分组（用 data/05_taxon_group.csv；不做就删掉那个文件）
CFG$taxon_group$enabled  <- TRUE
CFG$taxon_group$file     <- "data/05_taxon_group.csv"
CFG$taxon_group$match_by <- "Genus"          # ASV | Phylum | Class | Order | Family | Genus | Species
# 82 类群共享 UpSet
CFG$taxon_share$enabled <- TRUE
CFG$taxon_share$factor  <- "FactorA"
# 83 差异类群热图（CLR + KW + Dunn + Cliff's delta）
CFG$taxon_diff$enabled <- TRUE
CFG$taxon_diff$rank    <- "Genus"

# ---- ★ 往已有模块里加的增量（**默认关**，打开会改变该模块的图与表）----
CFG$composition$indicator$enabled <- FALSE   # 10 号：指示种 IndVal 面板
CFG$env$rf$enabled                <- FALSE   # 40 号：随机森林重要性面板
CFG$env$rf$response               <- "Shannon"
CFG$response_models$pattern_test  <- "kruskal_dunn"   # 或 "scheirer_ray_hare"
CFG$response_models$linkage       <- c("alpha", "ordination", "pairwise")
# 可选追加："procrustes"（群落排序一致性）、"rf"（群落组成 -> 响应的重要变量）

# ---- ★ 图件 ----
CFG$figure$base_pt <- 9
CFG$figure$cjk_font <- NULL                  # 中文图设 "Microsoft YaHei"
CFG$figure$palette$FactorA <- c(L1 = "#2B7A78", L2 = "#D9A05B")
CFG$figure$palette$FactorB <- c(G1 = "#8FAE5B", G2 = "#3E7C4A")

# ---- ★ 运行 ----
CFG$opts$seed <- 123                         # 整篇论文固定一个值
CFG$opts$strict <- TRUE
