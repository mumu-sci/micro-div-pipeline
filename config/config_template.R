# =============================================================================
# config/config_template.R  —— 全参数模板
# -----------------------------------------------------------------------------
# 用法：复制本文件为 <论文项目>/config/config.R，只改本文件，不动 R/ 里的模块代码。
#
# 写法约定（**请勿破坏**）：
#   每一项都写成「值 + 注释里的【默认 / 可选项 / 推荐 / 理由】」。
#   开工前必须与用户逐项过一遍（见 SKILL.md 第 3–4 步）；
#   标注了 【不许代选】 的项，用户没明确拍板就不许执行。
#
# seed 全包唯一：CFG$opts$seed。任何模块不许出现字面量种子。
# =============================================================================

# 小工具：本包只用这一个运算符，避免为了一行代码引入整个 rlang
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

# 定位包根（config/ 的上一级）；供本文件与各 config_*.R 使用
md_pkg_root <- function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    p <- dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE))
    if (basename(p) == "config") return(dirname(p))
    return(p)
  }
  wd <- normalizePath(getwd(), mustWork = FALSE)
  if (dir.exists(file.path(wd, "config")) && dir.exists(file.path(wd, "R"))) return(wd)
  if (file.exists(file.path(wd, "config", "config_template.R"))) return(wd)
  if (basename(wd) == "config") return(dirname(wd))
  wd
}

CFG <- list(

  # ---------------------------------------------------------------------------
  # 项目
  # ---------------------------------------------------------------------------
  project = list(
    dir   = ".",              # 项目根目录（绝对路径或相对当前工作目录）
    label = "unnamed-paper"   # 只用于日志与图注前缀
  ),

  # ---------------------------------------------------------------------------
  # 路径（相对 project$dir）
  # ---------------------------------------------------------------------------
  paths = list(
    data    = "data",
    cache   = "cache",
    results = "results",
    figures = "figures",
    logs    = "logs",
    code    = "code"
  ),

  # ---------------------------------------------------------------------------
  # 输入契约：五个纯文本文件
  # ---------------------------------------------------------------------------
  data = list(
    asv        = "data/01_asv_table.csv",      # 首列 asv_id，其余列 = 样本丰度（原始 reads）
    taxonomy   = "data/02_taxonomy.csv",       # 首列 asv_id + Kingdom..Species（可少不可多）
    meta       = "data/03_sample_meta.csv",    # 首列 sample_id + 分组因子 + 协变量
    env        = "data/04_env_response.csv",   # 首列 sample_id + 环境变量 + 响应变量（可无）
    tree       = "data/99_tree.nwk",           # Newick；无树写 NULL（系统发育模块自动降级）
    # 可选第 6 文件：**类群 → 分组** 映射表（模块 81 用）。两列：taxon, group
    #   分组含义不限（功能群 / 生态位 / R-K 对策 / 任意自定义）；
    #   按哪个分类层级匹配由 taxon_group$match_by 声明。
    #   不提供（写 NULL 或文件不存在）时模块 81 自动跳过，**五文件契约不变**。
    taxon_group = "data/05_taxon_group.csv",
    id_col     = "asv_id",
    sample_col = "sample_id",
    # 【默认】七级；少一级就删一级，顺序不能乱
    tax_ranks  = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ),

  # ---------------------------------------------------------------------------
  # 实验设计
  # ---------------------------------------------------------------------------
  design = list(
    # 【不许代选】分组因子（1–3 个），必须是 meta 里真实存在的列
    factors = c("Zone", "Vegetation"),
    # 【不许代选】因子水平顺序（决定箱线图顺序与字母标注顺序）；NULL = 按字母序
    ordered_levels = list(
      Zone       = c("Coastal", "Nearshore", "Xisha"),
      Vegetation = c("Bare", "Herb", "Tree")
    ),
    # 【默认 independent】independent(各样点独立) | fixed_site(固定站点重复采样) | fixed_plot(固定样方跟踪)
    #  【推荐】有重复采样/跟踪设计时必须改成对应值，否则违反独立性假设
    sampling = "independent",
    # 站点列（随机效应、站点级聚合、伪重复控制用）；没有就写 NULL
    site_var = "Site",
    replicate_var = NULL      # 可选的重复编号列
  ),

  # ---------------------------------------------------------------------------
  # 过滤（在标准化之前）
  # ---------------------------------------------------------------------------
  filter = list(
    min_asv_reads     = 10,   # 【默认 10】ASV 总 reads 下限 | 【推荐】10，去掉测序噪声
    min_prev_samples  = 2,    # 【默认 2】ASV 至少在几个样本中出现 | 【推荐】2；稀释曲线/Chao1 想保留 singleton 时改 1
    min_sample_reads  = 1000, # 【默认 1000】样本 reads 下限 | 【推荐】按数据量级定，低于此值的样本会被剔除
    # 类群子集过滤（通用机制）：按某个分类层级的关键词保留/剔除
    #   例：只留固氮相关类群 -> list(rank="Class", keep=c("Alphaproteobacteria"), drop=NULL)
    #   例：剔除叶绿体/线粒体 -> list(rank="Order", keep=NULL, drop=c("Chloroplast","Mitochondria"))
    # 【默认 NULL = 不过滤】【推荐】16S 务必 drop 叶绿体与线粒体
    taxon_subset = NULL
  ),

  # ---------------------------------------------------------------------------
  # 标准化
  # ---------------------------------------------------------------------------
  normalize = list(
    # 【默认 rarefy】rarefy(抽平到统一深度) | relative(相对丰度) | css(累积和标度) | none(不处理) | subsample_avg(多次亚抽样取均值)
    #  【推荐】多样性指数与距离用 rarefy（与 nifH 论文口径一致、可复现）；
    #         组成堆叠图与稀释曲线用 none/relative（保留全部 reads，图形更平滑）
    method = "rarefy",
    # 【默认 NULL】NULL = 最小正样本 reads（与 nifH 论文一致）| 也可显式给整数
    #  【推荐】NULL；但必须在结果表脚注写明实际深度。ITS 因长度变异大，常需显式指定
    depth = NULL,
    subsample_reps = 20,      # method = "subsample_avg" 时的重复次数
    # 逐模块覆盖（key 见 references/03-module-specs.md）
    by_module = list(
      composition = "none",   # 门/属组成用未抽平
      rarefaction = "none",   # 稀释曲线必须用原始 reads
      alpha       = "rarefy",
      distances   = "rarefy",
      phylo       = "rarefy",
      assembly    = "rarefy",
      partition   = "rarefy",
      env         = "rarefy",
      model       = "rarefy"
    )
  ),

  # ---------------------------------------------------------------------------
  # 响应变量（生态系统功能等）。不填 enabled=FALSE 时，70/71 模块自动跳过
  # ---------------------------------------------------------------------------
  response = list(
    # 【不许代选】哪一列是响应变量；声明"本论文没有响应变量"就写 FALSE
    enabled = TRUE,
    var     = "EF",
    label   = "Ecosystem function",
    unit    = "",             # 写进坐标轴与表头
    # 【默认 none】none | log | sqrt | z
    trans   = "none"
  ),

  # ---------------------------------------------------------------------------
  # 系统发育
  # ---------------------------------------------------------------------------
  phylo = list(
    enabled     = TRUE,       # 【不许代选】无树写 FALSE，相关模块自动降级并在日志里说明
    prev_min    = 2,          # 系统发育 β 的最少出现样本数（betapart 用）
    root_method = "midpoint", # midpoint | outgroup | none
    # 系统发育信号要检验的性状/指数列（来自 meta 或 env）；NULL = 默认用 alpha 指数
    signal_traits = NULL
  ),

  # ---------------------------------------------------------------------------
  # α 多样性
  # ---------------------------------------------------------------------------
  alpha = list(
    # 【默认全部】Observed | Chao1 | ACE | Shannon | Simpson | Pielou | PD
    #  【推荐】Observed + Shannon + Pielou + PD（论文最常用；Chao1/ACE 对测序深度敏感）
    indices = c("Observed", "Chao1", "Shannon", "Simpson", "Pielou", "PD"),
    # 【默认 kruskal_dunn】anova_tukey | kruskal_dunn | srh_nemenyi | lmm | none
    #  【推荐】先看 01_data_health 的正态性诊断：正态且方差齐 -> anova_tukey；否则 kruskal_dunn
    test = "kruskal_dunn",
    p_adjust = "BH",
    # 【默认 auto】auto(>4 组自动字母，<=4 组可用括号) | letters | brackets | none
    sig_label = "auto",
    rounding  = 3,
    rarefaction = list(
      curve = "vegan",        # vegan(rarecurve) | inext(iNEXT 稀疏外推) | both
      step  = 200,
      ci    = TRUE,           # 是否画 95% 置信带（需样本量 >= 3）
      nboot = 50
    )
  ),

  # ---------------------------------------------------------------------------
  # β 多样性
  # ---------------------------------------------------------------------------
  distances = list(
    # 【默认 bray+jaccard+wunifrac】bray | jaccard | unifrac | wunifrac
    #  【推荐】Bray（丰度）+ Jaccard（有无）+ wUniFrac（系统发育，需树）
    metrics = c("bray", "jaccard", "wunifrac"),
    binary_jaccard = TRUE
  ),

  ordination = list(
    # 【默认 pcoa】pcoa | nmds | both
    #  【推荐】pcoa（与 PERMANOVA 同源、可报轴解释率）；NMDS 只作稳健性补充
    method = "pcoa",
    nmds_trymax = 200,
    ellipse = TRUE,
    ellipse_level = 0.95,
    ellipse_min_n = 4        # 每组少于这个数不画椭圆
  ),

  beta_tests = list(
    permanova_terms = TRUE,   # by = "terms"（注意：本机 vegan 的 by="margin" 只返末项）
    anosim          = TRUE,
    permdisp        = TRUE,   # PERMDISP / betadisper（检验离散度差异）
    mrpp            = FALSE,
    pairwise        = TRUE,
    permutations    = 999,    # 【默认 999】正式写论文用 999；快速试跑用 199
    p_adjust        = "fdr",
    min_group_n     = 3       # 每组少于这个数不参与两两比较
  ),

  # ---------------------------------------------------------------------------
  # 多样性分区（α / β / γ）
  # ---------------------------------------------------------------------------
  partition = list(
    # 【默认 iNEXT】iNEXT(分类学 Hill 数) | iNEXT.3D(可做系统发育 PD) | both
    engine = "iNEXT",
    # 【默认 0,1,2】Hill 阶；【推荐】论文报 q=0/1/2 三档
    q      = c(0, 1, 2),
    nboot  = 50,
    conf   = 0.95,
    t_alpha = 1,              # α 的稀释水平（采样单元数）
    # 【默认 NULL = 自动】γ 的稀释水平；NULL 时取「各组采样单元数的中位数」并对不足的组标注
    t_gamma = NULL,
    coverage_aware = TRUE     # 是否同时报样本覆盖度 SC 并在表注提示外推风险
  ),

  # ---------------------------------------------------------------------------
  # 群落组装零模型（βNTI / RCbray）
  # ---------------------------------------------------------------------------
  assembly = list(
    # 【不许代选】随机化次数：正式 999，快速试跑 99–199。**写进论文前必须确认是 999**
    runs       = 999,
    null_model = "independentswap",  # picante::randomizeMatrix 的零模型
    fixedmar   = "columns",          # vegan::permatfull 的固定边际
    bnti_thr   = 2,                  # |βNTI| > 2 判定确定性
    rcbray_thr = 0.95,               # |RCbray| > 0.95 判定扩散过程
    # 大样本提速（16S 数百样本时 RCbray 的 n^2 x runs 数组会爆内存）
    ncores     = 4,                  # 【默认 4】并行核数
    chunk      = NULL,               # NULL = 不分块；给整数则按行分块计算再合并
    subset_top = NULL                # NULL = 用全部 ASV；给整数则只用丰度前 N 个 ASV（**会改变结论，须在表注说明**）
  ),

  # ---------------------------------------------------------------------------
  # 环境变量
  # ---------------------------------------------------------------------------
  env = list(
    vars = c("MAT", "MAP", "SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe"),
    # 变量分组（用于 VPA 分块与热图分组）
    groups = list(
      Climate = c("MAT", "MAP"),
      Soil    = c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe")
    ),
    pca_vars = c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe"),
    # 【默认 FALSE】是否对 PCA 输入做 log10(x+1e-4)；【推荐】量纲差异大时 TRUE，
    #   pH 属于对数刻度，log 变换时应排除（见 pca_pH_exempt_log）
    pca_log = FALSE,
    pca_pH_exempt_log = TRUE,
    pca_show = "sites",       # 【默认 sites】sites | species | both  ← 【不许代选】样方点/变量点二选一
    vif_thr = 10,
    cor_method = "spearman",  # spearman | pearson
    cor_p_adjust = "BH",      # 只对上三角独立配对校正，再镜像
    # 变差分解分块；键名会出现在图例里
    vpa_blocks = list(
      Geo     = c("Latitude", "Longitude"),
      Soil    = c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe"),
      Climate = c("MAT", "MAP")
    ),
    vpa_transform = "hellinger",  # hellinger | none
    vpa_permutations = 999,
    vpa_p_adjust = "holm",
    # ---- 随机森林变量重要性（模块 40 的 D 面板；**默认关**）----
    # 通用化：response 可以是任意连续列，predictors 可以是任意解释变量集 ——
    # 不限于"环境 → 群落"。毕业论文 Fig3b-B 用的是 ranger + permutation importance。
    rf = list(
      enabled = FALSE,           # 【默认 FALSE】打开后 40 号多一个 D 面板 + RF_importance 表
      response = "Shannon",      # 任意连续列（α 指数 / 排序轴 / 环境量 / 响应变量）
      predictors = NULL,         # NULL = 用 env$vars；也可显式给一列向量
      engine = "ranger",         # ranger（推荐，论文口径）| randomForest | cforest
      num_trees = 2000,          # 【默认 2000】论文口径；importance 在 ~2000 树已稳定
      importance = "permutation",# permutation（推荐）| impurity
      top_n = 10,                # 条形图显示的变量数
      pdp = FALSE,               # 是否加偏依赖（PDP）面板（需 pdp 包）
      pdp_top = 3,
      style = list(label_r2 = TRUE, fill = "#3C5488")
    )
  ),
  # （env 块结束）

  # ---------------------------------------------------------------------------
  # 建模（LMM / GLMM / GAMM 与条件效应）
  # ---------------------------------------------------------------------------
  model = list(
    # 随机效应（lme4 语法）；无随机效应写 NULL 退化为普通回归
    random = "(1|Site)",
    reml   = FALSE,           # 【默认 FALSE】便于用似然比做模型比较
    p_method = "satterthwaite",  # satterthwaite | kenward_roger | wald
    p_adjust = "BH",
    robust_se = FALSE,        # 是否同时报聚类稳健标准误（站点级整簇）
    # GAMM 非线性判定：edf 大于该值才认为非线性、才画平滑曲线
    edf_nonlin = 1.1,
    grid_n     = 200,         # 条件效应预测网格点数
    standardize = TRUE,       # 是否把预测变量 z 化后建模（系数可横向比较）
    # 【不许代选】焦点变量（画条件效应曲线）与协变量集合
    focus       = c("Shannon"),
    covariates  = c("MAT", "Soil_PC1"),
    # 由环境变量派生的中间变量（模块 40 生成后可在此引用）
    derived = c("Soil_PC1")
  ),

  # ---------------------------------------------------------------------------
  # SEM：声明式边表（**结构由用户声明，skill 不自动探索**）
  # ---------------------------------------------------------------------------
  sem = list(
    # 【不许代选】每行一条路径。layer 决定路径图纵向分层；model 用 lm 或 lmer
    #   例：data.frame(from="MAT", to="Soil_PC1", layer="soil", model="lm")
    #       data.frame(from="Shannon", to="EF", layer="response", model="lmer")
    edges = data.frame(
      from  = character(0), to = character(0),
      layer = character(0), model = character(0),
      stringsAsFactors = FALSE
    ),
    # 每一层在路径图上的纵向位置（从上到下），层名要与 edges$layer 一致
    layers = c("climate", "soil", "diversity", "response"),
    # 变量的显示名（图上标签）；没写的用原名
    labels = NULL,
    dsep = TRUE,              # 是否做 d-separation / Fisher's C 检验
    p_adjust = "BH"           # 对同一模型内全部路径校正；写 NULL 则不校正
  ),

  # ---------------------------------------------------------------------------
  # 敏感性分析（模块 53）：剔除特定组/站点/样本后重跑同一模型并并列对比
  # 审稿人一定会问「结论是不是被某几个样本或某个特殊组驱动的」
  # ---------------------------------------------------------------------------
  sensitivity = list(
    enabled = TRUE,
    # 按因子水平剔除：【推荐】把有特殊成因的组列出来（人工岛 / 极端污染点 / 单一样地）
    #   例：list(Zone = "Nansha")  或  list(Habitat = c("Urban","Mine"))
    drop_levels = list(),
    # 按站点剔除（用 design$site_var 匹配）
    drop_sites = character(0),
    # 按样本剔除（直接写 sample_id）
    drop_samples = character(0),
    # 重跑哪类模型：lmm（混合模型系数）/ sem（整条路径）
    target = c("lmm"),
    max_subsets = 6,          # 最多展示几个子集
    min_n = 12                # 子集样本量低于此值不跑
  ),

  # ---------------------------------------------------------------------------
  # 空间
  # ---------------------------------------------------------------------------
  spatial = list(
    # 【默认 auto】auto(坐标有重复 -> 站点级；全唯一 -> 样本级) | sample | site
    #  【推荐】auto。站内距离恒为 0 时做样本级空间自相关是伪重复
    level = "auto",
    moran_k = NULL,           # NULL = max(3, min(n/4, 10)) 的 k 近邻反距离权重
    permutations = 999,
    dist_classes = c(0, 10, 100, 600, 1600),  # correlogram 距离分段（km）
    min_links = 3             # 分段内连接数少于该值不出点
  ),

  # ---------------------------------------------------------------------------
  # 研究区地图（需要联网取 DEM；无网时自动降级为无 DEM 散点图）
  # ---------------------------------------------------------------------------
  map = list(
    enabled    = FALSE,       # 【默认 FALSE】有经纬度且要出研究区图时打开
    bbox       = c(103, 2, 123, 26),      # xmin, ymin, xmax, ymax
    zoom_bbox  = c(108, 17, 112, 22),
    dem_zoom   = 6,
    dem_fact   = 4,
    dem_limits = c(0, 2000),
    highlight  = NULL,        # 主图里要单独框出的区域 bbox
    point_label = TRUE,
    label_abbr = NULL,        # 站点名 -> 图面缩写；NULL = 用原站点名
    countries  = NULL,        # 是否画世界 inset（NULL = 不画）。注意：原脚本 inset 里**没有**
                              #   国家高亮层，所以这里不做高亮，只作"画不画 inset"的开关
    offline    = FALSE        # TRUE = 只用缓存 DEM，不联网
  ),

  # ---------------------------------------------------------------------------
  # 响应变量关系（模块 70 / 71）
  # ---------------------------------------------------------------------------
  response_models = list(
    # 响应变量的因子格局检验方法：
    #   【默认 kruskal_dunn】单因子 KW + Dunn + CLD 字母
    #   scheirer_ray_hare = 双因子非参数交互检验（源：毕业论文 8_ARA_function.R §2）；
    #     需要 rcompanion；缺包时回退到 ART + 双因素 ANOVA。
    #   ⚠ 换口径会改变 70 号图上的 P 值与字母，必须事前与用户确认（硬门）。
    pattern_test = "kruskal_dunn",
    # 跑哪几类关联：
    #   alpha      响应 vs α 多样性（单回归 + 偏回归 + LMM/GAMM）
    #   ordination 响应 vs 排序轴的条件效应
    #   pairwise   成对级 |Δ响应| vs 群落距离（MRM/Mantel + bootstrap）
    #   procrustes 群落排序 vs 响应/环境排序的 Procrustes（源：8_ARA_function.R §D）
    #   rf         群落组成 -> 响应的随机森林重要性（源：8_ARA_function.R §D）
    linkage = c("alpha", "ordination", "pairwise"),
    alpha_focus = c("Shannon", "PD"),
    ordination_axes = c("PCoA1"),    # 用哪个排序轴做条件效应
    pair_metric = "bray",            # 成对级用哪个距离
    bootstrap_B = 1000,              # 效应量 bootstrap 重复数
    # Procrustes（linkage 含 "procrustes" 时生效）
    procrustes = list(axes = 5, permutations = 999),
    # 群落 -> 响应的随机森林（linkage 含 "rf" 时生效）
    #   源：毕业论文 8_ARA_function.R:1299-1360
    #   特征 = PCoA 前 k 轴 + 某一分类层级的相对丰度
    rf = list(engine = "randomForest", num_trees = 3000, importance = "permutation",
              pcoa_axes = 5, rank = "Phylum", top_n = 15)
  ),

  # ---------------------------------------------------------------------------
  # 物种组成（模块 10）
  # ---------------------------------------------------------------------------
  composition = list(
    rank        = "Phylum",   # 主面板用的分类层级（taxonomy 表的列名）
    # 【重要】额外层级：一个层级一个堆叠柱面板（**只出图、不进表**）。
    # 你论文的 Fig.1B（门组成）与 Fig.S3（属组成）是两张图，所以默认再出一个属层级。
    # 想只出单张就设 character(0)。
    extra_ranks = c("Genus"),
    group_by    = NULL,       # 组成图的分组变量；NULL = **前两个设计因子的交互组**
                              #   （如 Zone_Vegetation = 3×3，与箱线图一致；原脚本同口径）。
                              #   想只按一个因子上色分组就写 "Zone"
    top_n       = 12,         # 每个层级显示前 N 个类群，其余归 "Other"
    # 分层级的 Top-N 覆盖（原脚本门用 12、属用 20；数量会直接决定属的配色分配）
    top_n_by_rank = list(Phylum = 12, Genus = 20),
    abund_t     = 0.1,        # 四象限：平均相对丰度阈值（%）
    freq_t      = 10,         # 四象限：出现率阈值（%）
    cell_colors = NULL,       # 四象限配色；NULL = 内置 quadrant 板（= nifH cols_2x2）
    # ---- 指示种 IndVal（模块 10 的 D 面板；**默认关**）----
    # 源：毕业论文 S4_supplementary_taxa.R §1（indicspecies::multipatt）
    # 通用化：按 design$factors 里声明的因子逐列出一张 lollipop（Top-N）。
    indicator = list(
      enabled  = FALSE,       # 【默认 FALSE】开启后 10 号多一个 D 面板 + Indicator_species 表
      rank     = "ASV",       # 在哪一层做指示种：ASV（原脚本口径）| Genus | Family | ...
      factors  = NULL,        # NULL = 用 design$factors（原脚本是 Zone 与 Vegetation 各一张）
      func     = "IndVal.g",  # IndVal.g（推荐，原脚本）| IndVal | r.g | r
      nperm    = 999,         # 【默认 999】置换次数
      adjust   = "BH",        # p 值校正；原脚本报表用 BH 的 fdr 列
      top_n    = 20,          # 每张图显示前 N 个（原脚本 20）
      fdr_thr  = 0.05,        # 只画 fdr < 该值的指示种
      style    = list(label_size_pt = 4.3, text_size = 2, point_size = 1.8)
    )
  ),

  # ---------------------------------------------------------------------------
  # 生态位宽度与稀有性（模块 12）
  #   源：毕业论文 1_community_composition.R §7（原 6_niche_rarity.R）
  # ---------------------------------------------------------------------------
  niche = list(
    enabled = TRUE,
    # 【默认 levins】levins（Levins' B，原脚本口径）| shannon | both
    index   = "levins",
    # 【默认 both】asv（ASV 级四象限 + 宽度）| sample（样本级生态策略组成）| both
    level   = "both",
    # 【默认 raw】raw（未抽平，原脚本用 ps_rare）| rarefied
    input   = "raw",
    # 四象限阈值（%，与 composition 的 abund_t / freq_t 独立）
    thr_abund = 0.1,
    thr_prev  = 10,
    # 阈值敏感性扫描（原脚本 SCI-8）：出现率阈值 x 丰度阈值
    sensitivity_occ   = c(0.1, 0.3, 0.5, 0.7, 0.9),
    sensitivity_abund = c(0.0001, 0.001, 0.01),
    # 样本级策略组成的横轴分组：
    #   NULL（推荐）= **前两个设计因子的交互组**（毕业论文 Fig1c 原图就是 Zone×Vegetation），
    #                 顺序为「第一因子为主序、第二因子为次序」
    #   给字符串     = 强制单因子（例如 "Zone"）
    groups = NULL,
    style = list(
      point_size = 0.4, point_alpha = 0.5,
      # A 面板（交互堆叠柱）的黑边与透明度 —— 逐字照抄 1_community_composition.R:701
      inter_alpha = 0.85, inter_lw = 0.3,
      # C 面板的 gam 平滑 —— 逐字照抄 1_community_composition.R:873-874
      gam = list(method = "gam", bs = "cs", colour = "#D55E00", lw = 0.6, alpha = 0.15)
    )
  ),

  # ---------------------------------------------------------------------------
  # 约束排序 db-RDA / RDA / CCA / capscale（模块 44）
  #   源：毕业论文 3_environmental_drivers.R:855-928（Fig3b-A 的 db-RDA triplot）
  # ---------------------------------------------------------------------------
  constrained_ord = list(
    enabled = TRUE,
    # 【默认 dbrda】dbrda（距离型，论文口径）| rda | cca | capscale
    method  = "dbrda",
    # 响应侧：community（用 ASV/类群丰度矩阵）| distance（用某一距离矩阵）
    #   【默认 community】换论文时因变量不一定是群落，这里也可以指向任意样本×变量矩阵
    response = "community",
    response_vars = NULL,     # response = "community" 时用哪些列；NULL = 用全部 ASV
    response_distance = NULL, # response = "distance" 时用哪个距离：bray | jaccard | wunifrac | unifrac
    # 解释变量：环境变量列名向量；NULL = 用 env$vars
    constraints = NULL,
    transform = "hellinger",  # hellinger | none（只对 response = "community" 生效）
    vif_thr   = 15,           # VIF 逐步剔除阈值（原脚本 15，不是 10）
    permutations = 999,
    # 箭头长度：原脚本按"最长箭头 = 站点云半径的 0.78"自动归一，这里保持同口径
    arrow_ratio = 0.78,
    style = list(arrow_col = "grey25", point_size = 1.4, label_size_pt = 6)
  ),

  # ---------------------------------------------------------------------------
  # GAM / GAMM 响应曲线（模块 45）
  #   源：毕业论文 3_environmental_drivers.R:1261-1345（Fig3b-C 的 GAM 响应）
  #   通用化：**因变量不必是群落** —— 任意连续列都行
  # ---------------------------------------------------------------------------
  response_curves = list(
    enabled = TRUE,
    # 因变量：任意连续列（排序轴 / α 指数 / 环境量 / 响应变量）
    y = "PCoA1",
    # 自变量：任意预测变量列名向量；NULL = 用 env$vars
    x = NULL,
    # 【默认 gam】gam（论文口径）| gamm（有随机效应时）| lm（线性对照）
    engine = "gam",
    k = 4,                    # 平滑基维数（原脚本 k = 4）
    random = NULL,            # 例如 "(1|Site)"；给出且 engine = "gam" 时自动升 GAMM
    grid_n = 100,             # 预测网格点数（原脚本 100）
    edf_nonlin = 1.1,         # edf 大于该值才算非线性
    # 【默认 0.1】只画 P 小于该值的平滑项（原脚本 3_environmental_drivers.R:1320 口径）；
    #   一个都没通过时**按原脚本不画曲线**，图里只给一句说明 + 表里列出全部 P。
    p_max = 0.1,
    # 【默认 3】最多画几个变量 —— 原脚本用随机森林的前 3 个
    max_panels = 3,
    # 预测变量从哪来：
    #   "rf_top"（**默认，推荐**）= 先用随机森林给 env$vars 排序、取前 max_panels 个
    #        （完全对应原脚本 `gam_preds <- names(sort(rf_A_pcoa1$variable.importance, decreasing=TRUE))[1:3]`）；
    #   NULL = 直接用 response_curves$x
    x_from = "rf_top",
    # x_from = "rf_top" 但读不到模块 40 / 71 的 RF 结果时，模块 45 自己跑一次用的参数
    rf = list(engine = "ranger", num_trees = 2000),
    style = list(line_col = "#0072B2", ribbon_fill = "grey80", point_alpha = 0.3)
  ),

  # ---------------------------------------------------------------------------
  # 共现网络（模块 80）—— 源：毕业论文 5_network_structure.R（Fig5a-e + FigS7）
  #   ⚠ 本模块会出 **4 张图**（节点属性 / 指标主图 / 指标附图 / 拓扑合图）
  # ---------------------------------------------------------------------------
  network = list(
    enabled = TRUE,
    # 建网口径：
    #   【默认 cor】相关系数矩阵 + 阈值（论文的"SparCC 近似"实际就是这一套）
    #   sparcc  = 真 SparCC（需 SpiecEasi，GitHub 安装；显著更慢）
    method = "cor",
    cor_method = "spearman",  # spearman（推荐，原脚本）| pearson
    min_cor = 0.5,            # 【默认 0.5】全局网络的 |r| 阈值（原脚本主网络用 0.5）
    sub_min_cor = 0.5,        # 丰富/稀有亚网络的 |r| 阈值（原脚本 0.5）
    min_prev = 0.02,          # 【默认 0.02】ASV 最少在多少比例的样本中出现
    # 节点最多保留多少个（按平均丰度降序）；NULL = 不限制。
    #   ⚠ 参与计算的 ASV 越多，相关矩阵 O(n^2) 与置换越慢；1000 个 ASV 约 ~30 s。
    max_asv = 1500,
    module_algo = "louvain",  # louvain（推荐，原脚本）| fast_greedy | walktrap
    zi_thr = 2.5,             # Zi 阈值（Deng et al. 2012）
    pi_thr = 0.62,            # Pi 阈值
    random_n = 100,           # 随机网络（Erdos-Renyi）模拟次数（原脚本 100）
    robustness_nsim = 30,     # 网络韧性模拟次数（原脚本 30）
    n_boot = 1000,            # 分组网络指标 bootstrap 均衡次数（原脚本 1000）
    permutations = 999,       # 样本级网络指标的 PERMANOVA 置换次数（原脚本 999）
    # 主图（原 Fig5b）与附图（原 FigS7）各 6 个指标；值 = 样本级指标列名
    #   逐字照抄毕业论文 5_network_structure.R:1443-1448 的 metric_order_main / _supp
    metrics_main = c("N_DetectedASV", "MeanNodeDegree", "MeanEigen",
                     "SubAvgPath", "SubNModules", "Robustness"),
    metrics_supp = c("SumNodeDegree", "SubDensity", "SubClustering",
                     "SubConnectedness", "SubModularity", "Cohesion_pos"),
    # 【默认 TRUE】附图（原 FigS7）是否出图；设 FALSE 可把总图数从 29 压到 26
    metrics_supp_figure = TRUE,
    # 分组网络（原 Fig5c）：用 design$factors 里的哪几个因子；NULL = 前 2 个
    group_factors = NULL,
    # 丰富/稀有亚网络（原 Fig5e）：用哪一列作为"细胞类型"分组
    #   abundance_rank = 按 mean_abund × occurrence 的 2x2 四象限（原脚本口径）
    sub_network_by = "abundance_rank",
    # 节点周转 Sankey（原 Fig5d-C）按哪个因子的相邻水平对比较；NULL = 第二个因子
    turnover_by = NULL,
    style = list(
      layout = "graphopt",    # 全局网络布局：graphopt（推荐，原脚本）| fr | nicely | circle
      layout_group = "fr",    # 分组网络布局：fr（原脚本，niter = 5000）| graphopt | nicely
      # 节点属性图（80d）版式：FALSE（默认）= 183 mm 宽三格纵向，与全包其余图同字号观感；
      #                          TRUE = 复刻原文 Fig5a 的 533.4 x 177.8 mm 单排
      wide_layout = FALSE,
      node_size = 1.6, node_size_hub = 2.6, edge_width = 0.18, edge_alpha = 0.45,
      node_alpha = 0.75, label_size_pt = 5, top_betweenness = 20
    )
  ),

  # ---------------------------------------------------------------------------
  # 通用类群分组（模块 81）—— 源：毕业论文 6_ecological_guilds.R
  #   分组含义不限：功能群 / 生态位 / R-K 对策 / 任意自定义。
  #   映射表 = data/05_taxon_group.csv（两列 taxon, group）；不提供则本模块跳过。
  # ---------------------------------------------------------------------------
  taxon_group = list(
    enabled = TRUE,
    file = "data/05_taxon_group.csv",  # NULL = 跳过；也可写别的路径
    taxon_col = "taxon",               # 映射表里的标识列
    group_col = "group",               # 映射表里的分组列
    # 按哪个层级匹配：ASV（用 asv_id 直接匹配）| Phylum | Class | Order | Family | Genus | Species
    match_by = "Genus",
    group_order = NULL,        # 分组显示顺序；NULL = 按映射表出现顺序
    group_label = "Group",     # 图例名（ASCII）
    drop_unmapped = FALSE,     # TRUE = 未匹配到的类群从组成图里剔除；FALSE = 归入 "Unmapped"
    # 类群指示种（IndVal）—— 与 composition$indicator 共用实现
    indicator = list(enabled = FALSE, func = "IndVal.g", nperm = 999, top_n = 20),
    style = list(palette = NULL, label_size_pt = 6, panel_heights = c(1, 1.2))
  ),

  # ---------------------------------------------------------------------------
  # 类群共享 UpSet（模块 82）—— 源：毕业论文 S4_supplementary_taxa.R §2
  # ---------------------------------------------------------------------------
  taxon_share = list(
    enabled = TRUE,
    factor = NULL,             # 按哪个因子分组交集；NULL = design$factors 的第一个
    rank = "ASV",              # ASV（原脚本口径）| Genus | Family | ...
    min_prev = 1,              # 类群至少在几个样本中出现才纳入（原脚本按有无划分）
    max_comb = 20,             # UpSet 最多显示多少个交集组合
    upset_top_n = 20,          # 顶部柱状图显示的交集个数
    show_venn = TRUE,          # 2/3 水平时是否附一张 Venn
    style = list(bar_color = "grey30", set_colors = NULL,
                 point_size = 1.6, line_width = 0.6)
  ),

  # ---------------------------------------------------------------------------
  # 差异类群热图（模块 83）—— 源：毕业论文 S4_supplementary_taxa.R §3
  # ---------------------------------------------------------------------------
  taxon_diff = list(
    enabled = TRUE,
    rank = "Genus",            # 在哪一层聚合：Genus（原脚本）| Family | Order | ASV
    clr = TRUE,                # 是否做 CLR 变换（原脚本口径）
    pseudocount = NULL,        # NULL = 最小正值的 1/2（原脚本口径）
    min_prev_frac = 0.05,      # 类群至少在多少比例的样本中出现（原脚本 0.05）
    min_n = 3,                 # 每组最少样本数（原脚本 3）
    test = "kruskal",          # kruskal（原脚本）| wilcox
    posthoc = "dunn",          # dunn（原脚本）| none
    adjust = "BH",
    effect = "cliff",          # cliff（原脚本）| none
    p_thr = 0.05,              # KW_padj 显著性阈值
    # 【默认 3】显著类群少于这个数时，热图回退为"按 KW_padj 排序的 Top N"，
    #   并在面板副标题里写明 "None/only N passed FDR ...; showing top N by KW p"。
    #   宁可如实标注，也不出一张 1 行、看不出任何模式的空壳热图。
    min_sig = 3,
    top_n = 25,                # 热图显示多少个差异类群（也是回退时的 Top-N）
    heatmap_engine = "ComplexHeatmap",  # ComplexHeatmap（推荐，原脚本）| pheatmap
    heatmap_scale = "row_z",   # row_z（推荐）| clr | rel
    show_box = FALSE,          # 是否在热图下补 Top 差异类群的箱线图
    style = list(row_label_pt = 5, col_label_pt = 7, heatmap_col = NULL)
  ),

  # ---------------------------------------------------------------------------
  # 图件
  # ---------------------------------------------------------------------------
  figure = list(
    width_single = 90,        # mm；单栏
    width_double = 183,       # mm；双栏（本包默认全用双栏）
    height_default = 110,     # mm；模块没指定高度时用
    base_pt = 9,              # 主图基准字号（pt）
    comp_ratio = 0.8,         # （已弃用）拼图字号现由 12_plot_theme.R 的 MD_PT_COMP 显式给出
    font     = "Arial",
    cjk_font = NULL,          # 【⚠ 只影响非图文字体】图中文字必须纯 ASCII，见 pitfalls B12
    save_singles = FALSE,     # 【默认 FALSE】是否额外落盘每个子面板
    # 【默认关图】只放「论文里没有对应图」的模块；其余模块各自出一张图。
    #  01 数据体检诊断图 / 05 过滤瀑布+覆盖度 / 13 距离矩阵 / 15 多元检验效应图
    #  22 系统发育 β（论文只有 Table S11）
    #  想让某个模块出图，就把它从这个向量里删掉。
    #  注：原 42_constrained_ord（RDA/CCA/db-RDA）已移出管线（本包自行扩展，论文未做）。
    skip_modules = c("01", "05", "13", "15", "22"),
    # 模块号 -> 论文图号（只影响图面标签与文件命名前缀）
    id_map = NULL,
    palette = list(
      # 只在**要覆盖内置 nifH 色板**时才填。键名 = 因子列名。
      # 内置色板（12_plot_theme.R 的 MD_BUILTIN_PAL）：
      #   zone     = Coastal #89C4D7 / Nearshore #4DAFB5 / Xisha #2B6AA0 / Nansha #3B2C6B
      #   veg      = Bare #989898 / Herb #959BB0 / Shrub #60859E / Tree #B98079
      #   assembly = 5 个 Stegen 过程
      #   quadrant = Core #D73027 / Permanently_rare #999999 /
      #              Conditionally_rare #E69F00 / Transient_rare #4575B4
      # 留 NULL 就用内置板
      Zone       = NULL,
      Vegetation = NULL,
      group      = NULL,
      # 地形高程渐变（地图用）；NULL = 内置 value，= nifH 07_S1_study_area_map.R §4 的两色
      dem        = NULL,
      diverging  = c("#2166AC", "#F5F5F5", "#B2182B"),
      sequential = c("#F7FBFF", "#4292C6", "#08306B"),
      qualitative = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                      "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")
    )
  ),

  # ---------------------------------------------------------------------------
  # 箱线图样式（默认完全复刻 nifH 的 plot_dual_factor_panel；换论文时逐项改）
  # ---------------------------------------------------------------------------
  plot = list(
    box = list(
      style          = "nifh",   # nifh | simple
      # 是否画 CLD 显著性字母（**默认 TRUE = 与原代码构建方法一致**）。
      # 以后写论文不想要字母就设 FALSE —— 只影响出图，字母本身仍写在 Descriptive 表里。
      letters        = TRUE,
      # 【默认 auto】auto(双因子用第二因子，单因子用第一因子) | f1 | f2 | none
      #   nifH 的做法是"用第二因子给箱体边框与抖动点上色"；
      #   如果你的设计不是双因子、或想统一上色，改成 "f1" 或 "none"
      color_by       = "auto",
      fill           = NA,       # NA = 空心箱体（nifH 口径）；给颜色则填充
      letters_at     = 1.08,     # 字母 y = 每组 max(y) * 该系数
      letters_face   = "bold",
      letters_colour = "grey20",
      show_n         = FALSE,    # 【默认 FALSE】不标样本量（nifH 没有这个标注）
      kw_annot       = TRUE,     # 右上角 KW 注释 "Zone: p = 0.033 *; Vegetation: p = 0.486 ns"
      y_from_zero    = TRUE,     # y 轴从 0 起
      y_expand       = 1.35,     # y 上限 = max * 该系数
      inward_ticks   = TRUE,     # 刻度朝内（nifH 用 unit(-0.12,"cm")）
      inter_sep      = "_",      # 双因子组合分隔符，**必须与生成 letters 的一致**
      x_labels       = "two_line"  # two_line(Zone\nVegetation) | one_line
    )
  ),

  # ---------------------------------------------------------------------------
  # 备用件（**默认全部关闭**，不影响论文的图与数值）
  # ---------------------------------------------------------------------------
  # 用途：把「论文里不常用、但换数据时可能要」的能力收在这里，别混进主流程。
  # 详见 references/07-spare.md。
  spare = list(
    # 备用面板 1 —— 任意**样本级指标** vs **群落周转**（原论文 Fig.3D 的通用化版本）
    #   放在模块 31_assembly_process 里，enabled=FALSE 时不出该面板
    process_vs_turnover = list(
      enabled = FALSE,
      x       = NULL,        # 横轴：任意样本级连续列（过程占比 / α 指数 / 环境变量）
      y       = "MeanBeta",  # 纵轴标签；实际值 = 该样本与其余样本的平均群落距离
      group   = NULL,        # 可选：分组着色
      covar   = NULL,        # 可选：加进 lm/gam 的协变量
      fit     = "lm"         # lm | gam
    ),
    # 备用面板 2 —— 一组**焦点变量 → 响应变量**的标准化系数森林图
    #   （原论文 Fig.S10 右半的通用化版本）；放在模块 50_lmm 里
    std_beta = list(
      enabled    = FALSE,
      response   = NULL,         # NULL = 用 CFG$response$var
      focal      = character(0), # 每个焦点变量各拟合一个模型
      covariates = character(0),
      random     = NULL,         # NULL = 用 CFG$model$random
      by         = NULL          # 可选：按该因子分组各拟合一套，并列呈现
    )
  ),

  # ---------------------------------------------------------------------------
  # 结果表
  # ---------------------------------------------------------------------------
  table = list(
    detail = FALSE,           # 【默认 FALSE】是否额外输出逐样本/逐对明细表（大样本时会出巨表）
    digits = 3,
    write_csv = FALSE         # 是否同时落一份 csv（便于 Python/Excel 读）
  ),

  # ---------------------------------------------------------------------------
  # 运行选项
  # ---------------------------------------------------------------------------
  opts = list(
    seed          = 123,      # ★ 全包唯一种子；所有随机过程都从这里派生
    force_recalc  = FALSE,    # TRUE = 忽略缓存重算全部
    figures_only  = FALSE,    # TRUE = 只重画图（走缓存），跳过全部重计算
    ncores        = 4,
    strict        = TRUE,     # TRUE = 数据体检不通过就停；FALSE = 只警告继续
    # 【默认 TRUE】缺包时**自动** install.packages 后再试；装不上才让该模块 SKIP。
    #   想完全离线/只报告不安装，就设 FALSE（此时行为退回"打印安装命令并跳过"）。
    auto_install  = TRUE,
    repos         = "https://cloud.r-project.org"
  )
)

# 允许用环境变量覆盖（例如临时快速试跑）
if (nzchar(Sys.getenv("MD_SEED")))    CFG$opts$seed <- as.integer(Sys.getenv("MD_SEED"))
if (nzchar(Sys.getenv("MD_FORCE")))   CFG$opts$force_recalc <- toupper(Sys.getenv("MD_FORCE")) %in% c("1","TRUE","YES")
if (nzchar(Sys.getenv("MD_FIGONLY"))) CFG$opts$figures_only <- toupper(Sys.getenv("MD_FIGONLY")) %in% c("1","TRUE","YES")
if (nzchar(Sys.getenv("MD_RUNS")))    CFG$assembly$runs     <- as.integer(Sys.getenv("MD_RUNS"))
