# 02 · 全参数总表（开工前逐项与用户确认）

> **用法**：按本表逐项问用户。每一项都给**推荐值 + 一句科学理由**。
> 标 **【不许代选】** 的七项，用户没明确拍板就**不许执行**；用户回「都行/你定/随便」不算授权。
> 「键」列是 `config/config.R` 里的路径；改配置只改这些，不动 `R/` 里的代码。

---

## A. 项目与路径

| 键 | 默认 | 可选项 | 推荐 | 理由 / 何时要改 |
|---|---|---|---|---|
| `project$label` | `unnamed-paper` | 任意短串 | 用论文短名 | 只用于日志与产物 `_meta`；换论文必须改 |
| `project$dir` | `.` | 绝对路径 | 由 runner 解析 | 跨项目跑时用 `--project=` 指定 |
| `paths$*` | `data/cache/results/figures/logs/code` | 任意 | 保持默认 | 目录名变化会影响所有下游路径 |

---

## B. 输入契约

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `data$asv` | `data/01_asv_table.csv` | 保持 | 见 `references/00-data-contract.md` |
| `data$taxonomy` | `data/02_taxonomy.csv` | 保持 | 无 taxonomy 时组成模块跳过 |
| `data$meta` | `data/03_sample_meta.csv` | 保持 | 必需 |
| `data$env` | `data/04_env_response.csv` | 保持 | 没有环境/响应变量时设 NULL |
| `data$tree` | `data/99_tree.nwk` | 有树就用 | 无树时设 **NULL**，20/21/22 与 wUniFrac 自动跳过 |
| `data$id_col` / `sample_col` | `asv_id` / `sample_id` | 保持 | 只有首列名不同才需要改 |
| `data$tax_ranks` | 7 级 | 按实际层级删减 | 顺序不能乱；少一级就删一级 |

---

## C. 实验设计

| 键 | 默认 | 可选项 | 推荐 | 理由 |
|---|---|---|---|---|
| `design$factors` | `c("Zone","Vegetation")` | 1–3 个 meta 实有列 | **【不许代选】** 按科学问题定 | 决定所有分组分析与图的分面 |
| `design$ordered_levels` | list | 每个因子的水平顺序 | **【不许代选】** 按梯度排（近→远 / 低→高） | 箱线图顺序与字母顺序；顺序错会让"趋势"读反 |
| `design$sampling` | `independent` | `independent` / `fixed_site` / `fixed_plot` | 有重复采样就改 | 固定站点/样方跟踪时，样点不独立，必须用随机效应吸收 |
| `design$site_var` | `Site` | meta 里的站点列；没有设 NULL | 尽量提供 | 随机效应、站点级聚合、伪重复控制都靠它 |
| `design$replicate_var` | NULL | meta 里的重复编号列 | 有就填 | 只用于描述性汇总 |

---

## D. 过滤（标准化之前）

| 键 | 默认 | 可选项 | 推荐 | 理由 |
|---|---|---|---|---|
| `filter$min_asv_reads` | 10 | 1–50 | **10** | 去掉测序噪声与嵌合体；低于 10 的 ASV 置信度低 |
| `filter$min_prev_samples` | 2 | 1–10 | **看数据**：常见 16S/ITS 用 **2–5**；极端单例化数据（如本包示范的 nifH）用 **1** | 出现率过滤去噪效果好，但对高周转数据会砍掉大半信息；先看 `01_data_health.xlsx` 的 ASV 出现率分布再定 |
| `filter$min_sample_reads` | 1000 | 依赖数据量级 | 取数据中最小可信深度 | 读长过低的样本抽平后会丢大量 ASV |
| `filter$taxon_subset` | NULL | `list(rank=, keep=, drop=)` | **16S 务必 drop 叶绿体与线粒体** | 见下方示例 |

**`taxon_subset` 的通用用法**（这就是把 原实现"Cluster IV 过滤"通用化后的机制）：
```r
# 16S：剔除叶绿体与线粒体
filter$taxon_subset <- list(rank = "Order", keep = NULL, drop = c("Chloroplast", "Mitochondria"))
# nifH：只保留固氮相关类群
filter$taxon_subset <- list(rank = "Class", keep = c("Alphaproteobacteria","Betaproteobacteria"), drop = NULL)
# 真菌：只留真菌界（若注释里混了非真菌）
filter$taxon_subset <- list(rank = "Kingdom", keep = "Fungi", drop = NULL)
```
匹配用 `grepl(..., ignore.case = TRUE)`，是**子串匹配**，写 `"Chloro"` 也能命中。

---

## E. 标准化（**不许代选**）

| 键 | 默认 | 可选项 | 推荐 | 理由 |
|---|---|---|---|---|
| `normalize$method` | `rarefy` | `rarefy` / `relative` / `css` / `none` / `subsample_avg` | **rarefy** | 抽平让 α 指数与距离不受测序深度干扰、可复现；`relative` 保留全部 reads 但对深度敏感；`css` 适合深度差异极大时；`subsample_avg` 降低单次抽平随机性但结果是小数 |
| `normalize$depth` | NULL | NULL 或整数 | **NULL = 最小正样本 reads** | 保留最多样本；显式指定会剔除低于该深度的样本 |
| `normalize$subsample_reps` | 20 | 10–50 | 20 | 仅 `subsample_avg` 用 |
| `normalize$by_module` | 见下 | 每个键都可覆盖 | **按默认** | 见下方口径表 |

**逐模块口径默认值**（改它需要科学理由，不是随手改）：

| 键 | 默认 | 理由 |
|---|---|---|
| `composition` | `none` | 组成堆叠图要保留全部 reads，抽平会压缩优势类群、图形失真 |
| `rarefaction` | `none` | 稀释/累积曲线**必须**用原始 reads，抽平后曲线没有意义 |
| `alpha` | `rarefy` | 丰富度类指数对深度极敏感 |
| `distances` | `rarefy` | 距离矩阵会被深度差异污染 |
| `phylo` | `rarefy` | 与 α 保持同口径 |
| `assembly` | `rarefy` | βNTI/RCbray 的零模型假设等深度 |
| `partition` | `rarefy` | iNEXT 用 incidence，但抽平能统一采样强度 |
| `env` | `rarefy` | RDA/VPA 对深度敏感 |
| `model` | `rarefy` | 与 α 保持同口径 |

---

## F. 响应变量（**【不许代选】**）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `response$enabled` | TRUE | 没有响应变量就 FALSE | FALSE 时 70/71 自动跳过 |
| `response$var` | `EF` | **用户指定** | 列必须在 `04_env_response.csv` 里 |
| `response$label` `/unit` | `Ecosystem function` / 空 | 写实际名称与单位 | 进坐标轴与表头，投稿前必须改对 |
| `response$trans` | `none` | `none` / `log` / `sqrt` / `z` | 看 `01_data_health` 的偏度与 Shapiro；右偏重时用 `log` |

---

## G. 系统发育（**【不许代选】** 是否有树）

| 键 | 默认 | 可选项 | 推荐 | 理由 |
|---|---|---|---|---|
| `phylo$enabled` | TRUE | TRUE / FALSE | **有树才 TRUE** | ITS 通常无树；确定无树时关掉可省大量时间并避免无意义结果 |
| `phylo$prev_min` | 2 | 1–5 | 2 | 系统发育 β 只保留出现 ≥ 该次数的类群，避免单例主导 |
| `phylo$root_method` | `midpoint` | `midpoint` / `outgroup` / `none` | `midpoint` | 多数 ASV 树无根；有明确外群时用 `outgroup` |

---

## H. α 多样性

| 键 | 默认 | 可选项 | 推荐 | 理由 |
|---|---|---|---|---|
| `alpha$indices` | 6 个 | `Observed/Chao1/ACE/Shannon/Simpson/Pielou/PD` | **Observed + Shannon + Pielou + PD** | Chao1/ACE 是外推估计，深度不够时虚高；PD 需树 |
| `alpha$test` | `kruskal_dunn` | `anova_tukey` / `kruskal_dunn` / `lmm` / `none` | **看 `01` 的诊断**：正态且方差齐 → `anova_tukey`；否则 `kruskal_dunn`；有重复采样 → `lmm` | 方法选错会直接影响 P 值 |
| `alpha$p_adjust` | `BH` | `BH` / `holm` / `fdr` / `bonferroni` / `none` | BH | 控制 FDR，生态学常用 |
| `alpha$sig_label` | `auto` | `auto` / `letters` / `brackets` / `none` | `auto` | **分组数 > 4 必须字母**，否则括号会挤成一团 |
| `alpha$rarefaction$curve` | `vegan` | `vegan` / `inext` / `both` | `vegan` | vegan 稀疏曲线够用；iNEXT 有外推但更慢 |

---

## I. β 多样性与排序

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `distances$metrics` | `bray,jaccard,wunifrac` | **三选**：Bray（丰度）+ Jaccard（有无）+ wUniFrac（系统发育，需树） | 三种口径互补；无树时 wUniFrac 自动跳过 |
| `ordination$method` | `pcoa` | `pcoa`；`nmds` 只作稳健性补充 | PCoA 与 PERMANOVA 同源、可报轴解释率 |
| `ordination$ellipse_min_n` | 4 | 4–9 | 每组太少时椭圆不可靠 |

## J. 多元差异检验

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `beta_tests$permanova_terms` | TRUE | TRUE | 用 `by="terms"` 才能得到各效应项；**本机 `by="margin"` 只返回末项** |
| `beta_tests$anosim` | TRUE | TRUE | 与 PERMANOVA 互补（基于秩） |
| `beta_tests$permdisp` | TRUE | **TRUE** | 显著 PERMANOVA 需要排除"只是离散度不同" |
| `beta_tests$pairwise` | TRUE | 因子水平 ≤ 6 时 TRUE | 水平太多时对数爆炸 |
| `beta_tests$permutations` | 999 | **999**（试跑可 199） | 论文标准；写进表注 |
| `beta_tests$p_adjust` | `fdr` | `fdr` | 两两比较必须校正 |
| `beta_tests$min_group_n` | 3 | 3 | 每组 < 3 无法做两两 |

## K. 多样性分区

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `partition$engine` | `iNEXT` | `iNEXT`（分类学）；要系统发育分区才用 `iNEXT.3D` | iNEXT.3D 慢很多且需树 |
| `partition$q` | `0,1,2` | **0,1,2 三档** | 分别对应"有没有/常见不常见/优势度"，审稿人常要求 |
| `partition$t_gamma` | NULL | **NULL（自动）** | 自动取各组采样单元数中位数，并对不足的组显式标注；手工剔除组会造成选择性报告 |
| `partition$coverage_aware` | TRUE | TRUE | 覆盖度低时 iNEXT 外推倍数会很大，**必须在表注披露** |

## L. 群落组装零模型（**【不许代选】** runs 与阈值）

| 键 | 默认 | 可选项 | 推荐 | 理由 |
|---|---|---|---|---|
| `assembly$runs` | 999 | 99 / 199 / 999 / 9999 | **999**；试跑可 99–199，**写论文前必须确认是 999** | 零模型精度的决定性参数 |
| `assembly$null_model` | `independentswap` | `independentswap` / `swap` / `richness` / `frequency` | `independentswap` | 与 Stegen 框架一致 |
| `assembly$fixedmar` | `columns` | `columns` / `rows` / `both` / `none` | `columns` | 固定样本丰富度，与 nifH 口径一致 |
| `assembly$bnti_thr` | 2 | 2（惯例） | **2** | Stegen 惯例，改它必须给理由 |
| `assembly$rcbray_thr` | 0.95 | 0.95（惯例） | **0.95** | 同上 |
| `assembly$ncores` | 4 | ≥1 | 物理核数 | 零模型是唯一的重计算瓶颈 |
| `assembly$chunk` | NULL | NULL 或整数 | 大样本（>200）时给分块大小 | RCbray 的 `n² × runs` 数组会爆内存 |
| `assembly$subset_top` | NULL | NULL 或整数 | **NULL**（用全部 ASV） | 只用丰度前 N 会造成选择性报告；非用不可时**必须在表注写明** |

## M. 环境变量

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `env$vars` | 13 个 | 按实际列名改 | 相关矩阵、VIF、PCA 都从这里取 |
| `env$groups` | Climate / Soil 两块 | 按科学假设分组 | 热图分组与解释用 |
| `env$pca_vars` | 11 个土壤指标 | 去掉共线极高的变量 | PCA 输入 |
| `env$pca_log` | FALSE | 量纲差异大时 TRUE | log 变换能压缩量纲；**pH 是对数刻度，log 时应豁免**（`pca_pH_exempt_log`） |
| `env$pca_show` | `sites` | **【不许代选】** sites / species / both | 样方多时同时画变量点会糊成一团 |
| `env$vif_thr` | 10 | 10（宽松）/ 5（严格） | 共线诊断门限 |
| `env$cor_p_adjust` | `BH` | BH | **只对上三角独立配对校正再镜像**（不能对含重复的整矩阵校正） |
| `env$vpa_blocks` | Geo / Soil / Climate | 按科学假设分块 | 键名会出现在图例里 |
| `env$vpa_transform` | `hellinger` | Hellinger | 群落数据的标准预处理 |

## N. 建模

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `model$random` | `(1|Site)` | 有站点就保留 | 站点级伪重复必须吸收 |
| `model$reml` | FALSE | FALSE | 便于用似然比做模型比较 |
| `model$p_method` | `satterthwaite` | Satterthwaite | 小样本下比 Wald 稳 |
| `model$robust_se` | FALSE | 站点少（<10）时 TRUE | 聚类稳健标准误能防过窄的置信区间 |
| `model$edf_nonlin` | 1.1 | 1.1 | GAMM 的有效自由度 > 该值才认为非线性 |
| `model$focus` | `Shannon` | **【不许代选】** 真正关心的变量 | 条件效应曲线的 x 轴 |
| `model$covariates` | `MAT, Soil_PC1` | 按科学假设 | 偏回归的协变量 |
| `model$standardize` | TRUE | TRUE | 系数可横向比较 |

## O. SEM（**【不许代选】** DAG 边表）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `sem$edges` | 空 | **用户必须声明** | 结构层不许自动探索；本包只做拟合与画图 |
| `sem$layers` | 4 层 | 与 `edges$layer` 一致，从上到下 | 决定路径图层次 |
| `sem$labels` | NULL | 给显示名 | 图上不用变量原名更易读 |
| `sem$dsep` | TRUE | TRUE | Fisher's C 是 SEM 的必备检验；不通过必须照实报告 |
| `sem$p_adjust` | `BH` | BH 或 NULL | 对同一模型内全部路径校正；**校正范围要写进表注** |

## P. 空间与地图

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `spatial$level` | `auto` | **auto** | 坐标重复（站内重复）时做样本级空间自相关是伪重复，auto 会自动降到站点级 |
| `spatial$moran_k` | NULL | NULL（自动 k=max(3,min(n/4,10))） | k 太大时权重矩阵被稀释 |
| `spatial$dist_classes` | `0,10,100,600,1600` | 按实际距离跨度改 | 分段要覆盖实际范围且每段有足够连接数 |
| `map$enabled` | FALSE | 有经纬度且要出图才 TRUE | 出图需联网取 DEM |
| `map$offline` | FALSE | 无网时 TRUE | 只用缓存 DEM；没有缓存则降级为散点图 |

## Q. 响应变量关系

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `response_models$linkage` | alpha+ordination+pairwise | 按科学问题挑 | 三类是不同尺度的证据，不必都做 |
| `response_models$alpha_focus` | `Shannon, PD` | 按假设挑 | 多做几个指数会引出多重比较问题 |
| `response_models$ordination_axes` | `PCoA1` | 只用第一轴 | 第二轴的生物学含义通常不清晰 |
| `response_models$pair_metric` | `bray` | bray | 与主分析同口径 |
| `response_models$bootstrap_B` | 1000 | 1000 | 效应量 CI 的精度 |

## R. 图件

### R1 尺寸与字体

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `figure$width_double` | 183 | 183 mm | ISME/Nature 双栏标准；单栏用 `width_single`=90 |
| `figure$base_pt` | 9 | 7–9 | 投稿字号；**必须用 pt 绝对值**，不要用相对字号 |
| `figure$font` | `Arial` | 期刊要求字体 | 只影响字体文件；**图中文字本身必须纯 ASCII**（见 pitfalls B12） |
| `figure$save_singles` | FALSE | FALSE | 每模块一张合成图；投稿要拆图时再开 |
| `figure$id_map` | NULL | 论文图号映射 | 只影响图面标签与文件名前缀，不影响分析 |

### R2 出哪些图（**默认关掉论文里没有的图**）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `figure$skip_modules` | `c("01","05","13","15","22")` | **保持默认**（除非你论文确实要这几张） | 这 5 个模块在 原始研究里**没有对应图**（诊断图 / 只出表的检验）。默认只出结果表；想让某个模块出图就把它从这个向量里删掉。完整对照见 `references/06-figure-map.md` |

### R3 配色（内置 nifH 那套，逐字照搬）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `figure$palette$Zone` | NULL → 内置 `zone` | 留 NULL | 内置 = `#89C4D7 / #4DAFB5 / #2B6AA0 / #3B2C6B`（nifH `cols_zone`） |
| `figure$palette$Vegetation` | NULL → 内置 `veg` | 留 NULL | 内置 = `#989898 / #959BB0 / #60859E / #B98079`（nifH `cols_veg`） |
| `figure$palette$<其他>` | NULL | 要覆盖才填 | 键名 = 因子列名；水平名要对得上，缺色会 `md_warn` 并补兜底色 |
| 组装过程 / 核心稀有 | 内置 `assembly` / `quadrant` | 留 NULL | 来自 nifH `cols_assembly` / `cols_2x2` |

### R4 箱线图样式（默认完全复刻 原实现 `plot_dual_factor_panel`）

| 键 | 默认 | 可选项 | 推荐 | 理由 |
|---|---|---|---|---|
| `plot$box$style` | `nifh` | `nifh` / `simple` | `nifh` | 与论文一致 |
| `plot$box$letters` | `TRUE` | TRUE / FALSE | `TRUE` | 是否画 CLD 显著性字母（**默认与原代码构建方法一致**）。设 `FALSE` 只影响出图，字母本身仍写在 `Descriptive` 表里 |
| `plot$box$color_by` | `auto` | `auto` / `f1` / `f2` / `none` | `auto` | auto = 双因子用**第二因子**、单因子用第一因子（nifH 口径）；**非双因子设计或想统一上色时改成 `none`** |
| `plot$box$fill` | `NA` | `NA` 或颜色 | `NA` | 空心箱体 + 边框上色是 原实现做法 |
| `plot$box$letters_at` | 1.08 | 1.0–1.2 | 1.08 | 字母 y = 每组 max(y) × 该系数 |
| `plot$box$show_n` | `FALSE` | TRUE/FALSE | **FALSE** | nifH 不标样本量；标了会和字母抢位置 |
| `plot$box$kw_annot` | TRUE | TRUE/FALSE | TRUE | 右上角 `Zone: p = 0.033 *; Vegetation: p = 0.486 ns` |
| `plot$box$y_from_zero` | TRUE | TRUE/FALSE | TRUE | 原实现 y 轴从 0 起 |
| `plot$box$y_expand` | 1.35 | 1.1–1.6 | 1.35 | y 上限 = max × 该系数 |
| `plot$box$inward_ticks` | TRUE | TRUE/FALSE | TRUE | nifH 用 `unit(-0.12,"cm")` 朝内刻度 |
| `plot$box$inter_sep` | `"_"` | 任意分隔符 | `"_"` | **必须与生成 CLD 字母时用的分隔符一致**，否则字母被丢弃（pitfalls C7） |
| `plot$box$x_labels` | `two_line` | `two_line` / `one_line` | `two_line` | nifH 是 `Zone\nVegetation` 两行标签；窄面板可改 `one_line` |

## S. 结果表与运行

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `table$detail` | FALSE | FALSE | TRUE 会输出逐样本/逐对明细，大样本时出巨表 |
| `table$write_csv` | FALSE | 下游用 Python 时 TRUE | openxlsx 产物在 Python 侧有已知坑（见 `04-pitfalls.md`） |
| `opts$seed` | **123** | **整篇论文固定一个值** | 全包唯一种子；换值会改变所有随机结果 |
| `opts$force_recalc` | FALSE | 只在改了分析逻辑时 TRUE | 平时走缓存，秒级重出图 |
| `opts$figures_only` | FALSE | 只调样式时 TRUE | 跳过重计算模块 |
| `opts$ncores` | 4 | 物理核数 | 只影响零模型并行 |
| `opts$strict` | TRUE | TRUE | 体检有 BLOCK 级问题就停，不带着坏数据往下跑 |
| `opts$auto_install` | TRUE | TRUE | 缺包先自动 `install.packages()` 再重试；只对 CRAN 有效（SpiecEasi 需 GitHub）。完全离线设 FALSE |
| `opts$repos` | CRAN 官方源 | 就近镜像 | 自动装包用的 CRAN 源 |

---

## T. 2026 补充轮（模块 12 / 44 / 45 / 80 / 81 / 82 / 83）

> 这一组全部是**可选**的：`enabled = FALSE` 时该模块直接 SKIP，不影响其余模块。
> 往已有模块里加的四条增量（`composition$indicator` / `env$rf` /
> `response_models$pattern_test` / `response_models$linkage`）**默认关**，
> 打开会改变该模块的图与表（见 `04-pitfalls.md` C25）。

### T1 `niche`（模块 12 生态位宽度）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `enabled` | TRUE | 做类群稀有性/生态位分析时 TRUE | 不做就设 FALSE |
| `index` | `levins` | `levins` | 原脚本口径（Levins' B）；`shannon` 只作稳健性补充，`both` 两列都进表 |
| `level` | `both` | `both` | ASV 级看类群策略，样本级看群落策略构成，两件事 |
| `input` | `raw` | `raw` | 生态位宽度对测序深度敏感，抽平会削平宽度（与原脚本一致） |
| `thr_abund` | 0.1 | 0.1 | 四象限的平均相对丰度阈值（%），原 `classify_2x2` 口径 |
| `thr_prev` | 10 | 10 | 四象限的出现率阈值（%），原脚本口径 |
| `groups` | NULL | 保持 NULL | NULL = **前两个设计因子的交互组**（原实现）；给字符串则强制单因子 |
| `sensitivity_occ` / `sensitivity_abund` | 原脚本 5×3 网格 | 保持 | 阈值敏感性扫描（原 SCI-8），验"核心/稀有划分是不是阈值造出来的" |

### T2 `constrained_ord`（模块 44 约束排序）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `enabled` | TRUE | 要"环境→群落"的约束排序时 TRUE | |
| `method` | `dbrda` | `dbrda` | 与 PERMANOVA 同源（都基于距离），可直接报 R²adj 与 P；`rda` 用欧氏、`cca` 用卡方、`capscale` 是 db-RDA 的通用写法 |
| `response` | `community` | `community` | 也可 `distance`（任意距离矩阵）；**换论文时因变量不必是群落** |
| `response_vars` | NULL | NULL | NULL = 用全部 ASV |
| `response_distance` | NULL | 有树时 `wunifrac` | `response="distance"` 时用哪个距离 |
| `constraints` | NULL | 按科学假设 | NULL = 用 `env$vars`；解释变量集合决定结论，必须事前定 |
| `transform` | `hellinger` | `hellinger` | 群落数据的标准预处理（只对 `response="community"` 生效） |
| `vif_thr` | 15 | 10（严格）/ 15（跟原脚本） | 原脚本用 15；审稿人常问共线，严格场景用 10 |
| `permutations` | 999 | 999 | 论文标准 |
| `arrow_ratio` | 0.78 | 0.78 | 最长 env 箭头 = 站点云半径 × 该系数（原脚本口径，替代硬编码 1.8） |

### T3 `response_curves`（模块 45 GAM/GAMM 响应曲线）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `enabled` | TRUE | 要响应形状时 TRUE | |
| `y` | `PCoA1` | 任意连续列 | **不必是群落**：排序轴 / α 指数 / 环境量 / 响应变量都行 |
| `x` | NULL | NULL | NULL = 用 `env$vars`；也可显式给列名向量 |
| `engine` | `gam` | `gam` | 论文口径（`s(x, k)` + REML）；`random` 非空时自动升 `gamm` |
| `k` | 4 | 4–5 | 平滑基维数；太小欠拟合、太大过拟合（原脚本 4） |
| `random` | NULL | 有重复采样时给 | 例如 `"(1\|Site)"`；给出则走 GAMM |
| `p_max` | 0.1 | 0.1 | 只画 P 小于该值的平滑项（原脚本口径）；全超阈值时保留 P 最小者 |
| `max_panels` | 3 | 3 | 原脚本用 RF top 3；面板太多会挤 |
| `edf_nonlin` | 1.1 | 1.1 | edf > 该值才认为非线性 |
| `x_from` | `"rf_top"` | `"rf_top"` | 先用随机森林给 `env$vars` 排序、取前 `max_panels` 个（对齐原脚本 `gam_preds <- names(sort(rf_A_pcoa1$variable.importance))[1:3]`）。读不到 40/71 的 RF 结果时，模块 45 **自己跑一次** `md_rf_fit()`；再不行才回退 `x` |
| `rf` | `list(engine="ranger", num_trees=2000)` | 保持 | `x_from="rf_top"` 且读不到 40/71 的 RF 结果时，模块 45 自己跑一次用的参数 |
| （已删除）`style$show_edf_bar` | —— | —— | 原实现的 **B 面板（EDF / -log10(P) 条形）是自造的**，原代码没有这张图，2026-09-26 已删除；本模块现在只有 GAM 曲线一组 |
| `p_max` | 0.1 | 0.1 | 只画 P 小于该值的平滑项（原脚本 口径）；**一个都没通过时按与模块 83 相同的回退口径只显示 P 最小的一个**，并在面板内标注 `p > 0.10 (not significant)` |

### T4 `network`（模块 80 共现网络 —— **会出 4 张图**）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `enabled` | TRUE | 做网络分析时 TRUE | 不做设 FALSE（这个模块最重） |
| `method` | `cor` | `cor` | 论文的"SparCC 近似"实际就是"相关 + 阈值"；`sparcc` 需 SpiecEasi（GitHub，且显著更慢），缺包会自动回退 |
| `cor_method` | `spearman` | `spearman` | 丰度数据非正态 |
| `min_cor` | 0.5 | 0.4–0.6 | 全局网络 `\|r\|` 阈值。**它决定节点数与模块数，必须在表注写明** |
| `sub_min_cor` | 0.5 | 同 `min_cor` | 丰富/稀有亚网络的阈值 |
| `min_prev` | 0.02 | 2–5 个样本 | ASV 至少在多少比例的样本中出现；数据太小会被全部过滤（见 pitfalls C20） |
| `max_asv` | 1500 | 500–2000 | 参与建网的最大 ASV 数（按平均丰度降序）。相关是 O(n²)，1000 个约 30 s |
| `module_algo` | `louvain` | `louvain` | Louvain 免调参、论文口径；`fast_greedy`/`walktrap` 作稳健性 |
| `zi_thr` / `pi_thr` | 2.5 / 0.62 | 保持 | Deng et al. 2012 的节点角色阈值；改它必须给理由 |
| `random_n` | 100 | 100 | Erdos-Renyi 随机网络次数（原脚本口径） |
| `robustness_nsim` | 30 | 30–100 | 网络韧性模拟次数；**这一项之前必须 `md_seed`**（原脚本漏了，见 pitfalls C14） |
| `n_boot` | 1000 | 1000 | 分组网络指标 bootstrap 均衡样本量的次数 |
| `metrics_main` / `metrics_supp` | 原脚本 `metric_order_main` / `_supp` | 保持 | 各 6 个样本级指标，主图 / 附图 |
| `metrics_supp_figure` | TRUE | TRUE | 设 FALSE 可把总图数 31 → 30 |
| `group_factors` | NULL | 前 2 个设计因子 | NULL = `design$factors[1:2]`（原脚本 Zone + Vegetation） |
| `sub_network_by` | `abundance_rank` | `abundance_rank` | 亚网络按丰度×出现率四象限分组（原脚本口径） |
| `turnover_by` | NULL | 第二个设计因子 | Sankey 节点周转按它的相邻水平对比较 |
| `style$layout` / `style$layout_group` | `graphopt` / `fr` | 保持 | **布局在出图段算**，属于样式键（不进缓存键） |
| `style$wide_layout` | FALSE | FALSE | 节点属性图（80d）：FALSE = 183 mm 三格纵向（与全包其余图同观感）；TRUE = 复刻原文 原图 533.4×177.8 mm 单排 |
| `style$*` | `MD_NET` | 一般不用改 | 节点/边尺寸、透明度、字号、top-N 都走 `style` |

### T5 `taxon_group`（模块 81 通用类群分组）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `enabled` | TRUE | 有映射表时 TRUE | 没表设 FALSE |
| `file` | `data/05_taxon_group.csv` | 保持 | NULL = 跳过本模块 |
| `taxon_col` / `group_col` | `taxon` / `group` | 保持 | 映射表的列名 |
| `match_by` | `Genus` | 按你的注释深度 | `ASV`（直接匹配 `asv_id`）或 `Phylum`…`Species`；**分组方式不限**（功能群 / 生态位 / R-K 对策 / 自定义） |
| `group_order` | NULL | 建议显式给 | NULL = 按映射表出现顺序；给向量可固定图例与堆叠顺序 |
| `group_label` | `Group` | 按语义写（ASCII） | 图例名，会出现在图上 |
| `drop_unmapped` | FALSE | FALSE | TRUE 会把映射表没覆盖的类群整个丢掉（**会改变总丰度构成，慎用**） |
| `indicator$enabled` | FALSE | 需要时 TRUE | 追加类群级 IndVal（多一个面板/表） |

### T6 `taxon_share`（模块 82 类群共享）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `enabled` | TRUE | 要比较各处理的类群重叠时 TRUE | |
| `factor` | NULL | 主设计因子 | NULL = `design$factors[1]` |
| `rank` | `ASV` | `ASV`（原脚本口径）或 `Genus` | 类群级共享更稳（ASV 级受单例影响极大） |
| `min_prev` | 1 | 2–5 | 类群至少在几个样本中出现才纳入；**1 会让单例把 UpSet 撑爆** |
| `max_comb` | 20 | 20 | 进入 UpSet 的交集组合上限（组合数指数增长，见 pitfalls C17） |
| `upset_top_n` | 20 | 10–20 | 顶部柱状图显示的交集个数 |
| `show_venn` | TRUE | **已弃用** | 键保留但不生效：2026-09-26 起本模块**只出 UpSet 一张图**，Venn / 共享条形面板已删除 |

### T7 `taxon_diff`（模块 83 差异类群）

| 键 | 默认 | 推荐 | 理由 |
|---|---|---|---|
| `enabled` | TRUE | 要做差异类群时 TRUE | |
| `rank` | `Genus` | `Genus` / `Family` | 聚合层级会影响"差异"的检出力；属级最常用（原脚本口径） |
| `clr` | TRUE | TRUE | CLR 是组成数据的标准变换（原脚本口径） |
| `pseudocount` | NULL | NULL | NULL = **最小正值的 1/2**（原脚本口径）；换值会改变所有 P，必须写表注 |
| `min_prev_frac` | 0.05 | 0.05–0.1 | 类群至少在多少比例的样本中出现（原脚本 0.05） |
| `min_n` | 3 | 3 | 每组最少样本数，低于则该类群不检验（原脚本口径） |
| `test` | `kruskal` | `kruskal` | 非参数，适合组成数据；`wilcox` 只在两组时用 |
| `posthoc` | `dunn` | `dunn` | `kruskal` 的事后标准选择；`none` 只报整体 P |
| `adjust` | `BH` | `BH` | 生态学惯例 |
| `effect` | `cliff` | `cliff` | 非参数效应量，与 KW/Dunn 同一族（原脚本口径） |
| `p_thr` | 0.05 | 0.05 | 显著阈值（用 BH 校正后的 p） |
| `top_n` | 25 | 15–30 | 热图显示多少个差异类群（再多行标签读不清） |
| `heatmap_engine` | `ComplexHeatmap` | `ComplexHeatmap` | 原脚本口径；缺包自动退 `pheatmap`，再退 `geom_tile` |
| `heatmap_scale` | `row_z` | `row_z` | 行 z 化后差异看得清；`rel` / `clr` 作对照 |
| `show_box` | FALSE | 需要时 TRUE | 追加 Top 差异类群箱线面板 |

### T8 已有模块的四条可选增量

| 模块 | 键 | 默认 | 打开后 |
|---|---|---|---|
| `10_composition` | `composition$indicator$enabled` | FALSE | 面板 D 指示种 IndVal（`func`/`nperm`/`top_n`/`fdr_thr`/`factors` 可调）+ `Indicator_species` sheet |
| `40_env_summary` | `env$rf$enabled` | FALSE | 填上**预留的 D 区**：随机森林重要性条形 + `RF_importance` sheet；`response` / `predictors` / `engine`(ranger/randomForest/cforest) / `num_trees` / `top_n` 可调；`pdp = TRUE` 再追加 E 面板 |
| `70_response_pattern` | `response_models$pattern_test` | `"kruskal_dunn"` | 改 `"scheirer_ray_hare"` 用双因子非参数检验（需 rcompanion，缺包退 ART + 双因素 ANOVA），额外写 `SRH_two_way` sheet |
| `71_response_linkage` | `response_models$linkage` | `c("alpha","ordination","pairwise")` | 追加 `"procrustes"`（C 面板 + 2 sheet）与 `"rf"`（D 面板 + 1 sheet）；`procrustes$axes/permutations`、`rf$engine/num_trees/pcoa_axes/rank/top_n` 可调 |
