# 07 · 备用件（spare）—— 默认全部关闭

> **备用件 = 默认不生效、不影响论文图件与数值的可选能力。**
> 主流程只保留「你论文里有的分析」；凡是**换数据时可能有用、但你论文没用/不常用**的东西，
> 一律收在这里，默认关着。开启方式见各条。

---

## A. 备用面板（出图用）

这两个面板都是**你论文里确实有、但不常用**的东西，因此做成**通用可配**而不是写死成一篇论文专用。

### A1 · 样本级指标 vs 群落周转（原论文 Fig.3D 的通用化）

| 项 | 内容 |
|---|---|
| 归属模块 | `31_assembly_process` |
| 原论文形态 | `lm(MeanBeta ~ DriftProp + Zone)`——漂变占比 vs 平均 β 多样性 |
| 通用化后 | **横轴 = 任意样本级连续列**（过程占比 / α 指数 / 环境变量…）；**纵轴 = 该样本与其余样本的平均群落距离**（"群落周转"）；可加分组着色与协变量；`lm` 或 `gam` |
| 开关 | `CFG$spare$process_vs_turnover$enabled`（默认 `FALSE`） |
| 关键参数 | `x`（横轴列名，必填才会出图）· `y`（纵轴标签，默认 `MeanBeta`）· `group` · `covar` · `fit` |

开启示例（复现你论文那张图）：
```r
CFG$spare$process_vs_turnover <- list(
  enabled = TRUE, x = "Drift", y = "MeanBeta", group = "Zone", covar = NULL, fit = "lm"
)
```
> 横轴列取自 `md_null_cache(CFG)$pair_df` 展平出的**样本级过程占比**，
> 或 `md_analysis_with_derived(CFG)` 里的任意连续列；纵轴由 `md_dist_get(CFG, <度量>)` 现算。

### A2 · 焦点变量的标准化系数森林图（原论文 Fig.S10 右半的通用化）

| 项 | 内容 |
|---|---|
| 归属模块 | `50_lmm` |
| 原论文形态 | 组装过程占比 → 响应变量的标准化 β 森林图 |
| 通用化后 | 对 `focal` 里**每个变量**各拟合 `response ~ focal_i + covariates [+ random]`，把**标准化系数 β**、SE、P 收成长表，用 `md_forest()` 画并列森林图；可选 `by` 分组并列 |
| 开关 | `CFG$spare$std_beta$enabled`（默认 `FALSE`） |
| 关键参数 | `response`（默认用 `CFG$response$var`）· `focal`（一组变量，必填才会出图）· `covariates` · `random`（默认用 `CFG$model$random`）· `by` |

开启示例：
```r
CFG$spare$std_beta <- list(
  enabled = TRUE, response = NULL,
  focal = c("Drift", "Homogenizing dispersal", "Variable selection"),
  covariates = c("MAT", "Soil_PC1"), random = "(1|Site)", by = NULL
)
```

---

## B. 备用口径开关（分析用）

这些开关**都住在各自的配置段里**（不是 `CFG$spare`），**默认值都等于你论文的口径**。
换数据时可能要，所以你论文的图与数值都不受影响。

| 开关 | 位置 | 默认 | 你论文用的 | 说明 |
|---|---|---|---|---|
| `normalize$method` | `CFG$normalize` | `"rarefy"` | `rarefy` | 可选 `relative` / `css` / `none` / `subsample_avg`；后两者是你论文没用过的口径 |
| `ordination$method` | `CFG$ordination` | `"pcoa"` | `pcoa` | 可选 `nmds` / `both`；你论文**只做 PCoA** |
| `beta_tests$mrpp` | `CFG$beta_tests` | `FALSE` | 未做 | MRPP 检验 |
| `assembly$subset_top` | `CFG$assembly` | `NULL` | 未做 | 只用丰度前 N 个 ASV 跑零模型（大数据降规模）；**用了必须在表注披露** |
| `assembly$chunk` | `CFG$assembly` | `NULL` | 未做 | 零模型分块（内存控制，**不改变统计口径**） |
| `distances$metrics` | `CFG$distances` | bray+jaccard+wunifrac | 同样三种 | 距离口径集合，按需增删 |
| `filter$taxon_subset` | `CFG$filter` | `NULL` | Cluster IV 过滤 | 通用化的**类群子集过滤**机制（你论文的 Cluster IV 过滤就是它的一个实例） |
| `composition$indicator$enabled` | `CFG$composition` | `FALSE` | 未做（毕业论文做过） | 指示种 IndVal 面板（2026 补充轮） |
| `env$rf$enabled` | `CFG$env` | `FALSE` | 未做（毕业论文做过） | 随机森林重要性面板（2026 补充轮） |
| `response_models$pattern_test` | `CFG$response_models` | `"kruskal_dunn"` | `kruskal_dunn` | 可切 `"scheirer_ray_hare"`（2026 补充轮） |
| `response_models$linkage` | `CFG$response_models` | alpha+ordination+pairwise | 同样三项 | 可加 `"procrustes"` / `"rf"`（2026 补充轮） |

> 这些开关的代码路径**都保留着**，但默认不触发。要删可以，删之前告知。

---

## C. 已移出管线的模块

| 模块 | 内容 | 状态 |
|---|---|---|
| ~~`42_constrained_ord`~~ | RDA / CCA / db-RDA + envfit | **已被 `44_constrained_ord` 取代**（按毕业论文 `3_environmental_drivers.R:855-1000` 的 db-RDA 样式重写，不再是自编版）。旧文件归档在**包外**构建记录 `analysis-packages/_build-notes/micro-div-pipeline/removed-code/42_constrained_ord.R`，仅作追溯（随本包单独分发时不存在）；`44` 与它的区别：支持 `community`/`distance` 两种响应侧、`method` 可选 dbrda/rda/cca/capscale、VIF 阈值跟原脚本用 15、箭头长度按站点云半径自动归一。 |

---

## D. 纪律

1. **备用件默认必须关着**，且开启后**不得改变主流程任何已有图件与数值**（只允许"多一个面板/多一张图"）。
2. 备用件**不许写死成某篇论文专用**：能用配置表达的，一律走配置。
3. 新增备用件时：写进本文件 + 在 `config_template.R` 的 `CFG$spare`（或对应段）里给默认值，
   并在模块的 `Note` sheet 里记明"该面板/口径为可选备用"。
