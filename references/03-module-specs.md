# 03 · 模块规格（签名 / 输入 / 输出 / 验收要点）

> 所有模块函数签名统一为 `md_xxx(CFG)`，返回 `md_result(module, status, reason, tables, figures, notes)`。
> status ∈ `OK` / `SKIP` / `PARTIAL` / `EMPTY`。
> 数据访问一律走 `md_core_*()` / `md_dist_get()` / `md_null_cache()` / `md_env_cache()` / `md_phylo_cache()`。
> **不要直接 `readRDS` 别人的缓存。**

---

## 组 0 · 数据（01–05）

### `01_data_health` — `md_data_health(CFG)`
| 项 | 内容 |
|---|---|
| 输入 | 五文件（`md_read_inputs`） |
| 产出 | `results/01_data_health.xlsx`（sheets：Overview / Files / Issues / Taxonomy / Factors / Factor_crosstab / Tree_check / Sample_reads / Response_diagnostics / Env_missing / `_meta`）<br>`figures/01_data_health.pdf`（A 每样本读长 B ASV 丰度-出现率 C 环境缺失率）<br>`cache/01_health.rds`（含 `issues`、`n_block`、`fingerprint`） |
| 硬门 | 有 BLOCK 级问题时 `stop()`（`opts$strict=TRUE` 时）；否则降为 WARN |
| 验收 | Issues sheet 里 BLOCK 数量与日志一致；Files sheet 有五个 md5 |
| 不做 | 一律不改数据（不插补、不删样本、不改名） |

### `05_assemble` — `md_assemble(CFG)`
| 项 | 内容 |
|---|---|
| 输入 | 五文件 |
| 顺序 | **过滤（ASV 读长 → 出现率 → 样本读长 → 类群子集）→ 再抽平**（顺序不能反，见 pitfalls A9） |
| 口径 | 按 `CFG$normalize$by_module` 逐模块定（composition/rarefaction 用 `none`，其余用 `rarefy`） |
| 产出 | `results/05_assemble.xlsx`（Filter_log / Taxon_subset / Normalization / Coverage / Factor_sizes / Depth_sensitivity / `_meta`）<br>`figures/05_assemble.pdf`<br>**`cache/05_core.rds`**（全下游的输入） |
| 验收 | `Filter_log` 逐步记录 ASV/样本数是单调下降；`Normalization` sheet 里有实际抽平深度 |

---

## 组 1 · 组成与多样性（10–17）

### `10_composition` — `md_composition(CFG)`
- **口径**：`md_core_asv(CFG, "composition")`（未抽平）
- **做**：任意层级组成堆叠柱 + Top-N 类群 + 核心/稀有类群四象限（丰度 × 出现率双阈值）+ ASV 丰度出现率表
- **产出**：`results/10_composition.xlsx`、`figures/10_composition.pdf`
- **验收**：图例不被裁（用 `md_guide_flex`）；核心类群为 0 时照常输出并在 Note 说明

### `11_alpha` — `md_alpha(CFG)`
- **口径**：`md_core_asv(CFG, "alpha")`（抽平）；稀释曲线用 `"rarefaction"`（未抽平）
- **做**：指数计算（Observed/Chao1/ACE/Shannon/Simpson/Pielou/PD）+ 组间检验（`CFG$alpha$test`）+ CLD 字母 + 稀释/累积曲线
- **产出**：`results/11_alpha.xlsx`（Index_values / Descriptive / Test / Posthoc / Letters / Rarefaction / `_meta`）、`figures/11_alpha.pdf`
- **降级**：无树时跳过 PD 并在日志与表里说明
- **验收**：字母与事后检验 p 矩阵自洽；每组样本量 < 2 时该因子跳过

### `13_distances` — `md_distances(CFG)`
- **口径**：`md_core_asv(CFG, "distances")`
- **做**：Bray / Jaccard(binary) / UniFrac / wUniFrac；无树自动跳过系统发育口径
- **产出**：**`cache/13_dist.rds`** = `list(dist=<具名 dist 列表>, meta=, params=)`、`results/13_distances.xlsx`、`figures/13_distances.pdf`
- **硬约束**：距离矩阵**只存 rds**，绝不并入 data.frame
- **验收**：下游 `md_dist_get(CFG,"bray")` 能取到对象；无树时 `dist` 里没有 unifrac 键且日志有说明

### `14_ordination` — `md_ordination(CFG)`
- **输入**：`md_dist_get(CFG, metric)`
- **做**：PCoA（Cailliez 校正，轴解释率 = 特征值 / 正特征值之和）与/或 NMDS；95% 椭圆（每组 ≥ `ellipse_min_n`）
- **产出**：`results/14_ordination.xlsx`（Axis_variance / Scores / Ellipse / `_meta`）、`figures/14_ordination.pdf`
- **验收**：轴解释率之和 ≤ 100%；椭圆样本不足时跳过并记 NA

### `15_beta_tests` — `md_beta_tests(CFG)`
- **做**：PERMANOVA（**`by="terms"`**）+ ANOSIM + PERMDISP + 可选 MRPP + 两两比较（每组 < `min_group_n` 记 NA，`p_adjust` 校正）
- **产出**：`results/15_beta_tests.xlsx`（PERMANOVA_terms / PERMANOVA_pairwise / ANOSIM / PERMDISP / MRPP / `_meta`）、`figures/15_beta_tests.pdf`
- **验收**：交互模型的效应项齐全（不是只有 Residual）；**必须用 `by="terms"`**（pitfalls A4）

### `16_beta_components` — `md_beta_components(CFG)`
- **做**：`betapart::beta.pair(index.family="sorensen")` 的 turnover/nestedness；可选 Bray 分解作敏感性；分类 vs 系统发育 β 的 Mantel
- **产出**：`results/16_beta_components.xlsx`、`figures/16_beta_components.pdf`
- **验收**：无树时跳过系统发育对比并说明

### `17_gamma_partition` — `md_gamma_partition(CFG)`
- **做**：iNEXT / iNEXT.3D 的 α/β/γ 分区，Hill 阶 `q`；`t_gamma = NULL` 时**自动判定**并对样本不足的组显式标注；报 AsyEst 与覆盖度 SC
- **产出**：`results/17_gamma_partition.xlsx`（Partition / AsyEst / Coverage / Sensitivity / `_meta`）、`figures/17_gamma_partition.pdf`
- **heavy**：TRUE；`figures_only` 模式跳过
- **验收**：覆盖度必须与分区结果同表出现；分组样本不足时该组记 NA 而非整模块失败

---

## 组 2 · 系统发育（20–22）

### `20_phylo_prep` — `md_phylo_prep(CFG)`
- **前置**：`md_precheck(..., need_phylo = TRUE)`；不满足则 SKIP
- **做**：tip 对齐检查、按 `CFG$phylo$root_method` 生根、`keep.tip`、多分支检测
- **产出**：**`cache/20_phylo_prep.rds`**（tree / tip 覆盖率 / 生根方式）、`results/20_phylo_prep.xlsx`、`figures/20_phylo_prep.pdf`（fan + rectangular）

### `21_phylo_diversity` — `md_phylo_diversity(CFG)`
- **做**：**Faith's PD + MPD/MNTD + NRI/NTI 合并**（你指定）。`ses.mpd/ses.mntd(null.model="taxa.labels", runs=CFG$assembly$runs, abundance.weighted=TRUE)`；NRI=−SES_MPD、NTI=−SES_MNTD
- **产出**：**`cache/21_phylo.rds`** = `data.frame(sample_id, PD, MPD, MNTD, NRI, NTI, 各因子列)`、`results/21_phylo_diversity.xlsx`、`figures/21_phylo_diversity.pdf`
- **heavy**：TRUE；并行用 `md_parallel()`，结束 `md_parallel_off()`
- **验收**：SES 非有限时记 NA；无树 SKIP

### `22_phylo_beta` — `md_phylo_beta(CFG)`
- **做**：`betapart::phylo.beta.pair(prev_min=CFG$phylo$prev_min)` 分解 + 与分类学 β 对比（Mantel / bootstrap ρ）
- **产出**：`results/22_phylo_beta.xlsx`、`figures/22_phylo_beta.pdf`

---

## 组 3 · 群落组装（30–34）

### `30_null_models` — `md_null_models(CFG)`
- **做**：βNTI（`picante::comdistnt` + `randomizeMatrix(null.model=)` × `runs`）与 RCbray（`vegan::permatfull(fixedmar=)` × `runs`）；**分块 + 并行 + rds 缓存**三件套
- **产出**：**`cache/30_null.rds`** = `list(bnti, rc, pair_df, process, params)`（`pair_df` 列名见 `MODULE_SPEC` §8）、`results/30_null_models.xlsx`、`figures/30_null_models.pdf`
- **heavy**：TRUE；必须打印内存与耗时预估
- **验收**：`pair_df` 与矩阵一致；`subset_top` 非 NULL 时表注写明代价

### `31_assembly_process` — `md_assembly_process(CFG)`
- **输入**：`md_null_cache(CFG)`
- **做**：Stegen 五过程判定（阈值 `bnti_thr`/`rcbray_thr`）+ 比例汇总（总体/因子1/因子2/交互）+ 四象限散点
- **产出**：`results/31_assembly_process.xlsx`、`figures/31_assembly_process.pdf`
- **注意**：**距离衰减不在本模块**（在 33）

### `33_distance_decay` — `md_distance_decay(CFG)`
- **做**：地理（Haversine）/ 环境（标准化欧氏）/ 群落（1 − Bray）三口径；`lm(log(Similarity) ~ log(Dist))`；最少对数门限参数化
- **产出**：`results/33_distance_decay.xlsx`、`figures/33_distance_decay.pdf`

### `34_variance_partition` — `md_variance_partition(CFG)`
- **做**：`vegan::varpart`，分块来自 `CFG$env$vpa_blocks`，Hellinger 转换，偏 RDA 显著性 + `holm`
- **产出**：`results/34_variance_partition.xlsx`（Partition / Significance / `_meta`）、`figures/34_variance_partition.pdf`
- **注意**：**校正 R² 为负要照实输出**（不截断为 0），并在表注解释

---

## 组 4 · 环境关联（40–43）

### `40_env_summary` — `md_env_summary(CFG)`
- **做**：描述统计 + VIF + 相关矩阵（上三角独立配对 BH 后镜像）+ 组均值热图 + 环境 PCA 双标图 + **`Soil_PC1`**
- **产出**：**`cache/40_env.rds`** = `list(env, vif, cor, pca, soil_pca, derived, ...)`（`derived` 至少含 `sample_id` 与 `Soil_PC1`）、`results/40_env_summary.xlsx`、`figures/40_env_summary.pdf`
- **硬约束**：`Soil_PC1` = `prcomp(scale(X))` 第一主轴得分，X = `CFG$env$pca_vars`，log 由 `CFG$env$pca_log` 定、**pH 在 log 时豁免**；全包唯一，不许别的模块重算

### `43_mantel_mrm` — `md_mantel_mrm(CFG)`
- **做**：Mantel / partial Mantel / MRM（样本级与成对级双口径）+ 效应量 bootstrap CI + **饱和与膨胀诊断**
- **产出**：`results/43_mantel_mrm.xlsx`、`figures/43_mantel_mrm.pdf`
- **硬约束**：置换层级必须在表里写清（对象级 n = 样本数）

---

## 组 5 · 建模（50–53）

### `50_lmm` — `md_lmm(CFG)`
- **做**：`lmerTest::lmer`（Satterthwaite P）+ `mgcv::gamm`（edf > `edf_nonlin` 判非线性）+ `performance::r2()` + 可选聚类稳健 SE + **站点级共线守卫**（`md_within_unique()`，站内唯一值=1 时降协变量并禁平滑）+ 条件效应曲线（z 化 → 网格 → 反变换）
- **输入**：`md_analysis_with_derived(CFG)`
- **产出**：`results/50_lmm.xlsx`（Coefficients / Fit / Nonlinearity / Collinearity_guard / `_meta`）、`figures/50_lmm.pdf`

### `52_sem` — `md_sem(CFG)`
- **做**：`CFG$sem$edges` 声明式边表 → `piecewiseSEM::psem` → `fisherC` d-separation → 路径系数（可选 BH）→ R² → 分层路径图（跨层绕行 + 标签避让；正 #2166AC / 负 #B2182B）
- **降级**：边表为空 → SKIP 并说明"用户未声明 DAG"
- **产出**：`results/52_sem.xlsx`（Path_coefficients / Model_fit / R2 / DSep / Note）、`figures/52_sem.pdf`
- **验收**：d-separation 被拒绝时必须照实报告，不许只报拟合好的部分

### `53_sensitivity` — `md_sensitivity(CFG)`
- **做**：按 `CFG$sensitivity$drop_levels` / `drop_sites` 剔组后重跑指定模型并与全量并列对比
- **产出**：`results/53_sensitivity.xlsx`、`figures/53_sensitivity.pdf`
- **注意**：未配置任何剔除项时 SKIP 并说明

---

## 组 6 · 空间（60–61）

### `60_spatial_autocorr` — `md_spatial_autocorr(CFG)`
- **做**：**Moran's I + correlogram 合并**（你指定）。按 `CFG$spatial$level`（auto/sample/site）决定层级；**auto 时坐标有重复则强制站点级**并登记"样本级不可检验"；k-NN 反距离权重；原始值 + 环境残差两种
- **产出**：`results/60_spatial_autocorr.xlsx`、`figures/60_spatial_autocorr.pdf`（Moran's I 条形 + correlogram）

### `61_map` — `md_map(CFG)`
- **前置**：`CFG$map$enabled = FALSE` 时 SKIP
- **做**：DEM + 站点 + 放大图 + 可选世界 inset；DEM 缓存 `cache/61_dem_z<N>.rds`；**无网/取不到 DEM 时降级为无 DEM 散点图**（不报错）
- **产出**：`results/61_map.xlsx`、`figures/61_map.pdf`

---

## 组 7 · 响应变量关系（70–71）

### `70_response_pattern` — `md_response_pattern(CFG)`
- **前置**：`md_precheck(..., need_response = TRUE)`
- **做**：双因子交互箱线（KW + Dunn + CLD）+ 单因子箱线 + 与连续变量的散点/回归
- **产出**：`results/70_response_pattern.xlsx`、`figures/70_response_pattern.pdf`

### `71_response_linkage` — `md_response_linkage(CFG)`
- **做**：按 `CFG$response_models$linkage` 三类：
  - `alpha`：响应 vs α 指数（单回归 + 偏回归 + LMM/GAMM）
  - `ordination`：响应 vs 排序轴条件效应
  - `pairwise`：成对 |Δ响应| vs 群落距离（MRM/Mantel + bootstrap CI 森林图 + 饱和提示）
- **输入**：`md_dist_get()`、`md_analysis_with_derived()`、`md_phylo_cache()`
- **产出**：`results/71_response_linkage.xlsx`、`figures/71_response_linkage.pdf`

---

## 组 8 · 运行与导出（90–95）

### `90_run_all` — runner（不是模块函数）
- **输入**：`run/manifest.R`（唯一清单）+ config
- **参数**：`--config=` `--project=` `--only=` `--figures-only` `--force`
- **行为**：拓扑排序执行；缺文件/缺包 → SKIP；每模块 `tryCatch` 不中断；warning 用 `withCallingHandlers` 不截断；有 ERROR → 非零退出
- **产出**：`logs/<label>_<run_id>.log`、`logs/summary_<run_id>.csv`、`cache/00_run_state.rds`
- **同时承担**原 `91_regenerate_figures` 的职能（`--figures-only`）

### `95_export_tables` — `md_export_tables(CFG)`
- **做**：清点 `results/*.xlsx` 与 `figures/*.pdf`，生成 `results/00_index.xlsx`（Index_tables / Index_figures / Validation / Run_summary / Skipped / Missing_packages / Warnings / Session）+ **后置校验**（PDF > 0 且可渲染、xlsx 按 sheet 名可读且非空、`data/` md5 未变）
- **位置**：manifest 里 `last = TRUE`，永远最后跑
- **产出**：`results/00_index.xlsx`
- **验收**：Validation sheet 里没有 FAIL

---

---

## 组 9 · 生态位 / 约束排序 / 响应曲线（12 / 44 / 45，2026 补充轮）

### `12_niche_breadth` — `md_niche_breadth(CFG)`
- **源**：原实现（原脚本）
- **输入**：`md_core(CFG)`（`$asv_raw` 或 `$asv_rare`，取决于 `niche$input`）
- **算**：`md_niche_table()`（Levins' B / Shannon B + 丰度×出现率四象限 + 出现样本数）+ `md_sample_strategy()`（样本级四类细胞占比）+ 阈值敏感性扫描
- **图**：`12_niche_breadth.pdf`（**200 × 150 mm**，逐项对齐原图 样本生态策略组成堆叠柱（x = **前两个设计因子的交互组**，alpha 0.85 + 黑边 0.3，标签 45° 两行）｜ B 丰度×出现率四象限散点（**从 10_composition 移入**，复刻原脚本 p_07sa：固定 0.1/10 阈值线 + 4 个角标 + log10 y + xlim 0-40，**无曲线**）｜ C Levins B vs 出现率 + **gam 平滑 + 95% CI**（复刻 p_07sb：`method='gam', formula=y~s(x,bs='cs'), colour='#D55E00', alpha=0.15`，k 自适应）
- **表**：`Niche_Breadth` \| `Cell_summary` \| `Sample_strategy` \| `Threshold_sensitivity` \| `Sample_stats` \| `Note`
- **开关**：`CFG$niche$enabled`；`index`(levins/shannon/both)、`level`(asv/sample/both)、`input`(raw/rarefied)、`thr_abund`(0.1)、`thr_prev`(10)、`groups`
- **验收**：角标文字纯 ASCII（`Core` / `CondRare` / `TransRare` / `PermRare`）；`index="both"` 时表里两列都有、图仍用 `Levins_B`

### `44_constrained_ord` — `md_constrained_ord(CFG)`
- **源**：原实现（Part A db-RDA triplot，原图）
- **输入**：`md_core_asv(CFG, "alpha")`（响应侧，**需转置为样本×变量**再传给 `md_cord_fit`）或 `md_dist_get()`；`md_core_analysis(CFG)` 取解释变量
- **算**：`md_cord_fit()`：VIF>15 逐步剔除 → `anova(perm=999)` → `anova(by="terms")` → `RsquareAdj` → `scores("sites"/"bp")` → 轴解释率
- **图**：`44_constrained_ord.pdf`（单面板，183 × 110 mm）triplot：样方点按 `factors[1]` 着色 / `factors[2]` 定形状 + env 箭头（长度按 `arrow_ratio` 归一）+ ggrepel 标签 + `R2adj`/P 副标题
- **表**：`Model_summary` \| `Terms` \| `VIF` \| `Site_scores` \| `Env_scores` \| `Note`
- **开关**：`constrained_ord$enabled`；`method`(dbrda/rda/cca/capscale)、`response`(community/distance)、`response_vars`、`response_distance`、`constraints`、`transform`(hellinger)、`vif_thr`(15)、`permutations`(999)、`arrow_ratio`(0.78)
- **降级**：`md_cord_fit()` 返回 NULL → SKIP（解释变量不足 / 样本 < 5 / vegan 拟合失败）

### `45_response_curves` — `md_response_curves(CFG)`
- **源**：原实现（§6 GAM 非线性响应，原图）
- **输入**：`md_core_analysis(CFG)`；`y` 任意连续列，`x` 任意列名向量（NULL → `env$vars`）
- **算**：每个 `x` 调一次 `md_gam_curve()`（`mgcv::gam(s(x, k), method="REML")`；`random` 非空时自动升 `gamm`）；只留 `p_smooth <= p_max`；全超阈值则留 P 最小者
- **图**：`45_response_curves.pdf`（183 × 70 mm）：**只有一组 GAM 曲线**（每变量一格，ncol 自适应）。预测变量默认取**随机森林前 `max_panels`（3）个**（对齐原脚本 `RF -> top3 -> GAM`，读不到模块 40/71 的 RF 结果时模块 45 自己跑一次 `md_rf_fit()`）；只画 `p_smooth <= p_max` 的平滑项，一个都没通过时按与 83 相同的回退口径只显示 P 最小的一个并在面板内标注 `p > 0.10 (not significant)`。原实现的 **B 面板（EDF / -log10(P) 条形）是自造的，2026-09-26 已删除**
- **表**：`Curve_fits` \| `Curve_grid` \| `Note`
- **开关**：`response_curves$enabled`；`y`、`x`、`engine`(gam/gamm/lm)、`k`(4)、`random`、`grid_n`(100)、`edf_nonlin`(1.1)、`p_max`(0.1)、`max_panels`(6)、`x_from`（NULL / `"rf_top"`）
- **验收**：**因变量不必是群落**（这是与原文最大的区别）；`random` 给出时 fit 表的 `engine` 为 `"gamm"`

## 组 10 · 网络与类群（80–83，2026 补充轮）

### `80_network` — `md_network(CFG)`（**默认图件清单里唯一出 4 张图的模块**）
- **源**：原实现（合并了 原脚本 + 原脚本）
- **输入**：`md_core_asv(CFG, "composition")`（相对丰度）+ `md_core(CFG)$tax` + `md_core_meta(CFG)`
- **算**：出现率过滤 → Spearman 相关 → `|r|` 阈值 → igraph 建图（**边权保留正负号**）→ 中心性（最大连通组分）→ 样本级指标（含 bootstrap 均衡）→ Louvain 模块 + Zi-Pi 角色 → 韧性模拟 → Erdos-Renyi 随机网络 → 分组子网络 → 丰富/稀有亚网络 → 节点周转
- **图（6 张，编号 80a-80f）**：`80a_network_global.pdf`（183×92，**左右各半**：A 全局网**无坐标轴** ｜ B vs 随机**保留**坐标轴）· `80b_network_by_group.pdf`（183×165，分组网**按 ceiling(sqrt(n)) 排布，6 个 = 2 行 × 3 列**、无坐标轴；下方图例带用 `heights = c(7, 1)` 压到约 12% 高度，**写出每个门名**）· `80c_network_turnover.pdf`（183×150，A 周转 Sankey ｜ B 丰富/稀有亚网络）· `80d_network_nodes.pdf`（183×220，A betweenness ｜ B 度-丰度 ｜ C Zi-Pi；`network$style$wide_layout = TRUE` 回到原文 533.4×177.8 单排）· `80e_network_metrics.pdf`（183×165，6 指标 × 2 因子 = 12 格）· `80f_network_metrics_supp.pdf`（同版式 12 格，`metrics_supp_figure` 控制）。
- **表**：`Node_centrality` \| `Sample_metrics` \| `PERMANOVA` \| `Group_topology` \| `Hub_nodes` \| `Global_topology` \| `Random_network` \| `Abundant_rare` \| `Node_turnover` \| `Note`
- **缓存**：`cache/80_analysis.rds`（分析段，`md_stage`）+ **`cache/80_network.rds`**（跨模块产物：`nodes` / `edges` / `params`）
- **开关**：`network$enabled`；`method`(cor/sparcc)、`cor_method`、`min_cor`(0.5)、`sub_min_cor`、`min_prev`(0.02)、`max_asv`(1500)、`module_algo`、`zi_thr`(2.5)、`pi_thr`(0.62)、`random_n`(100)、`robustness_nsim`(30)、`n_boot`(1000)、`metrics_main`、`metrics_supp`、`metrics_supp_figure`、`group_factors`、`sub_network_by`、`turnover_by`、`style$*`
- **降级**：过滤后 ASV < 3 → SKIP 并写明原因；`method="sparcc"` 但 SpiecEasi 缺失 → 回退 `cor` + notes
- **验收**：`robustness_R()` 之前必须 `md_seed`（修掉原脚本的 seed 错位）；布局在**出图段**算（`style$layout` 不进缓存键）；图上文字纯 ASCII（原文的中文指标名已换英文，见 `ROUND5_MODULE_SPEC.md` §5）

### `81_taxon_group` — `md_taxon_group(CFG)`
- **源**：原实现；**硬编码固氮属表已删除**，改读 `data/05_taxon_group.csv`
- **输入**：可选第 6 文件（`taxon, group`）+ `md_core_asv`；网络角色部分读 `cache/80_network.rds`（**软依赖**，不用会 stop 的访问器）
- **算**：组成（组×样本相对丰度）→ 梯度检验 `md_group_test()` → 网络角色（组×Zi-Pi role）→ 响应变量 Spearman → 跨组边 Fisher → 组指示种 `md_indval()`
- **图**：`81_taxon_group.pdf`（183 × 150 mm，heights `c(1, 1.2)`）：A 组×网络角色堆叠 ｜ B 组组成堆叠（**缺 80 缓存时退化为只有 B**）
- **表**：`Mapping` \| `Composition` \| `Gradient_test` \| `Network_centrality` \| `Response_correlation` \| `Edge_test` \| `Indicator_species` \| `Note`
- **开关**：`taxon_group$enabled`；`file`、`taxon_col`、`group_col`、`match_by`（`ASV`/`Phylum`…/`Species`）、`group_order`、`group_label`、`drop_unmapped`、`indicator$*`、`style$*`
- **降级**：文件不存在 / `file=NULL` / `enabled=FALSE` → SKIP；80 缓存缺失 → 跳过 A 面板 + notes；`response$enabled=FALSE` → 跳过功能关联
- **验收**：未匹配类群归 `Unmapped` 并在 notes 报出未匹配 ASV 数与占比；覆盖率 < 50% 时 `md_warn`

### `82_taxon_share` — `md_taxon_share(CFG)`
- **源**：原实现（原脚本）
- **输入**：`md_core_asv` + `md_core_meta`；分组取 `taxon_share$factor`（NULL → `design$factors[1]`）
- **算**：各水平"至少在 `min_prev` 个样本出现"的类群集合 → 组合交集统计
- **图**：`82_taxon_share.pdf`（183 × 110 mm）：**单面板 UpSet**（原脚本 附图b 就是一张 UpSet）。原 B 面板（Venn / 共享-独有条形）已于 2026-09-26 按用户裁定删除，`show_venn` 键保留但不再生效；单面板不进 patchwork，因此图是**全矢量**的（0 个位图 XObject）
- **表**：`Set_size` \| `Intersection` \| `Sharing_stats` \| `Note`
- **开关**：`taxon_share$enabled`；`factor`、`rank`(ASV)、`min_prev`、`max_comb`(20)、`upset_top_n`(20)、`show_venn`、`style$*`
- **降级**：因子水平 < 2 → SKIP；`ComplexHeatmap` 缺失 → 用条形/geom_tile 兜底
- **验收**：组合数受 `max_comb` 限制，其余组合仍完整写进 `Intersection`；水平 ≥ 4 时不画 Venn
- **⚠ 已知取舍**：Venn/UpSet 面板走 `annotation_raster()` 栅格化后再进 patchwork（理由见模块头注释），
  所以这两个面板**不是矢量**。若期刊硬性要求全矢量，需要另做一版。

### `83_taxon_diff` — `md_taxon_diff(CFG)`
- **源**：原实现（原脚本）
- **输入**：`md_core_asv(CFG, "composition")` + `md_core(CFG)$tax` + `md_core_meta`
- **算**：出现率过滤（`min_prev_frac` 0.05）→ `md_glom()` 到 `rank` → `md_rel()`/`md_clr()` → `md_diff_taxa_kw()`（KW+BH+Cliff's δ+富集组）→ `md_taxon_dunn()`
- **图**：`83_taxon_diff.pdf`（183 × 150 mm）：A 差异类群热图（Top `top_n`，`heatmap_engine` ComplexHeatmap/pheatmap，`heatmap_scale` row_z/rel/clr）｜ B（`show_box`）Top 差异类群箱线
- **表**：`Top_taxa` \| `Kruskal` \| `Dunn` \| `Cliff_delta` \| `Note`
- **开关**：`taxon_diff$enabled`；`rank`(Genus)、`clr`、`pseudocount`、`min_prev_frac`、`min_n`(3)、`test`(kruskal)、`posthoc`(dunn)、`adjust`(BH)、`effect`(cliff)、`p_thr`(0.05)、`top_n`(25)、`heatmap_engine`、`heatmap_scale`、`show_box`、`style$*`
- **降级**：显著类群 < `min_sig`（默认 3）时，热图**回退为按 KW_padj 排序的 Top-N**，副标题写明 `None/only N passed FDR < 0.05; showing top N by KW p`（与模块 45 的回退同一口径）；列数 < 2 或矩阵为空 → SKIP；热图包全缺 → `geom_tile()` 兜底（离散轴，不会再出现 0.00-1.00 的数值轴）
- **验收**：CLR 伪计数缺省 = 最小正值的 1/2（原脚本口径，见 pitfalls C18）

## 组 11 · 已有模块的**可选增量**（默认全关，2026 补充轮）

| 模块 | 新增路径 | 开关 | 默认 | 影响 |
|---|---|---|---|---|
| `10_composition` | 指示种 IndVal 面板 D（+ `Indicator_species` sheet） | `composition$indicator$enabled` | `FALSE` | 图每多一个因子高 70 mm；表多 1 sheet |
| `40_env_summary` | 随机森林重要性填上**设计串里预留的 D 区**（+ `RF_importance` sheet；可选 PDP 面板 E） | `env$rf$enabled` | `FALSE` | 图高度 190 mm 不变；开 PDP 时 245 mm |
| `70_response_pattern` | `pattern_test = "scheirer_ray_hare"`（+ `SRH_two_way` sheet） | `response_models$pattern_test` | `"kruskal_dunn"` | 面板 A 右上角注释改为 SRH 的 P；字母仍来自 KW + Dunn（SRH 不做事后） |
| `71_response_linkage` | `linkage` 追加 `"procrustes"` / `"rf"`（+ `Procrustes` / `Procrustes_summary` / `RF_importance` sheet 与 C/D 面板） | `response_models$linkage` | 不含这两项 | 面板数 2 → 4，图高 95 → 115 mm |

> ⚠ 这四条在默认配置下**完全不执行**，已有图件与数据表逐像素一致（验收见
> `references/04-pitfalls.md` C25）。

---

## 附录 · 实际产物清单（2026-09-26 实核）

> 上面各模块条目里的 sheet 名是**设计时**写的；下表是从实际跑出的文件里读出来的，**以本表为准**。

### 结果表与 sheet

| 文件 | sheet |
|---|---|
| `00_index.xlsx` | Index_tables \| Index_figures \| Validation \| Run_summary \| Skipped \| Session \| `_meta` |
| `01_data_health.xlsx` | Overview \| Files \| Issues \| Taxonomy \| Factors \| Factor_crosstab \| Tree_check \| Sample_reads \| Response_diagnostics \| Env_missing \| `_meta` |
| `05_assemble.xlsx` | Filter_log \| Normalization \| Coverage \| Factor_sizes \| Depth_sensitivity \| `_meta` |
| `10_composition.xlsx` | Taxa_top \| Group_means \| ASV_quadrant \| Quadrant_counts \| Indicator_species \| Notes \| `_meta`（`Indicator_species` 只在 `composition$indicator$enabled = TRUE` 时存在） |
| `11_alpha.xlsx` | Descriptive \| Tests \| Posthoc \| Rarefaction_curve \| Notes \| `_meta` |
| `13_distances.xlsx` | Distance_summary \| Parameters \| `_meta` |
| `14_ordination.xlsx` | Axis_variance \| Ellipse_diagnostics \| Site_scores \| Notes \| `_meta` |
| `15_beta_tests.xlsx` | Effect_sizes \| Pairwise \| PERMDISP \| Note \| `_meta` |
| `16_beta_components.xlsx` | Pair_summary \| By_levelpair \| Abundance \| Mantel \| Note \| `_meta` |
| `17_gamma_partition.xlsx` | Partition \| Coverage \| Note \| `_meta` |
| `20_phylo_prep.xlsx` | Tree_summary \| Tip_alignment \| Taxon_ring \| Note \| `_meta` |
| `21_phylo_diversity.xlsx` | Phylo_indices \| Group_desc \| Group_test \| Note \| `_meta` |
| `22_phylo_beta.xlsx` | Region_pair \| Phylo_means \| Mantel \| Pairwise_rho \| Note \| `_meta` |
| `30_null_models.xlsx` | Process_summary \| BetaNTI_stats \| RCbray_stats \| By_factor1 \| By_factor2 \| Note \| `_meta` |
| `31_assembly_process.xlsx` | Process_overall \| By_factor1 \| By_factor2 \| By_interaction \| BetaNTI_RCbray \| Quadrant \| Data_features \| Note \| `_meta` |
| `33_distance_decay.xlsx` | DDR_fits \| Distance_summary \| Note \| `_meta` |
| `34_variance_partition.xlsx` | **Fractions \| Blocks \| Total** \| Note \| `_meta` |
| `40_env_summary.xlsx` | Descriptive | VIF | Correlations | Cor_matrix | P_adjusted | Soil_PCA | Env_PCA | Group_means | Derived | RF_importance | Note | `_meta`（`RF_importance` 只在 `env$rf$enabled = TRUE` 时存在） |
| ~~42_constrained_ord.xlsx~~ | **已作废**：2026 补充轮以 `44_constrained_ord` 重新入管线（见 `07-spare.md` C 节） |
| `43_mantel_mrm.xlsx` | Mantel \| MRM \| Effect_sizes \| Saturation \| Pair_levels \| Variables \| Note \| `_meta` |
| `50_lmm.xlsx` | Coefficients \| Fit \| Nonlinearity \| Collinearity_guard \| Variables \| Note \| `_meta` |
| `52_sem.xlsx` | Path_coefficients \| Fit \| Fit_Interpretation \| dSeparation \| Variables \| Note \| `_meta` |
| `53_sensitivity.xlsx` | Robustness_summary \| Coefficients \| Model_fit \| Design \| Note \| `_meta` |
| `60_spatial_autocorr.xlsx` | Moran_I \| Correlogram \| Coordinates \| Notes \| `_meta` |
| `70_response_pattern.xlsx` | Descriptive | Tests | Posthoc | Interaction_model | Continuity | SRH_two_way | Notes | `_meta`（`SRH_two_way` 只在 `pattern_test = "scheirer_ray_hare"` 时存在） |
| `71_response_linkage.xlsx` | Alpha_marginal | Alpha_partial | Ordination_axes | Pairwise_MRM | Pairwise_Mantel | Pairwise_effect | Pairwise_saturation | Procrustes | Procrustes_summary | RF_importance | Notes | `_meta` |

| `12_niche_breadth.xlsx` | Niche_Breadth \| Cell_summary \| Sample_strategy \| Threshold_sensitivity \| Sample_stats \| Gradient_test \| Note \| `_meta` |
| `44_constrained_ord.xlsx` | Model_summary \| Terms \| VIF \| Site_scores \| Env_scores \| Note \| `_meta` |
| `45_response_curves.xlsx` | Curve_fits \| Curve_grid \| Note \| `_meta` |
| `80_network.xlsx` | Node_centrality \| Sample_metrics \| PERMANOVA \| Group_topology \| Hub_nodes \| Global_topology \| Random_network \| Abundant_rare \| Node_turnover \| Note \| `_meta`（另含 Kruskal_Wallis \| Bootstrap_metrics \| Jaccard_overlap \| Node_persistence \| ASV_cells \| Module_size） |
| `81_taxon_group.xlsx` | Mapping \| Composition \| Gradient_test \| Network_centrality \| Response_correlation \| Edge_test \| Indicator_species \| Mapping_coverage \| Note \| `_meta` |
| `82_taxon_share.xlsx` | Set_size \| Intersection \| Sharing_stats \| Note \| `_meta` |
| `83_taxon_diff.xlsx` | Top_taxa \| Kruskal \| Dunn \| Cliff_delta \| Note \| `_meta` |
**每个结果表都带一个 `_meta` sheet**（包版本 / R 版本 / seed / run_id / 项目 / 数据 md5 / 用到的包版本），用于回溯。

### 图件（示范数据下 **31 张**；`61_map` 要联网取 DEM，默认关）

```
本包口径（示范数据实产 19 张）：
  10_composition  11_alpha  14_ordination  16_beta_components  17_gamma_partition
  20_phylo_prep  21_phylo_diversity  30_null_models  31_assembly_process
  33_distance_decay  34_variance_partition  40_env_summary  43_mantel_mrm
  50_lmm  52_sem  53_sensitivity  60_spatial_autocorr  70_response_pattern  71_response_linkage

2026 参考实现补充（12 张）：
  12_niche_breadth  44_constrained_ord  45_response_curves
  80a_network_global  80b_network_by_group  80c_network_turnover
  80d_network_nodes  80e_network_metrics  80f_network_metrics_supp
  81_taxon_group  82_taxon_share  83_taxon_diff
```

- `61_map` 默认不出图（`map$enabled = FALSE`，需联网取 DEM）。
- `80` 是唯一出 4 张图的模块；把 `network$metrics_supp_figure` 设为 `FALSE` 可少出 1 张。
- 另有 4 条**可选面板**不占新图号（`10` 指示种 / `40` 随机森林 / `70` SRH 口径 / `71` Procrustes+RF），
  默认关；全部打开时由 `config/config_demo_extras.R` 演示。

### 空 sheet 的行为

`md_write_xlsx()` 会**丢弃零行 sheet** 并把被丢的名字记进 `_meta` 的 `empty_sheets_dropped`，
避免 `95_export_tables` 把"本论文没有这类结果"误报成产物校验 FAIL。
所以某张表少了一个 sheet 不一定是 bug —— 先看它是否真的没有结果。

> ⚠ 一个例外已经改掉：模块 81 的 `Indicator_species` 以前在"一个都没通过 FDR"时
> 整张 sheet 被丢掉，于是"跑过了没结果"与"根本没跑"分不清。现在改为
> **按 IndVal 降序保留 top_n 并加 `fdr_sig` 列**，一个都不过时也写表 + 记 note。

