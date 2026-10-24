# 06 · 图件清单与对照表（按模块）

> 单位 = **模块**。论文的图常是跨模块拼的（例如 Fig.1 = 系统发育树 + 物种组成 + α 多样性），
> 本包按模块拆开出图。
>
> **两份来源**：
> 1. 自有项目 —— 全套图件，全部有归属模块；
> 2. 2026 参考实现（`本地工作目录 参考实现\13-数据分析\03-生态学数据分析`）——
>    补进 7 个新模块 + 4 条可选路径，共新增 10 张图（见第 1 节末尾）。

## 1. 出图清单（31 张）

### 1.1 自有项目（全套图件）

| # | 图件 | 模块 | 面板构成 | 对应论文的 |
|---|---|---|---|---|
| 1 | `10_composition.pdf` | 物种组成 | A 门组成 ｜ B 属组成（x 轴 = **设计因子交互 交互组**，与箱线图一致）；打开 `composition$indicator$enabled` 后追加 C = 指示种 IndVal。⚠ 原「丰度×出现率四象限」面板已于 2026-09-26 **移入 `12_niche_breadth`**（对应原图 | Fig.1B、Fig.S3 |
| 2 | `11_alpha.pdf` | α 多样性 | A–D 各 α 指数箱线（Shannon / Pielou / Richness / …）｜ E 稀释曲线 | Fig.1D、Fig.S4A、Fig.S4B、Fig.S2 |
| 3 | `20_phylo_prep.pdf` | 系统发育树 | 树（fan 布局 + 分类色环） | Fig.1A |
| 4 | `21_phylo_diversity.pdf` | 系统发育多样性 | A Faith's PD ｜ B NRI ｜ C NTI | Fig.1E、Fig.S4C、Fig.3A |
| 5 | `14_ordination.pdf` | 排序 | A PCoA(Bray) ｜ B PCoA(wUniFrac) | Fig.2B、Fig.2C |
| 6 | `16_beta_components.pdf` | β 周转分解 | 单面板：按因子水平对的周转/嵌套分解（竖向短标签 + 45°，底部两行图例） | Fig.2D |
| 7 | `17_gamma_partition.pdf` | α-β-γ 分区 | 单面板：γ 分解（Hill 刻度；原脚本 §E 只有这一格，**没有覆盖度面板**） | Fig.2E |
| 8 | `30_null_models.pdf` | 零模型 | A βNTI 分布 ｜ B RCbray 分布 | Fig.3C、Fig.S5A、Fig.S5B |
| 9 | `31_assembly_process.pdf` | 组装过程 | A 过程比例 ｜ B 随机性 vs β ｜ C βNTI×RCbray 四象限 | Fig.3B、Fig.3D、Fig.S5C |
| 10 | `34_variance_partition.pdf` | 变差分解 | A 各块唯一解释力 ｜ B **三圆 Venn**（`draw_venn3` 逐行移植；负交叠项按原脚本截 0） | Fig.3E |
| 11 | `33_distance_decay.pdf` | 距离衰减 | A/B/C 三格 log-log DDR（**原脚本无斜率森林图**，已移除） | Fig.S6 |
| 12 | `40_env_summary.pdf` | 环境表征 | A 土壤 PCA 双标图 ｜ B 组均值热图 ｜ C 相关矩阵 ｜ **D（可选）= 随机森林重要性**（`env$rf$enabled`，填上设计串里预留的 D 区） | Fig.S7、Fig.S8A、Fig.S8B |
| 13 | `50_lmm.pdf` | 混合模型 | A 系数森林图 ｜ B–E 条件效应（MAT/Vegetation/Shannon/PD） | Fig.4C、Fig.4D、Fig.5B–E |
| 14 | `52_sem.pdf` | SEM | 路径图（主模型 / 系统发育版） | Fig.4A、Fig.4B |
| 15 | `53_sensitivity.pdf` | 敏感性 | 剔组后系数并列 | Fig.S11 |
| 16 | `61_map.pdf` | 地图 | 研究区主图 + 放大图 | Fig.S1 |
| 17 | `70_response_pattern.pdf` | 响应变量格局 | 双因子交互箱线 + 连续变量散点（`pattern_test` 可切 Scheirer-Ray-Hare，图右上角注释随之改变） | Fig.5A |
| 18 | `71_response_linkage.pdf` | 响应变量关联 | A 响应 vs α（单回归+偏回归）｜ B PCoA1 条件效应 ｜ **C / D（可选）= Procrustes 残差 / 群落→响应 RF**（`linkage` 含 `"procrustes"` / `"rf"`） | Fig.S9、Fig.S12A |
| 19 | `43_mantel_mrm.pdf` | Mantel / MRM | 单面板：\*\*hexbin\*\* \|Δ响应\| vs 群落距离 + MRM/Mantel/饱和度多行标注（**原脚本已把效应量森林图从图里移除**，数字只在 `Effect_sizes` 表） | Fig.S12B |
| 20 | `60_spatial_autocorr.pdf` | 空间自相关 | A Moran's I ｜ B correlogram | Fig.S13 |

> 1.1 表里是 20 行（`61_map` 需联网取 DEM，demo 默认关 → 实际产出 19 张），编号沿用历史顺序。

### 1.2 2026 参考实现补充（12 张；demo 配置下 19 → 31）

| # | 图件 | 模块 | 面板构成 | 对应参考实现的 |
|---|---|---|---|---|
| 21 | `12_niche_breadth.pdf` | 生态位宽度 | **A** 样本生态策略组成堆叠柱（x = **前两个设计因子的交互组**）｜**B** 丰度×出现率四象限散点（**从 10_composition 移入**）｜**C** Levins' B vs 出现率 + gam 平滑 + 95% CI。版式与尺寸（200×150 mm）逐项对齐原图 | 原图 |
| 22 | `44_constrained_ord.pdf` | 约束排序 | 单面板 triplot：样方点（按因子着色/定形状）+ env 箭头 + ggrepel 标签；副标题 `R2adj` 与 P | 原图（db-RDA） |
| 23 | `45_response_curves.pdf` | GAM 响应曲线 | 一排 GAM 曲线（预测变量 = **随机森林前 3 个**，与原脚本 `RF → top3 → GAM` 同流水线）；只画 `p_smooth <= p_max` 的平滑项 | 原图 |
| 24 | `80a_network_global.pdf` | 共现网络·全局 | **左右各占一半版幅**：A 全局网络（按门着色 + 正负边，**无坐标轴**）｜ B 观测 vs 随机网络（柱 + 随机均值±SD，**保留**坐标轴）；183×92 mm | 原图 |
| 25 | `80b_network_by_group.pdf` | 共现网络·分组 | 分组子网络按 `ceiling(sqrt(n))` 自动排布（**6 个网络 = 2 行 × 3 列**），每格标题 `组名 (n=, N=, E=)`，**全部无坐标轴**；下方一条**窄**图例带（`heights = c(7, 1)`，只占约 12% 高度 —— 门色 + **每个门写出名字** + 相关方向 Positive/Negative）；183×165 mm | 原图 |
| 26 | `80c_network_turnover.pdf` | 共现网络·周转 | A 节点周转 Sankey ｜ B 丰富/稀有亚网络指标；183×150 mm | 原图、原图 |
| 27 | `80d_network_nodes.pdf` | 共现网络·节点 | A betweenness top-N ｜ B 度-丰度 ｜ C Zi-Pi 四区（**183 mm 三格纵向**，`network$style$wide_layout = TRUE` 可回到原文的 533×178 mm 单排） | 原图 |
| 28 | `80e_network_metrics.pdf` | 共现网络·指标主图 | `metrics_main` **6 指标 × 设计因子（前 2 个）= 12 格** | 原图 |
| 29 | `80f_network_metrics_supp.pdf` | 共现网络·指标附图 | `metrics_supp` 另外 6 指标，同版式 12 格（`network$metrics_supp_figure` 可关） | 附图 |
| 30 | `81_taxon_group.pdf` | 通用类群分组 | A 分组 × 网络角色堆叠（需 80 的缓存，缺则退化） ｜ B 分组组成堆叠 | 原图 |
| 31 | `82_taxon_share.pdf` | 类群共享 | **单面板 UpSet**（原脚本 附图b 就是一张 UpSet）。原 B 面板（Venn / 共享条形）已于 2026-09-26 按用户裁定删除 | 附图b |
| 32 | `83_taxon_diff.pdf` | 差异类群 | A 差异类群热图（CLR + KW + Dunn + Cliff's δ）；显著类群 < `taxon_diff$min_sig`（默认 3）时回退为 Top-N 并在副标题写明 ｜ B（可选）Top 差异类群箱线 | 附图c |

**总数 = 19（demo 实产） + 12 = 31 张**（把 `network$metrics_supp_figure` 设为 `FALSE` 则为 30 张）。

另有 4 条**可选面板**不新增图号，只往已有图里加一格（默认关）：
`10` 指示种 C ｜ `40` 随机森林 D ｜ `70` Scheirer-Ray-Hare 口径 ｜ `71` Procrustes C / RF D。

## 2. 不出图（只出结果表）—— 你论文里没有对应图

由 `CFG$figure$skip_modules = c("01","05","13","15","22")` 控制。

| 模块 | 它会算的分析 | 你论文里的对应物 |
|---|---|---|
| `01_data_health` | 五文件体检（读长/缺失/因子交叉表/响应变量诊断） | 无（数据质检） |
| `05_assemble` | 过滤 + 标准化 + Good's coverage | 无（方法学中间产物） |
| `13_distances` | Bray / Jaccard / UniFrac / wUniFrac 距离矩阵 | 无（距离只在 PCoA / PERMANOVA 里用） |
| `15_beta_tests` | PERMANOVA / ANOSIM / PERMDISP | Table S4、Table 01（**只有表**） |
| `22_phylo_beta` | 系统发育 β 分解 | Table S11（**只有表**） |

## 3. 已移出管线 / 已被取代的模块

| 模块 | 内容 | 现状 |
|---|---|---|
| `42_constrained_ord` | RDA / CCA / db-RDA + envfit | **已被 `44_constrained_ord` 取代**（按原实现 的 db-RDA 样式重写，不再是早期的自编版）。`44` 与它的差别：支持 `community` / `distance` 两种响应侧，`method` 可选 `dbrda` / `rda` / `cca` / `capscale`，VIF 阈值跟原脚本用 15，箭头长度按站点云半径自动归一。旧文件归档在**包外**构建记录 `analysis-packages/_build-notes/micro-div-pipeline/removed-code/`，仅作追溯（随本包单独分发时不存在）。 |

## 4. 视觉口径来源（**唯一真身**）

所有图的外观只由 `R/shared/12_plot_theme.R` 决定，且**逐项对齐原实现**：

| 项 | 原代码 | 本包 |
|---|---|---|
| 主题 | `theme_isme()` | `theme_md()` |
| 字号（主图） | `S_TICK/S_TTL/S_LEG/S_ANN/S_TAG` | `MD_PT` |
| 字号（拼图） | `S_*_COMP`（显式值，不是统一比例） | `MD_PT_COMP` |
| 散点 / 抖动点 / 线宽 | `S_PT`(0.65) / `PT_SIZE_STD`(0.9) / `PT_ALPHA_STD`(0.4) / `JIT_SIZE_STD`(1.8) / `JIT_ALPHA_STD`(0.67) / `JIT_WIDTH_STD`(0.15) / `LW_AXIS`(0.35) / `BOX_LW_STD`(0.4) | `MD_PT_SIZE` / `MD_LW` |
| 回归线与置信带 | `REG_LINE_COL`(#1E72B0) / `REG_LINE_LW`(0.4) / `REG_RIBBON_A`(0.1) | `MD_REG` |
| Zone 配色 | `cols_zone` | `MD_BUILTIN_PAL$zone` |
| Vegetation 配色 | `cols_veg` | `MD_BUILTIN_PAL$veg` |
| 组装过程配色 | `cols_assembly` | `MD_BUILTIN_PAL$assembly` |
| 核心/稀有配色 | `cols_2x2` | `MD_BUILTIN_PAL$quadrant` |
| 箱线图 | `plot_dual_factor_panel()` | `md_box_two_factor()` / `md_box_group()`（走 `CFG$plot$box`） |

**改外观只改这一处**；`CFG$figure$palette` 与 `CFG$plot$box` 可逐项覆盖。

**图例默认位置 = 图下方、横向排列**（`theme_md()` / `theme_md_comp()` 都设了 `legend.position = "bottom"`）。
`md_guide_flex()` 的默认 `position` 也是 `"bottom"`。
需要放右边或隐藏的模块显式覆盖即可（例如 `70` 用 `md_compose(legend = "none")`）。

## 5. 备用件（默认全关）

**备用件已单独成篇 → `references/07-spare.md`。** 概要：

| 类别 | 内容 | 默认 |
|---|---|---|
| 备用**面板** | A1 样本级指标 vs 群落周转（你论文 Fig.3D 的**通用化**版本，在 `31`）；A2 焦点变量的标准化系数森林图（你论文 Fig.S10 右半的**通用化**版本，在 `50`） | `CFG$spare$*$enabled = FALSE` |
| 备用**口径** | `normalize` 的 CSS / 多次亚抽样、`ordination` 的 NMDS、`beta_tests` 的 MRPP、`assembly` 的 `subset_top`/`chunk` | 各自段里的默认值**等于你论文的口径** |
| 已被取代 | `42_constrained_ord`（RDA/CCA/db-RDA） | 能力由 **`44_constrained_ord`** 承接（见 §3）；旧文件归档在**包外** `analysis-packages/_build-notes/micro-div-pipeline/removed-code/` |

**备用件开启后只允许"多一个面板/多一张图"，不得改变任何已有图件与数值。**
