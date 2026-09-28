# 快速上手 · micro-div-pipeline 扩增子微生物多样性分析代码包

> 本文件原名 `README.md`：本机的 skill 校验器（`skill-creator/scripts/validate_skill.py:80-81`）
> 把 skill 目录根下的 `README.md` 判为 **ERROR**（"put docs in SKILL.md or references/"），
> 所以这份**面向人的介绍**放在 `references/` 里。只是位置变化，内容与用法都不受影响。

**一句话**：把微生物多样性论文的分析流程标准化——整理一次数据，说清科学问题，勾选方法与参数，生成代码，确认后执行，交付结果表与 pdf 图。

> 适用范围：从一篇 nifH 扩增子论文（147 样本 / 23,873 ASV / 4 分区 × 4 植被型）的 20 个 R 脚本
> 通用化而来。**换数据不改代码**：16S、ITS、nifH 及其他 marker 都可用，只要按下面五文件交数据。
>
> **交付物**：`results/*.xlsx`（每模块一张多 sheet 结果表）+ `figures/*.pdf`（每模块一张合成图）
> + `results/00_index.xlsx`（产物索引 + 校验 + 运行汇总 + 会话信息）。

---

## 1. 怎么用（三步）

### 第 1 步：整理数据（只做一次）
在新论文项目里建目录，把数据整理成 **5 个纯文本文件**（详见 `references/00-data-contract.md`）：

```
<论文项目>/
├─ data/
│  ├─ 01_asv_table.csv      ASV × 样本 丰度矩阵（首列 asv_id，其余列 = 样本，原始 reads）
│  ├─ 02_taxonomy.csv       asv_id + Kingdom..Species
│  ├─ 03_sample_meta.csv    首列 sample_id + 任意个分组因子/协变量
│  ├─ 04_env_response.csv   首列 sample_id + 环境变量 + 响应变量（可无响应变量）
│  └─ 99_tree.nwk           Newick 树（tip = asv_id）；没有就写 NULL
├─ config/config.R          从本包 config/config_template.R 复制修改
├─ code/  cache/  results/  figures/  logs/
```

### 第 2 步：改配置
只改 `config/config.R`，不动 `R/` 里的模块代码。要改的通常是：
- `design$factors` / `design$ordered_levels` — 分组因子与其水平顺序
- `design$sampling` / `design$site_var` — 采样设计与站点列（随机效应用）
- `normalize$method` / `normalize$depth` / `normalize$by_module` — 标准化口径
- `response$var` — 哪一列是响应变量（没有就 `enabled = FALSE`）
- `filter$*` — 读长与出现率阈值
- `env$vars` / `env$groups` / `env$vpa_blocks` — 环境变量与其分组
- `model$*` — 随机效应、焦点变量、协变量
- `sem$edges` — SEM 的 DAG（真的是一张边表，不是公式）
- `assembly$runs` — 零模型随机化次数（正式写论文用 999）
- `opts$seed` — 全包唯一种子

`references/02-params.md` 有全部参数的默认值、可选项、**推荐值与理由**。

### 第 3 步：跑
```powershell
# 先自检（只读！不装包、不写 results/，只落 logs/preflight_*.csv）
D:/R-4.5.2/bin/Rscript.exe --vanilla R/90_run_all.R --config=config/config.R --check

# 用示范配置（36 样本 × 327 ASV，随包交付）
D:/R-4.5.2/bin/Rscript.exe --vanilla R/90_run_all.R --config=config/config_demo.R

# 用你自己的配置
D:/R-4.5.2/bin/Rscript.exe --vanilla R/90_run_all.R --config=config/config.R --project=D:/path/to/paper

# 只跑部分模块（PowerShell 里逗号必须加引号）
D:/R-4.5.2/bin/Rscript.exe --vanilla R/90_run_all.R '--only=01,05,10,11,13'

# 只重出图（走缓存，跳过全部重计算）
D:/R-4.5.2/bin/Rscript.exe --vanilla R/90_run_all.R --figures-only

# 忽略缓存全部重算
D:/R-4.5.2/bin/Rscript.exe --vanilla R/90_run_all.R --force
```
结果 → `results/`、图 → `figures/`、日志 → `logs/`。
**任一模块 ERROR 会以非零退出码结束**（不会静默通过）。

`--check` = `SKILL.md` 第 6 步「自检（不执行）」的可执行版本，跑四项：语法 `parse()`、
依赖是否齐（**只报告不安装**，"缺 X → 该模块将 SKIP"）、配置声明的因子/环境变量/响应变量/`sem$edges`
是否都在数据列名里、以及五文件的输入指纹（路径 + 行数 + ASV 数 + md5）。
有 BLOCK（语法错 / 必需列缺失）才返回非零；缺包只算 WARN。可与 `--only=` 组合。

临时快速试跑可用环境变量覆盖：`MD_RUNS=99`、`MD_SEED=1`、`MD_FORCE=1`、`MD_FIGONLY=1`。

---

## 2. 模块清单（5 个共享库 + 32 个分析模块）

| 组 | 编号 | 模块 | 一句话 |
|---|---|---|---|
| 数据 | 01 / 05 | 数据体检 / 装配与标准化 | 五文件体检 + 对齐 + 过滤 + 逐模块标准化 + Good's coverage |
| 组成与多样性 | 10 | 物种组成 | 任意层级组成堆叠柱 + Top-N + 核心/稀有四象限表（+ 可选指示种） |
| | 11 | α 多样性 | 6 个指数 + 组间检验 + CLD 字母 + 稀释/累积曲线 |
| | 12 | 生态位宽度 | 样本策略组成（双因子交互堆叠）+ 四象限散点 + Levins' B 与 gam 平滑（对应毕业论文 Fig1c） |
| | 13 | 距离矩阵 | Bray / Jaccard / UniFrac / wUniFrac，单独存 rds |
| | 14 | 排序 | PCoA / NMDS + 95% 椭圆 + 轴解释率 |
| | 15 | 多元差异检验 | PERMANOVA（`by="terms"`）+ ANOSIM + PERMDISP + MRPP + 两两比较 |
| | 16 | β 周转分解 | turnover / nestedness（betapart）+ 分类 vs 系统发育 β |
| | 17 | α-β-γ 分区 | iNEXT / iNEXT.3D，Hill 数 0/1/2 + 覆盖度 |
| 系统发育 | 20 | 树准备与可视化 | 校验 / 生根 / 子集 + fan 与 rectangular 出图 |
| | 21 | 系统发育多样性 | Faith's PD + MPD/MNTD + NRI/NTI（SES 零模型） |
| | 22 | 系统发育 β | phylo Sørensen 分解 + 与分类学 β 对比 |
| 群落组装 | 30 | 零模型引擎 | βNTI + RCbray（分块 + 并行 + 缓存） |
| | 31 | 组装过程判定 | Stegen 五过程 + 比例汇总 + 四象限散点 |
| | 33 | 距离衰减 | DDR（地理 / 环境 / 群落三口径，log-log） |
| | 34 | 变差分解 | VPA/varpart，任意分块，Hellinger 转换 |
| 环境关联 | 40 | 环境表征 | 描述统计 + VIF + 相关矩阵 + 组均值热图 + 环境 PCA + `Soil_PC1`（+ 可选随机森林） |
| | 43 | Mantel / MRM | 样本级与成对级双口径 + 饱和与膨胀诊断 |
| | 44 | 约束排序 | db-RDA / RDA / CCA / capscale + envfit + VIF（取代早期的 `42`） |
| | 45 | GAM/GAMM 响应曲线 | 预测变量取随机森林前 N 个，逐个拟合平滑并画 95% CI |
| 建模 | 50 | 混合模型 | LMM / GLMM / GAMM + 条件效应曲线 + 共线守卫 + 森林图 |
| | 52 | SEM | 声明式边表 → psem → d-separation → 分层路径图 |
| | 53 | 敏感性分析 | 剔因子水平 / 站点 / 离群点后重跑对比 |
| 空间 | 60 | 空间自相关 | Moran's I + correlogram（站点级/样本级自动守卫） |
| | 61 | 研究区地图 | DEM + 站点 + inset；无网自动降级 |
| 响应变量 | 70 | 响应变量格局 | 双因子交互箱线 + 连续变量散点（可切 Scheirer-Ray-Hare） |
| | 71 | 响应变量关联 | 响应 vs α / 排序轴 / 成对距离（MRM + bootstrap 森林图）（+ 可选 Procrustes / RF） |
| 网络与类群 | 80 | 共现网络 | 全局网 / 分组网 / 节点属性 / 拓扑指标 / 随机网络 / 丰富-稀有亚网络（出 6 张图 `80a`–`80f`） |
| | 81 | 通用类群分组 | 任意分组方案（固氮菌型、生活史策略、自定义表）× 网络角色 / 组成 |
| | 82 | 类群共享 | UpSet 交集统计（单面板；原 Venn 面板已按用户裁定删除） |
| | 83 | 差异类群 | KW + BH + Dunn + Cliff's δ + CLR 热图（显著类群不足时回退 Top-N 并标注） |
| 运行导出 | 90 / 95 | runner / 索引与校验 | manifest 驱动、日志、非零退出、`--check` 只读自检；产物索引 + 后置校验 |

> 曾经的 `42_constrained_ord` 已移出管线，其能力由 **`44_constrained_ord`** 承接
> （按毕业论文 `3_environmental_drivers.R:855-1000` 的 db-RDA 重写）。旧文件归档在包外的
> `_build-notes/micro-div-pipeline/removed-code/`，只作追溯用。

**不做**：宏基因组 / KO-KEGG 通路、差异丰度（DESeq2 / ANCOM-BC / LEfSe）、
功能预测（FAPROTAX / PICRUSt2）、中性群落模型 NCM、iCAMP、从序列建树、docx 三线表。

---

## 3. 六条使用纪律

1. **参数不写死**：所有分析参数集中在 `config/`，模块代码不硬编码；`references/02-params.md` 是参数总表。
2. **缺数据只报告、不动手**：不自动建树、不自动插补、不自动改类群名。
3. **口径不锁死**：每个模块保留可选方法，由用户按数据与科学问题拍板（包会给有依据的推荐）。
4. **结果可回溯**：每个统计量都能指回结果表的具体 sheet 与行；每个结果表都带 `_meta` sheet（版本 + seed + 时间）。
5. **输出精简**：每模块一张图 + 一张表（只出核心表）；距离矩阵单独存 rds；不出汇总总表。
6. **数据只读**：`data/` 永不被改写；`01` 记录 md5，`95` 交付前比对。

---

## 4. 装包

核心包（缺了跑不动）：`ggplot2` `patchwork` `dplyr` `tidyr` `tibble` `purrr` `stringr` `openxlsx` `digest`。

按模块需要的包组（缺哪组，对应模块**自动跳过并打印安装命令**，不影响其它模块）：

```r
# 生态与系统发育（13-22, 33, 34, 40-43, 60, 70, 71）
install.packages(c("vegan", "ape", "picante", "betapart", "GUniFrac", "phangorn",
                   "ecodist", "geosphere", "spdep"))
# 分区（17）
install.packages(c("iNEXT", "iNEXT.3D"))
# 统计检验（11, 15, 70）
install.packages(c("dunn.test", "rstatix", "multcompView", "FSA", "car"))
# 建模与 SEM（50, 52, 53, 71）
install.packages(c("lme4", "lmerTest", "performance", "mgcv", "piecewiseSEM"))
# 环境表征（40）
install.packages(c("FactoMineR", "factoextra"))
# 绘图（全部出图模块）
install.packages(c("cowplot", "scales", "RColorBrewer", "ggrepel"))
# 地图（61）
install.packages(c("sf", "terra", "rnaturalearth", "ggspatial"))
# 并行（21, 30）
install.packages(c("future", "future.apply"))
# 可选：CSS 标准化
BiocManager::install("metagenomeSeq")
```

已在本机（R 4.5.2）验证过的版本可直接跑；`95_export_tables` 会把实际用到的包版本写进索引表。

---

## 5. 文档索引

| 文件 | 内容 |
|---|---|
| `SKILL.md` | 给 AI 的工作流与提问协议（11 步 + 七道硬门） |
| `references/00-data-contract.md` | 五文件输入契约与命名规范 |
| `references/01-method-menu.md` | 科学问题 → 模块；模块 → 方法选项与推荐规则 |
| `references/02-params.md` | **全参数总表**（参数 / 默认 / 可选项 / 推荐 / 理由 / 何时要改） |
| `references/03-module-specs.md` | 32 个模块的签名、输入输出、消费/产出的缓存、验收要点 |
| `references/04-pitfalls.md` | 踩坑库（含从 nifH 论文移植的工程坑 + 本包实测） |
| `references/06-figure-map.md` | 哪张图来自哪个模块、对应论文哪张图；视觉口径的"唯一真身"对照表 |
| `references/07-spare.md` | 默认关着的备用能力有哪些、怎么开 |
| `references/08-maintenance.md` | **改这个包本身**的流程：`--check` 自检、回归基线与比对、验收硬线、禁止项 |
| `config/config_template.R` | 配置模板（全参数、每项带推荐注释） |
| `config/config_demo.R` | 示范数据配置 |
| `assets/make_demo_data.R` | 示范数据生成脚本（从 nifH 论文数据只读抽取，固定 seed） |
| `assets/project_skeleton/` | 新论文项目骨架 |
| `assets/compare_figures.R` | 两张图逐像素比对（改包回归用） |
| `scripts/audit_figures.R` | 图件体检：非 ASCII 审计 + 位图 XObject 计数 + 页尺寸 |
| `scripts/snapshot_results.R` | 结果表逐 sheet 元素指纹（改包回归用） |
| `run/manifest.R` | 模块清单（依赖、产物、是否重计算）—— 唯一真身 |
| `run/_TEMPLATE_run_all.R` | 调用脚本模板 |
| `data/` | 示范数据（**模板专用，不得用于发表**） |

---

## 6. 示范数据是什么

`data/` 里是 36 样本 × 327 ASV 的示范子集，由 `assets/make_demo_data.R` 从一篇已定稿的
nifH 扩增子论文数据中**只读抽取**而来（固定 seed = 123）：

- **样本**：3 分区（Coastal / Nearshore / Xisha）× 3 植被型（Bare / Herb / Tree）× 4 重复 = 36，
  跨 11 个站点（`Site`）
- **ASV**：327 个（每个样本丰度前 10 个 ASV 的并集），抽平深度 9968
- **响应变量**：`EF`（生态系统功能，原论文为固氮酶活性），已做术语中性化
- **环境**：`MAT MAP SWC NH4 NO3 pH OM TN TP EC TK Ca Fe` + `Latitude Longitude`

> ⚠ `data/` 只是**模板与自检用**，不得用于发表。要用自己的数据，请按 `references/00-data-contract.md` 组织。
>
> ⚠ 这批数据**极端单例化**（原数据 23,251 / 23,873 个 ASV 在 147 样本中只出现 1 次），
> 所以示范结果里"核心类群几乎为 0、β 周转极高、漂变占绝对主导、环境解释力只有百分之几"
> 都是**数据的真实特征**（与那篇论文的结论一致），不是程序错误。
