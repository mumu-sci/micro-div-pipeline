# micro-div-pipeline

**扩增子微生物多样性论文的分析代码包** · 16S / ITS / nifH / 其他 marker 通用

> 整理一次数据 → 改一个 `config` → 一条命令跑完 30+ 个分析模块 → 拿到 `results/*.xlsx` 与 `figures/*.pdf`。
> 换数据不改代码：`R/` 里任何模块都不需要动。

---

## 它是什么

由一整套扩增子分析流程沉淀为模块化包。

| 维度 | 做法 |
|---|---|
| 数据无关 | 只改 `config/config.R`，不动 `R/` 模块代码 |
| 交付形态 | 每模块一张 patchwork 合成图（`pdf`，期刊物理尺寸 90 / 183 mm）+ 一张多 sheet `xlsx` |
| 产物索引 | `results/00_index.xlsx`：有哪些表、哪些图、校验结果、运行汇总、会话信息 |
| 缺数据 | 只报告不动手（不建树、不插补、不改名） |
| 缺 R 包 | `CFG$opts$auto_install = TRUE` 时自动装，装不上才让该模块 SKIP |
| 可复现 | 全包唯一种子 `CFG$opts$seed`；统计量只搬运不发明 |
| 缓存 | 模块级 `.rds`；样式参数不进缓存键 → 改配色/字号可秒级重出图 |

## 三步上手

### 1. 准备数据（5 个纯文本文件）

```
<论文项目>/
├─ data/
│  ├─ 01_asv_table.csv     ASV × 样本 丰度矩阵（首列 asv_id，其余列 = 样本，原始 reads）
│  ├─ 02_taxonomy.csv      asv_id + Kingdom..Species
│  ├─ 03_sample_meta.csv   首列 sample_id + 任意个分组因子 / 协变量
│  ├─ 04_env_response.csv  首列 sample_id + 环境变量 + 响应变量（可以没有响应变量）
│  └─ 99_tree.nwk          Newick 树（tip = asv_id）；没有就写 NULL
├─ config/config.R         从 config/config_template.R 复制修改
└─ code/  cache/  results/  figures/  logs/
```

### 2. 改配置

只改 `config/config.R`。常改的几组：`design$factors`（分组因子）、`normalize$*`（标准化口径）、`response$var`（响应变量）、`filter$*`（读长 / 出现率阈值）、`env$vpa_blocks`（变差分解分块）、`model$*`（随机效应与焦点变量）、`sem$edges`（SEM 的 DAG，是一张边表而不是公式）、`assembly$runs`（零模型随机化次数，论文正式用 999）、`opts$seed`。

### 3. 跑

```powershell
# 只自检（只读：不装包、不写 results/，只落 logs/）
Rscript --vanilla R/90_run_all.R --config=config/config.R --check

# 跑随包示范数据（36 样本 × 327 ASV）
Rscript --vanilla R/90_run_all.R --config=config/config_demo.R

# 跑自己的项目
Rscript --vanilla R/90_run_all.R --config=config/config.R --project=D:/path/to/paper

# 只跑部分模块（PowerShell 下逗号必须加引号）
Rscript --vanilla R/90_run_all.R --config=config/config.R --only="11_alpha,15_beta_tests"
```

## 模块一览

组号沿用 `R/` 的文件名前缀，模块详细规格见 `references/03-module-specs.md`。

| 组 | 模块 | 做什么 |
|---|---|---|
| 数据 | `01_data_health` · `05_assemble` | 五文件体检（维度、读长、ID 交集、分类完整度、树覆盖、设计平衡、响应变量正态性），BLOCK 级问题直接停；装配与标准化各模块统一口径的 ASV 表 |
| 组成与多样性 | `10_composition` · `11_alpha` · `13_distances` · `14_ordination` · `15_beta_tests` · `16_beta_components` · `17_gamma_partition` | 组成堆叠柱 + 核心/稀有类群四象限；α 指数与组间检验 + CLD 字母 + 稀释曲线；Bray/Jaccard/UniFrac；PCoA/NMDS + 95% 椭圆；PERMANOVA + ANOSIM + PERMDISP；β 周转/嵌套分解；iNEXT.3D 的 α/β/γ 分区 |
| 系统发育 | `20_phylo_prep` · `21_phylo_diversity` · `22_phylo_beta` | tip 对齐 / 生根 / 多分支检测；Faith's PD + MPD/MNTD + NRI/NTI；系统发育 β 分解并与分类学 β 对比 |
| 群落组装 | `30_null_models` · `31_assembly_process` · `33_distance_decay` · `34_variance_partition` | βNTI 与 RCbray；Stegen 五过程判定 + 比例汇总 + 四象限；地理/环境/群落三口径距离衰减；varpart 变差分解 + 偏 RDA 显著性 |
| 环境关联 | `40_env_summary` · `43_mantel_mrm` · `44_constrained_ord` | 描述统计 + VIF + 相关矩阵 + 环境 PCA；Mantel / partial Mantel / MRM + bootstrap CI + 膨胀诊断；db-RDA / RDA / CCA |
| 建模 | `50_lmm` · `52_sem` · `53_sensitivity` | LMM（Satterthwaite）+ GAMM（edf 判非线性）+ R²；声明式边表 → piecewiseSEM → d-separation 检验 → 路径系数；剔组/剔站重跑的敏感性分析 |
| 空间 | `60_spatial_autocorr` · `61_map` | Moran's I + correlogram（sample / site / auto 三档）；DEM + 站点 + 放大图（无网自动降级为散点图） |
| 响应变量 | `45_response_curves` · `70_response_pattern` · `71_response_linkage` | GAM 响应曲线；双因子交互箱线（KW + Dunn + CLD）与回归；响应变量与群落指标的 linkage 模型 |
| 网络与类群 | `80_network` · `81_taxon_group` · `82_taxon_share` · `83_taxon_diff` | 共现网络（拓扑 / 中心性 / Zi-Pi / 随机网络 / 丰富-稀有亚网络，唯一出 4 张图的模块）；类群分组；类群共享 UpSet；差异类群热图与指示种 |
| 生态位 | `12_niche_breadth` | Levins' B 生态位宽度 |

## 产物

```
results/
├─ 00_index.xlsx        产物索引（Index_tables / Index_figures / 校验 / 运行汇总 / sessionInfo）
├─ 01_data_health.xlsx  每模块一张多 sheet 结果表
└─ ...
figures/
├─ 11_alpha.pdf         每模块一张 patchwork 合成图（单栏 90 mm / 双栏 183 mm，pt 绝对字号，cairo_pdf）
└─ ...
logs/                   运行日志与 preflight 自检记录
cache/                  模块级 .rds 缓存
```

## 运行环境

- R 4.x（开发与交付用 4.5.2），Windows / Linux 均可
- 主要依赖按模块按需加载：`vegan`、`picante`、`betapart`、`iNEXT` / `iNEXT.3D`、`lmerTest`、`mgcv`、`performance`、`piecewiseSEM`、`betapart`、`randomForest`、`terra` / `raster`（DEM）、`patchwork` / `ggplot2`；`SpiecEasi` 需从 GitHub 安装
- `--check` 只报告不安装；正式运行缺包时按 `CFG$opts$auto_install` 自动装，装不上该模块 SKIP

## 目录结构

```
R/              分析模块（文件名即模块号）+ R/shared/ 五个共享库
                （00_setup / 10_io / 11_stats / 12_plot_theme / 13_plot_types）
config/         config_template.R（全参数与默认值）、config_demo*.R（随包示范）
references/     数据契约、方法菜单、参数表、模块规格、踩坑、快速上手、图件地图、维护
data/           示范数据（36 样本 × 327 ASV + 分类 + 元数据 + 环境 + 树）
assets/         示范数据生成脚本与项目骨架
run/            manifest.R 与 run_all 模板
scripts/        图件审计、结果快照等辅助脚本
```

## 文档索引

| 文件 | 内容 |
|---|---|
| `SKILL.md` | 包定位、11 步工作流、验收口径、常见问题 |
| `references/00-data-contract.md` | 五个输入文件的字段契约与命名规则 |
| `references/01-method-menu.md` | 每个科学问题可选的统计方法与其前提 |
| `references/02-params.md` | 全部参数的默认值、可选项、推荐值与理由 |
| `references/03-module-specs.md` | 逐模块的输入口径、产出、降级与验收 |
| `references/04-pitfalls.md` | 实测踩坑与规避 |
| `references/05-quickstart.md` | 上手流程与命令 |
| `references/06-figure-map.md` | 模块 → 图件清单对照 |

## 边界

- 只接收现成 Newick 树，**不生成树**；系统发育树构建见 [`amplicon-phylo-tree`](https://github.com/mumu-sci/amplicon-phylo-tree)
- 不覆盖宏基因组 / 功能预测（FAPROTAX、PICRUSt2）/ 中性群落模型 / iCAMP
- 姊妹包：[`plant-div-pipeline`](https://github.com/mumu-sci/plant-div-pipeline)（植物多维度多样性，6 组 23 个模块）