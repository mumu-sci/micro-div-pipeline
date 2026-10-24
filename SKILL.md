---
name: micro-div-pipeline
description: 扩增子微生物多样性论文的标准分析代码包（16S / ITS / nifH 等 marker 通用），覆盖装配与标准化、物种组成、α/β 多样性、排序与多元检验、β 周转分解、α-β-γ 分区、系统发育（PD / MPD-MNTD / NRI-NTI）、群落组装零模型（βNTI / RCbray）、距离衰减、变差分解、环境表征、约束排序（db-RDA / RDA / CCA）、GAM/GAMM 响应曲线、随机森林、Mantel/MRM、LMM/SEM、Scheirer-Ray-Hare、Procrustes、空间自相关、研究区地图、共现网络（拓扑 / 中心性 / Zi-Pi / 随机网络 / 丰富-稀有亚网络）、类群分组、指示种 IndVal、类群共享 UpSet、差异类群热图、生态位宽度 Levins' B、响应变量关系；32 个模块 + 5 个库 + 全参数 config。Use for 微生物多样性论文的数据分析、16S/ITS/nifH 扩增子分析、ASV 表分析、α 多样性、PCoA / NMDS / PERMANOVA、Faith's PD / NRI-NTI、βNTI / RCbray 群落组装、距离衰减、RDA / db-RDA / Mantel / MRM、SEM 结构方程、LMM / GAMM、随机森林重要性、共现网络 / 关键类群 / Zi-Pi、指示种 IndVal、UpSet / 韦恩图、差异类群热图、生态位宽度，或为新论文搭标准化分析代码。Do NOT use for 宏基因组 / KO-KEGG 通路、功能预测（FAPROTAX / PICRUSt2）、中性群落模型 NCM、iCAMP、从序列建树（只接收现成 Newick）、docx 三线表排版；也不用于「只把已有图美化一遍」或接单交付包装。
---

# micro-div-pipeline · 扩增子微生物多样性分析代码包

> 一次整理数据 → 说清科学问题 → 勾选方法与参数 → 生成代码 → 自检 → **用户确认** → 执行 → 交付结果表与 pdf 图。

## 包的定位

| 项 | 设定 |
|---|---|
| 存放 | `<包目录>`（**不接入自动加载**，由用户在对话中指定使用） |
| 数据无关 | 换论文只改 `config/`，不动 `R/` 里任何模块代码 |
| 覆盖范围 | 32 个分析模块 + 5 个共享库，由一整套扩增子分析流程沉淀而来 |
| 交付形态 | 生成代码 → **自检** → **给用户确认** → 才执行 → 交付 `results/*.xlsx` + `figures/*.pdf` |
| 代码口径 | **不锁死**：每个模块保留可选方法，由用户拍板，我方须给**有科学依据的推荐** |
| 缺数据 | **只报告，不动手**（不建树、不插补、不改名） |
| 缺包 | **只报告前先自动装**：`CFG$opts$auto_install = TRUE`（默认）时 `md_need()` 会 `install.packages()` 后重试；装不上才让该模块 SKIP |
| seed | **全包唯一**：`CFG$opts$seed`（示范配置为 123）。任何模块不许出现字面量种子 |
| 输出 | 每模块**至少一张 patchwork 合成图（pdf）** + **一张多 sheet xlsx**（只出核心表）；子面板单独落盘默认关 |
| 图规格 | 期刊物理尺寸（单栏 90mm / 双栏 183mm）+ **pt 绝对字号** + cairo_pdf |
| 缓存 | 模块级 `.rds` 缓存（`cache/`）；**绘图样式参数不进缓存键** → 改配色/字号秒级重出图 |
| 署名口径 | 统计量只搬运不发明；口径变化必须写进表注或 `_meta` sheet |

---

## 工作流（第 0–10 步，共 11 步）

### 第 0 步 · 定位项目
确认论文项目目录。若无 `data/` 五文件，进入数据整理（见 `references/00-data-contract.md`）。

### 第 1 步 · 数据体检（`R/01_data_health.R`）
读 `config/config.R` → 输出 `results/01_data_health.xlsx`：
- 五文件维度/列名/读长分布/零 reads 样本
- `sample_id` 与丰度表列名的交集与差集
- taxonomy 各层级完整度（未分类比例）
- 树 tip 覆盖率与生根状态
- 每因子水平数 + 每组合样本量 + **因子交叉表（自动发现缺格）**
- 响应变量诊断：Shapiro-Wilk、Bartlett/Levene、偏度、离群点

**体检不通过（BLOCK 级）→ 停止，报告缺口，等用户处置。**
例外：`CFG$opts$strict = FALSE` 时只警告继续。

### 第 2 步 · 科学问题录入
用户以**科学问题**为入口（"分区和植被如何影响群落组成"），而不是直接点模块。
按 `references/01-method-menu.md` 的映射表，输出**推荐模块组合**。

### 第 3 步 · 一次性方法勾选
把「模块 × 可选方法 + 图形/设计选项」列成**一张表**给用户一次勾完（**一次问完，不挤牙膏**）。
**每一行必须带推荐标记 + 一句推荐理由（依据来自第 1 步的诊断结果）。** 禁止只罗列不推荐。

### 第 4 步 · 参数逐项协商（**硬门**）
按 `references/02-params.md` 的全参数总表，逐项请用户拍板，落进 `config/config.R`。
**每一项都必须给出推荐值 + 理由**（用户可能不懂怎么选）。

### 第 5 步 · 生成配置
从 `config/config_template.R` 复制，只改本文件，不动 `R/` 里的模块代码。

### 第 6 步 · 自检（不执行）
```
Rscript --vanilla R/90_run_all.R --config=config/config.R --check
```
一条命令跑完四项，**只读**（不装包、不跑分析、不写 `results/`），报告落 `logs/preflight_<run_id>.csv`：
- 语法检查（`parse()` 全部通过）
- 依赖检查（缺包 → 标注"该模块将降级跳过" + 给安装命令；`--check` 下**只报告不安装**）
- 列名检查（`config` 声明的因子 / 环境变量 / 响应变量是否真实存在）
- 输入指纹（文件路径 + 行数 + ASV 数 + md5）

**有 BLOCK（语法错 / 数据读不进来 / 必需列缺失）必须停下来改配置或数据**，不许带着 BLOCK 进入第 7 步。
缺包只算 WARN，可以继续（对应模块会 SKIP）。

### 第 7 步 · 交用户确认
把「将执行什么 + 预期产出什么」摆出来，等确认。

**什么算确认**：明确的开始执行语 —— 「按计划执行」「开始」「确认执行」「跑吧」「执行」。
**不算确认**：「我看看」「有道理」「先这样」「嗯」「都行」—— **继续停在第 7 步**，不要自行解释成已确认。

### 第 8 步 · 执行
```
Rscript --vanilla R/90_run_all.R --config=config/config.R
```
逐步运行，日志写 `logs/`。**有任一模块 ERROR 则整体非零退出。**
长任务（置换 ≥999、零模型、iNEXT 外推）**执行前先给替代方案与时间预估**。

### 第 9 步 · 交付
- `results/NN_模块_*.xlsx`（每模块一个，多 sheet）
- `figures/NN_模块_*.pdf`（每模块一张合成图）
- `results/00_index.xlsx`（产物索引 + 校验 + 运行汇总 + 会话信息）
- 结果摘要（关键统计量 + 与预期/既有文献是否一致 + 哪些模块降级跳过及原因）

### 第 10 步 · 复核
对每个统计量标注其在结果表中的位置，供用户回溯。

---

## 七道硬门（`config` 里为空 / 缺失就停下问用户，**不许代选**）

1. **分组因子及其水平顺序**（`design$factors` / `design$ordered_levels`），以及采样设计
   （`independent` / `fixed_site` / `fixed_plot`）
2. **标准化策略**（rarefy / relative / CSS / none / subsample_avg）与**抽平深度**
3. **响应变量**是哪一列（`response$var`），或声明"本论文没有响应变量"
4. **是否使用系统发育树**（`phylo$enabled`），无树则哪些模块降级跳过
5. **出哪些模块**（模块勾选；不是全跑）
6. **零模型 runs 与判定阈值**（`assembly$runs` / `bnti_thr` / `rcbray_thr`）
7. **SEM 的 DAG 边表**（`sem$edges`）—— 结构必须由用户声明，**不许自动探索**

**规则**：每一项都要带推荐 + 理由；用户回「都行 / 你定 / 随便」**不算授权** —— 把推荐再问一次并说清后果。

---

## 标准化默认口径表

写进 `references/02-params.md`，每次开工复述并请用户确认。

| 分析 | 默认输入口径 | 理由 |
|---|---|---|
| 稀释/累积曲线、门/属水平组成 | **未抽平**（相对丰度） | 抽平会削平曲线、压缩优势类群 |
| α 多样性、β 距离、PERMANOVA、系统发育、组装、多样性分区、环境关联、建模 | **抽平后** | 与 原始研究口径一致，指数不受测序深度干扰 |
| 抽平深度 | 默认 = 最小正样本 reads（`depth = NULL`） | 保留最多样本；也可显式指定 |
| 抽平随机种子 | `CFG$opts$seed`（示范 123） | 全包唯一，可复现 |
| 抽平前先做的过滤 | ASV 读长 ≥ 10、出现率 ≥ 1（可调）、样本读长 ≥ 1000（可调） | 去测序噪声；**阈值不是普适的，须按数据类型调** |

> **重要例外**：本包示范数据为规范化后的示例数据，**极端单例化**（绝大多数 ASV 在多样本中只出现 1 次），
> 所以示范配置把出现率下限设成 1。换成常见的 16S/ITS 数据时，出现率下限通常应设 2 甚至更高。
> 详见 `references/02-params.md` 的「按数据类型的推荐」。

---

## 模块清单（5 个共享库 + 32 个分析模块）

### 共享库（`R/shared/`，无分析语义，只提供函数）

| 文件 | 内容 |
|---|---|
| `00_setup.R` | 配置解析、包管理（缺包降级）、路径、日志、缓存引擎、**seed 唯一来源**、`md_init()` |
| `10_io.R` | 五文件读写与对齐、核心对象访问器、跨模块产物访问器、结果表写出、格式化器、数据护栏 |
| `11_stats.R` | 组间比较与 CLD 字母、距离矩阵工具、Haversine、Mantel/MRM、bootstrap、VIF、抽平 |
| `12_plot_theme.R` | 主题、pt 字号双轨（`pt2gg`）、版式 mm 计算、配色、`md_save_figure`、`md_compose` |
| `13_plot_types.R` | 箱线、双因子面板、散点+拟合、条件效应曲线、热图、堆叠柱、森林图、直方图、四象限、NRI/NTI 面板 |

### 分析模块（`R/`，32 个）

| 组 | 模块 |
|---|---|
| 数据 | `01_data_health` · `05_assemble` |
| 组成与多样性 | `10_composition`（+ 可选指示种）· `11_alpha` · `12_niche_breadth` · `13_distances` · `14_ordination` · `15_beta_tests` · `16_beta_components` · `17_gamma_partition` |
| 系统发育 | `20_phylo_prep` · `21_phylo_diversity`（PD + MPD/MNTD + NRI/NTI）· `22_phylo_beta` |
| 群落组装 | `30_null_models`（βNTI + RCbray）· `31_assembly_process` · `33_distance_decay` · `34_variance_partition` |
| 环境关联与建模 | `40_env_summary`（+ 可选随机森林）· `43_mantel_mrm` · `44_constrained_ord`（db-RDA/RDA/CCA/capscale）· `45_response_curves`（GAM/GAMM）· `50_lmm` · `52_sem` · `53_sensitivity` |
| 空间 | `60_spatial_autocorr` · `61_map` |
| 响应变量关系 | `70_response_pattern`（可切 Scheirer-Ray-Hare）· `71_response_linkage`（可选 Procrustes / RF） |
| 网络与类群 | `80_network` · `81_taxon_group` · `82_taxon_share` · `83_taxon_diff` |
| 运行与导出 | `90_run_all`（runner）· `95_export_tables`（索引 + 产物校验） |

每个模块的**签名 / 输入 / 输出 / 消费的缓存 / 验收要点**见 `references/03-module-specs.md`。

> **一个模块可以出多张图**：`md_result(figures = c(...))` 接受字符向量，文件名用 `NNx_` 前缀（`md_module_of("80a_...")` 仍返回 `"80"`）。
> 目前 `80_network` 出满 6 张（`80a` 全局网+vs随机 / `80b` 分组网 / `80c` 周转+亚网络 / `80d` 节点属性 / `80e` 指标主图 / `80f` 指标附图）。
> 拆图原则：**版式差异大才拆**（网络快照 vs 箱线 vs 热图），同版式的多面板仍拼一张；画幅统一 183 mm（宽幅会让字看起来变小，见 pitfalls C34）。

**缺包时的行为**：该模块**不报错退出**，打印「缺少包 X，本模块跳过」+ 安装命令，继续后续模块。
安装命令与逐包用途 → `references/05-quickstart.md`。

---

## 硬约束（违者返工）

1. **参数一律从 `CFG` 读**，模块里不许写死分析参数；`CFG` 没有的键用 `%||%` 给默认值。
2. **seed 唯一**：所有随机过程前 `md_seed(CFG, offset)` 或传 `CFG$opts$seed`。模块里不许出现字面量种子。
3. **检查项不通过就停**，不"顺手"补数据。
4. 输出的每个统计量必须能在结果表中定位。
5. **距离矩阵单独存 rds**，绝不并入结果 data.frame。
6. 图中的 `geom_text/annotate` 的 `size` **必须** `pt2gg(pt)`；`element_text(size=)` 才用 pt。
7. **显著性标注**：分组数 > 4 一律用字母；≤4 也可用字母（本包统一字母）。
8. 拼图已包含的散图 = 冗余（`figure$save_singles=FALSE` 默认）；其他情况一律不删。
9. **`data/` 永不改写** —— `md_assert_not_data()` 拦写出口，`01` 记录 md5、`95` 交付前比对。
10. **产物必须过校验**：`95_export_tables` 检查 PDF > 0 且可渲染、xlsx 按 sheet 名可读且非空、**图中文字纯 ASCII**。
11. **图中文字必须纯 ASCII**（码点 ≤ 0x7E）。中文只允许出现在日志、结果表与 Note sheet 里。
    `95_export_tables` 会用 `pdftools::pdf_text()` 抽出每张图的文字层做审计，出现任何非 ASCII 可见字符即 **FAIL**。
    替换表见 `references/04-pitfalls.md` B12。
12. **视觉口径唯一真身在 `R/shared/12_plot_theme.R`**，且必须与 原实现 对齐
    （主题 / 字号 / 抖动点 / 回归线 / `cols_*` 色板）。**不要在任何模块里自造主题或配色**；
    要改就用 `CFG$figure$palette` 与 `CFG$plot$box` 覆盖。见 `references/06-figure-map.md` §4。
13. **分析必须可缓存**：分析段用 `md_cache(CFG, "NN_x", {...})` 或 `md_stage(CFG, "NN_analysis", {...})` 包起来，
    图与结果表放在缓存之外。这样 `--figures-only` 重出图时**不会重跑分析**（绘图样式不进缓存键）。
    **一个模块的图里混进分析计算，就等于每次改配色都要重算一遍。**
14. **箱线图默认照 nifH**：空心箱体 + `color_by` 因子给边框与抖动点上色 + 字母画在每组 `max(y)×letters_at` +
    右上角 KW 注释 + y 轴从 0 起 + 刻度朝内 + **不标样本量**。逐项可用 `CFG$plot$box` 改。
    画字母前必须保证字母的 `group` 与图中 x 水平**同名**（`inter_sep` 一致），否则字母会被丢弃。
15. **往已有模块加能力，默认必须关**。本轮给 `10`（指示种）、`40`（随机森林）、`70`（Scheirer-Ray-Hare）、
    `71`（Procrustes / RF）各加了一条可选路径，全部由 config 开关控制且**默认 `FALSE`** —— 只有这样，
    默认配置跑出来的图与表才与加这批能力之前**逐像素一致**。验收：`--force` 后与基线指纹比对，
    新增能力关闭时必须是 IDENTICAL（见 `references/04-pitfalls.md` C25）。
16. **每种分析的多方法一律做成 config 选项，并在流程第 3 步请用户一次勾完**。
    不许在模块里替用户选方法；`references/01-method-menu.md` 的每一行都要带**推荐 + 理由**。
17. **风格键统一收进各配置块的 `style` 子列表**。`md_cfg_for_hash()` 会把
    `network$style` / `taxon_group$style` / `taxon_diff$style` / `niche$style` /
    `constrained_ord$style` / `response_curves$style` / `composition$indicator$style` /
    `env$rf$style` 整体剔除 —— 于是改配色、字号、**网络布局**都只重出图，不触发重算。
    新加样式参数时**必须放进 `style`**，否则改一次配色就要重跑一次分析。

---

## 文件索引（按需读取，不必一次全读）

| 何时读 | 文件 |
|---|---|
| 要交数据 / 命名 / 五文件格式 | `references/00-data-contract.md` |
| 科学问题 → 选哪些模块/方法 | `references/01-method-menu.md` |
| **要协商参数（每项带推荐）** | `references/02-params.md` |
| 要查某模块的签名/输入输出/验收点 | `references/03-module-specs.md` |
| 报错 / 结果怪 / 口径存疑 | `references/04-pitfalls.md` |
| 填配置 | `config/config_template.R`、`config/config_demo.R` |
| 跑起来 | `R/90_run_all.R`、`run/manifest.R`、`run/_TEMPLATE_run_all.R` |
| **改这个包本身**（自检 / 回归基线 / 验收硬线 / 禁止项） | `references/08-maintenance.md` |
| 图件体检、结果表指纹（维护用） | `scripts/audit_figures.R`、`scripts/snapshot_results.R`、`assets/compare_figures.R` |
| 读函数实现 | `R/shared/*`、`R/NN_*.R` |
| 建新项目 / 要示范数据 | `assets/project_skeleton/`、`assets/make_demo_data.R`、`data/` |
| 给人看的介绍（三步上手 / 模块清单 / 装包 / 纪律 / 示范数据说明） | `references/05-quickstart.md` |
| 查「哪张图来自哪个模块、对应论文哪张图」 | `references/06-figure-map.md` |
| 查「默认关着的备用能力有哪些、怎么开」 | `references/07-spare.md` |

---

## 常见错误与处置（只列最高频的 6 条；完整踩坑库见 `references/04-pitfalls.md`）

| 现象 | 处置 |
|---|---|
| 一堆模块 SKIP，日志说"缺包 X" | `CFG$opts$auto_install = TRUE`（默认）会自动装；装不上（如 SpiecEasi 只在 GitHub）按日志里的 `install.packages(...)` 手工装，或用 `--only=` 只跑具备条件的模块 |
| 某模块 SKIP，说"缺文件 `data/05_taxon_group.csv`" | 做类群分组就按 `references/00-data-contract.md` §6 补上；不做就把 `CFG$taxon_group$file` 设为 `NULL` |
| `--figures-only` 还是跑了很久 | 该模块的分析没包进 `md_stage()` / `md_cache()` → pitfalls B14 |
| 改了配色/字号/网络布局，结果仍重算了一遍 | 样式参数没放进该配置块的 `style` 子列表 → 硬约束 17 |
| 图每次跑都不一样（节点位置 / 抖动点） | 布局要在**出图段**算并先 `md_seed()`；抖动点用 `position_jitter(seed = MD_JITTER_SEED)` → pitfalls C14 / C33 |
| 95 报 `图中文字纯 ASCII ... FAIL` | 图里混进了中文或非 ASCII 符号 → pitfalls B12 的替换表 |

其余（UpSet 挤成一条线、类群分组全是 `Unmapped`、图是 0 字节、图上出现「4」而不是字母等）
一律先查 `references/04-pitfalls.md`。
