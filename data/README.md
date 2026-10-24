# data/ —— 示范数据（模板与自检专用）

> ⚠ **这里的文件只能用于演示与自检，不得用于发表。**

## 这是什么

36 样本 × 327 ASV 的示范子集，由 `../assets/make_demo_data.R` 从一篇已定稿的
nifH 扩增子论文数据中**只读抽取**而来（固定 seed = 123）。

| 文件 | 内容 |
|---|---|
| `01_asv_table.csv` | 327 ASV × 36 样本丰度矩阵（每样本丰度前 10 个 ASV 的并集） |
| `02_taxonomy.csv` | `asv_id` + Kingdom..Species（已去 `k__/p__` 前缀，未分类写 `Unclassified`） |
| `03_sample_meta.csv` | `sample_id` + `Site` + `Zone`(3 水平) + `Vegetation`(3 水平) |
| `04_env_response.csv` | `sample_id` + `Latitude` `Longitude` + `MAT` `MAP` + 11 项土壤指标 + `EF`（响应变量） |
| `05_taxon_group.csv` | **可选第 6 文件**：`taxon,group` 两列，示范里把 45 个属分成 Symbiotic / Associative / Free_living / Other 四组（模块 81 用） |
| `99_tree.nwk` | Newick 树，327 tip，与 ASV 一一对应 |

## 采样设计

```
3 Zone（Coastal / Nearshore / Xisha，离大陆由近到远）
  × 3 Vegetation（Bare / Herb / Tree，离海岸线由近到远）
  × 4 重复（优先取自不同 Site）
  = 36 样本，跨 11 个站点
抽平深度 = 9968 reads
```

**为什么不含 Nansha、不含 Shrub**：Nansha 是人工岛（原论文为它专门做过剔除敏感性分析，
放进示范数据会混淆"示范结论"与"论文结论"）；Shrub 在 Coastal / Nearshore 只有 3 个样本，
凑不齐 4 重复。详见**包外**构建记录 `analysis-packages/_build-notes/micro-div-pipeline/DECISIONS.md` 的 D-03
（该目录不在本包内，随本包单独分发时不存在）。

## 这批数据的真实特征（不是程序错误）

- **极端单例化**：原数据中大量 ASV 在多样本中只出现 1 次。
  所以示范结果里会出现：**核心类群几乎为 0**（原论文实测 Core 0.06%、Transient_rare 98%）、
  **β 周转极高**、**漂变占绝对主导**（原论文实测 85.13%）、**环境解释力只有百分之几**。
- 这些正是那篇论文的核心发现，不是 bug。

## 怎么重新生成

```powershell
cd <micro-div-pipeline>
Rscript --vanilla assets/make_demo_data.R
```

默认从 `path/to/source/data` 读取源数据（**只读**）。
换源目录用环境变量：`$env:MD_SRC='<路径>'`。

## 数据护栏

- `data/` 目录**永不被本包改写**。`md_assert_not_data()` 会拦下所有写向此目录的操作。
- `01_data_health` 会记录这五个文件的 md5；`95_export_tables` 在交付前比对，
  不一致会在 `results/00_index.xlsx` 的 `Validation` sheet 里报 **FAIL**。
  （`05_taxon_group.csv` 是**可选**输入，不参与 md5 比对 —— 换分组表不需要重跑全流程。）
- 所以：**改了数据就要重跑**（用 `--force` 或先清空 `results/` 与 `figures/`）。

## 用你自己的数据

按 `../references/00-data-contract.md` 把数据整理成同样的五个文件（做类群分组分析时
再加可选的 `05_taxon_group.csv`），再照 `../assets/project_skeleton/` 建新项目、
改 `config/config.R`。
