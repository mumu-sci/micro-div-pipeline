# 00 · 数据契约（五文件 + 一个可选文件）

> 本包只吃 **5 个必需纯文本文件** + **1 个可选纯文本文件**。不读 phyloseq RDS、不读 QIIME2 `.qza`、不读 Excel。
> 理由：纯文本跨平台（QIIME2 / mothur / DADA2 / USEARCH 都能导出）、无编码坑、可 diff、可版本管理。

---

## 1. 目录与文件

```
<论文项目>/
├─ data/
│  ├─ 01_asv_table.csv      必需
│  ├─ 02_taxonomy.csv       必需
│  ├─ 03_sample_meta.csv    必需
│  ├─ 04_env_response.csv   可选（没有环境/响应变量时留空文件或不提供）
│  ├─ 05_taxon_group.csv    **可选**（做「类群 -> 分组」分析时提供；模块 81 用）
│  └─ 99_tree.nwk           可选（没有树就把 config 的 data$tree 设为 NULL）
```

文件名可以在 `config/config.R` 的 `data$*` 里改，但**建议保持默认**，方便复用。

---

## 2. `01_asv_table.csv` — ASV × 样本丰度矩阵

| 要求 | 说明 |
|---|---|
| 编码 | UTF-8（无 BOM） |
| 分隔符 | 逗号 |
| 首行 | 表头。**第一列必须是 ASV 标识列**（默认名 `asv_id`，可在 `data$id_col` 改） |
| 其余列 | 每列一个样本，列名 = 该样本的 `sample_id` |
| 取值 | **原始 reads（整数）**，不要预先抽平、不要预先转相对丰度、不要预先做任何标准化 |
| 缺失 | 写 `0`，不要留空 |

示例：
```csv
asv_id,QH1-L,QH1-H,QH1-M,XW1-L
ASV0001,1532,88,4,990
ASV0002,0,421,77,0
```

**注意**
- `asv_id` 必须唯一，且在 `02_taxonomy.csv` 与 `99_tree.nwk` 的 tip.label 里一致。
- 列名不要用 R 的特殊字符开头的形式（如 `X1`、`2A`），R 的 `read.csv(check.names=FALSE)` 会原样保留，但下游匹配容易出错。
- 不要包含 `taxonomy` 之类的注释列在这个文件里 —— 注释放 `02_taxonomy.csv`。

---

## 3. `02_taxonomy.csv` — 分类注释

| 要求 | 说明 |
|---|---|
| 首列 | 与丰度表同名的 ASV 标识列（`asv_id`） |
| 其余列 | 分类层级，**列名必须与 `data$tax_ranks` 一致**（默认 `Kingdom, Phylum, Class, Order, Family, Genus, Species`） |
| 取值 | 纯类群名，**不带前缀**（不要写 `p__Proteobacteria`，写 `Proteobacteria`） |
| 未分类 | 写 `Unclassified`（也可写空或 `NA`，包会归一） |
| 可少不可多 | 只有到属也可以；多出来的列会被忽略，但**顺序不能乱** |

示例：
```csv
asv_id,Kingdom,Phylum,Class,Order,Family,Genus,Species
ASV0001,Bacteria,Proteobacteria,Alphaproteobacteria,Rhizobiales,Bradyrhizobiaceae,Bradyrhizobium,Unclassified
```

**如果你的原始注释带 `k__/p__/c__` 前缀**（QIIME2 常见），先自己拆列：
```r
tax <- read.csv("raw_taxonomy.csv")
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
sp <- strsplit(tax$Taxon, "; ?")
tax[ranks] <- do.call(rbind, lapply(sp, function(x) {
  x <- sub("^[a-z]__", "", x); length(x) <- 7; x[is.na(x)] <- "Unclassified"; x
}))
```
（示范数据的生成脚本 `assets/make_demo_data.R` §5 里就是这么做的，可以直接照抄。）

---

## 4. `03_sample_meta.csv` — 样本元数据

| 要求 | 说明 |
|---|---|
| 首列 | `sample_id`，必须与丰度表列名**一一对应** |
| 其余列 | 任意个**分组因子**与协变量。因子列建议用字符串，水平顺序在 `config` 的 `design$ordered_levels` 里声明 |
| 必备 | 在 `design$factors` 里声明的列必须存在 |
| 推荐 | 有一列**站点标识**（`design$site_var`，如 `Site`），用于随机效应与站点级聚合；有站点跟踪设计时必需 |

示例：
```csv
sample_id,Site,Zone,Vegetation,Replicate
QH1-L,QH,Bare,Coastal,1
QH1-H,QH,Herb,Coastal,1
```

**注意**
- `sample_id` 里的连字符、点号都行，但**不要有前后空格**（包会 `trimws`，但列名对不上时会静默丢样本）。
- 不要在这里放环境数值 —— 环境与响应变量放 `04_env_response.csv`。
- **不要在元数据里放地理坐标和响应变量**（放 04），否则列名冲突时可能被覆盖。

---

## 5. `04_env_response.csv` — 环境变量与响应变量

| 要求 | 说明 |
|---|---|
| 首列 | `sample_id`，与丰度表列名对应 |
| 其余列 | 环境变量（连续量）+ 可选的地理坐标 + 可选的响应变量 |
| 坐标列名 | 若要用空间模块，需要 `Latitude` 与 `Longitude` 两列（可在 `config$env` / `md_geo_dist()` 的参数里改） |
| 响应变量 | 一列连续量，列名写进 `config` 的 `response$var` |
| 缺失 | 包不做插补，**缺失会被 `complete.cases` 静默剔除**，请自己先处理或如实保留并在表注说明 |

示例：
```csv
sample_id,Latitude,Longitude,MAT,MAP,pH,SWC,EF
QH1-L,20.25,109.93,24.1,1646,7.69,0.127,14.62
```

**响应变量是什么？** 任何"你想解释的生态系统/过程量"，例如：酶活、生产力、呼吸速率、
养分转化速率、病原抑制率、降解率。示范数据里叫 `EF`（ecosystem function）。
**没有响应变量就不提供这个文件**，并把 `config` 的 `response$enabled` 设为 `FALSE`（70/71 会自动跳过）。

---

## 6. `05_taxon_group.csv` — 类群 → 分组映射（**可选**，模块 81 用）

| 要求 | 说明 |
|---|---|
| 是否必需 | **否**。不提供（或 `config` 的 `data$taxon_group = NULL`、`taxon_group$file = NULL`）时模块 81 自动跳过，**其余模块完全不受影响** |
| 列 | 两列：`taxon` 与 `group`（列名可由 `config$taxon_group$taxon_col` / `$group_col` 改） |
| `taxon` | 与哪一个字段匹配由 `config$taxon_group$match_by` 声明：`ASV`（用 `asv_id` 直接匹配）或 `Phylum` / `Class` / `Order` / `Family` / `Genus` / `Species`（用 `02_taxonomy.csv` 的对应列匹配） |
| `group` | **分组标签，含义不限** —— 功能群 / 生态位 / R-K 对策 / 生活型 / 任何自定义分类 |
| 组数 | 任意（1 个组没有意义，建议 2–8 个；配色超出内置板时会用兜底定性色） |
| 未收录的类群 | `config$taxon_group$drop_unmapped = FALSE`（默认）时归入 `Unmapped`，结果表里单列一行；`TRUE` 时从分析中剔除 |

示例（示范数据里就是这个，按固氮生态类群分）：

```csv
taxon,group
Frankia,Symbiotic
Bradyrhizobium,Symbiotic
Azotobacter,Associative
Pseudomonas,Free_living
Proteobacteria,Other
```

**为什么单独一个文件**：分组方式因论文而异（16S 常按门/纲，ITS 常按功能群或营养型，
nifH 常按共生/联合/自生），而且同一篇论文也可能同时想看好几种分组 ——
所以映射表是**数据**不是**代码**，换论文只换这张表，模块代码一行不动。

---

## 7. `99_tree.nwk` — 系统发育树

| 要求 | 说明 |
|---|---|
| 格式 | Newick（`.nwk` / `.tre` / `.tree`） |
| 有根性 | 有根无根都行；无根时包按 `config$phylo$root_method` 自动生根（默认 midpoint） |
| tip.label | **必须与 `asv_id` 完全一致**（含大小写；引号会被自动剥离） |
| tip 数与 ASV 数的关系 | 不必一一对应。包会检查覆盖率：< 50% 会**阻断**并建议关闭系统发育模块；50–95% 会**警告**并以子集口径计算（须在表注披露） |

**上游怎么建树**：本包**不建树**，只接收现成 Newick。常见做法：
- QIIME2：`qiime phylogeny align-to-tree-mafft-fasttree`
- 命令行：`mafft rep_seqs.fasta > aln.fasta` → `iqtree -s aln.fasta -m MFP`
- 已有参考树：用 `ape::keep.tip(tree, asv_ids)` 剪到你的 ASV 集

---

## 8. 命名规范

| 对象 | 规范 | 示例 |
|---|---|---|
| ASV/OTU 标识 | 字母数字 + 下划线，全局唯一，**不要含逗号/引号/空格** | `ASV0001` |
| 样本标识 | 字母数字 + `-` / `_` / `.`，全局唯一，**不要含逗号/空格** | `QH1-L` |
| 因子水平 | 短字符串，**不要含 `_` 与 `|`**（`_` 会被用于交互标签拼接） | `Coastal` `Bare` |
| 环境变量 | 短字符串，与 `config$env$vars` 一致 | `MAT` `pH` `SWC` |
| 响应变量 | 一列，短字符串 | `EF` |

---

## 9. 数据护栏（硬约束）

- **`data/` 目录永不被改名**。`md_assert_not_data()` 拦下所有写向 `data/` 的操作。
- `01_data_health` 记录五文件的 md5 到 `results/01_data_health.xlsx`（Files sheet）与
  `cache/01_health.rds`；`95_export_tables` 交付前比对，不一致会在 Validation sheet 里报 **FAIL**。
- 所以：**跑完一轮之后如果改了数据，必须重跑**（缓存签名会自动失效，但旧结果还在 `results/` 里，
  请用 `--force` 重跑或先清空 `results/` `figures/`）。

---

## 10. 换数据类型时最容易踩的三个坑

1. **注释列的前缀与层级数**。QIIME2 的 `k__/p__/...`、SILVA 的 7 级、UNITE 的 ITS 注释（`k__Fungi;p__...`）
   都不一样。统一规则：拆成独立列、去前缀、未分类写 `Unclassified`。
2. **ITS 通常没有可靠的系统发育树**。把 `data$tree` 设为 `NULL`、`phylo$enabled = FALSE`；
   `13_distances` 会自动只算 Bray/Jaccard，`14_ordination` 会跳过 wUniFrac，`20–22` 全部跳过。
3. **出现率过滤阈值不能照搬**。示范数据（nifH）里 97% 的 ASV 是单例，
   所以示范配置用 `min_prev_samples = 1`。常见 16S/ITS 数据一般用 2–5，
   甚至可以先按"至少在 10% 样本中出现"过滤掉噪声。**看 `results/01_data_health.xlsx` 的
   `Taxonomy` 与 `Sample_reads` 再决定**，不要盲抄。
