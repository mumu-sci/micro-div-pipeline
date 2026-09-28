# 新论文项目骨架

把本目录复制成你的论文项目，然后按下面三步填。

```
<论文项目>/
├─ config/config.R           ← 从本骨架复制，或从包的 config/config_template.R 复制
├─ data/                     ← 放你的五个文件（见 references/00-data-contract.md）
│  ├─ 01_asv_table.csv
│  ├─ 02_taxonomy.csv
│  ├─ 03_sample_meta.csv
│  ├─ 04_env_response.csv
│  ├─ 05_taxon_group.csv     ← **可选第 6 文件**：类群 -> 分组映射（模块 81 用；不做就删掉）
│  └─ 99_tree.nwk
├─ run/run_all.R             ← 从包的 run/_TEMPLATE_run_all.R 复制
└─ cache/  results/  figures/  logs/   （会自动创建，不用手工建）
```

## 三步

1. **填数据**：按 `references/00-data-contract.md` 把数据整理成五个纯文本文件放进 `data/`。
   （`data/` 下的 `*.csv` 模板只给了表头，照格式填。）

2. **改配置**：`config/config.R` 里必须由你拍板的七项（`SKILL.md` 的「七道硬门」）：
   - `design$factors` 与 `design$ordered_levels`
   - `design$sampling` 与 `design$site_var`
   - `normalize$method` 与 `normalize$depth`
   - `response$var`（或 `response$enabled = FALSE`）
   - `phylo$enabled` 与 `data$tree`
   - `assembly$runs` 与阈值
   - `sem$edges`（要跑 SEM 就填边表，否则留空 → 52 会跳过）
   其余参数见 `references/02-params.md`。

3. **跑**：
   ```powershell
   cd <论文项目>
   D:/R-4.5.2/bin/Rscript.exe --vanilla D:/01-Academic/00-学术文件/analysis-packages/micro-div-pipeline/R/90_run_all.R --config=config/config.R --project=.
   ```

## 首次运行建议顺序

```powershell
# ① 体检（会拦截 BLOCK 级问题）
... --only=01
# ② 装配与标准化（产生 cache/05_core.rds，所有下游都依赖它）
... --only=01,05
# ③ 主力分析
... '--only=01,05,10,11,13,14,15,17,20,21,40,50,70'
# ④ 其余
... （不带 --only，全跑）
```
每一步都看 `logs/summary_*.csv` 与 `results/00_index.xlsx`。

## 目录说明

| 目录 | 内容 | 能否删 |
|---|---|---|
| `data/` | 源数据五文件 | **不能删，且永不改写** |
| `config/` | 配置 | 不能删 |
| `cache/` | 计算缓存（含 `05_core.rds`） | 可删，但删了要重跑耗时计算 |
| `results/` | 结果表 xlsx | 可删（重跑即重建） |
| `figures/` | 图件 pdf | 可删（重跑即重建） |
| `logs/` | 运行日志与汇总 | 建议留档 |
