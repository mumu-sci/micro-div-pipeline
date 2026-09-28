# 08 维护与回归 —— 「改这个包本身」的流程

> 这份文档面向**改包的人**，不是面向跑分析的论文作者。
> 论文作者只需要看 `references/05-quickstart.md`。

---

## 1. 唯一真身清单（改之前先确认改的是真身）

| 事项 | 唯一真身 | 说明 |
|---|---|---|
| 模块清单 / 执行顺序 / 依赖 / 是否重计算 | `run/manifest.R` | 加模块、改顺序、改 `needs` 只动这里 |
| 视觉口径（主题、字号、抖动点、回归线、色板） | `R/shared/12_plot_theme.R` | **任何模块都不许自造主题或配色**（硬约束 12） |
| 图件归属与画幅 | `references/06-figure-map.md` | 哪张图来自哪个模块、对应论文哪张图 |
| 参数默认值与推荐 | `config/config_template.R` + `references/02-params.md` | 两处必须同步 |
| 缓存键的构成 | `md_cfg_for_hash()`（`R/shared/00_setup.R`） | 新增样式参数必须放进对应配置块的 `style` 子列表（硬约束 17） |

---

## 2. 改包的标准循环（6 步）

```text
1. Rscript --vanilla R/90_run_all.R --check          # 只读自检：语法 / 依赖 / 列名 / 输入指纹
2. Rscript --vanilla scripts/snapshot_results.R --out=logs/fp_before.csv
   （图件基线：先 cp figures/*.pdf 到一个临时目录）
3. 改代码
4. Rscript --vanilla R/90_run_all.R --force          # 全量重算
5. 逐像素比对 + 图件体检：
     Rscript --vanilla assets/compare_figures.R <基线目录> figures <out.csv>
     Rscript --vanilla scripts/audit_figures.R
6. 数据表比对：
     Rscript --vanilla scripts/snapshot_results.R --out=logs/fp_after.csv
     然后按脚本头部那行 R 代码 diff 两份 csv
```

`--check` 是 `SKILL.md` 第 6 步「自检（不执行）」的可执行版本，四项：
语法 `parse()`、依赖是否齐（**只报告不安装**）、配置声明的列名是否都存在、
五文件输入指纹（路径 + 行数 + ASV 数 + md5）。有 BLOCK 才非零退出；缺包只算 WARN。

---

## 3. 验收硬线（**这一条最容易漏**）

| 改了什么 | 验收线 |
|---|---|
| 只加**默认关**的可选能力 | 关闭时，默认配置跑出来的图与表必须与改动前 **逐像素 / 逐 sheet 一致**（IDENTICAL） |
| 只改某个模块 | **未触及的模块必须 IDENTICAL**；触及的模块按预期变化，并说明变了什么 |
| 改共享库（`R/shared/`） | 影响面 = 全部模块，必须全量比对 |

**允许的例外**：`max_abs_diff ≤ 1/255`（约 0.004）且 `frac_gt5pct = 0` 的**亚像素差**。
已知成因：抖动点定种子（`MD_JITTER_SEED`）会让箱线图上的点位置有极微小位移，
渲染后灰度差 < 1 级，肉眼与打印都不可辨。超过这个量级就不算"没变"。

**PDF 的 MD5 不能用作基线**：`cairo_pdf` 会在文件里写入时间戳，同一张图每次导出的字节都不同。
图件必须**逐像素**比；结果表必须**逐 sheet**比（xlsx 的 zip 也带时间戳）。

---

## 4. 基线放哪

| 资产 | 位置 | 说明 |
|---|---|---|
| 图件指纹基线 | 包外 `analysis-packages/_build-notes/micro-div-pipeline/round5-baseline/fingerprint.csv` | **不在本包内**，随本包单独分发时不存在 |
| 构建记录 | 包外 `_build-notes/micro-div-pipeline/` 的 `VALIDATION.md` / `DECISIONS.md` / `BUILD_LOG.md` | 每轮的验证记录、决策与沿革 |
| 已归档代码 | 包外 `_build-notes/micro-div-pipeline/removed-code/` | 例如被 `44` 取代的 `42_constrained_ord.R` |
| 论文目录只读快照 | 包外 `_build-notes/micro-div-pipeline/thesis-readonly-baseline.csv` | 证明改包过程**没有碰过论文原始文件** |

> 也就是说：**本包不自足于维护场景**。要完整跑一遍回归，需要 `analysis-packages/` 这一层，
> 而不是只有 `micro-div-pipeline/`。分发单包时请一并说明这一点。

---

## 5. 三条固定命令

```powershell
# 只读自检（语法 / 依赖 / 列名 / 输入指纹）
Rscript --vanilla R/90_run_all.R --check

# 图件体检：非 ASCII 审计 + 位图 XObject 计数 + 页尺寸 + 文字层是否为空
Rscript --vanilla scripts/audit_figures.R

# 结果表逐 sheet 元素指纹
Rscript --vanilla scripts/snapshot_results.R --out=logs/fp_after.csv

# 图件逐像素比对（改动前后各渲染一份 PNG 目录，或用 PDF）
Rscript --vanilla assets/compare_figures.R <基线目录> figures <out.csv>
```

`audit_figures.R` 的退出码：出现非 ASCII 字符或某张图文字层为空 → `1`。
"含位图 XObject"只提示不判失败（既有的 `40_env_summary` 与 `43_mantel_mrm` 本来就含位图）。

---

## 6. 禁止项

1. **不许改 `data/`** —— 源数据只读。`md_assert_not_data()` 拦写出口，`01` 记 md5、`95` 交付前比对。
2. **不许在模块里写死分析参数或字面量种子** —— 一律从 `CFG` 读、一律走 `md_seed()`（硬约束 1、2）。
3. **不许把历史决策只写进对话** —— 每轮的决定要回写包外的 `DECISIONS.md`，否则下一个人会用错口径。
4. **不许把样式参数漏在 `style` 子列表之外** —— 否则改一次配色就要重跑一遍分析（硬约束 17）。
5. **不许为了"让图好看"改统计口径** —— 视觉与数值必须分开改，改完按第 3 节的验收线比。

---

## 7. 版本控制现状

**本包目前未纳入版本控制**（目录下没有 `.git`）。`.gitignore` 是**为将来预留**的，
它把 `cache/ results/ figures/ logs/` 视为可重建产物 —— 这几类确实都能由
`Rscript --vanilla R/90_run_all.R --force` 重建。

在补上版本控制之前，唯一的变更历史就是包外 `_build-notes/micro-div-pipeline/` 里的
`BUILD_LOG.md` / `DECISIONS.md` / `VALIDATION.md`。**改动前请先在这三份里记一笔。**
