# 04 · 踩坑库

> 分两部分：**A 段**是从 原始研究原项目移植过来的工程坑（那些坑在那个项目里都真实发生过）；
> **B 段**是本包建设与验证过程中新踩到的。改代码前先扫一遍。

---

## A. 从 原始研究移植的工程坑

### A1 · Rscript `--vanilla` 下中文路径读不到
`Rscript --vanilla` 用 C locale，`read_excel("data/环境数据.xlsx")` 会报 "path does not exist"。
**处置**：每个入口脚本开头 `Sys.setlocale("LC_CTYPE", "")`（本包 `md_init()` 已内置）。
更稳的做法是**输入文件一律用纯 ASCII 文件名**（本包五文件契约就是这么定的）。

### A2 · `library()` 必须在 `unit()` / `theme_set()` 之前
原脚本（nifH）的 §0 里 `grid::unit` 在 source 时就求值，而 ggplot2 对 grid 只是 Imports，
不自动挂搜索路径 → source 直接失败。
**处置**：本包 `md_init()` 里先 `library(ggplot2); library(patchwork)`，其余按需 `md_load()`。

### A3 · `tryCatch(warning=)` 会截断 source
用 `tryCatch(..., warning = function(w) ...)` 捕获 warning 会 unwind，
`source()` 只跑到第一条 warning 就结束，后面的代码根本不执行（而且不报错）。
**处置**：本包 runner 用 `withCallingHandlers(..., warning = function(w) { log; invokeRestart("muffleWarning") })`。

### A4 · `vegan::adonis2(by = "margin")` 在本机只返回末项
实测：`by="margin"` 只返回最后一项 + Residual + Total，导致 PERMANOVA 表里只有交互项一行。
**处置**：**一律用 `by = "terms"`**（`15_beta_tests` 已内置）。

### A5 · 直读 xlsx 不能按 `sheetN.xml` 序号对表名
`xl/worksheets/sheet1.xml` 的编号**不等于** `workbook.xml` 里的 sheet 顺序。
按序号对齐会**静默读到别的工作表**且不报错（nifH 项目据此误判过一张表"名实错置"）。
**处置**：一律 `openxlsx::read.xlsx(path, sheet = "表名")`（本包 `md_read_xlsx()`）。

### A6 · openpyxl 读不了 openxlsx 生成的 xlsx
openxlsx 会写出**悬空 drawing 引用**（`xl/drawings/drawing1.xml` 不存在），
`openpyxl.load_workbook()` 直接 `KeyError`。
**处置**：下游要用 Python 处理时，把 `CFG$table$write_csv` 设为 `TRUE` 同步落一份 csv；
或改用 zipfile + `workbook.xml` / `rels` 自写最小读取器。

### A7 · 改图后的验收不能用 MD5
PDF 内嵌 `CreationDate`，重跑后**未改动的图 MD5 也会变**，无法判断有没有误伤。
**处置**：渲染成位图逐像素比对（`pdftools::pdf_render_page` + `magick`）。实测未改动图 `meanAbsDiff = 0.000000`。

### A8 · 拼图版式的三个具体坑
- **图例溢出被裁**：图例项多（>10）时在 183mm 下并排会溢出，最右一列被裁成半截词，
  PDF 不报错只少半个词。处置：`legend.box = "vertical"` + 图例字号压到 2/3 + `md_guide_flex()`。
- **直方图两侧标注被裁**：标注锚在 `x = ∓2.5`（数据范围之外）时会被裁掉。
  处置：锚在面板内部 `x = ∓Inf` + `hjust` 朝内（`md_hist()` 已内置）。
- **45° 双行标签只在整页宽面板安全**：约 79mm 装 16 组必然叠字。
  处置：窄面板改单行（`x_labels_one_line = TRUE`）。

### A9 · 先抽平后删类群 → 文库量不再相等
nifH 项目里 `noC4` 对象是先抽平再删类群，导致文库量变成 27,495–30,553（不等），
α 多样性口径因此不干净。
**处置**：本包 `05_assemble` 的顺序固定为 **先过滤（含类群子集）→ 再抽平**，并在 `Filter_log` sheet 里
逐步记录每一步的 ASV / 样本数。

### A10 · 清单在两个 runner 之间漂移
nifH 项目 `_run_all.R` 与 `_regenerate_figures.R` 各复制了一份脚本清单，后来漂移了
（曾漏 12/14/15/16，导致"重出全部图表"其实没刷新几张表）。
**处置**：本包只有一份 `run/manifest.R`，两个模式（全跑 / `--figures-only`）都由它生成。

### A11 · 全局环境互相污染
原实现 runner 在同一 R 会话里顺序 `source`，顶层变量互相覆盖（02/05/08 都定义过 `results`）。
**处置**：本包 runner 自己的变量一律 `.run_` 前缀；模块只导出函数，数据只走 `cache/`。
**并且**：本包不再用"load(RData) 往全局注入几十个对象"的模式 —— 全部走 `md_core_*()` 访问器。

---

## B. 本包建设与验证中踩到的坑

### B1 · 示范数据按"平均相对丰度取前 300 ASV"会做不出可分析的群落
实测：nifH 全数据里 **大量 ASV 在多样本中**只出现 1 次**（论文的极端周转）。
按平均丰度取前 300，抽到的几乎全是"单样本爆发型"ASV —— 300 个里 240 个只出现在 1 个样本中；
按"36 样本内出现率 ≥ 2"筛只剩 104 个，按 ≥3 筛只剩 33 个；还会出现 0 reads 样本、抽平深度只剩个位数。
**处置**：改成 **「每个样本丰度前 k 个 ASV 的并集」**（k 取使并集 ≥ 目标数的最小值）。
示范数据最终 k = 10 → 327 个 ASV，覆盖 84% reads，每样本最低仍有 9968 reads。
**教训**：极端高周转的数据集上，"全局丰度排序"和"出现率过滤"都不是好的子集选择策略。

### B2 · 出现率过滤阈值不能跨数据集照搬
承 B1：示范配置里 `min_prev_samples = 1`，而常见 16S/ITS 数据一般用 2–5。
这个差异**必须写在配置注释与文档里**，否则用户换数据后会莫名其妙丢掉大半 ASV。

### B3 · PowerShell 会把 `--only=10,11` 的逗号拆成多个参数
结果 R 只收到 `--only=10`，其余模块静默不跑。
**处置**：命令行参数加引号（`'--only=10,11'`）；runner 里也做了兜底：
分隔符接受 `,` `;` `+` 空格，并支持 `--only 10 11` 的形式。

### B4 · PowerShell `Set-Content` 回写会破坏 UTF-8 中文源文件
用 `(Get-Content x) -replace ... | Set-Content x -Encoding UTF8` 处理含中文的 R 文件后，
文件里的中文注释全部变成乱码（`缺少` → `缂哄皯`），且不报错。
**处置**：**不要用 PowerShell 回写 UTF-8 源文件**；用编辑器/专用工具改。
（本包文件全部为 UTF-8 无 BOM，不要在 Windows 上用默认编码工具批量改。）

### B5 · `phyloseq::sample_data()` 需要非零列
`sample_data(data.frame(row.names = samples))` 会报
"Sample Data must have non-zero dimensions"。
**处置**：本包 `05_assemble` 干脆**不依赖 phyloseq**，抽平用自实现的
`md_rarefy()`（无放回超几何抽样，与 `rarefy_even_depth` 结果口径等价，更快、依赖更少、种子可控）。

### B6 · `md_need(c("phyloseq","ecology"))` 传的是包组名不是包名
早期实现的 `md_need()` 只处理 `length(what)==1` 的组名，传向量时会去找一个叫 `"ecology"` 的包，
导致模块被无理由 SKIP。
**处置**：`md_expand_pkgs()` 先把每个元素展开成真实包名再检查。

### B7 · 站点级坐标不能做样本级空间自相关
nifH 数据 147 行只有 14 个唯一 `(Lat, Lon)`，站内距离恒为 0。
**处置**：`60_spatial_autocorr` 在 `spatial$level = "auto"` 时检测唯一坐标数，
有重复就强制站点级聚合并如实登记"样本级不可检验"。示范数据 36 样本 / 11 唯一坐标，同样触发。

### B8 · `by = "terms"` 与 `md_vif()` 的共线交互
环境变量高度共线时（如 pH 与 Ca、MAT 与站点），VPA 的校正 R² 会变成负值、SEM 会不可识别。
**处置**：不是 bug，**不要截断为 0、不要删变量了事**；照实输出 + 在表注写明"校正 R² 为负表示解释力
不优于随机"，`50_lmm` 的 `md_within_unique()` 守卫会自动降协变量并记录。

### B9 · `geom_text(size=)` 的单位是 mm 不是 pt
忘了套 `pt2gg()` 时字会大 2.8 倍；而且 `geom_text` 的 size 在**构造时固化**，
事后 `+ theme()` 改不动 —— 窄面板必须在构造时传入小字号参数。
**处置**：`13_plot_types.R` 里所有图形基元都把字号作为参数暴露（`pt =`），并统一走 `pt2gg()`。

### B10 · 抽平深度自动取值会随过滤结果变化
`depth = NULL` 时会取"过滤之后"的最小正样本 reads。若中途改了过滤阈值，
深度会静默变化，历史结果的 α 多样性不可比。
**处置**：实际使用的深度写进 `results/05_assemble.xlsx` 的 `Normalization` sheet 与每个结果表的 `_meta` sheet。
**换过滤阈值后请用 `--force` 重跑全部下游。**

---

## C. 共享库的已知陷阱（子代理实测报告 + 已修）

### C1 · `theme_md(CFG)` 位置传参会报 `non-numeric argument to binary operator`
`theme_md()` 的签名是 `(base_size, base_family, comp, CFG)`，CFG 在第 4 位。
写 `theme_md(CFG)` 会把列表绑到 `base_size` 上。
**已修**：函数内加了容错（`is.list(base_size)` 时自动把 CFG 挪正）+ 类型守卫给出可读报错。
新代码仍然建议写 `theme_md(CFG = CFG)`。

### C2 · `md_site_aggregate()` 的 `vars` 语义是"白名单"
显式传 `vars` 时**只保留这些列**，其余列（包括残差协变量）一律丢弃；不传才是"全部数值列"。
做站点级聚合时如果要用残差协变量，**必须把它们并进 `vars`**。

### C3 · `md_box_two_factor()` 的 `letters` 表必须含 `group` / `letter` 两列
列名不对时报的是 `arguments imply differing number of rows: 0, 9, 1`，极难定位。
**已修**：加了列名校验与可读报错。

### C4 · `md_add_p_annot(g, labels, ...)` 的 `labels` 必须具名传
写成 `g + md_add_p_annot("...")` 会把字符串当 `g` 传进去，报 `argument "labels" is missing`。

### C5 · `md_compose(legend = "bottom")` 会覆盖子面板自己的 `legend.position`
默认值就是 `"bottom"`，所以会把你已经设好的图例位置顶掉。
要保留子面板设置请显式传 `legend = NULL`。

### C6 · 本机 OpenBLAS 与 sf / terra 冲突会导致进程 abort（退出码 134）
现象：runner 里加载 `sf` / `terra`（即 `61_map`）时
`OpenBLAS error: Memory allocation still failed after 10 retries`，R 直接 abort，
看起来像"61 崩了"。
**已修**：`md_init()` 里默认设 `OMP_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1` / `MKL_NUM_THREADS=1`
（本包瓶颈在置换与 R 循环，不在 BLAS，几乎无性能损失）。
需要多线程 BLAS 的用户可在调用 `md_init()` 前自行设置这三个环境变量。

### C7 · 箱线图字母"消失"，箱子上只剩样本量数字（2026-09-26 实测）
现象：`11_alpha.pdf` 每个箱子上方只有一个"4"（每组样本量），一个显著性字母都没有，且**不报错**。
根因：`11_alpha` 用 `local_group_var(sep = "_")` 生成 `Coastal_Bare` 这类分组名去算 CLD 字母，
而 `md_box_two_factor()` 内部用 `sep = " | "` 生成 `Coastal | Bare` 作为 x 轴水平；
`geom_text()` 收到的是**图上级别里不存在的 x 值**，ggplot 静默丢弃这些点。
处置（两处都要）：
1. 画图与生成字母必须用**同一个分隔符**（本包统一 `inter_sep = "_"`，在 `CFG$plot$box` 里配）；
2. `local_letters_df()` 在画字母前**校验** `letters$group` 是否都在 x 水平里，不匹配就 `md_warn` 报出来（不再静默）。
另外 `CFG$plot$box$show_n` 默认 `FALSE` —— 原实现箱线图本来就不标样本量。

### C8 · 图上出现中文（2026-09-26 实测）
现象：`05_assemble.pdf` 的 x 轴是"ASV 总读长 >= 10 / 出现样本数 >= 1 / 样本读长 >= 1000"。
根因：该模块把**中文的过滤日志 `step` 字段**直接当成了柱状图的 x 轴标签。
处置：过滤日志同时保留 `step`（中文，给人看表）与 `step_label`（**纯 ASCII**），**图只用 `step_label`**。
护栏：`95_export_tables` 用 `pdftools::pdf_text()` 抽出每张图的文字层，出现任何非 ASCII 可见字符即 FAIL。
**替换表**（图上文字一律用右列）：

| 原 | 改为 | 原 | 改为 |
|---|---|---|---|
| `²` | `R2` | `≥ ≤` | `>= <=` |
| `Δ` | `delta` | `×` | `x` |
| `ρ` | `rho` | `°` | `deg` |
| `β` | `beta` | `ø`（Sørensen） | `o`（Sorensen） |
| `→` | `->` | `– —` | `-` |
| `∩` | `x` | `±` | `+/-` |
| `≈` | `~` | 弯引号 `“ ” ’` | 直引号 `" '` |

注意：**plotmath 的 `R^2` 也会在 PDF 文字层留下 `²`**（`expression(R^2)` 实测会），
所以严格 ASCII 要求下请直接用 `R2`。
另：批量替换弯引号时**小心把 R 字符串截断**（`"…"中文"…"` 会报 `unexpected symbol`）——
中文串里的引号请用单引号 `'`。

### C9 · 不要在模块里自造主题或配色（2026-09-26 实测）
现象：早期版本在 `config_demo.R` 里自造了 `Zone = #2B7A78/#3AAFA9/#8AC6D1`、
`Vegetation = #D9A05B/#8FAE5B/#3E7C4A`，与 原实现 `cols_zone` / `cols_veg` 完全不同；
主题的 border 线宽、axis.line、图例标题字重、strip 背景也与 `theme_isme()` 不一致。
处置：`R/shared/12_plot_theme.R` 是**视觉唯一真身**，色板逐字照搬 nifH；config 只做**覆盖**不做**另立一套**。

### C10 · `--figures-only` 曾经会跳过重计算模块，导致图表不刷新（2026-09-26 修）
早期实现里 figures_only 直接跳过 `heavy = TRUE` 的模块 —— 那些模块的图**根本没被重画**。
现在改为：figures_only **照跑所有模块**，重计算部分由 `md_cache()` / `md_stage()` 命中缓存而跳过。
因此**每个模块的分析都必须包进缓存**，否则 figures_only 仍会重算它。
用法：把分析段包起来，图与表留在外面。
```r
md_stage(CFG, "11_analysis", {
  idx_m <- ...      # 这一段内**新建**的变量会被缓存
  tab_test <- ...
})
# 下面照常画图/写表，用的还是 idx_m / tab_test，每次都用新样式重画
```
缓存键 = 数据指纹 + 分析相关配置（**绘图样式不进键**）→ 只改配色/字号时必然命中。
已验证：连续两次 `--figures-only`，135 个数据 sheet 的指纹**完全一致**（数字没被缓存改动）。

### C11 · 把变量定义挪进缓存块后，块外仍在引用它
把 `n_used <- md_parallel(...)` 挪进 `md_cache(...)` 之后，缓存命中时块不执行，
块外引用 `n_used` 就报 `object 'n_used' not found`。
**处置**：凡是块外还要用的变量，**定义在缓存块之外**，块内只做昂贵计算。

### C12 · `ggtreeExtra::geom_fruit(axis.params = <变量>)` 会静默丢掉整条环轴（2026-09-26 实测）
现象：树的多环图里，某一环的**轴文字全部消失**（环名、丰度刻度一个都不出），但**不报任何错**。
根因：`ggtreeExtra` 1.20.1 的 `axis.params` 经过 `enquo`/`reset_params` 处理后，
传**已经存进变量的 list** 会失效；传**就地写的字面量 list** 才生效。
**处置**：照 nifH 原脚本的写法，`axis.params = list(axis = "y", text.size = ..., nbreak = 3)` **就地写**，
不要先赋值给一个变量再传进来。已在 `R/20_phylo_prep.R` 内注释标注规避。

### C13 · `--figures-only` 不会让图自动更新（已修）
早期实现里 `figures_only` 会**直接跳过** `heavy = TRUE` 的模块 —— 那些模块的图**根本没被重画**，
用户以为"改了配色重出图"生效了，其实拿到的还是旧图。
**已修**：`--figures-only` 现在照跑所有模块，重计算部分由 `md_cache()` / `md_stage()` 命中缓存而跳过。
**代价**：每个模块的分析都必须包进缓存，否则 figures_only 仍会重算它（见 C10）。

---

## C-2. 2026 补充轮（模块 12 / 44 / 45 / 80 / 81 / 82 / 83）新踩到的

### C14 · igraph 的布局是随机的 —— 不设种子，图每次都不一样
`layout_with_graphopt()` / `layout_with_fr()` / `layout_nicely()` 都从随机起点出发。
论文原脚本里 `<local>` 块先 `set.seed(...)` 再算坐标，而且**布局是在出图段算的**。
**处置**：布局放在**出图段**（`md_stage()` 之外），并在这之前
`md_seed(CFG, md_offset(CFG, "network_layout", 8230L))`。
放进分析段会导致"改布局要重算整个网络"。

### C15 · base R 的 `plot.igraph()` 不能直接进 patchwork
`plot.igraph()` 画在图形设备上，返回的 `NULL` 不能参与 `+` 拼接。
**处置**：`ggplotify::as.ggplot(function() plot(g, ...))` 包成 ggplot 对象再拼。
注意 `as.ggplot()` 里访问的变量必须在**当前求值环境**可见（原脚本用"提升到全局环境"
绕过了这件事；本包统一在同一个模块函数体内构造，不需要那一步）。

### C16 · ComplexHeatmap / pheatmap / UpSet 不是 ggplot —— 尺寸必须显式给
它们返回 gtable / 网格对象。`ggplotify::as.ggplot()` 能包装，但**不会**自动继承
`md_save_figure()` 的 mm 尺寸；行标签、列标签的字号也不会跟随 `theme_md()`。
**处置**：字号走各自的参数（`row_names_gp` / `fontsize_row`），并在包装后
用 `md_save_figure(..., width, height)` 明确给尺寸。

### C17 · UpSet 的组合数会指数爆炸
`ComplexHeatmap::make_comb_mat()` 在 8 个以上集合时会产出几百个组合，图被压成一条线。
**处置**：`CFG$taxon_share$max_comb`（默认 20）限制进入 UpSet 的组合数，
其余组合仍完整写进 `Intersection` sheet。集合数 >= 4 时不再画 Venn（只画 UpSet）。

### C18 · CLR 的伪计数会改变所有下游 p 值
`s4_supplementary_taxa.R` 的口径是"**最小正值 / 2**"，不是常见的 +1 或 +0.5。
**处置**：`md_clr(mat, pseudocount = NULL)` 就是原脚本口径；
要换口径必须 `CFG$taxon_diff$pseudocount` 显式给，并写进表注（口径变化要披露）。

### C19 · 类群映射表匹配不到时要"报告"而不是"静默丢"
`05_taxon_group.csv` 与 `taxonomy` 的某一层级对不上（拼写、大小写、前缀）时，
匹配失败的 ASV 会全部落空。
**处置**：模块 81 在 `drop_unmapped = FALSE`（默认）下把未匹配的归入 `Unmapped`
并在图上与表里都保留这一组；同时在 notes 里报出**未匹配的 ASV 数与占比**。
映射表覆盖不到 50% 的 ASV 时额外 `md_warn`。

### C20 · 网络在"样本少 + 出现率阈值高"时会被全部过滤掉
`min_prev = 0.02` 在 36 个样本上 = 至少出现在 1 个样本里才保留；
但在 6 个样本上同样的 0.02 也是 1，而 `min_prev = 0.2` 就要求 2 个样本，
样本量再小会把节点全砍光，`cor()` 直接报错。
**处置**：模块 80 在 `nrow(mat) < 3` 或 `ncol(mat) < 4` 时
**明确 SKIP 并写出原因**（当前样本数、过滤后 ASV 数、所用阈值），不要抛 R 的原始错误。

### C21 · SpiecEasi 不在 CRAN
`CFG$network$method = "sparcc"` 需要 `SpiecEasi`，它只在 GitHub（还依赖 JuliaCall）。
**处置**：自动装包会失败并留日志；模块 80 必须**回退到 `method = "cor"`** 并在 notes 说明，
不能让整条流水线挂掉。手工安装命令：
`remotes::install_github("zdk123/SpiecEasi")`。

### C22 · 自动装包的两条纪律
`md_need()` 现在默认会自动 `install.packages()`（`CFG$opts$auto_install`）。
两条纪律：
1. **同一个包一次运行只尝试安装一次**（否则每个模块各试一次会反复卡几分钟）；
2. 装包失败**不算 ERROR**，只让该模块 SKIP。完全离线请设 `CFG$opts$auto_install = FALSE`。

### C23 · patchwork 的 design 里"用不到的字母"= 空白区域，也可能是**预留位**
`40_env_summary` 的 design 是 `"ADD\nBCC"`，而默认只放 A/B/C 三个面板 ——
实测 **D 区是一整块空白**（渲染后该象限非白像素占比 0.0000）。
这不是渲染 bug：那个 D 就是给"可选面板"预留的位置。
往已有图里加面板时，**优先去填已有的预留区**（`all_panels$D <- p_rf`），
而不是改 design 字符串 —— 改字符串会挪动已有面板的位置，违反"不动已有图"的约束。

### C24 · `md_stage()` 只记住"段内新建的变量"
`md_stage(CFG, "NN_analysis", { ... })` 命中缓存时，把段内**新建**的变量灌回环境；
修改已有变量不会被记住。
**处置**：分析段里一律用新名字；模块外还要用的开关/参数在段外先定义好（同 C11）。

### C25 · 往已有模块加可选能力时，默认必须关
本轮往 `10` / `40` / `70` / `71` 各加了一条可选路径（指示种 / 随机森林 /
Scheirer-Ray-Hare / Procrustes+RF）。全部由 config 开关控制且**默认 FALSE** ——
只有这样，默认配置跑出来的图与表才与加这批能力之前**逐字节一致**。
验收方式：跑 `--force` 后与基线 MD5/指纹比对（新增能力关 = 必须 IDENTICAL）。

### C26 · `40` 的随机森林因变量必须在 `md_core_analysis()` 里，或是一个 α 指数
`md_core_analysis()` = meta + env + 派生变量（`soil`/`env` 的 PC 列），**不含** α 指数 ——
所以 `env$rf$response = "Shannon"` 会报"不在分析总表里"。
**处置**：模块 40 现在会**现算** `Observed / Shannon / Simpson / Pielou` 补齐
（`local_rf_alpha()`，口径同模块 11）；`PD`、排序轴等不在其中，要用就写总表里真有的列
（如示范数据里的 `EF`）。取不到时该面板跳过并记 notes，不报错。

### C27 · `pdp::partial()` 的输出列名随预测变量而变，不能直接 `rbind`
每个变量返回的列是"该预测变量名 + `yhat`"，列名各不相同；
`do.call(rbind, list_of_partials)` 会报 **`names do not match previous names`**。
**处置**：先统一成 `(Variable, x, yhat)` 三列再合并（`R/shared/11_stats.R::md_rf_fit_impl`）。
另外 `md_rf_fit()` 外面套了一层 `tryCatch` —— 它**永不抛错**，任何异常只降级为
`md_warn()` + 返回 NULL，避免一个可选面板把整个模块打成 ERROR。

### C28 · `any()` 遇到 NA 会返回 NA，`if (!any(...))` 直接报错
`indicator$fdr < thr` 在 `fdr` 为 NA 时给出 NA；`any(NA)` = NA，`!NA` = NA，
`if (NA)` 报 **`missing value where TRUE/FALSE needed`**。
**处置**：先显式消 NA：`sig <- !is.na(x) & x < thr`，再 `if (!any(sig))`。
本包所有"阈值判定"都要按这个写法。

### C29 · 一个模块出多张图时，`figures` 与 `figures_stale` 都要按向量处理
`md_result(figures = c("a.pdf","b.pdf"))` 是允许的（模块 80 出 6 张）。
`95_export_tables` 的 `Figures_stale` 检查必须 `rep(..., length(stale))`，
否则会写出长度不等的 data.frame 而报错。

### C30 · `facet_wrap` 的 facet 变量在「数据层」与「文本层」水平不一致 → 幽灵面板
现象：`80_network_metrics` 本该 6 格，实际出来 **7 格** —— 第 7 格没有箱体、没有 strip 文字，
只有 CLD 字母与 `R2 = ...` 注释，y 轴范围也和数据层对不上。
根因：`ggplot()` 的数据层 `MetricLabel` 是 6 个水平的因子，而传给 `geom_text(data = stats)`
与 `geom_text(data = lab)` 的**文本层**，其 `MetricLabel` 被 `factor(..., levels = <查表得到的 NA 向量>)`
压成了**别的水平集合**；`facet_wrap` 取所有层的并集，于是多出一格。
本例的具体写法错误是"把 display 名又当键查了一次标签表"：
`unname(labels[metric_levels])`，而 `metric_levels` 传进来时已经是 display 名。
**处置**：① facet 变量在**所有层**用同一份 `levels`；② 传参前想清楚是"键"还是"显示名"，
不要二次查表；③ 图出来后用 `pdftools::pdf_data()` 数一遍 strip 标题的出现次数
（每个度量应对应"因子个数"次），这是最省事的幽灵面板检测。

### C31 · `ggplotify::as.ggplot(<ComplexHeatmap::Heatmap>)` 会返回**空图**
现象：`83_taxon_diff` 的热图在 PDF 里是一片空白，还带着 `x` / `y` 数值轴与 0.00–1.00 刻度，
但行名、列名、`Z-score` 图例都在。
根因：`as.ggplot()` 对未 `draw()` 的 `Heatmap` S4 对象无法取到绘制内容，返回一个"只有空面板"
的 ggplot；后面再被主题化就长出默认数值轴。
**处置**：一律
```r
gr <- md_grid_grab(ComplexHeatmap::draw(ht, newpage = TRUE))   # 见 C39
g  <- ggplotify::as.ggplot(gr)
```
并加一道**有效性护栏**：`length(ggplot2::ggplot_build(g)$data) > 0`，不满足就换下一档渲染
（`pheatmap` → `geom_tile`）。同一条路也在 `82_taxon_share` 的 UpSet 上用。

### C32 · `md_compose()` 会把任何面板重新主题化 —— 传给它的图**必须自带数据**
`md_compose()` 结尾是 `pw & theme_md_comp(CFG)`，它会把 `axis.text` / `axis.title` 重新打开。
对普通 ggplot 无害（轴本来就在），但对两类图是灾难：
- `theme_void()` 的图（Venn / 网络快照）→ 长出数值轴（`34_variance_partition`、`82` 都踩过）；
- `ggplotify` 包装的 grob 图（`0-1` 坐标系 + 空面板）→ 长出 `0.00–1.00` 的 x/y 轴与轴题
  （`83_taxon_diff` 的热图踩过）。
**处置**：单面板图**不要**为了统一而走 `md_compose()`，直接 `md_save_figure(CFG, p, ...)`；
必须留在 patchwork 里的空图/装饰图，显式 `scale_x_continuous(breaks = NULL)` +
`scale_y_continuous(breaks = NULL)`。

### C33 · `geom_jitter()` **没有** `seed` 参数 —— 抖动点定种子要用 `position_jitter()`
`geom_jitter(seed = 123)` 不报错，只会在日志里留一句 **"Ignoring unknown parameters: `seed`"**，
然后**照旧从全局 RNG 抽位置**。于是"同一张图在全流程里跑与单跑不一样"，
而且定点种子的修复看起来"生效了"（因为大多数模块在画图前已经 `md_seed()` 过）。
正解：
```r
ggplot2::geom_jitter(size = ..., alpha = ...,
  position = ggplot2::position_jitter(width = w, seed = MD_JITTER_SEED))
```
`MD_JITTER_SEED` 定义在 `R/shared/12_plot_theme.R`，是**渲染常量**（不承载统计信息），
与统计用的 `md_seed(CFG, md_offset(...))` 两回事，不要混。

### C34 · 宽幅图在屏幕上"看起来字小"——统一到 183 mm
`80_network_nodes` 原本按原文出 **533.4 × 177.8 mm**（21×7 in）单排三格。
字号是绝对 pt，与 183 mm 的图**一模一样**；但看图时按屏宽缩放，整幅被压到 1/3 大小，
于是"字和点都比别的图小"。这不是字号配置错，是**画幅**的问题。
**处置**：默认改成 183 mm 宽、三格纵向排列（`design = "A\nB\nC"`），
把原文的宽幅单排留成开关 `CFG$network$style$wide_layout`。

### C35 · 拼图里"某一个面板不要坐标轴"时，`pw & theme(...)` 会把它又打开
网络快照用 `element_blank()` 关掉了 `axis.text` / `axis.ticks` / `panel.grid`，
但 `local80_compose()` 结尾统一执行 `pw & theme_md_comp(CFG)` —— 这条**后来居上**，
把轴又画了回来（`80a` 的 A 面板长了 `-1.0…1.0` 的坐标轴）。
`pw & theme()` 无法"只作用于某一个子图"，而 `[[<-` 改子图的写法依赖 patchwork 版本。
**处置**：给 `local80_compose()` 加 `theme_comp` 开关。
整张图里有面板不需要坐标轴时置 `FALSE`，然后**只给真正需要缩字号的那个面板**
单独 `p + theme_md_comp(CFG)`（`80a` 的观测-随机面板、`80c` 的两个面板这样做）。
判据很简单：**看 PDF 文字层里有没有孤零零的数字刻度**。

### C36 · 自绘图例带必须把"每一类"写出来
`80b` 的门色图例原来只画了一排彩色圆点、**没有写门名**（读者只看到色块）；
相关方向那行也没标。自绘图例（`theme_void()` + `geom_*`）不像 `scale_*_manual()`
那样自动带标签，**每个类别都要显式 `geom_text()` 画名字**，类别多时折行。
修好后 `80b` 的文字层能直接抽到
`Phylum Actinobacteria Bacteria Cyanobacteria Fibrobacteres Firmicutes Proteobacteria`，
这就是验收判据。

### C37 · 子图格数别写死 —— 用 `ceiling(sqrt(n))`
`local80_panel_groups()` 原来写死 `ncol = 4`。演示数据有 6 个分组网络，
于是排成 4 + 2，第二行只有两格、右边一大片空白。
**处置**：`ncol <- max(1L, min(4L, ceiling(sqrt(length(keys)))))` ——
4/6/8/9/12 个网络分别得到 2/3/3/3/4 列，都接近方形。

### C38 · `plot_layout(design = )` 的行**默认等分** —— 图例带会吃掉一半高度
`80b` 用 `design = "AAAA\nLLLL"`（上面是网络区，下面是图例带），
patchwork 把两行**等分**，于是一条本来就又扁又宽的图例带占掉了 50% 的画幅，
把 2×3 的网络格挤成一半高（用户报的"图例版幅过大"）。
**处置**：`plot_layout(design = ..., heights = c(7, 1))` 显式给行权重。
实测图例带从 50% 降到约 **12%**（与 `1/8` 的比例吻合）。
同一条也适用于 `widths`。判据：把 PDF 里"图例文字的最小 y"除以页面高度。

### C39 · `grid::grid.grabExpr()` 会在包根掉一个 `Rplots.pdf`
现象：`--force` 跑完，包根多出一个 3 KB 的 `Rplots.pdf`（mtime = 运行结束那一刻）。
根因：`grid.grabExpr()` 在**当前没有图形设备**时会打开 R 的**默认设备**；
Rscript 下默认设备是 `pdf()`，文件名就是当前目录的 `Rplots.pdf`。
`82`（`ComplexHeatmap::UpSet`）与 `83`（热图）都走这条路。
实测两种"看起来对"的修法**都无效**：
- `grid.grabExpr(expr, device = function(...) pdf(NULL, ...))` —— 该参数在此 R 版本下不生效；
- 只写 `pdf(NULL)` 而后不加清理 —— 设备泄漏。

**处置**：统一走 `md_grid_grab()`（`R/shared/12_plot_theme.R`，紧挨 `md_save_figure`）：
```r
gg <- ggplotify::as.ggplot(md_grid_grab(ComplexHeatmap::draw(ht, newpage = TRUE)))
```
它先判断有没有设备，没有就开一个 `pdf(NULL)`（**不写文件**），抓完把自己开的关掉；
调用方本来就开着设备（例如 `md_save_figure` 的 cairo_pdf）时则完全不碰。
**验收**：`--force` 之后 `Test-Path Rplots.pdf` 必须为假。
`.gitignore` 里也已加上 `Rplots.pdf` 兜底。

---

## D. 排查顺序（照这个顺序查，别乱试）

1. `logs/<label>_<run_id>.log` 里的 `[ERROR]` / `[SKIP]` / `[WARN]`
2. `logs/summary_<run_id>.csv`：哪个模块什么状态
3. `results/00_index.xlsx`：`Validation`（产物校验）/ `Run_summary` / `Skipped` / `Missing_packages`
4. `results/01_data_health.xlsx`：`Issues` sheet 列出全部 BLOCK/WARN 级问题
5. `results/05_assemble.xlsx`：`Filter_log` + `Normalization` + `Coverage`
6. 还没头绪 → `--only=<该模块>` 单跑，看完整输出
