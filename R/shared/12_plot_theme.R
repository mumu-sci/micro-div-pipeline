# =============================================================================
# R/shared/12_plot_theme.R —— 主题 / 字号 / 版式 / 出图 / 配色（**视觉唯一真身**）
# -----------------------------------------------------------------------------
# 库文件。本包所有图件必须经 md_save_figure() 落盘（统一 pdf、统一校验）。
#
# ⚠ 视觉口径来自 nifH 论文的 `R代码/11_shared.R` §0，**逐项对齐，不要自造**：
#     theme_md()  ↔ theme_isme()        （border 0.5 / axis.line / legend.title 加粗 / strip 透明）
#     MD_PT       ↔ S_TICK/S_TTL/S_LEG/S_ANN/S_TAG
#     MD_PT_COMP  ↔ S_*_COMP            （显式值，**不是**统一 0.8 倍）
#     MD_PT_SIZE  ↔ JIT_* / PT_* / S_PT
#     MD_BUILTIN_PAL ↔ cols_zone / cols_veg / cols_assembly / cols_2x2
#     回归线     ↔ REG_LINE_COL / REG_LINE_LW / REG_RIBBON_A
#
# 两条字号轨道（**nifH 论文踩过的坑，务必遵守**）：
#   element_text(size = )        单位是 pt   -> 直接用 MD_PT$*
#   geom_text / annotate(size = ) 单位是 mm   -> **必须** pt2gg(pt)
# 且 geom_text 的 size 在构造时就固化，事后 + theme() 改不动；
# 小面板必须把字号作为函数参数传进去。
#
# ⚠ 图中文字必须**纯 ASCII**（码点 <= 0x7E）。中文/特殊符号一律换成 ASCII 写法，
#   见 references/04-pitfalls.md B12；95_export_tables 有自动审计。
# =============================================================================

# ---- 字号（pt）—— 主图，对应 nifH 的 S_TICK / S_TTL / S_LEG / S_ANN / S_TAG ----
MD_PT <- list(
  tick   = 8,      # 坐标轴刻度        S_TICK
  title  = 9,      # 坐标轴标题        S_TTL
  legend = 7.5,    # 图例              S_LEG
  ann    = 7.5,    # 图内标注          S_ANN
  tag    = 10,     # (a)(b) 面板标签   S_TAG
  strip  = 7.5     # 分面标签
)

# ---- 字号（pt）—— 拼图专用，对应 nifH 的 S_*_COMP（**显式值，不是统一比例**）----
MD_PT_COMP <- list(
  tick   = 6.5,    # S_TICK_COMP
  title  = 7.5,    # S_TTL_COMP
  legend = 6,      # S_LEG_COMP
  ann    = 5.5,    # S_ANN_COMP
  tag    = 8,      # S_TAG_COMP
  strip  = 6       # S_LEG_COMP（nifH 的 strip.text 在拼图里用图例字号）
)

# 兼容旧字段（早期版本用统一比例）；新代码请用 MD_PT_COMP
MD_PT_COMP_RATIO <- 0.8

# ---- 线宽 / 点大小 —— 对应 nifH 的 LW_AXIS / BOX_LW_STD / JIT_* / PT_* / S_PT ----
MD_LW <- list(
  axis = 0.35,     # LW_AXIS
  box  = 0.4,      # BOX_LW_STD
  line = 0.4,      # REG_LINE_LW
  ref  = 0.3,
  bar  = 0.6
)

MD_PT_SIZE <- list(
  point        = 0.9,    # PT_SIZE_STD  回归图散点
  scatter      = 0.65,   # S_PT         一般散点（nifH：全稿散点缩小为原 1/2）
  pcoa         = 1.4667, # 2.2 * 2/3    PCoA 排序散点（02_Fig02.R §B/C 的写法）
  jitter       = 1.8,    # JIT_SIZE_STD 箱线图抖动点
  alpha        = 0.67,   # JIT_ALPHA_STD
  point_alpha  = 0.4,    # PT_ALPHA_STD 回归图散点透明度
  jitter_width = 0.15,   # JIT_WIDTH_STD
  stroke       = 0,      # 抖动点描边
  errorbar     = 0.35
)

# ---- 回归线与置信带 —— 对应 REG_LINE_COL / REG_LINE_LW / REG_RIBBON_A ----
MD_REG <- list(col = "#1E72B0", lw = 0.4, ribbon_alpha = 0.1)

MD_FONT <- "Arial"

# 抖动点的**渲染**种子（不是统计种子）。
# 为什么需要它：`geom_jitter()` 默认从全局 RNG 抽位置，于是**同一张图在全流程里跑
# 与单跑会不一样**（前面的模块消耗掉的随机数不同）。抖动位置纯属呈现，不承载任何
# 统计信息，所以固定成常量最合理 —— 这样"同一份数据 + 同一份配置"画出来的图逐像素一致。
# 统计随机过程仍然必须走 `md_seed(CFG, md_offset(...))`，两者不要混。
MD_JITTER_SEED <- 20260926L

#' pt -> ggplot geom_text 的 size（mm）
pt2gg <- function(pt) pt * 0.352778

#' 取字号表。comp = TRUE 时返回拼图用的显式字号。
md_pt <- function(CFG, comp = FALSE) {
  if (comp) return(MD_PT_COMP)
  p <- MD_PT
  if (!is.null(CFG) && !is.null(CFG$figure$base_pt)) {
    # base_pt 只按比例微调整套主图字号（默认 9 与 MD_PT$title 一致，等于不变）
    k <- CFG$figure$base_pt / MD_PT$title
    if (is.finite(k) && k > 0) p <- lapply(p, function(x) x * k)
  }
  p
}

# ---------------------------------------------------------------------------
# 主题（逐项对齐 nifH theme_isme）
# ---------------------------------------------------------------------------
#' @param base_size 基准字号（pt）
#' @param comp      是否用拼图字号
#' @param CFG       传进来可让 config 覆盖字体（figure$font / figure$cjk_font）
theme_md <- function(base_size = MD_PT$title, base_family = MD_FONT,
                     comp = FALSE, CFG = NULL) {
  # 容错：允许 theme_md(CFG) 这种位置调用（把整个 CFG 当 base_size 传进来）
  if (is.list(base_size)) { CFG <- base_size; base_size <- MD_PT$title }
  if (!is.numeric(base_size)) {
    stop("theme_md(): base_size 必须是数值（pt）。",
         "如果你想把 CFG 传进来，请用 theme_md(CFG = CFG)。", call. = FALSE)
  }
  if (!is.null(CFG)) {
    base_family <- CFG$figure$font %||% base_family
    if (!is.null(CFG$figure$cjk_font)) base_family <- CFG$figure$cjk_font
  }
  p <- md_pt(CFG, comp)
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # —— 对齐 theme_isme() ——
      text              = ggplot2::element_text(family = base_family),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(colour = "black", linewidth = 0.5, fill = NA),
      panel.background  = ggplot2::element_rect(fill = "white"),
      axis.line         = ggplot2::element_line(linewidth = MD_LW$axis, colour = "black"),
      axis.text         = ggplot2::element_text(size = p$tick, colour = "black"),
      axis.title        = ggplot2::element_text(size = p$title, colour = "black"),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = MD_LW$axis),
      axis.ticks.length = grid::unit(0.15, "cm"),
      legend.title      = ggplot2::element_text(size = p$legend, colour = "black", face = "bold"),
      legend.text       = ggplot2::element_text(size = p$legend, colour = "black"),
      legend.key.size   = grid::unit(0.35, "cm"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      # **图例默认放图下面、横向排列**（用户口径；原稿多数面板也是 bottom）
      # 需要放右边/隐藏的模块请显式覆盖（如 70 用 md_compose(legend = "none")）
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.box        = "vertical",
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(0.1, "cm"),
      legend.spacing    = grid::unit(0.05, "cm"),
      strip.text        = ggplot2::element_text(size = p$strip, face = "bold"),
      strip.background  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(size = p$title, face = "bold", hjust = 0.5),
      plot.tag          = ggplot2::element_text(size = p$tag, face = "bold"),
      plot.margin       = ggplot2::margin(4, 4, 4, 4, unit = "pt")
    )
}

#' 拼图专用：套用显式拼图字号（对应 nifH theme_comp()）
theme_md_comp <- function(CFG = NULL) {
  p <- MD_PT_COMP
  ggplot2::theme(
    axis.text    = ggplot2::element_text(size = p$tick),
    axis.title   = ggplot2::element_text(size = p$title),
    legend.title = ggplot2::element_text(size = p$legend, colour = "black", face = "bold"),
    legend.text  = ggplot2::element_text(size = p$legend, colour = "black"),
    legend.key.size = grid::unit(0.28, "cm"),
    # 拼图层同样默认底部 —— patchwork 的 guides="collect" 是按**合成层主题**放图例的，
    # 这里不设就会落到 ggplot 默认的 right，把子面板自己设的 bottom 顶掉。
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box      = "vertical",
    legend.margin   = ggplot2::margin(0, 0, 0, 0),
    strip.text   = ggplot2::element_text(size = p$strip, face = "bold"),
    plot.tag     = ggplot2::element_text(size = p$tag, face = "bold"),
    plot.margin  = ggplot2::margin(2, 2, 2, 2)
  )
}

#' 图例列数自适应（防溢出 —— nifH 的 FigS3 图例被裁过）
#' @param nrow 给了就按行排（nifH 的 PCoA 图常用 `guide_legend(nrow = 1)`）；
#'   给了 nrow 时忽略 ncol。
md_guide_flex <- function(n_items, position = "bottom", bottom_max = 4, nrow = NULL) {
  if (!is.null(nrow)) return(ggplot2::guide_legend(nrow = max(1L, as.integer(nrow))))
  if (identical(position, "bottom")) {
    ncol <- if (n_items <= bottom_max) n_items else bottom_max
  } else {
    ncol <- if (n_items <= 10) 1 else 2
  }
  ggplot2::guide_legend(ncol = max(1L, ncol), byrow = TRUE)
}

# ---------------------------------------------------------------------------
# 配色 —— 内置具名色板，逐字照搬 nifH 11_shared.R
# ---------------------------------------------------------------------------
MD_BUILTIN_PAL <- list(
  # cols_zone（11_shared.R:162-167）
  zone = c(Coastal = "#89C4D7", Nearshore = "#4DAFB5",
           Xisha = "#2B6AA0", Nansha = "#3B2C6B"),
  # cols_veg（11_shared.R:172-177）
  veg = c(Bare = "#989898", Herb = "#959BB0",
          Shrub = "#60859E", Tree = "#B98079"),
  # cols_assembly（11_shared.R:259-265）—— 逐字照搬，注意 nifH 本身有两组同色
  assembly = c("Homogeneous selection"  = "#80B1D2",
               "Variable selection"     = "#80B1D2",
               "Dispersal limitation"   = "#C7CEEA",
               "Homogenizing dispersal" = "#F6AAA1",
               "Drift"                  = "#C7CEEA"),
  # cols_2x2（11_shared.R:823-826）
  quadrant = c(Core = "#D73027", Permanently_rare = "#999999",
               Conditionally_rare = "#E69F00", Transient_rare = "#4575B4"),
  # 02_Fig02.R §D 周转 / 嵌套两色（原稿写死常量）
  beta = c(turnover = "#5B93C4", nestedness = "#E0877C"),
  # 02_Fig02.R §E 的 cols_ep4（γ 分解：分类学浅色族 / 系统发育深色族）
  gamma = c(taxonomic_alpha = "#AECFE6", taxonomic_turnover = "#F7CCC5",
            phylogenetic_alpha = "#5B93C4", phylogenetic_turnover = "#E0877C"),
  # 03_Fig03_assembly.R §E 变差分解的分块色
  vpa = c(geo = "#80B1D2", soil = "#F6AAA1", climate = "#C7CEEA")
)

# 08_Supp_diversity.R:277-282 的属配色（20 色硬编码 + Others 灰）
MD_TAXON_GENUS20 <- c(
  "#FAF5B5", "#F9B063", "#F7CBDF", "#FF9AA2", "#FFB7B2", "#FAD0C4", "#FDE2E4", "#FBE7C6",
  "#FFDFD3", "#E6C9A8", "#F8E9A1", "#E2F0CB", "#B5EAD7", "#A8E6CF", "#D4F0F0", "#A0E7E5",
  "#B9CCED", "#DCD6F7", "#E0BBE4", "#957DAD"
)
MD_TAXON_OTHER   <- "#D7D7D5"   # 属的 "Others"
MD_TAXON_UNCLASS <- "#D3D3D3"   # 门的 "Unclassified"
MD_TAXON_PHYLUM_OTHER <- "#A9A9A9"  # 门的 "Other"

#' 分类单元配色 —— **逐项照抄 nifH 原脚本**，不是自造色板
#'
#' · `Phylum`（`01_Fig01.R:254-256`）：
#'   `colorRampPalette(brewer.pal(12,"Set3"))(n)`，names 按**丰度降序**的 top_phyla 给，
#'   再加 `Unclassified = #D3D3D3`、`Other = #A9A9A9`。
#' · `Genus`（`08_Supp_diversity.R:277-283`）：
#'   20 色硬编码，names = `c(rev(top20_genera), "Others")` —— 即**升序**分配，
#'   丰度最低的那个属拿第一色；其余归 "Others"（`#D7D7D5`）。
#' · 其它层级原稿没有 → 退回 `md_pal()` 定性板并 `md_warn`，**不自造**。
#'
#' @param taxa 已按丰度**降序**排好的 top 类群名（本包 `top_taxa` 就是这个序）
#' @param color_order 用于**给色板命名**的类群顺序。门的原脚本有个关键细节：
#'   `phylum_colors` 的名字来自 `top_phyla`（**按 ASV 数量降序**，`01_Fig01.R:217`），
#'   而面板的堆叠顺序来自 `top12_phyla`（**按平均丰度降序**，`01_Fig01.R:369`）；
#'   `scale_fill_manual` 是**按名字**取色的，所以实际效果是「数量序配色 + 丰度序堆叠」。
#'   要逐色复刻就必须把 `color_order` 传成 ASV 数量降序的类群名。
md_taxon_pal <- function(CFG, rank, taxa, other_name = "Other", color_order = NULL) {
  ov <- CFG$figure$palette[[paste0("taxon_", rank)]]
  if (!is.null(ov)) return(ov)
  taxa <- as.character(taxa)
  n <- length(taxa)
  if (identical(rank, "Phylum")) {
    ramp <- if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
    } else {
      md_warn(NULL, "缺 RColorBrewer，门的 Set3 色板退回定性板"); NULL
    }
    if (is.null(ramp)) {
      pal <- stats::setNames(md_pal(NULL, taxa, name = "taxon"), taxa)
    } else {
      ord <- if (!is.null(color_order)) intersect(as.character(color_order), taxa) else taxa
      if (!length(ord)) ord <- taxa
      ramp_v <- stats::setNames(ramp(length(ord)), ord)
      pal <- ramp_v[taxa]
      if (anyNA(pal)) {   # 不在数量序里的类群兜底色（原脚本这类会是 NA -> 默认灰）
        fb <- md_pal(NULL, taxa, name = "taxon")
        pal[is.na(pal)] <- fb[is.na(pal)]
      }
    }
    pal[taxa == "Unclassified"] <- MD_TAXON_UNCLASS
    pal[taxa == "Other"]        <- MD_TAXON_PHYLUM_OTHER
    if (other_name %in% names(pal)) pal[other_name] <- MD_TAXON_PHYLUM_OTHER
    return(pal)
  }
  if (identical(rank, "Genus")) {
    # 原脚本：names(cols_genus) <- c(rev(top20_genera), "Others")  → 升序分配
    idx <- seq_len(min(n, length(MD_TAXON_GENUS20)))
    pal <- stats::setNames(MD_TAXON_GENUS20[idx], rev(taxa)[idx])
    if (other_name %in% taxa || TRUE) pal[other_name] <- MD_TAXON_OTHER
    return(pal)
  }
  md_warn(NULL, "分类层级 '", rank, "' 在原代码里没有对应色板，退回定性板（可用 ",
          "CFG$figure$palette$taxon_", rank, " 指定）")
  md_pal(NULL, taxa, name = paste0("taxon_", rank))
}

# 通用兜底（非 nifH 具名因子的场合）
MD_PAL_DEFAULT <- list(
  diverging  = c("#2166AC", "#F5F5F5", "#B2182B"),
  sequential = c("#F7FBFF", "#4292C6", "#08306B"),
  qualitative = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
                  "#86BCB6", "#D37295", "#A0CBE8", "#FFBE7D"),
  # 地形高程渐变（nifH 07_S1_study_area_map.R §4 逐字两色）
  elevation = c("#B9D7A8", "#C9AE7E"),
  sem_pos = "#2166AC", sem_neg = "#B2182B",
  # ---- 2026 毕业论文补充轮 ----
  # 共现网络正/负边（毕业论文 5_network_structure.R:169 注释："正=steelblue, 负=indianred"）
  net_pos = "#4682B4", net_neg = "#CD5C5C",
  # Zi-Pi 四类节点角色（毕业论文 5_network_structure.R:1064-1065 逐字）
  zi_pi = c("Peripheral" = "#BDBDBD", "Connector" = "#4DAFB5",
            "Module Hub" = "#D73027", "Network Hub" = "#B2182B")
)

# ---------------------------------------------------------------------------
# 共现网络的样式常量（模块 80）
#   逐项照抄毕业论文 5_network_structure.R；模块不得自造网络配色/尺寸。
#   布局（layout）走 CFG$network$style$layout —— 它只影响节点坐标，不影响数值。
# ---------------------------------------------------------------------------
MD_NET <- list(
  node_size      = 1.6,      # 一般节点
  node_size_hub  = 2.6,      # hub / keystone 节点
  node_alpha     = 0.75,
  edge_width     = 0.18,
  edge_alpha     = 0.45,
  label_size_pt  = 5,        # 节点标签字号（pt）
  edge_col_pos   = "#4682B4",  # steelblue
  edge_col_neg   = "#CD5C5C",  # indianred
  layout         = "graphopt", # 全局网络布局（原脚本 layout_with_graphopt）
  layout_group   = "fr",       # 分组网络布局（原脚本 layout_with_fr, niter = 5000）
  fr_niter       = 5000,
  # 节点属性图（80d）的版式：FALSE（默认）= 183 mm 宽、三格纵向排列，
  # 与全包其余模块的"一图一号字梯"一致；TRUE = 复刻原文 Fig5a 的 533.4 x 177.8 mm 单排。
  wide_layout    = FALSE,
  role_colors    = c("Peripheral" = "#BDBDBD", "Connector" = "#4DAFB5",
                     "Module Hub" = "#D73027", "Network Hub" = "#B2182B"),
  zi_thr         = 2.5,      # Deng et al. 2012 阈值
  pi_thr         = 0.62,
  # 随机网络比对的配色（**逐字照抄 5_network_structure.R:1750-1752**）
  random_cols    = c("Higher than random" = "#E31A1C",
                     "Lower than random"  = "#2166AC",
                     "Not significant"    = "#999999"),
  # 节点周转 Sankey（ggalluvial）各"细胞类型"配色沿用 2x2 四象限
  cell_labels    = c(Core               = "Core",
                     Permanently_rare   = "Permanently rare",
                     Conditionally_rare = "Conditionally rare",
                     Transient_rare     = "Transient rare")
)

#' 取网络样式（CFG 覆盖优先）
md_net_cfg <- function(CFG) {
  st <- CFG$network$style %||% list()
  utils::modifyList(MD_NET, st)
}

#' 按因子名 / 水平名自动挑选内置具名色板；挑不到返回 NULL
#'   name 命中 zone/region -> zone 板；命中 veg/habitat -> veg 板
#'   否则看 levels 是否被某个板完全覆盖
md_pick_board <- function(name = NULL, levels = NULL) {
  nm <- if (is.null(name)) "" else tolower(as.character(name)[1])
  if (nzchar(nm)) {
    if (grepl("zone|region|island|site_group|landform", nm)) return(MD_BUILTIN_PAL$zone)
    if (grepl("veg|habitat|plant|cover", nm))              return(MD_BUILTIN_PAL$veg)
  }
  lv <- as.character(levels %||% character(0))
  if (length(lv)) {
    for (b in MD_BUILTIN_PAL) {
      if (all(lv %in% names(b))) return(b)
    }
  }
  NULL
}

#' 取具名配色向量。查找顺序：
#'   1) CFG$figure$palette[[name]]   （用户显式覆盖）
#'   2) 内置具名板（按 name 或 levels 命中）
#'   3) 定性色板兜底
#' 缺色时 **warning 提示**并补兜底色，绝不静默换整套。
md_pal <- function(CFG, levels, name = NULL, n = NULL) {
  levels <- as.character(levels)
  # 1) config 覆盖
  if (!is.null(name) && !is.null(CFG$figure$palette[[name]])) {
    p <- CFG$figure$palette[[name]]
    if (!is.null(names(p))) {
      out <- p[levels]
      if (anyNA(out)) {
        fb <- md_pick_board(name, levels) %||% MD_PAL_DEFAULT$qualitative
        miss <- levels[is.na(out)]
        out[is.na(out)] <- unname(fb)[seq_along(miss)][seq_along(miss)]
        md_warn(CFG, "md_pal(): 自定义色板缺少水平 ", paste(miss, collapse = ", "),
                " —— 已用兜底色补齐")
      }
      return(out)
    }
    return(stats::setNames(rep(p, length.out = length(levels)), levels))
  }
  # 2) 内置具名板
  board <- md_pick_board(name, levels)
  if (!is.null(board)) {
    out <- board[levels]
    if (anyNA(out)) {
      miss <- levels[is.na(out)]
      out[is.na(out)] <- unname(MD_PAL_DEFAULT$qualitative)[seq_along(miss)]
      md_warn(CFG, "md_pal(): 内置色板 ", paste(names(board)[1], "…"), " 不含水平 ",
              paste(miss, collapse = ", "), " —— 已用兜底色补齐")
    }
    return(out)
  }
  # 3) 兜底
  n <- n %||% length(levels)
  stats::setNames(unname(MD_PAL_DEFAULT$qualitative)[seq_len(max(n, 1))][seq_along(levels)],
                  levels)
}

#' 连续色标
md_pal_cont <- function(CFG, which = "sequential") {
  CFG$figure$palette[[which]] %||% MD_PAL_DEFAULT[[which]]
}

# ---------------------------------------------------------------------------
# 箱线图样式：默认完全复刻 nifH 的 plot_dual_factor_panel，但每项可覆盖
# ---------------------------------------------------------------------------
MD_BOX_DEFAULT <- list(
  style          = "nifh",      # nifh | simple
  color_by       = "auto",      # auto(双因子=f2, 单因子=f1) | f1 | f2 | none
  fill           = NA,          # NA = 空心箱体（nifH 用 fill=NA + 边框上色）
  # 是否画 CLD 显著性字母。**默认 TRUE = 与原代码构建方法一致**；
  # 以后写论文不想要字母时，设 CFG$plot$box$letters = FALSE 即可（只影响出图，不影响表）。
  letters        = TRUE,
  letters_at     = 1.08,        # 字母 y = 每组 max(y) * 该系数（nifH 口径）
  letters_face   = "bold",
  letters_colour = "grey20",
  show_n         = FALSE,       # 默认不标样本量（nifH 没有这个标注）
  y_mult         = 0.05,        # y scale 的 expansion(mult=)：nifH 用 ggplot2 默认 0.05
                                #   （05_Fig05_decoupling.R §A 显式写 c(0.05,0.05)）
  kw_annot       = TRUE,        # 右上角 KW 注释
  y_from_zero    = TRUE,        # y 轴从 0 起
  y_expand       = 1.35,        # y 上限 = max * 该系数
  inward_ticks   = TRUE,        # 刻度朝内（nifH 用 unit(-0.12,"cm")）
  inter_sep      = "_",         # 双因子组合分隔符，务必与生成 letters 的一致
  x_labels       = "two_line"   # two_line | one_line
)

#' 取箱线图样式（默认值 <- CFG$plot$box <- 本次调用的显式覆盖）
md_box_cfg <- function(...) {
  b <- utils::modifyList(MD_BOX_DEFAULT, .md_state$box %||% list())
  ov <- list(...)
  ov <- ov[!vapply(ov, is.null, logical(1))]
  if (length(ov)) b <- utils::modifyList(b, ov)
  b
}

# ---- 显著性标记与 P 值文案（对齐 nifH 的 sig_label / kw_fmt，**小写 p + 星号**）----
#' 星号标记；p<0.05 -> *, <0.01 -> **, <0.001 -> ***, 其余 ns
md_sig_label <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

#' "p = 0.033 ***" / "p < 0.001 ***"（**纯 ASCII**）
md_kw_fmt <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("p = NA")
  paste0(if (p < 0.001) "p < 0.001" else paste0("p = ", format(p, digits = 3)),
         " ", md_sig_label(p))
}

#' 组装面板右上角的 KW 注释文本，如 "Zone: p = 0.033 *; Vegetation: p = 0.486 ns"
#' @param ps 具名数值向量或具名 list（名字 = 因子名，值 = P）
md_kw_text <- function(ps) {
  if (is.null(ps) || length(ps) == 0) return(NULL)
  ps <- ps[!vapply(ps, function(x) is.null(x) || (length(x) == 1 && is.na(x)), logical(1))]
  if (length(ps) == 0) return(NULL)
  paste(vapply(names(ps), function(nm) paste0(nm, ": ", md_kw_fmt(ps[[nm]])),
               character(1)), collapse = "; ")
}

# ---------------------------------------------------------------------------
# 版式：design 字符串 -> 每个面板的 mm 尺寸
# ---------------------------------------------------------------------------
#' @param design 形如 "AABB\nAACC" 的字符矩阵；或 NULL（用 n_row/n_col 均分）
#' @param widths,heights 面板宽/高的相对比例
#' @return data.frame(panel, row, col, width_mm, height_mm)
md_panel_dims <- function(design = NULL, total_w = 183, total_h = 110,
                          widths = NULL, heights = NULL, n_row = NULL, n_col = NULL,
                          gap_mm = 1.5) {
  if (!is.null(design)) {
    rows <- strsplit(design, "\n")[[1]]
    nr <- length(rows); nc <- max(nchar(rows))
    mat <- do.call(rbind, lapply(rows, function(r) {
      s <- strsplit(r, "")[[1]]
      length(s) <- nc; s
    }))
    panels <- sort(unique(as.vector(mat[!is.na(mat)])))
    col_of <- function(ch) range(which(mat == ch, arr.ind = TRUE)[, "col"])
    row_of <- function(ch) range(which(mat == ch, arr.ind = TRUE)[, "row"])
  } else {
    nr <- n_row %||% 1L; nc <- n_col %||% 1L
    letters_all <- LETTERS[seq_len(nr * nc)]
    mat <- matrix(letters_all, nr, nc, byrow = TRUE)
    panels <- letters_all
    col_of <- function(ch) { w <- which(mat == ch, arr.ind = TRUE); range(w[, "col"]) }
    row_of <- function(ch) { w <- which(mat == ch, arr.ind = TRUE); range(w[, "row"]) }
  }
  w_rel <- if (is.null(widths)) rep(1, nc) else rep(widths, length.out = nc)
  h_rel <- if (is.null(heights)) rep(1, nr) else rep(heights, length.out = nr)
  col_mm <- (total_w - gap_mm * (nc - 1)) * w_rel / sum(w_rel)
  row_mm <- (total_h - gap_mm * (nr - 1)) * h_rel / sum(h_rel)
  do.call(rbind, lapply(panels, function(ch) {
    cc <- col_of(ch); rr <- row_of(ch)
    data.frame(panel = ch, row = rr[1], col = cc[1],
               width_mm = sum(col_mm[cc[1]:cc[2]]),
               height_mm = sum(row_mm[rr[1]:rr[2]]),
               stringsAsFactors = FALSE)
  }))
}

#' 取某个面板的 mm 尺寸
md_get_panel_dim <- function(CFG, fig_id, panel) {
  spec <- CFG$figure$panels[[fig_id]]
  if (is.null(spec)) return(c(CFG$figure$width_double, CFG$figure$height_default))
  d <- md_panel_dims(spec$design, spec$total_w %||% CFG$figure$width_double,
                     spec$total_h %||% CFG$figure$height_default,
                     spec$widths, spec$heights, spec$n_row, spec$n_col)
  r <- d[d$panel == panel, ]
  if (nrow(r) == 0) return(c(CFG$figure$width_double, CFG$figure$height_default))
  c(r$width_mm[1], r$height_mm[1])
}

# ---------------------------------------------------------------------------
# 出图开关：哪些模块默认不出图
# ---------------------------------------------------------------------------
#' 模块是否允许出图。CFG$figure$skip_modules 里的模块只出结果表、不出图。
md_figure_on <- function(CFG, module) {
  skip <- CFG$figure$skip_modules %||% character(0)
  !(as.character(module) %in% as.character(skip))
}

#' 从图件名（"11_alpha"）里取模块号（"11"）
md_module_of <- function(name) sub("^([0-9]+).*$", "\\1", as.character(name))

# ---------------------------------------------------------------------------
# 出图
# ---------------------------------------------------------------------------
#' 在**不落盘**的图形设备上执行一段"会绘图"的表达式。
#'
#' 背景：Rscript 下**没有**图形设备，此时任何 base-R / grid 的绘图动作都会打开
#' R 的**默认设备**，也就是当前目录的 `Rplots.pdf`（3 KB 的散件，每次 `--force`
#' 都会掉一个）。已知两类：
#'   * `grid::grid.grabExpr()`（82 的 UpSet、83 的热图）
#'   * **`vegan::rarecurve()` 即使 `tidy = TRUE` 也会开设备**（11 的稀释曲线）
#' 还有一些只是为了"取数据"而不得不触发的绘图副作用。
#'
#' 本函数：当前没设备就先开一个 `pdf(NULL)`（**不写任何文件**），跑完把
#' **自己开的**设备关掉；调用方本来就开着设备时完全不碰。
#' 见 `references/04-pitfalls.md` C39。
#'
#' @param expr 要执行的表达式（其返回值就是本函数的返回值）
md_offscreen <- function(expr) {
  had <- length(grDevices::dev.list()) > 0L
  if (!had) grDevices::pdf(NULL)
  on.exit({
    if (!had) while (length(grDevices::dev.list()))
      try(grDevices::dev.off(), silent = TRUE)
  }, add = TRUE)
  expr
}

#' 抓取一段 grid 绘图表达式，且**保证不在包根掉 Rplots.pdf**。
#'
#' `grid::grid.grabExpr()` 在**当前没有图形设备**时会打开 R 的**默认设备**。
#' 实测 `grid.grabExpr(expr, device = function(...) pdf(NULL, ...))` **无效**
#' （该参数在此 R 版本下不生效）；先开一个 null 设备则有效 —— 所以走 `md_offscreen()`。
#'
#' @param expr 一段 grid 绘图表达式（如 `ComplexHeatmap::draw(ht, newpage = TRUE)`）
#' @return `grid::grid.grabExpr()` 的返回值（grob / gTree）
md_grid_grab <- function(expr) {
  md_offscreen(grid::grid.grabExpr(expr))
}

#' 统一落盘。写 figures/<name>.pdf 并登记，供 95_export_tables 后置校验。
#' 若该模块在 CFG$figure$skip_modules 里，则跳过错图（只打日志，不算失败）。
md_save_figure <- function(CFG, plot, name, width = NULL, height = NULL,
                           singles = NULL) {
  mod <- md_module_of(name)
  if (!md_figure_on(CFG, mod)) {
    md_log(CFG, "SKIP-FIG", name, " —— 模块 ", mod, " 在 figure$skip_modules 中（只出结果表）")
    return(invisible(NULL))
  }
  if (!md_need(CFG, "core", "figure")) return(invisible(NULL))
  width  <- width  %||% CFG$figure$width_double
  height <- height %||% CFG$figure$height_default
  singles <- singles %||% CFG$figure$save_singles
  p <- md_path(CFG, "figures", paste0(name, ".pdf"))
  md_assert_not_data(CFG, p)
  ok <- tryCatch({
    ggplot2::ggsave(p, plot = plot, width = width, height = height, units = "mm",
                    device = grDevices::cairo_pdf, limitsize = FALSE)
    TRUE
  }, error = function(e) {
    md_warn(CFG, "出图失败 ", name, ": ", conditionMessage(e)); FALSE
  })
  if (ok) {
    sz <- file.info(p)$size
    if (is.na(sz) || sz == 0) md_warn(CFG, "图件为 0 字节: ", basename(p))
    else md_log(CFG, "OK", "出图 ", basename(p), " (", width, "x", height, " mm, ",
                round(sz / 1024), " KB)")
    .md_state$figures <- c(.md_state$figures %||% character(0), p)
  }
  invisible(p)
}

#' 一组面板拼成一张图。plots 为具名 list。design 里的字母 = 面板顺序。
#' @param legend 传给所有子面板的 legend.position。
#'   **默认 NULL = 不动子面板自己的设置**；非 NULL 时会覆盖子图已设的 legend.position。
md_compose <- function(CFG, plots, design = NULL, ncol = NULL, nrow = NULL,
                       guides = "collect", legend = NULL,
                       tag_levels = "A") {
  plots <- Filter(function(x) inherits(x, "gg") || inherits(x, "patchwork"), plots)
  if (length(plots) == 0) stop("md_compose(): 没有可拼的面板", call. = FALSE)
  pw <- if (!is.null(design)) {
    Reduce(`+`, plots) + patchwork::plot_layout(design = design, guides = guides)
  } else {
    patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow, guides = guides)
  }
  pw <- pw + patchwork::plot_annotation(tag_levels = tag_levels)
  pw <- pw & theme_md_comp(CFG)
  if (!is.null(legend)) pw <- pw & ggplot2::theme(legend.position = legend)
  pw
}

#' 单面板落盘（受 figure$save_singles 控制，默认关）
md_maybe_single <- function(CFG, plot, name, fig_id = NULL, panel = NULL) {
  if (!isTRUE(CFG$figure$save_singles)) return(invisible(NULL))
  w <- h <- NULL
  if (!is.null(fig_id) && !is.null(panel)) {
    d <- md_get_panel_dim(CFG, fig_id, panel); w <- d[1]; h <- d[2]
  }
  md_save_figure(CFG, plot, name, w, h)
}

md_figures_written <- function() .md_state$figures %||% character(0)
