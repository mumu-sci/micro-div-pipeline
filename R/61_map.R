# =============================================================================
# R/61_map.R ---- 研究区地图（主图 + 放大图 + 可选世界 inset）
# -----------------------------------------------------------------------------
# 契约见 _build-notes/micro-div-pipeline/MODULE_SPEC.md
#
# 做什么
#   面板 A  研究区主图：DEM 栅格（可选）+ 海岸线 + 采样点 + 站名标注
#   面板 B  局部放大图（CFG$map$zoom_bbox）
#   inset   世界小图（仅当 CFG$map$countries 非 NULL 时）
#
# 画法来源：07_S1_study_area_map.R §4（DEM 配色）/ §6（主图）/ §7（世界小图）/
#   §8（放大图）/ §9（拼图与尺寸）。原脚本只画 ne_coastline，没有国界图层，
#   故本模块也不取 ne_countries、不画国界。
#
# 关键设计：**不允许因联网失败而中断**
#   - DEM 走 cache/61_dem_z<zoom>.rds 缓存；已缓存则直接读，完全不联网
#   - CFG$map$offline = TRUE 时只用缓存，绝不联网
#   - 联网取 DEM 有超时保护（getOption("md.map.dem_timeout", 30) 秒）；
#     超时 / 报错 / 返回空 一律**降级为无 DEM 的散点图**，并把原因写进表与日志
#   - 海岸线数据（rnaturalearth）取不到时同样降级（只画点 + 网格）
#
# 参数来源（全部 CFG，未给的走 %||% 默认）
#   CFG$map$enabled / bbox / zoom_bbox / dem_zoom / dem_fact / dem_limits /
#   point_label / label_abbr / countries / offline / highlight
#
# 依赖：05（cache/05_core.rds）。坐标取自环境表（env），站点聚合后标点。
# =============================================================================

md_map <- function(CFG) {
  module <- "61"

  if (!isTRUE(CFG$map$enabled)) {
    md_skip(CFG, module, "map$enabled = FALSE（未开启研究区地图）")
    return(md_result(module, "SKIP", "map$enabled = FALSE"))
  }
  if (!md_need(CFG, "map", module))
    return(md_result(module, "SKIP", "缺少地图依赖包（sf/terra/rnaturalearth/ggspatial）"))
  if (!md_precheck(CFG, module))
    return(md_result(module, "SKIP", "前置条件不满足"))

  d <- tryCatch(md_core_analysis(CFG), error = function(e) NULL)
  if (is.null(d)) {
    md_skip(CFG, module, "取不到分析总表（cache/05_core.rds 缺失？）")
    return(md_result(module, "SKIP", "取不到分析总表"))
  }
  lat <- "Latitude"; lon <- "Longitude"
  if (!all(c(lat, lon) %in% names(d))) {
    md_skip(CFG, module, "数据里没有经纬度列（Latitude / Longitude）")
    return(md_result(module, "SKIP", "缺经纬度列"))
  }

  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  site_var <- CFG$design$site_var
  fac <- CFG$design$factors %||% character(0)
  fac <- fac[fac %in% names(d)]

  # ---------------------------------------------------------------------------
  # 1. 站点表（坐标去重；同一坐标的样本合并为一个点）
  # ---------------------------------------------------------------------------
  dd <- d[!is.na(d[[lat]]) & !is.na(d[[lon]]), , drop = FALSE]
  if (nrow(dd) == 0) {
    md_skip(CFG, module, "没有任何带经纬度的样本")
    return(md_result(module, "SKIP", "无坐标"))
  }
  if (is.null(site_var) || !site_var %in% names(dd)) {
    site_var <- ".site"
    dd[[site_var]] <- paste(round(dd[[lon]], 6), round(dd[[lat]], 6), sep = "_")
    add_note("design$site_var 不可用，站点按唯一坐标自动编号。")
  }
  grp <- as.character(dd[[site_var]])
  ord <- unique(grp)
  st <- data.frame(site = ord, stringsAsFactors = FALSE)
  st[[lat]] <- vapply(ord, function(g) dd[[lat]][grp == g][1], numeric(1))
  st[[lon]] <- vapply(ord, function(g) dd[[lon]][grp == g][1], numeric(1))
  for (f in fac) st[[f]] <- vapply(ord, function(g) as.character(dd[[f]][grp == g][1]), character(1))
  st$n_sample <- vapply(ord, function(g) sum(grp == g), integer(1))

  # ---------------------------------------------------------------------------
  # 2. 图面标签
  # ---------------------------------------------------------------------------
  lab_col <- "label"
  if (isTRUE(CFG$map$point_label)) {
    ab <- CFG$map$label_abbr
    st[[lab_col]] <- if (!is.null(ab) && !is.null(names(ab))) {
      out <- unname(ab[st$site]); out[is.na(out)] <- st$site[is.na(out)]; out
    } else st$site
  } else {
    st[[lab_col]] <- ""
  }

  # ---------------------------------------------------------------------------
  # 3. 边界 / 海岸线（取不到就降级）
  # ---------------------------------------------------------------------------
  bbox <- CFG$map$bbox %||% c(min(st[[lon]]) - 2, min(st[[lat]]) - 2,
                              max(st[[lon]]) + 2, max(st[[lat]]) + 2)
  bbox <- as.numeric(bbox)
  zoom_bbox <- CFG$map$zoom_bbox %||% bbox
  zoom_bbox <- as.numeric(zoom_bbox)

  # 原脚本 07_S1_study_area_map.R §6 只 geom_sf(ne_coastline)；它取过的
  # ne_countries（map_sf / map_coastal）是死代码，画面上没有国界图层 ——
  # 故这里同样只取海岸线，不再取国界（逐像素对齐原图）。
  coast_sf <- NULL
  ctry_arg <- CFG$map$countries
  coast_sf <- tryCatch(rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf"),
                       error = function(e) NULL)
  if (is.null(coast_sf)) {
    add_note("海岸线数据取不到（rnaturalearth 需要联网或本地数据缺失），",
             "主图降级为纯散点 + 坐标网格。")
    md_warn(CFG, "[61] 海岸线不可用，降级为无底图散点")
  }

  # ---------------------------------------------------------------------------
  # 4. DEM（缓存 + 超时 + 失败降级）
  # ---------------------------------------------------------------------------
  dem_df <- NULL
  dem_zoom <- CFG$map$dem_zoom %||% 6
  dem_fact <- CFG$map$dem_fact %||% 4
  dem_lim  <- CFG$map$dem_limits %||% c(0, 2000)
  dem_cache <- md_path(CFG, "cache", paste0("61_dem_z", dem_zoom, ".rds"))
  dem_reason <- NULL
  dem_src <- NULL

  if (isTRUE(CFG$map$offline)) {
    if (file.exists(dem_cache)) {
      dem_df <- tryCatch(readRDS(dem_cache), error = function(e) NULL)
      dem_src <- "cache"
    } else {
      dem_reason <- paste0("map$offline = TRUE 且缓存不存在（", basename(dem_cache), "）")
    }
  } else if (file.exists(dem_cache) && !md_recalc(CFG)) {
    dem_df <- tryCatch(readRDS(dem_cache), error = function(e) NULL)
    dem_src <- "cache"
    if (is.null(dem_df)) dem_reason <- "DEM 缓存存在但读取失败"
  } else {
    dem_timeout <- getOption("md.map.dem_timeout", 30)
    md_log(CFG, "INFO", "[61] 联网取 DEM（bbox 外扩，z = ", dem_zoom,
           "，超时保护 ", dem_timeout, " s）…")
    # 超时保护只覆盖【下载】这一步：setTimeLimit 一旦设成 transient=FALSE 就常驻，
    # 会把后续 terra 聚合也一起掐断（曾导致 DEM 已到手却出图失败）。
    # 因此下载完立刻复位，再用 elapsed 时间做判定。
    dem_elapsed <- NA_real_
    r <- tryCatch({
      bb <- sf::st_bbox(c(xmin = min(bbox[1], zoom_bbox[1]) - 1,
                          ymin = min(bbox[2], zoom_bbox[2]) - 1,
                          xmax = max(bbox[3], zoom_bbox[3]) + 1,
                          ymax = max(bbox[4], zoom_bbox[4]) + 1),
                        crs = sf::st_crs(4326))
      poly <- sf::st_as_sf(sf::st_as_sfc(bb))
      t0 <- Sys.time()
      setTimeLimit(cpu = Inf, elapsed = dem_timeout, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
      out <- elevatr::get_elev_raster(locations = poly, z = dem_zoom, clip = "locations")
      dem_elapsed <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      out
    }, error = function(e) { dem_reason <<- conditionMessage(e); NULL })

    if (!is.null(r)) {
      if (is.finite(dem_elapsed) && dem_elapsed > dem_timeout) {
        dem_reason <- paste0("取 DEM 用时 ", sprintf("%.1f", dem_elapsed),
                             " s，超过超时门限 ", dem_timeout, " s")
        r <- NULL
      } else {
        dem_df <- tryCatch({
          ra <- terra::rast(r)
          if (dem_fact > 1) ra <- terra::aggregate(ra, fact = dem_fact, fun = mean)
          xy <- terra::xyFromCell(ra, seq_len(terra::ncell(ra)))
          vv <- as.numeric(terra::values(ra))
          df <- data.frame(x = xy[, 1], y = xy[, 2], elev = vv)
          df <- df[!is.na(df$elev) & df$elev >= 0, , drop = FALSE]
          if (nrow(df) == 0) stop("DEM 无有效像元", call. = FALSE)
          saveRDS(df, dem_cache)
          df
        }, error = function(e) { dem_reason <<- conditionMessage(e); NULL })
        if (!is.null(dem_df)) dem_src <- "download"
      }
    }
  }

  if (is.null(dem_df)) {
    add_note("DEM 不可用（", dem_reason %||% "未知原因",
             "）-> 主图/放大图自动降级为**无 DEM 的散点图**（不报错）。")
    md_warn(CFG, "[61] DEM 不可用，降级为无 DEM 散点图：", dem_reason %||% "未知")
  } else {
    add_note("DEM 来源 = ", dem_src, "（z = ", dem_zoom, "，", dem_fact,
             "x", dem_fact, " 聚合，", nrow(dem_df), " 个有效像元，缓存 ",
             basename(dem_cache), "）。")
    if (identical(dem_src, "cache")) md_log(CFG, "INFO", "[61] 使用 DEM 缓存")
  }

  # ---------------------------------------------------------------------------
  # 5. 绘图（画法逐项对齐 07_S1_study_area_map.R §6–§9）
  # ---------------------------------------------------------------------------
  pt <- md_pt(CFG)
  # 暖砂褐地形配色（原脚本 §4 逐字：低地 #B9D7A8 -> 高地 #C9AE7E，两色线性渐变）
  dem_cols <- CFG$figure$palette$dem %||% MD_PAL_DEFAULT$elevation
  fill_var <- fac[1]
  if (is.null(fill_var)) {
    st$.grp <- "all"
    fill_var <- ".grp"
  }
  pal <- md_pal(CFG, unique(st[[fill_var]]), fill_var)
  # 点大小照抄原脚本：主图 2.5 * 2/3、放大图 3 * 2/3（原脚本 §6 / §8 的 geom_point(size=)）
  pt_main <- 2.5 * 2 / 3
  pt_zoom <- 3 * 2 / 3
  # DEM 色条刻度：原脚本 seq(0, 2000, by = 1000)（本包 dem_limits 可配置，按其跨度三等分）
  dem_breaks <- seq(dem_lim[1], dem_lim[2], by = (dem_lim[2] - dem_lim[1]) / 2)

  base_layers <- function(g, keep_legend = TRUE) {
    if (!is.null(dem_df)) {
      g <- g + ggplot2::geom_tile(
        data = dem_df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$elev),
        inherit.aes = FALSE) +
        ggplot2::scale_fill_gradientn(
          colours = dem_cols, na.value = "transparent", name = "Elevation (m)",
          limits = dem_lim, breaks = dem_breaks,
          guide = if (keep_legend)
            ggplot2::guide_colourbar(barwidth = 10, barheight = 0.8,
                                     title.position = "top") else "none")
    }
    # 海岸线（原脚本 §6：不裁剪，覆盖整个区域；grey30 / 0.2）
    if (!is.null(coast_sf)) {
      g <- g + ggplot2::geom_sf(data = coast_sf, colour = "grey30",
                                linewidth = 0.2, inherit.aes = FALSE)
    }
    g
  }

  mk_map <- function(bb, show_lab, tag, show_y = TRUE, keep_dem_legend = TRUE,
                     grid_brk = NULL, point_size = pt_main, zoom = FALSE) {
    g <- ggplot2::ggplot()
    g <- base_layers(g, keep_legend = keep_dem_legend)
    g <- g + ggplot2::geom_point(data = st,
                                 ggplot2::aes(x = .data[[lon]], y = .data[[lat]],
                                              colour = .data[[fill_var]]),
                                 size = point_size, shape = 16)
    # Zone 图例（原脚本 §6 的 guide_legend：实心圆、放大的图例键）
    g <- g + ggplot2::scale_colour_manual(
      values = pal, name = fill_var,
      guide = if (identical(tag, "A"))
        ggplot2::guide_legend(override.aes = list(shape = 16, size = 3),
                              label.theme = ggplot2::element_text(size = pt$tick),
                              keywidth = 0.5, keyheight = 0.5) else "none")
    if (show_lab && any(nzchar(st[[lab_col]]))) {
      # ggrepel 参数照抄原脚本（主图 §6 force = 10；放大图 §8 force = 8）；
      # seed 用全包唯一种子 CFG$opts$seed（原脚本写死 42，本包禁用字面量种子）
      g <- g + ggrepel::geom_text_repel(
        data = st, ggplot2::aes(x = .data[[lon]], y = .data[[lat]],
                                label = .data[[lab_col]]),
        size = pt2gg(pt$ann),
        box.padding = if (zoom) 0.4 else 0.3,
        point.padding = 0.3,
        max.overlaps = 20, segment.size = 0.2,
        force = if (zoom) 8 else 10,
        min.segment.length = 0.1, seed = CFG$opts$seed)
    }
    if (!is.null(CFG$map$highlight) && length(CFG$map$highlight) == 4) {
      hb <- as.numeric(CFG$map$highlight)
      g <- g + ggplot2::annotate("rect", xmin = hb[1], xmax = hb[3],
                                 ymin = hb[2], ymax = hb[4],
                                 fill = NA, colour = "grey25",
                                 linewidth = 0.3, linetype = "dashed")
    }
    g <- g + ggplot2::coord_sf(xlim = c(bb[1], bb[3]), ylim = c(bb[2], bb[4]),
                               crs = 4326)
    if (!is.null(grid_brk)) {
      # 刻度标签强制纯 ASCII：coord_sf(crs = 4326) 默认会写 "108°E"，
      # 而本包硬规则是图中文字一律 ASCII（见 FIGURE_ALIGN_SPEC §1.6 / pitfalls B12），
      # 与轴题 "Longitude (deg E)" 的写法保持一致。
      lab_xy <- function(suffix) function(v)
        paste0(format(v, trim = TRUE, scientific = FALSE), suffix)
      g <- g + ggplot2::scale_x_continuous(breaks = grid_brk$x,
                                           labels = lab_xy(" E")) +
        ggplot2::scale_y_continuous(breaks = grid_brk$y,
                                    labels = lab_xy(" N"))
    }
    g + ggplot2::labs(x = "Longitude (deg E)",
                      y = if (show_y) "Latitude (deg N)" else NULL, tag = tag) +
      theme_md(CFG = CFG) +
      # 主题照抄原脚本 §6/§8：grey85 网格 / aliceblue 底色 / 图例无框；
      # 放大图不显示图例、x 刻度 45 度
      ggplot2::theme(panel.grid = ggplot2::element_line(colour = "grey85", linewidth = 0.3),
                     panel.background = ggplot2::element_rect(fill = "aliceblue", colour = NA),
                     legend.key = ggplot2::element_blank(),
                     legend.background = ggplot2::element_blank(),
                     legend.box.background = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(5, 10, 5, 5)) +
      (if (zoom)
         ggplot2::theme(legend.position = "none",
                        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                           size = pt$tick),
                        axis.text.y = ggplot2::element_text(size = pt$tick))
       else ggplot2::theme())
  }

  brk_fun <- function(v) {
    r <- diff(range(v))
    step <- signif(r / 4, 1)
    seq(floor(min(v) / step) * step, ceiling(max(v) / step) * step, by = step)
  }
  pA <- mk_map(bbox, isTRUE(CFG$map$point_label) && is.null(CFG$map$countries),
               "A", TRUE, TRUE,
               list(x = brk_fun(bbox[c(1, 3)]), y = brk_fun(bbox[c(2, 4)])),
               point_size = pt_main, zoom = FALSE)
  pB <- mk_map(zoom_bbox, isTRUE(CFG$map$point_label), "B", FALSE, FALSE,
               list(x = brk_fun(zoom_bbox[c(1, 3)]), y = brk_fun(zoom_bbox[c(2, 4)])),
               point_size = pt_zoom, zoom = TRUE)

  panels <- list(A = pA, B = pB)
  # 原脚本 §9：plot_layout(widths = c(1, 0.8))，即 A 占 5/9、B 占 4/9
  design <- "AAAAABBBB"

  # ---- 可选世界 inset（原脚本 §7）----
  if (!is.null(ctry_arg) && length(ctry_arg)) {
    wld <- tryCatch(rnaturalearth::ne_countries(scale = "small", returnclass = "sf"),
                    error = function(e) NULL)
    if (!is.null(wld)) {
      pC <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = wld, fill = "grey88", colour = NA) +
        ggplot2::annotate("rect", xmin = bbox[1], xmax = bbox[3],
                          ymin = bbox[2], ymax = bbox[4],
                          fill = NA, colour = "red", linewidth = 0.7,
                          linetype = "dashed") +
        ggplot2::coord_sf(xlim = c(-180, 180), ylim = c(-60, 85)) +
        ggplot2::labs(tag = "C") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white",
                                                               colour = "grey70",
                                                               linewidth = 0.4),
                       plot.margin = ggplot2::margin(0, 0, 0, 0),
                       plot.tag = ggplot2::element_text(size = pt$tag, face = "bold"))
      panels$C <- pC
      design <- "AB\nCC"
      add_note("世界 inset 已绘出（灰色陆地 + 红色虚线研究区框，画法同原脚本 07 §7）；",
               "map$countries = ", paste(ctry_arg, collapse = ", "), "。")
    } else {
      add_note("世界 inset 数据取不到（需联网），已省略 inset 面板。")
    }
  } else {
    add_note("map$countries = NULL，按配置不画世界 inset。")
  }

  fig <- md_compose(CFG, panels, design = design, guides = "collect", legend = "bottom")
  # 原脚本 FIGURE_PANELS$FigS1 = 183 x 120 mm
  md_save_figure(CFG, fig, "61_map", height = if (length(panels) > 2) 130 else 120)

  # ---------------------------------------------------------------------------
  # 6. 出表
  # ---------------------------------------------------------------------------
  note_tab <- data.frame(
    Item = c("map$enabled", "Main bbox", "Zoom bbox", "DEM available", "DEM source",
             "DEM cache", "DEM z / aggregate", "DEM limits", "Basemap",
             "Point labels", "Countries (world inset)", "offline"),
    Value = c("TRUE",
              paste(bbox, collapse = ", "),
              paste(zoom_bbox, collapse = ", "),
              !is.null(dem_df),
              if (is.null(dem_df)) paste0("none (", dem_reason %||% "unknown", ")")
              else dem_src,
              basename(dem_cache),
              paste0(dem_zoom, " / ", dem_fact, "x", dem_fact),
              paste(dem_lim, collapse = "-"),
              if (is.null(coast_sf)) "none (points only)" else "coastline",
              isTRUE(CFG$map$point_label),
              if (is.null(ctry_arg)) "none" else paste(ctry_arg, collapse = ", "),
              isTRUE(CFG$map$offline)),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "61_map", list(
    Sites = st,
    Settings = note_tab,
    Notes = data.frame(Note = if (length(notes)) notes else "none",
                       stringsAsFactors = FALSE)
  ), meta_extra = c(dem_available = !is.null(dem_df),
                    dem_source = if (is.null(dem_df)) "none" else dem_src))

  md_result(module, "OK",
            tables = basename(tp), figures = "figures/61_map.pdf",
            notes = notes)
}
