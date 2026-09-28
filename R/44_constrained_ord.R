# =============================================================================
# R/44_constrained_ord.R ---- 模块 44：约束排序（db-RDA / RDA / CCA / capscale）
# -----------------------------------------------------------------------------
# 源：毕业论文 3_environmental_drivers.R:855-1000（Part A: db-RDA triplot，Fig3b-A）
#     契约见 _build-notes/micro-div-pipeline/ROUND5_MODULE_SPEC.md
#
# 导出：md_constrained_ord(CFG)
#
# 计算：全部走共享 md_cord_fit() ——
#       VIF 逐步剔除（阈值 15）→ anova(permutations = 999) → anova(by = "terms")
#       → RsquareAdj → scores(display = "sites" / "bp") → 特征值解释率。
#  ⚠ md_cord_fit(CFG, Y, X, ...) 的 **Y 是"样本 × 变量"**（它用 rownames(Y) 与样本
#    交集取 ids），所以村落矩阵在这里先 t() 成 样本 × ASV 再传进去。
#     任务书里"ASV×样本 直接传"的说法与共享函数签名不一致，以共享函数为准。
#
#  响应侧（CFG$constrained_ord$response）：
#    community —— md_core_asv(CFG, "alpha") 的 t()；response_vars 可只取一部分行
#    distance  —— md_dist_get(CFG, response_distance)，取不到就退回 community
#  解释侧：md_core_analysis(CFG) 里取 constraints（默认 env$vars），
#          模块只做"字符列转数值"的预处理，其余（complete.cases / 数值筛选 /
#          VIF 逐步剔除）由 md_cord_fit() 内部完成。
#
# 图：figures/44_constrained_ord.pdf（单面板 triplot，183 × 110 mm）
#     样方点按 design$factors[1] 着色、factors[2] 定形状；env 箭头 +
#     ggrepel 排斥标签（白色描边）；箭头长度按"最长箭头 = 站点云半径 × arrow_ratio"
#     归一；副标题 "R2adj = ..., P ..."（md_fmt_p_fig()，纯 ASCII）。
# 表：results/44_constrained_ord.xlsx（Model_summary / Terms / VIF / Site_scores /
#     Env_scores / Note）
# =============================================================================

# ---------------------------------------------------------------------------
# local_* 辅助（只在本模块内使用）
# ---------------------------------------------------------------------------

#' 图上文字一律纯 ASCII（95_export_tables 会对 PDF 里的可见字符做审计）
local_cord_ascii <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- suppressWarnings(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"))
  x[is.na(x)] <- ""
  gsub("[^ -~]", "_", x)
}

#' 轴名前缀：照抄原图 dbRDA1 / dbRDA2（rda 时写 RDA1/RDA2）
local_cord_prefix <- function(method) {
  switch(tolower(as.character(method)),
         dbrda = "dbRDA", capscale = "CAP", rda = "RDA", cca = "CCA", "dbRDA")
}

#' 面板标题里的方法名（原图写 "db-RDA"）
local_cord_label <- function(method) {
  switch(tolower(as.character(method)),
         dbrda = "db-RDA", capscale = "capscale", rda = "RDA", cca = "CCA", "db-RDA")
}

#' 形状值：与 14_ordination 同源（02_Fig02.R §B/C 的 16/17/15/18）
local_cord_shape <- function(lv) {
  builtin <- c(Bare = 16, Herb = 17, Shrub = 15, Tree = 18)
  out <- unname(builtin[as.character(lv)])
  if (anyNA(out)) {
    pool <- setdiff(c(16, 17, 15, 18, 8, 3, 7, 4, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14), out)
    out[is.na(out)] <- pool[seq_len(sum(is.na(out)))]
  }
  out
}

# ---------------------------------------------------------------------------
# 模块主体
# ---------------------------------------------------------------------------
md_constrained_ord <- function(CFG) {
  module <- "44"
  if (!isTRUE(CFG$constrained_ord$enabled %||% TRUE))
    return(md_result(module, "SKIP", "constrained_ord$enabled = FALSE"))
  if (!md_need(CFG, "ecology", module))
    return(md_result(module, "SKIP", "缺少依赖包（vegan）"))

  ana <- tryCatch(md_core_analysis(CFG), error = function(e) NULL)
  smp <- CFG$data$sample_col
  if (is.null(ana) || nrow(ana) < 5 || !smp %in% names(ana)) {
    msg <- "分析总表不可用（先跑模块 05；或样本数 < 5）"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }

  notes <- character(0)

  # ---- 分析段：全部计算都在这里，命中缓存时整段跳过 ----
  md_stage(CFG, "44_analysis", {
    md_seed(CFG, md_offset(CFG, "constrained_ord_44", 8400L))
    co_notes <- character(0)

    co_method <- tolower(as.character(CFG$constrained_ord$method %||% "dbrda"))
    if (!co_method %in% c("dbrda", "rda", "cca", "capscale")) {
      co_notes <- c(co_notes, paste0("method '", co_method, "' 不可识别 —— 已退回 dbrda"))
      co_method <- "dbrda"
    }
    co_response  <- tolower(as.character(CFG$constrained_ord$response %||% "community"))
    co_transform <- as.character(CFG$constrained_ord$transform %||% "hellinger")
    co_vif_thr   <- CFG$constrained_ord$vif_thr %||% 15
    co_perm      <- as.integer(CFG$constrained_ord$permutations %||% 999L)
    co_arrow     <- CFG$constrained_ord$arrow_ratio %||% 0.78
    co_metric    <- as.character(CFG$constrained_ord$response_distance %||% "bray")

    # ---- 解释变量 ----
    co_constr <- as.character(CFG$constrained_ord$constraints %||%
                                CFG$env$vars %||% character(0))
    co_constr <- intersect(co_constr, names(ana))
    miss_c <- setdiff(as.character(CFG$constrained_ord$constraints %||%
                                     CFG$env$vars %||% character(0)), names(ana))
    if (length(miss_c))
      co_notes <- c(co_notes, paste0("CFG 里点名但数据中缺失的解释变量：",
                                     paste(miss_c, collapse = ", ")))
    if (!length(co_constr))
      co_notes <- c(co_notes, "没有任何解释变量可用（constraints / env$vars）")
    # 字符列先尝试转数值（CSV 读进来常常是 character），转不动就交给 md_cord_fit() 丢掉
    for (v in co_constr) {
      if (!is.numeric(ana[[v]])) {
        z <- suppressWarnings(as.numeric(ana[[v]]))
        if (any(is.finite(z))) ana[[v]] <- z
      }
    }
    co_constr <- co_constr[vapply(ana[co_constr], is.numeric, logical(1))]

    # ---- 响应侧 ----
    co_Y <- NULL
    if (identical(co_response, "distance")) {
      co_Y <- tryCatch(md_dist_get(CFG, co_metric), error = function(e) NULL)
      if (is.null(co_Y))
        co_notes <- c(co_notes, sprintf(
          "response = 'distance' 但取不到距离 '%s'（缺 13 号缓存）—— 已退回 community", co_metric))
    }
    if (is.null(co_Y)) {
      co_response <- "community"
      co_M <- tryCatch(md_core_asv(CFG, "alpha"), error = function(e) NULL)
      if (!is.null(co_M) && nrow(co_M) >= 2 && ncol(co_M) >= 5) {
        rv <- CFG$constrained_ord$response_vars
        if (!is.null(rv)) {
          rv <- intersect(as.character(rv), rownames(co_M))
          if (length(rv)) co_M <- co_M[rv, , drop = FALSE]
          else co_notes <- c(co_notes, "response_vars 与矩阵行名无交集 —— 已改用全部行")
        }
        co_M <- co_M[rowSums(co_M) > 0, , drop = FALSE]
        if (nrow(co_M) >= 2) co_Y <- t(co_M)   # md_cord_fit() 要 样本 × 变量
      }
      if (is.null(co_Y))
        co_notes <- c(co_notes, "响应矩阵不可用（05 核心对象缺失或小于 5 样本）")
    }

    # ---- 拟合 ----
    co_X <- NULL
    if (!is.null(co_Y) && length(co_constr) >= 1) {
      co_X <- as.data.frame(ana[, co_constr, drop = FALSE], check.names = FALSE)
      rownames(co_X) <- as.character(ana[[smp]])
    }
    co_fit <- NULL
    if (!is.null(co_X)) {
      co_fit <- tryCatch(
        md_cord_fit(CFG, co_Y, co_X, method = co_method, transform = co_transform,
                    vif_thr = co_vif_thr, permutations = co_perm),
        error = function(e) {
          co_notes <<- c(co_notes, paste0("md_cord_fit() 报错：", conditionMessage(e)))
          NULL
        })
      if (is.null(co_fit))
        co_notes <- c(co_notes, "md_cord_fit() 返回 NULL（样本不足 / 解释变量全被剔除 / 拟合失败）")
    }

    # ---- 把得分、箭头、轴题解出来（供段外画图与写表用）----
    co_sites <- NULL
    co_envs  <- NULL
    co_arrow_mul <- NA_real_
    co_meta  <- NULL
    if (!is.null(co_fit)) {
      sc  <- co_fit$sites
      sid <- as.character(sc$sample_id)
      if (!"sample_id" %in% names(sc)) { sc$sample_id <- rownames(sc); sid <- rownames(sc) }
      ax <- setdiff(names(sc), "sample_id")
      if (length(ax) >= 2) {
        sc$.a1 <- as.numeric(sc[[ax[1]]])
        sc$.a2 <- as.numeric(sc[[ax[2]]])
        facs <- intersect(CFG$design$factors %||% character(0), names(ana))
        if (length(facs)) {
          sc <- merge(sc, ana[, c(smp, facs), drop = FALSE],
                      by.x = "sample_id", by.y = smp, all.x = TRUE, sort = FALSE)
          sc <- sc[match(sid, as.character(sc$sample_id)), , drop = FALSE]
        }
        co_sites <- sc
      } else {
        co_notes <- c(co_notes, "scores(display = 'sites') 少于 2 根轴 —— 无法出 triplot")
      }

      es <- co_fit$env_scores
      if (!is.null(es) && nrow(es) > 0 && length(ax) >= 2) {
        if (!"env_var" %in% names(es)) es$env_var <- rownames(es)
        eax <- setdiff(names(es), "env_var")
        if (length(eax) >= 2) {
          es$.e1 <- as.numeric(es[[eax[1]]])
          es$.e2 <- as.numeric(es[[eax[2]]])
          # 箭头长度归一：最长箭头 = 站点云半径 × arrow_ratio（原脚本 0.78）
          sr <- max(sqrt(co_sites$.a1^2 + co_sites$.a2^2), na.rm = TRUE)
          br <- max(sqrt(es$.e1^2 + es$.e2^2), na.rm = TRUE)
          co_arrow_mul <- if (is.finite(sr) && is.finite(br) && sr > 0 && br > 0)
            co_arrow * sr / br else 0
          es$.x <- es$.e1 * co_arrow_mul
          es$.y <- es$.e2 * co_arrow_mul
          co_envs <- es
        }
      }
      if (is.null(co_envs))
        co_notes <- c(co_notes, "scores(display = 'bp') 不可用 —— 只出样方点，不出环境箭头")

      # 轴题（dbRDA1 (xx.x%) / dbRDA2 (xx.x%)）
      ve <- as.numeric(co_fit$var_exp)
      ve[!is.finite(ve)] <- NA_real_
      co_pfx  <- local_cord_prefix(co_method)
      co_xlab <- if (is.na(ve[1])) paste0(co_pfx, "1") else sprintf("%s1 (%.1f%%)", co_pfx, ve[1])
      co_ylab <- if (is.na(ve[2])) paste0(co_pfx, "2") else sprintf("%s2 (%.1f%%)", co_pfx, ve[2])

      # 模型整体检验（anova(permutations) 里第一行有 Pr(>F) 的）
      co_p <- NA_real_; co_F <- NA_real_
      an <- co_fit$anova
      if (!is.null(an) && "Pr(>F)" %in% names(an)) {
        ii <- which(!is.na(an[["Pr(>F)"]]))
        if (length(ii)) {
          co_p <- as.numeric(an[["Pr(>F)"]][ii[1]])
          if ("F" %in% names(an)) co_F <- as.numeric(an[["F"]][ii[1]])
        }
      }
      co_sub <- sprintf("R2adj = %s, %s",
                        if (is.finite(co_fit$R2adj)) formatC(co_fit$R2adj, format = "f", digits = 3) else "NA",
                        md_fmt_p_fig(co_p))
      resp_lab <- if (identical(co_fit$response, "distance")) {
        if (tolower(co_metric) == "bray") "Bray-Curtis" else co_metric
      } else if (identical(co_transform, "hellinger")) "Hellinger" else "raw community"
      co_title <- sprintf("%s (%s)", local_cord_label(co_method), resp_lab)

      co_meta <- list(method = co_method, response = co_fit$response,
                      transform = co_transform, distance = co_metric,
                      n = co_fit$n, R2 = co_fit$R2, R2adj = co_fit$R2adj,
                      F = co_F, p = co_p, permutations = co_perm,
                      var_exp = ve, variables = co_fit$variables,
                      arrow_ratio = co_arrow, arrow_mul = co_arrow_mul)
    }
  })
  if (length(co_notes)) notes <- c(notes, co_notes)

  if (is.null(co_fit) || is.null(co_sites) || !all(c(".a1", ".a2") %in% names(co_sites))) {
    msg <- if (length(co_notes)) paste(co_notes, collapse = " | ") else "约束排序不可用"
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", "约束排序不可用（见 notes）", notes = co_notes))
  }

  # ---------------------------------------------------------------------------
  # 出图：单面板 triplot（原文 Fig3b-A 的美化逐项复刻）
  #   ① 白底无网格（theme_md，不用 ggplot 默认灰底）
  #   ② ggrepel 排斥标签 + 白色描边
  #   ③ 箭头按站点云半径归一（最长箭头 = 半径 × arrow_ratio）
  #   ④ 过原点浅灰参考线
  #   ⑤ R2adj / P 写成副标题（原文由面板内注释上移到副标题）
  # ---------------------------------------------------------------------------
  pt    <- md_pt(CFG)
  style <- CFG$constrained_ord$style %||% list()
  facs  <- intersect(CFG$design$factors %||% character(0), names(co_sites))
  colour_var <- if (length(facs) >= 1) facs[1] else NULL
  shape_var  <- if (length(facs) >= 2) facs[2] else NULL
  if (!is.null(shape_var) &&
      length(unique(stats::na.omit(as.character(co_sites[[shape_var]])))) < 2)
    shape_var <- NULL
  if (!is.null(colour_var) &&
      length(unique(stats::na.omit(as.character(co_sites[[colour_var]])))) < 2)
    colour_var <- NULL

  local_cord_panel <- function() {
    d <- co_sites
    d <- d[!is.na(d$.a1) & !is.na(d$.a2), , drop = FALSE]
    g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$.a1, y = .data$.a2)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey90") +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey90")

    if (!is.null(colour_var)) {
      d[[colour_var]] <- factor(as.character(d[[colour_var]]))
      pal <- md_pal(CFG, levels(d[[colour_var]]), colour_var)
      if (!is.null(shape_var)) {
        d[[shape_var]] <- factor(as.character(d[[shape_var]]))
        g <- g + ggplot2::geom_point(
          data = d, ggplot2::aes(colour = .data[[colour_var]], shape = .data[[shape_var]]),
          size = style$point_size %||% MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha,
          stroke = 0)
        g <- g + ggplot2::scale_shape_manual(
          values = stats::setNames(local_cord_shape(levels(d[[shape_var]])), levels(d[[shape_var]])),
          name = shape_var)
      } else {
        g <- g + ggplot2::geom_point(
          data = d, ggplot2::aes(colour = .data[[colour_var]]),
          size = style$point_size %||% MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha,
          stroke = 0)
      }
      g <- g + ggplot2::scale_colour_manual(values = pal, name = colour_var, drop = FALSE) +
        ggplot2::guides(colour = md_guide_flex(length(pal), "bottom", 4))
    } else {
      g <- g + ggplot2::geom_point(
        data = d, size = style$point_size %||% MD_PT_SIZE$point,
        alpha = MD_PT_SIZE$point_alpha, stroke = 0, colour = "grey25")
    }

    if (!is.null(co_envs) && nrow(co_envs) > 0 && is.finite(co_arrow_mul) && co_arrow_mul > 0) {
      g <- g + ggplot2::geom_segment(
        data = co_envs,
        ggplot2::aes(x = 0, y = 0, xend = .data$.x, yend = .data$.y),
        inherit.aes = FALSE,
        arrow = ggplot2::arrow(length = grid::unit(0.09, "cm"), type = "closed"),
        colour = style$arrow_col %||% "#3182BD", linewidth = 0.45)
      lab_df <- data.frame(x = co_envs$.x, y = co_envs$.y,
                           lab = local_cord_ascii(co_envs$env_var),
                           stringsAsFactors = FALSE)
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        g <- g + ggrepel::geom_text_repel(
          data = lab_df, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
          inherit.aes = FALSE,
          size = pt2gg(style$label_size_pt %||% pt$ann), colour = "black",
          box.padding = 0.26, point.padding = 0.16, min.segment.length = 0,
          segment.size = 0.2, segment.colour = "grey55",
          max.overlaps = Inf, seed = CFG$opts$seed, bg.colour = "white", bg.r = 0.12)
      } else {
        g <- g + ggplot2::geom_text(
          data = lab_df, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
          inherit.aes = FALSE, size = pt2gg(style$label_size_pt %||% pt$ann),
          colour = "black", vjust = -0.25)
      }
    }

    g + ggplot2::coord_fixed() +
      ggplot2::labs(x = co_xlab, y = co_ylab, title = co_title, subtitle = co_sub) +
      theme_md(CFG = CFG) +
      ggplot2::theme(legend.position = "bottom",
                     legend.box = "vertical",
                     legend.box.just = "center",
                     legend.margin = ggplot2::margin(0, 0, 0, 0))
  }

  pA <- tryCatch(local_cord_panel(), error = function(e) {
    md_warn(CFG, "44: triplot 绘制失败 —— ", conditionMessage(e))
    md_empty_plot(CFG, "Triplot unavailable")
  })
  md_save_figure(CFG, pA, "44_constrained_ord", height = 110)
  md_maybe_single(CFG, pA, "44_constrained_ord_A")

  # ---------------------------------------------------------------------------
  # 结果表
  # ---------------------------------------------------------------------------
  co_pfx2 <- local_cord_prefix(co_meta$method)
  term_tab <- NULL
  if (!is.null(co_fit$terms)) {
    t <- as.data.frame(co_fit$terms, check.names = FALSE)
    t$Term <- rownames(t)
    rownames(t) <- NULL
    t <- t[!t$Term %in% c("Residual", "Total"), , drop = FALSE]
    if (nrow(t)) {
      t$Var_pct <- if ("SumOfSqs" %in% names(t)) t$SumOfSqs / sum(t$SumOfSqs, na.rm = TRUE) * 100 else NA_real_
      pcol <- intersect(c("Pr(>F)", "Pr..F."), names(t))[1]
      t$sig <- if (!is.na(pcol)) md_sig_label(t[[pcol]]) else NA_character_
      t <- t[, c("Term", setdiff(names(t), "Term")), drop = FALSE]
      term_tab <- t
    }
  }
  if (!is.null(term_tab))
    notes <- c(notes, sprintf("Terms sheet 的 sig 列用 md_sig_label()（* ** *** ns），共 %d 个边际项。",
                              nrow(term_tab)))
  else
    notes <- c(notes, "anova(by = 'terms') 不可用 —— Terms sheet 为空。")

  vif_tab <- NULL
  if (!is.null(co_fit$vif) && nrow(co_fit$vif)) {
    vif_tab <- co_fit$vif
    vif_tab$threshold <- co_vif_thr
    vif_tab$retained <- vif_tab$variable %in% co_meta$variables
    vif_tab$flag <- ifelse(is.na(vif_tab$vif), "NA",
                           ifelse(vif_tab$vif > co_vif_thr, "collinear", "ok"))
  }
  notes <- c(notes, sprintf(
    "VIF sheet 是**逐步剔除后**的最终 VIF（原文 886-888 行还打印了初始 VIF；md_cord_fit() 只返回最终值）。阈值 %s，保留 %d / %d 个解释变量。",
    as.character(co_vif_thr), length(co_meta$variables), length(co_constr)))

  site_tab <- data.frame(sample_id = as.character(co_sites$sample_id),
                         Axis1 = co_sites$.a1, Axis2 = co_sites$.a2,
                         stringsAsFactors = FALSE)
  names(site_tab)[2:3] <- c(paste0(co_pfx2, "1"), paste0(co_pfx2, "2"))
  for (f in facs) if (f %in% names(co_sites)) site_tab[[f]] <- as.character(co_sites[[f]])

  env_tab <- NULL
  if (!is.null(co_envs) && nrow(co_envs)) {
    env_tab <- data.frame(
      Variable      = local_cord_ascii(co_envs$env_var),
      Axis1         = co_envs$.e1,
      Axis2         = co_envs$.e2,
      Axis1_scaled  = co_envs$.x,
      Axis2_scaled  = co_envs$.y,
      stringsAsFactors = FALSE)
    names(env_tab)[2:5] <- c(paste0(co_pfx2, "1"), paste0(co_pfx2, "2"),
                             paste0(co_pfx2, "1_scaled"), paste0(co_pfx2, "2_scaled"))
  }

  ms <- data.frame(
    item = c("module", "method", "response", "transform", "distance", "n_samples",
             "n_constraints_retained", "R2", "R2_adj", "F", "p", "permutations",
             "axis1_explained_pct", "axis2_explained_pct", "arrow_ratio",
             "arrow_multiplier", "variables_retained"),
    value = c("44", co_meta$method, co_meta$response, co_meta$transform, co_meta$distance,
              as.character(co_meta$n), as.character(length(co_meta$variables)),
              as.character(co_meta$R2), as.character(co_meta$R2adj),
              as.character(co_meta$F), as.character(co_meta$p),
              as.character(co_meta$permutations),
              as.character(co_meta$var_exp[1]), as.character(co_meta$var_exp[2]),
              as.character(co_meta$arrow_ratio), as.character(co_meta$arrow_mul),
              paste(co_meta$variables, collapse = ", ")),
    stringsAsFactors = FALSE)

  note_df <- data.frame(Note = c(
    sprintf("口径源：毕业论文 3_environmental_drivers.R:855-1000（db-RDA triplot，Fig3b-A）。"),
    sprintf("VIF 逐步剔除阈值 = %s（原脚本 while (any(vif > 15))）；置换次数 = %d。",
            as.character(co_vif_thr), as.integer(co_perm)),
    sprintf("解释率 = CCA$eig[1:2] / tot.chi * 100（原脚本 928 行口径）。"),
    sprintf("箭头长度：最长箭头 = 站点云半径 × %s（原脚本 0.78；这里由 CFG$constrained_ord$arrow_ratio 控制）。",
            as.character(co_meta$arrow_ratio)),
    sprintf("样方点着色 = %s；形状 = %s（原文 Fig3b-A 只按 Vegetation 着色，本包按 design$factors 顺序取）。",
            colour_var %||% "none", shape_var %||% "none"),
    notes))

  tp <- md_write_xlsx(CFG, "44_constrained_ord", list(
    Model_summary = ms,
    Terms         = term_tab,
    VIF           = vif_tab,
    Site_scores   = site_tab,
    Env_scores    = env_tab,
    Note          = note_df),
    meta_extra = list(n = co_meta$n, method = co_meta$method,
                      response = co_meta$response, R2adj = co_meta$R2adj))

  md_result(module, "OK", tables = basename(tp),
            figures = "figures/44_constrained_ord.pdf",
            notes = paste(notes, collapse = " | "))
}
