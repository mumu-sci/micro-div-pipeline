# =============================================================================
# R/81_taxon_group.R ---- 模块 81：通用类群分组（组成 / 网络角色 / 功能关联）
# -----------------------------------------------------------------------------
# 源：毕业论文 6_ecological_guilds.R（1492 行）
#   参考段落：84-130（映射表构造）、1239-1333（组成与网络角色绘图段）、
#             1460-1472（出图规格：heights = c(1, 1.2)）
#
# ★ 与原文的关键区别（用户明确要求）：
#   **删除硬编码的固氮属 -> 共生/联合/自生 映射表**，改为读通用映射表
#   `data/05_taxon_group.csv`（两列 taxon, group；列名可由
#   CFG$taxon_group$taxon_col / $group_col 指定）。
#   分组含义不限（功能群 / 生态位 / R-K 对策 / 任意自定义）；
#   按哪个分类层级匹配由 CFG$taxon_group$match_by 决定：
#     "ASV" -> 用 asv_id 直接匹配；其它（Phylum..Species）-> 用 md_core(CFG)$tax 对应列匹配。
#   未匹配到的 ASV：drop_unmapped = FALSE 归入 "Unmapped"，TRUE 则剔除。
#   文件不存在 / file = NULL / enabled = FALSE -> md_skip() + SKIP，不报错。
#
# 做什么
#   1 组成：每个样本内各 group 的相对丰度（%）
#   2 梯度检验：group ~ 各设计因子 的 KW + Dunn + FDR（直接用 md_group_test）
#   3 网络角色（软依赖）：读 cache/80_network.rds；不存在则跳过 A 面板并在 notes 说明
#   4 响应变量关联：group 相对丰度 vs CFG$response$var 的 Spearman
#   5 跨组互作强度：80 边表上 组内边 vs 组间边的 Fisher 精确检验
#   6 类群指示种（CFG$taxon_group$indicator$enabled = TRUE）：md_indval()
#
# 图 81_taxon_group.pdf（183 x 150 mm，heights = style$panel_heights = c(1, 1.2)）
#   A group x 网络角色 堆叠柱（角色配色 = MD_NET$role_colors）；缺 80 产物则整图退化为只有 B
#   B group 组成堆叠柱（横轴 = 样本分组，fill = group）
#
# 依赖共享库：md_core / md_core_asv / md_rel / md_glom / md_group_test / md_indval /
#             md_pal / md_compose / md_save_figure / md_write_xlsx / md_path /
#             MD_NET / md_empty_plot / pt2gg / md_result
# =============================================================================

#' 模块 81 入口
#' @param CFG config 列表（见 config/config_template.R 的 taxon_group 块）
md_taxon_group <- function(CFG) {
  module <- "81"

  if (!isTRUE(CFG$taxon_group$enabled %||% TRUE))
    return(md_result(module, "SKIP", "taxon_group$enabled = FALSE"))
  if (!md_need(CFG, c("phyloseq"), module))
    return(md_result(module, "SKIP", "缺少依赖包（phyloseq）"))

  tg      <- CFG$taxon_group
  f_cfg   <- tg$file %||% CFG$data$taxon_group
  if (is.null(f_cfg) || !nzchar(as.character(f_cfg))) {
    md_skip(CFG, module, "taxon_group$file = NULL（未提供类群 -> 分组映射表）")
    return(md_result(module, "SKIP", "未提供映射表"))
  }
  map_path <- if (grepl("^(/|[A-Za-z]:)", f_cfg)) f_cfg else file.path(CFG$project$dir, f_cfg)
  if (!file.exists(map_path)) {
    md_skip(CFG, module, paste0("类群映射表不存在: ", map_path))
    return(md_result(module, "SKIP", "映射表不存在"))
  }

  core <- tryCatch(md_core(CFG), error = function(e) NULL)
  if (is.null(core) || is.null(core$asv_raw)) {
    md_skip(CFG, module, "取不到 cache/05_core.rds")
    return(md_result(module, "SKIP", "缺 05 产物"))
  }

  tcol     <- tg$taxon_col %||% "taxon"
  gcol     <- tg$group_col %||% "group"
  match_by <- tg$match_by  %||% "Genus"
  label    <- tg$group_label %||% "Group"
  drop_un  <- isTRUE(tg$drop_unmapped %||% FALSE)
  grp_ord  <- tg$group_order
  st       <- tg$style %||% list()
  pt       <- md_pt(CFG, comp = FALSE)
  digits   <- CFG$table$digits %||% 3
  id_col   <- CFG$data$id_col %||% "asv_id"
  smp      <- CFG$data$sample_col %||% "sample_id"

  # ---- 读映射表（列为空 / 列名不对 -> 跳过，不报错）----
  map <- tryCatch(utils::read.csv(map_path, stringsAsFactors = FALSE,
                                  fileEncoding = "UTF-8-BOM"),
                  error = function(e) NULL)
  if (is.null(map)) map <- tryCatch(utils::read.csv(map_path, stringsAsFactors = FALSE),
                                    error = function(e) NULL)
  if (is.null(map) || !all(c(tcol, gcol) %in% names(map))) {
    msg <- paste0("映射表缺少列 '", tcol, "' / '", gcol, "'（实际列: ",
                  paste(names(map), collapse = ", "), "）")
    md_skip(CFG, module, msg)
    return(md_result(module, "SKIP", msg))
  }
  map <- map[!is.na(map[[tcol]]) & !is.na(map[[gcol]]), c(tcol, gcol), drop = FALSE]
  map[[tcol]] <- trimws(as.character(map[[tcol]]))
  map[[gcol]] <- trimws(as.character(map[[gcol]]))
  map <- map[nzchar(map[[tcol]]) & nzchar(map[[gcol]]), , drop = FALSE]
  dup <- duplicated(map[[tcol]])
  if (any(dup))
    md_warn(CFG, "映射表里 ", sum(dup), " 个 taxon 有重复条目，只保留第一条")
  map <- map[!dup, , drop = FALSE]
  if (nrow(map) == 0) {
    md_skip(CFG, module, "映射表没有有效行")
    return(md_result(module, "SKIP", "映射表为空"))
  }
  map_lk  <- stats::setNames(map[[gcol]], map[[tcol]])
  map_lv  <- unique(map[[gcol]])
  if (!is.null(grp_ord)) {
    grp_lv <- c(intersect(as.character(grp_ord), map_lv),
                setdiff(map_lv, as.character(grp_ord)))
  } else {
    grp_lv <- map_lv
  }
  unmapped_lab <- tg$unmapped_label %||% "Unmapped"

  # ---- 组成口径矩阵（默认走 composition 口径 = 未抽平，与 10_composition 一致）----
  m <- md_core_asv(CFG, tg$normalize_key %||% "composition")
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
    md_skip(CFG, module, "类群组成矩阵为空")
    return(md_result(module, "SKIP", "矩阵为空"))
  }
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  meta <- core$meta
  if (!is.null(meta)) meta <- meta[match(colnames(m), meta[[smp]]), , drop = FALSE]

  fac <- intersect(CFG$design$factors %||% character(0),
                   if (is.null(meta)) character(0) else names(meta))
  fac1 <- fac[1] %||% NULL
  grp_levels <- NULL
  if (!is.null(fac1)) {
    ol  <- CFG$design$ordered_levels[[fac1]]
    obs <- unique(as.character(meta[[fac1]]))
    obs <- obs[!is.na(obs)]
    grp_levels <- if (is.null(ol)) obs else intersect(as.character(ol), obs)
    if (length(grp_levels) == 0) grp_levels <- obs
  }

  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  # ===========================================================================
  # 分析段 —— 整体包进 md_stage()，图与表在段外
  # ===========================================================================
  md_stage(CFG, "81_analysis", {

    # ---- 1. 把每个 ASV 映射到 group ----
    # ⚠ core$tax 的**行数可能多于**过滤后的 ASV（原始注释表里还有没通过丰度/出现率
    #   过滤的 ASV）。所以必须先 `ti <- match(core$tax[[id_col]], asv_ids)` 再
    #   `core$tax[[rank]][ti]` —— 两步缺一不可：
    #     · 只写 `core$tax[[rank]][match(...)]` 会拿**整数下标**去索引数据框的列，
    #       长度仍是 327（不是 321），构造 data.frame 时直接报行数不一致；
    #     · 不做 match 会带上多余的 ASV。
    asv_ids <- rownames(m)
    n_asv <- length(asv_ids)
    tax_key <- rep(NA_character_, n_asv)
    if (!is.null(core$tax) && match_by %in% names(core$tax)) {
      ti <- match(as.character(core$tax[[id_col]]), asv_ids)
      hit <- which(!is.na(ti))
      tax_key[ti[hit]] <- as.character(core$tax[[match_by]][hit])
    }
    if (identical(toupper(match_by), "ASV")) {
      group_of <- unname(map_lk[asv_ids])
      match_level <- "ASV"
      tax_key <- asv_ids
    } else if (is.null(core$tax) || !match_by %in% names(core$tax)) {
      group_of <- rep(NA_character_, n_asv)
      match_level <- "none"
      tax_key <- rep(NA_character_, n_asv)
    } else {
      group_of <- unname(map_lk[tax_key])
      match_level <- match_by
    }
    n_mapped <- sum(!is.na(group_of))

    mapping_tab <- data.frame(
      ASV_ID = asv_ids,
      Match_by = match_level,
      Taxon = tax_key,
      Group = ifelse(is.na(group_of), unmapped_lab, group_of),
      Mapped = !is.na(group_of),
      stringsAsFactors = FALSE)
    if (!is.null(core$tax)) {
      ti <- match(as.character(core$tax[[id_col]]), asv_ids)
      hit <- which(!is.na(ti))
      for (rk in intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                           setdiff(names(core$tax), id_col))) {
        v <- rep(NA_character_, n_asv)
        v[ti[hit]] <- as.character(core$tax[[rk]][hit])
        mapping_tab[[rk]] <- v
      }
    }

    keep_asv <- if (drop_un) !is.na(group_of) else rep(TRUE, length(asv_ids))
    m_keep   <- m[keep_asv, , drop = FALSE]
    g_keep   <- ifelse(is.na(group_of[keep_asv]), unmapped_lab, group_of[keep_asv])
    names(g_keep) <- rownames(m_keep)
    grp_all  <- c(grp_lv, if (!drop_un) unmapped_lab else NULL)

    # ---- 2. 组成：每个样本内各 group 的相对丰度（%）----
    glom <- md_glom(m_keep, data.frame(asv_id = rownames(m_keep), .g = g_keep,
                                       stringsAsFactors = FALSE),
                    rank = ".g", id_col = "asv_id")
    rel_g <- md_rel(glom$mat) * 100
    rel_g[!is.finite(rel_g)] <- 0
    rel_g <- rel_g[intersect(grp_all, rownames(rel_g)), , drop = FALSE]

    comp_long <- data.frame(
      Sample  = rep(colnames(rel_g), each = nrow(rel_g)),
      Group   = rep(rownames(rel_g), times = ncol(rel_g)),
      RelAbund = as.numeric(rel_g),
      stringsAsFactors = FALSE)
    comp_long$SampleGroup <- if (!is.null(fac1))
      as.character(meta[[fac1]][match(comp_long$Sample, meta[[smp]])]) else "All samples"
    for (fv in setdiff(fac, fac1))
      comp_long[[fv]] <- as.character(meta[[fv]][match(comp_long$Sample, meta[[smp]])])

    # 分组 x group 的组平均（堆叠柱数据）
    comp_mean <- if (nrow(comp_long)) {
      agg <- stats::aggregate(RelAbund ~ SampleGroup + Group, data = comp_long, FUN = mean)
      agg$n_sample <- as.integer(table(comp_long$SampleGroup)[agg$SampleGroup])
      agg$Group <- factor(agg$Group, levels = intersect(grp_all, unique(agg$Group)))
      agg$SampleGroup <- factor(agg$SampleGroup,
                                levels = if (length(grp_levels)) grp_levels else unique(agg$SampleGroup))
      agg[order(agg$SampleGroup, agg$Group), , drop = FALSE]
    } else NULL

    # ---- 3. 梯度检验：group ~ 各设计因子（KW + Dunn + FDR，走 md_group_test）----
    gradient_test <- list()
    present_groups <- unique(stats::na.omit(comp_long$Group))
    for (fv in fac) {
      for (gv in present_groups) {
        sub <- comp_long[comp_long$Group == gv, , drop = FALSE]
        if (length(unique(stats::na.omit(sub[[fv]]))) < 2) next
        gt <- tryCatch(md_group_test(CFG, sub, "RelAbund", fv,
                                     test = CFG$alpha$test %||% "kruskal_dunn",
                                     p_adjust = CFG$alpha$p_adjust %||% "BH",
                                     min_n = 3, order_levels = CFG$design$ordered_levels[[fv]]),
                       error = function(e) NULL)
        if (is.null(gt) || is.null(gt$test) || nrow(gt$test) == 0) next
        tt <- gt$test
        tt$Group <- gv; tt$Factor <- fv; tt$n <- nrow(sub)
        gradient_test[[paste(fv, gv, sep = "::")]] <- tt
      }
    }
    gradient_test <- if (length(gradient_test)) {
      gt <- do.call(rbind, gradient_test)
      gt$p_adj_fdr <- md_p_adjust(gt$p, "BH")
      gt <- gt[, c("Group", "Factor", "test", "statistic", "df", "p", "method",
                   "effect", "p_adj_fdr", "n"), drop = FALSE]
      gt[order(gt$p_adj_fdr), , drop = FALSE]
    } else NULL

    # ---- 4. 网络角色（软依赖：80_network.rds）----
    net_path <- md_path(CFG, "cache", "80_network.rds")
    net <- if (file.exists(net_path)) tryCatch(readRDS(net_path), error = function(e) NULL) else NULL
    if (is.null(net)) add_note("未找到 cache/80_network.rds（或读取失败），跳过网络角色面板 A，仅出组成面板 B。")

    net_nodes <- if (!is.null(net) && !is.null(net$nodes)) as.data.frame(net$nodes) else NULL
    net_edges <- if (!is.null(net) && !is.null(net$edges)) as.data.frame(net$edges) else NULL
    if (!is.null(net_nodes) && !id_col %in% names(net_nodes) && "ASV_ID" %in% names(net_nodes))
      net_nodes[[id_col]] <- net_nodes$ASV_ID

    # 节点表列名归一（80 号产物字段固定，这里仍做一次容错）
    pick_col <- function(df, cands) {
      hit <- intersect(cands, names(df))
      if (length(hit)) hit[1] else NULL
    }
    role_tab <- NULL; centrality_by_group <- NULL
    if (!is.null(net_nodes) && id_col %in% names(net_nodes)) {
      grp_nodes <- ifelse(is.na(group_of[match(net_nodes[[id_col]], asv_ids)]),
                          unmapped_lab, group_of[match(net_nodes[[id_col]], asv_ids)])
      net_nodes$.Group <- grp_nodes
      role_col <- NULL
      if ("role" %in% names(net_nodes)) role_col <- "role"
      else {
        zc <- pick_col(net_nodes, c("Zi", "zi")); pc <- pick_col(net_nodes, c("Pi", "pi"))
        if (!is.null(zc) && !is.null(pc)) {
          net_nodes$role <- ifelse(net_nodes[[zc]] > MD_NET$zi_thr & net_nodes[[pc]] > MD_NET$pi_thr, "Network Hub",
                            ifelse(net_nodes[[zc]] > MD_NET$zi_thr, "Module Hub",
                            ifelse(net_nodes[[pc]] > MD_NET$pi_thr, "Connector", "Peripheral")))
          role_col <- "role"
        }
      }
      if (!is.null(role_col)) {
        keep_role <- !is.na(net_nodes$.Group)
        role_tab <- as.data.frame(table(Group = net_nodes$.Group[keep_role],
                                        Role = net_nodes[[role_col]][keep_role]),
                                  stringsAsFactors = FALSE)
        names(role_tab) <- c("Group", "Role", "n_node")
        tot_g <- tapply(role_tab$n_node, role_tab$Group, sum)
        role_tab$n_group_nodes <- as.integer(tot_g[role_tab$Group])
        role_tab$pct_node <- 100 * role_tab$n_node / pmax(role_tab$n_group_nodes, 1)
        role_tab$Role <- factor(role_tab$Role,
                                levels = intersect(names(MD_NET$role_colors),
                                                   unique(role_tab$Role)))
        role_tab <- role_tab[!is.na(role_tab$Role), , drop = FALSE]
        role_tab <- role_tab[order(role_tab$Group, role_tab$Role), , drop = FALSE]
      }
      # 各类群中心性汇总（只保留节点表里真实存在的列）
      num_cols <- intersect(c("degree", "Degree", "betweenness", "Betweenness",
                              "closeness", "Closeness", "eigen", "Eigen_centrality",
                              "Zi", "Pi", "module", "Module"), names(net_nodes))
      if (length(num_cols)) {
        centrality_by_group <- do.call(rbind, lapply(sort(unique(net_nodes$.Group)), function(g) {
          sub <- net_nodes[net_nodes$.Group == g, , drop = FALSE]
          row <- data.frame(Group = g, n_node = nrow(sub), stringsAsFactors = FALSE)
          for (cc in num_cols) row[[paste0("mean_", cc)]] <- mean(as.numeric(sub[[cc]]), na.rm = TRUE)
          row
        }))
      }
    }

    # ---- 5. 响应变量关联：group 相对丰度 vs response$var 的 Spearman ----
    response_cor <- NULL
    if (isTRUE(CFG$response$enabled) && !is.null(CFG$response$var)) {
      rv <- CFG$response$var
      env <- core$env
      if (!is.null(env) && rv %in% names(env)) {
        # ⚠ 用**局部副本**算相关：comp_long 是结果表 Composition 的数据源，
        #   往里加 resp 列会让结果表多出一列中间变量。
        cor_df <- comp_long
        cor_df$resp <- suppressWarnings(as.numeric(env[[rv]][match(cor_df$Sample, env[[smp]])]))
        response_cor <- do.call(rbind, lapply(unique(cor_df$Group), function(g) {
          sub <- cor_df[cor_df$Group == g, , drop = FALSE]
          sub <- sub[is.finite(sub$resp) & is.finite(sub$RelAbund), , drop = FALSE]
          if (nrow(sub) < 5) return(NULL)
          ct <- tryCatch(stats::cor.test(sub$RelAbund, sub$resp, method = "spearman",
                                         exact = FALSE), error = function(e) NULL)
          if (is.null(ct)) return(NULL)
          data.frame(Group = g, Response = rv, n = nrow(sub),
                     Spearman_rho = unname(ct$estimate), p_value = ct$p.value,
                     stringsAsFactors = FALSE)
        }))
        if (!is.null(response_cor) && nrow(response_cor)) {
          response_cor$p_adj_fdr <- md_p_adjust(response_cor$p_value, "BH")
          response_cor <- response_cor[order(response_cor$p_adj_fdr), , drop = FALSE]
        } else response_cor <- NULL
      } else {
        add_note("response$var = ", CFG$response$var, " 不在 env 表里，跳过响应变量关联。")
      }
    }

    # ---- 6. 跨组互作强度：组内边 vs 组间边（Fisher 精确检验）----
    edge_test <- NULL
    edge_long <- NULL
    if (!is.null(net_edges) && nrow(net_edges) > 0) {
      ec_from <- pick_col(net_edges, c("from", "From", "source", "ASV_1", "asv1"))
      ec_to   <- pick_col(net_edges, c("to", "To", "target", "ASV_2", "asv2"))
      if (!is.null(ec_from) && !is.null(ec_to)) {
        g1 <- ifelse(is.na(group_of[match(net_edges[[ec_from]], asv_ids)]), unmapped_lab,
                     group_of[match(net_edges[[ec_from]], asv_ids)])
        g2 <- ifelse(is.na(group_of[match(net_edges[[ec_to]], asv_ids)]), unmapped_lab,
                     group_of[match(net_edges[[ec_to]], asv_ids)])
        keep_e <- !is.na(g1) & !is.na(g2)
        same <- g1[keep_e] == g2[keep_e]
        obs_within <- sum(same); obs_between <- sum(!same)
        # 零模型期望 = 按**网络节点的分组频率**随机连线（与原文 C9 同口径）
        net_asvs <- unique(c(as.character(net_edges[[ec_from]][keep_e]),
                             as.character(net_edges[[ec_to]][keep_e])))
        ng <- ifelse(is.na(group_of[match(net_asvs, asv_ids)]), unmapped_lab,
                     group_of[match(net_asvs, asv_ids)])
        freq <- table(ng) / length(ng)
        p_within <- sum(freq^2)
        n_e <- sum(keep_e)
        exp_within <- n_e * p_within
        exp_between <- n_e - exp_within
        ft <- tryCatch(stats::fisher.test(matrix(c(obs_within, obs_between,
                                                   round(exp_within), round(exp_between)),
                                                 nrow = 2, byrow = TRUE)),
                       error = function(e) NULL)
        edge_test <- data.frame(
          Metric = c("Total_edges", "Classifiable_edges", "Within_group_edges",
                     "Between_group_edges", "Expected_within_group", "Expected_between_group"),
          Value  = c(nrow(net_edges), n_e, obs_within, obs_between,
                     round(exp_within, 1), round(exp_between, 1)),
          stringsAsFactors = FALSE)
        edge_test <- rbind(edge_test, data.frame(
          Metric = c("Fisher_odds_ratio", "Fisher_p_value", "Test"),
          Value  = c(if (is.null(ft)) NA_real_ else unname(ft$estimate),
                     if (is.null(ft)) NA_real_ else ft$p.value,
                     if (is.null(ft)) NA else "Fisher exact (observed vs expected within/between group edges, 2x2)"),
          stringsAsFactors = FALSE))
        edge_long <- data.frame(from = net_edges[[ec_from]][keep_e],
                                to = net_edges[[ec_to]][keep_e],
                                Group_from = g1[keep_e], Group_to = g2[keep_e],
                                Type = ifelse(same, "Within_group", "Between_group"),
                                stringsAsFactors = FALSE)
        if (!is.null(pick_col(net_edges, c("weight", "Weight")))) {
          wc <- pick_col(net_edges, c("weight", "Weight"))
          edge_long$weight <- as.numeric(net_edges[[wc]][keep_e])
        }
      }
    }

    # ---- 7. 类群指示种（IndVal）----
    indicator <- NULL
    ind_cfg <- tg$indicator %||% list()
    if (isTRUE(ind_cfg$enabled) && !is.null(fac1) && length(grp_levels) >= 2) {
      comm_rel <- md_rel(m_keep)
      grp_vec  <- ifelse(is.na(group_of), unmapped_lab, group_of)
      names(grp_vec) <- asv_ids
      if (identical(toupper(match_by), "ASV") || !is.null(core$tax)) {
        comm_species <- if (identical(toupper(match_by), "ASV")) m_keep else {
          # 非 ASV 层级：在该层级上聚合，再按该层级的 group 归属做指示种
          rank_vec <- tax_key
          rank_vec[is.na(rank_vec) | !nzchar(rank_vec)] <- "Unclassified"
          rowsum(m_keep, group = rank_vec[keep_asv], reorder = TRUE)
        }
        # 注意：聚合发生在**类群层**，列仍是样本；所以两种情形都用 m_keep 的列名
        grp_for_ind <- as.character(meta[[fac1]][match(colnames(m_keep), meta[[smp]])])
        if (length(grp_for_ind) == ncol(comm_species)) {
          iv <- md_indval(CFG, comm_species, grp_for_ind,
                          func = ind_cfg$func %||% "IndVal.g",
                          nperm = as.integer(ind_cfg$nperm %||% 999),
                          adjust = ind_cfg$adjust %||% "BH")
          if (!is.null(iv) && nrow(iv)) {
            indicator <- iv
            top_n <- as.integer(ind_cfg$top_n %||% 20)
            fdr_thr <- as.numeric(ind_cfg$fdr_thr %||% 0.05)
            # 先按 IndVal 降序取 top_n（**不过滤**）落到表里，让用户能看到未达标的类群；
            # 同时用 fdr_sig 标出哪些通过阈值。一个都不过时仍然写表 + 记 note，
            # 而不是把整张 sheet 丢掉（否则"跑过了没结果"与"根本没跑"分不清）。
            indicator <- indicator[order(-indicator$IndVal), , drop = FALSE]
            indicator$fdr_sig <- !is.na(indicator$fdr) & indicator$fdr < fdr_thr
            if (nrow(indicator) > top_n) indicator <- utils::head(indicator, top_n)
            if (!any(indicator$fdr_sig)) {
              add_note("IndVal 没有通过 fdr < ", fdr_thr, " 的类群；",
                       "表里仍保留按 IndVal 降序的前 ", nrow(indicator),
                       " 个（fdr_sig = FALSE）供查阅。")
            }
          } else {
            add_note("IndVal 无结果（indicspecies 缺失或分组不足）。")
          }
        } else {
          add_note("IndVal 跳过：样本对齐失败。")
        }
      }
    }

    nf <- function(x) as.character(round(as.numeric(x), digits))
    coverage <- data.frame(
      item = c("n_ASV_total", "n_ASV_mapped", "pct_ASV_mapped",
               "pct_reads_mapped", "match_by",
               "reads_mapped", "reads_total"),
      value = c(nf(length(asv_ids)), nf(n_mapped),
                nf(100 * n_mapped / max(1L, length(asv_ids))),
                nf(if (sum(m) > 0) 100 * sum(m_keep) / sum(m) else NA_real_),
                as.character(match_level),
                nf(sum(m_keep)), nf(sum(m))),
      stringsAsFactors = FALSE)
  })

  if (n_mapped == 0) {
    md_skip(CFG, module, paste0("映射表与数据没有任何匹配（match_by = ",
                                match_by, "），检查匹配层级或表内容"))
    return(md_result(module, "SKIP", "映射表零匹配"))
  }
  if (n_mapped < length(asv_ids))
    add_note("有 ", length(asv_ids) - n_mapped, " / ", length(asv_ids), " 个 ASV 未匹配到分组",
             if (drop_un) "，按 drop_unmapped = TRUE 已从组成中剔除。"
             else paste0("，已归入 \"", unmapped_lab, "\"。"))

  # ---- 配色：走 md_pal()（CFG$figure$palette[[label]] -> 内置板 -> 定性板兜底）----
  group_cols <- st$palette %||% md_pal(CFG, levels(comp_mean$Group) %||% grp_all, name = label)
  if (!is.null(st$palette) && is.null(names(st$palette)))
    group_cols <- stats::setNames(rep(st$palette, length.out = length(grp_all)), grp_all)
  role_cols <- MD_NET$role_colors

  # ---------------------------------------------------------------------------
  # 出图段（每次运行都执行；样式不进缓存键）
  # ---------------------------------------------------------------------------
  panels <- list()
  if (!is.null(role_tab) && nrow(role_tab) > 0) {
    role_tab$Group <- factor(role_tab$Group, levels = unique(as.character(role_tab$Group)))
    pA <- ggplot2::ggplot(role_tab,
                          ggplot2::aes(x = .data$Group, y = .data$pct_node, fill = .data$Role)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.65,
                        colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = role_cols, name = "Node role", drop = FALSE,
                                 guide = md_guide_flex(length(role_cols), "bottom")) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 101),
                                  labels = function(x) paste0(x, "%")) +
      ggplot2::scale_x_discrete(labels = function(v) gsub("_", "\n", v, fixed = TRUE)) +
      ggplot2::labs(x = label, y = "% of network nodes", tag = "A") +
      theme_md(CFG = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = pt$tick),
                     legend.position = "bottom", legend.box = "vertical")
    panels$A <- pA
  } else {
    add_note("网络角色数据不可用（80 未运行或节点表无角色列），面板 A 缺失。")
  }

  if (is.null(comp_mean) || nrow(comp_mean) == 0) {
    pB <- md_empty_plot(CFG, "No group composition")
  } else {
    pB <- ggplot2::ggplot(comp_mean,
                          ggplot2::aes(x = .data$SampleGroup, y = .data$RelAbund,
                                       fill = .data$Group)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.65,
                        colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = group_cols, name = label, drop = FALSE,
                                 guide = md_guide_flex(length(group_cols), "bottom", nrow = 2)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 101),
                                  labels = function(x) paste0(x, "%")) +
      ggplot2::scale_x_discrete(labels = function(v) gsub("_", "\n", v, fixed = TRUE)) +
      ggplot2::labs(x = fac1 %||% NULL, y = "Relative abundance (%)",
                    tag = if (length(panels)) "B" else "A") +
      theme_md(CFG = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = pt$tick),
                     legend.position = "bottom", legend.box = "vertical")
  }
  panels[[if (length(panels)) "B" else "A"]] <- pB

  ph <- as.numeric(st$panel_heights %||% c(1, 1.2))
  # guides = "keep"：两个面板都用了 fill，collect 会把「Node role」与「Group」的图例
  # 合并成同一个图例（6+ 个条目混在一起，读者分不清哪条属于哪个面板）；原文也是各自独立。
  fig <- if (length(panels) >= 2) {
    md_compose(CFG, panels, design = "AABB", guides = "keep", tag_levels = NULL)
  } else {
    md_compose(CFG, panels, design = "A", guides = "keep", tag_levels = NULL)
  }
  md_save_figure(CFG, fig, "81_taxon_group",
                 width = 183, height = if (length(panels) >= 2) 150 else 100)
  for (nm in names(panels)) md_maybe_single(CFG, panels[[nm]], paste0("81_taxon_group_", nm))

  # ---------------------------------------------------------------------------
  # 写表段
  # ---------------------------------------------------------------------------
  comp_out <- comp_long[, c("Sample", "SampleGroup", "Group", "RelAbund",
                            setdiff(fac, fac1)), drop = FALSE]
  comp_out$RelAbund <- round(comp_out$RelAbund, digits)

  note <- data.frame(
    item = c("Mapping file", "taxon_col", "group_col", "match_by", "group_label",
             "Group levels", "drop_unmapped", "Unmapped label", "Composition matrix",
             "Normalization (composition)", "Panel A", "Panel B", "Edge test",
             "Original script difference"),
    value = c(map_path, tcol, gcol, match_level, label,
              paste(levels(comp_mean$Group) %||% grp_all, collapse = ", "),
              as.character(drop_un), unmapped_lab,
              paste0(nrow(m_keep), " ASV x ", ncol(m_keep), " samples"),
              as.character(CFG$normalize$by_module$composition %||% CFG$normalize$method),
              if (!is.null(role_tab)) "Percentage of network nodes per group by Zi-Pi role" else "skipped (no 80_network.rds)",
              "Mean per-sample relative abundance (%) of each group, by sample grouping factor",
              "Fisher exact test on observed within-group vs between-group edges against a random-wiring expectation built from the node group frequencies",
              paste0("The hard-coded diazotroph genus -> Symbiotic/Associative/Free-living table ",
                     "of 6_ecological_guilds.R:102-268 was REMOVED; grouping now comes entirely ",
                     "from the user-supplied mapping table. Group meaning is arbitrary ",
                     "(functional guild / niche / R-K strategy / any custom scheme).")),
    stringsAsFactors = FALSE)

  tp <- md_write_xlsx(CFG, "81_taxon_group", list(
    Mapping            = mapping_tab,
    Composition        = comp_out,
    Gradient_test      = gradient_test,
    Network_centrality = centrality_by_group,
    Response_correlation = response_cor,
    Edge_test          = edge_test,
    Indicator_species  = indicator,
    Mapping_coverage   = coverage,
    Note               = note))

  md_result(module, "OK", NULL,
            tables  = "81_taxon_group.xlsx",
            figures = "figures/81_taxon_group.pdf",
            notes = c(notes,
                      paste0("match_by=", match_level, "; groups=", length(grp_all),
                             "; mapped ASV=", n_mapped, "/", length(asv_ids)),
                      if (is.null(role_tab)) "面板 A（网络角色）缺失，整图只有组成面板 B"
                      else "面板 A + B 齐全"))
}
