# =============================================================================
# R/80_network.R ---- 模块 80：共现网络（全局拓扑 + 节点角色 + 分组网络 + 亚网络）
# -----------------------------------------------------------------------------
# 源：毕业论文 `5_network_structure.R`（物理合并了 12_network_structure.R 与
#     13_abundant_rare_network.R 两个 local() 成员块）—— 图 Fig5a / Fig5b /
#     FigS7 / Fig5c / Fig5d / Fig5e + 原脚本后来并入的 Sankey 面板。
# 契约：_build-notes/micro-div-pipeline/ROUND5_MODULE_SPEC.md
#
# 口径（逐项照抄原脚本，不做"顺手优化"）：
#   建网         = 相对丰度(%) -> 出现率过滤 -> Spearman 相关 -> |r| 阈值稀疏化
#                  -> graph_from_adjacency_matrix(加权、无向) -> 边权保留正负号
#   节点中心性   = get_centrality()：先 simplify，再取**最大连通组分**后算
#                  betweenness / closeness / degree(归一) / eigen / transitivity
#   样本级指标   = 每个样本检出的 ASV 在全局网络上诱导子图的拓扑 + 正负内聚力
#   网络韧性     = robustness_R()：随机逐次去除节点直到最大组分 < 50%，R = 去除比例
#   模块 / Zi-Pi = Louvain(权重 |w|) + Deng et al. 2012 的 Zi/Pi 与四类角色
#   随机网络     = Erdos-Renyi sample_gnm(n, m)，Z > 1.96 / < -1.96 判高/低于随机
#   亚网络       = 丰度 x 出现率 2x2 四象限（Core / Permanently rare /
#                  Conditionally rare / Transient rare）x 分组因子
#   节点周转     = 各分组网络节点集合的 Jaccard 重叠 + 相邻水平间 共有/丢失/新增
#
# 出图（4 张）：
#   figures/80_network_nodes.pdf      A betweenness top-N | B 度-丰度 | C Zi-Pi
#   figures/80_network_metrics.pdf    metrics_main 指标 x 前 1~2 个因子
#   figures/80_network_metrics_supp.pdf  metrics_supp（仅当 metrics_supp_figure = TRUE）
#   figures/80_network_topology.pdf   A 全局网络 | B 分组网络 | C vs 随机 |
#                                     D 节点周转 Sankey | E 丰富/稀有亚网络
# 出表：results/80_network.xlsx
# 跨模块缓存：cache/80_network.rds（nodes / edges / params），供模块 81 读
#
# 与原文的差异（每一处都在此登记，理由见对应代码处的注释）：
#   1. **修 bug**：原脚本 `set.seed(2026)` 排在 robustness_R() 的 30 次模拟**之后**，
#      导致样本级 Robustness 每次运行都不同。本模块在调用 robustness_R() 之前
#      先 `md_seed(CFG, md_offset(CFG, "network_robustness", 8220L))`。
#   2. 网络布局（图 80_network_topology 的 A / B 面板）是随机过程，出图段先
#      `md_seed(CFG, md_offset(CFG, "network_layout", 8230L))` 再算坐标。
#   3. 所有图上文字纯 ASCII；度量名按 ROUND5_MODULE_SPEC 第 5 节的替换表转英文。
#   4. 原脚本靠 `00_environment.RData` 注入 cols_zone / cols_veg / cols_2x2 /
#      theme_isme() / S_* 字号 / stat_cor 等；本模块一律改走 R/shared/ 的
#      theme_md() / md_pt() / md_pal() / md_pt_size 常量，不复制那些函数。
#   5. 原脚本的 6 张图在 183 mm 版心下字号只剩约 5 pt 可读；本模块统一走
#      主图字号轨道（md_pt(CFG)），与全包其余模块的"一图一号字梯"保持一致
#      （原图 533.4 mm 版幅的 80_network_nodes.pdf 仍按原文 9 pt 主图轨道）。
# =============================================================================

# -----------------------------------------------------------------------------
# 模块内工具（local_ 前缀；不改共享库）
# -----------------------------------------------------------------------------

#' 因子 x 标签度量名（ASCII）。键 = sample_metrics 列名。
local80_metric_labels <- function() {
  c(
    N_DetectedASV     = "Network Size",
    MeanNodeDegree    = "Mean Node Degree",
    SumNodeDegree     = "Sum Node Degree",
    MeanBetweenness   = "Betweenness",
    MeanCloseness     = "Closeness",
    MeanEigen         = "Eigenvector Centrality",
    MeanWithinCorr    = "Mean Within-Corr",
    Density_sample    = "Density",
    Clustering_sample = "Clustering Coef.",
    Connectedness_sample = "Connectedness",
    Modularity_sample = "Modularity",
    SubDensity        = "Density",
    SubClustering     = "Clustering Coef.",
    SubConnectedness  = "Connectedness",
    SubAvgPath        = "Avg Path Length",
    SubHeterogeneity  = "Heterogeneity",
    SubModularity     = "Modularity",
    SubNModules       = "N Modules",
    Cohesion_pos      = "Positive Cohesion",
    Cohesion_neg      = "Negative Cohesion",
    Robustness        = "Robustness"
  )
}

#' 度量名 -> 显示名（表里没登记的键原样返回）
#' Sub* / *_sample 两套度量语义不同（前者 = 样本内诱导子图，后者 = 整张全局网络），
#' 但原脚本给它们的显示名是一样的（都叫 Density / Clustering / ...）。若直接照抄，
#' 同一张图的 facet 会出现两个同名面板、读者无法分辨。故在**整张网**口径的名字后
#' 补 "(network)" 加以区分 —— 这是对原文显示名的唯一改动，数值口径未动。
local80_lab <- function(m, labels = local80_metric_labels()) {
  v <- unname(labels[as.character(m)])
  v[is.na(v)] <- as.character(m)[is.na(v)]
  v <- ifelse(as.character(m) %in% c("Density_sample", "Clustering_sample",
                                     "Connectedness_sample", "Modularity_sample"),
              paste0(v, " (network)"), v)
  v
}

#' 取某因子的有序水平（与 10/11 模块同口径）
local80_levels <- function(CFG, meta, var) {
  ol  <- CFG$design$ordered_levels[[var]]
  obs <- unique(as.character(meta[[var]]))
  obs <- obs[!is.na(obs)]
  if (!is.null(ol)) {
    hit <- intersect(as.character(ol), obs)
    if (length(hit)) return(hit)
  }
  obs
}

#' 分组配色：优先内置具名板（zone / veg），否则 md_pal() 兜底
local80_group_cols <- function(CFG, levels, name) {
  p <- CFG$figure$palette[[name]]
  if (!is.null(p) && !is.null(names(p)) && all(levels %in% names(p))) return(p[levels])
  b <- md_pick_board(name, levels)
  if (!is.null(b) && all(levels %in% names(b))) return(b[levels])
  md_pal(CFG, levels, name)
}

# -----------------------------------------------------------------------------
# 建网 / 拓扑 / 节点角色（逐行照抄原脚本）
# -----------------------------------------------------------------------------

#' 全局网络构建（源：5_network_structure.R:142-175）
#' @param otu_mat ASV x 样本 **相对丰度（%）**
local80_calc_network <- function(otu_mat, min_cor = 0.5, min_prev = 0.02) {
  prev <- rowSums(otu_mat > 0.01) / ncol(otu_mat)
  keep <- prev >= min_prev
  mat  <- otu_mat[keep, , drop = FALSE]
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (nrow(mat) < 3) return(NULL)
  cor_mat <- stats::cor(t(mat), method = "spearman")
  cor_mat[abs(cor_mat) < min_cor] <- 0
  diag(cor_mat) <- 0
  g <- igraph::graph_from_adjacency_matrix(abs(cor_mat) * (abs(cor_mat) > 0),
                                           mode = "undirected", weighted = TRUE,
                                           diag = FALSE)
  if (igraph::ecount(g) == 0) return(NULL)
  el <- igraph::as_edgelist(g, names = FALSE)
  edge_weights <- vapply(seq_len(igraph::ecount(g)),
                         function(i) cor_mat[el[i, 1], el[i, 2]], numeric(1))
  igraph::E(g)$weight <- edge_weights
  igraph::E(g)$sign   <- ifelse(edge_weights > 0, "pos", "neg")
  list(g = g, cor_mat = cor_mat, kept_asv = rownames(mat))
}

#' 分组/亚网络构建（源：5_network_structure.R:790-822 与 1910-1928）
#' 与全局网络的唯一差别：先剔除零方差 ASV、样本数 < 4 直接 NULL
local80_group_net <- function(otu_mat, min_cor = 0.5, min_prev = 0.02) {
  if (is.null(otu_mat) || ncol(otu_mat) < 4) return(NULL)
  prev <- rowSums(otu_mat > 0.01) / ncol(otu_mat)
  mat  <- otu_mat[prev >= min_prev, , drop = FALSE]
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (nrow(mat) < 3) return(NULL)
  keep_var <- apply(mat, 1, stats::var, na.rm = TRUE) > 0
  keep_var[is.na(keep_var)] <- FALSE
  mat <- mat[keep_var, , drop = FALSE]
  if (nrow(mat) < 3) return(NULL)
  cor_mat <- stats::cor(t(mat), method = "spearman")
  cor_mat[abs(cor_mat) < min_cor] <- 0
  diag(cor_mat) <- 0
  g <- igraph::graph_from_adjacency_matrix(abs(cor_mat) * (abs(cor_mat) > 0),
                                           mode = "undirected", weighted = TRUE,
                                           diag = FALSE)
  if (igraph::ecount(g) == 0) return(NULL)
  el <- igraph::as_edgelist(g, names = FALSE)
  w  <- vapply(seq_len(igraph::ecount(g)),
               function(i) cor_mat[el[i, 1], el[i, 2]], numeric(1))
  igraph::E(g)$weight <- w
  igraph::E(g)$sign   <- ifelse(w > 0, "pos", "neg")
  g
}

#' 节点中心性（源：5_network_structure.R:205-241）
local80_centrality <- function(g) {
  empty <- data.frame(ASV_ID = character(0), Betweenness = numeric(0),
                      Closeness = numeric(0), Degree = numeric(0),
                      Eigen_centrality = numeric(0), Transitivity = numeric(0),
                      stringsAsFactors = FALSE)
  if (is.null(g) || igraph::vcount(g) < 3) return(empty)
  g <- igraph::simplify(g, remove.loops = TRUE, remove.multiple = FALSE)
  if (igraph::vcount(g) < 3) return(empty)
  if (!igraph::is_connected(g)) {
    comp <- igraph::components(g)
    gcc_idx <- which(comp$membership == which.max(comp$csize))
    g_cc <- igraph::induced_subgraph(g, gcc_idx)
  } else {
    g_cc <- g
  }
  btw <- tryCatch(igraph::betweenness(g_cc, directed = FALSE, normalized = TRUE),
                  error = function(e) rep(NA_real_, igraph::vcount(g_cc)))
  clo <- tryCatch(igraph::closeness(g_cc, directed = FALSE, normalized = TRUE),
                  error = function(e) rep(NA_real_, igraph::vcount(g_cc)))
  deg <- tryCatch(igraph::degree(g_cc, normalized = TRUE),
                  error = function(e) rep(NA_real_, igraph::vcount(g_cc)))
  eig <- tryCatch(igraph::eigen_centrality(g_cc, directed = FALSE)$vector,
                  error = function(e) rep(NA_real_, igraph::vcount(g_cc)))
  cc  <- tryCatch(igraph::transitivity(g_cc, type = "global"), error = function(e) NA_real_)
  data.frame(ASV_ID = igraph::V(g_cc)$name,
             Betweenness = as.numeric(btw),
             Closeness = as.numeric(clo),
             Degree = as.numeric(deg),
             Eigen_centrality = as.numeric(eig),
             Transitivity = rep(cc, igraph::vcount(g_cc)),
             stringsAsFactors = FALSE)
}

#' 网络韧性 R（源：5_network_structure.R:268-281）
#' ⚠ 原脚本的 set.seed(2026) 排在本函数**之外、30 次模拟之后**（L401），
#'   造成 Fig5b 不可复现。本模块要求调用方在进入样本循环**之前**先 md_seed()。
local80_robustness <- function(g, nsim = 30L) {
  n <- igraph::vcount(g)
  if (n < 3) return(NA_real_)
  outs <- numeric(nsim)
  for (s in seq_len(nsim)) {
    idx <- sample(seq_len(n), n)
    target <- 0.5 * n
    for (r in seq_len(n - 1)) {
      gg <- igraph::induced_subgraph(g, idx[-(1:r)])
      if (igraph::vcount(gg) == 0) { outs[s] <- r / n; break }
      if (max(igraph::components(gg)$csize) < target) { outs[s] <- r / n; break }
      if (r == n - 1) outs[s] <- 1
    }
  }
  mean(outs)
}

#' 样本级网络指标（源：5_network_structure.R:283-384）
local80_sample_metrics <- function(g, otu_mat, meta, sample_col,
                                   factors, robustness_nsim = 30L,
                                   group_net = FALSE) {
  ids  <- colnames(otu_mat)
  node_names <- igraph::V(g)$name
  out <- data.frame(SampleID = ids, stringsAsFactors = FALSE)
  for (f in factors) {
    v <- as.character(meta[[f]][match(ids, as.character(meta[[sample_col]]))])
    out[[f]] <- v
  }
  out$MeanNodeDegree  <- NA_real_
  out$SumNodeDegree   <- NA_real_
  out$MeanBetweenness <- NA_real_
  out$MeanCloseness   <- NA_real_
  out$MeanEigen       <- NA_real_
  out$MeanWithinCorr  <- NA_real_
  out$N_DetectedASV   <- NA_integer_
  out$SubDensity      <- NA_real_
  out$SubClustering   <- NA_real_
  out$SubConnectedness <- NA_real_
  out$SubAvgPath      <- NA_real_
  out$SubHeterogeneity <- NA_real_
  out$SubModularity   <- NA_real_
  out$SubNModules     <- NA_integer_
  out$Cohesion_pos    <- NA_real_
  out$Cohesion_neg    <- NA_real_
  out$Robustness      <- NA_real_

  for (i in seq_along(ids)) {
    present_in_net <- intersect(names(which(otu_mat[, i] > 0)), node_names)
    n_det <- length(present_in_net)
    out$N_DetectedASV[i] <- n_det
    if (n_det < 2) next

    v_idx <- which(node_names %in% present_in_net)
    sub_g <- igraph::induced_subgraph(g, v_idx)

    out$MeanNodeDegree[i] <- mean(igraph::degree(sub_g, normalized = TRUE), na.rm = TRUE)
    out$SumNodeDegree[i]  <- sum(igraph::degree(sub_g, normalized = FALSE), na.rm = TRUE)

    btw_sub <- tryCatch(igraph::betweenness(sub_g, directed = FALSE, normalized = TRUE),
                        error = function(e) NA_real_)
    clo_sub <- tryCatch(igraph::closeness(sub_g, directed = FALSE, normalized = TRUE),
                        error = function(e) NA_real_)
    eig_sub <- tryCatch(igraph::eigen_centrality(sub_g, directed = FALSE)$vector,
                        error = function(e) NA_real_)
    out$MeanBetweenness[i] <- mean(btw_sub, na.rm = TRUE)
    out$MeanCloseness[i]   <- mean(clo_sub, na.rm = TRUE)
    out$MeanEigen[i]       <- mean(eig_sub, na.rm = TRUE)

    adj_sub <- igraph::as_adjacency_matrix(sub_g, type = "upper", attr = "weight",
                                           sparse = FALSE)
    diag(adj_sub) <- 0
    wts <- adj_sub[upper.tri(adj_sub)]
    wts <- wts[wts != 0]
    out$MeanWithinCorr[i] <- if (length(wts) > 0) mean(wts) else NA_real_

    e_sign <- igraph::E(sub_g)$sign
    e_wt   <- igraph::E(sub_g)$weight
    pos_w  <- e_wt[e_sign == "pos"]
    neg_w  <- e_wt[e_sign == "neg"]
    out$Cohesion_pos[i] <- if (length(pos_w) > 0) mean(pos_w, na.rm = TRUE) else NA_real_
    out$Cohesion_neg[i] <- if (length(neg_w) > 0) mean(abs(neg_w), na.rm = TRUE) else NA_real_
    if (n_det >= 3) out$Robustness[i] <- local80_robustness(sub_g, nsim = robustness_nsim)

    out$SubDensity[i] <- tryCatch(igraph::edge_density(sub_g), error = function(e) NA_real_)
    out$SubClustering[i] <- tryCatch(igraph::transitivity(sub_g, type = "global"),
                                     error = function(e) NA_real_)
    comp_sub <- tryCatch(igraph::components(sub_g), error = function(e) NULL)
    out$SubConnectedness[i] <- if (!is.null(comp_sub))
      max(comp_sub$csize) / igraph::vcount(sub_g) else NA_real_
    out$SubAvgPath[i] <- tryCatch(igraph::mean_distance(sub_g), error = function(e) NA_real_)
    deg_r <- igraph::degree(sub_g, normalized = FALSE)
    out$SubHeterogeneity[i] <- if (mean(deg_r) > 0) stats::sd(deg_r) / mean(deg_r) else NA_real_
    louv_sub <- tryCatch(igraph::cluster_louvain(sub_g), error = function(e) NULL)
    out$SubModularity[i] <- if (!is.null(louv_sub)) igraph::modularity(louv_sub) else NA_real_
    out$SubNModules[i]   <- if (!is.null(louv_sub))
      length(unique(igraph::membership(louv_sub))) else NA_integer_
  }
  # 全局网络口径下，把"整张网的等价指标"另存一份前缀列，便于与分组网络对照
  if (isTRUE(group_net)) {
    out$Density_sample       <- igraph::edge_density(g)
    out$Clustering_sample    <- igraph::transitivity(g, type = "global")
    cmp <- igraph::components(g)
    out$Connectedness_sample <- if (cmp$no == 1) 1 else
      sum(vapply(cmp$csize, function(s) s * (s - 1) / 2, numeric(1))) /
      (igraph::vcount(g) * (igraph::vcount(g) - 1) / 2)
    out$Modularity_sample    <- NA_real_
  }
  out
}

#' 单组网络拓扑指标（源：5_network_structure.R:825-885）
local80_group_topology <- function(g, group_name, group_type, n_samples) {
  empty_row <- data.frame(
    Group = group_name, Type = group_type, N_Samples = n_samples,
    Nodes = 0, Edges = 0, Avg_Degree = NA_real_, Density = NA_real_,
    Clustering = NA_real_, Connectedness = NA_real_, Avg_Path_Length = NA_real_,
    Modularity = NA_real_, N_Modules = NA_integer_, Heterogeneity = NA_real_,
    Natural_Conn = NA_real_, N_Pos = 0, N_Neg = 0, stringsAsFactors = FALSE)
  if (is.null(g) || igraph::vcount(g) < 2) return(empty_row)

  deg_all <- igraph::degree(g, normalized = FALSE)
  n_pos <- sum(igraph::E(g)$sign == "pos", na.rm = TRUE)
  n_neg <- sum(igraph::E(g)$sign == "neg", na.rm = TRUE)
  n_nodes <- igraph::vcount(g)

  connectedness <- tryCatch({
    comp <- igraph::components(g)
    if (comp$no == 1) 1.0 else
      sum(vapply(comp$csize, function(s) s * (s - 1) / 2, numeric(1))) /
      (n_nodes * (n_nodes - 1) / 2)
  }, error = function(e) NA_real_)
  avg_path <- tryCatch(igraph::mean_distance(g, directed = FALSE),
                       error = function(e) NA_real_)
  mod_result <- tryCatch(igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight)),
                         error = function(e) NULL)
  modularity_val <- if (!is.null(mod_result)) round(igraph::modularity(mod_result), 4) else NA_real_
  n_modules <- if (!is.null(mod_result))
    length(unique(igraph::membership(mod_result))) else NA_integer_
  heterogeneity <- if (mean(deg_all) > 0) round(stats::sd(deg_all) / mean(deg_all), 4) else NA_real_
  natural_conn <- tryCatch({
    adj <- igraph::as_adjacency_matrix(g, sparse = FALSE)
    ev <- Re(eigen(adj, only.values = TRUE)$values)
    ev <- ev[ev > 0]
    if (length(ev) > 0) log(mean(exp(ev))) else 0
  }, error = function(e) NA_real_)

  data.frame(
    Group = group_name, Type = group_type, N_Samples = n_samples,
    Nodes = n_nodes,
    Edges = igraph::ecount(g),
    Avg_Degree = round(mean(deg_all), 2),
    Density = round(igraph::edge_density(g), 4),
    Clustering = round(igraph::transitivity(g, type = "global"), 4),
    Connectedness = round(connectedness, 4),
    Avg_Path_Length = round(avg_path, 3),
    Modularity = modularity_val,
    N_Modules = n_modules,
    Heterogeneity = heterogeneity,
    Natural_Conn = round(natural_conn, 4),
    N_Pos = n_pos, N_Neg = n_neg,
    stringsAsFactors = FALSE)
}

#' Zi-Pi 节点角色（源：5_network_structure.R:1024-1058）
local80_zi_pi <- function(g, zi_thr = 2.5, pi_thr = 0.62) {
  modules <- igraph::V(g)$module
  deg_all <- igraph::degree(g)
  n <- igraph::vcount(g)
  zi_vec <- numeric(n); pi_vec <- numeric(n)
  for (i in seq_len(n)) {
    m_i <- modules[i]
    nbrs <- igraph::neighbors(g, i)
    all_in_m <- which(modules == m_i)
    deg_in_m <- vapply(all_in_m, function(v)
      sum(modules[igraph::neighbors(g, v)] == m_i), numeric(1))
    mu_m  <- mean(deg_in_m)
    sig_m <- stats::sd(deg_in_m)
    zi_vec[i] <- if (is.na(sig_m) || sig_m == 0) 0 else
      (sum(modules[nbrs] == m_i) - mu_m) / sig_m
    k_i <- length(nbrs)
    pi_vec[i] <- if (k_i > 0) { md <- table(modules[nbrs]); 1 - sum((md / k_i)^2) } else 0
  }
  out <- data.frame(ASV_ID = igraph::V(g)$name, degree = as.numeric(deg_all),
                    module = as.character(modules),
                    Zi = round(zi_vec, 4), Pi = round(pi_vec, 4),
                    stringsAsFactors = FALSE)
  out$role <- ifelse(out$Zi > zi_thr & out$Pi > pi_thr, "Network Hub",
              ifelse(out$Zi > zi_thr, "Module Hub",
              ifelse(out$Pi > pi_thr, "Connector", "Peripheral")))
  out
}

#' 2x2 丰度 x 出现率四象限（源：5_network_structure.R:1884-1893 的 classify_2x2）
#' 复用共享库 md_niche_table()（同一口径，含 Levins' B 便于附表），阈值走 CFG。
local80_cells <- function(rel, a_thr, p_thr) {
  nt <- md_niche_table(rel, a_thr = a_thr, p_thr = p_thr, index = "levins")
  data.frame(ASV_ID = nt$ASV_ID, Mean_RA = nt$Mean_Abund,
             Occurrence = nt$Occurrence / 100, Biosphere = nt$Cell,
             stringsAsFactors = FALSE)
}

#' 与随机网络比对（源：5_network_structure.R:1109-1149）
local80_random <- function(g, n_random = 100L) {
  n_nodes <- igraph::vcount(g); n_edges <- igraph::ecount(g)
  rs <- data.frame(avg_degree = numeric(n_random), density = numeric(n_random),
                   avg_path = numeric(n_random), clust_co = numeric(n_random),
                   modularity = numeric(n_random))
  for (i in seq_len(n_random)) {
    gr <- igraph::sample_gnm(n_nodes, n_edges)
    rs$avg_degree[i] <- mean(igraph::degree(gr, normalized = FALSE))
    rs$density[i]    <- igraph::edge_density(gr)
    rs$avg_path[i]   <- tryCatch(igraph::mean_distance(gr), error = function(e) NA_real_)
    rs$clust_co[i]   <- tryCatch(igraph::transitivity(gr, type = "global"),
                                 error = function(e) NA_real_)
    rs$modularity[i] <- igraph::modularity(igraph::cluster_louvain(gr))
  }
  cmp <- data.frame(
    Metric = c("Avg Degree", "Density", "Avg Path", "Clustering Coef", "Modularity"),
    Real = c(mean(igraph::degree(g, normalized = FALSE)), igraph::edge_density(g),
             tryCatch(igraph::mean_distance(g, directed = FALSE), error = function(e) NA_real_),
             igraph::transitivity(g, type = "global"),
             igraph::modularity(igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight)))),
    Rand_Mean = c(mean(rs$avg_degree), mean(rs$density),
                  mean(rs$avg_path, na.rm = TRUE),
                  mean(rs$clust_co, na.rm = TRUE), mean(rs$modularity)),
    Rand_SD = c(stats::sd(rs$avg_degree), stats::sd(rs$density),
                stats::sd(rs$avg_path, na.rm = TRUE),
                stats::sd(rs$clust_co, na.rm = TRUE), stats::sd(rs$modularity)),
    stringsAsFactors = FALSE)
  cmp$Z_score <- round((cmp$Real - cmp$Rand_Mean) / cmp$Rand_SD, 3)
  # ⚠ 原脚本把 NaN（Rand_SD = 0 的确定性零模型，如 fixed n/m 的 Avg Degree 与
  #   Density）经 case_when 落进 "Not significant"；这里保留同一分类结果，
  #   但把 NA 显式留成 NA 并在 Note sheet 说明"检验退化"，不再把 NA 说成不显著。
  cmp$Interpretation <- ifelse(is.na(cmp$Z_score), NA_character_,
                        ifelse(cmp$Z_score > 1.96, "Higher than random",
                        ifelse(cmp$Z_score < -1.96, "Lower than random",
                               "Not significant")))
  list(stats = rs, comparison = cmp)
}

# -----------------------------------------------------------------------------
# 出图用：ggplot 版 igraph 网络面板（样式走 MD_NET）
# -----------------------------------------------------------------------------

#' igraph -> ggplot 网络面板（按门着色、正负边着色）
#' @param layout 已算好的布局矩阵；NULL 时按 st$layout 现算（调用方必须先 md_seed）
#' @param legend 是否在面板内显示 Phylum 图例。默认 FALSE —— 网络的 Phylum /
#'   Correlation 图例由专门的 local80_net_legend() 面板承载，避免 8 张分组网络
#'   各画一份重复图例（原图 Fig5c 的做法也是把图例单独放在第 3 行）。
#' @param guides_fill 是否保留 fill 图例（"keep" / "none"）。分组网络面板用
#'   "none"：它们与图例带同名同义，收集时会合成一坨重复图例。global 面板用
#'   "keep"：它含全部 6 个门，图例带上缺的门（如 Actinobacteria）由它补全。
local80_net_gg <- function(CFG, g, st, phylum_of, phylo_colors, lay = NULL,
                           title = NULL, subtitle = NULL, tag = NULL,
                           legend = FALSE, guides_fill = "keep") {
  if (is.null(g) || igraph::ecount(g) == 0) return(md_empty_plot(CFG, "No valid network"))
  vn <- igraph::V(g)$name
  if (is.null(lay)) {
    lay <- switch(st$layout %||% "graphopt",
                  circle  = igraph::layout_in_circle(g),
                  nicely  = igraph::layout_nicely(g),
                  fr      = igraph::layout_with_fr(g),
                  igraph::layout_with_graphopt(g))
    lay <- igraph::norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  }
  lay <- as.matrix(lay)
  el  <- igraph::as_edgelist(g, names = FALSE)
  deg <- as.numeric(igraph::degree(g, normalized = TRUE))
  nd  <- data.frame(x = lay[, 1], y = lay[, 2], ASV_ID = vn,
                    degree = deg, stringsAsFactors = FALSE)
  nd$phylum <- phylum_of(nd$ASV_ID)
  ed <- data.frame(x = lay[el[, 1], 1], y = lay[el[, 1], 2],
                   xend = lay[el[, 2], 1], yend = lay[el[, 2], 2],
                   sign = igraph::E(g)$sign, stringsAsFactors = FALSE)
  ed$edge_col <- ifelse(ed$sign == "pos", st$edge_col_pos, st$edge_col_neg)

  lv <- names(phylo_colors)
  nd$phylum <- ifelse(nd$phylum %in% lv, nd$phylum, "Unknown")
  if (!"Unknown" %in% lv) {
    phylo_colors <- c(phylo_colors, Unknown = "#BDBDBD")
  }
  g <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed,
                          ggplot2::aes(x = .data$x, y = .data$y,
                                       xend = .data$xend, yend = .data$yend),
                          colour = ed$edge_col, linewidth = st$edge_width,
                          alpha = st$edge_alpha) +
    ggplot2::geom_point(data = nd,
                        ggplot2::aes(x = .data$x, y = .data$y,
                                     fill = .data$phylum,
                                     size = .data$degree),
                        shape = 21, colour = "black", stroke = 0.15,
                        alpha = st$node_alpha) +
    ggplot2::scale_fill_manual(values = phylo_colors, drop = FALSE, name = "Phylum") +
    ggplot2::scale_size_continuous(range = c(st$node_size, st$node_size_hub),
                                   guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, tag = tag, x = NULL, y = NULL) +
    ggplot2::coord_equal() +
    theme_md(CFG = CFG) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = if (isTRUE(legend)) "bottom" else "none",
      legend.key.size = grid::unit(0.3, "cm"),
      plot.title = ggplot2::element_text(size = md_pt(CFG, TRUE)$title,
                                         face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = md_pt(CFG, TRUE)$ann,
                                            colour = "grey40", hjust = 0.5, vjust = 1.6),
      plot.margin = ggplot2::margin(2, 2, 2, 2))
  if (identical(guides_fill, "none")) g <- g + ggplot2::guides(fill = "none")
  g
}

#' 网络图例（门色 + 正负边）单独一张 ggplot，供拼图时占一整条
#'
#' ⚠ 2026-09-26（第六轮）重做：原版只画了一排彩色圆点、**没有写门名**，
#'   读者只能看到 6 个色块却不知道哪个颜色是哪个门（用户报的"Phylum 没有注释"）。
#'   现在每个圆点正下方写门名；门多时按每行 8 个折行。下方单独一行讲相关方向。
#'   本函数用 `theme_void()`，调用方必须用 `local80_compose(..., theme_comp = FALSE)`，
#'   否则 `theme_md_comp()` 会把坐标轴重新打开。
local80_net_legend <- function(CFG, phylo_colors, pt) {
  lv <- names(phylo_colors)
  n  <- length(lv)
  if (n == 0) return(md_empty_plot(CFG, "No phylum legend"))
  per <- 8L
  nrow_p <- ceiling(n / per)
  idx  <- seq_len(n)
  pos <- data.frame(
    col   = ((idx - 1L) %% per) + 1L,
    row   = ((idx - 1L) %/% per),
    phylum = factor(lv, levels = lv),
    lab    = lv,
    stringsAsFactors = FALSE)
  pos$x  <- pos$col
  pos$y  <- -pos$row          # 圆点所在行
  pos$yl <- pos$y - 0.40      # 门名写在该行圆点正下方
  y_corr <- -nrow_p - 0.95    # 相关方向那一行
  ed <- data.frame(x = c(1.3, 5.3), y = y_corr,
                   lab = c("Positive", "Negative"),
                   col = c(MD_NET$edge_col_pos, MD_NET$edge_col_neg))
  ggplot2::ggplot() +
    ggplot2::geom_point(data = pos,
                        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$phylum),
                        shape = 21, size = 2.1, colour = "black", stroke = 0.2) +
    ggplot2::geom_text(data = pos,
                       ggplot2::aes(x = .data$x, y = .data$yl, label = .data$lab),
                       size = pt2gg(pt$ann), colour = "grey15") +
    ggplot2::geom_text(data = data.frame(x = 1, y = 0.75, lab = "Phylum"),
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                       size = pt2gg(pt$ann), fontface = "bold", hjust = 0.6) +
    ggplot2::geom_text(data = data.frame(x = 1, y = y_corr, lab = "Correlation"),
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                       size = pt2gg(pt$ann), fontface = "bold", hjust = 0.6) +
    ggplot2::geom_segment(data = ed,
                          ggplot2::aes(x = .data$x - 0.45, y = .data$y,
                                       xend = .data$x + 0.45, yend = .data$y),
                          colour = ed$col, linewidth = 1.2) +
    ggplot2::geom_text(data = ed,
                       ggplot2::aes(x = .data$x + 0.6, y = .data$y, label = .data$lab),
                       size = pt2gg(pt$ann), colour = "grey15", hjust = 0) +
    ggplot2::scale_fill_manual(values = phylo_colors, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0.2, per + 1.9), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(y_corr - 0.6, 1.15), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void()
}

#' 径向展开 + 归一化（原脚本 layout 后处理，逐行照抄）
local80_radial <- function(l) {
  l <- igraph::norm_coords(l, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  rc <- colMeans(l)
  lc <- sweep(l, 2, rc)
  r  <- sqrt(rowSums(lc^2))
  rmax <- max(r)
  if (!is.finite(rmax) || rmax <= 0) return(l)
  rnew <- rmax * (r / rmax)^0.5
  sc <- rnew / r
  sc[!is.finite(sc)] <- 0
  igraph::norm_coords(lc * sc, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
}

#' 图 B 用的 fr 布局（原脚本：niter 5000、start.temp = sqrt(n) * 2）
local80_layout_fr <- function(g) {
  local80_radial(igraph::layout_with_fr(g, niter = 5000,
                                        start.temp = sqrt(igraph::vcount(g)) * 2))
}

# -----------------------------------------------------------------------------
# 度量箱线图（源：5_network_structure.R:1481-1597 的 fig5b_* 三函数）
# -----------------------------------------------------------------------------

#' 每个度量在某个因子上的统计：主检验 P、rank-based R2、CLD 字母与字母基线
#' 统计链与全包一致：md_group_test()（KW + Dunn(BH) + multcompView 紧凑字母）；
#' rank-based R2 口径照抄原脚本（1 - SS_within / SS_total，秩空间）。
local80_metric_stats <- function(CFG, df, grp, metric_levels, labels) {
  rows <- lapply(metric_levels, function(ml) {
    d <- df[!is.na(df$Value) & as.character(df$MetricLabel) == ml, , drop = FALSE]
    if (nrow(d) < 4) {
      return(data.frame(MetricLabel = ml, p = NA_real_, R2 = NA_real_,
                        Group = character(0), Letter = character(0),
                        letter_y = numeric(0), stringsAsFactors = FALSE))
    }
    gt <- tryCatch(md_group_test(CFG, d, "Value", grp, test = "kruskal_dunn",
                                 p_adjust = "BH", min_n = 2,
                                 order_levels = unique(as.character(d[[grp]]))),
                   error = function(e) NULL)
    if (is.null(gt) || is.null(gt$desc) || nrow(gt$desc) == 0) {
      return(data.frame(MetricLabel = ml, p = NA_real_, R2 = NA_real_,
                        Group = character(0), Letter = character(0),
                        letter_y = numeric(0), stringsAsFactors = FALSE))
    }
    p_kw <- if (!is.null(gt$test) && "p" %in% names(gt$test))
      suppressWarnings(min(gt$test$p, na.rm = TRUE)) else NA_real_
    if (!is.finite(p_kw)) p_kw <- NA_real_
    r  <- rank(d$Value)
    ss_t <- sum((r - mean(r))^2)
    ss_w <- sum(vapply(split(r, as.character(d[[grp]])),
                       function(x) sum((x - mean(x))^2), numeric(1)))
    rng <- diff(range(d$Value, na.rm = TRUE))
    off <- if (is.finite(rng) && rng > 0) 0.05 * rng else 0
    lett <- if (!is.null(gt$letters)) gt$letters else NULL
    data.frame(
      MetricLabel = ml,
      p    = p_kw,
      R2   = if (ss_t > 0) round(1 - ss_w / ss_t, 3) else NA_real_,
      Group    = gt$desc$group,
      Letter   = if (is.null(lett)) "" else
        as.character(lett$letter[match(gt$desc$group, lett$group)]),
      letter_y = gt$desc$max + off,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(NULL)
  out$Letter[is.na(out$Letter)] <- ""
  # ⚠ 显示名可能重复（Density_sample 与 SubDensity 都叫 "Density"），
  #   factor 的 levels 必须唯一，否则报 "factor level [12] is duplicated"。
  out$MetricLabel <- factor(out$MetricLabel, levels = unique(unname(labels)))
  out
}

#' 每个度量一行的面板注释："R2 = 0.21, p < 0.001 ***"
local80_metric_lab <- function(stats)
  data.frame(MetricLabel = stats$MetricLabel,
             label = sprintf("R2 = %s, %s%s",
                             ifelse(is.na(stats$R2), "NA",
                                    formatC(stats$R2, format = "f", digits = 2)),
                             ifelse(is.na(stats$p), "p = NA",
                                    ifelse(stats$p < 0.001, "p < 0.001",
                                           paste0("p = ", format(stats$p, digits = 3)))),
                             ifelse(is.na(stats$p), "", paste0(" ", md_sig_label(stats$p)))),
             stringsAsFactors = FALSE)

#' 一个因子一格的度量箱线面板（照抄 make_fig5b_panel 的版式）
local80_fig5b_panel <- function(CFG, df, grp, stats, metric_levels, pal,
                                title = NULL, ncol = 3, pt = NULL) {
  pt <- pt %||% md_pt(CFG)
  lab <- local80_metric_lab(stats)
  stats$MetricLabel <- factor(stats$MetricLabel, levels = metric_levels)
  lab$MetricLabel   <- factor(lab$MetricLabel,   levels = metric_levels)
  d <- df
  d$.g <- as.character(d[[grp]])
  ggplot2::ggplot(d, ggplot2::aes(x = .data$.g, y = .data$Value)) +
    ggplot2::geom_jitter(ggplot2::aes(colour = .data$.g), size = MD_PT_SIZE$jitter,
                         alpha = MD_PT_SIZE$alpha,
                         position = ggplot2::position_jitter(width = MD_PT_SIZE$jitter_width,
                                                             seed = MD_JITTER_SEED),
                         stroke = MD_PT_SIZE$stroke, show.legend = FALSE) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = .data$.g), fill = NA, width = 0.55,
                          linewidth = MD_LW$box, outlier.shape = NA) +
    ggplot2::geom_text(data = stats,
                       ggplot2::aes(x = .data$Group, y = .data$letter_y,
                                    label = .data$Letter),
                       inherit.aes = FALSE, vjust = 0, size = pt2gg(pt$ann),
                       fontface = "bold", colour = "black") +
    ggplot2::geom_text(data = lab,
                       ggplot2::aes(x = Inf, y = Inf, label = .data$label),
                       inherit.aes = FALSE, hjust = 1.05, vjust = 1.5,
                       size = pt2gg(pt$ann), colour = "grey30") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30))) +
    ggplot2::facet_wrap(~ MetricLabel, scales = "free_y", ncol = ncol) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    theme_md(CFG = CFG) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = pt$strip, face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = pt$tick),
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold",
                                                      size = pt$title))
}

#' 一张度量图 = 每个因子一横行（与 make_fig5b_figure 同版式）
#' 因子只有 1 个时只出 1 行；2 个时上下两行。
local80_metrics_figure <- function(CFG, df, factors, labels, metric_levels,
                                   pt_comp = TRUE) {
  pt <- md_pt(CFG, comp = pt_comp)
  plist <- list()
  for (f in factors) {
    st <- local80_metric_stats(CFG, df, f, metric_levels, labels)
    if (is.null(st) || nrow(st) == 0) next
    cols <- local80_group_cols(CFG, unique(as.character(st$Group)), f)
    # ⚠ 2026-09-26（第六轮）修 bug：原来这里传的是 `unname(labels[metric_levels])` ——
    #   而 `metric_levels` 传进来时**已经是 display 名**（local80_fig5b_panel 的调用方
    #   已经查过一次表），二次查表得到 `c(NA,...,"Robustness")`，于是 stats / lab 两个
    #   **文本层**的 MetricLabel 被压成 1 个水平，而数据层仍是 6 个水平 ——
    #   `facet_wrap` 取并集后多出一个"只有字母与 R2 注释、没有箱体"的幽灵面板（实测 7 格）。
    #   直接把 metric_levels 透传即可，现在是 **6 + 6 = 12 格**。
    plist[[f]] <- local80_fig5b_panel(
      CFG, df, f, st, metric_levels, cols,
      title = paste0("Sample-level Network Metrics by ", f),
      ncol = 3, pt = pt)
  }
  if (length(plist) == 0) return(NULL)
  if (length(plist) == 1) return(plist[[1]])
  Reduce(`/`, plist)
}

# =============================================================================
# 模块入口
# =============================================================================
md_network <- function(CFG) {
  module <- "80"
  t_start <- Sys.time()
  notes <- character(0)
  add_note <- function(...) notes <<- c(notes, paste0(...))

  if (!isTRUE(CFG$network$enabled %||% TRUE))
    return(md_result(module, "SKIP", "network$enabled = FALSE"))

  # SpiecEasi（network_alt 组）是 GitHub 包，缺了必须**回退**而不是让模块 SKIP
  method <- as.character(CFG$network$method %||% "cor")[1]
  method_note <- NULL
  if (identical(method, "sparcc")) {
    ok_alt <- tryCatch(md_need(CFG, "network_alt", module), error = function(e) FALSE)
    if (!isTRUE(ok_alt)) {
      md_warn(CFG, "network$method = 'sparcc' 但 SpiecEasi 不可用 —— 回退到 'cor'（Spearman 相关近似）")
      method <- "cor"
      method_note <- "network$method = 'sparcc' 但 SpiecEasi 不可用，已回退到 'cor'"
      add_note(method_note)
    }
  }

  if (!md_need(CFG, c("network", "ecology"), module))
    return(md_result(module, "SKIP", "缺少依赖包（igraph / tidygraph / vegan）"))

  core <- tryCatch(md_core(CFG), error = function(e) NULL)
  if (is.null(core)) {
    md_skip(CFG, module, "取不到 cache/05_core.rds（请先跑模块 05）")
    return(md_result(module, "SKIP", "缺 05 产物"))
  }

  md_load(c("igraph", "tidyr"))
  if (requireNamespace("tidygraph", quietly = TRUE))
    suppressPackageStartupMessages(library(tidygraph, quietly = TRUE,
                                           warn.conflicts = FALSE))

  # ---- 参数 ----
  factors    <- as.character(CFG$design$factors %||% character(0))
  min_cor    <- as.numeric(CFG$network$min_cor %||% 0.5)
  sub_min_cor <- as.numeric(CFG$network$sub_min_cor %||% min_cor)
  min_prev   <- as.numeric(CFG$network$min_prev %||% 0.02)
  max_asv    <- CFG$network$max_asv
  algo       <- as.character(CFG$network$module_algo %||% "louvain")[1]
  zi_thr     <- as.numeric(CFG$network$zi_thr %||% MD_NET$zi_thr)
  pi_thr     <- as.numeric(CFG$network$pi_thr %||% MD_NET$pi_thr)
  random_n   <- as.integer(CFG$network$random_n %||% 100L)
  nsim       <- as.integer(CFG$network$robustness_nsim %||% 30L)
  n_boot     <- as.integer(CFG$network$n_boot %||% 1000L)
  sub_by     <- as.character(CFG$network$sub_network_by %||% "abundance_rank")[1]
  turnover_by <- if (is.null(CFG$network$turnover_by)) NULL else
    as.character(CFG$network$turnover_by)[1]
  gf_cfg     <- CFG$network$group_factors
  group_factors <- if (is.null(gf_cfg)) factors else intersect(as.character(gf_cfg), factors)
  if (length(group_factors) == 0) group_factors <- utils::head(factors, 2)
  metrics_main <- as.character(CFG$network$metrics_main %||%
                                 c("MeanNodeDegree", "SumNodeDegree", "MeanBetweenness",
                                   "MeanEigen", "MeanWithinCorr", "Cohesion_pos"))
  metrics_supp <- as.character(CFG$network$metrics_supp %||%
                                 c("MeanCloseness", "Density_sample", "Clustering_sample",
                                   "Connectedness_sample", "Modularity_sample",
                                   "Cohesion_neg"))
  supp_fig   <- isTRUE(CFG$network$metrics_supp_figure %||% TRUE)
  a_thr <- CFG$composition$abund_t %||% CFG$niche$thr_abund %||% 0.1
  p_thr <- CFG$composition$freq_t  %||% CFG$niche$thr_prev  %||% 10
  labels <- local80_metric_labels()

  # ---- 分析段（全部计算；图与表在段外）----
  md_stage(CFG, "80_analysis", {
    otu_raw <- core$asv_raw
    meta    <- as.data.frame(core$meta, stringsAsFactors = FALSE)
    smp     <- CFG$data$sample_col
    tax     <- as.data.frame(core$tax, stringsAsFactors = FALSE)
    id_col  <- names(tax)[1]
    if (!id_col %in% names(tax)) id_col <- CFG$data$id_col
    if (!"ASV_ID" %in% names(tax)) tax$ASV_ID <- tax[[id_col]]

    meta <- meta[match(colnames(otu_raw), as.character(meta[[smp]])), , drop = FALSE]
    otu_rel <- md_rel(otu_raw) * 100
    n_sample <- ncol(otu_rel)

    # 1) 全局网络（照抄 L117-182：出现率过滤 -> 丰度上限截断走 CFG$network$max_asv）
    prev  <- rowSums(otu_rel > 0.01) / n_sample
    keep0 <- which(prev >= min_prev)
    if (!is.null(max_asv) && length(keep0) > as.integer(max_asv)) {
      ord <- order(rowMeans(otu_rel[keep0, , drop = FALSE]), decreasing = TRUE)
      keep0 <- keep0[ord[seq_len(as.integer(max_asv))]]
    }
    otu_filt <- otu_rel[keep0, , drop = FALSE]
    if (!is.matrix(otu_filt)) otu_filt <- as.matrix(otu_filt)

    net_built <- FALSE
    if (nrow(otu_filt) >= 3L) {
      res <- tryCatch(local80_calc_network(otu_filt, min_cor, min_prev),
                      error = function(e) NULL)
      if (!is.null(res)) {
        net_built <- TRUE
        g_global  <- res$g
        cor_mat   <- res$cor_mat
        kept_asv  <- res$kept_asv
      }
    }
    if (!net_built) {
      g_global <- NULL; cor_mat <- NULL; kept_asv <- character(0)
      skip_reason <- sprintf("过滤后 ASV = %d（阈值 min_prev = %g），或 |r| >= %g 阈值下无任何边",
                             nrow(otu_filt), min_prev, min_cor)
    } else {
      skip_reason <- NULL
    }

    # 1a) 门的配色（照抄 L184-198：npg nrc <= 10 色，> 10 时后半段 viridis）
    tax_filt <- tax[tax$ASV_ID %in% rownames(otu_filt), , drop = FALSE]
    tax_filt$Phylum <- as.character(tax_filt$Phylum)
    tax_filt$Phylum[is.na(tax_filt$Phylum) | !nzchar(tax_filt$Phylum)] <- "Unknown"
    tax_filt$Genus <- as.character(tax_filt$Genus)
    all_phyla <- sort(unique(tax_filt$Phylum))
    n_phyla   <- length(all_phyla)
    if (n_phyla > 10 && !requireNamespace("viridis", quietly = TRUE)) {
      md_need(CFG, "viridis", "80")
    }
    phylo_colors <- local80_phylum_colors(all_phyla)
    phylum_of <- function(ids) {
      v <- tax_filt$Phylum[match(ids, tax_filt$ASV_ID)]
      v[is.na(v) | !nzchar(v)] <- "Unknown"
      v
    }

    if (net_built) {
      # 2) 节点中心性 + 度数
      cent_df <- local80_centrality(g_global)
      deg_raw <- igraph::degree(g_global, normalized = FALSE)
      names(deg_raw) <- igraph::V(g_global)$name
      e_df    <- igraph::as_edgelist(g_global)
      edge_wt <- igraph::E(g_global)$weight
      edge_sign <- igraph::E(g_global)$sign

      # 3) 样本级指标
      # ⚠ 本行是相对原文的**唯一 bug 修复**：原脚本的 set.seed 在 robustness_R()
      #   的 30 次模拟之后（L401），导致 Fig5b 每次都不一样。这里提前播种。
      md_seed(CFG, md_offset(CFG, "network_robustness", 8220L))
      sample_metrics <- local80_sample_metrics(
        g_global, otu_filt, meta, smp, factors,
        robustness_nsim = nsim, group_net = TRUE)
      long_cols <- intersect(c(
        "MeanNodeDegree", "SumNodeDegree", "MeanBetweenness", "MeanCloseness",
        "MeanEigen", "MeanWithinCorr", "N_DetectedASV",
        "SubDensity", "SubClustering", "SubConnectedness", "SubAvgPath",
        "SubHeterogeneity", "SubModularity", "SubNModules",
        "Cohesion_pos", "Cohesion_neg", "Robustness",
        "Density_sample", "Clustering_sample", "Connectedness_sample",
        "Modularity_sample"), names(sample_metrics))
      sample_metrics_long <- tidyr::pivot_longer(
        sample_metrics, cols = tidyselect::all_of(long_cols),
        names_to = "Metric", values_to = "Value")
      sample_metrics_long$MetricLabel <- factor(local80_lab(sample_metrics_long$Metric, labels),
                                                levels = unique(unname(labels)))

      # 4) 组间检验（KW + BH）与 PERMANOVA
      metrics_vec <- unique(c(metrics_main, metrics_supp))
      kw_results <- local80_kw(CFG, sample_metrics_long, metrics_vec, factors)
      permanova  <- local80_permanova(CFG, sample_metrics, metrics_main, factors)

      # 5) 分组网络 + 拓扑
      gnet <- local80_group_networks(CFG, otu_filt, meta, group_factors,
                                     min_cor, min_prev)
      group_topo  <- gnet$topo
      group_nets  <- gnet$nets
      group_ns    <- gnet$ns
      group_labels <- gnet$labels

      # 6) 模块检测 + Zi-Pi + 全局拓扑 + Hub
      md_seed(CFG, md_offset(CFG, "network_module", 8210L))
      comm_global <- local80_communities(g_global, algo)
      igraph::V(g_global)$module <- as.character(igraph::membership(comm_global))
      modularity_global <- igraph::modularity(comm_global)
      n_modules_global  <- length(unique(igraph::V(g_global)$module))
      module_size <- as.data.frame(table(Module = igraph::V(g_global)$module),
                                   stringsAsFactors = FALSE)
      names(module_size) <- c("Module", "N_Nodes")

      zi_pi_df <- local80_zi_pi(g_global, zi_thr, pi_thr)
      zi_pi_df <- merge(zi_pi_df,
                        tax_filt[, c("ASV_ID", "Phylum", "Genus")],
                        by = "ASV_ID", all.x = TRUE, sort = FALSE)
      zi_pi_df <- zi_pi_df[match(igraph::V(g_global)$name, zi_pi_df$ASV_ID), , drop = FALSE]
      zi_pi_df$Phylum[is.na(zi_pi_df$Phylum) | !nzchar(zi_pi_df$Phylum)] <- "Unknown"
      role_levels <- names(MD_NET$role_colors)
      zi_pi_df$role <- factor(zi_pi_df$role, levels = role_levels)
      hub_df <- zi_pi_df[as.character(zi_pi_df$role) %in%
                           c("Network Hub", "Module Hub", "Connector"), , drop = FALSE]
      hub_df <- hub_df[order(as.character(hub_df$role), -hub_df$Zi), , drop = FALSE]

      # 节点表（跨模块契约）：中心性 + 模块 + Zi-Pi + 分类学
      nodes_tbl <- local80_nodes_table(g_global, cent_df, zi_pi_df, tax_filt, deg_raw)
      edges_tbl <- data.frame(from = e_df[, 1], to = e_df[, 2],
                              weight = as.numeric(edge_wt),
                              sign = as.character(edge_sign),
                              stringsAsFactors = FALSE)

      topo_summary <- local80_global_topo(g_global, modularity_global,
                                          n_modules_global)

      # 7) 随机网络比对
      md_seed(CFG, md_offset(CFG, "network_random", 8250L))
      rnd <- local80_random(g_global, random_n)
      random_stats <- rnd$stats
      comparison_df <- rnd$comparison

      # 8) 丰富/稀有亚网络 + 节点周转
      rel_for_cells <- otu_raw
      cells <- local80_cells(md_rel(rel_for_cells) * 100, a_thr, p_thr)
      cells$Cell <- factor(cells$Biosphere, levels = names(MD_NET$cell_labels))
      abr <- local80_abundant_rare(CFG, otu_filt, meta, smp, cells, sub_by,
                                   sub_min_cor, min_prev)
      turnover <- local80_turnover(CFG, otu_filt, meta, smp, cells, turnover_by,
                                   group_factors, sub_min_cor, min_prev)

      # 9) bootstrap 均衡（照抄 L391-498：分组内抽到最小样本量，报均值 + 95% CI）
      md_seed(CFG, md_offset(CFG, "network_bootstrap", 8260L))
      boot_summary <- local80_bootstrap(sample_metrics, factors, metrics_main, n_boot)
    }
  })

  # ---- 网络建不起来：优雅降级 ----
  if (!isTRUE(net_built)) {
    md_skip(CFG, module, skip_reason)
    note <- data.frame(
      item = c("Status", "Reason", "Filter", "Parameters"),
      value = c("SKIP - global co-occurrence network could not be built",
                skip_reason,
                sprintf("ASV x sample = %d x %d after occurrence filter (min_prev = %g)",
                        nrow(otu_filt), ncol(otu_filt), min_prev),
                sprintf("method = %s, cor = %s, min_cor = %g",
                        method, CFG$network$cor_method %||% "spearman", min_cor)),
      stringsAsFactors = FALSE)
    md_write_xlsx(CFG, "80_network", list(Note = note))
    return(md_result(module, "SKIP", skip_reason))
  }

  st <- md_net_cfg(CFG)
  pt <- md_pt(CFG)
  ptc <- md_pt(CFG, comp = TRUE)

  # ---- 跨模块契约：cache/80_network.rds（模块 81 读网络角色）----
  cache_file <- md_path(CFG, "cache", "80_network.rds")
  tryCatch({
    params <- list(
      method = method, cor_method = CFG$network$cor_method %||% "spearman",
      min_cor = min_cor, sub_min_cor = sub_min_cor, min_prev = min_prev,
      max_asv = max_asv, module_algo = algo, zi_thr = zi_thr, pi_thr = pi_thr,
      n_nodes = igraph::vcount(g_global), n_edges = igraph::ecount(g_global),
      n_modules = n_modules_global, modularity = modularity_global,
      sample_col = CFG$data$sample_col, id_col = CFG$data$id_col,
      factors = CFG$design$factors, seed = CFG$opts$seed)
    saveRDS(list(nodes = nodes_tbl, edges = edges_tbl, params = params), cache_file)
    md_log(CFG, "OK", "写出 ", basename(cache_file), "（",
           nrow(nodes_tbl), " 节点 / ", nrow(edges_tbl), " 边）")
  }, error = function(e) {
    md_warn(CFG, "写 cache/80_network.rds 失败：", conditionMessage(e))
    add_note(paste0("cache/80_network.rds 写出失败：", conditionMessage(e)))
  })

  figures <- character(0)

  # ===========================================================================
  # 图：80a - 80f（2026-09-26 第六轮按用户裁定重排）
  #   80a  全局网络 + 观测 vs 随机
  #   80b  分组子网络 2x4（+ 图例带）
  #   80c  节点周转 Sankey + 丰富/稀有亚网络
  #   80d  节点属性（betweenness | 度-丰度 | Zi-Pi）
  #   80e  指标主图（6 指标 x 因子）
  #   80f  指标附图（另外 6 指标）
  #   注：`md_module_of("80a_...")` 仍返回 "80"，figure$skip_modules 与 95 的归属不受影响。
  # ===========================================================================
  # ⚠ 布局是随机过程 —— 在**出图段**先固定种子（不在 md_stage 里）
  md_seed(CFG, md_offset(CFG, "network_layout", 8230L))

  # ggalluvial 的 stat_stratum / stat_alluvium 只在包内注册 —— 必须 library() 进
  # 搜索路径，geom_text(stat = "stratum") 才找得到（否则报 Can't find stat "stratum"）。
  if (requireNamespace("ggalluvial", quietly = TRUE))
    suppressPackageStartupMessages(library(ggalluvial, quietly = TRUE,
                                           warn.conflicts = FALSE))

  pT_A   <- local80_panel_global(CFG, g_global, st, phylum_of, phylo_colors, pt)
  pT_B   <- local80_panel_groups(CFG, group_nets, group_labels, group_ns, st,
                                 phylum_of, phylo_colors)
  pT_Bleg <- local80_net_legend(CFG, phylo_colors, ptc)
  pT_C   <- local80_panel_random(CFG, comparison_df, random_n, ptc)
  pT_D   <- local80_panel_sankey(CFG, turnover, ptc)
  pT_E   <- local80_panel_abundant_rare(CFG, abr, ptc)

  # ---- 80a：全局网络 + vs 随机 ----
  #   `theme_comp = FALSE`：**不能**给整张图套 theme_md_comp() —— 它会把 A 面板
  #   （网络快照，已 element_blank() 掉坐标轴）的 x/y 轴重新打开（用户报的第 1 条）。
  #   需要缩字号的只有 B 面板，单独加。
  pT_C <- pT_C + theme_md_comp(CFG)
  # 左右各占一半版幅（"AABB"）。画幅配成 183 x 92 mm，使每个面板接近方形 ——
  # 网络是 coord_equal 的，格子太宽只会在左右留白。
  f_a <- md_save_figure(CFG, local80_compose(CFG, list(A = pT_A, B = pT_C),
                                             design = "AABB",
                                             guides = "collect", legend = "bottom",
                                             theme_comp = FALSE),
                        "80a_network_global", width = 183, height = 92)
  if (!is.null(f_a)) figures <- c(figures, "figures/80a_network_global.pdf")

  # ---- 80b：分组子网络（6 个 -> 2 行 x 3 列）+ 图例带 ----
  #   pT_B 本身是若干个网络的 patchwork，必须先用 wrap_elements() 把它当**一个**面板，
  #   否则 design 的格子不够、patchwork 会报 "Too few patch areas to hold all plots" 并丢图。
  #   同样 `theme_comp = FALSE`：这些网络图与图例带都不该有坐标轴。
  #   ⚠ `heights = c(7, 1)` 是关键：design 的行**默认等分**，不加权重的话
  #     一条又扁又宽的图例带会和整个 2x3 网络区各吃一半高度（用户报的"图例版幅过大"）。
  f_b <- md_save_figure(CFG,
                        local80_compose(CFG,
                                        list(A = patchwork::wrap_elements(pT_B), L = pT_Bleg),
                                        design = "AAAA\nLLLL",
                                        guides = "collect", legend = "bottom",
                                        theme_comp = FALSE,
                                        heights = c(7, 1)),
                        "80b_network_by_group", width = 183, height = 165)
  if (!is.null(f_b)) figures <- c(figures, "figures/80b_network_by_group.pdf")

  # ---- 80c：节点周转 Sankey + 丰富/稀有亚网络 ----
  f_c <- md_save_figure(CFG, md_compose(CFG, list(A = pT_D, B = pT_E),
                                        design = "AABB", guides = "collect",
                                        legend = "bottom"),
                        "80c_network_turnover", width = 183, height = 150)
  if (!is.null(f_c)) figures <- c(figures, "figures/80c_network_turnover.pdf")

  # ===========================================================================
  # 图 80d：节点属性（原文 Fig5a：A betweenness | B 度-丰度 | C Zi-Pi）
  #   ⚠ 2026-09-26（第六轮）：原实现按原文的 533.4 x 177.8 mm 单排三格出图，
  #     在屏幕上按宽度缩放显示时字与点比其余 183 mm 的图明显小。用户要求"统一配置"，
  #     故改为 **183 mm 宽、三格纵向排列**，字号/点尺寸与全包其余模块完全一致。
  #     想回到原文的宽幅单排：把 CFG$network$style$wide_layout 设为 TRUE。
  # ===========================================================================
  pA <- local80_plot_betweenness(CFG, nodes_tbl, phylo_colors, st, pt)
  pB <- local80_plot_degree_abund(CFG, otu_filt, kept_asv, g_global, tax_filt,
                                 phylo_colors, pt)
  pC <- local80_plot_zipi(CFG, zi_pi_df, pt, zi_thr, pi_thr)
  if (isTRUE(st$wide_layout)) {
    p_nodes <- md_compose(CFG, list(A = pA, B = pB, C = pC), design = "ABC",
                          guides = "collect", legend = "bottom")
    f_d <- md_save_figure(CFG, p_nodes, "80d_network_nodes",
                          width = 533.4, height = 177.8)
  } else {
    p_nodes <- md_compose(CFG, list(A = pA, B = pB, C = pC), design = "A\nB\nC",
                          guides = "collect", legend = "bottom")
    f_d <- md_save_figure(CFG, p_nodes, "80d_network_nodes",
                          width = 183, height = 220)
  }
  if (!is.null(f_d)) figures <- c(figures, "figures/80d_network_nodes.pdf")

  # ===========================================================================
  # 图 2 / 3：80_network_metrics(.pdf / _supp.pdf)
  # ===========================================================================
  mf <- utils::head(factors, 2)
  df_main <- sample_metrics_long[as.character(sample_metrics_long$Metric) %in% metrics_main, ,
                                 drop = FALSE]
  ok_main <- local80_ok_metrics(sample_metrics, metrics_main)
  add_note("Sample_metrics 表里 Density_sample / Clustering_sample / ",
           "Connectedness_sample / Modularity_sample 是**整张全局网络**的拓扑量",
           "（原脚本 sample_metrics 里没有这四列的逐样本版本），",
           "逐样本的等价量请用 SubDensity / SubClustering / SubConnectedness / ",
           "SubModularity 四列")
  if (length(ok_main) == 0) {
    add_note("metrics_main 全部为 NA（无有效样本级指标），跳过 80_network_metrics.pdf")
  } else {
    if (length(ok_main) < length(metrics_main))
      add_note("metrics_main 中全为 NA 的度量已从面板剔除：",
               paste(setdiff(metrics_main, ok_main), collapse = ", "))
    dfm <- df_main[as.character(df_main$Metric) %in% ok_main, , drop = FALSE]
    dfm <- local80_reorder_metric(dfm, ok_main, labels)
    pm <- local80_metrics_figure(CFG, dfm, mf, labels, unname(labels[ok_main]))
    if (!is.null(pm)) {
      f2 <- md_save_figure(CFG, pm, "80e_network_metrics", width = 183, height = 165)
      if (!is.null(f2)) figures <- c(figures, "figures/80e_network_metrics.pdf")
    }
  }

  if (isTRUE(supp_fig)) {
    ok_supp <- local80_ok_metrics(sample_metrics, metrics_supp)
    if (length(ok_supp) == 0) {
      add_note("metrics_supp 全部为 NA，跳过 80_network_metrics_supp.pdf")
    } else {
      if (length(ok_supp) < length(metrics_supp))
        add_note("metrics_supp 中全为 NA 的度量已从面板剔除：",
                 paste(setdiff(metrics_supp, ok_supp), collapse = ", "))
      dfs <- sample_metrics_long[as.character(sample_metrics_long$Metric) %in% ok_supp, ,
                                 drop = FALSE]
      dfs <- local80_reorder_metric(dfs, ok_supp, labels)
      ps <- local80_metrics_figure(CFG, dfs, mf, labels, unname(labels[ok_supp]))
      if (!is.null(ps)) {
        f3 <- md_save_figure(CFG, ps, "80f_network_metrics_supp",
                             width = 183, height = 165)
        if (!is.null(f3)) figures <- c(figures, "figures/80f_network_metrics_supp.pdf")
      }
    }
  } else {
    add_note("CFG$network$metrics_supp_figure = FALSE：未出 80_network_metrics_supp.pdf")
  }
  # 负内聚力在原始 |r| 阈值下常常无负边 —— 全 NA 时把原因写清楚，
  # 不要让它看起来像"算错了"（原脚本 Fig5b 的注释也记录了同一现象）。
  if ("Cohesion_neg" %in% names(sample_metrics) &&
      !any(is.finite(sample_metrics$Cohesion_neg))) {
    add_note("Cohesion_neg 全为 NA：原始 |r| 阈值下该网络没有负相关边")
  }

  # ===========================================================================
  # 图 80a / 80b / 80c 已在上面（图件重排时提前输出）
  #   原实现是一张 80_network_topology.pdf，含 6 个面板（全局网 / 分组网 / 图例带 /
  #   vs 随机 / Sankey / 亚网络）。2026-09-26 第六轮按用户裁定拆成 3 张，
  #   编号统一成 80a-80f。
  # ===========================================================================

  # ===========================================================================
  # 表：results/80_network.xlsx
  # ===========================================================================
  note <- local80_note(CFG, method, method_note, min_cor, sub_min_cor, min_prev,
                       max_asv, algo, zi_thr, pi_thr, random_n, nsim, n_boot,
                       metrics_main, metrics_supp, supp_fig, mf, otu_filt,
                       g_global, modularity_global, n_modules_global,
                       comparison_df, notes)
  sheets <- list(
    Node_centrality = nodes_tbl,
    Sample_metrics  = sample_metrics,
    PERMANOVA       = permanova,
    Group_topology  = group_topo,
    Hub_nodes       = hub_df,
    Global_topology = topo_summary,
    Random_network  = comparison_df,
    Abundant_rare   = abr$topo,
    Node_turnover   = turnover$flows,
    Note            = note)
  # 额外的完整结果（不改变上面 10 个契约 sheet 存在的语义）
  sheets$Kruskal_Wallis      <- kw_results
  sheets$Bootstrap_metrics   <- boot_summary
  sheets$Jaccard_overlap     <- turnover$jaccard
  sheets$Node_persistence    <- turnover$persistence
  sheets$ASV_cells           <- turnover$cells
  sheets$Module_size         <- module_size
  tp <- md_write_xlsx(CFG, "80_network", sheets,
                      meta_extra = list(network_method = method,
                                        min_cor = min_cor, min_prev = min_prev,
                                        nodes = igraph::vcount(g_global),
                                        edges = igraph::ecount(g_global)))

  md_log(CFG, "OK", "模块 80 用时 ",
         format(round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)),
         " s")
  md_result(module, "OK",
            tables  = basename(tp %||% "80_network.xlsx"),
            figures = figures,
            notes   = paste(c(sprintf("nodes=%d edges=%d modules=%d",
                                      igraph::vcount(g_global),
                                      igraph::ecount(g_global), n_modules_global),
                              notes), collapse = " | "))
}

# =============================================================================
# local80_ 辅助函数（分析）
# =============================================================================

#' 门配色（源：5_network_structure.R:184-198）
local80_phylum_colors <- function(all_phyla) {
  n <- length(all_phyla)
  if (n == 0) return(stats::setNames(character(0), character(0)))
  npg <- if (requireNamespace("ggsci", quietly = TRUE))
    ggsci::pal_npg("nrc")(10) else
    unname(MD_PAL_DEFAULT$qualitative)
  pal <- if (n <= 10) npg[seq_len(n)] else {
    tail <- if (requireNamespace("viridis", quietly = TRUE))
      viridis::viridis(n - 10) else
      grDevices::colorRampPalette(MD_PAL_DEFAULT$sequential)(n - 10)
    c(npg, tail)
  }
  stats::setNames(as.character(pal), all_phyla)
}

#' 社区检测（CFG$network$module_algo；louvain 是原文口径）
local80_communities <- function(g, algo = "louvain") {
  w <- abs(igraph::E(g)$weight)
  res <- switch(as.character(algo)[1],
                fast_greedy = tryCatch(igraph::cluster_fast_greedy(g, weights = w),
                                       error = function(e) NULL),
                walktrap   = tryCatch(igraph::cluster_walktrap(g, weights = w),
                                      error = function(e) NULL),
                tryCatch(igraph::cluster_louvain(g, weights = w),
                         error = function(e) NULL))
  if (is.null(res)) {
    res <- tryCatch(igraph::cluster_louvain(g, weights = w), error = function(e) NULL)
  }
  if (is.null(res)) igraph::cluster_fast_greedy(g, weights = w) else res
}

#' 节点表（跨模块契约 nodes 的载体）
local80_nodes_table <- function(g, cent_df, zi_pi_df, tax_filt, deg_raw) {
  vn <- igraph::V(g)$name
  out <- data.frame(
    ASV_ID = vn,
    degree = as.numeric(deg_raw[vn]),
    betweenness = cent_df$Betweenness[match(vn, cent_df$ASV_ID)],
    closeness = cent_df$Closeness[match(vn, cent_df$ASV_ID)],
    eigen = cent_df$Eigen_centrality[match(vn, cent_df$ASV_ID)],
    transitivity = cent_df$Transitivity[match(vn, cent_df$ASV_ID)],
    module = zi_pi_df$module[match(vn, as.character(zi_pi_df$ASV_ID))],
    Zi = zi_pi_df$Zi[match(vn, as.character(zi_pi_df$ASV_ID))],
    Pi = zi_pi_df$Pi[match(vn, as.character(zi_pi_df$ASV_ID))],
    role = as.character(zi_pi_df$role[match(vn, as.character(zi_pi_df$ASV_ID))]),
    Phylum = zi_pi_df$Phylum[match(vn, as.character(zi_pi_df$ASV_ID))],
    Genus = zi_pi_df$Genus[match(vn, as.character(zi_pi_df$ASV_ID))],
    stringsAsFactors = FALSE)
  out$Phylum[is.na(out$Phylum) | !nzchar(out$Phylum)] <- "Unknown"
  out
}

#' 全局拓扑汇总（源：5_network_structure.R:1078-1098）
local80_global_topo <- function(g, modularity_global, n_modules_global) {
  avg_path <- tryCatch(igraph::mean_distance(g, directed = FALSE),
                       error = function(e) NA_real_)
  data.frame(
    Metric = c("Nodes", "Edges", "Positive_Edges", "Negative_Edges",
               "Avg_Degree", "Density", "Avg_Path_Length",
               "Clustering_Coefficient", "Modularity", "Modules"),
    Value = c(igraph::vcount(g), igraph::ecount(g),
              sum(igraph::E(g)$sign == "pos", na.rm = TRUE),
              sum(igraph::E(g)$sign == "neg", na.rm = TRUE),
              round(mean(igraph::degree(g, normalized = FALSE)), 3),
              round(igraph::edge_density(g), 5),
              round(avg_path, 3),
              round(igraph::transitivity(g, type = "global"), 4),
              round(modularity_global, 4), n_modules_global),
    stringsAsFactors = FALSE)
}

#' KW 检验 + BH（源：5_network_structure.R:527-581）
local80_kw <- function(CFG, long_df, metrics_vec, factors) {
  rows <- list()
  for (f in factors) {
    for (m in metrics_vec) {
      sub <- long_df[as.character(long_df$Metric) == m, , drop = FALSE]
      p <- NA_real_; stat <- NA_real_
      ok <- nrow(sub) > 0 && length(unique(stats::na.omit(sub[[f]]))) >= 2 &&
        length(unique(stats::na.omit(sub$Value))) >= 2
      if (ok) {
        kt <- tryCatch(stats::kruskal.test(sub$Value ~ sub[[f]]),
                       error = function(e) NULL)
        if (!is.null(kt)) { stat <- unname(kt$statistic); p <- kt$p.value }
      }
      rows[[length(rows) + 1L]] <- data.frame(
        Metric = m, MetricLabel = local80_lab(m), Test = "Kruskal-Wallis",
        Chi_sq = stat, p_value = p, Group = f, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out$FDR <- md_p_adjust(out$p_value, "BH")
  out$sig <- ifelse(is.na(out$FDR), "",
             ifelse(out$FDR < 0.001, "***", ifelse(out$FDR < 0.01, "**",
             ifelse(out$FDR < 0.05, "*", "ns"))))
  out
}

#' PERMANOVA（源：5_network_structure.R:583-604，欧氏距离 + adonis2 置换）
#' 置换次数走 `CFG$network$permutations` -> `CFG$beta_tests$permutations` -> 999。
local80_permanova <- function(CFG, sample_metrics, metrics_main, factors) {
  n_perm <- as.integer(CFG$network$permutations %||%
                         CFG$beta_tests$permutations %||% 999)
  md_seed(CFG, md_offset(CFG, "network_permanova", 8240L))
  mm <- intersect(metrics_main, names(sample_metrics))
  mm <- mm[vapply(mm, function(z) any(is.finite(sample_metrics[[z]])), logical(1))]
  if (length(mm) < 2 || nrow(sample_metrics) < 4) return(NULL)
  d <- sample_metrics[, mm, drop = FALSE]
  ok <- stats::complete.cases(d)
  if (sum(ok) < 4) return(NULL)
  dm <- tryCatch(vegan::vegdist(d[ok, , drop = FALSE], method = "euclidean"),
                 error = function(e) NULL)
  if (is.null(dm)) return(NULL)
  md_tbl <- sample_metrics[ok, , drop = FALSE]
  out <- list()
  for (f in factors) {
    if (!f %in% names(md_tbl)) next
    if (length(unique(stats::na.omit(md_tbl[[f]]))) < 2) next
    ad <- tryCatch(vegan::adonis2(dm ~ md_tbl[[f]], permutations = n_perm),
                   error = function(e) NULL)
    if (is.null(ad)) next
    df <- as.data.frame(ad)
    df$term <- rownames(df)
    df$Group <- f
    df$Metric <- "All_metrics"
    out[[f]] <- df
  }
  if (length(out) == 0) return(NULL)
  res <- do.call(rbind, out)
  nm <- names(res)
  if ("Pr(>F)" %in% nm) names(res)[nm == "Pr(>F)"] <- "p_value"
  if ("R2" %in% nm)     names(res)[nm == "R2"] <- "R2"
  first <- intersect(c("Group", "Metric", "term", "Df", "SumOfSqs", "R2",
                       "F", "p_value"), names(res))
  cbind(res[, first, drop = FALSE], res[, setdiff(names(res), first), drop = FALSE])
}

#' 分组网络 + 拓扑（源：5_network_structure.R:928-985）
local80_group_networks <- function(CFG, otu_filt, meta, group_factors,
                                   min_cor, min_prev) {
  nets <- list(); topo <- list(); ns <- integer(0); labs <- character(0)
  for (f in group_factors) {
    lv <- local80_levels(CFG, meta, f)
    lv <- lv[vapply(lv, function(x)
      sum(as.character(meta[[f]]) == x, na.rm = TRUE) >= 4, logical(1))]
    for (x in lv) {
      idx <- which(as.character(meta[[f]]) == x)
      sub <- otu_filt[, idx, drop = FALSE]
      g   <- tryCatch(local80_group_net(sub, min_cor, min_prev),
                      error = function(e) NULL)
      key <- paste(f, x, sep = "|")
      nets[[key]] <- g
      topo[[key]] <- local80_group_topology(g, x, f, length(idx))
      ns[key] <- length(idx); labs[key] <- x
    }
  }
  tp <- if (length(topo)) do.call(rbind, topo) else NULL
  if (!is.null(tp)) rownames(tp) <- NULL
  list(nets = nets, topo = tp, ns = ns, labels = labs)
}

#' 丰富/稀有亚网络（源：5_network_structure.R:1952-1990 与 2200-2246）
#' @param sub_by CFG$network$sub_network_by；只实现 "abundance_rank"（2x2 四象限，
#'   原文口径）。其它取值原脚本没有对应实现 —— 直接在表里说明并降级，不猜语义。
local80_abundant_rare <- function(CFG, otu_filt, meta, smp, cells, sub_by,
                                  min_cor, min_prev) {
  bio_levels <- names(MD_NET$cell_labels)
  if (!identical(sub_by, "abundance_rank")) {
    return(list(topo = NULL, grid = NULL, nets = NULL, xvar = NULL,
                note = paste0("sub_network_by = '", sub_by,
                              "' 未实现（原脚本只按 mean abundance x occurrence 的",
                              " 2x2 四象限分类），亚网络面板与表已降级跳过")))
  }
  grp_vec <- stats::setNames(as.character(cells$Cell), cells$ASV_ID)
  g_levels <- bio_levels
  # 横轴：设计因子 1（没有因子时退化为单一水平 "All"）
  xvar <- (CFG$design$factors %||% character(0))[1]
  if (is.null(xvar) || !xvar %in% names(meta)) xvar <- NULL
  xlv <- if (is.null(xvar)) "All" else local80_levels(CFG, meta, xvar)

  rows <- list(); nets <- list()
  for (x in xlv) {
    idx <- if (is.null(xvar)) seq_len(ncol(otu_filt)) else
      which(as.character(meta[[xvar]]) == x)
    if (length(idx) < 4) next
    sub_otu <- otu_filt[, idx, drop = FALSE]
    for (cl in g_levels) {
      asvs <- intersect(names(grp_vec)[!is.na(grp_vec) & grp_vec == cl],
                        rownames(sub_otu))
      if (length(asvs) < 3) {
        rows[[length(rows) + 1L]] <- data.frame(
          Vegetation = x, Cell = cl, Nodes = 0, Edges = 0, Avg_Degree = NA_real_,
          Density = NA_real_, Clustering = NA_real_, Connectedness = NA_real_,
          Modularity = NA_real_, N_Modules = NA_integer_, Heterogeneity = NA_real_,
          stringsAsFactors = FALSE)
        next
      }
      g <- tryCatch(local80_group_net(sub_otu[asvs, , drop = FALSE], min_cor, min_prev),
                    error = function(e) NULL)
      nets[[paste(x, cl, sep = "|")]] <- g
      rows[[length(rows) + 1L]] <- local80_calc_topo(g, x, cl)
    }
  }
  topo <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(topo)) {
    topo$Vegetation <- factor(topo$Vegetation, levels = xlv)
    topo$Cell <- factor(topo$Cell, levels = g_levels)
  }
  grid <- NULL
  if (!is.null(topo) && any(topo$Nodes > 0)) {
    metrics <- c(Nodes = "Network Size", Edges = "Total Links",
                 Avg_Degree = "Average Degree", Clustering = "Clustering",
                 Connectedness = "Connectedness", Modularity = "Modularity")
    long <- do.call(rbind, lapply(names(metrics), function(m) {
      data.frame(Vegetation = topo$Vegetation, Cell = topo$Cell,
                 Metric = m, MetricLabel = metrics[[m]],
                 Value = topo[[m]], stringsAsFactors = FALSE)
    }))
    long <- long[is.finite(long$Value), , drop = FALSE]
    long$MetricLabel <- factor(long$MetricLabel, levels = unname(metrics))
    grid <- long
  }
  list(topo = topo, grid = grid, nets = nets, note = NULL, xvar = xvar,
       gkey = "abundance_rank")
}

#' 亚网络拓扑（源：5_network_structure.R:1931-1950 的 calc_topo）
local80_calc_topo <- function(g, xv, cl) {
  if (is.null(g) || igraph::vcount(g) < 2) {
    return(data.frame(Vegetation = xv, Cell = cl, Nodes = 0, Edges = 0,
                      Avg_Degree = NA_real_, Density = NA_real_,
                      Clustering = NA_real_, Connectedness = NA_real_,
                      Modularity = NA_real_, N_Modules = NA_integer_,
                      Heterogeneity = NA_real_, stringsAsFactors = FALSE))
  }
  deg  <- igraph::degree(g)
  comp <- igraph::components(g)
  connectedness <- if (comp$no == 1) 1.0 else
    sum(vapply(comp$csize, function(s) s * (s - 1) / 2, numeric(1))) /
    (igraph::vcount(g) * (igraph::vcount(g) - 1) / 2)
  mr <- tryCatch(igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight)),
                 error = function(e) NULL)
  data.frame(
    Vegetation = xv, Cell = cl,
    Nodes = igraph::vcount(g), Edges = igraph::ecount(g),
    Avg_Degree = mean(deg),
    Density = igraph::edge_density(g),
    Clustering = igraph::transitivity(g, type = "global"),
    Connectedness = connectedness,
    Modularity = if (!is.null(mr)) igraph::modularity(mr) else NA_real_,
    N_Modules = if (!is.null(mr)) length(unique(igraph::membership(mr))) else NA_integer_,
    Heterogeneity = if (mean(deg) > 0) stats::sd(deg) / mean(deg) else NA_real_,
    stringsAsFactors = FALSE)
}

#' 节点周转（源：5_network_structure.R:2261-2348）
local80_turnover <- function(CFG, otu_filt, meta, smp, cells, turnover_by,
                             group_factors, min_cor, min_prev) {
  f <- turnover_by %||% (group_factors[2] %||% (CFG$design$factors %||% character(0))[2])
  if (is.null(f) || !f %in% names(meta)) {
    f <- (CFG$design$factors %||% character(0))[1]
  }
  empty <- list(flows = NULL, jaccard = NULL, persistence = NULL, cells = cells,
                sankey = NULL, factor = f, factor2 = NA_character_,
                levels = character(0))
  if (is.null(f) || !f %in% names(meta)) return(empty)

  lv <- local80_levels(CFG, meta, f)
  lv <- lv[vapply(lv, function(x)
    sum(as.character(meta[[f]]) == x, na.rm = TRUE) >= 4, logical(1))]
  if (length(lv) < 2) return(empty)

  net_nodes <- list()
  for (x in lv) {
    idx <- which(as.character(meta[[f]]) == x)
    g <- tryCatch(local80_group_net(otu_filt[, idx, drop = FALSE], min_cor, min_prev),
                  error = function(e) NULL)
    net_nodes[[x]] <- if (is.null(g)) character(0) else igraph::V(g)$name
  }

  # 相邻水平之间的 共有 / 丢失 / 新增
  flows <- do.call(rbind, lapply(seq_len(length(lv) - 1), function(i) {
    a <- net_nodes[[lv[i]]]; b <- net_nodes[[lv[i + 1]]]
    data.frame(From = lv[i], To = lv[i + 1],
               Category = c("Shared", "Lost", "Gained"),
               N = c(length(intersect(a, b)), length(setdiff(a, b)),
                     length(setdiff(b, a))),
               stringsAsFactors = FALSE)
  }))

  # Jaccard 重叠矩阵
  jac <- matrix(NA_real_, length(lv), length(lv), dimnames = list(lv, lv))
  for (i in seq_along(lv)) for (j in seq_along(lv)) {
    a <- net_nodes[[lv[i]]]; b <- net_nodes[[lv[j]]]
    jac[i, j] <- if (length(a) == 0 || length(b) == 0) 0 else
      length(intersect(a, b)) / length(union(a, b))
  }

  # 节点在多少个网络里出现
  all_nodes <- unique(unlist(net_nodes))
  if (length(all_nodes) == 0) return(empty)
  np <- vapply(all_nodes, function(nd) sum(vapply(net_nodes, function(x) nd %in% x,
                                                  logical(1))), numeric(1))
  pers <- data.frame(ASV_ID = all_nodes, N_Networks = as.numeric(np),
                     stringsAsFactors = FALSE)
  pers <- merge(pers, cells[, c("ASV_ID", "Biosphere", "Mean_RA")],
                by = "ASV_ID", all.x = TRUE, sort = FALSE)

  # Sankey 三轴：因子（Vegetation）x 细胞类型（Cell）x 设计因子 2（Landform）
  f2 <- (CFG$design$factors %||% character(0))[2]
  if (is.null(f2) || !f2 %in% names(meta) || identical(f2, f)) 
    f2 <- (CFG$design$factors %||% character(0))[1]
  lv2 <- if (!is.null(f2) && f2 %in% names(meta)) local80_levels(CFG, meta, f2) else character(0)
  sk <- NULL
  if (length(lv2) > 0) {
    rel_all <- md_rel(md_core(CFG)$asv_raw) * 100
    rows <- list()
    for (x in lv) {
      nodes_v <- net_nodes[[x]]
      if (length(nodes_v) == 0) next
      for (z in lv2) {
        sam <- which(as.character(meta[[f]]) == x & as.character(meta[[f2]]) == z)
        if (length(sam) == 0) next
        det <- rownames(rel_all)[rowSums(rel_all[, sam, drop = FALSE] > 0.01) > 0]
        det_nodes <- intersect(nodes_v, det)
        if (length(det_nodes) == 0) next
        tb <- table(as.character(cells$Biosphere[match(det_nodes, cells$ASV_ID)]))
        for (cn in names(tb)) {
          rows[[length(rows) + 1L]] <- data.frame(
            Axis1 = x, Axis2 = cn, Axis3 = z, Count = as.numeric(tb[[cn]]),
            stringsAsFactors = FALSE)
        }
      }
    }
    if (length(rows)) {
      sk <- do.call(rbind, rows)
      sk$Axis1 <- factor(sk$Axis1, levels = lv)
      sk$Axis2 <- factor(sk$Axis2, levels = names(MD_NET$cell_labels),
                         labels = unname(MD_NET$cell_labels))
      sk$Axis3 <- factor(sk$Axis3, levels = lv2)
    }
  }
  list(flows = flows, jaccard = local80_jaccard_df(jac), persistence = pers,
       cells = cells, sankey = sk, factor = f, factor2 = f2, levels = lv)
}

#' 矩阵 -> 长表（Jaccard 写表用）
local80_jaccard_df <- function(m) {
  d <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(d) <- c("From", "To", "Jaccard")
  d
}

#' bootstrap 均衡后的组均值 ± 95% CI（源：5_network_structure.R:413-469）
local80_bootstrap <- function(sample_metrics, factors, metrics_main, n_boot) {
  if (length(factors) == 0 || length(metrics_main) == 0) return(NULL)
  rows <- list()
  for (f in factors) {
    if (!f %in% names(sample_metrics)) next
    tgt <- min(table(sample_metrics[[f]]))
    if (!is.finite(tgt) || tgt < 2) next
    for (m in metrics_main) {
      if (!m %in% names(sample_metrics)) next
      vals <- sample_metrics[[m]]; grp <- as.character(sample_metrics[[f]])
      lvs <- unique(grp[!is.na(grp)])
      mat <- matrix(NA_real_, n_boot, length(lvs), dimnames = list(NULL, lvs))
      for (g in lvs) {
        gv <- vals[grp == g & !is.na(grp)]
        gv <- gv[is.finite(gv)]
        if (length(gv) < 2) next
        for (b in seq_len(n_boot)) {
          idx <- if (length(gv) >= tgt) sample(seq_along(gv), tgt, replace = FALSE) else
            sample(seq_along(gv), tgt, replace = TRUE)
          mat[b, g] <- mean(gv[idx], na.rm = TRUE)
        }
      }
      rows[[length(rows) + 1L]] <- data.frame(
        Metric = m, Factor = f, Group = lvs,
        BootMean = colMeans(mat, na.rm = TRUE),
        BootSD   = apply(mat, 2, stats::sd, na.rm = TRUE),
        CI_lower = apply(mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
        CI_upper = apply(mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
        Target_n = tgt, stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' 度量是否有任何有限值（用于剔除全 NA 面板）
local80_ok_metrics <- function(sample_metrics, metrics) {
  metrics[vapply(metrics, function(m)
    m %in% names(sample_metrics) && any(is.finite(sample_metrics[[m]])), logical(1))]
}

#' 把 long 表的 MetricLabel 重排到给定度量顺序（facet 出场顺序）
local80_reorder_metric <- function(df, metrics, labels) {
  df$MetricLabel <- factor(local80_lab(as.character(df$Metric), labels),
                           levels = unname(labels[metrics]))
  df[!is.na(df$MetricLabel), , drop = FALSE]
}

# =============================================================================
# local80_ 辅助函数（出图）
# =============================================================================

#' 图 1-A：betweenness top-N 条形（源：5_network_structure.R:710-722 + 1369-1388）
local80_plot_betweenness <- function(CFG, nodes_tbl, phylo_colors, st, pt) {
  top_n <- as.integer(st$top_betweenness %||% 20L)
  d <- nodes_tbl[order(-nodes_tbl$betweenness), , drop = FALSE]
  d <- utils::head(d, top_n)
  if (nrow(d) == 0) return(md_empty_plot(CFG, "No node centrality"))
  lab <- ifelse(is.na(d$Genus) | !nzchar(d$Genus), d$ASV_ID, d$Genus)
  dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
  lab[dup] <- paste0(lab[dup], " [", d$ASV_ID[dup], "]")
  d$Label <- lab
  d$Label <- factor(d$Label, levels = rev(d$Label))
  # 门的 factor 只含**本面板实际出现**的水平 —— scale_fill_manual(values = 全套,
  # drop = FALSE) 会为缺席的门生成一个 NA 图例项，与 B/C 面板的 Phylum 图例对不上
  d$Phylum <- factor(d$Phylum, levels = intersect(names(phylo_colors),
                                                  unique(d$Phylum)))
  d <- d[!is.na(d$Phylum), , drop = FALSE]
  pal <- phylo_colors[levels(d$Phylum)]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$Label, y = .data$betweenness,
                                  fill = .data$Phylum)) +
    ggplot2::geom_col(width = 0.7, colour = "white") +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = "Phylum") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Betweenness centrality (normalized)",
                  title = sprintf("Top %d ASVs by betweenness centrality", nrow(d)),
                  tag = "A") +
    theme_md(CFG = CFG) +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = pt$legend),
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold",
                                                      size = pt$title))
}

#' 图 1-B：节点度 vs 平均相对丰度（源：5_network_structure.R:1390-1410）
#' 原文的 ggpubr::stat_cor() 换成本地 Pearson r / p 注释（同样的数值，纯 ASCII）
local80_plot_degree_abund <- function(CFG, otu_filt, kept_asv, g_global, tax_filt,
                                      phylo_colors, pt) {
  node_abund <- rowMeans(otu_filt[kept_asv, , drop = FALSE])
  deg_norm   <- igraph::degree(g_global, normalized = TRUE)
  d <- data.frame(ASV_ID = kept_asv, MeanAbundance = node_abund,
                  Degree = as.numeric(deg_norm[kept_asv]),
                  stringsAsFactors = FALSE)
  d$Phylum <- tax_filt$Phylum[match(d$ASV_ID, tax_filt$ASV_ID)]
  d$Phylum[is.na(d$Phylum) | !nzchar(d$Phylum)] <- "Unknown"
  d$x <- log10(d$MeanAbundance + 1e-6)
  d <- d[is.finite(d$x) & is.finite(d$Degree), , drop = FALSE]
  if (nrow(d) < 3) return(md_empty_plot(CFG, "Degree-abundance unavailable"))
  d$Phylum <- factor(d$Phylum, levels = intersect(names(phylo_colors),
                                                  unique(d$Phylum)))
  d <- d[!is.na(d$Phylum), , drop = FALSE]
  if (nrow(d) < 3) return(md_empty_plot(CFG, "Degree-abundance unavailable"))
  ct <- tryCatch(stats::cor.test(d$x, d$Degree, method = "pearson"),
                 error = function(e) NULL)
  lab <- if (is.null(ct)) "r = NA" else
    sprintf("r = %.2f, p %s", unname(ct$estimate),
            ifelse(ct$p.value < 0.001, "< 0.001",
                   paste0("= ", format(ct$p.value, digits = 2))))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$Degree,
                                  colour = .data$Phylum)) +
    ggplot2::geom_point(size = MD_PT_SIZE$point, alpha = MD_PT_SIZE$point_alpha,
                         stroke = 0) +
    # 每门一条 lm：门内点数 < 4 时不画（点太少时 qt() 自由度异常 → "NaNs produced"，
    # 且空 CI ribbon 会报 "no non-missing arguments to max"；两条都只是噪声告警）。
    ggplot2::geom_smooth(
      data = d[ave(seq_len(nrow(d)), d$Phylum, FUN = length) >= 4, , drop = FALSE],
      ggplot2::aes(colour = .data$Phylum), method = "lm",
      formula = y ~ x, se = TRUE, linewidth = 0.6, alpha = 0.1,
      show.legend = FALSE) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         colour = "grey20", linewidth = 0.8, alpha = 0.15,
                         inherit.aes = FALSE,
                         ggplot2::aes(x = .data$x, y = .data$Degree)) +
    ggplot2::scale_colour_manual(values = phylo_colors[levels(d$Phylum)],
                                 drop = FALSE, name = "Phylum") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6,
                      label = lab, size = pt2gg(pt$ann), colour = "grey25") +
    ggplot2::labs(x = "log10(Mean relative abundance %)", y = "Normalized degree",
                  title = "Node connectivity vs abundance",
                  tag = "B") +
    theme_md(CFG = CFG) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold",
                                                      size = pt$title)) +
    ggplot2::guides(colour = md_guide_flex(length(phylo_colors), "bottom"))
}

#' 图 1-C：Zi-Pi 四区散点（源：5_network_structure.R:1681-1713）
local80_plot_zipi <- function(CFG, zi_pi_df, pt, zi_thr, pi_thr) {
  d <- zi_pi_df[is.finite(zi_pi_df$Zi) & is.finite(zi_pi_df$Pi), , drop = FALSE]
  if (nrow(d) < 3) return(md_empty_plot(CFG, "Zi-Pi unavailable"))
  d$lab <- ifelse(!is.na(d$Genus) & nzchar(d$Genus), d$Genus, as.character(d$ASV_ID))
  role_cols <- MD_NET$role_colors
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Pi, y = .data$Zi,
                                       colour = .data$role)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.5) +
    ggplot2::geom_hline(yintercept = zi_thr, linetype = "dotted",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = pi_thr, linetype = "dotted",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = role_cols, name = "Node role",
                                 drop = FALSE) +
    ggplot2::annotate("text", x = 0.02, y = Inf, label = "Module Hub",
                      colour = role_cols[["Module Hub"]], size = pt2gg(pt$ann),
                      fontface = "bold", hjust = 0, vjust = 1.6) +
    ggplot2::annotate("text", x = 0.98, y = Inf, label = "Network Hub",
                      colour = role_cols[["Network Hub"]], size = pt2gg(pt$ann),
                      fontface = "bold", hjust = 1, vjust = 1.6) +
    ggplot2::annotate("text", x = 0.98, y = -Inf, label = "Connector",
                      colour = role_cols[["Connector"]], size = pt2gg(pt$ann),
                      fontface = "bold", hjust = 1, vjust = -1.0) +
    ggplot2::annotate("text", x = 0.02, y = -Inf, label = "Peripheral",
                      colour = role_cols[["Peripheral"]], size = pt2gg(pt$ann),
                      fontface = "bold", hjust = 0, vjust = -1.0) +
    ggplot2::labs(x = expression(italic("P")[i] ~ "(among-module connectivity)"),
                  y = expression(italic("Z")[i] ~ "(within-module connectivity)"),
                  title = "Node role classification",
                  tag = "C") +
    theme_md(CFG = CFG) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold",
                                                      size = pt$title))
  # ggrepel 可用时给非 Peripheral 节点加标签（原文口径）
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    g <- g + ggrepel::geom_text_repel(
      data = d[d$role != "Peripheral", , drop = FALSE],
      ggplot2::aes(label = .data$lab), size = pt2gg(pt$ann), max.overlaps = 30,
      box.padding = 0.35, point.padding = 0.25, show.legend = FALSE,
      colour = "grey20", min.segment.length = 0.4)
  }
  g
}

#' 一组面板拼成一张图。
#' 与共享库 md_compose() 的唯一差别：**不自动挂 A/B/C... 面板标签**
#' （tag_levels = NULL）—— 模块 80 的拓扑合图里 B 面板内部本身就是 8 张子图的
#' patchwork，自动编号会把 8 个子图与图例各编一个号，导致主面板编号与图注对不上。
#' 其余（design 布局、guides 收集、拼图字号、图例位置）逐项与 md_compose() 一致。
local80_compose <- function(CFG, plots, design, guides = "collect", legend = "bottom",
                            theme_comp = TRUE, heights = NULL, widths = NULL) {
  plots <- Filter(function(x) inherits(x, "gg") || inherits(x, "patchwork"), plots)
  if (length(plots) == 0) return(md_empty_plot(CFG, "no panels"))
  # heights / widths：design 里各**行 / 列**的权重。**默认是等分** —— 所以一条又扁又宽的
  #   图例带如果和网络格子各占一行，会白白吃掉一半高度。要显式给权重（见 80b）。
  pw <- Reduce(`+`, plots) +
    patchwork::plot_layout(design = design, guides = guides,
                           heights = heights, widths = widths) +
    patchwork::plot_annotation(tag_levels = NULL)
  # ⚠ theme_comp = FALSE 时**不要**统一套 theme_md_comp()：
  #   它会把 axis.text / axis.ticks / panel.grid 重新打开，把网络快照与图例带
  #   （local80_net_gg / local80_net_legend 里已经 element_blank() / theme_void()）
  #   的坐标轴又画回来 —— 见 references/04-pitfalls.md 的 C32。
  #   这类图请自己在**需要缩小字号的**那一个面板上单独加 theme_md_comp(CFG)。
  if (isTRUE(theme_comp)) pw <- pw & theme_md_comp(CFG)
  if (!is.null(legend)) pw <- pw & ggplot2::theme(legend.position = legend)
  pw
}

#' 图 4-A：全局网络快照（按门着色、正负边着色）
local80_panel_global <- function(CFG, g_global, st, phylum_of, phylo_colors, pt) {
  lay <- switch(st$layout %||% "graphopt",
                circle = igraph::layout_in_circle(g_global),
                nicely = igraph::layout_nicely(g_global),
                fr     = local80_layout_fr(g_global),
                igraph::layout_with_graphopt(g_global))
  lay <- local80_radial(lay)
  sub <- sprintf("%d ASVs | %d edges | |r| threshold applied", igraph::vcount(g_global),
                 igraph::ecount(g_global))
  local80_net_gg(CFG, g_global, st, phylum_of, phylo_colors, lay = lay,
                 title = "Global ASV co-occurrence network", subtitle = sub,
                 tag = "A")
}

#' 图 4-B：分组网络 2x4（每个网络一张 ggplot 面板拼成）
#' @param ncol 每行几个网络（原图 2 行 x 4 列）
#' @param tag  子面板标签（原图 Fig5c 的 8 个网络均无独立标签，默认 NULL）
local80_panel_groups <- function(CFG, nets, labels, ns, st, phylum_of,
                                 phylo_colors, ncol = NULL, tag = NULL) {
  keys <- names(nets)
  if (length(keys) == 0) return(md_empty_plot(CFG, "No group network"))
  # 格子数自适应：6 个网络 = 2 行 x 3 列（原来写死 4 列，6 个会排成 4 + 2 很难看）；
  # ceil(sqrt(n)) 对 4/6/8/9/12 分别给出 2/3/3/3/4，都是接近方形的排布。
  if (is.null(ncol)) ncol <- max(1L, min(4L, ceiling(sqrt(length(keys)))))
  plist <- lapply(keys, function(k) {
    g <- nets[[k]]
    if (is.null(g) || igraph::ecount(g) < 2) {
      return(md_empty_plot(CFG, paste0(labels[[k]], "  (no valid network)")))
    }
    lay <- local80_layout_fr(g)
    ttl <- sprintf("%s (n=%d, N=%d, E=%d)", labels[[k]], ns[[k]],
                   igraph::vcount(g), igraph::ecount(g))
    local80_net_gg(CFG, g, st, phylum_of, phylo_colors, lay = lay,
                   title = ttl, tag = tag, guides_fill = "none")
  })
  patchwork::wrap_plots(plist, ncol = ncol)
}

#' 图 4-C：观测 vs 随机网络（源：5_network_structure.R:1741-1757）
local80_panel_random <- function(CFG, comparison_df, random_n, pt) {
  if (is.null(comparison_df) || nrow(comparison_df) == 0)
    return(md_empty_plot(CFG, "Random comparison unavailable"))
  d <- comparison_df
  # 原脚本把 NaN（Rand_SD = 0 的确定性零模型，如 fixed n/m 的 Avg Degree 与
  # Density）经 case_when 落进 "Not significant" —— 这里**保持同一文案与色值**
  # （`MD_NET$random_cols` 逐字照抄 5_network_structure.R:1750-1752），
  # 但把"哪些度量的 Z 检验退化"显式写进 Note sheet，不让读者把它误读成"检验过、不显著"。
  d$Interpretation[is.na(d$Interpretation)] <- "Not significant"
  d$Metric <- factor(d$Metric,
                     levels = c("Avg Degree", "Density", "Avg Path",
                                "Clustering Coef", "Modularity"))
  cols <- MD_NET$random_cols
  ggplot2::ggplot(d, ggplot2::aes(x = .data$Metric, y = .data$Real,
                                  fill = .data$Interpretation)) +
    ggplot2::geom_bar(stat = "identity", width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Rand_Mean - .data$Rand_SD,
                                        ymax = .data$Rand_Mean + .data$Rand_SD),
                           width = 0.2, colour = "grey40", linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = .data$Rand_Mean), shape = 21, size = 3,
                        fill = "grey50") +
    ggplot2::scale_fill_manual(values = cols, name = "Comparison") +
    ggplot2::labs(x = NULL, y = "Value",
                  title = sprintf("Observed vs random networks (n = %d)", random_n),
                  subtitle = "Points with error bars = random network mean +/- SD",
                  tag = "B") +
    theme_md(CFG = CFG) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                   legend.position = "bottom",
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold",
                                                      size = pt$title),
                   plot.subtitle = ggplot2::element_text(size = pt$ann,
                                                         colour = "grey40"))
}

#' 图 4-D：节点周转 Sankey（源：5_network_structure.R:2321-2376）
local80_panel_sankey <- function(CFG, turnover, pt) {
  sk <- turnover$sankey
  if (is.null(sk) || nrow(sk) == 0 ||
      !requireNamespace("ggalluvial", quietly = TRUE))
    return(md_empty_plot(CFG, "Node turnover Sankey unavailable"))
  bio_cols <- md_pal(CFG, names(MD_NET$cell_labels), "quadrant")
  names(bio_cols) <- unname(MD_NET$cell_labels)
  x_lab <- c(turnover$factor, "Cell", turnover$factor2)
  ggplot2::ggplot(sk, ggplot2::aes(axis1 = .data$Axis1, axis2 = .data$Axis2,
                                   axis3 = .data$Axis3, y = .data$Count)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$Axis2), width = 1 / 4,
                              alpha = 0.7, knot.pos = 0.4) +
    ggalluvial::geom_stratum(width = 1 / 4, colour = "black", linewidth = MD_LW$axis) +
    ggplot2::geom_text(stat = "stratum",
                       ggplot2::aes(label = after_stat(
                         ifelse(prop >= 0.02, as.character(stratum), ""))),
                       size = pt2gg(pt$ann)) +
    ggplot2::scale_x_discrete(limits = x_lab, expand = c(0.12, 0.08)) +
    ggplot2::scale_fill_manual(values = bio_cols, name = "Cell", drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Number of network nodes",
                  title = "Network node composition flow", tag = "A") +
    theme_md(CFG = CFG) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold",
                                                      size = pt$title))
}

#' 图 4-E：丰富/稀有亚网络指标对比（源：5_network_structure.R:2200-2246）
local80_panel_abundant_rare <- function(CFG, abr, pt) {
  long <- abr$grid
  if (is.null(long) || nrow(long) == 0)
    return(md_empty_plot(CFG, "Abundant/rare sub-networks unavailable"))
  bio_cols <- md_pal(CFG, names(MD_NET$cell_labels), "quadrant")
  names(bio_cols) <- unname(MD_NET$cell_labels)
  long$CellDisp <- factor(unname(MD_NET$cell_labels[as.character(long$Cell)]),
                          levels = unname(MD_NET$cell_labels))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$Vegetation, y = .data$Value,
                                     fill = .data$CellDisp)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6,
                      colour = "black", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = bio_cols, name = "Cell", drop = FALSE) +
    ggplot2::facet_wrap(~ MetricLabel, scales = "free_y", ncol = 3) +
    ggplot2::labs(x = abr$xvar %||% NULL, y = NULL,
                  title = "Abundant / rare sub-network topology", tag = "B") +
    theme_md(CFG = CFG) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = pt$strip, face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1,
                                                       size = pt$tick),
                   legend.position = "bottom",
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold",
                                                      size = pt$title))
}

#' Note sheet
local80_note <- function(CFG, method, method_note, min_cor, sub_min_cor, min_prev,
                         max_asv, algo, zi_thr, pi_thr, random_n, nsim, n_boot,
                         metrics_main, metrics_supp, supp_fig, mf, otu_filt,
                         g_global, modularity_global, n_modules_global,
                         comparison_df, notes) {
  degenerate <- if (!is.null(comparison_df))
    paste(comparison_df$Metric[is.na(comparison_df$Z_score)], collapse = ", ") else ""
  base <- data.frame(
    item = c("Network method", "Correlation", "|r| threshold (global)",
             "|r| threshold (sub-networks)", "Min occurrence", "Max ASVs",
             "Module algorithm", "Zi / Pi thresholds", "Random networks",
             "Robustness simulations", "Bootstrap replicates",
             "Group factors", "Metrics main", "Metrics supp",
             "Supp figure", "Nodes", "Edges", "Modules", "Modularity",
             "Degenerate Z tests", "Seed", "Row counts"),
    value = c(method,
              as.character(CFG$network$cor_method %||% "spearman"),
              as.character(min_cor), as.character(sub_min_cor),
              as.character(min_prev),
              as.character(max_asv %||% "unlimited"),
              as.character(algo),
              sprintf("Zi > %g / Pi > %g (Deng et al. 2012)", zi_thr, pi_thr),
              sprintf("Erdos-Renyi sample_gnm(n, m), n = %d", random_n),
              sprintf("node removal until largest component < 50%%, n = %d", nsim),
              sprintf("group means resampled to the smallest group size, B = %d", n_boot),
              paste(mf, collapse = " | "),
              paste(metrics_main, collapse = ", "),
              paste(metrics_supp, collapse = ", "),
              as.character(supp_fig),
              as.character(igraph::vcount(g_global)),
              as.character(igraph::ecount(g_global)),
              as.character(n_modules_global),
              as.character(round(modularity_global, 4)),
              if (nzchar(degenerate)) paste0(degenerate,
                " (zero SD between random replicates: the Z test is degenerate, not 'not significant')") else "none",
              as.character(CFG$opts$seed),
              sprintf("ASV for network = %d after occurrence filter", nrow(otu_filt))),
    stringsAsFactors = FALSE)
  extra <- c(
    if (!is.null(method_note)) method_note,
    "Seed fix: robustness_R() simulations are seeded BEFORE the 30 node-removal runs (the source script seeded after them, which made the metric non-reproducible).",
    "Layout: network layout coordinates are seeded in the plotting section (offset network_layout), so figures are reproducible without recomputing the network.",
    "Reproducibility: this module writes its own sheets; it does not depend on the upstream 5-1..5-8 tables.",
    "Cross-module: cache/80_network.rds (nodes / edges / params) is consumed by module 81.",
    notes)
  if (length(extra)) {
    base <- rbind(base, data.frame(item = paste0("Note", seq_along(extra)),
                                   value = extra, stringsAsFactors = FALSE))
  }
  base
}
