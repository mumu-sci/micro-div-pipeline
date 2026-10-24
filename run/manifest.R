# =============================================================================
# run/manifest.R —— 单一模块清单（**唯一真身**）
# -----------------------------------------------------------------------------
# 90_run_all.R 与 91_regenerate_figures.R 都由本清单生成，不再各写一份
# （原始研究的 _run_all.R 与 _regenerate_figures.R 各复制了一份清单，后来漂移了）。
#
# 每个条目：
#   id      模块号（同时也是文件名前缀与输出前缀）
#   file    R/ 下的文件名
#   fun     该文件导出的函数名（签名统一 md_xxx(CFG)）
#   after   依赖的模块 id；只对**清单里真实存在的 id** 构成约束（拓扑排序用）
#   heavy   TRUE = 重计算，figures_only 模式下跳过
#   needs   需要的包组
#   last    TRUE = 永远排在最后（收尾模块）
# =============================================================================

MD_MANIFEST <- list(
  list(id = "01", file = "01_data_health.R",        fun = "md_data_health",
       after = character(0),        heavy = FALSE, needs = c("ecology")),
  list(id = "05", file = "05_assemble.R",           fun = "md_assemble",
       after = "01",                heavy = FALSE, needs = c("phyloseq", "ecology")),

  # ---- 组成与多样性 ----
  list(id = "10", file = "10_composition.R",        fun = "md_composition",
       after = "05",                heavy = FALSE, needs = c("phyloseq")),
  list(id = "11", file = "11_alpha.R",              fun = "md_alpha",
       after = "05",                heavy = FALSE, needs = c("phyloseq", "ecology")),
  # 12 = 2026 补充轮新增：生态位宽度（Levins' B）与生态策略组成
  #      参考实现（原脚本）
  list(id = "12", file = "12_niche_breadth.R",      fun = "md_niche_breadth",
       after = "05",                heavy = FALSE, needs = c("phyloseq")),
  list(id = "13", file = "13_distances.R",          fun = "md_distances",
       after = "05",                heavy = FALSE, needs = c("ecology", "unifrac")),
  list(id = "14", file = "14_ordination.R",         fun = "md_ordination",
       after = "13",                heavy = FALSE, needs = c("ecology")),
  list(id = "15", file = "15_beta_tests.R",         fun = "md_beta_tests",
       after = "13",                heavy = FALSE, needs = c("ecology")),
  list(id = "16", file = "16_beta_components.R",    fun = "md_beta_components",
       after = "13",                heavy = FALSE, needs = c("ecology")),
  list(id = "17", file = "17_gamma_partition.R",    fun = "md_gamma_partition",
       after = "05",                heavy = TRUE,  needs = c("partition")),

  # ---- 系统发育（无树时模块自行降级跳过）----
  list(id = "20", file = "20_phylo_prep.R",         fun = "md_phylo_prep",
       after = "05",                heavy = FALSE, needs = c("ecology")),
  list(id = "21", file = "21_phylo_diversity.R",    fun = "md_phylo_diversity",
       after = "20",                heavy = TRUE,  needs = c("ecology", "parallel")),
  list(id = "22", file = "22_phylo_beta.R",         fun = "md_phylo_beta",
       after = "20",                heavy = FALSE, needs = c("ecology")),

  # ---- 群落组装 ----
  list(id = "30", file = "30_null_models.R",        fun = "md_null_models",
       after = "05",                heavy = TRUE,  needs = c("ecology", "parallel")),
  list(id = "31", file = "31_assembly_process.R",   fun = "md_assembly_process",
       after = "30",                heavy = FALSE, needs = c("ecology")),
  list(id = "33", file = "33_distance_decay.R",     fun = "md_distance_decay",
       after = c("13", "05"),       heavy = FALSE, needs = c("ecology", "mantel")),
  list(id = "34", file = "34_variance_partition.R", fun = "md_variance_partition",
       after = "05",                heavy = FALSE, needs = c("ecology", "venn")),

  # ---- 环境关联 ----
  list(id = "40", file = "40_env_summary.R",        fun = "md_env_summary",
       after = "05",                heavy = FALSE, needs = c("ecology")),
  # 注：原 42_constrained_ord（RDA/CCA/db-RDA）已移出管线 —— 那是**本包自行扩展**的能力，
  #     原始研究里没有做约束排序。文件归档在 _build-notes/micro-div-pipeline/removed-code/。
  list(id = "43", file = "43_mantel_mrm.R",         fun = "md_mantel_mrm",
       after = "13",                heavy = FALSE, needs = c("mantel")),
  # 44 / 45 = 2026 补充轮新增
  #   44 约束排序（db-RDA / RDA / CCA / capscale）—— 参考实现
  #   45 GAM / GAMM 响应曲线                —— 参考实现
  list(id = "44", file = "44_constrained_ord.R",    fun = "md_constrained_ord",
       after = "05",                heavy = FALSE, needs = c("ecology")),
  list(id = "45", file = "45_response_curves.R",    fun = "md_response_curves",
       after = "05",                heavy = FALSE, needs = c("gam")),

  # ---- 建模 ----
  list(id = "50", file = "50_lmm.R",                fun = "md_lmm",
       after = "40",                heavy = FALSE, needs = c("model", "gam")),
  list(id = "52", file = "52_sem.R",                fun = "md_sem",
       after = "40",                heavy = FALSE, needs = c("sem", "model")),
  list(id = "53", file = "53_sensitivity.R",        fun = "md_sensitivity",
       after = "50",                heavy = FALSE, needs = c("model")),

  # ---- 空间 ----
  list(id = "60", file = "60_spatial_autocorr.R",   fun = "md_spatial_autocorr",
       after = "05",                heavy = TRUE,  needs = c("mantel")),
  list(id = "61", file = "61_map.R",                fun = "md_map",
       after = "05",                heavy = FALSE, needs = c("map")),

  # ---- 响应变量关系（response$enabled=FALSE 时模块自行跳过）----
  list(id = "70", file = "70_response_pattern.R",   fun = "md_response_pattern",
       after = "05",                heavy = FALSE, needs = c("ecology")),
  list(id = "71", file = "71_response_linkage.R",   fun = "md_response_linkage",
       after = c("13", "40"),       heavy = FALSE, needs = c("model", "mantel")),

  # ---- 网络与类群补充（2026 补充轮）----
  #   80 共现网络（含丰富/稀有亚网络）—— 参考实现
  #   81 通用类群分组（组成 / 网络角色 / 功能关联 / 跨组边）—— 参考实现
  #   82 ASV/类群共享 UpSet                       —— 参考实现
  #   83 差异类群热图（CLR + KW + Dunn + Cliff's delta）—— 参考实现
  #   注 81 对 80 是**软依赖**：80 跳过/缺缓存时 81 只出组成与表，并在 notes 说明。
  list(id = "80", file = "80_network.R",            fun = "md_network",
       after = "05",                heavy = TRUE,  needs = c("network", "ecology")),
  list(id = "81", file = "81_taxon_group.R",        fun = "md_taxon_group",
       after = c("05", "80"),       heavy = FALSE, needs = c("phyloseq")),
  list(id = "82", file = "82_taxon_share.R",        fun = "md_taxon_share",
       after = "05",                heavy = FALSE, needs = c("upset")),
  list(id = "83", file = "83_taxon_diff.R",         fun = "md_taxon_diff",
       after = "05",                heavy = FALSE, needs = c("taxon_heat")),

  # ---- 收尾 ----
  list(id = "95", file = "95_export_tables.R",      fun = "md_export_tables",
       after = character(0),        heavy = FALSE, needs = c("io"), last = TRUE)
)

MD_MANIFEST_IDS <- function(manifest = MD_MANIFEST)
  vapply(manifest, function(x) x$id, character(1))

#' 拓扑排序后的执行顺序（生产者先于消费者；last=TRUE 的排最后）
md_manifest_order <- function(manifest = MD_MANIFEST) {
  ids  <- MD_MANIFEST_IDS(manifest)
  last <- vapply(manifest, function(x) isTRUE(x$last), logical(1))
  done <- character(0); out <- character(0); guard <- 0L
  todo <- manifest[!last]
  while (length(out) < length(todo) && guard < 1000L) {
    guard <- guard + 1L
    progressed <- FALSE
    for (m in todo) {
      if (m$id %in% done) next
      deps <- intersect(m$after, ids)
      if (all(deps %in% done)) {
        out <- c(out, m$id); done <- c(done, m$id); progressed <- TRUE
      }
    }
    if (!progressed) {
      rest <- setdiff(vapply(todo, function(x) x$id, character(1)), done)
      warning("manifest 依赖存在环，按声明顺序追加: ", paste(rest, collapse = ", "))
      out <- c(out, rest); break
    }
  }
  c(out, MD_MANIFEST_IDS(manifest)[last])
}

md_manifest_get <- function(id, manifest = MD_MANIFEST) {
  for (m in manifest) if (m$id == id) return(m)
  NULL
}

#' 模块文件的绝对路径
md_module_file <- function(CFG, m) file.path(CFG$project$dir, "R", m$file)
