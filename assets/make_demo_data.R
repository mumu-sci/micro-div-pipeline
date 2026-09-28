# =============================================================================
# assets/make_demo_data.R
# 从 nifH 论文源数据抽取示范子集：3 Zone x 4 Site x 3 Vegetation = 36 样本 x 300 ASV
# -----------------------------------------------------------------------------
# 用途：本 skill 的示范/自检数据。**只读源数据，只写本包 data/ 目录。**
#
# 抽样口径（详见 skill 外 _build-notes/micro-div-pipeline/DECISIONS.md）：
#   Zone      : Coastal / Nearshore / Xisha（保留「离大陆由近到远」的空间梯度；
#               不用 Nansha —— 它是人工岛，论文里专门做过剔除敏感性分析）
#   Vegetation: Bare / Herb / Tree（不用 Shrub —— Coastal/Nearshore 的 Shrub 只有 3 个
#               样本，凑不齐 4 重复）
#   每格样本  : 4 个，优先取自**不同 Site**（轮转取样，seed 固定）
#   ASV      : 按 36 个样本内的平均相对丰度取前 300
#   树       : 从 tree_full.nwk 用 keep.tip 取这 300 个 ASV
#
# 运行：
#   Rscript --vanilla assets/make_demo_data.R
# 可用环境变量 MD_SRC 覆盖源数据目录（默认指向论文）。
# =============================================================================

.libPaths(c("D:/R-4.5.2/library", .libPaths()))
Sys.setlocale("LC_CTYPE", "")

suppressPackageStartupMessages({
  library(readxl); library(ape)
})

# ---- 参数（集中在此，改这里即可）----
SEED      <- 123L
ZONES     <- c("Coastal", "Nearshore", "Xisha")
VEGS      <- c("Bare", "Herb", "Tree")
N_REP     <- 4L          # 每 Zone x Vegetation 格抽几个样本
N_ASV     <- 300L        # ASV 数目标（实际由「每样本 top-k 并集」决定，见 [3a]）
SOIL_VARS <- c("SWC", "NH4", "NO3", "pH", "OM", "TN", "TP", "EC", "TK", "Ca", "Fe")
BIO_VARS  <- c("bio_1", "bio_12")   # 年均温 / 年降水

SRC <- Sys.getenv("MD_SRC", "D:/01-Academic/06-nifH_BEF/data")

# 定位包根：优先命令行 --pkg=，其次脚本自身路径，最后工作目录
arg_pkg <- grep("^--pkg=", commandArgs(TRUE), value = TRUE)
if (length(arg_pkg)) {
  PKG <- normalizePath(sub("^--pkg=", "", arg_pkg[1]), mustWork = FALSE)
} else {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) {
    PKG <- normalizePath(file.path(dirname(sub("^--file=", "", fa[1])), ".."), mustWork = FALSE)
  } else {
    PKG <- getwd()
    if (!dir.exists(file.path(PKG, "assets"))) PKG <- dirname(PKG)
  }
}
OUT <- file.path(PKG, "data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

set.seed(SEED)
cat("===== make_demo_data.R =====\n")
cat("源数据:", SRC, "\n输出  :", OUT, "\nseed  :", SEED, "\n\n")

stopifnot(dir.exists(SRC))
need <- c("allsample_ASV_taxonomy.txt", "样本顺序及其分组.xlsx",
          "环境数据.xlsx", "固氮酶活性数据.xlsx", "tree_full.nwk")
miss <- need[!file.exists(file.path(SRC, need))]
if (length(miss)) stop("源数据缺少: ", paste(miss, collapse = ", "))

# ---- 1. 样本元数据 ----
meta <- as.data.frame(suppressMessages(
  read_excel(file.path(SRC, "样本顺序及其分组.xlsx"), sheet = "sample_metadata")))
names(meta) <- c("Order", "sample_id", "Site", "Vegetation", "Zone")
meta <- meta[, c("sample_id", "Site", "Vegetation", "Zone")]
meta$sample_id  <- trimws(as.character(meta$sample_id))
meta$Site       <- trimws(as.character(meta$Site))
meta$Vegetation <- trimws(as.character(meta$Vegetation))
meta$Zone       <- trimws(as.character(meta$Zone))
stopifnot(!anyDuplicated(meta$sample_id))

keep <- meta$Zone %in% ZONES & meta$Vegetation %in% VEGS
pool <- meta[keep, , drop = FALSE]
cat("[1] 候选样本池:", nrow(pool), "\n")
print(table(pool$Zone, pool$Vegetation))

# ---- 2. 每格轮转取样（优先不同 Site）----
pick_cell <- function(df, n) {
  # 按 Site 分组轮转，保证尽量分散到不同 Site
  spl <- split(df, df$Site)
  spl <- spl[order(names(spl))]
  out <- list(); k <- 1L
  repeat {
    added <- FALSE
    for (nm in names(spl)) {
      if (nrow(spl[[nm]]) == 0) next
      out[[length(out) + 1L]] <- spl[[nm]][1, , drop = FALSE]
      spl[[nm]] <- spl[[nm]][-1, , drop = FALSE]
      added <- TRUE
      if (length(out) >= n) break
    }
    if (length(out) >= n || !added) break
  }
  do.call(rbind, out[seq_len(min(n, length(out)))])
}

sel <- do.call(rbind, lapply(ZONES, function(z) {
  do.call(rbind, lapply(VEGS, function(v) {
    cell <- pool[pool$Zone == z & pool$Vegetation == v, , drop = FALSE]
    if (nrow(cell) < N_REP) stop("格不足: ", z, " / ", v, " 只有 ", nrow(cell))
    cell <- cell[order(cell$Site, cell$sample_id), , drop = FALSE]
    pick_cell(cell, N_REP)
  }))
}))
sel <- sel[order(sel$Zone, sel$Site, sel$Vegetation), , drop = FALSE]
rownames(sel) <- NULL
cat("\n[2] 抽中样本:", nrow(sel), " 站点数:", length(unique(sel$Site)), "\n")
print(table(sel$Zone, sel$Vegetation))
cat("\n每站点样本数:\n"); print(table(sel$Site, sel$Zone))

# ---- 3. ASV 表 ----
f1 <- file.path(SRC, "allsample_ASV_taxonomy.txt")
d1 <- read.table(f1, skip = 1, header = TRUE, sep = "\t", check.names = FALSE,
                 quote = "", comment.char = "")
asv_id <- trimws(as.character(d1[[1]]))
tax_raw <- as.character(d1[[ncol(d1)]])
mat <- as.matrix(d1[, 2:(ncol(d1) - 1)])
rownames(mat) <- asv_id
storage.mode(mat) <- "double"

stopifnot(all(sel$sample_id %in% colnames(mat)))
sub <- mat[, sel$sample_id, drop = FALSE]

# ---- 4. 取 ASV：**每个样本丰度前 k 个 ASV 的并集** ----
# 为什么不用"按平均相对丰度取前 300"（实测，详见 DECISIONS.md D-04）：
#   这批 nifH 数据里 23,251/23,873 个 ASV 在全部 147 样本中只出现 1 次。
#   纯按平均丰度取前 300，抽到的基本是"单样本爆发型"ASV —— 300 个里 240 个只出现
#   在 1 个样本中，还会出现 0 reads 样本、抽平深度只剩个位数，群落完全不可分析。
#   改成"每样本 top-k 并集"后，每个样本都有自己的优势类群被保留，覆盖度与
#   抽平深度都正常，ASV 规模也仍然很小。
# k 的选定：从 1 递增，取第一个使并集 >= N_ASV 的 k。
rel <- sweep(sub, 2, pmax(colSums(sub), 1), "/")
k <- 1L; keep_asv <- rownames(sub)
repeat {
  keep_asv <- unique(unlist(lapply(seq_len(ncol(rel)), function(j)
    order(rel[, j], decreasing = TRUE)[seq_len(min(k, nrow(rel)))])))
  if (length(keep_asv) >= N_ASV || k >= nrow(sub)) break
  k <- k + 1L
}
keep_asv <- rownames(sub)[keep_asv]
cat("\n[3a] 每样本 top-", k, " ASV 的并集 = ", length(keep_asv), " 个 ASV",
    "（目标 ~", N_ASV, "）\n", sep = "")
sub <- sub[keep_asv, , drop = FALSE]
cat("\n[3] 保留 ASV:", nrow(sub), " 占子集 reads:", 
    sprintf("%.1f%%", 100 * sum(sub) / sum(mat[, sel$sample_id])), "\n")
cat("    样本 reads 范围:", min(colSums(sub)), "-", max(colSums(sub)), "\n")
cat("    抽平深度(最小正样本 reads):", min(colSums(sub)[colSums(sub) > 0]), "\n")
cat("    零 reads 样本数:", sum(colSums(sub) == 0), "\n")
cat("    出现样本数 分布: min", min(rowSums(sub > 0)), "中位",
    stats::median(rowSums(sub > 0)), "max", max(rowSums(sub > 0)), "\n")

# ---- 5. taxonomy 解析 ----
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
pre   <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
tax <- do.call(rbind, strsplit(tax_raw[match(keep_asv, asv_id)], "; ", fixed = TRUE))
tax <- tax[, seq_along(ranks), drop = FALSE]
for (i in seq_along(ranks)) {
  x <- tax[, i]
  x <- sub(paste0("^", pre[i]), "", x)
  x <- gsub("^unclassified_", "", x)
  x <- trimws(x)
  x[x == "" | x == "unclassified" | is.na(x)] <- "Unclassified"
  tax[, i] <- x
}
colnames(tax) <- ranks
taxdf <- data.frame(asv_id = keep_asv, tax, stringsAsFactors = FALSE, check.names = FALSE)

# ---- 6. 环境 + 响应变量 ----
soil <- as.data.frame(suppressMessages(
  read_excel(file.path(SRC, "环境数据.xlsx"), sheet = "土壤理化性质")))
geo <- as.data.frame(suppressMessages(
  read_excel(file.path(SRC, "环境数据.xlsx"), sheet = "地理和气候")))
ara <- as.data.frame(suppressMessages(
  read_excel(file.path(SRC, "固氮酶活性数据.xlsx"), sheet = 1)))

norm_id <- function(x) trimws(as.character(x))
soil$MatchID <- norm_id(soil$MatchID); geo$MatchID <- norm_id(geo$MatchID)
ara$MatchID  <- norm_id(ara$MatchID)

env <- data.frame(sample_id = sel$sample_id, stringsAsFactors = FALSE)
env$Lat <- geo$Lat[match(env$sample_id, geo$MatchID)]
env$Lon <- geo$Lon[match(env$sample_id, geo$MatchID)]
for (v in BIO_VARS) env[[if (v == "bio_1") "MAT" else "MAP"]] <-
  geo[[v]][match(env$sample_id, geo$MatchID)]
for (v in SOIL_VARS) env[[v]] <- soil[[v]][match(env$sample_id, soil$MatchID)]
env$EF <- ara$ARA[match(env$sample_id, ara$MatchID)]   # 生态系统功能（原 ARA）

nm <- c(sample_id = "sample_id", Lat = "Latitude", Lon = "Longitude",
        MAT = "MAT", MAP = "MAP")
names(env)[match(names(nm), names(env))] <- nm
stopifnot(!anyNA(env))

# ---- 7. 树 ----
tr <- read.tree(file.path(SRC, "tree_full.nwk"))
cat("\n[4] 树 tip 覆盖保留 ASV:", sprintf("%.1f%%", 100 * mean(keep_asv %in% tr$tip.label)), "\n")
missing_tip <- setdiff(keep_asv, tr$tip.label)
if (length(missing_tip)) {
  cat("    树中缺失", length(missing_tip), "个 ASV，从丰度表与注释表中同步剔除\n")
  keep_asv <- intersect(keep_asv, tr$tip.label)
  sub <- sub[keep_asv, , drop = FALSE]
  taxdf <- taxdf[taxdf$asv_id %in% keep_asv, , drop = FALSE]
  tr <- keep.tip(tr, keep_asv)
} else {
  tr <- keep.tip(tr, keep_asv)
}
tr$node.label <- NULL

# ---- 8. 类群 -> 分组映射表（可选第 6 文件；模块 81 用）----
#   契约：两列 taxon, group。分组含义不限（功能群 / 生态位 / R-K 对策 / 任意自定义），
#   按哪一层级匹配由 config 的 taxon_group$match_by 声明（示范用 Genus）。
#   示范数据是 nifH（固氮）数据，这里给一份「固氮生态类群」的通用示例：
#   共生 / 联合 / 自生 / 其他。换成 16S/ITS 时把这张表整个替换成你的分组即可。
GROUP_MAP <- list(
  Symbiotic = c("Bradyrhizobium", "Frankia", "Mesorhizobium", "Rhizobiales",
                "Azospirillum"),
  Associative = c("Azotobacter", "Azohydromonas", "Methylocystis", "Geobacter",
                  "Pelosinus", "Anaerobacillus", "Clostridium"),
  Free_living = c("Pseudomonas", "Paenibacillus", "Gloeothece", "Cyanobacterium",
                  "Cyanobacteria", "Oscillatoriales", "Oscillatoriophycideae",
                  "Chroococcales", "Ectothiorhodospira", "Acetonema",
                  "Fibrobacter", "Zehria")
)
grp_df <- do.call(rbind, lapply(names(GROUP_MAP), function(g)
  data.frame(taxon = GROUP_MAP[[g]], group = g, stringsAsFactors = FALSE)))
# 未在映射表里出现的 Genus 一律落到 Other（示范口径：drop_unmapped = FALSE 时进组成图）
other_genus <- setdiff(unique(taxdf$Genus), grp_df$taxon)
if (length(other_genus))
  grp_df <- rbind(grp_df, data.frame(taxon = other_genus, group = "Other",
                                     stringsAsFactors = FALSE))
grp_df <- grp_df[order(grp_df$group, grp_df$taxon), , drop = FALSE]
rownames(grp_df) <- NULL

# ---- 9. 落盘（纯 ASCII 内容，UTF-8 无 BOM）----
asv_out <- data.frame(asv_id = rownames(sub), sub, check.names = FALSE, stringsAsFactors = FALSE)
write.csv(asv_out, file.path(OUT, "01_asv_table.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(taxdf,   file.path(OUT, "02_taxonomy.csv"),  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sel,     file.path(OUT, "03_sample_meta.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(env,     file.path(OUT, "04_env_response.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(grp_df,  file.path(OUT, "05_taxon_group.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.tree(tr,     file.path(OUT, "99_tree.nwk"))

cat("\n===== 完成 =====\n")
cat("01_asv_table.csv     ", nrow(asv_out), "x", ncol(asv_out), "\n")
cat("02_taxonomy.csv      ", nrow(taxdf), "x", ncol(taxdf), "\n")
cat("03_sample_meta.csv   ", nrow(sel), "x", ncol(sel), "\n")
cat("04_env_response.csv  ", nrow(env), "x", ncol(env), "\n")
cat("05_taxon_group.csv   ", nrow(grp_df), "x", ncol(grp_df),
    "（", length(unique(grp_df$group)), " 个分组）\n", sep = "")
cat("99_tree.nwk          Ntip =", Ntip(tr), "\n")
cat("\n因子水平:\n"); print(table(sel$Zone, sel$Vegetation))
cat("\nEF 范围:", sprintf("%.2f - %.2f", min(env$EF), max(env$EF)), "\n")
