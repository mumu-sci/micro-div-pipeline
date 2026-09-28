# =============================================================================
# R/shared/00_setup.R —— 基础设施：配置解析 / 包管理 / 日志 / 缓存 / 种子
# -----------------------------------------------------------------------------
# 本文件是**库**，不含任何分析逻辑。所有模块通过 md_*(CFG, ...) 调用它。
#
# 设计要点（相对 nifH 论文原代码的改进）：
#  1. 显式 CFG 对象，不再靠 load(RData) 往全局环境注入几十个对象
#  2. source 期**无副作用**：不 setwd、不 theme_set；一切通过 md_init(CFG) 显式初始化
#  3. 所有 IO 走绝对路径 md_path(CFG, ...)，不再依赖工作目录
#  4. 缺包 -> 打印安装命令 + 登记跳过，模块不报错退出
#  5. 缓存键 = 数据指纹 + 参数签名；**绘图样式参数不进键**（改配色不触发重算）
#  6. seed 全包唯一来源 CFG$opts$seed
# =============================================================================

MD_VERSION <- "0.1.0"

# ---- 包清单：组名 -> 包名。模块用 md_need(CFG, "ecology") 声明依赖 ----
MD_PKGS <- list(
  core      = c("dplyr", "tidyr", "tibble", "purrr", "stringr", "ggplot2", "patchwork"),
  io        = c("openxlsx"),
  digest    = c("digest"),
  phyloseq  = c("phyloseq"),
  ecology   = c("vegan", "ape", "picante", "betapart"),
  unifrac   = c("GUniFrac", "phangorn"),
  partition = c("iNEXT", "iNEXT.3D"),
  stats     = c("dunn.test", "rstatix", "multcompView", "FSA"),
  model     = c("lme4", "lmerTest", "performance"),
  gam       = c("mgcv"),
  sem       = c("piecewiseSEM"),
  mantel    = c("ecodist", "geosphere"),
  spatial   = c("spdep"),
  pca       = c("FactoMineR", "factoextra"),
  figure    = c("cowplot", "scales", "RColorBrewer", "ggrepel"),
  venn      = c("ggforce"),
  tree      = c("ggtree", "ggtreeExtra", "ggnewscale", "treeio", "tidytree"),
  table     = c("flextable", "officer"),
  parallel  = c("future", "future.apply"),
  map       = c("sf", "terra", "rnaturalearth", "ggspatial"),
  # ---- 以下为 2026 毕业论文补充轮新增（对应模块 12 / 44 / 45 / 80 / 81 / 82 / 83）----
  # 共现网络（80）：默认口径 = 相关矩阵 + 阈值建图
  network   = c("igraph", "tidygraph"),
  # 共现网络备选口径（network$method = "sparcc"）；**默认不装**，用户显式选择时才需要
  network_alt = c("SpiecEasi"),
  # 指示种 IndVal（10_composition$indicator / 81_taxon_group）
  indicator = c("indicspecies", "permute"),
  # ASV/类群共享 UpSet（82）
  upset     = c("ComplexHeatmap", "circlize", "ggplotify"),
  # 差异类群热图（83）；heatmap 组是备选引擎（heatmap_engine = "pheatmap"）
  taxon_heat = c("pheatmap"),
  # 随机森林重要性（40_env_summary$rf / 71_response_linkage 的 "rf" 口径）
  rf        = c("ranger", "randomForest", "party", "pdp"),
  # 节点周转 Sankey（80_network_topology 的 D 面板）
  alluvial  = c("ggalluvial", "ggplotify")
)

# ---- 运行时状态（不放进 CFG，避免污染缓存签名）----
.md_state <- new.env(parent = emptyenv())
.md_state$log_file  <- NULL
.md_state$skipped   <- list()
.md_state$missing   <- list()
.md_state$warned    <- character(0)
.md_state$install_tried <- character(0)   # 自动装包的"试过一次"名单
.md_state$initialized <- FALSE

# ---------------------------------------------------------------------------
# 初始化
# ---------------------------------------------------------------------------
#' 初始化一个分析项目。返回补全后的 CFG。
#' @param CFG      config 文件里定义的列表
#' @param project_dir 项目根目录；NULL 时用 CFG$project$dir
md_init <- function(CFG, project_dir = NULL) {
  stopifnot(is.list(CFG))

  # --- 1. 项目根目录 ---
  root <- project_dir %||% CFG$project$dir %||% "."
  if (identical(root, ".")) {
    fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(fa)) root <- dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE))
    if (basename(root) %in% c("R", "run")) root <- dirname(root)
  }
  root <- normalizePath(root, mustWork = FALSE)
  CFG$project$dir <- root

  # --- 2. 目录 ---
  CFG$paths$abs <- lapply(CFG$paths, function(p) {
    a <- if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(root, p)
    dir.create(a, showWarnings = FALSE, recursive = TRUE); a
  })
  CFG$paths$abs$project <- root

  # --- 3. locale（Rscript --vanilla 下中文路径必需）---
  suppressWarnings(Sys.setlocale("LC_CTYPE", ""))

  # --- 3b. BLAS 线程数 ---
  # 本机实测：runner 进程内加载 sf / terra 时，多线程 OpenBLAS 会
  # `OpenBLAS error: Memory allocation still failed after 10 retries` 直接 abort（退出码 134）。
  # 单线程 BLAS 稳定；本包的分析瓶颈在置换与 R 循环，不在 BLAS，几乎无性能损失。
  # 需要多线程 BLAS 的用户可在调用 md_init() 之前自己设好环境变量。
  if (!nzchar(Sys.getenv("OMP_NUM_THREADS")))      Sys.setenv(OMP_NUM_THREADS = "1")
  if (!nzchar(Sys.getenv("OPENBLAS_NUM_THREADS"))) Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  if (!nzchar(Sys.getenv("MKL_NUM_THREADS")))      Sys.setenv(MKL_NUM_THREADS = "1")

  # --- 4. 核心包（缺了就跑不动）---
  core <- MD_PKGS$core
  miss_core <- core[!vapply(core, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_core)) {
    stop("核心包缺失，无法继续：", paste(miss_core, collapse = ", "),
         "\n请先运行：install.packages(c(",
         paste0('"', miss_core, '"', collapse = ", "), "))", call. = FALSE)
  }
  suppressPackageStartupMessages({
    library(ggplot2); library(patchwork)
  })

  # --- 5. 日志 ---
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  CFG$opts$run_id <- ts
  .md_state$log_file <- file.path(CFG$paths$abs$logs, paste0(CFG$project$label, "_", ts, ".log"))
  .md_state$skipped <- list(); .md_state$missing <- list(); .md_state$warned <- character(0)
  .md_state$install_tried <- character(0)

  # --- 6. 种子（全包唯一来源）---
  set.seed(CFG$opts$seed)

  # --- 6b. 绘图样式选项挂到运行时状态，供 md_box_cfg() 等读取 ---
  .md_state$box <- CFG$plot$box %||% list()

  .md_state$initialized <- TRUE
  CFG$opts$r_version <- paste(R.version$major, R.version$minor, sep = ".")
  CFG$opts$md_version <- MD_VERSION

  md_log(CFG, "INFO", "===== micro-div-pipeline ", MD_VERSION, " =====")
  md_log(CFG, "INFO", "项目: ", CFG$project$label, " @ ", root)
  md_log(CFG, "INFO", "seed: ", CFG$opts$seed,
         " | 置换: ", CFG$assembly$permutations %||% CFG$beta_tests$permutations,
         " | 核数: ", CFG$opts$ncores)
  md_log(CFG, "INFO", "R: ", CFG$opts$r_version, " | 工作目录: ", getwd())

  # 响应变量与树的一次性硬门检查
  if (isTRUE(CFG$phylo$enabled) && is.null(CFG$data$tree))
    md_log(CFG, "WARN", "phylo$enabled=TRUE 但 data$tree 为 NULL —— 系统发育模块将降级跳过")
  CFG
}

# ---------------------------------------------------------------------------
# 日志
# ---------------------------------------------------------------------------
md_log <- function(CFG, level = "INFO", ...) {
  msg <- paste0(vapply(list(...), function(x) paste(as.character(x), collapse = " "),
                       character(1)), collapse = "")
  line <- sprintf("[%s] [%-5s] %s", format(Sys.time(), "%H:%M:%S"), level, msg)
  cat(line, "\n", sep = "")
  f <- .md_state$log_file
  if (!is.null(f)) {
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    cat(line, "\n", sep = "", file = f, append = TRUE)
  }
  invisible(line)
}

md_warn <- function(CFG, ...) {
  key <- paste0(...)
  .md_state$warned <- unique(c(.md_state$warned, key))
  md_log(CFG, "WARN", key)
}

#' 登记模块跳过（不报错），供 runner 汇总
md_skip <- function(CFG, module, reason) {
  .md_state$skipped[[module]] <- reason
  md_log(CFG, "SKIP", module, " —— ", reason)
  invisible(FALSE)
}

md_skipped <- function() .md_state$skipped
md_warnings <- function() .md_state$warned

# ---------------------------------------------------------------------------
# 包管理
# ---------------------------------------------------------------------------
#' 把「组名或包名」的混合向量展开成真实包名
md_expand_pkgs <- function(what) {
  unique(unlist(lapply(what, function(w)
    if (is.character(w) && length(w) == 1 && w %in% names(MD_PKGS)) MD_PKGS[[w]] else w),
    use.names = FALSE))
}

#' 检查一组包是否可用。缺包 -> **先自动安装**，装完仍缺才登记跳过并返回 FALSE。
#'
#' 自动安装由 `CFG$opts$auto_install` 控制（默认 TRUE）。CRAN 源可用
#' `CFG$opts$repos` 覆盖。同一个包在一次运行里只尝试安装一次（避免反复卡住）。
#' 装包失败 **不** 让整个运行报错 —— 该模块 SKIP，其余模块继续。
#' @param what MD_PKGS 的组名，或包名字符串，或它们的混合向量
md_need <- function(CFG, what, module = NULL) {
  pkgs <- md_expand_pkgs(what)
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (all(ok)) return(invisible(TRUE))
  miss <- pkgs[!ok]

  auto <- !identical(CFG$opts$auto_install, FALSE)
  tried <- .md_state$install_tried %||% character(0)
  todo <- setdiff(miss, tried)
  if (auto && length(todo)) {
    .md_state$install_tried <- unique(c(tried, todo))
    repos <- CFG$opts$repos %||% "https://cloud.r-project.org"
    md_log(CFG, "INFO", "缺少包: ", paste(todo, collapse = ", "),
           " —— 尝试自动安装（repos = ", repos, "）")
    ok_ins <- tryCatch({
      utils::install.packages(todo, repos = repos, quiet = TRUE)
      TRUE
    }, error = function(e) {
      md_log(CFG, "WARN", "自动安装失败: ", conditionMessage(e)); FALSE
    })
    if (ok_ins) {
      still <- todo[!vapply(todo, requireNamespace, logical(1), quietly = TRUE)]
      if (!length(still)) {
        md_log(CFG, "OK", "已自动安装: ", paste(todo, collapse = ", "))
      } else {
        md_log(CFG, "WARN", "自动安装后仍不可用: ", paste(still, collapse = ", "),
               "（可能不在 CRAN，例如 SpiecEasi 需从 GitHub 安装）")
      }
      miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
      if (!length(miss)) return(invisible(TRUE))
    }
  }

  .md_state$missing[[module %||% paste(what, collapse = "+")]] <- miss
  md_log(CFG, "WARN", "缺少包: ", paste(miss, collapse = ", "),
         " | 手工安装: install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))")
  invisible(FALSE)
}

md_missing_pkgs <- function() .md_state$missing

#' 载入一组包（已确认存在时）
md_load <- function(what) {
  for (p in md_expand_pkgs(what)) suppressWarnings(suppressPackageStartupMessages(
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 路径
# ---------------------------------------------------------------------------
md_path <- function(CFG, which, ..., must_exist = FALSE) {
  base <- CFG$paths$abs[[which]]
  if (is.null(base)) stop("未知路径键: ", which, call. = FALSE)
  p <- do.call(file.path, c(list(base), list(...)))
  if (must_exist && !file.exists(p)) stop("文件不存在: ", p, call. = FALSE)
  p
}

#' 输入数据文件的绝对路径
md_data_path <- function(CFG, key) {
  p <- CFG$data[[key]]
  if (is.null(p)) return(NULL)
  if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(CFG$project$dir, p)
}

# ---------------------------------------------------------------------------
# 缓存
# ---------------------------------------------------------------------------
#' 参与缓存签名的配置：**排除绘图样式与运行开关**（改配色不触发重算）
md_cfg_for_hash <- function(CFG) {
  x <- CFG
  x$figure <- NULL            # 样式不进键
  x$table  <- NULL
  x$opts$force_recalc <- NULL
  x$opts$figures_only <- NULL
  x$opts$run_id <- NULL
  x$opts$ncores <- NULL       # 核数不影响结果
  x$opts$auto_install <- NULL # 运行开关，不影响结果
  x$opts$repos <- NULL
  x$map    <- NULL            # 地图样式
  x$sem$layers <- NULL
  x$sem$labels <- NULL
  # 2026 补充轮：新增模块一律把**样式**收进各自的 `style` 子列表，
  # 这里统一剔除 —— 于是改配色/字号/布局只重出图，不触发重算。
  # 布局（network$layout）也在 style 里：它只影响节点坐标，不影响任何拓扑数值。
  for (blk in c("network", "taxon_group", "taxon_share", "taxon_diff",
                "niche", "constrained_ord", "response_curves")) {
    if (!is.null(x[[blk]])) x[[blk]]$style <- NULL
  }
  if (!is.null(x$composition$indicator)) x$composition$indicator$style <- NULL
  if (!is.null(x$env$rf))                x$env$rf$style <- NULL
  x
}

md_signature <- function(CFG, deps = NULL) {
  fp <- tryCatch(md_data_md5(CFG), error = function(e) "no-fingerprint")
  digest::digest(list(cfg = md_cfg_for_hash(CFG), data = fp, deps = deps),
                 algo = "xxhash64")
}

#' 模块级缓存。key 建议写 "模块号_用途"，deps 写上游产物的签名或值。
#' figures_only=TRUE 时隐含走缓存；缓存缺失则照常计算。
md_cache <- function(CFG, key, expr, deps = NULL, env = parent.frame()) {
  f <- md_path(CFG, "cache", paste0(key, ".rds"))
  sig <- md_signature(CFG, deps)
  force_calc <- isTRUE(CFG$opts$force_recalc)
  if (!force_calc && file.exists(f)) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(obj) && identical(attr(obj, "md_sig"), sig)) {
      md_log(CFG, "CACHE", key)
      return(obj)
    }
    if (isTRUE(CFG$opts$figures_only)) {
      md_log(CFG, "WARN", "缓存 ", key, " 签名不匹配但 figures_only=TRUE，仍使用旧缓存")
      if (!is.null(obj)) return(obj)
    }
  }
  md_log(CFG, "CALC", key)
  v <- eval(substitute(expr), env)
  attr(v, "md_sig") <- sig
  attr(v, "md_time") <- Sys.time()
  saveRDS(v, f)
  v
}

#' 供 deps 使用的轻量签名（对任意对象）
md_dep <- function(...) digest::digest(list(...), algo = "xxhash64")

#' **分段缓存**：把一段"分析"产生的**全部局部变量**缓存起来。
#'
#' 命中缓存时把这些变量直接灌回调用环境，从而**整段分析被跳过**；
#' 图与结果表在该段之外，所以**每次都会用新样式重画/重写**。
#' 这样「改配色/字号后重出图」不会重跑任何分析。
#'
#' 用法（在模块函数体内，把分析段包起来即可）：
#' ```r
#' md_stage(CFG, "11_analysis", {
#'   idx_m <- ...
#'   tab_test <- ...
#' })
#' # 下面照常画图 / 写表，用的还是 idx_m / tab_test
#' ```
#' 注意：
#'  - 缓存键 = 数据指纹 + 分析相关配置（**绘图样式不进键**），与 md_cache 一致
#'  - `opts$force_recalc = TRUE` 时强制重算
#'  - 该段内**新建**的变量才会被缓存；修改已有变量不会被记住（如需记住请新建变量）
md_stage <- function(CFG, key, expr, env = parent.frame(), only = NULL) {
  f <- md_path(CFG, "cache", paste0(key, ".rds"))
  sig <- md_signature(CFG)
  if (!isTRUE(CFG$opts$force_recalc) && file.exists(f)) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(obj) && identical(attr(obj, "md_sig"), sig)) {
      md_log(CFG, "CACHE", key, " —— 跳过该段分析")
      for (nm in names(obj$vars)) assign(nm, obj$vars[[nm]], envir = env)
      return(invisible(TRUE))
    }
  }
  md_log(CFG, "CALC", key)
  before <- ls(env, all.names = TRUE)
  eval(substitute(expr), envir = env)
  after <- ls(env, all.names = TRUE)
  new <- setdiff(after, before)
  if (!is.null(only)) new <- intersect(new, only)
  vars <- if (length(new)) mget(new, envir = env, ifnotfound = list(NULL)) else list()
  out <- list(vars = vars, key = key, n = length(vars))
  attr(out, "md_sig") <- sig
  attr(out, "md_time") <- Sys.time()
  saveRDS(out, f)
  invisible(FALSE)
}

md_recalc <- function(CFG) isTRUE(CFG$opts$force_recalc)
md_figonly <- function(CFG) isTRUE(CFG$opts$figures_only)

# ---------------------------------------------------------------------------
# 种子
# ---------------------------------------------------------------------------
#' 从全包唯一种子派生子种子。模块里凡有随机过程都必须先调用它。
md_seed <- function(CFG, offset = 0L) {
  s <- as.integer(CFG$opts$seed) + as.integer(offset)
  set.seed(s)
  invisible(s)
}

# ---------------------------------------------------------------------------
# 并行
# ---------------------------------------------------------------------------
#' 设置并行后端（future）。返回实际使用的核数。
md_parallel <- function(CFG, ncores = NULL, seed = NULL) {
  n <- ncores %||% CFG$opts$ncores %||% 1L
  n <- max(1L, min(as.integer(n), parallel::detectCores(logical = FALSE)))
  if (n <= 1L || !requireNamespace("future", quietly = TRUE)) return(1L)
  future::plan(future::multisession, workers = n)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  if (!is.null(seed)) set.seed(seed)
  md_log(CFG, "INFO", "并行后端: multisession x", n)
  n
}

md_parallel_off <- function() {
  if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 会话信息（可复现性）
# ---------------------------------------------------------------------------
md_session_info <- function(CFG) {
  ip <- utils::installed.packages()[, "Version"]
  used <- unique(unlist(MD_PKGS))
  used <- intersect(used, names(ip))
  data.frame(
    item = c("md_version", "R_version", "platform", "seed", "run_id",
             "project", "project_dir", paste0("pkg:", used)),
    value = c(MD_VERSION, paste(R.version$major, R.version$minor, sep = "."),
              R.version$platform, as.character(CFG$opts$seed), CFG$opts$run_id,
              CFG$project$label, CFG$project$dir, as.character(ip[used])),
    stringsAsFactors = FALSE
  )
}

#' 把会话信息与数据指纹写进 results/
md_write_session <- function(CFG) {
  si <- md_session_info(CFG)
  fp <- tryCatch(md_data_md5(CFG), error = function(e) NULL)
  if (!is.null(fp)) {
    si <- rbind(si, data.frame(item = paste0("md5:", names(fp)),
                               value = as.character(fp), stringsAsFactors = FALSE))
  }
  p <- md_path(CFG, "results", "00_session_info.csv")
  utils::write.csv(si, p, row.names = FALSE, fileEncoding = "UTF-8")
  p
}

# ---------------------------------------------------------------------------
# 模块返回契约
# ---------------------------------------------------------------------------
#' 每个分析模块必须返回本函数的产物。status: OK | SKIP | PARTIAL | EMPTY
md_result <- function(module, status = "OK", reason = NULL,
                      tables = character(0), figures = character(0),
                      notes = NULL) {
  list(module = module, status = status, reason = reason,
       tables = tables, figures = figures, notes = notes)
}

#' 模块统一开头：检查响应变量 / 系统发育树等前置条件，返回 TRUE 表示可以继续
md_precheck <- function(CFG, module, need_response = FALSE, need_phylo = FALSE) {
  if (need_response && !isTRUE(CFG$response$enabled)) {
    md_skip(CFG, module, "response$enabled = FALSE（本论文没有响应变量）")
    return(FALSE)
  }
  if (need_phylo) {
    if (!isTRUE(CFG$phylo$enabled)) {
      md_skip(CFG, module, "phylo$enabled = FALSE")
      return(FALSE)
    }
    if (is.null(CFG$data$tree)) {
      md_skip(CFG, module, "未提供系统发育树（data$tree = NULL）")
      return(FALSE)
    }
    p <- md_data_path(CFG, "tree")
    if (is.null(p) || !file.exists(p)) {
      md_skip(CFG, module, paste0("树文件不存在: ", p %||% "NULL"))
      return(FALSE)
    }
  }
  TRUE
}
