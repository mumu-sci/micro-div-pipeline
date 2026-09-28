# =============================================================================
# R/90_run_all.R —— 统一 runner（由 run/manifest.R 驱动）
# -----------------------------------------------------------------------------
# 用法：
#   Rscript --vanilla R/90_run_all.R                                  # 用 config/config_demo.R
#   Rscript --vanilla R/90_run_all.R --config=config/config.R          # 指定配置
#   Rscript --vanilla R/90_run_all.R --project=D:/path/to/paper        # 指定论文项目目录
#   Rscript --vanilla R/90_run_all.R --only=01,05,10                   # 只跑指定模块
#   Rscript --vanilla R/90_run_all.R --figures-only                    # 只重出图（走缓存）
#   Rscript --vanilla R/90_run_all.R --force                           # 忽略缓存全部重算
#
# 行为约定：
#   - 任一模块 ERROR -> 继续跑其余模块，最后以非零退出码结束（**不静默通过**）
#   - warning 不截断 source（用 withCallingHandlers + muffleWarning，
#     不能用 tryCatch(warning=)，那会 unwind 导致脚本只跑到第一条 warning）
#   - 缺文件 / 缺包 -> SKIP，不算 ERROR
#   - 每步写日志到 logs/<label>_<run_id>.log
# =============================================================================

# ---- 参数 ----
.args <- commandArgs(TRUE)
.get_arg <- function(name, default = NULL) {
  hits <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hits)) return(sub(paste0("^--", name, "="), "", hits[1]))
  # 兜底：PowerShell 把 "a,b,c" 拆成多参数时，把裸值重新粘回来
  i <- which(.args == paste0("--", name))
  if (length(i) && i[1] < length(.args)) {
    tail_args <- .args[(i[1] + 1):length(.args)]
    tail_args <- tail_args[!grepl("^--", tail_args)]
    if (length(tail_args)) return(paste(tail_args, collapse = ","))
  }
  default
}
.has_flag <- function(name) any(.args == paste0("--", name))

cfg_arg     <- .get_arg("config", "config/config_demo.R")
proj_arg    <- .get_arg("project", NULL)
only_arg    <- .get_arg("only", NULL)
fig_only    <- .has_flag("figures-only")
force_calc  <- .has_flag("force")
# --check / --dry-run：只读自检。**保持为局部变量，不写进 CFG$opts**，
# 因此完全不用动 md_cfg_for_hash()，缓存键不受影响（见文件末的 check 分支）。
check_only  <- .has_flag("check") || .has_flag("dry-run")

# ---- 定位包根 ----
.self <- tryCatch({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE) else NA_character_
}, error = function(e) NA_character_)
PKG_DIR <- if (!is.na(.self)) dirname(dirname(.self)) else normalizePath(getwd())
if (!dir.exists(file.path(PKG_DIR, "R", "shared"))) {
  if (dir.exists(file.path(getwd(), "R", "shared"))) PKG_DIR <- getwd()
}
PROJECT_DIR <- proj_arg %||% PKG_DIR

suppressWarnings(Sys.setlocale("LC_CTYPE", ""))

# ---- 载入共享库（顺序固定）----
for (f in c("00_setup.R", "10_io.R", "11_stats.R", "12_plot_theme.R", "13_plot_types.R")) {
  p <- file.path(PKG_DIR, "R", "shared", f)
  if (!file.exists(p)) stop("共享库缺失: ", p, call. = FALSE)
  source(p, local = FALSE)
}
source(file.path(PKG_DIR, "run", "manifest.R"), local = FALSE)

# ---- 载入配置 ----
cfg_path <- if (grepl("^(/|[A-Za-z]:)", cfg_arg)) cfg_arg else file.path(PROJECT_DIR, cfg_arg)
if (!file.exists(cfg_path)) stop("配置文件不存在: ", cfg_path, call. = FALSE)
# 先 parse 再 source：语法错时给一句干净的中文报错（附行号），
# 而不是 R 的裸 "unexpected symbol" —— 非 --check 运行同样受益。
.cfg_syn <- tryCatch({ parse(cfg_path); NA_character_ },
                     error = function(e) conditionMessage(e))
if (!is.na(.cfg_syn))
  stop("配置文件语法错误: ", cfg_path, "\n  ", gsub("[\r\n]+", " ", .cfg_syn),
       call. = FALSE)
CFG <- NULL
source(cfg_path, local = FALSE)
if (!exists("CFG") || !is.list(CFG)) stop("配置文件未定义 CFG 列表: ", cfg_path, call. = FALSE)

if (fig_only)   CFG$opts$figures_only <- TRUE
if (force_calc) CFG$opts$force_recalc <- TRUE

CFG <- md_init(CFG, project_dir = PROJECT_DIR)
md_log(CFG, "INFO", "配置: ", cfg_path)

# ---- 执行顺序 ----
order_ids <- md_manifest_order()
if (!is.null(only_arg)) {
  # 分隔符接受 , ; + 空格（PowerShell 会把未加引号的逗号拆成多个参数，故一并兜住）
  want <- trimws(unlist(strsplit(only_arg, "[,;+[:space:]]+")))
  want <- want[nzchar(want)]
  order_ids <- order_ids[order_ids %in% want]
  md_log(CFG, "INFO", "只跑: ", paste(order_ids, collapse = ", "))
}
if (fig_only) md_log(CFG, "INFO", "figures_only = TRUE（跳过 heavy 模块，其余走缓存）")

# =============================================================================
# --check / --dry-run：只读自检（SKILL.md 第 6 步「自检（不执行）」的可执行版本）
# -----------------------------------------------------------------------------
# 约定：
#   * **不装包**（强制 auto_install = FALSE，md_need() 会遵守）
#   * **不写 results/**、**不跑任何分析**、**不碰 cache/**
#   * 只落一份 logs/preflight_<run_id>.csv
#   * 语法错 / 数据不可读 / 必需列缺失 -> BLOCK -> 退出码 1
#     缺包 -> WARN -> 退出码 0（缺包只让对应模块 SKIP，不算失败）
# =============================================================================
if (check_only) {
  CFG$opts$auto_install <- FALSE
  .pf <- list()
  .pf_add <- function(section, item, status, detail) {
    .pf[[length(.pf) + 1L]] <<- data.frame(
      section = section, item = item, status = status, detail = detail,
      stringsAsFactors = FALSE)
  }
  .pf_add("meta", "config", "OK", cfg_path)
  .pf_add("meta", "project_dir", "OK", PROJECT_DIR)
  .pf_add("meta", "modules", "OK", paste(order_ids, collapse = ", "))

  # ---- 1. 语法：R/（含 R/shared/）+ run/ ----
  #   本配置不进这个扫描：它在上面已被 parse 守卫过（语法错会直接停在那一步），
  #   放进来只会永远返回 OK，属于误导性的死代码。
  syn <- unique(c(list.files(file.path(PKG_DIR, "R"), "[.]R$",
                             recursive = TRUE, full.names = TRUE),
                  list.files(file.path(PKG_DIR, "run"), "[.]R$", full.names = TRUE)))
  n_bad_syn <- 0L
  for (f in syn) {
    msg <- tryCatch({ parse(f); NA_character_ },
                    error = function(e) conditionMessage(e))
    if (!is.na(msg)) {
      n_bad_syn <- n_bad_syn + 1L
      .pf_add("syntax", basename(f), "BLOCK", gsub("[\r\n]+", " ", msg))
    }
  }
  .pf_add("syntax", sprintf("%d 个 R 文件", length(syn)),
          if (n_bad_syn) "BLOCK" else "OK",
          sprintf("%d 个解析失败", n_bad_syn))

  # ---- 2. 依赖：只查不装 ----
  dep <- do.call(rbind, lapply(order_ids, function(mid) {
    m <- md_manifest_get(mid)
    if (is.null(m$needs) || !length(m$needs)) return(NULL)
    ok <- isTRUE(tryCatch(md_need(CFG, m$needs, mid), error = function(e) FALSE))
    data.frame(id = mid, needs = paste(m$needs, collapse = ", "), ok = ok,
               stringsAsFactors = FALSE)
  }))
  miss <- md_missing_pkgs()
  n_miss_mod <- 0L
  if (!is.null(dep) && nrow(dep)) {
    for (i in seq_len(nrow(dep))) {
      if (isTRUE(dep$ok[i])) next
      n_miss_mod <- n_miss_mod + 1L
      .pf_add("deps", dep$id[i], "WARN",
              paste0("缺 ", paste(miss[[dep$id[i]]] %||% dep$needs[i], collapse = ", "),
                     " —— 该模块将 SKIP（装上即可）"))
    }
    .pf_add("deps", sprintf("%d 个模块声明了包依赖", nrow(dep)),
            if (n_miss_mod) "WARN" else "OK",
            sprintf("缺包模块 %d 个", n_miss_mod))
  } else {
    .pf_add("deps", "所有模块", "OK", "manifest 里没有声明包依赖")
  }

  # ---- 3. 列名 ----
  #   分两档，避免误判：
  #   A) **严格列**（BLOCK）：必须真的存在于 03_sample_meta.csv / 04_env_response.csv 里
  #      —— design$factors / design$ordered_levels / design$site_var / env$vars /
  #         env$groups 的**成员变量**（注意：env$groups 的 list **名字**是分组标签，不是列名）
  #   B) **可派生名**（WARN）：响应变量 / SEM 边表 / model$focus·covariates·derived。
  #      它们可能由上游模块在运行时生成（Soil_PC1 来自 40、Shannon 来自 α 多样性、
  #      PCoA1 来自 13），所以解析链放宽为：
  #        原始列 ∪ md_core_analysis() 列 ∪ md_derived() 名 ∪ α 指数名 ∪ PCoA 轴名
  #      对不上只提示"可能是拼写错误，也可能要先把上游模块跑一遍"。
  inp <- tryCatch(md_read_inputs(CFG, quiet = TRUE), error = function(e) e)
  if (inherits(inp, "error")) {
    .pf_add("columns", "五文件", "BLOCK",
            paste0("读不进来：", gsub("[\r\n]+", " ", conditionMessage(inp))))
  } else {
    na_meta <- names(inp$meta); na_env <- names(inp$env); na_tax <- names(inp$tax)
    raw_cols <- unique(c(na_meta, na_env, na_tax))

    # ---- A) 严格列 ----
    strict <- list()
    .push_s <- function(where, name) {
      if (is.null(name) || !length(name) || !nzchar(as.character(name[1]))) return(invisible())
      strict[[length(strict) + 1L]] <<- list(where = where, name = as.character(name[1]))
    }
    for (f in as.character(CFG$design$factors %||% character(0)))
      .push_s("design$factors", f)
    for (f in names(CFG$design$ordered_levels %||% list()))
      .push_s("design$ordered_levels", f)
    .push_s("design$site_var", CFG$design$site_var %||% "")
    for (v in as.character(CFG$env$vars %||% character(0)))
      .push_s("env$vars", v)
    for (gname in names(CFG$env$groups %||% list()))
      for (v in as.character(CFG$env$groups[[gname]]))
        .push_s(paste0("env$groups$", gname), v)
    n_bad_col <- 0L
    for (z in strict) {
      if (z$name %in% raw_cols) next
      n_bad_col <- n_bad_col + 1L
      .pf_add("columns", paste0(z$where, " = ", z$name), "BLOCK",
              paste0("原始数据里没有这一列（可用列：",
                     paste(utils::head(raw_cols, 30), collapse = ", "),
                     if (length(raw_cols) > 30) " ..." else "", "）"))
    }
    .pf_add("columns", sprintf("严格列 %d 个（设计因子 / 站点列 / 环境变量 / env$groups 成员）",
                               length(strict)),
            if (n_bad_col) "BLOCK" else "OK",
            sprintf("缺失 %d 个；meta %d 列 / env %d 列 / tax %d 列",
                    n_bad_col, length(na_meta), length(na_env), length(na_tax)))

    # ---- B) 可派生名：放宽解析链 ----
    pool <- raw_cols
    ana <- tryCatch(md_core_analysis(CFG), error = function(e) NULL)
    if (!is.null(ana)) pool <- unique(c(pool, names(ana)))
    ec <- tryCatch(md_env_cache(CFG), error = function(e) NULL)
    if (!is.null(ec) && !is.null(ec$derived)) pool <- unique(c(pool, names(ec$derived)))
    pool <- unique(c(pool, "Observed", "Shannon", "Simpson", "Pielou",
                     paste0("PCoA", 1:3),
                     as.character(CFG$model$derived %||% character(0))))

    lv <- list()
    .push_l <- function(where, name) {
      if (is.null(name) || !length(name) || !nzchar(as.character(name[1]))) return(invisible())
      lv[[length(lv) + 1L]] <<- list(where = where, name = as.character(name[1]))
    }
    .push_l("response$var", CFG$response$var %||% "")
    for (v in as.character(CFG$model$focus %||% character(0)))      .push_l("model$focus", v)
    for (v in as.character(CFG$model$covariates %||% character(0))) .push_l("model$covariates", v)
    for (v in as.character(CFG$model$derived %||% character(0)))    .push_l("model$derived", v)
    ed <- CFG$sem$edges
    if (is.data.frame(ed) && all(c("from", "to") %in% names(ed)) && nrow(ed)) {
      for (v in unique(c(as.character(ed$from), as.character(ed$to))))
        .push_l("sem$edges", v)
    }
    n_warn_col <- 0L
    for (z in lv) {
      if (z$name %in% pool) next
      n_warn_col <- n_warn_col + 1L
      .pf_add("derived", paste0(z$where, " = ", z$name), "WARN",
              paste0("在原始列与已知派生名里都没找到 —— 可能是拼写错误；",
                     "也可能它要由上游模块运行时生成（Soil_PC1←40 / Shannon←α / PCoA1←13），",
                     "那种情况请先跑上游模块再 --check"))
    }
    .pf_add("derived", sprintf("可派生名 %d 个（响应变量 / SEM 边表 / model 焦点·协变量·派生）",
                               length(lv)),
            if (n_warn_col) "WARN" else "OK",
            sprintf("对不上 %d 个（当前可解析名 %d 个%s）",
                    n_warn_col, length(pool),
                    if (is.null(ana)) "；分析总表不可用，已退化为原始列 + 已知派生名" else ""))
  }

  # ---- 4. 输入指纹：路径 + 行数 + md5 ----
  md5 <- md_data_md5(CFG)
  for (k in names(md5)) {
    p <- md_data_path(CFG, k)
    nrow_k <- NA_integer_
    if (!is.null(p) && file.exists(p) && !identical(md5[[k]], "MISSING")) {
      nrow_k <- tryCatch(length(readLines(p, warn = FALSE)) - 1L,
                         error = function(e) NA_integer_)
    }
    .pf_add("fingerprint", k,
            if (identical(md5[[k]], "MISSING")) "BLOCK" else "OK",
            sprintf("%s | 行数(不含表头) %s | md5 %s",
                    if (is.null(p)) "(未指定)" else p,
                    if (is.na(nrow_k)) "NA" else as.character(nrow_k), md5[[k]]))
  }
  if (!inherits(inp, "error"))
    .pf_add("fingerprint", "ASV x 样本",
            "OK", sprintf("%d ASV x %d 样本（交集后）",
                          nrow(inp$asv), ncol(inp$asv)))

  # ---- 汇总 ----
  pf <- do.call(rbind, .pf)
  pf_file <- file.path(CFG$paths$abs$logs,
                       paste0("preflight_", CFG$opts$run_id, ".csv"))
  utils::write.csv(pf, pf_file, row.names = FALSE, fileEncoding = "UTF-8")
  n_block <- sum(pf$status == "BLOCK"); n_warn <- sum(pf$status == "WARN")
  print(pf[pf$status != "OK", c("section", "item", "status", "detail")],
        row.names = FALSE)
  cat("\n")
  md_log(CFG, "INFO", "预检 ==============================================")
  md_log(CFG, if (n_block) "ERROR" else "INFO",
         sprintf("BLOCK %d | WARN %d | OK %d | 共 %d 项",
                 n_block, n_warn, sum(pf$status == "OK"), nrow(pf)))
  md_log(CFG, "INFO", "报告: ", pf_file)
  md_log(CFG, "INFO", "（--check 是只读预检：未装包、未跑分析、未写 results/）")
  cat("\n")
  quit(status = if (n_block) 1L else 0L, save = "no")
}

md_write_session(CFG)

# ---- 逐个执行 ----
.summary <- vector("list", length(order_ids))
t_all <- Sys.time()

for (i in seq_along(order_ids)) {
  mid <- order_ids[i]
  m   <- md_manifest_get(mid)
  t0  <- Sys.time()
  md_log(CFG, "INFO", sprintf("[%d/%d] %s %s", i, length(order_ids), mid, m$fun))

  # 1) 文件存在？
  mod_file <- md_module_file(CFG, m)
  if (!file.exists(mod_file)) {
    md_skip(CFG, mid, paste0("模块文件尚未实现: ", m$file))
    .summary[[i]] <- list(id = mid, status = "SKIP", reason = "not implemented",
                          secs = 0, tables = "", figures = "")
    next
  }
  # 2) figures_only：**不再跳过重计算模块** —— 各模块的分析都走 md_cache/md_stage，
  #    只要数据与分析方法没变（绘图样式不进缓存键），缓存就会命中而不重算；
  #    图与结果表照常按新样式重画/重写。
  #    （早期版本在这里直接跳过 heavy 模块，会让"只重出图"反而漏刷图表。）
  if (fig_only && isTRUE(m$heavy)) {
    md_log(CFG, "INFO", mid, "：figures_only —— 分析段走缓存（缓存缺失时才算一次并落缓存）")
  }
  # 3) 包检查（缺则 SKIP，不 ERROR）
  if (!is.null(m$needs) && !md_need(CFG, m$needs, mid)) {
    md_skip(CFG, mid, paste0("缺少依赖包: ", paste(md_missing_pkgs()[[mid]], collapse = ", ")))
    .summary[[i]] <- list(id = mid, status = "SKIP", reason = "missing packages",
                          secs = 0, tables = "", figures = "")
    next
  }

  # 4) 执行
  res <- NULL
  err <- NULL
  t1 <- Sys.time()
  tryCatch({
    source(mod_file, local = FALSE)
    if (!exists(m$fun, mode = "function")) stop("模块未导出函数 ", m$fun)
    res <- withCallingHandlers(
      do.call(m$fun, list(CFG)),
      warning = function(w) {
        # 只记录再 muffle，绝不 unwind（否则 source 会中途截断）
        md_log(CFG, "WARN", mid, " : ", conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  }, error = function(e) { err <<- e })

  secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  if (!is.null(err)) {
    md_log(CFG, "ERROR", mid, " : ", conditionMessage(err))
    .summary[[i]] <- list(id = mid, status = "ERROR", reason = conditionMessage(err),
                          secs = secs, tables = "", figures = "")
  } else {
    st <- res$status %||% "OK"
    md_log(CFG, if (st == "OK") "DONE" else st, mid,
           sprintf(" (%.1f s)", secs),
           if (!is.null(res$reason)) paste0(" —— ", res$reason) else "")
    .summary[[i]] <- list(id = mid, status = st, reason = res$reason %||% "",
                          secs = secs,
                          tables = paste(res$tables %||% character(0), collapse = "; "),
                          figures = paste(res$figures %||% character(0), collapse = "; "))
  }

  # 增量落盘运行状态：95（manifest 里 last=TRUE）必须读到**本轮**的结果，
  # 而不是上一轮的（否则 00_index 的 Run_summary 永远落后一轮）。
  .done <- Filter(Negate(is.null), .summary)
  if (length(.done)) {
    .cur <- do.call(rbind, lapply(.done, function(x)
      data.frame(id = x$id, status = x$status, secs = round(x$secs, 1),
                 reason = x$reason, tables = x$tables, figures = x$figures,
                 stringsAsFactors = FALSE)))
    saveRDS(list(summary = .cur, skipped = md_skipped(), warnings = md_warnings(),
                 missing = md_missing_pkgs(), figures = md_figures_written(),
                 run_id = CFG$opts$run_id),
            file.path(CFG$paths$abs$cache, "00_run_state.rds"))
  }
}

# ---- 汇总 ----
sumdf <- do.call(rbind, lapply(.summary, function(x)
  data.frame(id = x$id, status = x$status, secs = round(x$secs, 1),
             reason = x$reason, tables = x$tables, figures = x$figures,
             stringsAsFactors = FALSE)))
sumdf <- sumdf[order(match(sumdf$id, order_ids)), , drop = FALSE]
rownames(sumdf) <- NULL

total_secs <- as.numeric(difftime(Sys.time(), t_all, units = "secs"))
md_log(CFG, "INFO", "==================================================")
md_log(CFG, "INFO", sprintf("总计 %.1f min | OK %d | SKIP %d | ERROR %d | 图 %d | 表 %d",
                            total_secs / 60,
                            sum(sumdf$status == "OK"), sum(sumdf$status == "SKIP"),
                            sum(sumdf$status == "ERROR"),
                            length(md_figures_written()),
                            length(list.files(CFG$paths$abs$results, "\\.xlsx$"))))
print(sumdf[, c("id", "status", "secs", "reason")], row.names = FALSE)

# 落盘汇总
utils::write.csv(sumdf, file.path(CFG$paths$abs$logs, paste0("summary_", CFG$opts$run_id, ".csv")),
                 row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(summary = sumdf, skipped = md_skipped(), warnings = md_warnings(),
             missing = md_missing_pkgs(), figures = md_figures_written(),
             run_id = CFG$opts$run_id),
        file.path(CFG$paths$abs$cache, "00_run_state.rds"))

if (any(sumdf$status == "ERROR")) {
  md_log(CFG, "ERROR", "有模块失败：", paste(sumdf$id[sumdf$status == "ERROR"], collapse = ", "))
  quit(status = 1L, save = "no")
}
md_log(CFG, "INFO", "全部模块执行完毕（无 ERROR）")
