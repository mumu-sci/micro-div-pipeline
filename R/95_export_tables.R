# =============================================================================
# R/95_export_tables.R —— 产物索引 + 后置校验（收尾模块，永远最后跑）
# -----------------------------------------------------------------------------
# 1. 清点 results/*.xlsx 与 figures/*.pdf，生成 results/00_index.xlsx
# 2. **后置校验**（封堵 nifH 那个 0 字节 PDF 类问题）：
#      - PDF: 文件大小 > 0 且能渲染出页
#      - xlsx: 能按 sheet 名读出，且行列数 > 0
#      - data/: md5 与 01 号体检时记录的一致（源数据未被改写）
# 3. 把 runner 的执行汇总并入索引
# =============================================================================

md_export_tables <- function(CFG) {
  module <- "95"
  if (!md_need(CFG, "io", module)) return(md_result(module, "SKIP", "缺少 openxlsx"))

  res_dir <- CFG$paths$abs$results
  fig_dir <- CFG$paths$abs$figures

  # ---- 1. 表索引 ----
  xs <- list.files(res_dir, "\\.xlsx$", full.names = TRUE)
  xs <- xs[!grepl("^00_index", basename(xs))]
  tab_rows <- list()
  val_rows <- list()
  for (f in xs) {
    sh <- tryCatch(openxlsx::getSheetNames(f), error = function(e) character(0))
    for (s in sh) {
      d <- tryCatch(md_read_xlsx(f, s), error = function(e) NULL)
      tab_rows[[length(tab_rows) + 1L]] <- data.frame(
        file = basename(f), sheet = s,
        n_row = if (is.null(d)) NA_integer_ else nrow(d),
        n_col = if (is.null(d)) NA_integer_ else ncol(d),
        stringsAsFactors = FALSE)
      if (is.null(d) || nrow(d) == 0 || ncol(d) == 0) {
        val_rows[[length(val_rows) + 1L]] <- data.frame(
          artifact = basename(f), type = "xlsx", check = paste0("sheet 可读且非空: ", s),
          status = "FAIL", detail = if (is.null(d)) "按 sheet 名读取失败" else "空表",
          stringsAsFactors = FALSE)
      }
    }
    if (length(sh) == 0) {
      val_rows[[length(val_rows) + 1L]] <- data.frame(
        artifact = basename(f), type = "xlsx", check = "工作簿可读", status = "FAIL",
        detail = "getSheetNames 失败（可能是 openxlsx 悬空 drawing）", stringsAsFactors = FALSE)
    }
  }
  tab_index <- if (length(tab_rows)) do.call(rbind, tab_rows) else
    data.frame(file = character(0), sheet = character(0), n_row = integer(0), n_col = integer(0))

  # ---- 2. 图索引 + 校验 ----
  pdfs <- list.files(fig_dir, "\\.pdf$", full.names = TRUE)
  fig_rows <- list()
  has_pdftools <- requireNamespace("pdftools", quietly = TRUE)
  ascii_rows <- list()
  for (f in pdfs) {
    sz <- file.info(f)$size
    np <- tryCatch({
      if (has_pdftools) length(pdftools::pdf_info(f)$pages) else NA_integer_
    }, error = function(e) NA_integer_)
    ok <- !is.na(sz) && sz > 0 && (is.na(np) || np > 0)
    fig_rows[[length(fig_rows) + 1L]] <- data.frame(
      file = basename(f), size_kb = round(sz / 1024, 1), pages = np,
      status = if (ok) "PASS" else "FAIL", stringsAsFactors = FALSE)
    if (!ok) val_rows[[length(val_rows) + 1L]] <- data.frame(
      artifact = basename(f), type = "pdf", check = "大小 > 0 且可渲染",
      status = "FAIL", detail = paste0("size=", sz, " bytes, pages=", np),
      stringsAsFactors = FALSE)

    # ---- 2b. 图中文字必须纯 ASCII（码点 <= 0x7E）----
    # 这是防复发的护栏：nifH 移植期曾把中文过滤日志当成了图的 x 轴标签。
    if (has_pdftools) {
      txt <- tryCatch(paste(pdftools::pdf_text(f), collapse = ""), error = function(e) "")
      ch <- unique(strsplit(txt, "")[[1]])
      # 只看**可见字符**：换行/制表等空白不算问题
      ch <- ch[!is.na(ch) & grepl("[^ -~]", ch) & !grepl("[[:space:]]", ch)]
      if (length(ch)) {
        labs <- vapply(ch, function(cc) sprintf("%s(U+%04X)", cc, utf8ToInt(cc)), character(1))
        ascii_rows[[length(ascii_rows) + 1L]] <- data.frame(
          file = basename(f), n_nonascii = length(ch),
          chars = paste(labs, collapse = " "), stringsAsFactors = FALSE)
        val_rows[[length(val_rows) + 1L]] <- data.frame(
          artifact = basename(f), type = "pdf", check = "图中文字纯 ASCII",
          status = "FAIL",
          detail = paste0("发现 ", length(ch), " 个非 ASCII 字符: ", paste(labs, collapse = " ")),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!has_pdftools) {
    val_rows[[length(val_rows) + 1L]] <- data.frame(
      artifact = "figures/", type = "pdf", check = "图中文字纯 ASCII",
      status = "WARN", detail = "未安装 pdftools，跳过该审计（install.packages(\"pdftools\")）",
      stringsAsFactors = FALSE)
  }
  fig_ascii <- if (length(ascii_rows)) do.call(rbind, ascii_rows) else
    data.frame(file = character(0), n_nonascii = integer(0), chars = character(0))

  # ---- 2c. 残留图告警：figures/ 里有本轮没写的 pdf ----
  # 典型情形：某模块本轮 SKIP（例如无树时的 21），它上一轮出的图还留在目录里，容易被误当成新图。
  run_state0 <- tryCatch(readRDS(md_path(CFG, "cache", "00_run_state.rds")),
                         error = function(e) NULL)
  written <- basename(run_state0$figures %||% character(0))
  stale <- if (length(pdfs)) setdiff(basename(pdfs), written) else character(0)
  if (length(stale)) {
    val_rows[[length(val_rows) + 1L]] <- data.frame(
      artifact = "figures/", type = "pdf", check = "无上一轮残留图",
      status = "WARN",
      detail = paste0(length(stale), " 张图本轮未写出（可能是该模块 SKIP 留下的旧图，",
                      "勿当作本轮结果）：", paste(stale, collapse = ", ")),
      stringsAsFactors = FALSE)
  }
  fig_stale <- data.frame(file = stale,
                          note = rep("本轮未写出（可能是上一轮残留）", length(stale)),
                          stringsAsFactors = FALSE)
  fig_index <- if (length(fig_rows)) do.call(rbind, fig_rows) else
    data.frame(file = character(0), size_kb = numeric(0), pages = integer(0), status = character(0))

  # ---- 3. 数据护栏：源数据 md5 是否与体检时一致 ----
  now_fp <- tryCatch(md_data_md5(CFG), error = function(e) NULL)
  old <- tryCatch(readRDS(md_path(CFG, "cache", "01_health.rds")), error = function(e) NULL)
  if (!is.null(now_fp) && !is.null(old$fingerprint)) {
    same <- identical(now_fp, old$fingerprint)
    val_rows[[length(val_rows) + 1L]] <- data.frame(
      artifact = "data/", type = "guard", check = "源数据 md5 未被改写",
      status = if (same) "PASS" else "FAIL",
      detail = if (same) "与 01 号体检时一致" else "md5 发生变化！",
      stringsAsFactors = FALSE)
  }
  # 兜底：results/ 与 figures/ 里不得出现任何指向 data/ 的写入痕迹
  val_rows[[length(val_rows) + 1L]] <- data.frame(
    artifact = "results/", type = "guard", check = "结果表不在 data/ 目录下",
    status = if (startsWith(normalizePath(res_dir, mustWork = FALSE),
                            normalizePath(CFG$paths$abs$data, mustWork = FALSE))) "FAIL" else "PASS",
    detail = res_dir, stringsAsFactors = FALSE)

  val_df <- if (length(val_rows)) do.call(rbind, val_rows) else
    data.frame(artifact = character(0), type = character(0), check = character(0),
               status = character(0), detail = character(0))
  n_fail <- sum(val_df$status == "FAIL")

  # ---- 4. runner 汇总 ----
  run_state <- tryCatch(readRDS(md_path(CFG, "cache", "00_run_state.rds")),
                        error = function(e) NULL)
  run_df <- if (!is.null(run_state$summary)) run_state$summary else
    data.frame(id = character(0), status = character(0), secs = numeric(0),
               reason = character(0), tables = character(0), figures = character(0))
  if (!is.null(run_state$missing) && length(run_state$missing)) {
    miss_df <- do.call(rbind, lapply(names(run_state$missing), function(k)
      data.frame(module = k, missing_pkgs = paste(run_state$missing[[k]], collapse = ", "),
                 stringsAsFactors = FALSE)))
  } else miss_df <- data.frame(module = character(0), missing_pkgs = character(0))

  tp <- md_write_xlsx(CFG, "00_index", list(
    Index_tables = tab_index,
    Index_figures = fig_index,
    Figure_text_ascii = fig_ascii,
    Figures_stale = fig_stale,
    Validation = val_df,
    Run_summary = run_df,
    Skipped = if (length(run_state$skipped))
      data.frame(module = names(run_state$skipped),
                 reason = unlist(run_state$skipped), stringsAsFactors = FALSE) else
      data.frame(module = character(0), reason = character(0)),
    Missing_packages = miss_df,
    Warnings = data.frame(warning = run_state$warnings %||% character(0), stringsAsFactors = FALSE),
    Session = md_session_info(CFG)))

  if (n_fail > 0) {
    md_log(CFG, "ERROR", "产物校验发现 ", n_fail, " 项 FAIL：")
    print(val_df[val_df$status == "FAIL", , drop = FALSE], row.names = FALSE)
    return(md_result(module, "PARTIAL", reason = paste0(n_fail, " 项产物校验 FAIL"),
                     tables = basename(tp)))
  }
  md_log(CFG, "INFO", "产物校验通过：", nrow(tab_index), " 个 sheet / ",
         nrow(fig_index), " 张图，全部可读且非空")
  md_result(module, "OK", tables = basename(tp))
}
