# =============================================================================
# run/_TEMPLATE_run_all.R —— 在新论文项目里调用本包的模板
# -----------------------------------------------------------------------------
# 用法：
#   1) 复制本文件到 <论文项目>/run/run_all.R
#   2) 改下面两处：PKG_DIR（本包位置）、CFG_FILE（配置文件名）
#   3) cd 到 <论文项目> 然后：
#        Rscript --vanilla run/run_all.R
#
# 更简单的替代做法（不用本模板）：
#   Rscript --vanilla <PKG_DIR>/R/90_run_all.R \
#       --config=config/config.R --project=<论文项目绝对路径>
# =============================================================================

PKG_DIR  <- "path/to/micro-div-pipeline"
CFG_FILE <- "config/config.R"        # 相对当前工作目录（= 论文项目根）

# 可选：只跑部分模块 / 只重出图 / 强制重算
ONLY       <- NULL                   # 例如 c("01","05","10")
FIG_ONLY   <- FALSE
FORCE      <- FALSE

# ---- 以下不用改 ----
args <- c(file.path(PKG_DIR, "R", "90_run_all.R"),
          paste0("--config=", CFG_FILE),
          paste0("--project=", normalizePath(getwd())))
if (!is.null(ONLY))  args <- c(args, paste0("--only=", paste(ONLY, collapse = ",")))
if (FIG_ONLY)        args <- c(args, "--figures-only")
if (FORCE)           args <- c(args, "--force")

status <- system2(file.path(dirname(R.home()), "bin", "Rscript.exe"),
                  c("--vanilla", shQuote(args)))
if (status != 0) stop("micro-div-pipeline 运行失败，退出码 ", status, call. = FALSE)
