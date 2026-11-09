# test_targets_pipeline.R -- end-to-end check of the ExpressOM `targets` workflow
# (mod_targets.R) with cheap stand-in stages: no data, no Ensembl package, no
# DESeq2. It runs the REAL pipeline definition (expressom_pipeline()), the REAL
# parameter groups / invalidation logic and the REAL run_expressom_targets() /
# expressom_stage_status() / expressom_targets_invalidate(); only the stage
# functions themselves are replaced (through expressom_pipeline(overrides = )).
#
# How to run (about 10 seconds), from the package source folder:
#     devtools::load_all()          # or library(ExpressOM)
#     source("test_targets_pipeline.R")
# Needs: targets, callr, withr, tidyselect. Prints one line per check and stops
# with an error if any check fails.

for (pkg in c("targets", "callr", "withr", "tidyselect")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install '", pkg, "'.")
}
ns <- function(name) getFromNamespace(name, "ExpressOM")

# the folder is created first so that normalizePath() really normalises it
# (forward slashes on Windows too: the paths are embedded in string literals)
root <- tempfile("xp_targets_test_")
dir.create(root, recursive = TRUE)
root      <- normalizePath(root, winslash = "/")
pipe_dir  <- file.path(root, "pipe")
out_dir   <- file.path(root, "out")
calls_log <- file.path(root, "calls.log")
flags     <- file.path(root, "flags")
for (d in c(pipe_dir, out_dir, flags)) dir.create(d, recursive = TRUE)

# ---- stand-in stages, written into the generated _targets.R ------------------
# Each stand-in logs its name (one line per execution) and returns a `files`
# manifest that changes every time the stage re-runs, like a real stage whose
# output files were rewritten. Touch <flags>/<stage> to make a stand-in fail.
# (wrapped in local() so that nothing leaks into the global environment; the
# value of the last expression of _targets.R is what targets reads)
script_lines <- c(
  'local({',
  sprintf('LOG <- "%s"; FLAGS <- "%s"', calls_log, flags),
  'xp_run <- function(stage) {',
  '  cat(stage, "\\n", file = LOG, append = TRUE, sep = "")',
  '  if (file.exists(file.path(FLAGS, stage))) stop("simulated failure of ", stage, call. = FALSE)',
  '  n <- sum(readLines(LOG) == stage)',
  '  list(n = n, files = data.frame(path = paste0(stage, ".txt"), size = n, mtime = n))',
  '}',
  'xp_file <- function(stage, ctx, name) {',
  '  xp_run(stage)',
  '  f <- file.path(ctx$p$out_dir, name); dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)',
  '  writeLines(paste(stage, sum(readLines(LOG) == stage)), f)',
  '  f',
  '}',
  'overrides <- list(',
  '  .pipeline_setup   = function(p, fp_counts, fp_labels) { xp_run("setup");',
  '    list(p = p, main_condition = "condition", model_eda = p$model, comp_name = paste0(p$level, "vs", p$base)) },',
  '  .stage_import     = function(ctx) xp_run("import"),',
  '  .stage_eda        = function(ctx, tx_data) xp_run("eda"),',
  '  .stage_dge        = function(ctx, tx_data) xp_run("dge"),',
  '  .stage_dge_plots  = function(ctx, tx_data, dge) xp_run("dge_plots"),',
  '  .stage_functional = function(ctx, dge) xp_run("functional"),',
  '  .stage_fgsea      = function(ctx, dge) xp_run("fgsea"),',
  '  .stage_report     = function(ctx, ...) xp_file("report", ctx, "Report/report.html"),',
  '  .stage_save_rdata = function(ctx, tx_data, dge, func, isoform) xp_file("rdata", ctx, "Save_rdata/results.RData")',
  ')',
  sprintf('ExpressOM::expressom_pipeline("%s/params.rds", overrides = overrides)', pipe_dir),
  '})'
)

# ---- helpers -------------------------------------------------------------------
n_seen <- 0L
make <- function(p, ...) {
  res <- withCallingHandlers(
    try(suppressMessages(ExpressOM::run_expressom_targets(
      p, dir = pipe_dir, script_lines = script_lines, callr = FALSE,
      reporter = "silent", ...
    )), silent = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  lines <- if (file.exists(calls_log)) readLines(calls_log) else character()
  new <- lines[seq_along(lines) > n_seen]
  n_seen <<- length(lines)
  list(res = res, ran = sort(unique(new)), failed = inherits(res, "try-error"))
}
warns <- character()
fails <- 0L
check <- function(label, ok) {
  cat(sprintf("%-78s %s\n", label, if (isTRUE(ok)) "ok" else "FAIL"))
  if (!isTRUE(ok)) fails <<- fails + 1L
}
same <- function(a, b) identical(sort(as.character(a)), sort(as.character(b)))
state_of <- function(p) {
  st <- ExpressOM::expressom_stage_status(p, dir = pipe_dir)
  stats::setNames(st$state, st$stage)
}

p <- ExpressOM::expressom_params(
  out_dir = out_dir, data_dir = file.path(root, "data"),
  sample_table = file.path(root, "samples.csv"),
  level = "treated", base = "control", model = "~ condition",
  run_dge = TRUE, run_isoform = FALSE, make_report = TRUE
)

# ---- 1. first run builds everything, once -------------------------------------
r <- make(p)
check("run 1: no error", !r$failed)
check("run 1: every stage ran once",
      same(r$ran, c("setup", "import", "eda", "dge", "dge_plots", "functional",
                    "fgsea", "report", "rdata")))
check("run 1: report and .RData written",
      file.exists(file.path(out_dir, "Report/report.html")) &&
        file.exists(file.path(out_dir, "Save_rdata/results.RData")))
check("run 1: all stages 'current'", all(state_of(p) == "current"))

# ---- 2. identical parameters: nothing re-runs ---------------------------------
r <- make(p)
check("run 2 (same parameters): nothing re-ran", length(r$ran) == 0L)

# ---- 3. a cosmetic parameter only touches its figure stage + the report -------
p3 <- p; p3$volcano_colors <- c(up = "orange", down = "purple", ns = "grey80")
r <- make(p3)
check("volcano_colors: only dge_plots + report re-ran (+ setup)",
      same(r$ran, c("setup", "dge_plots", "report")))

# ---- 4. padj_cutoff re-fits DESeq2 and everything after it, not import / EDA ---
p4 <- p3; p4$padj_cutoff <- 0.05
r <- make(p4)
check("padj_cutoff: import and eda were NOT re-run",
      !any(c("import", "eda") %in% r$ran))
check("padj_cutoff: dge, plots, functional, fgsea, report, rdata re-ran",
      all(c("dge", "dge_plots", "functional", "fgsea", "report", "rdata") %in% r$ran))

# ---- 5. a failing stage: the error surfaces, finished stages are recorded -------
# change one parameter of the `import` group (count_type) and one of the `dge` group
p5 <- p4; p5$count_type <- "kallisto"; p5$padj_cutoff <- 0.03
file.create(file.path(flags, "dge"))
r <- make(p5)
check("failing dge: run_expressom_targets() raised an error", r$failed)
st <- state_of(p5)
check("failing dge: dge is 'errored'", identical(unname(st[["dge"]]), "errored"))
hashes <- readRDS(file.path(pipe_dir, "stage_hashes.rds"))
check("failing dge: import (finished before the failure) is recorded as current",
      identical(hashes$import, ns(".stage_hashes")(p5)[["import"]]))
check("failing dge: dge hash is still the OLD one (stays stale)",
      !identical(hashes$dge, ns(".stage_hashes")(p5)[["dge"]]))
unlink(file.path(flags, "dge"))

# ---- 6. recovery after fixing the cause -----------------------------------------
r <- make(p5)
check("recovery: no error", !r$failed)
check("recovery: dge ran again", "dge" %in% r$ran)
check("recovery: all stages 'current'", all(state_of(p5) == "current"))

# ---- 7. a failing report must not stop the .RData --------------------------------
p7 <- p5; p7$padj_cutoff <- 0.1
file.create(file.path(flags, "report"))
warns <- character()
r <- make(p7)
check("failing report: no error (it is a leaf)", !r$failed)
check("failing report: .RData still written", "rdata" %in% r$ran)
check("failing report: a warning names the failed stage",
      any(grepl("Interactive report", warns, fixed = TRUE)))
check("failing report: its state is 'errored'",
      identical(unname(state_of(p7)[["report"]]), "errored"))
unlink(file.path(flags, "report"))
r <- make(p7)
check("after the fix the report re-runs on its own", "report" %in% r$ran)

# ---- 8. building only up to a stage ----------------------------------------------
ExpressOM::expressom_targets_invalidate(pipe_dir, c("dge", "dge_plots", "functional", "fgsea", "report", "rdata"))
r <- make(p7, stages = "dge")
check("stages = 'dge': dge ran, report / rdata did not",
      "dge" %in% r$ran && !any(c("report", "rdata") %in% r$ran))

# ---- 9. forced re-run ------------------------------------------------------------
r <- make(p7)                                   # finish the pipeline
ExpressOM::expressom_targets_invalidate(pipe_dir, "dge")
r <- make(p7)
check("invalidate('dge'): dge + downstream re-ran, import / eda did not",
      all(c("dge", "dge_plots", "functional", "fgsea", "report", "rdata") %in% r$ran) &&
        !any(c("import", "eda") %in% r$ran))

cat(sprintf("\n%s\n", if (fails == 0L) "ALL TARGETS CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", fails)))
unlink(root, recursive = TRUE)
if (fails > 0L) stop(fails, " targets check(s) failed.", call. = FALSE)
