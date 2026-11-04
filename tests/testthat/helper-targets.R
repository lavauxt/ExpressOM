# Helpers to run the real `targets` pipeline definition with cheap stage
# stand-ins (no EnsDb, no DESeq2). Every stand-in appends its name to
# <out_dir>/calls.log, so a test can see exactly which stages re-ran.

mock_pipeline_script <- function(params_file) {
  c(
    'log_call <- function(od, stage) cat(stage, "\\n", file = file.path(od, "calls.log"), append = TRUE)',
    'mock <- function(stage) {',
    '  force(stage)',
    '  function(ctx, ...) {',
    '    od <- ctx$p$out_dir',
    '    log_call(od, stage)',
    '    if (file.exists(file.path(od, paste0("FAIL_", stage)))) stop("boom in ", stage)',
    '    n <- sum(trimws(readLines(file.path(od, "calls.log"))) == stage)',
    '    d <- file.path(od, "Plots"); dir.create(d, showWarnings = FALSE)',
    '    writeLines("x", file.path(d, paste0(stage, "_", n, ".txt")))',
    '    list(files = ExpressOM:::.out_manifest(od, "Plots"), n = n)',
    '  }',
    '}',
    'ov <- list(',
    '  .pipeline_setup = function(p, fp_counts, fp_labels) {',
    '    dir.create(p$out_dir, showWarnings = FALSE, recursive = TRUE)',
    '    log_call(p$out_dir, "setup")',
    '    list(p = p, main_condition = "condition", model_eda = p$model, comp_name = paste0(p$level, "vs", p$base))',
    '  },',
    '  .stage_import = function(ctx) { log_call(ctx$p$out_dir, "import"); list(x = unname(tools::md5sum(ctx$p$matrix_file))) },',
    '  .stage_eda = mock("eda"), .stage_dge = mock("dge"), .stage_dge_plots = mock("dge_plots"),',
    '  .stage_functional = mock("functional"), .stage_fgsea = mock("fgsea"), .stage_isoform = mock("isoform"),',
    '  .stage_report = function(ctx, ...) {',
    '    log_call(ctx$p$out_dir, "report")',
    '    if (file.exists(file.path(ctx$p$out_dir, "FAIL_report"))) stop("report boom")',
    '    f <- file.path(ctx$p$out_dir, "report.html"); writeLines(format(Sys.time(), "%H:%M:%OS6"), f); f',
    '  },',
    '  .stage_save_rdata = function(ctx, ...) {',
    '    log_call(ctx$p$out_dir, "rdata")',
    '    f <- file.path(ctx$p$out_dir, "res.RData"); writeLines("rd", f); f',
    '  }',
    ')',
    sprintf('ExpressOM::expressom_pipeline("%s", overrides = ov)', params_file)
  )
}

mock_params <- function(root, ...) {
  inputs <- make_inputs(file.path(root, "inputs"))
  args <- list(
    out_dir = file.path(root, "out"),
    data_dir = inputs$dir,
    sample_table = inputs$sample_table,
    matrix_file = inputs$matrix_file,
    count_type = "matrix",
    level = "treated", base = "control", model = "~ condition"
  )
  # later values replace the defaults above (a duplicated name would not)
  over <- list(...)
  for (nm in names(over)) {
    args[nm] <- list(over[[nm]])
  }
  do.call(expressom_params, args)
}

read_calls <- function(out_dir) {
  f <- file.path(out_dir, "calls.log")
  if (file.exists(f)) trimws(readLines(f, warn = FALSE)) else character()
}

# Run the pipeline in-process; returns the stages that ran during this call
# (in order of first appearance) and the status table.
run_mock <- function(params, ...) {
  dir <- expressom_targets_dir(params$out_dir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  pf <- file.path(normalizePath(dir, winslash = "/"), "params.rds")
  before <- length(read_calls(params$out_dir))
  status <- run_expressom_targets(
    params, script_lines = mock_pipeline_script(pf),
    callr = FALSE, reporter = "silent", ...
  )
  calls <- read_calls(params$out_dir)
  ran <- if (length(calls) > before) calls[(before + 1L):length(calls)] else character()
  list(ran = unique(ran), status = status)
}
