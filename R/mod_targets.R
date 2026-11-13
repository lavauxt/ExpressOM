# mod_targets.R -- the ExpressOM pipeline as a `targets` workflow.
#
# The same stage functions that expressom() runs one after the other
# (mod_stages.R) become targets here. What `targets` adds:
#
#  * caching: only the stages whose inputs changed are re-run;
#  * fine-grained invalidation: each stage depends on a small *subset* of the
#    parameters (see .stage_param_groups()), so changing the volcano colours or# mod_targets.R -- the ExpressOM pipeline as a `targets` workflow.
#
# The same stage functions that expressom() runs one after the other
# (mod_stages.R) become targets here. What `targets` adds:
#
#  * caching: only the stages whose inputs changed are re-run;
#  * fine-grained invalidation: each stage depends on a small *subset* of the
#    parameters (see .stage_param_groups()), so changing the volcano colours or
#    `top_genes` never re-imports the counts or re-fits DESeq2 -- targets
#    compares the subset's value, not the whole parameter list. Note that
#    `padj_cutoff` / `deg_lfc_cutoff` ARE part of the `dge` group (they set the
#    DESeq2 `alpha` and the DEG lists), so changing them re-runs the DESeq2
#    stage and everything downstream of it -- but not the import or the EDA;
#    `deg_lfc_cutoff` is also part of the `decouple` group, so the per-TF
#    regulon plots (volcano / leading edge) use the same |log2FC| floor as
#    the rest of the pipeline;
#  * code changes are only tracked for the stage *wrapper* functions of
#    mod_stages.R (they are copied into the targets environment, see
#    .pipeline_env()). Targets does not follow the functions they call
#    (run_deseq2_analysis(), run_functional_analysis(), ...), so after
#    upgrading ExpressOM use expressom_targets_invalidate() on the stages
#    whose code changed;
#  * stage-level re-runs and status, used by the Shiny app (mod_shiny.R).
#
# Layout of a pipeline directory (default `<out_dir>/.expressom_targets/`):
#   params.rds   the parameter list (written by write_expressom_targets())
#   _targets.R   generated, one line: ExpressOM::expressom_pipeline(<params.rds>)
#   store/       the targets data store
#   stage_hashes.rds  parameter hashes of the last successful run (stale flags)

# ---- parameters ----------------------------------------------------------------------------------

#' Default parameters of the pipeline
#'
#' The default value of every argument of [expressom()], as a named list.
#'
#' @return A named list.
#' @export
#' @examples
#' p <- expressom_default_params()
#' p$padj_cutoff
expressom_default_params <- function() {
  f <- formals(expressom)
  out <- lapply(names(f), function(nm) {
    e <- f[[nm]]
    if (nm == "bpparam") {
      return(NULL) # resolved at run time (a live BiocParallel object does not serialise)
    }
    if (is.symbol(e) && identical(as.character(e), "")) {
      return(NULL)
    }
    eval(e, envir = baseenv())
  })
  names(out) <- names(f)
  out$execution_order <- out$execution_order[[1L]]
  out$isoform_test_engine <- out$isoform_test_engine[[1L]]
  out
}

# parameters that are file / directory paths (made absolute so the pipeline
# does not depend on the working directory it is launched from)
.path_param_names <- function() {
  c(
    "data_dir", "out_dir", "sample_table", "matrix_file", "custom_tx2gene",
    "custom_gene_map", "isoform_fasta", "isoform_gff", "resume_isoform_from",
    "custom_transcript_id_map", "sample_labels"
  )
}

# normalizePath() leaves a not-yet-existing relative path relative, so anchor
# it at the working directory first
.abs_path <- function(v) {
  v <- path.expand(v)
  rel <- !grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", v)
  v[rel] <- file.path(getwd(), v[rel])
  normalizePath(v, winslash = "/", mustWork = FALSE)
}

.absolutize_paths <- function(p) {
  for (nm in intersect(.path_param_names(), names(p))) {
    v <- p[[nm]]
    if (is.character(v) && length(v) == 1L && nzchar(v) && !is.na(v)) {
      # sample_labels / custom_transcript_id_map may also be inline values
      if (nm %in% c("sample_labels", "custom_transcript_id_map") &&
          !file.exists(v)) {
        next
      }
      p[[nm]] <- .abs_path(v)
    }
  }
  if (is.character(p$gmt_file)) {
    # MSigDB collection codes ("C2:CP:KEGG") stay as they are; paths become absolute
    is_file <- file.exists(p$gmt_file)
    p$gmt_file[is_file] <- .abs_path(p$gmt_file[is_file])
  }
  p
}

#' Build a parameter list for the pipeline
#'
#' Starts from [expressom_default_params()] and overrides the arguments you
#' give. Unknown names are rejected, so a typo cannot silently be ignored.
#' File and directory arguments are made absolute.
#'
#' @param ... Arguments of [expressom()] to override (e.g. `level = "treated"`).
#' @param defaults Parameter list to start from (default
#'   [expressom_default_params()]). Use [expressom_params_from_run()] to start
#'   from the parameters of a previous run.
#' @return A named list accepted by [run_expressom_targets()] and
#'   [run_expressom_app()].
#' @export
#' @examples
#' p <- expressom_params(out_dir = "results", level = "treated", base = "control")
#' p$padj_cutoff
expressom_params <- function(..., defaults = expressom_default_params()) {
  over <- list(...)
  if (length(over) && (is.null(names(over)) || any(!nzchar(names(over))))) {
    stop("All arguments to `expressom_params()` must be named.", call. = FALSE)
  }
  unknown <- setdiff(names(over), names(formals(expressom)))
  if (length(unknown)) {
    stop(
      "Unknown ExpressOM parameter(s): ", paste(unknown, collapse = ", "),
      ". See ?expressom for the available arguments.",
      call. = FALSE
    )
  }
  p <- defaults
  for (nm in names(over)) {
    p[nm] <- list(over[[nm]])
  }
  .absolutize_paths(p)
}

#' Recover the parameters of a previous run
#'
#' Reads the parameter snapshot (`log/expressom_run.json`) that [expressom()]
#' and the pipeline write next to the results, so a run can be reproduced or
#' tweaked ("same as before, but another cut-off").
#'
#' @param out_dir Output directory of a previous run.
#' @return A parameter list (see [expressom_params()]).
#' @export
expressom_params_from_run <- function(out_dir) {
  cfg <- .rep_read_json(.expressom_log_file(out_dir, "expressom_run.json"))
  if (is.null(cfg)) {
    stop(
      "No run snapshot (log/expressom_run.json) found in '", out_dir, "'.",
      call. = FALSE
    )
  }
  cfg$written_at <- NULL
  cfg$expressom_version <- NULL
  defaults <- expressom_default_params()
  keep <- intersect(names(cfg), names(defaults))
  for (nm in keep) {
    v <- cfg[[nm]]
    d <- defaults[[nm]]
    # named vectors come back from JSON as named lists
    if (is.list(v) && is.atomic(d) && !is.null(d)) {
      v <- unlist(v)
    } else if (is.list(v) && is.null(d) && all(vapply(v, is.atomic, logical(1)))) {
      v <- if (is.null(names(v))) unlist(v) else v
    }
    defaults[nm] <- list(v)
  }
  defaults$out_dir <- .abs_path(out_dir)
  defaults
}

# readRDS() that returns `default` for a missing / unreadable file, silently
.read_rds_safe <- function(path, default = NULL) {
  if (!file.exists(path)) {
    return(default)
  }
  out <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(out)) default else out
}

# ---- which stage uses which parameter ----------------------------------------------------------------

.stage_param_groups <- function() {
  common <- c("out_dir", "ensembl_package_name")
  list(
    import = c(
      common, "data_dir", "sample_table", "count_type", "matrix_file",
      "subset_sample", "remove_sample", "custom_tx2gene", "custom_gene_map"
    ),
    eda = c(
      common, "run_dge", "level", "base", "replicate_col", "group_col",
      "batch_col", "pca_ntop", "pca_colors", "sample_labels"
    ),
    dge = c(
      common, "run_dge", "level", "base", "model", "replicate_col",
      "shrink_method", "padj_cutoff", "test", "reduced", "deg_lfc_cutoff",
      "deg_padj_cutoffs"
    ),
    dge_plots = c(
      common, "level", "base", "top_genes", "padj_cutoff", "highlight_genes",
      "batch_col", "pca_ntop", "pca_colors", "sample_labels", "volcano_colors",
      "deg_lfc_cutoff", "zscore_genes", "gene_sets_zscore"
    ),
    functional = c(
      common, "level", "base", "top_genes", "padj_cutoff", "go_pvalue_cutoff",
      "go_qvalue_cutoff", "gsea_metric", "test", "ora_padj_cutoff",
      "ora_min_genes", "ora_lfc_cutoff", "run_spia"
    ),
    fgsea = c(common, "gmt_file", "padj_cutoff", "gsea_metric", "test"),
    decouple = c(
      common, "level", "base", "padj_cutoff",
      "run_decouple", "decouple_resources", "decouple_method",
      "decouple_min_size",
      # NEW: |log2FC| floor used by the per-TF regulon volcano and the
      # "Activé / Réprimé / Autre" categorisation in tf_regulon_summary.csv.
      # Adding it here means changing deg_lfc_cutoff invalidates the
      # decouple stage and regenerates those plots with the new floor.
      "deg_lfc_cutoff"
    ),
    isoform = c(
      common, "data_dir", "sample_table", "count_type", "matrix_file",
      "subset_sample", "remove_sample", "custom_tx2gene", "custom_gene_map",
      "custom_transcript_id_map", "isoform_fasta", "isoform_gff", "model",
      "level", "base", "test", "reduced", "padj_cutoff", "run_isoform",
      "run_dexseq", "run_isoform_enrichment", "ora_padj_cutoff",
      "ora_min_genes", "ora_lfc_cutoff", "gsea_metric", "gmt_file",
      "top_genes", "go_pvalue_cutoff", "go_qvalue_cutoff", "run_predictors",
      "predictor_cpu", "use_wsl", "wsl_distro", "run_spia",
      "isoform_test_engine", "resume_isoform_from", "batch_col", "pca_ntop",
      "pca_colors", "sample_labels", "skip_fasta_filter", "plot_topology",
      "isoform_report_genes", "isoform_plot_top_n", "bpparam", "make_report"
    ),
    report = c(
      "out_dir", "run_dge", "level", "base", "report_settings", "make_report"
    ),
    rdata = c(
      common, "level", "base", "model", "test", "batch_col", "run_dge",
      "run_isoform"
    )
  )
}

# Values .expressom_prepare() derives from several parameters. A stage keeps
# only the ones it reads: otherwise changing `level` would change the context of
# the import stage (comp_name) and needlessly re-import the counts.
.stage_derived_values <- function() {
  list(
    import   = character(),
    eda      = c("main_condition", "model_eda"),
    dge      = character(),
    dge_plots = "main_condition",
    functional = character(),
    fgsea    = "comp_name",
    decouple = "comp_name",
    isoform  = "main_condition",
    report   = character(),
    rdata    = c("main_condition", "comp_name")
  )
}

#' Restrict a stage context to the parameters its stage uses
#' @keywords internal
.stage_ctx <- function(ctx, stage) {
  keep <- .stage_param_groups()[[stage]]
  ctx$p <- ctx$p[intersect(keep, names(ctx$p))]
  derived <- c("main_condition", "model_eda", "comp_name")
  for (nm in setdiff(derived, .stage_derived_values()[[stage]])) {
    ctx[nm] <- list(NULL)
  }
  ctx
}

#' Pipeline stages: id, label and the targets they correspond to
#'
#' @return A data.frame with columns `stage`, `label`, `target` and `upstream`
#'   (space-separated stages it reads from).
#' @export
expressom_pipeline_stages <- function() {
  data.frame(
    stage = c("setup", "import", "eda", "dge", "dge_plots", "functional",
              "fgsea", "decouple", "isoform", "report", "rdata"),
    label = c("Validate parameters", "Import counts", "EDA / PCA",
              "DESeq2 differential expression", "DGE figures",
              "ORA / GSEA / SPIA", "FGSEA",
              "TF activity (decoupleR)", "Isoforms (DTE / DTU / switch)",
              "Interactive report", "Save .RData"),
    target = c("setup", "import", "eda", "dge", "dge_plots", "functional",
               "fgsea", "decouple", "isoform", "report", "rdata"),
    upstream = c("", "setup", "import", "import", "dge", "dge", "dge", "dge",
                 "setup", "eda dge_plots functional fgsea decouple isoform",
                 "import dge functional decouple isoform"),
    stringsAsFactors = FALSE
  )
}

# which optional stages exist for these parameters (mirrors .expressom_prepare())
.pipeline_flags <- function(p) {
  full <- !isTRUE(p$eda_only) && !is.null(p$model) && !is.null(p$level) &&
    !is.null(p$base)
  list(
    run_dge      = full && isTRUE(p$run_dge),
    run_decouple = full && isTRUE(p$run_dge) && isTRUE(p$run_decouple),
    run_isoform  = full && isTRUE(p$run_isoform),
    make_report  = isTRUE(p$make_report)
  )
}

# one hash per stage, over the parameters that stage uses (sorted by name)
.stage_hashes <- function(p) {
  vapply(.stage_param_groups(), function(keep) {
    used <- sort(intersect(keep, names(p)))
    .object_md5(p[used])
  }, character(1))
}

# ---- targets definition -------------------------------------------------------------------------------

.fp_slim <- function(fp) {
  # path + size + md5: a touched-but-identical file must not invalidate anything
  if (is.data.frame(fp)) fp[, c("path", "size", "md5"), drop = FALSE] else fp
}

.fp_counts <- function(p) {
  .fp_slim(.file_fingerprint(c(
    p$sample_table, p$matrix_file, p$custom_tx2gene, p$custom_gene_map,
    .quant_input_files(p, gene_level = TRUE)
  )))
}

.fp_labels <- function(p) {
  v <- p$sample_labels
  .fp_slim(if (is.character(v) && length(v) == 1L) .file_fingerprint(v) else character())
}

.fp_isoform <- function(p) {
  paths <- c(p$isoform_fasta, p$isoform_gff, .quant_input_files(p))
  if (is.character(p$custom_transcript_id_map) &&
      length(p$custom_transcript_id_map) == 1L) {
    paths <- c(paths, p$custom_transcript_id_map)
  }
  .fp_slim(.file_fingerprint(paths))
}

.pipeline_setup <- function(p, fp_counts, fp_labels) {
  ctx <- .expressom_prepare(p)
  .write_run_config(ctx$p$out_dir, ctx$p)
  ctx
}

# The functions the target commands call. They are copied into the environment
# targets evaluates in, so the generated _targets.R only needs
# ExpressOM::expressom_pipeline() -- no `:::` and no library(ExpressOM) -- and
# targets hashes their code: reinstalling a package with an edited stage
# *wrapper* re-runs that stage (the code those wrappers call is not tracked, see
# the header of this file). `overrides` (a named list of functions) replaces any
# of them; the test-suite uses it to run the DAG with cheap stand-ins.
.pipeline_functions <- function() {
  c(
    ".stage_ctx", ".fp_counts", ".fp_labels", ".fp_isoform", ".pipeline_setup",
    ".predictor_environment_fingerprint",
    ".stage_import", ".stage_eda", ".stage_dge", ".stage_dge_plots",
    ".stage_functional", ".stage_fgsea", ".stage_decouple",
    ".stage_isoform", ".stage_report", ".stage_save_rdata"
  )
}

.pipeline_env <- function(overrides = list()) {
  ns <- environment(.pipeline_env)
  env <- new.env(parent = globalenv())
  for (nm in .pipeline_functions()) {
    assign(nm, get(nm, envir = ns), envir = env)
  }
  for (nm in names(overrides)) {
    assign(nm, overrides[[nm]], envir = env)
  }
  env
}

#' The ExpressOM `targets` pipeline
#'
#' Returns the list of targets for the parameters stored in `params_file`.
#' It is what the generated `_targets.R` calls; you normally do not call it
#' yourself -- use [run_expressom_targets()] or [run_expressom_app()].
#'
#' The stages form a DAG: `setup` -> `import` -> `eda` / `dge` -> `dge_plots`,
#' `functional`, `fgsea`, `decouple`; `isoform` is independent of `dge`;
#' `report` and `rdata` close the graph. Each stage target depends on a small
#' context target (`ctx_<stage>`) holding only the parameters it uses, so a
#' parameter that concerns one stage leaves the cached results of the others
#' untouched. When `run_predictors = TRUE`, the `isoform` target additionally
#' depends on `fp_predictors` (the external predictor environment fingerprint),
#' so it re-runs when a predictor tool or database is installed or removed even
#' though no parameter or input file changed.
#'
#' @param params_file Path to the `params.rds` written by
#'   [write_expressom_targets()].
#' @param overrides Named list of functions replacing internal pipeline
#'   functions (for tests).
#' @return A list of `targets` target objects.
#' @export
expressom_pipeline <- function(params_file, overrides = list()) {
  .rep_require("targets", "to run the pipeline with targets")
  p0 <- readRDS(params_file)
  flags <- .pipeline_flags(p0)

  # Every WSL/predictor log must land under <out_dir>/log/WSL/ instead of the
  # targets working directory. Set it here (in the child process that
  # tar_make() spins up) so all downstream .wsl_exec_script() calls -- target
  # code AND .predictor_environment_fingerprint() probes -- write to log/WSL.
  if (!is.null(p0$out_dir) && nzchar(p0$out_dir)) {
    wsl_log_dir <- file.path(.expressom_log_dir(p0$out_dir, create = TRUE), "WSL")
    if (!dir.exists(wsl_log_dir)) {
      dir.create(wsl_log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    options(expressom.wsl_log_dir = wsl_log_dir)
  }

  targets::tar_option_set(
    envir = .pipeline_env(overrides),
    packages = character(0),
    format = "rds",
    memory = "transient",
    garbage_collection = TRUE,
    error = "stop"
  )

  tgt <- function(name, command, ...) {
    targets::tar_target_raw(name, command, ...)
  }
  ctx_target <- function(stage) {
    tgt(paste0("ctx_", stage), bquote(.stage_ctx(setup, .(stage))))
  }
  always <- targets::tar_cue(mode = "always")

  out <- list(
    tgt("params_file", params_file, format = "file"),
    tgt("params", quote(readRDS(params_file))),
    # fingerprints of the *content* of the input files, re-computed each run
    tgt("fp_counts", quote(.fp_counts(params)), cue = always),
    tgt("fp_labels", quote(.fp_labels(params)), cue = always),
    tgt("setup", quote(.pipeline_setup(params, fp_counts, fp_labels))),
    ctx_target("import"),
    ctx_target("eda"),
    tgt("import", quote({
      fp_counts
      .stage_import(ctx_import)
    })),
    tgt("eda", quote(.stage_eda(ctx_eda, import))),
    tgt("m_eda", quote(eda$files))
  )
  manifests <- "m_eda"

  if (flags$run_dge) {
    out <- c(out, list(
      ctx_target("dge"),
      ctx_target("dge_plots"),
      ctx_target("functional"),
      ctx_target("fgsea"),
      tgt("dge", quote(.stage_dge(ctx_dge, import))),
      tgt("dge_plots", quote(.stage_dge_plots(ctx_dge_plots, import, dge))),
      tgt("functional", quote(.stage_functional(ctx_functional, dge))),
      tgt("fgsea", quote(.stage_fgsea(ctx_fgsea, dge))),
      tgt("m_dge", quote(dge$files)),
      tgt("m_dge_plots", quote(dge_plots$files)),
      tgt("m_functional", quote(functional$files)),
      tgt("m_fgsea", quote(fgsea$files))
    ))
    manifests <- c(manifests, "m_dge", "m_dge_plots", "m_functional", "m_fgsea")
  }

  if (isTRUE(flags$run_decouple)) {
    out <- c(out, list(
      ctx_target("decouple"),
      tgt("decouple", quote(.stage_decouple(ctx_decouple, dge))),
      tgt("m_decouple", quote(decouple$files))
    ))
    manifests <- c(manifests, "m_decouple")
  }

  if (flags$run_isoform) {
    isoform_deps <- list(
      tgt("fp_isoform", quote(.fp_isoform(params)), cue = always),
      ctx_target("isoform")
    )
    # Re-probe the predictor environment on every run if predictors are
    # enabled (~1-2s WSL round-trip). Without this, "isoform" stays cached
    # even after a new predictor tool / database was installed, and the
    # pipeline reuses the old, predictor-empty checkpoints. Combined with
    # `predictor_env` now being part of `.stage_isoform()`'s checkpoint
    # signature, this is what makes "install a tool → re-run → predictor
    # results appear" work end-to-end.
    #
    # The explicit log_dir makes the WSL commands it issues land under
    # log/WSL/ (matching the option set at the top of expressom_pipeline()),
    # so the fingerprint probe shows up in the report's WSL log section.
    if (isTRUE(p0$run_predictors)) {
      isoform_deps <- c(isoform_deps, list(
        tgt(
          "fp_predictors",
          quote(.predictor_environment_fingerprint(
            wsl_distro = params$wsl_distro,
            use_wsl    = params$use_wsl,
            log_dir    = file.path(.expressom_log_dir(params$out_dir, create = TRUE), "WSL")
          )),
          cue = always
        )
      ))
    }
    isoform_body <- if (isTRUE(p0$run_predictors)) {
      quote({
        fp_counts
        fp_isoform
        fp_predictors
        .stage_isoform(ctx_isoform)
      })
    } else {
      quote({
        fp_counts
        fp_isoform
        .stage_isoform(ctx_isoform)
      })
    }
    out <- c(out, c(isoform_deps, list(
      tgt("isoform", isoform_body),
      tgt("m_isoform", quote(isoform$files))
    )))
    manifests <- c(manifests, "m_isoform")
  }

  if (flags$make_report) {
    report_cmd <- as.call(c(
      list(quote(.stage_report), quote(ctx_report)),
      lapply(manifests, as.name)
    ))
    # a failing report must not keep the .RData from being written (nor the
    # reverse): both are leaves of the graph
    out <- c(out, list(
      ctx_target("report"),
      tgt("report", report_cmd, format = "file", error = "continue")
    ))
  }

  if (flags$run_dge || flags$run_isoform) {
    arg <- function(target, present) if (present) as.name(target) else quote(NULL)
    rdata_cmd <- call(
      ".stage_save_rdata", quote(ctx_rdata), quote(import),
      arg("dge",      flags$run_dge),
      arg("functional", flags$run_dge),
      arg("isoform",  flags$run_isoform),
      arg("decouple", flags$run_decouple)
    )
    out <- c(out, list(
      ctx_target("rdata"),
      tgt("rdata", rdata_cmd, format = "file", error = "continue")
    ))
  }
  out
}

# ---- writing / running -----------------------------------------------------------------------------------

#' Default directory of the pipeline for an output directory
#' @param out_dir Output directory of the run.
#' @return `<out_dir>/.expressom_targets`.
#' @export
expressom_targets_dir <- function(out_dir) {
  file.path(out_dir, ".expressom_targets")
}

#' Write the `targets` pipeline for a set of parameters
#'
#' Stores the parameters (`params.rds`) and generates the one-line
#' `_targets.R`. Re-writing identical parameters leaves the store untouched, so
#' nothing is invalidated; changed parameters invalidate exactly the stages
#' that use them.
#'
#' @param params Parameter list from [expressom_params()].
#' @param dir Pipeline directory (default [expressom_targets_dir()] of
#'   `params$out_dir`).
#' @param script_lines Lines of `_targets.R`; by default the one-liner calling
#'   [expressom_pipeline()] (a custom script is used by the tests).
#' @return A list with `dir`, `script`, `store` and `params_file` (invisibly).
#' @export
write_expressom_targets <- function(params, dir = NULL, script_lines = NULL) {
  if (!is.list(params) || is.null(params$out_dir)) {
    stop("`params` must be a list from expressom_params() with an `out_dir`.", call. = FALSE)
  }
  params <- .absolutize_paths(params)
  dir <- dir %||% expressom_targets_dir(params$out_dir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir, winslash = "/", mustWork = TRUE)

  params_file <- file.path(dir, "params.rds")
  # only touch the file when the content changed: its hash is a target, but a
  # stable mtime also keeps `tar_make()` from having to re-read it
  old <- .read_rds_safe(params_file)
  if (is.null(old) || !identical(.object_md5(old), .object_md5(params))) {
    saveRDS(params, params_file)
  }

  script <- file.path(dir, "_targets.R")
  if (is.null(script_lines)) {
    script_lines <- c(
      "# Generated by ExpressOM::write_expressom_targets() -- do not edit.",
      "# The parameters live in params.rds next to this file; change them with",
      "# expressom_params() + run_expressom_targets(), or in the Shiny app.",
      sprintf('ExpressOM::expressom_pipeline("%s")', params_file)
    )
  }
  writeLines(script_lines, script)

  invisible(list(
    dir = dir,
    script = script,
    store = file.path(dir, "store"),
    params_file = params_file
  ))
}

# targets (and their upstream) needed to build `stages`
.stage_targets <- function(stages) {
  tbl <- expressom_pipeline_stages()
  bad <- setdiff(stages, tbl$stage)
  if (length(bad)) {
    stop(
      "Unknown stage(s): ", paste(bad, collapse = ", "),
      ". Choose among: ", paste(tbl$stage, collapse = ", "), ".",
      call. = FALSE
    )
  }
  tbl$target[match(stages, tbl$stage)]
}

#' Run the pipeline with `targets`
#'
#' Writes the pipeline for `params` and runs `targets::tar_make()`. Stages whose
#' inputs did not change are skipped; the interactive report is rebuilt when
#' (and only when) one of the stages it reads from re-ran.
#'
#' @details **Re-run only one stage.**
#'
#' `stages = "dge"` builds the DESeq2 stage and its ancestors, nothing more.
#' To force a stage that is already up to date, invalidate it first:
#'
#' ```r
#' pipe <- expressom_targets_dir(out_dir)
#' expressom_targets_invalidate(pipe, "decouple")
#' run_expressom_targets(stages = c("decouple", "report"), dir = pipe)
#' ```
#'
#' The convenience wrapper `rerun_decouple_only(out_dir)` does exactly that
#' for the decoupleR stage (the one most likely to need a manual re-run:
#' OmniPath unreachable at first attempt, static fallback used for
#' CollecTRI/DoRothEA, PROGENy skipped, …). Passing `stages = "report"`
#' alone also works -- `targets` walks the ancestor graph, so `decouple` is
#' rebuilt if it was invalidated, and everything else is skipped when
#' nothing changed.
#'
#' @param params Parameter list from [expressom_params()]. `NULL` re-uses the
#'   parameters already stored in `dir`.
#' @param dir Pipeline directory (default `<out_dir>/.expressom_targets`).
#' @param stages Optional character vector of stages to build, among
#'   [expressom_pipeline_stages()]`$stage`. Each is built together with the
#'   stages it needs; nothing else runs (`stages = "dge"` stops after the
#'   DESeq2 stage, without functional analysis or report). `NULL` builds
#'   everything.
#' @param callr Run in a clean background R process (default; recommended).
#' @param reporter `targets` console reporter.
#' @param script_lines Custom `_targets.R` lines (see
#'   [write_expressom_targets()]); mostly for tests.
#' @param ... Further arguments passed to [targets::tar_make()].
#' @return The stage status table ([expressom_stage_status()]), invisibly. A
#'   stage that stops the pipeline (import, DESeq2 ...) raises an error; the
#'   report and the `.RData` are leaves that do not block each other, so if one
#'   of them fails the run completes with a warning naming it.
#' @export
#' @examples
#' \dontrun{
#' p <- expressom_params(
#'   data_dir = "data", out_dir = "results", sample_table = "samples.csv",
#'   level = "treated", base = "control", model = "~ condition"
#' )
#' run_expressom_targets(p)                      # first run: everything
#' p$padj_cutoff <- 0.01
#' run_expressom_targets(p)                      # only the stages using padj_cutoff
#' run_expressom_targets(p, stages = "dge")      # up to the DESeq2 stage only
#' }
run_expressom_targets <- function(
  params = NULL,
  dir = NULL,
  stages = NULL,
  callr = TRUE,
  reporter = "summary",
  script_lines = NULL,
  ...
) {
  .rep_require("targets", "to run the pipeline with targets")
  if (is.null(params)) {
    if (is.null(dir)) {
      stop("Give `params`, or the `dir` of an existing pipeline.", call. = FALSE)
    }
    params <- readRDS(file.path(dir, "params.rds"))
  }
  wanted <- if (is.null(stages)) NULL else .stage_targets(stages)

  paths <- write_expressom_targets(params, dir, script_lines)
  params <- readRDS(paths$params_file)

  # The child process that runs tar_make() inherits this option only if it
  # was set inside expressom_pipeline() (which is what actually matters for
  # the targets DAG). Setting it here too covers edge cases where a target
  # shells out before expressom_pipeline() has run (e.g. custom _targets.R).
  old_wsl_log <- getOption("expressom.wsl_log_dir")
  options(expressom.wsl_log_dir = file.path(.expressom_log_dir(params$out_dir, create = TRUE), "WSL"))
  on.exit(options(expressom.wsl_log_dir = old_wsl_log), add = TRUE)

  args <- list(
    script = paths$script,
    store = paths$store,
    reporter = reporter,
    callr_function = if (isTRUE(callr)) callr::r else NULL
  )
  if (!is.null(wanted)) {
    args$names <- bquote(tidyselect::any_of(.(wanted)))
  }
  args <- c(args, list(...))

  # A stage that stops the pipeline (error = "stop") makes tar_make() throw --
  # but the stages that completed BEFORE it are up to date in the store. Record
  # their parameter hashes either way, otherwise the app keeps flagging them
  # "stale" after a failed run and offers to re-run work that is already done.
  make_error <- NULL
  tryCatch(
    withr::with_dir(
      paths$dir,
      do.call(targets::tar_make, args)
    ),
    error = function(e) make_error <<- e
  )

  .record_stage_hashes(params, paths$dir)
  if (!is.null(make_error)) {
    stop(make_error)
  }
  status <- expressom_stage_status(params, dir = paths$dir)

  # stages that may fail without stopping the pipeline (the report and the
  # .RData are leaves): say so instead of finishing silently
  failed <- status[status$state == "errored", , drop = FALSE]
  if (nrow(failed)) {
    meta <- tryCatch(
      targets::tar_meta(store = paths$store, fields = "error"),
      error = function(e) NULL
    )
    tg <- expressom_pipeline_stages()$target[match(failed$stage, expressom_pipeline_stages()$stage)]
    msgs <- vapply(tg, function(t) {
      m <- if (is.null(meta)) NA_character_ else as.character(meta$error[match(t, meta$name)])
      if (length(m) && !is.na(m) && nzchar(m)) paste0(": ", m) else ""
    }, character(1))
    warning(
      "These stages failed: ", paste0(failed$label, msgs, collapse = "; "),
      ". The other stages completed; fix the cause and run again.",
      call. = FALSE
    )
  }
  invisible(status)
}

# After a make: remember, per stage, the parameter hash it is now consistent
# with. Only stages the make touched (completed or skipped = up to date) get
# their hash refreshed; the others stay as they were, so they keep showing as
# stale in the app.
.record_stage_hashes <- function(params, dir) {
  hashes_file <- file.path(dir, "stage_hashes.rds")
  old <- .read_rds_safe(hashes_file, default = list())
  store <- file.path(dir, "store")
  prog <- tryCatch(targets::tar_progress(store = store), error = function(e) NULL)
  if (is.null(prog) || !nrow(prog)) {
    return(invisible(old))
  }
  tbl <- expressom_pipeline_stages()
  cur <- .stage_hashes(params)
  for (st in names(cur)) {
    tg <- tbl$target[tbl$stage == st]
    ok <- prog$progress[prog$name == tg]
    if (length(ok) && ok[[1L]] %in% c("completed", "skipped")) {
      old[[st]] <- cur[[st]]
    }
  }
  saveRDS(old, hashes_file)
  invisible(old)
}

#' Re-run only the decoupleR stage (and refresh the report)
#'
#' Shortcut for `expressom_targets_invalidate(pipe, "decouple")` followed by
#' a run of `decouple` (and, optionally, `report`). Nothing else recomputes:
#' import, EDA, DESeq2, functional analysis, FGSEA and isoform analysis are
#' all skipped because their inputs did not change.
#'
#' Use this when the decoupleR stage silently degraded -- typically because
#' OmniPath was unreachable on the first run (CollecTRI/DoRothEA fell back
#' to older static tables, PROGENy was skipped) -- and you want it to run
#' again now that the network or proxy is fixed. The per-resource fetch
#' status is written to
#' `<out_dir>/DecoupleR/decoupleR_network_status_<comp>.csv`.
#'
#' @param out_dir Output directory of a run (its `.expressom_targets/` pipeline
#'   directory must already exist -- run the pipeline once before calling this).
#' @param rebuild_report Also refresh the interactive HTML report once the
#'   stage has re-run (default `TRUE`). Set `FALSE` to only rewrite the
#'   `DecoupleR/` tables.
#' @param ... Passed to [run_expressom_targets()].
#' @return The stage status table ([expressom_stage_status()]), invisibly.
#' @export
#' @examples
#' \dontrun{
#' # First run succeeded but PROGENy was skipped (offline)
#' rerun_decouple_only("results")
#'
#' # Fix proxy, then retry only the decoupleR stage without touching the report
#' Sys.setenv(https_proxy = "http://proxy:8080")
#' rerun_decouple_only("results", rebuild_report = FALSE)
#' }
rerun_decouple_only <- function(out_dir, rebuild_report = TRUE, ...) {
  if (!is.character(out_dir) || length(out_dir) != 1L || !dir.exists(out_dir)) {
    stop("`out_dir` must be an existing ExpressOM output directory.",
         call. = FALSE)
  }
  pipe <- expressom_targets_dir(out_dir)
  if (!dir.exists(pipe) || !file.exists(file.path(pipe, "params.rds"))) {
    stop(
      "No pipeline found in '", pipe, "'.\n",
      "Run the pipeline once with run_expressom_targets() or the Shiny app ",
      "before calling rerun_decouple_only().",
      call. = FALSE
    )
  }
  expressom_targets_invalidate(pipe, "decouple")
  stages <- c("decouple", if (isTRUE(rebuild_report)) "report")
  run_expressom_targets(stages = stages, dir = pipe, ...)
}

#' Status of every pipeline stage
#'
#' For each stage: whether it ran, how long it took, and whether the current
#' parameters differ from the ones it last ran with. Staleness is propagated
#' downstream: a stage is stale when its own parameters changed or any stage it
#' reads from is stale.
#'
#' @param params Current parameter list (to compute staleness); `NULL` uses the
#'   parameters stored in `dir`.
#' @param dir Pipeline directory.
#' @param running Is a run in progress right now? `NULL` (default) trusts the
#'   `targets` progress file; `FALSE` treats stages left "running" by a killed
#'   process as errored (used by the app).
#' @return A data.frame with `stage`, `label`, `progress`, `seconds`,
#'   `finished` and `state` ("never run", "current", "stale", "running",
#'   "errored"). Only the stages that exist for these parameters are listed.
#' @export
expressom_stage_status <- function(params = NULL, dir, running = NULL) {
  stages <- expressom_pipeline_stages()
  store <- file.path(dir, "store")

  if (is.null(params)) {
    pf <- file.path(dir, "params.rds")
    params <- if (file.exists(pf)) readRDS(pf) else list()
  }
  flags <- .pipeline_flags(params)

  prog <- NULL
  meta <- NULL
  if (dir.exists(store) && requireNamespace("targets", quietly = TRUE)) {
    prog <- tryCatch(targets::tar_progress(store = store), error = function(e) NULL)
    meta <- tryCatch(
      targets::tar_meta(store = store, fields = c("time", "seconds", "error")),
      error = function(e) NULL
    )
  }
  pick <- function(df, col, tg, default) {
    if (is.null(df) || !tg %in% df$name) default else df[[col]][match(tg, df$name)]
  }

  stages$progress <- vapply(
    stages$target, function(t) as.character(pick(prog, "progress", t, "none")), character(1)
  )
  stages$seconds <- vapply(
    stages$target, function(t) as.numeric(pick(meta, "seconds", t, NA_real_)), numeric(1)
  )
  stages$finished <- vapply(stages$target, function(t) {
    tm <- pick(meta, "time", t, NA)
    if (length(tm) != 1L || is.na(tm)) NA_character_ else format(tm, "%Y-%m-%d %H:%M:%S")
  }, character(1))
  stages$errored <- vapply(stages$target, function(t) {
    e <- pick(meta, "error", t, NA_character_)
    !is.na(e) && nzchar(e)
  }, logical(1))
  stages$built <- stages$target %in% (if (is.null(meta)) character() else meta$name[is.na(meta$error)])

  # stages that do not exist for these parameters
  exists <- stats::setNames(rep(TRUE, nrow(stages)), stages$stage)
  exists[c("dge", "dge_plots", "functional", "fgsea")] <- flags$run_dge
  exists["decouple"] <- flags$run_decouple
  exists["isoform"]  <- flags$run_isoform
  exists["report"]   <- flags$make_report
  exists["rdata"]    <- flags$run_dge || flags$run_isoform
  stages <- stages[exists[stages$stage], , drop = FALSE]

  last <- .read_rds_safe(file.path(dir, "stage_hashes.rds"), default = list())
  now <- .stage_hashes(params)
  # "setup" has no parameter group of its own: it re-validates on any change
  # and is never reported as stale by itself
  stale <- stats::setNames(
    vapply(stages$stage, function(s) {
      s %in% names(now) && !identical(last[[s]], now[[s]])
    }, logical(1)),
    stages$stage
  )
  for (i in seq_len(nrow(stages))) {
    ups <- intersect(strsplit(stages$upstream[[i]], " ", fixed = TRUE)[[1L]], stages$stage)
    if (length(ups) && any(stale[ups])) {
      stale[[stages$stage[[i]]]] <- TRUE
    }
  }

  state <- unname(ifelse(stale[stages$stage], "stale", "current"))
  state[!stages$built] <- "never run"
  state[stages$errored | stages$progress %in% c("errored", "canceled")] <- "errored"
  active <- stages$progress %in% c("dispatched", "running")
  state[active] <- if (isFALSE(running)) "errored" else "running"
  stages$state <- state
  rownames(stages) <- NULL
  stages[, c("stage", "label", "progress", "seconds", "finished", "state")]
}

#' Force stages to re-run
#'
#' Marks the targets of `stages` as outdated so the next
#' [run_expressom_targets()] rebuilds them (and the report) even though their
#' inputs did not change -- for instance after editing an input file in place
#' or deleting figures.
#'
#' @param dir Pipeline directory.
#' @param stages Stages to invalidate (see [expressom_pipeline_stages()]).
#' @return `NULL` (invisibly).
#' @export
expressom_targets_invalidate <- function(dir, stages) {
  .rep_require("targets", "to run the pipeline with targets")
  tg <- .stage_targets(stages)
  tg <- c(tg, paste0("m_", stages))
  store <- file.path(dir, "store")
  if (dir.exists(store)) {
    targets::tar_invalidate(
      names = tidyselect::any_of(tg),
      store = store
    )
  }
  hashes_file <- file.path(dir, "stage_hashes.rds")
  if (file.exists(hashes_file)) {
    h <- readRDS(hashes_file)
    h[stages] <- NULL
    saveRDS(h, hashes_file)
  }
  invisible(NULL)
}
# ---- parameters ----------------------------------------------------------------------------------

#' Default parameters of the pipeline
#'
#' The default value of every argument of [expressom()], as a named list.
#'
#' @return A named list.
#' @export
#' @examples
#' p <- expressom_default_params()
#' p$padj_cutoff
expressom_default_params <- function() {
  f <- formals(expressom)
  out <- lapply(names(f), function(nm) {
    e <- f[[nm]]
    if (nm == "bpparam") {
      return(NULL) # resolved at run time (a live BiocParallel object does not serialise)
    }
    if (is.symbol(e) && identical(as.character(e), "")) {
      return(NULL)
    }
    eval(e, envir = baseenv())
  })
  names(out) <- names(f)
  out$execution_order <- out$execution_order[[1L]]
  out$isoform_test_engine <- out$isoform_test_engine[[1L]]
  out
}

# parameters that are file / directory paths (made absolute so the pipeline
# does not depend on the working directory it is launched from)
.path_param_names <- function() {
  c(
    "data_dir", "out_dir", "sample_table", "matrix_file", "custom_tx2gene",
    "custom_gene_map", "isoform_fasta", "isoform_gff", "resume_isoform_from",
    "custom_transcript_id_map", "sample_labels"
  )
}

# normalizePath() leaves a not-yet-existing relative path relative, so anchor
# it at the working directory first
.abs_path <- function(v) {
  v <- path.expand(v)
  rel <- !grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", v)
  v[rel] <- file.path(getwd(), v[rel])
  normalizePath(v, winslash = "/", mustWork = FALSE)
}

.absolutize_paths <- function(p) {
  for (nm in intersect(.path_param_names(), names(p))) {
    v <- p[[nm]]
    if (is.character(v) && length(v) == 1L && nzchar(v) && !is.na(v)) {
      # sample_labels / custom_transcript_id_map may also be inline values
      if (nm %in% c("sample_labels", "custom_transcript_id_map") &&
          !file.exists(v)) {
        next
      }
      p[[nm]] <- .abs_path(v)
    }
  }
  if (is.character(p$gmt_file)) {
    # MSigDB collection codes ("C2:CP:KEGG") stay as they are; paths become absolute
    is_file <- file.exists(p$gmt_file)
    p$gmt_file[is_file] <- .abs_path(p$gmt_file[is_file])
  }
  p
}

#' Build a parameter list for the pipeline
#'
#' Starts from [expressom_default_params()] and overrides the arguments you
#' give. Unknown names are rejected, so a typo cannot silently be ignored.
#' File and directory arguments are made absolute.
#'
#' @param ... Arguments of [expressom()] to override (e.g. `level = "treated"`).
#' @param defaults Parameter list to start from (default
#'   [expressom_default_params()]). Use [expressom_params_from_run()] to start
#'   from the parameters of a previous run.
#' @return A named list accepted by [run_expressom_targets()] and
#'   [run_expressom_app()].
#' @export
#' @examples
#' p <- expressom_params(out_dir = "results", level = "treated", base = "control")
#' p$padj_cutoff
expressom_params <- function(..., defaults = expressom_default_params()) {
  over <- list(...)
  if (length(over) && (is.null(names(over)) || any(!nzchar(names(over))))) {
    stop("All arguments to `expressom_params()` must be named.", call. = FALSE)
  }
  unknown <- setdiff(names(over), names(formals(expressom)))
  if (length(unknown)) {
    stop(
      "Unknown ExpressOM parameter(s): ", paste(unknown, collapse = ", "),
      ". See ?expressom for the available arguments.",
      call. = FALSE
    )
  }
  p <- defaults
  for (nm in names(over)) {
    p[nm] <- list(over[[nm]])
  }
  .absolutize_paths(p)
}

#' Recover the parameters of a previous run
#'
#' Reads the parameter snapshot (`log/expressom_run.json`) that [expressom()]
#' and the pipeline write next to the results, so a run can be reproduced or
#' tweaked ("same as before, but another cut-off").
#'
#' @param out_dir Output directory of a previous run.
#' @return A parameter list (see [expressom_params()]).
#' @export
expressom_params_from_run <- function(out_dir) {
  cfg <- .rep_read_json(.expressom_log_file(out_dir, "expressom_run.json"))
  if (is.null(cfg)) {
    stop(
      "No run snapshot (log/expressom_run.json) found in '", out_dir, "'.",
      call. = FALSE
    )
  }
  cfg$written_at <- NULL
  cfg$expressom_version <- NULL
  defaults <- expressom_default_params()
  keep <- intersect(names(cfg), names(defaults))
  for (nm in keep) {
    v <- cfg[[nm]]
    d <- defaults[[nm]]
    # named vectors come back from JSON as named lists
    if (is.list(v) && is.atomic(d) && !is.null(d)) {
      v <- unlist(v)
    } else if (is.list(v) && is.null(d) && all(vapply(v, is.atomic, logical(1)))) {
      v <- if (is.null(names(v))) unlist(v) else v
    }
    defaults[nm] <- list(v)
  }
  defaults$out_dir <- .abs_path(out_dir)
  defaults
}

# readRDS() that returns `default` for a missing / unreadable file, silently
.read_rds_safe <- function(path, default = NULL) {
  if (!file.exists(path)) {
    return(default)
  }
  out <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(out)) default else out
}

# ---- which stage uses which parameter ----------------------------------------------------------------

.stage_param_groups <- function() {
  common <- c("out_dir", "ensembl_package_name")
  list(
    import = c(
      common, "data_dir", "sample_table", "count_type", "matrix_file",
      "subset_sample", "remove_sample", "custom_tx2gene", "custom_gene_map"
    ),
    eda = c(
      common, "run_dge", "level", "base", "replicate_col", "group_col",
      "batch_col", "pca_ntop", "pca_colors", "sample_labels"
    ),
    dge = c(
      common, "run_dge", "level", "base", "model", "replicate_col",
      "shrink_method", "padj_cutoff", "test", "reduced", "deg_lfc_cutoff",
      "deg_padj_cutoffs"
    ),
    dge_plots = c(
      common, "level", "base", "top_genes", "padj_cutoff", "highlight_genes",
      "batch_col", "pca_ntop", "pca_colors", "sample_labels", "volcano_colors",
      "deg_lfc_cutoff", "zscore_genes", "gene_sets_zscore"
    ),
    functional = c(
      common, "level", "base", "top_genes", "padj_cutoff", "go_pvalue_cutoff",
      "go_qvalue_cutoff", "gsea_metric", "test", "ora_padj_cutoff",
      "ora_min_genes", "ora_lfc_cutoff", "run_spia"
    ),
    fgsea = c(common, "gmt_file", "padj_cutoff", "gsea_metric", "test"),
    decouple = c(
      common, "level", "base", "padj_cutoff",
      "run_decouple", "decouple_resources", "decouple_method",
      "decouple_min_size",
      # NEW: |log2FC| floor used by the per-TF regulon volcano and the
      # "Activé / Réprimé / Autre" categorisation in tf_regulon_summary.csv.
      # Adding it here means changing deg_lfc_cutoff invalidates the
      # decouple stage and regenerates those plots with the new floor.
      "deg_lfc_cutoff"
    ),
    isoform = c(
      common, "data_dir", "sample_table", "count_type", "matrix_file",
      "subset_sample", "remove_sample", "custom_tx2gene", "custom_gene_map",
      "custom_transcript_id_map", "isoform_fasta", "isoform_gff", "model",
      "level", "base", "test", "reduced", "padj_cutoff", "run_isoform",
      "run_dexseq", "run_isoform_enrichment", "ora_padj_cutoff",
      "ora_min_genes", "ora_lfc_cutoff", "gsea_metric", "gmt_file",
      "top_genes", "go_pvalue_cutoff", "go_qvalue_cutoff", "run_predictors",
      "predictor_cpu", "use_wsl", "wsl_distro", "run_spia",
      "isoform_test_engine", "resume_isoform_from", "batch_col", "pca_ntop",
      "pca_colors", "sample_labels", "skip_fasta_filter", "plot_topology",
      "isoform_report_genes", "isoform_plot_top_n", "bpparam", "make_report"
    ),
    report = c(
      "out_dir", "run_dge", "level", "base", "report_settings", "make_report"
    ),
    rdata = c(
      common, "level", "base", "model", "test", "batch_col", "run_dge",
      "run_isoform"
    )
  )
}

# Values .expressom_prepare() derives from several parameters. A stage keeps
# only the ones it reads: otherwise changing `level` would change the context of
# the import stage (comp_name) and needlessly re-import the counts.
.stage_derived_values <- function() {
  list(
    import   = character(),
    eda      = c("main_condition", "model_eda"),
    dge      = character(),
    dge_plots = "main_condition",
    functional = character(),
    fgsea    = "comp_name",
    decouple = "comp_name",
    isoform  = "main_condition",
    report   = character(),
    rdata    = c("main_condition", "comp_name")
  )
}

#' Restrict a stage context to the parameters its stage uses
#' @keywords internal
.stage_ctx <- function(ctx, stage) {
  keep <- .stage_param_groups()[[stage]]
  ctx$p <- ctx$p[intersect(keep, names(ctx$p))]
  derived <- c("main_condition", "model_eda", "comp_name")
  for (nm in setdiff(derived, .stage_derived_values()[[stage]])) {
    ctx[nm] <- list(NULL)
  }
  ctx
}

#' Pipeline stages: id, label and the targets they correspond to
#'
#' @return A data.frame with columns `stage`, `label`, `target` and `upstream`
#'   (space-separated stages it reads from).
#' @export
expressom_pipeline_stages <- function() {
  data.frame(
    stage = c("setup", "import", "eda", "dge", "dge_plots", "functional",
              "fgsea", "decouple", "isoform", "report", "rdata"),
    label = c("Validate parameters", "Import counts", "EDA / PCA",
              "DESeq2 differential expression", "DGE figures",
              "ORA / GSEA / SPIA", "FGSEA",
              "TF activity (decoupleR)", "Isoforms (DTE / DTU / switch)",
              "Interactive report", "Save .RData"),
    target = c("setup", "import", "eda", "dge", "dge_plots", "functional",
               "fgsea", "decouple", "isoform", "report", "rdata"),
    upstream = c("", "setup", "import", "import", "dge", "dge", "dge", "dge",
                 "setup", "eda dge_plots functional fgsea decouple isoform",
                 "import dge functional decouple isoform"),
    stringsAsFactors = FALSE
  )
}

# which optional stages exist for these parameters (mirrors .expressom_prepare())
.pipeline_flags <- function(p) {
  full <- !isTRUE(p$eda_only) && !is.null(p$model) && !is.null(p$level) &&
    !is.null(p$base)
  list(
    run_dge      = full && isTRUE(p$run_dge),
    run_decouple = full && isTRUE(p$run_dge) && isTRUE(p$run_decouple),
    run_isoform  = full && isTRUE(p$run_isoform),
    make_report  = isTRUE(p$make_report)
  )
}

# one hash per stage, over the parameters that stage uses (sorted by name)
.stage_hashes <- function(p) {
  vapply(.stage_param_groups(), function(keep) {
    used <- sort(intersect(keep, names(p)))
    .object_md5(p[used])
  }, character(1))
}

# ---- targets definition -------------------------------------------------------------------------------

.fp_slim <- function(fp) {
  # path + size + md5: a touched-but-identical file must not invalidate anything
  if (is.data.frame(fp)) fp[, c("path", "size", "md5"), drop = FALSE] else fp
}

.fp_counts <- function(p) {
  .fp_slim(.file_fingerprint(c(
    p$sample_table, p$matrix_file, p$custom_tx2gene, p$custom_gene_map,
    .quant_input_files(p, gene_level = TRUE)
  )))
}

.fp_labels <- function(p) {
  v <- p$sample_labels
  .fp_slim(if (is.character(v) && length(v) == 1L) .file_fingerprint(v) else character())
}

.fp_isoform <- function(p) {
  paths <- c(p$isoform_fasta, p$isoform_gff, .quant_input_files(p))
  if (is.character(p$custom_transcript_id_map) &&
      length(p$custom_transcript_id_map) == 1L) {
    paths <- c(paths, p$custom_transcript_id_map)
  }
  .fp_slim(.file_fingerprint(paths))
}

.pipeline_setup <- function(p, fp_counts, fp_labels) {
  ctx <- .expressom_prepare(p)
  .write_run_config(ctx$p$out_dir, ctx$p)
  ctx
}

# The functions the target commands call. They are copied into the environment
# targets evaluates in, so the generated _targets.R only needs
# ExpressOM::expressom_pipeline() -- no `:::` and no library(ExpressOM) -- and
# targets hashes their code: reinstalling a package with an edited stage
# *wrapper* re-runs that stage (the code those wrappers call is not tracked, see
# the header of this file). `overrides` (a named list of functions) replaces any
# of them; the test-suite uses it to run the DAG with cheap stand-ins.
.pipeline_functions <- function() {
  c(
    ".stage_ctx", ".fp_counts", ".fp_labels", ".fp_isoform", ".pipeline_setup",
    ".predictor_environment_fingerprint",
    ".stage_import", ".stage_eda", ".stage_dge", ".stage_dge_plots",
    ".stage_functional", ".stage_fgsea", ".stage_decouple",
    ".stage_isoform", ".stage_report", ".stage_save_rdata"
  )
}

.pipeline_env <- function(overrides = list()) {
  ns <- environment(.pipeline_env)
  env <- new.env(parent = globalenv())
  for (nm in .pipeline_functions()) {
    assign(nm, get(nm, envir = ns), envir = env)
  }
  for (nm in names(overrides)) {
    assign(nm, overrides[[nm]], envir = env)
  }
  env
}

#' The ExpressOM `targets` pipeline
#'
#' Returns the list of targets for the parameters stored in `params_file`.
#' It is what the generated `_targets.R` calls; you normally do not call it
#' yourself -- use [run_expressom_targets()] or [run_expressom_app()].
#'
#' The stages form a DAG: `setup` -> `import` -> `eda` / `dge` -> `dge_plots`,
#' `functional`, `fgsea`, `decouple`; `isoform` is independent of `dge`;
#' `report` and `rdata` close the graph. Each stage target depends on a small
#' context target (`ctx_<stage>`) holding only the parameters it uses, so a
#' parameter that concerns one stage leaves the cached results of the others
#' untouched. When `run_predictors = TRUE`, the `isoform` target additionally
#' depends on `fp_predictors` (the external predictor environment fingerprint),
#' so it re-runs when a predictor tool or database is installed or removed even
#' though no parameter or input file changed.
#'
#' @param params_file Path to the `params.rds` written by
#'   [write_expressom_targets()].
#' @param overrides Named list of functions replacing internal pipeline
#'   functions (for tests).
#' @return A list of `targets` target objects.
#' @export
expressom_pipeline <- function(params_file, overrides = list()) {
  .rep_require("targets", "to run the pipeline with targets")
  p0 <- readRDS(params_file)
  flags <- .pipeline_flags(p0)

  # Every WSL/predictor log must land under <out_dir>/log/WSL/ instead of the
  # targets working directory. Set it here (in the child process that
  # tar_make() spins up) so all downstream .wsl_exec_script() calls -- target
  # code AND .predictor_environment_fingerprint() probes -- write to log/WSL.
  if (!is.null(p0$out_dir) && nzchar(p0$out_dir)) {
    wsl_log_dir <- file.path(.expressom_log_dir(p0$out_dir, create = TRUE), "WSL")
    if (!dir.exists(wsl_log_dir)) {
      dir.create(wsl_log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    options(expressom.wsl_log_dir = wsl_log_dir)
  }

  targets::tar_option_set(
    envir = .pipeline_env(overrides),
    packages = character(0),
    format = "rds",
    memory = "transient",
    garbage_collection = TRUE,
    error = "stop"
  )

  tgt <- function(name, command, ...) {
    targets::tar_target_raw(name, command, ...)
  }
  ctx_target <- function(stage) {
    tgt(paste0("ctx_", stage), bquote(.stage_ctx(setup, .(stage))))
  }
  always <- targets::tar_cue(mode = "always")

  out <- list(
    tgt("params_file", params_file, format = "file"),
    tgt("params", quote(readRDS(params_file))),
    # fingerprints of the *content* of the input files, re-computed each run
    tgt("fp_counts", quote(.fp_counts(params)), cue = always),
    tgt("fp_labels", quote(.fp_labels(params)), cue = always),
    tgt("setup", quote(.pipeline_setup(params, fp_counts, fp_labels))),
    ctx_target("import"),
    ctx_target("eda"),
    tgt("import", quote({
      fp_counts
      .stage_import(ctx_import)
    })),
    tgt("eda", quote(.stage_eda(ctx_eda, import))),
    tgt("m_eda", quote(eda$files))
  )
  manifests <- "m_eda"

  if (flags$run_dge) {
    out <- c(out, list(
      ctx_target("dge"),
      ctx_target("dge_plots"),
      ctx_target("functional"),
      ctx_target("fgsea"),
      tgt("dge", quote(.stage_dge(ctx_dge, import))),
      tgt("dge_plots", quote(.stage_dge_plots(ctx_dge_plots, import, dge))),
      tgt("functional", quote(.stage_functional(ctx_functional, dge))),
      tgt("fgsea", quote(.stage_fgsea(ctx_fgsea, dge))),
      tgt("m_dge", quote(dge$files)),
      tgt("m_dge_plots", quote(dge_plots$files)),
      tgt("m_functional", quote(functional$files)),
      tgt("m_fgsea", quote(fgsea$files))
    ))
    manifests <- c(manifests, "m_dge", "m_dge_plots", "m_functional", "m_fgsea")
  }

  if (isTRUE(flags$run_decouple)) {
    out <- c(out, list(
      ctx_target("decouple"),
      tgt("decouple", quote(.stage_decouple(ctx_decouple, dge))),
      tgt("m_decouple", quote(decouple$files))
    ))
    manifests <- c(manifests, "m_decouple")
  }

  if (flags$run_isoform) {
    isoform_deps <- list(
      tgt("fp_isoform", quote(.fp_isoform(params)), cue = always),
      ctx_target("isoform")
    )
    # Re-probe the predictor environment on every run if predictors are
    # enabled (~1-2s WSL round-trip). Without this, "isoform" stays cached
    # even after a new predictor tool / database was installed, and the
    # pipeline reuses the old, predictor-empty checkpoints. Combined with
    # `predictor_env` now being part of `.stage_isoform()`'s checkpoint
    # signature, this is what makes "install a tool → re-run → predictor
    # results appear" work end-to-end.
    #
    # The explicit log_dir makes the WSL commands it issues land under
    # log/WSL/ (matching the option set at the top of expressom_pipeline()),
    # so the fingerprint probe shows up in the report's WSL log section.
    if (isTRUE(p0$run_predictors)) {
      isoform_deps <- c(isoform_deps, list(
        tgt(
          "fp_predictors",
          quote(.predictor_environment_fingerprint(
            wsl_distro = params$wsl_distro,
            use_wsl    = params$use_wsl,
            log_dir    = file.path(.expressom_log_dir(params$out_dir, create = TRUE), "WSL")
          )),
          cue = always
        )
      ))
    }
    isoform_body <- if (isTRUE(p0$run_predictors)) {
      quote({
        fp_counts
        fp_isoform
        fp_predictors
        .stage_isoform(ctx_isoform)
      })
    } else {
      quote({
        fp_counts
        fp_isoform
        .stage_isoform(ctx_isoform)
      })
    }
    out <- c(out, c(isoform_deps, list(
      tgt("isoform", isoform_body),
      tgt("m_isoform", quote(isoform$files))
    )))
    manifests <- c(manifests, "m_isoform")
  }

  if (flags$make_report) {
    report_cmd <- as.call(c(
      list(quote(.stage_report), quote(ctx_report)),
      lapply(manifests, as.name)
    ))
    # a failing report must not keep the .RData from being written (nor the
    # reverse): both are leaves of the graph
    out <- c(out, list(
      ctx_target("report"),
      tgt("report", report_cmd, format = "file", error = "continue")
    ))
  }

  if (flags$run_dge || flags$run_isoform) {
    arg <- function(target, present) if (present) as.name(target) else quote(NULL)
    rdata_cmd <- call(
      ".stage_save_rdata", quote(ctx_rdata), quote(import),
      arg("dge",      flags$run_dge),
      arg("functional", flags$run_dge),
      arg("isoform",  flags$run_isoform),
      arg("decouple", flags$run_decouple)
    )
    out <- c(out, list(
      ctx_target("rdata"),
      tgt("rdata", rdata_cmd, format = "file", error = "continue")
    ))
  }
  out
}

# ---- writing / running -----------------------------------------------------------------------------------

#' Default directory of the pipeline for an output directory
#' @param out_dir Output directory of the run.
#' @return `<out_dir>/.expressom_targets`.
#' @export
expressom_targets_dir <- function(out_dir) {
  file.path(out_dir, ".expressom_targets")
}

#' Write the `targets` pipeline for a set of parameters
#'
#' Stores the parameters (`params.rds`) and generates the one-line
#' `_targets.R`. Re-writing identical parameters leaves the store untouched, so
#' nothing is invalidated; changed parameters invalidate exactly the stages
#' that use them.
#'
#' @param params Parameter list from [expressom_params()].
#' @param dir Pipeline directory (default [expressom_targets_dir()] of
#'   `params$out_dir`).
#' @param script_lines Lines of `_targets.R`; by default the one-liner calling
#'   [expressom_pipeline()] (a custom script is used by the tests).
#' @return A list with `dir`, `script`, `store` and `params_file` (invisibly).
#' @export
write_expressom_targets <- function(params, dir = NULL, script_lines = NULL) {
  if (!is.list(params) || is.null(params$out_dir)) {
    stop("`params` must be a list from expressom_params() with an `out_dir`.", call. = FALSE)
  }
  params <- .absolutize_paths(params)
  dir <- dir %||% expressom_targets_dir(params$out_dir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir, winslash = "/", mustWork = TRUE)

  params_file <- file.path(dir, "params.rds")
  # only touch the file when the content changed: its hash is a target, but a
  # stable mtime also keeps `tar_make()` from having to re-read it
  old <- .read_rds_safe(params_file)
  if (is.null(old) || !identical(.object_md5(old), .object_md5(params))) {
    saveRDS(params, params_file)
  }

  script <- file.path(dir, "_targets.R")
  if (is.null(script_lines)) {
    script_lines <- c(
      "# Generated by ExpressOM::write_expressom_targets() -- do not edit.",
      "# The parameters live in params.rds next to this file; change them with",
      "# expressom_params() + run_expressom_targets(), or in the Shiny app.",
      sprintf('ExpressOM::expressom_pipeline("%s")', params_file)
    )
  }
  writeLines(script_lines, script)

  invisible(list(
    dir = dir,
    script = script,
    store = file.path(dir, "store"),
    params_file = params_file
  ))
}

# targets (and their upstream) needed to build `stages`
.stage_targets <- function(stages) {
  tbl <- expressom_pipeline_stages()
  bad <- setdiff(stages, tbl$stage)
  if (length(bad)) {
    stop(
      "Unknown stage(s): ", paste(bad, collapse = ", "),
      ". Choose among: ", paste(tbl$stage, collapse = ", "), ".",
      call. = FALSE
    )
  }
  tbl$target[match(stages, tbl$stage)]
}

#' Run the pipeline with `targets`
#'
#' Writes the pipeline for `params` and runs `targets::tar_make()`. Stages whose
#' inputs did not change are skipped; the interactive report is rebuilt when
#' (and only when) one of the stages it reads from re-ran.
#'
#' @param params Parameter list from [expressom_params()]. `NULL` re-uses the
#'   parameters already stored in `dir`.
#' @param dir Pipeline directory (default `<out_dir>/.expressom_targets`).
#' @param stages Optional character vector of stages to build, among
#'   [expressom_pipeline_stages()]`$stage`. Each is built together with the
#'   stages it needs; nothing else runs (`stages = "dge"` stops after the
#'   DESeq2 stage, without functional analysis or report). `NULL` builds
#'   everything.
#' @param callr Run in a clean background R process (default; recommended).
#' @param reporter `targets` console reporter.
#' @param script_lines Custom `_targets.R` lines (see
#'   [write_expressom_targets()]); mostly for tests.
#' @param ... Further arguments passed to [targets::tar_make()].
#' @return The stage status table ([expressom_stage_status()]), invisibly. A
#'   stage that stops the pipeline (import, DESeq2 ...) raises an error; the
#'   report and the `.RData` are leaves that do not block each other, so if one
#'   of them fails the run completes with a warning naming it.
#' @export
#' @examples
#' \dontrun{
#' p <- expressom_params(
#'   data_dir = "data", out_dir = "results", sample_table = "samples.csv",
#'   level = "treated", base = "control", model = "~ condition"
#' )
#' run_expressom_targets(p)                      # first run: everything
#' p$padj_cutoff <- 0.01
#' run_expressom_targets(p)                      # only the stages using padj_cutoff
#' run_expressom_targets(p, stages = "dge")      # up to the DESeq2 stage only
#' }
run_expressom_targets <- function(
  params = NULL,
  dir = NULL,
  stages = NULL,
  callr = TRUE,
  reporter = "summary",
  script_lines = NULL,
  ...
) {
  .rep_require("targets", "to run the pipeline with targets")
  if (is.null(params)) {
    if (is.null(dir)) {
      stop("Give `params`, or the `dir` of an existing pipeline.", call. = FALSE)
    }
    params <- readRDS(file.path(dir, "params.rds"))
  }
  wanted <- if (is.null(stages)) NULL else .stage_targets(stages)

  paths <- write_expressom_targets(params, dir, script_lines)
  params <- readRDS(paths$params_file)

  # The child process that runs tar_make() inherits this option only if it
  # was set inside expressom_pipeline() (which is what actually matters for
  # the targets DAG). Setting it here too covers edge cases where a target
  # shells out before expressom_pipeline() has run (e.g. custom _targets.R).
  old_wsl_log <- getOption("expressom.wsl_log_dir")
  options(expressom.wsl_log_dir = file.path(.expressom_log_dir(params$out_dir, create = TRUE), "WSL"))
  on.exit(options(expressom.wsl_log_dir = old_wsl_log), add = TRUE)

  args <- list(
    script = paths$script,
    store = paths$store,
    reporter = reporter,
    callr_function = if (isTRUE(callr)) callr::r else NULL
  )
  if (!is.null(wanted)) {
    args$names <- bquote(tidyselect::any_of(.(wanted)))
  }
  args <- c(args, list(...))

  # A stage that stops the pipeline (error = "stop") makes tar_make() throw --
  # but the stages that completed BEFORE it are up to date in the store. Record
  # their parameter hashes either way, otherwise the app keeps flagging them
  # "stale" after a failed run and offers to re-run work that is already done.
  make_error <- NULL
  tryCatch(
    withr::with_dir(
      paths$dir,
      do.call(targets::tar_make, args)
    ),
    error = function(e) make_error <<- e
  )

  .record_stage_hashes(params, paths$dir)
  if (!is.null(make_error)) {
    stop(make_error)
  }
  status <- expressom_stage_status(params, dir = paths$dir)

  # stages that may fail without stopping the pipeline (the report and the
  # .RData are leaves): say so instead of finishing silently
  failed <- status[status$state == "errored", , drop = FALSE]
  if (nrow(failed)) {
    meta <- tryCatch(
      targets::tar_meta(store = paths$store, fields = "error"),
      error = function(e) NULL
    )
    tg <- expressom_pipeline_stages()$target[match(failed$stage, expressom_pipeline_stages()$stage)]
    msgs <- vapply(tg, function(t) {
      m <- if (is.null(meta)) NA_character_ else as.character(meta$error[match(t, meta$name)])
      if (length(m) && !is.na(m) && nzchar(m)) paste0(": ", m) else ""
    }, character(1))
    warning(
      "These stages failed: ", paste0(failed$label, msgs, collapse = "; "),
      ". The other stages completed; fix the cause and run again.",
      call. = FALSE
    )
  }
  invisible(status)
}

# After a make: remember, per stage, the parameter hash it is now consistent
# with. Only stages the make touched (completed or skipped = up to date) get
# their hash refreshed; the others stay as they were, so they keep showing as
# stale in the app.
.record_stage_hashes <- function(params, dir) {
  hashes_file <- file.path(dir, "stage_hashes.rds")
  old <- .read_rds_safe(hashes_file, default = list())
  store <- file.path(dir, "store")
  prog <- tryCatch(targets::tar_progress(store = store), error = function(e) NULL)
  if (is.null(prog) || !nrow(prog)) {
    return(invisible(old))
  }
  tbl <- expressom_pipeline_stages()
  cur <- .stage_hashes(params)
  for (st in names(cur)) {
    tg <- tbl$target[tbl$stage == st]
    ok <- prog$progress[prog$name == tg]
    if (length(ok) && ok[[1L]] %in% c("completed", "skipped")) {
      old[[st]] <- cur[[st]]
    }
  }
  saveRDS(old, hashes_file)
  invisible(old)
}

#' Status of every pipeline stage
#'
#' For each stage: whether it ran, how long it took, and whether the current
#' parameters differ from the ones it last ran with. Staleness is propagated
#' downstream: a stage is stale when its own parameters changed or any stage it
#' reads from is stale.
#'
#' @param params Current parameter list (to compute staleness); `NULL` uses the
#'   parameters stored in `dir`.
#' @param dir Pipeline directory.
#' @param running Is a run in progress right now? `NULL` (default) trusts the
#'   `targets` progress file; `FALSE` treats stages left "running" by a killed
#'   process as errored (used by the app).
#' @return A data.frame with `stage`, `label`, `progress`, `seconds`,
#'   `finished` and `state` ("never run", "current", "stale", "running",
#'   "errored"). Only the stages that exist for these parameters are listed.
#' @export
expressom_stage_status <- function(params = NULL, dir, running = NULL) {
  stages <- expressom_pipeline_stages()
  store <- file.path(dir, "store")

  if (is.null(params)) {
    pf <- file.path(dir, "params.rds")
    params <- if (file.exists(pf)) readRDS(pf) else list()
  }
  flags <- .pipeline_flags(params)

  prog <- NULL
  meta <- NULL
  if (dir.exists(store) && requireNamespace("targets", quietly = TRUE)) {
    prog <- tryCatch(targets::tar_progress(store = store), error = function(e) NULL)
    meta <- tryCatch(
      targets::tar_meta(store = store, fields = c("time", "seconds", "error")),
      error = function(e) NULL
    )
  }
  pick <- function(df, col, tg, default) {
    if (is.null(df) || !tg %in% df$name) default else df[[col]][match(tg, df$name)]
  }

  stages$progress <- vapply(
    stages$target, function(t) as.character(pick(prog, "progress", t, "none")), character(1)
  )
  stages$seconds <- vapply(
    stages$target, function(t) as.numeric(pick(meta, "seconds", t, NA_real_)), numeric(1)
  )
  stages$finished <- vapply(stages$target, function(t) {
    tm <- pick(meta, "time", t, NA)
    if (length(tm) != 1L || is.na(tm)) NA_character_ else format(tm, "%Y-%m-%d %H:%M:%S")
  }, character(1))
  stages$errored <- vapply(stages$target, function(t) {
    e <- pick(meta, "error", t, NA_character_)
    !is.na(e) && nzchar(e)
  }, logical(1))
  stages$built <- stages$target %in% (if (is.null(meta)) character() else meta$name[is.na(meta$error)])

  # stages that do not exist for these parameters
  exists <- stats::setNames(rep(TRUE, nrow(stages)), stages$stage)
  exists[c("dge", "dge_plots", "functional", "fgsea")] <- flags$run_dge
  exists["decouple"] <- flags$run_decouple
  exists["isoform"]  <- flags$run_isoform
  exists["report"]   <- flags$make_report
  exists["rdata"]    <- flags$run_dge || flags$run_isoform
  stages <- stages[exists[stages$stage], , drop = FALSE]

  last <- .read_rds_safe(file.path(dir, "stage_hashes.rds"), default = list())
  now <- .stage_hashes(params)
  # "setup" has no parameter group of its own: it re-validates on any change
  # and is never reported as stale by itself
  stale <- stats::setNames(
    vapply(stages$stage, function(s) {
      s %in% names(now) && !identical(last[[s]], now[[s]])
    }, logical(1)),
    stages$stage
  )
  for (i in seq_len(nrow(stages))) {
    ups <- intersect(strsplit(stages$upstream[[i]], " ", fixed = TRUE)[[1L]], stages$stage)
    if (length(ups) && any(stale[ups])) {
      stale[[stages$stage[[i]]]] <- TRUE
    }
  }

  state <- unname(ifelse(stale[stages$stage], "stale", "current"))
  state[!stages$built] <- "never run"
  state[stages$errored | stages$progress %in% c("errored", "canceled")] <- "errored"
  active <- stages$progress %in% c("dispatched", "running")
  state[active] <- if (isFALSE(running)) "errored" else "running"
  stages$state <- state
  rownames(stages) <- NULL
  stages[, c("stage", "label", "progress", "seconds", "finished", "state")]
}

#' Force stages to re-run
#'
#' Marks the targets of `stages` as outdated so the next
#' [run_expressom_targets()] rebuilds them (and the report) even though their
#' inputs did not change -- for instance after editing an input file in place
#' or deleting figures.
#'
#' @param dir Pipeline directory.
#' @param stages Stages to invalidate (see [expressom_pipeline_stages()]).
#' @return `NULL` (invisibly).
#' @export
expressom_targets_invalidate <- function(dir, stages) {
  .rep_require("targets", "to run the pipeline with targets")
  tg <- .stage_targets(stages)
  tg <- c(tg, paste0("m_", stages))
  store <- file.path(dir, "store")
  if (dir.exists(store)) {
    targets::tar_invalidate(
      names = tidyselect::any_of(tg),
      store = store
    )
  }
  hashes_file <- file.path(dir, "stage_hashes.rds")
  if (file.exists(hashes_file)) {
    h <- readRDS(hashes_file)
    h[stages] <- NULL
    saveRDS(h, hashes_file)
  }
  invisible(NULL)
}