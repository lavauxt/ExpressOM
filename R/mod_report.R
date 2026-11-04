# mod_report.R -- the ONE unified, self-contained, interactive HTML report.
#
# build_expressom_report() reads an ExpressOM output directory through the
# shared loader (mod_report_data.R), builds the tab tree (mod_report_tabs.R)
# and renders inst/rmd/expressom_report.Rmd to a single HTML file with plotly
# figures and DT tables in tabs. It replaces the RegionReport report and the
# separate DTE/DTU report.

#' Build the unified interactive HTML report
#'
#' Regroups every result of an ExpressOM run -- sample QC / PCA, differential
#' expression (volcano, MA, heatmap, searchable tables), functional enrichment
#' (ORA, GSEA, FGSEA, SPIA, transcription factors), isoform analysis (DTE, DTU,
#' DEXSeq, isoform-switch plots), static figures and run logs -- in **one
#' self-contained HTML file** with plotly figures and DT tables in tabs. The
#' report only *reads* what the pipeline already wrote; it never recomputes an
#' analysis, so it can be rebuilt in seconds (for instance with other cut-offs).
#'
#' Tabs for stages that were not run (no isoform analysis, EDA-only run ...)
#' are left out.
#'
#' @param out_dir Output directory of [expressom()] (its `out_dir` argument).
#' @param output_file Path of the HTML file to write. Default:
#'   `<out_dir>/Report/ExpressOM_report_<comparison>.html`.
#' @param comparison Comparison to report (`<level>_vs_<base>`); see
#'   [list_expressom_comparisons()]. `NULL` uses the most recent one.
#' @param settings List from [expressom_report_settings()] (title, number of
#'   points / genes / terms drawn, static-figure embedding, cut-off overrides).
#' @param results Optional pre-loaded [load_expressom_results()] object (skips
#'   re-reading the files).
#' @param quiet Suppress knitr / pandoc console output.
#' @return The path of the HTML report (invisibly).
#' @export
#' @examples
#' \dontrun{
#' expressom(out_dir = "results", ...)
#' build_expressom_report("results")
#' build_expressom_report(
#'   "results",
#'   settings = expressom_report_settings(padj_cutoff = 0.01, highlight_genes = c("TP53", "MYC"))
#' )
#' }
build_expressom_report <- function(
  out_dir,
  output_file = NULL,
  comparison = NULL,
  settings = expressom_report_settings(),
  results = NULL,
  quiet = TRUE
) {
  for (pkg in c("rmarkdown", "knitr", "plotly", "DT", "htmltools", "htmlwidgets")) {
    .rep_require(pkg, "to build the interactive report")
  }
  if (!isTRUE(rmarkdown::pandoc_available())) {
    stop(
      "Pandoc was not found, but it is needed to write a self-contained HTML ",
      "report. It is bundled with RStudio / Positron; otherwise install it ",
      "from https://pandoc.org/installing.html",
      call. = FALSE
    )
  }

  if (is.null(results)) {
    results <- load_expressom_results(out_dir, comparison = comparison)
  } else if (!inherits(results, "expressom_results")) {
    stop("`results` must come from load_expressom_results().", call. = FALSE)
  }
  out_dir <- results$out_dir

  if (is.null(output_file)) {
    stem <- if (is.na(results$comparison)) {
      "ExpressOM_report"
    } else {
      paste0("ExpressOM_report_", results$comparison)
    }
    output_file <- file.path(out_dir, "Report", paste0(stem, ".html"))
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  work <- tempfile("expressom_report_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  template <- file.path(work, "expressom_report.Rmd")
  file.copy(.expressom_rmd_path("expressom_report.Rmd"), template)
  results_rds <- file.path(work, "results.rds")
  settings_rds <- file.path(work, "settings.rds")
  saveRDS(results, results_rds)
  saveRDS(settings, settings_rds)

  # Render in a fresh environment whose parent is the package namespace, so the
  # template can call the (internal) tab builders without `:::`.
  render_env <- new.env(parent = environment(build_expressom_report))

  rmarkdown::render(
    input = template,
    output_file = basename(output_file),
    output_dir = normalizePath(dirname(output_file), winslash = "/"),
    intermediates_dir = work,
    knit_root_dir = work,
    params = list(results_rds = results_rds, settings_rds = settings_rds),
    envir = render_env,
    quiet = quiet
  )

  size_mb <- file.size(output_file) / 1024^2
  message(
    "Report written: ", output_file, " (", format(round(size_mb, 1), nsmall = 1), " MB)"
  )
  invisible(normalizePath(output_file, winslash = "/"))
}
