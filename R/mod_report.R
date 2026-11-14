#' Build the unified interactive HTML report
#'
#' Regroups every result of an ExpressOM run -- sample QC / PCA, differential
#' expression (volcano, MA, heatmap, searchable tables), functional enrichment
#' (ORA, GSEA, FGSEA, SPIA, transcription factors, decoupleR), isoform analysis
#' (DTE, DTU, DEXSeq, isoform-switch plots), static figures and run logs -- in
#' one interactive HTML file with plotly figures and DT tables in tabs.
#'
#' By default the report is rendered with `rmarkdown` + pandoc as a single
#' self-contained `.html` file written **at the root of `<out_dir>`** (i.e.
#' `<out_dir>/ExpressOM_report_<comparison>.html`), so it sits right next to
#' `log/`, `Plots/`, `DE_raw_results/` and every other result folder, and the
#' relative links it contains (`log/…`, `Plots/…`, `DE_raw_results/…`) resolve
#' without any further setup when the file is opened from its own directory.
#'
#' Pass `settings = expressom_report_settings(self_contained = FALSE)` for the
#' legacy mode that writes a folder with `lib/` and `figures/` next to the
#' HTML.
#'
#' @param out_dir Output directory of [expressom()] (its `out_dir` argument).
#' @param output_file Path of the HTML file to write. Default:
#'   `<out_dir>/ExpressOM_report_<comparison>.html`.
#' @param comparison Comparison to report (`<level>_vs_<base>`); see
#'   [list_expressom_comparisons()]. `NULL` uses the most recent one.
#' @param settings List from [expressom_report_settings()].
#' @param results Optional pre-loaded [load_expressom_results()] object.
#' @param engine Report engine: `"htmltools"` (default; fast, dependency-free,
#'   needs plotly + DT), `"pandoc"` (the legacy rmarkdown self-contained
#'   single file) or `"quarto"` (Quarto dashboard with Observable JS
#'   explorers, DuckDB-Wasm filtering and float32 binary payloads; needs the
#'   quarto CLI, see `vignette` and `engine = "quarto"` in README).
#' @param quiet Suppress console output.
#' @return The path of the HTML report (invisibly).
#' @export
build_expressom_report <- function(
  out_dir,
  output_file = NULL,
  comparison = NULL,
  settings = expressom_report_settings(),
  results = NULL,
  engine = c("auto", "htmltools", "pandoc", "quarto"),
  quiet = TRUE
) {
  engine <- match.arg(engine)
  if (identical(engine, "auto")) {
    engine <- if (isTRUE(settings$self_contained)) "pandoc" else "htmltools"
  }
  if (identical(engine, "quarto")) {
    .rep_require("htmltools", "to build the interactive report")
    .rep_require("jsonlite", "to build the Quarto dashboard report")
    if (is.null(results)) {
      results <- load_expressom_results(out_dir, comparison = comparison)
    } else if (!inherits(results, "expressom_results")) {
      stop("`results` must come from load_expressom_results().", call. = FALSE)
    }
    out_dir <- results$out_dir
    if (is.null(output_file)) {
      stem <- if (is.na(results$comparison)) "ExpressOM_report" else {
        paste0("ExpressOM_report_", results$comparison)
      }
      output_file <- file.path(out_dir, paste0(stem, "_quarto.html"))
    }
    path <- .xq_render_dashboard(results, settings, output_file, quiet = quiet)
    size_mb <- file.size(path) / 1024^2
    message(
      "Quarto dashboard report written: ", path, " (",
      format(round(size_mb, 1), nsmall = 1), " MB, self-contained)"
    )
    return(invisible(path))
  }
  for (pkg in c("htmltools", "htmlwidgets")) {
    .rep_require(pkg, "to build the interactive report")
  }
  pandoc_mode <- identical(engine, "pandoc")
  if (pandoc_mode) {
    for (pkg in c("rmarkdown", "knitr", "plotly", "DT")) {
      .rep_require(pkg, "to build the interactive report")
    }
    if (!isTRUE(rmarkdown::pandoc_available())) {
      stop(
        "Pandoc is required for the self-contained single-file report ",
        "(the default). Install pandoc (RStudio / Positron bundle it), or ",
        "pass settings = expressom_report_settings(self_contained = FALSE) ",
        "to produce a folder with the HTML + lib/ + figures/ instead.",
        call. = FALSE
      )
    }
  } else {
    for (pkg in c("plotly", "DT")) {
      .rep_require(pkg, "to build the interactive report")
    }
  }

  if (is.null(results)) {
    results <- load_expressom_results(out_dir, comparison = comparison)
  } else if (!inherits(results, "expressom_results")) {
    stop("`results` must come from load_expressom_results().", call. = FALSE)
  }
  out_dir <- results$out_dir

  # ---- report path: HTML now lives at the ROOT of out_dir ---------------
  if (is.null(output_file)) {
    stem <- if (is.na(results$comparison)) {
      "ExpressOM_report"
    } else {
      paste0("ExpressOM_report_", results$comparison)
    }
    output_file <- file.path(out_dir, paste0(stem, ".html"))
  }
  report_dir <- dirname(output_file)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Configuration for the leaves: figures go to `figures/` next to the HTML in
  # external mode, or are inlined as base64 (pandoc mode).
  # ---------------------------------------------------------------------------
  st <- settings
  st$figures_dir <- if (pandoc_mode) NULL else file.path(report_dir, "figures")
  if (!is.null(st$figures_dir)) {
    unlink(st$figures_dir, recursive = TRUE)
    dir.create(st$figures_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (pandoc_mode) {
    work <- tempfile("expressom_report_")
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)

    template <- file.path(work, "expressom_report.Rmd")
    writeLines(.embedded_template_text("expressom_report.Rmd"), template, useBytes = TRUE)

    results_rds  <- file.path(work, "results.rds")
    settings_rds <- file.path(work, "settings.rds")
    saveRDS(results, results_rds, compress = FALSE)
    saveRDS(settings, settings_rds, compress = FALSE)

    render_env <- new.env(parent = environment(build_expressom_report))
    rmarkdown::render(
      input = template,
      output_file = basename(output_file),
      output_dir = normalizePath(report_dir, winslash = "/"),
      intermediates_dir = work,
      knit_root_dir = work,
      params = list(results_rds = results_rds, settings_rds = settings_rds),
      envir = render_env,
      quiet = quiet
    )
  } else {
    tabs <- expressom_report_tabs(results, st)
    ui <- .xp_build_html(tabs, results, st)
    htmltools::save_html(ui, file = output_file, libdir = "lib", lang = "en")
  }

  size_mb <- file.size(output_file) / 1024^2
  message(
    "Report written: ", output_file, " (", format(round(size_mb, 1), nsmall = 1),
    " MB",
    if (!pandoc_mode) " + figures/ and lib/" else ", self-contained",
    ")"
  )
  invisible(normalizePath(output_file, winslash = "/"))
}