# Report open-speed behaviour: bounded table payloads, slimmed columns,
# precomputed box-plot statistics and the plotly.js bundle selection.

skip_if_not_installed("plotly")
skip_if_not_installed("DT")

local_bundle_state <- function(ok = NULL) {
  # isolate the session-level bundle memo between tests
  env <- ExpressOM:::.rep_plotly_bundle_env
  old <- env$cartesian_ok
  if (!is.null(ok)) env$cartesian_ok <- ok
  withr::defer(assign("cartesian_ok", old, envir = env), envir = parent.frame())
  env
}
test_that(".xp_prune_table records the columns it hides", {
  df <- data.frame(
    gene = c("g1", "g2"), baseMean = 1:2, log2FoldChange = c(1, -1),
    lfcSE = c(.1, .1), stat = c(9, -9), padj = c(.01, .02),
    geneID = c("A/B/C", "D/E/F"), core_enrichment = c("A;B;C", "D;E;F"),
    stringsAsFactors = FALSE
  )
  slim <- ExpressOM:::.xp_prune_table(df, list(slim_tables = TRUE))
  expect_false(any(c("lfcSE", "stat", "geneID", "core_enrichment") %in% colnames(slim)))
  expect_setequal(
    attr(slim, "xp_dropped_cols"),
    c("lfcSE", "stat", "geneID", "core_enrichment")
  )
  # the essential result columns survive
  expect_true(all(c("gene", "log2FoldChange", "padj") %in% colnames(slim)))

  # slimming disabled -> the table is returned untouched, without attribute
  full <- ExpressOM:::.xp_prune_table(df, list(slim_tables = FALSE))
  expect_identical(colnames(full), colnames(df))
  expect_null(attr(full, "xp_dropped_cols"))

  # a settings list that predates slim_tables keeps slimming on ( %||% TRUE )
  legacy <- ExpressOM:::.xp_prune_table(df, list())
  expect_false("geneID" %in% colnames(legacy))
})

test_that("itable documents hidden columns and capped rows in the caption", {
  df <- data.frame(
    gene = paste0("g", 1:30), padj = sort(stats::runif(30)),
    geneID = rep("A/B/C", 30), stringsAsFactors = FALSE
  )
  attr(df, "xp_dropped_cols") <- "geneID"
  w <- itable(df, caption = "DE results", max_rows = 10)
  expect_s3_class(w, "datatables")
  expect_match(w$x$caption, "DE results", fixed = TRUE)
  expect_match(w$x$caption, "columns geneID are hidden", fixed = TRUE)
  expect_match(w$x$caption, "first 10 of 30", fixed = TRUE)
  expect_equal(nrow(w$x$data), 10L)
  # the attribute must not leak into the widget payload
  expect_null(attr(w$x$data, "xp_dropped_cols"))
})

test_that("the default table cap keeps the biggest tables browser-sized", {
  res <- fixture_results()
  # 400-gene fixture: default settings embed everything ...
  expect_identical(expressom_report_settings()$table_rows, 5000L)
  # ... but a table larger than the cap is truncated by itable, not by the loader
  big <- res$dge$results[rep(seq_len(nrow(res$dge$results)), length.out = 6000), ]
  w <- itable(big, max_rows = expressom_report_settings()$table_rows)
  expect_equal(nrow(w$x$data), 5000L)
})

test_that("box-plot statistics are computed in R and travel, not the raw genes", {
  set.seed(1)
  v <- c(rexp(1900, 1 / 8), rep(-40, 80))          # 1980 values, 80 outliers
  bs <- ExpressOM:::.rep_box_stats(v)
  expect_equal(bs$median, stats::median(v))
  # boxplot.stats() hinges are fivenum() quartiles, not type-7 quantiles
  fn <- stats::fivenum(v)
  expect_equal(bs$q1, fn[[2L]])
  expect_equal(bs$q3, fn[[4L]])
  expect_true(all(bs$outliers < bs$lowerfence | bs$outliers > bs$upperfence))
  expect_lte(length(bs$outliers), 250L)

  # tiny samples fall back to NULL so the trace is skipped, not broken
  expect_null(ExpressOM:::.rep_box_stats(c(1, NA, 3)))
})

test_that("the p-value histogram ships 20 bins, not every p-value", {
  p <- c(stats::runif(60000), rep(0.0001, 5000))    # 65k values, small-p excess
  w <- plotly::plotly_build(iplot_pvalue_hist(p))
  tr <- w$x$data[[1L]]
  expect_identical(tr$type, "bar")
  expect_equal(length(tr$y), 20L)                     # 20 fixed bins over [0, 1]
  expect_equal(sum(unlist(tr$y)), 65000L)             # nothing lost outside the bins
  # the spike lands in the first bin
  expect_true(tr$y[[1L]] > 5000L)
  expect_match(tr$text[[1L]], "^\\[0.00, 0.05", fixed = FALSE)
})

test_that("the p-value Q-Q plot subsamples instead of shipping every test", {
  p <- stats::runif(60000)
  w <- plotly::plotly_build(ExpressOM:::.rep_iplot_pvalue_qq(p, max_points = 2000L))
  tr <- w$x$data[[1L]]
  expect_equal(length(tr$x), 2000L)
  # endpoints survive: the smallest and largest p-values are kept
  expect_equal(min(unlist(tr$y)), -log10(max(p)), tolerance = 1e-8)
})

test_that("QC box plots send five-number summaries, not every gene", {
  res <- fixture_results()
  mat <- res$dge$expression$distribution_matrix
  n_vals <- nrow(mat) * ncol(mat)

  p <- plotly::plotly_build(.rep_iplot_expression_boxplot(mat))
  traces <- p$x$data
  expect_gte(length(traces), 2L)
  expect_true(all(vapply(traces, function(t) !is.null(t$lowerfence), logical(1))))
  expect_true(all(vapply(traces, function(t) !is.null(t$median), logical(1))))
  # per-sample payload: summary + capped outliers, never the raw n_genes values
  n_points <- sum(vapply(traces, function(t) {
    y <- t$y
    if (is.null(y)) 0L else if (is.list(y)) length(unlist(y)) else length(y)
  }, numeric(1)))
  expect_lt(n_points, 0.2 * n_vals)

  # RLE keeps its zero reference line and gains the same summaries
  r <- plotly::plotly_build(.rep_iplot_rle(mat))
  expect_true(all(vapply(r$x$data, function(t) !is.null(t$lowerfence), logical(1))))
  expect_true(any(vapply(r$x$layout$shapes, function(s) identical(s$y0, 0),
                         logical(1))))
})

test_that("no report trace requires a plotly.js bundle larger than cartesian", {
  res <- fixture_results()
  inv <- res$functional$inventory
  ora_i <- which(inv$kind == "ora")[[1L]]
  gsea_i <- which(inv$kind == "gsea")[[1L]]
  std_ora <- .rep_standardise_enrichment(res$functional$tables[[inv$path[[ora_i]]]], "ora")
  std_gsea <- .rep_standardise_enrichment(res$functional$tables[[inv$path[[gsea_i]]]], "gsea")
  widgets <- list(
    iplot_volcano(res$dge$results, 0.05, 1),
    iplot_ma(res$dge$results, 0.05, 1),
    iplot_heatmap(res$dge$expression$matrix, max_genes = 20),
    iplot_pvalue_hist(res$dge$results$pvalue),
    .rep_iplot_mean_variance(res$dge$expression$qc_matrix),
    .rep_iplot_expression_density(res$dge$expression$distribution_matrix),
    iplot_dtu(res$isoform$dtu),
    iplot_enrichment(std_ora, top_n = 10),
    iplot_enrichment(std_gsea, top_n = 10),
    iplot_dte_volcano(res$isoform$dte)
  )
  types <- unique(unlist(lapply(widgets, function(w) {
    vapply(plotly::plotly_build(w)$x$data, function(t) t$type %||% "scatter",
           character(1))
  })))
  # the cartesian bundle (the smaller plotly.js build the reports embed)
  # implements exactly: scatter, bar, pie, box, heatmap, histogram,
  # histogram2d, histogram2dcontour, contour, scatterternary, violin
  cartesian <- c("scatter", "bar", "pie", "box", "heatmap", "histogram",
                 "histogram2d", "histogram2dcontour", "contour",
                 "scatterternary", "violin")
  expect_true(all(types %in% cartesian), info = paste(types, collapse = ", "))
})

test_that("the plotly.js bundle choice degrades safely without network", {
  env <- local_bundle_state(ok = FALSE)              # pretend the CDN is unreachable
  withr::local_options(expressom.plotly_bundle = "auto")
  expect_identical(ExpressOM:::.rep_plotly_bundle_mode(), "full")
  withr::local_options(expressom.plotly_bundle = "full")
  expect_identical(ExpressOM:::.rep_plotly_bundle_mode(), "full")

  # forced cartesian: one clear message, then silence, and still a report
  withr::local_options(expressom.plotly_bundle = "cartesian")
  env$notified <- NULL
  expect_message(
    expect_identical(ExpressOM:::.rep_plotly_bundle_mode(), "full"),
    "full bundle instead"
  )
  expect_identical(suppressMessages(ExpressOM:::.rep_plotly_bundle_mode()), "full")
})

test_that("a cached cartesian bundle swaps into every widget dependency", {
  skip_if_not_installed("htmltools")
  path <- ExpressOM:::.rep_cartesian_bundle_path()
  skip_if(is.null(path), "plotly bundle version cannot be read")
  # a local stand-in for the CDN download: the swap only needs the file
  writeLines("// plotly-cartesian stand-in", path)
  withr::defer(unlink(path), envir = parent.frame())
  local_bundle_state(ok = TRUE)
  withr::local_options(expressom.plotly_bundle = "cartesian")

  res <- fixture_results()
  p <- .rep_iplot_mean_variance(res$dge$expression$qc_matrix)
  deps <- p$dependencies
  expect_true(any(vapply(deps, function(d) identical(d$name, "plotly-cartesian"),
                         logical(1))))
  expect_false(any(vapply(deps, function(d) identical(d$name, "plotly-main"),
                          logical(1))))
  # the swapped dependency points at the cached local file
  swapped <- deps[[which(vapply(deps, function(d) identical(d$name, "plotly-cartesian"),
                                logical(1)))[[1L]]]]
  expect_identical(swapped$script, basename(path))
  expect_true(file.exists(file.path(swapped$src$file, swapped$script)))
})
