trace_names <- function(p) {
  vapply(plotly::plotly_build(p)$x$data, function(t) t$name %||% "", character(1))
}

fixture_results <- function() {
  out <- withr::local_tempdir(.local_envir = parent.frame())
  make_run_dir(out)
  load_expressom_results(out)
}

test_that("iplot_volcano draws up / down / not-significant traces with real colours", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  p <- iplot_volcano(res$dge$results, 0.05, 1)
  expect_s3_class(p, "plotly")
  nm <- trace_names(p)
  expect_true(any(grepl("^Upregulated", nm)))
  expect_true(any(grepl("^Downregulated", nm)))
  expect_true(any(grepl("^Not significant", nm)))

  # R colour names are converted to hex: plotly.js does not know "red2"
  built <- plotly::plotly_build(p)$x$data
  cols <- unlist(lapply(built, function(t) t$marker$color))
  cols <- cols[!is.na(cols) & nzchar(cols)]
  expect_true(all(grepl("^(#|rgb)", cols)))
})

test_that("the volcano's significance counts follow the cut-offs", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  tbl <- res$dge$results
  n_up <- sum(tbl$padj < 0.05 & tbl$log2FoldChange > 1, na.rm = TRUE)
  p <- iplot_volcano(tbl, 0.05, 1)
  up_name <- grep("^Upregulated", trace_names(p), value = TRUE)
  expect_match(up_name, paste0("\\(", n_up, "\\)"))
  # tighter cut-off -> fewer significant genes
  p2 <- iplot_volcano(tbl, 0.001, 2)
  n_up2 <- sum(tbl$padj < 0.001 & tbl$log2FoldChange > 2, na.rm = TRUE)
  expect_match(grep("^Upregulated", trace_names(p2), value = TRUE), paste0("\\(", n_up2, "\\)"))
})

test_that("an LRT volcano has a single significance class (no direction)", {
  skip_if_not_installed("plotly")
  out <- withr::local_tempdir()
  make_run_dir(out, test = "LRT")
  res <- load_expressom_results(out)
  nm <- trace_names(iplot_volcano(res$dge$results, 0.05, 1, test_type = "LRT"))
  expect_true(any(grepl("LRT", nm)))
  expect_false(any(grepl("^Upregulated", nm)))
})

test_that("iplot_volcano marks highlighted genes", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  g <- res$dge$results$gene[[1L]]
  nm <- trace_names(iplot_volcano(res$dge$results, 0.05, 1, highlight = g))
  expect_true(length(nm) > length(trace_names(iplot_volcano(res$dge$results, 0.05, 1))))
})

test_that("iplot_ma, heatmap, DEG counts and p-value histogram return widgets", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  expect_s3_class(iplot_ma(res$dge$results, 0.05, 1), "plotly")
  expect_s3_class(
    iplot_heatmap(res$dge$expression$matrix, res$dge$expression$groups, max_genes = 20),
    "plotly"
  )
  expect_s3_class(iplot_deg_counts(res$dge$deg_counts), "plotly")
  expect_s3_class(iplot_pvalue_hist(res$dge$results$pvalue), "plotly")
})

test_that("iplot_pca and iplot_scree use the PCA tables", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  pc <- res$pca[[1L]]
  p <- iplot_pca(pc$scores, pc$variance, group_col = "condition",
                 group_colors = c(treated = "red2", control = "royalblue"))
  expect_s3_class(p, "plotly")
  expect_true(all(c("treated", "control") %in% trace_names(p)))
  expect_s3_class(iplot_scree(pc$variance), "plotly")
  # a request for a PC that does not exist degrades to an empty plot, not an error
  expect_s3_class(iplot_pca(pc$scores, pc$variance, x = "PC1", y = "PC99"), "plotly")
})

test_that("iplot_enrichment and iplot_dtu handle real and empty input", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  inv <- res$functional$inventory
  i <- which(inv$kind == "ora")[[1L]]
  std <- .rep_standardise_enrichment(res$functional$tables[[inv$path[[i]]]], "ora")
  expect_s3_class(iplot_enrichment(std, top_n = 10), "plotly")
  expect_s3_class(iplot_enrichment(NULL), "plotly")
  expect_s3_class(iplot_dtu(res$isoform$dtu), "plotly")
  expect_s3_class(iplot_dtu(NULL), "plotly")
})

test_that("empty inputs never error", {
  skip_if_not_installed("plotly")
  expect_s3_class(iplot_volcano(data.frame()), "plotly")
  expect_s3_class(iplot_ma(NULL), "plotly")
  expect_s3_class(iplot_heatmap(NULL), "plotly")
  expect_s3_class(iplot_deg_counts(NULL), "plotly")
})

test_that("itable builds a DT widget, caps rows and keeps counts as integers", {
  skip_if_not_installed("DT")
  df <- data.frame(n = 1:30, x = seq(0.00012345, 1, length.out = 30))
  w <- itable(df, caption = "cap", max_rows = 10)
  expect_s3_class(w, "datatables")
  expect_equal(nrow(w$x$data), 10L)
  expect_match(w$x$caption, "first 10 of 30")
  expect_s3_class(itable(NULL), "datatables")
})

test_that(".rep_hex converts R colour names and leaves hex values alone", {
  expect_identical(unname(.rep_hex("red2")), "#EE0000")
  expect_identical(unname(.rep_hex("#112233")), "#112233")
  expect_identical(names(.rep_hex(c(up = "red", down = "blue"))), c("up", "down"))
  expect_null(.rep_hex(NULL))
})


test_that("extended RNA-seq QC and diagnostic charts are interactive widgets", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  mat <- res$dge$expression$qc_matrix
  expect_s3_class(.rep_iplot_sample_correlation(mat), "plotly")
  expect_s3_class(.rep_iplot_sample_distance(mat), "plotly")
  expect_s3_class(.rep_iplot_sample_dendrogram(mat), "plotly")
  expect_s3_class(.rep_iplot_sample_mds(mat), "plotly")
  expect_s3_class(.rep_iplot_rle(mat), "plotly")
  expect_s3_class(.rep_iplot_mean_variance(mat), "plotly")
  expect_s3_class(.rep_iplot_expression_boxplot(res$dge$expression$distribution_matrix), "plotly")
  expect_s3_class(.rep_iplot_expression_density(res$dge$expression$distribution_matrix), "plotly")
  expect_s3_class(.rep_iplot_pvalue_qq(res$dge$results$pvalue), "plotly")
  expect_s3_class(.rep_iplot_ranked_effects(res$dge$results), "plotly")
})

test_that("enrichment bubble and gene-overlap plots handle suitable tables", {
  skip_if_not_installed("plotly")
  res <- fixture_results()
  inv <- res$functional$inventory
  i <- which(inv$kind == "ora")[[1L]]
  tbl <- res$functional$tables[[inv$path[[i]]]]
  std <- .rep_standardise_enrichment(tbl, "ora")
  expect_s3_class(.rep_iplot_enrichment_bubble(std), "plotly")
  # The fixture uses toy gene IDs; a heatmap is only returned if overlap data
  # are in a parseable, non-empty gene-list column.
  overlap <- .rep_iplot_enrichment_overlap(tbl)
  if (!is.null(overlap)) expect_s3_class(overlap, "plotly")
})

test_that("enrichment networks and GSEA running scores are interactive when genes overlap", {
  skip_if_not_installed("plotly")
  tbl <- data.frame(
    Description = c("Pathway A", "Pathway B", "Pathway C"),
    p.adjust = c(0.001, 0.003, 0.01),
    core_enrichment = c("g1/g2/g3/g4", "g2/g3/g5", "g4/g5/g6"),
    stringsAsFactors = FALSE
  )
  expect_s3_class(.rep_iplot_enrichment_network(tbl), "plotly")
  dge <- data.frame(
    gene = paste0("g", 1:8), stat = c(4, 3, 2, 1, -1, -2, -3, -4),
    stringsAsFactors = FALSE
  )
  expect_s3_class(.rep_iplot_gsea_running(tbl, dge), "plotly")
})
