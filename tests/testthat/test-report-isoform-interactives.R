test_that("DTE interactive diagnostics use transcript-level columns", {
  skip_if_not_installed("plotly")
  dte <- data.frame(
    transcript_id = paste0("ENST", seq_len(50)),
    gene_symbol = rep(c("GENE_A", "GENE_B"), 25),
    baseMean = seq(10, 500, length.out = 50),
    log2FoldChange = seq(-3, 3, length.out = 50),
    pvalue = pmax(1e-8, seq(0.0001, 0.9, length.out = 50)),
    padj = pmax(1e-7, seq(0.0002, 0.95, length.out = 50)),
    stringsAsFactors = FALSE
  )
  expect_s3_class(iplot_dte_volcano(dte, padj_cutoff = 0.05), "plotly")
  expect_s3_class(iplot_dte_ma(dte, padj_cutoff = 0.05), "plotly")
  expect_s3_class(iplot_dte_ranked(dte, top_n = 15), "plotly")
  expect_s3_class(iplot_dte_abundance_significance(dte), "plotly")
  expect_s3_class(iplot_dte_volcano(data.frame()), "plotly")
})

test_that("DTU diagnostics use LR and adjusted p-values without inventing fold change", {
  skip_if_not_installed("plotly")
  dtu <- data.frame(
    gene = rep(c("GENE_A", "GENE_B", "GENE_C"), 20),
    feature_id = paste0("ENST", seq_len(60)),
    lr = seq(0.1, 80, length.out = 60),
    df = rep(c(1, 2, 3), 20),
    pvalue = seq(0.0001, 0.9, length.out = 60),
    adj_pvalue = seq(0.0002, 0.95, length.out = 60),
    stringsAsFactors = FALSE
  )
  expect_s3_class(iplot_dtu_lr_q(dtu, padj_cutoff = 0.05), "plotly")
  expect_s3_class(iplot_dtu_lr_hist(dtu), "plotly")
  expect_s3_class(iplot_dtu_lr_df(dtu), "plotly")
  expect_s3_class(iplot_dtu_gene_summary(dtu), "plotly")
  expect_s3_class(iplot_dtu_lr_q(data.frame()), "plotly")
})

test_that("DT tables default to all rows with persistent easy filters", {
  skip_if_not_installed("DT")
  d <- data.frame(id = seq_len(1200), group = rep(c("A", "B", "C"), 400))
  w <- itable(d)
  expect_s3_class(w, "datatables")
  expect_equal(nrow(w$x$data), 1200L)
  expect_identical(w$x$filter, "top")
  expect_true(isTRUE(w$x$options$stateSave))
  expect_equal(w$x$options$searchDelay, 250)
})
