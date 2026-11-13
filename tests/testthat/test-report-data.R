test_that("list_expressom_comparisons finds every comparison in a run", {
  out <- withr::local_tempdir()
  make_run_dir(out, level = "treated", base = "control")
  make_run_dir(out, level = "drugB", base = "control", functional = FALSE, isoform = FALSE)
  cmp <- list_expressom_comparisons(out)
  expect_setequal(cmp, c("treated_vs_control", "drugB_vs_control"))
})

test_that("list_expressom_comparisons is empty for an empty / EDA-only directory", {
  out <- withr::local_tempdir()
  expect_identical(length(list_expressom_comparisons(out)), 0L)
})

test_that("load_expressom_results returns the typed result contract", {
  out <- withr::local_tempdir()
  make_run_dir(out)
  res <- load_expressom_results(out)

  expect_s3_class(res, "expressom_results")
  expect_named(
    res,
    c("out_dir", "comparison", "comparisons", "config", "dge", "pca",
      "functional", "isoform", "figures", "logs", "files", "notes"),
    ignore.order = TRUE
  )
  expect_identical(res$comparison, "treated_vs_control")

  # configuration comes from the run snapshot + DGE log
  expect_identical(res$config$level, "treated")
  expect_identical(res$config$base, "control")
  expect_identical(res$config$test, "Wald")
  expect_identical(res$config$model, "~ condition")
  expect_equal(res$config$padj_cutoff, 0.05)
  expect_equal(res$config$lfc_cutoff, 1)
  expect_identical(res$config$main_condition, "condition")

  # differential expression
  expect_s3_class(res$dge$results, "data.frame")
  expect_equal(nrow(res$dge$results), 400L)
  expect_true(all(c("gene", "log2FoldChange", "padj") %in% colnames(res$dge$results)))
  expect_type(res$dge$results$padj, "double")
  expect_s3_class(res$dge$deg_counts, "data.frame")
  expect_equal(nrow(res$dge$deg_counts), 2L)
  expect_true(length(res$dge$sig_lists) >= 1L)
  expect_true(is.matrix(res$dge$expression$matrix))
  expect_true(is.matrix(res$dge$expression$qc_matrix))
  expect_true(is.matrix(res$dge$expression$distribution_matrix))
  expect_equal(ncol(res$dge$expression$qc_matrix), 6L)
  expect_equal(ncol(res$dge$expression$distribution_matrix), 6L)
  expect_identical(res$dge$expression$samples, paste0("S", 1:6))
  # heatmap columns are labelled with the PCA's grouping factor
  expect_identical(unname(res$dge$expression$groups), rep(c("control", "treated"), each = 3))

  # PCA
  expect_length(res$pca, 1L)
  expect_identical(res$pca[[1L]]$pcs, paste0("PC", 1:6))
  expect_equal(unname(res$pca[[1L]]$variance[["PC1"]]), 55)

  # figures: PDFs/PNGs of the run, but not the isoform plots (shown in their own tab)
  expect_true(all(c("pdf", "png") %in% res$figures$type))
  expect_false(any(grepl("DTE_volcano", res$figures$name)))
  expect_true(nrow(res$files) > 20L)
  expect_true(is.character(res$logs$session_info))
})

test_that("the loader discovers every enrichment table and labels its source", {
  out <- withr::local_tempdir()
  make_run_dir(out)
  inv <- load_expressom_results(out)$functional$inventory

  expect_setequal(unique(inv$source), c("DGE", "DTE", "DTU"))
  expect_setequal(unique(inv$kind), c("ora", "gsea", "spia", "tf"))
  # the "_significant" FGSEA companion duplicates the full table and is skipped
  expect_false(any(grepl("_significant", inv$path)))
  expect_true("FGSEA" %in% inv$method)
  expect_setequal(inv$collection[inv$method == "FGSEA"], c("C2_CP_REACTOME", "H"))
  expect_true(all(inv$source[grepl("DTE_Enrichment", inv$path)] == "DTE"))
})

test_that("enrichment tables are standardised to term / padj / size / score", {
  out <- withr::local_tempdir()
  make_run_dir(out)
  res <- load_expressom_results(out)
  inv <- res$functional$inventory
  for (i in seq_len(nrow(inv))) {
    std <- .rep_standardise_enrichment(
      res$functional$tables[[inv$path[[i]]]], inv$kind[[i]]
    )
    expect_named(std, c("term", "padj", "size", "score"), info = inv$path[[i]])
    expect_false(is.unsorted(std$padj), info = inv$path[[i]])
  }
  expect_null(.rep_standardise_enrichment(NULL))
  expect_null(.rep_standardise_enrichment(data.frame(x = 1)))
})

test_that("an LRT run has no directional counts", {
  out <- withr::local_tempdir()
  make_run_dir(out, test = "LRT")
  res <- load_expressom_results(out)
  expect_identical(res$config$test, "LRT")
  expect_identical(res$config$reduced, "~ 1")
  expect_true(all(is.na(res$dge$deg_counts$n_up)))
})

test_that("stages that were not run leave their elements empty", {
  out <- withr::local_tempdir()
  make_run_dir(out, functional = FALSE, isoform = FALSE, figures = FALSE)
  res <- load_expressom_results(out)
  expect_equal(nrow(res$functional$inventory), 0L)
  expect_null(res$isoform)
  expect_equal(nrow(res$figures), 0L)
  expect_false(is.null(res$dge))
})

test_that("an EDA-only directory loads (no comparison, no DGE)", {
  out <- withr::local_tempdir()
  make_run_dir(out)
  unlink(file.path(out, c("DE_raw_results", "ORA", "GSEA", "SPIA", "Transcription_Factors", "IsoformSwitch", "log")), recursive = TRUE)
  file.rename(
    file.path(out, "Plots", "PCA_treated_vs_control_sample_scores.tsv"),
    file.path(out, "Plots", "PCA_EDA_sample_scores.tsv")
  )
  file.rename(
    file.path(out, "Plots", "PCA_treated_vs_control_variance_explained.tsv"),
    file.path(out, "Plots", "PCA_EDA_variance_explained.tsv")
  )
  res <- load_expressom_results(out)
  expect_true(is.na(res$comparison))
  expect_null(res$dge)
  expect_length(res$pca, 1L)
})

test_that("loader errors are explicit", {
  expect_error(load_expressom_results(file.path(tempdir(), "no-such-dir")), "existing ExpressOM output directory")
  out <- withr::local_tempdir()
  make_run_dir(out)
  expect_error(load_expressom_results(out, comparison = "x_vs_y"), "not found")
})

test_that("the table reader is tolerant", {
  d <- withr::local_tempdir()
  expect_null(.rep_read_table(file.path(d, "missing.csv")))
  empty <- file.path(d, "empty.csv")
  file.create(empty)
  expect_null(.rep_read_table(empty))
  tsv <- file.path(d, "a.tsv")
  writeLines(c("a\tb", "1\t2", "3\t4"), tsv)
  expect_equal(nrow(.rep_read_table(tsv)), 2L)
  expect_equal(nrow(.rep_read_table(tsv, max_rows = 1)), 1L)
})

test_that("print.expressom_results summarises the run", {
  out <- withr::local_tempdir()
  make_run_dir(out)
  expect_output(print(load_expressom_results(out)), "treated_vs_control")
})

test_that("the run snapshot keeps simple values and drops what cannot be serialised", {
  out <- withr::local_tempdir()
  .write_run_config(out, list(
    level = "a", base = "b", padj_cutoff = 0.01,
    flag = TRUE,
    named = c(up = "red", down = "blue"),
    fun = function(x) x,
    df = data.frame(x = 1:3),
    long = seq_len(500),
    nothing = NULL
  ))
  cfg <- jsonlite::read_json(file.path(out, "log", "expressom_run.json"), simplifyVector = TRUE)
  expect_identical(cfg$level, "a")
  expect_equal(cfg$padj_cutoff, 0.01)
  expect_true(cfg$flag)
  expect_false("fun" %in% names(cfg))
  expect_false("df" %in% names(cfg))
  expect_false("long" %in% names(cfg))
  expect_true(!is.null(cfg$written_at))
  # named vectors are stored as objects so their names survive
  expect_identical(unlist(cfg$named), c(up = "red", down = "blue"))
})

test_that("the colours chosen for a run are read back by the loader", {
  out <- withr::local_tempdir()
  make_run_dir(out, snapshot = FALSE)
  .write_run_config(out, list(
    level = "treated", base = "control", test = "Wald",
    volcano_colors = c(up = "orange", down = "navy", ns = "grey70"),
    pca_colors = c(level = "darkgreen", base = "purple")
  ))
  res <- load_expressom_results(out)
  expect_identical(unname(res$config$volcano_colors[["up"]]), "orange")
  expect_identical(unname(res$config$volcano_colors[["ns"]]), "grey70")
  expect_identical(unname(res$config$pca_colors[["base"]]), "purple")
})
