test_that("expressom_default_params mirrors the arguments of expressom()", {
  p <- expressom_default_params()
  expect_setequal(names(p), names(formals(expressom)))
  expect_identical(p$level, "treated")
  expect_identical(p$padj_cutoff, 0.01)
  # choice arguments resolve to their first choice, the live backend object is left out
  expect_identical(p$execution_order, "dge_first")
  expect_identical(p$isoform_test_engine, "DEXSeq")
  expect_null(p$bpparam)
  expect_true(is.character(p$gmt_file) && length(p$gmt_file) > 1L)
})

test_that("expressom_params overrides, validates names and makes paths absolute", {
  p <- expressom_params(level = "a", base = "b", out_dir = "some/rel/dir", padj_cutoff = 0.2)
  expect_identical(p$level, "a")
  expect_equal(p$padj_cutoff, 0.2)
  expect_true(grepl("^(/|[A-Za-z]:)", p$out_dir))
  expect_match(p$out_dir, "some/rel/dir$")
  # NULL is a legitimate override (e.g. "no batch column")
  expect_null(expressom_params(batch_col = NULL)$batch_col)

  expect_error(expressom_params(not_a_param = 1), "Unknown ExpressOM parameter")
  expect_error(expressom_params(1), "must be named")
  # a typo is reported, not silently ignored
  expect_error(expressom_params(padj_cutof = 0.1), "padj_cutof")
})

test_that("expressom_params starts from another parameter set", {
  base <- expressom_params(level = "x", base = "y")
  p <- expressom_params(top_genes = 5, defaults = base)
  expect_identical(p$level, "x")
  expect_equal(p$top_genes, 5)
})

test_that("expressom_params_from_run restores a run from its snapshot", {
  out <- withr::local_tempdir()
  make_run_dir(out, level = "drugB")
  p <- expressom_params_from_run(out)
  expect_identical(p$level, "drugB")
  expect_identical(p$test, "Wald")
  expect_equal(p$padj_cutoff, 0.05)
  expect_identical(p$count_type, "matrix")
  expect_equal(normalizePath(p$out_dir), normalizePath(out))
  # named colour vectors survive the JSON round trip
  expect_setequal(names(p$volcano_colors), c("up", "down", "ns"))
  # arguments that were not recorded keep their defaults
  expect_identical(p$shrink_method, "ashr")
  expect_identical(p$top_genes, expressom_default_params()$top_genes)

  expect_error(expressom_params_from_run(withr::local_tempdir()), "No run snapshot")
})

test_that("the stage table describes the DAG", {
  st <- expressom_pipeline_stages()
  expect_named(st, c("stage", "label", "target", "upstream"))
  expect_true(all(c("import", "eda", "dge", "dge_plots", "functional", "fgsea", "isoform", "report", "rdata") %in% st$stage))
  ups <- strsplit(st$upstream, " ")
  names(ups) <- st$stage
  expect_true(all(unlist(ups) %in% st$stage))
  # the isoform stage does not depend on the gene-level DESeq2 stage
  expect_false("dge" %in% ups$isoform)
  expect_true(all(c("eda", "dge_plots", "functional", "fgsea", "isoform") %in% ups$report))
})

test_that("write_expressom_targets stores the parameters and a one-line script", {
  skip_if_not_installed("targets")
  root <- withr::local_tempdir()
  p <- mock_params(root)
  paths <- write_expressom_targets(p)
  expect_true(file.exists(paths$params_file))
  expect_true(file.exists(paths$script))
  expect_identical(readRDS(paths$params_file), p)
  script <- readLines(paths$script)
  expect_true(any(grepl("ExpressOM::expressom_pipeline", script, fixed = TRUE)))
  expect_match(paths$dir, "\\.expressom_targets$")

  # identical parameters leave the file untouched (no spurious invalidation)
  m1 <- file.mtime(paths$params_file)
  Sys.sleep(1.1)
  write_expressom_targets(p)
  expect_identical(file.mtime(paths$params_file), m1)
  write_expressom_targets(utils::modifyList(p, list(padj_cutoff = 0.2)))
  expect_gt(as.numeric(file.mtime(paths$params_file)), as.numeric(m1))

  expect_error(write_expressom_targets(list(a = 1)), "out_dir")
})

test_that("the pipeline contains the stages that exist for the parameters", {
  skip_if_not_installed("targets")
  target_names <- function(...) {
    root <- withr::local_tempdir(.local_envir = parent.frame())
    p <- mock_params(root, ...)
    pf <- write_expressom_targets(p)$params_file
    vapply(expressom_pipeline(pf), function(t) t$settings$name, character(1))
  }

  full <- target_names()
  expect_true(all(c("params", "setup", "import", "eda", "dge", "dge_plots", "functional", "fgsea", "report", "rdata") %in% full))
  expect_false("isoform" %in% full)

  iso <- target_names(run_isoform = TRUE)
  expect_true(all(c("isoform", "m_isoform", "fp_isoform", "ctx_isoform") %in% iso))

  no_report <- target_names(make_report = FALSE)
  expect_false("report" %in% no_report)

  eda <- target_names(eda_only = TRUE)
  expect_true(all(c("import", "eda") %in% eda))
  expect_false(any(c("dge", "dge_plots", "functional", "fgsea", "rdata") %in% eda))

  no_dge <- target_names(run_dge = FALSE)
  expect_false("dge" %in% no_dge)
})

test_that("stage status before any run says 'never run' and follows the parameters", {
  skip_if_not_installed("targets")
  root <- withr::local_tempdir()
  p <- mock_params(root, run_isoform = TRUE)
  st <- expressom_stage_status(p, dir = expressom_targets_dir(p$out_dir))
  expect_true(all(st$state == "never run"))
  expect_true("isoform" %in% st$stage)
  expect_false("isoform" %in% expressom_stage_status(
    utils::modifyList(p, list(run_isoform = FALSE)), dir = expressom_targets_dir(p$out_dir)
  )$stage)
})

test_that("unknown stages are rejected", {
  skip_if_not_installed("targets")
  expect_error(.stage_targets("nope"), "Unknown stage")
  expect_error(expressom_targets_invalidate(withr::local_tempdir(), "nope"), "Unknown stage")
})
