# These tests run the real pipeline definition (expressom_pipeline()) with stage
# stand-ins, in-process, and look at which stages re-ran.

skip_if_not_installed("targets")
skip_if_not_installed("withr")
skip_on_cran()

test_that("the first run builds every stage once, in dependency order", {
  root <- withr::local_tempdir()
  p <- mock_params(root, run_isoform = TRUE)
  r <- run_mock(p)

  expect_setequal(r$ran, c("setup", "import", "eda", "dge", "dge_plots", "functional", "fgsea", "isoform", "report", "rdata"))
  expect_true(all(table(read_calls(p$out_dir)) == 1L))
  expect_lt(match("import", r$ran), match("eda", r$ran))
  expect_lt(match("dge", r$ran), match("dge_plots", r$ran))
  expect_lt(match("dge", r$ran), match("functional", r$ran))
  expect_lt(match("dge", r$ran), match("fgsea", r$ran))
  expect_lt(match("eda", r$ran), match("report", r$ran))
  expect_lt(match("isoform", r$ran), match("report", r$ran))
  expect_true(all(r$status$state == "current"))
  expect_true(file.exists(file.path(p$out_dir, "report.html")))
})

test_that("re-running with identical parameters does nothing", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  run_mock(p)
  expect_identical(run_mock(p)$ran, character())
})

test_that("a cut-off change re-runs only the stages that use it", {
  root <- withr::local_tempdir()
  p <- mock_params(root, run_isoform = TRUE)
  run_mock(p)
  r <- run_mock(utils::modifyList(p, list(padj_cutoff = 0.2)))
  expect_false(any(c("import", "eda") %in% r$ran))
  expect_true(all(c("dge", "dge_plots", "functional", "fgsea", "isoform", "report") %in% r$ran))
})

test_that("a figure-only parameter re-runs the figures and the report, nothing else", {
  root <- withr::local_tempdir()
  p <- mock_params(root, run_isoform = TRUE)
  run_mock(p)
  r <- run_mock(utils::modifyList(p, list(volcano_colors = c(up = "orange", down = "navy", ns = "grey70"))))
  expect_setequal(setdiff(r$ran, "setup"), c("dge_plots", "report"))
})

test_that("a report setting re-runs only the report", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  run_mock(p)
  r <- run_mock(utils::modifyList(p, list(report_settings = list(title = "Other title"))))
  expect_setequal(setdiff(r$ran, "setup"), "report")
})

test_that("changing the contrast never re-imports the counts", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  run_mock(p)
  r <- run_mock(utils::modifyList(p, list(level = "drugB")))
  expect_false("import" %in% r$ran)
  expect_true(all(c("eda", "dge", "dge_plots", "functional", "fgsea", "report") %in% r$ran))
})

test_that("an isoform-only parameter leaves the gene-level stages alone", {
  root <- withr::local_tempdir()
  p <- mock_params(root, run_isoform = TRUE)
  run_mock(p)
  r <- run_mock(utils::modifyList(p, list(isoform_plot_top_n = 3)))
  expect_setequal(setdiff(r$ran, "setup"), c("isoform", "rdata", "report"))
})

test_that("editing an input file re-imports; touching it without changing it does not", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  run_mock(p)

  Sys.setFileTime(p$matrix_file, Sys.time() + 60)
  expect_false("import" %in% run_mock(p)$ran)

  writeLines(c("gene,S1,S2,S3,S4,S5,S6", "g1,9,9,9,9,9,9"), p$matrix_file)
  r <- run_mock(p)
  expect_true("import" %in% r$ran)
  expect_true(all(c("eda", "dge") %in% r$ran))
})

test_that("`stages` builds only what is asked for, and what it needs", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  r <- run_mock(p, stages = "dge")
  expect_setequal(setdiff(r$ran, "setup"), c("import", "dge"))
  expect_equal(r$status$state[r$status$stage == "report"], "never run")
  expect_equal(r$status$state[r$status$stage == "dge"], "current")

  r2 <- run_mock(p)
  expect_setequal(setdiff(r2$ran, "setup"), c("eda", "dge_plots", "functional", "fgsea", "report", "rdata"))
})

test_that("stale stages are reported before they are re-run, and propagate downstream", {
  root <- withr::local_tempdir()
  p <- mock_params(root, run_isoform = TRUE)
  run_mock(p)
  dir <- expressom_targets_dir(p$out_dir)

  st <- expressom_stage_status(utils::modifyList(p, list(volcano_colors = c(up = "orange", down = "navy", ns = "grey"))), dir = dir)
  stale <- st$stage[st$state == "stale"]
  expect_setequal(stale, c("dge_plots", "report"))

  st2 <- expressom_stage_status(utils::modifyList(p, list(padj_cutoff = 0.3)), dir = dir)
  expect_setequal(
    st2$stage[st2$state == "stale"],
    c("dge", "dge_plots", "functional", "fgsea", "isoform", "report", "rdata")
  )
  expect_equal(st2$state[st2$stage == "import"], "current")

  # saving the same parameters again does not make anything stale
  expect_true(all(expressom_stage_status(p, dir = dir)$state == "current"))
})

test_that("forcing a stage rebuilds it and what depends on it", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  run_mock(p)
  dir <- expressom_targets_dir(p$out_dir)
  expressom_targets_invalidate(dir, "dge")
  st <- expressom_stage_status(p, dir = dir)
  expect_true(st$state[st$stage == "dge"] %in% c("never run", "stale"))
  r <- run_mock(p)
  expect_true(all(c("dge", "report") %in% r$ran))
  expect_false(any(c("import", "eda") %in% r$ran))
})

test_that("a failing stage stops the run, is reported as errored and can be resumed", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  dir.create(p$out_dir, recursive = TRUE)
  file.create(file.path(p$out_dir, "FAIL_dge"))

  expect_error(run_mock(p))
  st <- expressom_stage_status(p, dir = expressom_targets_dir(p$out_dir))
  expect_equal(st$state[st$stage == "dge"], "errored")
  expect_false(st$state[st$stage == "report"] %in% c("current", "stale"))
  calls_before <- read_calls(p$out_dir)
  expect_equal(sum(calls_before == "import"), 1L)

  unlink(file.path(p$out_dir, "FAIL_dge"))
  r <- run_mock(p)
  expect_true(all(c("dge", "dge_plots", "report") %in% r$ran))
  # the import that had succeeded before the failure is not repeated
  expect_equal(sum(read_calls(p$out_dir) == "import"), 1L)
  expect_true(all(r$status$state == "current"))
})

test_that("a failing report does not prevent the other leaves, and is reported", {
  root <- withr::local_tempdir()
  p <- mock_params(root)
  dir.create(p$out_dir, recursive = TRUE)
  file.create(file.path(p$out_dir, "FAIL_report"))

  expect_warning(r <- run_mock(p), "Interactive report: report boom")
  expect_true("rdata" %in% r$ran)
  expect_equal(r$status$state[r$status$stage == "report"], "errored")
  expect_equal(r$status$state[r$status$stage == "rdata"], "current")

  # once the cause is fixed only the report is rebuilt
  unlink(file.path(p$out_dir, "FAIL_report"))
  r2 <- run_mock(p)
  expect_setequal(setdiff(r2$ran, "setup"), "report")
  expect_true(all(r2$status$state == "current"))
})
