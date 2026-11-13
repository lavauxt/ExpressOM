test_that("new ExpressOM logs use lowercase log/ and legacy Log/ remains readable", {
  out <- withr::local_tempdir()
  current <- .expressom_log_dir(out, create = TRUE)
  expect_identical(basename(current), "log")
  expect_true(dir.exists(current))

  legacy_out <- withr::local_tempdir()
  dir.create(file.path(legacy_out, "Log", "DGE"), recursive = TRUE)
  writeLines("{}", file.path(legacy_out, "Log", "DGE", "DGE_params_A_vs_B.json"))
  expect_identical(.expressom_log_dir(legacy_out), file.path(legacy_out, "Log"))
  expect_true(file.exists(.expressom_log_file(
    legacy_out, "DGE", "DGE_params_A_vs_B.json"
  )))
})

test_that("report log loader discovers OmniPathR and targets logs", {
  out <- withr::local_tempdir()
  dir.create(file.path(out, "log", "omnipathr"), recursive = TRUE)
  dir.create(file.path(out, "log", "targets"), recursive = TRUE)
  writeLines("OmniPath request log", file.path(out, "log", "omnipathr", "omnipathr.log"))
  writeLines("targets log", file.path(out, "log", "targets", "last_run.log"))
  logs <- .rep_load_logs(out, dge_log = NULL)
  expect_length(logs$omnipath_logs, 1L)
  expect_length(logs$targets_logs, 1L)
})

test_that("OmniPath cache recovery is opt-in and checks the configured cache", {
  skip_if_not_installed("OmnipathR")
  out <- withr::local_tempdir()
  cache_dir <- file.path(out, "cache")
  dir.create(cache_dir)
  writeLines("{}", file.path(cache_dir, "cache.json"))
  lock <- file.path(cache_dir, "cache.lock")
  writeLines("stale lock", lock)
  Sys.setFileTime(lock, Sys.time() - 3600)

  # Even a clearly old lock is left alone unless the operator explicitly opts
  # in after stopping other workers that could use this shared cache.
  withr::local_options(list(
    ExpressOM.omnipath_unlock_cache = FALSE,
    ExpressOM.omnipath_cache_recovery_checked = FALSE,
    omnipathr.cachedir = cache_dir
  ))
  expect_false(.decouple_unlock_stale_cache_db(cache_dir = cache_dir,
                                                stale_after = 600))
  expect_true(file.exists(lock))

  withr::local_options(list(
    ExpressOM.omnipath_unlock_cache = TRUE,
    ExpressOM.omnipath_cache_recovery_checked = FALSE,
    omnipathr.cachedir = cache_dir
  ))
  expect_true(.decouple_unlock_stale_cache_db(cache_dir = cache_dir,
                                               stale_after = 600))
  expect_false(file.exists(lock))
  expect_identical(getOption("omnipathr.cachedir"), cache_dir)
})

test_that("OmniPathR uses current option prefixes and writes its log inside log/", {
  skip_if_not_installed("OmnipathR")
  out <- withr::local_tempdir()
  withr::local_options(list(
    omnipathr.logdir = NULL,
    omnipathr.logfile = NULL,
    omnipathr.retry_downloads = NULL,
    omnipathr.retry_downloads_in_seconds = NULL,
    omnipathr.loglevel = NULL,
    omnipathr.console_loglevel = NULL,
    omnipathr.http_timeout = NULL,
    omnipathr.connect_timeout = NULL,
    omnipath.max_retries = NULL,
    ExpressOM.omnipathr_log_identity = NULL,
    ExpressOM.omnipathr_logfile = NULL,
    ExpressOM.omnipathr_log_run_id = NULL,
    ExpressOM.omnipathr_configured = NULL,
    # Avoid cache mutation in this option-only test.
    ExpressOM.omnipath_cache_recovery_checked = TRUE
  ))
  .configure_decouple_omnipathr(out_dir = out, comp_name = "case_vs_ctrl")
  expect_true(dir.exists(file.path(out, "log", "omnipathr")))
  expect_match(getOption("omnipathr.logfile"),
               "[/\\\\]log[/\\\\]omnipathr[/\\\\]")
  expect_equal(getOption("omnipathr.retry_downloads"), 3L)
  expect_equal(getOption("omnipathr.retry_downloads_in_seconds"), 2L)
  expect_null(getOption("omnipath.max_retries"))
})
