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

# ---- OmniPath reachability probe ---------------------------------------------

test_that("the OmniPath probe is a plain HTTP check and reports the real verdict", {
  skip_on_cran()
  skip_if_not_installed("OmnipathR")
  skip_on_os("windows")  # file:// download method differs on Windows
  out <- withr::local_tempdir()

  live_file <- file.path(withr::local_tempdir(), "queries.json")
  writeLines('{"query_type":"annotations;complexes;enzsub;interactions;intercell"}',
             live_file)
  ok <- .decouple_omnipath_probe(
    api_url = paste0("file://", live_file),
    out_dir = out, comp_name = "case_vs_ctrl"
  )
  expect_true(ok$available)
  expect_null(ok$reason)

  fail <- .decouple_omnipath_probe(
    api_url = "http://127.0.0.1:9/queries",  # port 9: connection refused
    timeout = 2,
    out_dir = out, comp_name = "case_vs_ctrl"
  )
  expect_false(fail$available)
  expect_match(fail$reason, "request failed")
  # regression guard: a code-level failure must never be reported as offline
  expect_false(grepl("not an exported object", fail$reason))
})

test_that("a broken HTTP configuration keeps the probe at NA instead of lying", {
  skip_on_cran()
  skip_if_not_installed("OmnipathR")
  out <- withr::local_tempdir()
  withr::local_options(list(download.file.method = "no-such-download-method"))
  probe <- .decouple_omnipath_probe(
    api_url = "https://omnipathdb.org/queries",
    out_dir = out, comp_name = "case_vs_ctrl"
  )
  expect_true(is.na(probe$available))
  expect_match(probe$reason, "method")
})

test_that("OmnipathR gets a log destination even without an out_dir", {
  withr::local_options(list(
    ExpressOM.omnipathr_log_identity = NULL,
    ExpressOM.omnipathr_logfile = NULL,
    ExpressOM.omnipathr_log_run_id = NULL
  ))
  # no out_dir / log_dir: the fallback must be under tempdir(), never the cwd
  .configure_decouple_omnipathr()
  lf <- getOption("omnipathr.logfile")
  expect_true(grepl(tempdir(), lf, fixed = TRUE))
  expect_false(grepl(getwd(), normalizePath(lf, winslash = "/"), fixed = TRUE))
  # a second call in the same session reuses the same file
  .configure_decouple_omnipathr()
  expect_identical(getOption("omnipathr.logfile"), lf)
})

test_that("stray omnipathr-log folders are swept into the results", {
  cwd <- withr::local_tempdir()
  out <- withr::local_tempdir()
  # folder form + flat form
  dir.create(file.path(cwd, "omnipathr-log"))
  writeLines("log1", file.path(cwd, "omnipathr-log", "omnipathr-log-run-x.log"))
  writeLines("log2", file.path(cwd, "omnipathr-log-flat.log"))
  moved <- .decouple_sweep_omnipathr_logs(out, cwd = cwd)
  expect_gte(moved, 2L)
  dest <- file.path(out, "log", "omnipathr")
  expect_true(file.exists(file.path(dest, "omnipathr-log-run-x.log")))
  expect_true(file.exists(file.path(dest, "omnipathr-log-flat.log")))
  # the emptied folder was removed; a later clash is resolved by suffixing,
  # the original file's content survives
  dir.create(file.path(cwd, "omnipathr-log"))
  writeLines("clash", file.path(cwd, "omnipathr-log", "omnipathr-log-flat.log"))
  writeLines("keep", file.path(cwd, "omnipathr-log", "omnipathr-log-keep.log"))
  moved2 <- .decouple_sweep_omnipathr_logs(out, cwd = cwd)
  expect_gte(moved2, 2L)
  expect_identical(readLines(file.path(dest, "omnipathr-log-run-x.log")), "log1")
  expect_true(length(list.files(dest, pattern = "flat")) >= 2L)
})

test_that("the network-status label follows the probe verdict only", {
  expect_identical(.decouple_net_status_label(TRUE), "network (live)")
  expect_identical(.decouple_net_status_label(FALSE),
                   "offline (cache / bundled fallback)")
  expect_identical(.decouple_net_status_label(NA), "network (OmniPath)")
})
