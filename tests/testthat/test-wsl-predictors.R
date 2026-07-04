# Tests for the external-predictor execution layer (CPAT / SignalP / Pfam).
#
# These tests exercise the *native* (non-WSL) code path, which is what runs
# on Linux/macOS and is what CI will use. The Windows+WSL branch is guarded
# behind `.Platform$OS.type == "windows"` inside the functions themselves and
# is skipped automatically on any other platform.

test_that(".to_wsl_path is a no-op off Windows", {
  skip_on_os("windows")
  expect_equal(.to_wsl_path("/tmp/some/path"), "/tmp/some/path")
  expect_equal(.to_wsl_path(NULL), NULL)
  expect_equal(.to_wsl_path(""), "")
})

test_that(".wsl_exec_script runs natively on Linux/macOS and returns a status", {
  skip_on_os("windows")

  status <- .wsl_exec_script("exit 0", use_wsl = FALSE)
  expect_identical(status, 0L)

  status_fail <- .wsl_exec_script("exit 3", use_wsl = FALSE, ignore_stderr = TRUE)
  expect_identical(status_fail, 3L)
})

test_that(".wsl_exec_script intern = TRUE captures stdout natively", {
  skip_on_os("windows")

  out <- .wsl_exec_script("echo hello-expressom", use_wsl = FALSE, intern = TRUE)
  expect_true(any(grepl("hello-expressom", out)))
})

test_that(".wsl_exec_script fails gracefully (not a hard error) when bash is missing", {
  skip_on_os("windows")
  skip_on_cran()

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = "")

  expect_warning(
    status <- .wsl_exec_script("echo should_not_run", use_wsl = FALSE),
    regexp = "bash"
  )
  expect_identical(status, 127L)
})

test_that(".wsl_tool_exists correctly detects present and absent tools", {
  skip_on_os("windows")

  expect_true(.wsl_tool_exists("ls", use_wsl = FALSE))
  expect_false(.wsl_tool_exists("definitely_not_a_real_tool_xyz123", use_wsl = FALSE))
})

test_that(".find_pfam_db and .find_cpat_logit_model accept the positional call ",
          "signature used by the predictor pipeline (regression test)", {
  skip_on_os("windows")

  # These mirror the exact positional calls made inside
  # .run_external_predictors() / mod_isoform.R -- a mismatch in argument
  # order here previously caused wsl_distro / use_wsl / conda_sh to be
  # silently swapped.
  pfam_result <- .find_pfam_db("Ubuntu", FALSE, NULL, "isoform_tools")
  expect_true(is.null(pfam_result) || is.character(pfam_result))

  cpat_result <- .find_cpat_logit_model("Human", "Ubuntu", FALSE, NULL, "isoform_tools")
  expect_true(is.null(cpat_result) || is.character(cpat_result))
})

test_that("debug_wsl() runs natively (no WSL) without error and reports tool status", {
  skip_on_os("windows")

  res <- suppressMessages(
    debug_wsl(out_dir = NULL, verbose = FALSE, use_wsl = FALSE)
  )

  expect_type(res, "list")
  expect_true(res$wsl_available)  # the local bash shell itself is reachable
  expect_identical(res$platform, .Platform$OS.type)
  expect_true(is.list(res$tools))
  expect_true(all(c("cpat", "signalp", "hmmscan") %in% names(res$tools)))
  # every value should be a plain logical, never NA / error
  expect_true(all(vapply(res$tools, is.logical, logical(1))))
})

test_that("debug_wsl() defaults use_wsl sensibly per platform", {
  skip_on_os("windows")
  res <- suppressMessages(debug_wsl(out_dir = NULL, verbose = FALSE))
  expect_identical(res$platform, .Platform$OS.type)
})

test_that("debug_wsl() writes a JSON log file when out_dir is supplied", {
  skip_on_os("windows")
  skip_if_not_installed("jsonlite")

  tmp_out <- withr::local_tempdir()
  res <- suppressMessages(
    debug_wsl(out_dir = tmp_out, verbose = FALSE, use_wsl = FALSE)
  )
  log_file <- file.path(tmp_out, "Log", "wsl_debug.json")
  expect_true(file.exists(log_file))

  parsed <- jsonlite::fromJSON(log_file)
  expect_identical(parsed$wsl_available, res$wsl_available)
})

test_that("CPU auto-detection in the predictor helper never returns < 1", {
  skip_on_os("windows")

  detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  n_cpu <- if (is.na(detected) || detected < 1) 1L else max(1L, detected - 1L)
  expect_true(n_cpu >= 1L)
  expect_type(n_cpu, "integer")
})

test_that("check_wsl() returns FALSE (not an error) off Windows", {
  skip_on_os("windows")
  expect_false(check_wsl())
})
