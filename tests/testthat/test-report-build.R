test_that("build_expressom_report writes one self-contained interactive HTML file", {
  skip_on_cran()
  skip_if_no_report_stack()
  skip_if_no_pandoc()
  out <- withr::local_tempdir()
  make_run_dir(out)

  path <- suppressMessages(build_expressom_report(out))
  expect_true(file.exists(path))
  expect_identical(basename(path), "ExpressOM_report_treated_vs_control.html")
  expect_identical(basename(dirname(path)), "Report")

  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  # one file: scripts / styles are inlined, nothing is fetched from a CDN
  expect_false(grepl("<script[^>]+src=\"http", html))
  expect_false(grepl("<link[^>]+href=\"http", html))
  # tabs, plotly figures and DT tables are in the page
  for (tab in c("Overview", "Differential expression", "Functional analysis", "Isoforms", "Run info")) {
    expect_match(html, tab, fixed = TRUE)
  }
  expect_match(html, "plotly", fixed = TRUE)
  expect_match(html, "datatables", fixed = TRUE)
  expect_match(html, "tabset", fixed = TRUE)
  expect_gt(length(gregexpr("html-widget", html, fixed = TRUE)[[1L]]), 5L)
})

test_that("the report has a light/dark switch and honours the theme setting", {
  skip_on_cran()
  skip_if_no_report_stack()
  skip_if_no_pandoc()
  out <- withr::local_tempdir()
  make_run_dir(out, functional = FALSE, isoform = FALSE, figures = FALSE)
  res <- load_expressom_results(out)

  read_html <- function(theme) {
    path <- suppressMessages(build_expressom_report(
      out, results = res, output_file = file.path(withr::local_tempdir(), "r.html"),
      settings = expressom_report_settings(theme = theme)
    ))
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }
  auto <- read_html("auto")
  expect_match(auto, 'id="xp-theme-toggle"', fixed = TRUE)
  expect_match(auto, "html[data-theme='dark']", fixed = TRUE)
  expect_match(auto, "prefers-color-scheme", fixed = TRUE)
  expect_match(auto, "var DEFAULT = 'auto';", fixed = TRUE)
  expect_false(grepl("__XP_DEFAULT__", auto, fixed = TRUE))
  expect_match(read_html("dark"), "var DEFAULT = 'dark';", fixed = TRUE)
})

test_that("the report honours output_file, comparison, title and a pre-loaded result", {
  skip_on_cran()
  skip_if_no_report_stack()
  skip_if_no_pandoc()
  out <- withr::local_tempdir()
  make_run_dir(out, functional = FALSE, isoform = FALSE, figures = FALSE)
  make_run_dir(out, level = "drugB", functional = FALSE, isoform = FALSE, figures = FALSE)
  target <- file.path(withr::local_tempdir(), "custom", "r.html")

  res <- load_expressom_results(out, comparison = "drugB_vs_control")
  path <- suppressMessages(build_expressom_report(
    out, output_file = target, results = res,
    settings = expressom_report_settings(title = "My study", max_points = 500)
  ))
  expect_identical(normalizePath(path), normalizePath(target))
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "My study", fixed = TRUE)
  expect_match(html, "drugB", fixed = TRUE)
})

test_that("an EDA-only directory still produces a report", {
  skip_on_cran()
  skip_if_no_report_stack()
  skip_if_no_pandoc()
  out <- withr::local_tempdir()
  make_run_dir(out, functional = FALSE, isoform = FALSE, figures = FALSE)
  unlink(file.path(out, c("DE_raw_results", "log")), recursive = TRUE)
  file.rename(
    file.path(out, "Plots", "PCA_treated_vs_control_sample_scores.tsv"),
    file.path(out, "Plots", "PCA_EDA_sample_scores.tsv")
  )
  path <- suppressMessages(build_expressom_report(out))
  expect_identical(basename(path), "ExpressOM_report.html")
})

test_that("build_expressom_report validates its inputs", {
  skip_if_no_report_stack()
  expect_error(build_expressom_report(file.path(tempdir(), "nope")), "existing ExpressOM output directory")
  out <- withr::local_tempdir()
  make_run_dir(out)
  expect_error(build_expressom_report(out, results = list()), "load_expressom_results")
})

test_that("a missing report package gives a clear message", {
  expect_error(.rep_require("definitelyNotAPackage123"), "definitelyNotAPackage123")
})
