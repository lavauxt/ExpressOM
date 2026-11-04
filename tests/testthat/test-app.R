skip_if_not_installed("shiny")

# ---- parameter <-> widget conversion -------------------------------------------------------

test_that("every widget value converts to a parameter and back", {
  spec <- .app_param_spec()
  d <- expressom_default_params()
  for (s in spec) {
    expect_true(s$name %in% names(d), info = s$name)
    widget <- .app_to_input(s, d[[s$name]])
    back <- .app_from_input(s, widget, d[[s$name]])
    expect_equal(back, d[[s$name]], info = s$name)
  }
})

test_that("widget values are cleaned on the way in", {
  txt <- .app_spec("x", "x", "text", "g")
  csv <- .app_spec("x", "x", "csv", "g")
  num <- .app_spec("x", "x", "number", "g")
  sel <- .app_spec("x", "x", "select", "g", c("a", "b"))
  expect_null(.app_from_input(txt, "   ", "default"))
  expect_identical(.app_from_input(txt, "  hello ", NULL), "hello")
  expect_identical(.app_from_input(csv, " TP53 , MYC;\nEGFR,, ", NULL), c("TP53", "MYC", "EGFR"))
  expect_null(.app_from_input(csv, " , ", "x"))
  expect_equal(.app_from_input(num, 0.05, NULL), 0.05)
  expect_null(.app_from_input(num, NA, 1))
  expect_identical(.app_from_input(sel, "b", "a"), "b")
  # a widget that does not exist yet keeps the parameter
  expect_identical(.app_from_input(txt, NULL, "kept"), "kept")
  # unknown select values fall back to the first choice
  expect_identical(.app_to_input(sel, "zzz"), "a")
})

test_that(".app_build_params applies widget values, report settings and JSON extras", {
  base <- expressom_params(out_dir = "x")
  vals <- list(
    padj_cutoff = 0.001, level = "drugB", highlight_genes = "TP53, MYC",
    r_title = "My study", r_max_points = 20000, r_top_terms = 20, r_heatmap_genes = 50,
    r_static_figures = TRUE
  )
  p <- .app_build_params(vals, base, '{"zscore_genes": ["A", "B"], "top_genes": 12}')
  expect_equal(p$padj_cutoff, 0.001)
  expect_identical(p$level, "drugB")
  expect_identical(p$highlight_genes, c("TP53", "MYC"))
  expect_identical(p$zscore_genes, c("A", "B"))
  expect_equal(p$top_genes, 12)
  # only settings that differ from the defaults are stored
  expect_identical(p$report_settings, list(title = "My study", max_points = 20000))

  expect_error(.app_build_params(vals, base, '{"nonsense": 1}'), "Unknown parameter")
  expect_error(.app_build_params(vals, base, "{not json"), "not valid JSON")
  expect_error(.app_build_params(vals, base, "[1, 2]"), "JSON object")
  # no JSON and default report widgets -> no report_settings
  p2 <- .app_build_params(list(r_title = "ExpressOM report"), base)
  expect_null(p2$report_settings)
})

test_that(".app_validate_params flags what would make the run fail", {
  root <- withr::local_tempdir()
  ok <- mock_params(root)
  expect_identical(.app_validate_params(ok), character(0))

  bad <- function(...) .app_validate_params(utils::modifyList(ok, list(...)))
  expect_match(bad(out_dir = ""), "Output directory")
  expect_match(bad(level = NULL), "Level")
  expect_match(bad(level = "x", base = "x"), "must differ")
  expect_match(bad(model = "~ ~ ~"), "cannot be parsed")
  expect_match(bad(test = "LRT", reduced = NULL), "reduced formula")
  expect_match(bad(padj_cutoff = 2), "cut-off")
  expect_match(bad(padj_cutoff = NULL), "cut-off")
  expect_match(bad(sample_table = file.path(root, "missing.csv")), "not found")
  expect_match(bad(count_type = "matrix", matrix_file = NULL), "matrix")
  # an EDA-only run needs no contrast
  expect_identical(bad(eda_only = TRUE, run_dge = FALSE, level = NULL, base = NULL, model = NULL), character(0))
})

test_that("initial parameters come from the defaults, a saved run, or both", {
  expect_identical(.app_initial_params()$level, "treated")

  out <- withr::local_tempdir()
  make_run_dir(out, level = "drugB")
  p <- .app_initial_params(out_dir = out)
  expect_identical(p$level, "drugB")
  expect_equal(normalizePath(p$out_dir), normalizePath(out))

  # explicit parameters win over the saved run
  p2 <- .app_initial_params(out_dir = out, params = list(level = "other", top_genes = 7))
  expect_identical(p2$level, "other")
  expect_equal(p2$top_genes, 7)
  expect_identical(p2$base, "control")

  # a folder with no saved run starts from the defaults, pointed at that folder
  empty <- withr::local_tempdir()
  expect_equal(normalizePath(.app_initial_params(out_dir = empty)$out_dir), normalizePath(empty))

  # the exact parameters of a targets run are preferred over the JSON snapshot
  root <- withr::local_tempdir()
  q <- mock_params(root, top_genes = 99)
  skip_if_not_installed("targets")
  write_expressom_targets(q)
  expect_equal(.app_initial_params(out_dir = q$out_dir)$top_genes, 99)
})

# ---- small renderers --------------------------------------------------------------------------

test_that("the stage table shows one pill per stage with its state", {
  st <- data.frame(
    stage = c("import", "dge"), label = c("Import counts", "DESeq2"),
    progress = c("completed", "none"), seconds = c(65, NA),
    finished = c("2026-01-01 10:00:00", NA), state = c("current", "never run"),
    stringsAsFactors = FALSE
  )
  html <- as.character(.app_stage_table(st))
  expect_match(html, "xp-current", fixed = TRUE)
  expect_match(html, "xp-never-run", fixed = TRUE)
  expect_match(html, "1.1 min", fixed = TRUE)
  expect_match(as.character(.app_stage_table(NULL)), "Nothing to show")
})

test_that("the log viewer tails the file and strips terminal colours", {
  f <- withr::local_tempfile()
  expect_match(.app_log_tail(f), "No run yet")
  writeLines(paste0("\033[32mline ", 1:600, "\033[0m"), f)
  out <- .app_log_tail(f, n = 5)
  expect_identical(strsplit(out, "\n")[[1L]], paste("line", 596:600))
})

# ---- the tab tree as shiny tabs ---------------------------------------------------------------------

test_that("the result tree becomes nested tabsets with one lazy output per leaf", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  tree <- expressom_report_tabs(load_expressom_results(out))
  tu <- .app_tree_ui(tree, prefix = "t1", ns = shiny::NS("res"))
  expect_equal(length(tu$leaves), length(tree_leaves(tree)))
  html <- as.character(tu$ui)
  expect_match(html, "Differential expression", fixed = TRUE)
  for (id in names(tu$leaves)) {
    expect_match(html, paste0('id="res-', id, '"'), fixed = TRUE)
  }
  # a leaf that fails shows an explanation instead of breaking the page
  safe <- as.character(.app_safe_leaf(function() stop("kaput")))
  expect_match(safe, "could not be built")
  expect_match(safe, "kaput")
})

# ---- modules (shiny::testServer) ---------------------------------------------------------------------

test_that("the parameter module turns widget values into a parameter list", {
  base <- expressom_params(out_dir = "x")
  shiny::testServer(.app_params_server, args = list(init = base), {
    session$setInputs(p_padj_cutoff = 0.05, p_level = "drugB", p_highlight_genes = "A, B", extra_json = "")
    p <- params()
    expect_equal(p$padj_cutoff, 0.05)
    expect_identical(p$level, "drugB")
    expect_identical(p$highlight_genes, c("A", "B"))

    session$setInputs(extra_json = '{"zzz": 1}')
    expect_s3_class(params(), "xp_param_error")
    expect_match(params()$message, "Unknown parameter")

    session$setInputs(extra_json = '{"top_genes": 5}')
    expect_equal(params()$top_genes, 5)
  })
})

test_that("the parameter module reloads a saved run into the widgets", {
  out <- withr::local_tempdir()
  make_run_dir(out, level = "drugB")
  shiny::testServer(.app_params_server, args = list(init = expressom_params(out_dir = out)), {
    session$setInputs(p_out_dir = out, load_run = 1)
    expect_identical(output$loaded, paste0("Parameters loaded from ", .abs_path(out)))
  })
})

test_that("the results module loads a run and builds one output per visible leaf", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  shiny::testServer(
    .app_results_server,
    args = list(out_dir = shiny::reactive(out), version = shiny::reactive(0L)),
    {
      session$setInputs(comparison = "treated_vs_control", padj = NA, lfc = NA, highlight = "", reload = 0)
      res <- results()
      expect_s3_class(res, "expressom_results")
      expect_identical(comparison(), "treated_vs_control")
      tu <- tree_ui()
      expect_gt(length(tu$leaves), 8L)
      expect_match(as.character(tu$ui), "Differential expression", fixed = TRUE)

      # display overrides rebuild the tree without touching the files
      session$setInputs(padj = 0.001, lfc = 2, highlight = "GENE1, GENE2")
      s <- settings()
      expect_equal(s$padj_cutoff, 0.001)
      expect_equal(s$lfc_cutoff, 2)
      expect_identical(s$highlight_genes, c("GENE1", "GENE2"))

      # a leaf renders through the module's own output
      lid <- names(tree_ui()$leaves)[[1L]]
      expect_false(is.null(output[[lid]]))
    }
  )
})

test_that("the results module waits politely when there is nothing to show", {
  skip_if_no_report_stack()
  empty <- withr::local_tempdir()
  shiny::testServer(
    .app_results_server,
    args = list(out_dir = shiny::reactive(file.path(empty, "nope")), version = shiny::reactive(0L)),
    {
      session$setInputs(comparison = "", padj = NA, lfc = NA, highlight = "", reload = 0)
      expect_error(results(), class = "shiny.silent.error")
    }
  )
})

test_that("the pipeline module launches a background run, follows it, and refreshes the status", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("targets")
  root <- withr::local_tempdir()
  p <- mock_params(root)
  marker <- file.path(root, "ran.txt")
  runner <- eval(substitute(
    function(params, dir, stages, force) {
      cat("hello from the runner\n")
      writeLines("done", MARKER)
      invisible(TRUE)
    },
    list(MARKER = marker)
  ))

  shiny::testServer(
    .app_pipeline_server,
    args = list(params = shiny::reactive(p), runner = runner, poll_ms = 100L),
    {
      session$flushReact()
      st <- status_val()
      expect_true(all(st$state == "never run"))
      expect_equal(version(), 0L)

      session$setInputs(run_all = 1)
      expect_false(is.null(proc()))
      for (i in 1:80) {
        Sys.sleep(0.25)
        session$elapse(300)
        if (is.null(proc())) break
      }
      expect_null(proc())
      expect_true(file.exists(marker))
      expect_equal(version(), 1L)
      expect_match(.app_log_tail(log_file()), "hello from the runner")
    }
  )
})

test_that("the pipeline module refuses to launch with invalid parameters", {
  skip_if_not_installed("callr")
  skip_if_not_installed("targets")
  root <- withr::local_tempdir()
  p <- mock_params(root, level = "same", base = "same")
  shiny::testServer(
    .app_pipeline_server,
    args = list(params = shiny::reactive(p), runner = function(...) stop("must not run"), poll_ms = 100L),
    {
      session$setInputs(run_all = 1)
      expect_null(proc())
    }
  )
})

test_that("the report module finds the newest report of a comparison", {
  out <- withr::local_tempdir()
  expect_null(.app_latest_report(out))
  dir.create(file.path(out, "Report"))
  a <- file.path(out, "Report", "ExpressOM_report_a_vs_b.html")
  b <- file.path(out, "Report", "ExpressOM_report_c_vs_d.html")
  writeLines("a", a)
  writeLines("b", b)
  Sys.setFileTime(a, Sys.time() - 100)
  expect_identical(basename(.app_latest_report(out)), basename(b))
  expect_identical(basename(.app_latest_report(out, "a_vs_b")), basename(a))
})

test_that("the app object assembles", {
  skip_if_not_installed("targets")
  skip_if_not_installed("callr")
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  app <- expressom_app(out_dir = out)
  expect_s3_class(app, "shiny.appobj")
  html <- as.character(app$httpHandler(list(PATH_INFO = "/", REQUEST_METHOD = "GET", rook.input = NULL, HTTP_ACCEPT = "text/html"))$content)
  expect_match(html, "ExpressOM", fixed = TRUE)
  expect_match(html, "Pipeline", fixed = TRUE)
})
