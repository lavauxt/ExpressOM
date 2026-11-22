skip_if_not_installed("htmltools")

test_that("the tab tree has the expected top-level tabs for a full run", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  tree <- expressom_report_tabs(load_expressom_results(out))
  expect_identical(
    tree_titles(tree),
    c("Overview", "Samples & QC", "Differential expression", "Functional analysis",
      "Isoforms", "Figures", "Run info")
  )
  run_info <- tree[[which(tree_titles(tree) == "Run info")]]
  expect_identical(
    vapply(run_info$children, `[[`, character(1), "title"),
    c("Parameters", "Logs", "Files")
  )
  qc_titles <- vapply(tree[[which(tree_titles(tree) == "Samples & QC")]]$children,
                      `[[`, character(1), "title")
  expect_true(all(c("Sample similarity", "Expression distributions") %in% qc_titles))
  de_titles <- vapply(tree[[which(tree_titles(tree) == "Differential expression")]]$children,
                      `[[`, character(1), "title")
  expect_true("Diagnostics" %in% de_titles)
})

test_that("tabs for missing stages are left out", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out, functional = FALSE, isoform = FALSE, figures = FALSE)
  titles <- tree_titles(expressom_report_tabs(load_expressom_results(out)))
  expect_false(any(c("Functional analysis", "Isoforms", "Figures") %in% titles))
  expect_true(all(c("Overview", "Differential expression", "Run info") %in% titles))
})

test_that("every leaf of every run variant builds without error", {
  skip_if_no_report_stack()
  variants <- list(
    wald = list(),
    lrt = list(test = "LRT"),
    no_isoform = list(isoform = FALSE),
    no_functional = list(functional = FALSE),
    no_snapshot = list(snapshot = FALSE)
  )
  for (nm in names(variants)) {
    out <- withr::local_tempdir()
    do.call(make_run_dir, c(list(dir = out), variants[[nm]]))
    tree <- expressom_report_tabs(load_expressom_results(out))
    leaves <- tree_leaves(tree)
    expect_gt(length(leaves), 8L)
    for (lf in leaves) {
      html <- tryCatch(
        htmltools::renderTags(lf$fn())$html,
        error = function(e) paste("ERROR:", conditionMessage(e))
      )
      expect_false(
        grepl("^ERROR:", paste(html, collapse = "")),
        info = paste(nm, "/", paste(lf$path, collapse = " > "), ":", paste(html, collapse = ""))
      )
    }
  }
})

test_that("leaves are lazy: building the tree does not build any figure", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  tree <- expressom_report_tabs(load_expressom_results(out))
  leaves <- tree_leaves(tree)
  expect_true(all(vapply(leaves, function(l) is.function(l$fn), logical(1))))
})

test_that("an LRT overview reports significant genes without direction", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out, test = "LRT")
  tree <- expressom_report_tabs(load_expressom_results(out))
  html <- paste(htmltools::renderTags(tree[[1L]]$content())$html, collapse = "")
  expect_match(html, "LRT")
  expect_false(grepl("UPREGULATED", html))
})

test_that("display cut-offs override the run's own values", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  res <- load_expressom_results(out)
  strict <- expressom_report_settings(padj_cutoff = 0.0001, lfc_cutoff = 3)
  cut <- .xp_cutoffs(res, strict)
  expect_equal(unname(unlist(cut[c("padj", "lfc")])), c(0.0001, 3))
  def <- .xp_cutoffs(res, expressom_report_settings())
  expect_equal(unname(unlist(def[c("padj", "lfc")])), c(0.05, 1))
})

test_that("report settings have sensible defaults and are overridable", {
  s <- expressom_report_settings()
  expect_identical(s$title, "ExpressOM report")
  expect_true(s$static_figures)
  s2 <- expressom_report_settings(title = "x", top_terms = 5L, highlight_genes = "TP53")
  expect_identical(s2$title, "x")
  expect_identical(s2$top_terms, 5L)
  expect_identical(s2$highlight_genes, "TP53")
})

test_that("the tab emitter produces one heading per node and nests tabsets", {
  skip_if_no_report_stack()
  tree <- list(
    .xp_node("A", .xp_node("A1", content = function() htmltools::p("one")), .xp_node("A2", content = function() htmltools::p("two"))),
    .xp_node("B", content = function() htmltools::p("three"))
  )
  # the emitter runs knit_child(), so it has to be called while knitting
  withr::local_dir(withr::local_tempdir())
  env <- list2env(list(tree = tree), parent = environment(.xp_emit_tabs))
  md <- knitr::knit(
    text = c("```{r, results='asis', echo=FALSE}", "cat(.xp_emit_tabs(tree, level = 2L), sep = '\\n')", "```"),
    envir = env, quiet = TRUE
  )
  md <- paste(md, collapse = "\n")
  expect_match(md, "## A \\{\\.tabset \\.tabset-pills\\}")
  expect_match(md, "### A1")
  expect_match(md, "### A2")
  expect_match(md, "## B")
  expect_match(md, "<p>one</p>", fixed = TRUE)
})

test_that("the theme setting is validated and defaults to the viewer's system", {
  expect_identical(expressom_report_settings()$theme, "auto")
  expect_identical(expressom_report_settings(theme = "dark")$theme, "dark")
  expect_error(expressom_report_settings(theme = "sepia"), "should be one of")
  # a hand-written or old settings list is tolerated
  expect_identical(.xp_theme_setting(list(theme = "light")), "light")
  expect_identical(.xp_theme_setting(list(theme = "sepia")), "auto")
  expect_identical(.xp_theme_setting(list()), "auto")
})

test_that("the header carries the theme switch, the dark stylesheet and the chosen default", {
  skip_if_not_installed("htmltools")
  res <- structure(list(comparison = NA_character_, config = list()), class = "expressom_results")
  dark <- as.character(.xp_header_html(res, expressom_report_settings(theme = "dark")))
  auto <- as.character(.xp_header_html(res, expressom_report_settings()))
  expect_match(dark, 'id="xp-theme-toggle"', fixed = TRUE)
  expect_match(dark, "html[data-theme='dark']", fixed = TRUE)
  expect_match(dark, "var DEFAULT = 'dark';", fixed = TRUE)
  expect_match(auto, "var DEFAULT = 'auto';", fixed = TRUE)
  expect_false(grepl("__XP_DEFAULT__", dark, fixed = TRUE))
  # every colour of the page goes through the variables the dark theme redefines
  css <- .xp_css()
  expect_match(css, "--xp-bg:#0f1720", fixed = TRUE)
  expect_match(css, "--xp-bg:#ffffff", fixed = TRUE)
  expect_false(grepl("color:\\s*#1f2933", css))
})

test_that("a failing tab shows its error and does not stop the other tabs", {
  skip_if_no_report_stack()
  tree <- list(
    .xp_node("Broken", content = function() stop("boom")),
    .xp_node("Fine", content = function() htmltools::p("still here"))
  )
  withr::local_dir(withr::local_tempdir())
  env <- list2env(list(tree = tree), parent = environment(.xp_emit_tabs))
  md <- knitr::knit(
    text = c("```{r, results='asis', echo=FALSE}", "cat(.xp_emit_tabs(tree, level = 2L), sep = '\\n')", "```"),
    envir = env, quiet = TRUE
  )
  md <- paste(md, collapse = "\n")
  expect_match(md, "This tab could not be built", fixed = TRUE)
  expect_match(md, "boom", fixed = TRUE)
  expect_match(md, "still here", fixed = TRUE)
})

test_that(".xp_node drops empty nodes and keeps leaves", {
  expect_null(.xp_node("empty"))
  expect_true(.xp_is_leaf(.xp_node("x", content = function() NULL)))
  n <- .xp_node("parent", NULL, .xp_node("kid", content = function() NULL))
  expect_length(n$children, 1L)
  expect_false(.xp_is_leaf(n))
})

# ---- decoupleR subtree --------------------------------------------------------

test_that("the decoupleR subtree has one tab per resource plus network status", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  tree <- expressom_report_tabs(load_expressom_results(out))
  func <- tree[[which(tree_titles(tree) == "Functional analysis")]]
  titles <- vapply(func$children, `[[`, character(1), "title")
  expect_true("decoupleR (TF / pathway activity)" %in% titles)
  dec <- func$children[[which(titles == "decoupleR (TF / pathway activity)")]]
  dec_titles <- vapply(dec$children, `[[`, character(1), "title")
  expect_setequal(dec_titles, c("COLLECTRI", "PROGENY", "DOROTHEA", "Network status"))
  # the network-status table is not listed as a resource tab of its own
  expect_false("DECOUPLER" %in% dec_titles)
})

test_that("every decoupleR leaf builds (including the network-status table)", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  res <- load_expressom_results(out)
  st <- expressom_report_settings()
  html <- paste(htmltools::renderTags(.xp_leaf_decouple_status(res, st)())$html, collapse = "")
  expect_match(html, "offline (cache / bundled fallback)", fixed = TRUE)
  expect_match(html, "decoupleR_network_status", fixed = TRUE)
  for (entry in res$decouple) {
    built <- .xp_leaf_decouple(res, st, entry)()
    expect_false(inherits(built, "error"))
    html <- paste(htmltools::renderTags(built)$html, collapse = "")
    expect_match(html, "How to read this plot", fixed = TRUE)
  }
})

test_that("a decoupleR table without score / p-value columns builds gracefully", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out, decouple = FALSE)
  res <- load_expressom_results(out)
  st <- expressom_report_settings()
  entry <- list(
    resource = "DecoupleR",
    path = "DecoupleR/decoupleR_network_status_x_vs_y.csv",
    table = data.frame(resource = "collectri", status = "network (live)",
                       n_edges = 10L, stringsAsFactors = FALSE)
  )
  html <- paste(htmltools::renderTags(.xp_leaf_decouple(res, st, entry)())$html, collapse = "")
  expect_false(grepl("could not be built", html))
  expect_match(html, "source(s) tested", fixed = TRUE)
})

# ---- merged enrichment tables -------------------------------------------------

test_that(".rep_rbind_fill unions columns of mismatched tables", {
  a <- data.frame(term = "x", padj = 0.1, extra_a = 1)
  b <- data.frame(term = "y", padj = 0.2, extra_b = "z", stringsAsFactors = FALSE)
  out <- .rep_rbind_fill(a, b, NULL, data.frame())
  expect_equal(nrow(out), 2L)
  expect_setequal(colnames(out), c("term", "padj", "extra_a", "extra_b"))
  expect_true(is.na(out$extra_a[[2L]]))
  expect_true(is.na(out$extra_b[[1L]]))
  expect_null(.rep_rbind_fill())
  # factors never require a shared level set
  f <- .rep_rbind_fill(
    data.frame(t = "a", g = factor("A")),
    data.frame(t = "b", g = factor("B"))
  )
  expect_identical(f$g, c("A", "B"))
})

test_that("the merged-collection leaf handles tables with different columns", {
  skip_if_no_report_stack()
  out <- withr::local_tempdir()
  make_run_dir(out)
  res <- load_expressom_results(out)
  st <- expressom_report_settings()
  inv <- res$functional$inventory
  ora <- inv[inv$source == "DGE" & grepl("ORA$", inv$method), , drop = FALSE]
  expect_equal(nrow(ora), 2L)
  # give the two ORA collections different columns to mimic mixed sources
  first <- res$functional$tables[[ora$path[[1L]]]]
  first$only_in_first <- seq_len(nrow(first))
  res$functional$tables[[ora$path[[1L]]]] <- first
  leaf <- .xp_leaf_enrichment_merged(res, st, ora$path, ora$collection, kind = "ora")
  html <- paste(htmltools::renderTags(leaf())$html, collapse = "")
  expect_false(grepl("could not be built", html))
  expect_match(html, "rows across 2 collections", fixed = TRUE)
})

# ---- per-tab help -------------------------------------------------------------

test_that("every known tab title has a one-line help and unknown titles get none", {
  known <- c("Samples & QC", "Differential expression", "Functional analysis",
             "Isoforms", "Figures", "Run info", "Volcano", "MA plot",
             "Top genes heatmap", "Parameters", "Logs", "decoupleR (TF / pathway activity)")
  for (t in known) expect_type(.xp_tab_help(t), "character")
  expect_null(.xp_tab_help("No such tab"))
  expect_null(.xp_tab_help(NA_character_))
  # helper columns / identities are dropped when slim_tables is on
  s <- expressom_report_settings()
  expect_identical(s$table_rows, 5000L)
  expect_true(s$slim_tables)
})

test_that("the htmltools renderer prepends the tab help", {
  skip_if_not_installed("htmltools")
  tree <- list(
    .xp_node("Differential expression",
             .xp_node("Volcano", content = function() htmltools::p("plot"))),
    .xp_node("No such tab", content = function() htmltools::p("quiet"))
  )
  out <- .xp_render_tabset(tree)
  html <- paste(htmltools::renderTags(htmltools::tagList(out$nav, out$content))$html, collapse = "")
  expect_match(html, "Gene-level differential expression", fixed = TRUE)
  expect_match(html, "points beyond the dashed lines", fixed = TRUE)
  expect_false(grepl("xp-note\">quiet", html))
})
