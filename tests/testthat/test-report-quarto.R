# Tests for the Quarto dashboard engine (engine = "quarto").

skip_if_no_quarto <- function() {
  bin <- Sys.which("quarto")
  if (nzchar(Sys.getenv("QUARTO_PATH"))) bin <- Sys.getenv("QUARTO_PATH")
  testthat::skip_if(!nzchar(bin), "the quarto CLI is not installed")
}

# decode a container written by .xq_encode_table() back into R vectors
xq_decode <- function(enc) {
  raw <- memDecompress(jsonlite::base64_dec(enc$b64), type = "gzip")
  hl <- readBin(raw[1:4], "integer", size = 4L, endian = "little")
  hdr <- jsonlite::fromJSON(rawToChar(raw[5:(4 + hl)]), simplifyVector = FALSE)
  base <- 8L * ceiling((4L + hl) / 8)
  cols <- lapply(hdr$cols, function(sp) {
    blk <- raw[(base + sp$o + 1L):(base + sp$o + sp$l)]
    switch(
      sp$t,
      f32 = , nlp = readBin(blk, "numeric", n = hdr$n, size = 4L, endian = "little"),
      f64 = readBin(blk, "numeric", n = hdr$n, size = 8L, endian = "little"),
      i32 = readBin(blk, "integer", n = hdr$n, size = 4L, endian = "little"),
      b = as.integer(blk),
      d = {
        codes <- readBin(blk, "integer", n = hdr$n, size = sp$w, signed = FALSE, endian = "little")
        lev <- unlist(sp$k)
        out <- lev[codes + 1L]
        out
      },
      s = strsplit(rawToChar(blk), "\n", fixed = TRUE)[[1L]]
    )
  })
  names(cols) <- vapply(hdr$cols, `[[`, character(1), "name")
  list(n = hdr$n, spec = hdr$cols, cols = cols)
}

# ---- data layer (R side) ------------------------------------------------------------

test_that("the columnar encoder picks the right type and keeps the values", {
  skip_on_cran()
  set.seed(3)
  df <- data.frame(
    gene = c(paste0("G", 1:11), "café"),
    baseMean = c(rexp(11, 1 / 100), NA),
    log2FoldChange = c(rnorm(11), -Inf),
    padj = c(1e-300, 10^-runif(8, 0, 40), 0.05, NA, 0),
    entrezid = c(1:11 * 1000000L, NA),
    big = c(1:11 * 987654321, NA),
    collection = rep(c("GO_BP", "Reactome"), 6),
    flag = c(rep(c(TRUE, FALSE), 5), NA, TRUE),
    stringsAsFactors = FALSE
  )
  enc <- ExpressOM:::.xq_encode_table(df)
  expect_identical(enc$meta$n, 12L)
  types <- stats::setNames(vapply(enc$meta$cols, `[[`, character(1), 2L),
                           vapply(enc$meta$cols, `[[`, character(1), 1L))
  expect_identical(types[["padj"]], "nlp")          # p-values keep 1e-300
  expect_identical(types[["baseMean"]], "f32")
  expect_identical(types[["entrezid"]], "i32")      # exact integers
  expect_identical(types[["big"]], "f64")           # beyond int32: stays exact
  expect_identical(types[["collection"]], "d")      # dictionary-coded
  expect_identical(types[["flag"]], "b")
  expect_identical(types[["gene"]], "s")

  dec <- xq_decode(enc)
  expect_equal(10^-dec$cols$padj[1], 1e-300, tolerance = 1e-5)
  expect_true(is.na(dec$cols$padj[11]))
  expect_equal(dec$cols$baseMean[1:11], df$baseMean[1:11], tolerance = 1e-6)
  expect_true(is.na(dec$cols$baseMean[12]))
  expect_identical(dec$cols$entrezid[c(1, 11)], c(1000000L, 11000000L))
  expect_true(is.na(dec$cols$entrezid[12]))
  expect_identical(dec$cols$big[11], 11 * 987654321)
  expect_identical(dec$cols$collection, df$collection)
  expect_identical(dec$cols$flag[c(1, 2, 11)], c(1L, 0L, 255L))
})

test_that("only p / q columns are stored as -log10 values", {
  expect_true(ExpressOM:::.xq_is_p_col("padj", c(0.5, 1e-10, NA)))
  expect_true(ExpressOM:::.xq_is_p_col("p.adjust", c(0.5, 0.01)))
  expect_true(ExpressOM:::.xq_is_p_col("isoform_switch_q_value", c(0.2, 0.01)))
  expect_false(ExpressOM:::.xq_is_p_col("padj_cutoff", c(0.05, 0.01)))   # a threshold
  expect_false(ExpressOM:::.xq_is_p_col("padj", c(0.5, 12)))             # not a probability
  expect_false(ExpressOM:::.xq_is_p_col("log2FoldChange", c(0.1, 0.2)))
  expect_false(ExpressOM:::.xq_is_p_col("gene", c("a", "b")))
})

test_that("embedded blocks can never close their own <script> element", {
  js <- ExpressOM:::.xq_json_script("k", list(a = "</script><!-- x", b = c(1, NA)))
  expect_false(grepl("</script><!--", js, fixed = TRUE))
  expect_equal(length(gregexpr("</script>", js, fixed = TRUE)[[1L]]), 1L)
  expect_match(js, 'id="xqj-k"', fixed = TRUE)
  tb <- ExpressOM:::.xq_table_script("de", ExpressOM:::.xq_encode_table(data.frame(x = 1:3 + 0.5)))
  expect_match(tb, 'type="application/x-xq-table" id="xqt-de"', fixed = TRUE)
  expect_null(ExpressOM:::.xq_encode_table(data.frame()))
})

test_that("the data layer round-trips through the JavaScript runtime", {
  skip_on_cran()
  node <- Sys.which("node")
  testthat::skip_if(!nzchar(node), "node is not installed")
  xq_dir <- system.file("xq", package = "ExpressOM")
  if (!nzchar(xq_dir)) xq_dir <- normalizePath(testthat::test_path("..", "..", "inst", "xq"))
  testthat::skip_if(!file.exists(file.path(xq_dir, "xq-core.js")), "runtime assets not found")
  set.seed(1)
  df <- data.frame(
    gene = c(paste0("G", 1:48), "Ünï\"cø,de", NA),
    baseMean = c(rexp(49, 1 / 100), NA),
    log2FoldChange = c(rnorm(49), -Inf),
    padj = c(10^-runif(46, 0, 50), 1e-300, 0, 0.05, NA),
    pvalue = c(runif(49), 1),
    entrezid = c(1:48 * 1000000L, 2000000000, NA),
    big = c(1:49 * 123456789, NA),
    collection = rep(c("GO_BP", "Reactome"), 25),
    flag = c(rep(c(TRUE, FALSE), 24), NA, TRUE),
    note = paste("line", 1:50, "\n x"),
    stringsAsFactors = FALSE
  )
  enc <- ExpressOM:::.xq_encode_table(df)
  payload <- withr::local_tempfile(fileext = ".json")
  jsonlite::write_json(list(b64 = enc$b64, meta = enc$meta, df = df), payload,
                       auto_unbox = TRUE, digits = NA, na = "null")
  js <- testthat::test_path("xq", "core_roundtrip.js")
  out <- suppressWarnings(system2(node, c(shQuote(js), shQuote(xq_dir), shQuote(payload)),
                                  stdout = TRUE, stderr = TRUE, timeout = 120))
  expect_false(any(grepl("^FAIL", out)), info = paste(out, collapse = "\n"))
  expect_true(any(grepl("ALL OK", out)), info = paste(out, collapse = "\n"))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
})

# ---- the "How to read this" notes ------------------------------------------------------

test_that("every help key yields a collapsible note and unknown keys are rejected", {
  ctx <- list(level = "treated", base = "control", padj = 0.05, lfc = 1, test = "Wald",
              lrt = FALSE, group = "condition")
  keys <- c("overview", "deg_counts", "contents", "pca", "similarity", "expr_qc", "volcano",
            "diagnostics", "heatmap", "de_table", "gene_lists", "sql_de", "enrichment",
            "decouple", "regulon", "switch", "dte", "dtu", "dexseq", "tx_pca", "iso_plots",
            "iso_methods", "figures", "params", "logs", "files")
  for (k in keys) {
    html <- ExpressOM:::.xq_help_html(k, ctx)
    expect_true(is.character(html) && length(html) == 1L && nzchar(html), info = k)
    expect_match(html, "<details class=\"xq-help\"", fixed = TRUE, info = k)
    expect_match(html, "<li>", fixed = TRUE, info = k)
  }
  # the notes name the contrast and the cut-offs of the run
  v <- ExpressOM:::.xq_help_html("volcano", ctx)
  expect_match(v, "treated", fixed = TRUE)
  expect_match(v, "control", fixed = TRUE)
  expect_match(v, "0.05", fixed = TRUE)
  # each enrichment kind gets its own wording
  expect_match(ExpressOM:::.xq_help_html("enrichment", ctx, kind = "gsea", method = "FGSEA"), "NES", fixed = TRUE)
  expect_match(ExpressOM:::.xq_help_html("enrichment", ctx, kind = "ora", method = "GO ORA"), "over-representation", ignore.case = TRUE)
})

test_that("the dashboard source gives every tab a note and covers every report page", {
  skip_on_cran()
  out <- withr::local_tempdir()
  make_run_dir(out)
  res <- load_expressom_results(out)
  st <- expressom_report_settings()
  col <- ExpressOM:::.xq_collect(res, st)
  src <- ExpressOM:::.xq_dashboard_source(res, st, col, "Quarto test")
  txt <- paste(src, collapse = "\n")
  has <- function(x, pat) grepl(pat, x, fixed = TRUE, useBytes = TRUE)
  for (page in c("# Overview", "# Samples & QC", "# Differential expression",
                 "# Functional analysis", "# Isoforms", "# Run info")) {
    expect_true(has(txt, page), info = page)
  }
  expect_true(has(txt, "embed-resources: true"))
  expect_true(has(txt, "format:\n  dashboard:"))

  # split into cards (### heading ... closing :::) and check each has a note,
  # either inline (static HTML) or registered for the view it mounts
  start <- grep("^### ", src, useBytes = TRUE)
  expect_gt(length(start), 20L)
  registered <- names(ExpressOM:::.xq_reg$help)
  for (i in seq_along(start)) {
    end <- if (i < length(start)) start[[i + 1L]] - 1L else length(src)
    block <- src[start[[i]]:end]
    title <- sub("^### ", "", block[[1L]])
    if (identical(title, "Notes")) next          # loader messages, not a data view
    mounted <- sub('^XQ\\.mount\\("([^"]+)".*$', "\\1",
                   grep('^XQ\\.mount\\(', block, value = TRUE, useBytes = TRUE))
    has_note <- any(grepl("xq-help", block, fixed = TRUE, useBytes = TRUE)) || any(mounted %in% registered)
    expect_true(has_note, info = paste("card without a note:", title))
  }
  # data never travels inside JavaScript source
  expect_false(grepl("JSON.parse\\(", txt))
})

test_that("the head payload carries styles, specs, notes, tables and runtime", {
  skip_on_cran()
  out <- withr::local_tempdir()
  make_run_dir(out)
  res <- load_expressom_results(out)
  st <- expressom_report_settings()
  col <- ExpressOM:::.xq_collect(res, st)
  src <- ExpressOM:::.xq_dashboard_source(res, st, col, "t")
  head <- paste(ExpressOM:::.xq_head_html(col, st), collapse = "\n")
  expect_match(head, 'id="xqt-de"', fixed = TRUE)
  expect_match(head, 'id="xqj-cfg"', fixed = TRUE)
  expect_match(head, 'id="xqh-', fixed = TRUE)
  expect_match(head, "globalThis.XQV", fixed = TRUE)
  expect_match(head, "XQ.mount", fixed = TRUE)
  expect_true(all(c("de", "dte", "dtu") %in% names(col$tables)))
})

# ---- end to end ------------------------------------------------------------------------------

test_that("build_expressom_report(engine = 'quarto') renders a self-contained dashboard", {
  skip_on_cran()
  skip_if_no_quarto()
  out <- withr::local_tempdir()
  make_run_dir(out)
  path <- suppressMessages(build_expressom_report(
    out, engine = "quarto",
    settings = expressom_report_settings(title = "Quarto test")
  ))
  expect_true(file.exists(path))
  html <- rawToChar(readBin(path, "raw", file.size(path)))
  Encoding(html) <- "bytes"
  # self-contained: no external scripts or stylesheets, nothing fetched on load
  expect_false(grepl('<script[^>]+src="http', html))
  expect_false(grepl('<link[^>]+href="http', html))
  expect_false(grepl("fonts.googleapis", html, fixed = TRUE))
  # title, with a real middle dot (not the literal text "\\u00b7")
  expect_match(html, "Quarto test · ", fixed = TRUE)
  expect_false(grepl("\\u00b7", html, fixed = TRUE))
  # binary payloads, specs and notes are in the page
  expect_match(html, 'id="xqt-de"', fixed = TRUE)
  expect_match(html, 'id="xqj-v_volcano"', fixed = TRUE)
  expect_gt(length(gregexpr('id="xqh-', html, fixed = TRUE)[[1L]]), 10L)
  # the file:// guards are neutralised so the report opens from a double-click
  expect_false(grepl("protocol === \"file:\"", html, fixed = TRUE))
  # every page of the htmltools report is represented
  for (tab in c("Overview", "Samples &amp; QC", "Differential expression",
                "Functional analysis", "Isoforms", "Run info")) {
    expect_match(html, tab, fixed = TRUE)
  }
})
