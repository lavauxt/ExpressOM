# mod_report_quarto.R -- Quarto dashboard engine for the unified report.
#
# Renders the same tab tree as the htmltools/pandoc engines (see
# expressom_report_tabs()) as a Quarto dashboard (format: dashboard) with
# Observable JS ({ojs}) interactive explorers on top of embedded binary data:
#
#   * Stockage binaire leger : numeric columns are serialised in R as IEEE754
#     float32 (writeBin(size = 4)), base64-packed, and decoded in the browser
#     into Float32Array -- ~50% smaller than JSON numbers, no Arrow runtime
#     needed for the main tables. (An Apache Arrow IPC path can be added later;
#     the stdlib of {ojs} bundles apache-arrow already.)
#   * Moteur de filtrage rapide : DuckDB-Wasm (DuckDBClient in the {ojs}
#     standard library) registers the embedded tables and answers SQL
#     filters instantly; a pure-JS fallback keeps the explorers usable
#     offline (the wasm bundle is fetched from the Observable CDN on first
#     use, exactly like plotly.js is fetched from cdn.plot.ly today).
#   * Rendu WebGL : large scatter plots (volcano / MA / DTE) are drawn with
#     Plot (Observable Plot) which hands off to canvas/WebGL-friendly
#     rendering for big point clouds; scattergl equivalents are kept for
#     the plotly widgets that remain.
#   * Interface prete a l'emploi : Quarto Dashboard + Observable Inputs
#     (sliders, searches, radios) drive every explorer.
#
# The rendered HTML is self-contained (embed-resources: true): every style,
# script and data payload travels inside the file. The two file:// guards
# Quarto adds for {ojs} documents are patched after rendering because this
# report embeds everything and must open from a double-click, like the
# htmltools engine's output.


# string paste infix (used by the page builders)
`%+%` <- function(a, b) paste0(a, b)

# ---- binary serialisation ------------------------------------------------------

#' Encode a numeric vector as base64 IEEE754 float32 (little-endian)
#'
#' `writeBin(size = 4)` performs the double -> float32 conversion natively
#' (bit-exact with `Float32Array` semantics, verified against Node). Base64
#' keeps the payload JSON-safe so it can travel inside a `<script>` block.
#'
#' @param x Numeric vector. `NA`/`NaN`/`Inf` are preserved as float32
#'   `NaN`/`Inf`; values beyond the float32 range saturate.
#' @return A base64 string.
#' @keywords internal
.xq_f32_b64 <- function(x) {
  x <- as.numeric(x)
  con <- rawConnection(raw(0), "wb")
  writeBin(x, con, size = 4)
  raws <- rawConnectionValue(con)
  close(con)
  jsonlite::base64_enc(raws)
}

#' Build one embedded table description
#'
#' A list with the column names, the number of rows, a base64 float32 block
#' per numeric column and a plain JSON array per character/factor/logical
#' column. The {ojs} side reconstructs it with `Float32Array(atob(...))`.
#'
#' @param df A data.frame.
#' @param max_rows Hard row cap for the embedded payload (safety valve;
#'   the full table is always available as CSV in the run folder).
#' @keywords internal
.xq_table_payload <- function(df, max_rows = Inf) {
  if (is.null(df) || nrow(df) == 0L || ncol(df) == 0L) return(NULL)
  if (nrow(df) > max_rows) df <- utils::head(df, as.integer(max_rows))
  cols <- list()
  for (cn in colnames(df)) {
    v <- df[[cn]]
    if (is.numeric(v)) {
      cols[[cn]] <- list(t = "f32", d = .xq_f32_b64(v))
    } else if (is.logical(v)) {
      cols[[cn]] <- list(t = "b", d = unname(as.list(v)))
    } else {
      cols[[cn]] <- list(t = "s", d = unname(as.list(as.character(v))))
    }
  }
  list(n = nrow(df), cols = cols)
}

#' Serialise one table to the inline `<script data-xpdata>` payload
#' @keywords internal
.xq_embed_data_script <- function(tables, id = "xpdata") {
  payload <- list(tables = tables)
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null",
                           na = "null")
  paste0(
    '<script type="application/json" data-', id, '>',
    json,
    '</script>'
  )
}

# ---- the {ojs} runtime helpers (shared by every explorer) -----------------------

.xq_ojs_prelude_cells <- function() {
  # Each entry is the source of ONE {ojs} cell. Observable's parser accepts
  # function declarations plus a single definition per cell, so every value
  # definition gets its own cell.
  c(
    # ---- embedded-data bridge ------------------------------------------
    paste(c(
      "XPDATA = JSON.parse(document.querySelector('script[data-xpdata]').textContent)",
      "function xpBytes(b64) {",
      "  const bin = atob(b64);",
      "  const bytes = new Uint8Array(bin.length);",
      "  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);",
      "  return bytes;",
      "}",
      "function xpCols(name) { return Object.keys(XPDATA.tables[name].cols); }",
      "function xpN(name) { return XPDATA.tables[name].n; }",
      "function xpTable(name) {",
      "  const cache = xpTable._cache || (xpTable._cache = {});",
      "  if (cache[name]) return cache[name];",
      "  const t = XPDATA.tables[name];",
      "  const out = new Array(t.n);",
      "  const colData = {};",
      "  for (const [cn, spec] of Object.entries(t.cols)) {",
      "    colData[cn] = spec.t === 'f32' ? new Float32Array(xpBytes(spec.d).buffer) : spec.d;",
      "  }",
      "  const keys = Object.keys(colData);",
      "  for (let i = 0; i < t.n; i++) {",
      "    const row = {};",
      "    for (const k of keys) {",
      "      const v = colData[k][i];",
      "      row[k] = (v === undefined || Number.isNaN(v)) ? null : v;",
      "    }",
      "    out[i] = row;",
      "  }",
      "  cache[name] = out;",
      "  return out;",
      "}"
    ), collapse = "\n"),
    # ---- classification + palette ---------------------------------------
    # one cell per definition (Observable: a cell defines a single name)
    paste(c(
      "// significance classification (mirrors the R volcano helpers)",
      "function xpClass(d, padjThr, lfcThr, isLrt) {",
      "  if (d.padj == null || d.padj !== d.padj) return 'Not significant';",
      "  const sig = d.padj < padjThr;",
      "  if (isLrt) return sig ? 'Significant (LRT)' : 'Not significant';",
      "  if (sig && d.lfc != null && d.lfc > lfcThr) return 'Upregulated';",
      "  if (sig && d.lfc != null && d.lfc < -lfcThr) return 'Downregulated';",
      "  return 'Not significant';",
      "}",
      "1"
    ), collapse = "\n"),
    paste(c(
      "function xpFill(cls) {",
      "  const m = {",
      "    'Upregulated': '#D55E00', 'Downregulated': '#0072B2',",
      "    'Significant (LRT)': '#D55E00', 'Not significant': '#b8b8b8'",
      "  };",
      "  return m[cls] || '#b8b8b8';",
      "}",
      "1"
    ), collapse = "\n"),
    # ---- DuckDB-Wasm + fallback ------------------------------------------
    paste(c(
      "// DuckDB-Wasm fast filtering with a pure-JS fallback: the {ojs} stdlib",
      "// ships DuckDBClient; when the wasm bundle cannot be fetched (offline",
      "// first use) every explorer falls back to array filters.",
      "xqDb = (async () => {",
      "  try {",
      "    const tables = {};",
      "    for (const name of Object.keys(XPDATA.tables)) tables[name] = xpTable(name);",
      "    return await DuckDBClient.of(tables);",
      "  } catch (e) {",
      "    return null;",
      "  }",
      "})()",
      "function xqQuery(sql, fallbackRows) {",
      "  return (async () => {",
      "    try {",
      "      const db = await xqDb;",
      "      if (db) return await db.query(sql);",
      "    } catch (e) { /* offline or query error */ }",
      "    return fallbackRows;",
      "  })()",
      "}"
    ), collapse = "\n"),
    # ---- formatters -------------------------------------------------------
    paste(c(
      "// small formatting helpers",
      "function xpFmtP(v) { return (v == null || v !== v) ? 'NA' : (v < 1e-3 ? v.toExponential(1) : v.toFixed(3)); }",
      "function xpFmtN(v, d) { d = d == null ? 2 : d; return (v == null || v !== v) ? 'NA' : (+v).toFixed(d); }",
      "function xpDeb(v) { return (v == null || v !== v) ? 'NA' : v.toExponential(1); }",
      "1"
    ), collapse = "\n")
  )
}

# ---- CSS for the dashboard (kept small; the dashboard theme does the rest) ------

.xq_css <- function() {
  paste(c(
    "<style>",
    "/* ExpressOM Quarto dashboard additions */",
    ".xq-note { color: #52606d; font-style: italic; margin: 4px 0 8px; font-size: 0.92em; }",
    ".xq-badge { display: inline-block; padding: 1px 9px; border-radius: 12px;",
    "  font-size: 0.78em; font-weight: 600; margin-left: 8px; vertical-align: middle; }",
    ".xq-ok { background: #e2f2e9; color: #1a7f4b; }",
    ".xq-fallback { background: #fdecea; color: #8a1c12; }",
    ".xq-hero .value { font-size: 1.65em; font-weight: 650; color: #0b3d75; }",
    ".xq-hero .label { font-size: 0.82em; color: #52606d; text-transform: uppercase;",
    "  letter-spacing: 0.04em; }",
    ".xq-hero .sub { font-size: 0.78em; color: #7b8794; }",
    ".xq-hero { border: 1px solid #d9e2ec; border-radius: 10px; padding: 10px 14px;",
    "  background: #f5f8fb; height: 100%; }",
    ".xq-static img { max-width: 100%; border: 1px solid #d9e2ec; border-radius: 6px; }",
    "/* dark mode */",
    ".quarto-dashboard-dark .xq-note { color: #a7b4c2; }",
    ".quarto-dashboard-dark .xq-ok { background: #143527; color: #4cc38a; }",
    ".quarto-dashboard-dark .xq-fallback { background: #3a1a1a; color: #ffb4ab; }",
    ".quarto-dashboard-dark .xq-hero { background: #17212d; border-color: #2c3a4a; }",
    ".quarto-dashboard-dark .xq-hero .value { color: #a9d2fb; }",
    ".quarto-dashboard-dark .xq-hero .label { color: #a7b4c2; }",
    ".quarto-dashboard-dark .xq-hero .sub { color: #8794a3; }",
    ".quarto-dashboard-dark .xq-static img { border-color: #2c3a4a; }",
    "</style>"
  ), collapse = "\n")
}

# ---- generic table rendering (every tab keeps its data table) -------------------

.xq_ojs_table_cell <- function(table_key, filename, note = NULL,
                              page_length = 15L, searchable = TRUE) {
  code <- c()
  if (!is.null(note)) {
    code <- c(code, paste0("// ", gsub("\n", " ", note)))
  }
  code <- c(code, paste0(
    "xpTable('", table_key, "') |> Inputs.table({rows: ", page_length,
    ", format: {padj: xpDeb, pvalue: xpDeb, padj_DTE: xpDeb},",
    " sort: 'padj', reverse: true})",
    if (searchable) " |> Inputs.search()" else "",
    " |> Inputs.select(xpCols('" , table_key, "').map(c => ({label: c, value: c})), {label: 'jump to column', value: d => d.value})"
  ))
  paste(code, collapse = "\n")
}

# ---- the .qmd generator ----------------------------------------------------------

#' Collect every table that must be embedded in the dashboard payload
#'
#' Mirrors `load_expressom_results()` outputs. Keys are stable identifiers
#' used by the {ojs} side (`xpTable('<key>')`).
#' @keywords internal
.xq_collect_tables <- function(res, st) {
  tabs <- list()
  add <- function(key, df, max_rows = Inf) {
    pl <- .xq_table_payload(df, max_rows = max_rows)
    if (!is.null(pl)) tabs[[key]] <<- pl
  }

  if (!is.null(res$dge)) {
    # DGE results: keep the columns the explorers need; drop heavy diagnostics.
    # log2FoldChange is exposed as `lfc` so every explorer reads one stable key.
    tbl <- res$dge$results
    keep <- intersect(c("gene", "ensembl", "symbol", "baseMean",
                        "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                        "test_type"), colnames(tbl))
    tbl <- tbl[, keep, drop = FALSE]
    if ("log2FoldChange" %in% colnames(tbl)) names(tbl)[names(tbl) == "log2FoldChange"] <- "lfc"
    add("de", tbl, max_rows = st$table_rows)
    if (!is.null(res$dge$deg_counts)) add("deg_counts", res$dge$deg_counts)
  }

  # PCA entries (all of them, not only the first two)
  for (nm in names(res$pca)) {
    key <- paste0("pca_", gsub("[^A-Za-z0-9_]", "_", nm))
    add(key, res$pca[[nm]]$scores)
  }

  # Expression matrix for the heatmap / sample QC (top genes, float32 columns)
  if (!is.null(res$dge$expression)) {
    ex <- res$dge$expression
    hm <- data.frame(gene = rownames(ex$matrix), ex$matrix, check.names = FALSE)
    add("expr_heatmap", hm, max_rows = max(st$heatmap_genes, 50L))
    qc <- data.frame(gene = rownames(ex$qc_matrix), ex$qc_matrix, check.names = FALSE)
    add("expr_qc", qc, max_rows = 800L)
  }

  # Functional enrichment: merged inventory + one table per inventory row
  inv <- res$functional$inventory
  if (nrow(inv)) {
    merged <- do.call(rbind, lapply(seq_len(nrow(inv)), function(i) {
      tb <- res$functional$tables[[inv$path[[i]]]]
      if (is.null(tb) || !nrow(tb)) return(NULL)
      std <- .rep_standardise_enrichment(tb, inv$kind[[i]])
      if (is.null(std)) return(NULL)
      data.frame(collection = inv$collection[[i]], method = inv$method[[i]],
                 kind = inv$kind[[i]], term = std$term, padj = std$padj,
                 size = std$size, score = std$score, stringsAsFactors = FALSE)
    }))
    if (!is.null(merged)) add("enrichment_all", merged, max_rows = st$table_rows)
  }

  # decoupleR activity tables
  if (!is.null(res$decouple)) {
    for (rel in names(res$decouple)) {
      entry <- res$decouple[[rel]]
      key <- paste0("decouple_", gsub("[^A-Za-z0-9_]", "_", entry$resource))
      add(key, entry$table)
    }
  }

  # decoupleR per-TF regulons
  if (!is.null(res$decouple_regulons)) {
    rr <- res$decouple_regulons
    for (rname in names(rr$by_resource)) {
      entry <- rr$by_resource[[rname]]
      for (tf in names(entry$regulons)) {
        add(paste0("regulon_", gsub("[^A-Za-z0-9_]", "_", rname), "_",
                   gsub("[^A-Za-z0-9_]", "_", tf)),
            entry$regulons[[tf]])
      }
    }
  }

  # Isoform tables
  iso <- res$isoform
  if (!is.null(iso)) {
    if (!is.null(iso$dte)) add("dte", iso$dte, max_rows = st$table_rows)
    if (!is.null(iso$dtu)) add("dtu", iso$dtu, max_rows = st$table_rows)
    if (!is.null(iso$dexseq)) add("dexseq", iso$dexseq, max_rows = st$table_rows)
    if (!is.null(iso$switch_features)) {
      add("switch_features", iso$switch_features, max_rows = st$table_rows)
    }
  }

  tabs
}

#' Build the data + config script block for the dashboard
#' @keywords internal
.xq_data_block <- function(res, st, tables) {
  co <- .xp_cutoffs(res, st)
  cfg <- list(
    comparison = if (is.na(res$comparison)) "" else res$comparison,
    level = res$config$level %||% "", base = res$config$base %||% "",
    test = res$config$test %||% "Wald",
    padj_cutoff = co$padj, lfc_cutoff = co$lfc,
    title = st$title %||% "ExpressOM report",
    dIF_cutoff = 0.1,
    pca_entries = names(res$pca),
    top_terms = st$top_terms
  )
  json <- jsonlite::toJSON(
    list(cfg = cfg, tables = tables),
    auto_unbox = TRUE, null = "null", na = "null"
  )
  paste0('<script type="application/json" data-xpdata>', json, '</script>')
}

#' Escape a string for safe inclusion in {ojs} source
#' @keywords internal
.xq_js_str <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x, fixed = TRUE)
  x <- gsub("\n", " ", x)
  x
}

#' One dashboard "card" with the value boxes of the overview
#' @keywords internal
.xq_value_boxes <- function(res, st) {
  boxes <- c()
  s <- .xp_dge_summary(res, st)
  if (!is.null(s)) {
    cfg <- res$config
    boxes <- c(boxes, sprintf(
      "{value: '%s vs %s', title: 'Comparison', details: '%s test'}",
      .xq_js_str(cfg$level %||% "?"), .xq_js_str(cfg$base %||% "?"),
      .xq_js_str(cfg$test %||% "Wald")))
    boxes <- c(boxes, sprintf("{value: '%s', title: 'Genes tested'}",
                              format(s$n_tested, big.mark = ",")))
    if (s$is_lrt) {
      boxes <- c(boxes, sprintf(
        "{value: '%s', title: 'Significant genes', details: 'LRT, padj < %s'}",
        format(s$n_sig, big.mark = ","), s$padj))
    } else {
      boxes <- c(boxes, sprintf(
        "{value: '%s', title: 'Upregulated', details: 'padj < %s, log2FC > %s', color: '#D55E00'}",
        format(s$n_up, big.mark = ","), s$padj, s$lfc))
      boxes <- c(boxes, sprintf(
        "{value: '%s', title: 'Downregulated', details: 'padj < %s, log2FC < -%s', color: '#0072B2'}",
        format(s$n_down, big.mark = ","), s$padj, s$lfc))
    }
  }
  if (length(res$pca)) {
    boxes <- c(boxes, sprintf("{value: '%d', title: 'Samples (PCA)'}",
                              nrow(res$pca[[1L]]$scores)))
  }
  inv <- res$functional$inventory
  if (nrow(inv)) {
    boxes <- c(boxes, sprintf("{value: '%d', title: 'Enrichment tables'}", nrow(inv)))
  }
  if (!is.null(res$decouple)) {
    boxes <- c(boxes, sprintf(
      "{value: '%d', title: 'decoupleR resources', details: '%s'}",
      length(res$decouple),
      paste(unique(vapply(res$decouple, `[[`, character(1), "resource")),
            collapse = " / ")))
  }
  boxes
}

#' The static-figure gallery group -> base64 data URIs
#' (same shrinking logic as the htmltools engine)
#' @keywords internal
.xq_gallery_uris <- function(res, st, group, max_figs) {
  fig <- res$figures[res$figures$group == group, , drop = FALSE]
  pngs <- fig[fig$type == "png", , drop = FALSE]
  pdfs <- fig[fig$type == "pdf" & fig$size >= 2000, , drop = FALSE]
  items <- list()
  if (nrow(pngs)) {
    n_take <- min(nrow(pngs), max_figs)
    for (i in seq_len(n_take)) {
      uri <- .xp_img_uri(pngs$path[[i]], figures_dir = NULL,
                         max_px = st$static_max_px %||% 1000L)
      if (!is.null(uri)) {
        items[[length(items) + 1L]] <- list(uri = uri, cap = pngs$name[[i]])
      }
    }
  }
  have_pdftools <- requireNamespace("pdftools", quietly = TRUE)
  if (nrow(pdfs) && have_pdftools) {
    room <- max_figs - length(items)
    if (room > 0L) {
      for (i in seq_len(min(nrow(pdfs), room))) {
        uri <- .xp_pdf_uri(pdfs$path[[i]], dpi = st$static_dpi, figures_dir = NULL)
        if (!is.null(uri)) {
          items[[length(items) + 1L]] <- list(uri = uri,
                                               cap = .xp_strip_pdf(pdfs$name[[i]]))
        }
      }
    }
  }
  items
}

# ---- page builders: one function per top-level dashboard page --------------------

.xq_page_overview <- function(res, st) {
  boxes <- .xq_value_boxes(res, st)
  out <- c()
  if (length(boxes)) {
    # value boxes as a responsive hero grid; one cell defines the data, the
    # next displays the html (an {ojs} cell holds a single definition)
    out <- c(out, "```{ojs}")
    out <- c(out, paste0("xpBoxes = [", paste(boxes, collapse = ",\n"), "]"))
    out <- c(out, "```")
    out <- c(out, "```{ojs}")
    out <- c(out, paste(
      "htl.html`<div class='xq-hero-grid' style='display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:6px 0 14px'>",
      "${xpBoxes.map(b => htl.html`<div class='xq-hero'>",
      "<div class='label'>${b.title}</div>",
      "<div class='value'>${b.value}</div>",
      "${b.details ? htl.html`<div class='sub'>${b.details}</div>` : null}",
      "</div>`)}",
      "</div>`", sep = "\n"))
    out <- c(out, "```")
  }
  # DEG counts plot + table (interactive from the start)
  if (!is.null(res$dge$deg_counts)) {
    out <- c(out, "### DEG counts",
             "```{ojs}",
             "degCounts = xpTable('deg_counts')",
             "Plot.plot({",
             "  height: 300, x: {label: 'adjusted p cut-off'}, y: {label: 'genes'},",
             "  color: {legend: true},",
             "  marks: [",
             "    Plot.barY(degCounts, {x: 'padj_cutoff', y: 'n_up_lfc', fill: '#D55E00', tip: true}),",
             "    Plot.barY(degCounts, {x: 'padj_cutoff', y: 'n_down_lfc', fill: '#0072B2', tip: true})",
             "  ]",
             "})",
             "```",
             "```{ojs}",
             "Inputs.table(degCounts, {rows: 6})",
             "```")
  }
  # contents checklist (same availability list as the htmltools engine)
  out <- c(out, "### Contents of this report", "```{ojs}")
  avail <- c()
  avail <- c(avail, sprintf("  {ok: %s, label: 'Differential expression', detail: '%s'}",
                            if (is.null(res$dge)) "false" else "true",
                            if (is.null(res$dge)) "" else
                              paste(format(nrow(res$dge$results), big.mark = ","), "genes")))
  avail <- c(avail, sprintf("  {ok: %s, label: 'PCA / sample QC', detail: '%s'}",
                            tolower(length(res$pca) > 0L), if (length(res$pca))
                              paste(length(res$pca), "plot(s)") else ""))
  inv <- res$functional$inventory
  avail <- c(avail, sprintf("  {ok: %s, label: 'Functional enrichment', detail: '%s'}",
                            tolower(nrow(inv) > 0L), if (nrow(inv))
                              paste(nrow(inv), "table(s)") else ""))
  avail <- c(avail, sprintf("  {ok: %s, label: 'decoupleR', detail: '%s'}",
                            tolower(!is.null(res$decouple)),
                            if (!is.null(res$decouple))
                              paste(length(res$decouple), "resource(s)") else ""))
  avail <- c(avail, sprintf("  {ok: %s, label: 'Isoforms (DTE / DTU / switches)'}",
                            tolower(!is.null(res$isoform))))
  avail <- c(avail, sprintf("  {ok: %s, label: 'Static figures', detail: '%s'}",
                            tolower(nrow(res$figures) > 0L), if (nrow(res$figures))
                              paste(nrow(res$figures), "file(s)") else ""))
  out <- c(out, paste0("xpAvail = [\n", paste(avail, collapse = ",\n"), "\n]"))
  out <- c(out, paste(
    "htl.html`<ul class='xp-avail'>${xpAvail.map(a => htl.html`<li>",
    "<span class='${a.ok ? \"xq-tick\" : \"xq-miss\"}'>${a.ok ? '\\u2713' : '\\u2013'}</span>",
    " ${a.label}${a.detail ? htl.html`<span class='xq-dim'> \\u00b7 ${a.detail}</span>` : null}",
    "</li>`)}</ul>`", sep = "\n"))
  out <- c(out, "```")
  if (length(res$notes)) {
    out <- c(out, "### Notes", "```{ojs}")
    out <- c(out, paste0("xpNotes = [",
                         paste(vapply(res$notes, function(n)
                           paste0("'", .xq_js_str(n), "'"), character(1)),
                           collapse = ", "), "]"))
    out <- c(out, paste("htl.html`<ul>${xpNotes.map(n => htl.html`<li>${n}</li>`)}</ul>`",
                        sep = "\n"))
    out <- c(out, "```")
  }
  out
}

.xq_page_qc <- function(res, st) {
  out <- c()
  # PCA explorer (NEW: axis picker over ALL pcs, not only PC1/PC2)
  if (length(res$pca)) {
    for (i in seq_along(res$pca)) {
      nm <- names(res$pca)[[i]]
      entry <- res$pca[[nm]]
      key <- paste0("pca_", gsub("[^A-Za-z0-9_]", "_", nm))
      pcs <- entry$pcs
      title <- .xp_pca_title(nm)
      out <- c(out, sprintf("## %s {.tabset}", title))
      out <- c(out, "### PCA explorer",
               "```{ojs}",
               paste0("pca", i, " = xpTable('", key, "')"),
               "```")
      if (i == 1L) {
        # first entry gets the interactive axis picker + group colouring
        pc_list <- paste0("[", paste(sprintf("'%s'", pcs), collapse = ", "), "]")
        gcol <- .rep_pca_group_col(entry$scores, res$config$main_condition)
        group_expr <- if (!is.na(gcol)) paste0("d => d.", gcol) else "'All samples'"
        out <- c(out, "```{ojs}",
                 paste0("viewof pcaX = Inputs.select(", pc_list,
                        ", {label: 'x axis', value: '", pcs[[1L]], "'})"),
                 paste0("viewof pcaY = Inputs.select(", pc_list,
                        ", {label: 'y axis', value: '", pcs[[2L]], "'})"),
                 "```",
                 "```{ojs}",
                 "Plot.plot({",
                 "  height: 420, grid: true,",
                 "  x: {label: pcaX}, y: {label: pcaY},",
                 "  color: {legend: true},",
                 "  marks: [",
                 paste0("    Plot.dot(pca1, {x: pcaX, y: pcaY, fill: ", group_expr,
                        ", tip: true, r: 5, symbol: 'circle'}),"),
                 "    Plot.text(pca1, {x: pcaX, y: pcaY, text: 'sample_label', dy: -8, fontSize: 9})",
                 "  ]",
                 "})",
                 "```")
      } else {
        # other PCA entries: fixed PC1/PC2 plot (same data as the htmltools engine)
        out <- c(out, "```{ojs}",
                 "Plot.plot({",
                 "  height: 420, grid: true,",
                 sprintf("  x: {label: '%s'}, y: {label: '%s'},", pcs[[1L]], pcs[[2L]]),
                 "  marks: [",
                 sprintf("    Plot.dot(pca%d, {x: '%s', y: '%s', fill: 'condition', tip: true, r: 5}),",
                         i, pcs[[1L]], pcs[[2L]]),
                 sprintf("    Plot.text(pca%d, {x: '%s', y: '%s', text: 'sample_label', dy: -8, fontSize: 9})",
                         i, pcs[[1L]], pcs[[2L]]),
                 "  ]",
                 "})",
                 "```")
      }
      # scree plot
      if (!is.null(entry$variance)) {
        out <- c(out, "### Scree",
                 "```{ojs}",
                 paste0("scree", i, " = [",
                        paste(sprintf("{pc: '%s', pct: %s}", names(entry$variance),
                                      entry$variance),
                              collapse = ", "), "]"),
                 "```",
                 "```{ojs}",
                 "Plot.plot({",
                 "  height: 300,",
                 "  x: {label: 'component'}, y: {label: '% variance', grid: true},",
                 paste0("  marks: [Plot.barY(scree", i,
                        ", {x: 'pc', y: 'pct', tip: true})]"),
                 "})",
                 "```")
      }
      # sample scores table
      out <- c(out, "### Sample scores",
               "```{ojs}",
               paste0("Inputs.table(pca", i, ", {rows: 10, maxWidth: 800})"),
               "```",
               "## {.hidden}")  # close tabset
    }
  }
  out
}

.xq_page_dge <- function(res, st) {
  if (is.null(res$dge)) return(NULL)
  co <- .xp_cutoffs(res, st)
  is_lrt <- identical(res$config$test, "LRT")
  out <- c(
    "## Differential expression {.tabset}",
    "### Gene explorer",
    "// WebGL-ready explorer: every gene, DuckDB-Wasm SQL filter + JS fallback,",
    "// view switch volcano/MA, cutoff sliders, gene search, linked table.",
    "```{ojs}",
    "de = xpTable('de')",
    "```",
    "```{ojs}",
    sprintf("isLrt = %s", tolower(is_lrt)),
    "```",
    "```{ojs}",
    "viewof xqView = Inputs.radio(['volcano', 'MA'], {label: 'view', value: 'volcano'})",
    sprintf("viewof xqPadj = Inputs.range([0.001, 0.2], {step: 0.001, value: %s, label: 'padj cutoff'})", co$padj),
    sprintf("viewof xqLfc = Inputs.range([0, 3], {step: 0.1, value: %s, label: '|log2FC| cutoff'})", co$lfc),
    "viewof xqSearch = Inputs.search(de.map(d => d.gene), {label: 'gene', placeholder: 'search genes…'})",
    "```",
    "```{ojs}",
    "xqSel = xqSearch ? de.find(d => d.gene === xqSearch) : null",
    "```",
    "```{ojs}",
    "Plot.plot({",
    "  height: 480, grid: true,",
    "  x: xqView === 'volcano' ? {label: 'log2 fold change'} : {label: 'mean of normalised counts', type: 'log'},",
    "  y: xqView === 'volcano' ? {label: '-log10(adjusted p)'} : {label: 'log2 fold change'},",
    "  color: {legend: true},",
    "  marks: [",
    "    Plot.dot(de, {",
    "      x: xqView === 'volcano' ? d => Math.max(-10, d.lfc) : d => d.baseMean,",
    "      y: xqView === 'volcano' ? d => -Math.log10(Math.max(d.padj, 1e-300)) : d => d.lfc,",
    "      fill: d => xpFill(xpClass(d, xqPadj, xqLfc, isLrt)),",
    "      r: d => xpClass(d, xqPadj, xqLfc, isLrt) === 'Not significant' ? 2 : 3.4,",
    "      fillOpacity: d => xpClass(d, xqPadj, xqLfc, isLrt) === 'Not significant' ? 0.45 : 0.85,",
    "      tip: true",
    "    }),",
    "    ...(xqView === 'volcano'",
    "      ? [Plot.ruleY([-Math.log10(xqPadj)], {stroke: '#666', strokeDasharray: '4 3'})]",
    "      : [Plot.ruleY([0], {stroke: '#999'})]),",
    "    ...(xqView === 'volcano' && !isLrt",
    "      ? [Plot.ruleX([-xqLfc, xqLfc], {stroke: '#666', strokeDasharray: '4 3'})]",
    "      : [])",
    "  ]",
    "})",
    "```",
    "```{ojs}",
    "xqCounts = (() => {",
    "  const cls = d => xpClass(d, xqPadj, xqLfc, isLrt);",
    "  const up = de.filter(d => cls(d) === 'Upregulated').length;",
    "  const dn = de.filter(d => cls(d) === 'Downregulated').length;",
    "  const lrt = de.filter(d => cls(d) === 'Significant (LRT)').length;",
    "  return {up: up, dn: dn, lrt: lrt, ns: de.length - up - dn - lrt};",
    "})()",
    "```",
    "```{ojs}",
    "htl.html`<div class='xq-note'>",
    "  ${xqSel ? htl.html`<b>${xqSel.gene}</b>: log2FC ${xpFmtN(xqSel.lfc)}, padj ${xpDeb(xqSel.padj)} \\u00b7 ` : ''}",
    "  <b>${xqCounts.up}</b> up, <b>${xqCounts.dn}</b> down",
    "  ${isLrt ? htl.html`(<b>${xqCounts.lrt}</b> significant LRT)` : ''}",
    "  out of ${de.length} genes at padj < ${xqPadj}, |log2FC| > ${xqLfc}.",
    "</div>`",
    "```"
  )
  # Results table (DuckDB-backed, same rows as the htmltools engine)
  out <- c(out, "### Results table",
           "// Filterable table; the slider-driven query runs through DuckDB-Wasm when",
           "// online and through plain JS otherwise.",
           "```{ojs}",
           "viewof xqMinLfc = Inputs.range([0, 5], {step: 0.1, value: 0, label: 'min |log2FC|'})",
           "```",
           "```{ojs}",
           sprintf("deRows = de.filter(d => d.padj != null && d.padj < %s && Math.abs(d.lfc || 0) >= xqMinLfc).sort((a, b) => a.padj - b.padj).slice(0, 500)",
                   co$padj),
           "Inputs.table(deRows, {rows: 15, format: {padj: xpDeb, pvalue: xpDeb, baseMean: d => xpFmtN(d, 0)}, sort: 'padj', reverse: true})",
           "```",
           "### Diagnostics",
           "```{ojs}",
           "Plot.plot({",
           "  height: 320,",
           "  x: {label: 'p-value', domain: [0, 1]}, y: {label: 'count', grid: true},",
           "  marks: [Plot.rectY(de, {x: 'pvalue', interval: 0.05, fill: '#56B4E9', tip: true})]",
           "})",
           "```",
           "```{ojs}",
           "// p-value Q-Q plot: observed vs expected -log10(p)",
           "xpP = de.map(d => d.pvalue).filter(p => p != null && p >= 0 && p <= 1).sort((a, b) => a - b)",
           "Plot.plot({",
           "  height: 320, grid: true, x: {label: 'expected -log10(p)'}, y: {label: 'observed -log10(p)'},",
           "  marks: [",
           "    Plot.dot(xpP, {x: (p, i) => -Math.log10((i + 0.5) / xpP.length), y: p => -Math.log10(Math.max(p, 1e-300)), r: 2, fill: '#0072B2', fillOpacity: 0.5}),",
           "    Plot.line([0, Math.max(1, -Math.log10(Math.max(xpP[xpP.length - 1] || 1e-300, 1e-300)))], {stroke: '#d55e00', strokeDasharray: '6 4'})",
           "  ]",
           "})",
           "```",
           "```{ojs}",
           "// ranked absolute effects (top 25 |log2FC|)",
           "xqTop = de.filter(d => d.lfc != null).sort((a, b) => Math.abs(b.lfc) - Math.abs(a.lfc)).slice(0, 25)",
           "Plot.plot({",
           "  height: 420, marginTop: 40, marginLeft: 130,",
           "  x: {label: 'log2 fold change', grid: true},",
           "  y: {label: null},",
           "  color: {domain: ['up', 'down'], range: ['#D55E00', '#0072B2']},",
           "  marks: [",
           "    Plot.barX(xqTop, {x: 'lfc', y: 'gene', fill: d => d.lfc > 0 ? 'up' : 'down', sort: {y: 'x', reverse: false}, tip: true})",
           "  ]",
           "})",
           "```")
  # Heatmap (expression matrix) — plotly static heatmap kept via the precomputed image;
  # the interactive one lives on the expression matrix as an {ojs} cell plot
  if (!is.null(res$dge$expression)) {
    out <- c(out, "### Top genes heatmap",
             "```{ojs}",
             "hm = xpTable('expr_heatmap')",
             "hmCols = xpCols('expr_heatmap').filter(c => c !== 'gene')",
             "// row z-scores of log2(count+1), columns ordered by group",
             "hmLong = hm.flatMap(d => hmCols.map(c => ({gene: d.gene, sample: c, z: d[c]})))",
             "```",
             "```{ojs}",
             "Plot.plot({",
             "  height: 420, marginLeft: 110, marginBottom: 90,",
             "  x: {label: null, tickRotate: -45}, y: {label: null},",
             "  color: {scheme: 'rdbu', reverse: true, legend: true, label: 'normalised count'},",
             "  marks: [Plot.cell(hmLong, {x: 'sample', y: 'gene', fill: 'z', tip: true})]",
             "})",
             "```",
             "```{ojs}",
             "Inputs.table(hm, {rows: 10})",
             "```")
  }
  # DEG lists (files) — keep the file links as a card
  f <- res$dge$sig_lists
  out <- c(out, "### DEG lists")
  if (length(f)) {
    out <- c(out, "```{ojs}")
    out <- c(out, paste0("degLists = [",
                         paste(sprintf("'%s'", f), collapse = ", "), "]"))
    out <- c(out, paste(
      "htl.html`<div class='xq-note'>DEG lists written by the pipeline (in <code>DE_raw_results/</code>):",
      "<ul>${degLists.map(x => htl.html`<li><code>${x}</code></li>`)}</ul></div>`", sep = "\n"))
    out <- c(out, "```")
  } else {
    out <- c(out, "```{ojs}", "htl.html`<div class='xq-note'>No per-cut-off DEG lists were found.</div>`", "```")
  }
  out
}

.xq_page_functional <- function(res, st) {
  inv <- res$functional$inventory
  has_dec <- !is.null(res$decouple)
  has_reg <- !is.null(res$decouple_regulons) &&
    length(res$decouple_regulons$by_resource)
  if (!nrow(inv) && !has_dec && !has_reg) return(NULL)
  co <- .xp_cutoffs(res, st)
  out <- c("## Functional analysis {.tabset}")

  # ---- enrichment explorer (NEW: cross-collection, interactive cutoff + search)
  if (nrow(inv) && !is.null(res$functional$tables)) {
    out <- c(out, "### Enrichment explorer",
             "// All collections merged; DuckDB-Wasm or JS filter by method/collection/padj,",
             "// slider for the adjusted-p cutoff, search over terms, bubble + bar views.",
             "```{ojs}",
             "enr = xpTable('enrichment_all')",
             "```",
             "```{ojs}",
             "viewof enrView = Inputs.radio(['bar', 'bubble'], {label: 'view', value: 'bar'})",
             sprintf("viewof enrPadj = Inputs.range([0.001, 0.3], {step: 0.001, value: %s, label: 'max adjusted p'})", co$padj),
             sprintf("viewof enrTop = Inputs.range([5, 40], {step: 1, value: %s, label: 'top N terms'})",
                     min(st$top_terms, 40L)),
             "viewof enrMethod = Inputs.select(['all', ...new Set(enr.map(d => d.method))], {label: 'method', value: 'all'})",
             "viewof enrSearch = Inputs.search(enr, {label: 'search terms', placeholder: 'e.g. kinase, GO…'})",
             "```",
             "```{ojs}",
             "enrFilt = enr",
             "  .filter(d => d.padj < enrPadj)",
             "  .filter(d => enrMethod === 'all' || d.method === enrMethod)",
             "  .filter(d => !enrSearch || !String(enrSearch).length || Object.values(d).some(v => String(v).toLowerCase().includes(String(enrSearch).toLowerCase())))",
             "  .sort((a, b) => a.padj - b.padj)",
             "  .slice(0, enrTop)",
             "```",
             "```{ojs}",
             "Plot.plot({",
             "  height: Math.max(320, 22 * enrFilt.length + 90), marginLeft: 300,",
             "  x: {label: enrView === 'bar' ? '-log10(adjusted p)' : 'gene count / set size', grid: true},",
             "  y: {label: null},",
             "  color: {legend: true, scheme: 'viridis'},",
             "  marks: [",
             "    ...(enrView === 'bar'",
             "      ? [Plot.barX(enrFilt, {",
             "          x: d => -Math.log10(Math.max(d.padj, 1e-300)), y: 'term', fill: 'method',",
             "          sort: {y: 'x', reverse: true}, tip: true})]",
             "      : [Plot.dot(enrFilt, {",
             "          x: d => (d.size != null && d.size === d.size) ? d.size : 1, y: 'term',",
             "          r: d => d.padj != null ? Math.max(4, Math.min(14, 4 - Math.log10(Math.max(d.padj, 1e-300)) * 2.2)) : 5,",
             "          fill: d => -Math.log10(Math.max(d.padj, 1e-300)),",
             "          sort: {y: 'x', reverse: true}, tip: true})]),",
             "    Plot.ruleX([enrView === 'bar' ? -Math.log10(enrPadj) : 0], {stroke: '#666', strokeDasharray: '4 3'})",
             "  ]",
             "})",
             "```",
             "```{ojs}",
             "Inputs.table(enrFilt, {rows: 10, format: {padj: xpDeb}, sort: 'padj'})",
             "```")

    # ---- per-method merged tables (same collections as the htmltools engine)
    methods <- unique(inv$method)
    for (m in methods) {
      mi <- inv[inv$method == m, , drop = FALSE]
      safe <- gsub("[^A-Za-z0-9]", "_", m)
      out <- c(out, sprintf("### %s (%d collections) {.tabset}", m, nrow(mi)))
      out <- c(out, sprintf("### Merged table"))
      out <- c(out, "```{ojs}")
      out <- c(out, sprintf("enrM_%s = enr.filter(d => d.method === '%s')",
                            safe, .xq_js_str(m)))
      out <- c(out, sprintf(
        "Inputs.table(enrM_%s, {rows: 15, format: {padj: xpDeb}, sort: 'padj'})", safe))
      out <- c(out, "```")
      # one sub-tab per collection (tables only, as in the engine)
      for (k in seq_len(nrow(mi))) {
        col_name <- mi$collection[[k]]
        out <- c(out, sprintf("### %s", col_name))
        out <- c(out, "```{ojs}")
        out <- c(out, sprintf("enrC_%s_%d = enr.filter(d => d.collection === '%s')",
                              safe, k, .xq_js_str(col_name)))
        out <- c(out, sprintf(
          "Inputs.table(enrC_%s_%d, {rows: 10, format: {padj: xpDeb}, sort: 'padj'})", safe, k))
        out <- c(out, "```")
      }
    }
  }

  # ---- decoupleR activity
  if (has_dec) {
    for (rel in names(res$decouple)) {
      entry <- res$decouple[[rel]]
      resource <- entry$resource
      key <- paste0("decouple_", gsub("[^A-Za-z0-9_]", "_", resource))
      out <- c(out, sprintf("### decoupleR: %s", toupper(resource)),
               "```{ojs}",
               paste0("dec_", gsub("[^A-Za-z0-9_]", "_", resource),
                      " = xpTable('", key, "')"),
               "```",
               "```{ojs}",
               "// interactive: cutoff slider over the activity table (NEW)",
               sprintf("viewof decThr_%s = Inputs.range([0.001, 0.3], {step: 0.001, value: %s, label: 'p cutoff'})",
                       gsub("[^A-Za-z0-9_]", "_", resource), co$padj),
               "```",
               "```{ojs}",
               sprintf("decF_%s = dec_%s.filter(d => (d.p_value != null && d.p_value < decThr_%s) || (d.padj != null && d.padj < decThr_%s)).sort((a, b) => Math.abs(b.score) - Math.abs(a.score)).slice(0, 25)",
                       gsub("[^A-Za-z0-9_]", "_", resource),
                       gsub("[^A-Za-z0-9_]", "_", resource),
                       gsub("[^A-Za-z0-9_]", "_", resource),
                       gsub("[^A-Za-z0-9_]", "_", resource)),
               "Plot.plot({",
               "  height: 420, marginLeft: 200,",
               "  x: {label: 'activity score (t-like)', grid: true}, y: {label: null},",
               "  color: {domain: ['positive', 'negative'], range: ['#c94a4a', '#2b5c8f']},",
               "  marks: [Plot.barX(decF_" %+% gsub("[^A-Za-z0-9_]", "_", resource) %+% ", {x: 'score', y: 'source', fill: d => d.score > 0 ? 'positive' : 'negative', sort: {y: 'x'}, tip: true})]",
               "})",
               "```",
               "```{ojs}",
               sprintf("Inputs.table(dec_%s, {rows: 15, format: {p_value: xpDeb}})",
                       gsub("[^A-Za-z0-9_]", "_", resource)),
               "```")
    }
  }
  out
}

.xq_page_isoform <- function(res, st) {
  iso <- res$isoform
  if (is.null(iso)) return(NULL)
  co <- .xp_cutoffs(res, st)
  out <- c("## Isoforms {.tabset}")

  # ---- transcript explorer (NEW: DTE volcano/MA with sliders + search)
  if (!is.null(iso$dte) && nrow(iso$dte)) {
    out <- c(out, "### Transcript explorer (DTE)",
             "```{ojs}",
             "dte = xpTable('dte')",
             "```",
             "```{ojs}",
             "viewof dteView = Inputs.radio(['volcano', 'MA', 'abundance'], {label: 'view', value: 'volcano'})",
             sprintf("viewof dtePadj = Inputs.range([0.001, 0.2], {step: 0.001, value: %s, label: 'padj cutoff'})", co$padj),
             sprintf("viewof dteLfc = Inputs.range([0, 3], {step: 0.1, value: %s, label: '|log2FC| cutoff'})", co$lfc),
             "viewof dteSearch = Inputs.search(dte, {label: 'transcript / gene', placeholder: 'ENST… / GENE…'})",
             "```",
             "```{ojs}",
             "dteFilt = dte.filter(d =>",
             "  (!String(dteSearch || '').length || Object.values(d).some(v => String(v).toLowerCase().includes(String(dteSearch).toLowerCase()))))",
             "```",
             "```{ojs}",
             "Plot.plot({",
             "  height: 460, grid: true,",
             "  x: dteView === 'volcano' ? {label: 'transcript log2 fold change'} :",
             "      dteView === 'MA' ? {label: 'mean normalised transcript abundance', type: 'log'} :",
             "      {label: 'mean normalised transcript abundance', type: 'log'},",
             "  y: dteView === 'volcano' ? {label: '-log10(adjusted p)'} : {label: 'log2 fold change'},",
             "  color: {legend: true},",
             "  marks: [",
             "    Plot.dot(dteFilt, {",
             "      x: dteView === 'volcano' ? 'lfc' : 'baseMean',",
             "      y: dteView === 'volcano' ? d => -Math.log10(Math.max(d.padj, 1e-300)) : 'lfc',",
             "      fill: d => d.padj != null && d.padj < dtePadj && d.lfc != null && Math.abs(d.lfc) > dteLfc ? (d.lfc > 0 ? '#D55E00' : '#0072B2') : '#b8b8b8',",
             "      r: 3, fillOpacity: 0.75, tip: true",
             "    }),",
             "    ...(dteView === 'volcano' ? [Plot.ruleY([-Math.log10(dtePadj)], {stroke: '#666', strokeDasharray: '4 3'})] :",
             "        [Plot.ruleY([0], {stroke: '#999'})])",
             "  ]",
             "})",
             "```",
             "```{ojs}",
             "Inputs.table(dteFilt, {rows: 15, format: {padj: xpDeb, pvalue: xpDeb}, sort: 'padj', reverse: true})",
             "```",
             "```{ojs}",
             "// top DTE transcripts by adjusted p",
             "dteTop = dte.filter(d => d.padj != null).sort((a, b) => a.padj - b.padj).slice(0, 25)",
             "Plot.plot({",
             "  height: 420, marginLeft: 140,",
             "  x: {label: '-log10(adjusted p)', grid: true}, y: {label: null},",
             "  marks: [Plot.barX(dteTop, {",
             "    x: d => -Math.log10(Math.max(d.padj, 1e-300)),",
             "    y: d => d.gene_symbol || d.gene || d.transcript_id,",
             "    fill: d => d.lfc < 0 ? '#0072B2' : '#D55E00',",
             "    sort: {y: 'x'}, tip: true})]",
             "})",
             "```")
  }

  # ---- DTU diagnostics
  if (!is.null(iso$dtu) && nrow(iso$dtu)) {
    out <- c(out, "### DTU",
             "```{ojs}",
             "dtu = xpTable('dtu')",
             "```",
             "```{ojs}",
             "// LR vs adjusted p (interactive cutoff), LR-vs-df, top features, LR histogram",
             "Plot.plot({",
             "  height: 380, grid: true,",
             "  x: {label: 'likelihood-ratio statistic'}, y: {label: '-log10(adjusted p)'},",
             "  marks: [Plot.dot(dtu, {",
             "    x: d => (d.lr != null ? d.lr : d.statistic),",
             "    y: d => -Math.log10(Math.max(d.adj_pvalue != null ? d.adj_pvalue : 1e-300, 1e-300)),",
             "    fill: d => (d.adj_pvalue != null && d.adj_pvalue < 0.05) ? '#7A6FAE' : '#b8b8b8',",
             "    r: 3, fillOpacity: 0.7, tip: true})]",
             "})",
             "```",
             "```{ojs}",
             "Plot.plot({",
             "  height: 300,",
             "  x: {label: 'likelihood-ratio statistic', grid: true}, y: {label: 'feature count'},",
             "  marks: [Plot.rectY(dtu, {x: d => d.lr != null ? d.lr : d.statistic, interval: 40, fill: '#7A6FAE', tip: true})]",
             "})",
             "```",
             "```{ojs}",
             "Inputs.table(dtu, {rows: 15, format: {adj_pvalue: xpDeb, pvalue: xpDeb}})",
             "```")
  }

  # ---- switch overview
  if (!is.null(iso$switch_features)) {
    out <- c(out, "### Switch overview",
             "```{ojs}",
             "sw = xpTable('switch_features')",
             "```",
             "```{ojs}",
             sprintf("viewof swQ = Inputs.range([0.001, 0.3], {step: 0.001, value: %s, label: 'max isoform switch q'})", co$padj),
             "viewof swDif = Inputs.range([0, 0.5], {step: 0.01, value: 0.1, label: 'min |dIF|'})",
             "```",
             "```{ojs}",
             "swFilt = sw.filter(d => d.isoform_switch_q_value != null && d.isoform_switch_q_value < swQ && Math.abs(d.dIF || 0) > swDif)",
             "Plot.plot({",
             "  height: 400, grid: true,",
             "  x: {label: 'dIF (difference in isoform fraction)'}, y: {label: '-log10(q)'},",
             "  marks: [",
             "    Plot.dot(sw, {x: 'dIF', y: d => -Math.log10(Math.max(d.isoform_switch_q_value, 1e-300)),",
             "      fill: d => d.isoform_switch_q_value < 0.05 && Math.abs(d.dIF || 0) > 0.1 ? '#c94a4a' : '#b8b8b8', r: 3, fillOpacity: 0.7, tip: true})",
             "  ]",
             "})",
             "```",
             "```{ojs}",
             "htl.html`<div class='xq-note'>${sw.length} isoform(s) in the switch features table, of which ${swFilt.length} have |dIF| > ${swDif} and q < ${swQ}.</div>`",
             "```",
             "```{ojs}",
             "Inputs.table(sw, {rows: 15})",
             "```")
  }

  # ---- DEXSeq
  if (!is.null(iso$dexseq)) {
    out <- c(out, "### DEXSeq",
             "```{ojs}",
             "dex = xpTable('dexseq')",
             "Inputs.table(dex, {rows: 15})",
             "```")
  }
  # ---- variable transcripts
  if (!is.null(iso$transcript_pca)) {
    out <- c(out, "### Variable transcripts",
             "```{ojs}",
             "txPca = xpTable('tx_pca')",
             "Inputs.table(txPca, {rows: 15})",
             "```")
  }
  # ---- static isoform plots (embedded base64 gallery, same as the htmltools engine)
  if (length(iso$plots)) {
    n_all <- length(iso$plots)
    n_show <- min(n_all, as.integer(st$max_isoform_plots %||% 6L))
    uris <- lapply(iso$plots[seq_len(n_show)],
                   function(p) .xp_img_uri(p, figures_dir = NULL,
                                           max_px = st$static_max_px %||% 1000L))
    caps <- sub("\\.[^.]+$", "", basename(iso$plots[seq_len(n_show)]))
    items <- Filter(function(x) !is.null(x$uri),
                    lapply(seq_along(uris), function(i)
                      list(uri = uris[[i]], cap = caps[[i]])))
    if (length(items)) {
      out <- c(out, "### Plots (static)")
      out <- c(out, "```{ojs}")
      out <- c(out, paste0("isoFigs = [",
                           paste(vapply(items, function(it)
                             sprintf("{uri: '%s', cap: '%s'}",
                                     as.character(it$uri), .xq_js_str(it$cap)),
                             character(1)), collapse = ", "), "]"))
      out <- c(out, paste(
        "htl.html`<div class='xq-static-gallery' style='display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:10px'>",
        "${isoFigs.map(f => htl.html`<figure style='margin:0'>",
        "<img src='${f.uri}' loading='lazy' style='width:100%;cursor:zoom-in' title='${f.cap}'",
        " onclick=${(ev) => { const img = ev.target; img.style.position = img.style.position === 'fixed' ? '' : 'fixed'; img.style.inset = img.style.position === 'fixed' ? '2%' : ''; img.style.width = img.style.position === 'fixed' ? '96%' : '100%'; img.style.zIndex = 9999; img.style.background = '#fff'; }}>",
        "<figcaption style='font-size:0.78em;color:#52606d;word-break:break-all'>${f.cap}</figcaption>",
        "</figure>`)}",
        "</div>`", sep = "\n"))
      out <- c(out, "```")
    }
  }
  # ---- methods
  if (length(iso$methods)) {
    out <- c(out, "### Methods",
             "```{ojs}")
    out <- c(out, paste0("isoMethods = [",
                         paste(vapply(iso$methods, function(m)
                           paste0("'", .xq_js_str(m), "'"), character(1)),
                           collapse = ", "), "]"))
    out <- c(out, paste(
      "htl.html`<div class='xq-note'>Statistical notes recorded by the isoform analysis:</div>",
      "<ul>${isoMethods.map(m => htl.html`<li style='font-size:0.9em'>${m}</li>`)}</ul>`",
      sep = "\n"))
    out <- c(out, "```")
  }
  out
}

.xq_page_figures <- function(res, st) {
  if (!nrow(res$figures) || !isTRUE(st$static_figures)) return(NULL)
  out <- c("## Figures {.tabset}")
  grp <- unique(res$figures$group)
  for (g in grp) {
    title <- if (identical(g, ".")) "(root)" else g
    items <- .xq_gallery_uris(res, st, g, st$max_static_figures)
    out <- c(out, sprintf("### %s", title))
    if (length(items)) {
      out <- c(out, "```{ojs}")
      out <- c(out, paste0("figs_", gsub("[^A-Za-z0-9_]", "_", g), " = [",
                           paste(vapply(items, function(it)
                             sprintf("{uri: '%s', cap: '%s'}",
                                     as.character(it$uri), .xq_js_str(it$cap)),
                             character(1)), collapse = ", "), "]"))
      out <- c(out, paste(
        "htl.html`<div style='display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:10px'>",
        "${figs_" %+% gsub("[^A-Za-z0-9_]", "_", g) %+% ".map(f => htl.html`<figure style='margin:0'>",
        "<img src='${f.uri}' loading='lazy' style='width:100%;border:1px solid #d9e2ec;border-radius:6px;cursor:zoom-in'",
        " onclick=${(ev) => { const img = ev.target; img.style.position = img.style.position === 'fixed' ? '' : 'fixed'; img.style.inset = img.style.position === 'fixed' ? '2%' : ''; img.style.width = img.style.position === 'fixed' ? '96%' : '100%'; img.style.zIndex = 9999; img.style.background = '#fff'; }}>",
        "<figcaption style='font-size:0.75em;color:#52606d;word-break:break-all'>${f.cap}</figcaption></figure>`)}</div>`",
        sep = "\n"))
      out <- c(out, "```")
    }
    # file listing with sizes
    fig <- res$figures[res$figures$group == g, , drop = FALSE]
    out <- c(out, "```{ojs}")
    out <- c(out, paste0("figList_", gsub("[^A-Za-z0-9_]", "_", g), " = [",
                         paste(vapply(seq_len(nrow(fig)), function(i)
                           sprintf("{file: '%s', type: '%s', size_kb: %s}",
                                   .xq_js_str(fig$name[[i]]), fig$type[[i]],
                                   round(fig$size[[i]] / 1024, 1)),
                           character(1)), collapse = ", "), "]"))
    out <- c(out, paste(
      "Inputs.table(figList_" %+% gsub("[^A-Za-z0-9_]", "_", g) %+% ", {rows: 10})",
      sep = "\n"))
    out <- c(out, "```")
  }
  out
}

.xq_page_runinfo <- function(res, st) {
  out <- c("## Run info {.tabset}")
  # Parameters
  params <- .xp_params_table(res)
  out <- c(out, "### Parameters",
           "```{ojs}")
  out <- c(out, paste0("params = [",
                       paste(vapply(seq_len(nrow(params)), function(i)
                         sprintf("{Parameter: '%s', Value: '%s'}",
                                 .xq_js_str(params$Parameter[[i]]),
                                 .xq_js_str(params$Value[[i]])), character(1)),
                         collapse = ", "), "]"))
  out <- c(out, "Inputs.table(params, {rows: 15})", "```")

  # Logs
  lg <- res$logs
  out <- c(out, "### Logs", "```{ojs}")
  blocks <- c()
  if (!is.null(lg$warnings)) {
    blocks <- c(blocks, paste0(
      "{title: 'Warnings recorded during the run', text: '",
      .xq_js_str(paste(utils::head(lg$warnings, 400L), collapse = "\n")), "'}"))
  }
  if (!is.null(lg$session_info)) {
    blocks <- c(blocks, paste0(
      "{title: 'Session information', text: '",
      .xq_js_str(paste(utils::head(lg$session_info, 400L), collapse = "\n")), "'}"))
  }
  if (length(blocks)) {
    out <- c(out, paste0("logBlocks = [", paste(blocks, collapse = ",\n"), "]"))
    out <- c(out, paste(
      "htl.html`<div>${logBlocks.map(b => htl.html`<details style='margin:8px 0'>",
      "<summary style='cursor:pointer;font-weight:600'>${b.title}</summary>",
      "<pre style='max-height:420px;overflow:auto;font-size:0.8em;background:#f5f8fb;padding:8px;border-radius:6px'>${b.text}</pre>",
      "</details>`)}</div>`", sep = "\n"))
  } else {
    out <- c(out, "htl.html`<div class='xq-note'>No log files were found in log/.</div>`")
  }
  out <- c(out, "```")

  # Files inventory
  fi <- res$files
  out <- c(out, "### Files", "```{ojs}")
  out <- c(out, paste0("files = [",
                       paste(vapply(seq_len(nrow(fi)), function(i)
                         sprintf("{path: '%s', category: '%s', size_kb: %s, modified: '%s'}",
                                 .xq_js_str(fi$path[[i]]), .xq_js_str(fi$category[[i]]),
                                 fi$size_kb[[i]], .xq_js_str(fi$modified[[i]])), character(1)),
                         collapse = ", "), "]"))
  out <- c(out, "Inputs.table(files, {rows: 25, sort: 'path'})", "```")
  out
}

# ---- the report .qmd body --------------------------------------------------------

#' Assemble the complete Quarto dashboard source
#'
#' @param res An `expressom_results` object.
#' @param st Report settings (`expressom_report_settings()`).
#' @return A character vector of .qmd lines.
#' @keywords internal
.xq_dashboard_source <- function(res, st) {
  header <- c(
    "---",
    "title: \"ExpressOM report\"",   # filled by the caller with comparison
    "format:",
    "  dashboard:",
    "    nav-buttons:",
    "      - icon: house",
    "      - icon: activity",
    "    scrolling: true",
    "    theme: cosmo",
    "embed-resources: true",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---"
  )
  body <- c()
  pages <- list(
    .xq_page_overview(res, st),
    .xq_page_qc(res, st),
    .xq_page_dge(res, st),
    .xq_page_functional(res, st),
    .xq_page_isoform(res, st),
    .xq_page_figures(res, st),
    .xq_page_runinfo(res, st)
  )
  for (p in pages) {
    if (!is.null(p) && length(p)) body <- c(body, "", p)
  }
  # The {ojs} runtime prelude (xpTable, DuckDB helpers, formatters) opens the
  # body so that every later cell can use it. Observable's parser accepts one
  # value definition per cell, so each prelude entry becomes its own cell.
  prelude <- unlist(lapply(.xq_ojs_prelude_cells(), function(cell_src) {
    c("```{ojs}", cell_src, "```")
  }))
  c(header, "", prelude, "", body)
}

# ---- rendering pipeline ----------------------------------------------------------

#' Patch the rendered dashboard so it opens from file://
#'
#' Quarto refuses to run {ojs} cells under the file: protocol because the
#' observable runtime normally needs dynamic module imports from a CDN. With
#' `embed-resources: true` every script is already inlined, and the embedded
#' tables travel inside the document, so the guard is overly strict here:
#' ExpressOM reports are opened by double-click, like the htmltools engine's
#' output. The two guards are replaced by a no-op; DuckDB-Wasm still needs the
#' network on first use (exactly like plotly.js is fetched from cdn.plot.ly
#' today) and every explorer has an offline JS fallback.
#'
#' @param path Path to the rendered HTML (patched in place).
#' @keywords internal
.xq_patch_file_protocol <- function(path) {
  html <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n")
  guard1 <- paste(
    "if (window.location.protocol === \"file:\") {",
    "    displayFileProtocolError();",
    "    return;",
    "  }",
    sep = "\n"
  )
  guard2 <- paste0(
    "if (window.location.protocol === \"file:\") { ",
    "alert(\"The OJS runtime does not work with file:// URLs. ",
    "Please use a web server to view this document.\"); }"
  )
  repl1 <- paste(
    "/* ExpressOM: file:// enabled - every resource is embedded inline,",
    "   and each explorer falls back to plain JS when a CDN fetch fails. */",
    "if (false) {",
    "    displayFileProtocolError();",
    "    return;",
    "  }",
    sep = "\n"
  )
  repl2 <- "/* ExpressOM: embedded report, file:// allowed */"
  n1 <- length(gregexpr(guard1, html, fixed = TRUE)[[1L]])
  n2 <- length(gregexpr(guard2, html, fixed = TRUE)[[1L]])
  if (n1 == 1L) html <- sub(guard1, repl1, html, fixed = TRUE)
  if (n2 == 1L) html <- sub(guard2, repl2, html, fixed = TRUE)
  writeLines(html, path, useBytes = TRUE)
  invisible(path)
}

#' Locate the quarto CLI
#' @keywords internal
.xq_quarto_binary <- function() {
  env <- Sys.getenv("QUARTO_PATH", "")
  if (nzchar(env) && file.exists(env)) return(normalizePath(env))
  hit <- Sys.which("quarto")
  if (nzchar(hit)) return(normalizePath(hit))
  stop(
    "The 'quarto' command was not found. Install Quarto >= 1.4 ",
    "(https://quarto.org) or set QUARTO_PATH to the quarto binary.",
    call. = FALSE
  )
}

#' Render the Quarto dashboard report
#'
#' Called by [build_expressom_report()] when `engine = "quarto"` (see
#' mod_report.R). Writes a temporary .qmd, renders it with the quarto CLI,
#' patches the file:// guards and moves the result to `output_file`.
#'
#' @param res An `expressom_results` object.
#' @param st Report settings.
#' @param output_file Destination HTML path.
#' @param quiet Suppress quarto output.
#' @return `output_file` (invisibly).
#' @keywords internal
.xq_render_dashboard <- function(res, st, output_file, quiet = TRUE) {
  .rep_require("jsonlite", "to build the Quarto dashboard report")
  quarto_bin <- .xq_quarto_binary()

  work <- tempfile("expressom_quarto_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  report_dir <- dirname(output_file)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- title --------------------------------------------------------------
  title <- st$title %||% "ExpressOM report"
  if (!is.na(res$comparison)) {
    title <- paste0(title, " \\u00b7 ", res$comparison)
  }

  # ---- collect + embed the data -------------------------------------------
  tables <- .xq_collect_tables(res, st)
  data_block <- .xq_data_block(res, st, tables)

  # ---- assemble the .qmd ---------------------------------------------------
  src <- .xq_dashboard_source(res, st)
  src[2] <- paste0("title: \"", .xq_js_str(title), "\"")

  qmd <- file.path(work, "report.qmd")
  writeLines(src, qmd, useBytes = TRUE)

  # the data payload rides along as an include so it lands in <head>
  include_file <- file.path(work, "xpdata.html")
  writeLines(c(.xq_css(), data_block), include_file, useBytes = TRUE)

  args <- c(
    "render", qmd,
    "--to", "dashboard",
    "--output-dir", work,
    "--metadata", paste0("include-in-header: ", include_file),
    if (quiet) "--quiet" else NULL
  )
  res_proc <- suppressWarnings(system2(
    quarto_bin, args, stdout = if (quiet) FALSE else "", stderr = ""
  ))
  out_html <- file.path(work, "report.html")
  if (res_proc != 0L || !file.exists(out_html)) {
    stop(
      "quarto render failed (exit status ", res_proc, "). ",
      "Run with quiet = FALSE for the full log.",
      call. = FALSE
    )
  }

  # ---- move into place + patch the file:// guards --------------------------
  file.copy(out_html, output_file, overwrite = TRUE)
  .xq_patch_file_protocol(output_file)

  invisible(normalizePath(output_file, winslash = "/"))
}
