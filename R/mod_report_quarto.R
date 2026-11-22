# mod_report_quarto.R -- the Quarto dashboard report engine
# (`build_expressom_report(engine = "quarto")`).
#
# Architecture, in one paragraph: R computes everything (see
# mod_report_quarto_collect.R) and embeds it as compressed binary columns and
# small JSON specs in <head>. The page body is a Quarto *dashboard*: pages are
# the navbar entries, tabs are cards, and every card is "a How-to-read note
# (static HTML) + one Observable JS cell that mounts a view" (inst/xq/*.js).
# Views build lazily, only when their tab is first shown, so the file opens
# fast no matter how many tables it contains. The JS libraries (Observable
# Plot / Inputs, regl-scatterplot for WebGL, fflate) are vendored inline: the
# report needs no network, except for DuckDB-Wasm on first free-SQL use.

# ---- assets ----------------------------------------------------------------------------

#' Locate a runtime asset (inst/xq/<name>)
#' @keywords internal
.xq_asset_path <- function(name) {
  cands <- c(
    file.path(getOption("expressom.xq_dir", ""), name),
    system.file("xq", name, package = "ExpressOM"),
    file.path("inst", "xq", name),
    file.path("..", "..", "inst", "xq", name)
  )
  cands <- cands[nzchar(cands) & file.exists(cands)]
  if (!length(cands)) {
    stop("ExpressOM runtime asset not found: xq/", name,
         " (reinstall the package, or set options(expressom.xq_dir = <inst/xq>)).",
         call. = FALSE)
  }
  cands[[1L]]
}

.xq_asset <- function(name) {
  path <- .xq_asset_path(name)
  txt <- readChar(path, file.size(path), useBytes = TRUE)
  Encoding(txt) <- "UTF-8"
  txt
}

# Per-card "How to read this" notes, filled while the qmd is written and
# embedded in <head> as <script type="text/html" id="xqh-<view id>">.
.xq_reg <- new.env(parent = emptyenv())
.xq_reg$help <- list()

# a script body must never contain a sequence that ends / confuses the element
.xq_safe_js <- function(js) {
  js <- gsub("</script", "<\\/script", js, fixed = TRUE)
  gsub("<!--", "<\\!--", js, fixed = TRUE)
}

#' The <head> include: styles, data payloads, runtime scripts
#' @keywords internal
.xq_head_html <- function(col, st) {
  blocks <- character(0)
  blocks <- c(blocks, paste0("<style id=\"xq-css\">\n", .xq_asset("xq.css"), "\n</style>"))
  blocks <- c(blocks, .xq_json_script("cfg", col$cfg))
  for (id in names(col$views)) {
    blocks <- c(blocks, .xq_json_script(paste0("v_", id), col$views[[id]]))
  }
  for (id in names(.xq_reg$help)) {
    blocks <- c(blocks, paste0("<script type=\"text/html\" id=\"xqh-", id, "\">",
                               .xq_safe_js(.xq_reg$help[[id]]), "</script>"))
  }
  keys <- names(col$tables)
  enc <- lapply(keys, function(k) {
    t <- col$tables[[k]]
    .xq_encode_table(t$df, max_rows = st$table_rows %||% Inf, no_p = t$no_p)
  })
  names(enc) <- keys
  for (k in keys) {
    if (!is.null(enc[[k]])) blocks <- c(blocks, .xq_table_script(.xq_key(k), enc[[k]]))
  }
  for (f in c("xq-vendor.min.js", "xq-core.js", "xq-widgets.js", "xq-scatter.js",
              "xq-charts.js", "xq-views.js")) {
    blocks <- c(blocks, paste0("<script>\n", .xq_safe_js(.xq_asset(f)), "\n</script>"))
  }
  attr(blocks, "bytes") <- sum(vapply(enc, function(e) if (is.null(e)) 0 else e$bytes, numeric(1)))
  blocks
}

# ---- small HTML / qmd emitters ----------------------------------------------------------------

.xq_raw <- function(html) c("```{=html}", html, "```", "")
.xq_ojs <- function(id) c("```{ojs}", sprintf("XQ.mount(\"%s\", invalidation)", id), "```", "")

#' One dashboard card = heading + how-to-read note + optional static HTML + view
#' @keywords internal
.xq_card <- function(title, help = NULL, view = NULL, html = NULL, extra_views = NULL) {
  # `::: {.card}` keeps everything of a tab in ONE card, stacked vertically
  # (a dashboard otherwise lays the blocks of a tab out side by side).
  # Quarto may also reorder raw HTML relative to code cells, so the note is
  # never a separate block: with a view it is registered and rendered by the
  # view itself (XQ.mount puts it on top); without one it leads the raw HTML.
  if (!is.null(view) && !is.null(help)) .xq_reg$help[[view]] <- help
  lead <- if (is.null(view)) help else NULL
  c(
    paste0("### ", .xp_clean_title(title)), "",
    "::: {.card}", "",
    if (!is.null(lead) || !is.null(html)) .xq_raw(paste0(lead %||% "", html %||% "")),
    if (!is.null(view)) .xq_ojs(view),
    unlist(lapply(extra_views, .xq_ojs)),
    ":::", ""
  )
}

.xq_page <- function(title, cards, tabset = TRUE) {
  cards <- cards[!vapply(cards, is.null, logical(1))]
  if (!length(cards)) return(NULL)
  c(paste0("# ", .xp_clean_title(title)), "",
    if (tabset) "## Row {.tabset}" else "## Row", "",
    unlist(cards))
}

.xq_tile <- function(label, value, sub = NULL, cls = NULL) {
  paste0(
    "<div class=\"xq-stat", if (!is.null(cls)) paste0(" ", cls) else "", "\">",
    "<div class=\"xq-stat-v\">", .xq_e(value), "</div>",
    "<div class=\"xq-stat-l\">", .xq_e(label), "</div>",
    if (!is.null(sub)) paste0("<div class=\"xq-stat-s\">", .xq_e(sub), "</div>") else "",
    "</div>"
  )
}

.xq_rel_link <- function(path, root, text = NULL) {
  rel <- .rep_rel(path, root)
  paste0("<a href=\"", .xq_e(rel), "\" target=\"_blank\" rel=\"noopener\">", .xq_e(text %||% rel), "</a>")
}

# NB: a <dl>, not a <table>: Quarto's dashboard script attaches a sticky header
# to every <table> and throws on one without <thead>, which blanks the page.
.xq_kv_html <- function(df) {
  paste0("<dl class=\"xq-kv\">",
         paste0("<dt>", .xq_e(df[[1L]]), "</dt><dd>", .xq_e(df[[2L]]), "</dd>", collapse = ""),
         "</dl>")
}

# ---- static pieces ---------------------------------------------------------------------------------

.xq_overview_html <- function(res, st, ctx) {
  s <- .xp_dge_summary(res, st)
  cfg <- res$config
  tiles <- character(0)
  add <- function(...) tiles <<- c(tiles, .xq_tile(...))
  if (!is.null(s)) {
    add("Comparison", paste(cfg$level %||% "?", "vs", cfg$base %||% "?"), paste(cfg$test, "test"))
    add("Genes tested", .xp_fmt(s$n_tested), "after independent filtering")
    if (s$is_lrt) {
      add("Significant genes", .xp_fmt(s$n_sig), paste0("LRT, padj < ", s$padj), "xq-up")
    } else {
      add("Upregulated", .xp_fmt(s$n_up), paste0("padj < ", s$padj, ", log2FC > ", s$lfc), "xq-up")
      add("Downregulated", .xp_fmt(s$n_down), paste0("padj < ", s$padj, ", log2FC < -", s$lfc), "xq-down")
    }
  }
  if (length(res$pca)) add("Samples", .xp_fmt(nrow(res$pca[[1L]]$scores)), "in the PCA")
  inv <- res$functional$inventory
  if (nrow(inv)) add("Enrichment tables", .xp_fmt(nrow(inv)), paste(unique(inv$source), collapse = " / "))
  if (!is.null(res$decouple)) {
    add("decoupleR resources", .xp_fmt(length(res$decouple)),
        paste(unique(vapply(res$decouple, `[[`, character(1), "resource")), collapse = " / "))
  }
  iso <- res$isoform
  if (!is.null(iso$dte)) {
    p <- .rep_num(iso$dte$padj)
    add("DTE transcripts", .xp_fmt(sum(!is.na(p) & p < ctx$padj)), "significant")
  }
  if (!is.null(iso$dtu)) {
    pc <- .rep_first_col(iso$dtu, c("adj_pvalue", "padj", "FDR"))
    if (!is.na(pc)) {
      p <- .rep_num(iso$dtu[[pc]])
      add("DTU features", .xp_fmt(sum(!is.na(p) & p < ctx$padj)), "significant")
    }
  }
  paste0("<div class=\"xq-stats\">", paste(tiles, collapse = ""), "</div>")
}

.xq_contents_html <- function(res) {
  inv <- res$functional$inventory
  iso <- res$isoform
  li <- function(ok, label, detail = NULL) {
    paste0("<li class=\"", if (ok) "xq-ok" else "xq-missing", "\"><span class=\"xq-tick\">",
           if (ok) "✓" else "–", "</span>", .xq_e(label),
           if (!is.null(detail)) paste0("<span class=\"xq-detail\">", .xq_e(detail), "</span>") else "", "</li>")
  }
  paste0("<ul class=\"xq-avail\">",
    li(!is.null(res$dge), "Differential expression", if (!is.null(res$dge)) paste(.xp_fmt(nrow(res$dge$results)), "genes")),
    li(length(res$pca) > 0L, "PCA / sample QC", if (length(res$pca)) paste(length(res$pca), "plot(s)")),
    li(nrow(inv) > 0L, "Functional enrichment (ORA / GSEA / FGSEA / SPIA / TF)", if (nrow(inv)) paste(nrow(inv), "table(s)")),
    li(!is.null(res$decouple), "decoupleR TF / pathway activity", if (!is.null(res$decouple)) paste(length(res$decouple), "resource(s)")),
    li(!is.null(res$decouple_regulons), "decoupleR TF regulon explorer", if (!is.null(res$decouple_regulons)) paste(length(res$decouple_regulons$by_resource), "resource(s)")),
    li(!is.null(iso), "Isoform analysis (DTE / DTU / switches)"),
    li(nrow(res$figures) > 0L, "Static figures", if (nrow(res$figures)) paste(nrow(res$figures), "file(s)")),
    "</ul>"
  )
}

.xq_gallery_items <- function(paths, st, cap_fun = function(p) sub("\\.[^.]+$", "", basename(p))) {
  if (!length(paths)) return(list())
  conv <- function(p) {
    if (grepl("\\.pdf$", p, ignore.case = TRUE)) .xp_pdf_uri(p, dpi = st$static_dpi %||% 60)
    else .xp_img_uri(p, max_px = st$static_max_px %||% 1000L)
  }
  uris <- if (.Platform$OS.type != "windows" && length(paths) > 1L) {
    parallel::mclapply(paths, conv, mc.cores = max(1L, min(parallel::detectCores() - 1L, 8L, length(paths))))
  } else lapply(paths, conv)
  items <- lapply(seq_along(paths), function(i) {
    u <- uris[[i]]
    if (is.null(u) || inherits(u, "try-error") || !is.character(u)) return(NULL)
    list(uri = u, cap = cap_fun(paths[[i]]))
  })
  items[!vapply(items, is.null, logical(1))]
}

.xq_gallery_html <- function(items) {
  if (!length(items)) return("<p class=\"xq-intro\">No previews available.</p>")
  paste0("<div class=\"xq-gallery\">", paste0(vapply(items, function(it) {
    paste0("<figure><img class=\"xq-zoom\" loading=\"lazy\" decoding=\"async\" src=\"", it$uri,
           "\" alt=\"", .xq_e(it$cap), "\" title=\"Click to enlarge\"><figcaption>", .xq_e(it$cap),
           "</figcaption></figure>")
  }, character(1)), collapse = ""), "</div>")
}

.xq_logs_html <- function(res) {
  out <- character(0)
  dl <- res$logs$dge_params
  if (!is.null(dl)) {
    out <- c(out, "<h4 class=\"xq-h4\">Differential expression parameters</h4>",
             paste0("<pre class=\"xq-pre\">", .xq_e(paste(utils::capture.output(print(unlist(dl))), collapse = "\n")), "</pre>"))
  }
  det <- function(title, lines) {
    if (is.null(lines) || !length(lines)) return(NULL)
    paste0("<details class=\"xq-sqlbox\"><summary>", .xq_e(title), "</summary><pre class=\"xq-pre\">",
           .xq_e(paste(utils::head(lines, 400L), collapse = "\n")), "</pre></details>")
  }
  out <- c(out, det("Warnings recorded during the run", res$logs$warnings), det("Session information", res$logs$session_info))
  links <- function(title, paths) {
    if (!length(paths)) return(NULL)
    c(paste0("<h4 class=\"xq-h4\">", .xq_e(title), "</h4>"),
      paste0("<ul class=\"xq-links\">", paste0("<li>", vapply(paths, function(p) .xq_rel_link(p, res$out_dir), character(1)), "</li>", collapse = ""), "</ul>"))
  }
  out <- c(out, links("Isoform logs", res$logs$isoform_logs), links("External-predictor (WSL) logs", res$logs$wsl_logs),
           links("OmnipathR logs", res$logs$omnipath_logs), links("targets pipeline logs", res$logs$targets_logs))
  if (!length(out)) return("<p class=\"xq-intro\">No log files were found in log/.</p>")
  paste(out, collapse = "\n")
}

# ---- pages ---------------------------------------------------------------------------------------------

.xq_page_overview <- function(res, st, col) {
  ctx <- col$ctx
  glance <- c(.xq_overview_html(res, st, ctx))
  row1 <- c("# Overview", "", "## Row", "",
            .xq_card("At a glance", help = .xq_help_html("overview", ctx), html = glance))
  tabs <- list(
    if (!is.null(col$views$deg_counts)) .xq_card("DEGs by cut-off", help = .xq_help_html("deg_counts", ctx),
                                                  view = "deg_counts", extra_views = "deg_counts_table"),
    .xq_card("Contents", help = .xq_help_html("contents", ctx), html = .xq_contents_html(res)),
    if (length(res$notes)) .xq_card("Notes", html = paste0("<ul>", paste0("<li>", .xq_e(res$notes), "</li>", collapse = ""), "</ul>"))
  )
  tabs <- tabs[!vapply(tabs, is.null, logical(1))]
  c(row1, "## Row {.tabset}", "", unlist(tabs))
}

.xq_page_qc <- function(res, st, col) {
  ctx <- col$ctx
  cards <- list()
  for (id in col$pca) {
    nm <- col$views[[id]]$filename
    cards[[length(cards) + 1L]] <- .xq_card(col$views[[id]]$title, help = .xq_help_html("pca", ctx, name = nm), view = id)
  }
  if (!is.null(col$views$similarity)) cards[[length(cards) + 1L]] <- .xq_card("Sample similarity", help = .xq_help_html("similarity", ctx), view = "similarity")
  if (!is.null(col$views$expr_qc)) cards[[length(cards) + 1L]] <- .xq_card("Expression distributions", help = .xq_help_html("expr_qc", ctx), view = "expr_qc")
  .xq_page("Samples & QC", cards)
}

.xq_page_dge <- function(res, st, col) {
  if (is.null(res$dge)) return(NULL)
  ctx <- col$ctx
  v <- col$views
  files <- res$dge$sig_lists
  lists_html <- if (length(files)) {
    paste0("<p class=\"xq-intro\">Files in <code>DE_raw_results/</code>:</p><ul class=\"xq-links\">",
           paste0("<li><a href=\"DE_raw_results/", .xq_e(files), "\" target=\"_blank\" rel=\"noopener\">", .xq_e(files), "</a></li>", collapse = ""), "</ul>")
  } else "<p class=\"xq-intro\">No per-cut-off DEG lists were found.</p>"
  cards <- list(
    .xq_card("Volcano / MA", help = .xq_help_html("volcano", ctx), view = "volcano"),
    .xq_card("Diagnostics", help = .xq_help_html("diagnostics", ctx), view = "diag"),
    if (!is.null(v$heat)) .xq_card("Top genes heatmap", help = .xq_help_html("heatmap", ctx), view = "heat"),
    .xq_card("Results table", help = .xq_help_html("de_table", ctx), view = "de_table"),
    .xq_card("DEG lists", help = .xq_help_html("gene_lists", ctx), html = lists_html),
    if (!is.null(v$sql)) .xq_card("SQL explorer", help = .xq_help_html("sql_de", ctx), view = "sql")
  )
  .xq_page("Differential expression", cards)
}

.xq_page_functional <- function(res, st, col) {
  ctx <- col$ctx
  cards <- list()
  multi_src <- length(unique(vapply(col$enrich, `[[`, character(1), "source"))) > 1L
  src_titles <- c(DGE = "Gene level", DTE = "DTE", DTU = "DTU")
  for (e in col$enrich) {
    title <- if (multi_src) paste0(src_titles[[e$source]] %||% e$source, " · ", e$method) else e$method
    cards[[length(cards) + 1L]] <- .xq_card(title, help = .xq_help_html("enrichment", ctx, kind = e$kind, method = e$method), view = e$id)
  }
  for (d in col$decouple) {
    cards[[length(cards) + 1L]] <- .xq_card(paste0("decoupleR · ", toupper(d$resource)),
                                            help = .xq_help_html("decouple", ctx, resource = d$resource), view = d$id)
  }
  for (r in col$regulon) {
    cards[[length(cards) + 1L]] <- .xq_card(paste0("TF regulons · ", toupper(r$resource)),
                                            help = .xq_help_html("regulon", ctx), view = r$id)
  }
  if (isTRUE(col$reg_summary)) {
    cards[[length(cards) + 1L]] <- .xq_card("TF regulon summary", help = .xq_help_html("regulon", ctx), view = "reg_summary")
  }
  if (isTRUE(col$dec_status)) {
    cards[[length(cards) + 1L]] <- .xq_card("decoupleR network status",
                                            help = .xq_help_html("decouple", ctx), view = "dec_status")
  }
  .xq_page("Functional analysis", cards)
}

.xq_page_isoform <- function(res, st, col) {
  if (is.null(res$isoform)) return(NULL)
  ctx <- col$ctx
  v <- col$views
  iso <- res$isoform
  cards <- list(
    if (!is.null(v$sw)) .xq_card("Switch overview", help = .xq_help_html("switch", ctx), view = "sw"),
    if (!is.null(v$dte)) .xq_card("DTE", help = .xq_help_html("dte", ctx), view = "dte"),
    if (!is.null(v$dte_table)) .xq_card("DTE table", help = .xq_help_html("de_table", ctx), view = "dte_table"),
    if (!is.null(v$dtu)) .xq_card("DTU", help = .xq_help_html("dtu", ctx), view = "dtu"),
    if (!is.null(v$dtu_table)) .xq_card("DTU table", help = .xq_help_html("de_table", ctx), view = "dtu_table"),
    if (!is.null(v$dex)) .xq_card("DEXSeq", help = .xq_help_html("dexseq", ctx), view = "dex"),
    if (!is.null(v$txpca)) .xq_card("Variable transcripts", help = .xq_help_html("tx_pca", ctx), view = "txpca")
  )
  if (length(iso$plots)) {
    n_show <- as.integer(st$max_isoform_plots %||% 24L)
    items <- .xq_gallery_items(utils::head(iso$plots, n_show), st)
    note <- if (length(iso$plots) > n_show) paste0("<p class=\"xq-intro\">Showing ", n_show, " of ", length(iso$plots), " plots.</p>") else ""
    cards[[length(cards) + 1L]] <- .xq_card("Plots", help = .xq_help_html("iso_plots", ctx), html = paste0(note, .xq_gallery_html(items)))
  }
  if (length(iso$methods)) {
    cards[[length(cards) + 1L]] <- .xq_card("Methods", help = .xq_help_html("iso_methods", ctx),
                                            html = paste0("<ul>", paste0("<li>", .xq_e(iso$methods), "</li>", collapse = ""), "</ul>"))
  }
  .xq_page("Isoforms", cards)
}

.xq_page_figures <- function(res, st, col) {
  if (!nrow(res$figures) || !isTRUE(st$static_figures)) return(NULL)
  ctx <- col$ctx
  groups <- unique(res$figures$group)
  have_pdf <- requireNamespace("pdftools", quietly = TRUE)
  cards <- lapply(groups, function(g) {
    fig <- res$figures[res$figures$group == g, , drop = FALSE]
    pngs <- fig[fig$type == "png", , drop = FALSE]
    pdfs <- fig[fig$type == "pdf" & fig$size >= 2000, , drop = FALSE]
    room <- as.integer(st$max_static_figures %||% 24L)
    take_png <- utils::head(pngs$path, room)
    take_pdf <- if (have_pdf) utils::head(pdfs$path, max(0L, room - length(take_png))) else character(0)
    items <- .xq_gallery_items(c(take_png, take_pdf), st)
    shown <- length(items)
    note <- paste0("<p class=\"xq-intro\">", nrow(fig), " file(s) in this folder; ", shown, " preview(s) shown.",
                   if (nrow(pdfs) && !have_pdf) " Install the 'pdftools' package to preview PDF figures." else "", "</p>")
    rest <- fig$path[!(fig$path %in% c(take_png, take_pdf))]
    links <- if (length(rest)) paste0("<details class=\"xq-sqlbox\"><summary>Other files (", length(rest), ")</summary><ul class=\"xq-links\">",
                                       paste0("<li>", vapply(utils::head(rest, 300L), function(p) .xq_rel_link(p, res$out_dir), character(1)), "</li>", collapse = ""),
                                       "</ul></details>") else ""
    .xq_card(if (identical(g, ".")) "(root)" else g, help = .xq_help_html("figures", ctx), html = paste0(note, .xq_gallery_html(items), links))
  })
  .xq_page("Figures", cards)
}

.xq_page_runinfo <- function(res, st, col) {
  ctx <- col$ctx
  cards <- list(
    .xq_card("Parameters", help = .xq_help_html("params", ctx), html = .xq_kv_html(.xp_params_table(res))),
    .xq_card("Logs", help = .xq_help_html("logs", ctx), html = .xq_logs_html(res)),
    if (!is.null(col$views$files)) .xq_card("Files", help = .xq_help_html("files", ctx), view = "files")
  )
  .xq_page("Run info", cards)
}

# ---- the document ------------------------------------------------------------------------------------------

.xq_yaml_str <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  paste0("\"", gsub("\"", "\\\"", x, fixed = TRUE), "\"")
}

#' Assemble the complete Quarto dashboard source
#' @keywords internal
.xq_dashboard_source <- function(res, st, col, title) {
  .xq_reg$help <- list()
  header <- c(
    "---",
    paste0("title: ", .xq_yaml_str(title)),
    "lang: en",
    "format:",
    "  dashboard:",
    "    theme: default",
    "    scrolling: true",
    "    expandable: false",
    "embed-resources: true",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---", ""
  )
  pages <- list(
    .xq_page_overview(res, st, col),
    .xq_page_qc(res, st, col),
    .xq_page_dge(res, st, col),
    .xq_page_functional(res, st, col),
    .xq_page_isoform(res, st, col),
    .xq_page_figures(res, st, col),
    .xq_page_runinfo(res, st, col)
  )
  body <- unlist(lapply(pages, function(p) if (length(p)) c("", p)))
  c(header, body)
}

# ---- rendering pipeline --------------------------------------------------------------------------------------

#' Patch the rendered dashboard so it opens from file://
#'
#' Quarto refuses to start the Observable runtime under the file: protocol
#' because it normally imports modules from a CDN. Here every script is already
#' inline, so the two guards are neutralised. Patterns are matched loosely
#' (quote style / spacing) so minor Quarto version changes do not break it.
#'
#' @param path Rendered HTML (patched in place).
#' @return The number of guards replaced (invisibly).
#' @keywords internal
.xq_patch_file_protocol <- function(path) {
  html <- rawToChar(readBin(path, "raw", file.size(path)))
  pat <- "window\\.location\\.protocol\\s*===?\\s*([\"'])file:\\1"
  m <- gregexpr(pat, html, perl = TRUE, useBytes = TRUE)[[1L]]
  n <- sum(m != -1L)
  if (n > 0L) {
    html <- gsub(pat, "false /* ExpressOM: embedded report, file:// allowed */", html,
                 perl = TRUE, useBytes = TRUE)
    writeBin(charToRaw(html), path)
  }
  invisible(n)
}

#' Inject the runtime (styles, data payloads, scripts) into the rendered page
#'
#' The payload is several MB of base64 and ~700 KB of JavaScript. It is added
#' to the finished HTML, just before `</head>`, rather than through Quarto's
#' `include-in-header`: that way pandoc and Quarto's Lua filters never scan it
#' (they warn on `<table` strings inside the vendored libraries) and the
#' render is much faster.
#'
#' @param path Rendered HTML (patched in place).
#' @param head_blocks Character vector from `.xq_head_html()`.
#' @return `TRUE` when the head was found (invisibly).
#' @keywords internal
.xq_inject_head <- function(path, head_blocks) {
  raw <- readBin(path, "raw", file.size(path))
  txt <- rawToChar(raw)
  pos <- regexpr("</head>", txt, fixed = TRUE, useBytes = TRUE)[[1L]]
  if (pos < 1L) stop("Could not find </head> in the rendered dashboard.", call. = FALSE)
  payload <- charToRaw(enc2utf8(paste0(paste(head_blocks, collapse = "\n"), "\n")))
  writeBin(c(raw[seq_len(pos - 1L)], payload, raw[pos:length(raw)]), path)
  invisible(TRUE)
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
#' Called by [build_expressom_report()] when `engine = "quarto"`.
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

  title <- st$title %||% "ExpressOM report"
  if (!is.na(res$comparison)) title <- paste0(title, " · ", res$comparison)

  col <- .xq_collect(res, st)
  src <- .xq_dashboard_source(res, st, col, title)   # registers the per-card notes
  head_blocks <- .xq_head_html(col, st)
  qmd <- file.path(work, "report.qmd")
  writeLines(enc2utf8(src), qmd, useBytes = TRUE)

  args <- c("render", qmd, "--to", "dashboard", "--output-dir", work, if (quiet) "--quiet" else NULL)
  status <- suppressWarnings(system2(quarto_bin, args, stdout = if (quiet) FALSE else "", stderr = ""))
  out_html <- file.path(work, "report.html")
  if (status != 0L || !file.exists(out_html)) {
    stop("quarto render failed (exit status ", status, "). Run with quiet = FALSE for the full log.", call. = FALSE)
  }

  .xq_inject_head(out_html, head_blocks)
  n_guard <- .xq_patch_file_protocol(out_html)
  if (n_guard < 1L) {
    warning("Could not find Quarto's file:// guard; the report may need a local web server to open.", call. = FALSE)
  }
  file.copy(out_html, output_file, overwrite = TRUE)
  invisible(normalizePath(output_file, winslash = "/"))
}
