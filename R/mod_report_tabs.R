# mod_report_tabs.R -- the tab tree shown by the unified HTML report and by
# the Shiny app ("same UI, precomputed backend").

#' Settings of the interactive report
#'
#' Defaults aim for a complete report: every figure, table and enrichment
#' plot is included. Lower any value to trade size/speed for content.
#'
#' @param title Report title.
#' @param max_points Points drawn in the volcano / MA plots.
#' @param heatmap_genes Genes in the expression heatmap.
#' @param top_terms Terms drawn in each enrichment plot.
#' @param table_rows Rows shown in each interactive table.
#' @param static_figures Embed the PDF figures of the run (needs `pdftools`).
#' @param max_static_figures Maximum number of embedded static figures per group.
#' @param static_dpi Resolution of the embedded static figures.
#' @param max_isoform_plots Maximum number of isoform-gallery PNGs.
#' @param plotly_enrichment_max Number of enrichment files that get an
#'   interactive plotly widget (the rest show the table only).
#' @param merge_enrichment_threshold Merge the tables of a method into a single
#'   interactive table when it has more than this many collections.
#' @param plotly_pca_max Number of PCA entries rendered as widgets (the extra
#'   ones are shown as tables only).
#' @param highlight_genes Genes to mark on the volcano plot.
#' @param padj_cutoff,lfc_cutoff Override the run's cut-offs.
#' @param theme "auto", "light" or "dark".
#' @param self_contained `TRUE` (default) renders with `rmarkdown` + pandoc
#'   for a single self-contained HTML file (every plotly / DT asset and every
#'   static figure inlined). `FALSE` falls back to `htmltools` and produces a
#'   folder with `lib/` and `figures/` next to the HTML (legacy).
#' @return A named list.
#' @export
expressom_report_settings <- function(
  title = "ExpressOM report",
  max_points = 2000L,
  heatmap_genes = 20L,
  top_terms = 15L,
  table_rows = .Machine$integer.max,
  static_figures = TRUE,
  max_static_figures = .Machine$integer.max,
  static_dpi = 60,
  max_isoform_plots = .Machine$integer.max,
  plotly_enrichment_max = .Machine$integer.max,
  merge_enrichment_threshold = 3L,
  plotly_pca_max = .Machine$integer.max,
  highlight_genes = NULL,
  padj_cutoff = NULL,
  lfc_cutoff = NULL,
  theme = c("auto", "light", "dark"),
  self_contained = TRUE
) {
  theme <- match.arg(theme)
  list(
    title = title,
    max_points = max_points,
    heatmap_genes = heatmap_genes,
    top_terms = top_terms,
    table_rows = table_rows,
    static_figures = static_figures,
    max_static_figures = max_static_figures,
    static_dpi = static_dpi,
    max_isoform_plots = max_isoform_plots,
    plotly_enrichment_max = plotly_enrichment_max,
    merge_enrichment_threshold = merge_enrichment_threshold,
    plotly_pca_max = plotly_pca_max,
    highlight_genes = highlight_genes,
    padj_cutoff = padj_cutoff,
    lfc_cutoff = lfc_cutoff,
    theme = theme,
    self_contained = self_contained
  )
}

# ---- tree constructors -------------------------------------------------------

.xp_node <- function(title, ..., content = NULL) {
  kids <- list(...)
  kids <- kids[!vapply(kids, is.null, logical(1))]
  if (is.null(content) && length(kids) == 0L) return(NULL)
  list(title = title, children = kids, content = content)
}

.xp_is_leaf <- function(node) is.function(node$content)

# ---- html building blocks ----------------------------------------------------

.xp_note <- function(...) htmltools::tags$p(class = "xp-note", ...)

.xp_box <- function(label, value, sub = NULL) {
  tg <- htmltools::tags
  tg$div(
    class = "xp-box",
    tg$div(class = "xp-box-label", label),
    tg$div(class = "xp-box-value", value),
    if (!is.null(sub)) tg$div(class = "xp-box-sub", sub)
  )
}

.xp_fmt <- function(x) {
  if (length(x) == 0L || is.na(x)) return("n/a")
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

.xp_kv_table <- function(kv) {
  tg <- htmltools::tags
  tg$table(
    class = "xp-kv",
    tg$tbody(lapply(seq_len(nrow(kv)), function(i) {
      tg$tr(tg$th(kv[[1L]][[i]]), tg$td(kv[[2L]][[i]]))
    }))
  )
}

# Returns either a data-URI string (figures_dir = NULL) or a list with a `src`
# pointing to "figures/<basename>" (figures_dir set, external mode).
.xp_img_uri <- function(path, figures_dir = NULL) {
  if (!is.null(figures_dir)) {
    dest <- file.path(figures_dir, basename(path))
    if (!file.exists(dest)) file.copy(path, dest, overwrite = TRUE)
    return(list(src = paste0("figures/", basename(path))))
  }
  tryCatch(knitr::image_uri(path), error = function(e) NULL)
}

# Same contract as .xp_img_uri.
.xp_pdf_uri <- function(pdf, dpi = 100, figures_dir = NULL) {
  if (!requireNamespace("pdftools", quietly = TRUE)) return(NULL)
  stem <- tempfile("xp_pdf_")
  tmp <- sprintf("%s_%d.%s", stem, 1L, "png")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch({
    pdftools::pdf_convert(
      pdf, format = "png", pages = 1L, dpi = dpi,
      filenames = paste0(gsub("%", "%%", stem, fixed = TRUE), "_%d.%s"),
      verbose = FALSE
    )
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(tmp)) return(NULL)
  if (!is.null(figures_dir)) {
    dest <- file.path(figures_dir, basename(sub("\\.pdf$", ".png", pdf)))
    file.copy(tmp, dest, overwrite = TRUE)
    return(list(src = paste0("figures/", basename(dest))))
  }
  .xp_img_uri(tmp)
}

# `uris` may contain data-URI strings, list(src = ...), or NULL entries.
.xp_gallery <- function(uris, captions) {
  tg <- htmltools::tags
  figures <- lapply(seq_along(uris), function(i) {
    u <- uris[[i]]
    if (is.null(u)) return(NULL)
    src <- if (is.list(u) && !is.null(u$src)) u$src else u
    tg$figure(
      tg$img(
        src = src, alt = captions[[i]], loading = "lazy",
        title = "Click to enlarge",
        onclick = "this.classList.toggle('xp-zoom')"
      ),
      tg$figcaption(captions[[i]])
    )
  })
  figures <- figures[!vapply(figures, is.null, logical(1))]
  tg$div(class = "xp-gallery", figures)
}

.xp_two_col <- function(...) htmltools::tags$div(class = "xp-grid", ...)

.xp_cutoffs <- function(res, st) {
  list(
    padj = st$padj_cutoff %||% res$config$padj_cutoff %||% 0.05,
    lfc  = st$lfc_cutoff  %||% res$config$lfc_cutoff  %||% 1
  )
}

# ---- summaries ---------------------------------------------------------------

.xp_dge_summary <- function(res, st) {
  tbl <- res$dge$results
  if (is.null(tbl)) return(NULL)
  co <- .xp_cutoffs(res, st)
  padj <- .rep_num(tbl$padj)
  lfc  <- .rep_num(tbl$log2FoldChange)
  sig  <- !is.na(padj) & padj < co$padj
  is_lrt <- identical(res$config$test, "LRT")
  list(
    n_tested = sum(!is.na(padj)),
    n_sig    = sum(sig),
    n_up     = if (is_lrt) NA_integer_ else sum(sig & !is.na(lfc) & lfc >  co$lfc),
    n_down   = if (is_lrt) NA_integer_ else sum(sig & !is.na(lfc) & lfc < -co$lfc),
    is_lrt   = is_lrt,
    padj     = co$padj,
    lfc      = co$lfc
  )
}

.xp_pca_title <- function(stub) {
  s <- sub("^PCA_", "", stub)
  pretty <- function(x) gsub("_", " ", x)
  if (grepl("^BatchCorrected_", s)) {
    paste0("Batch-corrected (", pretty(sub("^BatchCorrected_", "", s)), ")")
  } else if (grepl("^ComparisonOnly_", s)) {
    paste0("Comparison samples (", pretty(sub("^ComparisonOnly_", "", s)), ")")
  } else if (identical(s, "EDA")) {
    "All samples (EDA)"
  } else {
    paste0("All samples (", pretty(s), ")")
  }
}

.xp_params_table <- function(res) {
  cfg <- res$config
  rows <- list(
    c("Comparison", if (is.na(res$comparison)) "none" else res$comparison),
    c("Contrast",   if (!is.null(cfg$level)) paste(cfg$level, "vs", cfg$base) else "n/a"),
    c("Design formula", cfg$model %||% "n/a"),
    c("Contrast factor", cfg$main_condition %||% "n/a"),
    c("Test", cfg$test %||% "n/a"),
    c("Reduced formula (LRT)", cfg$reduced %||% "n/a"),
    c("LFC shrinkage", cfg$shrink_method %||% "n/a"),
    c("Adjusted p-value cut-off", cfg$padj_cutoff),
    c("|log2FC| cut-off", cfg$lfc_cutoff),
    c("Count type", cfg$count_type %||% "n/a"),
    c("Output directory", res$out_dir),
    c("Run written at", cfg$run$written_at %||% "n/a"),
    c("ExpressOM version", cfg$run$expressom_version %||% "n/a"),
    c("Report generated", format(Sys.time(), "%Y-%m-%d %H:%M"))
  )
  data.frame(
    Parameter = vapply(rows, function(r) r[[1L]], character(1)),
    Value     = vapply(rows, function(r) as.character(r[[2L]]), character(1)),
    stringsAsFactors = FALSE
  )
}

# ---- leaves ------------------------------------------------------------------

.xp_leaf_overview <- function(res, st) {
  function() {
    tg <- htmltools::tags
    boxes <- list()
    add <- function(...) boxes[[length(boxes) + 1L]] <<- .xp_box(...)
    cfg <- res$config

    s <- .xp_dge_summary(res, st)
    if (!is.null(s)) {
      add("Comparison", paste(cfg$level %||% "?", "vs", cfg$base %||% "?"),
          paste(cfg$test, "test"))
      add("Genes tested", .xp_fmt(s$n_tested), "after independent filtering")
      if (s$is_lrt) {
        add("Significant genes", .xp_fmt(s$n_sig), paste0("LRT, padj < ", s$padj))
      } else {
        add("Upregulated",   .xp_fmt(s$n_up),   paste0("padj < ", s$padj, ", log2FC > ", s$lfc))
        add("Downregulated", .xp_fmt(s$n_down), paste0("padj < ", s$padj, ", log2FC < -", s$lfc))
      }
    }
    if (length(res$pca)) add("Samples", .xp_fmt(nrow(res$pca[[1L]]$scores)), "in the PCA")
    inv <- res$functional$inventory
    if (nrow(inv)) add("Enrichment tables", .xp_fmt(nrow(inv)), paste(unique(inv$source), collapse = " / "))
    iso <- res$isoform
    if (!is.null(iso$dte)) {
      p <- .rep_num(iso$dte$padj)
      add("DTE transcripts", .xp_fmt(sum(!is.na(p) & p < .xp_cutoffs(res, st)$padj)), "significant")
    }
    if (!is.null(iso$dtu)) {
      pc <- .rep_first_col(iso$dtu, c("adj_pvalue", "padj", "FDR"))
      if (!is.na(pc)) {
        p <- .rep_num(iso$dtu[[pc]])
        add("DTU features", .xp_fmt(sum(!is.na(p) & p < .xp_cutoffs(res, st)$padj)), "significant")
      }
    }

    avail <- function(ok, label, detail = NULL) {
      tg$li(
        class = if (ok) "xp-ok" else "xp-missing",
        tg$span(class = "xp-tick", if (ok) "\u2713" else "\u2013"),
        label,
        if (!is.null(detail)) tg$span(class = "xp-detail", detail)
      )
    }

    tg_out <- list(
      tg$div(class = "xp-boxes", boxes),
      tg$h4("Contents of this report"),
      tg$ul(
        class = "xp-avail",
        avail(!is.null(res$dge), "Differential expression",
              if (!is.null(res$dge)) paste(.xp_fmt(nrow(res$dge$results)), "genes")),
        avail(length(res$pca) > 0L, "PCA / sample QC",
              if (length(res$pca)) paste(length(res$pca), "plot(s)")),
        avail(nrow(inv) > 0L, "Functional enrichment (ORA / GSEA / FGSEA / SPIA / TF)",
              if (nrow(inv)) paste(nrow(inv), "table(s)")),
        avail(!is.null(iso), "Isoform analysis (DTE / DTU / switches)"),
        avail(nrow(res$figures) > 0L, "Static figures",
              if (nrow(res$figures)) paste(nrow(res$figures), "file(s)"))
      )
    )

    if (!is.null(res$dge$deg_counts)) {
      tg_out <- c(tg_out, list(
        tg$h4("Differentially expressed genes by cut-off"),
        .xp_two_col(
          iplot_deg_counts(res$dge$deg_counts, cfg$volcano_colors, cfg$test),
          itable(res$dge$deg_counts, page_length = 10L, filename = "DEG_counts")
        )
      ))
    }
    if (length(res$notes)) {
      tg_out <- c(tg_out, list(tg$h4("Notes"), tg$ul(lapply(res$notes, tg$li))))
    }
    htmltools::tagList(tg_out)
  }
}

.xp_leaf_pca <- function(res, st, entry) {
  function() {
    tg <- htmltools::tags
    cfg <- res$config
    gcol <- .rep_pca_group_col(entry$scores, cfg$main_condition)
    gcols <- NULL
    if (!is.na(gcol)) {
      gcols <- character(0)
      if (!is.null(cfg$level)) gcols[as.character(cfg$level)] <- unname(cfg$pca_colors[["level"]])
      if (!is.null(cfg$base))  gcols[as.character(cfg$base)]  <- unname(cfg$pca_colors[["base"]])
    }
    pcs <- entry$pcs
    plots <- list(
      iplot_pca(entry$scores, entry$variance, pcs[[1L]], pcs[[2L]], gcol, gcols,
                title = paste(pcs[[1L]], "vs", pcs[[2L]]))
    )
    if (!is.null(entry$variance)) {
      plots <- c(plots, list(iplot_scree(entry$variance)))
    }
    htmltools::tagList(
      do.call(.xp_two_col, plots),
      tg$h4("Sample scores"),
      itable(entry$scores, page_length = 10L, filename = entry$name)
    )
  }
}

.xp_leaf_volcano <- function(res, st) {
  function() {
    co <- .xp_cutoffs(res, st)
    iplot_volcano(
      res$dge$results, padj_cutoff = co$padj, lfc_cutoff = co$lfc,
      colors = res$config$volcano_colors, test_type = res$config$test,
      max_points = st$max_points, highlight = st$highlight_genes,
      title = paste0("Volcano plot: ", res$config$level %||% "", " vs ", res$config$base %||% "")
    )
  }
}

.xp_leaf_ma <- function(res, st) {
  function() {
    co <- .xp_cutoffs(res, st)
    iplot_ma(
      res$dge$results, padj_cutoff = co$padj, lfc_cutoff = co$lfc,
      colors = res$config$volcano_colors, test_type = res$config$test,
      max_points = st$max_points
    )
  }
}

.xp_leaf_heatmap <- function(res, st) {
  function() {
    ex <- res$dge$expression
    iplot_heatmap(
      ex$matrix, groups = ex$groups, colors = res$config$volcano_colors,
      max_genes = st$heatmap_genes
    )
  }
}

.xp_leaf_de_table <- function(res, st) {
  function() {
    tbl <- res$dge$results
    tbl <- tbl[order(tbl$padj, na.last = TRUE), , drop = FALSE]
    htmltools::tagList(
      .xp_note("Genes are sorted by adjusted p-value. Use the column filters to ",
               "select genes of interest; the buttons export what is currently ",
               "displayed."),
      itable(tbl, page_length = 15L, max_rows = st$table_rows,
             filename = paste0("DEgenes_", res$comparison))
    )
  }
}

.xp_leaf_gene_lists <- function(res) {
  function() {
    tg <- htmltools::tags
    f <- res$dge$sig_lists
    if (!length(f)) return(.xp_note("No per-cut-off DEG lists were found."))
    tg$div(
      .xp_note("DEG lists written by the pipeline (in `DE_raw_results/`):"),
      tg$ul(lapply(f, function(x) tg$li(tg$code(x))))
    )
  }
}

.xp_leaf_enrichment <- function(res, st, path, kind, show_plot = TRUE) {
  function() {
    tb <- res$functional$tables[[path]]
    co <- .xp_cutoffs(res, st)
    std <- .rep_standardise_enrichment(tb, kind)
    n_sig <- if (is.null(std)) NA_integer_ else sum(std$padj < co$padj, na.rm = TRUE)
    htmltools::tagList(
      .xp_note(
        paste0(
          .xp_fmt(nrow(tb)), " row(s) in the table",
          if (!is.na(n_sig)) paste0("; ", .xp_fmt(n_sig), " with adjusted p < ", co$padj),
          ". Source file: "
        ),
        htmltools::tags$code(path)
      ),
      if (isTRUE(show_plot) && !is.null(std)) {
        iplot_enrichment(
          std, top_n = st$top_terms, padj_cutoff = co$padj,
          colors = res$config$volcano_colors,
          title = paste0("Top terms \u2013 ", basename(path))
        )
      } else if (is.null(std)) {
        .xp_note("This table has no recognised term / adjusted p-value columns; showing the raw table only.")
      } else {
        .xp_note("Interactive plot omitted here to keep the report light \u2014 see the table below.")
      },
      itable(tb, page_length = 10L, max_rows = st$table_rows,
             filename = sub("\\.[^.]+$", "", basename(path)))
    )
  }
}

# Merge every collection of a method into one table. `paths`/`collections`
# are parallel vectors; a `collection` column is added.
.xp_leaf_enrichment_merged <- function(res, st, paths, collections, kind) {
  function() {
    merged <- do.call(rbind, lapply(seq_along(paths), function(i) {
      tb <- res$functional$tables[[paths[[i]]]]
      if (is.null(tb) || !nrow(tb)) return(NULL)
      tb$collection <- collections[[i]]
      tb
    }))
    if (is.null(merged) || !nrow(merged)) {
      return(.xp_note("No rows to show across this method."))
    }
    htmltools::tagList(
      .xp_note(sprintf("%s rows across %d collections. Filter by 'collection'.",
                       .xp_fmt(nrow(merged)), length(paths))),
      itable(merged, page_length = 15L, max_rows = st$table_rows,
             filename = paste0("merged_", kind))
    )
  }
}

.xp_leaf_dte <- function(res, st) {
  function() {
    co <- .xp_cutoffs(res, st)
    dte <- res$isoform$dte
    lab <- .rep_first_col(dte, c("gene", "transcript_id"))
    htmltools::tagList(
      iplot_volcano(
        dte, padj_cutoff = co$padj, lfc_cutoff = co$lfc,
        colors = res$config$volcano_colors,
        label_col = if (is.na(lab)) "gene" else lab,
        max_points = st$max_points,
        title = "DTE volcano plot (transcript level)"
      ),
      itable(dte[order(.rep_num(dte$padj), na.last = TRUE), , drop = FALSE],
             page_length = 15L, max_rows = st$table_rows, filename = "DTE_results")
    )
  }
}

.xp_leaf_dtu <- function(res, st) {
  function() {
    co <- .xp_cutoffs(res, st)
    dtu <- res$isoform$dtu
    pc <- .rep_first_col(dtu, c("adj_pvalue", "padj", "FDR", "pvalue"))
    raw_p <- .rep_first_col(dtu, c("pvalue", "p_value"))
    ord <- if (is.na(pc)) seq_len(nrow(dtu)) else order(.rep_num(dtu[[pc]]), na.last = TRUE)
    htmltools::tagList(
      .xp_two_col(
        iplot_dtu(dtu, top_n = 30L, padj_cutoff = co$padj),
        if (!is.na(raw_p)) iplot_pvalue_hist(dtu[[raw_p]], "DTU p-value distribution")
      ),
      itable(dtu[ord, , drop = FALSE], page_length = 15L,
             max_rows = st$table_rows, filename = "DTU_results")
    )
  }
}

.xp_leaf_table <- function(df, filename, note = NULL, st) {
  function() {
    htmltools::tagList(
      if (!is.null(note)) .xp_note(note),
      itable(df, page_length = 15L, max_rows = st$table_rows, filename = filename)
    )
  }
}

.xp_leaf_switch_plots <- function(res, st) {
  function() {
    plots <- res$isoform$plots
    if (!length(plots)) {
      return(.xp_note("No isoform-switch / DTU / DTE PNG plots were found."))
    }
    n_all  <- length(plots)
    n_show <- as.integer(st$max_isoform_plots %||% 6L)
    plots  <- utils::head(plots, n_show)
    uris <- if (.Platform$OS.type != "windows" && length(plots) > 1L) {
      parallel::mclapply(plots, .xp_img_uri, figures_dir = st$figures_dir,
                         mc.cores = min(parallel::detectCores() - 1L, 8L, length(plots)))
    } else {
      lapply(plots, .xp_img_uri, figures_dir = st$figures_dir)
    }
    cap <- sub("\\.[^.]+$", "", basename(plots))
    htmltools::tagList(
      .xp_note(
        "Static plots produced by the isoform analysis (click an image to enlarge). ",
        if (n_all > length(plots)) paste0("Showing ", length(plots), " of ", n_all, ".")
      ),
      .xp_gallery(uris, cap)
    )
  }
}

.xp_leaf_methods <- function(res) {
  function() {
    txt <- res$isoform$methods
    htmltools::tagList(
      .xp_note("Statistical notes recorded by the isoform analysis:"),
      htmltools::tags$ul(lapply(txt, htmltools::tags$li))
    )
  }
}

.xp_leaf_figures <- function(res, st, group) {
  function() {
    fig  <- res$figures[res$figures$group == group, , drop = FALSE]
    pdfs <- fig[fig$type == "pdf" & fig$size >= 2000, , drop = FALSE]
    pngs <- fig[fig$type == "png", , drop = FALSE]
    tg <- htmltools::tags
    if (!isTRUE(st$static_figures)) {
      return(.xp_note("Static figures are disabled for this report."))
    }

    items <- list()
    if (nrow(pngs)) {
      n_take <- min(nrow(pngs), st$max_static_figures)
      uris <- if (.Platform$OS.type != "windows" && n_take > 1L) {
        parallel::mclapply(pngs$path[seq_len(n_take)], .xp_img_uri,
                           figures_dir = st$figures_dir,
                           mc.cores = min(parallel::detectCores() - 1L, 8L, n_take))
      } else {
        lapply(pngs$path[seq_len(n_take)], .xp_img_uri,
               figures_dir = st$figures_dir)
      }
      for (i in seq_along(uris)) {
        if (!is.null(uris[[i]])) {
          items[[length(items) + 1L]] <- list(uri = uris[[i]], cap = pngs$name[[i]])
        }
      }
    }
    have_pdftools <- requireNamespace("pdftools", quietly = TRUE)
    if (nrow(pdfs) && have_pdftools) {
      room <- st$max_static_figures - length(items)
      n_take <- max(0L, min(nrow(pdfs), room))
      if (n_take > 0L) {
        batch <- pdfs$path[seq_len(n_take)]
        uris <- if (.Platform$OS.type != "windows" && n_take > 1L) {
          parallel::mclapply(batch, .xp_pdf_uri, dpi = st$static_dpi,
                             figures_dir = st$figures_dir,
                             mc.cores = min(parallel::detectCores() - 1L, 8L, n_take))
        } else {
          lapply(batch, .xp_pdf_uri, dpi = st$static_dpi,
                 figures_dir = st$figures_dir)
        }
        for (i in seq_along(uris)) {
          if (!is.null(uris[[i]])) {
            items[[length(items) + 1L]] <- list(uri = uris[[i]],
                                                 cap = .xp_strip_pdf(pdfs$name[[i]]))
          }
        }
      }
    }
    listing <- itable(
      data.frame(file = .rep_rel(fig$path, res$out_dir), type = fig$type,
                 size_kb = round(fig$size / 1024, 1), stringsAsFactors = FALSE),
      page_length = 10L, filename = paste0("figures_", group)
    )
    htmltools::tagList(
      if (length(items)) {
        .xp_gallery(lapply(items, `[[`, "uri"),
                    vapply(items, `[[`, character(1), "cap"))
      } else if (nrow(pdfs) && !have_pdftools) {
        .xp_note("Install the 'pdftools' package to preview the PDF figures here. ",
                 "The files are listed below.")
      },
      tg$h4("Files"), listing
    )
  }
}

.xp_strip_pdf <- function(x) sub("\\.pdf$", "", x, ignore.case = TRUE)

.xp_leaf_logs <- function(res) {
  function() {
    tg <- htmltools::tags
    blocks <- list()
    dl <- res$logs$dge_params
    if (!is.null(dl)) {
      blocks <- c(blocks, list(
        tg$h4("Differential expression parameters"),
        tg$pre(paste(utils::capture.output(print(unlist(dl))), collapse = "\n"))
      ))
    }
    det <- function(title, lines) {
      if (is.null(lines) || !length(lines)) return(NULL)
      tg$details(tg$summary(title),
                 tg$pre(paste(utils::head(lines, 400L), collapse = "\n")))
    }
    blocks <- c(blocks, list(
      det("Warnings recorded during the run", res$logs$warnings),
      det("Session information", res$logs$session_info)
    ))
    if (length(res$logs$isoform_logs)) {
      blocks <- c(blocks, list(
        tg$h4("Isoform logs"),
        tg$ul(lapply(.rep_rel(res$logs$isoform_logs, res$out_dir),
                     function(x) tg$li(tg$code(x))))
      ))
    }
    if (length(blocks) == 0L) return(.xp_note("No log files were found in Log/."))
    htmltools::tagList(blocks)
  }
}

# ---- the tree ----------------------------------------------------------------

#' Build the tab tree of the report / Shiny viewer
#'
#' @param res An `expressom_results` object ([load_expressom_results()]).
#' @param settings A list from [expressom_report_settings()].
#' @return A list of nodes (see the file header of `mod_report_tabs.R`).
#' @export
expressom_report_tabs <- function(res, settings = expressom_report_settings()) {
  stopifnot(inherits(res, "expressom_results"))
  st <- settings

  plotly_budget   <- as.integer(st$plotly_enrichment_max %||% 3L)
  merge_threshold <- as.integer(st$merge_enrichment_threshold %||% 3L)
  pca_budget      <- as.integer(st$plotly_pca_max %||% 2L)
  plotly_used     <- 0L

  take_plotly <- function() {
    if (plotly_used < plotly_budget) {
      plotly_used <<- plotly_used + 1L
      TRUE
    } else FALSE
  }

  # -- QC / PCA
  pca_keep <- utils::head(names(res$pca), pca_budget)
  pca_drop <- setdiff(names(res$pca), pca_keep)
  pca_nodes <- lapply(pca_keep, function(nm) {
    .xp_node(.xp_pca_title(nm), content = .xp_leaf_pca(res, st, res$pca[[nm]]))
  })
  if (length(pca_drop)) {
    pca_nodes <- c(pca_nodes, list(
      .xp_node("Other PCAs (tables only)", content = function() {
        htmltools::tagList(lapply(pca_drop, function(nm) {
          htmltools::tagList(
            htmltools::tags$h4(.xp_pca_title(nm)),
            itable(res$pca[[nm]]$scores, page_length = 10L,
                   max_rows = st$table_rows, filename = nm)
          )
        }))
      })
    ))
  }
  qc <- if (length(pca_nodes)) {
    do.call(.xp_node, c(list("Samples & QC"), pca_nodes))
  } else NULL

  # -- DGE
  dge <- NULL
  if (!is.null(res$dge)) {
    dge <- .xp_node(
      "Differential expression",
      .xp_node("Volcano", content = .xp_leaf_volcano(res, st)),
      .xp_node("MA plot", content = .xp_leaf_ma(res, st)),
      if (!is.null(res$dge$expression)) {
        .xp_node("Top genes heatmap", content = .xp_leaf_heatmap(res, st))
      },
      .xp_node("Results table", content = .xp_leaf_de_table(res, st)),
      .xp_node("DEG lists", content = .xp_leaf_gene_lists(res))
    )
  }

  # -- functional
  inv <- res$functional$inventory
  func <- NULL
  if (nrow(inv)) {
    by_source <- lapply(split(inv, inv$source), function(si) {
      method_nodes <- lapply(split(si, si$method), function(mi) {
        n_col <- nrow(mi)
        if (n_col > merge_threshold) {
          children <- list(
            .xp_node(sprintf("Merged table (%d collections)", n_col),
                     content = .xp_leaf_enrichment_merged(
                       res, st, mi$path, mi$collection, mi$kind[[1]]))
          )
          for (k in seq_len(n_col)) {
            children <- c(children, list(
              .xp_node(mi$collection[[k]],
                       content = .xp_leaf_enrichment(res, st, mi$path[[k]],
                                                     mi$kind[[k]],
                                                     show_plot = take_plotly()))
            ))
          }
          do.call(.xp_node, c(list(unique(mi$method)), children))
        } else {
          leaves <- lapply(seq_len(n_col), function(k) {
            .xp_node(mi$collection[[k]],
                     content = .xp_leaf_enrichment(res, st, mi$path[[k]],
                                                   mi$kind[[k]],
                                                   show_plot = take_plotly()))
          })
          if (length(leaves) == 1L) {
            leaves[[1L]]$title <- unique(mi$method)
            leaves[[1L]]
          } else {
            do.call(.xp_node, c(list(unique(mi$method)), leaves))
          }
        }
      })
      do.call(.xp_node, c(list(unique(si$source)), unname(method_nodes)))
    })
    func <- if (length(by_source) == 1L) {
      do.call(.xp_node, c(list("Functional analysis"), unname(by_source[[1L]]$children)))
    } else {
      src_titles <- c(DGE = "Gene level", DTE = "DTE", DTU = "DTU")
      by_source <- unname(lapply(by_source, function(n) {
        n$title <- src_titles[[n$title]] %||% n$title
        n
      }))
      do.call(.xp_node, c(list("Functional analysis"), by_source))
    }
  }

  # -- isoform
  iso <- res$isoform
  isoform <- NULL
  if (!is.null(iso)) {
    isoform <- .xp_node(
      "Isoforms",
      if (!is.null(iso$dte)) .xp_node("DTE", content = .xp_leaf_dte(res, st)),
      if (!is.null(iso$dtu)) .xp_node("DTU", content = .xp_leaf_dtu(res, st)),
      if (!is.null(iso$dexseq)) {
        .xp_node("DEXSeq", content = .xp_leaf_table(
          iso$dexseq, "DEXSeq_results",
          "DEXSeq (transcript features as bins) results.", st))
      },
      if (!is.null(iso$transcript_pca)) {
        .xp_node("Variable transcripts", content = .xp_leaf_table(
          iso$transcript_pca, "PCA_transcripts",
          "Most variable transcripts used for the transcript-level PCA.", st))
      },
      if (length(iso$plots)) .xp_node("Plots", content = .xp_leaf_switch_plots(res, st)),
      if (length(iso$methods)) .xp_node("Methods", content = .xp_leaf_methods(res))
    )
  }

  # -- static figures
  figs <- NULL
  if (nrow(res$figures) && isTRUE(st$static_figures)) {
    grp <- unique(res$figures$group)
    figs <- do.call(.xp_node, c(
      list("Figures"),
      lapply(grp, function(g) {
        .xp_node(if (identical(g, ".")) "(root)" else g,
                 content = .xp_leaf_figures(res, st, g))
      })
    ))
  }

  # -- run info
  info <- .xp_node(
    "Run info",
    .xp_node("Parameters", content = .xp_leaf_table(
      .xp_params_table(res), "parameters", NULL, st)),
    .xp_node("Logs", content = .xp_leaf_logs(res)),
    .xp_node("Files", content = .xp_leaf_table(
      res$files, "file_inventory",
      "Every file written to the output directory.", st))
  )

  tree <- list(
    .xp_node("Overview", content = .xp_leaf_overview(res, st)),
    qc, dge, func, isoform, figs, info
  )
  tree[!vapply(tree, is.null, logical(1))]
}

# ---- rendering: htmltools engine ---------------------------------------------

.xp_clean_title <- function(x) {
  x <- gsub("[{}#\\[\\]]", "", x, perl = TRUE)
  trimws(x)
}

.xp_safe_leaf <- function(fn) {
  force(fn)
  function() {
    tryCatch(fn(), error = function(e) {
      htmltools::tags$div(class = "xp-problems",
                          "This tab could not be built: ", conditionMessage(e))
    })
  }
}

# Build one tabset (nav + content) from a list of nodes; returns a list of
# `nav` and `content` tag objects. Ids are made unique per call.
.xp_render_tabset <- function(nodes, classes = "nav-pills", counter = NULL) {
  tg <- htmltools::tags
  if (is.null(counter)) counter <- new.env(parent = emptyenv())
  if (is.null(counter$n)) counter$n <- 0L
  next_id <- function() {
    counter$n <- counter$n + 1L
    paste0("xp-tab-", counter$n)
  }
  navs  <- list()
  panes <- list()
  for (i in seq_along(nodes)) {
    nd      <- nodes[[i]]
    title   <- .xp_clean_title(nd$title)
    pane_id <- next_id()
    body <- if (is.function(nd$content)) {
      .xp_safe_leaf(nd$content)()
    } else {
      nested <- .xp_render_tabset(nd$children, classes, counter)
      tg$div(
        nested$nav,
        nested$content
      )
    }
    navs[[i]] <- tg$li(
      role = "presentation",
      class = if (i == 1L) "active" else NULL,
      tg$a(
        href = paste0("#", pane_id),
        `aria-controls` = pane_id,
        role = "tab",
        `data-toggle` = "tab",
        title
      )
    )
    panes[[i]] <- tg$div(
      role = "tabpanel",
      class = paste("tab-pane fade", if (i == 1L) "in active" else ""),
      id = pane_id,
      body
    )
  }
  list(
    nav     = tg$ul(class = paste("nav", classes), role = "tablist", navs),
    content = tg$div(class = "tab-content", panes)
  )
}

# Full standalone HTML page for the htmltools engine. Bootstrap 3 is pulled
# from CDN, plotly/DT widget dependencies are handled by htmltools::save_html().
.xp_build_html <- function(nodes, res, st) {
  tg <- htmltools::tags
  counter <- new.env(parent = emptyenv()); counter$n <- 0L

  # Top-level tabsets are pills (matches the pandoc template)
  top <- .xp_render_tabset(nodes, classes = "nav-pills xp-top", counter = counter)

  tg$html(
    lang = "en",
    tg$head(
      tg$meta(charset = "utf-8"),
      tg$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tg$title(st$title %||% "ExpressOM report"),
      tg$link(
        rel = "stylesheet",
        href = "https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css"
      ),
      tg$style(htmltools::HTML(paste(.xp_css(), .xp_dark_css_text(), sep = "\n"))),
      tg$script(src = "https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js"),
      tg$script(src = "https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js")
    ),
    tg$body(
      tg$div(
        class = "container-fluid",
        .xp_header_html(res, st),
        top$nav,
        top$content
      )
    )
  )
}

# ---- rendering: pandoc engine (legacy, used when self_contained = TRUE) ------

.xp_emit_tabs <- function(nodes, level = 2L) {
  counter <- 0L
  leaf_fns <- list()
  emit <- function(nodes, level) {
    out <- character(0)
    for (nd in nodes) {
      title <- .xp_clean_title(nd$title)
      if (.xp_is_leaf(nd)) {
        counter <<- counter + 1L
        id <- paste0("leaf_", counter)
        leaf_fns[[id]] <<- .xp_safe_leaf(nd$content)
        child <- paste0(
          paste0(strrep("#", level), " ", title, "\n\n"),
          "```{r, echo=FALSE, results='markup'}\n",
          "leaf_fns[['", id, "']]()\n",
          "```\n"
        )
        out <- c(out, knitr::knit_child(
          text = child, envir = environment(), quiet = TRUE
        ))
      } else {
        out <- c(
          out,
          paste0("\n\n", strrep("#", level), " ", title, " {.tabset .tabset-pills}\n\n"),
          emit(nd$children, level + 1L)
        )
      }
    }
    out
  }
  emit(nodes, level)
}

# ---- CSS + theme -------------------------------------------------------------

.xp_css_light_vars <- function() {
  paste(
    "--xp-bg:#ffffff; --xp-fg:#1f2933; --xp-muted:#52606d; --xp-faint:#7b8794;",
    "--xp-accent:#0b5cad; --xp-heading:#0b3d75; --xp-border:#d9e2ec;",
    "--xp-card:#f5f8fb; --xp-card2:#e8eef4; --xp-ok:#1a7f4b; --xp-missing:#9aa5b1;",
    "--xp-bad-bg:#fdecea; --xp-bad-border:#f5c2bd; --xp-bad-fg:#8a1c12;",
    "--xp-img-bg:#ffffff; --xp-shadow:rgba(0,0,0,.45);",
    "--xp-plot-bg:#ffffff; --xp-plot-grid:#e5eaf0; --xp-plot-zero:#6b7a8a;",
    "--xp-ref-line:#666666; --xp-contrast:#000000;"
  )
}

.xp_css_dark_vars <- function() {
  paste(
    "--xp-bg:#0f1720; --xp-fg:#e4e9ef; --xp-muted:#a7b4c2; --xp-faint:#8794a3;",
    "--xp-accent:#6cb2ff; --xp-heading:#a9d2fb; --xp-border:#2c3a4a;",
    "--xp-card:#17212d; --xp-card2:#1f2c3b; --xp-ok:#4cc38a; --xp-missing:#6b7a8a;",
    "--xp-bad-bg:#3a1a1a; --xp-bad-border:#7a2e2a; --xp-bad-fg:#ffb4ab;",
    "--xp-img-bg:#ffffff; --xp-shadow:rgba(0,0,0,.8);",
    "--xp-plot-bg:#0f1720; --xp-plot-grid:#263343; --xp-plot-zero:#6f8091;",
    "--xp-ref-line:#9fb3c8; --xp-contrast:#f5f7fa;"
  )
}

.xp_css <- function() {
  paste0(
    ":root { ", .xp_css_light_vars(), " }\n",
    "html[data-theme='dark'] { ", .xp_css_dark_vars(), " }\n",
    "@media print { html[data-theme='dark'] { ", .xp_css_light_vars(), " } }\n",
    "
.main-container { max-width: 1480px !important; }
body { font-family: -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; color: var(--xp-fg); }
.container-fluid { max-width: 1480px; padding: 0 16px; }
.section.level1.xp-top > h1 { display: none; }
.xp-header { display:flex; flex-wrap:wrap; align-items:center; gap:6px 18px; padding: 16px 0 8px; border-bottom: 2px solid var(--xp-accent); margin-bottom: 14px; }
.xp-header h1 { margin: 0; font-size: 28px; font-weight: 650; color: var(--xp-heading); }
.xp-header .xp-sub { color: var(--xp-muted); font-size: 14px; }
.xp-theme-toggle { margin-left:auto; display:inline-flex; align-items:center; gap:6px; padding:5px 12px; font-size:13px; line-height:1.2; color: var(--xp-fg); background: var(--xp-card); border:1px solid var(--xp-border); border-radius:16px; cursor:pointer; }
.xp-theme-toggle:hover { background: var(--xp-card2); }
.xp-theme-toggle svg { flex:none; }
.xp-theme-toggle .xp-ico-sun { display:none; }
html[data-theme='dark'] .xp-theme-toggle .xp-ico-sun { display:inline; }
html[data-theme='dark'] .xp-theme-toggle .xp-ico-moon { display:none; }
@media print { .xp-theme-toggle { display:none; } }
.nav-pills > li > a { padding: 6px 14px; border-radius: 16px; font-size: 14px; }
.nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus { background: var(--xp-accent); }
.nav-pills { display: flex; flex-wrap: wrap; gap: 4px; }
.tab-content { padding-top: 12px; }
.xp-problems { background: var(--xp-bad-bg); border:1px solid var(--xp-bad-border); color: var(--xp-bad-fg); border-radius:8px; padding:8px 12px; margin-top:10px; font-size:13px; }
.xp-boxes { display:grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 12px; margin: 6px 0 18px; }
.xp-box { border: 1px solid var(--xp-border); border-radius: 10px; padding: 12px 16px; background: var(--xp-card); }
.xp-box-label { font-size: 12px; text-transform: uppercase; letter-spacing: .04em; color: var(--xp-muted); }
.xp-box-value { font-size: 26px; font-weight: 650; color: var(--xp-heading); line-height: 1.25; }
.xp-box-sub { font-size: 12px; color: var(--xp-faint); }
.xp-grid { display:grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 14px; align-items: start; margin-bottom: 14px; }
.xp-note { color: var(--xp-muted); font-style: italic; margin: 6px 0 10px; }
.xp-avail { list-style:none; padding-left:0; }
.xp-avail li { padding: 3px 0; }
.xp-tick { display:inline-block; width: 22px; font-weight:700; }
.xp-ok .xp-tick { color: var(--xp-ok); }
.xp-missing { color: var(--xp-missing); }
.xp-detail { color: var(--xp-faint); margin-left:10px; font-size: 12px; }
.xp-kv th { text-align:left; padding: 3px 18px 3px 0; color: var(--xp-muted); font-weight:600; }
.xp-gallery { display:grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 12px; }
.xp-gallery figure { margin:0; }
.xp-gallery img { width:100%; border:1px solid var(--xp-border); border-radius:6px; cursor: zoom-in; background: var(--xp-img-bg); }
.xp-gallery figcaption { font-size: 12px; color: var(--xp-muted); word-break: break-all; padding-top: 2px; }
img.xp-zoom { position: fixed; inset: 2%; width: 96%; height: 96%; object-fit: contain; z-index: 9999; background: var(--xp-img-bg); box-shadow: 0 0 60px var(--xp-shadow); cursor: zoom-out; }
details { margin: 8px 0; } details > summary { cursor: pointer; font-weight: 600; }
pre { max-height: 420px; overflow:auto; font-size: 12px; }
table.dataTable { font-size: 13px; }
"
  )
}

# Embedded stylesheet/script (no disk lookup needed). Legacy fallback to
# .expressom_rmd_path() is kept only for the two .Rmd templates, which do exist
# on disk if the caller overrides them.
.xp_asset_text <- function(file) {
  if (identical(file, "expressom_dark.css"))  return(.xp_dark_css_text())
  if (identical(file, "expressom_theme.js")) return(.xp_theme_js_text())
  path <- tryCatch(.expressom_rmd_path(file), error = function(e) NULL)
  if (is.null(path)) return("")
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.xp_theme_setting <- function(st) {
  th <- st$theme
  if (is.character(th) && length(th) == 1L && th %in% c("auto", "light", "dark")) th else "auto"
}

.xp_theme_toggle <- function() {
  htmltools::HTML(paste0(
    '<button type="button" id="xp-theme-toggle" class="xp-theme-toggle" ',
    'title="Switch to dark mode" aria-pressed="false">',
    '<svg class="xp-ico-moon" viewBox="0 0 24 24" width="15" height="15" aria-hidden="true">',
    '<path fill="currentColor" d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/></svg>',
    '<svg class="xp-ico-sun" viewBox="0 0 24 24" width="15" height="15" aria-hidden="true">',
    '<circle cx="12" cy="12" r="4.2" fill="currentColor"/>',
    '<path stroke="currentColor" stroke-width="2" stroke-linecap="round" ',
    'd="M12 2v2.5M12 19.5V22M2 12h2.5M19.5 12H22M4.9 4.9l1.8 1.8M17.3 17.3l1.8 1.8M4.9 19.1l1.8-1.8M17.3 6.7l1.8-1.8"/></svg>',
    '<span class="xp-theme-label">Dark</span></button>'
  ))
}

.xp_header_html <- function(res, st) {
  tg <- htmltools::tags
  sub <- c(
    if (!is.na(res$comparison)) paste0(res$config$level %||% "", " vs ", res$config$base %||% ""),
    if (!is.null(res$config$model)) paste0("design ", res$config$model),
    if (!is.null(res$config$test))  paste0(res$config$test, " test"),
    format(Sys.time(), "%Y-%m-%d %H:%M")
  )
  htmltools::tagList(
    tg$style(htmltools::HTML(paste(.xp_css(), .xp_dark_css_text(), sep = "\n"))),
    tg$div(
      class = "xp-header",
      tg$h1(st$title %||% "ExpressOM report"),
      tg$span(class = "xp-sub", paste(sub, collapse = "  \u00b7  ")),
      .xp_theme_toggle()
    ),
    tg$script(htmltools::HTML(.xp_theme_js_text(.xp_theme_setting(st))))
  )
}