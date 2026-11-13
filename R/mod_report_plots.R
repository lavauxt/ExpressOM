# mod_report_plots.R -- interactive (plotly / DT) counterparts of the PDF
# figures, built from the result tables returned by load_expressom_results().
#
# Every function takes plain data.frames, returns an htmlwidget and is used
# both by the self-contained report and by the Shiny app, so the two always
# show the same figures. plotly and DT are only needed at report time.

.rep_plotly_config <- function(p, filename = "expressom_plot") {
  p$sizingPolicy <- htmlwidgets::sizingPolicy(
    defaultWidth = "100%",
    defaultHeight = 520,
    knitr.figure = FALSE,
    knitr.defaultWidth = "100%",
    knitr.defaultHeight = "520px",
    browser.fill = FALSE,
    viewer.fill = FALSE,
    padding = 0
  )
  plotly::config(
    p,
    displaylogo = FALSE,
    responsive = TRUE,
    toImageButtonOptions = list(format = "png", filename = filename, scale = 2)
  )
}

# plotly.js only understands CSS colour names / hex codes, while R colours such
# as "red2" or "grey70" are R-specific: convert everything to #RRGGBB.
.rep_hex <- function(x) {
  if (is.null(x)) return(x)
  nm <- names(x)
  out <- vapply(as.character(x), function(col) {
    tryCatch(
      {
        m <- grDevices::col2rgb(col)
        grDevices::rgb(m[1L], m[2L], m[3L], maxColorValue = 255)
      },
      error = function(e) "#999999"
    )
  }, character(1), USE.NAMES = FALSE)
  names(out) <- nm
  out
}

.rep_hline <- function(y, color = "#666666") {
  list(
    type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
    line = list(color = color, width = 1, dash = "dash")
  )
}

.rep_vline <- function(x, color = "#666666") {
  list(
    type = "line", x0 = x, x1 = x, y0 = 0, y1 = 1, yref = "paper",
    line = list(color = color, width = 1, dash = "dash")
  )
}

.rep_empty_plot <- function(msg = "No data to display") {
  .rep_require("plotly")
  plotly::plot_ly(type = "scatter", mode = "markers") |>
    plotly::layout(
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE),
      annotations = list(list(
        text = msg, x = 0.5, y = 0.5, xref = "paper", yref = "paper",
        showarrow = FALSE, font = list(size = 14, color = "#666666")
      ))
    ) |>
    .rep_plotly_config()
}

.rep_trunc <- function(x, n = 60L) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 1L), "…"), x)
}

.rep_palette <- function(n) {
  if (n <= 0L) return(character(0))
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
            "#56B4E9", "#F0E442", "#999999")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

# ---- volcano / MA -------------------------------------------------------------

.rep_de_frame <- function(
  res_tbl,
  padj_cutoff,
  lfc_cutoff,
  test_type = "Wald",
  label_col = "gene",
  max_points = 2000L
) {
  lfc <- .rep_num(res_tbl$log2FoldChange)
  padj <- .rep_num(res_tbl$padj)
  base_mean <- if ("baseMean" %in% colnames(res_tbl)) {
    .rep_num(res_tbl$baseMean)
  } else {
    rep(NA_real_, length(lfc))
  }
  label <- if (label_col %in% colnames(res_tbl)) {
    as.character(res_tbl[[label_col]])
  } else {
    as.character(seq_along(lfc))
  }

  keep <- is.finite(lfc) & !is.na(padj)
  df <- data.frame(
    label = label[keep],
    lfc = lfc[keep],
    padj = padj[keep],
    base_mean = base_mean[keep],
    stringsAsFactors = FALSE
  )
  if (nrow(df) == 0L) return(df)
  df$mlog10 <- -log10(pmax(df$padj, .Machine$double.xmin))

  is_lrt <- identical(test_type, "LRT")
  sig <- df$padj < padj_cutoff
  if (is_lrt) {
    df$class <- ifelse(sig, "Significant (LRT)", "Not significant")
  } else {
    sig_lfc <- sig & abs(df$lfc) > lfc_cutoff
    df$class <- "Not significant"
    df$class[sig_lfc & df$lfc > 0] <- "Upregulated"
    df$class[sig_lfc & df$lfc < 0] <- "Downregulated"
  }

  if (nrow(df) > max_points) {
    ns <- which(df$class == "Not significant")
    n_keep <- max(0L, max_points - (nrow(df) - length(ns)))
    if (length(ns) > n_keep) {
      keep_idx <- if (requireNamespace("withr", quietly = TRUE)) {
        withr::with_seed(1L, sample.int(length(ns), n_keep))
      } else {
        sample.int(length(ns), n_keep)
      }
      drop <- ns[-keep_idx]
      df <- df[-drop, , drop = FALSE]
    }
  }

  df$lfc <- round(df$lfc, 3)
  df$padj <- signif(df$padj, 3)
  df$base_mean <- signif(df$base_mean, 4)
  df
}

.rep_de_hover <- function(df) {
  paste0(
    "<b>", df$label, "</b><br>log2FC: ", signif(df$lfc, 3),
    "<br>adj. p: ", signif(df$padj, 3),
    ifelse(is.na(df$base_mean), "", paste0("<br>baseMean: ", signif(df$base_mean, 4)))
  )
}

.rep_de_traces <- function(p, df, colors, test_type) {
  colors <- .rep_hex(colors)
  is_lrt <- identical(test_type, "LRT")
  cls <- if (is_lrt) {
    c("Not significant", "Significant (LRT)")
  } else {
    c("Not significant", "Downregulated", "Upregulated")
  }
  cols <- if (is_lrt) {
    c("Not significant" = unname(colors[["ns"]]),
      "Significant (LRT)" = unname(colors[["up"]]))
  } else {
    c(
      "Not significant" = unname(colors[["ns"]]),
      "Downregulated" = unname(colors[["down"]]),
      "Upregulated" = unname(colors[["up"]])
    )
  }
  trace_type <- if (nrow(df) > 500L) "scattergl" else "scatter"
  for (k in cls) {
    d <- df[df$class == k, , drop = FALSE]
    if (nrow(d) == 0L) next
    p <- plotly::add_trace(
      p, data = d, x = ~x_val, y = ~y_val,
      type = trace_type, mode = "markers",
      name = paste0(k, " (", nrow(d), ")"),
      text = .rep_de_hover(d), hoverinfo = "text",
      marker = list(
        color = cols[[k]],
        size = if (k == "Not significant") 4 else 6,
        opacity = if (k == "Not significant") 0.45 else 0.85
      )
    )
  }
  p
}

#' Interactive volcano plot
#' @export
iplot_volcano <- function(
  res_tbl, padj_cutoff = 0.05, lfc_cutoff = 1,
  colors = .de_direction_colors(), test_type = "Wald",
  label_col = "gene", max_points = 5000L, highlight = NULL, title = NULL
) {
  .rep_require("plotly")
  df <- .rep_de_frame(res_tbl, padj_cutoff, lfc_cutoff, test_type, label_col, max_points)
  if (nrow(df) == 0L) {
    return(.rep_empty_plot("No genes with a finite log2 fold change and adjusted p-value"))
  }
  df$x_val <- df$lfc
  df$y_val <- df$mlog10

  p <- plotly::plot_ly()
  p <- .rep_de_traces(p, df, colors, test_type)

  if (length(highlight)) {
    hl <- df[df$label %in% highlight, , drop = FALSE]
    if (nrow(hl)) {
      p <- plotly::add_trace(
        p, data = hl, x = ~x_val, y = ~y_val, type = "scatter",
        mode = "markers+text", text = hl$label, textposition = "top center",
        hovertext = .rep_de_hover(hl), hoverinfo = "text",
        name = paste0("Highlighted (", nrow(hl), ")"),
        meta = "xp-contrast",
        marker = list(color = "black", size = 9, symbol = "diamond-open",
                      line = list(width = 2))
      )
    }
  }

  shapes <- list(.rep_hline(-log10(padj_cutoff)))
  if (!identical(test_type, "LRT") && lfc_cutoff > 0) {
    shapes <- c(shapes, list(.rep_vline(-lfc_cutoff), .rep_vline(lfc_cutoff)))
  }

  p |>
    plotly::layout(
      title = list(text = title %||% "Volcano plot", x = 0.02),
      xaxis = list(title = "log2 fold change", zeroline = TRUE),
      yaxis = list(title = "-log10(adjusted p-value)"),
      shapes = shapes,
      legend = list(orientation = "h", y = -0.18),
      hovermode = "closest"
    ) |>
    .rep_plotly_config("volcano")
}

#' Interactive MA plot
#' @inheritParams iplot_volcano
#' @export
iplot_ma <- function(
  res_tbl, padj_cutoff = 0.05, lfc_cutoff = 1,
  colors = .de_direction_colors(), test_type = "Wald",
  label_col = "gene", max_points = 5000L, title = NULL
) {
  .rep_require("plotly")
  df <- .rep_de_frame(res_tbl, padj_cutoff, lfc_cutoff, test_type, label_col, max_points)
  df <- df[is.finite(df$base_mean) & df$base_mean > 0, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(.rep_empty_plot("No baseMean column available for the MA plot"))
  }
  df$x_val <- df$base_mean
  df$y_val <- df$lfc

  p <- .rep_de_traces(plotly::plot_ly(), df, colors, test_type)
  shapes <- list(.rep_hline(0, "#999999"))
  if (!identical(test_type, "LRT") && lfc_cutoff > 0) {
    shapes <- c(shapes, list(.rep_hline(lfc_cutoff), .rep_hline(-lfc_cutoff)))
  }
  p |>
    plotly::layout(
      title = list(text = title %||% "MA plot", x = 0.02),
      xaxis = list(title = "mean of normalised counts", type = "log"),
      yaxis = list(title = "log2 fold change"),
      shapes = shapes,
      legend = list(orientation = "h", y = -0.18)
    ) |>
    .rep_plotly_config("ma_plot")
}

# ---- PCA ------------------------------------------------------------------------

#' Interactive PCA
#' @export
iplot_pca <- function(
  scores, variance = NULL, x = "PC1", y = "PC2",
  group_col = NA_character_, group_colors = NULL, title = NULL
) {
  .rep_require("plotly")
  group_colors <- .rep_hex(group_colors)
  if (is.null(scores) || !all(c(x, y) %in% colnames(scores))) {
    return(.rep_empty_plot("PCA scores not available"))
  }
  df <- scores
  df$.id <- as.character(if ("plot_label" %in% colnames(df))
                         df$plot_label else df$sample_label)
  df$.grp <- if (!is.na(group_col) && group_col %in% colnames(df)) {
    as.character(df[[group_col]])
  } else {
    "All samples"
  }
  groups <- unique(df$.grp)
  cols <- stats::setNames(.rep_palette(length(groups)), groups)
  if (length(group_colors)) {
    hit <- intersect(names(group_colors), groups)
    cols[hit] <- group_colors[hit]
  }

  axis_title <- function(pc) {
    v <- if (!is.null(variance) && pc %in% names(variance))
      paste0(": ", variance[[pc]], "% variance") else ""
    paste0(pc, v)
  }

  p <- plotly::plot_ly()
  for (g in groups) {
    d <- df[df$.grp == g, , drop = FALSE]
    p <- plotly::add_trace(
      p,
      x = d[[x]], y = d[[y]], type = "scatter", mode = "markers+text",
      text = d$.id, textposition = "top center",
      textfont = list(size = 10),
      hovertext = paste0(
        "<b>", d$.id, "</b><br>",
        if (!is.na(group_col)) paste0(group_col, ": ", g, "<br>") else "",
        x, ": ", signif(d[[x]], 4), "<br>", y, ": ", signif(d[[y]], 4)
      ),
      hoverinfo = "text", name = g,
      marker = list(size = 11, color = cols[[g]], opacity = 0.85,
                    line = list(width = 1, color = "white"))
    )
  }
  p |>
    plotly::layout(
      title = list(text = title %||% paste(x, "vs", y), x = 0.02),
      xaxis = list(title = axis_title(x), zeroline = TRUE),
      yaxis = list(title = axis_title(y), zeroline = TRUE),
      legend = list(title = list(text = if (is.na(group_col)) "" else group_col))
    ) |>
    .rep_plotly_config("pca")
}

#' Interactive scree plot
#' @export
iplot_scree <- function(variance, max_pcs = 10L) {
  .rep_require("plotly")
  if (is.null(variance) || length(variance) == 0L) {
    return(.rep_empty_plot("Variance explained not available"))
  }
  v <- utils::head(variance, max_pcs)
  plotly::plot_ly(
    x = factor(names(v), levels = names(v)), y = as.numeric(v),
    type = "bar", marker = list(color = "#0072B2"),
    hovertemplate = "%{x}: %{y}%<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = "Variance explained", x = 0.02),
      xaxis = list(title = ""), yaxis = list(title = "% variance")
    ) |>
    .rep_plotly_config("scree")
}

# ---- DEG counts -----------------------------------------------------------------

#' Interactive DEG-count bar chart
#' @export
iplot_deg_counts <- function(
  deg_counts, colors = .de_direction_colors(),
  test_type = "Wald", use_lfc = TRUE
) {
  .rep_require("plotly")
  colors <- .rep_hex(colors)
  if (is.null(deg_counts) || nrow(deg_counts) == 0L ||
      !"padj_cutoff" %in% colnames(deg_counts)) {
    return(.rep_empty_plot("DEG counts not available"))
  }
  lab <- paste0("padj < ", deg_counts$padj_cutoff)
  lab <- factor(lab, levels = unique(lab))
  up_col <- if (use_lfc && "n_up_lfc" %in% colnames(deg_counts)) "n_up_lfc" else "n_up"
  dn_col <- if (use_lfc && "n_down_lfc" %in% colnames(deg_counts)) "n_down_lfc" else "n_down"
  tot_col <- if (use_lfc && "n_sig_lfc" %in% colnames(deg_counts)) "n_sig_lfc" else "n_sig"

  no_dir <- identical(test_type, "LRT") ||
    !all(c(up_col, dn_col) %in% colnames(deg_counts)) ||
    all(is.na(deg_counts[[up_col]]))

  p <- plotly::plot_ly()
  if (no_dir) {
    if (!tot_col %in% colnames(deg_counts)) {
      return(.rep_empty_plot("DEG counts not available"))
    }
    p <- plotly::add_bars(
      p, x = lab, y = deg_counts[[tot_col]], name = "Significant",
      marker = list(color = unname(colors[["up"]]))
    )
  } else {
    p <- plotly::add_bars(
      p, x = lab, y = deg_counts[[up_col]], name = "Upregulated",
      marker = list(color = unname(colors[["up"]]))
    )
    p <- plotly::add_bars(
      p, x = lab, y = deg_counts[[dn_col]], name = "Downregulated",
      marker = list(color = unname(colors[["down"]]))
    )
  }
  p |>
    plotly::layout(
      barmode = "group",
      title = list(
        text = paste0(
          "Differentially expressed genes",
          if (use_lfc && "log2FC_cutoff" %in% colnames(deg_counts))
            paste0(" (|log2FC| > ", deg_counts$log2FC_cutoff[[1L]], ")") else ""
        ),
        x = 0.02
      ),
      xaxis = list(title = ""), yaxis = list(title = "number of genes")
    ) |>
    .rep_plotly_config("deg_counts")
}

# ---- expression heatmap -----------------------------------------------------------

#' Interactive heatmap of normalised counts (row z-scores)
#' @export
iplot_heatmap <- function(
  mat, groups = NULL, colors = .de_direction_colors(),
  max_genes = 50L, title = NULL, cluster_columns = TRUE
) {
  .rep_require("plotly")
  colors <- .rep_hex(colors)
  if (is.null(mat) || nrow(mat) < 2L || ncol(mat) < 2L) {
    return(.rep_empty_plot("Normalised counts not available for the heatmap"))
  }
  mat <- mat[seq_len(min(nrow(mat), max_genes)), , drop = FALSE]
  z <- t(scale(t(log2(mat + 1))))
  z[!is.finite(z)] <- 0

  if (nrow(z) >= 3L) {
    ord <- tryCatch(stats::hclust(stats::dist(z))$order,
                    error = function(e) seq_len(nrow(z)))
    z <- z[ord, , drop = FALSE]
  }
  col_ord <- seq_len(ncol(z))
  if (isTRUE(cluster_columns) && ncol(z) >= 3L) {
    col_ord <- tryCatch(
      stats::hclust(stats::dist(t(z)))$order,
      error = function(e) seq_len(ncol(z))
    )
  } else if (!is.null(groups) && length(groups) == ncol(z)) {
    col_ord <- order(factor(groups, levels = unique(groups)), seq_len(ncol(z)))
  }
  z <- z[, col_ord, drop = FALSE]
  x_lab <- colnames(z)
  if (!is.null(groups) && length(groups) == length(col_ord)) {
    x_lab <- paste0(groups[col_ord], " · ", colnames(z))
  }

  zmax <- max(abs(z), na.rm = TRUE)
  plotly::plot_ly(
    x = x_lab, y = rownames(z), z = z, type = "heatmap",
    zmin = -zmax, zmax = zmax,
    colorscale = list(
      c(0, unname(colors[["down"]])), c(0.5, "#FFFFFF"),
      c(1, unname(colors[["up"]]))
    ),
    colorbar = list(title = "z-score"),
    hovertemplate = "%{y}<br>%{x}<br>z = %{z:.2f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(
        text = title %||% "Top differentially expressed genes (row z-score of log2 normalised counts)",
        x = 0.02
      ),
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "", autorange = "reversed", tickfont = list(size = 9)),
      margin = list(l = 110, b = 120)
    ) |>
    .rep_plotly_config("heatmap")
}

# ---- enrichment ---------------------------------------------------------------------

#' Interactive enrichment plot
#' @export
iplot_enrichment <- function(
  std, top_n = 20L, padj_cutoff = 0.05,
  colors = .de_direction_colors(), title = NULL
) {
  .rep_require("plotly")
  colors <- .rep_hex(colors)
  if (is.null(std) || nrow(std) == 0L) {
    return(.rep_empty_plot("No enrichment results"))
  }
  d <- utils::head(std[order(std$padj), , drop = FALSE], top_n)
  d$label <- make.unique(.rep_trunc(d$term, 55L))
  height <- max(320L, 24L * nrow(d) + 140L)
  is_gsea <- all(!is.na(d$score))

  fmt_hover <- function(d) {
    paste0(
      "<b>", d$term, "</b><br>adj. p: ", signif(d$padj, 3),
      ifelse(is.na(d$size), "", paste0("<br>size/count: ", d$size)),
      ifelse(is.na(d$score), "", paste0("<br>score: ", signif(d$score, 3)))
    )
  }

  if (is_gsea) {
    d <- d[order(d$score), , drop = FALSE]
    col <- ifelse(d$score >= 0, unname(colors[["up"]]), unname(colors[["down"]]))
    p <- plotly::plot_ly(
      x = d$score, y = factor(d$label, levels = d$label), type = "bar",
      orientation = "h", marker = list(color = col),
      hovertext = fmt_hover(d), hoverinfo = "text", height = height
    ) |>
      plotly::layout(xaxis = list(title = "normalised enrichment score (NES)"))
    ref_x <- 0
  } else {
    d <- d[order(-d$padj), , drop = FALSE]
    p <- plotly::plot_ly(
      x = -log10(pmax(d$padj, .Machine$double.xmin)),
      y = factor(d$label, levels = d$label), type = "bar", orientation = "h",
      marker = if (all(is.na(d$size))) {
        list(color = "#0072B2")
      } else {
        list(
          color = d$size, colorscale = "Viridis", reversescale = TRUE,
          showscale = TRUE, colorbar = list(title = "count")
        )
      },
      hovertext = fmt_hover(d), hoverinfo = "text", height = height
    ) |>
      plotly::layout(xaxis = list(title = "-log10(adjusted p-value)"))
    ref_x <- -log10(padj_cutoff)
  }
  p |>
    plotly::layout(
      title = list(text = title %||% "Top enriched terms", x = 0.02),
      yaxis = list(title = "", tickfont = list(size = 10)),
      margin = list(l = 330),
      shapes = list(list(
        type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = ref_x, x1 = ref_x,
        line = list(color = "#666666", width = 1, dash = "dot")
      ))
    ) |>
    .rep_plotly_config("enrichment")
}

# ---- isoform ---------------------------------------------------------------------

#' Interactive DTU bar chart: strongest differential transcript usage features
#' @export
iplot_dtu <- function(dtu, top_n = 30L, padj_cutoff = 0.05) {
  .rep_require("plotly")
  if (is.null(dtu) || nrow(dtu) == 0L) {
    return(.rep_empty_plot("No DTU results"))
  }
  pcol <- .rep_first_col(dtu, c("adj_pvalue", "padj", "FDR", "pvalue"))
  if (is.na(pcol)) {
    return(.rep_empty_plot("No p-value column in the DTU table"))
  }
  lab_col <- .rep_first_col(dtu, c("gene_label", "gene_symbol", "gene_id", "feature_id"))
  d <- dtu
  d$.p <- .rep_num(d[[pcol]])
  d <- d[!is.na(d$.p), , drop = FALSE]
  if (nrow(d) == 0L) {
    return(.rep_empty_plot("No finite p-values in the DTU table"))
  }
  d <- utils::head(d[order(d$.p), , drop = FALSE], top_n)
  d <- d[order(-d$.p), , drop = FALSE]
  labels <- if (is.na(lab_col)) as.character(seq_len(nrow(d)))
            else as.character(d[[lab_col]])
  labels <- make.unique(.rep_trunc(labels, 50L))
  height <- max(320L, 22L * nrow(d) + 140L)

  plotly::plot_ly(
    x = -log10(pmax(d$.p, .Machine$double.xmin)),
    y = factor(labels, levels = labels), type = "bar", orientation = "h",
    marker = list(color = "#009E73"),
    hovertext = paste0("<b>", labels, "</b><br>", pcol, ": ", signif(d$.p, 3)),
    hoverinfo = "text", height = height
  ) |>
    plotly::layout(
      title = list(text = paste0("Top DTU features (", pcol, ")"), x = 0.02),
      xaxis = list(title = paste0("-log10(", pcol, ")")),
      yaxis = list(title = "", tickfont = list(size = 10)),
      margin = list(l = 260),
      shapes = list(list(
        type = "line", y0 = 0, y1 = 1, yref = "paper",
        x0 = -log10(padj_cutoff), x1 = -log10(padj_cutoff),
        line = list(color = "#666666", width = 1, dash = "dot")
      ))
    ) |>
    .rep_plotly_config("dtu_top")
}

#' Interactive p-value histogram
#' @export
iplot_pvalue_hist <- function(p, title = "p-value distribution") {
  .rep_require("plotly")
  p <- .rep_num(p)
  p <- p[is.finite(p)]
  if (length(p) == 0L) {
    return(.rep_empty_plot("No p-values"))
  }
  plotly::plot_ly(
    x = p, type = "histogram", xbins = list(start = 0, end = 1, size = 0.05),
    marker = list(color = "#56B4E9", line = list(color = "white", width = 1))
  ) |>
    plotly::layout(
      title = list(text = title, x = 0.02),
      xaxis = list(title = "p-value", range = c(0, 1)),
      yaxis = list(title = "count")
    ) |>
    .rep_plotly_config("pvalue_hist")
}

# ---- DTE / DTU interactive diagnostics --------------------------------------------

#' Interactive transcript-level DTE volcano plot
#' @keywords internal
iplot_dte_volcano <- function(dte, padj_cutoff = 0.05, lfc_cutoff = 1,
                              label_col = NULL, max_points = 5000L) {
  .rep_require("plotly")
  if (is.null(dte) || !nrow(dte)) return(.rep_empty_plot("No DTE results"))
  lfc_col <- .rep_first_col(dte, c("log2FoldChange", "log2FC", "logFC", "lfc"))
  q_col <- .rep_first_col(dte, c("padj", "adj_pvalue", "FDR", "qvalue", "p.adjust"))
  if (is.na(lfc_col) || is.na(q_col)) {
    return(.rep_empty_plot("DTE volcano requires a transcript log2 fold-change and adjusted p-value"))
  }
  if (is.null(label_col) || !label_col %in% colnames(dte)) {
    label_col <- .rep_first_col(dte, c("gene_label", "gene_symbol", "gene", "transcript_id", "feature_id"))
  }
  label <- if (is.na(label_col)) as.character(seq_len(nrow(dte))) else as.character(dte[[label_col]])
  d <- data.frame(
    label = label,
    lfc = .rep_num(dte[[lfc_col]]),
    q = .rep_num(dte[[q_col]]),
    stringsAsFactors = FALSE
  )
  d <- d[is.finite(d$lfc) & is.finite(d$q) & d$q >= 0 & d$q <= 1, , drop = FALSE]
  if (!nrow(d)) return(.rep_empty_plot("No finite DTE fold-changes and adjusted p-values"))
  d$mlog10 <- -log10(pmax(d$q, .Machine$double.xmin))
  d$class <- "Not significant"
  d$class[d$q < padj_cutoff & d$lfc > lfc_cutoff] <- "Upregulated"
  d$class[d$q < padj_cutoff & d$lfc < -lfc_cutoff] <- "Downregulated"
  if (nrow(d) > max_points) {
    keep <- order(d$q, -abs(d$lfc))[seq_len(max_points)]
    d <- d[keep, , drop = FALSE]
  }
  palette <- c("Not significant" = "#999999", "Downregulated" = "#0072B2", "Upregulated" = "#D55E00")
  p <- plotly::plot_ly()
  for (nm in names(palette)) {
    dd <- d[d$class == nm, , drop = FALSE]
    if (!nrow(dd)) next
    p <- plotly::add_trace(
      p, data = dd, x = ~lfc, y = ~mlog10,
      type = "scattergl", mode = "markers", name = paste0(nm, " (", nrow(dd), ")"),
      text = paste0("<b>", dd$label, "</b><br>", lfc_col, ": ", signif(dd$lfc, 3),
                    "<br>", q_col, ": ", signif(dd$q, 3)),
      hoverinfo = "text",
      marker = list(color = unname(palette[[nm]]), size = 6, opacity = 0.72)
    )
  }
  p |>
    plotly::layout(
      title = list(text = "DTE volcano (transcript level)", x = 0.02),
      xaxis = list(title = "Transcript log2 fold change", zeroline = TRUE),
      yaxis = list(title = "-log10(adjusted p-value)"),
      shapes = list(
        list(type = "line", x0 = -lfc_cutoff, x1 = -lfc_cutoff, y0 = 0, y1 = 1,
             yref = "paper", line = list(color = "#666666", dash = "dot")),
        list(type = "line", x0 = lfc_cutoff, x1 = lfc_cutoff, y0 = 0, y1 = 1,
             yref = "paper", line = list(color = "#666666", dash = "dot")),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = -log10(padj_cutoff),
             y1 = -log10(padj_cutoff), line = list(color = "#666666", dash = "dot"))
      ),
      legend = list(orientation = "h", y = -0.18), hovermode = "closest"
    ) |>
    .rep_plotly_config("dte_volcano")
}

#' Interactive DTE MA plot
#' @keywords internal
iplot_dte_ma <- function(dte, padj_cutoff = 0.05, lfc_cutoff = 1,
                         label_col = NULL, max_points = 5000L) {
  .rep_require("plotly")
  if (is.null(dte) || !nrow(dte)) return(.rep_empty_plot("No DTE results"))
  lfc_col <- .rep_first_col(dte, c("log2FoldChange", "log2FC", "logFC", "lfc"))
  mean_col <- .rep_first_col(dte, c("baseMean", "base_mean", "mean", "mean_expression"))
  q_col <- .rep_first_col(dte, c("padj", "adj_pvalue", "FDR", "qvalue", "p.adjust"))
  if (is.na(lfc_col) || is.na(mean_col)) {
    return(.rep_empty_plot("DTE MA plot requires mean transcript abundance and log2 fold-change"))
  }
  if (is.null(label_col) || !label_col %in% colnames(dte)) {
    label_col <- .rep_first_col(dte, c("gene_label", "gene_symbol", "gene", "transcript_id", "feature_id"))
  }
  label <- if (is.na(label_col)) as.character(seq_len(nrow(dte))) else as.character(dte[[label_col]])
  q <- if (is.na(q_col)) rep(NA_real_, nrow(dte)) else .rep_num(dte[[q_col]])
  d <- data.frame(label = label, mean = .rep_num(dte[[mean_col]]),
                  lfc = .rep_num(dte[[lfc_col]]), q = q, stringsAsFactors = FALSE)
  d <- d[is.finite(d$mean) & d$mean > 0 & is.finite(d$lfc), , drop = FALSE]
  if (!nrow(d)) return(.rep_empty_plot("No finite positive DTE abundance values"))
  d$class <- ifelse(!is.na(d$q) & d$q < padj_cutoff & abs(d$lfc) > lfc_cutoff,
                    ifelse(d$lfc > 0, "Upregulated", "Downregulated"), "Not significant")
  if (nrow(d) > max_points) d <- d[order(d$q, -abs(d$lfc), na.last = TRUE)[seq_len(max_points)], , drop = FALSE]
  palette <- c("Not significant" = "#999999", "Downregulated" = "#0072B2", "Upregulated" = "#D55E00")
  p <- plotly::plot_ly()
  for (nm in names(palette)) {
    dd <- d[d$class == nm, , drop = FALSE]
    if (!nrow(dd)) next
    p <- plotly::add_trace(
      p, data = dd, x = ~mean, y = ~lfc, type = "scattergl", mode = "markers",
      name = paste0(nm, " (", nrow(dd), ")"),
      text = paste0("<b>", dd$label, "</b><br>mean: ", signif(dd$mean, 4),
                    "<br>log2FC: ", signif(dd$lfc, 3),
                    ifelse(is.na(dd$q), "", paste0("<br>adjusted p: ", signif(dd$q, 3)))),
      hoverinfo = "text", marker = list(color = unname(palette[[nm]]), size = 6, opacity = 0.72)
    )
  }
  p |>
    plotly::layout(
      title = list(text = "DTE MA plot", x = 0.02),
      xaxis = list(title = "Mean normalised transcript abundance", type = "log"),
      yaxis = list(title = "Transcript log2 fold change"),
      shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0,
                         line = list(color = "#666666", dash = "dot"))),
      legend = list(orientation = "h", y = -0.18)
    ) |>
    .rep_plotly_config("dte_ma")
}

#' Interactive ranked DTE transcript bar chart
#' @keywords internal
iplot_dte_ranked <- function(dte, top_n = 25L) {
  .rep_require("plotly")
  if (is.null(dte) || !nrow(dte)) return(.rep_empty_plot("No DTE results"))
  q_col <- .rep_first_col(dte, c("padj", "adj_pvalue", "FDR", "qvalue", "p.adjust"))
  id_col <- .rep_first_col(dte, c("gene_label", "gene_symbol", "gene", "transcript_id", "feature_id"))
  lfc_col <- .rep_first_col(dte, c("log2FoldChange", "log2FC", "logFC", "lfc"))
  if (is.na(q_col)) return(.rep_empty_plot("No adjusted p-value column in DTE table"))
  q <- .rep_num(dte[[q_col]])
  lfc <- if (is.na(lfc_col)) rep(NA_real_, nrow(dte)) else .rep_num(dte[[lfc_col]])
  label <- if (is.na(id_col)) as.character(seq_len(nrow(dte))) else as.character(dte[[id_col]])
  ok <- is.finite(q) & q >= 0 & q <= 1
  if (!any(ok)) return(.rep_empty_plot("No finite DTE adjusted p-values"))
  idx <- which(ok)[order(q[ok])]
  idx <- utils::head(idx, as.integer(top_n))
  idx <- rev(idx)
  labs <- make.unique(.rep_trunc(label[idx], 55L))
  cols <- ifelse(is.finite(lfc[idx]) & lfc[idx] < 0, "#0072B2", "#D55E00")
  plotly::plot_ly(
    x = -log10(pmax(q[idx], .Machine$double.xmin)),
    y = factor(labs, levels = labs), type = "bar", orientation = "h",
    marker = list(color = cols),
    text = paste0("<b>", labs, "</b><br>adjusted p: ", signif(q[idx], 3),
                  ifelse(is.finite(lfc[idx]), paste0("<br>log2FC: ", signif(lfc[idx], 3)), "")),
    hoverinfo = "text"
  ) |>
    plotly::layout(
      title = list(text = paste0("Top ", length(idx), " DTE transcripts by adjusted p-value"), x = 0.02),
      xaxis = list(title = "-log10(adjusted p-value)"),
      yaxis = list(title = "Transcript / gene", automargin = TRUE), margin = list(l = 260)
    ) |>
    .rep_plotly_config("dte_ranked")
}

#' Interactive DTU likelihood-ratio versus adjusted significance scatter
#' @keywords internal
iplot_dtu_lr_q <- function(dtu, padj_cutoff = 0.05, max_points = 5000L) {
  .rep_require("plotly")
  if (is.null(dtu) || !nrow(dtu)) return(.rep_empty_plot("No DTU results"))
  lr_col <- .rep_first_col(dtu, c("lr", "LR", "likelihood_ratio", "statistic"))
  q_col <- .rep_first_col(dtu, c("adj_pvalue", "padj", "FDR", "qvalue", "p.adjust"))
  id_col <- .rep_first_col(dtu, c("gene_label", "gene_symbol", "gene", "feature_id", "transcript_id"))
  if (is.na(lr_col) || is.na(q_col)) {
    return(.rep_empty_plot("DTU LR/significance view requires LR and adjusted p-value columns"))
  }
  label <- if (is.na(id_col)) as.character(seq_len(nrow(dtu))) else as.character(dtu[[id_col]])
  d <- data.frame(label = label, lr = .rep_num(dtu[[lr_col]]), q = .rep_num(dtu[[q_col]]),
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$lr) & is.finite(d$q) & d$q >= 0 & d$q <= 1, , drop = FALSE]
  if (!nrow(d)) return(.rep_empty_plot("No finite DTU LR and adjusted p-values"))
  d$mlog10 <- -log10(pmax(d$q, .Machine$double.xmin))
  d$class <- ifelse(d$q < padj_cutoff, "Adjusted p below cutoff", "Not significant")
  if (nrow(d) > max_points) d <- d[order(d$q)[seq_len(max_points)], , drop = FALSE]
  cols <- c("Not significant" = "#999999", "Adjusted p below cutoff" = "#7A6FAE")
  p <- plotly::plot_ly()
  for (nm in names(cols)) {
    dd <- d[d$class == nm, , drop = FALSE]
    if (!nrow(dd)) next
    p <- plotly::add_trace(
      p, data = dd, x = ~lr, y = ~mlog10, type = "scattergl", mode = "markers",
      name = paste0(nm, " (", nrow(dd), ")"),
      text = paste0("<b>", dd$label, "</b><br>LR: ", signif(dd$lr, 4),
                    "<br>adjusted p: ", signif(dd$q, 3)),
      hoverinfo = "text", marker = list(color = unname(cols[[nm]]), size = 6, opacity = 0.72)
    )
  }
  p |>
    plotly::layout(
      title = list(text = "DTU likelihood-ratio vs adjusted p-value", x = 0.02),
      xaxis = list(title = "Likelihood-ratio statistic"),
      yaxis = list(title = "-log10(adjusted p-value)"),
      shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper",
                         y0 = -log10(padj_cutoff), y1 = -log10(padj_cutoff),
                         line = list(color = "#666666", dash = "dot"))),
      legend = list(orientation = "h", y = -0.18)
    ) |>
    .rep_plotly_config("dtu_lr_q")
}

#' Interactive histogram of DTU likelihood-ratio statistics
#' @keywords internal
iplot_dtu_lr_hist <- function(dtu) {
  .rep_require("plotly")
  if (is.null(dtu) || !nrow(dtu)) return(.rep_empty_plot("No DTU results"))
  lr_col <- .rep_first_col(dtu, c("lr", "LR", "likelihood_ratio", "statistic"))
  if (is.na(lr_col)) return(.rep_empty_plot("No likelihood-ratio column in DTU table"))
  x <- .rep_num(dtu[[lr_col]])
  x <- x[is.finite(x)]
  if (!length(x)) return(.rep_empty_plot("No finite DTU likelihood-ratio statistics"))
  plotly::plot_ly(x = x, type = "histogram", nbinsx = 40,
                  marker = list(color = "#7A6FAE")) |>
    plotly::layout(
      title = list(text = "DTU likelihood-ratio distribution", x = 0.02),
      xaxis = list(title = "Likelihood-ratio statistic"),
      yaxis = list(title = "Feature count")
    ) |>
    .rep_plotly_config("dtu_lr_hist")
}

#' Interactive DTE transcript abundance versus adjusted significance
#' @keywords internal
iplot_dte_abundance_significance <- function(dte, padj_cutoff = 0.05, max_points = 5000L) {
  .rep_require("plotly")
  if (is.null(dte) || !nrow(dte)) return(.rep_empty_plot("No DTE results"))
  mean_col <- .rep_first_col(dte, c("baseMean", "base_mean", "mean", "mean_expression"))
  q_col <- .rep_first_col(dte, c("padj", "adj_pvalue", "FDR", "qvalue", "p.adjust"))
  id_col <- .rep_first_col(dte, c("gene_label", "gene_symbol", "gene", "transcript_id", "feature_id"))
  if (is.na(mean_col) || is.na(q_col)) {
    return(.rep_empty_plot("DTE abundance-significance view requires mean abundance and adjusted p-value"))
  }
  label <- if (is.na(id_col)) as.character(seq_len(nrow(dte))) else as.character(dte[[id_col]])
  d <- data.frame(label = label, mean = .rep_num(dte[[mean_col]]),
                  q = .rep_num(dte[[q_col]]), stringsAsFactors = FALSE)
  d <- d[is.finite(d$mean) & d$mean > 0 & is.finite(d$q) & d$q >= 0 & d$q <= 1, , drop = FALSE]
  if (!nrow(d)) return(.rep_empty_plot("No finite positive abundance and adjusted p-values"))
  d$mlog10 <- -log10(pmax(d$q, .Machine$double.xmin))
  d$class <- ifelse(d$q < padj_cutoff, "Adjusted p below cutoff", "Not significant")
  if (nrow(d) > max_points) d <- d[order(d$q, -d$mlog10)[seq_len(max_points)], , drop = FALSE]
  cols <- c("Not significant" = "#999999", "Adjusted p below cutoff" = "#7A6FAE")
  p <- plotly::plot_ly()
  for (nm in names(cols)) {
    dd <- d[d$class == nm, , drop = FALSE]
    if (!nrow(dd)) next
    p <- plotly::add_trace(
      p, data = dd, x = ~mean, y = ~mlog10, type = "scattergl", mode = "markers",
      name = paste0(nm, " (", nrow(dd), ")"),
      text = paste0("<b>", dd$label, "</b><br>mean abundance: ", signif(dd$mean, 4),
                    "<br>adjusted p: ", signif(dd$q, 3)),
      hoverinfo = "text", marker = list(color = unname(cols[[nm]]), size = 6, opacity = 0.72)
    )
  }
  p |>
    plotly::layout(
      title = list(text = "DTE transcript abundance vs adjusted significance", x = 0.02),
      xaxis = list(title = "Mean normalised transcript abundance", type = "log"),
      yaxis = list(title = "-log10(adjusted p-value)"),
      shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper",
                         y0 = -log10(padj_cutoff), y1 = -log10(padj_cutoff),
                         line = list(color = "#666666", dash = "dot"))),
      legend = list(orientation = "h", y = -0.18)
    ) |>
    .rep_plotly_config("dte_abundance_significance")
}

#' Interactive DTU likelihood-ratio versus degrees-of-freedom diagnostic
#' @keywords internal
iplot_dtu_lr_df <- function(dtu, max_points = 5000L) {
  .rep_require("plotly")
  if (is.null(dtu) || !nrow(dtu)) return(.rep_empty_plot("No DTU results"))
  lr_col <- .rep_first_col(dtu, c("lr", "LR", "likelihood_ratio", "statistic"))
  df_col <- .rep_first_col(dtu, c("df", "Df", "degrees_of_freedom"))
  id_col <- .rep_first_col(dtu, c("gene_label", "gene_symbol", "gene", "feature_id", "transcript_id"))
  if (is.na(lr_col) || is.na(df_col)) return(.rep_empty_plot("DTU LR-vs-df view requires LR and df columns"))
  label <- if (is.na(id_col)) as.character(seq_len(nrow(dtu))) else as.character(dtu[[id_col]])
  d <- data.frame(label = label, lr = .rep_num(dtu[[lr_col]]), df = .rep_num(dtu[[df_col]]),
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$lr) & is.finite(d$df), , drop = FALSE]
  if (!nrow(d)) return(.rep_empty_plot("No finite DTU LR and df values"))
  if (nrow(d) > max_points) d <- d[order(d$lr, decreasing = TRUE)[seq_len(max_points)], , drop = FALSE]
  plotly::plot_ly(
    d, x = ~df, y = ~lr, type = "scattergl", mode = "markers",
    text = ~paste0("<b>", label, "</b><br>df: ", df, "<br>LR: ", signif(lr, 4)),
    hoverinfo = "text", marker = list(color = "#4C78A8", size = 6, opacity = 0.7)
  ) |>
    plotly::layout(
      title = list(text = "DTU likelihood-ratio vs degrees of freedom", x = 0.02),
      xaxis = list(title = "Degrees of freedom"), yaxis = list(title = "Likelihood-ratio statistic")
    ) |>
    .rep_plotly_config("dtu_lr_df")
}

#' Interactive DTU feature counts by gene (descriptive, not gene-level FDR)
#' @keywords internal
iplot_dtu_gene_summary <- function(dtu, top_n = 30L) {
  .rep_require("plotly")
  if (is.null(dtu) || !nrow(dtu)) return(.rep_empty_plot("No DTU results"))
  gene_col <- .rep_first_col(dtu, c("gene_symbol", "gene", "gene_id", "gene_label"))
  q_col <- .rep_first_col(dtu, c("adj_pvalue", "padj", "FDR", "qvalue", "p.adjust"))
  if (is.na(gene_col)) return(.rep_empty_plot("DTU gene summary requires a gene identifier"))
  gene <- as.character(dtu[[gene_col]])
  gene[is.na(gene) | !nzchar(gene)] <- "Unannotated"
  groups <- split(seq_along(gene), gene)
  q <- if (is.na(q_col)) rep(NA_real_, nrow(dtu)) else .rep_num(dtu[[q_col]])
  summary <- data.frame(
    gene = names(groups),
    n_features = vapply(groups, length, integer(1)),
    best_feature_q = vapply(groups, function(i) {
      z <- q[i][is.finite(q[i])]
      if (!length(z)) NA_real_ else min(z)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  summary <- summary[order(summary$n_features, summary$best_feature_q, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  summary <- utils::head(summary, as.integer(top_n))
  summary$gene <- factor(summary$gene, levels = rev(summary$gene))
  plotly::plot_ly(
    summary, x = ~n_features, y = ~gene, type = "bar", orientation = "h",
    marker = list(color = "#7A6FAE"),
    text = ~paste0("<b>", gene, "</b><br>DTU features: ", n_features,
                   ifelse(is.na(best_feature_q), "", paste0("<br>best feature-level adjusted p: ", signif(best_feature_q, 3)))),
    hoverinfo = "text"
  ) |>
    plotly::layout(
      title = list(text = "Genes with the most tested DTU features", x = 0.02),
      xaxis = list(title = "Number of DTU features"), yaxis = list(title = "Gene", automargin = TRUE),
      margin = list(l = 220)
    ) |>
    .rep_plotly_config("dtu_gene_summary")
}

# ---- isoform switch ---------------------------------------------------------------

#' Interactive dIF vs q-value scatter of isoform switches
#' @export
iplot_isoform_switch <- function(
  feat, padj_cutoff = 0.05, dIF_cutoff = 0.1, max_points = 5000L,
  positive_col = "#c94a4a", neutral_col = "#999999"
) {
  .rep_require("plotly")
  if (is.null(feat) || !nrow(feat)) {
    return(.rep_empty_plot("No isoform-switch features"))
  }
  if (!all(c("dIF", "isoform_switch_q_value") %in% colnames(feat))) {
    return(.rep_empty_plot("dIF / isoform_switch_q_value not found"))
  }
  d <- feat
  d$dIF <- .rep_num(d$dIF)
  d$.q  <- .rep_num(d$isoform_switch_q_value)
  d <- d[is.finite(d$dIF) & !is.na(d$.q), , drop = FALSE]
  if (!nrow(d)) return(.rep_empty_plot("No finite switch statistics"))
  d$.mlog10 <- -log10(pmax(d$.q, .Machine$double.xmin))
  d$.sig <- d$.q < padj_cutoff & abs(d$dIF) > dIF_cutoff
  d$.label <- if ("gene_name" %in% colnames(d)) d$gene_name
              else if ("gene_id" %in% colnames(d)) d$gene_id
              else rep("", nrow(d))

  if (nrow(d) > max_points) {
    keep_sig  <- which(d$.sig)
    keep_rest <- setdiff(seq_len(nrow(d)), keep_sig)
    n_rest <- max(0L, max_points - length(keep_sig))
    if (length(keep_rest) > n_rest) {
      keep_rest <- if (requireNamespace("withr", quietly = TRUE)) {
        withr::with_seed(1L, sample(keep_rest, n_rest))
      } else {
        sample(keep_rest, n_rest)
      }
    }
    d <- d[sort(c(keep_sig, keep_rest)), , drop = FALSE]
  }

  p <- plotly::plot_ly()
  for (k in c(FALSE, TRUE)) {
    dd <- d[d$.sig == k, , drop = FALSE]
    if (!nrow(dd)) next
    p <- plotly::add_trace(
      p, data = dd, x = ~dIF, y = ~.mlog10,
      type = "scatter", mode = "markers",
      name = if (k) paste0("Switch (", nrow(dd), ")")
             else paste0("Not significant (", nrow(dd), ")"),
      marker = list(
        color   = if (k) positive_col else neutral_col,
        size    = if (k) 7 else 4,
        opacity = if (k) 0.85 else 0.45
      ),
      text = paste0("<b>", dd$.label, "</b><br>dIF: ", signif(dd$dIF, 3),
                    "<br>q: ", signif(dd$.q, 3)),
      hoverinfo = "text"
    )
  }
  p |>
    plotly::layout(
      title = list(text = "Isoform switch overview", x = 0.02),
      xaxis = list(title = "dIF (isoform fraction change)", zeroline = TRUE),
      yaxis = list(title = "-log10(q)"),
      shapes = list(
        .rep_vline(-dIF_cutoff), .rep_vline(dIF_cutoff),
        .rep_hline(-log10(padj_cutoff))
      ),
      legend = list(orientation = "h", y = -0.18)
    ) |>
    .rep_plotly_config("isoform_switch")
}

# ---- per-gene transcript proportions ----------------------------------------------

#' Interactive stacked bar of transcript proportions per sample for one gene
#' @export
iplot_transcript_proportions <- function(isoform_obj, gene_symbol,
                                        condition = "condition") {
  .rep_require("plotly")
  if (is.null(isoform_obj) || is.null(gene_symbol) || !nzchar(gene_symbol)) {
    return(.rep_empty_plot("No gene requested"))
  }
  counts <- if (identical(isoform_obj$type, "tximport")) isoform_obj$txi$counts
            else isoform_obj$counts
  gene_ids <- isoform_obj$gene_map$ensembl[
    !is.na(isoform_obj$gene_map$symbol) & isoform_obj$gene_map$symbol == gene_symbol
  ]
  if (!length(gene_ids)) {
    return(.rep_empty_plot(paste0("Gene not found: ", gene_symbol)))
  }
  tx_ids <- isoform_obj$tx2gene$tx_id[isoform_obj$tx2gene$gene_id %in% gene_ids]
  tx_ids <- intersect(tx_ids, rownames(counts))
  if (length(tx_ids) < 2L) {
    return(.rep_empty_plot("Fewer than 2 quantified transcripts"))
  }
  sub    <- counts[tx_ids, , drop = FALSE]
  totals <- colSums(sub)
  props  <- sweep(sub, 2, ifelse(totals > 0, totals, NA_real_), "/")
  props[is.na(props)] <- 0

  long <- data.frame(
    sample     = rep(colnames(props), each = nrow(props)),
    transcript = rep(rownames(props), times = ncol(props)),
    proportion = as.vector(props),
    stringsAsFactors = FALSE
  )
  if (condition %in% colnames(isoform_obj$meta)) {
    long$group <- as.character(
      isoform_obj$meta[[condition]][match(long$sample, rownames(isoform_obj$meta))]
    )
  } else {
    long$group <- "All"
  }
  long$group  <- factor(long$group, levels = unique(long$group))
  long$sample <- factor(long$sample, levels = colnames(props))

  plotly::plot_ly(
    x = long$sample, y = long$proportion, color = long$transcript,
    type = "bar",
    hovertext = paste0(long$sample, "<br>", long$transcript, ": ",
                       scales::percent(long$proportion, accuracy = 0.1)),
    hoverinfo = "text"
  ) |>
    plotly::layout(
      barmode = "stack",
      title = list(text = paste0("Transcript proportions – ", gene_symbol), x = 0.02),
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "proportion", tickformat = ".0%"),
      legend = list(title = list(text = "transcript"))
    ) |>
    .rep_plotly_config("transcript_proportions")
}

# ---- decoupleR activity ------------------------------------------------------------

#' Interactive barplot of decoupleR inferred TF/pathway activities
#' @export
iplot_decouple_activity <- function(
  activity_df, top_n = 25L, padj_cutoff = 0.05,
  title = "TF activity (decoupleR)",
  positive_col = "#c94a4a", negative_col = "#2b5c8f"
) {
  .rep_require("plotly")
  if (is.null(activity_df) || !nrow(activity_df)) {
    return(.rep_empty_plot("No decoupleR activity results"))
  }
  pcol <- .rep_first_col(activity_df, c("p_value", "padj", "pval", "p.value"))
  d <- activity_df
  if (!is.na(pcol)) {
    d$.p <- .rep_num(d[[pcol]])
    d <- d[!is.na(d$.p) & d$.p < padj_cutoff, , drop = FALSE]
  } else {
    d$.p <- NA_real_
  }
  if (!nrow(d)) {
    return(.rep_empty_plot(paste0("No source with p < ", padj_cutoff)))
  }
  d <- utils::head(d[order(-abs(d$score)), , drop = FALSE], top_n)
  d <- d[order(d$score), , drop = FALSE]
  d$label <- make.unique(.rep_trunc(d$source, 40L))
  cols <- ifelse(d$score > 0, positive_col, negative_col)
  height <- max(320L, 22L * nrow(d) + 140L)

  plotly::plot_ly(
    x = d$score, y = factor(d$label, levels = d$label),
    type = "bar", orientation = "h", height = height,
    marker = list(color = cols),
    hovertext = paste0(
      "<b>", d$source, "</b><br>score: ", signif(d$score, 3),
      ifelse(is.na(d$.p), "", paste0("<br>p: ", signif(d$.p, 3)))
    ),
    hoverinfo = "text"
  ) |>
    plotly::layout(
      title = list(text = title, x = 0.02),
      xaxis = list(title = "activity score (t-like)"),
      yaxis = list(title = "", tickfont = list(size = 10)),
      margin = list(l = 240),
      shapes = list(.rep_vline(0, "#999999"))
    ) |>
    .rep_plotly_config("decouple")
}

# ---- TF regulon interactive plots -------------------------------------------------

#' Interactive targeted volcano for one TF's regulon
#'
#' Points are colored by whether the gene behaves as predicted by the network
#' (mor direction), exactly as the static PDF: "Activated (expected)",
#' "Repressed (expected)", "Other change", "Not significant".
#'
#' @param regulon_df Regulon table (`target`, `log2FoldChange`, `padj`, `mor`,
#'   `Regulation`) written by `.run_decouple_tf_regulon_plots()`.
#' @param tf          TF name (for the title).
#' @param padj_cutoff,lfc_cutoff Thresholds drawn as dashed lines.
#' @export
iplot_tf_regulon_volcano <- function(
  regulon_df, tf, padj_cutoff = 0.05, lfc_cutoff = 1
) {
  .rep_require("plotly")
  if (is.null(regulon_df) || !nrow(regulon_df)) {
    return(.rep_empty_plot(paste("No regulon targets for", tf)))
  }
  d <- regulon_df
  d$.mlog10 <- -log10(pmax(.rep_num(d$padj), .Machine$double.xmin))
  d$.lfc    <- .rep_num(d$log2FoldChange)
  d$.mor    <- .rep_num(d$mor)

  cols <- c(
    "Not significant"      = "#BDBDBD",
    "Other change"         = "#ff7f00",
    "Activated (expected)" = "#1f78b4",
    "Repressed (expected)" = "#e31a1c"
  )
  ord <- c("Not significant", "Other change",
           "Activated (expected)", "Repressed (expected)")

  p <- plotly::plot_ly()
  for (cat in ord) {
    dd <- d[as.character(d$Regulation) == cat, , drop = FALSE]
    if (!nrow(dd)) next
    p <- plotly::add_trace(
      p, data = dd, x = ~.lfc, y = ~.mlog10,
      type = "scatter", mode = "markers",
      name = paste0(cat, " (", nrow(dd), ")"),
      marker = list(
        color   = cols[[cat]],
        size    = if (cat == "Not significant") 6 else 9,
        opacity = if (cat == "Not significant") 0.55 else 0.9
      ),
      text = paste0(
        "<b>", dd$target, "</b><br>log2FC: ", signif(dd$.lfc, 3),
        "<br>padj: ",   signif(.rep_num(dd$padj), 3),
        "<br>mor: ",    signif(dd$.mor, 3)
      ),
      hoverinfo = "text"
    )
  }
  shapes <- list(
    .rep_vline(-lfc_cutoff), .rep_vline(lfc_cutoff),
    .rep_hline(-log10(padj_cutoff))
  )
  p |>
    plotly::layout(
      title = list(text = paste("Regulon:", tf), x = 0.02),
      xaxis = list(title = "log2 fold change", zeroline = TRUE),
      yaxis = list(title = "-log10(padj)"),
      shapes = shapes,
      legend = list(orientation = "h", y = -0.18)
    ) |>
    .rep_plotly_config("tf_regulon_volcano")
}

#' Interactive leading-edge plot for one TF's regulon
#'
#' Bar chart of the running enrichment score over all genes ranked by stat,
#' with vertical markers at the rank of each of the TF's targets. The
#' running ES curve is subsampled to ~1000 points (keeps the JSON payload
#' light) while all targets are kept as markers.
#'
#' @param le_df     Leading-edge table (`rank`, `target`, `stat`) written by
#'   `.run_decouple_tf_regulon_plots()`.
#' @param dge_results Full DGE results table (`gene`, `stat` / `contrast_stat`).
#' @param tf        TF name (for the title).
#' @param test_type "Wald" or "LRT" -- determines which stat column ranks.
#' @export
iplot_tf_leading_edge <- function(le_df, dge_results, tf, test_type = "Wald") {
  .rep_require("plotly")
  if (is.null(le_df) || !nrow(le_df) || is.null(dge_results) || !nrow(dge_results)) {
    return(.rep_empty_plot(paste("No leading-edge data for", tf)))
  }
  is_lrt <- identical(test_type, "LRT") && "contrast_stat" %in% colnames(dge_results)
  stat_col <- if (is_lrt) "contrast_stat" else "stat"
  if (!stat_col %in% colnames(dge_results)) {
    return(.rep_empty_plot("No stat column available"))
  }

  keep <- !is.na(dge_results[[stat_col]]) &
          !is.na(dge_results$gene) & dge_results$gene != ""
  r <- dge_results[keep, c("gene", stat_col), drop = FALSE]
  r <- r[order(-r[[stat_col]]), , drop = FALSE]
  r <- r[!duplicated(r$gene), , drop = FALSE]
  n_genes <- nrow(r)
  if (!n_genes) return(.rep_empty_plot("No ranked genes"))

  positions <- le_df$rank[order(le_df$rank)]
  positions <- positions[!is.na(positions) & positions >= 1L & positions <= n_genes]
  n_hits <- length(positions)
  if (!n_hits) return(.rep_empty_plot("No hits within rank range"))

  is_hit <- rep(FALSE, n_genes)
  is_hit[positions] <- TRUE
  hit_w  <- 1 / n_hits
  miss_w <- -1 / max(1L, n_genes - n_hits)
  es_vec <- cumsum(ifelse(is_hit, hit_w, miss_w))
  peak_idx <- which.max(abs(es_vec))

  # Downsample curve to <= 1000 points
  idx <- if (n_genes > 1000) {
    round(seq(1, n_genes, length.out = 1000))
  } else seq_len(n_genes)

  plotly::plot_ly() |>
    plotly::add_trace(
      x = idx, y = es_vec[idx],
      type = "scatter", mode = "lines",
      name = "Running ES",
      line = list(color = "#1f78b4", width = 2),
      hoverinfo = "skip"
    ) |>
    plotly::add_trace(
      x = positions,
      y = rep(min(es_vec) - 0.02, n_hits),
      type = "scatter", mode = "markers",
      name = paste0("Targets (", n_hits, ")"),
      marker = list(
        color = "#e31a1c", symbol = "line-ns",
        size = 12, line = list(width = 2)
      ),
      text = paste0(
        "<b>", r$gene[positions], "</b><br>rank: ", positions,
        "<br>stat: ", signif(r[[stat_col]][positions], 3)
      ),
      hoverinfo = "text"
    ) |>
    plotly::layout(
      title = list(text = paste("Leading edge:", tf), x = 0.02),
      xaxis = list(title = "Rank (by stat, decreasing)"),
      yaxis = list(title = "Running enrichment score"),
      shapes = list(.rep_vline(peak_idx, "#888888")),
      legend = list(orientation = "h", y = -0.18)
    ) |>
    .rep_plotly_config("tf_leading_edge")
}


# ---- additional sample-level QC and differential diagnostics -----------------

.rep_matrix_for_qc <- function(mat, min_rows = 2L) {
  if (is.null(mat)) return(NULL)
  mat <- as.matrix(mat)
  if (nrow(mat) < min_rows || ncol(mat) < 2L) return(NULL)
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- NA_real_
  mat
}

#' Interactive sample-to-sample correlation heatmap
#' @keywords internal
.rep_iplot_sample_correlation <- function(mat, method = c("pearson", "spearman")) {
  .rep_require("plotly")
  method <- match.arg(method)
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat)) return(.rep_empty_plot("At least two samples and two genes are needed"))
  x <- log2(mat + 1)
  cm <- suppressWarnings(stats::cor(x, method = method, use = "pairwise.complete.obs"))
  cm[!is.finite(cm)] <- 0
  diag(cm) <- 1
  if (nrow(cm) >= 3L) {
    ord <- tryCatch(stats::hclust(stats::as.dist(1 - cm))$order,
                    error = function(e) seq_len(nrow(cm)))
    cm <- cm[ord, ord, drop = FALSE]
  }
  plotly::plot_ly(
    x = colnames(cm), y = rownames(cm), z = cm, type = "heatmap",
    zmin = -1, zmax = 1,
    colorscale = list(c(0, "#2b5c8f"), c(0.5, "#f7f7f7"), c(1, "#c94a4a")),
    colorbar = list(title = "r"),
    hovertemplate = "%{y} vs %{x}<br>correlation = %{z:.3f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = paste("Sample correlation (", method, ")", sep = ""), x = 0.02),
      xaxis = list(title = "", tickangle = -35),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 100, b = 100)
    ) |>
    .rep_plotly_config("sample_correlation")
}

#' Interactive sample-to-sample Euclidean distance heatmap
#' @keywords internal
.rep_iplot_sample_distance <- function(mat) {
  .rep_require("plotly")
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat)) return(.rep_empty_plot("Sample distance matrix not available"))
  x <- t(log2(mat + 1))
  d <- as.matrix(stats::dist(x, method = "euclidean"))
  if (nrow(d) >= 3L) {
    ord <- tryCatch(stats::hclust(stats::as.dist(d))$order,
                    error = function(e) seq_len(nrow(d)))
    d <- d[ord, ord, drop = FALSE]
  }
  plotly::plot_ly(
    x = colnames(d), y = rownames(d), z = d, type = "heatmap",
    colorscale = "Viridis", colorbar = list(title = "distance"),
    hovertemplate = "%{y} vs %{x}<br>Euclidean distance = %{z:.2f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = "Sample-to-sample distance", x = 0.02),
      xaxis = list(title = "", tickangle = -35),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 100, b = 100)
    ) |>
    .rep_plotly_config("sample_distance")
}

#' Interactive hierarchical sample dendrogram
#' @keywords internal
.rep_iplot_sample_dendrogram <- function(mat, method = "complete") {
  .rep_require("plotly")
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat) || ncol(mat) < 3L) {
    return(.rep_empty_plot("At least three samples are needed for a dendrogram"))
  }
  x <- t(log2(mat + 1))
  hc <- tryCatch(stats::hclust(stats::dist(x), method = method),
                 error = function(e) NULL)
  if (is.null(hc)) return(.rep_empty_plot("Hierarchical clustering could not be computed"))

  n <- length(hc$order)
  xpos <- rep(NA_real_, 2L * n - 1L)
  ypos <- rep(0, 2L * n - 1L)
  for (i in seq_len(n)) xpos[i] <- match(i, hc$order)
  seg_x <- seg_y <- numeric(0)
  for (k in seq_len(n - 1L)) {
    child_pos <- function(v) if (v < 0L) -v else n + v
    a <- child_pos(hc$merge[k, 1L])
    b <- child_pos(hc$merge[k, 2L])
    parent <- n + k
    xpos[parent] <- mean(c(xpos[a], xpos[b]))
    ypos[parent] <- hc$height[k]
    seg_x <- c(seg_x, xpos[a], xpos[a], NA, xpos[b], xpos[b], NA, xpos[a], xpos[b], NA)
    seg_y <- c(seg_y, ypos[a], hc$height[k], NA, ypos[b], hc$height[k], NA,
               hc$height[k], hc$height[k], NA)
  }
  labels <- colnames(mat)[hc$order]
  plotly::plot_ly() |>
    plotly::add_trace(
      x = seg_x, y = seg_y, type = "scatter", mode = "lines",
      line = list(color = "#496a81", width = 2), hoverinfo = "skip",
      showlegend = FALSE
    ) |>
    plotly::layout(
      title = list(text = paste("Hierarchical clustering (", method, ")", sep = ""), x = 0.02),
      xaxis = list(title = "", tickmode = "array", tickvals = seq_len(n),
                   ticktext = labels, tickangle = -45, range = c(0.5, n + 0.5)),
      yaxis = list(title = "Euclidean distance", rangemode = "tozero"),
      margin = list(l = 60, b = 120)
    ) |>
    .rep_plotly_config("sample_dendrogram")
}

#' Interactive per-sample boxplots of normalised expression
#' @keywords internal
.rep_iplot_expression_boxplot <- function(mat) {
  .rep_require("plotly")
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat)) return(.rep_empty_plot("Expression matrix not available"))
  long <- data.frame(
    sample = rep(colnames(mat), each = nrow(mat)),
    expression = as.vector(log2(mat + 1)),
    stringsAsFactors = FALSE
  )
  long <- long[is.finite(long$expression), , drop = FALSE]
  plotly::plot_ly(
    long, x = ~sample, y = ~expression, type = "box", color = ~sample,
    boxpoints = FALSE, jitter = 0.25, pointpos = 0,
    hovertemplate = "%{x}<br>log2(normalised count + 1) = %{y:.2f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = "Expression distribution by sample", x = 0.02),
      showlegend = FALSE,
      xaxis = list(title = "Sample", tickangle = -35),
      yaxis = list(title = "log2(normalised count + 1)")
    ) |>
    .rep_plotly_config("expression_boxplot")
}

#' Interactive per-sample expression density curves
#' @keywords internal
.rep_iplot_expression_density <- function(mat) {
  .rep_require("plotly")
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat)) return(.rep_empty_plot("Expression matrix not available"))
  p <- plotly::plot_ly()
  for (sample in colnames(mat)) {
    v <- log2(mat[, sample] + 1)
    v <- v[is.finite(v)]
    if (length(v) < 2L) next
    de <- stats::density(v, na.rm = TRUE, n = 256L)
    p <- plotly::add_lines(
      p, x = de$x, y = de$y, name = sample,
      hovertemplate = paste0(sample, "<br>expression = %{x:.2f}<br>density = %{y:.3f}<extra></extra>")
    )
  }
  p |>
    plotly::layout(
      title = list(text = "Expression density by sample", x = 0.02),
      xaxis = list(title = "log2(normalised count + 1)"),
      yaxis = list(title = "Density"),
      legend = list(orientation = "h", y = -0.2)
    ) |>
    .rep_plotly_config("expression_density")
}

#' Interactive multidimensional scaling of samples
#' @keywords internal
.rep_iplot_sample_mds <- function(mat) {
  .rep_require("plotly")
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat) || ncol(mat) < 3L) {
    return(.rep_empty_plot("At least three samples are needed for MDS"))
  }
  x <- t(log2(mat + 1))
  d <- stats::dist(x, method = "euclidean")
  coords <- tryCatch(stats::cmdscale(d, k = 2L), error = function(e) NULL)
  if (is.null(coords) || ncol(coords) < 2L) {
    return(.rep_empty_plot("MDS coordinates could not be computed"))
  }
  labs <- rownames(coords) %||% colnames(mat)
  plotly::plot_ly(
    x = coords[, 1L], y = coords[, 2L], type = "scatter", mode = "markers+text",
    text = labs, textposition = "top center",
    marker = list(size = 11, color = "#159a8c", line = list(color = "#ffffff", width = 1)),
    hovertemplate = "%{text}<br>MDS1 = %{x:.3f}<br>MDS2 = %{y:.3f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = "Multidimensional scaling of samples", x = 0.02),
      xaxis = list(title = "MDS1"), yaxis = list(title = "MDS2")
    ) |>
    .rep_plotly_config("sample_mds")
}

#' Interactive relative log expression (RLE) boxplots
#' @keywords internal
.rep_iplot_rle <- function(mat) {
  .rep_require("plotly")
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat)) return(.rep_empty_plot("Expression matrix not available for RLE"))
  x <- log2(mat + 1)
  gene_medians <- apply(x, 1L, stats::median, na.rm = TRUE)
  rle <- sweep(x, 1L, gene_medians, FUN = "-")
  long <- data.frame(
    sample = rep(colnames(rle), each = nrow(rle)),
    rle = as.vector(rle), stringsAsFactors = FALSE
  )
  long <- long[is.finite(long$rle), , drop = FALSE]
  plotly::plot_ly(
    long, x = ~sample, y = ~rle, type = "box", color = ~sample,
    boxpoints = FALSE,
    hovertemplate = "%{x}<br>relative log expression = %{y:.3f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = "Relative log expression by sample", x = 0.02),
      showlegend = FALSE,
      xaxis = list(title = "Sample", tickangle = -35),
      yaxis = list(title = "log2 count + 1 minus gene median"),
      shapes = list(.rep_hline(0, "#777777"))
    ) |>
    .rep_plotly_config("expression_rle")
}

#' Interactive mean-variance diagnostic
#' @keywords internal
.rep_iplot_mean_variance <- function(mat) {
  .rep_require("plotly")
  mat <- .rep_matrix_for_qc(mat)
  if (is.null(mat)) return(.rep_empty_plot("Expression matrix not available for mean-variance plot"))
  means <- rowMeans(mat, na.rm = TRUE)
  vars <- apply(mat, 1L, stats::var, na.rm = TRUE)
  keep <- is.finite(means) & is.finite(vars)
  if (!any(keep)) return(.rep_empty_plot("No finite mean-variance pairs"))
  labs <- rownames(mat)
  if (is.null(labs)) labs <- paste0("feature_", seq_len(nrow(mat)))
  plotly::plot_ly(
    x = log10(means[keep] + 1), y = log10(vars[keep] + 1),
    type = "scattergl", mode = "markers",
    text = labs[keep],
    marker = list(size = 5, color = "#8b6fc1", opacity = 0.6),
    hovertemplate = "%{text}<br>log10(mean + 1) = %{x:.2f}<br>log10(variance + 1) = %{y:.2f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = "Mean-variance relationship", x = 0.02),
      xaxis = list(title = "log10(mean normalised count + 1)"),
      yaxis = list(title = "log10(variance + 1)")
    ) |>
    .rep_plotly_config("mean_variance")
}

#' Interactive p-value Q-Q plot
#' @keywords internal
.rep_iplot_pvalue_qq <- function(pvalues) {
  .rep_require("plotly")
  p <- .rep_num(pvalues)
  p <- sort(p[is.finite(p) & p >= 0 & p <= 1])
  if (length(p) < 3L) return(.rep_empty_plot("At least three finite p-values are needed"))
  n <- length(p)
  expected <- -log10((seq_len(n) - 0.5) / n)
  observed <- -log10(pmax(p, .Machine$double.xmin))
  lim <- max(c(expected, observed), na.rm = TRUE)
  plotly::plot_ly(x = expected, y = observed, type = "scattergl", mode = "markers",
                  marker = list(size = 5, color = "#0072B2", opacity = 0.55),
                  hovertemplate = "Expected = %{x:.2f}<br>Observed = %{y:.2f}<extra></extra>") |>
    plotly::add_lines(x = c(0, lim), y = c(0, lim), name = "Null expectation",
                      line = list(color = "#d55e00", dash = "dash"), hoverinfo = "skip") |>
    plotly::layout(
      title = list(text = "P-value Q-Q plot", x = 0.02),
      xaxis = list(title = "Expected -log10(p)"),
      yaxis = list(title = "Observed -log10(p)")
    ) |>
    .rep_plotly_config("pvalue_qq")
}

#' Interactive ranked absolute fold changes
#' @keywords internal
.rep_iplot_ranked_effects <- function(res_tbl, top_n = 25L) {
  .rep_require("plotly")
  if (is.null(res_tbl) || !nrow(res_tbl) ||
      !all(c("gene", "log2FoldChange") %in% colnames(res_tbl))) {
    return(.rep_empty_plot("Fold-change results are not available"))
  }
  d <- res_tbl
  d$log2FoldChange <- .rep_num(d$log2FoldChange)
  d$padj <- if ("padj" %in% colnames(d)) .rep_num(d$padj) else NA_real_
  d <- d[is.finite(d$log2FoldChange), , drop = FALSE]
  d <- utils::head(d[order(-abs(d$log2FoldChange)), , drop = FALSE], top_n)
  d <- d[order(d$log2FoldChange), , drop = FALSE]
  d$direction <- ifelse(d$log2FoldChange > 0, "Positive", "Negative")
  plotly::plot_ly(
    d, x = ~log2FoldChange, y = ~factor(gene, levels = gene),
    type = "bar", orientation = "h", color = ~direction,
    colors = c("Negative" = "#2b5c8f", "Positive" = "#c94a4a"),
    hovertemplate = paste0("%{y}<br>log2FC = %{x:.2f}<extra></extra>")
  ) |>
    plotly::layout(
      title = list(text = paste0("Top ", nrow(d), " genes by |log2FC|"), x = 0.02),
      xaxis = list(title = "log2 fold change", zeroline = TRUE),
      yaxis = list(title = ""),
      margin = list(l = 140),
      shapes = list(.rep_vline(0, "#777777"))
    ) |>
    .rep_plotly_config("ranked_effects")
}

#' Interactive enrichment bubble plot
#' @keywords internal
.rep_iplot_enrichment_bubble <- function(std, top_n = 20L, title = "Enrichment overview") {
  .rep_require("plotly")
  if (is.null(std) || !nrow(std)) return(.rep_empty_plot("No enrichment results"))
  d <- utils::head(std[order(std$padj), , drop = FALSE], top_n)
  d$label <- make.unique(.rep_trunc(d$term, 52L))
  d$mlog10 <- -log10(pmax(d$padj, .Machine$double.xmin))
  is_gsea <- all(is.finite(d$score))
  d$x <- if (is_gsea) d$score else ifelse(is.finite(d$size), d$size, seq_len(nrow(d)))
  d$size_plot <- ifelse(is.finite(d$size), pmax(6, pmin(28, sqrt(d$size) * 3)), 10)
  xtitle <- if (is_gsea) "Normalised enrichment score (NES)" else "Gene count / overlap size"
  plotly::plot_ly(
    d, x = ~x, y = ~factor(label, levels = label),
    type = "scatter", mode = "markers",
    marker = list(size = ~size_plot, color = ~mlog10, colorscale = "Viridis",
                  showscale = TRUE, colorbar = list(title = "-log10(FDR)"),
                  sizemode = "diameter", line = list(color = "#ffffff", width = 1)),
    text = ~paste0("<b>", term, "</b><br>FDR: ", signif(padj, 3),
                   "<br>size: ", ifelse(is.finite(size), size, "n/a"),
                   "<br>score: ", ifelse(is.finite(score), signif(score, 3), "n/a")),
    hoverinfo = "text"
  ) |>
    plotly::layout(
      title = list(text = title, x = 0.02),
      xaxis = list(title = xtitle, zeroline = TRUE),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 270)
    ) |>
    .rep_plotly_config("enrichment_bubble")
}

#' Interactive enrichment-map network based on shared genes
#' @keywords internal
.rep_iplot_enrichment_network <- function(tbl, top_n = 14L, min_jaccard = 0.08) {
  .rep_require("plotly")
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  member_col <- .rep_first_col(tbl, c("core_enrichment", "geneID", "gene_id",
                                      "geneIDs", "leadingEdge", "leading_edge", "Genes"))
  term_col <- .rep_first_col(tbl, c("Description", "pathway", "Name", "Term", "ID", "source"))
  padj_col <- .rep_first_col(tbl, c("p.adjust", "padj", "pGFdr", "Adjusted.P.value",
                                     "qvalue", "qvalues", "FDR"))
  if (is.na(member_col) || is.na(term_col) || is.na(padj_col)) return(NULL)
  d <- tbl
  d$.padj <- .rep_num(d[[padj_col]])
  d <- d[is.finite(d$.padj), , drop = FALSE]
  d <- utils::head(d[order(d$.padj), , drop = FALSE], top_n)
  split_members <- strsplit(as.character(d[[member_col]]), "[/;,| ]+", perl = TRUE)
  sets <- lapply(split_members, function(x) unique(x[nzchar(x)]))
  keep <- lengths(sets) >= 2L
  d <- d[keep, , drop = FALSE]
  sets <- sets[keep]
  n <- length(sets)
  if (n < 3L) return(NULL)

  j <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) for (k in seq_len(n)) {
    den <- length(union(sets[[i]], sets[[k]]))
    j[i, k] <- if (den) length(intersect(sets[[i]], sets[[k]])) / den else 0
  }
  labs <- make.unique(.rep_trunc(as.character(d[[term_col]]), 46L))
  coords <- tryCatch(
    stats::cmdscale(stats::as.dist(1 - j), k = 2L),
    error = function(e) NULL
  )
  if (is.null(coords) || ncol(coords) < 2L || any(!is.finite(coords))) {
    ang <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
    coords <- cbind(cos(ang), sin(ang))
  }

  p <- plotly::plot_ly()
  for (i in seq_len(n - 1L)) for (k in (i + 1L):n) {
    if (j[i, k] >= min_jaccard) {
      p <- plotly::add_trace(
        p, x = coords[c(i, k), 1L], y = coords[c(i, k), 2L],
        type = "scatter", mode = "lines", showlegend = FALSE,
        line = list(color = "rgba(105,125,145,0.35)",
                    width = 1 + 5 * j[i, k]),
        hoverinfo = "skip"
      )
    }
  }
  p <- plotly::add_trace(
    p, x = coords[, 1L], y = coords[, 2L],
    type = "scatter", mode = "markers+text", text = labs,
    textposition = "top center", showlegend = FALSE,
    marker = list(
      size = pmax(12, pmin(34, sqrt(lengths(sets)) * 5)),
      color = -log10(pmax(d$.padj, .Machine$double.xmin)),
      colorscale = "Viridis", showscale = TRUE,
      colorbar = list(title = "-log10(FDR)"),
      line = list(color = "#ffffff", width = 1)
    ),
    hovertext = paste0(
      "<b>", as.character(d[[term_col]]), "</b><br>FDR: ",
      signif(d$.padj, 3), "<br>genes: ", lengths(sets),
      "<br>connected terms: ", rowSums(j >= min_jaccard) - 1L
    ),
    hoverinfo = "text"
  )
  p |>
    plotly::layout(
      title = list(text = "Enrichment map: shared-gene network", x = 0.02),
      xaxis = list(title = "", visible = FALSE),
      yaxis = list(title = "", visible = FALSE),
      showlegend = FALSE, margin = list(l = 30, r = 30, b = 30, t = 65)
    ) |>
    .rep_plotly_config("enrichment_network")
}

#' Interactive running enrichment score for a ranked GSEA gene set
#' @keywords internal
.rep_iplot_gsea_running <- function(enrichment_tbl, dge_results) {
  .rep_require("plotly")
  if (is.null(enrichment_tbl) || !nrow(enrichment_tbl) ||
      is.null(dge_results) || !nrow(dge_results)) return(NULL)
  member_col <- .rep_first_col(enrichment_tbl, c("core_enrichment", "leading_edge",
                                                  "leadingEdge", "geneID", "gene_id"))
  if (is.na(member_col)) return(NULL)
  term_col <- .rep_first_col(enrichment_tbl, c("Description", "pathway", "Name", "Term", "ID"))
  stat_col <- .rep_first_col(dge_results, c("contrast_stat", "stat"))
  gene_col <- .rep_first_col(dge_results, c("gene", "symbol", "gene_id"))
  if (is.na(stat_col) || is.na(gene_col)) return(NULL)
  pcol <- .rep_first_col(enrichment_tbl, c("p.adjust", "padj", "qvalue", "FDR"))
  d <- enrichment_tbl
  if (!is.na(pcol)) {
    d$.padj <- .rep_num(d[[pcol]])
    d <- d[order(d$.padj, na.last = TRUE), , drop = FALSE]
  }
  if (!nrow(d)) return(NULL)
  set_genes <- unique(strsplit(as.character(d[[member_col]][[1L]]),
                               "[/;,| ]+", perl = TRUE)[[1L]])
  set_genes <- set_genes[nzchar(set_genes)]
  ranked <- data.frame(
    gene = as.character(dge_results[[gene_col]]),
    stat = .rep_num(dge_results[[stat_col]]),
    stringsAsFactors = FALSE
  )
  ranked <- ranked[is.finite(ranked$stat) & !is.na(ranked$gene) &
                     nzchar(ranked$gene), , drop = FALSE]
  ranked <- ranked[order(-ranked$stat), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$gene), , drop = FALSE]
  hit <- ranked$gene %in% set_genes
  nh <- sum(hit)
  if (nh < 2L || nh >= nrow(ranked)) return(NULL)
  weights <- abs(ranked$stat[hit])
  if (!sum(weights)) weights <- rep(1, nh)
  hit_step <- rep(0, nrow(ranked))
  hit_step[hit] <- weights / sum(weights)
  miss_step <- 1 / (nrow(ranked) - nh)
  es <- cumsum(ifelse(hit, hit_step, -miss_step))
  title_term <- if (!is.na(term_col)) as.character(d[[term_col]][[1L]]) else "Selected gene set"
  plotly::plot_ly(
    x = seq_along(es), y = es, type = "scatter", mode = "lines",
    line = list(color = "#2b5c8f", width = 2.5),
    hovertemplate = "rank = %{x}<br>running ES = %{y:.3f}<extra></extra>"
  ) |>
    plotly::add_trace(
      x = which(hit), y = es[hit], type = "scatter", mode = "markers",
      marker = list(color = "#c94a4a", size = 5),
      text = ranked$gene[hit], name = "Leading-edge genes",
      hovertemplate = "%{text}<br>rank = %{x}<br>ES = %{y:.3f}<extra></extra>"
    ) |>
    plotly::layout(
      title = list(text = paste0("GSEA running score: ", .rep_trunc(title_term, 70L)), x = 0.02),
      xaxis = list(title = "Genes ranked by signed test statistic"),
      yaxis = list(title = "Running enrichment score"),
      shapes = list(.rep_hline(0, "#777777"))
    ) |>
    .rep_plotly_config("gsea_running_score")
}

#' Interactive pairwise gene-overlap heatmap between enriched terms
#' @keywords internal
.rep_iplot_enrichment_overlap <- function(tbl, top_n = 20L) {
  .rep_require("plotly")
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  member_col <- .rep_first_col(tbl, c("core_enrichment", "geneID", "gene_id",
                                       "geneIDs", "leadingEdge", "leading_edge", "Genes"))
  term_col <- .rep_first_col(tbl, c("Description", "pathway", "Name", "Term", "ID", "source"))
  padj_col <- .rep_first_col(tbl, c("p.adjust", "padj", "pGFdr", "Adjusted.P.value",
                                     "qvalue", "qvalues", "FDR"))
  if (is.na(member_col) || is.na(term_col) || is.na(padj_col)) return(NULL)
  d <- tbl
  d$.padj <- .rep_num(d[[padj_col]])
  d <- d[is.finite(d$.padj), , drop = FALSE]
  d <- utils::head(d[order(d$.padj), , drop = FALSE], top_n)
  if (nrow(d) < 3L) return(NULL)
  split_members <- strsplit(as.character(d[[member_col]]), "[/;,| ]+", perl = TRUE)
  sets <- lapply(split_members, function(x) unique(x[nzchar(x)]))
  if (sum(lengths(sets) > 0L) < 3L) return(NULL)
  labs <- make.unique(.rep_trunc(as.character(d[[term_col]]), 42L))
  j <- matrix(0, nrow = length(sets), ncol = length(sets), dimnames = list(labs, labs))
  for (i in seq_along(sets)) for (k in seq_along(sets)) {
    den <- length(union(sets[[i]], sets[[k]]))
    j[i, k] <- if (den) length(intersect(sets[[i]], sets[[k]])) / den else 0
  }
  if (max(j, na.rm = TRUE) <= 0) return(NULL)
  plotly::plot_ly(
    x = labs, y = labs, z = j, type = "heatmap", zmin = 0, zmax = 1,
    colorscale = "Blues", colorbar = list(title = "Jaccard"),
    hovertemplate = "%{y}<br>%{x}<br>shared-gene Jaccard = %{z:.2f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = "Gene overlap between enriched terms", x = 0.02),
      xaxis = list(title = "", tickangle = -40),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 260, b = 150)
    ) |>
    .rep_plotly_config("enrichment_overlap")
}

# ---- tables ------------------------------------------------------------------

#' Interactive result table
#'
#' A DT table with per-column filters and copy / CSV / Excel buttons.
#'
#' @param df A data.frame.
#' @param caption Optional caption.
#' @param page_length Rows per page.
#' @param max_rows Rows beyond this limit are dropped (a note is added to the
#'   caption); the complete table is always available as a file in `out_dir`.
#' @param digits Significant digits for numeric columns.
#' @param filename Base name of the exported CSV / Excel files.
#' @param round_payload `TRUE` (default) rounds real-valued columns to
#'   `digits` significant digits in R BEFORE serialisation.
#' @param escape HTML-escape cell contents (default TRUE). Set to `FALSE`
#'   when a column holds pre-built `<a href=...>` markup (e.g. clickable
#'   file links built by `.xp_file_link()`).
#' @return A DT htmlwidget.
#' @export
itable <- function(
  df, caption = NULL, page_length = 25L, max_rows = Inf, digits = 3L,
  filename = "expressom_table", round_payload = TRUE, escape = TRUE
) {
  .rep_require("DT", "to display interactive tables")
  if (is.null(df) || nrow(df) == 0L) {
    df <- data.frame(Note = "No rows to display", stringsAsFactors = FALSE)
  }
  if (nrow(df) > max_rows) {
    caption <- paste0(
      caption %||% "",
      " (showing the first ", format(max_rows, big.mark = ","), " of ",
      format(nrow(df), big.mark = ","), " rows; see the exported file for all rows)"
    )
    df <- df[seq_len(max_rows), , drop = FALSE]
  }
  # significant digits only for real-valued columns: counts stay integers
  num_cols <- colnames(df)[vapply(df, function(x) {
    is.numeric(x) && any(is.finite(x) & x != round(x))
  }, logical(1))]

  if (round_payload && length(num_cols)) {
    for (cn in num_cols) {
      df[[cn]] <- signif(df[[cn]], digits)
    }
  }

  dt <- DT::datatable(
    df,
    caption = caption,
    rownames = FALSE,
    filter = "top",
    extensions = "Buttons",
    class = "compact stripe hover",
    escape = escape,
    options = list(
      pageLength = page_length,
      lengthMenu = list(c(10, 25, 50, 100, 250, 500), c(10, 25, 50, 100, 250, 500)),
      searchDelay = 250,
      stateSave = TRUE,
      orderCellsTop = TRUE,
      autoWidth = FALSE,
      scrollX = TRUE,
      deferRender = TRUE,
      dom = "Bfltip",
      buttons = list(
        "copy",
        list(extend = "csv", filename = filename),
        list(extend = "excel", filename = filename)
      )
    )
  )
  if (length(num_cols)) {
    dt <- DT::formatSignif(dt, columns = num_cols, digits = digits)
  }
  dt
}