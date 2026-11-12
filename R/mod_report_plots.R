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
  max_genes = 50L, title = NULL
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
  if (!is.null(groups) && length(groups) == ncol(z)) {
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
#' (mor direction), exactly as the static PDF: "Activé (attendu)",
#' "Réprimé (attendu)", "Autre variation", "Non significatif".
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
    "Non significatif"  = "#BDBDBD",
    "Autre variation"   = "#ff7f00",
    "Activé (attendu)"  = "#1f78b4",
    "Réprimé (attendu)" = "#e31a1c"
  )
  ord <- c("Non significatif", "Autre variation",
           "Activé (attendu)", "Réprimé (attendu)")

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
        size    = if (cat == "Non significatif") 6 else 9,
        opacity = if (cat == "Non significatif") 0.55 else 0.9
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
  df, caption = NULL, page_length = 10L, max_rows = 800L, digits = 3L,
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
      scrollX = TRUE,
      deferRender = TRUE,
      dom = "Bfrtip",
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