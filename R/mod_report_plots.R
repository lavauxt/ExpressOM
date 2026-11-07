# mod_report_plots.R -- interactive (plotly / DT) counterparts of the PDF
# figures, built from the result tables returned by load_expressom_results().
#
# Every function takes plain data.frames, returns an htmlwidget and is used
# both by the self-contained report and by the Shiny app, so the two always
# show the same figures. plotly and DT are only needed at report time.

.rep_plotly_config <- function(p, filename = "expressom_plot") {
  # Fill the container (a tab / grid cell) instead of knitr's fixed 672 px
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
  if (is.null(x)) {
    return(x)
  }
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

# A colour-blind friendly default palette for categorical groups
.rep_palette <- function(n) {
  if (n <= 0L) {
    return(character(0))
  }
  base <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#F0E442", "#999999"
  )
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

# ---- volcano / MA -------------------------------------------------------------

# Shared preparation for the volcano and MA plots: classify every gene as
# up / down / not significant and thin the non-significant cloud so that the
# self-contained HTML stays light even for 60k-gene tables.
#
# max_points is a hard cap on the total number of points drawn: all
# significant genes are kept and the non-significant cloud is thinned to
# fill the remaining budget (the sampling is seeded, with or without withr).
# The serialised numbers are rounded here because htmlwidgets serialises at
# up to 16 significant digits, which dominates the size of the embedded JSON.
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
  if (nrow(df) == 0L) {
    return(df)
  }
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

  # shrink the serialised payload (full-precision doubles -> 3-4 digits);
  # invisible at plot scale and typically halves the embedded JSON
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
    c("Not significant" = unname(colors[["ns"]]), "Significant (LRT)" = unname(colors[["up"]]))
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
      p,
      data = d,
      x = ~x_val, y = ~y_val,
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
#'
#' @param res_tbl DE result table (`log2FoldChange`, `padj`, a label column).
#' @param padj_cutoff,lfc_cutoff Significance thresholds (drawn as dashed lines).
#' @param colors Named vector with `up`, `down`, `ns` colours.
#' @param test_type "Wald" or "LRT" (LRT p-values are omnibus: no direction).
#' @param label_col Column used for hover labels and highlighting.
#' @param max_points Cap on the points drawn; all significant genes are kept,
#'   only the non-significant cloud is thinned.
#' @param highlight Optional gene labels to mark and name on the plot.
#' @param title Plot title.
#' @return A plotly htmlwidget.
#' @export
iplot_volcano <- function(
  res_tbl,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  colors = .de_direction_colors(),
  test_type = "Wald",
  label_col = "gene",
  max_points = 5000L,
  highlight = NULL,
  title = NULL
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
        meta = "xp-contrast", # drawn near-white in the dark theme (see expressom_theme.js)
        marker = list(color = "black", size = 9, symbol = "diamond-open", line = list(width = 2))
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
  res_tbl,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  colors = .de_direction_colors(),
  test_type = "Wald",
  label_col = "gene",
  max_points = 5000L,
  title = NULL
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
#'
#' @param scores PCA score table (`sample_label`, grouping columns, `PC1..PCn`).
#' @param variance Named numeric vector of percent variance per PC (or NULL).
#' @param x,y Principal components to plot.
#' @param group_col Column used to colour the samples (`NA` = single colour).
#' @param group_colors Named colour vector (values of `group_col`); missing
#'   groups get a colour-blind friendly palette.
#' @param title Plot title.
#' @return A plotly htmlwidget.
#' @export
iplot_pca <- function(
  scores,
  variance = NULL,
  x = "PC1",
  y = "PC2",
  group_col = NA_character_,
  group_colors = NULL,
  title = NULL
) {
  .rep_require("plotly")
  group_colors <- .rep_hex(group_colors)
  if (is.null(scores) || !all(c(x, y) %in% colnames(scores))) {
    return(.rep_empty_plot("PCA scores not available"))
  }
  df <- scores
  df$.id <- as.character(if ("plot_label" %in% colnames(df)) df$plot_label else df$sample_label)
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
    v <- if (!is.null(variance) && pc %in% names(variance)) {
      paste0(": ", variance[[pc]], "% variance")
    } else {
      ""
    }
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
        "<b>", d$.id, "</b><br>", if (!is.na(group_col)) paste0(group_col, ": ", g, "<br>") else "",
        x, ": ", signif(d[[x]], 4), "<br>", y, ": ", signif(d[[y]], 4)
      ),
      hoverinfo = "text", name = g,
      marker = list(size = 11, color = cols[[g]], opacity = 0.85, line = list(width = 1, color = "white"))
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

#' Interactive scree plot (percent variance per principal component)
#' @param variance Named numeric vector of percent variance per PC.
#' @param max_pcs Number of components shown.
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
#'
#' @param deg_counts Table written by `export_significant_results()`
#'   (`DEG_counts_<comparison>.txt`).
#' @param colors Named vector with `up` / `down` colours.
#' @param test_type "Wald" or "LRT".
#' @param use_lfc Count genes passing the log2FC filter (default) or not.
#' @export
iplot_deg_counts <- function(
  deg_counts,
  colors = .de_direction_colors(),
  test_type = "Wald",
  use_lfc = TRUE
) {
  .rep_require("plotly")
  colors <- .rep_hex(colors)
  if (is.null(deg_counts) || nrow(deg_counts) == 0L || !"padj_cutoff" %in% colnames(deg_counts)) {
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
          if (use_lfc && "log2FC_cutoff" %in% colnames(deg_counts)) {
            paste0(" (|log2FC| > ", deg_counts$log2FC_cutoff[[1L]], ")")
          } else {
            ""
          }
        ),
        x = 0.02
      ),
      xaxis = list(title = ""), yaxis = list(title = "number of genes")
    ) |>
    .rep_plotly_config("deg_counts")
}

# ---- expression heatmap -----------------------------------------------------------

#' Interactive heatmap of normalised counts (row z-scores)
#'
#' @param mat Matrix of normalised counts (genes x samples).
#' @param groups Optional group of each column (used to order and label samples).
#' @param colors Named vector with `up` / `down` colours for the diverging scale.
#' @param max_genes Number of genes shown.
#' @param title Plot title.
#' @export
iplot_heatmap <- function(
  mat,
  groups = NULL,
  colors = .de_direction_colors(),
  max_genes = 50L,
  title = NULL
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
    ord <- tryCatch(stats::hclust(stats::dist(z))$order, error = function(e) seq_len(nrow(z)))
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
      c(0, unname(colors[["down"]])), c(0.5, "#FFFFFF"), c(1, unname(colors[["up"]]))
    ),
    colorbar = list(title = "z-score"),
    hovertemplate = "%{y}<br>%{x}<br>z = %{z:.2f}<extra></extra>"
  ) |>
    plotly::layout(
      title = list(text = title %||% "Top differentially expressed genes (row z-score of log2 normalised counts)", x = 0.02),
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "", autorange = "reversed", tickfont = list(size = 9)),
      margin = list(l = 110, b = 120)
    ) |>
    .rep_plotly_config("heatmap")
}

# ---- enrichment ---------------------------------------------------------------------

#' Interactive enrichment plot
#'
#' Bar chart of the top terms of an ORA / GSEA / SPIA / TF table. GSEA tables
#' (with a normalised enrichment score) are drawn as NES bars coloured by
#' direction; the others as -log10(adjusted p) bars coloured by gene count.
#'
#' @param std Output of the enrichment standardiser (`term`, `padj`, `size`, `score`).
#' @param top_n Number of terms shown.
#' @param padj_cutoff Terms above this adjusted p-value are greyed out in the
#'   hover but still listed (so near-misses stay visible).
#' @param colors Named vector with `up` / `down` colours.
#' @param title Plot title.
#' @export
iplot_enrichment <- function(
  std,
  top_n = 20L,
  padj_cutoff = 0.05,
  colors = .de_direction_colors(),
  title = NULL
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

# ---- isoform ---------------------------------------------------------------------------

#' Interactive DTU bar chart: strongest differential transcript usage features
#'
#' @param dtu DTU result table (annotated CSV of the isoform report).
#' @param top_n Number of features shown.
#' @param padj_cutoff Adjusted p-value cut-off (dotted line).
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
  labels <- if (is.na(lab_col)) as.character(seq_len(nrow(d))) else as.character(d[[lab_col]])
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
#' @param p Numeric vector of p-values.
#' @param title Plot title.
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

# ---- tables ----------------------------------------------------------------------------

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
#'   `digits` significant digits in R BEFORE serialisation: DT::formatSignif
#'   only changes the rendered text, while the data itself would be embedded
#'   at full double precision (16 significant digits), which is what makes
#'   big tables heavy.
#' @return A DT htmlwidget.
#' @export
itable <- function(
  df,
  caption = NULL,
  page_length = 10L,
  max_rows = 800L,
  digits = 3L,
  filename = "expressom_table",
  round_payload = TRUE
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
