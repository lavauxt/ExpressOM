#' Exploratory Data Analysis for DESeq2
#'
#' @param dds DESeqDataSet object
#' @param edb Ensembl Database
#' @param out_dir Output directory
#' @param level Level comparison (e.g. "Treated")
#' @param base Base Level (e.g. "Control")
#' @param main_condition The primary modeled condition
#' @export
run_eda <- function(dds, edb, out_dir, level, base, main_condition) {
  rld <- DESeq2::rlog(dds, blind = TRUE)

  org_info <- get_organism_info(edb)
  org_obj <- .load_org_db(org_info$org_db)
  message("Mapping IDs for EDA plots...")
  symbols <- tryCatch({
    suppressMessages(AnnotationDbi::mapIds(
      org_obj, keys = rownames(rld), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
    ))
  }, error = function(e) {
    message("   -> Warning: Could not map keys. They may already be symbols. Using original row names.")
    stats::setNames(rownames(rld), rownames(rld))
  })

  symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
  rownames(rld) <- symbols

  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  if (requireNamespace("pcaExplorer", quietly = TRUE)) {
    pca_plot <- pcaExplorer::pcaplot(rld, intgroup = main_condition, ntop = 500)
    pdf(file = file.path(plot_dir, paste0("PCA_", level, "_vs_", base, ".pdf")), width = 9, height = 7)
    try({ print(pca_plot) })
    dev.off()
  } else {
    message("Skipping PCA plot: 'pcaExplorer' is not installed.")
  }

  rld_cor <- cor(SummarizedExperiment::assay(rld))
  anno <- as.data.frame(SummarizedExperiment::colData(dds)[, main_condition, drop = FALSE])

  safe_pdf(file.path(plot_dir, paste0("HeatMap_", level, "_vs_", base, ".pdf")), expr = {
    pheatmap::pheatmap(rld_cor, annotation_col = anno, main = "Sample Correlation (Gene Symbols)")
  })
}

#' Generate Bulk Visualizations
#'
#' This function generates a series of visualizations to summarize the differential expression results.
#' @export
#' @param dds DESeq object
#' @param edb Ensembl database object
#' @param res_shrunken Shrunken DE results
#' @param res_unshrunken Unshrunken DE results
#' @param results_data Object from export_significant_results
#' @param out_dir Output directory
#' @param level Level
#' @param base Base
#' @param main_condition Extracted primary condition factor
#' @param top_genes N top genes
#' @param padj_cutoff Adjusted p-value significance cutoff
#' @param highlight_genes Optional character vector of gene names to highlight in the Volcano plot
#' @return NULL
generate_bulk_visualizations <- function(dds, edb, res_shrunken, res_unshrunken, results_data, out_dir, level, base, main_condition, top_genes, padj_cutoff, highlight_genes = NULL) {
  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  org_info <- get_organism_info(edb)
  org_obj <- .load_org_db(org_info$org_db)
  safe_pdf(file.path(plot_dir, paste0("MAplot_unshrunken_", level, "_vs_", base, ".pdf")), expr = {
    DESeq2::plotMA(res_unshrunken, ylim = c(-2, 2))
    graphics::abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)
  })

  safe_pdf(file.path(plot_dir, paste0("MAplot_shrunken_", level, "_vs_", base, ".pdf")), expr = {
    DESeq2::plotMA(res_shrunken, ylim = c(-2, 2))
    graphics::abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)
  })

  topx_sigOE_genes <- head(results_data$sig_res[order(results_data$sig_res$padj), "gene"], top_genes)
  topx_sigOE_norm <- results_data$normalized_counts[results_data$normalized_counts$gene %in% topx_sigOE_genes, ]

  gathered_top <- tidyr::gather(topx_sigOE_norm, key = "sample", value = "normalized_counts", -c(gene, ensembl, symbol, entrezid))
  gathered_top$normalized_counts <- as.numeric(gathered_top$normalized_counts)

  meta_df <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_df$sample <- rownames(meta_df)
  topx_final <- dplyr::inner_join(meta_df, gathered_top, by = "sample")

  safe_pdf(file.path(plot_dir, paste0("Top", top_genes, "_DE_Genes_", level, "_vs_", base, ".pdf")), expr = {
    p_top <- ggplot2::ggplot(topx_final, ggplot2::aes(x = gene, y = normalized_counts, fill = .data[[main_condition]])) +
      ggplot2::geom_boxplot() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ggtitle(paste("Top", top_genes, "Significant DE Genes"))
    print(p_top)
  })
  safe_pdf(file.path(plot_dir, paste0("DE_Volcanoplot_", level, "_vs_", base, ".pdf")), expr = {
    suppressWarnings({
      vp <- EnhancedVolcano::EnhancedVolcano(
        results_data$res_tbl, lab = results_data$res_tbl$gene,
        selectLab = highlight_genes,
        drawConnectors = !is.null(highlight_genes),
        x = 'log2FoldChange', y = 'padj',
        title = paste(level, "vs", base),
        pCutoff = padj_cutoff, FCcutoff = 1.0, pointSize = 2.0, labSize = 4.0,
        col = c('grey', 'grey', 'grey', 'red2')
      )
      print(vp)
    })
  })
}


#' Custom PCA Plot (ggplot2, fully customisable)
#' @param vsd VST-transformed DESeqDataSet or matrix
#' @param condition Character column name in colData for colour grouping
#' @param batch Optional batch column for shape grouping
#' @param title Plot title
#' @param ellipse Logical, whether to add 95% confidence ellipses
#' @param return_plot If TRUE, returns ggplot object; otherwise prints
#' @export
plot_custom_pca <- function(vsd, condition, batch = NULL, title = "PCA", ellipse = TRUE, return_plot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  if (inherits(vsd, "SummarizedExperiment")) {
    mat     <- SummarizedExperiment::assay(vsd)
    coldata <- as.data.frame(SummarizedExperiment::colData(vsd))
  } else if (is.matrix(vsd)) {
    stop("If vsd is a matrix, you must provide coldata separately. This function currently expects a SummarizedExperiment object.")
  } else {
    stop("vsd must be a SummarizedExperiment object (e.g., DESeqDataSet or DESeqTransform) or a matrix")
  }

  pca        <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pca_df     <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], coldata)

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
                                             color = .data[[condition]]))
  if (ellipse && length(unique(pca_df[[condition]])) > 1) {
    p <- p + ggplot2::stat_ellipse(level = 0.95, linetype = 2)
  }
  p <- p + ggplot2::geom_point(size = 3.5, alpha = 0.8)
  if (!is.null(batch) && batch %in% colnames(pca_df)) {
    p <- p + ggplot2::aes(shape = .data[[batch]])
  }
  p <- p +
    ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = title)
  if (return_plot) return(p) else print(p)
}

#' Sample Correlation Heatmap
#' @param vsd VST-transformed DESeqDataSet
#' @param condition_col Column name for sample annotation
#' @param plot_dir Directory to save PDF
#' @param comp_name Comparison name for file naming
#' @export
plot_sample_correlation <- function(vsd, condition_col, plot_dir, comp_name) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed, skipping correlation plot")
    return(invisible(NULL))
  }
  expr_mat   <- SummarizedExperiment::assay(vsd)
  cor_mat    <- cor(expr_mat, method = "pearson")
  annotation <- as.data.frame(SummarizedExperiment::colData(vsd)[, condition_col, drop = FALSE])
  rownames(annotation) <- colnames(cor_mat)
  colnames(annotation) <- condition_col

  pdf(file.path(plot_dir, paste0("SampleCorrelation_", comp_name, ".pdf")), width = 8, height = 7)
  pheatmap::pheatmap(cor_mat,
                     annotation_col  = annotation,
                     main            = paste0("Sample Correlation Matrix - ", comp_name),
                     color           = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                     cluster_rows    = TRUE, cluster_cols = TRUE,
                     display_numbers = FALSE, fontsize_row = 8)
  dev.off()
  message("   -> Sample correlation heatmap saved to: ", plot_dir)
}

#' Top Genes Expression Heatmap
#'
#' Selects the top-N significant DE genes by padj, retrieves their normalised
#' counts from the DESeqDataSet using Ensembl IDs (the actual rownames of the
#' count matrix), then labels rows with gene symbols for readability.
#'
#' @param dds DESeqDataSet
#' @param results_data List containing res_tbl (data.frame with gene, ensembl, log2FoldChange, padj)
#' @param condition_col Column name for condition grouping
#' @param level Treatment group
#' @param base Control group
#' @param top_n Number of top DE genes to plot (by padj)
#' @param padj_cutoff Adjusted p-value threshold
#' @param plot_dir Output directory
#' @param batch_col Optional batch column for expression correction
#' @export
plot_top_genes_heatmap <- function(dds, results_data, condition_col, level, base,
                                   top_n = 30, padj_cutoff = 0.01,
                                   plot_dir, batch_col = NULL) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed, skipping top genes heatmap")
    return(invisible(NULL))
  }

  res_tbl <- results_data$res_tbl
  sig     <- res_tbl[which(!is.na(res_tbl$padj) & res_tbl$padj < padj_cutoff), ]
  if (nrow(sig) == 0) {
    message("No significant genes at padj < ", padj_cutoff, " – skipping top genes heatmap")
    return(invisible(NULL))
  }
  top_sig <- head(sig[order(sig$padj), ], n = top_n)

  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  # BUG FIX: dds rownames are Ensembl IDs; res_tbl$gene is symbols.
  # Match via the 'ensembl' column, then label rows with symbols.
  if ("ensembl" %in% colnames(top_sig) && any(!is.na(top_sig$ensembl))) {
    candidate_ids <- top_sig$ensembl[!is.na(top_sig$ensembl)]
    present_ids   <- intersect(candidate_ids, rownames(norm_counts))
    if (length(present_ids) == 0) {
      message("None of the top genes' Ensembl IDs found in counts matrix – skipping heatmap")
      return(invisible(NULL))
    }
    mat         <- norm_counts[present_ids, , drop = FALSE]
    sym_labels  <- top_sig$gene[match(present_ids, top_sig$ensembl)]
    sym_labels[is.na(sym_labels)] <- present_ids[is.na(sym_labels)]
    rownames(mat) <- make.unique(sym_labels)
  } else {
    # Fallback: rownames already gene symbols (e.g. dds_rep was passed)
    present_ids <- intersect(top_sig$gene, rownames(norm_counts))
    if (length(present_ids) == 0) {
      message("None of the top genes found in counts matrix – skipping heatmap")
      return(invisible(NULL))
    }
    mat <- norm_counts[present_ids, , drop = FALSE]
  }

  # Optional batch correction
  if (!is.null(batch_col) && batch_col %in% colnames(SummarizedExperiment::colData(dds)) &&
      requireNamespace("limma", quietly = TRUE)) {
    batch_vec <- SummarizedExperiment::colData(dds)[[batch_col]]
    if (length(unique(batch_vec)) > 1) {
      design_mat <- model.matrix(as.formula(paste0("~ ", condition_col)),
                                 data = as.data.frame(SummarizedExperiment::colData(dds)))
      mat <- tryCatch(
        limma::removeBatchEffect(mat, batch = batch_vec, design = design_mat),
        error = function(e) mat
      )
    }
  }

  mat_z            <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0

  coldata        <- as.data.frame(SummarizedExperiment::colData(dds))
  annotation_col <- coldata[, condition_col, drop = FALSE]
  rownames(annotation_col) <- colnames(mat)
  colnames(annotation_col) <- condition_col

  col_pal <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

  pdf(file.path(plot_dir, paste0("TopGenes_Heatmap_", level, "_vs_", base, ".pdf")),
      width = 10, height = max(6, nrow(mat_z) * 0.3))
  pheatmap::pheatmap(mat_z,
                     annotation_col = annotation_col,
                     main           = paste0("Top ", nrow(mat_z), " DE genes (", level, " vs ", base, ")"),
                     color          = col_pal,
                     cluster_rows   = TRUE, cluster_cols = TRUE,
                     scale          = "none",
                     fontsize_row   = 8,
                     show_rownames  = TRUE)
  dev.off()
  message("   -> Top genes heatmap saved to: ", plot_dir)
}

# ---------------------- Unchanged functions below ----------------------

#' Plot Volcano
#' @keywords internal
plot_volcano <- function(res_tbl, padj_cutoff, highlight_genes = NULL, title = "") {
  EnhancedVolcano::EnhancedVolcano(
    res_tbl, lab = res_tbl$gene,
    selectLab      = highlight_genes,
    drawConnectors = !is.null(highlight_genes),
    x = "log2FoldChange", y = "padj",
    title     = title,
    pCutoff   = padj_cutoff, FCcutoff = 1.0, pointSize = 2.0, labSize = 4.0,
    col = c("grey", "grey", "grey", "red2")
  )
}

#' Plot Individual Sample Z-Score Heatmap
#' @param dds A DESeqDataSet object
#' @param selected_genes Character vector of gene symbols to plot
#' @param condition_col Character string representing the design metadata column
#' @param level String representing the foreground group
#' @param base String representing the background reference group
#' @param plot_dir String directory path where the PDF will be saved
#' @export
plot_sample_zscore <- function(dds, selected_genes, condition_col, level, base, plot_dir) {
  norm_counts   <- DESeq2::counts(dds, normalized = TRUE)
  present_genes <- intersect(selected_genes, rownames(norm_counts))
  if (length(present_genes) == 0) stop("None of the specified genes were found.")
  meta_data          <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_data$sample   <- rownames(meta_data)
  valid_samples      <- meta_data$sample[meta_data[[condition_col]] %in% c(base, level)]
  expr_mat           <- norm_counts[present_genes, valid_samples, drop = FALSE]
  row_means          <- rowMeans(expr_mat)
  row_sds            <- apply(expr_mat, 1, stats::sd)
  row_sds[row_sds == 0] <- 1
  z_mat              <- (expr_mat - row_means) / row_sds
  df_long            <- as.data.frame(z_mat)
  df_long$gene       <- rownames(df_long)
  plot_df            <- tidyr::pivot_longer(df_long, cols = -gene, names_to = "sample", values_to = "z_score")
  plot_df            <- merge(plot_df, meta_data, by = "sample")
  plot_df$Condition  <- factor(plot_df[[condition_col]], levels = c(base, level))
  plot_df$gene       <- factor(plot_df$gene, levels = rev(present_genes))
  comp_name          <- paste0(level, "_vs_", base)
  pdf_path           <- file.path(plot_dir, paste0("Sample_Zscore_Heatmap_", comp_name, ".pdf"))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = gene, fill = z_score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::facet_grid(~Condition, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "dodgerblue4", mid = "white", high = "red3", midpoint = 0,
      name = "z-score", limits = c(-1.5, 1.5), oob = scales::squish,
      guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(z_score, 1)), color = "black", size = 2.5) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 9),
      axis.text.y = ggplot2::element_text(face = "italic", size = 10),
      panel.grid  = ggplot2::element_blank(),
      strip.text  = ggplot2::element_text(face = "bold", size = 12, margin = ggplot2::margin(b = 10)),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::labs(title = "Sample Expression Heatmap", x = NULL, y = NULL)
  calc_height <- max(4, length(present_genes) * 0.4)
  calc_width  <- max(5, length(valid_samples) * 0.6 + 2)
  ggplot2::ggsave(filename = pdf_path, plot = p, width = calc_width, height = calc_height, device = "pdf")
  message("   -> Sample Heatmap successfully exported to: ", pdf_path)
}

#' Plot Log2 Fold Change Heatmap (1-Column Format)
#' @param dds A DESeqDataSet object
#' @param selected_genes Character vector of gene symbols to plot
#' @param condition_col Character string representing the design metadata column
#' @param level String representing the foreground group
#' @param base String representing the background reference group
#' @param plot_dir String directory path where the PDF will be saved
#' @export
plot_l2fc_heatmap <- function(dds, selected_genes, condition_col, level, base, plot_dir) {
  res         <- DESeq2::results(dds, contrast = c(condition_col, level, base))
  res_df      <- as.data.frame(res)
  present_genes <- intersect(selected_genes, rownames(res_df))
  if (length(present_genes) == 0) stop("None of the specified genes were found in the results.")
  plot_df     <- data.frame(gene = present_genes, L2FC = res_df[present_genes, "log2FoldChange"])
  plot_df$gene <- factor(plot_df$gene, levels = rev(present_genes))
  plot_df$Comparison <- paste0(level, " vs ", base)
  comp_name   <- paste0(level, "_vs_", base)
  pdf_path    <- file.path(plot_dir, paste0("L2FC_Heatmap_", comp_name, ".pdf"))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Comparison, y = gene, fill = L2FC)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "dodgerblue4", mid = "white", high = "red3", midpoint = 0,
      name = "Log2 FC",
      guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(L2FC, 2)), color = "black", size = 3.5) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "italic", size = 10),
      panel.grid  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold", margin = ggplot2::margin(b = 15))
    ) +
    ggplot2::labs(title = "Log2 Fold Change Heatmap", x = NULL, y = NULL)
  calc_height <- max(4, length(present_genes) * 0.4)
  ggplot2::ggsave(filename = pdf_path, plot = p, width = 4, height = calc_height, device = "pdf")
  message("   -> L2FC Heatmap successfully exported to: ", pdf_path)
}

#' SPIA two-way evidence plot (ggplot2 + ggrepel)
#' @param x A SPIA results data.frame
#' @param threshold Threshold for significance
#' @export
plotP_fork <- function(x, threshold = 0.01) {
  if (!inherits(x, "data.frame") | dim(x)[1] < 1 |
      !all(c("ID", "pNDE", "pPERT", "pG", "pGFdr", "pGFWER") %in% names(x))) {
    stop("SPIA graph can be applied only to a dataframe produced by SPIA function")
  }
  if (threshold < x[1, "pGFdr"]) {
    message("The threshold value was corrected to be equal to ", x[1, "pGFdr"])
    threshold <- x[1, "pGFdr"]
  }
  df  <- x
  pb  <- df$pPERT
  ph  <- df$pNDE
  # .combfunc and .getP2 are defined once in utils_core.R — no duplicate needed here
  combinemethod <- ifelse(
    sum(.combfunc(pb, ph, "fisher") == df$pG) > sum(.combfunc(pb, ph, "norminv") == df$pG),
    "fisher", "norminv"
  )
  okx <- (ph < 1e-6)
  oky <- (pb < 1e-6)
  ph[ph < 1e-6] <- 1e-6
  pb[pb < 1e-6] <- 1e-6
  df$x_val  <- -log(ph)
  df$y_val  <- -log(pb)
  df$Group  <- "Not Significant"
  df$Group[df$pGFdr  <= threshold] <- "FDR"
  df$Group[df$pGFWER <= threshold] <- "FWER"
  df$Group  <- factor(df$Group, levels = c("Not Significant", "FDR", "FWER"))
  df$Label  <- ""
  sig_idx   <- df$Group %in% c("FDR", "FWER")
  if (any(sig_idx)) df$Label[sig_idx] <- as.character(df$ID[sig_idx])

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_val, y = y_val)) +
    ggplot2::geom_point(ggplot2::aes(color = Group), size = 2.5) +
    ggplot2::scale_color_manual(values = c("Not Significant" = "black", "FDR" = "blue", "FWER" = "red"))

  tr_red <- threshold / nrow(na.omit(x))
  if (combinemethod == "fisher") {
    val_red  <- .getP2(tr_red, "fisher") / 2
    line_red <- data.frame(x = c(0, val_red), y = c(val_red, 0))
    p <- p + ggplot2::geom_path(data = line_red, ggplot2::aes(x = x, y = y), color = "red", linewidth = 1)
  } else {
    somep1 <- exp(seq(from = min(log(ph)), to = max(log(ph)), length = 200))
    somep2 <- pnorm(qnorm(tr_red) * sqrt(2) - qnorm(somep1))
    p <- p + ggplot2::geom_line(data = data.frame(x = -log(somep1), y = -log(somep2)),
                                ggplot2::aes(x = x, y = y), color = "red", linewidth = 1)
  }

  tr_blue_old <- tr_red
  tr_blue     <- suppressWarnings(max(df$pG[df$pGFdr <= threshold], na.rm = TRUE))
  if (is.infinite(tr_blue) || tr_blue <= tr_blue_old) tr_blue <- tr_blue_old * 1.03
  if (combinemethod == "fisher") {
    val_blue  <- .getP2(tr_blue, "fisher") / 2
    line_blue <- data.frame(x = c(0, val_blue), y = c(val_blue, 0))
    p <- p + ggplot2::geom_path(data = line_blue, ggplot2::aes(x = x, y = y), color = "blue", linewidth = 1)
  } else {
    somep1 <- exp(seq(from = min(log(ph)), to = max(log(ph)), length = 200))
    somep2 <- pnorm(qnorm(tr_blue) * sqrt(2) - qnorm(somep1))
    p <- p + ggplot2::geom_line(data = data.frame(x = -log(somep1), y = -log(somep2)),
                                ggplot2::aes(x = x, y = y), color = "blue", linewidth = 1)
  }

  p <- p + ggrepel::geom_text_repel(
    ggplot2::aes(label = Label, color = Group),
    size = 3.5, box.padding = 0.6, max.overlaps = Inf, show.legend = FALSE
  )
  if (any(okx)) {
    p <- p + ggplot2::geom_text(data = df[okx, ], ggplot2::aes(x = x_val - 0.15, y = y_val),
                                label = "|", size = 4, color = "black")
  }
  if (any(oky)) {
    p <- p + ggplot2::geom_text(data = df[oky, ], ggplot2::aes(x = x_val, y = y_val - 0.15),
                                label = "_", size = 4, color = "black", vjust = 1)
  }
  max_val <- max(c(df$x_val, df$y_val) + 1, na.rm = TRUE)
  p <- p +
    ggplot2::coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    ggplot2::labs(title = "SPIA two-way evidence plot",
                  x = "-log(P NDE)", y = "-log(P PERT)", color = "Significance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  print(p)
  return(invisible(p))
}
# NOTE: .combfunc and .getP2 are defined in utils_core.R.
# Duplicate definitions have been removed from this file to avoid divergence.