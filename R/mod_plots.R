#' Exploratory Data Analysis for DESeq2
#'
#' @param dds DESeqDataSet object
#' @param edb Ensembl Database
#' @param out_dir Output directory
#' @param level Level comparison (e.g. "Treated")
#' @param base Base Level (e.g. "Control")
#' @param main_condition The primary modeled condition
#' @param batch_col Optional batch column name for limma correction and before/after PCA
#' @export
run_eda <- function(dds, edb, out_dir, level, base, main_condition, batch_col = NULL) {
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

  # --- PCA before batch correction ---
  p_orig <- plot_custom_pca(rld, condition = main_condition, batch = batch_col,
                            title = paste0("PCA Before Correction (", level, " vs ", base, ")"),
                            return_plot = TRUE)
  pdf(file.path(plot_dir, paste0("PCA_", level, "_vs_", base, ".pdf")), width = 9, height = 7)
  print(p_orig)
  dev.off()

  # --- Batch correction + PCA after ---
  if (!is.null(batch_col) && batch_col %in% colnames(SummarizedExperiment::colData(dds)) &&
      requireNamespace("limma", quietly = TRUE)) {
    batch_vec <- SummarizedExperiment::colData(dds)[[batch_col]]
    if (length(unique(batch_vec)) > 1) {
      message("Applying limma batch correction for PCA visualisation...")
      design_mat <- model.matrix(as.formula(paste0("~ ", main_condition)),
                                 data = as.data.frame(SummarizedExperiment::colData(rld)))
      corrected_mat <- tryCatch(
        limma::removeBatchEffect(SummarizedExperiment::assay(rld), batch = batch_vec, design = design_mat),
        error = function(e) { warning("Batch correction failed: ", e$message); NULL }
      )
      if (!is.null(corrected_mat)) {
        rld_corrected <- rld
        SummarizedExperiment::assay(rld_corrected) <- corrected_mat
        p_corr <- plot_custom_pca(rld_corrected, condition = main_condition, batch = batch_col,
                                  title = "PCA After Batch Correction (limma)",
                                  return_plot = TRUE)
        pdf(file.path(plot_dir, paste0("PCA_BatchCorrected_", level, "_vs_", base, ".pdf")), width = 9, height = 7)
        print(p_corr)
        dev.off()
        message("   -> Batch-corrected PCA saved to: ", plot_dir)
      }
    }
  }

  # Check correlation matrix for NA/NaN
  rld_cor <- cor(SummarizedExperiment::assay(rld))
  if (anyNA(rld_cor)) {
    warning("Correlation matrix contains NA/NaN values. ",
            "Number of NA: ", sum(is.na(rld_cor)))
    # Optionally replace NA with 0 or remove affected samples
    rld_cor[is.na(rld_cor)] <- 0
  }

  # Verify annotation rownames match
  anno <- as.data.frame(SummarizedExperiment::colData(dds)[, main_condition, drop = FALSE])
  if (!all(rownames(anno) == colnames(rld_cor))) {
    warning("Annotation rownames do not match correlation matrix columns.")
  }

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
generate_bulk_visualizations <- function(dds, edb, res_shrunken, res_unshrunken, results_data, out_dir, level, base, main_condition, top_genes, padj_cutoff, highlight_genes = NULL, batch_col = NULL) {
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

  # pivot_longer replaces the deprecated tidyr::gather()
  id_cols <- intersect(c("gene", "ensembl", "symbol", "entrezid"), colnames(topx_sigOE_norm))
  gathered_top <- tidyr::pivot_longer(
    topx_sigOE_norm,
    cols      = -dplyr::all_of(id_cols),
    names_to  = "sample",
    values_to = "normalized_counts"
  )
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
  pca_df$sample_label <- rownames(pca_df)

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
                                             color = .data[[condition]]))
  # stat_ellipse needs >= 3 points per group; also requires > 1 group
  min_pts_per_group <- min(table(pca_df[[condition]]))
  if (ellipse && length(unique(pca_df[[condition]])) > 1 && min_pts_per_group >= 3) {
    p <- p + ggplot2::stat_ellipse(level = 0.95, linetype = 2)
  } else if (ellipse) {
    message("plot_custom_pca: ellipse skipped (need >= 3 samples per group and >= 2 groups)")
  }
  p <- p + ggplot2::geom_point(size = 3.5, alpha = 0.8)
  if (!is.null(batch) && batch %in% colnames(pca_df)) {
    p <- p + ggplot2::aes(shape = .data[[batch]])
  }
  p <- p +
    ggrepel::geom_text_repel(ggplot2::aes(label = .data[["sample_label"]]),
                             size = 3, show.legend = FALSE,
                             box.padding = 0.4, max.overlaps = Inf) +
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

  # Drop genes with zero variance before z-scoring (scale() produces NaN otherwise)
  row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  zero_var <- row_sds == 0 | is.na(row_sds)
  if (any(zero_var)) {
    message("   -> plot_top_genes_heatmap: removing ", sum(zero_var),
            " zero-variance gene(s) before z-score scaling")
    mat <- mat[!zero_var, , drop = FALSE]
  }
  if (nrow(mat) == 0) {
    message("   -> No genes with non-zero variance – skipping heatmap")
    return(invisible(NULL))
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

#' Plot Average Z-Score by Condition for Gene Set(s)
#'
#' Computes a per‑gene z‑score across samples, then averages the z‑scores across
#' all genes within each gene set for every sample, producing one "module score"
#' per sample per set. If `set_name` is supplied, only one set is expected and
#' the plot is saved as a single PDF (without faceting) using that name in the
#' title and filename.
#'
#' @param dds A DESeqDataSet object (rownames should already be gene symbols)
#' @param gene_sets Named list of character vectors, e.g.
#'   \code{list(Tightness = c("Cdh5","Pdgfa"), "Lipid Scavengers" = c("Cd36","Stab1"))}.
#'   List names are used as facet/panel titles unless `set_name` is provided.
#' @param condition_col Character string representing the design metadata column
#' @param level String representing the foreground group
#' @param base String representing the background reference group
#' @param plot_dir String directory path where the PDF will be saved
#' @param include_global Logical: append a pooled "Global (All Genes)" panel
#'   combining every gene across all supplied sets (default TRUE, ignored if
#'   `set_name` is not NULL)
#' @param set_name Optional character string. If provided, `gene_sets` must be a
#'   list of length 1; this name is used in the output filename and plot title,
#'   and faceting is suppressed.
#' @return Invisibly returns the ggplot object
#' @export
plot_geneset_zscore_avg <- function(dds, gene_sets, condition_col, level, base,
                                    plot_dir, include_global = TRUE,
                                    set_name = NULL) {
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("`gene_sets` must be a named list, e.g. list(SetA = c('Gene1','Gene2'), SetB = c('Gene3')).")
  }

  # If set_name is given, we treat gene_sets as a single set (list of length 1)
  if (!is.null(set_name)) {
    if (length(gene_sets) != 1) {
      warning("set_name provided but gene_sets has length != 1. Using only the first set.")
      gene_sets <- gene_sets[1]
    }
    # Override include_global: no global panel when plotting a single set separately
    include_global <- FALSE
  }

  norm_counts    <- DESeq2::counts(dds, normalized = TRUE)
  meta_data      <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_data$sample <- rownames(meta_data)
  valid_samples  <- meta_data$sample[meta_data[[condition_col]] %in% c(base, level)]

  # Helper: given a character vector of genes, return a per-sample average
  # z-score (NULL if none of the genes are present) plus which genes were missing.
  .avg_zscore_for_genes <- function(genes) {
    present_genes <- intersect(genes, rownames(norm_counts))
    missing_genes <- setdiff(genes, rownames(norm_counts))
    if (length(present_genes) == 0) {
      return(list(avg_z = NULL, present = present_genes, missing = missing_genes))
    }
    expr_mat              <- norm_counts[present_genes, valid_samples, drop = FALSE]
    row_means             <- rowMeans(expr_mat)
    row_sds               <- apply(expr_mat, 1, stats::sd)
    row_sds[row_sds == 0] <- 1
    z_mat                 <- (expr_mat - row_means) / row_sds
    list(avg_z = colMeans(z_mat, na.rm = TRUE), present = present_genes, missing = missing_genes)
  }

  set_results    <- list()
  missing_report <- list()

  for (set_name_ in names(gene_sets)) {
    res <- .avg_zscore_for_genes(gene_sets[[set_name_]])
    if (length(res$missing) > 0) missing_report[[set_name_]] <- res$missing
    if (is.null(res$avg_z)) {
      warning("None of the genes in gene set '", set_name_, "' were found in the data. Skipping this set.")
      next
    }
    set_results[[set_name_]] <- data.frame(
      sample     = names(res$avg_z),
      avg_zscore = as.numeric(res$avg_z),
      gene_set   = paste0(set_name_, " (n=", length(res$present), ")"),
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(include_global) && length(gene_sets) > 1) {
    all_genes_pooled <- unique(unlist(gene_sets, use.names = FALSE))
    res_global <- .avg_zscore_for_genes(all_genes_pooled)
    if (!is.null(res_global$avg_z)) {
      set_results[["__global__"]] <- data.frame(
        sample     = names(res_global$avg_z),
        avg_zscore = as.numeric(res_global$avg_z),
        gene_set   = paste0("Global (All Genes) (n=", length(res_global$present), ")"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(missing_report) > 0) {
    for (nm in names(missing_report)) {
      message("   -> Note: genes not found in '", nm, "' set (check spelling/case): ",
              paste(missing_report[[nm]], collapse = ", "))
    }
  }

  if (length(set_results) == 0) stop("None of the specified gene sets contained any matching genes.")

  plot_df           <- do.call(rbind, set_results)
  plot_df           <- merge(plot_df, meta_data[, c("sample", condition_col)], by = "sample")
  plot_df$Condition <- factor(plot_df[[condition_col]], levels = c(base, level))

  # For single set, no need to set factor levels for faceting
  if (is.null(set_name)) {
    set_levels <- vapply(names(set_results), function(nm) unique(set_results[[nm]]$gene_set), character(1))
    plot_df$gene_set <- factor(plot_df$gene_set, levels = set_levels)
  }

  summary_df <- do.call(rbind, lapply(split(plot_df, list(plot_df$gene_set, plot_df$Condition)), function(d) {
    if (nrow(d) == 0) return(NULL)
    data.frame(
      gene_set  = d$gene_set[1],
      Condition = d$Condition[1],
      mean_z    = mean(d$avg_zscore),
      sem       = if (nrow(d) > 1) stats::sd(d$avg_zscore) / sqrt(nrow(d)) else 0
    )
  }))

  comp_name <- paste0(level, "_vs_", base)

  # Build base plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = summary_df,
      ggplot2::aes(x = Condition, y = mean_z, fill = Condition),
      stat = "identity", alpha = 0.6, width = 0.6
    ) +
    ggplot2::geom_errorbar(
      data = summary_df,
      ggplot2::aes(x = Condition, ymin = mean_z - sem, ymax = mean_z + sem),
      width = 0.2
    ) +
    ggplot2::geom_jitter(
      data = plot_df,
      ggplot2::aes(x = Condition, y = avg_zscore, fill = Condition),
      width = 0.15, size = 2, shape = 21, color = "black"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = stats::setNames(c("dodgerblue4", "red3"), c(base, level))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text      = ggplot2::element_text(face = "bold", size = 11),
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  # Add faceting or single‑set title
  if (!is.null(set_name)) {
    p <- p + ggplot2::labs(
      title = paste0("Average Z‑score: ", set_name, " (", level, " vs ", base, ")"),
      x = NULL, y = "Mean Z‑score (\u00b1 SEM)"
    )
    pdf_path <- file.path(plot_dir, paste0("GeneSet_Zscore_Average_", comp_name, "_", set_name, ".pdf"))
  } else {
    p <- p + ggplot2::facet_wrap(~gene_set, scales = "free_y") +
      ggplot2::labs(
        title = "Average Z‑score by Condition (Gene Set Module Score)",
        x = NULL, y = "Mean Z‑score (\u00b1 SEM)"
      )
    pdf_path <- file.path(plot_dir, paste0("GeneSet_Zscore_Average_", comp_name, ".pdf"))
  }

  # Save with appropriate dimensions
  if (!is.null(set_name)) {
    width <- 6
  } else {
    width <- 4.5 * length(set_results) + 1
  }
  ggplot2::ggsave(filename = pdf_path, plot = p, width = width, height = 5, device = "pdf")
  message("   -> Gene set average Z‑score plot saved to: ", pdf_path)
  return(invisible(p))
}

#' Plot Z-Score by Condition for Individual Genes
#'
#' Companion to \code{plot_geneset_zscore_avg}: computes the same per‑gene
#' z‑score across samples and plots a grouped bar chart (mean ± SEM) with
#' every gene along the x‑axis, faceted by gene set if multiple sets are
#' supplied. When `set_name` is given, a single set is plotted without faceting,
#' and the filename includes the set name.
#'
#' @param dds A DESeqDataSet object (rownames should already be gene symbols)
#' @param gene_sets Named list of character vectors.
#' @param condition_col Character column for grouping
#' @param level Foreground group
#' @param base Background group
#' @param plot_dir Output directory
#' @param show_points Logical: overlay individual sample‑level jittered points (default TRUE)
#' @param set_name Optional character string. If provided, `gene_sets` must be a
#'   list of length 1; this name is used in the output filename and plot title,
#'   and faceting is suppressed.
#' @return Invisibly returns the ggplot object
#' @export
plot_gene_zscore_individual <- function(dds, gene_sets, condition_col, level, base,
                                        plot_dir, show_points = TRUE,
                                        set_name = NULL) {
  if (is.character(gene_sets)) gene_sets <- list("Genes" = gene_sets)
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("`gene_sets` must be a named list, e.g. list(SetA = c('Gene1','Gene2'), SetB = c('Gene3')).")
  }

  if (!is.null(set_name)) {
    if (length(gene_sets) != 1) {
      warning("set_name provided but gene_sets has length != 1. Using only the first set.")
      gene_sets <- gene_sets[1]
    }
  }

  norm_counts      <- DESeq2::counts(dds, normalized = TRUE)
  meta_data        <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_data$sample <- rownames(meta_data)
  valid_samples    <- meta_data$sample[meta_data[[condition_col]] %in% c(base, level)]

  # Map each gene to the (first) set it belongs to, preserving list order
  gene_lookup <- unlist(lapply(names(gene_sets), function(nm) {
    stats::setNames(rep(nm, length(gene_sets[[nm]])), gene_sets[[nm]])
  }))
  gene_lookup <- gene_lookup[!duplicated(names(gene_lookup))]

  all_genes     <- names(gene_lookup)
  present_genes <- intersect(all_genes, rownames(norm_counts))
  missing_genes <- setdiff(all_genes, rownames(norm_counts))
  if (length(present_genes) == 0) stop("None of the specified genes were found in the data.")
  if (length(missing_genes) > 0) {
    message("   -> Note: genes not found (check spelling/case): ", paste(missing_genes, collapse = ", "))
  }

  expr_mat              <- norm_counts[present_genes, valid_samples, drop = FALSE]
  row_means             <- rowMeans(expr_mat)
  row_sds               <- apply(expr_mat, 1, stats::sd)
  row_sds[row_sds == 0] <- 1
  z_mat                 <- (expr_mat - row_means) / row_sds

  df_long      <- as.data.frame(z_mat)
  df_long$gene <- rownames(df_long)
  plot_df      <- tidyr::pivot_longer(df_long, cols = -gene, names_to = "sample", values_to = "z_score")
  plot_df      <- merge(plot_df, meta_data[, c("sample", condition_col)], by = "sample")
  plot_df$Condition <- factor(plot_df[[condition_col]], levels = c(base, level))

  # Preserve gene order within each set, and keep genes grouped along the x‑axis
  gene_order       <- present_genes[order(match(gene_lookup[present_genes], names(gene_sets)))]
  plot_df$gene_set <- factor(gene_lookup[plot_df$gene], levels = names(gene_sets))
  plot_df$gene     <- factor(plot_df$gene, levels = gene_order)

  summary_df <- do.call(rbind, lapply(split(plot_df, list(plot_df$gene, plot_df$Condition)), function(d) {
    if (nrow(d) == 0) return(NULL)
    data.frame(
      gene      = d$gene[1],
      gene_set  = d$gene_set[1],
      Condition = d$Condition[1],
      mean_z    = mean(d$z_score),
      sem       = if (nrow(d) > 1) stats::sd(d$z_score) / sqrt(nrow(d)) else 0
    )
  }))

  comp_name  <- paste0(level, "_vs_", base)
  dodge_pos  <- ggplot2::position_dodge(width = 0.75)

  p <- ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = summary_df,
      ggplot2::aes(x = gene, y = mean_z, fill = Condition),
      stat = "identity", position = dodge_pos, alpha = 0.75, width = 0.7, color = "black", linewidth = 0.3
    ) +
    ggplot2::geom_errorbar(
      data = summary_df,
      ggplot2::aes(x = gene, ymin = mean_z - sem, ymax = mean_z + sem, group = Condition),
      position = dodge_pos, width = 0.25
    )

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = gene, y = z_score, group = Condition),
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      size = 1.3, alpha = 0.6, color = "black"
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = stats::setNames(c("dodgerblue4", "red3"), c(base, level))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text      = ggplot2::element_text(face = "bold", size = 11),
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title    = ggplot2::element_blank()
    )

  if (!is.null(set_name)) {
    p <- p + ggplot2::labs(
      title = paste0("Individual Gene Z‑scores: ", set_name, " (", level, " vs ", base, ")"),
      x = "Features", y = "Z‑score (\u00b1 SEM)"
    )
    pdf_path <- file.path(plot_dir, paste0("Gene_Zscore_Individual_", comp_name, "_", set_name, ".pdf"))
  } else {
    n_sets <- length(unique(plot_df$gene_set))
    if (n_sets > 1) {
      p <- p + ggplot2::facet_wrap(~gene_set, scales = "free_x", nrow = 1)
    }
    p <- p + ggplot2::labs(
      title = "Z‑score by Condition, Individual Genes",
      x = "Features", y = "Z‑score (\u00b1 SEM)"
    )
    pdf_path <- file.path(plot_dir, paste0("Gene_Zscore_Individual_", comp_name, ".pdf"))
  }

  plot_width <- max(6, 0.9 * length(present_genes) + 2)
  ggplot2::ggsave(filename = pdf_path, plot = p,
                  width = plot_width, height = 5, device = "pdf", limitsize = FALSE)
  message("   -> Individual gene Z‑score plot saved to: ", pdf_path)
  return(invisible(p))
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