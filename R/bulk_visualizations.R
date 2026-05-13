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
  # 1. Clean background devices
  graphics.off()
  
  # 2. Transformation
  rld <- DESeq2::rlog(dds, blind=TRUE)
  
   # 3. Detect Species from edb
  org_info <- get_organism_info(edb)
  org_db <- org_info$org_db

  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop("Package '", org_db, "' is required for identifier mapping in EDA. Please install it.")
  }
  
  # Extract the actual database object
  org_obj <- getExportedValue(org_db, org_db)

  # 4. Map Ensembl to Symbols for Plot Labels
  message("Mapping IDs for EDA plots...")
  symbols <- tryCatch({
    suppressMessages(AnnotationDbi::mapIds(
      org_obj, keys = rownames(rld), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
    ))
  }, error = function(e) {
    message("   -> Warning: Could not map keys. They may already be symbols. Using original row names.")
    stats::setNames(rownames(rld), rownames(rld))
  })
  
  # Safeguard: Keep original ID if mapping fails
  symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
  rownames(rld) <- symbols

  # 5. PCA Plot
  if (requireNamespace("pcaExplorer", quietly = TRUE)) {
    plot_dir <- file.path(out_dir, "Plots")
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    
    pca_plot <- pcaExplorer::pcaplot(rld, intgroup = main_condition, ntop = 500)
    pdf(file = file.path(plot_dir, paste0("PCA_", level, "_vs_", base, ".pdf")), width = 9, height = 7)
    try({
      print(pca_plot)
    })
    dev.off()
  } else {
    message("Skipping PCA plot: 'pcaExplorer' is not installed.")
  }
  
  # 6. Correlation Heatmap
  rld_cor <- cor(SummarizedExperiment::assay(rld))
  anno <- as.data.frame(SummarizedExperiment::colData(dds)[, main_condition, drop=FALSE])
  
  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  pdf(file = file.path(plot_dir, paste0("HeatMap_", level, "_vs_", base, ".pdf")))
  pheatmap::pheatmap(rld_cor, annotation_col = anno, main = "Sample Correlation (Gene Symbols)")
  dev.off()
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
  # Setup Directories
  plot_dir <- file.path(out_dir, "Plots")
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  org_info <- get_organism_info(edb)
  org_db <- org_info$org_db
  org_obj <- getExportedValue(org_db, org_db)

  # 1. MA Plots (QC mapped into Plots)
  pdf(file = file.path(plot_dir, paste0("MAplot_unshrunken_",level,"_vs_",base,".pdf")))
  DESeq2::plotMA(res_unshrunken, ylim=c(-2,2))
  graphics::abline(h=c(-1,1), col="dodgerblue", lwd=2)
  dev.off()
  
  pdf(file = file.path(plot_dir, paste0("MAplot_shrunken_",level,"_vs_",base,".pdf")))
  DESeq2::plotMA(res_shrunken, ylim=c(-2,2))
  graphics::abline(h=c(-1,1), col="dodgerblue", lwd=2)
  dev.off()
  
  # 2. Top Genes Normalized Counts Plot
  topx_sigOE_genes <- head(results_data$sig_res[order(results_data$sig_res$padj), "gene"], top_genes)
  topx_sigOE_norm <- results_data$normalized_counts[results_data$normalized_counts$gene %in% topx_sigOE_genes, ]
  
  # Ensure we subtract all identifier columns from the gather/pivot
  gathered_top <- tidyr::gather(topx_sigOE_norm, key = "sample", value = "normalized_counts", -c(gene, ensembl, symbol))
  gathered_top$normalized_counts <- as.numeric(gathered_top$normalized_counts)
  
  meta_df <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_df$sample <- rownames(meta_df)
  topx_final <- dplyr::inner_join(meta_df, gathered_top, by="sample")
  
  pdf(file= file.path(plot_dir, paste0("Top",top_genes,"_DE_Genes_",level,"_vs_",base,".pdf")))
  p_top <- ggplot2::ggplot(topx_final, ggplot2::aes(x = gene, y = normalized_counts, fill = !!rlang::sym(main_condition))) +
    ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle(paste("Top", top_genes, "Significant DE Genes"))
  print(p_top)
  dev.off()

  # 3. Volcano Plot
pdf(file= file.path(plot_dir, paste0("DE_Volcanoplot_",level,"_vs_",base,".pdf")))
suppressWarnings({
    vp <- EnhancedVolcano::EnhancedVolcano(
      results_data$res_tbl, lab = results_data$res_tbl$gene, 
      selectLab = highlight_genes, 
      drawConnectors = !is.null(highlight_genes), 
      x = 'log2FoldChange', y = 'padj', 
      title = paste(level, "vs", base), 
      pCutoff = padj_cutoff, FCcutoff = 1.0, pointSize = 2.0, labSize = 4.0, 
      col=c('grey','grey','grey', 'red2')
    )
    print(vp)
  })
  dev.off()
}