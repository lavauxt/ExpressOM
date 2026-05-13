#' Import RNA-seq Counts and Prepare Metadata
#' 
#' @param data_dir Folder where the input data is stored.
#' @param sample_table Path to the sample table.
#' @param ensembl_package_name Name of the installed Ensembl database package.
#' @param count_type Type of RNA count (e.g., "salmon" or "matrix").
#' @param out_dir Output directory.
#' @param matrix_file Path to raw counts file if count_type is "matrix".
#' @return A list containing txi (or counts), meta, edb, gene_map, and type.
#' @export
import_counts <- function(data_dir, sample_table, ensembl_package_name, count_type = "salmon", out_dir, matrix_file = NULL) {
  
  # 1. Load Sample Table
  if (!file.exists(sample_table)) stop("Sample table not found: ", sample_table)
  sample_df <- data.table::fread(sample_table, header = TRUE, data.table = FALSE)
  # Handle varying sample ID column name conventions
  sample_col <- if ("Sample" %in% colnames(sample_df)) "Sample" else "sample_id"
  rownames(sample_df) <- sample_df[[sample_col]]
  
  # 2. Access the Homemade Database
  edb <- getExportedValue(ensembl_package_name, ensembl_package_name)
  
  # 3. Build tx2gene map
  tx2gene <- ensembldb::transcripts(edb, columns = c("tx_id", "gene_id"), return.type = "DataFrame")
  tx2gene <- as.data.frame(tx2gene)
  
  # 4. Build Detailed Gene Map
  gene_map <- ensembldb::genes(edb, columns = c("gene_id", "gene_name"), return.type = "DataFrame")
  gene_map <- as.data.frame(gene_map)
  colnames(gene_map) <- c("ensembl", "symbol")
  
  # Fetch organism info to map Entrez IDs using org.db
  org_info <- get_organism_info(edb)
  org_db <- org_info$org_db
  if (requireNamespace(org_db, quietly = TRUE)) {
    org_obj <- getExportedValue(org_db, org_db)
    mapped_entrez <- suppressMessages(AnnotationDbi::mapIds(
      org_obj, 
      keys = gene_map$ensembl, 
      column = "ENTREZID", 
      keytype = "ENSEMBL", 
      multiVals = "first"
    ))
    gene_map$entrezid <- as.character(mapped_entrez)
  } else {
    gene_map$entrezid <- NA
  }

  if (count_type == "salmon") {
    # 5. Import via tximport
    tximport_file_list <- file.path(data_dir, sample_df[[sample_col]], 
                                   paste0(sample_df[[sample_col]], ".", count_type, ".quant"), "quant.sf")
    names(tximport_file_list) <- sample_df[[sample_col]]

    txi <- tximport::tximport(
      tximport_file_list, 
      type = count_type, 
      tx2gene = tx2gene, 
      ignoreTxVersion = TRUE,
      countsFromAbundance = "lengthScaledTPM"
    )

    # 6. Align Metadata
    meta <- sample_df[colnames(txi$counts), , drop = FALSE]
    return(list(txi = txi, meta = meta, edb = edb, gene_map = gene_map, type = "salmon"))
    
  } else if (count_type == "matrix") {
    # 5. Process Raw Matrix
    if (is.null(matrix_file) || !file.exists(matrix_file)) stop("Matrix file not specified or found.")
    
    counts_df <- data.table::fread(matrix_file, data.table = FALSE)
    # Assume first column is Gene ID if no header, else use first column
    rownames(counts_df) <- counts_df[, 1]
    
    # Filter columns to exactly match metadata sample IDs
    valid_samples <- intersect(colnames(counts_df), rownames(sample_df))
    if (length(valid_samples) == 0) stop("No matching sample names between matrix and sample table.")
    
    count_mat <- as.matrix(counts_df[, valid_samples, drop = FALSE])
    # Ensure mode is numeric
    mode(count_mat) <- "numeric"
    # Overwrite NA/NaN with 0 if necessary
    count_mat[is.na(count_mat)] <- 0
    
    meta <- sample_df[valid_samples, , drop = FALSE]
    return(list(counts = count_mat, meta = meta, edb = edb, gene_map = gene_map, type = "matrix"))
  } else {
    stop("Unsupported count_type: ", count_type)
  }
}

#' Create DESeqDataSet object
#' 
#' @export
#' @param tx_data List returned by import_counts
#' @param level Target level
#' @param base Reference level
#' @param model Design formula as a string
#' @param replicate_col Column containing replicate information
#' @return DESeqDataSet object
create_dds_object <- function(tx_data, level, base, model, replicate_col) {
  
  meta <- tx_data$meta
  main_condition <- tail(all.vars(as.formula(model)), 1)
  meta[[main_condition]] <- as.factor(meta[[main_condition]])
  
  if (tx_data$type == "salmon") {
    dds <- DESeq2::DESeqDataSetFromTximport(tx_data$txi, colData = meta, design = as.formula(model))
  } else if (tx_data$type == "matrix") {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(tx_data$counts), colData = meta, design = as.formula(model))
  }
  
  dds[[main_condition]] <- relevel(dds[[main_condition]], base)
  
  keep <- rowSums(DESeq2::counts(dds)) >= 1
  dds <- dds[keep,]
  
  if (!is.null(replicate_col) && replicate_col %in% colnames(meta)) {
    dds <- DESeq2::collapseReplicates(dds, groupby = dds[[replicate_col]])
  }
  
  return(dds)
}

#' Run DESeq2 Analysis
#' @export
#' @param dds DESeqDataSet object
#' @param model Design formula constraint
#' @param level the active group
#' @param base the reference group
#' @param shrink_method Shrinkage string
#' @param out_dir Output path
#' @param padj_cutoff Adjusted p-value cutoff
#' @param test Type of statistical test ("Wald" or "LRT")
#' @param reduced Reduced formula for LRT test
#' @return List containing res_unshrunken, res_shrunken, and the computed dds
run_deseq2_analysis <- function(dds, model, level, base, shrink_method, out_dir, padj_cutoff, test = "Wald", reduced = NULL) {
  
  if (test == "LRT") {
    if (is.null(reduced)) stop("You must provide a reduced model for LRT.")
    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = as.formula(reduced))
  } else {
    dds <- DESeq2::DESeq(dds)
  }
  
  main_condition <- tail(all.vars(as.formula(model)), 1)
  
  log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  comparison_raw <- DESeq2::resultsNames(dds)[-1]
  write.table(comparison_raw, file.path(log_dir, "Possible_Models.txt"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  
  if (test == "LRT") {
    write.table("LRT Test", file.path(log_dir, "Final_Model.txt"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
    res_unshrunken <- DESeq2::results(dds, alpha = padj_cutoff)
    res_shrunken <- suppressMessages(DESeq2::lfcShrink(dds, res = res_unshrunken, type = "ashr"))
  } else {
    comparison <- c(main_condition, level, base)
    write.table(as.character(comparison), file.path(log_dir, "Final_Model.txt"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
    res_unshrunken <- DESeq2::results(dds, contrast = comparison, alpha = padj_cutoff)
    res_shrunken <- DESeq2::lfcShrink(dds, coef = 2, type = shrink_method)
  }
  
  return(list(dds = dds, res_unshrunken = res_unshrunken, res_shrunken = res_shrunken))
}

#' Export Significant Results
#' @export
#' @param res_shrunken Shrunken DESeq results
#' @param res_unshrunken Unshrunken DESeq results
#' @param dds DESeqDataSet object
#' @param out_dir Output directory
#' @param level Comparison level
#' @param base Base level
#' @param gene_map Gene symbol mapping dataframe
#' @param padj_cutoff Adjusted p-value cutoff
#' @return List of processed dataframes
export_significant_results <- function(res_shrunken, res_unshrunken, dds, out_dir, level, base, gene_map, padj_cutoff) {
  fc_dir <- file.path(out_dir, "DE_Results")
  if (!dir.exists(fc_dir)) dir.create(fc_dir, recursive = TRUE)

  res_tbl <- as.data.frame(res_shrunken)
  res_tbl$ensembl <- rownames(res_tbl)
  res_tbl <- merge(res_tbl, gene_map, by = "ensembl", all.x = TRUE)
  res_tbl$gene <- ifelse(is.na(res_tbl$symbol) | res_tbl$symbol == "", res_tbl$ensembl, res_tbl$symbol)

  # Normalized counts
  counts_df <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  counts_df$ensembl <- rownames(counts_df)
  counts_df <- merge(counts_df, gene_map, by = "ensembl", all.x = TRUE)
  counts_df$gene <- ifelse(is.na(counts_df$symbol) | counts_df$symbol == "", counts_df$ensembl, counts_df$symbol)
  utils::write.table(counts_df, file.path(fc_dir, paste0("DESeq_Normalized_Counts_", level, "_vs_", base, ".txt")), sep="\t", quote=FALSE, row.names=FALSE)
  
  # Raw counts
  raw_counts <- as.data.frame(DESeq2::counts(dds, normalized = FALSE))
  raw_counts$ensembl <- rownames(raw_counts)
  raw_counts <- merge(raw_counts, gene_map, by = "ensembl", all.x = TRUE)
  raw_counts$gene <- ifelse(is.na(raw_counts$symbol) | raw_counts$symbol == "", raw_counts$ensembl, raw_counts$symbol)
  utils::write.table(raw_counts, file.path(fc_dir, paste0("DESeq_Raw_Counts_", level, "_vs_", base, ".txt")), sep="\t", quote=FALSE, row.names=FALSE)

  res_cols <- c("gene", "ensembl", "entrezid", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
  res_tbl <- res_tbl[, res_cols]
  sig_res <- res_tbl[!is.na(res_tbl$padj) & res_tbl$padj < padj_cutoff, ]
  
  utils::write.table(res_tbl, file.path(fc_dir, paste0("DEgenes_raw_", level, "_vs_", base, ".txt")), sep="\t", quote=FALSE, row.names = FALSE)
  
  write_filtered <- function(df, tag) {
    prefix <- file.path(fc_dir, paste0("DEgenes_", tag, "_", level, "_vs_", base))
    utils::write.table(df, paste0(prefix, ".txt"), sep="\t", quote=FALSE, row.names = FALSE)
    utils::write.table(df[df$log2FoldChange > 1, ], paste0(prefix, "_LFC1.txt"), sep="\t", quote=FALSE, row.names = FALSE)
    utils::write.table(df[df$log2FoldChange < -1, ], paste0(prefix, "_LFC-1.txt"), sep="\t", quote=FALSE, row.names = FALSE)
  }
  
  padj_tag <- paste0("pval_", gsub("\\.", "_", as.character(padj_cutoff)))
  write_filtered(sig_res, padj_tag)
  
  return(list(res_tbl = res_tbl, sig_res = sig_res, normalized_counts = counts_df, raw_counts = raw_counts))
}