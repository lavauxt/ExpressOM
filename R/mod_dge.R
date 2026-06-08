#' Write structured DGE parameters log
#' @keywords internal
.write_dge_log <- function(out_dir, comp_name, model, test, reduced, 
                           level, base, shrink_method, padj_cutoff,
                           n_genes_input, n_genes_after_filter,
                           n_de_genes_up, n_de_genes_down,
                           filter_criterion = "rowSums(counts) >= 1") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not installed; writing log as text file.")
    log_dir <- file.path(out_dir, "Log")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, paste0("DGE_params_", comp_name, ".txt"))
    sink(log_file)
    cat(paste("Timestamp:", Sys.time(), "\n"))
    cat(paste("Comparison:", comp_name, "\n"))
    cat(paste("Design formula:", model, "\n"))
    cat(paste("Test type:", test, "\n"))
    if (test == "LRT") cat(paste("Reduced formula:", reduced, "\n"))
    cat(paste("Contrast:", level, "vs", base, "\n"))
    cat(paste("Shrinkage method:", shrink_method, "\n"))
    cat(paste("padj cutoff:", padj_cutoff, "\n"))
    cat(paste("Filter criterion:", filter_criterion, "\n"))
    cat(paste("Genes before filter:", n_genes_input, "\n"))
    cat(paste("Genes after filter:", n_genes_after_filter, "\n"))
    cat(paste("DE genes up (log2FC>1):", n_de_genes_up, "\n"))
    cat(paste("DE genes down (log2FC<-1):", n_de_genes_down, "\n"))
    sink()
    return()
  }
  
  main_condition <- tail(all.vars(as.formula(model)), 1)
  params <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    comparison = comp_name,
    design_formula = model,
    test_type = test,
    reduced_formula = if (!is.null(reduced)) reduced else "none",
    contrast = list(factor = main_condition, level = level, base = base),
    shrinkage_method = shrink_method,
    padj_cutoff = padj_cutoff,
    filter_criterion = filter_criterion,
    genes_before_filter = n_genes_input,
    genes_after_filter = n_genes_after_filter,
    de_genes_up = n_de_genes_up,
    de_genes_down = n_de_genes_down,
    de_genes_total = n_de_genes_up + n_de_genes_down
  )
  log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  log_file <- file.path(log_dir, paste0("DGE_params_", comp_name, ".json"))
  jsonlite::write_json(params, log_file, pretty = TRUE, auto_unbox = TRUE)
  invisible(params)
}

#' Import RNA-seq Counts and Prepare Metadata
#' 
#' @param data_dir Folder where the input data is stored.
#' @param sample_table Path to the sample table.
#' @param ensembl_package_name Name of the installed Ensembl database package.
#' @param count_type Type of RNA count (e.g., "salmon" or "matrix").
#' @param out_dir Output directory.
#' @param matrix_file Path to raw counts file if count_type is "matrix".
#' @param subset_sample Optional string to filter the sample table.
#' @param remove_sample Optional character vector of sample IDs to exclude
#' @return A list containing txi (or counts), meta, edb, gene_map, and type.
#' @export
import_counts <- function(data_dir, sample_table, ensembl_package_name, count_type = "salmon", out_dir, matrix_file = NULL, subset_sample = NULL, remove_sample = NULL) {
  
  if (!file.exists(sample_table)) stop("Sample table not found: ", sample_table)
  sample_df <- data.table::fread(sample_table, header = TRUE, data.table = FALSE)

  sample_col <- if ("Sample" %in% colnames(sample_df)) "Sample" else "sample_id"

  if (!is.null(remove_sample)) {
    message("   -> Excluding requested samples: ", paste(remove_sample, collapse = ", "))
    keep_indices <- !(sample_df[[sample_col]] %in% remove_sample)
    sample_df <- sample_df[keep_indices, , drop = FALSE]
    if (nrow(sample_df) == 0) stop("The remove_sample constraint removed all available samples from your metadata!")
  }

  if (!is.null(subset_sample)) {
    message("   -> Applying subset condition: ", subset_sample)
    sample_df <- tryCatch({
      filter_expr <- rlang::parse_expr(subset_sample)
      subset_indices <- eval(filter_expr, envir = sample_df)
      sample_df[subset_indices, , drop = FALSE]
    }, error = function(e) {
      stop("Failed to evaluate subset_sample condition. Error: ", e$message)
    })
  }

  rownames(sample_df) <- sample_df[[sample_col]]
  edb <- getExportedValue(ensembl_package_name, ensembl_package_name)
  tx2gene <- ensembldb::transcripts(edb, columns = c("tx_id", "gene_id"), return.type = "DataFrame")
  tx2gene <- as.data.frame(tx2gene)
  gene_map <- ensembldb::genes(edb, columns = c("gene_id", "gene_name"), return.type = "DataFrame")
  gene_map <- as.data.frame(gene_map)
  colnames(gene_map) <- c("ensembl", "symbol")
  
  org_info <- get_organism_info(edb)
  org_db <- org_info$org_db
  if (requireNamespace(org_db, quietly = TRUE)) {
    org_obj <- .load_org_db(org_db)
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

  if (count_type != "matrix") {
    count_file_name <- switch(count_type,
      "salmon"    = "quant.sf",
      "kallisto"  = "abundance.tsv",
      "rsem"      = "quant.genes.results",
      "stringtie" = "t_data.ctab",
      stop("Unsupported count_type for tximport: ", count_type)
    )
    
    tximport_file_list <- sapply(sample_df[[sample_col]], function(sid) {
      p_nested <- file.path(data_dir, sid, paste0(sid, ".", count_type, ".quant"), count_file_name)
      p_direct <- file.path(data_dir, sid, count_file_name)
      if (file.exists(p_nested)) return(p_nested)
      return(p_direct)
    })
    
    names(tximport_file_list) <- sample_df[[sample_col]]

    missing_files <- tximport_file_list[!file.exists(tximport_file_list)]
    if (length(missing_files) > 0) {
      stop("Missing quantification files for samples: ", paste(names(missing_files), collapse = ", "),
           "\nExpected file like: ", count_file_name)
    }

    message("Importing ", count_type, " files via tximport...")
    txi <- tximport::tximport(
      tximport_file_list, 
      type = count_type, 
      tx2gene = tx2gene, 
      ignoreTxVersion = TRUE,
      countsFromAbundance = "lengthScaledTPM"
    )

    meta <- sample_df[colnames(txi$counts), , drop = FALSE]
    return(list(txi = txi, meta = meta, edb = edb, gene_map = gene_map, type = "tximport"))
    
  } else if (count_type == "matrix") {
    if (is.null(matrix_file) || !file.exists(matrix_file)) stop("Matrix file not specified or found.")
    
    counts_df <- data.table::fread(matrix_file, data.table = FALSE)
    rownames(counts_df) <- counts_df[, 1]
    
    valid_samples <- intersect(colnames(counts_df), rownames(sample_df))
    if (length(valid_samples) == 0) stop("No matching sample names between matrix and sample table.")
    
    count_mat <- as.matrix(counts_df[, valid_samples, drop = FALSE])
    mode(count_mat) <- "numeric"
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

  design_vars <- all.vars(as.formula(model))
    for (var in design_vars) {
      if (var %in% colnames(meta) && is.character(meta[[var]])) {
        meta[[var]] <- as.factor(meta[[var]])
      }
    }

  main_condition <- tail(all.vars(as.formula(model)), 1)
  meta[[main_condition]] <- as.factor(meta[[main_condition]])
  
  if (tx_data$type == "tximport") {
    dds <- DESeq2::DESeqDataSetFromTximport(tx_data$txi, colData = meta, design = as.formula(model))
  } else if (tx_data$type == "matrix") {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(tx_data$counts), colData = meta, design = as.formula(model))
  } else {
    stop("Invalid tx_data type: ", tx_data$type)
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
  
  # Store number of genes before any filtering (for logging)
  n_genes_input <- nrow(dds)
  
  if (test == "LRT") {
    if (is.null(reduced)) stop("You must provide a reduced model for LRT.")
    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = as.formula(reduced))
  } else {
    dds <- DESeq2::DESeq(dds)
  }
  
  main_condition <- tail(all.vars(as.formula(model)), 1)

  if (test == "LRT") {
    res_unshrunken <- DESeq2::results(dds, alpha = padj_cutoff)
    res_shrunken   <- suppressMessages(DESeq2::lfcShrink(dds, res = res_unshrunken, type = shrink_method))
    res_shrunken$stat <- res_unshrunken$stat
  } else {
    comparison     <- c(main_condition, level, base)
    res_unshrunken <- DESeq2::results(dds, contrast = comparison, alpha = padj_cutoff)
    res_shrunken <- DESeq2::lfcShrink(dds, contrast = comparison, type = shrink_method)
    res_shrunken$stat <- res_unshrunken$stat
  }
  
  # After obtaining results, compute DE counts for logging
  res_df <- as.data.frame(res_shrunken)
  sig_up <- sum(res_df$padj < padj_cutoff & res_df$log2FoldChange > 1, na.rm = TRUE)
  sig_down <- sum(res_df$padj < padj_cutoff & res_df$log2FoldChange < -1, na.rm = TRUE)
  n_genes_after_filter <- nrow(dds)
  
  # Write structured log (replaces old text logs)
  comp_name <- paste0(level, "_vs_", base)

base_dir <- file.path(out_dir, "DE_raw_results")
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

.write_dge_log(base_dir, comp_name, model, test, reduced, level, base, 
               shrink_method, padj_cutoff, n_genes_input, n_genes_after_filter,
               sig_up, sig_down)
  
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
  fc_dir <- file.path(out_dir, "DE_raw_results")
  if (!dir.exists(fc_dir)) dir.create(fc_dir, recursive = TRUE)

  log_dir <- file.path(fc_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  res_tbl <- as.data.frame(res_shrunken)
  res_tbl$ensembl <- rownames(res_tbl)
  res_tbl <- merge(res_tbl, gene_map, by = "ensembl", all.x = TRUE)
  res_tbl$gene <- ifelse(is.na(res_tbl$symbol) | res_tbl$symbol == "", res_tbl$ensembl, res_tbl$symbol)

  counts_df <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  counts_df$ensembl <- rownames(counts_df)
  counts_df <- merge(counts_df, gene_map, by = "ensembl", all.x = TRUE)
  counts_df$gene <- ifelse(is.na(counts_df$symbol) | counts_df$symbol == "", counts_df$ensembl, counts_df$symbol)
  utils::write.table(counts_df, file.path(fc_dir, paste0("DESeq_Normalized_Counts_", level, "_vs_", base, ".txt")), sep="\t", quote=FALSE, row.names=FALSE)
  
  raw_counts <- as.data.frame(DESeq2::counts(dds, normalized = FALSE))
  raw_counts$ensembl <- rownames(raw_counts)
  raw_counts <- merge(raw_counts, gene_map, by = "ensembl", all.x = TRUE)
  raw_counts$gene <- ifelse(is.na(raw_counts$symbol) | raw_counts$symbol == "", raw_counts$ensembl, raw_counts$symbol)
  utils::write.table(raw_counts, file.path(fc_dir, paste0("DESeq_Raw_Counts_", level, "_vs_", base, ".txt")), sep="\t", quote=FALSE, row.names=FALSE)

  res_cols <- c("gene", "ensembl", "entrezid", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "stat")
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