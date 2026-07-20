# ==============================================================================
# mod_dge.R - Gene-level DGE functions with enhanced Entrez mapping
# ==============================================================================

#' Write structured DGE parameters log
#' @keywords internal
.write_dge_log <- function(out_dir, comp_name, model, test, reduced,
                           level, base, shrink_method, padj_cutoff,
                           n_genes_input, n_genes_after_filter,
                           n_de_genes_up, n_de_genes_down,
                           filter_criterion = "rowSums(counts) >= 1") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not installed; writing log as text file.")
    log_dir <- file.path(out_dir, "Log", "DGE")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, paste0("DGE_params_", comp_name, ".txt"))
    con <- file(log_file, open = "wt")
    on.exit(close(con), add = TRUE)
    cat(paste("Timestamp:", Sys.time(), "\n"), file = con)
    cat(paste("Comparison:", comp_name, "\n"), file = con)
    cat(paste("Design formula:", model, "\n"), file = con)
    cat(paste("Test type:", test, "\n"), file = con)
    if (test == "LRT") cat(paste("Reduced formula:", reduced, "\n"), file = con)
    cat(paste("Contrast:", level, "vs", base, "\n"), file = con)
    cat(paste("Shrinkage method:", shrink_method, "\n"), file = con)
    cat(paste("padj cutoff:", padj_cutoff, "\n"), file = con)
    cat(paste("Filter criterion:", filter_criterion, "\n"), file = con)
    cat(paste("Genes before filter:", n_genes_input, "\n"), file = con)
    cat(paste("Genes after filter:", n_genes_after_filter, "\n"), file = con)
    cat(paste("DE genes up (log2FC>1):", n_de_genes_up, "\n"), file = con)
    cat(paste("DE genes down (log2FC<-1):", n_de_genes_down, "\n"), file = con)
    return()
  }

  main_condition <- tail(all.vars(as.formula(model)), 1)
  params <- list(
    timestamp          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    comparison         = comp_name,
    design_formula     = model,
    test_type          = test,
    reduced_formula    = if (!is.null(reduced)) reduced else "none",
    contrast           = list(factor = main_condition, level = level, base = base),
    shrinkage_method   = shrink_method,
    padj_cutoff        = padj_cutoff,
    filter_criterion   = filter_criterion,
    genes_before_filter = n_genes_input,
    genes_after_filter  = n_genes_after_filter,
    de_genes_up        = n_de_genes_up,
    de_genes_down      = n_de_genes_down,
    de_genes_total     = n_de_genes_up + n_de_genes_down
  )
  log_dir  <- file.path(out_dir, "Log", "DGE")
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
#' @param custom_tx2gene Optional path to a custom tx2gene file (TSV with columns 'tx_id' and 'gene_id')
#' @param custom_gene_map Optional path to a custom gene annotation file (TSV with columns 'gene_id', 'symbol', and optionally 'entrezid')
#' @return A list containing txi (or counts), meta, edb, gene_map, and type.
#' @export
import_counts <- function(data_dir, sample_table, ensembl_package_name, count_type = "salmon",
                          out_dir, matrix_file = NULL, subset_sample = NULL, remove_sample = NULL,
                          custom_tx2gene = NULL, custom_gene_map = NULL) {

  if (!requireNamespace("ensembldb", quietly = TRUE)) {
    stop("Package 'ensembldb' is required but not installed.")
  }

  if (!file.exists(sample_table)) stop("Sample table not found: ", sample_table)
  sample_df  <- data.table::fread(sample_table, header = TRUE, data.table = FALSE)
  sample_col <- if ("Sample" %in% colnames(sample_df)) "Sample" else "sample_id"

  sample_df <- .apply_sample_filters(sample_df, sample_col, remove_sample, subset_sample)

  rownames(sample_df) <- sample_df[[sample_col]]
  edb     <- getExportedValue(ensembl_package_name, ensembl_package_name)

  # ---- Build tx2gene (custom or from Ensembl) ----
  if (!is.null(custom_tx2gene) && file.exists(custom_tx2gene)) {
    message("Using custom tx2gene file: ", custom_tx2gene)
    tx2gene <- data.table::fread(custom_tx2gene, header = TRUE, data.table = FALSE)
    if (!all(c("tx_id", "gene_id") %in% colnames(tx2gene))) {
      stop("Custom tx2gene must contain columns 'tx_id' and 'gene_id'")
    }
    tx2gene <- tx2gene[, c("tx_id", "gene_id")]
    tx2gene$tx_id <- strip_ensembl_version(tx2gene$tx_id)
    tx2gene$gene_id <- strip_ensembl_version(tx2gene$gene_id)
  } else {
    tx2gene <- ensembldb::transcripts(edb, columns = c("tx_id", "gene_id"), return.type = "DataFrame")
    tx2gene <- as.data.frame(tx2gene)
    tx2gene$tx_id <- strip_ensembl_version(tx2gene$tx_id)
    tx2gene$gene_id <- strip_ensembl_version(tx2gene$gene_id)
  }

  # ---- Build gene_map (custom or from Ensembl) ----
  org_info <- get_organism_info(edb)
  org_db   <- org_info$org_db
  org_obj  <- if (requireNamespace(org_db, quietly = TRUE)) .load_org_db(org_db) else NULL

  if (!is.null(custom_gene_map) && file.exists(custom_gene_map)) {
    message("Using custom gene annotation file: ", custom_gene_map)
    gene_map <- data.table::fread(custom_gene_map, header = TRUE, data.table = FALSE)
    # Accept either 'gene_id' or 'ensembl' as ID column
    if (!("gene_id" %in% colnames(gene_map)) && "ensembl" %in% colnames(gene_map)) {
      colnames(gene_map)[colnames(gene_map) == "ensembl"] <- "gene_id"
    }
    if (!all(c("gene_id", "symbol") %in% colnames(gene_map))) {
      stop("Custom gene map must contain columns 'gene_id' (or 'ensembl') and 'symbol'")
    }
    # Strip version and rename
    gene_map$gene_id <- strip_ensembl_version(gene_map$gene_id)
    colnames(gene_map)[colnames(gene_map) == "gene_id"] <- "ensembl"
    if (!"entrezid" %in% colnames(gene_map)) {
      gene_map$entrezid <- NA_character_
    }
    # Keep only required columns
    gene_map <- gene_map[, c("ensembl", "symbol", "entrezid")]
    # Remove duplicates
    gene_map <- gene_map[!duplicated(gene_map$ensembl), ]
    # Fill missing symbol with ensembl
    gene_map$symbol[is.na(gene_map$symbol) | gene_map$symbol == ""] <- gene_map$ensembl[is.na(gene_map$symbol) | gene_map$symbol == ""]

    # ---- ENHANCED FALLBACK: fill missing Entrez using bitr ----
    if (!is.null(org_obj)) {
      gene_map <- .fill_entrez_with_bitr(gene_map, org_obj, id_col = "ensembl", symbol_col = "symbol")
    }
    
    # ---- FINAL CHECK: report Entrez coverage ----
    entrez_present <- sum(!is.na(gene_map$entrezid) & gene_map$entrezid != "")
    message("  Gene map loaded: ", nrow(gene_map), " genes, ", entrez_present, " with Entrez IDs")
    
  } else {
    # No custom gene map: build from Ensembl
    gene_map <- ensembldb::genes(edb, columns = c("gene_id", "gene_name"), return.type = "DataFrame")
    gene_map <- as.data.frame(gene_map)
    colnames(gene_map) <- c("ensembl", "symbol")
    gene_map$ensembl <- strip_ensembl_version(gene_map$ensembl)
    if (!is.null(org_obj)) {
      mapped_entrez <- suppressMessages(
        AnnotationDbi::mapIds(org_obj,
                              keys = gene_map$ensembl,
                              column = "ENTREZID",
                              keytype = "ENSEMBL",
                              multiVals = "first")
      )
      gene_map$entrezid <- as.character(mapped_entrez)
    } else {
      gene_map$entrezid <- NA_character_
    }
  }

  # --- Continue with tximport or matrix import ---
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

    if (count_type == "rsem")
      message("NOTE: rsem mode imports gene-level counts (quant.genes.results). For isoform analysis use isoforms.results.")
    missing_files <- tximport_file_list[!file.exists(tximport_file_list)]
    if (length(missing_files) > 0) {
      stop("Missing quantification files for samples: ", paste(names(missing_files), collapse = ", "),
           "\nExpected file like: ", count_file_name)
    }

    message("Importing ", count_type, " files via tximport...")
    txi <- tximport::tximport(
      tximport_file_list,
      type                 = count_type,
      tx2gene              = tx2gene,
      ignoreTxVersion      = TRUE,
      countsFromAbundance  = "lengthScaledTPM"
    )

    meta <- sample_df[colnames(txi$counts), , drop = FALSE]
    return(list(txi = txi, meta = meta, edb = edb, gene_map = gene_map, type = "tximport"))

  } else if (count_type == "matrix") {
    if (is.null(matrix_file) || !file.exists(matrix_file)) stop("Matrix file not specified or found.")

    counts_df <- data.table::fread(matrix_file, data.table = FALSE)
    rownames(counts_df) <- counts_df[, 1]
    counts_df <- counts_df[, -1, drop = FALSE]

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
#' @param tx_data List returned by import_counts
#' @param level Target level (can be NULL if model == "~1")
#' @param base Reference level (can be NULL if model == "~1")
#' @param model Design formula as a string (e.g., "~1" for EDA only)
#' @param replicate_col Column containing replicate information
#' @return DESeqDataSet object
#' @export
create_dds_object <- function(tx_data, level, base, model, replicate_col) {

  meta <- tx_data$meta
  design_formula <- as.formula(model)

  # If design is ~1, we skip factor conversions and releveling
  if (identical(design_formula, ~1)) {
    dds <- if (tx_data$type == "tximport") {
      DESeq2::DESeqDataSetFromTximport(tx_data$txi, colData = meta, design = ~1)
    } else {
      DESeq2::DESeqDataSetFromMatrix(countData = round(tx_data$counts),
                                     colData = meta, design = ~1)
    }
    keep <- rowSums(DESeq2::counts(dds)) >= 1
    dds <- dds[keep, ]
    return(dds)
  }

  # Normal path (model has factors)
  design_vars <- all.vars(design_formula)
  for (var in design_vars) {
    if (var %in% colnames(meta) && is.character(meta[[var]])) {
      meta[[var]] <- as.factor(meta[[var]])
    }
  }
  main_condition <- tail(design_vars, 1)
  meta[[main_condition]] <- as.factor(meta[[main_condition]])

  # Build dds
  dds <- if (tx_data$type == "tximport") {
    DESeq2::DESeqDataSetFromTximport(tx_data$txi, colData = meta, design = design_formula)
  } else {
    DESeq2::DESeqDataSetFromMatrix(countData = round(tx_data$counts),
                                   colData = meta, design = design_formula)
  }

  # Only relevel if level and base are provided (non-NULL)
  if (!is.null(level) && !is.null(base)) {
    dds[[main_condition]] <- relevel(dds[[main_condition]], base)
  } else {
    message("Note: level/base not provided; skipping releveling.")
  }

  # Minimal filter
  keep <- rowSums(DESeq2::counts(dds)) >= 1
  dds <- dds[keep, ]

  if (!is.null(replicate_col) && replicate_col %in% colnames(meta)) {
    dds <- DESeq2::collapseReplicates(dds, groupby = dds[[replicate_col]])
  }

  return(dds)
}

#' Run DESeq2 Analysis
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
#' @export
run_deseq2_analysis <- function(dds, model, level, base, shrink_method, out_dir,
                                padj_cutoff, test = "Wald", reduced = NULL) {

  n_genes_input <- nrow(dds)

  if (test == "LRT") {
    if (is.null(reduced)) stop("You must provide a reduced model for LRT.")
    if (shrink_method != "ashr") {
      message("NOTE: LRT is only compatible with shrink_method = 'ashr'. Overriding to 'ashr'.")
      shrink_method <- "ashr"
    }
    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = as.formula(reduced))
  } else {
    dds <- DESeq2::DESeq(dds)
  }

  main_condition <- tail(all.vars(as.formula(model)), 1)

  if (test == "LRT") {
    res_unshrunken <- DESeq2::results(dds, alpha = padj_cutoff)
    res_shrunken   <- suppressMessages(DESeq2::lfcShrink(dds, res = res_unshrunken, type = "ashr"))
    res_shrunken$stat <- res_unshrunken$stat
  } else {
    comparison     <- c(main_condition, level, base)
    res_unshrunken <- DESeq2::results(dds, contrast = comparison, alpha = padj_cutoff)
    res_shrunken   <- DESeq2::lfcShrink(dds, contrast = comparison, type = shrink_method)
    res_shrunken$stat <- res_unshrunken$stat
  }

  res_df      <- as.data.frame(res_shrunken)
  sig_up      <- sum(res_df$padj < padj_cutoff & res_df$log2FoldChange > 1,  na.rm = TRUE)
  sig_down    <- sum(res_df$padj < padj_cutoff & res_df$log2FoldChange < -1, na.rm = TRUE)
  n_genes_after_filter <- nrow(dds)

  comp_name <- paste0(level, "_vs_", base)

  base_dir <- file.path(out_dir, "DE_raw_results")
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  .write_dge_log(out_dir, comp_name, model, test, reduced, level, base,
                 shrink_method, padj_cutoff, n_genes_input, n_genes_after_filter,
                 sig_up, sig_down)

  return(list(dds = dds, res_unshrunken = res_unshrunken, res_shrunken = res_shrunken))
}

#' Export Significant Results
#' @param res_shrunken Shrunken DESeq results
#' @param res_unshrunken Unshrunken DESeq results
#' @param dds DESeqDataSet object
#' @param out_dir Output directory
#' @param level Comparison level
#' @param base Base level
#' @param gene_map Gene symbol mapping dataframe
#' @param padj_cutoff Adjusted p-value cutoff
#' @return List of processed dataframes
#' @export
export_significant_results <- function(res_shrunken, res_unshrunken, dds, out_dir,
                                       level, base, gene_map, padj_cutoff) {
  fc_dir <- file.path(out_dir, "DE_raw_results")
  if (!dir.exists(fc_dir)) dir.create(fc_dir, recursive = TRUE)

  # --- MAIN RESULTS TABLE ---
  res_tbl <- as.data.frame(res_shrunken)
  # Strip version suffixes from Ensembl IDs BEFORE merging (MSTRG/other custom
  # IDs containing dots are left untouched, see strip_ensembl_version())
  res_tbl$ensembl <- strip_ensembl_version(rownames(res_tbl))
  res_tbl <- merge(res_tbl, gene_map, by = "ensembl", all.x = TRUE)
  
  # --- RECOVER ENTREZ IDs FROM GENE MAP (fallback for any that were missed) ---
  if (any(is.na(res_tbl$entrezid) | res_tbl$entrezid == "")) {
    message("  Recovering missing Entrez IDs from gene_map...")
    
    # Clean gene_map for matching
    gene_map_clean <- gene_map[!is.na(gene_map$entrezid) & gene_map$entrezid != "", ]
    gene_map_clean$ensembl_clean <- strip_ensembl_version(gene_map_clean$ensembl)
    
    # First pass: match by gene symbol already resolved on res_tbl (i.e. the
    # primary merge succeeded on "ensembl" but entrezid itself was blank in
    # gene_map for that row)
    missing_idx <- is.na(res_tbl$entrezid) | res_tbl$entrezid == ""
    if (any(missing_idx)) {
      for (i in which(missing_idx)) {
        sym <- res_tbl$symbol[i]
        if (!is.na(sym) && sym != "" && sym != res_tbl$ensembl[i]) {
          match_idx <- which(gene_map_clean$symbol == sym)
          if (length(match_idx) > 0) {
            res_tbl$entrezid[i] <- gene_map_clean$entrezid[match_idx[1]]
          }
        }
      }
    }
    
    # Second pass: match by Ensembl ID (version-stripped). Mostly redundant
    # with the primary merge, kept as a safety net.
    missing_idx <- is.na(res_tbl$entrezid) | res_tbl$entrezid == ""
    if (any(missing_idx)) {
      for (i in which(missing_idx)) {
        ens <- res_tbl$ensembl[i]
        if (!is.na(ens) && ens != "") {
          match_idx <- which(gene_map_clean$ensembl_clean == ens)
          if (length(match_idx) > 0) {
            res_tbl$entrezid[i] <- gene_map_clean$entrezid[match_idx[1]]
          }
        }
      }
    }

    # Third pass (the one that was MISSING): res_tbl$ensembl may itself be a
    # bare gene SYMBOL rather than a gene_id, because a custom_tx2gene file
    # can use gene symbols as its aggregation key (e.g. tximport's tx2gene
    # maps "ENST00000832824" -> "DDX11L16" directly, instead of an
    # ENSG/MSTRG id). In that case the two passes above never fire because
    # res_tbl$symbol is still NA (the primary "by = ensembl" merge found no
    # match at all) and res_tbl$ensembl doesn't look like an Ensembl ID
    # either. Try matching res_tbl$ensembl against gene_map's symbol column.
    missing_idx <- is.na(res_tbl$entrezid) | res_tbl$entrezid == ""
    if (any(missing_idx)) {
      gene_map_by_symbol <- gene_map_clean[!duplicated(gene_map_clean$symbol), ]
      for (i in which(missing_idx)) {
        ens <- res_tbl$ensembl[i]
        if (!is.na(ens) && ens != "") {
          match_idx <- which(gene_map_by_symbol$symbol == ens)
          if (length(match_idx) > 0) {
            res_tbl$entrezid[i] <- gene_map_by_symbol$entrezid[match_idx[1]]
            if (is.na(res_tbl$symbol[i]) || res_tbl$symbol[i] == "") {
              res_tbl$symbol[i] <- ens
            }
          }
        }
      }
    }
    
    recovered <- sum(!is.na(res_tbl$entrezid) & res_tbl$entrezid != "")
    total <- nrow(res_tbl)
    message("    Recovered ", recovered, "/", total, " Entrez IDs")
  }

  # Set gene column for display
  res_tbl$gene <- ifelse(is.na(res_tbl$symbol) | res_tbl$symbol == "", res_tbl$ensembl, res_tbl$symbol)

  # --- NORMALIZED COUNTS ---
  counts_df <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  counts_df$ensembl <- strip_ensembl_version(rownames(counts_df))
  counts_df <- merge(counts_df, gene_map, by = "ensembl", all.x = TRUE)
  counts_df$gene <- ifelse(is.na(counts_df$symbol) | counts_df$symbol == "", counts_df$ensembl, counts_df$symbol)
  
  tryCatch(
    utils::write.table(counts_df, file.path(fc_dir, paste0("DESeq_Normalized_Counts_", level, "_vs_", base, ".txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE),
    error = function(e) warning("Could not write normalized counts: ", e$message)
  )

  # --- RAW COUNTS ---
  raw_counts <- as.data.frame(DESeq2::counts(dds, normalized = FALSE))
  raw_counts$ensembl <- strip_ensembl_version(rownames(raw_counts))
  raw_counts <- merge(raw_counts, gene_map, by = "ensembl", all.x = TRUE)
  raw_counts$gene <- ifelse(is.na(raw_counts$symbol) | raw_counts$symbol == "", raw_counts$ensembl, raw_counts$symbol)
  
  tryCatch(
    utils::write.table(raw_counts, file.path(fc_dir, paste0("DESeq_Raw_Counts_", level, "_vs_", base, ".txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE),
    error = function(e) warning("Could not write raw counts: ", e$message)
  )

  # --- FILTER RESULTS TABLE ---
  desired_cols <- c("gene", "ensembl", "entrezid", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "stat")
  res_cols <- intersect(desired_cols, colnames(res_tbl))
  res_tbl <- res_tbl[, res_cols, drop = FALSE]

  # --- SIGNIFICANT RESULTS ---
  sig_res <- res_tbl[!is.na(res_tbl$padj) & res_tbl$padj < padj_cutoff, ]

  # --- WRITE FILES ---
  utils::write.table(res_tbl, file.path(fc_dir, paste0("DEgenes_raw_", level, "_vs_", base, ".txt")),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  write_filtered <- function(df, tag) {
    prefix <- file.path(fc_dir, paste0("DEgenes_", tag, "_", level, "_vs_", base))
    utils::write.table(df,                                   paste0(prefix, ".txt"),      sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(df[df$log2FoldChange > 1,  ],         paste0(prefix, "_LFC1.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(df[df$log2FoldChange < -1, ],         paste0(prefix, "_LFC-1.txt"),sep = "\t", quote = FALSE, row.names = FALSE)
  }

  padj_tag <- paste0("pval_", gsub("\\.", "_", as.character(padj_cutoff)))
  write_filtered(sig_res, padj_tag)

  return(list(res_tbl = res_tbl, sig_res = sig_res, normalized_counts = counts_df, raw_counts = raw_counts))
}