
convert_pdf_to_png <- function(pdf_file, dpi = 200) {

  if (!file.exists(pdf_file)) {
    warning("PDF not found: ", pdf_file)
    return(NULL)
  }

  png_file <- sub("\\.pdf$", ".png", pdf_file)

  pdftools::pdf_convert(
    pdf = pdf_file,
    filenames = png_file,
    dpi = dpi
  )

  if (file.exists(png_file)) {
    return(normalizePath(
      png_file,
      winslash = "/",
      mustWork = FALSE
    ))
  }

  return(NULL)
}


# ==============================================================================
# mod_isoform.R - DTE, DTU, and IsoformSwitchAnalyzeR integration
# ==============================================================================

#' Import transcript-level counts for isoform analysis
#'
#' @param data_dir Folder with quantification files
#' @param sample_table Path to sample metadata
#' @param ensembl_package_name EnsDb package name
#' @param count_type Quantification type ("salmon", "kallisto", etc.)
#' @param matrix_file Path to raw counts matrix (if count_type = "matrix")
#' @param subset_sample Optional filter expression
#' @param remove_sample Optional sample IDs to exclude
#' @return List with txi (transcript counts), meta, tx2gene, gene_map
#' @export
import_transcript_counts <- function(data_dir, sample_table, ensembl_package_name,
                                      count_type = "salmon", matrix_file = NULL,
                                      subset_sample = NULL, remove_sample = NULL) {
  
  if (!file.exists(sample_table)) stop("Sample table not found: ", sample_table)
  sample_df <- data.table::fread(sample_table, header = TRUE, data.table = FALSE)
  sample_col <- if ("Sample" %in% colnames(sample_df)) "Sample" else "sample_id"
  
  if (!is.null(remove_sample)) {
    sample_df <- sample_df[!(sample_df[[sample_col]] %in% remove_sample), , drop = FALSE]
    if (nrow(sample_df) == 0) stop("All samples removed by remove_sample")
  }
  
  if (!is.null(subset_sample)) {
    filter_expr <- rlang::parse_expr(subset_sample)
    subset_indices <- eval(filter_expr, envir = sample_df)
    sample_df <- sample_df[subset_indices, , drop = FALSE]
  }
  
  rownames(sample_df) <- sample_df[[sample_col]]
  edb <- getExportedValue(ensembl_package_name, ensembl_package_name)
  tx2gene <- ensembldb::transcripts(edb, columns = c("tx_id", "gene_id"), return.type = "DataFrame")
  tx2gene <- as.data.frame(tx2gene)
  # --- CRITICAL: strip version suffixes from tx2gene ONCE ---
  tx2gene$tx_id <- sub("\\..*$", "", tx2gene$tx_id)
  
  gene_map <- ensembldb::genes(edb, columns = c("gene_id", "gene_name"), return.type = "DataFrame")
  gene_map <- as.data.frame(gene_map)
  colnames(gene_map) <- c("ensembl", "symbol")

  if (count_type != "matrix") {
    count_file_name <- switch(count_type,
                              "salmon" = "quant.sf",
                              "kallisto" = "abundance.tsv",
                              "rsem" = "quant.genes.results",
                              "stringtie" = "t_data.ctab",
                              stop("Unsupported count_type for tximport: ", count_type))
    
    file_list <- sapply(sample_df[[sample_col]], function(sid) {
      p_nested <- file.path(data_dir, sid, paste0(sid, ".", count_type, ".quant"), count_file_name)
      p_direct <- file.path(data_dir, sid, count_file_name)
      if (file.exists(p_nested)) return(p_nested)
      return(p_direct)
    })
    names(file_list) <- sample_df[[sample_col]]
    
    txi <- tximport::tximport(file_list, type = count_type, txOut = TRUE, countsFromAbundance = "lengthScaledTPM")
    # --- CRITICAL: strip version suffixes from transcript IDs in the count matrix ---
    rownames(txi$counts) <- sub("\\..*$", "", rownames(txi$counts))
    meta <- sample_df[colnames(txi$counts), , drop = FALSE]
    return(list(txi = txi, meta = meta, tx2gene = tx2gene, gene_map = gene_map, type = "tximport"))
    
  } else {
    if (is.null(matrix_file)) stop("matrix_file required for count_type='matrix'")
    counts_df <- data.table::fread(matrix_file, data.table = FALSE)
    rownames(counts_df) <- counts_df[, 1]   
    counts_df <- counts_df[, -1, drop = FALSE]
    valid_samples <- intersect(colnames(counts_df), rownames(sample_df))
    count_mat <- as.matrix(counts_df[, valid_samples, drop = FALSE])
    mode(count_mat) <- "numeric"
    count_mat[is.na(count_mat)] <- 0
    # --- CRITICAL: strip version suffixes from matrix row names ---
    rownames(count_mat) <- sub("\\..*$", "", rownames(count_mat))
    meta <- sample_df[valid_samples, , drop = FALSE]
    return(list(counts = count_mat, meta = meta, tx2gene = tx2gene, gene_map = gene_map, type = "matrix"))
  }
}

#' Run DTE (Differential Transcript Expression) using DESeq2
#'
#' @param isoform_obj Output from import_transcript_counts()
#' @param condition Column name in metadata containing groups
#' @param level Foreground level for contrast
#' @param base Reference level
#' @param padj_cutoff Adjusted p-value threshold
#' @return Data frame of DTE results
#' @export
run_dte <- function(isoform_obj, condition, level, base, padj_cutoff = 0.05) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for DTE analysis. Please install it.")
  }
  if (isoform_obj$type == "tximport") {
    dds <- DESeq2::DESeqDataSetFromTximport(isoform_obj$txi, colData = isoform_obj$meta,
                                            design = as.formula(paste0("~", condition)))
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(isoform_obj$counts),
                                          colData = isoform_obj$meta,
                                          design = as.formula(paste0("~", condition)))
  }
  dds[[condition]] <- relevel(dds[[condition]], base)
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds, test = "Wald")
  res <- DESeq2::results(dds, contrast = c(condition, level, base))
  
  res_df <- as.data.frame(res)
  res_df$transcript_id <- rownames(res_df)
  # Direct merge using already version‑free tx2gene (no stripping needed)
  res_df <- merge(res_df, isoform_obj$tx2gene, by.x = "transcript_id", by.y = "tx_id", all.x = TRUE)
  res_df$gene_symbol <- isoform_obj$gene_map$symbol[match(res_df$gene_id, isoform_obj$gene_map$ensembl)]
  res_df$signif <- !is.na(res_df$padj) & res_df$padj < padj_cutoff &
                   !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) > 1
  return(res_df)
}

#' Run DTU (Differential Transcript Usage) using DRIMSeq with improved filtering
#'
#' @param isoform_obj Output from import_transcript_counts()
#' @param condition Column name in metadata
#' @param level Foreground level for contrast
#' @param base Reference level
#' @param min_gene_expr Minimum total expression per gene (default: 10)
#' @param min_transcript_expr Minimum transcript proportion per gene (default: 0.05)
#' @param min_samps_gene_expr Min samples for gene expression (default: 50% of samples)
#' @param min_samps_feature_expr Min samples for transcript expression (default: 3)
#' @param chunk_size Number of genes to process at once (default: 5000)
#' @param max_transcripts Maximum number of transcripts per gene to keep (default: 300)
#' @param min_transcript_total Minimum total counts across all samples for a transcript (default: 10)
#' @return List with dtu_results
#' @export
run_dtu <- function(
  isoform_obj,
  condition,
  level,
  base,
  min_gene_expr = 10,
  min_transcript_expr = 0.05,
  min_samps_gene_expr = NULL,
  min_samps_feature_expr = 3,
  chunk_size = 5000,
  max_transcripts = 300,
  min_transcript_total = 10
) {

  if (!requireNamespace("DRIMSeq", quietly = TRUE)) {
    stop("DRIMSeq is required for DTU analysis.")
  }

  BiocParallel::register(BiocParallel::SerialParam())
  bp_param <- BiocParallel::SerialParam()

  counts <- if (isoform_obj$type == "tximport") {
    isoform_obj$txi$counts
  } else {
    isoform_obj$counts
  }

  sample_data <- isoform_obj$meta
  sample_data$condition <- factor(sample_data[[condition]], levels = c(base, level))

  num_samples <- ncol(counts)

  if (is.null(min_samps_gene_expr)) {
    min_samps_gene_expr <- ceiling(num_samples * 0.5)
  }

  min_samps_gene_expr <- min(min_samps_gene_expr, num_samples)
  min_samps_feature_expr <- min(min_samps_feature_expr, num_samples)

  tx2gene <- isoform_obj$tx2gene

  # Direct mapping: both row names and tx2gene$tx_id are already version‑free
  gene_id_map <- tx2gene$gene_id[match(rownames(counts), tx2gene$tx_id)]

  # Filter out transcripts with no gene mapping
  keep_mapped <- !is.na(gene_id_map)
  counts <- counts[keep_mapped, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_mapped]
  tx_ids_original <- rownames(counts)   # original (version‑free) IDs

  if (nrow(counts) == 0) {
    stop("No transcripts could be mapped to genes after version handling.")
  }

  message("Mapped transcripts: ", nrow(counts))

  keep_tx <- rowSums(counts) >= min_transcript_total
  counts <- counts[keep_tx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_tx]
  tx_ids_original <- tx_ids_original[keep_tx]

  if (nrow(counts) == 0) {
    stop("No transcripts passed min_transcript_total filter.")
  }

  keep_tx <- rowSums(counts > min_transcript_expr) >= min_samps_feature_expr
  counts <- counts[keep_tx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_tx]
  tx_ids_original <- tx_ids_original[keep_tx]

  if (nrow(counts) == 0) {
    stop("No transcripts passed expression filtering.")
  }

  gene_split <- split(seq_len(nrow(counts)), gene_id_map)

  keep_idx <- unlist(lapply(gene_split, function(idx) {
    if (length(idx) <= max_transcripts) return(idx)
    gene_expr <- rowMeans(counts[idx, , drop = FALSE])
    idx[order(gene_expr, decreasing = TRUE)[seq_len(max_transcripts)]]
  }))

  counts <- counts[keep_idx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_idx]
  tx_ids_original <- tx_ids_original[keep_idx]

  message("Final transcripts after capping: ", nrow(counts))
  message("Genes retained: ", length(unique(gene_id_map)))

  unique_genes <- unique(gene_id_map)
  gene_chunks <- split(unique_genes,
                       ceiling(seq_along(unique_genes) / chunk_size))

  all_results <- list()

  for (i in seq_along(gene_chunks)) {
    message("Processing chunk ", i, " / ", length(gene_chunks))

    current_genes <- gene_chunks[[i]]
    idx <- which(gene_id_map %in% current_genes)

    curr_counts <- counts[idx, , drop = FALSE]
    curr_gene_ids <- gene_id_map[idx]
    curr_tx_ids <- tx_ids_original[idx]

    curr_df <- data.frame(
      gene_id = curr_gene_ids,
      feature_id = curr_tx_ids,
      curr_counts,
      check.names = FALSE
    )

    d <- DRIMSeq::dmDSdata(counts = curr_df, samples = sample_data)

    d <- DRIMSeq::dmFilter(
      d,
      min_samps_gene_expr = min_samps_gene_expr,
      min_samps_feature_expr = min_samps_feature_expr,
      min_gene_expr = min_gene_expr,
      min_feature_expr = min_transcript_expr
    )

    if (nrow(DRIMSeq::counts(d)) == 0) {
      gc()
      next
    }

    design <- model.matrix(~ condition, data = DRIMSeq::samples(d))

    res <- tryCatch({
      d <- DRIMSeq::dmPrecision(d, design = design, BPPARAM = bp_param)
      d <- DRIMSeq::dmFit(d, design = design, BPPARAM = bp_param)
      d <- DRIMSeq::dmTest(d, coef = paste0("condition", level), BPPARAM = bp_param)
      DRIMSeq::results(d)
    }, error = function(e) {
      message("Chunk ", i, " failed: ", e$message)
      NULL
    })

    if (!is.null(res)) {
      all_results[[length(all_results) + 1]] <- res
    }

    rm(d, curr_counts, curr_df)
    gc()
  }

  dtu_results <- do.call(rbind, all_results)

  if (nrow(dtu_results) == 0) {
    stop("No DTU results were produced.")
  }

  pcol <- intersect(c("pvalue", "p_value"), colnames(dtu_results))[1]
  if (is.na(pcol)) {
    stop("No p-value column found in results.")
  }
  dtu_results$adj_pvalue <- p.adjust(dtu_results[[pcol]], method = "BH")

  return(list(dtu_results = dtu_results))
}

#' Run IsoformSwitchAnalyzeR analysis using isoformSwitchAnalysisCombined
#'
#' @param dte_results (unused, kept for compatibility)
#' @param dtu_results (unused, kept for compatibility)
#' @param isoform_obj Imported isoform data (from import_transcript_counts)
#' @param condition Column name in metadata containing groups
#' @param level Foreground level for contrast
#' @param base Reference level
#' @param fasta_file Path to transcript FASTA file
#' @param gff_file Path to GFF/GTF annotation
#' @param out_dir Output directory
#' @param use_wsl Logical: use WSL for external tools (Windows only)
#' @param wsl_distro WSL distribution name (default "Ubuntu-22.04")
#' @param save_dir Optional directory to save RDS files
#' @param resume_from Optional path to saved switch_list.rds
#' @param bsgenome_name Optional BSgenome package name (e.g. "BSgenome.Hsapiens.UCSC.hg38")
#' @return IsoformSwitchAnalyzeR results object
#' @export
run_isoform_switch <- function(dte_results = NULL, dtu_results = NULL,
                               isoform_obj, condition, level, base,
                               fasta_file, gff_file, out_dir,
                               use_wsl = FALSE, wsl_distro = "Ubuntu-22.04",
                               save_dir = NULL, resume_from = NULL,
                               bsgenome_name = NULL) {

  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE))
    stop("Please install IsoformSwitchAnalyzeR: BiocManager::install('IsoformSwitchAnalyzeR')")
  if (!requireNamespace("Biostrings", quietly = TRUE))
    install.packages("Biostrings")
  if (!requireNamespace("BSgenome", quietly = TRUE))
    BiocManager::install("BSgenome")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # --------------------------------------------------------------------------
  # 1. Load or build the SwitchList
  # --------------------------------------------------------------------------
  if (!is.null(resume_from) && file.exists(file.path(resume_from, "switch_list.rds"))) {
    message("Resuming from saved SwitchList: ", resume_from)
    switch_list <- readRDS(file.path(resume_from, "switch_list.rds"))
    if (!is.null(switch_list$isoformSwitchAnalysis)) {
      message("Full analysis already present – skipping.")
      return(switch_list)
    }
    run_analysis <- TRUE
  } else {
    message("Building SwitchList from raw data...")

    # Extract counts and metadata
    if (isoform_obj$type == "tximport") {
      count_matrix <- isoform_obj$txi$counts
    } else {
      count_matrix <- isoform_obj$counts
    }

    sample_col <- if ("Sample" %in% colnames(isoform_obj$meta)) "Sample" else "sample_id"
    sample_vector <- isoform_obj$meta[[sample_col]]
    if (is.null(sample_vector) || length(sample_vector) == 0) {
      sample_vector <- rownames(isoform_obj$meta)
    }
    if (!all(sample_vector %in% colnames(count_matrix))) {
      stop("Sample IDs in metadata do not match count matrix column names.")
    }

    # ----------------------------------------------------------------------
    # PRE-FILTERING: keep only transcripts that exist in the FASTA file
    # ----------------------------------------------------------------------
    message("Pre‑filtering transcripts to match FASTA file...")
    fasta_seqs <- Biostrings::fasta.seqlengths(fasta_file)
    fasta_ids <- names(fasta_seqs)

    # Clean transcript IDs (remove version suffixes, pipe, spaces)
    clean_id <- function(x) {
      x <- sub("\\..*$", "", x)  # remove version number
      x <- sub("\\|.*$", "", x)  # remove pipe
      x <- sub(" .*$", "", x)    # remove space
      x
    }
    clean_fasta_ids <- clean_id(fasta_ids)
    clean_rownames <- clean_id(rownames(count_matrix))

    # Keep only transcripts present in FASTA
    keep_in_fasta <- clean_rownames %in% clean_fasta_ids
    count_matrix <- count_matrix[keep_in_fasta, , drop = FALSE]
    message("  Kept ", nrow(count_matrix), " transcripts out of ", length(clean_rownames),
            " that have a match in the FASTA file.")

    if (nrow(count_matrix) == 0) {
      stop("No transcript IDs match between count matrix and FASTA file. ",
           "Check ID formats and consider using ignoreAfterPeriod = TRUE in importRdata().")
    }

    # Remove genes that end up with only one transcript (DRIMSeq requirement)
    tx2gene <- isoform_obj$tx2gene
    tx2gene$tx_clean <- clean_id(tx2gene$tx_id)
    keep_tx <- rownames(count_matrix) %in% tx2gene$tx_clean
    count_matrix <- count_matrix[keep_tx, , drop = FALSE]
    tx2gene <- tx2gene[match(rownames(count_matrix), tx2gene$tx_clean), ]
    gene_counts <- table(tx2gene$gene_id)
    genes_with_multiple <- names(gene_counts[gene_counts > 1])
    keep_genes <- tx2gene$gene_id %in% genes_with_multiple
    count_matrix <- count_matrix[keep_genes, , drop = FALSE]
    tx2gene <- tx2gene[keep_genes, ]
    message("  Kept ", nrow(count_matrix), " transcripts from ", length(unique(tx2gene$gene_id)),
            " genes after removing single‑transcript genes.")

    # Re‑build isoform_count_matrix with cleaned row names
    rownames(count_matrix) <- clean_id(rownames(count_matrix))
    isoform_count_matrix <- round(count_matrix)

    # Design matrix with correct factor levels
    design_matrix <- data.frame(
      sampleID = sample_vector,
      condition = factor(isoform_obj$meta[[condition]], levels = c(base, level)),
      stringsAsFactors = FALSE
    )
    rownames(design_matrix) <- design_matrix$sampleID

    # Ensure dplyr is loaded
    if (!"package:dplyr" %in% search()) attachNamespace("dplyr")

    # Import data
    switch_list <- IsoformSwitchAnalyzeR::importRdata(
      isoformCountMatrix   = isoform_count_matrix,
      isoformRepExpression = isoform_count_matrix,
      designMatrix         = design_matrix,
      isoformExonAnnoation = gff_file,
      isoformNtFasta       = fasta_file,
      ignoreAfterPeriod    = TRUE,
      showProgress         = TRUE
    )
    message("Isoform data import completed.")
    run_analysis <- TRUE
  }

  # --------------------------------------------------------------------------
  # 2. Run combined analysis (Part1 + Part2)
  # --------------------------------------------------------------------------
  if (run_analysis) {
    message("Running isoformSwitchAnalysisCombined...")

    # Determine if we need a BSgenome fallback
    if (!is.null(switch_list$ntSequence) && length(switch_list$ntSequence) > 0) {
      genome_object <- NULL
    } else if (!is.null(bsgenome_name) && requireNamespace(bsgenome_name, quietly = TRUE)) {
      genome_object <- getExportedValue(bsgenome_name, bsgenome_name)
      message("Using BSgenome: ", bsgenome_name)
    } else if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      # Use the default human BSgenome – note the quoted package name
      genome_object <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
      message("Using default BSgenome.Hsapiens.UCSC.hg38")
    } else {
      message("No BSgenome available – running analysis WITHOUT ORF prediction.")
      genome_object <- NULL
    }
    # Run the full workflow
    switch_list <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisCombined(
      switchAnalyzeRlist   = switch_list,
      genomeObject         = genome_object,
      pathToOutput         = out_dir,
      n                    = 50
    )
    message("Combined analysis completed.")
  }

  # --------------------------------------------------------------------------
  # 3. Save results
  # --------------------------------------------------------------------------
  if (!is.null(save_dir)) {
    saveRDS(switch_list, file.path(save_dir, "switch_list.rds"))
    message("Saved final SwitchList to ", save_dir)
  }

  message("Isoform switch analysis completed. Results saved in: ", out_dir)
  return(switch_list)
}

# Internal: run external predictors (Pfam, CPAT, SignalP) via system or WSL
.run_external_predictors <- function(switch_list, fasta_file, out_dir, use_wsl, wsl_distro, isoform_obj, save_dir) {
  is_windows <- .Platform$OS.type == "windows"
  
  if (is_windows && use_wsl) {
    base_prefix <- sprintf('wsl -d %s --', wsl_distro)
  } else {
    base_prefix <- ""
  }
  
  has_conda_env <- FALSE
  conda_sh <- NULL
  
  if (is_windows && use_wsl) {
    # 1. Find conda executable
    find_conda <- sprintf('wsl -d %s bash -c "command -v conda"', wsl_distro)
    conda_path <- system(find_conda, intern = TRUE, ignore.stderr = TRUE)
    if (length(conda_path) > 0 && nchar(conda_path[1]) > 0) {
      test_file <- sprintf('wsl -d %s test -f %s', wsl_distro, conda_path[1])
      if (system(test_file, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
        conda_base <- dirname(dirname(conda_path[1]))
        conda_sh <- file.path(conda_base, "etc/profile.d/conda.sh")
        test_sh <- sprintf('wsl -d %s test -f %s', wsl_distro, conda_sh)
        if (system(test_sh, ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) conda_sh <- NULL
      }
    }
    # 2. Fallback: search for conda.sh in /home and /opt
    if (is.null(conda_sh)) {
      find_sh <- sprintf(
        'wsl -d %s bash -c "find /home /opt -maxdepth 5 -type f -path \\"*/etc/profile.d/conda.sh\\" 2>/dev/null | head -1"',
        wsl_distro
      )
      conda_sh <- system(find_sh, intern = TRUE, ignore.stderr = TRUE)
      if (length(conda_sh) == 0 || nchar(conda_sh[1]) == 0) conda_sh <- NULL
    }
    # 3. Verify isoform_tools environment
    if (!is.null(conda_sh)) {
      check_env <- sprintf(
        'wsl -d %s bash -c "source %s && conda env list | grep -q isoform_tools"',
        wsl_distro, conda_sh
      )
      has_conda_env <- (system(check_env, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0)
    }
    
    if (has_conda_env) {
      message("Detected conda environment 'isoform_tools'")
      cmd_prefix <- sprintf(
        'wsl -d %s bash -c "source %s && conda run -n isoform_tools ',
        wsl_distro, conda_sh
      )
    } else {
      message("Conda environment 'isoform_tools' not found. Falling back to system PATH.")
      cmd_prefix <- base_prefix
    }
  } else {
    # Native Linux
    if (system("command -v conda", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0 &&
        system("conda env list | grep -q isoform_tools", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
      has_conda_env <- TRUE
      cmd_prefix <- "conda run -n isoform_tools "
      message("Detected conda environment 'isoform_tools' on native Linux.")
    } else {
      cmd_prefix <- ""
      message("No conda environment 'isoform_tools' found. Relying on system PATH.")
    }
  }
  
  # Hexamer file for CPAT
  first_symbol <- isoform_obj$gene_map$symbol[1]
  organism <- if (grepl("^[A-Z]+$", first_symbol) && nchar(first_symbol) > 2) "Human" else "Mouse"
  hexamer_file <- system.file("extdata", paste0(organism, "_Hexamer.tsv"), package = "IsoformSwitchAnalyzeR")
  if (!file.exists(hexamer_file)) {
    message("Hexamer file for ", organism, " not found. Falling back to Human_Hexamer.tsv")
    hexamer_file <- system.file("extdata", "Human_Hexamer.tsv", package = "IsoformSwitchAnalyzeR")
  }
  
  # ---- CPAT ----
  message("Running CPAT to assess coding potential...")
  cpat_out <- file.path(out_dir, "cpat_results.txt")
  check_cpat_cmd <- if (has_conda_env) {
    paste0(cmd_prefix, "run_cpat.py --help\"")
  } else {
    paste(cmd_prefix, "run_cpat.py --help")
  }
  check_cpat <- system(check_cpat_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (check_cpat == 0) {
    cpat_cmd <- if (has_conda_env) {
      sprintf('%s run_cpat.py -x %s -d %s -o %s"', cmd_prefix, fasta_file, hexamer_file, cpat_out)
    } else {
      sprintf('%s run_cpat.py -x %s -d %s -o %s', cmd_prefix, fasta_file, hexamer_file, cpat_out)
    }
    status <- system(cpat_cmd, wait = TRUE)
    if (status == 0) {
      switch_list <- IsoformSwitchAnalyzeR::addCPATanalysis(switch_list, cpat_out)
    } else {
      message("CPAT execution failed with status ", status)
    }
  } else {
    message("CPAT not found or not executable. Skipping coding potential prediction.")
  }
  
  # ---- SignalP ----
  message("Running SignalP (signal peptide prediction)...")
  signalp_out <- file.path(out_dir, "signalp_results")
  check_signalp_cmd <- if (has_conda_env) {
    paste0(cmd_prefix, "signalp -h\"")
  } else {
    paste(cmd_prefix, "signalp -h")
  }
  check_signalp <- system(check_signalp_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (check_signalp == 0) {
    signalp_cmd <- if (has_conda_env) {
      sprintf('%s signalp -fasta %s -output %s -format short"', cmd_prefix, fasta_file, signalp_out)
    } else {
      sprintf('%s signalp -fasta %s -output %s -format short', cmd_prefix, fasta_file, signalp_out)
    }
    status <- system(signalp_cmd, wait = TRUE)
    if (status == 0) {
      switch_list <- IsoformSwitchAnalyzeR::addSignalPanalysis(switch_list, signalp_out)
    } else {
      message("SignalP execution failed with status ", status)
    }
  } else {
    message("SignalP not found. Skipping signal peptide analysis.")
  }
  
  # ---- Pfam domain search ----
  message("Running Pfam domain search (via InterProScan or hmmscan)...")
  pfam_out <- file.path(out_dir, "pfam_results.xml")
  check_interpro_cmd <- if (has_conda_env) {
    paste0(cmd_prefix, "interproscan.sh -h\"")
  } else {
    paste(cmd_prefix, "interproscan.sh -h")
  }
  check_interpro <- system(check_interpro_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (check_interpro == 0) {
    interpro_cmd <- if (has_conda_env) {
      sprintf('%s interproscan.sh -i %s -f XML -o %s -dp -goterms -iprlookup"', cmd_prefix, fasta_file, pfam_out)
    } else {
      sprintf('%s interproscan.sh -i %s -f XML -o %s -dp -goterms -iprlookup', cmd_prefix, fasta_file, pfam_out)
    }
    status <- system(interpro_cmd, wait = TRUE)
    if (status == 0) {
      switch_list <- IsoformSwitchAnalyzeR::addPfamAnalysis(switch_list, pfam_out, "InterProScan")
    } else {
      message("InterProScan failed. Trying hmmscan...")
      pfam_db <- system.file("extdata", "Pfam-A.hmm", package = "IsoformSwitchAnalyzeR")
      if (file.exists(pfam_db)) {
        check_hmmscan_cmd <- if (has_conda_env) {
          paste0(cmd_prefix, "hmmscan -h\"")
        } else {
          paste(cmd_prefix, "hmmscan -h")
        }
        if (system(check_hmmscan_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
          hmmscan_cmd <- if (has_conda_env) {
            sprintf('%s hmmscan --tblout %s %s %s"', cmd_prefix, pfam_out, pfam_db, fasta_file)
          } else {
            sprintf('%s hmmscan --tblout %s %s %s', cmd_prefix, pfam_out, pfam_db, fasta_file)
          }
          status2 <- system(hmmscan_cmd, wait = TRUE)
          if (status2 == 0) {
            switch_list <- IsoformSwitchAnalyzeR::addPfamAnalysis(switch_list, pfam_out, "hmmscan")
          } else {
            message("hmmscan failed.")
          }
        } else {
          message("hmmscan not found.")
        }
      } else {
        message("Pfam database not found.")
      }
    }
  } else {
    message("InterProScan not found. Trying hmmscan...")
    pfam_db <- system.file("extdata", "Pfam-A.hmm", package = "IsoformSwitchAnalyzeR")
    if (file.exists(pfam_db)) {
      check_hmmscan_cmd <- if (has_conda_env) {
        paste0(cmd_prefix, "hmmscan -h\"")
      } else {
        paste(cmd_prefix, "hmmscan -h")
      }
      if (system(check_hmmscan_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
        hmmscan_cmd <- if (has_conda_env) {
          sprintf('%s hmmscan --tblout %s %s %s"', cmd_prefix, pfam_out, pfam_db, fasta_file)
        } else {
          sprintf('%s hmmscan --tblout %s %s %s', cmd_prefix, pfam_out, pfam_db, fasta_file)
        }
        status <- system(hmmscan_cmd, wait = TRUE)
        if (status == 0) {
          switch_list <- IsoformSwitchAnalyzeR::addPfamAnalysis(switch_list, pfam_out, "hmmscan")
        } else {
          message("hmmscan failed.")
        }
      } else {
        message("hmmscan not found.")
      }
    } else {
      message("Pfam database not found.")
    }
  }
  
  if (!is.null(save_dir)) {
    saveRDS(switch_list, file.path(save_dir, "switch_list.rds"))
    message("Saved final SwitchList to ", save_dir)
  }
  
  return(switch_list)
}

#' Generate DTE/DTU report with plots and interactive HTML table
#'
#' @param dte_results Data frame from run_dte()
#' @param dtu_results List from run_dtu() (contains dtu_results data frame)
#' @param isoform_obj Isoform import object (for gene mapping and counts)
#' @param out_dir Output directory (results root)
#' @param condition Condition column name
#' @param level Treatment level
#' @param base Control level
#' @param genes_of_interest Character vector of gene symbols to highlight/plot
#' @param top_n Number of top features to show in barplots
#' @return Invisible NULL, creates report and plots in IsoformSwitch/DTU_DTE_report/
#' @export
generate_dte_dtu_report <- function(dte_results, dtu_results, isoform_obj,
                                    out_dir, condition, level, base,
                                    genes_of_interest = NULL, top_n = 15) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("DT is required for interactive tables")
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown is required to generate the report")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required for data manipulation")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("tidyr is required for data reshaping")
  }
  
  report_dir <- file.path(out_dir, "IsoformSwitch", "DTU_DTE_report")
  if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
  
  # ---- 1. Annotate DTE results ----
  dte <- dte_results
  dte$gene_label <- paste0(dte$gene_symbol, " (", dte$transcript_id, ")")
  dte$signif_label <- ifelse(dte$signif, "Significant", "Not significant")
  
  # ---- 2. Annotate DTU results with gene symbols ----
  dtu <- dtu_results$dtu_results
  gene_map <- isoform_obj$gene_map[, c("ensembl", "symbol")]
  colnames(gene_map) <- c("gene_id", "gene_symbol")
  dtu <- merge(dtu, gene_map, by = "gene_id", all.x = TRUE)
  dtu$gene_label <- paste0(dtu$gene_symbol, " (", dtu$gene_id, ")")
  
  # Save annotated tables as CSV
  write.csv(dte, file.path(report_dir, "DTE_results_annotated.csv"), row.names = FALSE)
  write.csv(dtu, file.path(report_dir, "DTU_results_annotated.csv"), row.names = FALSE)
  
  # ---- 3. Generate plots ----
  plot_dir <- file.path(report_dir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir)
  
  # 3a. DTE Volcano plot
  dte_sig <- dte[!is.na(dte$padj), ]
  p_volcano <- ggplot2::ggplot(dte_sig, ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = signif_label)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "red")) +
    ggplot2::geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggplot2::labs(title = paste("DTE Volcano plot:", level, "vs", base),
         x = "log2 Fold Change", y = "-log10(adj. p-value)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  ggplot2::ggsave(file.path(plot_dir, "DTE_volcano.pdf"), p_volcano, width = 8, height = 6)
  
  # 3b. DTE MA plot (filter zero baseMean to avoid log10 infinite values)
  dte_ma <- dte[!is.na(dte$baseMean) & dte$baseMean > 0, ]
  p_ma <- ggplot2::ggplot(dte_ma, ggplot2::aes(x = baseMean, y = log2FoldChange, color = signif_label)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "red")) +
    ggplot2::geom_hline(yintercept = c(-1,1), linetype = "dashed") +
    ggplot2::labs(title = "DTE MA plot", x = "Mean of normalized counts", y = "log2 Fold Change") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  ggplot2::ggsave(file.path(plot_dir, "DTE_MA.pdf"), p_ma, width = 8, height = 6)
  
  # 3c. Top significant transcripts barplot (absolute log2FC) - FIXED: removed pipe
  dte_filtered <- dte[dte$signif, ]
  if (nrow(dte_filtered) > 0) {
    dte_ordered <- dte_filtered[order(dte_filtered$padj), ]
    top_dte <- head(dte_ordered, top_n)
    p_bar <- ggplot2::ggplot(top_dte, ggplot2::aes(x = reorder(gene_label, abs(log2FoldChange)), y = log2FoldChange)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste("Top", top_n, "significant DE transcripts"),
           x = "", y = "log2 Fold Change") +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(plot_dir, "DTE_top_barplot.pdf"), p_bar, width = 10, height = max(6, top_n*0.3))
  } else {
    message("No significant DTE transcripts found; skipping top barplot.")
  }
  
  # 3d. DTU p-value histogram
  p_hist <- ggplot2::ggplot(dtu, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
    ggplot2::labs(title = "DTU p-value distribution", x = "p-value", y = "Count") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(plot_dir, "DTU_pvalue_hist.pdf"), p_hist, width = 6, height = 4)
  
  # 3e. For each gene of interest, plot transcript proportions
  if (!is.null(genes_of_interest)) {
    counts <- if (isoform_obj$type == "tximport") isoform_obj$txi$counts else isoform_obj$counts
    meta <- isoform_obj$meta
    tx2gene <- isoform_obj$tx2gene
    
    for (gene_sym in genes_of_interest) {
      gene_ens <- isoform_obj$gene_map$ensembl[isoform_obj$gene_map$symbol == gene_sym]
      if (length(gene_ens) == 0) {
        message("Gene ", gene_sym, " not found in gene map. Skipping.")
        next
      }
      tx_ids <- tx2gene$tx_id[tx2gene$gene_id == gene_ens]
      if (length(tx_ids) == 0) next
      tx_counts <- counts[rownames(counts) %in% tx_ids, , drop = FALSE]
      # Convert to proportions per sample
      prop_df <- as.data.frame(t(tx_counts))  # samples as rows
      prop_df <- prop_df / rowSums(prop_df)
      prop_df$sample <- rownames(prop_df)
      prop_long <- tidyr::pivot_longer(prop_df, cols = -sample, names_to = "transcript", values_to = "proportion")
      prop_long <- merge(prop_long, meta, by.x = "sample", by.y = 0)
      
      p_prop <- ggplot2::ggplot(prop_long, ggplot2::aes(x = sample, y = proportion, fill = transcript)) +
        ggplot2::geom_bar(stat = "identity", position = "fill") +
        ggplot2::facet_grid(stats::as.formula(paste("~", condition)), scales = "free_x", space = "free_x") +
        ggplot2::labs(title = paste("Transcript proportions for", gene_sym),
             y = "Proportion", x = "Sample") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      ggplot2::ggsave(file.path(plot_dir, paste0("proportions_", gene_sym, ".pdf")), p_prop, width = 10, height = 6)
    }
  }
  
  # ---- 4. Create R Markdown report using absolute paths for CSV files ----
  rmd_file <- tempfile(fileext = ".Rmd")
  on.exit(unlink(rmd_file), add = TRUE)
  
  # Absolute paths for the CSV files (so the Rmd is self-contained).
  # winslash = "/" ensures forward slashes on Windows: backslashes in embedded
  # R string literals trigger the '\U used without hex digits' parse error
  # (e.g. C:\Users\... → \U is treated as a Unicode escape sequence).
  dte_abs_path <- normalizePath(file.path(report_dir, "DTE_results_annotated.csv"), winslash = "/")
  dtu_abs_path <- normalizePath(file.path(report_dir, "DTU_results_annotated.csv"), winslash = "/")
  
  rmd_lines <- c(
    "---",
    paste0("title: \"DTE/DTU Analysis Report (", level, " vs ", base, ")\""),
    "author: \"ExpressOM\"",
    paste0("date: \"`r Sys.Date()`\""),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_depth: 2",
    "    theme: flatly",
    "    highlight: tango",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(DT)",
    "library(ggplot2)",
    "```",
    "",
    "## Overview",
    "",
    paste0("- **Comparison**: ", level, " vs ", base),
    paste0("- **Date**: `r Sys.time()`"),
    paste0("- **Top N shown**: ", top_n),
    paste0("- **Genes of interest**: ", paste(genes_of_interest, collapse = ", ")),
    "",
    "## DTE Results (Differential Transcript Expression)",
    "",
    "### Volcano Plot",
    "![Volcano](plots/DTE_volcano.pdf)",
    "",
    "### MA Plot",
    "![MA](plots/DTE_MA.pdf)",
    "",
    "### Top Significant Transcripts",
    "![Top Barplot](plots/DTE_top_barplot.pdf)",
    "",
    "### Interactive Table (top 1000 by padj)",
    "```{r dt_dte}",
    paste0("dte_data <- read.csv(\"", dte_abs_path, "\", check.names = FALSE)"),
    "dte_data <- dte_data[order(dte_data$padj), ][1:min(1000, nrow(dte_data)), ]",
    "DT::datatable(dte_data, filter = \"top\", options = list(scrollX = TRUE))",
    "```",
    "",
    "## DTU Results (Differential Transcript Usage)",
    "",
    "### P-value Distribution",
    "![P-value histogram](plots/DTU_pvalue_hist.pdf)",
    "",
    "### Interactive Table",
    "```{r dt_dtu}",
    paste0("dtu_data <- read.csv(\"", dtu_abs_path, "\", check.names = FALSE)"),
    "DT::datatable(dtu_data, filter = \"top\", options = list(scrollX = TRUE))",
    "```",
    "",
    "## Gene‑specific Plots"
  )
  
  # Add gene-specific sections if requested
  if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
    for (g in genes_of_interest) {
      # Only add if the plot file exists
      if (file.exists(file.path(plot_dir, paste0("proportions_", g, ".pdf")))) {
        rmd_lines <- c(rmd_lines,
          paste0("\n### ", g),
          paste0("![Proportions](plots/proportions_", g, ".pdf)")
        )
      }
    }
  } else {
    rmd_lines <- c(rmd_lines, "\nNo custom genes requested.")
  }
  
  # Write to temporary Rmd file
  writeLines(rmd_lines, rmd_file)
  
  # Render the report
  rmarkdown::render(rmd_file, output_file = "report.html", output_dir = report_dir, quiet = TRUE)
  
  message("DTE/DTU report generated in: ", report_dir)
  invisible(NULL)
}