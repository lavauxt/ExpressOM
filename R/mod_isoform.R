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
#' @param bpparam Optional BiocParallel parameter object for DRIMSeq operations
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
  min_transcript_total = 10,
  bpparam = NULL
) {

  if (!requireNamespace("DRIMSeq", quietly = TRUE)) {
    stop("DRIMSeq is required for DTU analysis.")
  }

  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam
  BiocParallel::register(bp_param)

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
#' @param run_predictors Logical: run external predictors (CPAT, SignalP, Pfam)
#' @param use_wsl Logical: use WSL for external tools (Windows only)
#' @param wsl_distro WSL distribution name (default "Ubuntu-22.04")
#' @param save_dir Optional directory to save per-step RDS checkpoints
#' @param resume_from Optional path to a save_dir with a prior switch_list.rds
#' @param bsgenome_name Optional BSgenome package name
#' @return IsoformSwitchAnalyzeR results object
#' @export
run_isoform_switch <- function(dte_results = NULL, dtu_results = NULL,
                               isoform_obj, condition, level, base,
                               fasta_file, gff_file, out_dir,
                               run_predictors = FALSE,
                               use_wsl = FALSE, wsl_distro = "Ubuntu-22.04",
                               save_dir = NULL, resume_from = NULL,
                               bsgenome_name = NULL) {

  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE))
    stop("Please install IsoformSwitchAnalyzeR: BiocManager::install('IsoformSwitchAnalyzeR')")
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Please install Biostrings: BiocManager::install('Biostrings')")

  if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
  if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # --------------------------------------------------------------------------
  # WSL debug and logging (if predictors requested)
  # --------------------------------------------------------------------------
  if (run_predictors && .Platform$OS.type == "windows" && use_wsl) {
    message("Running WSL environment check...")
    debug_info <- debug_wsl(distro = wsl_distro, out_dir = out_dir,
                            conda_env = "isoform_tools", verbose = TRUE)
    if (!debug_info$wsl_available) {
      stop("WSL is not available. Please install WSL or set use_wsl=FALSE.")
    }
    if (run_predictors && !debug_info$conda_env_exists) {
      warning("Conda environment 'isoform_tools' not found. Predictors may fail.\n",
              "Run install_wsl_isoform_tools() to set up the environment.")
    }
    # Write summary to log
    .log_wsl_command("debug_wsl()", exit_code = 0,
                     stdout = capture.output(print(debug_info)),
                     log_dir = log_dir)
  }

  # --------------------------------------------------------------------------
  # Checkpoint helpers (scoped to this call)
  # --------------------------------------------------------------------------
  .ckpt_path   <- function(nm) if (!is.null(save_dir)) file.path(save_dir, nm) else ""
  .ckpt_exists <- function(nm) { p <- .ckpt_path(nm); nzchar(p) && file.exists(p) }
  .ckpt_save   <- function(obj, nm) {
    p <- .ckpt_path(nm)
    if (nzchar(p)) {
      saveRDS(obj, p)
      message("Checkpoint saved -> ", nm)
    }
    invisible(obj)
  }
  .ckpt_load   <- function(nm) {
    p <- .ckpt_path(nm)
    message("Resuming from checkpoint: ", nm)
    readRDS(p)
  }

  # --------------------------------------------------------------------------
  # Step 1 – Load or build the SwitchList
  # --------------------------------------------------------------------------
  switch_list      <- NULL
  already_analyzed <- FALSE

  # Priority: explicit resume_from > step2 checkpoint > step1 checkpoint > build fresh

  if (!is.null(resume_from) && file.exists(file.path(resume_from, "switch_list.rds"))) {
    message("Resuming from saved SwitchList: ", resume_from)
    switch_list <- readRDS(file.path(resume_from, "switch_list.rds"))
    if (!is.null(switch_list$isoformSwitchAnalysis)) {
      message("Full analysis already present - skipping combined analysis.")
      already_analyzed <- TRUE
    }

  } else if (.ckpt_exists("step2_analyzed.rds")) {
    switch_list      <- .ckpt_load("step2_analyzed.rds")
    already_analyzed <- TRUE

  } else if (.ckpt_exists("step1_imported.rds")) {
    switch_list <- .ckpt_load("step1_imported.rds")

  } else {
    # ---- Build SwitchList from raw data ----
    message("Building SwitchList from raw data...")

    if (isoform_obj$type == "tximport") {
      count_matrix <- isoform_obj$txi$counts
    } else {
      count_matrix <- isoform_obj$counts
    }

    sample_col    <- if ("Sample" %in% colnames(isoform_obj$meta)) "Sample" else "sample_id"
    sample_vector <- isoform_obj$meta[[sample_col]]
    if (is.null(sample_vector) || length(sample_vector) == 0)
      sample_vector <- rownames(isoform_obj$meta)
    if (!all(sample_vector %in% colnames(count_matrix)))
      stop("Sample IDs in metadata do not match count matrix column names.")

    # Pre-filter: keep only transcripts present in the FASTA file
    message("Pre-filtering transcripts to match FASTA file...")
    fasta_seqs    <- Biostrings::fasta.seqlengths(fasta_file)
    clean_id      <- function(x) {
      x <- sub("\\..*$", "", x)
      x <- sub("\\|.*$", "", x)
      x <- sub(" .*$",   "", x)
      x
    }
    clean_fasta_ids  <- clean_id(names(fasta_seqs))
    clean_rownames   <- clean_id(rownames(count_matrix))
    keep_in_fasta    <- clean_rownames %in% clean_fasta_ids
    count_matrix     <- count_matrix[keep_in_fasta, , drop = FALSE]
    message("  Kept ", nrow(count_matrix), " / ", length(clean_rownames),
            " transcripts matching the FASTA file.")
    if (nrow(count_matrix) == 0)
      stop("No transcript IDs match between count matrix and FASTA file. ",
           "Check ID formats or set ignoreAfterPeriod = TRUE in importRdata().")

    # Remove genes with only one transcript (DRIMSeq requirement)
    tx2gene      <- isoform_obj$tx2gene
    tx2gene$tx_clean <- clean_id(tx2gene$tx_id)
    keep_tx      <- rownames(count_matrix) %in% tx2gene$tx_clean
    count_matrix <- count_matrix[keep_tx, , drop = FALSE]
    tx2gene      <- tx2gene[match(rownames(count_matrix), tx2gene$tx_clean), ]
    gene_counts  <- table(tx2gene$gene_id)
    genes_multi  <- names(gene_counts[gene_counts > 1])
    keep_genes   <- tx2gene$gene_id %in% genes_multi
    count_matrix <- count_matrix[keep_genes, , drop = FALSE]
    tx2gene      <- tx2gene[keep_genes, ]
    message("  Kept ", nrow(count_matrix), " transcripts from ",
            length(unique(tx2gene$gene_id)), " genes after removing single-transcript genes.")

    rownames(count_matrix)  <- clean_id(rownames(count_matrix))
    isoform_count_matrix    <- round(count_matrix)

    design_matrix <- data.frame(
      sampleID  = sample_vector,
      condition = factor(isoform_obj$meta[[condition]], levels = c(base, level)),
      stringsAsFactors = FALSE
    )
    rownames(design_matrix) <- design_matrix$sampleID

    if (!"package:dplyr" %in% search()) attachNamespace("dplyr")

    # ---- importRdata: ALL PARAMETERS UNCHANGED ----
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

    # Save step-1 checkpoint immediately
    .ckpt_save(switch_list, "step1_imported.rds")
  }

  # --------------------------------------------------------------------------
  # Step 2 – Run combined analysis (Part1 + Part2)
  # --------------------------------------------------------------------------
  if (!already_analyzed) {
    message("Running isoformSwitchAnalysisCombined...")

    if (!is.null(switch_list$ntSequence) && length(switch_list$ntSequence) > 0) {
      genome_object <- NULL
    } else if (!is.null(bsgenome_name) && requireNamespace(bsgenome_name, quietly = TRUE)) {
      genome_object <- getExportedValue(bsgenome_name, bsgenome_name)
      message("Using BSgenome: ", bsgenome_name)
    } else if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      genome_object <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
      message("Using default BSgenome.Hsapiens.UCSC.hg38")
    } else {
      message("No BSgenome available - running WITHOUT ORF prediction.")
      genome_object <- NULL
    }

    # ---- isoformSwitchAnalysisCombined: ALL PARAMETERS UNCHANGED ----
    switch_list <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisCombined(
      switchAnalyzeRlist = switch_list,
      genomeObject       = genome_object,
      pathToOutput       = out_dir,
      n                  = 50
    )
    message("Combined analysis completed.")

    # Save step-2 checkpoint immediately
    .ckpt_save(switch_list, "step2_analyzed.rds")
  }

  # --------------------------------------------------------------------------
  # Step 3 – External predictors (CPAT, SignalP, Pfam)
  # --------------------------------------------------------------------------
  if (isTRUE(run_predictors)) {
    if (.ckpt_exists("step3_predictors.rds")) {
      message("Loading external predictor results from step-3 checkpoint.")
      switch_list <- .ckpt_load("step3_predictors.rds")
    } else {
      message("Running external isoform predictors...")
      switch_list <- .run_external_predictors(
        switch_list = switch_list,
        fasta_file  = fasta_file,
        out_dir     = out_dir,
        use_wsl     = use_wsl,
        wsl_distro  = wsl_distro,
        isoform_obj = isoform_obj,
        save_dir    = save_dir
      )
      .ckpt_save(switch_list, "step3_predictors.rds")
    }
  }

  # --------------------------------------------------------------------------
  # Final save
  # --------------------------------------------------------------------------
  if (!is.null(save_dir)) {
    saveRDS(switch_list, file.path(save_dir, "switch_list.rds"))
    message("Saved final SwitchList to ", save_dir)
  }

  message("Isoform switch analysis completed. Results saved in: ", out_dir)
  return(switch_list)
}

# ==============================================================================
# Internal: run external predictors (CPAT, SignalP, Pfam) via WSL or natively
# ==============================================================================
#
# Key design principles:
#   * All system calls use .wsl_exec_script() which writes a temp bash script
#     and executes it - completely avoids shell-quoting hell on Windows.
#   * Every Windows path is converted to a WSL Unix path via .to_wsl_path()
#     before being embedded in bash commands.
#   * CPAT: correct flag order -g (fasta) -x (hexamer) -d (logit model) -o (prefix)
#           supports CPAT3 (`cpat`) and CPAT2 (`run_cpat.py`)
#   * SignalP: supports v5 (`signalp`) and v6 (`signalp6`)
#   * Pfam: InterProScan XML first; falls back to hmmscan --domtblout
#   * ISA functions: analyzeCPAT / analyzeSignalIP / analyzePFAM / analyzeInterProScan
#
.run_external_predictors <- function(switch_list, fasta_file, out_dir,
                                      use_wsl, wsl_distro, isoform_obj, save_dir) {
  is_windows <- .Platform$OS.type == "windows"
  via_wsl    <- is_windows && use_wsl

  log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # --------------------------------------------------------------------------
  # 1. Discover conda.sh and verify 'isoform_tools' environment
  # --------------------------------------------------------------------------
  message("Detecting conda environment in execution context...")
  conda_sh      <- NULL
  has_conda_env <- FALSE

  find_out <- .wsl_exec_script(
    bash_body = paste(
      "find /home /opt /root -maxdepth 8",
      "-name 'conda.sh' -path '*/profile.d/*'",
      "2>/dev/null | head -5"
    ),
    wsl_distro = wsl_distro, use_wsl = via_wsl,
    conda_sh = NULL, intern = TRUE, ignore_stderr = TRUE,
    log_dir = log_dir
  )
  find_out <- trimws(find_out[nzchar(trimws(find_out))])

  for (csh in find_out) {
    ok <- .wsl_exec_script(
      bash_body = c(
        sprintf('. "%s"', csh),
        'conda env list 2>/dev/null | grep -q "isoform_tools"'
      ),
      wsl_distro = wsl_distro, use_wsl = via_wsl,
      conda_sh = NULL, intern = FALSE, ignore_stderr = TRUE,
      log_dir = log_dir
    )
    if (isTRUE(ok == 0L)) { conda_sh <- csh; has_conda_env <- TRUE; break }
  }

  active_conda <- if (has_conda_env) conda_sh else NULL

  if (has_conda_env)
    message("  Using conda env 'isoform_tools' (", conda_sh, ")")
  else
    message("  No 'isoform_tools' env found. Using system PATH.")

  # --------------------------------------------------------------------------
  # 2. Convenience wrappers
  # --------------------------------------------------------------------------

  # Run a single bash command string (conda activation is handled automatically)
  run_tool <- function(cmd_str, show_stderr = TRUE) {
    .wsl_exec_script(
      bash_body     = cmd_str,
      wsl_distro    = wsl_distro,
      use_wsl       = via_wsl,
      conda_sh      = active_conda,
      conda_env     = "isoform_tools",
      intern        = FALSE,
      ignore_stderr = !show_stderr,
      log_dir       = log_dir
    )
  }

  # Check tool availability (conda-aware)
  tool_ok <- function(t) {
    .wsl_tool_exists(t, wsl_distro, via_wsl, active_conda, "isoform_tools")
  }

  # Convert a Windows path to a WSL path (identity on Linux/Mac)
  w2l <- function(p) {
    if (is.null(p) || !nzchar(p)) return(p)
    if (is_windows && use_wsl) .to_wsl_path(p, wsl_distro) else p
  }

  # --------------------------------------------------------------------------
  # 3. Determine organism (for CPAT hexamer / logit model selection)
  # --------------------------------------------------------------------------
  syms     <- isoform_obj$gene_map$symbol[!is.na(isoform_obj$gene_map$symbol)]
  first_sym <- if (length(syms) > 0) syms[1] else ""
  organism <- if (grepl("^[A-Z]+$", first_sym) && nchar(first_sym) > 2) "Human" else "Mouse"
  message("Organism detected for CPAT: ", organism)

  # --------------------------------------------------------------------------
  # 4. Extract NT and AA sequences from the switch_list
  #    CPAT needs NT; SignalP and Pfam need AA (protein) sequences.
  # --------------------------------------------------------------------------
  seq_dir <- file.path(out_dir, "sequences")
  dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

  nt_fa_local <- file.path(seq_dir, "isoform_NT.fa")
  aa_fa_local <- file.path(seq_dir, "isoform_AA.fa")

  if (!file.exists(nt_fa_local) || !file.exists(aa_fa_local)) {
    message("Extracting NT / AA sequences from SwitchList for external tools...")
    tryCatch({
      IsoformSwitchAnalyzeR::extractSequence(
        switch_list,
        onlySwitchingGenes       = FALSE,
        extractNTseq             = TRUE,
        extractAAseq             = TRUE,
        writeToFile              = TRUE,
        pathToOutput             = seq_dir,
        addToSwitchAnalyzeRlist  = FALSE,
        quiet                    = TRUE
      )
      # ISA writes NT.fasta / AA.fasta; rename to stable names
      for (pair in list(c("NT.fasta", nt_fa_local), c("AA.fasta", aa_fa_local))) {
        src <- file.path(seq_dir, pair[1])
        if (file.exists(src) && normalizePath(src) != normalizePath(pair[2]))
          file.rename(src, pair[2])
      }
    }, error = function(e) {
      message("  Could not extract sequences from SwitchList: ", e$message)
      message("  CPAT will use the provided fasta_file; SignalP / Pfam may be skipped.")
    })
  }

  # Input selection per tool
  cpat_fa_local <- if (file.exists(nt_fa_local)) nt_fa_local else fasta_file
  sp_fa_local   <- if (file.exists(aa_fa_local)) aa_fa_local else NULL
  pfam_fa_local <- if (file.exists(aa_fa_local)) aa_fa_local else NULL

  # Convert to WSL paths
  cpat_fa_w  <- w2l(cpat_fa_local)
  sp_fa_w    <- w2l(sp_fa_local)
  pfam_fa_w  <- w2l(pfam_fa_local)
  out_dir_w  <- w2l(out_dir)

  # --------------------------------------------------------------------------
  # 5. CPAT – coding potential prediction
  # --------------------------------------------------------------------------
  message("\n--- CPAT (coding potential) ---")

  hexamer_local <- system.file("extdata", paste0(organism, "_Hexamer.tsv"),
                                package = "IsoformSwitchAnalyzeR")
  logit_local   <- .find_cpat_logit_model(organism, wsl_distro, via_wsl, active_conda)

  if (!nzchar(hexamer_local) || !file.exists(hexamer_local)) {
    message("  Hexamer table not found in ISA package. Skipping CPAT.")
  } else if (is.null(logit_local)) {
    message("  Logit model not found. Run install_isoform_databases() first. Skipping CPAT.")
  } else {
    hexamer_w      <- w2l(hexamer_local)
    logit_w        <- logit_local  # already a valid path in execution env
    cpat_prefix_w  <- paste0(out_dir_w, "/cpat_out")
    cpat_prefix_l  <- file.path(out_dir, "cpat_out")

    has_cpat3 <- tool_ok("cpat")
    has_cpat2 <- tool_ok("run_cpat.py")

    if (!has_cpat3 && !has_cpat2) {
      message("  Neither 'cpat' nor 'run_cpat.py' found. Skipping.")
    } else {
      cpat_exe    <- if (has_cpat3) "cpat" else "run_cpat.py"
      # Correct CPAT flag order: -g <NT fasta>  -x <hexamer>  -d <logit model>  -o <prefix>
      cpat_cmd <- sprintf(
        "%s -g %s -x %s -d %s -o %s",
        cpat_exe,
        shQuote(cpat_fa_w,     type = "sh"),
        shQuote(hexamer_w,     type = "sh"),
        shQuote(logit_w,       type = "sh"),
        shQuote(cpat_prefix_w, type = "sh")
      )
      message("  Running: ", cpat_exe)
      cpat_status <- run_tool(cpat_cmd)

      # CPAT3 main result: <prefix>.ORF_prob.best.tsv; CPAT2: <prefix>.r
      cpat_result <- if (has_cpat3) paste0(cpat_prefix_l, ".ORF_prob.best.tsv")
                     else           paste0(cpat_prefix_l, ".r")

      if (cpat_status == 0L && file.exists(cpat_result)) {
        switch_list <- tryCatch(
          IsoformSwitchAnalyzeR::analyzeCPAT(
            switchAnalyzeRlist       = switch_list,
            pathToAllCPATresultFiles = cpat_result,
            codingCutoff             = NULL,
            removeNoncodinORFs       = TRUE,
            quiet                    = FALSE
          ),
          error = function(e) { message("  analyzeCPAT(): ", e$message); switch_list }
        )
        message("  CPAT results imported successfully.")
      } else {
        message("  CPAT failed or result file not found (", basename(cpat_result), "). Skipping.")
      }
    }
  }

  # --------------------------------------------------------------------------
  # 6. SignalP – signal peptide prediction (requires AA sequences)
  # --------------------------------------------------------------------------
  message("\n--- SignalP (signal peptide prediction) ---")

  if (is.null(sp_fa_w)) {
    message("  No AA FASTA available (ORF analysis may not have run). Skipping SignalP.")
  } else {
    sp_out_dir_l <- file.path(out_dir, "signalp_out")
    dir.create(sp_out_dir_l, recursive = TRUE, showWarnings = FALSE)
    sp_out_dir_w <- w2l(sp_out_dir_l)

    has_sp6 <- tool_ok("signalp6")
    has_sp5 <- tool_ok("signalp")

    sp_result_file <- NULL
    sp_status      <- 1L

    if (has_sp6) {
      # SignalP 6 syntax
      sp_cmd <- sprintf(
        "signalp6 --fastafile %s --organism eukarya --output_dir %s --format txt --mode fast",
        shQuote(sp_fa_w, type = "sh"), shQuote(sp_out_dir_w, type = "sh")
      )
      message("  Running SignalP 6...")
      sp_status <- run_tool(sp_cmd)
      if (sp_status == 0L) {
        candidates <- c(
          file.path(sp_out_dir_l, "prediction_results.txt"),
          file.path(sp_out_dir_l, "output.gff3")
        )
        found <- Filter(file.exists, candidates)
        if (length(found) > 0) sp_result_file <- found[1]
      }

    } else if (has_sp5) {
      # SignalP 5 syntax
      sp_cmd <- sprintf(
        "signalp -fasta %s -org euk -format short -output %s",
        shQuote(sp_fa_w, type = "sh"), shQuote(sp_out_dir_w, type = "sh")
      )
      message("  Running SignalP 5...")
      sp_status <- run_tool(sp_cmd)
      if (sp_status == 0L) {
        sp_files <- list.files(sp_out_dir_l,
                               pattern    = "_summary\\.signalp5$",
                               full.names = TRUE)
        if (length(sp_files) > 0) sp_result_file <- sp_files[1]
      }

    } else {
      message("  Neither signalp6 nor signalp found in PATH / conda env. Skipping.")
    }

    if (!is.null(sp_result_file)) {
      switch_list <- tryCatch(
        IsoformSwitchAnalyzeR::analyzeSignalIP(
          switchAnalyzeRlist      = switch_list,
          pathToSignalPresultFile = sp_result_file,
          quiet                   = FALSE
        ),
        error = function(e) { message("  analyzeSignalIP(): ", e$message); switch_list }
      )
      message("  SignalP results imported successfully.")
    } else if (sp_status != 0L) {
      message("  SignalP execution failed (exit ", sp_status, "). Skipping.")
    }
  }

  # --------------------------------------------------------------------------
  # 7. Pfam domain annotation (InterProScan -> hmmscan fallback)
  #    Both require AA / protein sequences.
  # --------------------------------------------------------------------------
  message("\n--- Pfam domain annotation ---")

  if (is.null(pfam_fa_w)) {
    message("  No AA FASTA available. Skipping Pfam.")
  } else {
    pfam_done <- FALSE

    # 7a. InterProScan (preferred – produces richer XML output)
    if (tool_ok("interproscan.sh")) {
      iprscan_xml_l <- file.path(out_dir, "interproscan.xml")
      iprscan_xml_w <- w2l(iprscan_xml_l)

      ips_cmd <- sprintf(
        "interproscan.sh -i %s -f XML -o %s -dp -appl Pfam -goterms -iprlookup",
        shQuote(pfam_fa_w,     type = "sh"),
        shQuote(iprscan_xml_w, type = "sh")
      )
      message("  Running InterProScan...")
      ips_status <- run_tool(ips_cmd)

      if (ips_status == 0L && file.exists(iprscan_xml_l)) {
        switch_list <- tryCatch(
          IsoformSwitchAnalyzeR::analyzeInterProScan(
            switchAnalyzeRlist           = switch_list,
            pathToInterProScanResultFile = iprscan_xml_l,
            quiet                        = FALSE
          ),
          error = function(e) { message("  analyzeInterProScan(): ", e$message); switch_list }
        )
        pfam_done <- TRUE
        message("  InterProScan Pfam results imported successfully.")
      } else {
        message("  InterProScan failed (exit ", ips_status, "). Trying hmmscan fallback...")
      }
    }

    # 7b. hmmscan fallback (requires Pfam-A.hmm indexed with hmmpress)
    if (!pfam_done) {
      if (!tool_ok("hmmscan")) {
        message("  hmmscan not found. Skipping Pfam analysis.")
      } else {
        pfam_db_w <- .find_pfam_db(wsl_distro, via_wsl, active_conda)

        if (is.null(pfam_db_w)) {
          message("  Pfam-A.hmm not found. Run install_isoform_databases() first. Skipping.")
        } else {
          pfam_tbl_l <- file.path(out_dir, "pfam_domtblout.txt")
          pfam_tbl_w <- w2l(pfam_tbl_l)

          # Use --domtblout (domain table) – the format ISA's analyzePFAM expects
          hm_cmd <- sprintf(
            "hmmscan --cpu 4 --domtblout %s %s %s",
            shQuote(pfam_tbl_w, type = "sh"),
            shQuote(pfam_db_w,  type = "sh"),
            shQuote(pfam_fa_w,  type = "sh")
          )
          message("  Running hmmscan...")
          hm_status <- run_tool(hm_cmd)

          if (hm_status == 0L && file.exists(pfam_tbl_l)) {
            switch_list <- tryCatch(
              IsoformSwitchAnalyzeR::analyzePFAM(
                switchAnalyzeRlist   = switch_list,
                pathToPFAMresultFile = pfam_tbl_l,
                quiet                = FALSE
              ),
              error = function(e) { message("  analyzePFAM(): ", e$message); switch_list }
            )
            message("  hmmscan Pfam results imported successfully.")
          } else {
            message("  hmmscan failed (exit ", hm_status, "). Skipping Pfam.")
          }
        }
      }
    }
  }

  return(switch_list)
}

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