# mod_isoform_dtu.R -- differential transcript usage (DEXSeq engine)
#

#' Shared count-matrix preparation and filtering for DTU testing
#' @keywords internal
.prepare_dtu_counts <- function(
  isoform_obj,
  condition,
  level,
  base,
  design = NULL,
  min_transcript_total,
  min_transcript_expr,
  min_samps_feature_expr,
  min_gene_expr,
  min_samps_gene_expr,
  max_transcripts,
  apply_gene_expr_filter = FALSE,
  require_multi_transcript = FALSE
) {
  counts <- if (isoform_obj$type == "tximport") {
    isoform_obj$txi$counts
  } else {
    isoform_obj$counts
  }

  design <- design %||% stats::reformulate(condition)
  design_info <- .validate_isoform_design(
    design,
    condition,
    isoform_obj$meta
  )
  sample_data <- design_info$metadata

  keep_samples <- sample_data[[condition]] %in% c(base, level)

  if (sum(keep_samples) < 2) {
    stop(
      "Fewer than 2 samples available for comparison. Need at least one sample in each condition."
    )
  }

  sample_data <- sample_data[keep_samples, , drop = FALSE]
  counts <- counts[, rownames(sample_data), drop = FALSE]

  sample_data <- .set_contrast_reference(
    sample_data,
    condition,
    base,
    level
  )

  sample_data$condition <- factor(
    sample_data[[condition]],
    levels = c(base, level)
  )

  num_samples <- ncol(counts)

  if (is.null(min_samps_gene_expr)) {
    min_samps_gene_expr <- ceiling(num_samples * 0.5)
  }

  min_samps_gene_expr <- min(min_samps_gene_expr, num_samples)
  min_samps_feature_expr <- min(min_samps_feature_expr, num_samples)

  tx2gene <- isoform_obj$tx2gene

  rownames(counts) <- clean_transcript_id(rownames(counts))
  tx2gene$tx_id <- clean_transcript_id(tx2gene$tx_id)

  tx_ids_original <- rownames(counts)

  gene_id_map <- tx2gene$gene_id[match(rownames(counts), tx2gene$tx_id)]
  keep_mapped <- !is.na(gene_id_map)

  n_total <- length(rownames(counts))
  n_mapped <- sum(keep_mapped)

  if (n_mapped == 0) {
    stop(
      "No transcripts could be mapped to genes after version handling.\n\n",
      "This usually means the transcript IDs in the count matrix and the transcript IDs in tx2gene are not using the same format or annotation source.\n\n",
      "Example count transcript IDs:\n  ",
      paste(utils::head(rownames(counts), 5), collapse = "\n  "),
      "\n\nExample tx2gene tx_id values:\n  ",
      paste(utils::head(tx2gene$tx_id, 5), collapse = "\n  "),
      "\n\nPossible fixes:\n",
      "  1. Update clean_transcript_id() to strip pipe/space/version suffixes.\n",
      "  2. Use a custom_tx2gene built from the exact same GTF/FASTA used for quantification.\n",
      "  3. Make sure the EnsDb package release matches the Ensembl/GENCODE release used for kallisto/Salmon indexing.\n",
      call. = FALSE
    )
  }

  if (n_mapped < 0.5 * n_total) {
    warning(
      "Only ",
      n_mapped,
      " / ",
      n_total,
      " transcripts (",
      round(100 * n_mapped / n_total, 1),
      "%) could be mapped to genes. This may indicate an annotation mismatch.",
      call. = FALSE
    )
  }

  counts <- counts[keep_mapped, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_mapped]
  tx_ids_original <- tx_ids_original[keep_mapped]

  if (nrow(counts) == 0) {
    stop("No transcripts remained after transcript-to-gene mapping.")
  }

  message("Mapped transcripts: ", nrow(counts))

  keep_tx <- rowSums(counts) >= min_transcript_total

  counts <- counts[keep_tx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_tx]
  tx_ids_original <- tx_ids_original[keep_tx]

  if (nrow(counts) == 0) {
    stop("No transcripts passed min_transcript_total filter.")
  }

  if (isTRUE(apply_gene_expr_filter)) {
    gene_expr_matrix <- rowsum(counts, group = gene_id_map)

    keep_genes_expr <- rownames(gene_expr_matrix)[
      rowSums(gene_expr_matrix >= min_gene_expr) >= min_samps_gene_expr
    ]

    keep_tx <- gene_id_map %in% keep_genes_expr

    counts <- counts[keep_tx, , drop = FALSE]
    gene_id_map <- gene_id_map[keep_tx]
    tx_ids_original <- tx_ids_original[keep_tx]

    if (nrow(counts) == 0) {
      stop("No genes passed the min_gene_expr / min_samps_gene_expr filter.")
    }
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
    if (length(idx) <= max_transcripts) {
      return(idx)
    }

    gene_expr <- rowMeans(counts[idx, , drop = FALSE])
    idx[order(gene_expr, decreasing = TRUE)[seq_len(max_transcripts)]]
  }))

  counts <- counts[keep_idx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_idx]
  tx_ids_original <- tx_ids_original[keep_idx]

  if (nrow(counts) == 0) {
    stop("No transcripts remained after limiting transcripts per gene.")
  }

  if (isTRUE(require_multi_transcript)) {
    gene_counts_tbl <- table(gene_id_map)
    multi_tx_genes <- names(gene_counts_tbl)[gene_counts_tbl > 1]
    keep_multi <- gene_id_map %in% multi_tx_genes

    counts <- counts[keep_multi, , drop = FALSE]
    gene_id_map <- gene_id_map[keep_multi]
    tx_ids_original <- tx_ids_original[keep_multi]

    if (nrow(counts) == 0) {
      stop("No multi-transcript genes remained after filtering.")
    }
  }

  message("Final transcripts after filtering: ", nrow(counts))
  message("Genes retained: ", length(unique(gene_id_map)))

  list(
    counts = counts,
    gene_id_map = gene_id_map,
    tx_ids_original = tx_ids_original,
    sample_data = sample_data,
    design = design_info$formula,
    min_samps_gene_expr = min_samps_gene_expr,
    min_samps_feature_expr = min_samps_feature_expr
  )
}

.adjust_dexseq_feature_pvalues <- function(results_df) {
  results_df$padj <- stats::p.adjust(results_df$pvalue, method = "BH")
  results_df
}

.integerize_dtu_counts <- function(counts) {
  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"
  counts
}

.normalize_dtu_pvalues <- function(results_df) {
  pcol <- intersect(c("pvalue", "p_value"), colnames(results_df))[1]
  if (is.na(pcol)) {
    stop("No p-value column found in results.", call. = FALSE)
  }
  results_df$pvalue <- suppressWarnings(as.numeric(as.character(results_df[[pcol]])))
  results_df
}

#' Run DTU using DRIMSeq
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
  bpparam = NULL,
  design = NULL
) {
  if (!requireNamespace("DRIMSeq", quietly = TRUE)) {
    stop("DRIMSeq is required for DTU analysis.")
  }

  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam
  BiocParallel::register(bp_param)

  prep <- .prepare_dtu_counts(
    isoform_obj,
    condition,
    level,
    base,
    design = design,
    min_transcript_total = min_transcript_total,
    min_transcript_expr = min_transcript_expr,
    min_samps_feature_expr = min_samps_feature_expr,
    min_gene_expr = min_gene_expr,
    min_samps_gene_expr = min_samps_gene_expr,
    max_transcripts = max_transcripts,
    apply_gene_expr_filter = FALSE,
    require_multi_transcript = FALSE
  )

  counts <- prep$counts
  gene_id_map <- prep$gene_id_map
  tx_ids_original <- prep$tx_ids_original
  sample_data <- prep$sample_data
  min_samps_gene_expr <- prep$min_samps_gene_expr
  min_samps_feature_expr <- prep$min_samps_feature_expr

  unique_genes <- unique(gene_id_map)
  gene_chunks <- split(
    unique_genes,
    ceiling(seq_along(unique_genes) / chunk_size)
  )

  all_results <- list()
  failed_chunks <- integer()

  rownames(sample_data) <- trimws(rownames(sample_data))
  colnames(counts) <- trimws(colnames(counts))

  for (i in seq_along(gene_chunks)) {
    message("Processing chunk ", i, " / ", length(gene_chunks))

    current_genes <- gene_chunks[[i]]
    idx <- which(gene_id_map %in% current_genes)

    curr_counts <- .integerize_dtu_counts(counts[idx, , drop = FALSE])
    curr_gene_ids <- gene_id_map[idx]
    curr_tx_ids <- tx_ids_original[idx]

    common_samples <- intersect(colnames(curr_counts), rownames(sample_data))

    if (length(common_samples) == 0) {
      warning(
        "Chunk ",
        i,
        ": No common samples between counts and metadata. Skipping."
      )
      failed_chunks <- c(failed_chunks, i)
      next
    }

    curr_counts <- curr_counts[, common_samples, drop = FALSE]
    sample_data_sub <- sample_data[common_samples, , drop = FALSE]

    curr_df <- data.frame(
      gene_id = curr_gene_ids,
      feature_id = curr_tx_ids,
      curr_counts,
      check.names = FALSE
    )

    d <- DRIMSeq::dmDSdata(counts = curr_df, samples = sample_data_sub)

    d <- DRIMSeq::dmFilter(
      d,
      min_samps_gene_expr = min_samps_gene_expr,
      min_samps_feature_expr = min_samps_feature_expr,
      min_gene_expr = min_gene_expr,
      min_feature_expr = min_transcript_expr
    )

    if (nrow(DRIMSeq::counts(d)) == 0) {
      # Filtering legitimately removes genes; this is not an engine failure.
      gc()
      next
    }

    counts_d <- DRIMSeq::counts(d)
    samples_d <- DRIMSeq::samples(d)

    if (!identical(rownames(samples_d), colnames(counts_d))) {
      sample_idx <- match(colnames(counts_d), rownames(samples_d))
      if (anyNA(sample_idx)) {
        warning(
          "Chunk ",
          i,
          ": count samples are missing from metadata. Skipping."
        )
        failed_chunks <- c(failed_chunks, i)
        next
      }
      samples_d <- samples_d[sample_idx, , drop = FALSE]
    }

    if (nrow(samples_d) != ncol(counts_d)) {
      warning(
        "Chunk ",
        i,
        ": Design matrix rows (",
        nrow(samples_d),
        ") do not match number of samples (",
        ncol(counts_d),
        "). Skipping."
      )
      failed_chunks <- c(failed_chunks, i)
      next
    }

    design_info <- .validate_isoform_design(
      design %||% stats::reformulate(condition),
      condition,
      samples_d
    )
    model_matrix <- stats::model.matrix(
      design_info$formula,
      data = design_info$metadata
    )
    coef_name <- .design_contrast_coef(
      design_info$formula,
      design_info$metadata,
      condition,
      level,
      base
    )

    res <- tryCatch(
      {
        d <- DRIMSeq::dmPrecision(d, design = model_matrix, BPPARAM = bp_param)
        d <- DRIMSeq::dmFit(d, design = model_matrix, BPPARAM = bp_param)
        d <- DRIMSeq::dmTest(d, coef = coef_name, BPPARAM = bp_param)
        DRIMSeq::results(d, level = "feature")
      },
      error = function(e) {
        message("Chunk ", i, " failed: ", e$message)
        failed_chunks <<- c(failed_chunks, i)
        NULL
      }
    )

    if (!is.null(res)) {
      all_results[[length(all_results) + 1]] <- res
    }

    rm(d, curr_counts, curr_df)
    gc()
  }

  .stop_on_failed_chunks(failed_chunks, "DTU")

  if (length(all_results) == 0) {
    stop(
      "No DTU results were produced. Check that your samples and metadata match, ",
      "and that there are genes with multiple transcripts meeting the filtering criteria."
    )
  }

  dtu_results <- do.call(rbind, all_results)
  rm(all_results)
  gc()

  dtu_results$test_unit <- "transcript-feature-level; gene-level FDR is not provided"

  dtu_results <- .normalize_dtu_pvalues(dtu_results)
  dtu_results$adj_pvalue <- stats::p.adjust(dtu_results$pvalue, method = "BH")

  dtu_results$gene_symbol <- isoform_obj$gene_map$symbol[
    match(dtu_results$gene_id, isoform_obj$gene_map$ensembl)
  ]

  dtu_results$entrezid <- isoform_obj$gene_map$entrezid[
    match(dtu_results$gene_id, isoform_obj$gene_map$ensembl)
  ]

  dtu_results$gene <- .coalesce_gene_label(
    dtu_results$gene_symbol,
    dtu_results$feature_id
  )

  return(list(dtu_results = dtu_results))
}

#' Run DEXSeq-based DTU
#' @export
run_dexseq_dtu <- function(
  isoform_obj,
  condition,
  level,
  base,
  min_gene_expr = 10,
  min_transcript_expr = 0.05,
  min_samps_gene_expr = NULL,
  min_samps_feature_expr = 3,
  chunk_size = 2000,
  max_transcripts = 300,
  min_transcript_total = 10,
  keep_dxr = TRUE,
  bpparam = NULL,
  design = NULL
) {
  if (!requireNamespace("DEXSeq", quietly = TRUE)) {
    stop(
      "DEXSeq is required for run_dexseq_dtu(). Install with BiocManager::install('DEXSeq')."
    )
  }

  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam
  BiocParallel::register(bp_param)

  prep <- .prepare_dtu_counts(
    isoform_obj,
    condition,
    level,
    base,
    design = design,
    min_transcript_total = min_transcript_total,
    min_transcript_expr = min_transcript_expr,
    min_samps_feature_expr = min_samps_feature_expr,
    min_gene_expr = min_gene_expr,
    min_samps_gene_expr = min_samps_gene_expr,
    max_transcripts = max_transcripts,
    apply_gene_expr_filter = TRUE,
    require_multi_transcript = TRUE
  )

  counts <- prep$counts
  gene_id_map <- prep$gene_id_map
  tx_ids_original <- prep$tx_ids_original
  sample_data <- prep$sample_data

  message(
    "DEXSeq DTU: ",
    nrow(counts),
    " transcripts across ",
    length(unique(gene_id_map)),
    " genes."
  )

  unique_genes <- unique(gene_id_map)
  gene_chunks <- split(
    unique_genes,
    ceiling(seq_along(unique_genes) / chunk_size)
  )

  sample_meta <- sample_data[colnames(counts), , drop = FALSE]
  dex_designs <- .dexseq_usage_designs(
    prep$design,
    condition
  )

  all_results <- list()
  dxr_list <- list()
  failed_chunks <- integer()

  for (i in seq_along(gene_chunks)) {
    message("DEXSeq DTU: processing chunk ", i, " / ", length(gene_chunks))

    current_genes <- gene_chunks[[i]]
    idx <- which(gene_id_map %in% current_genes)

    curr_counts <- round(as.matrix(counts[idx, , drop = FALSE]))
    curr_gene_ids <- gene_id_map[idx]
    curr_tx_ids <- tx_ids_original[idx]

    dxd <- tryCatch(
      {
        DEXSeq::DEXSeqDataSet(
          countData = curr_counts,
          sampleData = sample_meta,
          design = dex_designs$full,
          featureID = curr_tx_ids,
          groupID = curr_gene_ids
        )
      },
      error = function(e) {
        message(
          "  Chunk ",
          i,
          " DEXSeqDataSet construction failed: ",
          e$message
        )
        NULL
      }
    )

    if (is.null(dxd)) {
      next
    }

    res <- tryCatch(
      {
        dxd <- DEXSeq::estimateSizeFactors(dxd)
        dxd <- DEXSeq::estimateDispersions(
          dxd,
          quiet = TRUE,
          BPPARAM = bp_param
        )
        dxd <- DEXSeq::testForDEU(
          dxd,
          fullModel = dex_designs$full,
          reducedModel = dex_designs$reduced,
          BPPARAM = bp_param
        )
        dxd <- DEXSeq::estimateExonFoldChanges(
          dxd,
          fitExpToVar = condition,
          BPPARAM = bp_param
        )
        DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
      },
      error = function(e) {
        message("  Chunk ", i, " DEXSeq testing failed: ", e$message)
        failed_chunks <<- c(failed_chunks, i)
        NULL
      }
    )

    if (!is.null(res)) {
      res_df <- as.data.frame(res)

      wanted <- intersect(
        c(
          "groupID",
          "featureID",
          "exonBaseMean",
          "dispersion",
          "stat",
          "pvalue",
          "padj",
          grep("^log2fold_", colnames(res_df), value = TRUE)
        ),
        colnames(res_df)
      )

      res_df <- res_df[, wanted, drop = FALSE]
      all_results[[length(all_results) + 1]] <- res_df

      if (isTRUE(keep_dxr)) dxr_list[[length(dxr_list) + 1]] <- res
    } else {
      failed_chunks <- c(failed_chunks, i)
    }

    rm(dxd, curr_counts)
    gc()
  }

  .stop_on_failed_chunks(failed_chunks, "DEXSeq DTU")

  if (length(all_results) == 0) {
    stop("No DEXSeq DTU results were produced.")
  }

  results_df <- do.call(rbind, all_results)
  rm(all_results)
  gc()

  # Chunks are an implementation detail: adjust feature-level p-values over
  # the full set of tested transcripts, not independently within each chunk.
  results_df <- .adjust_dexseq_feature_pvalues(results_df)
  results_df$test_unit <- "transcript-feature rows; DEXSeq feature-level results"

  pval_by_gene <- split(results_df$pvalue, results_df$groupID)

  gene_pvals <- vapply(
    pval_by_gene,
    function(p) {
      p <- p[!is.na(p)]
      if (length(p) == 0) {
        return(NA_real_)
      }
      n <- length(p)
      min(sort(p) * n / seq_len(n))
    },
    numeric(1)
  )

  gene_qvals <- stats::p.adjust(gene_pvals, method = "BH")

  results_df$gene_pvalue_simes <- gene_pvals[results_df$groupID]
  results_df$gene_qvalue <- gene_qvals[results_df$groupID]

  colnames(results_df)[colnames(results_df) == "featureID"] <- "transcript_id"
  colnames(results_df)[colnames(results_df) == "groupID"] <- "gene_id"

  results_df$gene_symbol <- isoform_obj$gene_map$symbol[
    match(results_df$gene_id, isoform_obj$gene_map$ensembl)
  ]

  results_df$gene <- .coalesce_gene_label(
    results_df$gene_symbol,
    results_df$transcript_id
  )

  n_sig_genes <- sum(gene_qvals < 0.05, na.rm = TRUE)

  message(
    "DEXSeq DTU complete: ",
    length(unique(results_df$gene_id)),
    " genes tested, ",
    n_sig_genes,
    " significant at gene q-value < 0.05."
  )

  list(results_df = results_df, dxr_list = dxr_list)
}

#' DEXSeq-style transcript usage plot for a single gene
#' @export
plot_dexseq_gene <- function(
  dxr_list,
  gene_id,
  plot_dir,
  gene_symbol = NULL,
  splicing = FALSE,
  condition = "condition"
) {
  if (!requireNamespace("DEXSeq", quietly = TRUE)) {
    message("DEXSeq is required for plot_dexseq_gene(). Skipping.")
    return(invisible(NULL))
  }

  if (is.null(dxr_list) || length(dxr_list) == 0) {
    message(
      "No DEXSeqResults available (run run_dexseq_dtu() with keep_dxr = TRUE). Skipping."
    )
    return(invisible(NULL))
  }

  label <- if (!is.null(gene_symbol) && nzchar(gene_symbol)) {
    gene_symbol
  } else {
    gene_id
  }
  pdf_path <- file.path(plot_dir, paste0("DEXSeq_", label, ".pdf"))

  plotted <- FALSE

  for (dxr in dxr_list) {
    if (is.null(dxr)) {
      next
    }
    if (!(gene_id %in% dxr$groupID)) {
      next
    }

    ok <- tryCatch(
      {
        .pdf_device()(pdf_path, width = 9, height = 6)

        DEXSeq::plotDEXSeq(
          dxr,
          geneID = gene_id,
          fitExpToVar = condition,
          expression = !splicing,
          splicing = splicing,
          legend = TRUE,
          cex.axis = 1.1,
          cex = 1.2,
          lwd = 2
        )

        dev.off()
        TRUE
      },
      error = function(e) {
        if (grDevices::dev.cur() > 1) {
          grDevices::dev.off()
        }
        message("   -> plotDEXSeq() failed for ", label, ": ", e$message)
        FALSE
      }
    )

    if (ok) {
      plotted <- TRUE
      break
    }
  }

  if (!plotted) {
    message(
      "   -> Gene '",
      label,
      "' not found in DEXSeq results (or not testable). Skipping plot."
    )
    return(invisible(NULL))
  }

  message("   -> DEXSeq transcript-usage plot saved to: ", pdf_path)

  invisible(pdf_path)
}

# mod_isoform_dte.R -- differential transcript expression (DESeq2) and transcript-level PCA
#

#' Run DTE using DESeq2
#' @param design Optional one-sided additive design formula. Defaults to
#'   `~ condition`; include the `condition` column and any adjustment covariates.
#' @export
run_dte <- function(isoform_obj, condition, level, base, padj_cutoff = 0.05,
                    bpparam = NULL, design = NULL) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for DTE analysis. Please install it.")
  }

  design <- design %||% stats::reformulate(condition)
  design_info <- .validate_isoform_design(
    design, condition, isoform_obj$meta
  )
  meta <- .set_contrast_reference(
    design_info$metadata, condition, base, level
  )

  if (isoform_obj$type == "tximport") {
    dds <- DESeq2::DESeqDataSetFromTximport(
      isoform_obj$txi,
      colData = meta,
      design = design_info$formula
    )
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(isoform_obj$counts),
      colData = meta,
      design = design_info$formula
    )
  }

  dds[[condition]] <- .safe_relevel_condition(dds[[condition]], base, label = condition)

  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]

  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam

  dds <- DESeq2::DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = bp_param)
  res <- DESeq2::results(
    dds,
    contrast = c(condition, level, base),
    alpha = padj_cutoff
  )
  res_df <- as.data.frame(res)

  res_df$transcript_id <- rownames(res_df)
  res_df$transcript_id <- clean_transcript_id(res_df$transcript_id)

  tx2gene_clean <- isoform_obj$tx2gene
  tx2gene_clean$tx_id <- clean_transcript_id(tx2gene_clean$tx_id)

  res_df <- merge(
    res_df,
    tx2gene_clean,
    by.x = "transcript_id",
    by.y = "tx_id",
    all.x = TRUE
  )

  res_df$gene_symbol <- isoform_obj$gene_map$symbol[
    match(res_df$gene_id, isoform_obj$gene_map$ensembl)
  ]

  res_df$entrezid <- isoform_obj$gene_map$entrezid[
    match(res_df$gene_id, isoform_obj$gene_map$ensembl)
  ]

  res_df$gene <- .coalesce_gene_label(res_df$gene_symbol, res_df$transcript_id)

  res_df$signif <- !is.na(res_df$padj) &
    res_df$padj < padj_cutoff &
    !is.na(res_df$log2FoldChange) &
    abs(res_df$log2FoldChange) > 1

  return(res_df)
}

#' Transcript-level PCA
#' @param pca_colors,sample_labels See `plot_custom_pca()`
#' @export
run_isoform_pca <- function(isoform_obj,
                            condition,
                            level = NULL,
                            base = NULL,
                            out_dir,
                            batch_col = NULL,
                            pca_ntop = 500,
                            save_dir = NULL,
                            pca_colors = NULL,
                            sample_labels = NULL) {

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for transcript-level PCA. Please install it.")
  }

  if (isoform_obj$type == "tximport") {
    dds <- DESeq2::DESeqDataSetFromTximport(
      isoform_obj$txi,
      colData = isoform_obj$meta,
      design = ~1
    )
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(isoform_obj$counts),
      colData = isoform_obj$meta,
      design = ~1
    )
  }

  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]

  if (nrow(dds) < 2) {
    message("   -> Fewer than 2 transcripts pass the expression filter; skipping transcript-level PCA.")
    return(invisible(NULL))
  }

  # vst() needs >= 1000 features; fall back like the gene-level PCA does
  vsd <- tryCatch(
    DESeq2::vst(dds, blind = TRUE),
    error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  )

  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  # PCA_transcripts_*.tsv/.rds data tables (as opposed to the PDF plot
  # itself, which stays in Plots/) go to save_dir when one is configured --
  # this was previously the one isoform-side intermediary left behind in
  # Plots/ instead of alongside the RDS checkpoints and cleaned_reference.
  data_dir <- if (!is.null(save_dir)) save_dir else plot_dir
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  comp_label <- if (!is.null(level) && !is.null(base)) paste0(level, "vs", base) else "AllSamples"
  cond_arg <- if (!is.null(condition) && condition %in% colnames(isoform_obj$meta)) condition else NULL

  res <- plot_custom_pca(
    vsd,
    condition = cond_arg,
    batch = batch_col,
    title = paste0(
      "Transcript-level PCA",
      if (!is.null(cond_arg)) paste0(" (", cond_arg, ")") else ""
    ),
    return_plot = TRUE,
    return_gene_list = TRUE,
    ntop = pca_ntop,
    level = level,
    base = base,
    pca_colors = pca_colors,
    sample_labels = sample_labels
  )

  .pdf_device()(file.path(plot_dir, paste0("PCA_transcripts_", comp_label, ".pdf")), width = 9, height = 7)
  print(res$plot)
  dev.off()

  tx2gene_clean <- isoform_obj$tx2gene
  tx2gene_clean$tx_id <- clean_transcript_id(tx2gene_clean$tx_id)

  if (!is.null(res$gene_info) && nrow(res$gene_info) > 0) {
    tx_info <- res$gene_info
    tx_info$gene <- clean_transcript_id(tx_info$gene)
    colnames(tx_info)[colnames(tx_info) == "gene"] <- "transcript_id"

    tx_info$gene_id <- tx2gene_clean$gene_id[
      match(tx_info$transcript_id, tx2gene_clean$tx_id)
    ]

    tx_info$gene_symbol <- isoform_obj$gene_map$symbol[
      match(tx_info$gene_id, isoform_obj$gene_map$ensembl)
    ]

    tx_info$gene_transcript <- ifelse(
      is.na(tx_info$gene_symbol) | tx_info$gene_symbol == "",
      tx_info$transcript_id,
      paste0(tx_info$gene_symbol, " (", tx_info$transcript_id, ")")
    )

    col_order <- c("transcript_id", "gene_id", "gene_symbol", "gene_transcript", "variance")
    tx_info <- tx_info[, intersect(col_order, colnames(tx_info))]

    write.table(
      tx_info,
      file = file.path(data_dir, paste0("PCA_transcripts_", comp_label, ".tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    message(
      "   -> Saved top variable transcripts (with gene names) used in PCA to: ",
      file.path(data_dir, paste0("PCA_transcripts_", comp_label, ".tsv"))
    )
  }

  tx_values <- res$gene_values

  if (!is.null(tx_values) && nrow(tx_values) > 0) {
    tx_values$gene <- clean_transcript_id(tx_values$gene)
    colnames(tx_values)[colnames(tx_values) == "gene"] <- "transcript_id"

    tx_values$gene_id <- tx2gene_clean$gene_id[
      match(tx_values$transcript_id, tx2gene_clean$tx_id)
    ]

    tx_values$gene_symbol <- isoform_obj$gene_map$symbol[
      match(tx_values$gene_id, isoform_obj$gene_map$ensembl)
    ]

    front_cols <- intersect(
      c("transcript_id", "gene_id", "gene_symbol"),
      colnames(tx_values)
    )

    tx_values <- tx_values[, c(front_cols, setdiff(colnames(tx_values), front_cols))]
  }

  .write_pca_scores(
    res$pca_scores,
    data_dir,
    paste0("PCA_transcripts_", comp_label),
    gene_values = tx_values
  )

  message("Transcript-level PCA completed. Plot saved in: ", plot_dir, "; data tables in: ", data_dir)

  invisible(res)
}

