# mod_isoform_dtu.R -- differential transcript usage (DEXSeq engine)
#

#' Shared count-matrix preparation and filtering for DTU testing
#' @keywords internal
.prepare_dtu_counts <- function(isoform_obj,
                                condition,
                                level,
                                base,
                                min_transcript_total,
                                min_transcript_expr,
                                min_samps_feature_expr,
                                min_gene_expr,
                                min_samps_gene_expr,
                                max_transcripts,
                                apply_gene_expr_filter = FALSE,
                                require_multi_transcript = FALSE) {

  counts <- if (isoform_obj$type == "tximport") {
    isoform_obj$txi$counts
  } else {
    isoform_obj$counts
  }

  sample_data <- isoform_obj$meta

  keep_samples <- sample_data[[condition]] %in% c(base, level)

  if (sum(keep_samples) < 2) {
    stop("Fewer than 2 samples available for comparison. Need at least one sample in each condition.")
  }

  sample_data <- sample_data[keep_samples, , drop = FALSE]
  counts <- counts[, rownames(sample_data), drop = FALSE]

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
      "Only ", n_mapped, " / ", n_total, " transcripts (",
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
    min_samps_gene_expr = min_samps_gene_expr,
    min_samps_feature_expr = min_samps_feature_expr
  )
}

#' Run DTU using DRIMSeq
#' @export
run_dtu <- function(isoform_obj,
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
                    bpparam = NULL) {

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
  gene_chunks <- split(unique_genes, ceiling(seq_along(unique_genes) / chunk_size))

  all_results <- list()

  rownames(sample_data) <- trimws(rownames(sample_data))
  colnames(counts) <- trimws(colnames(counts))

  for (i in seq_along(gene_chunks)) {
    message("Processing chunk ", i, " / ", length(gene_chunks))

    current_genes <- gene_chunks[[i]]
    idx <- which(gene_id_map %in% current_genes)

    curr_counts <- counts[idx, , drop = FALSE]
    curr_gene_ids <- gene_id_map[idx]
    curr_tx_ids <- tx_ids_original[idx]

    common_samples <- intersect(colnames(curr_counts), rownames(sample_data))

    if (length(common_samples) == 0) {
      warning("Chunk ", i, ": No common samples between counts and metadata. Skipping.")
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
      gc()
      next
    }

    counts_d <- DRIMSeq::counts(d)
    samples_d <- DRIMSeq::samples(d)

    if (!identical(rownames(samples_d), colnames(counts_d))) {
      samples_d <- samples_d[colnames(counts_d), , drop = FALSE]
    }

    if (nrow(samples_d) != ncol(counts_d)) {
      warning(
        "Chunk ", i, ": Design matrix rows (", nrow(samples_d),
        ") do not match number of samples (", ncol(counts_d), "). Skipping."
      )
      next
    }

    design <- model.matrix(~ condition, data = samples_d)

    res <- tryCatch(
      {
        d <- DRIMSeq::dmPrecision(d, design = design, BPPARAM = bp_param)
        d <- DRIMSeq::dmFit(d, design = design, BPPARAM = bp_param)
        d <- DRIMSeq::dmTest(d, coef = paste0("condition", level), BPPARAM = bp_param)
        DRIMSeq::results(d, level = "feature")
      },
      error = function(e) {
        message("Chunk ", i, " failed: ", e$message)
        NULL
      }
    )

    if (!is.null(res)) {
      all_results[[length(all_results) + 1]] <- res
    }

    rm(d, curr_counts, curr_df)
    gc()
  }

  if (length(all_results) == 0) {
    stop(
      "No DTU results were produced. Check that your samples and metadata match, ",
      "and that there are genes with multiple transcripts meeting the filtering criteria."
    )
  }

  dtu_results <- do.call(rbind, all_results)
  rm(all_results)
  gc()

  pcol <- intersect(c("pvalue", "p_value"), colnames(dtu_results))[1]

  if (is.na(pcol)) {
    stop("No p-value column found in results.")
  }

  dtu_results$adj_pvalue <- p.adjust(dtu_results[[pcol]], method = "BH")

  dtu_results$gene_symbol <- isoform_obj$gene_map$symbol[
    match(dtu_results$gene_id, isoform_obj$gene_map$ensembl)
  ]

  dtu_results$entrezid <- isoform_obj$gene_map$entrezid[
    match(dtu_results$gene_id, isoform_obj$gene_map$ensembl)
  ]

  dtu_results$gene <- .coalesce_gene_label(dtu_results$gene_symbol, dtu_results$feature_id)

  return(list(dtu_results = dtu_results))
}

#' Run DEXSeq-based DTU
#' @export
run_dexseq_dtu <- function(isoform_obj,
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
                           bpparam = NULL) {

  if (!requireNamespace("DEXSeq", quietly = TRUE)) {
    stop("DEXSeq is required for run_dexseq_dtu(). Install with BiocManager::install('DEXSeq').")
  }

  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam
  BiocParallel::register(bp_param)

  prep <- .prepare_dtu_counts(
    isoform_obj,
    condition,
    level,
    base,
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
    "DEXSeq DTU: ", nrow(counts), " transcripts across ",
    length(unique(gene_id_map)), " genes."
  )

  unique_genes <- unique(gene_id_map)
  gene_chunks <- split(unique_genes, ceiling(seq_along(unique_genes) / chunk_size))

  sample_meta <- data.frame(
    condition = sample_data[colnames(counts), "condition"],
    row.names = colnames(counts)
  )

  all_results <- list()
  dxr_list <- list()

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
          design = ~sample + exon + condition:exon,
          featureID = curr_tx_ids,
          groupID = curr_gene_ids
        )
      },
      error = function(e) {
        message("  Chunk ", i, " DEXSeqDataSet construction failed: ", e$message)
        NULL
      }
    )

    if (is.null(dxd)) next

    res <- tryCatch(
      {
        dxd <- DEXSeq::estimateSizeFactors(dxd)
        dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE, BPPARAM = bp_param)
        dxd <- DEXSeq::testForDEU(dxd, BPPARAM = bp_param)
        dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = bp_param)
        DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
      },
      error = function(e) {
        message("  Chunk ", i, " DEXSeq testing failed: ", e$message)
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
    }

    rm(dxd, curr_counts)
    gc()
  }

  if (length(all_results) == 0) stop("No DEXSeq DTU results were produced.")

  results_df <- do.call(rbind, all_results)
  rm(all_results)
  gc()

  pval_by_gene <- split(results_df$pvalue, results_df$groupID)

  gene_pvals <- vapply(pval_by_gene, function(p) {
    p <- p[!is.na(p)]
    if (length(p) == 0) return(NA_real_)
    n <- length(p)
    min(sort(p) * n / seq_len(n))
  }, numeric(1))

  gene_qvals <- stats::p.adjust(gene_pvals, method = "BH")

  results_df$gene_pvalue_simes <- gene_pvals[results_df$groupID]
  results_df$gene_qvalue <- gene_qvals[results_df$groupID]

  colnames(results_df)[colnames(results_df) == "featureID"] <- "transcript_id"
  colnames(results_df)[colnames(results_df) == "groupID"] <- "gene_id"

  results_df$gene_symbol <- isoform_obj$gene_map$symbol[
    match(results_df$gene_id, isoform_obj$gene_map$ensembl)
  ]

  results_df$gene <- .coalesce_gene_label(results_df$gene_symbol, results_df$transcript_id)

  n_sig_genes <- sum(gene_qvals < 0.05, na.rm = TRUE)

  message(
    "DEXSeq DTU complete: ", length(unique(results_df$gene_id)), " genes tested, ",
    n_sig_genes, " significant at gene q-value < 0.05."
  )

  list(results_df = results_df, dxr_list = dxr_list)
}

#' DEXSeq-style transcript usage plot for a single gene
#' @export
plot_dexseq_gene <- function(dxr_list, gene_id, plot_dir, gene_symbol = NULL, splicing = FALSE) {
  if (!requireNamespace("DEXSeq", quietly = TRUE)) {
    message("DEXSeq is required for plot_dexseq_gene(). Skipping.")
    return(invisible(NULL))
  }

  if (is.null(dxr_list) || length(dxr_list) == 0) {
    message("No DEXSeqResults available (run run_dexseq_dtu() with keep_dxr = TRUE). Skipping.")
    return(invisible(NULL))
  }

  label <- if (!is.null(gene_symbol) && nzchar(gene_symbol)) gene_symbol else gene_id
  pdf_path <- file.path(plot_dir, paste0("DEXSeq_", label, ".pdf"))

  plotted <- FALSE

  for (dxr in dxr_list) {
    if (is.null(dxr)) next
    if (!(gene_id %in% dxr$groupID)) next

    ok <- tryCatch(
      {
        .pdf_device()(pdf_path, width = 9, height = 6)

        DEXSeq::plotDEXSeq(
          dxr,
          geneID = gene_id,
          fitExpToVar = "condition",
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
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
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
    message("   -> Gene '", label, "' not found in DEXSeq results (or not testable). Skipping plot.")
    return(invisible(NULL))
  }

  message("   -> DEXSeq transcript-usage plot saved to: ", pdf_path)

  invisible(pdf_path)
}

