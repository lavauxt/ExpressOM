# mod_isoform_dte.R -- differential transcript expression (DESeq2) and transcript-level PCA
#

#' Run DTE using DESeq2
#' @export
run_dte <- function(isoform_obj, condition, level, base, padj_cutoff = 0.05, bpparam = NULL) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for DTE analysis. Please install it.")
  }

  if (isoform_obj$type == "tximport") {
    dds <- DESeq2::DESeqDataSetFromTximport(
      isoform_obj$txi,
      colData = isoform_obj$meta,
      design = as.formula(paste0("~ ", condition))
    )
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(isoform_obj$counts),
      colData = isoform_obj$meta,
      design = as.formula(paste0("~ ", condition))
    )
  }

  dds[[condition]] <- .safe_relevel_condition(dds[[condition]], base, label = condition)

  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]

  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam

  dds <- DESeq2::DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = bp_param)
  res <- DESeq2::results(dds, contrast = c(condition, level, base))
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
#' @export
run_isoform_pca <- function(isoform_obj,
                            condition,
                            level = NULL,
                            base = NULL,
                            out_dir,
                            batch_col = NULL,
                            pca_ntop = 500) {

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

  vsd <- DESeq2::vst(dds, blind = TRUE)

  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

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
    ntop = pca_ntop
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
      file = file.path(plot_dir, paste0("PCA_transcripts_", comp_label, ".tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    message(
      "   -> Saved top variable transcripts (with gene names) used in PCA to: ",
      file.path(plot_dir, paste0("PCA_transcripts_", comp_label, ".tsv"))
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
    plot_dir,
    paste0("PCA_transcripts_", comp_label),
    gene_values = tx_values
  )

  message("Transcript-level PCA completed. Plots saved in: ", plot_dir)

  invisible(res)
}

