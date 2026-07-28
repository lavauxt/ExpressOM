# ==============================================================================
# mod_isoform.R - DTE, DTU, and IsoformSwitchAnalyzeR integration
# ==============================================================================
# NOTE: .fill_entrez_with_bitr() previously had a duplicate, drifting copy
# here (marked "copied from mod_dge.R"). It now lives solely in utils_core.R
# and is used by both the DGE and isoform import paths.

#' Convert a PDF plot to a PNG alongside it
#'
#' Used by \code{generate_dte_dtu_report()} so plots can be embedded in the
#' HTML report via \code{knitr::include_graphics()}: browsers cannot render a
#' PDF through an \code{<img>} tag, so every figure needs a PNG rendition to
#' actually display in \code{report.html} (PDFs remain the primary output
#' format for the standalone plot files themselves).
#'
#' @param pdf_file Path to the source PDF.
#' @param dpi Resolution for the rendered PNG.
#' @return Path to the created PNG, or NULL if conversion wasn't possible.
#' @keywords internal
convert_pdf_to_png <- function(pdf_file, dpi = 200) {

  if (!file.exists(pdf_file)) {
    warning("PDF not found: ", pdf_file)
    return(NULL)
  }
  if (!requireNamespace("pdftools", quietly = TRUE)) {
    warning("Package 'pdftools' not installed; cannot convert ", pdf_file,
            " to PNG for HTML embedding (it will fall back to the unrendered PDF path).")
    return(NULL)
  }

  png_file <- sub("\\.pdf$", ".png", pdf_file)

  # pdftools::pdf_convert() silently ignores the exact `filenames` value
  # whenever the source PDF has more than one page: for a single-page PDF it
  # writes exactly `png_file`, but for a multi-page one (e.g. switchPlot()'s
  # multi-panel output) it instead writes "<name>_1.png", "<name>_2.png", ...
  # and returns the vector of paths it actually created. The previous code
  # ignored that return value and only ever checked file.exists(png_file) --
  # always FALSE for a multi-page source, so every such figure was silently
  # treated as "couldn't convert" and never embedded, regardless of whether
  # pdftools was installed or the conversion itself worked fine.
  produced <- tryCatch({
    pdftools::pdf_convert(
      pdf = pdf_file,
      filenames = png_file,
      dpi = dpi,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("Failed to convert ", pdf_file, " to PNG: ", conditionMessage(e))
    character(0)
  })

  existing <- produced[file.exists(produced)]
  if (length(existing) == 0) return(NULL)
  if (length(existing) > 1) {
    message("   -> ", basename(pdf_file), " has ", length(existing),
            " pages; using page 1 (", basename(existing[1]), ") for report embedding.")
  }

  normalizePath(existing[1], winslash = "/", mustWork = FALSE)
}

#' Import transcript-level counts for isoform analysis
#'
#' @param data_dir Folder with quantification files
#' @param sample_table Path to sample metadata
#' @param ensembl_package_name EnsDb package name
#' @param count_type Quantification type ("salmon", "kallisto", etc.)
#' @param matrix_file Path to raw counts matrix (if count_type = "matrix")
#' @param subset_sample Optional filter expression
#' @param remove_sample Optional sample IDs to exclude
#' @param custom_tx2gene Optional path to a custom tx2gene file (TSV with columns 'tx_id' and 'gene_id')
#' @param custom_gene_map Optional path to a custom gene annotation file (TSV with columns 'gene_id', 'symbol', and optionally 'entrezid')
#' @return List with txi (transcript counts), meta, tx2gene, gene_map
#' @export
import_transcript_counts <- function(data_dir, sample_table, ensembl_package_name,
                                      count_type = "salmon", matrix_file = NULL,
                                      subset_sample = NULL, remove_sample = NULL,
                                      custom_tx2gene = NULL, custom_gene_map = NULL) {

  if (!file.exists(sample_table)) stop("Sample table not found: ", sample_table)
  sample_df <- data.table::fread(sample_table, header = TRUE, data.table = FALSE)
  sample_col <- if ("Sample" %in% colnames(sample_df)) "Sample" else "sample_id"

  sample_df <- .apply_sample_filters(sample_df, sample_col, remove_sample, subset_sample)

  rownames(sample_df) <- sample_df[[sample_col]]
  edb <- getExportedValue(ensembl_package_name, ensembl_package_name)

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
    if (!("gene_id" %in% colnames(gene_map)) && "ensembl" %in% colnames(gene_map)) {
      colnames(gene_map)[colnames(gene_map) == "ensembl"] <- "gene_id"
    }
    if (!all(c("gene_id", "symbol") %in% colnames(gene_map))) {
      stop("Custom gene map must contain columns 'gene_id' (or 'ensembl') and 'symbol'")
    }
    gene_map$gene_id <- strip_ensembl_version(gene_map$gene_id)
    colnames(gene_map)[colnames(gene_map) == "gene_id"] <- "ensembl"
    if (!"entrezid" %in% colnames(gene_map)) {
      gene_map$entrezid <- NA_character_
    }
    gene_map <- gene_map[, c("ensembl", "symbol", "entrezid")]
    gene_map <- gene_map[!duplicated(gene_map$ensembl), ]
    gene_map$symbol[is.na(gene_map$symbol) | gene_map$symbol == ""] <- gene_map$ensembl[is.na(gene_map$symbol) | gene_map$symbol == ""]

    # ---- ENHANCED FALLBACK: fill missing Entrez using bitr ----
    if (!is.null(org_obj)) {
      gene_map <- .fill_entrez_with_bitr(gene_map, org_obj, id_col = "ensembl", symbol_col = "symbol")
    }
    
    # ---- FINAL CHECK: report Entrez coverage ----
    entrez_present <- sum(!is.na(gene_map$entrezid) & gene_map$entrezid != "")
    message("  Gene map loaded: ", nrow(gene_map), " genes, ", entrez_present, " with Entrez IDs")
    
  } else {
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

  # ---- Continue with tximport or matrix ----
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
    # clean_transcript_id() (not the bare strip_ensembl_version()) -- kallisto/
    # salmon target_ids can be GENCODE-style pipe-delimited compound headers
    # ("ENST00000456328.2|ENSG00000223972.5|..."), and strip_ensembl_version()
    # alone is a silent no-op on those (see clean_transcript_id() in
    # utils_core.R). Without this, rownames(txi$counts) never matches
    # tx2gene$tx_id and .prepare_dtu_counts() fails with "No transcripts
    # could be mapped to genes after version handling."
    rownames(txi$counts) <- clean_transcript_id(rownames(txi$counts))
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
    rownames(count_mat) <- clean_transcript_id(rownames(count_mat))
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
run_dte <- function(isoform_obj, condition, level, base, padj_cutoff = 0.05, bpparam = NULL) {
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
  # parallel=TRUE is a no-op unless BPPARAM specifies >1 worker, so this is
  # always safe to pass; on large transcript-level datasets (far more rows
  # than gene-level DGE) DESeq2's per-row GLM fitting is the dominant cost
  # and benefits from the same multi-core BiocParallel param already
  # threaded through run_dtu() / run_dexseq_dtu() for consistency.
  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam
  dds <- DESeq2::DESeq(dds, test = "Wald", parallel = TRUE, BPPARAM = bp_param)
  res <- DESeq2::results(dds, contrast = c(condition, level, base))
  
  res_df <- as.data.frame(res)
  res_df$transcript_id <- rownames(res_df)

  # Defensively re-clean both sides of the join, mirroring
  # .prepare_dtu_counts()'s identical guard (see the comment there): IDs are
  # already cleaned once at import_transcript_counts() time, but this makes
  # run_dte() idempotent-safe against the same failure modes DTU already
  # guards against -- a hand-built/resumed isoform_obj, or a checkpointed
  # isoform_import.rds predating a cleaning fix -- rather than being the one
  # analysis path that silently trusts un-re-verified IDs. Without this,
  # rownames(dds) can still carry a raw GENCODE pipe-delimited FASTA header
  # (tximport takes salmon's quant.sf "Name" field verbatim, which is the
  # *entire* header when the index was built from a standard GENCODE cDNA
  # FASTA), the merge below matches nothing, gene_id/gene_symbol come out
  # NA for every affected row, and every plot/table built from dte_results
  # downstream (DTE_top_barplot, the DTE volcano, DTE_results_annotated.csv)
  # displays that raw header instead of a transcript ID.
  res_df$transcript_id   <- clean_transcript_id(res_df$transcript_id)
  tx2gene_clean           <- isoform_obj$tx2gene
  tx2gene_clean$tx_id     <- clean_transcript_id(tx2gene_clean$tx_id)
  res_df <- merge(res_df, tx2gene_clean, by.x = "transcript_id", by.y = "tx_id", all.x = TRUE)
  
  # --- Add gene symbol and Entrez ID ---
  res_df$gene_symbol <- isoform_obj$gene_map$symbol[match(res_df$gene_id, isoform_obj$gene_map$ensembl)]
  res_df$entrezid <- isoform_obj$gene_map$entrezid[match(res_df$gene_id, isoform_obj$gene_map$ensembl)]
  # Display-safe column: never bare "NA" (paste0(NA, ...) would otherwise
  # render literally as "NA") for transcripts with no resolvable gene
  # symbol -- falls back to the (now-clean) transcript_id. gene_symbol
  # itself is left as a real NA, since that's the semantically accurate
  # annotation for a genuinely novel/unannotated locus.
  res_df$gene <- .coalesce_gene_label(res_df$gene_symbol, res_df$transcript_id)
  
  res_df$signif <- !is.na(res_df$padj) & res_df$padj < padj_cutoff &
                   !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) > 1
  return(res_df)
}

#' Transcript-level PCA (companion to run_eda()'s gene-level PCA)
#'
#' Runs the same PCA plot + per-sample-scores + gene(transcript)-loadings
#' table as \code{run_eda()}, but on transcript-level counts, by reusing
#' \code{plot_custom_pca()} / \code{.write_pca_scores()} directly rather
#' than duplicating that logic for a second feature type.
#'
#' Uses \code{DESeq2::vst()} rather than \code{rlog()} for the
#' variance-stabilizing transform, unlike the gene-level PCA in
#' \code{run_eda()}: a transcript-level count matrix is typically several
#' times larger (tens of thousands of transcripts vs. ~20k genes), and
#' \code{rlog()}'s per-feature Bayesian shrinkage fit does not scale well to
#' that size, whereas \code{vst()}'s closed-form dispersion-trend fit gives
#' near-identical results for datasets this size in a fraction of the time
#' (this is DESeq2's own documented recommendation once sample/feature count
#' is large -- \code{rlog()} is kept mainly for small-n compatibility).
#'
#' @param isoform_obj Output from import_transcript_counts()
#' @param condition Column name in isoform_obj$meta to colour PCA points by;
#'   pass NULL for an ungrouped PCA
#' @param level,base Optional comparison labels, used only for output file
#'   naming (e.g. "PCA_transcripts_Treated_vs_Control.pdf"); pass NULL for
#'   an all-samples PCA not tied to one comparison
#' @param out_dir Output directory (a "Plots" subdirectory is created if needed)
#' @param batch_col Optional column in isoform_obj$meta for a shape aesthetic.
#'   No batch *correction* is attempted at the transcript level (unlike the
#'   optional limma step in run_eda()) -- pass a pre-corrected isoform_obj if
#'   you need that.
#' @param pca_ntop Number of most variable transcripts to use (default 500;
#'   NULL uses all transcripts that passed the expression filter)
#' @return Invisibly, the list returned by plot_custom_pca() (plot, gene_info,
#'   pca_scores)
#' @export
run_isoform_pca <- function(isoform_obj, condition, level = NULL, base = NULL,
                            out_dir, batch_col = NULL, pca_ntop = 500) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for transcript-level PCA. Please install it.")
  }
  if (isoform_obj$type == "tximport") {
    dds <- DESeq2::DESeqDataSetFromTximport(isoform_obj$txi, colData = isoform_obj$meta, design = ~1)
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(isoform_obj$counts),
                                          colData = isoform_obj$meta, design = ~1)
  }
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds  <- dds[keep, ]
  if (nrow(dds) < 2) {
    message("   -> Fewer than 2 transcripts pass the expression filter; skipping transcript-level PCA.")
    return(invisible(NULL))
  }

  vsd <- DESeq2::vst(dds, blind = TRUE)

  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  comp_label <- if (!is.null(level) && !is.null(base)) paste0(level, "_vs_", base) else "AllSamples"
  cond_arg   <- if (!is.null(condition) && condition %in% colnames(isoform_obj$meta)) condition else NULL

  res <- plot_custom_pca(vsd, condition = cond_arg, batch = batch_col,
                        title = paste0("Transcript-level PCA",
                                       if (!is.null(cond_arg)) paste0(" (", cond_arg, ")") else ""),
                        return_plot = TRUE, return_gene_list = TRUE, ntop = pca_ntop)

  pdf(file.path(plot_dir, paste0("PCA_transcripts_", comp_label, ".pdf")), width = 9, height = 7)
  print(res$plot)
  dev.off()

  if (!is.null(res$gene_info) && nrow(res$gene_info) > 0) {
    # plot_custom_pca()'s loadings table names its identifier column "gene"
    # regardless of feature type; rename for clarity since these rows are
    # transcripts, and clean the ID the same way every other transcript_id
    # in this pipeline is cleaned (rownames(dds) here are already clean --
    # inherited from isoform_obj -- but re-cleaning is idempotent and cheap,
    # the same defensive guard used throughout the isoform module).
    tx_info <- res$gene_info
    tx_info$gene <- clean_transcript_id(tx_info$gene)
    colnames(tx_info)[colnames(tx_info) == "gene"] <- "transcript_id"
    write.table(tx_info,
                file = file.path(plot_dir, paste0("PCA_transcripts_", comp_label, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    message("   -> Saved top variable transcripts used in PCA to: ",
            file.path(plot_dir, paste0("PCA_transcripts_", comp_label, ".tsv")))
  }
  tx_values <- res$gene_values
  if (!is.null(tx_values) && nrow(tx_values) > 0) {
    tx_values$gene <- clean_transcript_id(tx_values$gene)
    colnames(tx_values)[colnames(tx_values) == "gene"] <- "transcript_id"
  }
  .write_pca_scores(res$pca_scores, plot_dir, paste0("PCA_transcripts_", comp_label), gene_values = tx_values)

  message("Transcript-level PCA completed. Plots saved in: ", plot_dir)
  invisible(res)
}

#' Shared count-matrix preparation & filtering for DTU testing (DRIMSeq/DEXSeq)
#'
#' \code{run_dtu()} (DRIMSeq engine) and \code{run_dexseq_dtu()} (DEXSeq
#' engine) both start from the same raw transcript count matrix and apply an
#' identical sequence of sample-restriction, gene-mapping, and
#' expression-based filters before handing off to their engine-specific
#' fitting step. This was previously ~60 duplicated lines per function (a
#' change to one filter silently wouldn't apply to the other engine) --
#' consolidated here as the single source of truth. Only the two genuinely
#' engine-specific filters -- DEXSeq's manual gene-level expression filter
#' (DRIMSeq's own \code{dmFilter()} applies its equivalent later, per-chunk,
#' so run_dtu() must NOT apply it here too) and DEXSeq's "genes need >= 2
#' features" requirement -- are optional steps controlled by arguments, so
#' each engine's prior filtering behaviour is preserved exactly.
#'
#' @param isoform_obj Output from import_transcript_counts()
#' @param condition,level,base As in run_dtu()/run_dexseq_dtu()
#' @param min_transcript_total,min_transcript_expr,min_samps_feature_expr,
#'   min_gene_expr,min_samps_gene_expr,max_transcripts As in
#'   run_dtu()/run_dexseq_dtu()
#' @param apply_gene_expr_filter Logical: apply the manual gene-level
#'   min_gene_expr/min_samps_gene_expr filter here (DEXSeq needs this;
#'   DRIMSeq does not, since dmFilter() does it internally)
#' @param require_multi_transcript Logical: drop genes left with only one
#'   transcript after filtering (DEXSeq needs >= 2 features per group)
#' @return List with \code{counts}, \code{gene_id_map}, \code{tx_ids_original},
#'   \code{sample_data}, and the (possibly num_samples-capped)
#'   \code{min_samps_gene_expr}/\code{min_samps_feature_expr} actually used
#' @keywords internal
.prepare_dtu_counts <- function(isoform_obj, condition, level, base,
                                min_transcript_total, min_transcript_expr,
                                min_samps_feature_expr, min_gene_expr,
                                min_samps_gene_expr, max_transcripts,
                                apply_gene_expr_filter = FALSE,
                                require_multi_transcript = FALSE) {

  counts <- if (isoform_obj$type == "tximport") isoform_obj$txi$counts else isoform_obj$counts
  sample_data <- isoform_obj$meta

  # --- Keep only samples belonging to the two conditions of interest ---
  keep_samples <- sample_data[[condition]] %in% c(base, level)
  if (sum(keep_samples) < 2) {
    stop("Fewer than 2 samples available for comparison. Need at least one sample in each condition.")
  }
  sample_data <- sample_data[keep_samples, , drop = FALSE]
  counts <- counts[, rownames(sample_data), drop = FALSE]
  sample_data$condition <- factor(sample_data[[condition]], levels = c(base, level))

  num_samples <- ncol(counts)
  if (is.null(min_samps_gene_expr)) min_samps_gene_expr <- ceiling(num_samples * 0.5)
  min_samps_gene_expr    <- min(min_samps_gene_expr, num_samples)
  min_samps_feature_expr <- min(min_samps_feature_expr, num_samples)

  tx2gene <- isoform_obj$tx2gene
  # Clean transcript IDs: remove version suffixes AND any pipe/space-
  # delimited compound-header text (already done at import time in
  # import_transcript_counts(), but re-applied here too since it's
  # idempotent and cheap -- guards against a hand-built isoform_obj that
  # bypassed that step). Must use clean_transcript_id(), not the bare
  # strip_ensembl_version() -- the latter is a silent no-op on GENCODE-style
  # compound IDs like "ENST00000456328.2|ENSG00000223972.5|...", which
  # previously left this match() 100% unmatched with no partial-mismatch
  # warning to flag it.
  rownames(counts) <- clean_transcript_id(rownames(counts))
  tx2gene$tx_id    <- clean_transcript_id(tx2gene$tx_id)
  gene_id_map <- tx2gene$gene_id[match(rownames(counts), tx2gene$tx_id)]

  # Filter out transcripts with no gene mapping
  keep_mapped <- !is.na(gene_id_map)
  counts <- counts[keep_mapped, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_mapped]
  tx_ids_original <- rownames(counts)
  if (nrow(counts) == 0) stop("No transcripts could be mapped to genes after version handling.")
  message("Mapped transcripts: ", nrow(counts))

  # Filter by total expression per transcript
  keep_tx <- rowSums(counts) >= min_transcript_total
  counts <- counts[keep_tx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_tx]
  tx_ids_original <- tx_ids_original[keep_tx]
  if (nrow(counts) == 0) stop("No transcripts passed min_transcript_total filter.")

  # Gene-level expression filter -- only for engines without their own
  # internal equivalent (mirrors DRIMSeq::dmFilter()'s min_gene_expr /
  # min_samps_gene_expr semantics: require a minimum number of samples where
  # the gene's TOTAL transcript expression reaches min_gene_expr).
  if (isTRUE(apply_gene_expr_filter)) {
    gene_expr_matrix <- rowsum(counts, group = gene_id_map)
    keep_genes_expr  <- rownames(gene_expr_matrix)[
      rowSums(gene_expr_matrix >= min_gene_expr) >= min_samps_gene_expr
    ]
    keep_tx <- gene_id_map %in% keep_genes_expr
    counts <- counts[keep_tx, , drop = FALSE]
    gene_id_map <- gene_id_map[keep_tx]
    tx_ids_original <- tx_ids_original[keep_tx]
    if (nrow(counts) == 0) stop("No genes passed the min_gene_expr / min_samps_gene_expr filter.")
  }

  # Filter by minimum proportion and number of samples with expression
  keep_tx <- rowSums(counts > min_transcript_expr) >= min_samps_feature_expr
  counts <- counts[keep_tx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_tx]
  tx_ids_original <- tx_ids_original[keep_tx]
  if (nrow(counts) == 0) stop("No transcripts passed expression filtering.")

  # Cap number of transcripts per gene
  gene_split <- split(seq_len(nrow(counts)), gene_id_map)
  keep_idx <- unlist(lapply(gene_split, function(idx) {
    if (length(idx) <= max_transcripts) return(idx)
    gene_expr <- rowMeans(counts[idx, , drop = FALSE])
    idx[order(gene_expr, decreasing = TRUE)[seq_len(max_transcripts)]]
  }))
  counts <- counts[keep_idx, , drop = FALSE]
  gene_id_map <- gene_id_map[keep_idx]
  tx_ids_original <- tx_ids_original[keep_idx]

  # DEXSeq needs >= 2 features per group to test differential usage; DRIMSeq
  # can in principle fit a single-feature gene (it's just uninformative), so
  # run_dtu() does not request this filter.
  if (isTRUE(require_multi_transcript)) {
    gene_counts_tbl <- table(gene_id_map)
    multi_tx_genes  <- names(gene_counts_tbl)[gene_counts_tbl > 1]
    keep_multi      <- gene_id_map %in% multi_tx_genes
    counts          <- counts[keep_multi, , drop = FALSE]
    gene_id_map     <- gene_id_map[keep_multi]
    tx_ids_original <- tx_ids_original[keep_multi]
    if (nrow(counts) == 0) stop("No multi-transcript genes remained after filtering.")
  }

  message("Final transcripts after filtering: ", nrow(counts))
  message("Genes retained: ", length(unique(gene_id_map)))

  list(counts = counts, gene_id_map = gene_id_map, tx_ids_original = tx_ids_original,
       sample_data = sample_data, min_samps_gene_expr = min_samps_gene_expr,
       min_samps_feature_expr = min_samps_feature_expr)
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
#' @return List with \code{dtu_results}, a TRANSCRIPT-level (feature-level)
#'   data.frame with columns \code{gene_id}, \code{feature_id}, \code{lr},
#'   \code{df}, \code{pvalue}, and a genome-wide-recomputed \code{adj_pvalue},
#'   plus \code{gene_symbol}/\code{entrezid}. Uses
#'   \code{DRIMSeq::results(d, level = "feature")} rather than DRIMSeq's
#'   default gene-level result: (1) that's the whole point of a
#'   *transcript*-usage test, and (2) \code{run_isoform_switch()}'s
#'   \code{test_engine = "DRIMSeq"} path requires a \code{feature_id} column
#'   that only the feature-level result has -- without this, that path could
#'   never succeed (see \code{required_cols} check there).
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

  prep <- .prepare_dtu_counts(
    isoform_obj, condition, level, base,
    min_transcript_total   = min_transcript_total,
    min_transcript_expr    = min_transcript_expr,
    min_samps_feature_expr = min_samps_feature_expr,
    min_gene_expr           = min_gene_expr,
    min_samps_gene_expr     = min_samps_gene_expr,
    max_transcripts         = max_transcripts,
    apply_gene_expr_filter   = FALSE,  # DRIMSeq's own dmFilter() applies this per-chunk below
    require_multi_transcript = FALSE
  )
  counts             <- prep$counts
  gene_id_map        <- prep$gene_id_map
  tx_ids_original    <- prep$tx_ids_original
  sample_data        <- prep$sample_data
  min_samps_gene_expr    <- prep$min_samps_gene_expr
  min_samps_feature_expr <- prep$min_samps_feature_expr

  unique_genes <- unique(gene_id_map)
  gene_chunks <- split(unique_genes,
                       ceiling(seq_along(unique_genes) / chunk_size))

  all_results <- list()

  # Clean sample IDs in metadata and counts
  rownames(sample_data) <- trimws(rownames(sample_data))
  colnames(counts) <- trimws(colnames(counts))

  for (i in seq_along(gene_chunks)) {
    message("Processing chunk ", i, " / ", length(gene_chunks))

    current_genes <- gene_chunks[[i]]
    idx <- which(gene_id_map %in% current_genes)

    curr_counts <- counts[idx, , drop = FALSE]
    curr_gene_ids <- gene_id_map[idx]
    curr_tx_ids <- tx_ids_original[idx]

    # Find common samples between current counts and metadata
    common_samples <- intersect(colnames(curr_counts), rownames(sample_data))
    if (length(common_samples) == 0) {
      warning("Chunk ", i, ": No common samples between counts and metadata. Skipping.")
      next
    }

    # Subset both count matrix and metadata to common samples
    curr_counts <- curr_counts[, common_samples, drop = FALSE]
    sample_data_sub <- sample_data[common_samples, , drop = FALSE]

    # Prepare data frame for dmDSdata
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

    # ----- Robust design matrix construction -----
    counts_d <- DRIMSeq::counts(d)
    samples_d <- DRIMSeq::samples(d)

    # Ensure samples_d rows are in the same order as the columns of counts_d
    if (!identical(rownames(samples_d), colnames(counts_d))) {
      samples_d <- samples_d[colnames(counts_d), , drop = FALSE]
    }

    # Check dimensions
    if (nrow(samples_d) != ncol(counts_d)) {
      warning("Chunk ", i, ": Design matrix rows (", nrow(samples_d), ") do not match number of samples (", ncol(counts_d), "). Skipping.")
      next
    }

    design <- model.matrix(~ condition, data = samples_d)

    res <- tryCatch({
      d <- DRIMSeq::dmPrecision(d, design = design, BPPARAM = bp_param)
      d <- DRIMSeq::dmFit(d, design = design, BPPARAM = bp_param)
      d <- DRIMSeq::dmTest(d, coef = paste0("condition", level), BPPARAM = bp_param)
      # level = "feature": DRIMSeq's default (level = "gene") returns ONE ROW
      # PER GENE with no feature_id column at all -- silently discarding the
      # transcript-level detail that is the entire point of a *transcript*
      # usage test, and leaving run_isoform_switch()'s test_engine = "DRIMSeq"
      # path (which requires feature_id, see required_cols check there)
      # unable to ever succeed. Confirmed against DRIMSeq's own
      # dmDStest-class docs: results_gene has gene_id/lr/df/pvalue/adj_pvalue;
      # results_feature adds feature_id and is only returned by
      # results(d, level = "feature").
      DRIMSeq::results(d, level = "feature")
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

  if (length(all_results) == 0) {
    stop("No DTU results were produced. Check that your samples and metadata match, and that there are genes with multiple transcripts meeting the filtering criteria.")
  }

  dtu_results <- do.call(rbind, all_results)
  rm(all_results)
  gc()

  pcol <- intersect(c("pvalue", "p_value"), colnames(dtu_results))[1]
  if (is.na(pcol)) {
    stop("No p-value column found in results.")
  }
  # Recomputed genome-wide (across all chunks) at the TRANSCRIPT level now,
  # rather than DRIMSeq's own per-chunk, gene-level adj_pvalue.
  dtu_results$adj_pvalue <- p.adjust(dtu_results[[pcol]], method = "BH")
  
  # --- Add gene symbol and Entrez ID to results ---
  dtu_results$gene_symbol <- isoform_obj$gene_map$symbol[match(dtu_results$gene_id, isoform_obj$gene_map$ensembl)]
  dtu_results$entrezid <- isoform_obj$gene_map$entrezid[match(dtu_results$gene_id, isoform_obj$gene_map$ensembl)]
  # Display-safe column -- see .coalesce_gene_label() / run_dte()'s $gene for
  # why this exists: falls back to feature_id (transcript) rather than
  # gene_id here since that's what generate_dte_dtu_report()'s DTU barplot
  # actually labels points by.
  dtu_results$gene <- .coalesce_gene_label(dtu_results$gene_symbol, dtu_results$feature_id)

  return(list(dtu_results = dtu_results))
}

# ==============================================================================
# DEXSeq-based DTU (complementary engine to run_dtu / DRIMSeq)
# ==============================================================================
#
# Follows the "transcript-level DEXSeq" workflow described by Soneson, Love &
# Robinson (F1000Research, 2016; "Swimming downstream: statistical analysis of
# differential transcript usage following Salmon quantification"): each
# transcript is treated as a DEXSeq exonic "bin" within its parent gene, and
# DEXSeq's negative-binomial GLM is used to test for differential usage. This
# is independent from the DRIMSeq engine in run_dtu() and is commonly run
# alongside it as a cross-validation of DTU calls. Crucially it also exposes
# the native DEXSeqResults object, which unlocks DEXSeq::plotDEXSeq() (wrapped
# below as plot_dexseq_gene()) for classic per-transcript expression/usage
# plots across conditions.

#' Run DEXSeq-based Differential Transcript Usage
#'
#' @param isoform_obj Output from import_transcript_counts()
#' @param condition Column name in metadata containing groups
#' @param level Foreground level for contrast
#' @param base Reference level
#' @param min_gene_expr Minimum total expression per gene (default: 10)
#' @param min_transcript_expr Minimum transcript proportion per gene (default: 0.05)
#' @param min_samps_gene_expr Min samples for gene expression (default: 50% of samples)
#' @param min_samps_feature_expr Min samples for transcript expression (default: 3)
#' @param chunk_size Number of genes to process at once (default: 2000)
#' @param max_transcripts Maximum number of transcripts per gene to keep (default: 300)
#' @param min_transcript_total Minimum total counts across all samples for a transcript (default: 10)
#' @param keep_dxr Logical: retain the per-chunk DEXSeqResults objects so
#'   individual genes can later be visualised with plot_dexseq_gene() (default TRUE).
#'   Set FALSE to reduce memory usage for very large datasets when only the flat
#'   results table is needed.
#' @param bpparam Optional BiocParallel parameter object for DEXSeq operations
#' @return List with `results_df` (annotated per-transcript data.frame including
#'   pvalue/padj plus a Simes-aggregated, BH-corrected gene-level q-value) and
#'   `dxr_list` (list of per-chunk DEXSeqResults objects, or an empty list if
#'   keep_dxr = FALSE)
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
  bpparam = NULL
) {

  if (!requireNamespace("DEXSeq", quietly = TRUE)) {
    stop("DEXSeq is required for run_dexseq_dtu(). Install with BiocManager::install('DEXSeq').")
  }

  bp_param <- if (is.null(bpparam)) BiocParallel::SerialParam() else bpparam
  BiocParallel::register(bp_param)

  # (mirrors run_dtu(); without restricting to base/level samples first, any
  # OTHER condition level gets an NA 'condition' factor below -- since
  # factor() only knows about c(base, level) -- and DEXSeqDataSet() then
  # aborts with "variables in design formula cannot contain NA: condition")
  prep <- .prepare_dtu_counts(
    isoform_obj, condition, level, base,
    min_transcript_total   = min_transcript_total,
    min_transcript_expr    = min_transcript_expr,
    min_samps_feature_expr = min_samps_feature_expr,
    min_gene_expr           = min_gene_expr,
    min_samps_gene_expr     = min_samps_gene_expr,
    max_transcripts         = max_transcripts,
    apply_gene_expr_filter   = TRUE,  # DEXSeq has no built-in equivalent of dmFilter()
    require_multi_transcript = TRUE   # DEXSeq needs >= 2 features per group to test
  )
  counts             <- prep$counts
  gene_id_map        <- prep$gene_id_map
  tx_ids_original    <- prep$tx_ids_original
  sample_data        <- prep$sample_data

  message("DEXSeq DTU: ", nrow(counts), " transcripts across ",
          length(unique(gene_id_map)), " genes.")

  unique_genes <- unique(gene_id_map)
  gene_chunks  <- split(unique_genes, ceiling(seq_along(unique_genes) / chunk_size))

  # colData for DEXSeqDataSet only needs the 'condition' column; the 'sample'
  # and 'exon' terms referenced in the design formula are generated internally.
  sample_meta <- data.frame(
    condition = sample_data[colnames(counts), "condition"],
    row.names = colnames(counts)
  )

  all_results <- list()
  dxr_list    <- list()

  for (i in seq_along(gene_chunks)) {
    message("DEXSeq DTU: processing chunk ", i, " / ", length(gene_chunks))

    current_genes <- gene_chunks[[i]]
    idx <- which(gene_id_map %in% current_genes)

    curr_counts   <- round(as.matrix(counts[idx, , drop = FALSE]))
    curr_gene_ids <- gene_id_map[idx]
    curr_tx_ids   <- tx_ids_original[idx]

    dxd <- tryCatch({
      DEXSeq::DEXSeqDataSet(
        countData  = curr_counts,
        sampleData = sample_meta,
        design     = ~sample + exon + condition:exon,
        featureID  = curr_tx_ids,
        groupID    = curr_gene_ids
      )
    }, error = function(e) {
      message("  Chunk ", i, " DEXSeqDataSet construction failed: ", e$message)
      NULL
    })

    if (is.null(dxd)) next

    res <- tryCatch({
      dxd <- DEXSeq::estimateSizeFactors(dxd)
      dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE, BPPARAM = bp_param)
      dxd <- DEXSeq::testForDEU(dxd, BPPARAM = bp_param)
      dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = bp_param)
      DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
    }, error = function(e) {
      message("  Chunk ", i, " DEXSeq testing failed: ", e$message)
      NULL
    })

    if (!is.null(res)) {
      res_df <- as.data.frame(res)
      wanted <- intersect(
        c("groupID", "featureID", "exonBaseMean", "dispersion", "stat", "pvalue", "padj",
          grep("^log2fold_", colnames(res_df), value = TRUE)),
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

  # Aggregate per-transcript p-values into a per-gene p-value via Simes' method,
  # then BH-adjust across all genes -- the same two-step approach used by
  # DEXSeq's own perGeneQValue(), reimplemented here so it stays correct when
  # results are combined across processing chunks.
  pval_by_gene <- split(results_df$pvalue, results_df$groupID)
  gene_pvals <- vapply(pval_by_gene, function(p) {
    p <- p[!is.na(p)]
    if (length(p) == 0) return(NA_real_)
    n <- length(p)
    min(sort(p) * n / seq_len(n))
  }, numeric(1))
  gene_qvals <- stats::p.adjust(gene_pvals, method = "BH")

  results_df$gene_pvalue_simes <- gene_pvals[results_df$groupID]
  results_df$gene_qvalue       <- gene_qvals[results_df$groupID]

  colnames(results_df)[colnames(results_df) == "featureID"] <- "transcript_id"
  colnames(results_df)[colnames(results_df) == "groupID"]   <- "gene_id"
  results_df$gene_symbol <- isoform_obj$gene_map$symbol[match(results_df$gene_id, isoform_obj$gene_map$ensembl)]
  # Display-safe column -- see .coalesce_gene_label() / run_dte()'s $gene.
  results_df$gene <- .coalesce_gene_label(results_df$gene_symbol, results_df$transcript_id)

  n_sig_genes <- sum(gene_qvals < 0.05, na.rm = TRUE)
  message("DEXSeq DTU complete: ", length(unique(results_df$gene_id)), " genes tested, ",
          n_sig_genes, " significant at gene q-value < 0.05.")

  list(results_df = results_df, dxr_list = dxr_list)
}

#' DEXSeq-style Transcript Usage Plot for a Single Gene
#'
#' Thin wrapper around \code{DEXSeq::plotDEXSeq()} producing the classic
#' per-feature (transcript-as-exon) expression-by-condition plot, with the
#' feature(s) showing significant differential usage highlighted. This is the
#' "DEXSeq plot" style requested for comparing transcript/exon usage across
#' conditions.
#'
#' @param dxr_list List of DEXSeqResults objects, as returned in the `dxr_list`
#'   element of run_dexseq_dtu() (requires keep_dxr = TRUE there)
#' @param gene_id Ensembl gene ID to plot (must be present in one of the DEXSeqResults chunks)
#' @param plot_dir Output directory for the PDF
#' @param gene_symbol Optional display name used for the output filename
#' @param splicing If TRUE, plots usage with the overall gene expression effect
#'   averaged out (isolates the usage/splicing signal); default FALSE shows
#'   fitted expression per transcript per condition
#' @return Invisibly, the path to the generated PDF (or NULL if the gene wasn't
#'   found in any chunk or plotting failed)
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

  label    <- if (!is.null(gene_symbol) && nzchar(gene_symbol)) gene_symbol else gene_id
  pdf_path <- file.path(plot_dir, paste0("DEXSeq_", label, ".pdf"))

  plotted <- FALSE
  for (dxr in dxr_list) {
    if (is.null(dxr)) next
    if (!(gene_id %in% dxr$groupID)) next
    ok <- tryCatch({
      pdf(pdf_path, width = 9, height = 6)
      DEXSeq::plotDEXSeq(dxr, geneID = gene_id, fitExpToVar = "condition",
                         expression = !splicing, splicing = splicing,
                         legend = TRUE, cex.axis = 1.1, cex = 1.2, lwd = 2)
      dev.off()
      TRUE
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      message("   -> plotDEXSeq() failed for ", label, ": ", e$message)
      FALSE
    })
    if (ok) { plotted <- TRUE; break }
  }

  if (!plotted) {
    message("   -> Gene '", label, "' not found in DEXSeq results (or not testable). Skipping plot.")
    return(invisible(NULL))
  }
  message("   -> DEXSeq transcript-usage plot saved to: ", pdf_path)
  invisible(pdf_path)
}

#' Run IsoformSwitchAnalyzeR analysis using isoformSwitchAnalysisCombined
#'
#' @param dte_results (unused, kept for compatibility)
#' @param dtu_results Required (and used) when \code{test_engine = "DRIMSeq"}:
#'   the \code{list(dtu_results = ...)} returned by \code{run_dtu()}, whose
#'   \code{feature_id} column is fed into IsoformSwitchAnalyzeR as
#'   pre-computed DTU p-values. Not used when \code{test_engine = "DEXSeq"}
#'   (that path calls \code{isoformSwitchTestDEXSeq()} itself instead).
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
#' @param predictor_cpu Integer: CPU threads for hmmscan / InterProScan
#'   (default NULL auto-detects with \code{parallel::detectCores() - 1})
#' @param log_dir Optional explicit directory for logs (WSL debug info, WSL
#'   command audit trail). Defaults to \code{file.path(out_dir, "Log")}. Set
#'   this when calling from a larger pipeline that already manages a single
#'   unified log tree elsewhere, so isoform logs land there instead of a
#'   second, separate location.
#' @param custom_transcript_id_map Path to a TSV/CSV with columns `count_id` and `fasta_id`.
#'   If provided, transcript IDs in the count matrix are renamed to match the FASTA file
#'   before any filtering; this resolves ID mismatches between custom references.
#' @param skip_fasta_filter Logical: if TRUE, skip the pre‑filtering that keeps only transcripts
#'   present in the FASTA file; rely on `importRdata()`'s own ID tolerance (e.g. `ignoreAfter*`).
#'   (Default: FALSE)
#' @param test_engine Which DTU test engine to use in IsoformSwitchAnalyzeR:
#'   `"DEXSeq"` (default, built‑in), `"DRIMSeq"` (uses pre‑computed results from `run_dtu()`),
#'   or `"satuRn"` (state‑of‑the‑art method). See `IsoformSwitchAnalyzeR::isoformSwitchTestSatuRn`.
#' @return IsoformSwitchAnalyzeR results object, or NULL if no genes were
#'   found to be switching (see the isoformSwitchAnalysisCombined() tryCatch
#'   below) -- callers should treat a NULL return as "no isoform switches
#'   detected", not as an error.
#' @export
run_isoform_switch <- function(dte_results = NULL, dtu_results = NULL,
                               isoform_obj, condition, level, base,
                               fasta_file, gff_file, out_dir,
                               run_predictors = FALSE,
                               use_wsl = (.Platform$OS.type == "windows"),
                               wsl_distro = "Ubuntu-22.04",
                               save_dir = NULL, resume_from = NULL,
                               bsgenome_name = NULL, predictor_cpu = NULL,
                               log_dir = NULL,
                               custom_transcript_id_map = NULL,
                               skip_fasta_filter = FALSE,
                               test_engine = c("DEXSeq", "DRIMSeq", "satuRn")) {

  test_engine <- match.arg(test_engine)

  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE))
    stop("Please install IsoformSwitchAnalyzeR: BiocManager::install('IsoformSwitchAnalyzeR')")
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Please install Biostrings: BiocManager::install('Biostrings')")

  if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
  if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  if (is.null(log_dir)) log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # --------------------------------------------------------------------------
  # Predictor environment check and logging (Windows+WSL or native Linux/macOS)
  # --------------------------------------------------------------------------
  if (isTRUE(run_predictors)) {
    message("Checking external predictor tool availability...")
    debug_info <- debug_wsl(distro = wsl_distro, out_dir = out_dir, log_dir = log_dir,
                            conda_env = "isoform_tools", verbose = TRUE,
                            use_wsl = use_wsl)
    if (!isTRUE(debug_info$wsl_available)) {
      stop(if (.Platform$OS.type == "windows" && use_wsl)
             "WSL is not available. Please install WSL or set use_wsl = FALSE."
           else
             paste0("Could not execute external predictor tools in this environment ",
                    "(no working bash shell found)."))
    }
    if (!isTRUE(debug_info$conda_env_exists)) {
      warning("Conda environment 'isoform_tools' not found. Predictors may still work if the ",
              "required tools (CPAT, SignalP, hmmscan/InterProScan) are already on PATH.\n",
              "On Windows/WSL run install_wsl_isoform_tools() to set up the environment; ",
              "on Linux/macOS install these tools manually or via conda/mamba.")
    }
    if (isTRUE(debug_info$pfam_db_found) && !isTRUE(debug_info$pfam_db_indexed)) {
      warning("Pfam-A.hmm was found but is not hmmpress-indexed -- hmmscan will fail against it. ",
              "Run install_isoform_databases() again (now conda-aware) or `hmmpress` it manually.")
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
    message("Building SwitchList from raw data...")

    if (isoform_obj$type == "tximport") {
      count_matrix <- isoform_obj$txi$counts
    } else {
      count_matrix <- isoform_obj$counts
    }

    sample_col <- if ("Sample" %in% colnames(isoform_obj$meta)) "Sample" else "sample_id"

    # --- Filter samples to only those in the comparison groups ---
    meta_sub <- isoform_obj$meta[isoform_obj$meta[[condition]] %in% c(base, level), , drop = FALSE]
    if (nrow(meta_sub) == 0) {
      stop("No samples found for the comparison groups: ", base, " and ", level)
    }
    # Ensure sample order matches count matrix columns
    keep_samples <- intersect(rownames(meta_sub), colnames(count_matrix))
    if (length(keep_samples) == 0) {
      stop("No matching sample IDs between metadata and count matrix after filtering.")
    }
    meta_sub <- meta_sub[keep_samples, , drop = FALSE]
    count_matrix <- count_matrix[, keep_samples, drop = FALSE]

    # Update sample vector
    sample_vector <- meta_sub[[sample_col]]
    # Ensure column names match the sample IDs
    colnames(count_matrix) <- sample_vector

    # Optionally warn if a group has fewer than 2 samples
    group_counts <- table(meta_sub[[condition]])
    if (any(group_counts < 2)) {
      warning("One or both groups have fewer than 2 samples. DTU/switch analysis may be unreliable.")
    }

    # ---- Custom ID mapping (NEW) ----
    if (!is.null(custom_transcript_id_map)) {
      if (is.character(custom_transcript_id_map) && file.exists(custom_transcript_id_map)) {
        map_df <- data.table::fread(custom_transcript_id_map, header = TRUE, data.table = FALSE)
      } else if (is.data.frame(custom_transcript_id_map)) {
        map_df <- custom_transcript_id_map
      } else {
        stop("custom_transcript_id_map must be a data.frame or a path to a TSV/CSV file.")
      }
      required_cols <- c("count_id", "fasta_id")
      if (!all(required_cols %in% colnames(map_df))) {
        stop("custom_transcript_id_map must contain columns: ", paste(required_cols, collapse = ", "))
      }
      # Rename count matrix rows
      old_rownames <- rownames(count_matrix)
      new_rownames <- map_df$fasta_id[match(old_rownames, map_df$count_id)]
      if (any(is.na(new_rownames))) {
        warning("Some count IDs were not found in the mapping file; they will be kept unchanged.")
        new_rownames[is.na(new_rownames)] <- old_rownames[is.na(new_rownames)]
      }
      rownames(count_matrix) <- new_rownames
      message("  Remapped ", sum(!is.na(new_rownames)), " transcript IDs using custom map.")
    }

    # ---- ID cleaning function (MUST BE DEFINED BEFORE USE) ----
    # Delegates to clean_transcript_id() in utils_core.R -- the single
    # source of truth for the bar/space-then-version-strip ordering (see that
    # function's docs for why the order matters and what silently broke
    # before: strip_ensembl_version() alone is a no-op on GENCODE-style
    # compound headers like "ENST00000456328.2|ENSG00000237683.5|...").
    # import_transcript_counts() and .prepare_dtu_counts() now use the same
    # helper, so IDs arriving here are typically already clean; this call
    # stays as a defensive no-op for any hand-built isoform_obj.
    clean_id <- clean_transcript_id

    # Persist the cleaned IDs onto the matrix itself -- not just into a local
    # variable used for the FASTA-membership test below. Every downstream
    # step in this function (the tx2gene/single-transcript-gene filter right
    # below, and the exact fasta/gtf ID match against
    # fasta_file_clean/gff_file_clean further down) reads
    # rownames(count_matrix) directly and assumes it's already clean.
    # Without this, that assumption silently fails whenever the incoming
    # isoform_obj isn't already clean -- e.g. an isoform_import.rds
    # checkpoint saved before clean_transcript_id() existed in
    # import_transcript_counts(). rownames(count_matrix) %in% tx2gene$tx_clean
    # then matches nothing, "Kept 0 transcripts from 0 genes" "succeeds"
    # with no error, and the resulting empty isoformRepExpression surfaces
    # many steps later as importRdata()'s unrelated-looking "contains
    # negative values" error (min()/max() of an empty vector return Inf/-Inf,
    # which is exactly the "no non-missing arguments" warning pair seen
    # right before that error).
    rownames(count_matrix) <- clean_id(rownames(count_matrix))

    # Pre-filter: keep only transcripts present in the FASTA file (unless skipped)
    if (!skip_fasta_filter) {
      message("Pre-filtering transcripts to match FASTA file...")
      fasta_seqs    <- Biostrings::fasta.seqlengths(fasta_file)
      clean_fasta_ids  <- clean_id(names(fasta_seqs))
      clean_rownames   <- clean_id(rownames(count_matrix))
      keep_in_fasta    <- clean_rownames %in% clean_fasta_ids
      count_matrix     <- count_matrix[keep_in_fasta, , drop = FALSE]
      message("  Kept ", nrow(count_matrix), " / ", length(clean_rownames),
              " transcripts matching the FASTA file.")
      if (nrow(count_matrix) == 0) {
        stop("No transcript IDs match between count matrix and FASTA file. ",
             "Check ID formats (pipe-delimited GENCODE headers, version suffixes, ",
             "or a genuinely different annotation/quantification source). ",
             "You may supply a custom_transcript_id_map or set skip_fasta_filter = TRUE.")
      }
    } else {
      message("Skipping FASTA pre-filtering (skip_fasta_filter = TRUE). ",
              "Will rely on importRdata()'s ignoreAfter* arguments.")
    }

    # Remove genes with only one transcript (DRIMSeq requirement)
    tx2gene      <- isoform_obj$tx2gene
    # Use cleaned IDs for consistency; if mapping was done, rownames already match
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

    isoform_count_matrix <- round(count_matrix)

    # Build design matrix from filtered metadata
    design_matrix <- data.frame(
      sampleID  = sample_vector,
      condition = factor(meta_sub[[condition]], levels = c(base, level)),
      stringsAsFactors = FALSE
    )
    rownames(design_matrix) <- design_matrix$sampleID

    if (!"package:dplyr" %in% search()) attachNamespace("dplyr")

    # ---- Sanity-check the fasta/gtf/count-matrix triplet BEFORE importRdata ----
    # importRdata()'s own error message for a scale mismatch here (a Jaccard
    # similarity check) is accurate but easy to misdiagnose as an ID-format
    # problem. Surface the same information earlier, in terms of the actual
    # input files, so a genuinely mismatched reference (e.g. a stale/different
    # SQANTI3 run) is obvious immediately rather than after minutes of DTU
    # fitting.
    .count_fasta_entries <- function(path) {
      tryCatch(length(Biostrings::fasta.seqlengths(path)), error = function(e) NA_integer_)
    }
    .count_gtf_transcripts <- function(path) {
      # Streamed line-by-line count of GTF "transcript" feature rows (column
      # 3) -- avoids loading a potentially very large (multi-GB, hundreds of
      # thousands of transcripts) GTF entirely into memory just to count it.
      tryCatch({
        con <- file(path, open = "rt")
        on.exit(close(con), add = TRUE)
        n <- 0L
        repeat {
          lines <- readLines(con, n = 200000L, warn = FALSE)
          if (length(lines) == 0L) break
          n <- n + sum(grepl("\ttranscript\t", lines, fixed = TRUE))
        }
        if (n == 0L) NA_integer_ else n
      }, error = function(e) NA_integer_)
    }
    n_fasta <- .count_fasta_entries(fasta_file[1])
    n_gtf   <- .count_gtf_transcripts(gff_file)
    n_count <- nrow(isoform_count_matrix)
    if (!is.na(n_fasta) && !is.na(n_gtf) && n_fasta > 0 && n_gtf > 0) {
      ratio <- min(n_fasta, n_gtf) / max(n_fasta, n_gtf)
      if (ratio < 0.5) {
        warning(
          "isoform_fasta and isoform_gff look like they describe very different ",
          "transcript sets (", n_fasta, " fasta sequences vs. ", n_gtf, " GTF transcripts; ",
          "count matrix after filtering has ", n_count, " transcripts). importRdata() will ",
          "likely fail its Jaccard-similarity check below. This almost always means the ",
          "fasta/gtf pair (and/or the tx2gene/count matrix) were built from different ",
          "pipeline runs -- e.g. a stale SQANTI3 output next to a newer one. Verify with, ",
          "e.g., `grep -c '^>' isoform_fasta` vs. `grep -c $'\\ttranscript\\t' isoform_gff`, ",
          "and make sure both come from the SAME sqanti3_rescue run as the kallisto/salmon ",
          "quantification tx2gene was built from.",
          immediate. = TRUE
        )
      }
    }

    # ---- importRdata ----
    # ignoreAfterBar/Space are safe to leave TRUE unconditionally: none of the
    # ID schemes seen here (Ensembl, StringTie2, SQANTI3) ever use '|' or ' '
    # *inside* a transcript ID itself -- only as description-field delimiters
    # -- so truncating at the first one only ever removes description text.
    #
    # ignoreAfterPeriod can't be handled with a single TRUE/FALSE flag here.
    # A SQANTI3 "rescued" fasta/gtf is typically a MIX of two ID styles in
    # the SAME file: transcripts that were rescued/matched to a known
    # reference keep that reference's (often versioned) Ensembl ID, while
    # genuinely novel transcripts get StringTie2/SQANTI3-style multi-dot IDs
    # (e.g. "PB.1.1", "MSTRG.6.1") where '.' is a structural separator, not a
    # version suffix.
    #   - ignoreAfterPeriod = TRUE truncates EVERY id at the first '.', so
    #     "PB.1.1" and "PB.1.2" (two different isoforms of the same locus)
    #     both collapse to "PB.1" -- importRdata()'s internal fixNames() then
    #     rejects the now-duplicated ID set ("...would cause IDs to be
    #     non-unique"). That's the crash this block used to hit.
    #   - ignoreAfterPeriod = FALSE avoids that collision but then silently
    #     stops stripping genuine Ensembl version suffixes too -- every
    #     rescued transcript that kept a versioned reference ID would fail to
    #     match the already version-stripped count matrix and get silently
    #     dropped from the analysis, with no error at all -- a worse outcome
    #     than the crash, since nothing would flag it.
    # Neither single setting is correct for a mixed file. The only correct
    # fix is to apply the SAME selective, per-ID logic as clean_id() /
    # strip_ensembl_version() above -- which strips a version suffix only
    # from IDs that actually look like Ensembl IDs, leaving PB./MSTRG-style
    # IDs untouched -- directly to the fasta headers and GTF transcript_id
    # values ourselves, so importRdata() receives IDs that already match
    # exactly and none of its own truncation is needed.
    clean_ref_dir <- file.path(out_dir, "cleaned_reference")
    if (!dir.exists(clean_ref_dir)) dir.create(clean_ref_dir, recursive = TRUE)
    fasta_file_clean <- file.path(clean_ref_dir, "isoform_fasta_ids_cleaned.fasta")
    gff_file_clean   <- file.path(clean_ref_dir, "isoform_gff_ids_cleaned.gtf")

    message("  Pre-cleaning transcript IDs in FASTA/GTF (same rule already used for the ",
            "count matrix above) so importRdata() gets exact matches instead of guessing...")

    # --- FASTA: rename headers only; Biostrings is already a hard dependency
    # of this function, and sequences themselves are untouched. ---
    seqs <- Biostrings::readDNAStringSet(fasta_file)
    names(seqs) <- clean_id(names(seqs))
    dup_fa <- unique(names(seqs)[duplicated(names(seqs))])
    if (length(dup_fa) > 0) {
      stop("Cleaning transcript IDs in isoform_fasta produced ", length(dup_fa),
           " duplicate ID(s) -- e.g. ", paste(utils::head(dup_fa, 5), collapse = ", "),
           ". The FASTA contains near-duplicate headers that only differ in the part ",
           "clean_id() strips ('|'/space-delimited description text, or an Ensembl version ",
           "suffix); this needs to be resolved in the source FASTA.")
    }
    Biostrings::writeXStringSet(seqs, filepath = fasta_file_clean)
    rm(seqs); gc()

    # --- GTF: stream line-by-line and rewrite only the transcript_id "..."
    # attribute (every exon/CDS/transcript row that carries one), leaving
    # coordinates, gene_id, and every other field untouched. Avoids pulling
    # in rtracklayer/GenomicFeatures just to rename one attribute on a file
    # that can run into the tens of millions of lines. ---
    .rewrite_gtf_transcript_ids <- function(path_in, path_out, clean_fn) {
      con_in  <- file(path_in,  open = "rt")
      on.exit(close(con_in), add = TRUE)
      con_out <- file(path_out, open = "wt")
      on.exit(close(con_out), add = TRUE)

      id_re      <- 'transcript_id "([^"]+)"'
      # Real Ensembl GTFs store gene_id/transcript_id/exon_id BARE and keep
      # the version in a SEPARATE "*_version" attribute, e.g.:
      #   gene_id "ENSG00000223972"; gene_version "5"; transcript_id
      #   "ENST00000456328"; transcript_version "2"; ...
      # (confirmed against actual Ensembl GTF output -- this is not the same
      # convention as the FASTA headers or kallisto/salmon target_ids, which
      # DO carry the version directly, e.g. "ENST00000456328.2".)
      # IsoformSwitchAnalyzeR's GTF import reconstructs the versioned
      # "transcript_id.transcript_version" form for matching -- the sensible
      # default for an unmodified Ensembl GTF paired with a versioned
      # Ensembl FASTA/quantification. But it means cleaning transcript_id
      # alone is a no-op here (it's already bare) while transcript_version
      # is left untouched, so a versioned ID gets reconstructed right back
      # on import -- one that our version-stripped count matrix can never
      # match (0 overlap, regardless of how clean transcript_id itself is).
      # Removing the *_version attributes entirely removes anything left
      # for that reconstruction to draw on.
      version_re  <- '\\s*(gene|transcript|exon)_version "[^"]*";'
      any_id_seen <- FALSE
      tx_ids_seen <- character(0)

      repeat {
        lines <- readLines(con_in, n = 200000L, warn = FALSE)
        if (length(lines) == 0L) break

        lines <- gsub(version_re, "", lines, perl = TRUE)

        m      <- regexpr(id_re, lines)
        has_id <- m != -1L
        if (any(has_id)) {
          any_id_seen <- TRUE
          raw_ids <- sub(id_re, "\\1", regmatches(lines, m))
          new_ids <- clean_fn(raw_ids)
          regmatches(lines, m) <- paste0('transcript_id "', new_ids, '"')

          tx_row <- has_id & grepl("\ttranscript\t", lines, fixed = TRUE)
          if (any(tx_row)) {
            m2 <- regexpr(id_re, lines[tx_row])
            tx_ids_seen <- c(tx_ids_seen, sub(id_re, "\\1", regmatches(lines[tx_row], m2)))
          }
        }
        writeLines(lines, con_out)
      }

      if (!any_id_seen) {
        stop("No `transcript_id \"...\"` attributes found in ", path_in,
             " -- expected standard GTF2-style attributes (if this is GFF3, convert to GTF first).")
      }
      dup <- unique(tx_ids_seen[duplicated(tx_ids_seen)])
      if (length(dup) > 0) {
        stop("Cleaning transcript IDs in isoform_gff produced ", length(dup),
             " duplicate transcript ID(s) among `transcript` feature rows -- e.g. ",
             paste(utils::head(dup, 5), collapse = ", "),
             ". The GTF contains near-duplicate transcript records that only differ in the ",
             "part clean_id() strips; this needs to be resolved in the source GTF before import.")
      }
      invisible(path_out)
    }
    .rewrite_gtf_transcript_ids(gff_file, gff_file_clean, clean_id)

    # estimateDifferentialGeneRange and detectAndCorrectUnwantedEffects both
    # trigger an internal SVA (surrogate-variable) fit inside importRdata()
    # itself (see IsoformSwitchAnalyzeR's NEWS: "importRdata() now
    # automatically corrects abundance and isoform fractions for unwanted
    # covariates...found via sva"). That internal fit has repeated,
    # longstanding edge-case bugs upstream where a discovered surrogate
    # variable's numeric value ends up baked into a design-matrix COLUMN NAME
    # (e.g. "sv1-0.127953077968942"), which is not a syntactically valid R
    # name and crashes limma::makeContrasts() deep inside importRdata()'s
    # internal estimateDifferentialRange() step.
    # ExpressOM's own design_matrix above is always just sampleID+condition
    # (no batch/cofactor column is exposed to run_isoform_switch() yet), so
    # detectAndCorrectUnwantedEffects would only ever be "correcting" for a
    # covariate discovered by that fragile automatic path -- never a
    # user-supplied one -- and estimateDifferentialGeneRange is, by
    # IsoformSwitchAnalyzeR's own documentation, just a rough pilot estimate
    # that "cannot be trusted" and that ExpressOM never reads back out of
    # switch_list anyway. Both are switched off: this sidesteps the crash
    # entirely and skips two computations whose result was never used, which
    # also speeds up import on larger datasets.
    # detectAndCorrectUnwantedEffects is only added if the installed
    # IsoformSwitchAnalyzeR actually exposes it (older versions don't have
    # the argument at all), so this stays compatible either way.
    importrdata_extra_args <- list(estimateDifferentialGeneRange = FALSE)
    if ("detectAndCorrectUnwantedEffects" %in% names(formals(IsoformSwitchAnalyzeR::importRdata))) {
      importrdata_extra_args$detectAndCorrectUnwantedEffects <- FALSE
    }

    # See .muffle_across_deprecation() in utils_core.R: importRdata() runs
    # its own internal dplyr::filter(across()) calls across the whole
    # genome-wide isoform/exon annotation, which on a large reference can
    # raise hundreds of thousands of dplyr-lifecycle deprecation warnings
    # that have nothing to do with ExpressOM's own code.
    switch_list <- .muffle_across_deprecation(do.call(IsoformSwitchAnalyzeR::importRdata, c(
      list(
        isoformCountMatrix   = isoform_count_matrix,
        isoformRepExpression = isoform_count_matrix,
        designMatrix         = design_matrix,
        isoformExonAnnoation = gff_file_clean,
        isoformNtFasta       = fasta_file_clean,
        ignoreAfterBar       = FALSE,
        ignoreAfterSpace     = FALSE,
        ignoreAfterPeriod    = FALSE,
        showProgress         = TRUE
      ),
      importrdata_extra_args
    )))
    # ========================================================================
    message("Isoform data import completed.")

    # Save step-1 checkpoint immediately
    .ckpt_save(switch_list, "step1_imported.rds")

    # Free the large count matrix now that switch_list is built
    if (!is.null(isoform_obj$txi$counts)) {
      isoform_obj$txi$counts <- NULL
      gc()
    }
    if (!is.null(isoform_obj$counts)) {
      isoform_obj$counts <- NULL
      gc()
    }
  }

  # --------------------------------------------------------------------------
  # Step 2 – Run DTU test according to chosen engine
  # --------------------------------------------------------------------------
  if (!already_analyzed) {
    message("Running DTU test using engine: ", test_engine)

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

    # Ensure dplyr is attached
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("dplyr is required but not installed.")
    }
    if (!"package:dplyr" %in% search()) {
      attachNamespace("dplyr")
    }

    # outputPlots = FALSE on every branch below: isoformSwitchAnalysisCombined()
    # / isoformSwitchAnalysisPart2() default to outputPlots = TRUE, which (with
    # pathToOutput = out_dir, as set here) silently writes ISA's own switch
    # plots directly under out_dir using its own comparison-named subfolder
    # (e.g. "IsoformSwitch/CONTROL_vs_NONMUT") -- a third, uncontrolled copy,
    # generated before predictor (CPAT/SignalP/Pfam) annotations even exist,
    # in addition to the two intentional copies this pipeline already
    # produces later: Step 3.5's predictor-refreshed set (out_dir/plots/
    # switch_plots_with_predictors/) and the report's own top_switches/
    # folder. Disabling it here removes that redundant/uncontrolled third
    # copy without losing anything -- both later copies supersede it.
    switch_list <- switch(
      test_engine,
      DEXSeq = {
        message("Using built-in DEXSeq via isoformSwitchAnalysisCombined...")
        tryCatch({
          sl <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisCombined(
            switchAnalyzeRlist = switch_list,
            genomeObject       = genome_object,
            pathToOutput       = out_dir,
            n                  = 50,
            outputPlots        = FALSE
          )
          message("Combined analysis (DEXSeq) completed.")
          sl
        }, error = function(e) {
          no_switches <- grepl("no genes were considered switching", conditionMessage(e), ignore.case = TRUE)
          if (no_switches) {
            message("\n=== No isoform switches detected ===")
            message("isoformSwitchAnalysisCombined() found zero genes meeting the switching ",
                    "cutoffs for '", level, "' vs '", base, "'. This usually means: the effect ",
                    "is genuinely weak/absent in this comparison, sample sizes are too small to ",
                    "reach significance, or the default cutoffs (alpha = 0.05, dIFcutoff = 0.1) ",
                    "are too strict for this dataset. Isoform switch analysis (predictors, ",
                    "consequence/splicing plots, switch plots) will be skipped for this ",
                    "comparison, but DTE and DTU results computed earlier are unaffected and ",
                    "have already been saved.")
          } else {
            message("\n=== isoformSwitchAnalysisCombined() failed ===")
            message("Error: ", conditionMessage(e))
            message("Isoform switch analysis will be skipped for this comparison; DTE and DTU ",
                    "results computed earlier are unaffected and have already been saved.")
          }
          NULL
        })
      },
      DRIMSeq = {
        message("Using pre‑computed DRIMSeq DTU results (from run_dtu)...")
        if (is.null(dtu_results)) {
          stop("test_engine = 'DRIMSeq' requires pre‑computed dtu_results. ",
               "Set run_dtu = TRUE in the main pipeline or pass dtu_results.")
        }
        # Ensure dtu_results is a data.frame with required columns
        drim_df <- dtu_results$dtu_results
        if (is.null(drim_df)) {
          stop("dtu_results does not contain a 'dtu_results' element.")
        }
        required_cols <- c("gene_id", "feature_id", "pvalue", "adj_pvalue")
        if (!all(required_cols %in% colnames(drim_df))) {
          stop("DRIMSeq results must contain columns: ", paste(required_cols, collapse = ", "))
        }
        # Rename 'feature_id' to 'feature_id' (already) and ensure 'adj_pvalue' exists
        # The function accepts dtuResults with columns: gene_id, feature_id, pvalue, adj_pvalue (and optionally log2FC)
        tryCatch({
          sl <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisCombined(
            switchAnalyzeRlist = switch_list,
            genomeObject       = genome_object,
            pathToOutput       = out_dir,
            n                  = 50,
            dtuResults         = drim_df,
            outputPlots        = FALSE
          )
          message("Combined analysis (DRIMSeq) completed.")
          sl
        }, error = function(e) {
          message("Error in isoformSwitchAnalysisCombined with DRIMSeq results: ", e$message)
          NULL
        })
      },
      satuRn = {
        message("Running satuRn DTU test...")
        if (!requireNamespace("satuRn", quietly = TRUE)) {
          stop("Package 'satuRn' is required for test_engine = 'satuRn'. Install with BiocManager::install('satuRn').")
        }
        # First, we need to run the satuRn test on the switchList
        # This function returns a switchList with the DTU results added
        sl_tested <- tryCatch({
          IsoformSwitchAnalyzeR::isoformSwitchTestSatuRn(
            switchAnalyzeRlist = switch_list,
            alpha = 0.05,               # could be made parameter
            dIFcutoff = 0.1,
            reduceToSwitchingGenes = FALSE
          )
        }, error = function(e) {
          message("isoformSwitchTestSatuRn failed: ", e$message)
          NULL
        })
        if (is.null(sl_tested)) {
          message("satuRn test did not return a switchList. Skipping further analysis.")
          NULL
        } else {
          # Now run Part2 (ORF prediction, consequence analysis, switch plots)
          message("Running isoformSwitchAnalysisPart2 (ORF, consequences, plots)...")
          tryCatch({
            sl_final <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisPart2(
              switchAnalyzeRlist = sl_tested,
              genomeObject       = genome_object,
              pathToOutput       = out_dir,
              n                  = 50,
              outputPlots        = FALSE
            )
            message("satuRn analysis completed.")
            sl_final
          }, error = function(e) {
            message("isoformSwitchAnalysisPart2 failed: ", e$message)
            # Return the tested list anyway, which already has the DTU results
            sl_tested
          })
        }
      }
    )

    if (is.null(switch_list)) {
      return(invisible(NULL))
    }

    switch_list <- tryCatch({
      message("Annotating alternative splicing events (exon skipping, intron retention, ",
              "alt. 5'/3' splice sites, alt. TSS/TES)...")
      sl <- .muffle_across_deprecation(IsoformSwitchAnalyzeR::analyzeAlternativeSplicing(
        switch_list,
        onlySwitchingGenes = TRUE,
        alpha     = 0.05,
        dIFcutoff = 0.1,
        quiet     = FALSE
      ))
      message("  Alternative splicing annotation completed.")
      sl
    }, error = function(e) {
      message("  Could not annotate alternative splicing events: ", e$message,
              "\n  Splicing-type summary/enrichment plots will be skipped; ",
              "all other results are unaffected.")
      switch_list
    })

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
        save_dir    = save_dir,
        n_cpu       = predictor_cpu,
        log_dir     = log_dir
      )
      .ckpt_save(switch_list, "step3_predictors.rds")
    }
  }

  # --------------------------------------------------------------------------
  # Step 3.5 - Refresh switch consequence analysis & switch plots now that
  # external predictor annotations (CPAT / SignalP / Pfam) are available.
  # --------------------------------------------------------------------------
  if (isTRUE(run_predictors)) {
    if (.ckpt_exists("step3_5_refreshed.rds")) {
      message("Loading refreshed switch consequences from step-3.5 checkpoint.")
      switch_list <- .ckpt_load("step3_5_refreshed.rds")
    } else {
      message("Refreshing switch consequence analysis and plots with predictor annotations...")
      plot_refresh_dir <- file.path(out_dir, "plots", "switch_plots_with_predictors")

      switch_list <- tryCatch({
        feat_cols <- colnames(switch_list$isoformFeatures)
        consequences_available <- c("intron_retention", "ORF_seq_similarity", "NMD_status")
        if (any(grepl("coding_potential|coding_prob", feat_cols, ignore.case = TRUE)))
          consequences_available <- c(consequences_available, "coding_potential")
        if (any(grepl("signal_peptide", feat_cols, ignore.case = TRUE)))
          consequences_available <- c(consequences_available, "signal_peptide_identified")
        if (any(grepl("^domain", feat_cols, ignore.case = TRUE)))
          consequences_available <- c(consequences_available, "domains_identified", "domain_isotype")

        sl <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
          switch_list,
          consequencesToAnalyze = consequences_available,
          alpha     = 0.05,
          dIFcutoff = 0.1
        )

        if (!dir.exists(plot_refresh_dir)) dir.create(plot_refresh_dir, recursive = TRUE)
        IsoformSwitchAnalyzeR::switchPlotTopSwitches(
          switchAnalyzeRlist          = sl,
          n                           = 50,
          filterForConsequences       = TRUE,
          splitFunctionalConsequences = TRUE,
          sortByQvals                 = TRUE,
          pathToOutput                = plot_refresh_dir,
          fileType                    = "pdf"
        )
        message("  Refreshed switch plots (now including predictor annotations) saved to: ",
                plot_refresh_dir)
        sl
      }, error = function(e) {
        message("  Could not refresh switch consequences/plots with predictor data: ", e$message,
                "\n  Falling back to the plots generated in Step 2 (without predictor annotations).")
        switch_list
      })

      .ckpt_save(switch_list, "step3_5_refreshed.rds")
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
                                      use_wsl, wsl_distro, isoform_obj, save_dir,
                                      n_cpu = NULL, log_dir = NULL) {
  is_windows <- .Platform$OS.type == "windows"
  via_wsl    <- is_windows && isTRUE(use_wsl)

  if (is.null(log_dir)) log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # --------------------------------------------------------------------------
  # Persistent, structured status log for this predictor run.
  #
  # .wsl_exec_script() already logs every actual SHELL command's exit code
  # and stdout/stderr to wsl_commands.log -- but the R-level glue *around*
  # those calls (sequence extraction, and each analyze*() import of a tool's
  # result file back into switch_list) previously had no persistent record
  # at all, only message()s that scroll away in a long/automated run. That
  # gap is exactly where "protein prediction didn't happen, and there's
  # nothing to check" reports come from: a tool can run fine as a shell
  # command (exit 0, result file present) while its R-level import into
  # switch_list silently fails (e.g. an ID mismatch between switch_list and
  # the result file) or never gets attempted. Every step below -- success or
  # failure -- is appended here so that's diagnosable after the fact instead
  # of only in real-time console output.
  # --------------------------------------------------------------------------
  predictor_log_path <- file.path(log_dir, "predictor_status.log")
  .log_predictor_status <- function(step, status, detail = "") {
    cat(sprintf("[%s] %-28s %-8s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step, status, detail),
        file = predictor_log_path, append = TRUE)
    invisible(NULL)
  }
  cat(sprintf("\n===== Predictor run started %s =====\n", format(Sys.time())),
      file = predictor_log_path, append = TRUE)
  message("Predictor status log: ", predictor_log_path)

  # --------------------------------------------------------------------------
  # 0. Decide how many CPU threads hmmscan / InterProScan may use. Auto-detect
  #    when not supplied, leaving one core free for the R session itself.
  # --------------------------------------------------------------------------
  if (is.null(n_cpu) || !is.finite(n_cpu) || n_cpu < 1) {
    detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
    n_cpu <- if (is.na(detected) || detected < 1) 1L else max(1L, detected - 1L)
  }
  n_cpu <- as.integer(n_cpu)
  message("Using ", n_cpu, " CPU thread(s) for hmmscan / InterProScan.")

  # --------------------------------------------------------------------------
  # 1. Discover conda.sh and verify 'isoform_tools' environment
  #    (cached per distro/execution-mode within the R session to avoid
  #    repeating a filesystem-wide `find` on every predictor run)
  # --------------------------------------------------------------------------
  message("Detecting conda environment in execution context...")
  cache_key <- paste0("conda_sh::", wsl_distro, "::", via_wsl)

  if (exists(cache_key, envir = .expressom_cache, inherits = FALSE)) {
    cached        <- get(cache_key, envir = .expressom_cache, inherits = FALSE)
    conda_sh      <- if (nzchar(cached)) cached else NULL
    has_conda_env <- nzchar(cached) > 0
    message("  Using cached conda lookup: ",
            if (has_conda_env) conda_sh else "no 'isoform_tools' env found")
  } else {
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

    assign(cache_key, if (has_conda_env) conda_sh else "", envir = .expressom_cache)

    if (has_conda_env)
      message("  Using conda env 'isoform_tools' (", conda_sh, ")")
    else
      message("  No 'isoform_tools' env found. Using system PATH.")
  }

  active_conda <- if (has_conda_env) conda_sh else NULL

  # --------------------------------------------------------------------------
  # 2. Convenience wrappers – now with enhanced logging
  # --------------------------------------------------------------------------

  # Run a single bash command string (conda activation is handled automatically)
  run_tool <- function(cmd_str, show_stderr = TRUE) {
    message("  Executing: ", cmd_str)
    res <- .wsl_exec_script(
      bash_body     = cmd_str,
      wsl_distro    = wsl_distro,
      use_wsl       = via_wsl,
      conda_sh      = active_conda,
      conda_env     = "isoform_tools",
      intern        = TRUE,              # capture stdout
      ignore_stderr = !show_stderr,
      log_dir       = log_dir
    )
    # Print stdout/stderr if any (res contains both if intern=TRUE)
    if (length(res) > 0) {
      msg <- paste(res, collapse = "\n")
      if (nzchar(msg)) message("  Output:\n", msg)
    }
    # Return exit code (the status attribute is stored on the output)
    attr(res, "status") %||% 0L
  }

  # Check tool availability (conda-aware)
  tool_ok <- function(t) {
    .wsl_tool_exists(t, wsl_distro, via_wsl, active_conda, "isoform_tools", log_dir = log_dir)
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
      .log_predictor_status("extractSequence", "FAILED", e$message)
    })
  }
  # nt/aa sequence counts -- a non-zero NT count with a zero/missing AA count
  # is the single clearest sign that ORF analysis never populated
  # switch_list$orfAnalysis (e.g. isoformSwitchAnalysisPart2() failed and
  # Step 2 fell back to sl_tested further up this file): CPAT can still run
  # off NT sequences alone, but SignalP/Pfam need AA and will be skipped.
  n_nt <- if (file.exists(nt_fa_local)) length(grep("^>", readLines(nt_fa_local, warn = FALSE))) else 0L
  n_aa <- if (file.exists(aa_fa_local)) length(grep("^>", readLines(aa_fa_local, warn = FALSE))) else 0L
  .log_predictor_status("extractSequence", if (n_aa > 0) "OK" else "PARTIAL",
                        sprintf("%d NT sequences, %d AA sequences (AA=0 means SignalP/Pfam will be skipped below)", n_nt, n_aa))

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
  logit_local   <- .find_cpat_logit_model(organism, wsl_distro, via_wsl, active_conda, log_dir = log_dir)

  if (!nzchar(hexamer_local) || !file.exists(hexamer_local)) {
    message("  Hexamer table not found in ISA package. Skipping CPAT.")
    .log_predictor_status("CPAT", "SKIPPED", "hexamer table not found in ISA package")
  } else if (is.null(logit_local)) {
    message("  Logit model not found. Run install_isoform_databases() first. Skipping CPAT.")
    .log_predictor_status("CPAT", "SKIPPED", "logit model not found (run install_isoform_databases())")
  } else {
    hexamer_w      <- w2l(hexamer_local)
    logit_w        <- logit_local  # already a valid path in execution env
    cpat_prefix_w  <- paste0(out_dir_w, "/cpat_out")
    cpat_prefix_l  <- file.path(out_dir, "cpat_out")

    has_cpat3 <- tool_ok("cpat")
    has_cpat2 <- tool_ok("run_cpat.py")

    if (!has_cpat3 && !has_cpat2) {
      message("  Neither 'cpat' nor 'run_cpat.py' found. Skipping.")
      .log_predictor_status("CPAT", "SKIPPED", "neither 'cpat' nor 'run_cpat.py' found on PATH/conda env")
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
        cpat_import_ok <- FALSE
        switch_list <- tryCatch({
          sl <- IsoformSwitchAnalyzeR::analyzeCPAT(
            switchAnalyzeRlist       = switch_list,
            pathToAllCPATresultFiles = cpat_result,
            codingCutoff             = NULL,
            removeNoncodinORFs       = TRUE,
            quiet                    = FALSE
          )
          cpat_import_ok <- TRUE
          sl
        }, error = function(e) {
          message("  analyzeCPAT(): ", e$message)
          switch_list
        })
        # This message previously printed unconditionally right after the
        # tryCatch above, even when analyzeCPAT() had just failed and been
        # caught by the error handler -- so a failed import (switch_list
        # left completely unchanged) still reported as a success right below
        # its own error message, with nothing to indicate they contradicted
        # each other. cpat_import_ok now tracks what actually happened.
        if (cpat_import_ok) {
          message("  CPAT results imported successfully.")
          .log_predictor_status("CPAT", "OK", paste("result file:", cpat_result))
        } else {
          message("  CPAT import failed; coding-potential annotation NOT added to switch_list.")
          .log_predictor_status("CPAT", "IMPORT_FAILED",
                                "shell command succeeded and result file exists, but analyzeCPAT() errored -- see message above")
        }
      } else {
        message("  CPAT failed or result file not found (", basename(cpat_result), "). Skipping.")
        .log_predictor_status("CPAT", "FAILED",
                              sprintf("exit code %s, result file exists: %s", cpat_status, file.exists(cpat_result)))
      }
    }
  }

  # --------------------------------------------------------------------------
  # 6. SignalP – signal peptide prediction (requires AA sequences)
  # --------------------------------------------------------------------------
  message("\n--- SignalP (signal peptide prediction) ---")

  if (is.null(sp_fa_w)) {
    message("  No AA FASTA available (ORF analysis may not have run). Skipping SignalP.")
    .log_predictor_status("SignalP", "SKIPPED", "no AA FASTA extracted (see extractSequence status above)")
  } else {
    sp_out_dir_l <- file.path(out_dir, "signalp_out")
    dir.create(sp_out_dir_l, recursive = TRUE, showWarnings = FALSE)
    sp_out_dir_w <- w2l(sp_out_dir_l)

    has_sp6 <- tool_ok("signalp6")
    has_sp5 <- tool_ok("signalp")

    sp_result_file <- NULL
    sp_attempted   <- FALSE
    sp_status      <- NA_integer_

    if (has_sp6) {
      sp_attempted <- TRUE
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
      sp_attempted <- TRUE

      sp_prefix_w <- paste0(sp_out_dir_w, "/signalp5")
      sp_cmd <- sprintf(
        "signalp -fasta %s -org euk -format short -prefix %s",
        shQuote(sp_fa_w, type = "sh"), shQuote(sp_prefix_w, type = "sh")
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
      .log_predictor_status("SignalP", "SKIPPED", "neither signalp6 nor signalp found on PATH/conda env")
    }

    if (!is.null(sp_result_file)) {
      sp_import_ok <- FALSE
      switch_list <- tryCatch({
        sl <- IsoformSwitchAnalyzeR::analyzeSignalIP(
          switchAnalyzeRlist      = switch_list,
          pathToSignalPresultFile = sp_result_file,
          quiet                   = FALSE
        )
        sp_import_ok <- TRUE
        sl
      }, error = function(e) {
        message("  analyzeSignalIP(): ", e$message)
        switch_list
      })
      # See the identical fix on the CPAT import above: this message used to
      # print unconditionally, even right after analyzeSignalIP() had just
      # failed and switch_list was left unchanged.
      if (sp_import_ok) {
        message("  SignalP results imported successfully.")
        .log_predictor_status("SignalP", "OK", paste("result file:", sp_result_file))
      } else {
        message("  SignalP import failed; signal-peptide annotation NOT added to switch_list.")
        .log_predictor_status("SignalP", "IMPORT_FAILED",
                              "shell command succeeded and result file exists, but analyzeSignalIP() errored -- see message above")
      }
    } else if (sp_attempted && !is.na(sp_status) && sp_status != 0L) {
      message("  SignalP execution failed (exit ", sp_status, "). Skipping.")
      .log_predictor_status("SignalP", "FAILED", paste("exit code", sp_status))
    } else if (sp_attempted) {
      message("  SignalP exited successfully but no expected result file was found. Skipping.")
      .log_predictor_status("SignalP", "FAILED", "exit code 0 but no expected output file found -- check output format/version")
    }
  }

  # --------------------------------------------------------------------------
  # 7. Pfam domain annotation (InterProScan -> hmmscan fallback)
  #    Both require AA / protein sequences.
  # --------------------------------------------------------------------------
  message("\n--- Pfam domain annotation ---")

  if (is.null(pfam_fa_w)) {
    message("  No AA FASTA available. Skipping Pfam.")
    .log_predictor_status("Pfam", "SKIPPED", "no AA FASTA extracted (see extractSequence status above)")
  } else {
    pfam_done <- FALSE

    # 7a. InterProScan (preferred – produces richer XML output)
    if (tool_ok("interproscan.sh")) {
      iprscan_xml_l <- file.path(out_dir, "interproscan.xml")
      iprscan_xml_w <- w2l(iprscan_xml_l)

      ips_cmd <- sprintf(
        "interproscan.sh -i %s -f XML -o %s -dp -appl Pfam -goterms -iprlookup -cpu %d",
        shQuote(pfam_fa_w,     type = "sh"),
        shQuote(iprscan_xml_w, type = "sh"),
        n_cpu
      )
      message("  Running InterProScan...")
      ips_status <- run_tool(ips_cmd)

      if (ips_status == 0L && file.exists(iprscan_xml_l)) {
        ips_import_ok <- FALSE
        switch_list <- tryCatch({
          sl <- IsoformSwitchAnalyzeR::analyzeInterProScan(
            switchAnalyzeRlist           = switch_list,
            pathToInterProScanResultFile = iprscan_xml_l,
            quiet                        = FALSE
          )
          ips_import_ok <- TRUE
          sl
        }, error = function(e) {
          message("  analyzeInterProScan(): ", e$message)
          switch_list
        })
        # pfam_done previously was set TRUE (and "imported successfully"
        # printed) unconditionally right after this tryCatch, even when
        # analyzeInterProScan() had just failed and been caught by the error
        # handler above -- switch_list left completely unchanged. Since
        # pfam_done gates the hmmscan fallback below (`if (!pfam_done)`), that
        # meant a failed InterProScan *import* (as opposed to a failed
        # InterProScan *run*, already handled by the ips_status/file.exists
        # check) silently skipped hmmscan too, leaving switch_list with NO
        # Pfam domain annotation at all while every message claimed success.
        pfam_done <- ips_import_ok
        if (ips_import_ok) {
          message("  InterProScan Pfam results imported successfully.")
          .log_predictor_status("Pfam-InterProScan", "OK", paste("result file:", iprscan_xml_l))
        } else {
          message("  InterProScan import failed; trying hmmscan fallback...")
          .log_predictor_status("Pfam-InterProScan", "IMPORT_FAILED",
                                "shell command succeeded and result file exists, but analyzeInterProScan() errored -- see message above; falling through to hmmscan")
        }
      } else {
        message("  InterProScan failed (exit ", ips_status, "). Trying hmmscan fallback...")
        .log_predictor_status("Pfam-InterProScan", "FAILED",
                              sprintf("exit code %s, result file exists: %s", ips_status, file.exists(iprscan_xml_l)))
      }
    } else {
      .log_predictor_status("Pfam-InterProScan", "SKIPPED", "interproscan.sh not found on PATH/conda env")
    }

    # 7b. hmmscan fallback (requires Pfam-A.hmm indexed with hmmpress)
    if (!pfam_done) {
      if (!tool_ok("hmmscan")) {
        message("  hmmscan not found. Skipping Pfam analysis.")
        .log_predictor_status("Pfam-hmmscan", "SKIPPED", "hmmscan not found on PATH/conda env")
      } else {
        pfam_db_w <- .find_pfam_db(wsl_distro, via_wsl, active_conda, log_dir = log_dir)

        if (is.null(pfam_db_w)) {
          message("  Pfam-A.hmm not found. Run install_isoform_databases() first. Skipping.")
          .log_predictor_status("Pfam-hmmscan", "SKIPPED", "Pfam-A.hmm not found (run install_isoform_databases())")
        } else {
          pfam_tbl_l <- file.path(out_dir, "pfam_domtblout.txt")
          pfam_tbl_w <- w2l(pfam_tbl_l)

          # Use --domtblout (domain table) – the format ISA's analyzePFAM expects
          hm_cmd <- sprintf(
            "hmmscan --cpu %d --domtblout %s %s %s",
            n_cpu,
            shQuote(pfam_tbl_w, type = "sh"),
            shQuote(pfam_db_w,  type = "sh"),
            shQuote(pfam_fa_w,  type = "sh")
          )
          message("  Running hmmscan...")
          hm_status <- run_tool(hm_cmd)

          if (hm_status == 0L && file.exists(pfam_tbl_l)) {
            hm_import_ok <- FALSE
            switch_list <- tryCatch({
              sl <- IsoformSwitchAnalyzeR::analyzePFAM(
                switchAnalyzeRlist   = switch_list,
                pathToPFAMresultFile = pfam_tbl_l,
                quiet                = FALSE
              )
              hm_import_ok <- TRUE
              sl
            }, error = function(e) {
              message("  analyzePFAM(): ", e$message)
              switch_list
            })
            # Same false-success fix as above: only report/log success when
            # analyzePFAM() actually succeeded.
            if (hm_import_ok) {
              message("  hmmscan Pfam results imported successfully.")
              .log_predictor_status("Pfam-hmmscan", "OK", paste("result file:", pfam_tbl_l))
            } else {
              message("  hmmscan import failed; Pfam domain annotation NOT added to switch_list.")
              .log_predictor_status("Pfam-hmmscan", "IMPORT_FAILED",
                                    "shell command succeeded and result file exists, but analyzePFAM() errored -- see message above")
            }
          } else {
            message("  hmmscan failed (exit ", hm_status, "). Skipping Pfam.")
            .log_predictor_status("Pfam-hmmscan", "FAILED",
                                  sprintf("exit code %s, result file exists: %s", hm_status, file.exists(pfam_tbl_l)))
          }
        }
      }
    }
  }

  # --------------------------------------------------------------------------
  # 8. Final summary: what actually ended up in switch_list, independent of
  #    what each step above *attempted*. This is the single place that
  #    answers "did protein prediction actually happen" -- both to the
  #    console and to predictor_status.log -- rather than needing to
  #    reconstruct it from scattered per-tool messages.
  # --------------------------------------------------------------------------
  n_orf    <- if (!is.null(switch_list$orfAnalysis))            nrow(switch_list$orfAnalysis) else 0L
  n_domain <- if (!is.null(switch_list$domainAnalysis))         length(unique(switch_list$domainAnalysis$isoform_id)) else 0L
  n_sigp   <- if (!is.null(switch_list$signalPeptideAnalysis))  length(unique(switch_list$signalPeptideAnalysis$isoform_id)) else 0L
  n_cds    <- if (!is.null(switch_list$isoformFeatures) && "codingPotential" %in% colnames(switch_list$isoformFeatures))
                sum(!is.na(switch_list$isoformFeatures$codingPotential)) else 0L
  summary_msg <- sprintf(
    "ORF analysis: %d transcripts with an ORF entry | Coding potential (CPAT): %d annotated | Domains (Pfam): %d transcripts annotated | Signal peptide: %d transcripts annotated",
    n_orf, n_cds, n_domain, n_sigp
  )
  message("\n--- Predictor summary ---\n  ", summary_msg)
  .log_predictor_status("SUMMARY", if (n_domain > 0 || n_sigp > 0 || n_cds > 0) "OK" else "EMPTY", summary_msg)

  return(switch_list)
}

# ==============================================================================
# Enhanced isoform-level visualization: switch summaries, sashimi-style
# junction plots, and exon-bin usage comparison
# ==============================================================================
#
# These functions build on the switchAnalyzeRlist produced by
# run_isoform_switch() to give richer graphic interpretation of isoform-level
# results:
#   - plot_isoform_switch_summary(): wraps IsoformSwitchAnalyzeR's own
#     switchPlot()/switchPlotTopSwitches() to render, per gene, the isoform
#     structure (exons/domains/ORF) alongside gene expression, isoform
#     expression, and isoform usage.
#   - plot_isoform_sashimi(): a custom, mirrored sashimi-style diagram of a
#     gene's exon structure with splice-junction arcs whose thickness reflects
#     isoform-fraction-weighted usage in each condition. NOTE this is a
#     transcript-quantification-based approximation (built from the per-isoform
#     IF1/IF2 values IsoformSwitchAnalyzeR already computes) and NOT derived
#     from spliced-read alignments/junction coverage -- true base-pair
#     resolution coverage requires BAM files, which are outside the scope of
#     this quantification-only (salmon/kallisto/rsem/stringtie) pipeline.
#   - plot_exon_usage_comparison(): collapses overlapping isoform exons into
#     non-overlapping genomic "exon bins" (the same disjoint-interval strategy
#     DEXSeq uses for its flattened exonic-part annotation) and compares the
#     isoform-expression-weighted signal assigned to each bin between conditions.

#' @keywords internal
.resolve_gene_id <- function(isoform_features, gene) {
  if (is.null(gene) || !nzchar(gene)) return(NULL)
  if (gene %in% isoform_features$gene_id) return(gene)
  if ("gene_name" %in% colnames(isoform_features)) {
    idx <- which(toupper(isoform_features$gene_name) == toupper(gene))
    if (length(idx) > 0) return(isoform_features$gene_id[idx[1]])
  }
  NULL
}

#' Build shared, non-overlapping exon bins for a gene's isoforms
#'
#' Uses \code{GenomicRanges::disjoin()} (the same strategy behind DEXSeq's
#' flattened exonic-part annotation) so that alternative splice sites produce
#' distinct adjacent bins rather than being silently merged together.
#' @keywords internal
.get_exon_bins <- function(switch_list, gene_id) {
  if (is.null(switch_list$exons)) {
    stop("switch_list has no 'exons' slot (unexpected switchAnalyzeRlist structure).")
  }
  exon_mcols <- S4Vectors::mcols(switch_list$exons)
  if (!"isoform_id" %in% colnames(exon_mcols)) {
    stop("switch_list$exons has no 'isoform_id' metadata column (unexpected switchAnalyzeRlist structure).")
  }

  feat    <- switch_list$isoformFeatures
  iso_ids <- unique(feat$isoform_id[feat$gene_id == gene_id])
  if (length(iso_ids) == 0) return(NULL)

  exon_gr <- switch_list$exons[exon_mcols$isoform_id %in% iso_ids]
  if (length(exon_gr) == 0) return(NULL)

  bins <- GenomicRanges::disjoin(exon_gr, ignore.strand = TRUE)
  bins <- sort(bins)
  S4Vectors::mcols(bins)$bin_id <- seq_along(bins)

  list(exon_gr = exon_gr, bins = bins, iso_ids = iso_ids)
}

#' Generate Isoform Switch Summary Plots (structure + expression + usage)
#'
#' Wraps \code{IsoformSwitchAnalyzeR::switchPlotTopSwitches()} to automatically
#' render the top-N most significant isoform switches, and additionally forces
#' individual \code{switchPlot()} figures for any explicitly requested genes
#' (even if not among the automatically selected top switches). Each figure is
#' a composite panel showing (1) isoform/transcript structure with any
#' annotated ORF, coding potential, protein domains, or signal peptides,
#' (2) gene expression, (3) isoform expression, and (4) isoform usage (IF),
#' including the switch-test result.
#'
#' @param switch_list Object returned by run_isoform_switch()
#' @param plot_dir Directory to save plots into (created if missing)
#' @param genes_of_interest Optional character vector of gene symbols to
#'   force-include regardless of significance
#' @param level Foreground condition label (used to resolve condition2 for
#'   forced gene plots; if NULL, taken from switch_list)
#' @param base Reference condition label (used to resolve condition1 for
#'   forced gene plots; if NULL, taken from switch_list)
#' @param top_n Number of top switching genes to auto-plot (default 10)
#' @param alpha Significance cutoff on the isoform switch q-value (default 0.05)
#' @param dIFcutoff Minimum |dIF| for a switch to be considered (default 0.1)
#' @return Invisibly, a character vector of generated PDF paths
#' @export
plot_isoform_switch_summary <- function(switch_list, plot_dir,
                                         genes_of_interest = NULL,
                                         level = NULL, base = NULL,
                                         top_n = 10,
                                         alpha = 0.05,
                                         dIFcutoff = 0.1) {
  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    message("IsoformSwitchAnalyzeR is required for switch summary plots. Skipping.")
    return(invisible(character(0)))
  }
  if (is.null(switch_list$isoformFeatures) || nrow(switch_list$isoformFeatures) == 0) {
    message("switch_list has no isoformFeatures - skipping switch summary plots.")
    return(invisible(character(0)))
  }
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  generated <- character(0)
  feat <- switch_list$isoformFeatures

  # ---- Automatic top-N plots (ranked by significance / effect size) --------
  top_dir <- file.path(plot_dir, "top_switches")
  auto_ok <- tryCatch({
    IsoformSwitchAnalyzeR::switchPlotTopSwitches(
      switchAnalyzeRlist          = switch_list,
      n                           = top_n,
      filterForConsequences       = FALSE,
      splitFunctionalConsequences = FALSE,
      sortByQvals                 = TRUE,
      alpha                       = alpha,
      dIFcutoff                   = dIFcutoff,
      pathToOutput                = top_dir,
      fileType                    = "pdf"
    )
    TRUE
  }, error = function(e) {
    message("   -> switchPlotTopSwitches() failed: ", e$message)
    FALSE
  })
  if (isTRUE(auto_ok) && dir.exists(top_dir)) {
    found <- list.files(top_dir, pattern = "\\.pdf$", full.names = TRUE, recursive = TRUE)
    generated <- c(generated, found)
    message("   -> Top ", top_n, " isoform switch plots saved to: ", top_dir,
            " (", length(found), " file(s))")
  }

  # ---- Force-include explicitly requested genes -----------------------------
  if (!is.null(genes_of_interest)) {
    cond1 <- if (!is.null(base))  base  else feat$condition_1[1]
    cond2 <- if (!is.null(level)) level else feat$condition_2[1]

    for (gene_sym in genes_of_interest) {
      gene_id <- .resolve_gene_id(feat, gene_sym)
      if (is.null(gene_id)) {
        message("   -> Gene '", gene_sym, "' not found in switch_list. Skipping switch plot.")
        next
      }
      pdf_path <- file.path(plot_dir, paste0("SwitchPlot_", gene_sym, ".pdf"))
      ok <- tryCatch({
        pdf(pdf_path, onefile = FALSE, width = 11, height = 9)
        IsoformSwitchAnalyzeR::switchPlot(switch_list, gene = gene_id,
                                          condition1 = cond1, condition2 = cond2)
        dev.off()
        TRUE
      }, error = function(e) {
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
        message("   -> switchPlot() failed for ", gene_sym, ": ", e$message)
        FALSE
      })
      if (ok) {
        generated <- c(generated, pdf_path)
        message("   -> Isoform switch plot saved to: ", pdf_path)
      }
    }
  }

  invisible(generated)
}

#' Pseudo-Sashimi Isoform Structure & Junction Usage Plot
#'
#' Draws a mirrored, sashimi-style diagram of a gene's exon structure with
#' splice-junction arcs whose thickness reflects isoform-fraction-weighted
#' usage in each condition: arcs above the exon track represent \code{level}
#' usage, arcs below represent \code{base} usage, so the two conditions can be
#' compared directly at a glance. See the section header above for an
#' important note on how this differs from a read-coverage sashimi plot.
#'
#' @param switch_list Object returned by run_isoform_switch()
#' @param gene Gene symbol or gene_id to plot
#' @param level Foreground condition label (arcs drawn above the exon track)
#' @param base Reference condition label (arcs drawn below the exon track)
#' @param plot_dir Output directory for the PDF
#' @param min_if Minimum isoform-fraction-weighted usage in either condition
#'   for a junction to be drawn (removes near-zero-usage clutter, default 0.01)
#' @return Invisibly, the path to the generated PDF (or NULL if skipped)
#' @export
plot_isoform_sashimi <- function(switch_list, gene, level, base, plot_dir, min_if = 0.01) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    message("GenomicRanges/S4Vectors are required for sashimi-style plots. Skipping.")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting")

  feat_all <- switch_list$isoformFeatures
  gene_id  <- .resolve_gene_id(feat_all, gene)
  if (is.null(gene_id)) {
    message("   -> Gene '", gene, "' not found in switch_list. Skipping sashimi plot.")
    return(invisible(NULL))
  }

  feat <- feat_all[feat_all$gene_id == gene_id &
                    feat_all$condition_1 == base &
                    feat_all$condition_2 == level, ]
  if (nrow(feat) == 0) {
    message("   -> No isoformFeatures rows for gene '", gene, "' with condition_1=", base,
            ", condition_2=", level, ". Skipping sashimi plot.")
    return(invisible(NULL))
  }

  ginfo <- tryCatch(.get_exon_bins(switch_list, gene_id), error = function(e) {
    message("   -> Could not build exon structure for '", gene, "': ", e$message)
    NULL
  })
  if (is.null(ginfo)) {
    message("   -> No exon structure found for '", gene, "'. Skipping sashimi plot.")
    return(invisible(NULL))
  }

  exon_gr    <- ginfo$exon_gr
  exon_mcols <- S4Vectors::mcols(exon_gr)
  exon_df <- data.frame(
    isoform_id = exon_mcols$isoform_id,
    start      = GenomicRanges::start(exon_gr),
    end        = GenomicRanges::end(exon_gr),
    stringsAsFactors = FALSE
  )

  # Build per-isoform junctions from consecutive exons in genomic order
  junction_list <- lapply(split(exon_df, exon_df$isoform_id), function(d) {
    d <- d[order(d$start), ]
    if (nrow(d) < 2) return(NULL)
    data.frame(
      isoform_id = d$isoform_id[1],
      jstart     = d$end[-nrow(d)],
      jend       = d$start[-1],
      stringsAsFactors = FALSE
    )
  })
  junctions <- do.call(rbind, junction_list)
  if (is.null(junctions) || nrow(junctions) == 0) {
    message("   -> Gene '", gene, "' has no multi-exon isoforms (no junctions to plot). Skipping.")
    return(invisible(NULL))
  }

  junctions <- merge(junctions, feat[, c("isoform_id", "IF1", "IF2")], by = "isoform_id")
  junc_summary <- stats::aggregate(cbind(IF1, IF2) ~ jstart + jend, data = junctions,
                                   FUN = function(x) sum(x, na.rm = TRUE))
  colnames(junc_summary)[colnames(junc_summary) == "IF1"] <- "base_weight"
  colnames(junc_summary)[colnames(junc_summary) == "IF2"] <- "level_weight"
  junc_summary <- junc_summary[junc_summary$base_weight >= min_if | junc_summary$level_weight >= min_if, ]
  if (nrow(junc_summary) == 0) {
    message("   -> No junctions passed the min_if = ", min_if, " threshold for '", gene, "'. Skipping.")
    return(invisible(NULL))
  }
  junc_summary$junction_id <- seq_len(nrow(junc_summary))

  .make_arc <- function(jstart, jend, height, junction_id, side, weight, n_points = 40) {
    t <- seq(0, 1, length.out = n_points)
    data.frame(
      x = jstart + t * (jend - jstart),
      y = 4 * height * t * (1 - t),
      junction_id = junction_id,
      side = side,
      weight = weight
    )
  }

  max_w <- max(c(junc_summary$base_weight, junc_summary$level_weight), na.rm = TRUE)
  max_w <- max(max_w, 1e-6)

  arc_rows <- list()
  for (i in seq_len(nrow(junc_summary))) {
    r <- junc_summary[i, ]
    if (r$level_weight >= min_if) {
      arc_rows[[length(arc_rows) + 1]] <- .make_arc(
        r$jstart, r$jend, height = r$level_weight / max_w,
        junction_id = r$junction_id, side = level, weight = r$level_weight
      )
    }
    if (r$base_weight >= min_if) {
      arc_rows[[length(arc_rows) + 1]] <- .make_arc(
        r$jstart, r$jend, height = -(r$base_weight / max_w),
        junction_id = r$junction_id, side = base, weight = r$base_weight
      )
    }
  }
  arc_df <- do.call(rbind, arc_rows)
  arc_df$side <- factor(arc_df$side, levels = c(base, level))

  bins_df <- as.data.frame(ginfo$bins)
  bins_df$bin_id <- S4Vectors::mcols(ginfo$bins)$bin_id
  chrom_label <- as.character(GenomicRanges::seqnames(ginfo$bins)[1])

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bins_df,
                       ggplot2::aes(xmin = start, xmax = end, ymin = -0.04, ymax = 0.04),
                       fill = "grey35", color = NA) +
    ggplot2::geom_path(data = arc_df,
                       ggplot2::aes(x = x, y = y, group = interaction(junction_id, side),
                                    linewidth = weight, color = side),
                       lineend = "round") +
    ggplot2::scale_color_manual(values = stats::setNames(c("dodgerblue4", "red3"), c(base, level)),
                                name = "Condition") +
    ggplot2::scale_linewidth(range = c(0.3, 3), name = "Usage (IF)") +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    ggplot2::labs(
      title = paste0("Isoform Structure & Junction Usage: ", gene),
      subtitle = paste0(level, " (above) vs ", base,
                        " (below) \u2014 arc thickness \u221d isoform-fraction-weighted junction usage"),
      x = paste0("Genomic position (", chrom_label, ")"),
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "grey30")
    )

  pdf_path <- file.path(plot_dir, paste0("Sashimi_", gene, ".pdf"))
  ggplot2::ggsave(filename = pdf_path, plot = p, width = 11, height = 5, device = "pdf")
  message("   -> Isoform sashimi-style junction plot saved to: ", pdf_path)
  invisible(pdf_path)
}

#' Compare Aggregated Exon-Bin Expression Between Conditions
#'
#' For a given gene, collapses all isoform exon structures into non-overlapping
#' genomic "exon bins" (via \code{GenomicRanges::disjoin()}, the same strategy
#' DEXSeq uses to build its flattened exonic-part annotation) and sums the
#' per-isoform mean expression already computed by IsoformSwitchAnalyzeR
#' (stored as \code{iso_value_1}/\code{iso_value_2} in isoformFeatures) across
#' every isoform overlapping each bin, per condition. This approximates
#' per-exon "coverage" from transcript-level quantification without requiring
#' aligned BAM files.
#'
#' @inheritParams plot_isoform_sashimi
#' @return Invisibly, the path to the generated PDF (or NULL if skipped)
#' @export
plot_exon_usage_comparison <- function(switch_list, gene, level, base, plot_dir) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    message("GenomicRanges/S4Vectors are required for exon usage plots. Skipping.")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required for plotting")

  feat_all <- switch_list$isoformFeatures
  gene_id  <- .resolve_gene_id(feat_all, gene)
  if (is.null(gene_id)) {
    message("   -> Gene '", gene, "' not found in switch_list. Skipping exon usage plot.")
    return(invisible(NULL))
  }

  feat <- feat_all[feat_all$gene_id == gene_id &
                    feat_all$condition_1 == base &
                    feat_all$condition_2 == level, ]
  if (nrow(feat) == 0) {
    message("   -> No isoformFeatures rows for gene '", gene, "' with condition_1=", base,
            ", condition_2=", level, ". Skipping exon usage plot.")
    return(invisible(NULL))
  }

  ginfo <- tryCatch(.get_exon_bins(switch_list, gene_id), error = function(e) {
    message("   -> Could not build exon structure for '", gene, "': ", e$message)
    NULL
  })
  if (is.null(ginfo)) {
    message("   -> No exon structure found for '", gene, "'. Skipping exon usage plot.")
    return(invisible(NULL))
  }

  exon_gr <- ginfo$exon_gr
  bins    <- ginfo$bins

  ov <- GenomicRanges::findOverlaps(exon_gr, bins, ignore.strand = TRUE)
  map_df <- data.frame(
    isoform_id = S4Vectors::mcols(exon_gr)$isoform_id[S4Vectors::queryHits(ov)],
    bin_id     = S4Vectors::mcols(bins)$bin_id[S4Vectors::subjectHits(ov)],
    stringsAsFactors = FALSE
  )
  map_df <- unique(map_df)
  map_df <- merge(map_df, feat[, c("isoform_id", "iso_value_1", "iso_value_2")], by = "isoform_id")
  if (nrow(map_df) == 0) {
    message("   -> No overlapping isoform/bin data for '", gene, "'. Skipping exon usage plot.")
    return(invisible(NULL))
  }

  bin_summary <- stats::aggregate(cbind(iso_value_1, iso_value_2) ~ bin_id, data = map_df, FUN = sum)
  colnames(bin_summary) <- c("bin_id", base, level)

  plot_df <- tidyr::pivot_longer(bin_summary, cols = -bin_id, names_to = "Condition", values_to = "expression")
  plot_df$Condition <- factor(plot_df$Condition, levels = c(base, level))
  plot_df$bin_id    <- factor(plot_df$bin_id, levels = sort(unique(plot_df$bin_id)))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = bin_id, y = expression, fill = Condition)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6,
                      color = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("dodgerblue4", "red3"), c(base, level))) +
    ggplot2::labs(
      title = paste0("Exon-Bin Expression Comparison: ", gene),
      subtitle = "Non-overlapping exon bins (DEXSeq-style flattened annotation); bars = summed isoform expression per bin",
      x = "Exon bin (5'\u2192 3', genomic order)", y = "Summed isoform expression"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "grey30"),
      legend.title = ggplot2::element_blank()
    )

  pdf_path <- file.path(plot_dir, paste0("ExonUsage_", gene, ".pdf"))
  ggplot2::ggsave(filename = pdf_path, plot = p,
                  width = max(7, 0.35 * nrow(bin_summary) + 2), height = 5, device = "pdf")
  message("   -> Exon-bin expression comparison plot saved to: ", pdf_path)
  invisible(pdf_path)
}

#' Generate DTE/DTU HTML Report with conditional image includes
#' (Now uses eval = file.exists(...) to avoid "Could not fetch resource" warnings)
generate_dte_dtu_report <- function(dte_results, dtu_results, isoform_obj,
                                    out_dir, condition, level, base,
                                    genes_of_interest = NULL, top_n = 15,
                                    switch_list = NULL,
                                    dexseq_results = NULL,
                                    plot_switch_summary = TRUE,
                                    switch_plot_top_n = 10,
                                    plot_sashimi = TRUE,
                                    plot_exon_usage = TRUE) {
  
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
  # transcript_id/gene here are already clean+coalesced by run_dte() itself;
  # re-deriving with .coalesce_gene_label() rather than trusting dte$gene
  # directly guards the same way .prepare_dtu_counts() does against a stale
  # checkpoint predating that fix (see run_dte()'s comment for the full
  # explanation of why an un-re-verified gene_symbol/transcript_id here was
  # previously showing up as literal "NA (raw|pipe|delimited|header)" in
  # DTE_top_barplot and the DTE volcano).
  dte$transcript_id <- clean_transcript_id(dte$transcript_id)
  dte$gene_label <- .format_gene_tx_label(
    .coalesce_gene_label(dte$gene_symbol, dte$transcript_id), dte$transcript_id
  )
  dte$signif_label <- ifelse(dte$signif, "Significant", "Not significant")
  
  # ---- 2. Annotate DTU results with gene symbols ----
  dtu <- dtu_results$dtu_results
  gene_map <- isoform_obj$gene_map[, c("ensembl", "symbol")]
  colnames(gene_map) <- c("gene_id", "gene_symbol")
  # If run_dtu()/run_dexseq_dtu() already resolved gene_symbol upstream,
  # re-merging here would silently blank it out again for anything gene_map
  # doesn't cover by Ensembl ID -- only fill in if not already present.
  if (!"gene_symbol" %in% colnames(dtu)) {
    dtu <- merge(dtu, gene_map, by = "gene_id", all.x = TRUE)
  }

  missing_sym <- is.na(dtu$gene_symbol) | dtu$gene_symbol == ""
  if (any(missing_sym)) {
    already_a_symbol <- dtu$gene_id[missing_sym] %in% gene_map$gene_symbol
    dtu$gene_symbol[missing_sym][already_a_symbol] <- dtu$gene_id[missing_sym][already_a_symbol]
  }
  # Display-safe label -- see .coalesce_gene_label(): a real gene_symbol NA
  # is expected for a genuinely novel/unannotated gene, but nothing
  # user-facing should ever print the literal text "NA" for it.
  dtu$gene_label <- .format_gene_tx_label(.coalesce_gene_label(dtu$gene_symbol, dtu$gene_id), dtu$gene_id)
  
  # Save annotated tables as CSV
  write.csv(dte, file.path(report_dir, "DTE_results_annotated.csv"), row.names = FALSE)
  write.csv(dtu, file.path(report_dir, "DTU_results_annotated.csv"), row.names = FALSE)

  # 2b. Save DEXSeq results as CSV, if the complementary DEXSeq DTU engine was run
  has_dexseq <- !is.null(dexseq_results) && !is.null(dexseq_results$results_df)
  if (has_dexseq) {
    write.csv(dexseq_results$results_df, file.path(report_dir, "DEXSeq_results_annotated.csv"),
              row.names = FALSE)
  }
  
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
  # Label the most extreme significant points with their gene so individual
  # points are identifiable -- previously this plot only had a
  # significant/not-significant colour legend, with no way to tell which
  # gene any given dot was.
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    dte_to_label <- head(dte_sig[dte_sig$signif, ][order(dte_sig[dte_sig$signif, ]$padj), ], 20)
    if (nrow(dte_to_label) > 0) {
      p_volcano <- p_volcano + ggrepel::geom_text_repel(
        data = dte_to_label,
        mapping = ggplot2::aes(x = log2FoldChange, y = -log10(padj), label = gene_label),
        inherit.aes = FALSE, size = 2.8, max.overlaps = 30,
        segment.size = 0.3, color = "black", seed = 42
      )
    }
  } else {
    message("   -> ggrepel not installed; DTE_volcano.pdf will show points without gene labels.")
  }
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
  
  # 3c. Top significant transcripts barplot (absolute log2FC)
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

  # 3d-bis. Top significant DTU transcripts barplot (by adjusted p-value)
  # Companion to the DTE top-barplot above. DRIMSeq's feature-level test
  # statistic ('lr', a likelihood-ratio stat) is unsigned -- there's no
  # DTE-style log2FoldChange to rank by -- so -log10(adj. p-value) stands in
  # as both the ranking measure and the bar height. Prefers the
  # transcript-level `feature_id` when present (the actual DTU unit being
  # tested); falls back to gene_id if it isn't (e.g. an older cached
  # checkpoint from before run_dtu() started requesting feature-level
  # DRIMSeq results) -- either way the label is built through the same
  # coalesce/format helpers as everywhere else, so a gene with no resolvable
  # symbol shows its (clean) id once instead of a bare/duplicated one.
  dtu_label_col <- if ("feature_id" %in% colnames(dtu)) "feature_id" else "gene_id"
  dtu_sig <- dtu[!is.na(dtu$adj_pvalue) & dtu$adj_pvalue < 0.05, ]
  if (nrow(dtu_sig) > 0) {
    dtu_ordered <- dtu_sig[order(dtu_sig$adj_pvalue), ]
    top_dtu <- head(dtu_ordered, top_n)
    top_dtu$dtu_label <- .format_gene_tx_label(
      .coalesce_gene_label(top_dtu$gene_symbol, top_dtu[[dtu_label_col]]),
      top_dtu[[dtu_label_col]]
    )
    top_dtu$neg_log10_padj <- -log10(top_dtu$adj_pvalue)
    p_dtu_bar <- ggplot2::ggplot(top_dtu, ggplot2::aes(x = reorder(dtu_label, neg_log10_padj), y = neg_log10_padj)) +
      ggplot2::geom_col(fill = "darkgreen") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste("Top", top_n, "significant DTU transcripts"),
           x = "", y = "-log10(adj. p-value)") +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(plot_dir, "DTU_top_barplot.pdf"), p_dtu_bar, width = 10, height = max(6, top_n*0.3))
  } else {
    message("No significant DTU transcripts found (adj_pvalue < 0.05); skipping DTU top barplot.")
  }

  # 3d-ter. Cross-engine significant-gene concordance.
  #
  # DTE, DTU, and (when run) DEXSeq each test a genuinely different
  # hypothesis (overall expression change; relative isoform usage change;
  # exon-bin usage change) and are shown throughout this report as separate,
  # independent sections -- there was previously no single place that showed
  # how much they actually agree on *which genes* are implicated. A gene
  # significant in all three is a much stronger candidate than one only
  # DEXSeq calls; this makes that visible without cross-referencing three
  # separate tables by hand.
  safe_run({
    sig_sets <- list(
      DTE = unique(dte$gene_id[!is.na(dte$gene_id) & dte$signif]),
      DTU = unique(dtu$gene_id[!is.na(dtu$gene_id) & !is.na(dtu$adj_pvalue) & dtu$adj_pvalue < 0.05])
    )
    if (has_dexseq && "gene_qvalue" %in% colnames(dexseq_results$results_df)) {
      dex <- dexseq_results$results_df
      sig_sets$DEXSeq <- unique(dex$gene_id[!is.na(dex$gene_id) & !is.na(dex$gene_qvalue) & dex$gene_qvalue < 0.05])
    }
    sig_sets <- sig_sets[lengths(sig_sets) > 0]

    if (length(sig_sets) >= 2) {
      all_genes  <- unique(unlist(sig_sets, use.names = FALSE))
      membership <- vapply(sig_sets, function(s) all_genes %in% s, logical(length(all_genes)))
      combo <- apply(membership, 1, function(row) paste(names(sig_sets)[row], collapse = " + "))
      combo_counts <- as.data.frame(table(combination = combo), stringsAsFactors = FALSE)
      colnames(combo_counts) <- c("combination", "n_genes")
      combo_counts <- combo_counts[order(-combo_counts$n_genes), ]
      combo_counts$combination <- factor(combo_counts$combination, levels = combo_counts$combination)

      p_concordance <- ggplot2::ggplot(combo_counts, ggplot2::aes(x = combination, y = n_genes)) +
        ggplot2::geom_col(fill = "darkslateblue") +
        ggplot2::geom_text(ggplot2::aes(label = n_genes), vjust = -0.4, size = 3.5) +
        ggplot2::labs(title = "Significant gene concordance across engines",
                      subtitle = paste(sprintf("%s: %d genes", names(sig_sets), lengths(sig_sets)), collapse = "   |   "),
                      x = NULL, y = "Number of genes") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
      ggplot2::ggsave(file.path(plot_dir, "Engine_concordance.pdf"), p_concordance,
                      width = max(7, 1.3 * nrow(combo_counts)), height = 6)
      message("   -> Engine concordance (significant genes): ",
              paste(sprintf("%s=%d", names(sig_sets), lengths(sig_sets)), collapse = ", "))
    } else {
      message("   -> Skipping engine concordance plot: fewer than 2 engines have any significant genes to compare.")
    }
  }, label = "Cross-engine (DTE/DTU/DEXSeq) concordance plot")

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

  # 3f. Isoform switch overview + summary plots (structure, expression, usage)
  has_switch <- !is.null(switch_list) && !is.null(switch_list$isoformFeatures) &&
                nrow(switch_list$isoformFeatures) > 0

  if (has_switch) {
    feat_sw <- switch_list$isoformFeatures
    if (all(c("dIF", "isoform_switch_q_value") %in% colnames(feat_sw))) {
      feat_sw_plot <- feat_sw[!is.na(feat_sw$isoform_switch_q_value), ]
      if (nrow(feat_sw_plot) > 0) {
        p_switch_overview <- ggplot2::ggplot(
          feat_sw_plot, ggplot2::aes(x = dIF, y = -log10(isoform_switch_q_value))
        ) +
          ggplot2::geom_point(
            ggplot2::aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
            size = 1, alpha = 0.6
          ) +
          ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
          ggplot2::geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
          ggplot2::scale_color_manual("Significant\nIsoform Switch", values = c("grey60", "red3")) +
          ggplot2::labs(title = "Isoform Switch Overview",
                        x = "dIF (Isoform Fraction change)",
                        y = "-log10(Isoform Switch q-value)") +
          ggplot2::theme_minimal()
        # Label the top switching genes -- previously this plot's only
        # legend was the significance colour key, with no way to tell which
        # gene any given point was.
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          gene_col <- if ("gene_name" %in% colnames(feat_sw_plot)) "gene_name" else "gene_id"
          feat_sw_plot$.label <- .coalesce_gene_label(feat_sw_plot[[gene_col]], feat_sw_plot$gene_id)
          feat_sig <- feat_sw_plot[!is.na(feat_sw_plot$isoform_switch_q_value) &
                                      feat_sw_plot$isoform_switch_q_value < 0.05 &
                                      abs(feat_sw_plot$dIF) > 0.1, ]
          feat_to_label <- head(feat_sig[order(feat_sig$isoform_switch_q_value), ], 20)
          if (nrow(feat_to_label) > 0) {
            p_switch_overview <- p_switch_overview + ggrepel::geom_text_repel(
              data = feat_to_label,
              mapping = ggplot2::aes(x = dIF, y = -log10(isoform_switch_q_value), label = .label),
              inherit.aes = FALSE, size = 2.8, max.overlaps = 30,
              segment.size = 0.3, color = "black", seed = 42
            )
          }
        } else {
          message("   -> ggrepel not installed; IsoformSwitch_overview.pdf will show points without gene labels.")
        }
        ggplot2::ggsave(file.path(plot_dir, "IsoformSwitch_overview.pdf"),
                        p_switch_overview, width = 7, height = 6)
      }
    }

    # 3f-bis. Genome-wide functional-consequence and alternative-splicing
    # summary plots.
    #
    # IsoformSwitchAnalyzeR's extractConsequenceSummary() / -Enrichment() /
    # extractSplicingSummary() / -Enrichment() do NOT return a ggplot object
    # to capture -- with returnResult = FALSE they return invisible(NULL),
    # and with returnResult = TRUE they return the underlying results
    # data.frame, never the plot itself. plot = TRUE instead makes them draw
    # directly (as a side effect, via print()) to whatever graphics device
    # is current at call time. The previous code assigned the return value
    # to p_*, checked inherits(p_*, "ggplot") -- always FALSE -- and so
    # never called ggsave(); meanwhile the actual plot silently went to
    # whatever device happened to be open (often nothing in a batch run, or
    # a leftover device from earlier code), so *_summary.pdf / *_enrichment.pdf
    # were never written, or written empty. Fixed the same way switchPlot()
    # is already handled a few hundred lines up in this file: open our own
    # pdf() device first via safe_pdf() so the function's internal print()
    # lands exactly where we want it, regardless of its return value.
    safe_run(
      safe_pdf(file.path(plot_dir, "ConsequenceSummary.pdf"), width = 9, height = 6, expr = {
        IsoformSwitchAnalyzeR::extractConsequenceSummary(
          switch_list, consequencesToAnalyze = "all",
          plot = TRUE, returnResult = FALSE
        )
      }),
      label = "Genome-wide consequence summary plot (extractConsequenceSummary)"
    )

    safe_run(
      safe_pdf(file.path(plot_dir, "ConsequenceEnrichment.pdf"), width = 9, height = 6, expr = {
        IsoformSwitchAnalyzeR::extractConsequenceEnrichment(
          switch_list, consequencesToAnalyze = "all",
          plot = TRUE, returnResult = FALSE
        )
      }),
      label = "Genome-wide consequence enrichment plot (extractConsequenceEnrichment)"
    )

    if (!is.null(switch_list[["AlternativeSplicingAnalysis"]])) {
      safe_run(
        safe_pdf(file.path(plot_dir, "SplicingSummary.pdf"), width = 9, height = 6, expr = {
          IsoformSwitchAnalyzeR::extractSplicingSummary(
            switch_list, splicingToAnalyze = "all",
            plot = TRUE, returnResult = FALSE
          )
        }),
        label = "Genome-wide alternative splicing summary plot (extractSplicingSummary)"
      )

      safe_run(
        safe_pdf(file.path(plot_dir, "SplicingEnrichment.pdf"), width = 9, height = 6, expr = {
          IsoformSwitchAnalyzeR::extractSplicingEnrichment(
            switch_list, splicingToAnalyze = "all",
            plot = TRUE, returnResult = FALSE
          )
        }),
        label = "Genome-wide alternative splicing enrichment plot (extractSplicingEnrichment)"
      )
    } else {
      message("  Skipping splicing summary/enrichment plots: no AlternativeSplicingAnalysis found ",
              "in switch_list (analyzeAlternativeSplicing() may have failed upstream -- see earlier messages).")
    }

    if (isTRUE(plot_switch_summary)) {
      safe_run(
        plot_isoform_switch_summary(switch_list, plot_dir = plot_dir,
                                    genes_of_interest = genes_of_interest,
                                    level = level, base = base,
                                    top_n = switch_plot_top_n),
        label = "Isoform switch summary plots"
      )
    }

    # Sashimi-style junction usage + exon-bin coverage comparison, per requested gene
    if (!is.null(genes_of_interest)) {
      for (gene_sym in genes_of_interest) {
        if (isTRUE(plot_sashimi)) {
          safe_run(
            plot_isoform_sashimi(switch_list, gene = gene_sym, level = level, base = base,
                                 plot_dir = plot_dir),
            label = paste0("Sashimi plot for ", gene_sym)
          )
        }
        if (isTRUE(plot_exon_usage)) {
          safe_run(
            plot_exon_usage_comparison(switch_list, gene = gene_sym, level = level, base = base,
                                       plot_dir = plot_dir),
            label = paste0("Exon usage plot for ", gene_sym)
          )
        }
      }
    }
  }

  # 3g. DEXSeq per-gene transcript-usage plots (if the DEXSeq DTU engine was run)
  if (has_dexseq && !is.null(dexseq_results$dxr_list) && !is.null(genes_of_interest)) {
    for (gene_sym in genes_of_interest) {
      gene_ens <- isoform_obj$gene_map$ensembl[isoform_obj$gene_map$symbol == gene_sym]
      if (length(gene_ens) == 0) next
      safe_run(
        plot_dexseq_gene(dexseq_results$dxr_list, gene_id = gene_ens[1],
                         plot_dir = plot_dir, gene_symbol = gene_sym),
        label = paste0("DEXSeq plot for ", gene_sym)
      )
    }
  }
  
  # ---- 4. Render the HTML report -------------------------------------
  # The report is a bundled, parameterized template (inst/rmd/dte_dtu_report.Rmd)
  # rather than an Rmd string built line-by-line here -- see
  # .expressom_rmd_path() / .render_placeholder_template() in utils_core.R for
  # the small sibling template used by the DGE RegionReport step.
  dte_abs_path <- normalizePath(file.path(report_dir, "DTE_results_annotated.csv"), winslash = "/")
  dtu_abs_path <- normalizePath(file.path(report_dir, "DTU_results_annotated.csv"), winslash = "/")
  dexseq_abs_path <- if (has_dexseq) {
    normalizePath(file.path(report_dir, "DEXSeq_results_annotated.csv"), winslash = "/")
  } else {
    NULL
  }

  # Convert every report plot from PDF to PNG up front (browsers can't
  # render a PDF through an <img> tag); the template's .fig() helper picks
  # the .png if present and falls back to .pdf otherwise.
  #
  # This block ALSO writes plot_manifest.csv, independent of anything the
  # Rmd template itself does with these files: a per-figure record of
  # whether the source PDF is suspiciously small/empty and whether PNG
  # conversion succeeded. This exists because "the report shows no
  # plots/data" can happen at several different, easily-confused layers --
  # (a) the figure was never generated, (b) it was generated but the PDF is
  # genuinely empty, (c) the PDF is fine but PNG conversion failed, or (d)
  # both PDF and PNG are fine but the *template*'s own embedding logic
  # doesn't pick them up -- and only this file distinguishes between them
  # without having to open the (multi-MB, mostly-base64) report.html by
  # hand. If every row here shows a reasonable size and "converted = yes"
  # but the report is still blank, the problem is inside dte_dtu_report.Rmd
  # itself, not in anything this function controls.
  plot_pdfs <- list.files(plot_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
  manifest <- data.frame(
    pdf_path       = plot_pdfs,
    size_bytes     = vapply(plot_pdfs, function(f) file.info(f)$size, numeric(1)),
    png_converted  = FALSE,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(plot_pdfs)) {
    if (!is.null(convert_pdf_to_png(plot_pdfs[i]))) manifest$png_converted[i] <- TRUE
  }
  manifest$likely_empty <- manifest$size_bytes < 2000  # a valid single-page PDF is essentially never this small
  manifest_path <- file.path(report_dir, "plot_manifest.csv")
  write.csv(manifest, manifest_path, row.names = FALSE)

  n_empty <- sum(manifest$likely_empty)
  n_unconverted <- sum(!manifest$png_converted)
  if (n_empty > 0) {
    message("   -> WARNING: ", n_empty, " plot PDF(s) are suspiciously small (<2KB, likely empty). See ",
            manifest_path, " for which ones.")
  }
  if (n_unconverted > 0 && nrow(manifest) > 0) {
    note_path <- file.path(report_dir, "README_IF_PLOTS_MISSING.txt")
    note_lines <- c(
      sprintf("%d of %d plot PDF(s) could not be converted to PNG.", n_unconverted, nrow(manifest)),
      "report.html embeds figures as PNG (browsers cannot display PDF through an <img> tag);",
      "any figure that didn't convert will not show up in the HTML report even though its",
      "PDF exists in plots/ and can be opened directly.",
      "",
      "Most common cause: the 'pdftools' R package is not installed. Run:",
      "  install.packages('pdftools')",
      "and re-generate the report (no need to re-run DTE/DTU/switch analysis).",
      "",
      "See plot_manifest.csv in this folder for the exact per-file status."
    )
    writeLines(note_lines, note_path)
    warning("Could not convert ", n_unconverted, " of ", nrow(manifest), " report plot(s) to PNG; ",
            "figures in report.html may not render. See ", note_path, " and ", manifest_path, ".")
  }

  rmd_file <- tempfile(fileext = ".Rmd")
  on.exit(unlink(rmd_file), add = TRUE)
  file.copy(.expressom_rmd_path("dte_dtu_report.Rmd"), rmd_file, overwrite = TRUE)

  rmarkdown::render(
    rmd_file,
    output_file  = "report.html",
    output_dir   = report_dir,
    quiet        = TRUE,
    # IMPORTANT: rmarkdown::render()'s default knit working directory is the
    # directory of the *input* .Rmd file -- since rmd_file is a tempfile()
    # (not inside report_dir), leaving this unset means every relative path
    # used inside the template (e.g. "plots/DTE_volcano.pdf") would resolve
    # against the wrong directory and silently fail its file.exists() guard,
    # so every figure would be skipped even though it was generated
    # correctly on disk. Pinning knit_root_dir to report_dir is what makes
    # the template's relative "plots/..." paths actually resolve.
    knit_root_dir = normalizePath(report_dir, winslash = "/"),
    envir        = new.env(parent = globalenv()),
    params = list(
      level              = level,
      base               = base,
      top_n              = top_n,
      genes_of_interest  = genes_of_interest,
      dte_csv            = dte_abs_path,
      dtu_csv            = dtu_abs_path,
      dexseq_csv         = dexseq_abs_path,
      has_dexseq         = has_dexseq,
      has_switch         = has_switch,
      plot_dir           = normalizePath(plot_dir, winslash = "/"),
      switch_plot_top_n  = switch_plot_top_n
    )
  )

  # Free large objects after report generation
  rm(dte, dtu, dte_sig, dte_ma, dte_filtered, top_dte)
  gc()

  message("DTE/DTU report generated in: ", report_dir)
  invisible(NULL)
}
