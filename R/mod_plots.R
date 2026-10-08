`%||%` <- function(x, y) if (is.null(x)) y else x

## ---------------------------------------------------------------------------
## Sample relabelling + colour helpers shared by the PCA and volcano plots
## ---------------------------------------------------------------------------

#' Normalise a user-supplied sample relabelling into a named character vector
#'
#' Accepts any of:
#' * `NULL` (no relabelling);
#' * a named character vector / named list: names are the sample IDs used in
#'   the sample table, values are the labels to print on the plot;
#' * a data.frame with a sample-ID column and a label column;
#' * the path to a CSV/TSV file holding such a table (delimiter auto-detected,
#'   so Excel's `;`-separated CSV works too).
#'
#' For tables, the ID column is the first one named sample / sample_id / id /
#' old / original / from ..., and the label column the first one named label /
#' new_label / new_name / display / to ...; when no name is recognised the
#' first two columns are used (ID, label). A header row is optional in files:
#' a first row that does not look like a header is read as data.
#'
#' @param sample_labels See above.
#' @return `NULL`, or a named character vector (names = original sample IDs).
#' @keywords internal
.read_sample_labels <- function(sample_labels) {
  if (is.null(sample_labels) || length(sample_labels) == 0) return(NULL)

  id_names    <- c("sample", "sample_id", "sampleid", "sample_name", "samplename",
                   "id", "old", "old_name", "old_label", "from", "original", "current")
  label_names <- c("label", "new_label", "newlabel", "new_name", "newname", "new",
                   "display", "display_name", "plot_label", "to", "rename",
                   "relabel", "alias", "name")
  norm_nm <- function(x) gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))

  ids <- labs <- NULL
  tbl <- NULL
  nm  <- character(0)

  if (is.character(sample_labels) && length(sample_labels) == 1L &&
      is.null(names(sample_labels))) {
    path <- path.expand(sample_labels)
    if (!file.exists(path)) {
      stop("sample_labels file not found: ", sample_labels, call. = FALSE)
    }

    raw <- tryCatch(
      data.table::fread(path, header = FALSE, colClasses = "character",
                        data.table = FALSE, blank.lines.skip = TRUE),
      error = function(e) {
        stop("Could not read sample_labels file '", sample_labels, "': ",
             conditionMessage(e), call. = FALSE)
      }
    )

    if (ncol(raw) < 2L || nrow(raw) == 0L) {
      stop("sample_labels file '", sample_labels, "' needs at least two columns ",
           "(sample ID, new label) and one row.", call. = FALSE)
    }

    first      <- norm_nm(unlist(raw[1, ], use.names = FALSE))
    has_header <- any(first %in% c(id_names, label_names))

    if (has_header) {
      tbl <- raw[-1, , drop = FALSE]
      nm  <- first
    } else {
      tbl <- raw
      nm  <- rep("", ncol(raw))
    }
  } else if (is.data.frame(sample_labels)) {
    if (ncol(sample_labels) < 2L) {
      stop("sample_labels data.frame needs at least two columns (sample ID, new label).",
           call. = FALSE)
    }
    tbl <- sample_labels
    nm  <- norm_nm(names(tbl))
  } else {
    v <- unlist(sample_labels, use.names = TRUE)
    if (is.null(names(v)) || any(is.na(names(v)) | !nzchar(names(v)))) {
      stop("`sample_labels` must be a named character vector (names = sample IDs, ",
           "values = new labels), a data.frame, or the path to a CSV/TSV file.",
           call. = FALSE)
    }
    ids  <- names(v)
    labs <- as.character(v)
  }

  if (is.null(ids)) {
    id_idx <- which(nm %in% id_names)[1]
    if (is.na(id_idx)) id_idx <- 1L

    lab_idx <- which(nm %in% label_names & seq_along(nm) != id_idx)[1]
    if (is.na(lab_idx)) {
      if (ncol(tbl) == 2L || !any(nzchar(nm))) {
        lab_idx <- setdiff(seq_len(ncol(tbl)), id_idx)[1]
      } else {
        stop("sample_labels table has ", ncol(tbl), " columns and none is named ",
             "'label' (or new_label, new_name, display, ...); please name the ",
             "column holding the new labels 'label'.", call. = FALSE)
      }
    }

    ids  <- as.character(tbl[[id_idx]])
    labs <- as.character(tbl[[lab_idx]])
  }

  ids  <- trimws(ids)
  labs <- trimws(labs)
  keep <- !is.na(ids) & nzchar(ids)
  ids  <- ids[keep]
  labs <- labs[keep]

  dup <- duplicated(ids)
  if (any(dup)) {
    warning("sample_labels: duplicated sample ID(s) ignored (first entry kept): ",
            paste(unique(ids[dup]), collapse = ", "), call. = FALSE)
    ids  <- ids[!dup]
    labs <- labs[!dup]
  }

  if (length(ids) == 0L) return(NULL)

  stats::setNames(labs, ids)
}

#' Display labels for a vector of sample IDs
#'
#' Samples without an entry (or with an empty/NA new label) keep their
#' original ID, so a partial mapping is fine.
#' @param sample_ids Character vector of sample IDs as found in the data
#' @param sample_labels Anything accepted by `.read_sample_labels()`
#' @return Character vector, same length as `sample_ids`
#' @keywords internal
.apply_sample_labels <- function(sample_ids, sample_labels = NULL) {
  sample_ids <- as.character(sample_ids)
  lab_map <- .read_sample_labels(sample_labels)
  if (is.null(lab_map)) return(sample_ids)

  out <- unname(lab_map[sample_ids])
  bad <- is.na(out) | !nzchar(out)
  out[bad] <- sample_ids[bad]
  out
}

#' Tell the user which relabelling entries / samples did not line up
#' @keywords internal
.report_unmatched_labels <- function(sample_labels, sample_ids) {
  lab_map <- .read_sample_labels(sample_labels)
  if (is.null(lab_map)) return(invisible(NULL))

  sample_ids <- as.character(sample_ids)
  unknown    <- setdiff(names(lab_map), sample_ids)
  unlabelled <- setdiff(sample_ids, names(lab_map))
  short <- function(x) {
    paste0(paste(utils::head(x, 8), collapse = ", "), if (length(x) > 8) ", ..." else "")
  }

  if (length(unknown) == length(lab_map)) {
    warning("sample_labels: none of the ", length(lab_map),
            " sample ID(s) match the samples in the data (data has e.g. ",
            short(sample_ids), "); PCA plots keep the original sample names.",
            call. = FALSE)
  } else {
    if (length(unknown) > 0) {
      message("   -> sample_labels: ", length(unknown),
              " ID(s) not found in the data and ignored: ", short(unknown))
    }
    if (length(unlabelled) > 0) {
      message("   -> sample_labels: ", length(unlabelled),
              " sample(s) without a new label keep their original name: ", short(unlabelled))
    }
  }

  invisible(NULL)
}

#' Resolve the PCA colours for the contrast levels
#'
#' Defaults: `level` (foreground/treated group) red, `base` (reference/control
#' group) blue. Override with `c(level = "...", base = "...")`; either name
#' may be omitted.
#' @keywords internal
.pca_colors <- function(pca_colors = NULL) {
  defaults <- c(level = "red2", base = "royalblue")
  if (is.null(pca_colors)) return(defaults)

  if (is.null(names(pca_colors)) || !all(names(pca_colors) %in% names(defaults))) {
    stop("`pca_colors` must be a named character vector with names among: level, base.",
         call. = FALSE)
  }

  defaults[names(pca_colors)] <- as.character(pca_colors)
  defaults
}

# Colours for groups that are neither `level` nor `base` (kept clear of red/blue)
.pca_extra_palette <- c("#2CA02C", "#FF7F0E", "#9467BD", "#8C564B",
                        "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")

#' Named colour vector for the groups of a PCA
#'
#' @param groups Factor / character vector of group memberships
#' @param level,base Contrast levels (may be NULL or absent from `groups`)
#' @param pca_colors See `.pca_colors()`
#' @return Character vector of colours named by group
#' @keywords internal
.pca_group_colors <- function(groups, level = NULL, base = NULL, pca_colors = NULL) {
  grp <- if (is.factor(groups)) {
    levels(droplevels(groups))
  } else {
    unique(as.character(groups[!is.na(groups)]))
  }

  cols <- .pca_colors(pca_colors)
  out  <- stats::setNames(rep(NA_character_, length(grp)), grp)

  if (length(level) == 1L && !is.na(level) && level %in% grp) out[[level]] <- cols[["level"]]
  if (length(base)  == 1L && !is.na(base)  && base  %in% grp) out[[base]]  <- cols[["base"]]

  rest <- is.na(out)
  if (any(rest)) {
    n   <- sum(rest)
    pal <- if (n <= length(.pca_extra_palette)) {
      .pca_extra_palette[seq_len(n)]
    } else {
      grDevices::colorRampPalette(.pca_extra_palette)(n)
    }
    out[rest] <- pal
  }

  out
}

#' Resolve the volcano colours (up / down / not significant)
#'
#' Defaults: up-regulated red, down-regulated blue, rest grey. Override with a
#' named vector such as `c(down = "darkgreen")`.
#' @keywords internal
.volcano_colors <- function(colors = NULL) {
  defaults <- .de_direction_colors()
  if (is.null(colors)) return(defaults)

  if (is.null(names(colors)) || !all(names(colors) %in% names(defaults))) {
    stop("`volcano_colors` must be a named character vector with names among: up, down, ns.",
         call. = FALSE)
  }

  defaults[names(colors)] <- as.character(colors)
  defaults
}

.plot_directional_ma <- function(res, title, padj_cutoff, lfc_cutoff) {
  plot_df <- as.data.frame(res)
  plot_df <- plot_df[
    !is.na(plot_df$baseMean) & plot_df$baseMean > 0 &
      is.finite(plot_df$log2FoldChange),
  ]
  plot_df$direction <- .de_direction_label(
    plot_df$log2FoldChange,
    !is.na(plot_df$padj) & plot_df$padj < padj_cutoff &
      abs(plot_df$log2FoldChange) > lfc_cutoff
  )
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = baseMean, y = log2FoldChange, color = direction)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = .de_direction_scale(), drop = FALSE) +
    ggplot2::geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff),
                        color = "grey40", linetype = "dashed") +
    ggplot2::coord_cartesian(ylim = c(-2, 2)) +
    ggplot2::labs(title = title, x = "Mean of normalized counts",
                  y = "log2 Fold Change", color = NULL) +
    ggplot2::theme_minimal()
}


#' Write PCA output tables: per-sample scores, percent variance, and
#' per-sample values for the genes/transcripts that went into the PCA
#'
#' Shared by both the gene-level PCA (run_eda(), via plot_custom_pca()) and
#' the transcript-level PCA (run_isoform_pca()) so the two stay consistent.
#' Previously the only PCA table written to disk was the gene-loadings/
#' variance table (which genes went into the PCA and their variance), with
#' no per-sample values of any kind -- neither each sample's PC coordinates
#' nor each gene's actual expression value in each sample.
#'
#' @param pca_scores data.frame from plot_custom_pca()'s `pca_scores` return
#'   element (or NULL, in which case that file is skipped), with a
#'   "percentVar" attribute set
#' @param gene_values data.frame from plot_custom_pca()'s `gene_values`
#'   return element (or NULL, in which case that file is skipped)
#' @param plot_dir Output directory
#' @param file_stub File name prefix, e.g. "PCA_Treated_vs_Control" -- writes
#'   "<file_stub>_sample_scores.tsv", "<file_stub>_variance_explained.tsv",
#'   and "<file_stub>_top_gene_values.tsv"
#' @keywords internal
.write_pca_scores <- function(pca_scores, plot_dir, file_stub, gene_values = NULL) {
  if (!is.null(pca_scores) && nrow(pca_scores) > 0) {
    scores_path <- file.path(plot_dir, paste0(file_stub, "_sample_scores.tsv"))
    write.table(pca_scores, file = scores_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("   -> Saved per-sample PCA scores (all PCs) to: ", scores_path)

    pct_var <- attr(pca_scores, "percentVar")
    if (!is.null(pct_var)) {
      var_path <- file.path(plot_dir, paste0(file_stub, "_variance_explained.tsv"))
      write.table(data.frame(PC = names(pct_var), PercentVariance = as.numeric(pct_var)),
                  file = var_path, sep = "\t", row.names = FALSE, quote = FALSE)
      message("   -> Saved percent variance explained per PC to: ", var_path)
    }
  }

  if (!is.null(gene_values) && nrow(gene_values) > 0) {
    values_path <- file.path(plot_dir, paste0(file_stub, "_top_gene_values.tsv"))
    write.table(gene_values, file = values_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("   -> Saved per-sample values for the ", nrow(gene_values),
            " gene(s)/transcript(s) used in this PCA to: ", values_path)
  }
  invisible(NULL)
}

#' Exploratory Data Analysis for DESeq2
#'
#' @param dds DESeqDataSet object
#' @param edb Ensembl Database
#' @param out_dir Output directory
#' @param level Level comparison (e.g. "Treated") – may be NULL for EDA-only
#' @param base Base Level (e.g. "Control") – may be NULL
#' @param main_condition The primary modeled condition (could be NULL)
#' @param group_col Optional column name to use for colouring/annotation (overrides main_condition if non-NULL)
#' @param batch_col Optional batch column name for limma correction and before/after PCA
#' @param pca_ntop Number of most variable genes to use for PCA (default 500; set NULL to use all)
#' @param pca_colors Optional named vector `c(level = ..., base = ...)` overriding the PCA colours
#'   (default: `level` red, `base` blue)
#' @param sample_labels Optional relabelling of the samples printed on the PCA plots: a named
#'   character vector (names = sample IDs), a data.frame, or the path to a CSV/TSV (sample ID + label)
#' @export
run_eda <- function(dds, edb, out_dir, level, base,
                    main_condition = NULL, group_col = NULL,
                    batch_col = NULL, pca_ntop = 500,
                    pca_colors = NULL, sample_labels = NULL) {
  cond_col <- group_col %||% main_condition

  sample_labels <- .read_sample_labels(sample_labels)
  .report_unmatched_labels(sample_labels, colnames(dds))

  rld <- DESeq2::rlog(dds, blind = TRUE)

  org_info <- get_organism_info(edb)
  org_obj <- .load_org_db(org_info$org_db)
  message("Mapping IDs for EDA plots...")
  symbols <- tryCatch({
    suppressMessages(AnnotationDbi::mapIds(
      org_obj, keys = rownames(rld), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
    ))
  }, error = function(e) {
    message("   -> Warning: Could not map keys. They may already be symbols. Using original row names.")
    stats::setNames(rownames(rld), rownames(rld))
  })

  symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
  rownames(rld) <- symbols

  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  if (!is.null(level) && !is.null(base)) {
    comp_label <- paste0(level, "_vs_", base)
  } else {
    comp_label <- "EDA"
  }

  if (!is.null(cond_col) && cond_col %in% colnames(SummarizedExperiment::colData(dds))) {
    res_orig <- plot_custom_pca(rld, condition = cond_col, batch = batch_col,
                                 title = paste0("PCA (", cond_col, ")", if (!is.null(level) && !is.null(base)) paste0(" - ", level, " vs ", base) else ""),
                                 return_plot = TRUE, return_gene_list = TRUE,
                                 ntop = pca_ntop,
                                 level = level, base = base,
                                 pca_colors = pca_colors, sample_labels = sample_labels)
    p_orig <- res_orig$plot
    gene_info_orig <- res_orig$gene_info
  } else {
    res_orig <- plot_custom_pca(rld, condition = NULL, batch = batch_col,
                                 title = "PCA (no condition grouping)",
                                 return_plot = TRUE, return_gene_list = TRUE,
                                 ntop = pca_ntop,
                                 level = level, base = base,
                                 pca_colors = pca_colors, sample_labels = sample_labels)
    p_orig <- res_orig$plot
    gene_info_orig <- res_orig$gene_info
  }

  .pdf_device()(file.path(plot_dir, paste0("PCA_", comp_label, ".pdf")), width = 9, height = 7)
  print(p_orig)
  dev.off()

  if (!is.null(gene_info_orig) && nrow(gene_info_orig) > 0) {
    write.table(gene_info_orig,
                file = file.path(plot_dir, paste0("PCA_", comp_label, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    message("   -> Saved top variable genes used in PCA to: ",
            file.path(plot_dir, paste0("PCA_", comp_label, ".tsv")))
  }
  .write_pca_scores(res_orig$pca_scores, plot_dir, paste0("PCA_", comp_label), gene_values = res_orig$gene_values)

  if (!is.null(batch_col) && batch_col %in% colnames(SummarizedExperiment::colData(dds)) &&
      requireNamespace("limma", quietly = TRUE)) {
    batch_vec <- SummarizedExperiment::colData(dds)[[batch_col]]
    if (length(unique(batch_vec)) > 1) {
      if (!is.null(cond_col) && cond_col %in% colnames(SummarizedExperiment::colData(dds))) {
        design_mat <- model.matrix(as.formula(paste0("~ ", cond_col)),
                                   data = as.data.frame(SummarizedExperiment::colData(dds)))
      } else {
        design_mat <- model.matrix(~1, data = as.data.frame(SummarizedExperiment::colData(dds)))
      }
      message("Applying limma batch correction for PCA visualisation...")
      corrected_mat <- tryCatch(
        limma::removeBatchEffect(SummarizedExperiment::assay(rld), batch = batch_vec, design = design_mat),
        error = function(e) { warning("Batch correction failed: ", e$message); NULL }
      )
      if (!is.null(corrected_mat)) {
        rld_corrected <- rld
        SummarizedExperiment::assay(rld_corrected) <- corrected_mat
        res_corr <- plot_custom_pca(rld_corrected, condition = cond_col, batch = batch_col,
                                     title = "PCA After Batch Correction (limma)",
                                     return_plot = TRUE, return_gene_list = TRUE,
                                     ntop = pca_ntop,
                                 level = level, base = base,
                                 pca_colors = pca_colors, sample_labels = sample_labels)
        p_corr <- res_corr$plot
        gene_info_corr <- res_corr$gene_info
        .pdf_device()(file.path(plot_dir, paste0("PCA_BatchCorrected_", comp_label, ".pdf")), width = 9, height = 7)
        print(p_corr)
        dev.off()
        if (!is.null(gene_info_corr) && nrow(gene_info_corr) > 0) {
          write.table(gene_info_corr,
                      file = file.path(plot_dir, paste0("PCA_BatchCorrected_", comp_label, ".tsv")),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          message("   -> Saved top variable genes used in batch-corrected PCA to: ",
                  file.path(plot_dir, paste0("PCA_BatchCorrected_", comp_label, ".tsv")))
        }
        .write_pca_scores(res_corr$pca_scores, plot_dir, paste0("PCA_BatchCorrected_", comp_label), gene_values = res_corr$gene_values)
      }
    }
  }

  rld_cor <- cor(SummarizedExperiment::assay(rld))
  if (anyNA(rld_cor)) {
    warning("Correlation matrix contains NA/NaN values. ",
            "Number of NA: ", sum(is.na(rld_cor)))
    rld_cor[is.na(rld_cor)] <- 0
  }

  if (!is.null(cond_col) && cond_col %in% colnames(SummarizedExperiment::colData(dds))) {
    anno <- as.data.frame(SummarizedExperiment::colData(dds)[, cond_col, drop = FALSE])
    if (!all(rownames(anno) == colnames(rld_cor))) {
      warning("Annotation rownames do not match correlation matrix columns.")
    }
  } else {
    anno <- NULL
  }

  safe_pdf(file.path(plot_dir, paste0("HeatMap_", comp_label, ".pdf")), expr = {
    pheatmap::pheatmap(rld_cor, annotation_col = anno, main = "Sample Correlation (Gene Symbols)")
  })

  message("EDA completed. Plots saved in: ", plot_dir)
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
#' @param pca_ntop Number of most variable genes to use for the supplementary
#'   comparison-only PCA (default 500; set NULL to use all). Matches the
#'   `pca_ntop` used for the main EDA PCA in run_eda().
#' @param pca_colors,sample_labels Passed to the PCA plots; see `plot_custom_pca()`
#' @param volcano_colors Optional named vector overriding the volcano colours (names: up, down, ns)
#' @param lfc_cutoff Absolute log2 fold-change threshold used by the volcano plot (default 1)
#' @return NULL
generate_bulk_visualizations <- function(dds, edb, res_shrunken, res_unshrunken, results_data, out_dir, level, base, main_condition, top_genes, padj_cutoff, highlight_genes = NULL, batch_col = NULL, pca_ntop = 500,
                                       pca_colors = NULL, sample_labels = NULL,
                                       volcano_colors = NULL, lfc_cutoff = 1) {
  plot_dir <- file.path(out_dir, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  org_info <- get_organism_info(edb)
  org_obj <- .load_org_db(org_info$org_db)


  if (!is.null(main_condition) && main_condition %in% colnames(SummarizedExperiment::colData(dds))) {
    vsd <- tryCatch(DESeq2::vst(dds, blind = TRUE),
                    error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
    plot_sample_correlation(vsd, main_condition, plot_dir, paste0(level, "_vs_", base))

    # Supplementary PCA restricted to just this comparison's samples (per
    # the sample table's main_condition column) -- see plot_comparison_pca().
    plot_comparison_pca(vsd, main_condition, level, base, batch = batch_col,
                        plot_dir = plot_dir, ntop = pca_ntop,
                        pca_colors = pca_colors, sample_labels = sample_labels)
  }

  safe_pdf(file.path(plot_dir, paste0("MAplot_unshrunken_", level, "_vs_", base, ".pdf")), expr = {
    print(.plot_directional_ma(res_unshrunken, "MA plot (unshrunken)",
                              padj_cutoff, lfc_cutoff))
  })

  safe_pdf(file.path(plot_dir, paste0("MAplot_shrunken_", level, "_vs_", base, ".pdf")), expr = {
    print(.plot_directional_ma(res_shrunken, "MA plot (shrunken)",
                              padj_cutoff, lfc_cutoff))
  })

  topx_sigOE_genes <- head(results_data$sig_res[order(results_data$sig_res$padj), "gene"], top_genes)
  topx_sigOE_norm <- results_data$normalized_counts[results_data$normalized_counts$gene %in% topx_sigOE_genes, ]

  id_cols <- intersect(c("gene", "ensembl", "symbol", "entrezid"), colnames(topx_sigOE_norm))
  gathered_top <- tidyr::pivot_longer(
    topx_sigOE_norm,
    cols      = -dplyr::all_of(id_cols),
    names_to  = "sample",
    values_to = "normalized_counts"
  )
  gathered_top$normalized_counts <- as.numeric(gathered_top$normalized_counts)

  meta_df <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_df$sample <- rownames(meta_df)
  topx_final <- dplyr::inner_join(meta_df, gathered_top, by = "sample")

  # Top-N DE gene boxplot. Colours follow the same convention as every other
  # plot in the pipeline (level = red, base = blue) via .pca_group_colors(),
  # which also handles the case where dds carries a third condition level
  # beyond level/base (those get a distinct palette colour instead of the
  # default ggplot hue). Previously this plot had no scale_fill_manual at
  # all, so it silently used ggplot's default palette.
  safe_pdf(file.path(plot_dir, paste0("Top", top_genes, "_DE_Genes_", level, "_vs_", base, ".pdf")), expr = {
    group_cols <- .pca_group_colors(
      factor(topx_final[[main_condition]]),
      level = level,
      base  = base
    )

    p_top <- ggplot2::ggplot(topx_final, ggplot2::aes(x = gene, y = normalized_counts, fill = .data[[main_condition]])) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_fill_manual(values = group_cols, name = main_condition) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ggtitle(paste("Top", top_genes, "Significant DE Genes"))
    print(p_top)
  })

  safe_pdf(file.path(plot_dir, paste0("DE_Volcanoplot_", level, "_vs_", base, ".pdf")), expr = {
    suppressWarnings(
      print(plot_volcano(results_data$res_tbl, padj_cutoff, highlight_genes,
                         title = paste(level, "vs", base),
                         lfc_cutoff = lfc_cutoff, colors = volcano_colors))
    )
  })

  plot_top_genes_heatmap(dds, results_data, main_condition, level, base,
                         top_n = top_genes, padj_cutoff = padj_cutoff,
                         plot_dir = plot_dir, batch_col = batch_col)
}


#' Custom PCA Plot (ggplot2, fully customisable)
#' @param vsd VST-transformed DESeqDataSet or matrix
#' @param condition Character column name in colData for colour grouping (can be NULL)
#' @param batch Optional batch column for shape grouping
#' @param title Plot title
#' @param ellipse Logical, whether to add 95% confidence ellipses
#' @param return_plot If TRUE, returns ggplot object; otherwise prints
#' @param ntop Number of top variable genes to use for PCA (NULL = all genes)
#' @param return_gene_list If TRUE, return a list with plot and gene_info; otherwise only the plot
#' @param level,base Optional contrast levels. When they are values of
#'   `condition`, the `level` group is drawn in `pca_colors[["level"]]` (red by
#'   default) and the `base` group in `pca_colors[["base"]]` (blue by default);
#'   any other group gets a distinct palette colour.
#' @param pca_colors Optional named character vector `c(level = ..., base = ...)`
#'   overriding the default colours (`"red2"` / `"royalblue"`).
#' @param sample_labels Optional relabelling of the sample names printed on the
#'   plot: a named character vector (names = sample IDs, values = new labels),
#'   a data.frame, or the path to a CSV/TSV with a sample-ID column and a label
#'   column. Samples without an entry keep their original name. The score
#'   tables keep the original ID in `sample_label` and add the printed label as
#'   `plot_label`.
#' @export
plot_custom_pca <- function(vsd, condition, batch = NULL, title = "PCA", ellipse = TRUE,
                            return_plot = TRUE, ntop = NULL, return_gene_list = FALSE,
                            level = NULL, base = NULL, pca_colors = NULL,
                            sample_labels = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  if (!inherits(vsd, "SummarizedExperiment")) {
    stop("vsd must be a SummarizedExperiment object (e.g., DESeqDataSet or DESeqTransform)")
  }

  mat <- SummarizedExperiment::assay(vsd)

  # optional = TRUE keeps non-syntactic column names ("cell type", "batch-id")
  # as they are. The default turns them into cell.type / batch.id, and the
  # `condition %in% colnames(...)` test below then silently fell back to an
  # ungrouped PCA for any such column.
  coldata <- as.data.frame(SummarizedExperiment::colData(vsd), optional = TRUE)

  gene_info <- NULL
  if (!is.null(ntop) && is.numeric(ntop) && ntop > 0 && ntop < nrow(mat)) {
    row_var <- apply(mat, 1, var, na.rm = TRUE)
    row_var[is.na(row_var)] <- -Inf
    top_idx <- order(row_var, decreasing = TRUE)[seq_len(min(ntop, nrow(mat)))]
    gene_names <- rownames(mat)[top_idx]
    gene_var   <- row_var[top_idx]
    mat <- mat[top_idx, , drop = FALSE]
    gene_info <- data.frame(gene = gene_names, variance = gene_var, stringsAsFactors = FALSE)
  } else {
    if (return_gene_list) {
      row_var <- apply(mat, 1, var, na.rm = TRUE)
      row_var[is.na(row_var)] <- NA_real_
      gene_info <- data.frame(gene = rownames(mat), variance = row_var, stringsAsFactors = FALSE)
    }
  }

  if (ncol(mat) < 2 || nrow(mat) < 2) {
    stop("PCA needs at least 2 samples and 2 genes/transcripts (got ",
         ncol(mat), " samples x ", nrow(mat), " features).")
  }

  pca        <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  sample_ids <- rownames(pca$x)

  pca_df <- coldata
  rownames(pca_df) <- sample_ids
  pca_df$PC1 <- pca$x[, 1]
  pca_df$PC2 <- pca$x[, 2]
  pca_df$sample_label <- .apply_sample_labels(sample_ids, sample_labels)

  score_front <- data.frame(sample_label = sample_ids, stringsAsFactors = FALSE)
  if (!is.null(sample_labels)) score_front$plot_label <- pca_df$sample_label
  pca_scores <- data.frame(score_front, coldata, pca$x, row.names = NULL, check.names = FALSE)
  attr(pca_scores, "percentVar") <- stats::setNames(percentVar, colnames(pca$x))

  gene_values <- data.frame(gene = rownames(mat), as.data.frame(mat, check.names = FALSE),
                            row.names = NULL, check.names = FALSE)

  group_cols <- NULL
  if (is.null(condition) || !(condition %in% colnames(pca_df))) {
    if (!is.null(condition)) {
      message("   -> PCA: column '", condition, "' not found in colData; plotting without grouping.")
    }
    pca_df$Group <- "All Samples"
    condition    <- "Group"
    ellipse      <- FALSE
    group_cols   <- c("All Samples" = "grey30")
  }

  discrete <- !is.numeric(pca_df[[condition]])
  if (discrete) {
    pca_df[[condition]] <- if (is.factor(pca_df[[condition]])) {
      droplevels(pca_df[[condition]])
    } else {
      factor(pca_df[[condition]])
    }
    if (is.null(group_cols)) {
      group_cols <- .pca_group_colors(pca_df[[condition]], level, base, pca_colors)
    }
  }

  # A numeric batch column cannot be mapped to shape, and ggplot silently
  # drops every point beyond the 6th level of a discrete shape scale.
  use_batch <- !is.null(batch) && batch %in% colnames(pca_df)
  if (use_batch) pca_df[[batch]] <- droplevels(as.factor(pca_df[[batch]]))

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
                                             color = .data[[condition]]))

  # stat_ellipse() needs at least 4 points per group (it computes a robust
  # covariance on n - 1 >= 3 degrees of freedom). With the usual n = 3 per
  # group the old `>= 3` guard still added the layer, which only produced
  # "Too few points to calculate an ellipse" messages and geom_path warnings.
  if (ellipse && discrete && nlevels(pca_df[[condition]]) > 1 &&
      min(table(pca_df[[condition]])) >= 4) {
    p <- p + ggplot2::stat_ellipse(level = 0.95, linetype = 2)
  }

  p <- p + ggplot2::geom_point(size = 3.5, alpha = 0.8)

  if (use_batch) {
    shapes <- c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10, 11, 12, 13, 14)
    p <- p +
      ggplot2::aes(shape = .data[[batch]]) +
      ggplot2::scale_shape_manual(values = rep_len(shapes, nlevels(pca_df[[batch]])))
  }

  if (discrete) {
    p <- p + ggplot2::scale_color_manual(values = group_cols, na.value = "grey60")
  }

  p <- p +
    ggrepel::geom_text_repel(ggplot2::aes(label = .data[["sample_label"]]),
                             size = 3, show.legend = FALSE,
                             box.padding = 0.4, max.overlaps = Inf) +
    ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = title)

  if (return_gene_list) {
    return(list(plot = p, gene_info = gene_info, pca_scores = pca_scores, gene_values = gene_values))
  } else {
    if (return_plot) return(p) else { print(p); invisible(p) }
  }
}

#' Supplementary PCA restricted to the samples in the current comparison
#'
#' run_eda()'s PCA is computed once, up front, across every sample in the
#' sample table -- including groups that aren't part of the current level
#' vs. base contrast (e.g. a third group sitting alongside a two-group
#' comparison). That's the right view for a global overview, but variance
#' coming from the excluded group(s) can dilute or obscure how well `level`
#' and `base` actually separate from each other. This subsets the
#' already-computed `vsd`/`rld` down to just those two groups (per the
#' sample table's condition column) and re-runs plot_custom_pca() on that
#' subset alone, so the comparison can be inspected on its own.
#'
#' Does not recompute the variance-stabilizing transform -- it subsets the
#' one already fit across the full cohort, which is more stable for small
#' comparison subsets than re-fitting vst()/rlog() dispersion from scratch
#' on just a handful of samples.
#'
#' @param vsd VST/rlog-transformed DESeqDataSet (already computed for ALL
#'   samples)
#' @param condition_col Column name in colData holding the contrasted
#'   factor -- i.e. main_condition, the column `level`/`base` are values of
#' @param level Foreground group of the comparison
#' @param base Reference group of the comparison
#' @param batch Optional batch column for shape grouping
#' @param plot_dir Output directory for the PDF/TSVs
#' @param ntop Number of most variable genes to use for PCA (NULL = all genes)
#' @param pca_colors,sample_labels See `plot_custom_pca()`
#' @return The ggplot object, invisibly (NULL if skipped)
#' @export
plot_comparison_pca <- function(vsd, condition_col, level, base, batch = NULL,
                                plot_dir, ntop = 500,
                                pca_colors = NULL, sample_labels = NULL) {
  coldata <- as.data.frame(SummarizedExperiment::colData(vsd))

  if (is.null(condition_col) || !condition_col %in% colnames(coldata)) {
    message("   -> Skipping comparison-only PCA: '", condition_col, "' not found in colData.")
    return(invisible(NULL))
  }

  keep <- coldata[[condition_col]] %in% c(level, base)

  if (sum(keep) < 2) {
    message("   -> Skipping comparison-only PCA: fewer than 2 samples match '",
            level, "' / '", base, "'.")
    return(invisible(NULL))
  }

  vsd_sub <- vsd[, keep]

  # Drop unused factor levels left over from the full cohort (e.g. a third
  # group not part of this contrast) -- otherwise plot_custom_pca()'s
  # min(table(condition)) ellipse check sees a phantom zero-count group and
  # silently disables the ellipse even when level/base both qualify.
  if (is.factor(SummarizedExperiment::colData(vsd_sub)[[condition_col]])) {
    SummarizedExperiment::colData(vsd_sub)[[condition_col]] <-
      droplevels(SummarizedExperiment::colData(vsd_sub)[[condition_col]])
  }

  n_groups <- length(unique(SummarizedExperiment::colData(vsd_sub)[[condition_col]]))
  if (n_groups < 2) {
    message("   -> Skipping comparison-only PCA: only one group present after ",
            "subsetting to '", level, "' / '", base, "'.")
    return(invisible(NULL))
  }

  res <- plot_custom_pca(vsd_sub, condition = condition_col, batch = batch,
                         title = paste0("PCA (Comparison Samples Only) - ", level, " vs ", base),
                         return_plot = TRUE, return_gene_list = TRUE,
                         ntop = ntop, level = level, base = base,
                         pca_colors = pca_colors, sample_labels = sample_labels)

  file_stub <- paste0("PCA_ComparisonOnly_", level, "_vs_", base)

  safe_pdf(file.path(plot_dir, paste0(file_stub, ".pdf")), width = 9, height = 7,
          expr = print(res$plot))

  if (!is.null(res$gene_info) && nrow(res$gene_info) > 0) {
    write.table(res$gene_info,
               file = file.path(plot_dir, paste0(file_stub, ".tsv")),
               sep = "\t", row.names = FALSE, quote = FALSE)
    message("   -> Saved top variable genes used in comparison-only PCA to: ",
           file.path(plot_dir, paste0(file_stub, ".tsv")))
  }

  .write_pca_scores(res$pca_scores, plot_dir, file_stub, gene_values = res$gene_values)

  message("   -> Comparison-only PCA (", level, " vs ", base, ", n=", sum(keep),
          " samples) saved to: ", plot_dir)

  invisible(res$plot)
}

#' Sample Correlation Heatmap
#' @param vsd VST-transformed DESeqDataSet
#' @param condition_col Column name for sample annotation
#' @param plot_dir Directory to save PDF
#' @param comp_name Comparison name for file naming
#' @export
plot_sample_correlation <- function(vsd, condition_col, plot_dir, comp_name) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed, skipping correlation plot")
    return(invisible(NULL))
  }
  expr_mat   <- SummarizedExperiment::assay(vsd)
  cor_mat    <- cor(expr_mat, method = "pearson")
  annotation <- as.data.frame(SummarizedExperiment::colData(vsd)[, condition_col, drop = FALSE])
  rownames(annotation) <- colnames(cor_mat)
  colnames(annotation) <- condition_col

  safe_pdf(file.path(plot_dir, paste0("SampleCorrelation_", comp_name, ".pdf")), width = 8, height = 7, expr = {
    pheatmap::pheatmap(cor_mat,
                       annotation_col  = annotation,
                       main            = paste0("Sample Correlation Matrix - ", comp_name),
                       color           = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                       cluster_rows    = TRUE, cluster_cols = TRUE,
                       display_numbers = FALSE, fontsize_row = 8)
  })
  message("   -> Sample correlation heatmap saved to: ", plot_dir)
}

#' Top Genes Expression Heatmap
#'
#' Selects the top-N significant DE genes by padj, retrieves their normalised
#' counts from the DESeqDataSet using Ensembl IDs (the actual rownames of the
#' count matrix), then labels rows with gene symbols for readability.
#'
#' @param dds DESeqDataSet
#' @param results_data List containing res_tbl (data.frame with gene, ensembl, log2FoldChange, padj)
#' @param condition_col Column name for condition grouping
#' @param level Treatment group
#' @param base Control group
#' @param top_n Number of top DE genes to plot (by padj)
#' @param padj_cutoff Adjusted p-value threshold
#' @param plot_dir Output directory
#' @param batch_col Optional batch column for expression correction
#' @export
plot_top_genes_heatmap <- function(dds, results_data, condition_col, level, base,
                                   top_n = 30, padj_cutoff = 0.01,
                                   plot_dir, batch_col = NULL) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed, skipping top genes heatmap")
    return(invisible(NULL))
  }

  res_tbl <- results_data$res_tbl
  sig     <- res_tbl[which(!is.na(res_tbl$padj) & res_tbl$padj < padj_cutoff), ]
  if (nrow(sig) == 0) {
    message("No significant genes at padj < ", padj_cutoff, " – skipping top genes heatmap")
    return(invisible(NULL))
  }
  top_sig <- head(sig[order(sig$padj), ], n = top_n)

  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  if ("ensembl" %in% colnames(top_sig) && any(!is.na(top_sig$ensembl))) {
    candidate_ids <- top_sig$ensembl[!is.na(top_sig$ensembl)]
    present_ids   <- intersect(candidate_ids, rownames(norm_counts))
    if (length(present_ids) == 0) {
      message("None of the top genes' Ensembl IDs found in counts matrix – skipping heatmap")
      return(invisible(NULL))
    }
    mat         <- norm_counts[present_ids, , drop = FALSE]
    sym_labels  <- top_sig$gene[match(present_ids, top_sig$ensembl)]
    sym_labels[is.na(sym_labels)] <- present_ids[is.na(sym_labels)]
    rownames(mat) <- make.unique(sym_labels)
  } else {
    present_ids <- intersect(top_sig$gene, rownames(norm_counts))
    if (length(present_ids) == 0) {
      message("None of the top genes found in counts matrix – skipping heatmap")
      return(invisible(NULL))
    }
    mat <- norm_counts[present_ids, , drop = FALSE]
  }

  if (!is.null(batch_col) && batch_col %in% colnames(SummarizedExperiment::colData(dds)) &&
      requireNamespace("limma", quietly = TRUE)) {
    batch_vec <- SummarizedExperiment::colData(dds)[[batch_col]]
    if (length(unique(batch_vec)) > 1) {
      design_mat <- model.matrix(as.formula(paste0("~ ", condition_col)),
                                 data = as.data.frame(SummarizedExperiment::colData(dds)))
      mat <- tryCatch(
        limma::removeBatchEffect(mat, batch = batch_vec, design = design_mat),
        error = function(e) mat
      )
    }
  }

  row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  zero_var <- row_sds == 0 | is.na(row_sds)
  if (any(zero_var)) {
    message("   -> plot_top_genes_heatmap: removing ", sum(zero_var),
            " zero-variance gene(s) before z-score scaling")
    mat <- mat[!zero_var, , drop = FALSE]
  }
  if (nrow(mat) == 0) {
    message("   -> No genes with non-zero variance – skipping heatmap")
    return(invisible(NULL))
  }
  mat_z            <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0

  coldata        <- as.data.frame(SummarizedExperiment::colData(dds))
  annotation_col <- coldata[, condition_col, drop = FALSE]
  rownames(annotation_col) <- colnames(mat)
  colnames(annotation_col) <- condition_col

  col_pal <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

  .pdf_device()(file.path(plot_dir, paste0("TopGenes_Heatmap_", level, "_vs_", base, ".pdf")),
      width = 10, height = max(6, nrow(mat_z) * 0.3))
  pheatmap::pheatmap(mat_z,
                     annotation_col = annotation_col,
                     main           = paste0("Top ", nrow(mat_z), " DE genes (", level, " vs ", base, ")"),
                     color          = col_pal,
                     cluster_rows   = TRUE, cluster_cols = TRUE,
                     scale          = "none",
                     fontsize_row   = 8,
                     show_rownames  = TRUE)
  dev.off()
  message("   -> Top genes heatmap saved to: ", plot_dir)
}

#' Plot Volcano
#'
#' Points are coloured by direction: up-regulated (padj < `padj_cutoff` and
#' log2FC > `lfc_cutoff`) in `colors[["up"]]` (red by default), down-regulated
#' (log2FC < -`lfc_cutoff`) in `colors[["down"]]` (blue by default), everything
#' else grey. The caption reports the number of up- and down-regulated genes
#' instead of EnhancedVolcano's default "total = N variables".
#'
#' @param res_tbl Results table with `gene`, `log2FoldChange` and `padj` columns
#' @param padj_cutoff Adjusted p-value threshold
#' @param highlight_genes Optional character vector of genes to label
#' @param title Plot title
#' @param lfc_cutoff Absolute log2 fold-change threshold (vertical lines,
#'   colouring and the counts in the caption)
#' @param colors Optional named character vector overriding the default colours;
#'   names among `up`, `down`, `ns`
#' @keywords internal
plot_volcano <- function(res_tbl, padj_cutoff, highlight_genes = NULL, title = "",
                         lfc_cutoff = 1, colors = NULL) {
  cols    <- .volcano_colors(colors)
  res_tbl <- as.data.frame(res_tbl)

  sig <- !is.na(res_tbl$padj) & !is.na(res_tbl$log2FoldChange) &
    res_tbl$padj < padj_cutoff & abs(res_tbl$log2FoldChange) > lfc_cutoff
  is_up   <- sig & res_tbl$log2FoldChange > 0
  is_down <- sig & res_tbl$log2FoldChange < 0
  n_up    <- sum(is_up)
  n_down  <- sum(is_down)

  # EnhancedVolcano colours each point by the per-row value of `colCustom` and
  # builds the legend from its names.
  key_col <- ifelse(is_up, cols[["up"]], ifelse(is_down, cols[["down"]], cols[["ns"]]))
  names(key_col) <- ifelse(is_up, "Upregulated",
                           ifelse(is_down, "Downregulated", "Not significant"))

  caption <- paste0("Upregulated: ", n_up, "  |  Downregulated: ", n_down,
                    "   (padj < ", format(padj_cutoff), ", |log2FC| > ", format(lfc_cutoff), ")")

  EnhancedVolcano::EnhancedVolcano(
    res_tbl, lab = res_tbl$gene,
    selectLab      = highlight_genes,
    drawConnectors = !is.null(highlight_genes),
    x = "log2FoldChange", y = "padj",
    title     = title,
    caption   = caption,
    pCutoff   = padj_cutoff, FCcutoff = lfc_cutoff, pointSize = 2.0, labSize = 4.0,
    colCustom = key_col
  )
}

#' Plot Individual Sample Z-Score Heatmap
#' @param dds A DESeqDataSet object
#' @param selected_genes Character vector of gene symbols to plot
#' @param condition_col Character string representing the design metadata column
#' @param level String representing the foreground group
#' @param base String representing the background reference group
#' @param plot_dir String directory path where the PDF will be saved
#' @export
plot_sample_zscore <- function(dds, selected_genes, condition_col, level, base, plot_dir) {
  norm_counts   <- DESeq2::counts(dds, normalized = TRUE)
  present_genes <- intersect(selected_genes, rownames(norm_counts))
  if (length(present_genes) == 0) stop("None of the specified genes were found.")
  meta_data          <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_data$sample   <- rownames(meta_data)
  valid_samples      <- meta_data$sample[meta_data[[condition_col]] %in% c(base, level)]
  expr_mat           <- norm_counts[present_genes, valid_samples, drop = FALSE]
  z_mat              <- .zscore_matrix(expr_mat)
  df_long            <- as.data.frame(z_mat)
  df_long$gene       <- rownames(df_long)
  plot_df            <- tidyr::pivot_longer(df_long, cols = -gene, names_to = "sample", values_to = "z_score")
  plot_df            <- merge(plot_df, meta_data, by = "sample")
  plot_df$Condition  <- factor(plot_df[[condition_col]], levels = c(base, level))
  plot_df$gene       <- factor(plot_df$gene, levels = rev(present_genes))
  comp_name          <- paste0(level, "_vs_", base)
  pdf_path           <- file.path(plot_dir, paste0("Sample_Zscore_Heatmap_", comp_name, ".pdf"))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = gene, fill = z_score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::facet_grid(~Condition, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "dodgerblue4", mid = "white", high = "red3", midpoint = 0,
      name = "z-score", limits = c(-1.5, 1.5), oob = scales::squish,
      guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(z_score, 1)), color = "black", size = 2.5) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 9),
      axis.text.y = ggplot2::element_text(face = "italic", size = 10),
      panel.grid  = ggplot2::element_blank(),
      strip.text  = ggplot2::element_text(face = "bold", size = 12, margin = ggplot2::margin(b = 10)),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::labs(title = "Sample Expression Heatmap", x = NULL, y = NULL)
  calc_height <- max(4, length(present_genes) * 0.4)
  calc_width  <- max(5, length(valid_samples) * 0.6 + 2)
  ggplot2::ggsave(filename = pdf_path, plot = p, width = calc_width, height = calc_height, device = .pdf_device())
  message("   -> Sample Heatmap successfully exported to: ", pdf_path)
}

#' Plot Log2 Fold Change Heatmap (1-Column Format)
#' @param dds A DESeqDataSet object
#' @param selected_genes Character vector of gene symbols to plot
#' @param condition_col Character string representing the design metadata column
#' @param level String representing the foreground group
#' @param base String representing the background reference group
#' @param plot_dir String directory path where the PDF will be saved
#' @export
plot_l2fc_heatmap <- function(dds, selected_genes, condition_col, level, base, plot_dir) {
  res         <- DESeq2::results(dds, contrast = c(condition_col, level, base))
  res_df      <- as.data.frame(res)
  present_genes <- intersect(selected_genes, rownames(res_df))
  if (length(present_genes) == 0) stop("None of the specified genes were found in the results.")
  plot_df     <- data.frame(gene = present_genes, L2FC = res_df[present_genes, "log2FoldChange"])
  plot_df$gene <- factor(plot_df$gene, levels = rev(present_genes))
  plot_df$Comparison <- paste0(level, " vs ", base)
  comp_name   <- paste0(level, "_vs_", base)
  pdf_path    <- file.path(plot_dir, paste0("L2FC_Heatmap_", comp_name, ".pdf"))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Comparison, y = gene, fill = L2FC)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = .de_direction_colors()[["down"]],
      mid = "white",
      high = .de_direction_colors()[["up"]],
      midpoint = 0,
      name = "Log2 FC",
      guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(L2FC, 2)), color = "black", size = 3.5) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "italic", size = 10),
      panel.grid  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold", margin = ggplot2::margin(b = 15))
    ) +
    ggplot2::labs(title = "Log2 Fold Change Heatmap", x = NULL, y = NULL)
  calc_height <- max(4, length(present_genes) * 0.4)
  ggplot2::ggsave(filename = pdf_path, plot = p, width = 4, height = calc_height, device = .pdf_device())
  message("   -> L2FC Heatmap successfully exported to: ", pdf_path)
}

#' Plot Average Z-Score by Condition for Gene Set(s)
#'
#' Computes a per‑gene z‑score across samples, then averages the z‑scores across
#' all genes within each gene set for every sample, producing one "module score"
#' per sample per set. If `set_name` is supplied, only one set is expected and
#' the plot is saved as a single PDF (without faceting) using that name in the
#' title and filename.
#'
#' @param dds A DESeqDataSet object (rownames should already be gene symbols)
#' @param gene_sets Named list of character vectors, e.g.
#'   \code{list(Tightness = c("Cdh5","Pdgfa"), "Lipid Scavengers" = c("Cd36","Stab1"))}.
#'   List names are used as facet/panel titles unless `set_name` is provided.
#' @param condition_col Character string representing the design metadata column
#' @param level String representing the foreground group
#' @param base String representing the background reference group
#' @param plot_dir String directory path where the PDF will be saved
#' @param include_global Logical: append a pooled "Global (All Genes)" panel
#'   combining every gene across all supplied sets (default TRUE, ignored if
#'   `set_name` is not NULL)
#' @param set_name Optional character string. If provided, `gene_sets` must be a
#'   list of length 1; this name is used in the output filename and plot title,
#'   and faceting is suppressed.
#' @return Invisibly returns the ggplot object
#' @export
plot_geneset_zscore_avg <- function(dds, gene_sets, condition_col, level, base,
                                    plot_dir, include_global = TRUE,
                                    set_name = NULL) {
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("`gene_sets` must be a named list, e.g. list(SetA = c('Gene1','Gene2'), SetB = c('Gene3')).")
  }

  if (!is.null(set_name)) {
    if (length(gene_sets) != 1) {
      warning("set_name provided but gene_sets has length != 1. Using only the first set.")
      gene_sets <- gene_sets[1]
    }
    include_global <- FALSE
  }

  norm_counts    <- DESeq2::counts(dds, normalized = TRUE)
  meta_data      <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_data$sample <- rownames(meta_data)
  valid_samples  <- meta_data$sample[meta_data[[condition_col]] %in% c(base, level)]

  .avg_zscore_for_genes <- function(genes) {
    present_genes <- intersect(genes, rownames(norm_counts))
    missing_genes <- setdiff(genes, rownames(norm_counts))
    if (length(present_genes) == 0) {
      return(list(avg_z = NULL, present = present_genes, missing = missing_genes))
    }
    expr_mat              <- norm_counts[present_genes, valid_samples, drop = FALSE]
    z_mat                 <- .zscore_matrix(expr_mat)
    list(avg_z = colMeans(z_mat, na.rm = TRUE), present = present_genes, missing = missing_genes)
  }

  set_results    <- list()
  missing_report <- list()

  for (set_name_ in names(gene_sets)) {
    res <- .avg_zscore_for_genes(gene_sets[[set_name_]])
    if (length(res$missing) > 0) missing_report[[set_name_]] <- res$missing
    if (is.null(res$avg_z)) {
      warning("None of the genes in gene set '", set_name_, "' were found in the data. Skipping this set.")
      next
    }
    set_results[[set_name_]] <- data.frame(
      sample     = names(res$avg_z),
      avg_zscore = as.numeric(res$avg_z),
      gene_set   = paste0(set_name_, " (n=", length(res$present), ")"),
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(include_global) && length(gene_sets) > 1) {
    all_genes_pooled <- unique(unlist(gene_sets, use.names = FALSE))
    res_global <- .avg_zscore_for_genes(all_genes_pooled)
    if (!is.null(res_global$avg_z)) {
      set_results[["__global__"]] <- data.frame(
        sample     = names(res_global$avg_z),
        avg_zscore = as.numeric(res_global$avg_z),
        gene_set   = paste0("Global (All Genes) (n=", length(res_global$present), ")"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(missing_report) > 0) {
    for (nm in names(missing_report)) {
      message("   -> Note: genes not found in '", nm, "' set (check spelling/case): ",
              paste(missing_report[[nm]], collapse = ", "))
    }
  }

  if (length(set_results) == 0) stop("None of the specified gene sets contained any matching genes.")

  plot_df           <- do.call(rbind, set_results)
  plot_df           <- merge(plot_df, meta_data[, c("sample", condition_col)], by = "sample")
  plot_df$Condition <- factor(plot_df[[condition_col]], levels = c(base, level))

  if (is.null(set_name)) {
    set_levels <- vapply(names(set_results), function(nm) unique(set_results[[nm]]$gene_set), character(1))
    plot_df$gene_set <- factor(plot_df$gene_set, levels = set_levels)
  }

  summary_df <- do.call(rbind, lapply(split(plot_df, list(plot_df$gene_set, plot_df$Condition)), function(d) {
    if (nrow(d) == 0) return(NULL)
    data.frame(
      gene_set  = d$gene_set[1],
      Condition = d$Condition[1],
      mean_z    = mean(d$avg_zscore),
      sem       = if (nrow(d) > 1) stats::sd(d$avg_zscore) / sqrt(nrow(d)) else 0
    )
  }))

  comp_name <- paste0(level, "_vs_", base)

  # Build base plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = summary_df,
      ggplot2::aes(x = Condition, y = mean_z, fill = Condition),
      stat = "identity", alpha = 0.6, width = 0.6
    ) +
    ggplot2::geom_errorbar(
      data = summary_df,
      ggplot2::aes(x = Condition, ymin = mean_z - sem, ymax = mean_z + sem),
      width = 0.2
    ) +
    ggplot2::geom_jitter(
      data = plot_df,
      ggplot2::aes(x = Condition, y = avg_zscore, fill = Condition),
      width = 0.15, size = 2, shape = 21, color = "black"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = stats::setNames(c("dodgerblue4", "red3"), c(base, level))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text      = ggplot2::element_text(face = "bold", size = 11),
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  # Add faceting or single‑set title
  if (!is.null(set_name)) {
    p <- p + ggplot2::labs(
      title = paste0("Average Z‑score: ", set_name, " (", level, " vs ", base, ")"),
      x = NULL, y = "Mean Z‑score (\u00b1 SEM)"
    )
    pdf_path <- file.path(plot_dir, paste0("GeneSet_Zscore_Average_", comp_name, "_", set_name, ".pdf"))
  } else {
    p <- p + ggplot2::facet_wrap(~gene_set, scales = "free_y") +
      ggplot2::labs(
        title = "Average Z‑score by Condition (Gene Set Module Score)",
        x = NULL, y = "Mean Z‑score (\u00b1 SEM)"
      )
    pdf_path <- file.path(plot_dir, paste0("GeneSet_Zscore_Average_", comp_name, ".pdf"))
  }

  # Save with appropriate dimensions
  if (!is.null(set_name)) {
    width <- 6
  } else {
    width <- 4.5 * length(set_results) + 1
  }
  ggplot2::ggsave(filename = pdf_path, plot = p, width = width, height = 5, device = .pdf_device())
  message("   -> Gene set average Z‑score plot saved to: ", pdf_path)
  return(invisible(p))
}

#' Plot Z-Score by Condition for Individual Genes
#'
#' Companion to \code{plot_geneset_zscore_avg}: computes the same per‑gene
#' z‑score across samples and plots a grouped bar chart (mean ± SEM) with
#' every gene along the x‑axis, faceted by gene set if multiple sets are
#' supplied. When `set_name` is given, a single set is plotted without faceting,
#' and the filename includes the set name.
#'
#' @param dds A DESeqDataSet object (rownames should already be gene symbols)
#' @param gene_sets Named list of character vectors.
#' @param condition_col Character column for grouping
#' @param level Foreground group
#' @param base Background group
#' @param plot_dir Output directory
#' @param show_points Logical: overlay individual sample‑level jittered points (default TRUE)
#' @param set_name Optional character string. If provided, `gene_sets` must be a
#'   list of length 1; this name is used in the output filename and plot title,
#'   and faceting is suppressed.
#' @return Invisibly returns the ggplot object
#' @export
plot_gene_zscore_individual <- function(dds, gene_sets, condition_col, level, base,
                                        plot_dir, show_points = TRUE,
                                        set_name = NULL) {
  if (is.character(gene_sets)) gene_sets <- list("Genes" = gene_sets)
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("`gene_sets` must be a named list, e.g. list(SetA = c('Gene1','Gene2'), SetB = c('Gene3')).")
  }

  if (!is.null(set_name)) {
    if (length(gene_sets) != 1) {
      warning("set_name provided but gene_sets has length != 1. Using only the first set.")
      gene_sets <- gene_sets[1]
    }
  }

  norm_counts      <- DESeq2::counts(dds, normalized = TRUE)
  meta_data        <- as.data.frame(SummarizedExperiment::colData(dds))
  meta_data$sample <- rownames(meta_data)
  valid_samples    <- meta_data$sample[meta_data[[condition_col]] %in% c(base, level)]

  gene_lookup <- unlist(lapply(names(gene_sets), function(nm) {
    stats::setNames(rep(nm, length(gene_sets[[nm]])), gene_sets[[nm]])
  }))
  gene_lookup <- gene_lookup[!duplicated(names(gene_lookup))]

  all_genes     <- names(gene_lookup)
  present_genes <- intersect(all_genes, rownames(norm_counts))
  missing_genes <- setdiff(all_genes, rownames(norm_counts))
  if (length(present_genes) == 0) stop("None of the specified genes were found in the data.")
  if (length(missing_genes) > 0) {
    message("   -> Note: genes not found (check spelling/case): ", paste(missing_genes, collapse = ", "))
  }

  expr_mat              <- norm_counts[present_genes, valid_samples, drop = FALSE]
  z_mat                 <- .zscore_matrix(expr_mat)

  df_long      <- as.data.frame(z_mat)
  df_long$gene <- rownames(df_long)
  plot_df      <- tidyr::pivot_longer(df_long, cols = -gene, names_to = "sample", values_to = "z_score")
  plot_df      <- merge(plot_df, meta_data[, c("sample", condition_col)], by = "sample")
  plot_df$Condition <- factor(plot_df[[condition_col]], levels = c(base, level))

  gene_order       <- present_genes[order(match(gene_lookup[present_genes], names(gene_sets)))]
  plot_df$gene_set <- factor(gene_lookup[plot_df$gene], levels = names(gene_sets))
  plot_df$gene     <- factor(plot_df$gene, levels = gene_order)

  summary_df <- do.call(rbind, lapply(split(plot_df, list(plot_df$gene, plot_df$Condition)), function(d) {
    if (nrow(d) == 0) return(NULL)
    data.frame(
      gene      = d$gene[1],
      gene_set  = d$gene_set[1],
      Condition = d$Condition[1],
      mean_z    = mean(d$z_score),
      sem       = if (nrow(d) > 1) stats::sd(d$z_score) / sqrt(nrow(d)) else 0
    )
  }))

  comp_name  <- paste0(level, "_vs_", base)
  dodge_pos  <- ggplot2::position_dodge(width = 0.75)

  p <- ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = summary_df,
      ggplot2::aes(x = gene, y = mean_z, fill = Condition),
      stat = "identity", position = dodge_pos, alpha = 0.75, width = 0.7, color = "black", linewidth = 0.3
    ) +
    ggplot2::geom_errorbar(
      data = summary_df,
      ggplot2::aes(x = gene, ymin = mean_z - sem, ymax = mean_z + sem, group = Condition),
      position = dodge_pos, width = 0.25
    )

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = gene, y = z_score, group = Condition),
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      size = 1.3, alpha = 0.6, color = "black"
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = stats::setNames(c("dodgerblue4", "red3"), c(base, level))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text      = ggplot2::element_text(face = "bold", size = 11),
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title    = ggplot2::element_blank()
    )

  if (!is.null(set_name)) {
    p <- p + ggplot2::labs(
      title = paste0("Individual Gene Z‑scores: ", set_name, " (", level, " vs ", base, ")"),
      x = "Features", y = "Z‑score (\u00b1 SEM)"
    )
    pdf_path <- file.path(plot_dir, paste0("Gene_Zscore_Individual_", comp_name, "_", set_name, ".pdf"))
  } else {
    n_sets <- length(unique(plot_df$gene_set))
    if (n_sets > 1) {
      p <- p + ggplot2::facet_wrap(~gene_set, scales = "free_x", nrow = 1)
    }
    p <- p + ggplot2::labs(
      title = "Z‑score by Condition, Individual Genes",
      x = "Features", y = "Z‑score (\u00b1 SEM)"
    )
    pdf_path <- file.path(plot_dir, paste0("Gene_Zscore_Individual_", comp_name, ".pdf"))
  }

  plot_width <- max(6, 0.9 * length(present_genes) + 2)
  ggplot2::ggsave(filename = pdf_path, plot = p,
                  width = plot_width, height = 5, device = .pdf_device(), limitsize = FALSE)
  message("   -> Individual gene Z‑score plot saved to: ", pdf_path)
  return(invisible(p))
}

#' SPIA two-way evidence plot (ggplot2 + ggrepel)
#' @param x A SPIA results data.frame
#' @param threshold Threshold for significance
#' @export
plotP_fork <- function(x, threshold = 0.01) {
  if (!inherits(x, "data.frame") | dim(x)[1] < 1 |
      !all(c("ID", "pNDE", "pPERT", "pG", "pGFdr", "pGFWER") %in% names(x))) {
    stop("SPIA graph can be applied only to a dataframe produced by SPIA function")
  }
  if (threshold < x[1, "pGFdr"]) {
    message("The threshold value was corrected to be equal to ", x[1, "pGFdr"])
    threshold <- x[1, "pGFdr"]
  }
  df  <- x
  pb  <- df$pPERT
  ph  <- df$pNDE

  combinemethod <- ifelse(
    sum(.combfunc(pb, ph, "fisher") == df$pG) > sum(.combfunc(pb, ph, "norminv") == df$pG),
    "fisher", "norminv"
  )
  okx <- (ph < 1e-6)
  oky <- (pb < 1e-6)
  ph[ph < 1e-6] <- 1e-6
  pb[pb < 1e-6] <- 1e-6
  df$x_val  <- -log(ph)
  df$y_val  <- -log(pb)
  df$Group  <- "Not Significant"
  df$Group[df$pGFdr  <= threshold] <- "FDR"
  df$Group[df$pGFWER <= threshold] <- "FWER"
  df$Group  <- factor(df$Group, levels = c("Not Significant", "FDR", "FWER"))
  df$Label  <- ""
  sig_idx   <- df$Group %in% c("FDR", "FWER")
  if (any(sig_idx)) df$Label[sig_idx] <- as.character(df$ID[sig_idx])

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_val, y = y_val)) +
    ggplot2::geom_point(ggplot2::aes(color = Group), size = 2.5) +
    ggplot2::scale_color_manual(values = c("Not Significant" = "black", "FDR" = "blue", "FWER" = "red"))

  tr_red <- threshold / nrow(na.omit(x))
  if (combinemethod == "fisher") {
    val_red  <- -log(.getP2(tr_red, "fisher")^2)
    if (is.finite(val_red)) {
      line_red <- data.frame(x = c(0, val_red), y = c(val_red, 0))
      p <- p + ggplot2::geom_path(data = line_red, ggplot2::aes(x = x, y = y), color = "red", linewidth = 1)
    } else {
      message("   -> plotP_fork: could not compute a finite FDR threshold value; skipping that line.")
    }
  } else {
    somep1  <- exp(seq(from = min(log(ph), na.rm = TRUE), to = max(log(ph), na.rm = TRUE), length = 200))
    somep2  <- pnorm(qnorm(tr_red) * sqrt(2) - qnorm(somep1))
    line_df <- data.frame(x = -log(somep1), y = -log(somep2))
    line_df <- line_df[is.finite(line_df$x) & is.finite(line_df$y), ]
    if (nrow(line_df) > 0) {
      p <- p + ggplot2::geom_line(data = line_df, ggplot2::aes(x = x, y = y), color = "red", linewidth = 1)
    } else {
      message("   -> plotP_fork: could not compute a finite FDR threshold curve; skipping that line.")
    }
  }

  tr_blue_old <- tr_red
  tr_blue     <- suppressWarnings(max(df$pG[df$pGFdr <= threshold], na.rm = TRUE))
  if (is.infinite(tr_blue) || tr_blue <= tr_blue_old) tr_blue <- tr_blue_old * 1.03
  if (combinemethod == "fisher") {
    val_blue  <- -log(.getP2(tr_blue, "fisher")^2)
    if (is.finite(val_blue)) {
      line_blue <- data.frame(x = c(0, val_blue), y = c(val_blue, 0))
      p <- p + ggplot2::geom_path(data = line_blue, ggplot2::aes(x = x, y = y), color = "blue", linewidth = 1)
    } else {
      message("   -> plotP_fork: could not compute a finite FWER threshold value; skipping that line.")
    }
  } else {
    somep1  <- exp(seq(from = min(log(ph), na.rm = TRUE), to = max(log(ph), na.rm = TRUE), length = 200))
    somep2  <- pnorm(qnorm(tr_blue) * sqrt(2) - qnorm(somep1))
    line_df <- data.frame(x = -log(somep1), y = -log(somep2))
    line_df <- line_df[is.finite(line_df$x) & is.finite(line_df$y), ]
    if (nrow(line_df) > 0) {
      p <- p + ggplot2::geom_line(data = line_df, ggplot2::aes(x = x, y = y), color = "blue", linewidth = 1)
    } else {
      message("   -> plotP_fork: could not compute a finite FWER threshold curve; skipping that line.")
    }
  }

  p <- p + ggrepel::geom_text_repel(
    ggplot2::aes(label = Label, color = Group),
    size = 3.5, box.padding = 0.6, max.overlaps = Inf, show.legend = FALSE
  )
  if (any(okx)) {
    p <- p + ggplot2::geom_text(data = df[okx, ], ggplot2::aes(x = x_val - 0.15, y = y_val),
                                label = "|", size = 4, color = "black")
  }
  if (any(oky)) {
    p <- p + ggplot2::geom_text(data = df[oky, ], ggplot2::aes(x = x_val, y = y_val - 0.15),
                                label = "_", size = 4, color = "black", vjust = 1)
  }
  max_val <- max(c(df$x_val, df$y_val) + 1, na.rm = TRUE)
  p <- p +
    ggplot2::coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    ggplot2::labs(title = "SPIA two-way evidence plot",
                  x = "-log(P NDE)", y = "-log(P PERT)", color = "Significance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  print(p)
  return(invisible(p))
}
