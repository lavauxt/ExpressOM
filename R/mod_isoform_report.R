# mod_isoform_report.R -- DTE/DTU HTML report generation
#

#' Generate DTE/DTU HTML report
#' @export
generate_dte_dtu_report <- function(dte_results,
                                    dtu_results,
                                    isoform_obj,
                                    out_dir,
                                    condition,
                                    level,
                                    base,
                                    genes_of_interest = NULL,
                                    top_n = 15,
                                    switch_list = NULL,
                                    dexseq_results = NULL,
                                    plot_switch_summary = TRUE,
                                    switch_plot_top_n = 10,
                                    plot_sashimi = TRUE,
                                    plot_exon_usage = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting")
  if (!requireNamespace("DT", quietly = TRUE)) stop("DT is required for interactive tables")
  if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("rmarkdown is required to generate the report")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required for data manipulation")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required for data reshaping")

  report_dir <- file.path(out_dir, "IsoformSwitch", "DTU_DTE_report")
  if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

  dte <- dte_results

  for (col in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
    if (col %in% colnames(dte)) {
      dte[[col]] <- suppressWarnings(as.numeric(as.character(dte[[col]])))
    }
  }

  if ("signif" %in% colnames(dte)) {
    dte$signif <- as.logical(dte$signif)
    dte$signif[is.na(dte$signif)] <- FALSE
  } else {
    dte$signif <- FALSE
  }

  dte$transcript_id <- clean_transcript_id(dte$transcript_id)

  dte$gene_label <- .format_gene_tx_label(
    .coalesce_gene_label(dte$gene_symbol, dte$transcript_id),
    dte$transcript_id
  )

  dte$signif_label <- ifelse(dte$signif, "Significant", "Not significant")

  dtu <- dtu_results$dtu_results

  for (col in c("pvalue", "adj_pvalue")) {
    if (col %in% colnames(dtu)) {
      dtu[[col]] <- suppressWarnings(as.numeric(as.character(dtu[[col]])))
    }
  }

  gene_map <- isoform_obj$gene_map[, c("ensembl", "symbol")]
  colnames(gene_map) <- c("gene_id", "gene_symbol")

  if (!"gene_symbol" %in% colnames(dtu)) {
    dtu <- merge(dtu, gene_map, by = "gene_id", all.x = TRUE)
  }

  missing_sym <- is.na(dtu$gene_symbol) | dtu$gene_symbol == ""

  if (any(missing_sym)) {
    already_a_symbol <- dtu$gene_id[missing_sym] %in% gene_map$gene_symbol
    dtu$gene_symbol[missing_sym][already_a_symbol] <- dtu$gene_id[missing_sym][already_a_symbol]
  }

  dtu$gene_label <- .format_gene_tx_label(
    .coalesce_gene_label(dtu$gene_symbol, dtu$gene_id),
    dtu$gene_id
  )

  write.csv(dte, file.path(report_dir, "DTE_results_annotated.csv"), row.names = FALSE)
  write.csv(dtu, file.path(report_dir, "DTU_results_annotated.csv"), row.names = FALSE)

  has_dexseq <- !is.null(dexseq_results) && !is.null(dexseq_results$results_df)

  if (has_dexseq) {
    write.csv(
      dexseq_results$results_df,
      file.path(report_dir, "DEXSeq_results_annotated.csv"),
      row.names = FALSE
    )
  }

  plot_dir <- file.path(report_dir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir)

  dte_sig <- dte[!is.na(dte$padj), ]

  p_volcano <- ggplot2::ggplot(
    dte_sig,
    ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = signif_label)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "red")) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggplot2::labs(
      title = paste("DTE Volcano plot:", level, "vs", base),
      x = "log2 Fold Change",
      y = "-log10(adj. p-value)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    dte_to_label <- head(
      dte_sig[dte_sig$signif, ][order(dte_sig[dte_sig$signif, ]$padj), ],
      20
    )

    if (nrow(dte_to_label) > 0) {
      p_volcano <- p_volcano + ggrepel::geom_text_repel(
        data = dte_to_label,
        mapping = ggplot2::aes(x = log2FoldChange, y = -log10(padj), label = gene_label),
        inherit.aes = FALSE,
        size = 2.8,
        max.overlaps = 30,
        segment.size = 0.3,
        color = "black",
        seed = 42
      )
    }
  } else {
    message("   -> ggrepel not installed; DTE_volcano.pdf will show points without gene labels.")
  }

  .safe_ggsave(filename = file.path(plot_dir, "DTE_volcano.pdf"), plot = p_volcano, width = 8, height = 6)

  dte_ma <- dte[!is.na(dte$baseMean) & dte$baseMean > 0, ]

  p_ma <- ggplot2::ggplot(
    dte_ma,
    ggplot2::aes(x = baseMean, y = log2FoldChange, color = signif_label)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "red")) +
    ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    ggplot2::labs(
      title = "DTE MA plot",
      x = "Mean of normalized counts",
      y = "log2 Fold Change"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  .safe_ggsave(filename = file.path(plot_dir, "DTE_MA.pdf"), plot = p_ma, width = 8, height = 6)

  dte_filtered <- dte[
    !is.na(dte$signif) &
      dte$signif &
      is.finite(dte$log2FoldChange) &
      is.finite(dte$padj),
  ]

  if (nrow(dte_filtered) > 0) {
    dte_ordered <- dte_filtered[order(dte_filtered$padj), ]
    top_dte <- head(dte_ordered, top_n)

    p_bar <- ggplot2::ggplot(
      top_dte,
      ggplot2::aes(x = reorder(gene_label, abs(log2FoldChange)), y = log2FoldChange)
    ) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste("Top", top_n, "significant DE transcripts"),
        x = "",
        y = "log2 Fold Change"
      ) +
      ggplot2::theme_minimal()

    .safe_ggsave(
      filename = file.path(plot_dir, "DTE_top_barplot.pdf"),
      plot = p_bar,
      width = 10,
      height = max(6, top_n * 0.3)
    )
  } else {
    message("No significant DTE transcripts found; skipping top barplot.")
  }

  p_hist <- ggplot2::ggplot(dtu, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
    ggplot2::labs(title = "DTU p-value distribution", x = "p-value", y = "Count") +
    ggplot2::theme_minimal()

  .safe_ggsave(filename = file.path(plot_dir, "DTU_pvalue_hist.pdf"), plot = p_hist, width = 6, height = 4)

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

    p_dtu_bar <- ggplot2::ggplot(
      top_dtu,
      ggplot2::aes(x = reorder(dtu_label, neg_log10_padj), y = neg_log10_padj)
    ) +
      ggplot2::geom_col(fill = "darkgreen") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste("Top", top_n, "significant DTU transcripts"),
        x = "",
        y = "-log10(adj. p-value)"
      ) +
      ggplot2::theme_minimal()

    .safe_ggsave(
      filename = file.path(plot_dir, "DTU_top_barplot.pdf"),
      plot = p_dtu_bar,
      width = 10,
      height = max(6, top_n * 0.3)
    )
  } else {
    message("No significant DTU transcripts found (adj_pvalue < 0.05); skipping DTU top barplot.")
  }

  safe_run(
    {
      sig_sets <- list(
        DTE = unique(dte$gene_id[!is.na(dte$gene_id) & dte$signif]),
        DTU = unique(dtu$gene_id[!is.na(dtu$gene_id) & !is.na(dtu$adj_pvalue) & dtu$adj_pvalue < 0.05])
      )

      if (has_dexseq && "gene_qvalue" %in% colnames(dexseq_results$results_df)) {
        dex <- dexseq_results$results_df
        sig_sets$DEXSeq <- unique(
          dex$gene_id[!is.na(dex$gene_id) & !is.na(dex$gene_qvalue) & dex$gene_qvalue < 0.05]
        )
      }

      sig_sets <- sig_sets[lengths(sig_sets) > 0]

      if (length(sig_sets) >= 2) {
        all_genes <- unique(unlist(sig_sets, use.names = FALSE))

        membership <- vapply(
          sig_sets,
          function(s) all_genes %in% s,
          logical(length(all_genes))
        )

        combo <- apply(
          membership,
          1,
          function(row) paste(names(sig_sets)[row], collapse = " + ")
        )

        combo_counts <- as.data.frame(table(combination = combo), stringsAsFactors = FALSE)
        colnames(combo_counts) <- c("combination", "n_genes")
        combo_counts <- combo_counts[order(-combo_counts$n_genes), ]

        combo_counts$combination <- factor(combo_counts$combination, levels = combo_counts$combination)

        p_concordance <- ggplot2::ggplot(
          combo_counts,
          ggplot2::aes(x = combination, y = n_genes)
        ) +
          ggplot2::geom_col(fill = "darkslateblue") +
          ggplot2::geom_text(ggplot2::aes(label = n_genes), vjust = -0.4, size = 3.5) +
          ggplot2::labs(
            title = "Significant gene concordance across engines",
            subtitle = paste(
              sprintf("%s: %d genes", names(sig_sets), lengths(sig_sets)),
              collapse = "   |   "
            ),
            x = NULL,
            y = "Number of genes"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

        .safe_ggsave(
          filename = file.path(plot_dir, "Engine_concordance.pdf"),
          plot = p_concordance,
          width = max(7, 1.3 * nrow(combo_counts)),
          height = 6
        )

        message(
          "   -> Engine concordance (significant genes): ",
          paste(sprintf("%s=%d", names(sig_sets), lengths(sig_sets)), collapse = ", ")
        )
      } else {
        message("   -> Skipping engine concordance plot: fewer than 2 engines have any significant genes to compare.")
      }
    },
    label = "Cross-engine (DTE/DTU/DEXSeq) concordance plot"
  )

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

      prop_df <- as.data.frame(t(tx_counts))
      prop_df <- prop_df / rowSums(prop_df)
      prop_df$sample <- rownames(prop_df)

      prop_long <- tidyr::pivot_longer(
        prop_df,
        cols = -sample,
        names_to = "transcript",
        values_to = "proportion"
      )

      prop_long <- merge(prop_long, meta, by.x = "sample", by.y = 0)

      p_prop <- ggplot2::ggplot(
        prop_long,
        ggplot2::aes(x = sample, y = proportion, fill = transcript)
      ) +
        ggplot2::geom_bar(stat = "identity", position = "fill") +
        ggplot2::facet_grid(
          stats::as.formula(paste("~", condition)),
          scales = "free_x",
          space = "free_x"
        ) +
        ggplot2::labs(
          title = paste("Transcript proportions for", gene_sym),
          y = "Proportion",
          x = "Sample"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      .safe_ggsave(
        filename = file.path(plot_dir, paste0("proportions_", gene_sym, ".pdf")),
        plot = p_prop,
        width = 10,
        height = 6
      )
    }
  }

  has_switch <- !is.null(switch_list) &&
    !is.null(switch_list$isoformFeatures) &&
    nrow(switch_list$isoformFeatures) > 0

  if (has_switch) {
    feat_sw <- switch_list$isoformFeatures

    if (all(c("dIF", "isoform_switch_q_value") %in% colnames(feat_sw))) {
      feat_sw_plot <- feat_sw

      feat_sw_plot$dIF <- suppressWarnings(as.numeric(as.character(feat_sw_plot$dIF)))
      feat_sw_plot$isoform_switch_q_value <- suppressWarnings(
        as.numeric(as.character(feat_sw_plot$isoform_switch_q_value))
      )

      feat_sw_plot <- feat_sw_plot[
        is.finite(feat_sw_plot$dIF) &
          is.finite(feat_sw_plot$isoform_switch_q_value),
      ]

      if (nrow(feat_sw_plot) > 0) {
        p_switch_overview <- ggplot2::ggplot(
          feat_sw_plot,
          ggplot2::aes(x = dIF, y = -log10(isoform_switch_q_value))
        ) +
          ggplot2::geom_point(
            ggplot2::aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
            size = 1,
            alpha = 0.6
          ) +
          ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
          ggplot2::geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
          ggplot2::scale_color_manual(
            "Significant\nIsoform Switch",
            values = c("grey60", "red3")
          ) +
          ggplot2::labs(
            title = "Isoform Switch Overview",
            x = "dIF (Isoform Fraction change)",
            y = "-log10(Isoform Switch q-value)"
          ) +
          ggplot2::theme_minimal()

        if (requireNamespace("ggrepel", quietly = TRUE)) {
          gene_col <- if ("gene_name" %in% colnames(feat_sw_plot)) "gene_name" else "gene_id"

          feat_sw_plot$.label <- .coalesce_gene_label(
            feat_sw_plot[[gene_col]],
            feat_sw_plot$gene_id
          )

          feat_sig <- feat_sw_plot[
            !is.na(feat_sw_plot$isoform_switch_q_value) &
              feat_sw_plot$isoform_switch_q_value < 0.05 &
              abs(feat_sw_plot$dIF) > 0.1,
          ]

          feat_to_label <- head(feat_sig[order(feat_sig$isoform_switch_q_value), ], 20)

          if (nrow(feat_to_label) > 0) {
            p_switch_overview <- p_switch_overview + ggrepel::geom_text_repel(
              data = feat_to_label,
              mapping = ggplot2::aes(
                x = dIF,
                y = -log10(isoform_switch_q_value),
                label = .label
              ),
              inherit.aes = FALSE,
              size = 2.8,
              max.overlaps = 30,
              segment.size = 0.3,
              color = "black",
              seed = 42
            )
          }
        } else {
          message("   -> ggrepel not installed; IsoformSwitch_overview.pdf will show points without gene labels.")
        }

        .safe_ggsave(
          filename = file.path(plot_dir, "IsoformSwitch_overview.pdf"),
          plot = p_switch_overview,
          width = 7,
          height = 6
        )
      }
    }

    .try_isa_plot <- function(pdf_path, width, height, isa_fn, isa_args, label) {
      .has_arg <- function(fn, nm) nm %in% names(formals(fn))

      skip_file <- sub("\\.pdf$", ".SKIP", pdf_path)

      check_args <- isa_args

      if (.has_arg(isa_fn, "plot")) check_args$plot <- FALSE
      if (.has_arg(isa_fn, "returnResult")) check_args$returnResult <- TRUE

      data_result <- tryCatch(
        do.call(isa_fn, check_args),
        error = function(e) {
          message("   -> ", label, ": data check failed: ", conditionMessage(e))
          NULL
        }
      )

      has_data <- !is.null(data_result) &&
        (
          (is.data.frame(data_result) && nrow(data_result) > 0) ||
            (is.list(data_result) && length(data_result) > 0) ||
            (is.atomic(data_result) && length(data_result) > 0)
        )

      if (!has_data) {
        message(
          "   -> Skipping ", basename(pdf_path),
          ": no data to plot (", label, " returned empty result)"
        )

        writeLines(
          c(
            paste("# Plot skipped:", label),
            "# Reason: no data available to summarize",
            paste("# Generated:", format(Sys.time()))
          ),
          skip_file
        )

        return(invisible(FALSE))
      }

      plot_args <- isa_args

      if (.has_arg(isa_fn, "plot")) plot_args$plot <- TRUE
      if (.has_arg(isa_fn, "returnResult")) plot_args$returnResult <- FALSE

      ok <- tryCatch(
        {
          grDevices::cairo_pdf(pdf_path, width = width, height = height)

          res <- do.call(isa_fn, plot_args)

          if (inherits(res, "ggplot")) {
            print(res)
          } else if (inherits(res, "grob")) {
            grid::grid.draw(res)
          } else if (inherits(res, c("arrangeGrob", "patchwork"))) {
            print(res)
          }

          dev.off()
          TRUE
        },
        error = function(e) {
          if (grDevices::dev.cur() > 1) dev.off()
          message("   -> ", label, " plotting failed: ", conditionMessage(e))
          FALSE
        }
      )

      sz <- if (file.exists(pdf_path)) file.info(pdf_path)$size else NA_real_

      if (ok && !is.na(sz) && sz >= 2000) {
        message("   -> Saved: ", pdf_path)
        return(invisible(TRUE))
      }

      if (file.exists(pdf_path)) file.remove(pdf_path)

      writeLines(
        c(
          paste("# Plot skipped:", label),
          "# Reason: PDF was empty/corrupt after plotting attempt",
          paste("# Generated:", format(Sys.time()))
        ),
        skip_file
      )

      invisible(FALSE)
    }

    .try_isa_plot(
      pdf_path = file.path(plot_dir, "ConsequenceSummary.pdf"),
      width = 9,
      height = 6,
      isa_fn = IsoformSwitchAnalyzeR::extractConsequenceSummary,
      isa_args = list(
        switchAnalyzeRlist = switch_list,
        consequencesToAnalyze = "all"
      ),
      label = "extractConsequenceSummary"
    )

    .try_isa_plot(
      pdf_path = file.path(plot_dir, "ConsequenceEnrichment.pdf"),
      width = 9,
      height = 6,
      isa_fn = IsoformSwitchAnalyzeR::extractConsequenceEnrichment,
      isa_args = list(
        switchAnalyzeRlist = switch_list,
        consequencesToAnalyze = "all"
      ),
      label = "extractConsequenceEnrichment"
    )

    if (!is.null(switch_list[["AlternativeSplicingAnalysis"]])) {
      .try_isa_plot(
        pdf_path = file.path(plot_dir, "SplicingSummary.pdf"),
        width = 9,
        height = 6,
        isa_fn = IsoformSwitchAnalyzeR::extractSplicingSummary,
        isa_args = list(
          switchAnalyzeRlist = switch_list,
          splicingToAnalyze = "all"
        ),
        label = "extractSplicingSummary"
      )

      .try_isa_plot(
        pdf_path = file.path(plot_dir, "SplicingEnrichment.pdf"),
        width = 9,
        height = 6,
        isa_fn = IsoformSwitchAnalyzeR::extractSplicingEnrichment,
        isa_args = list(
          switchAnalyzeRlist = switch_list,
          splicingToAnalyze = "all"
        ),
        label = "extractSplicingEnrichment"
      )
    } else {
      message(
        "  Skipping splicing summary/enrichment plots: no AlternativeSplicingAnalysis found ",
        "in switch_list (analyzeAlternativeSplicing() may have failed upstream)."
      )

      for (fn in c("SplicingSummary", "SplicingEnrichment")) {
        writeLines(
          c(
            paste("# Plot skipped:", fn),
            "# Reason: no AlternativeSplicingAnalysis in switch_list",
            paste("# Generated:", format(Sys.time()))
          ),
          file.path(plot_dir, paste0(fn, ".pdf.SKIP"))
        )
      }
    }

    if (isTRUE(plot_switch_summary)) {
      safe_run(
        plot_isoform_switch_summary(
          switch_list,
          plot_dir = plot_dir,
          genes_of_interest = genes_of_interest,
          level = level,
          base = base,
          top_n = switch_plot_top_n
        ),
        label = "Isoform switch summary plots"
      )
    }

    if (!is.null(genes_of_interest)) {
      for (gene_sym in genes_of_interest) {
        if (isTRUE(plot_sashimi)) {
          safe_run(
            plot_isoform_sashimi(
              switch_list,
              gene = gene_sym,
              level = level,
              base = base,
              plot_dir = plot_dir
            ),
            label = paste0("Sashimi plot for ", gene_sym)
          )
        }

        if (isTRUE(plot_exon_usage)) {
          safe_run(
            plot_exon_usage_comparison(
              switch_list,
              gene = gene_sym,
              level = level,
              base = base,
              plot_dir = plot_dir
            ),
            label = paste0("Exon usage plot for ", gene_sym)
          )
        }
      }
    }
  }

  if (has_dexseq && !is.null(dexseq_results$dxr_list) && !is.null(genes_of_interest)) {
    for (gene_sym in genes_of_interest) {
      gene_ens <- isoform_obj$gene_map$ensembl[isoform_obj$gene_map$symbol == gene_sym]
      if (length(gene_ens) == 0) next

      safe_run(
        plot_dexseq_gene(
          dexseq_results$dxr_list,
          gene_id = gene_ens[1],
          plot_dir = plot_dir,
          gene_symbol = gene_sym
        ),
        label = paste0("DEXSeq plot for ", gene_sym)
      )
    }
  }

  all_plot_pdfs <- list.files(
    plot_dir,
    pattern = "\\.pdf$",
    recursive = TRUE,
    full.names = TRUE
  )

  for (f in all_plot_pdfs) {
    sz <- file.info(f)$size

    if (is.na(sz) || sz < 2000) {
      message("   -> Removing empty/corrupt PDF: ", f)

      file.remove(f)

      writeLines(
        c(
          "# Plot removed: empty/corrupt PDF",
          paste("# File:", basename(f)),
          paste("# Generated:", format(Sys.time()))
        ),
        sub("\\.pdf$", ".SKIP", f)
      )
    }
  }

  dte_abs_path <- normalizePath(
    file.path(report_dir, "DTE_results_annotated.csv"),
    winslash = "/"
  )

  dtu_abs_path <- normalizePath(
    file.path(report_dir, "DTU_results_annotated.csv"),
    winslash = "/"
  )

  dexseq_abs_path <- if (has_dexseq) {
    normalizePath(
      file.path(report_dir, "DEXSeq_results_annotated.csv"),
      winslash = "/"
    )
  } else {
    NULL
  }

  plot_pdfs <- list.files(
    plot_dir,
    pattern = "\\.pdf$",
    recursive = TRUE,
    full.names = TRUE
  )

  pdf_sizes <- vapply(
    plot_pdfs,
    function(f) file.info(f)$size,
    numeric(1)
  )

  valid_pdfs <- plot_pdfs[pdf_sizes >= 2000]

  if (length(valid_pdfs) < length(plot_pdfs)) {
    skipped <- setdiff(basename(plot_pdfs), basename(valid_pdfs))

    message(
      "   -> WARNING: ", length(skipped),
      " PDF(s) are <2KB (likely empty) and will NOT be converted to PNG: ",
      paste(skipped, collapse = ", ")
    )
  }

  manifest <- data.frame(
    pdf_path = plot_pdfs,
    size_bytes = pdf_sizes,
    likely_empty = pdf_sizes < 2000,
    png_converted = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(plot_pdfs)) {
    if (manifest$likely_empty[i]) next

    png_path <- convert_pdf_to_png(plot_pdfs[i])

    if (!is.null(png_path)) {
      manifest$png_converted[i] <- TRUE
    }
  }

  manifest_path <- file.path(report_dir, "plot_manifest.csv")
  write.csv(manifest, manifest_path, row.names = FALSE)

  n_empty <- sum(manifest$likely_empty)
  n_unconverted <- sum(!manifest$png_converted & !manifest$likely_empty)

  if (n_empty > 0) {
    message(
      "   -> WARNING: ", n_empty, " plot PDF(s) are likely empty (<2KB). ",
      "These will not appear in the HTML report. See ", manifest_path, "."
    )
  }

  if (n_unconverted > 0) {
    note_path <- file.path(report_dir, "README_IF_PLOTS_MISSING.txt")

    writeLines(
      c(
        sprintf(
          "%d of %d non-empty plot PDF(s) could not be converted to PNG.",
          n_unconverted,
          sum(!manifest$likely_empty)
        ),
        "DTU_DTE_report.html embeds figures as PNG because browsers cannot display PDF via <img>.",
        "Any figure that did not convert will not show in the HTML report, but its",
        "PDF exists in plots/ and can be opened directly.",
        "",
        "Recommended fix:",
        "  install.packages('pdftools')",
        "",
        "Optional fallback:",
        "  install.packages('magick')",
        "  plus Ghostscript on your system.",
        "",
        "See plot_manifest.csv for exact per-file status."
      ),
      note_path
    )

    warning(
      "Could not convert ", n_unconverted, " report plot(s) to PNG. See ",
      note_path, " and ", manifest_path, "."
    )
  }

  skip_files <- list.files(
    plot_dir,
    pattern = "\\.SKIP$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(skip_files) > 0) {
    skip_summary <- data.frame(
      plot_name = sub("\\.pdf\\.SKIP$", "", basename(skip_files)),
      reason = vapply(
        skip_files,
        function(f) {
          lines <- readLines(f, n = 2, warn = FALSE)
          if (length(lines) >= 2) lines[2] else "unknown"
        },
        character(1)
      ),
      stringsAsFactors = FALSE
    )

    write.csv(
      skip_summary,
      file.path(report_dir, "skipped_plots.csv"),
      row.names = FALSE
    )

    message(
      "   -> ", nrow(skip_summary),
      " plot(s) were intentionally skipped or removed. See skipped_plots.csv."
    )
  }

  template_path <- tryCatch(
    .expressom_rmd_path("dte_dtu_report.Rmd"),
    error = function(e) NULL
  )

  if (!is.null(template_path)) {
    rmd_file <- tempfile(fileext = ".Rmd")
    on.exit(unlink(rmd_file), add = TRUE)

    file.copy(template_path, rmd_file, overwrite = TRUE)
  } else {
    warning(
      "Template 'dte_dtu_report.Rmd' not found. Using minimal fallback report.",
      call. = FALSE
    )

    rmd_file <- file.path(report_dir, "report_fallback.Rmd")

    writeLines(
      c(
        "---",
        "title: 'DTE / DTU fallback report'",
        "output: html_document",
        "params:",
        "  level: ''",
        "  base: ''",
        "  top_n: 15",
        "  genes_of_interest: ''",
        "  dte_csv: ''",
        "  dtu_csv: ''",
        "  dexseq_csv: ''",
        "  has_dexseq: no",
        "  has_switch: no",
        "  plot_dir: ''",
        "  switch_plot_top_n: 10",
        "---",
        "",
        "```{r setup, include=FALSE}",
        "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
        "```",
        "",
        "# DTE results",
        "",
        "```{r}",
        "if (nzchar(params$dte_csv) && file.exists(params$dte_csv)) {",
        "  DT::datatable(utils::read.csv(params$dte_csv), filter = 'top')",
        "} else {",
        "  knitr::asis_output('DTE result file not found.')",
        "}",
        "```",
        "",
        "# DTU results",
        "",
        "```{r}",
        "if (nzchar(params$dtu_csv) && file.exists(params$dtu_csv)) {",
        "  DT::datatable(utils::read.csv(params$dtu_csv), filter = 'top')",
        "} else {",
        "  knitr::asis_output('DTU result file not found.')",
        "}",
        "```",
        "",
        "# DEXSeq results",
        "",
        "```{r}",
        "if (isTRUE(params$has_dexseq) && nzchar(params$dexseq_csv) && file.exists(params$dexseq_csv)) {",
        "  DT::datatable(utils::read.csv(params$dexseq_csv), filter = 'top')",
        "} else {",
        "  knitr::asis_output('No DEXSeq results available.')",
        "}",
        "```"
      ),
      rmd_file
    )
  }

  rmarkdown::render(
    rmd_file,
    output_file = "DTU_DTE_report.html",
    output_dir = report_dir,
    quiet = TRUE,
    knit_root_dir = normalizePath(report_dir, winslash = "/"),
    envir = new.env(parent = globalenv()),
    params = list(
      level = level,
      base = base,
      top_n = top_n,
      genes_of_interest = genes_of_interest,
      dte_csv = dte_abs_path,
      dtu_csv = dtu_abs_path,
      dexseq_csv = dexseq_abs_path,
      has_dexseq = has_dexseq,
      has_switch = has_switch,
      plot_dir = normalizePath(plot_dir, winslash = "/"),
      switch_plot_top_n = switch_plot_top_n
    )
  )

  gc()

  message("DTE/DTU report generated in: ", report_dir)

  invisible(NULL)
}
