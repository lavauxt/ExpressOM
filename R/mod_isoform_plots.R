# mod_isoform_plots.R -- isoform-level plotting helpers (switch summary, sashimi,
# exon usage) shared by run_isoform_switch() and generate_dte_dtu_report()


.resolve_gene_id <- function(isoform_features, gene) {
  if (is.null(gene) || !nzchar(gene)) return(NULL)
  if (gene %in% isoform_features$gene_id) return(gene)

  if ("gene_name" %in% colnames(isoform_features)) {
    idx <- which(toupper(isoform_features$gene_name) == toupper(gene))
    if (length(idx) > 0) return(isoform_features$gene_id[idx[1]])
  }

  NULL
}

.get_exon_bins <- function(switch_list, gene_id) {
  if (is.null(switch_list$exons)) {
    stop("switch_list has no 'exons' slot (unexpected switchAnalyzeRlist structure).")
  }

  exon_mcols <- S4Vectors::mcols(switch_list$exons)

  if (!"isoform_id" %in% colnames(exon_mcols)) {
    stop("switch_list$exons has no 'isoform_id' metadata column.")
  }

  feat <- switch_list$isoformFeatures
  iso_ids <- unique(feat$isoform_id[feat$gene_id == gene_id])

  if (length(iso_ids) == 0) return(NULL)

  exon_gr <- switch_list$exons[exon_mcols$isoform_id %in% iso_ids]

  if (length(exon_gr) == 0) return(NULL)

  bins <- GenomicRanges::disjoin(exon_gr, ignore.strand = TRUE)
  bins <- sort(bins)

  S4Vectors::mcols(bins)$bin_id <- seq_along(bins)

  list(exon_gr = exon_gr, bins = bins, iso_ids = iso_ids)
}

#' @export
plot_isoform_switch_summary <- function(switch_list,
                                        plot_dir,
                                        genes_of_interest = NULL,
                                        level = NULL,
                                        base = NULL,
                                        top_n = 10,
                                        alpha = 0.05,
                                        dIFcutoff = 0.1,
                                        plot_topology = TRUE) {

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

  top_dir <- file.path(plot_dir, "top_switches")

  auto_ok <- tryCatch(
    {
      run_switch_plots <- function() {
        IsoformSwitchAnalyzeR::switchPlotTopSwitches(
          switchAnalyzeRlist = switch_list,
          n = top_n,
          filterForConsequences = FALSE,
          splitComparison = FALSE,
          splitFunctionalConsequences = FALSE,
          sortByQvals = TRUE,
          alpha = alpha,
          dIFcutoff = dIFcutoff,
          pathToOutput = top_dir,
          fileType = "pdf"
        )
      }

      if (isTRUE(plot_topology)) run_switch_plots() else suppressMessages(run_switch_plots())

      TRUE
    },
    error = function(e) {
      message("   -> switchPlotTopSwitches() failed: ", e$message)
      FALSE
    }
  )

  if (isTRUE(auto_ok) && dir.exists(top_dir)) {
    found <- list.files(top_dir, pattern = "\\.pdf$", full.names = TRUE, recursive = TRUE)
    generated <- c(generated, found)

    message(
      "   -> Top ", top_n, " isoform switch plots saved to: ", top_dir,
      " (", length(found), " file(s))"
    )
  }

  if (!is.null(genes_of_interest)) {
    cond1 <- if (!is.null(base)) base else feat$condition_1[1]
    cond2 <- if (!is.null(level)) level else feat$condition_2[1]

    for (gene_sym in genes_of_interest) {
      gene_id <- .resolve_gene_id(feat, gene_sym)

      if (is.null(gene_id)) {
        message("   -> Gene '", gene_sym, "' not found in switch_list. Skipping switch plot.")
        next
      }

      pdf_path <- file.path(plot_dir, paste0("SwitchPlot_", gene_sym, ".pdf"))

      n_domain_labels <- tryCatch(
        {
          da <- switch_list$domainAnalysis
          if (is.null(da) || !all(c("isoform_id", "hmm_name") %in% colnames(da))) {
            0L
          } else {
            iso_ids <- feat$isoform_id[feat$gene_id == gene_id]
            da_gene <- da[da$isoform_id %in% iso_ids & !is.na(da$hmm_name), ]

            labels <- if ("domain_isotype_simple" %in% colnames(da)) {
              paste0(
                da_gene$hmm_name,
                ifelse(
                  is.na(da_gene$domain_isotype_simple) | da_gene$domain_isotype_simple == "Reference",
                  "",
                  " (Non-ref Isotype)"
                )
              )
            } else {
              da_gene$hmm_name
            }

            length(unique(labels))
          }
        },
        error = function(e) 0L
      )

      # Baseline matches the previous fixed 11x9; grows once there are more
      # labels than comfortably fit at that size (empirically ~10-12).
      plot_height <- 9 + 0.22 * max(0, n_domain_labels - 10)
      legend_theme <- ggplot2::theme_bw(base_size = 8) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = 5),
          legend.title = ggplot2::element_text(size = 6),
          legend.key.size = ggplot2::unit(0.3, "lines"),
          legend.spacing.y = ggplot2::unit(0.1, "lines")
        )

      ok <- tryCatch(
        {
          .pdf_device()(pdf_path, onefile = FALSE, width = 11, height = plot_height)

          IsoformSwitchAnalyzeR::switchPlot(
            switch_list,
            gene = gene_id,
            condition1 = cond1,
            condition2 = cond2,
            localTheme = legend_theme,
            plotTopology = plot_topology
          )

          dev.off()
          TRUE
        },
        error = function(e) {
          if (grDevices::dev.cur() > 1) grDevices::dev.off()
          message("   -> switchPlot() failed for ", gene_sym, ": ", e$message)
          FALSE
        }
      )

      if (ok) {
        generated <- c(generated, pdf_path)
        message("   -> Isoform switch plot saved to: ", pdf_path)
      }
    }
  }

  invisible(generated)
}

#' @export
plot_isoform_sashimi <- function(switch_list, gene, level, base, plot_dir, min_if = 0.01) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    message("GenomicRanges/S4Vectors are required for sashimi-style plots. Skipping.")
    return(invisible(NULL))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting")

  feat_all <- switch_list$isoformFeatures
  gene_id <- .resolve_gene_id(feat_all, gene)

  if (is.null(gene_id)) {
    message("   -> Gene '", gene, "' not found in switch_list. Skipping sashimi plot.")
    return(invisible(NULL))
  }

  feat <- feat_all[feat_all$gene_id == gene_id &
                     feat_all$condition_1 == base &
                     feat_all$condition_2 == level, ]

  if (nrow(feat) == 0) {
    message(
      "   -> No isoformFeatures rows for gene '", gene, "' with condition_1=", base,
      ", condition_2=", level, ". Skipping sashimi plot."
    )
    return(invisible(NULL))
  }

  ginfo <- tryCatch(
    .get_exon_bins(switch_list, gene_id),
    error = function(e) {
      message("   -> Could not build exon structure for '", gene, "': ", e$message)
      NULL
    }
  )

  if (is.null(ginfo)) {
    message("   -> No exon structure found for '", gene, "'. Skipping sashimi plot.")
    return(invisible(NULL))
  }

  exon_gr <- ginfo$exon_gr
  exon_mcols <- S4Vectors::mcols(exon_gr)

  exon_df <- data.frame(
    isoform_id = exon_mcols$isoform_id,
    start = GenomicRanges::start(exon_gr),
    end = GenomicRanges::end(exon_gr),
    stringsAsFactors = FALSE
  )

  junction_list <- lapply(split(exon_df, exon_df$isoform_id), function(d) {
    d <- d[order(d$start), ]
    if (nrow(d) < 2) return(NULL)

    data.frame(
      isoform_id = d$isoform_id[1],
      jstart = d$end[-nrow(d)],
      jend = d$start[-1],
      stringsAsFactors = FALSE
    )
  })

  junctions <- do.call(rbind, junction_list)

  if (is.null(junctions) || nrow(junctions) == 0) {
    message("   -> Gene '", gene, "' has no multi-exon isoforms (no junctions to plot). Skipping.")
    return(invisible(NULL))
  }

  junctions <- merge(junctions, feat[, c("isoform_id", "IF1", "IF2")], by = "isoform_id")

  junc_summary <- stats::aggregate(
    cbind(IF1, IF2) ~ jstart + jend,
    data = junctions,
    FUN = function(x) sum(x, na.rm = TRUE)
  )

  colnames(junc_summary)[colnames(junc_summary) == "IF1"] <- "base_weight"
  colnames(junc_summary)[colnames(junc_summary) == "IF2"] <- "level_weight"

  junc_summary <- junc_summary[junc_summary$base_weight >= min_if |
                                 junc_summary$level_weight >= min_if, ]

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
        r$jstart,
        r$jend,
        height = r$level_weight / max_w,
        junction_id = r$junction_id,
        side = level,
        weight = r$level_weight
      )
    }

    if (r$base_weight >= min_if) {
      arc_rows[[length(arc_rows) + 1]] <- .make_arc(
        r$jstart,
        r$jend,
        height = -(r$base_weight / max_w),
        junction_id = r$junction_id,
        side = base,
        weight = r$base_weight
      )
    }
  }

  arc_df <- do.call(rbind, arc_rows)
  arc_df$side <- factor(arc_df$side, levels = c(base, level))

  bins_df <- as.data.frame(ginfo$bins)
  bins_df$bin_id <- S4Vectors::mcols(ginfo$bins)$bin_id

  chrom_label <- as.character(GenomicRanges::seqnames(ginfo$bins)[1])

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = bins_df,
      ggplot2::aes(xmin = start, xmax = end, ymin = -0.04, ymax = 0.04),
      fill = "grey35",
      color = NA
    ) +
    ggplot2::geom_path(
      data = arc_df,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(junction_id, side),
        linewidth = weight,
        color = side
      ),
      lineend = "round"
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(c("dodgerblue4", "red3"), c(base, level)),
      name = "Condition"
    ) +
    ggplot2::scale_linewidth(range = c(0.3, 3), name = "Usage (IF)") +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    ggplot2::labs(
      title = paste0("Isoform Structure & Junction Usage: ", gene),
      subtitle = paste0(
        level,
        " (above) vs ",
        base,
        " (below) — arc thickness ∝ isoform-fraction-weighted junction usage"
      ),
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

  ggplot2::ggsave(
    filename = pdf_path,
    plot = p,
    width = 11,
    height = 5,
    device = .pdf_device()
  )

  message("   -> Isoform sashimi-style junction plot saved to: ", pdf_path)

  invisible(pdf_path)
}

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
  gene_id <- .resolve_gene_id(feat_all, gene)

  if (is.null(gene_id)) {
    message("   -> Gene '", gene, "' not found in switch_list. Skipping exon usage plot.")
    return(invisible(NULL))
  }

  feat <- feat_all[feat_all$gene_id == gene_id &
                     feat_all$condition_1 == base &
                     feat_all$condition_2 == level, ]

  if (nrow(feat) == 0) {
    message(
      "   -> No isoformFeatures rows for gene '", gene, "' with condition_1=", base,
      ", condition_2=", level, ". Skipping exon usage plot."
    )
    return(invisible(NULL))
  }

  ginfo <- tryCatch(
    .get_exon_bins(switch_list, gene_id),
    error = function(e) {
      message("   -> Could not build exon structure for '", gene, "': ", e$message)
      NULL
    }
  )

  if (is.null(ginfo)) {
    message("   -> No exon structure found for '", gene, "'. Skipping exon usage plot.")
    return(invisible(NULL))
  }

  exon_gr <- ginfo$exon_gr
  bins <- ginfo$bins

  ov <- GenomicRanges::findOverlaps(exon_gr, bins, ignore.strand = TRUE)

  map_df <- data.frame(
    isoform_id = S4Vectors::mcols(exon_gr)$isoform_id[S4Vectors::queryHits(ov)],
    bin_id = S4Vectors::mcols(bins)$bin_id[S4Vectors::subjectHits(ov)],
    stringsAsFactors = FALSE
  )

  map_df <- unique(map_df)
  map_df <- merge(map_df, feat[, c("isoform_id", "iso_value_1", "iso_value_2")], by = "isoform_id")

  if (nrow(map_df) == 0) {
    message("   -> No overlapping isoform/bin data for '", gene, "'. Skipping exon usage plot.")
    return(invisible(NULL))
  }

  bin_summary <- stats::aggregate(
    cbind(iso_value_1, iso_value_2) ~ bin_id,
    data = map_df,
    FUN = sum
  )

  colnames(bin_summary) <- c("bin_id", base, level)

  plot_df <- tidyr::pivot_longer(
    bin_summary,
    cols = -bin_id,
    names_to = "Condition",
    values_to = "expression"
  )

  plot_df$Condition <- factor(plot_df$Condition, levels = c(base, level))
  plot_df$bin_id <- factor(plot_df$bin_id, levels = sort(unique(plot_df$bin_id)))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = bin_id, y = expression, fill = Condition)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.6,
      color = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(c("dodgerblue4", "red3"), c(base, level))
    ) +
    ggplot2::labs(
      title = paste0("Exon-Bin Expression Comparison: ", gene),
      subtitle = "Non-overlapping exon bins; bars = summed isoform expression per bin",
      x = "Exon bin (5'→3', genomic order)",
      y = "Summed isoform expression"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "grey30"),
      legend.title = ggplot2::element_blank()
    )

  pdf_path <- file.path(plot_dir, paste0("ExonUsage_", gene, ".pdf"))

  ggplot2::ggsave(
    filename = pdf_path,
    plot = p,
    width = max(7, 0.35 * nrow(bin_summary) + 2),
    height = 5,
    device = .pdf_device()
  )

  message("   -> Exon-bin expression comparison plot saved to: ", pdf_path)

  invisible(pdf_path)
}

