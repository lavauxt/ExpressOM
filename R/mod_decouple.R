# mod_decouple.R -- decoupleR-based TF/pathway activity inference.
#
# Adds three "human resources" (or their mouse counterparts) to the pipeline:
#   * CollecTRI  - TF -> target regulons with sign (mor)
#   * PROGENy    - pathway footprints (source = pathway, target = gene, mor = weight)
#   * DoRothEA   - TF regulons filtered to confidence A/B/C
#
# The stat column of the DGE table (DESeq2 Wald stat, or signed -log10 padj)
# is used as the gene-level score. Results are written to
# <out_dir>/DecoupleR/<resource>/DecoupleR_<resource>_<comp>.csv and
# rendered as interactive plotly barplots in the unified report.
#
# NEW: For every source with p_value < padj_cutoff in a resource, a targeted
# volcano plot and a leading-edge (barcode) plot are also generated:
#   * regulon_<TF>.csv         -- matched targets + mor + padj + LFC + category
#   * leading_edge_<TF>.csv    -- rank + target + stat (sorted by rank)
#   * tf_volcano_<TF>.pdf      -- volcano limited to the TF's regulon
#   * tf_leading_edge_<TF>.pdf -- barcode plot over all genes ranked by stat
#   * tf_regulon_summary.csv   -- per-source counts (targets / up / down / ns)
# The interactive counterparts are exposed via iplot_tf_regulon_volcano() and
# iplot_tf_leading_edge() in mod_report_plots.R, and displayed per-resource /
# per-TF in the "TF regulon plots" tab of the unified report.

#' Fetch a decoupleR network for one of the supported resources.
#' @keywords internal
.decouple_get_net <- function(resource, organism) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required. Install with BiocManager::install('decoupleR').")
  }
  org <- if (grepl("mus|mouse", organism, ignore.case = TRUE)) "mouse" else "human"
  net <- switch(
    resource,
    collectri = {
      n <- decoupleR::get_collectri(organism = org, split_complexes = FALSE)
      n <- as.data.frame(n)
      if (!"mor" %in% colnames(n)) n$mor <- 1
      n[, c("source", "target", "mor")]
    },
    progeny = {
      n <- decoupleR::get_progeny(organism = org)
      n <- as.data.frame(n)
      data.frame(
        source = as.character(n$source),
        target = as.character(n$target),
        mor    = as.numeric(n$weight),
        stringsAsFactors = FALSE
      )
    },
    dorothea = {
      n <- decoupleR::get_dorothea(organism = org, levels = c("A", "B", "C"))
      n <- as.data.frame(n)
      if (!"mor" %in% colnames(n)) n$mor <- 1
      n$mor[n$mor == 0] <- 1
      n[, c("source", "target", "mor")]
    },
    stop("Unknown decoupleR resource: ", resource, call. = FALSE)
  )
  net <- net[!is.na(net$mor) & net$mor != 0 &
             nzchar(net$source) & nzchar(net$target), , drop = FALSE]
  unique(net)
}

#' Run decoupleR TF/pathway activity inference
#'
#' @param res_tbl DGE results table (`gene`, `stat`, `log2FoldChange`, `padj`).
#' @param out_dir Output directory.
#' @param comp_name Comparison name used in file names.
#' @param organism "Homo sapiens" / "Mus musculus" (or plain "human"/"mouse").
#' @param resources decoupleR resources to run.
#' @param method One of `"ulm"`, `"mlm"`, `"wsum"`.
#' @param min_size Minimum number of targets per source.
#' @param padj_cutoff Adjusted p-value cutoff used downstream for filtering
#'   (the full tables are always written).
#' @param stat_col Name of the signed statistic column. When absent and both
#'   `log2FoldChange` and `padj` are present, `sign(log2FC) * -log10(padj)`
#'   is used.
#' @param lfc_cutoff Absolute log2 fold-change floor used to categorize regulon
#'   targets (activated/repressed/other) in the targeted volcano and the summary
#'   table. Default 1.
#' @return list(results = named list per resource, files = manifest)
#' @export
run_decouple_analysis <- function(
  res_tbl,
  out_dir,
  comp_name,
  organism = "human",
  resources = c("collectri", "progeny", "dorothea"),
  method = c("ulm", "mlm", "wsum"),
  min_size = 5,
  padj_cutoff = 0.05,
  stat_col = "stat",
  lfc_cutoff = 1
) {
  method <- match.arg(method)
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required. Install with BiocManager::install('decoupleR').")
  }
  res_tbl <- as.data.frame(res_tbl)
  if (nrow(res_tbl) == 0L) return(NULL)

  # Resolve the score column ------------------------------------------------
  if (!stat_col %in% colnames(res_tbl)) {
    if (all(c("log2FoldChange", "padj") %in% colnames(res_tbl))) {
      res_tbl[[stat_col]] <- sign(res_tbl$log2FoldChange) *
        -log10(pmax(res_tbl$padj, .Machine$double.xmin))
    } else {
      stop("No '", stat_col, "' column and cannot derive one from log2FoldChange/padj.")
    }
  }

  # ---- FIX: LRT fallback ----------------------------------------------------
  # Under test = "LRT", run_deseq2_analysis() stores the OMNIBUS LRT
  # statistic in `res_shrunken$stat`. That statistic is unsigned (a
  # likelihood-ratio chi-square: always >= 0), so feeding it to decoupleR
  # would silently kill all directional information: every TF / pathway
  # would look "activated" (positive score) regardless of whether it is
  # actually activated or repressed. The contrast-specific Wald statistic
  # is stored in `contrast_stat`; use it whenever it is available.
  #
  # Wald runs (the default) are unaffected: test_type == "Wald", so the
  # block is skipped.
  if (identical(stat_col, "stat") &&
      "test_type" %in% colnames(res_tbl) &&
      any(toupper(as.character(res_tbl$test_type)) == "LRT", na.rm = TRUE) &&
      "contrast_stat" %in% colnames(res_tbl)) {
    message(
      "   -> LRT detected: using the contrast-specific 'contrast_stat' ",
      "column for decoupleR (the omnibus 'stat' is unsigned)."
    )
    stat_col <- "contrast_stat"
  }

  gene_col <- intersect(c("gene", "symbol", "gene_symbol"), colnames(res_tbl))
  if (length(gene_col) == 0L) stop("No gene/symbol column in res_tbl.")
  gene_col <- gene_col[1]

  keep <- !is.na(res_tbl[[stat_col]]) & res_tbl[[stat_col]] != 0 &
          !is.na(res_tbl[[gene_col]]) & res_tbl[[gene_col]] != ""
  d <- res_tbl[keep, c(gene_col, stat_col), drop = FALSE]
  if (nrow(d) == 0L) return(NULL)

  # One row per symbol, keeping the most extreme |stat| (matches pipeline
  # convention elsewhere: GSEA and ORA input use the same rule).
  d$.abs <- abs(d[[stat_col]])
  d <- d[order(-d$.abs), , drop = FALSE]
  d <- d[!duplicated(d[[gene_col]]), , drop = FALSE]

  mat <- matrix(d[[stat_col]], ncol = 1L,
                dimnames = list(as.character(d[[gene_col]]), comp_name))

  out <- list()
  for (resource in resources) {
    message("   -> decoupleR / ", resource, " (", method, ")")
    net <- tryCatch(.decouple_get_net(resource, organism),
                    error = function(e) {
                      message("      ", conditionMessage(e))
                      NULL
                    })
    if (is.null(net) || nrow(net) == 0L) next

    net <- net[net$target %in% rownames(mat), , drop = FALSE]
    if (nrow(net) < min_size) {
      message("      ", resource, ": only ", nrow(net),
              " overlapping gene(s). Skipping.")
      next
    }

    act <- tryCatch(
      switch(
        method,
        ulm  = decoupleR::run_ulm(mat, net, .source = "source",
                                  .target = "target", .mor = "mor",
                                  minsize = min_size),
        mlm  = decoupleR::run_mlm(mat, net, .source = "source",
                                  .target = "target", .mor = "mor",
                                  minsize = min_size),
        wsum = decoupleR::run_wsum(mat, net, .source = "source",
                                   .target = "target", .mor = "mor",
                                   minsize = min_size)
      ),
      error = function(e) {
        message("      ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(act) || nrow(act) == 0L) next

    act_df <- as.data.frame(act)
    act_df$resource   <- resource
    act_df$comparison <- comp_name

    dir_res <- safe_dir(file.path(out_dir, "DecoupleR", resource))
    utils::write.csv(
      act_df,
      file.path(dir_res, paste0("DecoupleR_", resource, "_", comp_name, ".csv")),
      row.names = FALSE
    )

    out[[resource]] <- list(
      activity    = act_df,
      network     = net,
      top_sources = utils::head(act_df$source[order(-abs(act_df$score))], 20L)
    )
    message("      ", resource, ": ", nrow(act_df), " source(s) tested.")

    # ---- NEW: regulon volcano + leading edge for the significant sources
    tryCatch(
      .run_decouple_tf_regulon_plots(
        activity_df = act_df,
        net         = net,
        res_tbl     = res_tbl,
        resource    = resource,
        out_dir     = out_dir,
        comp_name   = comp_name,
        padj_cutoff = padj_cutoff,
        lfc_cutoff  = lfc_cutoff,
        top_n       = 20L
      ),
      error = function(e) message("      [", resource,
                                  "] regulon plots failed: ", e$message)
    )
  }

  if (length(out) == 0L) return(NULL)
  list(results = out,
       files   = .out_manifest(out_dir, "DecoupleR"))
}

# =============================================================================
# TF regulon plots: targeted volcano + leading edge, per significant source
# =============================================================================
#
# Called from run_decouple_analysis() once per resource, AFTER the activity
# table is written. For every source with p < padj_cutoff (top N by |score|):
#
#   * regulon_<TF>.csv       : matched targets + mor + padj + LFC + category
#   * leading_edge_<TF>.csv  : rank + target + stat (sorted by rank)
#   * tf_volcano_<TF>.pdf    : volcano limited to the TF's regulon
#   * tf_leading_edge_<TF>.pdf : barcode plot over all genes ranked by stat
#   * tf_regulon_summary.csv : per-source counts (targets / up / down / ns)
#
# The category assigned to each target follows the user's rule:
#   mor > 0 & LFC >  lfc_cutoff & padj < padj_cutoff  -> "Activated (expected)"
#   mor < 0 & LFC < -lfc_cutoff & padj < padj_cutoff  -> "Repressed (expected)"
#   padj < padj_cutoff & |LFC| > lfc_cutoff           -> "Other change"
#   else                                              -> "Not significant"
#
# Under test = "LRT", the DGE results table's `stat` column is the unsigned
# OMNIBUS statistic; the signed contrast Wald statistic lives in
# `contrast_stat` and is used for ranking instead (see run_decouple_analysis()
# for the same fallback on the activity table side).

#' Build per-TF regulon plots for the significant sources of one resource
#'
#' @param activity_df Activity table returned by decoupleR (source / score / p_value)
#' @param net          The resource's network (source / target / mor)
#' @param res_tbl      The DGE results table (gene / stat / contrast_stat /
#'                     log2FoldChange / padj / test_type)
#' @param resource     Resource name ("collectri", "progeny", "dorothea")
#' @param out_dir      Pipeline output directory
#' @param comp_name    Comparison name for messages
#' @param padj_cutoff  Adjusted p-value cutoff for regulon categories
#' @param lfc_cutoff   |log2FC| floor for regulon categories
#' @param top_n        Cap on the number of sources plotted (by |score|)
#' @return Invisible list of per-source summary rows
#' @keywords internal
.run_decouple_tf_regulon_plots <- function(
  activity_df,
  net,
  res_tbl,
  resource,
  out_dir,
  comp_name,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  top_n = 20L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (is.null(activity_df) || !nrow(activity_df)) return(invisible(NULL))
  if (is.null(net) || !nrow(net)) return(invisible(NULL))
  if (is.null(res_tbl) || !nrow(res_tbl)) return(invisible(NULL))

  pcol <- .rep_first_col(activity_df, c("p_value", "padj", "pval", "p.value"))
  if (is.na(pcol)) return(invisible(NULL))

  sig <- activity_df[!is.na(activity_df[[pcol]]) &
                     activity_df[[pcol]] < padj_cutoff, , drop = FALSE]
  if (!nrow(sig)) {
    message("      [", resource, "] no source with p < ", padj_cutoff,
            " -- skipping regulon plots.")
    return(invisible(NULL))
  }
  sig <- sig[order(-abs(sig$score)), , drop = FALSE]
  sig <- utils::head(sig, top_n)

  # Signed statistic for the leading-edge ranking (LRT-safe).
  is_lrt <- "test_type" %in% colnames(res_tbl) &&
    any(toupper(as.character(res_tbl$test_type)) == "LRT", na.rm = TRUE)
  stat_col <- if (is_lrt && "contrast_stat" %in% colnames(res_tbl)) {
    "contrast_stat"
  } else if ("stat" %in% colnames(res_tbl)) {
    "stat"
  } else {
    "log2FoldChange"
  }

  keep <- !is.na(res_tbl[[stat_col]]) &
          !is.na(res_tbl$gene) & res_tbl$gene != ""
  ranks_df <- res_tbl[keep, c("gene", stat_col), drop = FALSE]
  ranks_df <- ranks_df[order(-ranks_df[[stat_col]]), , drop = FALSE]
  ranks_df <- ranks_df[!duplicated(ranks_df$gene), , drop = FALSE]
  ranks <- stats::setNames(ranks_df[[stat_col]], ranks_df$gene)
  ranks <- ranks[order(ranks, decreasing = TRUE)]

  plots_dir <- safe_dir(file.path(out_dir, "DecoupleR", resource, "plots"))

  summary_rows <- list()

  for (i in seq_len(nrow(sig))) {
    tf     <- as.character(sig$source[i])
    tf_s   <- sig$score[i]
    tf_p   <- sig[[pcol]][i]

    tf_net <- net[net$source == tf, , drop = FALSE]
    if (!nrow(tf_net)) next

    regulon <- merge(
      tf_net[, c("target", "mor"), drop = FALSE],
      res_tbl[, intersect(
        c("gene", "log2FoldChange", "padj", stat_col, "test_type"),
        colnames(res_tbl)
      ), drop = FALSE],
      by.x = "target", by.y = "gene", all = FALSE
    )
    if (!nrow(regulon)) {
      message("      ", tf, ": no target overlaps the DGE table.")
      next
    }

    # Categorise each target -- user's rule, exact.
    reg <- rep("Not significant", nrow(regulon))
    up_expected <- !is.na(regulon$mor) & regulon$mor > 0 &
                   !is.na(regulon$log2FoldChange) &
                   regulon$log2FoldChange >  lfc_cutoff &
                   !is.na(regulon$padj) & regulon$padj < padj_cutoff
    dn_expected <- !is.na(regulon$mor) & regulon$mor < 0 &
                   !is.na(regulon$log2FoldChange) &
                   regulon$log2FoldChange < -lfc_cutoff &
                   !is.na(regulon$padj) & regulon$padj < padj_cutoff
    reg[up_expected] <- "Activated (expected)"
    reg[dn_expected] <- "Repressed (expected)"
    other <- !is.na(regulon$padj) & regulon$padj < padj_cutoff &
             !is.na(regulon$log2FoldChange) &
             abs(regulon$log2FoldChange) > lfc_cutoff &
             reg == "Not significant"
    reg[other] <- "Other change"
    regulon$Regulation <- factor(
      reg,
      levels = c("Not significant", "Other change",
                 "Activated (expected)", "Repressed (expected)")
    )

    safe_tf <- gsub("[^A-Za-z0-9_-]", "_", tf)

    # ---- save regulon table ----
    utils::write.csv(
      regulon,
      file.path(plots_dir, paste0("regulon_", safe_tf, ".csv")),
      row.names = FALSE
    )

    # ---- save leading-edge table (rank of each target in the global ordering)
    le <- data.frame(
      rank   = match(regulon$target, names(ranks)),
      target = regulon$target,
      stat   = unname(ranks[regulon$target]),
      stringsAsFactors = FALSE
    )
    le <- le[!is.na(le$rank), , drop = FALSE]
    le <- le[order(le$rank), , drop = FALSE]
    if (nrow(le)) {
      utils::write.csv(
        le,
        file.path(plots_dir, paste0("leading_edge_", safe_tf, ".csv")),
        row.names = FALSE
      )
    }

    # ---- static PDFs ----
    p_vol <- .tf_regulon_volcano_plot(regulon, tf, tf_s, tf_p,
                                       padj_cutoff, lfc_cutoff)
    if (!is.null(p_vol)) {
      .safe_ggsave(
        file.path(plots_dir, paste0("tf_volcano_", safe_tf, ".pdf")),
        p_vol, width = 9, height = 7
      )
    }
    if (nrow(le) && length(ranks) > 0) {
      p_le <- .tf_leading_edge_plot(le, ranks, tf, tf_s, tf_p)
      if (!is.null(p_le)) {
        .safe_ggsave(
          file.path(plots_dir, paste0("tf_leading_edge_", safe_tf, ".pdf")),
          p_le, width = 10, height = 5
        )
      }
    }

    summary_rows[[tf]] <- data.frame(
      resource          = resource,
      source            = tf,
      score             = tf_s,
      p_value           = tf_p,
      n_targets_total   = nrow(tf_net),
      n_targets_matched = nrow(regulon),
      n_activated       = sum(regulon$Regulation == "Activated (expected)"),
      n_repressed       = sum(regulon$Regulation == "Repressed (expected)"),
      n_other_var       = sum(regulon$Regulation == "Other change"),
      n_ns              = sum(regulon$Regulation == "Not significant"),
      stringsAsFactors  = FALSE
    )
    message("      ", resource, " / ", tf,
            ": ", nrow(regulon), " targets (",
            sum(regulon$Regulation == "Activated (expected)"), " up, ",
            sum(regulon$Regulation == "Repressed (expected)"), " down)")
  }

  if (length(summary_rows)) {
    summary_df <- do.call(rbind, summary_rows)
    summary_df <- summary_df[order(summary_df$resource,
                                    -abs(summary_df$score)), , drop = FALSE]
    utils::write.csv(
      summary_df,
      file.path(out_dir, "DecoupleR", "tf_regulon_summary.csv"),
      row.names = FALSE
    )
  }

  invisible(summary_rows)
}

#' Targeted volcano plot for one TF's regulon
#' @keywords internal
.tf_regulon_volcano_plot <- function(
  regulon, tf, tf_score, tf_p, padj_cutoff, lfc_cutoff
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  regulon$padj_plot <- pmax(regulon$padj, .Machine$double.xmin)
  colors <- c(
    "Not significant"      = "grey80",
    "Other change"         = "#ff7f00",
    "Activated (expected)" = "#1f78b4",
    "Repressed (expected)" = "#e31a1c"
  )
  p <- ggplot2::ggplot(
    regulon,
    ggplot2::aes(x = log2FoldChange, y = -log10(padj_plot), color = Regulation)
  ) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                        linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                        linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste("Transcriptional targets of", tf),
      subtitle = sprintf(
        "Score: %+.2f  |  p: %.2e  |  %d targets in the DGE (CollecTRI/DoRothEA)",
        tf_score, tf_p, nrow(regulon)
      ),
      x = "log2 Fold Change",
      y = "-log10(padj)",
      color = NULL
    )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    lbl <- regulon[regulon$Regulation %in% c("Activated (expected)",
                                             "Repressed (expected)"),
                   , drop = FALSE]
    if (nrow(lbl)) {
      lbl <- utils::head(lbl[order(lbl$padj), , drop = FALSE], 20)
      p <- p + ggrepel::geom_text_repel(
        data = lbl,
        mapping = ggplot2::aes(label = target),
        size = 3.2, max.overlaps = 25, seed = 42, inherit.aes = TRUE
      )
    }
  }
  p
}

#' Barcode / leading-edge plot for one TF's regulon
#' @keywords internal
.tf_leading_edge_plot <- function(le, ranks, tf, tf_score, tf_p) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  n_genes <- length(ranks)
  positions <- le$rank[order(le$rank)]
  n_hits <- length(positions)
  if (!n_genes || !n_hits) return(NULL)

  is_hit <- rep(FALSE, n_genes)
  is_hit[positions] <- TRUE
  hit_w  <- 1 / n_hits
  miss_w <- -1 / max(1L, n_genes - n_hits)
  es_vec <- cumsum(ifelse(is_hit, hit_w, miss_w))
  peak_idx <- which.max(abs(es_vec))

  es_df   <- data.frame(rank = seq_len(n_genes), es = es_vec)
  hits_df <- data.frame(rank = positions)

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = hits_df,
      ggplot2::aes(xmin = rank - 0.5, xmax = rank + 0.5,
                   ymin = min(es_vec) - 0.02,
                   ymax = min(es_vec) - 0.005),
      fill = "#e31a1c", alpha = 0.8
    ) +
    ggplot2::geom_line(data = es_df,
                       ggplot2::aes(x = rank, y = es),
                       color = "#1f78b4", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = peak_idx,
                        color = "grey40", linetype = "dotted") +
    ggplot2::labs(
      title = paste("Leading edge of", tf, "targets"),
      subtitle = sprintf(
        "Score: %+.2f  |  p: %.2e  |  %d targets  |  max ES = %.3f at rank %d",
        tf_score, tf_p, n_hits, es_vec[peak_idx], peak_idx
      ),
      x = "Gene rank (by statistic, decreasing)",
      y = "Cumulative enrichment score (running ES)"
    ) +
    ggplot2::theme_minimal()
}