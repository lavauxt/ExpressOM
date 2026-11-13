# Test fixtures: a synthetic ExpressOM output directory.
#
# It mirrors, file for file, what the pipeline's writers produce (DE tables,
# PCA scores, run snapshots, ORA / GSEA / FGSEA / SPIA / TF tables, isoform
# tables and plots), without running DESeq2 or needing an Ensembl database, so
# the report / app / targets tests run anywhere.

.fx_write <- function(df, ..., sep = ",") {
  path <- file.path(...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (identical(sep, ",")) {
    utils::write.csv(df, path, row.names = FALSE)
  } else {
    utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  invisible(path)
}

.fx_png <- function(path, label = "plot") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 320, height = 220)
  graphics::plot(1:5, main = label)
  grDevices::dev.off()
  invisible(path)
}

.fx_ora <- function(n, prefix) {
  p <- sort(stats::runif(n)^3)
  data.frame(
    ID = paste0("ID:", seq_len(n)),
    Description = paste(prefix, seq_len(n), "pathway"),
    GeneRatio = paste0(sample(3:20, n, TRUE), "/200"),
    BgRatio = paste0(sample(30:300, n, TRUE), "/18000"),
    pvalue = p, p.adjust = pmin(1, p * 3), qvalue = pmin(1, p * 2.5),
    geneID = "A/B/C", Count = sample(3:20, n, TRUE),
    stringsAsFactors = FALSE
  )
}

.fx_gsea <- function(n, prefix) {
  p <- sort(stats::runif(n)^2)
  nes <- stats::rnorm(n, 0, 1.8)
  data.frame(
    ID = paste0("ID:", seq_len(n)),
    Description = paste(prefix, seq_len(n), "pathway"),
    setSize = sample(15:300, n, TRUE), enrichmentScore = nes / 3, NES = nes,
    pvalue = p, p.adjust = pmin(1, p * 2), qvalue = pmin(1, p * 2),
    stringsAsFactors = FALSE
  )
}

.fx_fgsea <- function(n, prefix) {
  p <- sort(stats::runif(n)^2)
  nes <- stats::rnorm(n, 0, 1.8)
  data.frame(
    pathway = paste0(prefix, "_", seq_len(n)), pval = p, padj = pmin(1, p * 2),
    log2err = 0.3, ES = nes / 3, NES = nes, size = sample(15:300, n, TRUE),
    stringsAsFactors = FALSE
  )
}

#' Create a synthetic ExpressOM output directory
#'
#' @param dir Directory to fill (created if needed).
#' @param level,base Contrast levels.
#' @param test "Wald" or "LRT".
#' @param n_genes Number of genes in the DE table.
#' @param functional,isoform,figures Which optional stages to include.
#' @param snapshot Write `log/expressom_run.json` (the run snapshot).
#' @return `dir`.
make_run_dir <- function(
  dir,
  level = "treated",
  base = "control",
  test = "Wald",
  n_genes = 400L,
  functional = TRUE,
  isoform = TRUE,
  figures = TRUE,
  snapshot = TRUE,
  seed = 11L
) {
  withr::local_seed(seed)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cmp <- paste0(level, "_vs_", base)
  samples <- paste0("S", 1:6)
  cond <- rep(c(base, level), each = 3)

  # ---- differential expression (writer layout of export_significant_results)
  lfc <- stats::rnorm(n_genes, 0, 1.4)
  pval <- 2 * stats::pnorm(-abs(lfc / 0.45))
  padj <- stats::p.adjust(pval, "BH")
  de <- data.frame(
    gene = paste0("GENE", seq_len(n_genes)),
    ensembl = sprintf("ENSG%011d", seq_len(n_genes)),
    entrezid = as.character(1000L + seq_len(n_genes)),
    baseMean = stats::rexp(n_genes, 1 / 150),
    log2FoldChange = lfc, lfcSE = 0.45, pvalue = pval, padj = padj,
    stat = lfc / 0.45, test_type = test,
    stringsAsFactors = FALSE
  )
  de <- de[order(de$padj), ]
  de_dir <- file.path(dir, "DE_raw_results")
  .fx_write(de, de_dir, paste0("DEgenes_raw_", cmp, ".txt"), sep = "\t")
  .fx_write(de[de$padj < 0.05, ], de_dir, paste0("DEgenes_pval_0_05_", cmp, ".txt"), sep = "\t")
  .fx_write(de[de$padj < 0.01, ], de_dir, paste0("DEgenes_pval_0_01_", cmp, ".txt"), sep = "\t")

  counts <- data.frame(
    padj_cutoff = c(0.01, 0.05),
    n_tested = n_genes,
    n_sig = c(sum(de$padj < 0.01), sum(de$padj < 0.05)),
    n_up = if (test == "LRT") NA_integer_ else c(sum(de$padj < 0.01 & de$log2FoldChange > 0), sum(de$padj < 0.05 & de$log2FoldChange > 0)),
    n_down = if (test == "LRT") NA_integer_ else c(sum(de$padj < 0.01 & de$log2FoldChange < 0), sum(de$padj < 0.05 & de$log2FoldChange < 0)),
    log2FC_cutoff = 1,
    n_sig_lfc = c(sum(de$padj < 0.01 & abs(de$log2FoldChange) > 1), sum(de$padj < 0.05 & abs(de$log2FoldChange) > 1)),
    n_up_lfc = if (test == "LRT") NA_integer_ else c(sum(de$padj < 0.01 & de$log2FoldChange > 1), sum(de$padj < 0.05 & de$log2FoldChange > 1)),
    n_down_lfc = if (test == "LRT") NA_integer_ else c(sum(de$padj < 0.01 & de$log2FoldChange < -1), sum(de$padj < 0.05 & de$log2FoldChange < -1)),
    test_type = test
  )
  .fx_write(counts, de_dir, paste0("DEG_counts_", cmp, ".txt"), sep = "\t")

  norm <- matrix(
    stats::rpois(n_genes * 6, lambda = rep(de$baseMean, 6)) + 1,
    nrow = n_genes, dimnames = list(NULL, samples)
  )
  norm <- data.frame(ensembl = de$ensembl, norm, symbol = de$gene, entrezid = de$entrezid, gene = de$gene, check.names = FALSE)
  .fx_write(norm, de_dir, paste0("DESeq_Normalized_Counts_", cmp, ".txt"), sep = "\t")

  # ---- PCA (writer layout of .write_pca_scores)
  pcs <- matrix(stats::rnorm(6 * 6, sd = c(20, 8, 4, 2, 1, 0.5)), nrow = 6, byrow = TRUE,
                dimnames = list(NULL, paste0("PC", 1:6)))
  scores <- data.frame(sample_label = samples, condition = cond, batch = rep(c("b1", "b2"), 3), pcs, check.names = FALSE)
  .fx_write(scores, dir, "Plots", paste0("PCA_", cmp, "_sample_scores.tsv"), sep = "\t")
  .fx_write(
    data.frame(PC = paste0("PC", 1:6), PercentVariance = c(55, 22, 10, 7, 4, 2)),
    dir, "Plots", paste0("PCA_", cmp, "_variance_explained.tsv"), sep = "\t"
  )

  # ---- logs / run snapshot
  dir.create(file.path(dir, "log", "DGE"), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    list(
      comparison = cmp, design_formula = "~ condition", test_type = test,
      reduced_formula = if (test == "LRT") "~ 1" else "none",
      contrast = list(factor = "condition", level = level, base = base),
      shrinkage_method = "ashr", padj_cutoff = 0.05, lfc_cutoff = 1,
      genes_before_filter = n_genes, genes_after_filter = n_genes
    ),
    file.path(dir, "log", "DGE", paste0("DGE_params_", cmp, ".json")),
    auto_unbox = TRUE, pretty = TRUE
  )
  if (snapshot) {
    ExpressOM:::.write_run_config(dir, list(
      level = level, base = base, model = "~ condition", test = test,
      padj_cutoff = 0.05, deg_lfc_cutoff = 1, shrink_method = "ashr",
      count_type = "matrix", out_dir = dir,
      volcano_colors = c(up = "red2", down = "royalblue", ns = "grey"),
      pca_colors = c(level = "red2", base = "royalblue")
    ))
  }
  writeLines("R version x.y.z", file.path(dir, "log", "SessionInfo.txt"))

  # ---- functional analysis
  if (functional) {
    .fx_write(.fx_ora(60, "GO BP"), dir, "ORA", "GO", paste0("GO_BP_", cmp, ".csv"))
    .fx_write(.fx_ora(30, "Reactome"), dir, "ORA", "Reactome", paste0("Reactome_ORA_", cmp, ".csv"))
    .fx_write(.fx_gsea(40, "GO BP"), dir, "GSEA", "GO", paste0("GO_GSEA_BP_", cmp, ".csv"))
    fg <- .fx_fgsea(80, "REACTOME")
    .fx_write(fg, dir, "GSEA", "FGSEA", "C2_CP_REACTOME", "FGSEA_Results_C2_CP_REACTOME.csv")
    .fx_write(fg[fg$padj < 0.05, ], dir, "GSEA", "FGSEA", "C2_CP_REACTOME", "FGSEA_Results_C2_CP_REACTOME_significant.csv")
    .fx_write(.fx_fgsea(30, "HALLMARK"), dir, "GSEA", "FGSEA", "H", "FGSEA_Results_H.csv")
    .fx_write(
      data.frame(
        Name = paste("KEGG", 1:15), ID = 4000 + 1:15, pSize = 40, NDE = 10,
        tA = stats::rnorm(15), pNDE = stats::runif(15), pPERT = stats::runif(15),
        pG = sort(stats::runif(15)^2), pGFdr = sort(stats::runif(15)^2), pGFWER = 1,
        Status = sample(c("Activated", "Inhibited"), 15, TRUE), KEGGLINK = "http://x"
      ),
      dir, "SPIA", paste0("SPIA_Results_", cmp, ".csv")
    )
    .fx_write(
      data.frame(
        Term = paste0("TF", 1:20), Overlap = "4/150", P.value = sort(stats::runif(20)^3),
        Adjusted.P.value = sort(stats::runif(20)^3) * 2, Odds.Ratio = stats::runif(20, 1, 9),
        Combined.Score = stats::runif(20, 5, 90), Genes = "A;B;C", Database = "ChEA_2022"
      ),
      dir, "Transcription_Factors", paste0("Enrichr_TF_", cmp, ".csv")
    )
  }

  # ---- isoforms
  if (isoform) {
    n <- 300L
    l <- stats::rnorm(n, 0, 1.3)
    p <- 2 * stats::pnorm(-abs(l / 0.5))
    dte <- data.frame(
      transcript_id = sprintf("ENST%011d", seq_len(n)), baseMean = stats::rexp(n, 1 / 200),
      log2FoldChange = l, lfcSE = 0.5, stat = l / 0.5, pvalue = p, padj = stats::p.adjust(p, "BH"),
      gene_id = sprintf("ENSG%011d", sample(1:100, n, TRUE)),
      gene_symbol = paste0("GENE", sample(1:100, n, TRUE)), entrezid = "1",
      stringsAsFactors = FALSE
    )
    dte$gene <- dte$gene_symbol
    dte$signif <- dte$padj < 0.05 & abs(dte$log2FoldChange) > 1
    dte$direction_label <- ifelse(dte$signif, ifelse(dte$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not significant")
    iso <- file.path(dir, "IsoformSwitch", "DTU_DTE_report")
    .fx_write(dte, iso, "DTE_results_annotated.csv")
    pv <- stats::runif(120)^2
    dtu <- data.frame(
      gene_id = sprintf("ENSG%011d", 1:120), feature_id = sprintf("ENST%011d", 1:120),
      pvalue = pv, adj_pvalue = stats::p.adjust(pv, "BH"), gene_symbol = paste0("GENE", 1:120),
      stringsAsFactors = FALSE
    )
    dtu$gene_label <- paste0(dtu$gene_symbol, " (", dtu$feature_id, ")")
    .fx_write(dtu, iso, "DTU_results_annotated.csv")
    writeLines("DTE adjusted p-values are transcript-level DESeq2 Wald results.", file.path(iso, "analysis_methods.txt"))
    .fx_png(file.path(iso, "plots", "DTE_volcano.png"), "DTE")
    .fx_write(.fx_ora(20, "DTE GO"), dir, "IsoformSwitch", "DTE_Enrichment", "ORA", "GO", paste0("GO_BP_", cmp, ".csv"))
    .fx_write(.fx_ora(15, "DTU GO"), dir, "IsoformSwitch", "DTU_Enrichment", "ORA", "GO", paste0("GO_BP_", cmp, ".csv"))
  }

  # ---- static figures
  if (figures) {
    dir.create(file.path(dir, "Plots"), showWarnings = FALSE)
    grDevices::pdf(file.path(dir, "Plots", paste0("DE_Volcanoplot_", cmp, ".pdf")), 5, 4)
    graphics::plot(1:10)
    grDevices::dev.off()
    .fx_png(file.path(dir, "Plots", paste0("Heatmap_", cmp, ".png")), "heatmap")
  }
  dir
}

# Minimal sample table + count matrix so file-fingerprinting has real inputs
make_inputs <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  st <- file.path(dir, "samples.csv")
  mx <- file.path(dir, "counts.csv")
  writeLines(c("sample,condition", paste0("S", 1:6, ",", rep(c("control", "treated"), each = 3))), st)
  writeLines(c("gene,S1,S2,S3,S4,S5,S6", "g1,1,2,3,4,5,6"), mx)
  list(sample_table = st, matrix_file = mx, dir = dir)
}

skip_if_no_report_stack <- function() {
  testthat::skip_if_not_installed("plotly")
  testthat::skip_if_not_installed("DT")
  testthat::skip_if_not_installed("rmarkdown")
  testthat::skip_if_not_installed("htmltools")
  testthat::skip_if_not_installed("htmlwidgets")
  testthat::skip_if_not_installed("knitr")
}

skip_if_no_pandoc <- function() {
  testthat::skip_if_not_installed("rmarkdown")
  testthat::skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
}
