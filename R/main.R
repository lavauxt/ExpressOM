#' Run Expressom Pipeline
#'
#' @description
#' Executes the complete bulk RNA-seq pipeline including data import, exploratory
#' data analysis, differential expression (DESeq2), functional analysis, and
#' optional isoform-level analysis (DTE, DTU, IsoformSwitchAnalyzeR). When it is
#' done it builds ONE interactive HTML report (plotly figures and searchable
#' tables in tabs) that regroups every result; see [build_expressom_report()].
#'
#' The pipeline is made of stages that are also driven by `targets`
#' ([run_expressom_targets()]) and by the Shiny app ([run_expressom_app()]),
#' which re-run only the stages affected by a parameter change.
#'
#' @export
#' @param count_type Type of RNA count (default: "salmon")
#' @param data_dir Folder where the input data is stored
#' @param out_dir Folder to save results
#' @param sample_table Path to the sample table metadata
#' @param level Foreground level for contrast
#' @param base Background level for contrast
#' @param model Design formula for DESeq2 as a string
#' @param batch_col Optional column specifying batch effects for correction/visualization (uses limma)
#' @param replicate_col Optional column specifying replicates
#' @param shrink_method Type of shrinkage (e.g., "ashr")
#' @param ensembl_package_name The Ensembl DB package. If not installed, will be built and installed automatically.
#' @param top_genes Number of top genes to plot
#' @param gmt_file Path to a local GMT file, or a list/vector of paths. If NULL, msigdbr downloads Hallmark.
#' @param padj_cutoff Adjusted p-value threshold for significance
#' @param test Type of statistical test ("Wald" or "LRT")
#' @param reduced Design formula for the reduced model (used if test = "LRT")
#' @param highlight_genes Optional character vector of gene names to highlight in the Volcano plot
#' @param go_pvalue_cutoff GO ORA p-value cutoff
#' @param go_qvalue_cutoff GO ORA q-value cutoff
#' @param ora_padj_cutoff Adjusted p-value cutoff for the gene list fed into
#'   ORA (GO/Reactome/DO). Kept separate from `padj_cutoff` since ORA is a
#'   hypergeometric test that's underpowered by construction when its input
#'   list is short. Default 0.05.
#' @param ora_min_genes Minimum gene count at `ora_padj_cutoff` before ORA's
#'   input list falls back to raw `pvalue < 0.05`. Default 10.
#' @param ora_lfc_cutoff Optional `abs(log2FoldChange) >=` floor applied
#'   uniformly to GO/Reactome/DO ORA input. `NULL` (default) applies none.
#' @param matrix_file Path to raw counts file if count_type = 'matrix'
#' @param custom_tx2gene Path to a custom transcript-to-gene mapping file (TSV with columns 'tx_id' and 'gene_id')
#' @param custom_gene_map Path to a custom gene annotation file (TSV with columns 'gene_id', 'symbol', and optionally 'entrezid')
#' @param gsea_metric Metric to rank genes for GSEA ("stat", "signed_pval", or "log2FoldChange")
#' @param subset_sample Optional string to filter the sample table (e.g., "cell_type == 'T_cells'")
#' @param remove_sample Optional character vector of sample IDs to exclude entirely
#' @param zscore_genes Optional character vector of gene names for targeted Z-score expression plotting
#' @param gene_sets_zscore Optional named list of character vectors of gene names.
#' @param run_dge Logical: perform standard gene-level differential expression (default: TRUE)
#' @param run_isoform Logical: perform isoform-level analysis (DTE, DTU, Isoform SwitchAnalyzeR)
#' @param run_predictors Logical: Run CPAT, Pfam, and SignalP via WSL/Conda during Isoform analysis
#' @param bpparam BiocParallel backend for multi-threading
#' @param execution_order String: "dge_first" or "isoform_first"
#' @param isoform_fasta Path to transcript FASTA file
#' @param isoform_gff Path to GFF/GTF annotation file
#' @param use_wsl Logical: route external predictor tools through WSL on Windows
#' @param wsl_distro Name of WSL distribution
#' @param predictor_cpu Integer: CPU threads for hmmscan / InterProScan
#' @param resume_isoform_from Path to directory with saved DTE/DTU RDS files to resume isoform analysis
#' @param isoform_report_genes Gene symbols for transcript-proportion plots and enhanced isoform visualizations
#' @param run_dexseq Logical: also run DEXSeq-based DTU engine
#' @param run_isoform_enrichment Logical: run functional enrichment (ORA, and
#'   FGSEA for DTE) on DTE/DTU results, using the same GO/Reactome/DO/MSigDB
#'   databases as the DGE-level analysis. Runs once DTE/DTU are available,
#'   before `run_isoform_switch()`. See `run_isoform_functional_analysis()`.
#'   Default TRUE.
#' @param isoform_plot_top_n Number of top isoform switches to render automatically
#' @param nBest Deprecated and ignored. It only controlled the RegionReport
#'   report, which was replaced by the unified interactive report.
#' @param eda_only Logical: if TRUE, run only import + EDA then stop
#' @param group_col Optional column name in metadata for colouring PCA/heatmap when no model is given
#' @param custom_transcript_id_map Path to a TSV/CSV with columns `count_id` and `fasta_id`
#' @param skip_fasta_filter Logical: if TRUE, skip FASTA pre-filtering in isoform import
#' @param isoform_test_engine Which DTU test engine to use in IsoformSwitchAnalyzeR
#' @param plot_topology Logical: include the topology (DeepTMHMM) track in switch plots.
#'   Set FALSE to silence "Omitting topology visualization..." when DeepTMHMM hasn't been run.
#' @param run_spia Logical: also run SPIA (pathway impact analysis). Off by default --
#'   its bootstrap permutation testing is slow and less commonly needed than GSEA/ORA.
#' @param pca_ntop Number of most variable genes used for the PCA plots (default 500; NULL uses all)
#' @param sample_labels Optional relabelling of the samples printed on the PCA plots: a named
#'   character vector (names = sample IDs as in the sample table, values = new labels), a
#'   data.frame, or the path to a CSV/TSV file with a sample-ID column and a label column
#'   (headers such as `sample`/`sample_id` and `label`/`new_label` are recognised, otherwise the
#'   first two columns are used). Samples without an entry keep their name. Default NULL.
#' @param pca_colors Named vector `c(level = ..., base = ...)` with the PCA colours of the `level`
#'   and `base` groups. Default `c(level = "red2", base = "royalblue")`; either name may be omitted.
#' @param volcano_colors Named vector with the volcano colours of `up`, `down` and `ns` (not
#'   significant) genes. Default `c(up = "red2", down = "royalblue", ns = "grey")`.
#' @param deg_padj_cutoffs Numeric vector of adjusted p-value cutoffs for which the number of DEGs
#'   (total / up / down, with and without the log2FC filter) is saved to
#'   `DE_raw_results/DEG_counts_<level>_vs_<base>.txt`, together with one DEG list per cutoff
#'   (`DEgenes_pval_<cutoff>_<level>_vs_<base>.txt`). `padj_cutoff` is always included.
#'   Default `c(0.01, 0.05)`; `NULL` writes only the list for `padj_cutoff`.
#' @param deg_lfc_cutoff Absolute log2 fold-change threshold used by the volcano plot (lines,
#'   colouring and the up/down counts in its caption) and by the DEG count table. Default 1.
#' @param make_report Logical: build the unified interactive HTML report in
#'   `<out_dir>/Report/` once the analyses are done (default `TRUE`). It needs
#'   plotly, DT, rmarkdown and Pandoc; if they are missing the analysis still
#'   completes and a warning explains how to build the report later with
#'   [build_expressom_report()]. When `FALSE`, the legacy stand-alone
#'   `DTU_DTE_report.html` of the isoform module is rendered instead.
#' @param report_settings Optional named list overriding defaults of
#'   [expressom_report_settings()] (e.g. `list(title = "My study",
#'   highlight_genes = c("TP53", "MYC"))`).
#' @return NULL invisibly
expressom <- function(
  count_type = "salmon",
  data_dir = "./data",
  out_dir = "./results",
  sample_table = "./sample_table.csv",
  matrix_file = NULL,
  custom_tx2gene = NULL,
  custom_gene_map = NULL,
  level = "treated",
  base = "control",
  model = "~ condition",
  batch_col = NULL,
  replicate_col = NULL,
  shrink_method = "ashr",
  pca_ntop = 500,
  ensembl_package_name = "EnsDb.Hsapiens.v107",
  top_genes = 30,
  gmt_file = c(
    "C2:CP:BIOCARTA",
    "C2:CP:KEGG_LEGACY",
    "C2:CP:PID",
    "C2:CP:REACTOME",
    "C2:CP:WIKIPATHWAYS",
    "C3:TFT",
    "C5:GO:BP",
    "C5:GO:CC",
    "C5:GO:MF",
    "C8"
  ),
  padj_cutoff = 0.01,
  go_pvalue_cutoff = 0.05,
  go_qvalue_cutoff = 0.2,
  ora_padj_cutoff = 0.05,
  ora_min_genes = 10,
  ora_lfc_cutoff = NULL,
  test = "Wald",
  reduced = NULL,
  highlight_genes = NULL,
  subset_sample = NULL,
  remove_sample = NULL,
  gsea_metric = "stat",
  nBest = NULL,
  zscore_genes = NULL,
  gene_sets_zscore = NULL,
  run_dge = TRUE,
  run_isoform = FALSE,
  run_predictors = FALSE,
  bpparam = BiocParallel::bpparam(),
  execution_order = c("dge_first", "isoform_first"),
  isoform_fasta = NULL,
  isoform_gff = NULL,
  use_wsl = (.Platform$OS.type == "windows"),
  wsl_distro = "Ubuntu-22.04",
  predictor_cpu = NULL,
  resume_isoform_from = NULL,
  isoform_report_genes = NULL,
  run_dexseq = FALSE,
  run_isoform_enrichment = TRUE,
  isoform_plot_top_n = 10,
  eda_only = FALSE,
  group_col = NULL,
  custom_transcript_id_map = NULL,
  skip_fasta_filter = FALSE,
  isoform_test_engine = c("DEXSeq", "DRIMSeq", "satuRn"),
  plot_topology = TRUE,
  run_spia = FALSE,
  sample_labels = NULL,
  pca_colors = c(level = "red2", base = "royalblue"),
  volcano_colors = c(up = "red2", down = "royalblue", ns = "grey"),
  deg_padj_cutoffs = c(0.01, 0.05),
  deg_lfc_cutoff = 1,
  make_report = TRUE,
  report_settings = NULL
) {
  # every argument, evaluated -- must stay the first statement so that no
  # local variable leaks into the parameter list
  p <- as.list(environment())

  ctx <- .expressom_prepare(p)
  p <- ctx$p
  out_dir <- p$out_dir

  # Normalise vector defaults that are single choices: the `p <- as.list(environment())`
  # above captures the evaluated length-2 default of execution_order (and of
  # isoform_test_engine), and a length-2 value turns
  # `p$execution_order == "dge_first"` into a length-2 logical -- an error
  # inside if() on R >= 4.2, a silent first-element pick on older versions.
  # Taking the first element keeps the documented default; if
  # .expressom_prepare() already normalised them this is a harmless no-op.
  p$execution_order <- if (length(p$execution_order)) p$execution_order[[1L]] else "dge_first"
  p$isoform_test_engine <- if (length(p$isoform_test_engine)) p$isoform_test_engine[[1L]] else "DEXSeq"

  old_warn <- getOption("nwarnings")
  options(nwarnings = 10000)

  on.exit(
    {
      options(nwarnings = old_warn)

      log_dir <- file.path(out_dir, "Log")
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

      if (dir.exists(log_dir)) {
        capture.output(
          sessionInfo(),
          file = file.path(log_dir, "SessionInfo.txt")
        )

        if (!is.null(warnings())) {
          capture.output(warnings(), file = file.path(log_dir, "Warnings.txt"))
        }
      }
    },
    add = TRUE
  )

  # parameter snapshot read back by the report and the Shiny app
  .write_run_config(out_dir, p)

  ## ------------------------------------------------------------------
  ## Always import counts and run EDA
  ## ------------------------------------------------------------------
  tx_data <- .stage_import(ctx)
  .stage_eda(ctx, tx_data)

  if (!p$run_dge && !p$run_isoform) {
    message(
      "EDA completed. Skipping differential expression and isoform analyses."
    )
    if (isTRUE(p$make_report)) {
      safe_run(.stage_report(ctx), label = "Interactive report")
    }
    return(invisible(NULL))
  }

  pipeline_steps <- if (p$execution_order == "dge_first") {
    c("dge", "isoform")
  } else {
    c("isoform", "dge")
  }

  dge <- NULL
  func <- NULL
  isoform <- NULL

  for (step in pipeline_steps) {
    if (step == "dge" && p$run_dge) {
      dge <- .stage_dge(ctx, tx_data)
      .stage_dge_plots(ctx, tx_data, dge)
      func <- .stage_functional(ctx, dge)
      .stage_fgsea(ctx, dge)
    }

    if (step == "isoform" && p$run_isoform) {
      isoform <- .stage_isoform(ctx)
    }
  }

  .stage_save_rdata(ctx, tx_data, dge, func, isoform)

  if (isTRUE(p$make_report)) {
    safe_run(.stage_report(ctx), label = "Interactive report")
  }

  message("=== Pipeline complete for: ", ctx$comp_name, " ===")

  invisible(NULL)
}

#' Run the full ExpressOM bulk RNA-seq pipeline
#' @rdname expressom
#' @export
run_bulk_pipeline <- expressom
