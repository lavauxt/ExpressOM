#' Run Bulk RNA-Seq Pipeline
#' 
#' @description
#' Executes the complete bulk RNA-seq pipeline including data import, exploratory
#' data analysis, differential expression (DESeq2), and functional analysis.
#' 
#' @export
#' @param count_type Type of RNA count (default: "salmon")
#' @param data_dir Folder where the input data is stored
#' @param out_dir Folder to save results
#' @param sample_table Path to the sample table metadata
#' @param level Foreground level for contrast
#' @param base Background level for contrast
#' @param model Design formula for DESeq2 as a string
#' @param replicate_col Optional column specifying replicates
#' @param shrink_method Type of shrinkage (e.g., "ashr")
#' @param ensembl_package_name The Ensembl DB package
#' @param top_genes Number of top genes to plot
#' @param gmt_file Path to local GMT file. If NULL, msigdbr downloads Hallmark.
#' @param padj_cutoff Adjusted p-value threshold for significance
#' @param test Type of statistical test ("Wald" or "LRT")
#' @param reduced Design formula for the reduced model (used if test = "LRT")
#' @param highlight_genes Optional character vector of gene names to highlight in the Volcano plot
#' @param go_pvalue_cutoff GO ORA p-value cutoff
#' @param go_qvalue_cutoff GO ORA q-value cutoff
#' @param go_pvalue_cutoff GO ORA p-value cutoff
#' @param go_qvalue_cutoff GO ORA q-value cutoff
#' @param matrix_file Path to raw counts file if count_type = 'matrix'
#' @return NULL
run_bulk_pipeline <- function(count_type = "salmon", 
                              data_dir = "./data", 
                              out_dir = "./results", 
                              sample_table = "./sample_table.csv",
                              matrix_file = NULL,
                              level = "treated", 
                              base = "control",
                              model = "~ condition", 
                              replicate_col = NULL,
                              shrink_method = "ashr", 
                              ensembl_package_name = "EnsDb.Hsapiens.v107", 
                              top_genes = 30, 
                              gmt_file = NULL,
                              padj_cutoff = 0.01,
                              go_pvalue_cutoff = 0.05,
                              go_qvalue_cutoff = 0.2,
                              test = "Wald",
                              reduced = NULL,
                              highlight_genes = NULL) {
  
  # Ensure we capture all warnings for the final log
  old_warn <- getOption("nwarnings")
  options(nwarnings = 10000)
  on.exit({
    options(nwarnings = old_warn)
    log_dir <- file.path(out_dir, "Log")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    if (dir.exists(log_dir)) {
      capture.output(sessionInfo(), file = file.path(log_dir, "SessionInfo.txt"))
      if (!is.null(warnings())) {
        capture.output(warnings(), file = file.path(log_dir, "Warnings.txt"))
      }
    }
  }, add = TRUE)
  
  comp_name <- paste0(level, "_vs_", base)
  # Extract the main condition variable from the model (assuming the last term is the primary group)
  main_condition <- tail(all.vars(as.formula(model)), 1)

  # ---------------------------------------------------------
  # 2. Import Data 
  # ---------------------------------------------------------
  tx_data <- import_counts(
    data_dir = data_dir, 
    sample_table = sample_table, 
    ensembl_package_name = ensembl_package_name, 
    count_type = count_type, 
    out_dir = out_dir,
    matrix_file = matrix_file
  )
  
  # ---------------------------------------------------------
  # 3. Create DESeq2 Object
  # ---------------------------------------------------------
  dds <- create_dds_object(tx_data, level, base, model, replicate_col)
  
  # ---------------------------------------------------------
  # 4. Exploratory Data Analysis & QC
  # ---------------------------------------------------------
  edb_obj <- getExportedValue(ensembl_package_name, ensembl_package_name)
  
  run_eda(
    dds = dds, 
    edb = edb_obj, 
    out_dir = out_dir, 
    level = level, 
    base = base, 
    main_condition = main_condition
  )
  
  # ---------------------------------------------------------
  # 5. Differential Expression
  # ---------------------------------------------------------
    res_list <- run_deseq2_analysis(
    dds, model, level, base, shrink_method, out_dir, padj_cutoff, test, reduced
  )
 
  # ---------------------------------------------------------
  # 6. Extract & Save Tables 
  # ---------------------------------------------------------
  results_data <- export_significant_results(
    res_shrunken = res_list$res_shrunken, 
    res_unshrunken = res_list$res_unshrunken, 
    dds = res_list$dds, 
    out_dir = out_dir, 
    level = level, 
    base = base,
    gene_map = tx_data$gene_map,
    padj_cutoff = padj_cutoff
  )
  
  # ---------------------------------------------------------
  # 7. Reports & Visualizations
  # ---------------------------------------------------------
  message("Converting identifiers for RegionReport...")
  dds_rep <- res_list$dds
  res_rep <- res_list$res_shrunken
  
  clean_ens <- gsub("\\..*$", "", rownames(dds_rep))
  sym_map <- tx_data$gene_map$symbol[match(clean_ens, tx_data$gene_map$ensembl)]
  sym_map[is.na(sym_map) | sym_map == ""] <- clean_ens[is.na(sym_map) | sym_map == ""]
  
  rownames(dds_rep) <- sym_map
  rownames(res_rep) <- sym_map
  
  if (requireNamespace("regionReport", quietly = TRUE)) {
    report_dir <- file.path(out_dir, "Reports", "RegionReport")
    if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
    regionReport::DESeq2Report(
      dds = dds_rep, 
      res = res_rep, 
      project = comp_name, 
      intgroup = main_condition, 
      outdir = report_dir, 
      output = paste0("RegionReport_", comp_name),
      nBest = 10000 
    )
  } else {
    message("Skipping DESeq2Report: 'regionReport' package is not installed.")
  }
  
  
  while (dev.cur() > 1) dev.off()
  
   
  generate_bulk_visualizations(
    dds = res_list$dds, 
    edb = edb_obj, 
    res_shrunken = res_list$res_shrunken, 
    res_unshrunken = res_list$res_unshrunken, 
    results_data = results_data, 
    out_dir = out_dir, 
    level = level, 
    base = base, 
    main_condition = main_condition, 
    top_genes = top_genes, 
    padj_cutoff = padj_cutoff,
    highlight_genes = highlight_genes
  )
  
  while (dev.cur() > 1) dev.off()

 # ---------------------------------------------------------
  # 8. Functional Analysis
  # ---------------------------------------------------------
  func_results <- run_functional_analysis(
    res_tbl = results_data$res_tbl, sig_res = results_data$sig_res,
    edb = edb_obj, out_dir = out_dir, level = level, base = base,
    top_genes = top_genes, padj_cutoff = padj_cutoff,
    go_pvalue_cutoff = go_pvalue_cutoff,
    go_qvalue_cutoff = go_qvalue_cutoff
  )

  # ---------------------------------------------------------
  # 9. FGSEA Analysis
  # ---------------------------------------------------------
  message("Running FGSEA Analysis...")
  run_fgsea_analysis(
    res_tbl = results_data$res_tbl, gmt_file = gmt_file,
    edb = edb_obj, out_dir = out_dir, comp_name = comp_name,
    padj_cutoff = padj_cutoff
  )

  # ---------------------------------------------------------
  # 10. Save Environment
  # ---------------------------------------------------------
  rdata_dir  <- safe_dir(file.path(out_dir, "Save_rdata"))
  rdata_path <- file.path(rdata_dir, paste0("Results_", comp_name, ".RData"))
  save(list = ls(all.names = TRUE), file = rdata_path)
  message("R environment saved to: ", rdata_path)

  # ---------------------------------------------------------
  # 11. Generate & Render Markdown Report
  # ---------------------------------------------------------
  reports_dir <- safe_dir(file.path(out_dir, "Reports"))
  rmd_file    <- file.path(reports_dir, paste0("Pipeline_Summary_", comp_name, ".Rmd"))

  rmd_lines <- c(
    "---",
    paste0("title: 'RNA-Seq Pipeline Summary: ", level, " vs ", base, "'"),
    "date: '`r Sys.Date()`'",
    "output: html_document",
    "---",
    "## Differential Expression Summary",
    "```{r deg-summary, echo=FALSE}",
    paste0("cat('Total DEGs (padj < ", padj_cutoff, "):', nrow(results_data$sig_res), '\\n')"),
    "cat('Upregulated (LFC > 0):', sum(results_data$sig_res$log2FoldChange > 0, na.rm=TRUE), '\\n')",
    "cat('Downregulated (LFC < 0):', sum(results_data$sig_res$log2FoldChange < 0, na.rm=TRUE), '\\n')",
    "```",
    "",
    "## Top 10 Significant DEGs",
    "```{r top-genes, echo=FALSE}",
    "if (nrow(results_data$sig_res) > 0) {",
    "  knitr::kable(head(results_data$sig_res[order(results_data$sig_res$padj), c('gene','log2FoldChange','padj')], 10),",
    "    row.names = FALSE, digits = 3, format.args = list(scientific = TRUE))",
    "} else { cat('No significant genes found.') }",
    "```"
  )
  writeLines(rmd_lines, rmd_file)

  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    safe_run(
      rmarkdown::render(rmd_file, output_dir = reports_dir, envir = environment(), quiet = TRUE),
      label = "Rmd render"
    )
  } else {
    message("Skipping Rmd render: 'rmarkdown' is not installed.")
  }

  message("=== Pipeline complete for: ", comp_name, " ===")
  invisible(NULL)
}