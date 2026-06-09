#' Run Expressom Pipeline
#'
#' @description
#' Executes the complete bulk RNA-seq pipeline including data import, exploratory
#' data analysis, differential expression (DESeq2), functional analysis, and
#' optional isoform-level analysis (DTE, DTU, IsoformSwitchAnalyzeR).
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
#' @param matrix_file Path to raw counts file if count_type = 'matrix'
#' @param gsea_metric Metric to rank genes for GSEA ("stat", "signed_pval", or "log2FoldChange")
#' @param subset_sample Optional string to filter the sample table (e.g., "cell_type == 'T_cells'")
#' @param remove_sample Optional character vector of sample IDs to exclude entirely
#' @param zscore_genes Optional character vector of gene names for targeted Z-score expression plotting
#' @param run_dge Logical: perform standard gene-level differential expression (default: TRUE)
#' @param run_isoform Logical: perform isoform-level analysis (DTE, DTU, IsoformSwitchAnalyzeR)
#' @param run_predictors Logical: Run CPAT, Pfam, and SignalP via WSL/Conda during Isoform analysis
#' @param bpparam BiocParallel backend for multi-threading (e.g. BiocParallel::MulticoreParam(4))
#' @param execution_order String: "dge_first" or "isoform_first" to prioritise which analysis runs first
#' @param isoform_fasta Path to transcript FASTA file (auto-downloads if NULL and run_isoform=TRUE)
#' @param isoform_gff Path to GFF/GTF annotation file (auto-downloads if NULL and run_isoform=TRUE)
#' @param use_wsl Logical: use WSL Ubuntu for external tools (Windows only)
#' @param wsl_distro Name of WSL distribution (default "Ubuntu-22.04")
#' @param resume_isoform_from Path to directory with saved DTE/DTU RDS files to resume isoform analysis
#' @param isoform_report_genes Gene symbols (e.g., c("TP53", "BCL2")) for transcript-proportion plots
#' @param nBest Number of top genes to include in RegionReport
#' @return NULL (invisibly)
expressom <- function(count_type        = "salmon",
                      data_dir          = "./data",
                      out_dir           = "./results",
                      sample_table      = "./sample_table.csv",
                      matrix_file       = NULL,
                      level             = "treated",
                      base              = "control",
                      model             = "~ condition",
                      batch_col         = NULL,
                      replicate_col     = NULL,
                      shrink_method     = "ashr",
                      ensembl_package_name = "EnsDb.Hsapiens.v107",
                      top_genes         = 30,
                      gmt_file          = c("C2", "C3", "C5", "C8"),
                      padj_cutoff       = 0.01,
                      go_pvalue_cutoff  = 0.05,
                      go_qvalue_cutoff  = 0.2,
                      test              = "Wald",
                      reduced           = NULL,
                      highlight_genes   = NULL,
                      subset_sample     = NULL,
                      remove_sample     = NULL,
                      gsea_metric       = "stat",
                      nBest             = 20000,
                      zscore_genes      = NULL,
                      run_dge           = TRUE,
                      run_isoform       = FALSE,
                      run_predictors    = FALSE,
                      bpparam           = BiocParallel::bpparam(),
                      execution_order   = c("dge_first", "isoform_first"),
                      isoform_fasta     = NULL,
                      isoform_gff       = NULL,
                      use_wsl           = FALSE,
                      wsl_distro        = "Ubuntu-22.04",
                      resume_isoform_from = NULL,
                      isoform_report_genes = NULL) {

  execution_order <- match.arg(execution_order)
  old_warn        <- getOption("nwarnings")
  options(nwarnings = 10000)
  on.exit({
    options(nwarnings = old_warn)
    log_dir <- file.path(out_dir, "DE_raw_results", "Log")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    if (dir.exists(log_dir)) {
      capture.output(sessionInfo(), file = file.path(log_dir, "SessionInfo.txt"))
      if (!is.null(warnings())) {
        capture.output(warnings(), file = file.path(log_dir, "Warnings.txt"))
      }
    }
  }, add = TRUE)

  # 1. Automated Ensembl DB Installation Check
  if (!requireNamespace(ensembl_package_name, quietly = TRUE)) {
    message("Package '", ensembl_package_name, "' not found. Attempting automatic build and installation...")
    matches     <- regexec("^EnsDb\\.(Hsapiens|Mmusculus)\\.v([0-9]+)$", ensembl_package_name)
    match_parts <- regmatches(ensembl_package_name, matches)[[1]]
    if (length(match_parts) == 3) {
      species <- if (match_parts[2] == "Hsapiens") "human" else "mouse"
      release <- match_parts[3]
      create_homemade_db(species = species, release = release)
      install_internal_db(pkg_name = ensembl_package_name)
    } else {
      stop("Could not parse ensembl_package_name '", ensembl_package_name,
           "' to build database automatically. Expected format like 'EnsDb.Hsapiens.v107'")
    }
  }

  comp_name      <- paste0(level, "_vs_", base)
  main_condition <- tail(all.vars(as.formula(model)), 1)
  edb_obj        <- getExportedValue(ensembl_package_name, ensembl_package_name)

  if (!is.null(batch_col) && !(batch_col %in% all.vars(as.formula(model)))) {
    message("WARNING: 'batch_col' [", batch_col, "] is specified for visualization, but missing from your design formula: '", model, "'.")
    message("Consider adding it to prevent confounding your DGE results (e.g., model = '~ ", batch_col, " + ", main_condition, "').")
  }

  pipeline_steps  <- if (execution_order == "dge_first") c("dge", "isoform") else c("isoform", "dge")
  dge_results     <- NULL
  isoform_results <- NULL

  for (step in pipeline_steps) {

    # =========================================================================
    # GENE-LEVEL DIFFERENTIAL EXPRESSION (DGE) BLOCK
    # =========================================================================
    if (step == "dge" && run_dge) {
      message("\n=== Running Gene-Level DGE Analysis ===")

      tx_data <- import_counts(
        data_dir             = data_dir,
        sample_table         = sample_table,
        ensembl_package_name = ensembl_package_name,
        count_type           = count_type,
        out_dir              = out_dir,
        matrix_file          = matrix_file,
        subset_sample        = subset_sample,
        remove_sample        = remove_sample
      )

      dds <- create_dds_object(tx_data, level, base, model, replicate_col)

      run_eda(dds            = dds,
              edb            = edb_obj,
              out_dir        = out_dir,
              level          = level,
              base           = base,
              main_condition = main_condition)

      res_list <- run_deseq2_analysis(dds, model, level, base, shrink_method, out_dir,
                                      padj_cutoff, test, reduced)

      results_data <- export_significant_results(
        res_shrunken   = res_list$res_shrunken,
        res_unshrunken = res_list$res_unshrunken,
        dds            = res_list$dds,
        out_dir        = out_dir,
        level          = level,
        base           = base,
        gene_map       = tx_data$gene_map,
        padj_cutoff    = padj_cutoff
      )

      # Build a symbol-keyed dds/res for RegionReport (which needs stable rownames)
      message("Converting identifiers for RegionReport...")
      dds_rep   <- res_list$dds
      res_rep   <- res_list$res_shrunken

      clean_ens <- gsub("\\..*$", "", rownames(dds_rep))
      sym_map   <- tx_data$gene_map$symbol[match(clean_ens, tx_data$gene_map$ensembl)]
      sym_map[is.na(sym_map) | sym_map == ""] <- clean_ens[is.na(sym_map) | sym_map == ""]
      sym_map   <- make.unique(sym_map)
      rownames(dds_rep) <- sym_map
      rownames(res_rep) <- sym_map

      generate_bulk_visualizations(
        dds            = res_list$dds,
        edb            = edb_obj,
        res_shrunken   = res_list$res_shrunken,
        res_unshrunken = res_list$res_unshrunken,
        results_data   = results_data,
        out_dir        = out_dir,
        level          = level,
        base           = base,
        main_condition = main_condition,
        top_genes      = top_genes,
        padj_cutoff    = padj_cutoff,
        highlight_genes = highlight_genes
      )
      while (grDevices::dev.cur() > 1) grDevices::dev.off()

      # --- RegionReport ---
      report_dir <- file.path(out_dir, "RegionReport")
      if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

      # Relative paths from the report output directory to the Plots directory
      plot_dir_rel  <- "../Plots"
      pca_file      <- file.path(plot_dir_rel, paste0("PCA_",              level, "_vs_", base, ".pdf"))
      pca_corr_file <- file.path(plot_dir_rel, paste0("PCA_BatchCorrected_", level, "_vs_", base, ".pdf"))
      heatmap_file  <- file.path(plot_dir_rel, paste0("SampleCorrelation_", level, "_vs_", base, ".pdf"))
      volcano_file  <- file.path(plot_dir_rel, paste0("DE_Volcanoplot_",   level, "_vs_", base, ".pdf"))
      ma_file       <- file.path(plot_dir_rel, paste0("MAplot_shrunken_",  level, "_vs_", base, ".pdf"))

      # BUG FIX 1: extension must be .Rmd, not .R — regionReport::DESeq2Report()
      #            expects an Rmd fragment for its customCode argument.
      # BUG FIX 2: each chunk header must be a single string; previously the code
      #            used c("...'", path_var, "'...") which split the header across
      #            three separate lines, producing invalid Rmd.  Use paste0() so
      #            each element of the writeLines vector is one complete line.
      custom_script <- tempfile(fileext = ".Rmd")
      writeLines(c(
        paste0("```{r pca_orig, echo=FALSE, fig.cap='PCA plot (Original)', eval=file.exists('", pca_file, "')}"),
        paste0("  knitr::include_graphics('", pca_file, "')"),
        "```",
        "",
        paste0("```{r pca_corr, echo=FALSE, fig.cap='PCA plot (Batch Corrected via limma)', eval=file.exists('", pca_corr_file, "')}"),
        paste0("  if (file.exists('", pca_corr_file, "')) knitr::include_graphics('", pca_corr_file, "')"),
        "```",
        "",
        paste0("```{r heatmap, echo=FALSE, fig.cap='Sample Correlation Heatmap', eval=file.exists('", heatmap_file, "')}"),
        paste0("  knitr::include_graphics('", heatmap_file, "')"),
        "```",
        "",
        paste0("```{r volcano, echo=FALSE, fig.cap='Volcano plot', eval=file.exists('", volcano_file, "')}"),
        paste0("  knitr::include_graphics('", volcano_file, "')"),
        "```",
        "",
        paste0("```{r maplot, echo=FALSE, fig.cap='MA plot (shrunken)', eval=file.exists('", ma_file, "')}"),
        paste0("  knitr::include_graphics('", ma_file, "')"),
        "```"
      ), custom_script)

      if (requireNamespace("regionReport", quietly = TRUE)) {
        regionReport::DESeq2Report(
          dds        = dds_rep,
          res        = res_rep,
          project    = comp_name,
          intgroup   = unique(c(main_condition, batch_col)),
          outdir     = report_dir,
          output     = paste0("RegionReport_", comp_name),
          nBest      = nBest,
          customCode = custom_script,
          echo       = FALSE
        )
      } else {
        message("Skipping DESeq2Report: 'regionReport' package is not installed.")
      }
      unlink(custom_script)

      while (grDevices::dev.cur() > 1) grDevices::dev.off()

      if (!is.null(zscore_genes)) {
        message("Generating targeted expression heatmaps (Sample Z-score & L2FC)...")
        plot_dir <- file.path(out_dir, "Plots")
        plot_sample_zscore(dds_rep, zscore_genes, main_condition, level, base, plot_dir)
        plot_l2fc_heatmap(dds_rep, zscore_genes, main_condition, level, base, plot_dir)
      }

      func_results <- run_functional_analysis(
        res_tbl          = results_data$res_tbl,
        sig_res          = results_data$sig_res,
        edb              = edb_obj,
        out_dir          = out_dir,
        level            = level,
        base             = base,
        top_genes        = top_genes,
        padj_cutoff      = padj_cutoff,
        go_pvalue_cutoff = go_pvalue_cutoff,
        go_qvalue_cutoff = go_qvalue_cutoff,
        gsea_metric      = gsea_metric,
        test_type        = test
      )

      message("Running FGSEA Analysis...")
      run_fgsea_analysis(
        res_tbl     = results_data$res_tbl,
        gmt_file    = gmt_file,
        edb         = edb_obj,
        out_dir     = out_dir,
        comp_name   = comp_name,
        padj_cutoff = padj_cutoff
      )

      dge_results <- list(
        tx_data      = tx_data,
        dds          = dds,
        res_list     = res_list,
        results_data = results_data,
        func_results = func_results
      )
    }

    # =========================================================================
    # ISOFORM-LEVEL ANALYSIS BLOCK
    # =========================================================================
    if (step == "isoform" && run_isoform) {
      message("\n=== Running Isoform-Level Analysis ===")

      dte_res        <- NULL
      dtu_res        <- NULL
      isoform_import <- NULL

      if (!is.null(resume_isoform_from) && dir.exists(resume_isoform_from)) {
        message("Resuming isoform analysis from: ", resume_isoform_from)
        loaded <- load_isoform_results(resume_isoform_from)
        if (!is.null(loaded$dte_results) && !is.null(loaded$dtu_results)) {
          dte_res <- loaded$dte_results
          dtu_res <- loaded$dtu_results
          message("Successfully loaded DTE and DTU results from saved RDS files.")
        } else {
          warning("Resume directory missing DTE or DTU RDS files. Proceeding with full isoform analysis.")
          resume_isoform_from <- NULL
        }
      }

      if (is.null(dte_res) || is.null(dtu_res)) {
        message("Running full isoform analysis (DTE, DTU, and Switch).")

        if (is.null(isoform_fasta) || is.null(isoform_gff)) {
          message("   -> Auto-fetching Ensembl FASTA and GTF references...")
          ref_dir       <- safe_dir(file.path(data_dir, "reference"))
          refs          <- download_ensembl_refs(ensembl_package_name = ensembl_package_name, out_dir = ref_dir)
          isoform_gff   <- refs$gtf
          isoform_fasta <- c(refs$cdna_fasta, refs$ncrna_fasta)
        }

        if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
          stop("Please install IsoformSwitchAnalyzeR: BiocManager::install('IsoformSwitchAnalyzeR')")
        }

        isoform_import <- import_transcript_counts(
          data_dir             = data_dir,
          sample_table         = sample_table,
          ensembl_package_name = ensembl_package_name,
          count_type           = count_type,
          matrix_file          = matrix_file,
          subset_sample        = subset_sample,
          remove_sample        = remove_sample
        )

        dte_res <- run_dte(isoform_import, main_condition, level, base, padj_cutoff)
        dtu_res <- run_dtu(isoform_import, main_condition, level, base, bpparam = bpparam)
      }

      if (is.null(isoform_import)) {
        message("Re-importing isoform counts for switch analysis (resume mode).")
        isoform_import <- import_transcript_counts(
          data_dir             = data_dir,
          sample_table         = sample_table,
          ensembl_package_name = ensembl_package_name,
          count_type           = count_type,
          matrix_file          = matrix_file,
          subset_sample        = subset_sample,
          remove_sample        = remove_sample
        )
        if (is.null(isoform_fasta) || is.null(isoform_gff)) {
          message("   -> Auto-fetching Ensembl FASTA and GTF references...")
          ref_dir       <- safe_dir(file.path(data_dir, "reference"))
          refs          <- download_ensembl_refs(ensembl_package_name = ensembl_package_name, out_dir = ref_dir)
          isoform_gff   <- refs$gtf
          isoform_fasta <- c(refs$cdna_fasta, refs$ncrna_fasta)
        }
      }

      if (!is.null(dte_res) && !is.null(dtu_res) && !is.null(isoform_import)) {
        generate_dte_dtu_report(
          dte_results      = dte_res,
          dtu_results      = dtu_res,
          isoform_obj      = isoform_import,
          out_dir          = out_dir,
          condition        = main_condition,
          level            = level,
          base             = base,
          genes_of_interest = isoform_report_genes,
          top_n            = 15
        )
      }

      switch_res <- run_isoform_switch(
        dte_results    = dte_res,
        dtu_results    = dtu_res,
        isoform_obj    = isoform_import,
        condition      = main_condition,
        level          = level,
        base           = base,
        fasta_file     = isoform_fasta,
        gff_file       = isoform_gff,
        out_dir        = file.path(out_dir, "IsoformSwitch"),
        run_predictors = run_predictors,
        use_wsl        = use_wsl,
        wsl_distro     = wsl_distro,
        save_dir       = file.path(out_dir, "IsoformSwitch", "saved"),
        resume_from    = resume_isoform_from
      )

      iso_dir <- file.path(out_dir, "IsoformSwitch")
      if (!dir.exists(iso_dir)) dir.create(iso_dir, recursive = TRUE)
      saveRDS(dte_res,    file.path(iso_dir, "dte_results.rds"))
      saveRDS(dtu_res,    file.path(iso_dir, "dtu_results.rds"))
      saveRDS(switch_res, file.path(iso_dir, "switch_list.rds"))

      message("Isoform analysis complete. Results saved in: ", iso_dir)

      isoform_results <- list(
        isoform_import = isoform_import,
        dte_res        = dte_res,
        dtu_res        = dtu_res,
        switch_res     = switch_res
      )
    }
  }

  # Save workspace
  rdata_dir  <- safe_dir(file.path(out_dir, "Save_rdata"))
  rdata_path <- file.path(rdata_dir, paste0("Results_", comp_name, ".RData"))

  to_save <- c("comp_name", "edb_obj", "main_condition", "level", "base", "model", "test", "batch_col")
  if (run_dge)      to_save <- c(to_save, "dge_results")
  if (run_isoform)  to_save <- c(to_save, "isoform_results")
  for (nm in c("tx_data", "dds", "res_list", "results_data", "func_results",
               "isoform_import", "dte_res", "dtu_res", "switch_res")) {
    if (exists(nm, envir = environment())) to_save <- c(to_save, nm)
  }

  save(list  = intersect(to_save, ls(envir = environment(), all.names = TRUE)),
       file  = rdata_path,
       envir = environment())
  message("\nR environment saved to: ", rdata_path)
  message("=== Pipeline complete for: ", comp_name, " ===")
  invisible(NULL)
}