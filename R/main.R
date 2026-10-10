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
#' @param nBest Number of top genes to include in RegionReport
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
  nBest = 20000,
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
  deg_lfc_cutoff = 1
) {
  execution_order <- match.arg(execution_order)
  isoform_test_engine <- match.arg(isoform_test_engine)
  requested_nbest <- nBest
  regionreport_nbest <- if (isTRUE(run_dge)) .validate_nbest(nBest) else NULL

  # Validate the plot / summary options up front, before any heavy computation
  sample_labels <- .read_sample_labels(sample_labels)
  deg_padj_cutoffs <- .validate_deg_cutoffs(deg_padj_cutoffs)
  invisible(.pca_colors(pca_colors))
  invisible(.volcano_colors(volcano_colors))

  if (
    !is.numeric(deg_lfc_cutoff) ||
      length(deg_lfc_cutoff) != 1 ||
      is.na(deg_lfc_cutoff) ||
      deg_lfc_cutoff < 0
  ) {
    stop("`deg_lfc_cutoff` must be a single non-negative number (e.g. 1).")
  }

  validate_environment(
    run_isoform = run_isoform,
    run_functional = run_dge
  )

  if (isTRUE(run_predictors) && .Platform$OS.type == "windows" && !use_wsl) {
    warning(
      "run_predictors = TRUE but use_wsl = FALSE on Windows. ",
      "External tools (CPAT, SignalP, Pfam) normally run inside WSL on Windows; ",
      "set use_wsl = TRUE, or ensure bash + the tools are directly reachable ",
      "(e.g. Git-Bash/MSYS2) if you intend to run them natively."
    )
  }

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

  if (!requireNamespace(ensembl_package_name, quietly = TRUE)) {
    message(
      "Package '",
      ensembl_package_name,
      "' not found. Attempting automatic build and installation..."
    )

    parsed_ensembl <- .parse_ensembl_package_name(ensembl_package_name)
    archive_dir <- tempfile("expressom_ensdb_")
    if (!dir.create(archive_dir)) {
      stop(
        "Could not create temporary Ensembl archive directory.",
        call. = FALSE
      )
    }
    on.exit(unlink(archive_dir, recursive = TRUE), add = TRUE)
    archive_path <- create_homemade_db(
      species = parsed_ensembl$species,
      release = parsed_ensembl$release,
      output_dir = archive_dir
    )
    install_internal_db(
      pkg_name = ensembl_package_name,
      archive_path = archive_path
    )
  }

  # Read the sample table once so .resolve_main_condition() can look at which
  # column actually contains level/base. Wrapped in tryCatch so a missing or
  # malformed file still surfaces later with a normal, clearer error from
  # import_counts() rather than here.
  sample_meta_for_cond <- tryCatch(
    data.table::fread(sample_table, header = TRUE, data.table = FALSE),
    error = function(e) NULL
  )

  if (isTRUE(eda_only) || is.null(model) || is.null(level) || is.null(base)) {
    message(
      "DESeq2 model not fully specified or eda_only = TRUE. ",
      "Running only data import and exploratory analysis (EDA)."
    )

    run_dge <- FALSE
    run_isoform <- FALSE
    model_eda <- "~1"

    sample_cols <- if (!is.null(sample_meta_for_cond)) {
      colnames(sample_meta_for_cond)
    } else {
      tryCatch(
        colnames(utils::read.csv(sample_table, nrows = 1, check.names = FALSE)),
        error = function(e) character(0)
      )
    }

    if (!is.null(group_col) && group_col %in% sample_cols) {
      main_condition <- group_col
    } else if (!is.null(model) && model != "~1") {
      main_condition <- .resolve_main_condition(
        stats::as.formula(model),
        sample_meta_for_cond,
        level,
        base
      )
    } else {
      main_condition <- NULL
    }
  } else {
    model_eda <- model
    main_condition <- .resolve_main_condition(
      stats::as.formula(model),
      sample_meta_for_cond,
      level,
      base
    )
  }

  # NOTE: the environment check itself (debug_wsl()) now lives solely in
  # run_isoform_switch() -- it used to also run here first, unconditionally
  # and with its result discarded, duplicating every WSL probe command
  # (command -v cpat, signalp6, hmmscan, interproscan.sh, ...) for no
  # benefit. run_isoform_switch() performs the same check plus validates
  # wsl_available/conda_env_exists/pfam_db_indexed against it and logs it,
  # so it's the single source of truth now.

  comp_name <- if (!is.null(level) && !is.null(base)) {
    paste0(level, "vs", base)
  } else {
    "EDA"
  }

  edb_obj <- getExportedValue(ensembl_package_name, ensembl_package_name)

  if (!is.null(batch_col) && !is.null(model)) {
    model_vars <- all.vars(stats::as.formula(model))

    if (!(batch_col %in% model_vars)) {
      message(
        "WARNING: 'batch_col' [",
        batch_col,
        "] is specified for visualization, but missing from your design formula: '",
        model,
        "'."
      )

      message(
        "Consider adding it to prevent confounding your DGE results (e.g., model = '~ ",
        batch_col,
        " + ",
        main_condition,
        "')."
      )
    }
  }

  ## ------------------------------------------------------------------
  ## Always import counts and run EDA
  ## ------------------------------------------------------------------
  tx_data <- import_counts(
    data_dir = data_dir,
    sample_table = sample_table,
    ensembl_package_name = ensembl_package_name,
    count_type = count_type,
    out_dir = out_dir,
    matrix_file = matrix_file,
    subset_sample = subset_sample,
    remove_sample = remove_sample,
    custom_tx2gene = custom_tx2gene,
    custom_gene_map = custom_gene_map
  )

  dds <- create_dds_object(
    tx_data,
    level = if (run_dge) level else NULL,
    base = if (run_dge) base else NULL,
    model = model_eda,
    replicate_col = replicate_col
  )

  run_eda(
    dds = dds,
    edb = edb_obj,
    out_dir = out_dir,
    level = level,
    base = base,
    main_condition = main_condition,
    group_col = group_col,
    batch_col = batch_col,
    pca_ntop = pca_ntop,
    pca_colors = pca_colors,
    sample_labels = sample_labels
  )

  if (!run_dge && !run_isoform) {
    message(
      "EDA completed. Skipping differential expression and isoform analyses."
    )
    return(invisible(NULL))
  }

  pipeline_steps <- if (execution_order == "dge_first") {
    c("dge", "isoform")
  } else {
    c("isoform", "dge")
  }

  dge_results <- NULL
  isoform_results <- NULL

  for (step in pipeline_steps) {
    if (step == "dge" && run_dge) {
      message("\n=== Running Gene-Level DGE Analysis ===")

      dds <- create_dds_object(
        tx_data,
        level,
        base,
        model,
        replicate_col
      )

      res_list <- run_deseq2_analysis(
        dds,
        model,
        level,
        base,
        shrink_method,
        out_dir,
        padj_cutoff,
        test,
        reduced,
        lfc_cutoff = deg_lfc_cutoff
      )

      # keep the fitted object (the one created above has no size factors /
      # dispersions / results), so it is what ends up in dge_results and the .RData
      dds <- res_list$dds

      results_data <- export_significant_results(
        res_shrunken = res_list$res_shrunken,
        res_unshrunken = res_list$res_unshrunken,
        dds = res_list$dds,
        out_dir = out_dir,
        level = level,
        base = base,
        gene_map = tx_data$gene_map,
        padj_cutoff = padj_cutoff,
        padj_cutoffs = deg_padj_cutoffs,
        lfc_cutoff = deg_lfc_cutoff,
        test_type = test
      )

      message("Converting identifiers for RegionReport...")

      dds_rep <- res_list$dds
      res_rep <- res_list$res_shrunken

      clean_ens <- strip_ensembl_version(rownames(dds_rep))

      sym_map <- tx_data$gene_map$symbol[
        match(clean_ens, tx_data$gene_map$ensembl)
      ]

      sym_map[is.na(sym_map) | sym_map == ""] <- clean_ens[
        is.na(sym_map) | sym_map == ""
      ]

      sym_map <- make.unique(sym_map)

      rownames(dds_rep) <- sym_map
      rownames(res_rep) <- sym_map

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
        highlight_genes = highlight_genes,
        batch_col = batch_col,
        pca_ntop = pca_ntop,
        pca_colors = pca_colors,
        sample_labels = sample_labels,
        volcano_colors = volcano_colors,
        lfc_cutoff = deg_lfc_cutoff
      )

      while (grDevices::dev.cur() > 1) {
        grDevices::dev.off()
      }

      report_dir <- file.path(out_dir, "RegionReport")

      if (!dir.exists(report_dir)) {
        dir.create(report_dir, recursive = TRUE)
      }

      plot_dir_abs <- file.path(out_dir, "Plots")

      .find_plot <- function(prefix) {
        if (!dir.exists(plot_dir_abs)) {
          return("")
        }

        files <- list.files(
          plot_dir_abs,
          pattern = "\\.pdf$",
          ignore.case = TRUE
        )

        if (length(files) == 0) {
          return("")
        }

        target <- gsub(
          "[^A-Za-z0-9]",
          "",
          paste0(prefix, "_", level, "vs", base)
        )

        clean_files <- gsub(
          "[^A-Za-z0-9]",
          "",
          tools::file_path_sans_ext(files)
        )

        idx <- which(clean_files == target)

        if (length(idx) == 0) {
          idx <- which(startsWith(clean_files, target))
        }

        if (length(idx) > 0) {
          return(normalizePath(
            file.path(plot_dir_abs, files[idx[1]]),
            winslash = "/",
            mustWork = FALSE
          ))
        }

        ""
      }

      pca_file <- .find_plot("PCA")
      pca_corr_file <- .find_plot("PCA_BatchCorrected")
      heatmap_file <- .find_plot("SampleCorrelation")

      if (!nzchar(heatmap_file)) {
        heatmap_file <- .find_plot("HeatMap")
      }

      volcano_file <- .find_plot("DE_Volcanoplot")
      ma_file <- .find_plot("MAplot_shrunken")

      if (exists("convert_pdf_to_png", mode = "function")) {
        for (v in c(
          "pca_file",
          "pca_corr_file",
          "heatmap_file",
          "volcano_file",
          "ma_file"
        )) {
          abs_pdf <- get(v)

          if (is.null(abs_pdf) || length(abs_pdf) == 0 || !nzchar(abs_pdf)) {
            next
          }

          if (file.exists(abs_pdf)) {
            png_abs <- convert_pdf_to_png(abs_pdf, dpi = 150)

            if (!is.null(png_abs)) {
              assign(v, png_abs)
            }
          }
        }
      }

      safe_template_value <- function(x) {
        if (is.null(x) || length(x) == 0) {
          return("")
        }

        x <- x[1]

        if (is.na(x)) {
          return("")
        }

        x <- as.character(x)

        x <- gsub("\\", "/", x, fixed = TRUE)
        x <- gsub('"', '\\\"', x, fixed = TRUE)

        x
      }

      custom_script <- tryCatch(
        .render_placeholder_template(
          "dge_regionreport_customcode.Rmd",
          values = list(
            PCA_FILE = safe_template_value(pca_file),
            PCA_CORR_FILE = safe_template_value(pca_corr_file),
            HEATMAP_FILE = safe_template_value(heatmap_file),
            VOLCANO_FILE = safe_template_value(volcano_file),
            MA_FILE = safe_template_value(ma_file)
          )
        ),
        error = function(e) {
          message(
            "Could not render RegionReport custom code template: ",
            conditionMessage(e)
          )
          NULL
        }
      )
      if (!is.null(custom_script)) {
        on.exit(unlink(custom_script), add = TRUE)
      }

      if (requireNamespace("regionReport", quietly = TRUE)) {
        if (requested_nbest > 1000L) {
          message(
            "   -> nBest = ",
            nBest,
            " is too large for RegionReport HTML output. ",
            "Capping RegionReport nBest to ",
            regionreport_nbest,
            ". ",
            "Use CSV result tables if you need the full gene list."
          )
        }

        intgroup <- unique(c(main_condition, batch_col))
        intgroup <- intgroup[!is.na(intgroup) & nzchar(intgroup)]

        tryCatch(
          {
            if (!is.null(custom_script)) {
              regionReport::DESeq2Report(
                dds = dds_rep,
                res = res_rep,
                project = comp_name,
                intgroup = intgroup,
                outdir = report_dir,
                output = paste0("RegionReport_", comp_name),
                nBest = regionreport_nbest,
                customCode = custom_script,
                device = "png",
                browse = FALSE,
                echo = FALSE
              )
            } else {
              regionReport::DESeq2Report(
                dds = dds_rep,
                res = res_rep,
                project = comp_name,
                intgroup = intgroup,
                outdir = report_dir,
                output = paste0("RegionReport_", comp_name),
                nBest = regionreport_nbest,
                device = "png",
                browse = FALSE,
                echo = FALSE
              )
            }

            report_html <- file.path(
              report_dir,
              paste0("RegionReport_", comp_name, ".html")
            )

            if (file.exists(report_html)) {
              report_size_mb <- file.info(report_html)$size / 1024^2

              message(
                "   -> RegionReport HTML written: ",
                report_html,
                " (",
                round(report_size_mb, 1),
                " MB)"
              )

              if (report_size_mb > 50) {
                message(
                  "   -> WARNING: RegionReport HTML is very large (>50 MB). ",
                  "Consider reducing nBest further, for example nBest = 200 or 500."
                )
              }
            }
          },
          error = function(e) {
            message("RegionReport generation failed: ", conditionMessage(e))
          }
        )
      } else {
        message(
          "Skipping DESeq2Report: 'regionReport' package is not installed."
        )
      }

      while (grDevices::dev.cur() > 1) {
        grDevices::dev.off()
      }

      if (!is.null(zscore_genes)) {
        message(
          "Generating targeted expression heatmaps (Sample Z-score & L2FC)..."
        )

        plot_dir <- file.path(out_dir, "Plots")

        plot_sample_zscore(
          dds_rep,
          zscore_genes,
          main_condition,
          level,
          base,
          plot_dir
        )

        plot_l2fc_heatmap(
          dds_rep,
          zscore_genes,
          main_condition,
          level,
          base,
          plot_dir
        )
      }

      if (!is.null(gene_sets_zscore)) {
        message("Generating separate Z-score plots for each gene set...")

        plot_dir <- file.path(out_dir, "Plots")
        if (!dir.exists(plot_dir)) {
          dir.create(plot_dir, recursive = TRUE)
        }

        for (set_name in names(gene_sets_zscore)) {
          single_set <- list(gene_sets_zscore[[set_name]])
          names(single_set) <- set_name

          message("   -> Plotting set: ", set_name)

          plot_geneset_zscore_avg(
            dds = dds_rep,
            gene_sets = single_set,
            condition_col = main_condition,
            level = level,
            base = base,
            plot_dir = plot_dir,
            set_name = set_name
          )

          plot_gene_zscore_individual(
            dds = dds_rep,
            gene_sets = single_set,
            condition_col = main_condition,
            level = level,
            base = base,
            plot_dir = plot_dir,
            set_name = set_name
          )
        }
      }

      func_results <- safe_run(
        run_functional_analysis(
          res_tbl = results_data$res_tbl,
          sig_res = results_data$sig_res,
          edb = edb_obj,
          out_dir = out_dir,
          level = level,
          base = base,
          top_genes = top_genes,
          padj_cutoff = padj_cutoff,
          go_pvalue_cutoff = go_pvalue_cutoff,
          go_qvalue_cutoff = go_qvalue_cutoff,
          gsea_metric = gsea_metric,
          test_type = test,
          ora_padj_cutoff = ora_padj_cutoff,
          ora_min_genes = ora_min_genes,
          ora_lfc_cutoff = ora_lfc_cutoff,
          run_spia = run_spia
        ),
        label = "Functional analysis"
      )

      message("Running FGSEA Analysis...")

      safe_run(
        run_fgsea_analysis(
          res_tbl = results_data$res_tbl,
          gmt_file = gmt_file,
          edb = edb_obj,
          out_dir = out_dir,
          comp_name = comp_name,
          padj_cutoff = padj_cutoff,
          gsea_metric = gsea_metric,
          test_type = test
        ),
        label = "FGSEA analysis"
      )

      dge_results <- list(
        tx_data = tx_data,
        dds = dds,
        res_list = res_list,
        results_data = results_data,
        func_results = func_results
      )
    }

    if (step == "isoform" && run_isoform) {
      message("\n=== Running Isoform-Level Analysis ===")

      iso_dir <- file.path(out_dir, "IsoformSwitch")
      iso_save_dir <- file.path(iso_dir, "saved")

      dir.create(iso_save_dir, recursive = TRUE, showWarnings = FALSE)

      .iso_ckpt_exists <- function(nm) {
        file.exists(file.path(iso_save_dir, nm))
      }

      .iso_ckpt_save <- function(obj, nm) {
        .checkpoint_save(
          obj,
          file.path(iso_save_dir, nm),
          iso_checkpoint_signature
        )
        message("Step checkpoint saved: ", nm)
        invisible(obj)
      }

      .iso_ckpt_load <- function(nm) {
        message("Resuming from step checkpoint: ", nm)
        .checkpoint_load(file.path(iso_save_dir, nm), iso_checkpoint_signature)
      }

      dte_res <- NULL
      dtu_res <- NULL
      dexseq_res <- NULL
      isoform_import <- NULL

      if (!is.null(resume_isoform_from) && dir.exists(resume_isoform_from)) {
        message("Resuming isoform analysis from: ", resume_isoform_from)

        loaded <- load_isoform_results(resume_isoform_from)

        isoform_import <- loaded$isoform_import
        dte_res <- loaded$dte_results
        dtu_res <- loaded$dtu_results
        dexseq_res <- loaded$dexseq_results

        if (!is.null(dte_res) && !is.null(dtu_res)) {
          message("Loaded DTE and DTU results from resume directory.")
        } else {
          message("Resume directory incomplete. Will run missing steps.")
        }
      }

      if (is.null(isoform_fasta) || is.null(isoform_gff)) {
        message("   -> Auto-fetching Ensembl FASTA and GTF references...")

        ref_dir <- safe_dir(file.path(data_dir, "reference"))

        refs <- download_ensembl_refs(
          ensembl_package_name = ensembl_package_name,
          out_dir = ref_dir
        )

        isoform_refs <- .resolve_isoform_references(
          isoform_fasta,
          isoform_gff,
          refs
        )
        isoform_fasta <- isoform_refs$fasta
        isoform_gff <- isoform_refs$gff
      }

      quant_files <- switch(
        count_type,
        salmon = "quant.sf",
        kallisto = "abundance.tsv",
        rsem = "quant.isoforms.results",
        stringtie = "t_data.ctab",
        character()
      )
      quant_paths <- if (length(quant_files) && dir.exists(data_dir)) {
        list.files(
          data_dir,
          pattern = paste0("^", quant_files, "$"),
          recursive = TRUE,
          full.names = TRUE
        )
      } else {
        character()
      }
      checkpoint_input_files <- c(
        sample_table,
        matrix_file,
        custom_tx2gene,
        custom_gene_map,
        if (!is.data.frame(custom_transcript_id_map)) custom_transcript_id_map,
        quant_paths,
        isoform_fasta,
        isoform_gff
      )
      iso_checkpoint_signature <- .object_md5(list(
        inputs = .file_fingerprint(checkpoint_input_files),
        custom_transcript_id_map = if (
          is.data.frame(custom_transcript_id_map)
        ) {
          custom_transcript_id_map
        } else {
          NULL
        },
        count_type = count_type,
        subset_sample = subset_sample,
        remove_sample = remove_sample,
        model = model,
        level = level,
        base = base,
        padj_cutoff = padj_cutoff,
        run_dexseq = run_dexseq,
        run_isoform_enrichment = run_isoform_enrichment,
        ora_padj_cutoff = ora_padj_cutoff,
        ora_min_genes = ora_min_genes,
        ora_lfc_cutoff = ora_lfc_cutoff,
        gsea_metric = gsea_metric,
        gmt_file = gmt_file,
        top_genes = top_genes,
        go_pvalue_cutoff = go_pvalue_cutoff,
        go_qvalue_cutoff = go_qvalue_cutoff,
        run_predictors = run_predictors,
        predictor_cpu = predictor_cpu,
        use_wsl = use_wsl,
        wsl_distro = wsl_distro,
        run_spia = run_spia,
        isoform_test_engine = isoform_test_engine
      ))

      if (is.null(isoform_import) && .iso_ckpt_exists("isoform_import.rds")) {
        isoform_import <- .iso_ckpt_load("isoform_import.rds")
      }
      if (is.null(dte_res) && .iso_ckpt_exists("dte_results.rds")) {
        dte_res <- .iso_ckpt_load("dte_results.rds")
      }
      if (is.null(dtu_res) && .iso_ckpt_exists("dtu_results.rds")) {
        dtu_res <- .iso_ckpt_load("dtu_results.rds")
      }
      if (is.null(dexseq_res) && .iso_ckpt_exists("dexseq_results.rds")) {
        dexseq_res <- .iso_ckpt_load("dexseq_results.rds")
      }

      if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
        stop(
          "Please install IsoformSwitchAnalyzeR: BiocManager::install('IsoformSwitchAnalyzeR')"
        )
      }

      if (is.null(isoform_import)) {
        message("Step A: Importing transcript-level counts...")

        isoform_import <- import_transcript_counts(
          data_dir = data_dir,
          sample_table = sample_table,
          ensembl_package_name = ensembl_package_name,
          count_type = count_type,
          matrix_file = matrix_file,
          subset_sample = subset_sample,
          remove_sample = remove_sample,
          custom_tx2gene = custom_tx2gene,
          custom_gene_map = custom_gene_map
        )

        .iso_ckpt_save(isoform_import, "isoform_import.rds")
      } else {
        message("Step A: Isoform import loaded from checkpoint. Skipping.")
      }

      safe_run(
        run_isoform_pca(
          isoform_import,
          main_condition,
          level,
          base,
          out_dir = iso_dir,
          batch_col = batch_col,
          pca_ntop = pca_ntop,
          save_dir = iso_save_dir,
          pca_colors = pca_colors,
          sample_labels = sample_labels
        ),
        label = "Transcript-level PCA"
      )

      if (is.null(dte_res)) {
        message("Step B: Running DTE (Differential Transcript Expression)...")

        dte_res <- run_dte(
          isoform_import,
          main_condition,
          level,
          base,
          padj_cutoff,
          bpparam = bpparam,
          design = model
        )

        .iso_ckpt_save(dte_res, "dte_results.rds")
      } else {
        message("Step B: DTE results loaded from checkpoint. Skipping.")
      }

      if (is.null(dtu_res)) {
        message("Step C: Running DTU (Differential Transcript Usage)...")

        dtu_res <- run_dtu(
          isoform_import,
          main_condition,
          level,
          base,
          bpparam = bpparam,
          design = model
        )

        .iso_ckpt_save(dtu_res, "dtu_results.rds")
      } else {
        message("Step C: DTU results loaded from checkpoint. Skipping.")
      }

      if (isTRUE(run_dexseq)) {
        if (is.null(dexseq_res)) {
          message(
            "Step C.5: Running DEXSeq-based DTU (complementary engine)..."
          )

          dexseq_res <- safe_run(
            run_dexseq_dtu(
              isoform_import,
              main_condition,
              level,
              base,
              bpparam = bpparam,
              design = model
            ),
            label = "DEXSeq DTU"
          )

          if (!is.null(dexseq_res)) {
            .iso_ckpt_save(dexseq_res, "dexseq_results.rds")
          }
        } else {
          message(
            "Step C.5: DEXSeq DTU results loaded from checkpoint. Skipping."
          )
        }
      }

      # Step D: functional enrichment on DTE/DTU results, using the exact
      # same ORA/FGSEA machinery and GO/Reactome/DO/MSigDB databases as the
      # DGE-level analysis (run_isoform_functional_analysis() in
      # mod_functional.R). Runs once DTE/DTU are available -- whether just
      # computed above or loaded from a checkpoint/resume directory -- and
      # deliberately before run_isoform_switch() below, since that's the
      # point in the pipeline the person asked for it.
      dte_func_res <- NULL
      dtu_func_res <- NULL

      if (isTRUE(run_isoform_enrichment)) {
        dte_func_checkpoint <- NULL
        if (.iso_ckpt_exists("dte_func_results.rds")) {
          dte_func_checkpoint <- .iso_ckpt_load("dte_func_results.rds")
        }
        if (!is.null(dte_func_checkpoint)) {
          dte_func_res <- dte_func_checkpoint
        } else if (!is.null(dte_res)) {
          message("Step D: Running functional enrichment on DTE results...")
          dte_func_res <- run_isoform_functional_analysis(
            res_tbl = dte_res,
            edb = edb_obj,
            out_dir = iso_dir,
            analysis_type = "DTE",
            level = level,
            base = base,
            top_genes = top_genes,
            padj_cutoff = padj_cutoff,
            go_pvalue_cutoff = go_pvalue_cutoff,
            go_qvalue_cutoff = go_qvalue_cutoff,
            gsea_metric = gsea_metric,
            ora_padj_cutoff = ora_padj_cutoff,
            ora_min_genes = ora_min_genes,
            ora_lfc_cutoff = ora_lfc_cutoff,
            gmt_file = gmt_file,
            run_spia = run_spia
          )
          .iso_ckpt_save(dte_func_res, "dte_func_results.rds")
        } else {
          message(
            "Step D: No DTE results available. Skipping DTE functional enrichment."
          )
        }

        dtu_func_checkpoint <- NULL
        if (.iso_ckpt_exists("dtu_func_results.rds")) {
          dtu_func_checkpoint <- .iso_ckpt_load("dtu_func_results.rds")
        }
        if (!is.null(dtu_func_checkpoint)) {
          dtu_func_res <- dtu_func_checkpoint
        } else if (!is.null(dtu_res) && !is.null(dtu_res$dtu_results)) {
          message("Step D: Running functional enrichment on DTU results...")
          dtu_func_res <- run_isoform_functional_analysis(
            res_tbl = dtu_res$dtu_results,
            edb = edb_obj,
            out_dir = iso_dir,
            analysis_type = "DTU",
            level = level,
            base = base,
            top_genes = top_genes,
            padj_cutoff = padj_cutoff,
            go_pvalue_cutoff = go_pvalue_cutoff,
            go_qvalue_cutoff = go_qvalue_cutoff,
            ora_padj_cutoff = ora_padj_cutoff,
            ora_min_genes = ora_min_genes,
            ora_lfc_cutoff = ora_lfc_cutoff
          )
          .iso_ckpt_save(dtu_func_res, "dtu_func_results.rds")
        } else {
          message(
            "Step D: No DTU results available. Skipping DTU functional enrichment."
          )
        }
      } else {
        message(
          "Step D: Skipping DTE/DTU functional enrichment (run_isoform_enrichment = FALSE)."
        )
      }

      .warn_isoform_switch_covariates(model, main_condition)

      switch_res <- run_isoform_switch(
        dte_results = dte_res,
        dtu_results = dtu_res,
        isoform_obj = isoform_import,
        condition = main_condition,
        level = level,
        base = base,
        fasta_file = isoform_fasta,
        gff_file = isoform_gff,
        out_dir = iso_dir,
        run_predictors = run_predictors,
        use_wsl = use_wsl,
        wsl_distro = wsl_distro,
        save_dir = iso_save_dir,
        resume_from = resume_isoform_from,
        predictor_cpu = predictor_cpu,
        log_dir = file.path(out_dir, "Log", "Isoform"),
        custom_transcript_id_map = custom_transcript_id_map,
        skip_fasta_filter = skip_fasta_filter,
        test_engine = isoform_test_engine,
        organism = get_organism_info(edb_obj)$name,
        plot_topology = plot_topology
      )

      if (!is.null(dte_res) && !is.null(dtu_res) && !is.null(isoform_import)) {
        generate_dte_dtu_report(
          dte_results = dte_res,
          dtu_results = dtu_res,
          isoform_obj = isoform_import,
          out_dir = out_dir,
          condition = main_condition,
          level = level,
          base = base,
          genes_of_interest = isoform_report_genes,
          top_n = 15,
          switch_list = switch_res,
          dexseq_results = dexseq_res,
          switch_plot_top_n = isoform_plot_top_n,
          plot_topology = plot_topology,
          padj_cutoff = padj_cutoff,
          design = model,
          dge_test = test,
          reduced = reduced
        )
      }

      if (!dir.exists(iso_dir)) {
        dir.create(iso_dir, recursive = TRUE)
      }

      message("Isoform analysis complete. Results saved in: ", iso_dir)

      while (grDevices::dev.cur() > 1) {
        grDevices::dev.off()
      }

      isoform_results <- list(
        isoform_import = isoform_import,
        dte_res = dte_res,
        dtu_res = dtu_res,
        dexseq_res = dexseq_res,
        dte_func_res = dte_func_res,
        dtu_func_res = dtu_func_res,
        switch_res = switch_res
      )
    }
  }

  rdata_dir <- safe_dir(file.path(out_dir, "Save_rdata"))
  rdata_path <- file.path(rdata_dir, paste0("Results_", comp_name, ".RData"))

  to_save <- c(
    "comp_name",
    "edb_obj",
    "main_condition",
    "level",
    "base",
    "model",
    "test",
    "batch_col"
  )

  if (run_dge) {
    to_save <- c(to_save, "dge_results")
  }

  if (run_isoform) {
    to_save <- c(to_save, "isoform_results")
  }

  for (nm in c(
    "tx_data",
    "dds",
    "res_list",
    "results_data",
    "func_results",
    "isoform_import",
    "dte_res",
    "dtu_res",
    "switch_res"
  )) {
    if (exists(nm, envir = environment())) {
      to_save <- c(to_save, nm)
    }
  }

  save(
    list = intersect(to_save, ls(envir = environment(), all.names = TRUE)),
    file = rdata_path,
    envir = environment()
  )

  message("\nR environment saved to: ", rdata_path)
  message("=== Pipeline complete for: ", comp_name, " ===")

  invisible(NULL)
}

#' Run the full ExpressOM bulk RNA-seq pipeline
#' @rdname expressom
#' @export
run_bulk_pipeline <- expressom
