# mod_stages.R -- the ExpressOM pipeline cut into stages.
#
# expressom() (sequential) and the `targets` pipeline (cached, re-runnable
# stage by stage; see mod_targets.R) both call exactly these functions, so a
# run launched from the console, from `targets` or from the Shiny app executes
# the same code.

# ---- helpers -----------------------------------------------------------------------------

.out_manifest <- function(out_dir, subdirs = NULL) {
  empty <- data.frame(
    path = character(0), size = numeric(0), mtime = numeric(0),
    stringsAsFactors = FALSE
  )
  roots <- if (is.null(subdirs)) out_dir else file.path(out_dir, subdirs)
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(empty)
  f <- unique(unlist(lapply(
    roots, list.files, recursive = TRUE, full.names = TRUE
  )))
  if (!length(f)) return(empty)
  info <- file.info(f)
  out <- data.frame(
    path = .rep_rel(f, out_dir),
    size = info$size,
    mtime = round(as.numeric(info$mtime)),
    stringsAsFactors = FALSE
  )
  out[order(out$path), , drop = FALSE]
}

.stage_edb <- function(ctx) {
  pkg <- ctx$p$ensembl_package_name
  getExportedValue(pkg, pkg)
}

.close_devices <- function() {
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  invisible(NULL)
}

.quant_input_files <- function(p, gene_level = FALSE) {
  quant_files <- switch(
    p$count_type,
    salmon = "quant.sf",
    kallisto = "abundance.tsv",
    rsem = if (gene_level) "quant.genes.results" else "quant.isoforms.results",
    stringtie = "t_data.ctab",
    character()
  )
  if (length(quant_files) && !is.null(p$data_dir) && dir.exists(p$data_dir)) {
    list.files(
      p$data_dir,
      pattern = paste0("^", quant_files, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
  } else {
    character()
  }
}

# ---- stage 0: validation, environment, contrast resolution -----------------------------------

#' Validate the arguments of expressom() and resolve the contrast
#' @keywords internal
.expressom_prepare <- function(p) {
  p$execution_order <- match.arg(
    p$execution_order, c("dge_first", "isoform_first")
  )
  p$isoform_test_engine <- match.arg(
    p$isoform_test_engine, c("DEXSeq", "DRIMSeq", "satuRn")
  )

  # All WSL/predictor tool logs land in <out_dir>/Log/WSL/ so Log/ remains
  # the single folder holding every run diagnostic. Restored on exit so we
  # don't pollute the caller's session.
  old_log_opt <- getOption("expressom.wsl_log_dir")
  options(expressom.wsl_log_dir = file.path(p$out_dir, "Log", "WSL"))
  on.exit(options(expressom.wsl_log_dir = old_log_opt), add = TRUE)

  p$sample_labels <- .read_sample_labels(p$sample_labels)
  p$deg_padj_cutoffs <- .validate_deg_cutoffs(p$deg_padj_cutoffs)
  invisible(.pca_colors(p$pca_colors))
  invisible(.volcano_colors(p$volcano_colors))

  if (
    !is.numeric(p$deg_lfc_cutoff) ||
      length(p$deg_lfc_cutoff) != 1 ||
      is.na(p$deg_lfc_cutoff) ||
      p$deg_lfc_cutoff < 0
  ) {
    stop("`deg_lfc_cutoff` must be a single non-negative number (e.g. 1).")
  }

  if (!is.null(p$nBest)) {
    warning(
      "`nBest` is deprecated and ignored: the RegionReport report was replaced ",
      "by the unified interactive report (see `report_settings`).",
      call. = FALSE
    )
    p$nBest <- NULL
  }

  validate_environment(
    run_isoform    = p$run_isoform,
    run_functional = p$run_dge,
    run_report     = isTRUE(p$make_report),
    # FIX: surface a missing decoupleR package before the pipeline starts
    # (see mod_decouple.R / utils_core.R).
    run_decouple   = isTRUE(p$run_dge) && isTRUE(p$run_decouple)
  )

  if (
    isTRUE(p$run_predictors) && .Platform$OS.type == "windows" && !p$use_wsl
  ) {
    warning(
      "run_predictors = TRUE but use_wsl = FALSE on Windows. ",
      "External tools (CPAT, SignalP, Pfam) normally run inside WSL on Windows; ",
      "set use_wsl = TRUE, or ensure bash + the tools are directly reachable ",
      "(e.g. Git-Bash/MSYS2) if you intend to run them natively."
    )
  }

  ensembl_package_name <- p$ensembl_package_name
  if (!requireNamespace(ensembl_package_name, quietly = TRUE)) {
    message(
      "Package '", ensembl_package_name,
      "' not found. Attempting automatic build and installation..."
    )
    parsed_ensembl <- .parse_ensembl_package_name(ensembl_package_name)
    archive_dir <- tempfile("expressom_ensdb_")
    if (!dir.create(archive_dir)) {
      stop("Could not create temporary Ensembl archive directory.", call. = FALSE)
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

  sample_meta_for_cond <- tryCatch(
    {
      meta <- data.table::fread(
        p$sample_table, header = TRUE, data.table = FALSE
      )
      .apply_sample_filters(
        meta,
        .sample_id_column(meta),
        p$remove_sample,
        p$subset_sample
      )
    },
    error = function(e) NULL
  )

  if (
    isTRUE(p$eda_only) || is.null(p$model) || is.null(p$level) || is.null(p$base)
  ) {
    message(
      "DESeq2 model not fully specified or eda_only = TRUE. ",
      "Running only data import and exploratory analysis (EDA)."
    )
    p$run_dge <- FALSE
    p$run_isoform <- FALSE
    model_eda <- "~1"

    sample_cols <- if (!is.null(sample_meta_for_cond)) {
      colnames(sample_meta_for_cond)
    } else {
      tryCatch(
        colnames(utils::read.csv(p$sample_table, nrows = 1, check.names = FALSE)),
        error = function(e) character(0)
      )
    }

    if (!is.null(p$group_col) && p$group_col %in% sample_cols) {
      main_condition <- p$group_col
    } else if (!is.null(p$model) && p$model != "~1") {
      main_condition <- .resolve_main_condition(
        stats::as.formula(p$model),
        sample_meta_for_cond,
        p$level,
        p$base
      )
    } else {
      main_condition <- NULL
    }
  } else {
    model_eda <- p$model
    main_condition <- .resolve_main_condition(
      stats::as.formula(p$model),
      sample_meta_for_cond,
      p$level,
      p$base
    )
  }

  comp_name <- if (!is.null(p$level) && !is.null(p$base)) {
    paste0(p$level, "vs", p$base)
  } else {
    "EDA"
  }

  if (!is.null(p$batch_col) && !is.null(p$model)) {
    model_vars <- all.vars(stats::as.formula(p$model))
    if (!(p$batch_col %in% model_vars)) {
      message(
        "WARNING: 'batch_col' [", p$batch_col,
        "] is specified for visualization, but missing from your design formula: '",
        p$model, "'."
      )
      message(
        "Consider adding it to prevent confounding your DGE results (e.g., model = '~ ",
        p$batch_col, " + ", main_condition, "')."
      )
    }
  }

  list(
    p = p,
    main_condition = main_condition,
    model_eda = model_eda,
    comp_name = comp_name
  )
}

# ---- stage 1: import + EDA ---------------------------------------------------------------------

.stage_import <- function(ctx) {
  p <- ctx$p
  import_counts(
    data_dir = p$data_dir,
    sample_table = p$sample_table,
    ensembl_package_name = p$ensembl_package_name,
    count_type = p$count_type,
    out_dir = p$out_dir,
    matrix_file = p$matrix_file,
    subset_sample = p$subset_sample,
    remove_sample = p$remove_sample,
    custom_tx2gene = p$custom_tx2gene,
    custom_gene_map = p$custom_gene_map
  )
}

.stage_eda <- function(ctx, tx_data) {
  p <- ctx$p
  edb_obj <- .stage_edb(ctx)

  dds <- create_dds_object(
    tx_data,
    level = if (p$run_dge) p$level else NULL,
    base = if (p$run_dge) p$base else NULL,
    model = ctx$model_eda,
    replicate_col = p$replicate_col
  )

  run_eda(
    dds = dds,
    edb = edb_obj,
    out_dir = p$out_dir,
    level = p$level,
    base = p$base,
    main_condition = ctx$main_condition,
    group_col = p$group_col,
    batch_col = p$batch_col,
    pca_ntop = p$pca_ntop,
    pca_colors = p$pca_colors,
    sample_labels = p$sample_labels
  )
  .close_devices()

  list(files = .out_manifest(p$out_dir, "Plots"))
}

# ---- stage 2: gene-level differential expression -----------------------------------------------

.stage_dge <- function(ctx, tx_data) {
  p <- ctx$p
  if (!isTRUE(p$run_dge)) return(NULL)
  message("\n=== Running Gene-Level DGE Analysis ===")

  dds <- create_dds_object(
    tx_data, p$level, p$base, p$model, p$replicate_col
  )

  res_list <- run_deseq2_analysis(
    dds, p$model, p$level, p$base, p$shrink_method, p$out_dir,
    p$padj_cutoff, p$test, p$reduced,
    lfc_cutoff = p$deg_lfc_cutoff
  )

  results_data <- export_significant_results(
    res_shrunken = res_list$res_shrunken,
    res_unshrunken = res_list$res_unshrunken,
    dds = res_list$dds,
    out_dir = p$out_dir,
    level = p$level,
    base = p$base,
    gene_map = tx_data$gene_map,
    padj_cutoff = p$padj_cutoff,
    padj_cutoffs = p$deg_padj_cutoffs,
    lfc_cutoff = p$deg_lfc_cutoff,
    test_type = p$test
  )

  list(
    dds = res_list$dds,
    res_list = res_list,
    results_data = results_data,
    files = .out_manifest(p$out_dir, "DE_raw_results")
  )
}

.stage_dge_plots <- function(ctx, tx_data, dge) {
  p <- ctx$p
  if (is.null(dge)) return(NULL)
  edb_obj <- .stage_edb(ctx)
  res_list <- dge$res_list
  plot_dir <- file.path(p$out_dir, "Plots")

  dds_rep <- res_list$dds
  clean_ens <- strip_ensembl_version(rownames(dds_rep))
  sym_map <- tx_data$gene_map$symbol[
    match(clean_ens, tx_data$gene_map$ensembl)
  ]
  sym_map[is.na(sym_map) | sym_map == ""] <- clean_ens[
    is.na(sym_map) | sym_map == ""
  ]
  sym_map <- make.unique(sym_map)
  rownames(dds_rep) <- sym_map

  generate_bulk_visualizations(
    dds = res_list$dds,
    edb = edb_obj,
    res_shrunken = res_list$res_shrunken,
    res_unshrunken = res_list$res_unshrunken,
    results_data = dge$results_data,
    out_dir = p$out_dir,
    level = p$level,
    base = p$base,
    main_condition = ctx$main_condition,
    top_genes = p$top_genes,
    padj_cutoff = p$padj_cutoff,
    highlight_genes = p$highlight_genes,
    batch_col = p$batch_col,
    pca_ntop = p$pca_ntop,
    pca_colors = p$pca_colors,
    sample_labels = p$sample_labels,
    volcano_colors = p$volcano_colors,
    lfc_cutoff = p$deg_lfc_cutoff
  )
  .close_devices()

  if (!is.null(p$zscore_genes)) {
    message("Generating targeted expression heatmaps (Sample Z-score & L2FC)...")
    plot_sample_zscore(
      dds_rep, p$zscore_genes, ctx$main_condition, p$level, p$base, plot_dir
    )
    plot_l2fc_heatmap(
      dds_rep, p$zscore_genes, ctx$main_condition, p$level, p$base, plot_dir
    )
  }

  if (!is.null(p$gene_sets_zscore)) {
    message("Generating separate Z-score plots for each gene set...")
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

    for (set_name in names(p$gene_sets_zscore)) {
      single_set <- list(p$gene_sets_zscore[[set_name]])
      names(single_set) <- set_name
      message("   -> Plotting set: ", set_name)

      plot_geneset_zscore_avg(
        dds = dds_rep, gene_sets = single_set,
        condition_col = ctx$main_condition,
        level = p$level, base = p$base,
        plot_dir = plot_dir, set_name = set_name
      )
      plot_gene_zscore_individual(
        dds = dds_rep, gene_sets = single_set,
        condition_col = ctx$main_condition,
        level = p$level, base = p$base,
        plot_dir = plot_dir, set_name = set_name
      )
    }
  }
  .close_devices()

  list(files = .out_manifest(p$out_dir, "Plots"))
}

# ---- stage 3: functional analysis --------------------------------------------------------------

.stage_functional <- function(ctx, dge) {
  p <- ctx$p
  if (is.null(dge)) return(NULL)
  edb_obj <- .stage_edb(ctx)

  func_results <- safe_run(
    run_functional_analysis(
      res_tbl = dge$results_data$res_tbl,
      sig_res = dge$results_data$sig_res,
      edb = edb_obj,
      out_dir = p$out_dir,
      level = p$level,
      base = p$base,
      top_genes = p$top_genes,
      padj_cutoff = p$padj_cutoff,
      go_pvalue_cutoff = p$go_pvalue_cutoff,
      go_qvalue_cutoff = p$go_qvalue_cutoff,
      gsea_metric = p$gsea_metric,
      test_type = p$test,
      ora_padj_cutoff = p$ora_padj_cutoff,
      ora_min_genes = p$ora_min_genes,
      ora_lfc_cutoff = p$ora_lfc_cutoff,
      run_spia = p$run_spia
    ),
    label = "Functional analysis"
  )
  .close_devices()

  list(
    func_results = func_results,
    files = .out_manifest(
      p$out_dir,
      c("ORA", "SPIA", "Transcription_Factors", file.path("GSEA", "GO"),
        file.path("GSEA", "Reactome"), file.path("GSEA", "KEGG"))
    )
  )
}

.stage_fgsea <- function(ctx, dge) {
  p <- ctx$p
  if (is.null(dge)) return(NULL)
  edb_obj <- .stage_edb(ctx)
  message("Running FGSEA Analysis...")

  safe_run(
    run_fgsea_analysis(
      res_tbl = dge$results_data$res_tbl,
      gmt_file = p$gmt_file,
      edb = edb_obj,
      out_dir = p$out_dir,
      comp_name = ctx$comp_name,
      padj_cutoff = p$padj_cutoff,
      gsea_metric = p$gsea_metric,
      test_type = p$test
    ),
    label = "FGSEA analysis"
  )
  .close_devices()

  list(files = .out_manifest(p$out_dir, file.path("GSEA", "FGSEA")))
}

# ---- stage 4: decoupleR TF/pathway activity ------------------------------------------------

.stage_decouple <- function(ctx, dge) {
  p <- ctx$p
  if (!isTRUE(p$run_dge) || !isTRUE(p$run_decouple) || is.null(dge)) return(NULL)
  message("\n=== Running decoupleR activity inference ===")

  res <- safe_run(
    run_decouple_analysis(
      res_tbl     = dge$results_data$res_tbl,
      out_dir     = p$out_dir,
      comp_name   = ctx$comp_name,
      organism    = get_organism_info(.stage_edb(ctx))$name,
      resources   = p$decouple_resources,
      method      = p$decouple_method,
      min_size    = p$decouple_min_size,
      padj_cutoff = p$padj_cutoff,
      # NEW: forwarded so the regulon volcano / leading edge categories
      # use the same |log2FC| floor as the rest of the pipeline.
      lfc_cutoff  = p$deg_lfc_cutoff
    ),
    label = "decoupleR activity inference"
  )
  .close_devices()

  list(
    decouple_results = res,
    files = .out_manifest(p$out_dir, "DecoupleR")
  )
}

# ---- stage 5: isoform analysis -------------------------------------------------------------------

.stage_isoform <- function(ctx) {
  p <- ctx$p
  if (!isTRUE(p$run_isoform)) return(NULL)
  message("\n=== Running Isoform-Level Analysis ===")

  out_dir <- p$out_dir
  data_dir <- p$data_dir
  ensembl_package_name <- p$ensembl_package_name
  count_type <- p$count_type
  matrix_file <- p$matrix_file
  custom_tx2gene <- p$custom_tx2gene
  custom_gene_map <- p$custom_gene_map
  custom_transcript_id_map <- p$custom_transcript_id_map
  subset_sample <- p$subset_sample
  remove_sample <- p$remove_sample
  sample_table <- p$sample_table
  isoform_fasta <- p$isoform_fasta
  isoform_gff <- p$isoform_gff
  model <- p$model
  level <- p$level
  base <- p$base
  test <- p$test
  reduced <- p$reduced
  padj_cutoff <- p$padj_cutoff
  run_dexseq <- p$run_dexseq
  run_isoform_enrichment <- p$run_isoform_enrichment
  ora_padj_cutoff <- p$ora_padj_cutoff
  ora_min_genes <- p$ora_min_genes
  ora_lfc_cutoff <- p$ora_lfc_cutoff
  gsea_metric <- p$gsea_metric
  gmt_file <- p$gmt_file
  top_genes <- p$top_genes
  go_pvalue_cutoff <- p$go_pvalue_cutoff
  go_qvalue_cutoff <- p$go_qvalue_cutoff
  run_predictors <- p$run_predictors
  predictor_cpu <- p$predictor_cpu
  use_wsl <- p$use_wsl
  wsl_distro <- p$wsl_distro
  run_spia <- p$run_spia
  isoform_test_engine <- p$isoform_test_engine
  resume_isoform_from <- p$resume_isoform_from
  main_condition <- ctx$main_condition
  batch_col <- p$batch_col
  pca_ntop <- p$pca_ntop
  pca_colors <- p$pca_colors
  sample_labels <- p$sample_labels
  skip_fasta_filter <- p$skip_fasta_filter
  plot_topology <- p$plot_topology
  isoform_report_genes <- p$isoform_report_genes
  isoform_plot_top_n <- p$isoform_plot_top_n
  bpparam <- p$bpparam %||% BiocParallel::bpparam()
  edb_obj <- .stage_edb(ctx)

  iso_dir <- file.path(out_dir, "IsoformSwitch")
  iso_save_dir <- file.path(iso_dir, "saved")
  dir.create(iso_save_dir, recursive = TRUE, showWarnings = FALSE)

  .iso_ckpt_exists <- function(nm) file.exists(file.path(iso_save_dir, nm))
  .iso_ckpt_save <- function(obj, nm) {
    .checkpoint_save(obj, file.path(iso_save_dir, nm), iso_checkpoint_signature)
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

  if (is.null(isoform_fasta) || is.null(isoform_gff)) {
    message("   -> Auto-fetching Ensembl FASTA and GTF references...")
    ref_dir <- safe_dir(file.path(data_dir, "reference"))
    refs <- download_ensembl_refs(
      ensembl_package_name = ensembl_package_name, out_dir = ref_dir
    )
    isoform_refs <- .resolve_isoform_references(isoform_fasta, isoform_gff, refs)
    isoform_fasta <- isoform_refs$fasta
    isoform_gff <- isoform_refs$gff
  }

  quant_paths <- .quant_input_files(p)
  checkpoint_input_files <- c(
    sample_table, matrix_file, custom_tx2gene, custom_gene_map,
    if (!is.data.frame(custom_transcript_id_map)) custom_transcript_id_map,
    quant_paths, isoform_fasta, isoform_gff
  )
  iso_checkpoint_signature <- .object_md5(list(
    inputs = .file_fingerprint(checkpoint_input_files),
    custom_transcript_id_map = if (is.data.frame(custom_transcript_id_map)) {
      custom_transcript_id_map
    } else NULL,
    count_type = count_type, subset_sample = subset_sample,
    remove_sample = remove_sample, model = model, level = level, base = base,
    padj_cutoff = padj_cutoff, run_dexseq = run_dexseq,
    run_isoform_enrichment = run_isoform_enrichment,
    ora_padj_cutoff = ora_padj_cutoff, ora_min_genes = ora_min_genes,
    ora_lfc_cutoff = ora_lfc_cutoff, gsea_metric = gsea_metric,
    gmt_file = gmt_file, top_genes = top_genes,
    go_pvalue_cutoff = go_pvalue_cutoff, go_qvalue_cutoff = go_qvalue_cutoff,
    run_predictors = run_predictors,
    predictor_env = if (isTRUE(run_predictors)) {
      .predictor_environment_fingerprint(
        wsl_distro = wsl_distro, use_wsl = use_wsl
      )
    } else NULL,
    predictor_cpu = predictor_cpu, use_wsl = use_wsl,
    wsl_distro = wsl_distro, run_spia = run_spia,
    isoform_test_engine = isoform_test_engine
  ))

  if (!is.null(resume_isoform_from) && dir.exists(resume_isoform_from)) {
    message("Resuming isoform analysis from: ", resume_isoform_from)
    loaded <- load_isoform_results(
      resume_isoform_from,
      signature = iso_checkpoint_signature,
      require_signature = TRUE
    )
    isoform_import <- loaded$isoform_import
    dte_res <- loaded$dte_results
    dtu_res <- loaded$dtu_results
    dexseq_res <- loaded$dexseq_results
    if (!is.null(dte_res) && !is.null(dtu_res)) {
      message("Loaded compatible DTE and DTU results from resume directory.")
    } else {
      message("No compatible complete resume results found. Missing steps will run.")
    }
  }

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
      "Please install IsoformSwitchAnalyzeR: ",
      "BiocManager::install('IsoformSwitchAnalyzeR')"
    )
  }

  if (is.null(isoform_import)) {
    message("Step A: Importing transcript-level counts...")
    isoform_import <- import_transcript_counts(
      data_dir = data_dir, sample_table = sample_table,
      ensembl_package_name = ensembl_package_name,
      count_type = count_type, matrix_file = matrix_file,
      subset_sample = subset_sample, remove_sample = remove_sample,
      custom_tx2gene = custom_tx2gene, custom_gene_map = custom_gene_map
    )
    .iso_ckpt_save(isoform_import, "isoform_import.rds")
  } else {
    message("Step A: Isoform import loaded from checkpoint. Skipping.")
  }

  safe_run(
    run_isoform_pca(
      isoform_import, main_condition, level, base,
      out_dir = iso_dir, batch_col = batch_col, pca_ntop = pca_ntop,
      save_dir = iso_save_dir, pca_colors = pca_colors,
      sample_labels = sample_labels
    ),
    label = "Transcript-level PCA"
  )

  if (is.null(dte_res)) {
    message("Step B: Running DTE (Differential Transcript Expression)...")
    dte_res <- run_dte(
      isoform_import, main_condition, level, base, padj_cutoff,
      bpparam = bpparam, design = model
    )
    .iso_ckpt_save(dte_res, "dte_results.rds")
  } else {
    message("Step B: DTE results loaded from checkpoint. Skipping.")
  }

  if (is.null(dtu_res)) {
    message("Step C: Running DTU (Differential Transcript Usage)...")
    dtu_res <- run_dtu(
      isoform_import, main_condition, level, base,
      bpparam = bpparam, design = model
    )
    .iso_ckpt_save(dtu_res, "dtu_results.rds")
  } else {
    message("Step C: DTU results loaded from checkpoint. Skipping.")
  }

  if (isTRUE(run_dexseq)) {
    if (is.null(dexseq_res)) {
      message("Step C.5: Running DEXSeq-based DTU (complementary engine)...")
      dexseq_res <- safe_run(
        run_dexseq_dtu(
          isoform_import, main_condition, level, base,
          bpparam = bpparam, design = model
        ),
        label = "DEXSeq DTU"
      )
      if (!is.null(dexseq_res)) .iso_ckpt_save(dexseq_res, "dexseq_results.rds")
    } else {
      message("Step C.5: DEXSeq DTU results loaded from checkpoint. Skipping.")
    }
  }

  dte_func_res <- NULL
  dtu_func_res <- NULL

  if (isTRUE(run_isoform_enrichment)) {
    dte_func_checkpoint <- if (.iso_ckpt_exists("dte_func_results.rds")) {
      .iso_ckpt_load("dte_func_results.rds")
    }
    if (!is.null(dte_func_checkpoint)) {
      dte_func_res <- dte_func_checkpoint
    } else if (!is.null(dte_res)) {
      message("Step D: Running functional enrichment on DTE results...")
      dte_func_res <- run_isoform_functional_analysis(
        res_tbl = dte_res, edb = edb_obj, out_dir = iso_dir,
        analysis_type = "DTE", level = level, base = base,
        top_genes = top_genes, padj_cutoff = padj_cutoff,
        go_pvalue_cutoff = go_pvalue_cutoff, go_qvalue_cutoff = go_qvalue_cutoff,
        gsea_metric = gsea_metric, ora_padj_cutoff = ora_padj_cutoff,
        ora_min_genes = ora_min_genes, ora_lfc_cutoff = ora_lfc_cutoff,
        gmt_file = gmt_file, run_spia = run_spia
      )
      .iso_ckpt_save(dte_func_res, "dte_func_results.rds")
    } else {
      message("Step D: No DTE results available. Skipping DTE functional enrichment.")
    }

    dtu_func_checkpoint <- if (.iso_ckpt_exists("dtu_func_results.rds")) {
      .iso_ckpt_load("dtu_func_results.rds")
    }
    if (!is.null(dtu_func_checkpoint)) {
      dtu_func_res <- dtu_func_checkpoint
    } else if (!is.null(dtu_res) && !is.null(dtu_res$dtu_results)) {
      message("Step D: Running functional enrichment on DTU results...")
      dtu_func_res <- run_isoform_functional_analysis(
        res_tbl = dtu_res$dtu_results, edb = edb_obj, out_dir = iso_dir,
        analysis_type = "DTU", level = level, base = base,
        top_genes = top_genes, padj_cutoff = padj_cutoff,
        go_pvalue_cutoff = go_pvalue_cutoff, go_qvalue_cutoff = go_qvalue_cutoff,
        ora_padj_cutoff = ora_padj_cutoff, ora_min_genes = ora_min_genes,
        ora_lfc_cutoff = ora_lfc_cutoff
      )
      .iso_ckpt_save(dtu_func_res, "dtu_func_results.rds")
    } else {
      message("Step D: No DTU results available. Skipping DTU functional enrichment.")
    }
  } else {
    message("Step D: Skipping DTE/DTU functional enrichment (run_isoform_enrichment = FALSE).")
  }

  .warn_isoform_switch_covariates(model, main_condition)

  switch_res <- run_isoform_switch(
    dte_results = dte_res, dtu_results = dtu_res,
    isoform_obj = isoform_import,
    condition = main_condition, level = level, base = base,
    fasta_file = isoform_fasta, gff_file = isoform_gff,
    out_dir = iso_dir, run_predictors = run_predictors,
    use_wsl = use_wsl, wsl_distro = wsl_distro,
    save_dir = iso_save_dir, resume_from = resume_isoform_from,
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
      dte_results = dte_res, dtu_results = dtu_res,
      isoform_obj = isoform_import, out_dir = out_dir,
      condition = main_condition, level = level, base = base,
      genes_of_interest = isoform_report_genes, top_n = 15,
      switch_list = switch_res, dexseq_results = dexseq_res,
      switch_plot_top_n = isoform_plot_top_n,
      plot_topology = plot_topology, padj_cutoff = padj_cutoff,
      design = model, dge_test = test, reduced = reduced,
      render_html = !isTRUE(p$make_report)
    )
  }

  if (!dir.exists(iso_dir)) dir.create(iso_dir, recursive = TRUE)
  message("Isoform analysis complete. Results saved in: ", iso_dir)
  .close_devices()

  list(
    isoform_import = isoform_import,
    dte_res = dte_res, dtu_res = dtu_res, dexseq_res = dexseq_res,
    dte_func_res = dte_func_res, dtu_func_res = dtu_func_res,
    switch_res = switch_res,
    files = .out_manifest(out_dir, "IsoformSwitch")
  )
}

# ---- stage 6: unified report + final .RData -----------------------------------------------------

.as_report_settings <- function(x) {
  base_settings <- expressom_report_settings()
  if (is.null(x)) return(base_settings)
  utils::modifyList(base_settings, as.list(x))
}

.stage_report <- function(ctx, ...) {
  p <- ctx$p
  comparison <- if (isTRUE(p$run_dge)) {
    paste0(p$level, "_vs_", p$base)
  } else {
    NULL
  }
  if (!is.null(comparison) &&
      !comparison %in% list_expressom_comparisons(p$out_dir)) {
    comparison <- NULL
  }
  build_expressom_report(
    out_dir = p$out_dir,
    comparison = comparison,
    settings = .as_report_settings(p$report_settings)
  )
}

.stage_save_rdata <- function(ctx, tx_data, dge, func, isoform, decouple = NULL) {
  p <- ctx$p
  rdata_dir <- safe_dir(file.path(p$out_dir, "Save_rdata"))
  rdata_path <- file.path(rdata_dir, paste0("Results_", ctx$comp_name, ".RData"))

  e <- new.env(parent = emptyenv())
  e$comp_name <- ctx$comp_name
  e$edb_obj <- .stage_edb(ctx)
  e$main_condition <- ctx$main_condition
  e$level <- p$level
  e$base <- p$base
  e$model <- p$model
  e$test <- p$test
  e$batch_col <- p$batch_col
  e$tx_data <- tx_data

  if (isTRUE(p$run_dge) && !is.null(dge)) {
    e$dds <- dge$dds
    e$res_list <- dge$res_list
    e$results_data <- dge$results_data
    e$func_results <- func$func_results
    e$decouple_results <- if (!is.null(decouple)) decouple$decouple_results else NULL
    e$dge_results <- list(
      tx_data = tx_data,
      dds = dge$dds,
      res_list = dge$res_list,
      results_data = dge$results_data,
      func_results = func$func_results,
      decouple_results = if (!is.null(decouple)) decouple$decouple_results else NULL
    )
  }

  if (isTRUE(p$run_isoform) && !is.null(isoform)) {
    e$isoform_import <- isoform$isoform_import
    e$dte_res <- isoform$dte_res
    e$dtu_res <- isoform$dtu_res
    e$switch_res <- isoform$switch_res
    e$isoform_results <- isoform[setdiff(names(isoform), "files")]
  }

  save(list = ls(e, all.names = TRUE), file = rdata_path, envir = e)
  message("\nR environment saved to: ", rdata_path)
  invisible(rdata_path)
}