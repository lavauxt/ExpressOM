# mod_isoform_switch.R -- IsoformSwitchAnalyzeR orchestration (run_isoform_switch)
#

#' Run IsoformSwitchAnalyzeR analysis
#' @export
run_isoform_switch <- function(dte_results = NULL,
                               dtu_results = NULL,
                               isoform_obj,
                               condition,
                               level,
                               base,
                               fasta_file,
                               gff_file,
                               out_dir,
                               run_predictors = FALSE,
                               use_wsl = (.Platform$OS.type == "windows"),
                               wsl_distro = "Ubuntu-22.04",
                               save_dir = NULL,
                               resume_from = NULL,
                               bsgenome_name = NULL,
                               predictor_cpu = NULL,
                               log_dir = NULL,
                               custom_transcript_id_map = NULL,
                               skip_fasta_filter = FALSE,
                               test_engine = c("DEXSeq", "DRIMSeq", "satuRn"),
                               organism = NULL,
                               plot_topology = TRUE) {

  test_engine <- match.arg(test_engine)

  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    stop("Please install IsoformSwitchAnalyzeR: BiocManager::install('IsoformSwitchAnalyzeR')")
  }

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  if (is.null(log_dir)) log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(run_predictors)) {
    message("Checking external predictor tool availability...")

    debug_info <- debug_wsl(
      distro = wsl_distro,
      out_dir = out_dir,
      log_dir = log_dir,
      conda_env = "isoform_tools",
      verbose = TRUE,
      use_wsl = use_wsl
    )

    if (!isTRUE(debug_info$wsl_available)) {
      stop(
        if (.Platform$OS.type == "windows" && use_wsl) {
          "WSL is not available. Please install WSL or set use_wsl = FALSE."
        } else {
          paste0(
            "Could not execute external predictor tools in this environment ",
            "(no working bash shell found)."
          )
        }
      )
    }

    if (!isTRUE(debug_info$conda_env_exists)) {
      warning(
        "Conda environment 'isoform_tools' not found. Predictors may still work if the ",
        "required tools (CPAT, SignalP, hmmscan/InterProScan) are already on PATH.\n",
        "On Windows/WSL run install_wsl_isoform_tools() to set up the environment; ",
        "on Linux/macOS install these tools manually or via conda/mamba."
      )
    }

    if (isTRUE(debug_info$pfam_db_found) && !isTRUE(debug_info$pfam_db_indexed)) {
      warning(
        "Pfam-A.hmm was found but is not hmmpress-indexed -- hmmscan will fail against it. ",
        "Run install_isoform_databases() again or `hmmpress` it manually."
      )
    }

    .log_wsl_command(
      "debug_wsl()",
      exit_code = 0,
      stdout = capture.output(print(debug_info)),
      log_dir = log_dir
    )
  } else {
    message(
      "run_predictors = FALSE: skipping debug_wsl() environment check and all ",
      "external predictor (CPAT / SignalP / Pfam) steps for the isoform analysis. ",
      "No WSL commands will be run, so no wsl_commands.log or wsl_debug.json will ",
      "be produced. Pass run_predictors = TRUE to run_isoform_switch() to enable this."
    )
  }

  .ckpt_path <- function(nm) if (!is.null(save_dir)) file.path(save_dir, nm) else ""
  .ckpt_exists <- function(nm) {
    p <- .ckpt_path(nm)
    nzchar(p) && file.exists(p)
  }

  .ckpt_save <- function(obj, nm) {
    p <- .ckpt_path(nm)
    if (nzchar(p)) {
      saveRDS(obj, p)
      message("Checkpoint saved -> ", nm)
    }
    invisible(obj)
  }

  .ckpt_load <- function(nm) {
    p <- .ckpt_path(nm)
    message("Resuming from checkpoint: ", nm)
    readRDS(p)
  }

  switch_list <- NULL
  already_analyzed <- FALSE

  if (!is.null(resume_from) && file.exists(file.path(resume_from, "switch_list.rds"))) {
    message("Resuming from saved SwitchList: ", resume_from)
    switch_list <- readRDS(file.path(resume_from, "switch_list.rds"))

    if (!is.null(switch_list$isoformSwitchAnalysis)) {
      message("Full analysis already present - skipping combined analysis.")
      already_analyzed <- TRUE
    }
  } else if (.ckpt_exists("step2_analyzed.rds")) {
    switch_list <- .ckpt_load("step2_analyzed.rds")
    already_analyzed <- TRUE
  } else if (.ckpt_exists("step1_imported.rds")) {
    switch_list <- .ckpt_load("step1_imported.rds")
  } else {
    message("Building SwitchList from raw data...")

    if (isoform_obj$type == "tximport") {
      count_matrix <- isoform_obj$txi$counts
    } else {
      count_matrix <- isoform_obj$counts
    }

    sample_col <- if ("Sample" %in% colnames(isoform_obj$meta)) "Sample" else "sample_id"

    meta_sub <- isoform_obj$meta[isoform_obj$meta[[condition]] %in% c(base, level), , drop = FALSE]

    if (nrow(meta_sub) == 0) {
      stop("No samples found for the comparison groups: ", base, " and ", level)
    }

    keep_samples <- intersect(rownames(meta_sub), colnames(count_matrix))

    if (length(keep_samples) == 0) {
      stop("No matching sample IDs between metadata and count matrix after filtering.")
    }

    meta_sub <- meta_sub[keep_samples, , drop = FALSE]
    count_matrix <- count_matrix[, keep_samples, drop = FALSE]

    sample_vector <- meta_sub[[sample_col]]
    colnames(count_matrix) <- sample_vector

    group_counts <- table(meta_sub[[condition]])

    if (any(group_counts < 2)) {
      warning("One or both groups have fewer than 2 samples. DTU/switch analysis may be unreliable.")
    }

    if (!is.null(custom_transcript_id_map)) {
      if (is.character(custom_transcript_id_map) && file.exists(custom_transcript_id_map)) {
        map_df <- data.table::fread(custom_transcript_id_map, header = TRUE, data.table = FALSE)
      } else if (is.data.frame(custom_transcript_id_map)) {
        map_df <- custom_transcript_id_map
      } else {
        stop("custom_transcript_id_map must be a data.frame or a path to a TSV/CSV file.")
      }

      required_cols <- c("count_id", "fasta_id")

      if (!all(required_cols %in% colnames(map_df))) {
        stop("custom_transcript_id_map must contain columns: ", paste(required_cols, collapse = ", "))
      }

      old_rownames <- rownames(count_matrix)
      new_rownames <- map_df$fasta_id[match(old_rownames, map_df$count_id)]

      if (any(is.na(new_rownames))) {
        warning("Some count IDs were not found in the mapping file; they will be kept unchanged.")
        new_rownames[is.na(new_rownames)] <- old_rownames[is.na(new_rownames)]
      }

      rownames(count_matrix) <- new_rownames
      message("  Remapped ", sum(!is.na(new_rownames)), " transcript IDs using custom map.")
    }

    clean_id <- clean_transcript_id
    rownames(count_matrix) <- clean_id(rownames(count_matrix))

    if (!skip_fasta_filter) {
      message("Pre-filtering transcripts to match FASTA file...")

      fasta_seqs <- Biostrings::fasta.seqlengths(fasta_file)
      clean_fasta_ids <- clean_id(names(fasta_seqs))
      clean_rownames <- clean_id(rownames(count_matrix))

      keep_in_fasta <- clean_rownames %in% clean_fasta_ids
      count_matrix <- count_matrix[keep_in_fasta, , drop = FALSE]

      message(
        "  Kept ", nrow(count_matrix), " / ", length(clean_rownames),
        " transcripts matching the FASTA file."
      )

      if (nrow(count_matrix) == 0) {
        stop(
          "No transcript IDs match between count matrix and FASTA file. ",
          "Check ID formats (pipe-delimited GENCODE headers, version suffixes, ",
          "or a genuinely different annotation/quantification source). ",
          "You may supply a custom_transcript_id_map or set skip_fasta_filter = TRUE."
        )
      }
    } else {
      message(
        "Skipping FASTA pre-filtering (skip_fasta_filter = TRUE). ",
        "Will rely on importRdata()'s ignoreAfter* arguments."
      )
    }

    tx2gene <- isoform_obj$tx2gene
    tx2gene$tx_clean <- clean_id(tx2gene$tx_id)

    keep_tx <- rownames(count_matrix) %in% tx2gene$tx_clean
    count_matrix <- count_matrix[keep_tx, , drop = FALSE]
    tx2gene <- tx2gene[match(rownames(count_matrix), tx2gene$tx_clean), ]

    gene_counts <- table(tx2gene$gene_id)
    genes_multi <- names(gene_counts[gene_counts > 1])
    keep_genes <- tx2gene$gene_id %in% genes_multi

    count_matrix <- count_matrix[keep_genes, , drop = FALSE]
    tx2gene <- tx2gene[keep_genes, ]

    message(
      "  Kept ", nrow(count_matrix), " transcripts from ",
      length(unique(tx2gene$gene_id)), " genes after removing single-transcript genes."
    )

    isoform_count_matrix <- round(count_matrix)

    design_matrix <- data.frame(
      sampleID = sample_vector,
      condition = factor(meta_sub[[condition]], levels = c(base, level)),
      stringsAsFactors = FALSE
    )

    rownames(design_matrix) <- design_matrix$sampleID

    if (!"package:dplyr" %in% search()) attachNamespace("dplyr")

    .count_fasta_entries <- function(path) {
      tryCatch(length(Biostrings::fasta.seqlengths(path)), error = function(e) NA_integer_)
    }

    .count_gtf_transcripts <- function(path) {
      tryCatch(
        {
          con <- file(path, open = "rt")
          on.exit(close(con), add = TRUE)

          n <- 0L

          repeat {
            lines <- readLines(con, n = 200000L, warn = FALSE)
            if (length(lines) == 0L) break
            n <- n + sum(grepl("\ttranscript\t", lines, fixed = TRUE))
          }

          if (n == 0L) NA_integer_ else n
        },
        error = function(e) NA_integer_
      )
    }

    n_fasta <- .count_fasta_entries(fasta_file[1])
    n_gtf <- .count_gtf_transcripts(gff_file)
    n_count <- nrow(isoform_count_matrix)

    if (!is.na(n_fasta) && !is.na(n_gtf) && n_fasta > 0 && n_gtf > 0) {
      ratio <- min(n_fasta, n_gtf) / max(n_fasta, n_gtf)

      if (ratio < 0.5) {
        warning(
          "isoform_fasta and isoform_gff look like they describe very different ",
          "transcript sets (", n_fasta, " fasta sequences vs. ", n_gtf, " GTF transcripts; ",
          "count matrix after filtering has ", n_count, " transcripts). importRdata() will ",
          "likely fail its Jaccard-similarity check below. Verify that both come from the ",
          "same annotation/quantification source.",
          immediate. = TRUE
        )
      }
    }

    clean_ref_dir <- if (!is.null(save_dir)) {
      file.path(save_dir, "cleaned_reference")
    } else {
      file.path(out_dir, "cleaned_reference")
    }
    if (!dir.exists(clean_ref_dir)) dir.create(clean_ref_dir, recursive = TRUE)

    fasta_file_clean <- file.path(clean_ref_dir, "isoform_fasta_ids_cleaned.fasta")
    gff_file_clean <- file.path(clean_ref_dir, "isoform_gff_ids_cleaned.gtf")

    message(
      "  Pre-cleaning transcript IDs in FASTA/GTF so importRdata() gets exact matches..."
    )

    seqs <- Biostrings::readDNAStringSet(fasta_file)
    names(seqs) <- clean_id(names(seqs))

    dup_fa <- unique(names(seqs)[duplicated(names(seqs))])

    if (length(dup_fa) > 0) {
      stop(
        "Cleaning transcript IDs in isoform_fasta produced ", length(dup_fa),
        " duplicate ID(s) -- e.g. ", paste(utils::head(dup_fa, 5), collapse = ", "),
        ". The FASTA contains near-duplicate headers that only differ in the part ",
        "clean_id() strips."
      )
    }

    Biostrings::writeXStringSet(seqs, filepath = fasta_file_clean)
    rm(seqs)
    gc()

    .rewrite_gtf_transcript_ids <- function(path_in, path_out, clean_fn) {
      con_in <- file(path_in, open = "rt")
      on.exit(close(con_in), add = TRUE)

      con_out <- file(path_out, open = "wt")
      on.exit(close(con_out), add = TRUE)

      id_re <- 'transcript_id "([^"]+)"'
      version_re <- '\\s*(gene|transcript|exon)_version "[^"]*";'

      any_id_seen <- FALSE
      tx_ids_seen <- character(0)

      repeat {
        lines <- readLines(con_in, n = 200000L, warn = FALSE)
        if (length(lines) == 0L) break

        lines <- gsub(version_re, "", lines, perl = TRUE)

        m <- regexpr(id_re, lines)
        has_id <- m != -1L

        if (any(has_id)) {
          any_id_seen <- TRUE

          raw_ids <- sub(id_re, "\\1", regmatches(lines, m))
          new_ids <- clean_fn(raw_ids)

          regmatches(lines, m) <- paste0('transcript_id "', new_ids, '"')

          tx_row <- has_id & grepl("\ttranscript\t", lines, fixed = TRUE)

          if (any(tx_row)) {
            m2 <- regexpr(id_re, lines[tx_row])
            tx_ids_seen <- c(tx_ids_seen, sub(id_re, "\\1", regmatches(lines[tx_row], m2)))
          }
        }

        writeLines(lines, con_out)
      }

      if (!any_id_seen) {
        stop(
          "No `transcript_id \"...\"` attributes found in ", path_in,
          " -- expected standard GTF2-style attributes."
        )
      }

      dup <- unique(tx_ids_seen[duplicated(tx_ids_seen)])

      if (length(dup) > 0) {
        stop(
          "Cleaning transcript IDs in isoform_gff produced ", length(dup),
          " duplicate transcript ID(s) among `transcript` feature rows -- e.g. ",
          paste(utils::head(dup, 5), collapse = ", "),
          ". The GTF contains near-duplicate transcript records that only differ in the ",
          "part clean_id() strips."
        )
      }

      invisible(path_out)
    }

    .rewrite_gtf_transcript_ids(gff_file, gff_file_clean, clean_id)

    importrdata_extra_args <- list(estimateDifferentialGeneRange = FALSE)

    if ("detectAndCorrectUnwantedEffects" %in% names(formals(IsoformSwitchAnalyzeR::importRdata))) {
      importrdata_extra_args$detectAndCorrectUnwantedEffects <- FALSE
    }

    switch_list <- .muffle_across_deprecation(
      do.call(
        IsoformSwitchAnalyzeR::importRdata,
        c(
          list(
            isoformCountMatrix = isoform_count_matrix,
            isoformRepExpression = isoform_count_matrix,
            designMatrix = design_matrix,
            isoformExonAnnoation = gff_file_clean,
            isoformNtFasta = fasta_file_clean,
            ignoreAfterBar = FALSE,
            ignoreAfterSpace = FALSE,
            ignoreAfterPeriod = FALSE,
            showProgress = TRUE
          ),
          importrdata_extra_args
        )
      )
    )

    message("Isoform data import completed.")
    .ckpt_save(switch_list, "step1_imported.rds")

    if (!is.null(isoform_obj$txi$counts)) {
      isoform_obj$txi$counts <- NULL
      gc()
    }

    if (!is.null(isoform_obj$counts)) {
      isoform_obj$counts <- NULL
      gc()
    }
  }

  if (!already_analyzed) {
    message("Running DTU test using engine: ", test_engine)

    if (!is.null(switch_list$ntSequence) && length(switch_list$ntSequence) > 0) {
      genome_object <- NULL
    } else if (!is.null(bsgenome_name) && requireNamespace(bsgenome_name, quietly = TRUE)) {
      genome_object <- getExportedValue(bsgenome_name, bsgenome_name)
      message("Using BSgenome: ", bsgenome_name)
    } else if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      genome_object <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
      message("Using default BSgenome.Hsapiens.UCSC.hg38")
    } else {
      message("No BSgenome available - running WITHOUT ORF prediction.")
      genome_object <- NULL
    }

    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("dplyr is required but not installed.")
    }

    if (!"package:dplyr" %in% search()) {
      attachNamespace("dplyr")
    }

    switch_list <- switch(
      test_engine,
      DEXSeq = {
        message("Using built-in DEXSeq via isoformSwitchAnalysisCombined...")

        tryCatch(
          {
            sl <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisCombined(
              switchAnalyzeRlist = switch_list,
              genomeObject = genome_object,
              pathToOutput = out_dir,
              n = 50,
              outputPlots = FALSE
            )

            message("Combined analysis (DEXSeq) completed.")
            sl
          },
          error = function(e) {
            no_switches <- grepl("no genes were considered switching", conditionMessage(e), ignore.case = TRUE)

            if (no_switches) {
              message("\n=== No isoform switches detected ===")
              message(
                "isoformSwitchAnalysisCombined() found zero genes meeting the switching ",
                "cutoffs for '", level, "' vs '", base, "'. Isoform switch analysis will be skipped, ",
                "but DTE and DTU results are unaffected."
              )
            } else {
              message("\n=== isoformSwitchAnalysisCombined() failed ===")
              message("Error: ", conditionMessage(e))
              message("Isoform switch analysis will be skipped; DTE and DTU results are unaffected.")
            }

            NULL
          }
        )
      },
      DRIMSeq = {
        message("Using pre-computed DRIMSeq DTU results (from run_dtu)...")

        if (is.null(dtu_results)) {
          stop(
            "test_engine = 'DRIMSeq' requires pre-computed dtu_results. ",
            "Run DTU first or pass dtu_results."
          )
        }

        drim_df <- dtu_results$dtu_results

        if (is.null(drim_df)) {
          stop("dtu_results does not contain a 'dtu_results' element.")
        }

        required_cols <- c("gene_id", "feature_id", "pvalue", "adj_pvalue")

        if (!all(required_cols %in% colnames(drim_df))) {
          stop("DRIMSeq results must contain columns: ", paste(required_cols, collapse = ", "))
        }

        tryCatch(
          {
            sl <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisCombined(
              switchAnalyzeRlist = switch_list,
              genomeObject = genome_object,
              pathToOutput = out_dir,
              n = 50,
              dtuResults = drim_df,
              outputPlots = FALSE
            )

            message("Combined analysis (DRIMSeq) completed.")
            sl
          },
          error = function(e) {
            message("Error in isoformSwitchAnalysisCombined with DRIMSeq results: ", e$message)
            NULL
          }
        )
      },
      satuRn = {
        message("Running satuRn DTU test...")

        if (!requireNamespace("satuRn", quietly = TRUE)) {
          stop("Package 'satuRn' is required for test_engine = 'satuRn'. Install with BiocManager::install('satuRn').")
        }

        sl_tested <- tryCatch(
          {
            IsoformSwitchAnalyzeR::isoformSwitchTestSatuRn(
              switchAnalyzeRlist = switch_list,
              alpha = 0.05,
              dIFcutoff = 0.1,
              reduceToSwitchingGenes = FALSE
            )
          },
          error = function(e) {
            message("isoformSwitchTestSatuRn failed: ", e$message)
            NULL
          }
        )

        if (is.null(sl_tested)) {
          message("satuRn test did not return a switchList. Skipping further analysis.")
          NULL
        } else {
          message("Running isoformSwitchAnalysisPart2...")

          tryCatch(
            {
              sl_final <- IsoformSwitchAnalyzeR::isoformSwitchAnalysisPart2(
                switchAnalyzeRlist = sl_tested,
                genomeObject = genome_object,
                pathToOutput = out_dir,
                n = 50,
                outputPlots = FALSE
              )

              message("satuRn analysis completed.")
              sl_final
            },
            error = function(e) {
              message("isoformSwitchAnalysisPart2 failed: ", e$message)
              sl_tested
            }
          )
        }
      }
    )

    if (is.null(switch_list)) {
      return(invisible(NULL))
    }

    switch_list <- tryCatch(
      {
        message("Annotating alternative splicing events...")

        sl <- .muffle_across_deprecation(
          IsoformSwitchAnalyzeR::analyzeAlternativeSplicing(
            switch_list,
            onlySwitchingGenes = TRUE,
            alpha = 0.05,
            dIFcutoff = 0.1,
            quiet = FALSE
          )
        )

        message("  Alternative splicing annotation completed.")
        sl
      },
      error = function(e) {
        message(
          "  Could not annotate alternative splicing events: ", e$message,
          "\n  Splicing-type summary/enrichment plots will be skipped; all other results are unaffected."
        )
        switch_list
      }
    )

    .ckpt_save(switch_list, "step2_analyzed.rds")
  }

  if (isTRUE(run_predictors)) {
    if (.ckpt_exists("step3_predictors.rds")) {
      message("Loading external predictor results from step-3 checkpoint.")
      switch_list <- .ckpt_load("step3_predictors.rds")
    } else {
      message("Running external isoform predictors...")

      switch_list <- .run_external_predictors(
        switch_list = switch_list,
        fasta_file = fasta_file,
        out_dir = out_dir,
        use_wsl = use_wsl,
        wsl_distro = wsl_distro,
        isoform_obj = isoform_obj,
        save_dir = save_dir,
        n_cpu = predictor_cpu,
        log_dir = log_dir,
        organism = organism
      )

      .ckpt_save(switch_list, "step3_predictors.rds")
    }
  } else {
    message(
      "Skipping external isoform predictors (CPAT / SignalP / Pfam) because ",
      "run_predictors = FALSE. No WSL commands will be run and no wsl_commands.log / ",
      "wsl_debug.json will be written for this step. Pass run_predictors = TRUE ",
      "(and use_wsl = TRUE if the tools live in WSL) to enable coding-potential, ",
      "signal-peptide, and domain predictions."
    )
  }

  if (isTRUE(run_predictors)) {
    if (.ckpt_exists("step3_5_refreshed.rds")) {
      message("Loading refreshed switch consequences from step-3.5 checkpoint.")
      switch_list <- .ckpt_load("step3_5_refreshed.rds")
    } else {
      message("Refreshing switch consequence analysis and plots with predictor annotations...")

      # Lives under DTU_DTE_report/plots (not out_dir/plots) so it's inside
      # the same tree generate_dte_dtu_report() already scans for its PDF
      # -> PNG conversion, and so dte_dtu_report.Rmd can list it (see the
      # switch_by_consequence_section chunk there) -- it used to sit in an
      # entirely separate out_dir/plots/switch_plots_with_predictors tree
      # that neither of those ever looked inside, making it invisible to
      # the HTML report and just a second, disconnected pile of PDFs next
      # to top_switches/. "with predictors" was also a misleading name: by
      # the time either this step or the report's top_switches step runs,
      # switch_list already carries predictor annotations either way -- the
      # actual difference is that this set is filtered to switches with an
      # identified consequence and split by consequence type.
      plot_refresh_dir <- file.path(out_dir, "DTU_DTE_report", "plots", "by_consequence")

      switch_list <- tryCatch(
        {
          feat_cols <- colnames(switch_list$isoformFeatures)
          consequences_available <- c("intron_retention", "ORF_seq_similarity", "NMD_status")

          if (any(grepl("coding_potential|coding_prob", feat_cols, ignore.case = TRUE))) {
            consequences_available <- c(consequences_available, "coding_potential")
          }

          if (any(grepl("signal_peptide", feat_cols, ignore.case = TRUE))) {
            consequences_available <- c(consequences_available, "signal_peptide_identified")
          }

          if (any(grepl("^domain", feat_cols, ignore.case = TRUE))) {
            consequences_available <- c(consequences_available, "domains_identified", "domain_isotype")
          }

          if (!is.null(switch_list$topologyAnalysis)) {
            consequences_available <- c(consequences_available, "isoform_topology")
          }

          sl <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
            switch_list,
            consequencesToAnalyze = consequences_available,
            alpha = 0.05,
            dIFcutoff = 0.1
          )

          if (!dir.exists(plot_refresh_dir)) dir.create(plot_refresh_dir, recursive = TRUE)

          # switchPlotTopSwitches() doesn't expose plotTopology as a
          # parameter (verified against its signature and the
          # additionalArguments whitelist it forwards to switchPlot() --
          # plotTopology isn't in either), so suppressMessages() is the
          # only way to silence its "Omitting topology visualization..."
          # message from here when there's no DeepTMHMM data to plot.
          run_switch_plots <- function() {
            IsoformSwitchAnalyzeR::switchPlotTopSwitches(
              switchAnalyzeRlist = sl,
              n = 50,
              filterForConsequences = TRUE,
              splitComparison = FALSE,
              splitFunctionalConsequences = FALSE,
              sortByQvals = TRUE,
              pathToOutput = plot_refresh_dir,
              fileType = "pdf"
            )
          }

          if (isTRUE(plot_topology)) run_switch_plots() else suppressMessages(run_switch_plots())

          message(
            "  Refreshed switch plots (now including predictor annotations) saved to: ",
            plot_refresh_dir
          )

          sl
        },
        error = function(e) {
          message(
            "  Could not refresh switch consequences/plots with predictor data: ", e$message,
            "\n  Falling back to the plots generated in Step 2 (without predictor annotations)."
          )
          switch_list
        }
      )

      .ckpt_save(switch_list, "step3_5_refreshed.rds")
    }
  }

  if (!is.null(save_dir)) {
    saveRDS(switch_list, file.path(save_dir, "switch_list.rds"))
    message("Saved final SwitchList to ", save_dir)
  }

  message("Isoform switch analysis completed. Results saved in: ", out_dir)

  return(switch_list)
}
