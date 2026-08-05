# mod_isoform_predictors.R -- external predictor tools for isoform switch analysis
# (CPAT, SignalP, DeepTMHMM, Pfam/hmmscan/InterProScan), run via WSL


.run_external_predictors <- function(switch_list,
                                     fasta_file,
                                     out_dir,
                                     use_wsl,
                                     wsl_distro,
                                     isoform_obj,
                                     save_dir,
                                     n_cpu = NULL,
                                     log_dir = NULL,
                                     organism = NULL) {

  is_windows <- .Platform$OS.type == "windows"
  via_wsl <- is_windows && isTRUE(use_wsl)

  if (is.null(log_dir)) log_dir <- file.path(out_dir, "Log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  predictor_log_path <- file.path(log_dir, "predictor_status.log")

  message(
    "Predictor run logging to: ", predictor_log_path, " (step-by-step status) and ",
    file.path(log_dir, "wsl_commands.log"), " (every shell command, in real time)."
  )

  .log_predictor_status <- function(step, status, detail = "") {
    cat(
      sprintf(
        "[%s] %-28s %-8s %s\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        step,
        status,
        detail
      ),
      file = predictor_log_path,
      append = TRUE
    )
    invisible(NULL)
  }

  cat("\n\n############################################################\n", file = predictor_log_path, append = TRUE)
  cat(sprintf("# NEW PREDICTOR RUN: %s\n", format(Sys.time())), file = predictor_log_path, append = TRUE)
  cat(sprintf("# Output dir: %s\n", out_dir), file = predictor_log_path, append = TRUE)
  cat("############################################################\n\n", file = predictor_log_path, append = TRUE)

  message("Predictor status log: ", predictor_log_path)

  if (is.null(n_cpu) || !is.finite(n_cpu) || n_cpu < 1) {
    detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
    n_cpu <- if (is.na(detected) || detected < 1) 1L else max(1L, detected - 1L)
  }

  n_cpu <- as.integer(n_cpu)
  message("Using ", n_cpu, " CPU thread(s) for hmmscan / InterProScan.")

  message("Detecting conda environment in execution context...")

  cache_key <- paste0("conda_sh::", wsl_distro, "::", via_wsl)

  if (exists(cache_key, envir = .expressom_cache, inherits = FALSE)) {
    cached <- get(cache_key, envir = .expressom_cache, inherits = FALSE)
    conda_sh <- if (nzchar(cached)) cached else NULL
    has_conda_env <- nzchar(cached) > 0

    message(
      "  Using cached conda lookup: ",
      if (has_conda_env) conda_sh else "no 'isoform_tools' env found"
    )
  } else {
    # Same discovery strategy as debug_wsl(): ask conda itself first
    # (`conda info --base`), then fall back to a short list of common
    # install prefixes. A blind filesystem `find` for any file named
    # conda.sh can land on a non-functional copy cached inside conda's own
    # pkgs/ directory -- `conda env list` still succeeds after sourcing it
    # (so the old code accepted it), but `conda activate` from that copy is
    # not guaranteed to put the env's tools on PATH, which was silently
    # failing every tool_ok() check below despite debug_wsl() reporting the
    # same tools as present.
    conda_sh <- .find_conda_sh(wsl_distro = wsl_distro, use_wsl = via_wsl, log_dir = log_dir)
    has_conda_env <- FALSE

    if (!is.null(conda_sh)) {
      env_check <- .wsl_exec_script(
        bash_body = c(
          sprintf('. %s 2>/dev/null || true', .dq(conda_sh)),
          'conda env list 2>/dev/null | grep -q "isoform_tools"'
        ),
        wsl_distro = wsl_distro,
        use_wsl = via_wsl,
        conda_sh = NULL,
        intern = FALSE,
        ignore_stderr = TRUE,
        log_dir = log_dir
      )

      has_conda_env <- isTRUE(env_check == 0L)
    }

    if (!has_conda_env) conda_sh <- NULL

    assign(cache_key, if (has_conda_env) conda_sh else "", envir = .expressom_cache)

    if (has_conda_env) {
      message("  Using conda env 'isoform_tools' (", conda_sh, ")")
    } else {
      message("  No 'isoform_tools' env found. Using system PATH.")
    }
  }

  active_conda <- if (has_conda_env) conda_sh else NULL

  run_tool <- function(cmd_str, show_stderr = TRUE, tool_name = "tool") {
    tool_log <- file.path(log_dir, paste0(tool_name, ".log"))

    cat(sprintf("\n===== %s — %s =====\n", tool_name, format(Sys.time())), file = tool_log, append = TRUE)
    cat(sprintf("Working dir: %s\n", getwd()), file = tool_log, append = TRUE)
    cat(sprintf("Command: %s\n\n", cmd_str), file = tool_log, append = TRUE)

    message("  [", tool_name, "] Executing: ", cmd_str)

    res <- .wsl_exec_script(
      bash_body = cmd_str,
      wsl_distro = wsl_distro,
      use_wsl = via_wsl,
      conda_sh = active_conda,
      conda_env = "isoform_tools",
      intern = TRUE,
      ignore_stderr = !show_stderr,
      log_dir = log_dir,
      verbose = FALSE
    )

    exit_code <- attr(res, "status") %||% 0L

    if (length(res) > 0) {
      cat("--- stdout/stderr ---\n", file = tool_log, append = TRUE)
      cat(paste(res, collapse = "\n"), "\n", file = tool_log, append = TRUE)
    }

    cat(sprintf("--- exit code: %d ---\n", as.integer(exit_code)), file = tool_log, append = TRUE)

    if (length(res) > 0) {
      msg <- paste(res, collapse = "\n")

      if (nzchar(msg)) {
        msg_lines <- strsplit(msg, "\n")[[1]]

        if (length(msg_lines) > 30) {
          message(
            "  [", tool_name, "] Output (first 15 + last 15 of ",
            length(msg_lines), " lines; full log: ", tool_log, "):"
          )
          message(paste(head(msg_lines, 15), collapse = "\n"))
          message("  ... (", length(msg_lines) - 30, " lines omitted) ...")
          message(paste(tail(msg_lines, 15), collapse = "\n"))
        } else {
          message("  [", tool_name, "] Output (full log: ", tool_log, "):")
          message(msg)
        }
      }
    }

    message(
      "  [", tool_name, "] Exit code: ", as.integer(exit_code),
      " — ", if (exit_code == 0L) "SUCCESS" else "FAILED"
    )

    .log_predictor_status(
      tool_name,
      if (exit_code == 0L) "EXEC_OK" else "EXEC_FAILED",
      sprintf("exit code %d (see %s)", as.integer(exit_code), tool_log)
    )

    as.integer(exit_code)
  }

  tool_ok <- function(t) {
    .wsl_tool_exists(
      t,
      wsl_distro,
      via_wsl,
      active_conda,
      "isoform_tools",
      log_dir = log_dir
    )
  }

  w2l <- function(p) {
    if (is.null(p) || length(p) == 0 || !any(nzchar(p))) return(p)
    if (length(p) > 1) {
      return(vapply(p, w2l, character(1), USE.NAMES = FALSE))
    }
    if (!nzchar(p)) return(p)
    if (is_windows && use_wsl) .to_wsl_path(p, wsl_distro) else p
  }

  syms <- isoform_obj$gene_map$symbol[!is.na(isoform_obj$gene_map$symbol)]
  first_sym <- if (length(syms) > 0) syms[1] else ""

  organism_cpat <- if (!is.null(organism) && grepl("Homo sapiens", organism, ignore.case = TRUE)) {
    "Human"
  } else if (!is.null(organism) && grepl("Mus musculus", organism, ignore.case = TRUE)) {
    "Mouse"
  } else {
    if (!is.null(organism)) {
      message(
        "  Organism '", organism, "' is not Human or Mouse (CPAT only ships prebuilt logit ",
        "models for those two); falling back to a gene-symbol-casing guess for the CPAT model choice."
      )
    }
    if (grepl("^[A-Z]+$", first_sym) && nchar(first_sym) > 2) "Human" else "Mouse"
  }

  message("Organism detected for CPAT: ", organism_cpat, if (!is.null(organism)) paste0(" (from EnsDb: ", organism, ")") else " (guessed from gene symbol casing -- pass organism = ... for a reliable result)")

  seq_dir <- file.path(out_dir, "sequences")
  dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

  nt_fa_local <- file.path(seq_dir, "isoform_NT.fa")
  aa_fa_local <- file.path(seq_dir, "isoform_AA.fa")

  if (!file.exists(nt_fa_local) || !file.exists(aa_fa_local)) {
    message("Extracting NT / AA sequences from SwitchList for external tools...")

    tryCatch(
      {
        IsoformSwitchAnalyzeR::extractSequence(
          switch_list,
          onlySwitchingGenes = FALSE,
          extractNTseq = TRUE,
          extractAAseq = TRUE,
          writeToFile = TRUE,
          pathToOutput = seq_dir,
          addToSwitchAnalyzeRlist = FALSE,
          quiet = TRUE
        )

        # extractSequence() names its output using its own outputPrefix
        # (e.g. "isoformSwitchAnalyzeR_isoform_AA.fasta" /
        # "..._nt.fasta"), not a literal "NT.fasta"/"AA.fasta" -- that
        # exact-name assumption never matched, so nt_fa_local/aa_fa_local
        # were never created and CPAT/SignalP/Pfam all silently ran with
        # no (or the wrong) sequence input. Find what was actually written
        # by pattern instead of assuming an exact name, so this survives
        # IsoformSwitchAnalyzeR version differences too.
        written <- list.files(seq_dir, pattern = "\\.fa(sta)?$", full.names = TRUE, ignore.case = TRUE)
        written <- setdiff(written, c(nt_fa_local, aa_fa_local))
        written <- written[!grepl("complete|split|subset", basename(written), ignore.case = TRUE)]

        nt_src <- written[grepl("(^|[_.-])nt([_.-]|$)", basename(written), ignore.case = TRUE)]
        aa_src <- written[grepl("(^|[_.-])aa([_.-]|$)", basename(written), ignore.case = TRUE)]

        if (length(nt_src) >= 1 && !file.exists(nt_fa_local)) file.rename(nt_src[1], nt_fa_local)
        if (length(aa_src) >= 1 && !file.exists(aa_fa_local)) file.rename(aa_src[1], aa_fa_local)

        if (!file.exists(nt_fa_local) || !file.exists(aa_fa_local)) {
          message(
            "  extractSequence() ran but expected output was not found by pattern in ",
            seq_dir, " (files present: ",
            paste(basename(list.files(seq_dir)), collapse = ", "), ")"
          )
        }
      },
      error = function(e) {
        message("  Could not extract sequences from SwitchList: ", e$message)
        message("  CPAT will use the provided fasta_file; SignalP / Pfam may be skipped.")
        .log_predictor_status("extractSequence", "FAILED", e$message)
      }
    )
  }

  n_nt <- if (file.exists(nt_fa_local)) length(grep("^>", readLines(nt_fa_local, warn = FALSE))) else 0L
  n_aa <- if (file.exists(aa_fa_local)) length(grep("^>", readLines(aa_fa_local, warn = FALSE))) else 0L

  .log_predictor_status(
    "extractSequence",
    if (n_aa > 0) "OK" else "PARTIAL",
    sprintf(
      "%d NT sequences, %d AA sequences (AA=0 means SignalP/Pfam will be skipped below)",
      n_nt,
      n_aa
    )
  )

  cpat_fa_local <- if (file.exists(nt_fa_local)) {
    nt_fa_local
  } else {
    if (length(fasta_file) > 1) {
      message(
        "  fasta_file has ", length(fasta_file), " entries (e.g. combined cDNA + ncRNA ",
        "references); CPAT takes a single FASTA, so only the first (", fasta_file[1],
        ") will be used as a fallback. Sequence extraction from the SwitchList failed above, ",
        "which is why this fallback path was needed at all -- CPAT results may be incomplete."
      )
    }
    fasta_file[1]
  }
  sp_fa_local <- if (file.exists(aa_fa_local)) aa_fa_local else NULL
  pfam_fa_local <- if (file.exists(aa_fa_local)) aa_fa_local else NULL

  cpat_fa_w <- w2l(cpat_fa_local)
  sp_fa_w <- w2l(sp_fa_local)
  pfam_fa_w <- w2l(pfam_fa_local)
  out_dir_w <- w2l(out_dir)

  message("\n--- CPAT (coding potential) ---")

  hexamer_local <- .find_cpat_hexamer(
    organism_cpat,
    wsl_distro,
    via_wsl,
    active_conda,
    log_dir = log_dir
  )

  logit_local <- .find_cpat_logit_model(
    organism_cpat,
    wsl_distro,
    via_wsl,
    active_conda,
    log_dir = log_dir
  )

  if (is.null(hexamer_local)) {
    message(
      "  Hexamer table (", organism_cpat, "_Hexamer.tsv) not found in $CPAT_DATA / $HOME/.cpat_data. ",
      "Run install_isoform_databases() first. Skipping CPAT."
    )
    .log_predictor_status("CPAT", "SKIPPED", "hexamer table not found (run install_isoform_databases())")
  } else if (is.null(logit_local)) {
    message("  Logit model not found. Run install_isoform_databases() first. Skipping CPAT.")
    .log_predictor_status("CPAT", "SKIPPED", "logit model not found (run install_isoform_databases())")
  } else {
    hexamer_w <- hexamer_local
    logit_w <- logit_local
    cpat_prefix_w <- paste0(out_dir_w, "/cpat_out")
    cpat_prefix_l <- file.path(out_dir, "cpat_out")

    has_cpat3 <- tool_ok("cpat")
    has_cpat2 <- if (!has_cpat3) tool_ok("run_cpat.py") else FALSE

    if (!has_cpat3 && !has_cpat2) {
      message("  Neither 'cpat' nor 'run_cpat.py' found. Skipping.")
      .log_predictor_status("CPAT", "SKIPPED", "neither 'cpat' nor 'run_cpat.py' found on PATH/conda env")
    } else {
      cpat_exe <- if (has_cpat3) "cpat" else "run_cpat.py"

      cpat_cmd <- sprintf(
        "%s -g %s -x %s -d %s -o %s",
        cpat_exe,
        shQuote(cpat_fa_w, type = "sh"),
        shQuote(hexamer_w, type = "sh"),
        shQuote(logit_w, type = "sh"),
        shQuote(cpat_prefix_w, type = "sh")
      )

      message("  Running: ", cpat_exe)
      cpat_status <- run_tool(cpat_cmd, tool_name = "CPAT")

      cpat_result <- if (has_cpat3) {
        paste0(cpat_prefix_l, ".ORF_prob.best.tsv")
      } else {
        paste0(cpat_prefix_l, ".r")
      }

      if (cpat_status == 0L && file.exists(cpat_result)) {
        cpat_import_ok <- FALSE

        switch_list <- tryCatch(
          {
            sl <- IsoformSwitchAnalyzeR::analyzeCPAT(
              switchAnalyzeRlist = switch_list,
              pathToAllCPATresultFiles = cpat_result,
              codingCutoff = NULL,
              removeNoncodinORFs = TRUE,
              quiet = FALSE
            )

            cpat_import_ok <- TRUE
            sl
          },
          error = function(e) {
            message("  analyzeCPAT(): ", e$message)
            switch_list
          }
        )

        if (cpat_import_ok) {
          message("  CPAT results imported successfully.")
          .log_predictor_status("CPAT", "OK", paste("result file:", cpat_result))
        } else {
          message("  CPAT import failed; coding-potential annotation NOT added to switch_list.")
          .log_predictor_status("CPAT", "IMPORT_FAILED", "analyzeCPAT() errored")
        }
      } else {
        message("  CPAT failed or result file not found (", basename(cpat_result), "). Skipping.")
        .log_predictor_status(
          "CPAT",
          "FAILED",
          sprintf("exit code %s, result file exists: %s", cpat_status, file.exists(cpat_result))
        )
      }
    }
  }

  message("\n--- SignalP (signal peptide prediction) ---")

  if (is.null(sp_fa_w)) {
    message("  No AA FASTA available. Skipping SignalP.")
    .log_predictor_status("SignalP", "SKIPPED", "no AA FASTA extracted")
  } else {
    sp_out_dir_l <- file.path(out_dir, "signalp_out")
    dir.create(sp_out_dir_l, recursive = TRUE, showWarnings = FALSE)

    sp_out_dir_w <- w2l(sp_out_dir_l)

    has_sp6 <- tool_ok("signalp6")

    sp_result_file <- NULL
    sp_attempted <- FALSE
    sp_status <- NA_integer_

    if (has_sp6) {
      sp_attempted <- TRUE

      sp_cmd <- sprintf(
        "signalp6 --fastafile %s --organism eukarya --output_dir %s --format txt --mode fast",
        shQuote(sp_fa_w, type = "sh"),
        shQuote(sp_out_dir_w, type = "sh")
      )

      message("  Running SignalP 6...")
      sp_status <- run_tool(sp_cmd, tool_name = "SignalP6")

      if (sp_status == 0L) {
        candidates <- c(
          file.path(sp_out_dir_l, "prediction_results.txt"),
          file.path(sp_out_dir_l, "output.gff3")
        )

        found <- Filter(file.exists, candidates)
        if (length(found) > 0) sp_result_file <- found[1]
      }
    } else {
      message("  signalp6 not found in PATH / conda env. Skipping.")
      .log_predictor_status("SignalP", "SKIPPED", "signalp6 not found on PATH/conda env")
    }

    if (!is.null(sp_result_file)) {
      sp_import_ok <- FALSE

      switch_list <- tryCatch(
        {
          sl <- IsoformSwitchAnalyzeR::analyzeSignalIP(
            switchAnalyzeRlist = switch_list,
            pathToSignalPresultFile = sp_result_file,
            quiet = FALSE
          )

          sp_import_ok <- TRUE
          sl
        },
        error = function(e) {
          message("  analyzeSignalIP(): ", e$message)
          switch_list
        }
      )

      if (sp_import_ok) {
        message("  SignalP results imported successfully.")
        .log_predictor_status("SignalP", "OK", paste("result file:", sp_result_file))
      } else {
        message("  SignalP import failed; signal-peptide annotation NOT added to switch_list.")
        .log_predictor_status("SignalP", "IMPORT_FAILED", "analyzeSignalIP() errored")
      }
    } else if (sp_attempted && !is.na(sp_status) && sp_status != 0L) {
      message("  SignalP execution failed (exit ", sp_status, "). Skipping.")
      .log_predictor_status("SignalP", "FAILED", paste("exit code", sp_status))
    } else if (sp_attempted) {
      message("  SignalP exited successfully but no expected result file was found. Skipping.")
      .log_predictor_status("SignalP", "FAILED", "exit code 0 but no expected output file found")
    }
  }

  message("\n--- DeepTMHMM (transmembrane topology) ---")

  if (is.null(sp_fa_w)) {
    message("  No AA FASTA available. Skipping DeepTMHMM.")
    .log_predictor_status("DeepTMHMM", "SKIPPED", "no AA FASTA extracted")
  } else if (!tool_ok("biolib")) {
    message("  'biolib' not found in PATH / conda env. Skipping DeepTMHMM.")
    message("  (install inside the isoform_tools env with: pip install pybiolib)")
    .log_predictor_status("DeepTMHMM", "SKIPPED", "biolib not found on PATH/conda env (pip install pybiolib)")
  } else {
    tmhmm_out_dir_l <- file.path(out_dir, "deeptmhmm_out")
    dir.create(tmhmm_out_dir_l, recursive = TRUE, showWarnings = FALSE)
    tmhmm_out_dir_w <- w2l(tmhmm_out_dir_l)

    tmhmm_cmd <- sprintf(
      "cd %s && biolib run DTU/DeepTMHMM --fasta %s",
      shQuote(tmhmm_out_dir_w, type = "sh"),
      shQuote(sp_fa_w, type = "sh")
    )

    message("  Running DeepTMHMM via biolib (sequences are uploaded to BioLib's cloud servers)...")
    tmhmm_status <- run_tool(tmhmm_cmd, tool_name = "DeepTMHMM")

    tmhmm_gff <- file.path(tmhmm_out_dir_l, "biolib_results", "TMRs.gff3")

    if (tmhmm_status == 0L && file.exists(tmhmm_gff)) {
      tmhmm_import_ok <- FALSE

      switch_list <- tryCatch(
        {
          sl <- IsoformSwitchAnalyzeR::analyzeDeepTMHMM(
            switchAnalyzeRlist = switch_list,
            pathToDeepTMHMMresultFile = tmhmm_gff,
            showProgress = FALSE,
            quiet = FALSE
          )

          tmhmm_import_ok <- TRUE
          sl
        },
        error = function(e) {
          message("  analyzeDeepTMHMM(): ", e$message)
          switch_list
        }
      )

      if (tmhmm_import_ok) {
        message("  DeepTMHMM results imported successfully.")
        .log_predictor_status("DeepTMHMM", "OK", paste("result file:", tmhmm_gff))
      } else {
        message("  DeepTMHMM import failed; topology annotation NOT added to switch_list.")
        .log_predictor_status("DeepTMHMM", "IMPORT_FAILED", "analyzeDeepTMHMM() errored")
      }
    } else {
      message("  DeepTMHMM failed or result file not found (", tmhmm_gff, "). Skipping.")
      .log_predictor_status(
        "DeepTMHMM",
        "FAILED",
        sprintf("exit code %s, result file exists: %s", tmhmm_status, file.exists(tmhmm_gff))
      )
    }
  }

  message("\n--- Pfam domain annotation ---")

  if (is.null(pfam_fa_w)) {
    message("  No AA FASTA available. Skipping Pfam.")
    .log_predictor_status("Pfam", "SKIPPED", "no AA FASTA extracted")
  } else {
    pfam_done <- FALSE

    if (tool_ok("interproscan.sh")) {
      iprscan_xml_l <- file.path(out_dir, "interproscan.xml")
      iprscan_xml_w <- w2l(iprscan_xml_l)

      ips_cmd <- sprintf(
        "interproscan.sh -i %s -f XML -o %s -dp -appl Pfam -goterms -iprlookup -cpu %d",
        shQuote(pfam_fa_w, type = "sh"),
        shQuote(iprscan_xml_w, type = "sh"),
        n_cpu
      )

      message("  Running InterProScan...")
      ips_status <- run_tool(ips_cmd, tool_name = "Pfam_InterProScan")

      if (ips_status == 0L && file.exists(iprscan_xml_l)) {
        ips_import_ok <- FALSE

        switch_list <- tryCatch(
          {
            sl <- IsoformSwitchAnalyzeR::analyzeInterProScan(
              switchAnalyzeRlist = switch_list,
              pathToInterProScanResultFile = iprscan_xml_l,
              quiet = FALSE
            )

            ips_import_ok <- TRUE
            sl
          },
          error = function(e) {
            message("  analyzeInterProScan(): ", e$message)
            switch_list
          }
        )

        pfam_done <- ips_import_ok

        if (ips_import_ok) {
          message("  InterProScan Pfam results imported successfully.")
          .log_predictor_status("Pfam-InterProScan", "OK", paste("result file:", iprscan_xml_l))
        } else {
          message("  InterProScan import failed; trying hmmscan fallback...")
          .log_predictor_status("Pfam-InterProScan", "IMPORT_FAILED", "analyzeInterProScan() errored")
        }
      } else {
        message("  InterProScan failed (exit ", ips_status, "). Trying hmmscan fallback...")
        .log_predictor_status(
          "Pfam-InterProScan",
          "FAILED",
          sprintf("exit code %s, result file exists: %s", ips_status, file.exists(iprscan_xml_l))
        )
      }
    } else {
      .log_predictor_status("Pfam-InterProScan", "SKIPPED", "interproscan.sh not found on PATH/conda env")
    }

    if (!pfam_done) {
      if (!tool_ok("hmmscan")) {
        message("  hmmscan not found. Skipping Pfam analysis.")
        .log_predictor_status("Pfam-hmmscan", "SKIPPED", "hmmscan not found on PATH/conda env")
      } else {
        pfam_db_w <- .find_pfam_db(wsl_distro, via_wsl, active_conda, log_dir = log_dir)

        if (is.null(pfam_db_w)) {
          message("  Pfam-A.hmm not found. Run install_isoform_databases() first. Skipping.")
          .log_predictor_status("Pfam-hmmscan", "SKIPPED", "Pfam-A.hmm not found")
        } else {
          pfam_tbl_l <- file.path(out_dir, "pfam_domtblout.txt")
          pfam_tbl_w <- w2l(pfam_tbl_l)

          hm_cmd <- sprintf(
            "hmmscan --cpu %d --domtblout %s %s %s",
            n_cpu,
            shQuote(pfam_tbl_w, type = "sh"),
            shQuote(pfam_db_w, type = "sh"),
            shQuote(pfam_fa_w, type = "sh")
          )

          message("  Running hmmscan...")
          hm_status <- run_tool(hm_cmd, tool_name = "Pfam_hmmscan")

          if (hm_status == 0L && file.exists(pfam_tbl_l)) {
            hm_import_ok <- FALSE

            switch_list <- tryCatch(
              {
                sl <- IsoformSwitchAnalyzeR::analyzePFAM(
                  switchAnalyzeRlist = switch_list,
                  pathToPFAMresultFile = pfam_tbl_l,
                  quiet = FALSE
                )

                hm_import_ok <- TRUE
                sl
              },
              error = function(e) {
                message("  analyzePFAM(): ", e$message)
                switch_list
              }
            )

            if (hm_import_ok) {
              message("  hmmscan Pfam results imported successfully.")
              .log_predictor_status("Pfam-hmmscan", "OK", paste("result file:", pfam_tbl_l))
            } else {
              message("  hmmscan import failed; Pfam domain annotation NOT added to switch_list.")
              .log_predictor_status("Pfam-hmmscan", "IMPORT_FAILED", "analyzePFAM() errored")
            }
          } else {
            message("  hmmscan failed (exit ", hm_status, "). Skipping Pfam.")
            .log_predictor_status(
              "Pfam-hmmscan",
              "FAILED",
              sprintf("exit code %s, result file exists: %s", hm_status, file.exists(pfam_tbl_l))
            )
          }
        }
      }
    }
  }

  n_orf <- if (!is.null(switch_list$orfAnalysis)) nrow(switch_list$orfAnalysis) else 0L

  n_domain <- if (!is.null(switch_list$domainAnalysis)) {
    length(unique(switch_list$domainAnalysis$isoform_id))
  } else {
    0L
  }

  n_sigp <- if (!is.null(switch_list$signalPeptideAnalysis)) {
    length(unique(switch_list$signalPeptideAnalysis$isoform_id))
  } else {
    0L
  }

  n_cds <- if (!is.null(switch_list$isoformFeatures) &&
               "codingPotential" %in% colnames(switch_list$isoformFeatures)) {
    sum(!is.na(switch_list$isoformFeatures$codingPotential))
  } else {
    0L
  }

  n_topology <- if (!is.null(switch_list$topologyAnalysis)) {
    length(unique(switch_list$topologyAnalysis$isoform_id))
  } else {
    0L
  }

  summary_msg <- sprintf(
    "ORF analysis: %d transcripts with an ORF entry | Coding potential (CPAT): %d annotated | Domains (Pfam): %d transcripts annotated | Signal peptide: %d transcripts annotated | Topology (DeepTMHMM): %d transcripts annotated",
    n_orf,
    n_cds,
    n_domain,
    n_sigp,
    n_topology
  )

  message("\n--- Predictor summary ---\n   ", summary_msg)

  .log_predictor_status(
    "SUMMARY",
    if (n_domain > 0 || n_sigp > 0 || n_cds > 0 || n_topology > 0) "OK" else "EMPTY",
    summary_msg
  )

  message("\n========== PREDICTOR RUN SUMMARY ==========")
  message("  Log directory: ", log_dir)
  message("  Per-tool logs:")

  for (lf in list.files(log_dir, pattern = "\\.log$", full.names = TRUE)) {
    sz <- file.info(lf)$size
    message("     ", basename(lf), " (", format(round(sz / 1024, 1)), " KB)")
  }

  message("  Status summary: ", summary_msg)
  message("===========================================\n")

  summary_json <- list(
    timestamp = format(Sys.time()),
    log_dir = log_dir,
    n_orf = n_orf,
    n_cds = n_cds,
    n_domain = n_domain,
    n_sigp = n_sigp,
    n_topology = n_topology,
    tool_logs = list.files(log_dir, pattern = "\\.log$", full.names = FALSE)
  )

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      summary_json,
      file.path(log_dir, "predictor_run_summary.json"),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }

  return(switch_list)
}
