# mod_isoform_predictors.R -- external predictor tools for isoform switch analysis
# (CPAT, SignalP, DeepTMHMM, Pfam/hmmscan/InterProScan), run via WSL

# ---------------------------------------------------------------------------
# Pfam result-format converters
#
# IsoformSwitchAnalyzeR::analyzePFAM() only parses the pfam_scan.pl 16-column
# layout (seq_id, alignment_start/end, envelope_start/end, hmm_acc, hmm_name,
# type, hmm_start/end, hmm_length, bit_score, E_value, significant, clan,
# [residue]) -- confirmed against the current source (R/analyze_external_
# sequence_analysis.R): it hard-checks that column 6 matches '^PF|^PB' and
# that column 15 (or 16) contains at least one '^CL...' clan accession,
# stopping with "not recogniced as a pfam output" otherwise. Neither
# interproscan.sh's TSV output nor hmmscan's --domtblout use that layout, so
# both need to be remapped before analyzePFAM() will accept them. There is
# also no analyzeInterProScan() anywhere in IsoformSwitchAnalyzeR (checked
# the NAMESPACE and full R/ source directly) -- calling it always errored
# with "could not find function", regardless of whether InterProScan itself
# ran fine, which is why this path funnels into analyzePFAM() too.
#
# Caveat: neither source gives us real Pfam clan membership (that's a
# separate Pfam-A.clans.tsv mapping we don't fetch), so clan is filled with
# a sentinel ("CL_unknown") that only exists to satisfy the '^CL' style
# check -- it is not a real clan id. hmmscan's domtblout carries every other
# column (coordinates, hmm length, bit score, E-value) through exactly;
# InterProScan's TSV doesn't report hmm-model coordinates, hmm length, or
# bit score at all, so those are left NA for that path. If exact domain
# boundaries/scores matter, prefer the hmmscan path.
# ---------------------------------------------------------------------------

#' @keywords internal
.write_pfam_scan_table <- function(df16, dest_path) {
  header <- paste(
    "# <seq id> <alignment start> <alignment end> <envelope start> <envelope end>",
    "<hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score>",
    "<E-value> <significance> <clan> <predicted_active_site_residues>"
  )

  con <- file(dest_path, "wb")
  writeLines(header, con, sep = "\n")
  close(con)

  utils::write.table(
    df16,
    dest_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    append = TRUE
  )

  dest_path
}

#' Convert an InterProScan '-f tsv' (Pfam-filtered) file into pfam_scan.pl format
#' @keywords internal
.pfam_from_interproscan_tsv <- function(ips_tsv_path, dest_path) {
  if (!file.exists(ips_tsv_path) || file.info(ips_tsv_path)$size == 0) return(NULL)

  raw <- tryCatch(
    utils::read.delim(
      ips_tsv_path,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      fill = TRUE,
      quote = "",
      col.names = paste0("V", 1:15)
    ),
    error = function(e) NULL
  )

  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # -appl Pfam already restricts the run to Pfam, and status "T" = successful
  # match; filter on both defensively rather than assuming.
  keep <- raw$V4 == "Pfam" & grepl("^PF|^PB", raw$V5) & toupper(raw$V10) == "T"
  raw <- raw[keep, , drop = FALSE]
  if (nrow(raw) == 0) return(NULL)

  clean_text <- function(x) {
    x <- gsub("[\t\r\n]+", " ", x)
    trimws(x)
  }

  hmm_name <- clean_text(raw$V6)
  hmm_name[!nzchar(hmm_name)] <- raw$V5[!nzchar(hmm_name)]

  out <- data.frame(
    V1  = raw$V1,
    V2  = raw$V7,
    V3  = raw$V8,
    V4  = raw$V7,
    V5  = raw$V8,
    V6  = raw$V5,
    V7  = hmm_name,
    V8  = "Domain",
    V9  = NA_character_,
    V10 = NA_character_,
    V11 = NA_character_,
    V12 = NA_character_,
    V13 = raw$V9,
    V14 = 1,
    V15 = "CL_unknown",
    V16 = "",
    stringsAsFactors = FALSE
  )

  .write_pfam_scan_table(out, dest_path)
}

#' Convert hmmscan '--domtblout' output into pfam_scan.pl format
#' @keywords internal
.pfam_from_hmmscan_domtblout <- function(domtblout_path, dest_path) {
  if (!file.exists(domtblout_path) || file.info(domtblout_path)$size == 0) return(NULL)

  lines <- readLines(domtblout_path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) return(NULL)

  # domtblout is whitespace-delimited with 23 fixed columns; column 23
  # ("description of target") can itself contain spaces, so cap the split.
  fields <- strsplit(trimws(lines), "\\s+")
  fields <- lapply(fields, function(f) {
    if (length(f) > 22) f[1:22] else f
  })

  ok <- vapply(fields, length, integer(1)) == 22
  fields <- fields[ok]
  if (length(fields) == 0) return(NULL)

  raw <- as.data.frame(
    do.call(rbind, fields),
    stringsAsFactors = FALSE
  )
  colnames(raw) <- c(
    "target_name", "target_acc", "tlen", "query_name", "query_acc", "qlen",
    "seq_evalue", "seq_score", "seq_bias", "dom_num", "dom_of",
    "c_evalue", "i_evalue", "dom_score", "dom_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"
  )

  raw <- raw[grepl("^PF|^PB", raw$target_acc), , drop = FALSE]
  if (nrow(raw) == 0) return(NULL)

  out <- data.frame(
    V1  = raw$query_name,
    V2  = raw$ali_from,
    V3  = raw$ali_to,
    V4  = raw$env_from,
    V5  = raw$env_to,
    V6  = raw$target_acc,
    V7  = raw$target_name,
    V8  = "Domain",
    V9  = raw$hmm_from,
    V10 = raw$hmm_to,
    V11 = raw$tlen,
    V12 = raw$dom_score,
    V13 = raw$i_evalue,
    V14 = 1,
    V15 = "CL_unknown",
    V16 = "",
    stringsAsFactors = FALSE
  )

  .write_pfam_scan_table(out, dest_path)
}

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

      # CPAT's --log-file defaults to a bare "CPAT_run_info.log" (see `cpat
      # -h`), written relative to the shell's cwd at invocation time -- NOT
      # next to the -o outputs. Since that cwd is wherever WSL starts the
      # script (not out_dir), the log was landing outside the results
      # folder entirely. Pointing it at out_dir_w explicitly keeps it next
      # to cpat_out.* like every other tool's log here. Only cpat (3.x)
      # is confirmed to support this flag, so it's gated to that binary --
      # the legacy run_cpat.py (2.x) command below is left as-is.
      cpat_log_w <- paste0(out_dir_w, "/CPAT_run_info.log")

      cpat_cmd <- if (has_cpat3) {
        sprintf(
          "%s -g %s -x %s -d %s -o %s --log-file %s",
          cpat_exe,
          shQuote(cpat_fa_w, type = "sh"),
          shQuote(hexamer_w, type = "sh"),
          shQuote(logit_w, type = "sh"),
          shQuote(cpat_prefix_w, type = "sh"),
          shQuote(cpat_log_w, type = "sh")
        )
      } else {
        sprintf(
          "%s -g %s -x %s -d %s -o %s",
          cpat_exe,
          shQuote(cpat_fa_w, type = "sh"),
          shQuote(hexamer_w, type = "sh"),
          shQuote(logit_w, type = "sh"),
          shQuote(cpat_prefix_w, type = "sh")
        )
      }

      message("  Running: ", cpat_exe)
      cpat_status <- run_tool(cpat_cmd, tool_name = "CPAT")

      cpat_result <- if (has_cpat3) {
        paste0(cpat_prefix_l, ".ORF_prob.best.tsv")
      } else {
        paste0(cpat_prefix_l, ".r")
      }

      if (cpat_status == 0L && file.exists(cpat_result)) {
        cpat_import_ok <- FALSE
        cpat_err <- NULL

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
            cpat_err <<- e$message
            message("  analyzeCPAT(): ", e$message)
            switch_list
          }
        )

        if (cpat_import_ok) {
          message("  CPAT results imported successfully.")
          .log_predictor_status("CPAT", "OK", paste("result file:", cpat_result))
        } else {
          message("  CPAT import failed; coding-potential annotation NOT added to switch_list.")
          .log_predictor_status("CPAT", "IMPORT_FAILED", paste("analyzeCPAT() errored:", cpat_err))
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
      sp_err <- NULL

      switch_list <- tryCatch(
        {
          # NB: the function is analyzeSignalP() (no "I") -- analyzeSignalIP()
          # has never existed as an exported IsoformSwitchAnalyzeR function;
          # calling it always errored with "could not find function",
          # independent of whether SignalP itself ran fine.
          sl <- IsoformSwitchAnalyzeR::analyzeSignalP(
            switchAnalyzeRlist = switch_list,
            pathToSignalPresultFile = sp_result_file,
            quiet = FALSE
          )

          sp_import_ok <- TRUE
          sl
        },
        error = function(e) {
          sp_err <<- e$message
          message("  analyzeSignalP(): ", e$message)
          switch_list
        }
      )

      if (sp_import_ok) {
        message("  SignalP results imported successfully.")
        .log_predictor_status("SignalP", "OK", paste("result file:", sp_result_file))
      } else {
        message("  SignalP import failed; signal-peptide annotation NOT added to switch_list.")
        .log_predictor_status("SignalP", "IMPORT_FAILED", paste("analyzeSignalP() errored:", sp_err))
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
      tmhmm_err <- NULL

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
          tmhmm_err <<- e$message
          message("  analyzeDeepTMHMM(): ", e$message)
          switch_list
        }
      )

      if (tmhmm_import_ok) {
        message("  DeepTMHMM results imported successfully.")
        .log_predictor_status("DeepTMHMM", "OK", paste("result file:", tmhmm_gff))
      } else {
        message("  DeepTMHMM import failed; topology annotation NOT added to switch_list.")
        .log_predictor_status("DeepTMHMM", "IMPORT_FAILED", paste("analyzeDeepTMHMM() errored:", tmhmm_err))
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
      # TSV instead of XML: IsoformSwitchAnalyzeR has no analyzeInterProScan()
      # (verified against the current NAMESPACE / R source), so the output
      # has to be reformatted into pfam_scan.pl's layout ourselves before
      # analyzePFAM() will read it -- TSV is a plain, documented column
      # format we can parse directly, XML would need an extra dependency
      # (xml2) for no benefit here.
      iprscan_tsv_l <- file.path(out_dir, "interproscan.tsv")
      iprscan_tsv_w <- w2l(iprscan_tsv_l)

      ips_cmd <- sprintf(
        "interproscan.sh -i %s -f tsv -o %s -dp -appl Pfam -goterms -iprlookup -cpu %d",
        shQuote(pfam_fa_w, type = "sh"),
        shQuote(iprscan_tsv_w, type = "sh"),
        n_cpu
      )

      message("  Running InterProScan...")
      ips_status <- run_tool(ips_cmd, tool_name = "Pfam_InterProScan")

      if (ips_status == 0L && file.exists(iprscan_tsv_l)) {
        ips_import_ok <- FALSE
        ips_err <- NULL

        ips_pfam_l <- file.path(out_dir, "interproscan_as_pfam_scan.txt")
        ips_converted <- .pfam_from_interproscan_tsv(iprscan_tsv_l, ips_pfam_l)

        if (is.null(ips_converted)) {
          message("  No Pfam hits could be parsed from InterProScan's output; trying hmmscan fallback...")
          .log_predictor_status(
            "Pfam-InterProScan",
            "IMPORT_FAILED",
            paste("no Pfam rows parsed from", iprscan_tsv_l)
          )
        } else {
          switch_list <- tryCatch(
            {
              sl <- IsoformSwitchAnalyzeR::analyzePFAM(
                switchAnalyzeRlist = switch_list,
                pathToPFAMresultFile = ips_converted,
                quiet = FALSE
              )

              ips_import_ok <- TRUE
              sl
            },
            error = function(e) {
              ips_err <<- e$message
              message("  analyzePFAM() [from InterProScan]: ", e$message)
              switch_list
            }
          )

          pfam_done <- ips_import_ok

          if (ips_import_ok) {
            message("  InterProScan Pfam results imported successfully.")
            .log_predictor_status("Pfam-InterProScan", "OK", paste("result file:", ips_converted))
          } else {
            message("  InterProScan import failed; trying hmmscan fallback...")
            .log_predictor_status(
              "Pfam-InterProScan",
              "IMPORT_FAILED",
              paste("analyzePFAM() errored:", ips_err)
            )
          }
        }
      } else {
        message("  InterProScan failed (exit ", ips_status, "). Trying hmmscan fallback...")
        .log_predictor_status(
          "Pfam-InterProScan",
          "FAILED",
          sprintf("exit code %s, result file exists: %s", ips_status, file.exists(iprscan_tsv_l))
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
            hm_err <- NULL

            # hmmscan's --domtblout is a different, incompatible column
            # layout from what analyzePFAM() parses (pfam_scan.pl's 16
            # columns) -- feeding it in raw is what was erroring here even
            # though hmmscan itself succeeded. Convert first.
            hm_pfam_l <- file.path(out_dir, "hmmscan_as_pfam_scan.txt")
            hm_converted <- .pfam_from_hmmscan_domtblout(pfam_tbl_l, hm_pfam_l)

            if (is.null(hm_converted)) {
              message("  No Pfam hits could be parsed from hmmscan's domtblout output.")
              .log_predictor_status(
                "Pfam-hmmscan",
                "IMPORT_FAILED",
                paste("no Pfam rows parsed from", pfam_tbl_l)
              )
            } else {
              switch_list <- tryCatch(
                {
                  sl <- IsoformSwitchAnalyzeR::analyzePFAM(
                    switchAnalyzeRlist = switch_list,
                    pathToPFAMresultFile = hm_converted,
                    quiet = FALSE
                  )

                  hm_import_ok <- TRUE
                  sl
                },
                error = function(e) {
                  hm_err <<- e$message
                  message("  analyzePFAM(): ", e$message)
                  switch_list
                }
              )

              if (hm_import_ok) {
                message("  hmmscan Pfam results imported successfully.")
                .log_predictor_status("Pfam-hmmscan", "OK", paste("result file:", hm_converted))
              } else {
                message("  hmmscan import failed; Pfam domain annotation NOT added to switch_list.")
                .log_predictor_status(
                  "Pfam-hmmscan",
                  "IMPORT_FAILED",
                  paste("analyzePFAM() errored:", hm_err)
                )
              }
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
