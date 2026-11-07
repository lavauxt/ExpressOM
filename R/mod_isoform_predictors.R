# mod_isoform_predictors.R -- external predictor tools for isoform switch analysis
# (CPAT, SignalP, DeepTMHMM, Pfam/hmmscan/InterProScan), run via WSL

# Pfam result-format converters: analyzePFAM() only parses the pfam_scan.pl
# 16-column layout (requires col 6 = hmm_acc matching '^PF|^PB' and a
# '^CL...' clan id in col 15/16). Neither interproscan.sh's TSV nor
# hmmscan's --domtblout use that layout, and analyzeInterProScan() doesn't
# exist in IsoformSwitchAnalyzeR at all, so both tools' output is remapped
# here before analyzePFAM() will accept either.
#
# Caveat: neither source gives real Pfam clan membership, so clan is filled
# with a sentinel ("CL_unknown") that only exists to satisfy the '^CL'
# check. hmmscan's domtblout carries every other column through exactly;
# InterProScan's TSV doesn't report hmm-model coordinates/length/bit score,
# so those are left NA for that path -- prefer hmmscan if exact domain
# boundaries/scores matter.

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
  if (!file.exists(ips_tsv_path) || file.info(ips_tsv_path)$size == 0) {
    return(NULL)
  }

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

  if (is.null(raw) || nrow(raw) == 0) {
    return(NULL)
  }

  # -appl Pfam already restricts the run to Pfam, and status "T" = successful
  # match; filter on both defensively rather than assuming.
  keep <- raw$V4 == "Pfam" & grepl("^PF|^PB", raw$V5) & toupper(raw$V10) == "T"
  raw <- raw[keep, , drop = FALSE]
  if (nrow(raw) == 0) {
    return(NULL)
  }

  clean_text <- function(x) {
    x <- gsub("[\t\r\n]+", " ", x)
    trimws(x)
  }

  hmm_name <- clean_text(raw$V6)
  hmm_name[!nzchar(hmm_name)] <- raw$V5[!nzchar(hmm_name)]

  out <- data.frame(
    V1 = raw$V1,
    V2 = raw$V7,
    V3 = raw$V8,
    V4 = raw$V7,
    V5 = raw$V8,
    V6 = raw$V5,
    V7 = hmm_name,
    V8 = "Domain",
    V9 = NA_character_,
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
  if (!file.exists(domtblout_path) || file.info(domtblout_path)$size == 0) {
    return(NULL)
  }

  lines <- readLines(domtblout_path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    return(NULL)
  }

  # domtblout is whitespace-delimited with 23 fixed columns; column 23
  # ("description of target") can itself contain spaces, so cap the split.
  fields <- strsplit(trimws(lines), "\\s+")
  fields <- lapply(fields, function(f) {
    if (length(f) > 22) f[1:22] else f
  })

  ok <- vapply(fields, length, integer(1)) == 22
  fields <- fields[ok]
  if (length(fields) == 0) {
    return(NULL)
  }

  raw <- as.data.frame(
    do.call(rbind, fields),
    stringsAsFactors = FALSE
  )
  colnames(raw) <- c(
    "target_name",
    "target_acc",
    "tlen",
    "query_name",
    "query_acc",
    "qlen",
    "seq_evalue",
    "seq_score",
    "seq_bias",
    "dom_num",
    "dom_of",
    "c_evalue",
    "i_evalue",
    "dom_score",
    "dom_bias",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to",
    "env_from",
    "env_to",
    "acc"
  )

  raw <- raw[grepl("^PF|^PB", raw$target_acc), , drop = FALSE]
  if (nrow(raw) == 0) {
    return(NULL)
  }

  out <- data.frame(
    V1 = raw$query_name,
    V2 = raw$ali_from,
    V3 = raw$ali_to,
    V4 = raw$env_from,
    V5 = raw$env_to,
    V6 = raw$target_acc,
    V7 = raw$target_name,
    V8 = "Domain",
    V9 = raw$hmm_from,
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

.run_external_predictors <- function(
  switch_list,
  fasta_file,
  out_dir,
  use_wsl,
  wsl_distro,
  isoform_obj,
  save_dir,
  n_cpu = NULL,
  log_dir = NULL,
  organism = NULL
) {
  is_windows <- .Platform$OS.type == "windows"
  via_wsl <- is_windows && isTRUE(use_wsl)

  if (is.null(log_dir)) {
    log_dir <- file.path(out_dir, "Log")
  }
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Each tool's own raw output (CPAT's ORF/probability tables, SignalP's
  # prediction_results.txt, the InterProScan/hmmscan Pfam tables and their
  # pfam_scan.pl-converted counterparts) goes into saved/<Tool>/ instead of
  # out_dir root, so results_.../IsoformSwitch stops accumulating a flat
  # pile of per-tool files next to Plots/, Log/, sequences/, etc.
  tool_save_dir <- function(family) {
    d <- if (!is.null(save_dir)) {
      file.path(save_dir, family)
    } else {
      file.path(out_dir, family)
    }
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
    d
  }

  predictor_log_path <- file.path(log_dir, "predictor_status.log")

  message(
    "Predictor run logging to: ",
    predictor_log_path,
    " (step-by-step status) and ",
    file.path(log_dir, "wsl_commands.log"),
    " (every shell command, in real time)."
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

  cat(
    "\n\n############################################################\n",
    file = predictor_log_path,
    append = TRUE
  )
  cat(
    sprintf("# NEW PREDICTOR RUN: %s\n", format(Sys.time())),
    file = predictor_log_path,
    append = TRUE
  )
  cat(
    sprintf("# Output dir: %s\n", out_dir),
    file = predictor_log_path,
    append = TRUE
  )
  cat(
    "############################################################\n\n",
    file = predictor_log_path,
    append = TRUE
  )

  message("Predictor status log: ", predictor_log_path)

  if (is.null(n_cpu) || !is.finite(n_cpu) || n_cpu < 1) {
    detected <- tryCatch(
      parallel::detectCores(logical = TRUE),
      error = function(e) NA_integer_
    )
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
    conda_sh <- .find_conda_sh(
      wsl_distro = wsl_distro,
      use_wsl = via_wsl,
      log_dir = log_dir
    )
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

    if (!has_conda_env) {
      conda_sh <- NULL
    }

    assign(
      cache_key,
      if (has_conda_env) conda_sh else "",
      envir = .expressom_cache
    )

    if (has_conda_env) {
      message("  Using conda env 'isoform_tools' (", conda_sh, ")")
    } else {
      message("  No 'isoform_tools' env found. Using system PATH.")
    }
  }

  active_conda <- if (has_conda_env) conda_sh else NULL

  run_tool <- function(cmd_str, show_stderr = TRUE, tool_name = "tool") {
    tool_family <- sub("_.*$", "", tool_name)
    tool_log_dir <- file.path(log_dir, tool_family)
    if (!dir.exists(tool_log_dir)) {
      dir.create(tool_log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    tool_log <- file.path(tool_log_dir, paste0(tool_name, ".log"))

    cat(
      sprintf("\n===== %s — %s =====\n", tool_name, format(Sys.time())),
      file = tool_log,
      append = TRUE
    )
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

    cat(
      sprintf("--- exit code: %d ---\n", as.integer(exit_code)),
      file = tool_log,
      append = TRUE
    )

    if (length(res) > 0) {
      msg <- paste(res, collapse = "\n")

      if (nzchar(msg)) {
        msg_lines <- strsplit(msg, "\n")[[1]]

        if (length(msg_lines) > 30) {
          message(
            "  [",
            tool_name,
            "] Output (first 15 + last 15 of ",
            length(msg_lines),
            " lines; full log: ",
            tool_log,
            "):"
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
      "  [",
      tool_name,
      "] Exit code: ",
      as.integer(exit_code),
      " — ",
      if (exit_code == 0L) "SUCCESS" else "FAILED"
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
    if (is.null(p) || length(p) == 0 || !any(nzchar(p))) {
      return(p)
    }
    if (length(p) > 1) {
      return(vapply(p, w2l, character(1), USE.NAMES = FALSE))
    }
    if (!nzchar(p)) {
      return(p)
    }
    if (is_windows && use_wsl) .to_wsl_path(p, wsl_distro) else p
  }

  syms <- isoform_obj$gene_map$symbol[!is.na(isoform_obj$gene_map$symbol)]
  first_sym <- if (length(syms) > 0) syms[1] else ""

  organism_cpat <- if (
    !is.null(organism) && grepl("Homo sapiens", organism, ignore.case = TRUE)
  ) {
    "Human"
  } else if (
    !is.null(organism) && grepl("Mus musculus", organism, ignore.case = TRUE)
  ) {
    "Mouse"
  } else {
    if (!is.null(organism)) {
      message(
        "  Organism '",
        organism,
        "' is not Human or Mouse (CPAT only ships prebuilt logit ",
        "models for those two); falling back to a gene-symbol-casing guess for the CPAT model choice."
      )
    }
    if (grepl("^[A-Z]+$", first_sym) && nchar(first_sym) > 2) {
      "Human"
    } else {
      "Mouse"
    }
  }

  message(
    "Organism detected for CPAT: ",
    organism_cpat,
    if (!is.null(organism)) {
      paste0(" (from EnsDb: ", organism, ")")
    } else {
      " (guessed from gene symbol casing -- pass organism = ... for a reliable result)"
    }
  )

  seq_dir <- file.path(out_dir, "sequences")
  dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

  nt_fa_local <- file.path(seq_dir, "isoform_NT.fa")
  aa_fa_local <- file.path(seq_dir, "isoform_AA.fa")

  if (!file.exists(nt_fa_local) || !file.exists(aa_fa_local)) {
    message(
      "Extracting NT / AA sequences from SwitchList for external tools..."
    )

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
        written <- list.files(
          seq_dir,
          pattern = "\\.fa(sta)?$",
          full.names = TRUE,
          ignore.case = TRUE
        )
        written <- setdiff(written, c(nt_fa_local, aa_fa_local))
        written <- written[
          !grepl("complete|split|subset", basename(written), ignore.case = TRUE)
        ]

        nt_src <- written[grepl(
          "(^|[_.-])nt([_.-]|$)",
          basename(written),
          ignore.case = TRUE
        )]
        aa_src <- written[grepl(
          "(^|[_.-])aa([_.-]|$)",
          basename(written),
          ignore.case = TRUE
        )]

        if (length(nt_src) >= 1 && !file.exists(nt_fa_local)) {
          file.rename(nt_src[1], nt_fa_local)
        }
        if (length(aa_src) >= 1 && !file.exists(aa_fa_local)) {
          file.rename(aa_src[1], aa_fa_local)
        }

        if (!file.exists(nt_fa_local) || !file.exists(aa_fa_local)) {
          message(
            "  extractSequence() ran but expected output was not found by pattern in ",
            seq_dir,
            " (files present: ",
            paste(basename(list.files(seq_dir)), collapse = ", "),
            ")"
          )
        }
      },
      error = function(e) {
        message("  Could not extract sequences from SwitchList: ", e$message)
        message(
          "  CPAT will use the provided fasta_file; SignalP / Pfam may be skipped."
        )
        .log_predictor_status("extractSequence", "FAILED", e$message)
      }
    )
  }

  n_nt <- if (file.exists(nt_fa_local)) {
    length(grep("^>", readLines(nt_fa_local, warn = FALSE)))
  } else {
    0L
  }
  n_aa <- if (file.exists(aa_fa_local)) {
    length(grep("^>", readLines(aa_fa_local, warn = FALSE)))
  } else {
    0L
  }

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
        "  fasta_file has ",
        length(fasta_file),
        " entries (e.g. combined cDNA + ncRNA ",
        "references); CPAT takes a single FASTA, so only the first (",
        fasta_file[1],
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
      "  Hexamer table (",
      organism_cpat,
      "_Hexamer.tsv) not found in $CPAT_DATA / $HOME/.cpat_data. ",
      "Run install_isoform_databases() first. Skipping CPAT."
    )
    .log_predictor_status(
      "CPAT",
      "SKIPPED",
      "hexamer table not found (run install_isoform_databases())"
    )
  } else if (is.null(logit_local)) {
    message(
      "  Logit model not found. Run install_isoform_databases() first. Skipping CPAT."
    )
    .log_predictor_status(
      "CPAT",
      "SKIPPED",
      "logit model not found (run install_isoform_databases())"
    )
  } else {
    hexamer_w <- hexamer_local
    logit_w <- logit_local
    cpat_prefix_l <- file.path(tool_save_dir("CPAT"), "cpat_out")
    cpat_prefix_w <- w2l(cpat_prefix_l)

    has_cpat3 <- tool_ok("cpat")
    has_cpat2 <- if (!has_cpat3) tool_ok("run_cpat.py") else FALSE

    if (!has_cpat3 && !has_cpat2) {
      message("  Neither 'cpat' nor 'run_cpat.py' found. Skipping.")
      .log_predictor_status(
        "CPAT",
        "SKIPPED",
        "neither 'cpat' nor 'run_cpat.py' found on PATH/conda env"
      )
    } else {
      cpat_exe <- if (has_cpat3) "cpat" else "run_cpat.py"

      cpat_log_dir_l <- file.path(log_dir, "CPAT")
      if (!dir.exists(cpat_log_dir_l)) {
        dir.create(cpat_log_dir_l, recursive = TRUE, showWarnings = FALSE)
      }

      # CPAT's --log-file defaults to a bare "CPAT_run_info.log" (see `cpat
      # -h`), written relative to the shell's cwd at invocation time -- NOT
      # next to the -o outputs. Since that cwd is wherever WSL starts the
      # script (not out_dir), the log was landing outside the results
      # folder entirely. Pointing it at Log/Isoform/CPAT/ explicitly keeps
      # it grouped with this tool's other logs. Only cpat (3.x) is
      # confirmed to support this flag, so it's gated to that binary -- the
      # legacy run_cpat.py (2.x) command below is left as-is.
      cpat_log_w <- paste0(w2l(cpat_log_dir_l), "/CPAT_run_info.log")

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

        # CPAT's own published optimum coding-probability cutoffs (from the
        # CPAT paper's TG-ROC analysis): 0.364 for human, 0.44 for mouse.
        coding_cutoff <- if (organism_cpat == "Human") 0.364 else 0.44

        switch_list <- tryCatch(
          {
            sl <- IsoformSwitchAnalyzeR::analyzeCPAT(
              switchAnalyzeRlist = switch_list,
              pathToCPATresultFile = cpat_result,
              codingCutoff = coding_cutoff,
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
          .log_predictor_status(
            "CPAT",
            "OK",
            paste("result file:", cpat_result)
          )
        } else {
          message(
            "  CPAT import failed; coding-potential annotation NOT added to switch_list."
          )
          .log_predictor_status(
            "CPAT",
            "IMPORT_FAILED",
            paste("analyzeCPAT() errored:", cpat_err)
          )
        }
      } else {
        message(
          "  CPAT failed or result file not found (",
          basename(cpat_result),
          "). Skipping."
        )
        .log_predictor_status(
          "CPAT",
          "FAILED",
          sprintf(
            "exit code %s, result file exists: %s",
            cpat_status,
            file.exists(cpat_result)
          )
        )
      }
    }
  }

  message("\n--- SignalP (signal peptide prediction) ---")

  if (is.null(sp_fa_w)) {
    message("  No AA FASTA available. Skipping SignalP.")
    .log_predictor_status("SignalP", "SKIPPED", "no AA FASTA extracted")
  } else {
    sp_out_dir_l <- file.path(tool_save_dir("SignalP"), "signalp_out")
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
      .log_predictor_status(
        "SignalP",
        "SKIPPED",
        "signalp6 not found on PATH/conda env"
      )
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
        .log_predictor_status(
          "SignalP",
          "OK",
          paste("result file:", sp_result_file)
        )
      } else {
        message(
          "  SignalP import failed; signal-peptide annotation NOT added to switch_list."
        )
        .log_predictor_status(
          "SignalP",
          "IMPORT_FAILED",
          paste("analyzeSignalP() errored:", sp_err)
        )
      }
    } else if (sp_attempted && !is.na(sp_status) && sp_status != 0L) {
      message("  SignalP execution failed (exit ", sp_status, "). Skipping.")
      .log_predictor_status("SignalP", "FAILED", paste("exit code", sp_status))
    } else if (sp_attempted) {
      message(
        "  SignalP exited successfully but no expected result file was found. Skipping."
      )
      .log_predictor_status(
        "SignalP",
        "FAILED",
        "exit code 0 but no expected output file found"
      )
    }
  }

  message("\n--- DeepTMHMM (transmembrane topology) ---")

  if (is.null(sp_fa_w)) {
    message("  No AA FASTA available. Skipping DeepTMHMM.")
    .log_predictor_status("DeepTMHMM", "SKIPPED", "no AA FASTA extracted")
  } else if (!tool_ok("biolib")) {
    message("  'biolib' not found in PATH / conda env. Skipping DeepTMHMM.")
    message(
      "  (install inside the isoform_tools env with: pip install pybiolib)"
    )
    .log_predictor_status(
      "DeepTMHMM",
      "SKIPPED",
      "biolib not found on PATH/conda env (pip install pybiolib)"
    )
  } else {
    tmhmm_out_dir_l <- file.path(out_dir, "deeptmhmm_out")
    dir.create(tmhmm_out_dir_l, recursive = TRUE, showWarnings = FALSE)
    tmhmm_out_dir_w <- w2l(tmhmm_out_dir_l)

    tmhmm_cmd <- sprintf(
      "cd %s && biolib run DTU/DeepTMHMM --fasta %s",
      shQuote(tmhmm_out_dir_w, type = "sh"),
      shQuote(sp_fa_w, type = "sh")
    )

    message(
      "  Running DeepTMHMM via biolib (sequences are uploaded to BioLib's cloud servers)..."
    )
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
        .log_predictor_status(
          "DeepTMHMM",
          "OK",
          paste("result file:", tmhmm_gff)
        )
      } else {
        message(
          "  DeepTMHMM import failed; topology annotation NOT added to switch_list."
        )
        .log_predictor_status(
          "DeepTMHMM",
          "IMPORT_FAILED",
          paste("analyzeDeepTMHMM() errored:", tmhmm_err)
        )
      }
    } else {
      message(
        "  DeepTMHMM failed or result file not found (",
        tmhmm_gff,
        "). Skipping."
      )
      .log_predictor_status(
        "DeepTMHMM",
        "FAILED",
        sprintf(
          "exit code %s, result file exists: %s",
          tmhmm_status,
          file.exists(tmhmm_gff)
        )
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
      iprscan_tsv_l <- file.path(tool_save_dir("Pfam"), "interproscan.tsv")
      iprscan_tsv_w <- w2l(iprscan_tsv_l)

      ips_temp_dir_w <- paste0(
        "$HOME/.isoform_tools_tmp/ips_",
        basename(tempfile(""))
      )

      ips_cmd <- sprintf(
        paste(
          "set +e; mkdir -p %1$s &&",
          "interproscan.sh -i %2$s -f tsv -o %3$s -dp -appl Pfam -goterms -iprlookup -cpu %4$d -T %1$s;",
          "_ips_status=$?; set -e; rm -rf %1$s; exit $_ips_status"
        ),
        .dq_home(ips_temp_dir_w),
        shQuote(pfam_fa_w, type = "sh"),
        shQuote(iprscan_tsv_w, type = "sh"),
        n_cpu
      )

      message("  Running InterProScan...")
      ips_status <- run_tool(ips_cmd, tool_name = "Pfam_InterProScan")

      # InterProScan can report success (exit 0) with a technically-existing
      # but still-empty output file, if the write to a WSL-mounted
      # (/mnt/c/...) destination hasn't fully synced by the time this check
      # runs. A short wait-and-recheck absorbs that without silently
      # mistaking it for "genuinely zero Pfam hits" or a parsing problem.
      ips_output_empty <- file.exists(iprscan_tsv_l) &&
        isTRUE(file.info(iprscan_tsv_l)$size == 0)

      if (ips_status == 0L && ips_output_empty) {
        Sys.sleep(2)
        ips_output_empty <- file.exists(iprscan_tsv_l) &&
          isTRUE(file.info(iprscan_tsv_l)$size == 0)
      }

      if (ips_status == 0L && ips_output_empty) {
        message(
          "  InterProScan reported success but its output file is still empty (0 bytes) after a ",
          "brief re-check -- this looks like a filesystem sync issue writing to a WSL-mounted path, ",
          "not a parsing problem. Trying hmmscan fallback..."
        )
        .log_predictor_status(
          "Pfam-InterProScan",
          "FAILED",
          paste(
            "exit code 0 but output file is empty (0 bytes):",
            iprscan_tsv_l
          )
        )
      } else if (ips_status == 0L && file.exists(iprscan_tsv_l)) {
        ips_import_ok <- FALSE
        ips_err <- NULL

        ips_pfam_l <- file.path(
          tool_save_dir("Pfam"),
          "interproscan_as_pfam_scan.txt"
        )
        ips_converted <- .pfam_from_interproscan_tsv(iprscan_tsv_l, ips_pfam_l)

        if (is.null(ips_converted)) {
          message(
            "  No Pfam hits could be parsed from InterProScan's output; trying hmmscan fallback..."
          )
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
            .log_predictor_status(
              "Pfam-InterProScan",
              "OK",
              paste("result file:", ips_converted)
            )
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
        message(
          "  InterProScan failed (exit ",
          ips_status,
          "). Trying hmmscan fallback..."
        )
        .log_predictor_status(
          "Pfam-InterProScan",
          "FAILED",
          sprintf(
            "exit code %s, result file exists: %s",
            ips_status,
            file.exists(iprscan_tsv_l)
          )
        )
      }
    } else {
      .log_predictor_status(
        "Pfam-InterProScan",
        "SKIPPED",
        "interproscan.sh not found on PATH/conda env"
      )
    }

    if (!pfam_done) {
      if (!tool_ok("hmmscan")) {
        message("  hmmscan not found. Skipping Pfam analysis.")
        .log_predictor_status(
          "Pfam-hmmscan",
          "SKIPPED",
          "hmmscan not found on PATH/conda env"
        )
      } else {
        pfam_db_w <- .find_pfam_db(
          wsl_distro,
          via_wsl,
          active_conda,
          log_dir = log_dir
        )

        if (is.null(pfam_db_w)) {
          message(
            "  Pfam-A.hmm not found. Run install_isoform_databases() first. Skipping."
          )
          .log_predictor_status(
            "Pfam-hmmscan",
            "SKIPPED",
            "Pfam-A.hmm not found"
          )
        } else {
          pfam_tbl_l <- file.path(tool_save_dir("Pfam"), "pfam_domtblout.txt")
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
            hm_pfam_l <- file.path(
              tool_save_dir("Pfam"),
              "hmmscan_as_pfam_scan.txt"
            )
            hm_converted <- .pfam_from_hmmscan_domtblout(pfam_tbl_l, hm_pfam_l)

            if (is.null(hm_converted)) {
              message(
                "  No Pfam hits could be parsed from hmmscan's domtblout output."
              )
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
                .log_predictor_status(
                  "Pfam-hmmscan",
                  "OK",
                  paste("result file:", hm_converted)
                )
              } else {
                message(
                  "  hmmscan import failed; Pfam domain annotation NOT added to switch_list."
                )
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
              sprintf(
                "exit code %s, result file exists: %s",
                hm_status,
                file.exists(pfam_tbl_l)
              )
            )
          }
        }
      }
    }
  }

  n_orf <- if (!is.null(switch_list$orfAnalysis)) {
    nrow(switch_list$orfAnalysis)
  } else {
    0L
  }

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

  n_cds <- if (
    !is.null(switch_list$isoformFeatures) &&
      "codingPotential" %in% colnames(switch_list$isoformFeatures)
  ) {
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
    if (n_domain > 0 || n_sigp > 0 || n_cds > 0 || n_topology > 0) {
      "OK"
    } else {
      "EMPTY"
    },
    summary_msg
  )

  message("\n========== PREDICTOR RUN SUMMARY ==========")
  message("  Log directory: ", log_dir)
  message("  Per-tool logs:")

  all_logs <- list.files(
    log_dir,
    pattern = "\\.log$",
    full.names = TRUE,
    recursive = TRUE
  )
  log_rel_paths <- sub(
    paste0("^", gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", log_dir), "[/\\\\]?"),
    "",
    all_logs
  )

  for (i in seq_along(all_logs)) {
    sz <- file.info(all_logs[i])$size
    message(
      "     ",
      log_rel_paths[i],
      " (",
      format(round(sz / 1024, 1)),
      " KB)"
    )
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
    tool_logs = log_rel_paths
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

# mod_isoform_switch.R -- IsoformSwitchAnalyzeR orchestration (run_isoform_switch)
#

#' Run IsoformSwitchAnalyzeR analysis
#' @param custom_transcript_id_map A data frame or path to a TSV/CSV with
#'   `count_id` and `fasta_id` columns. Count IDs map to FASTA identifiers, and
#'   count IDs must correspond to transcript IDs in the imported annotation.
#' @export
run_isoform_switch <- function(
  dte_results = NULL,
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
  plot_topology = TRUE
) {
  test_engine <- match.arg(test_engine)

  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    stop(
      "Please install IsoformSwitchAnalyzeR: BiocManager::install('IsoformSwitchAnalyzeR')"
    )
  }

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  if (!is.null(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  switch_checkpoint_signature <- .object_md5(list(
    counts = if (identical(isoform_obj$type, "tximport")) {
      isoform_obj$txi$counts
    } else {
      isoform_obj$counts
    },
    metadata = isoform_obj$meta,
    tx2gene = isoform_obj$tx2gene,
    precomputed_results = .switch_results_fingerprint(dte_results, dtu_results),
    condition = condition,
    level = level,
    base = base,
    fasta = .file_fingerprint(fasta_file),
    gff = .file_fingerprint(gff_file),
    custom_transcript_id_map = if (is.data.frame(custom_transcript_id_map)) {
      custom_transcript_id_map
    } else {
      .file_fingerprint(custom_transcript_id_map)
    },
    test_engine = test_engine,
    skip_fasta_filter = skip_fasta_filter,
    plot_topology = plot_topology,
    run_predictors = run_predictors,
    predictor_cpu = predictor_cpu,
    use_wsl = use_wsl,
    wsl_distro = wsl_distro,
    organism = organism
  ))

  if (is.null(log_dir)) {
    log_dir <- file.path(out_dir, "Log")
  }
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

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

    if (
      isTRUE(debug_info$pfam_db_found) && !isTRUE(debug_info$pfam_db_indexed)
    ) {
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

  .ckpt_path <- function(nm) {
    if (!is.null(save_dir)) file.path(save_dir, nm) else ""
  }
  .ckpt_exists <- function(nm) {
    p <- .ckpt_path(nm)
    nzchar(p) && file.exists(p)
  }

  .ckpt_save <- function(obj, nm) {
    p <- .ckpt_path(nm)
    if (nzchar(p)) {
      .checkpoint_save(obj, p, switch_checkpoint_signature)
      message("Checkpoint saved -> ", nm)
    }
    invisible(obj)
  }

  .ckpt_load <- function(nm) {
    p <- .ckpt_path(nm)
    message("Resuming from checkpoint: ", nm)
    .checkpoint_load(p, switch_checkpoint_signature)
  }

  switch_list <- NULL
  already_analyzed <- FALSE

  if (
    !is.null(resume_from) &&
      file.exists(file.path(resume_from, "switch_list.rds"))
  ) {
    signature_path <- file.path(resume_from, "switch_list.signature.rds")
    if (
      file.exists(signature_path) &&
        identical(readRDS(signature_path), switch_checkpoint_signature)
    ) {
      message("Resuming from compatible saved SwitchList: ", resume_from)
      switch_list <- readRDS(file.path(resume_from, "switch_list.rds"))
      if (!is.null(switch_list$isoformSwitchAnalysis)) {
        message("Full analysis already present - skipping combined analysis.")
        already_analyzed <- TRUE
      }
    } else {
      message(
        "Ignoring incompatible or legacy saved SwitchList in: ",
        resume_from
      )
    }
  }

  if (is.null(switch_list)) {
    if (.ckpt_exists("step2_analyzed.rds")) {
      switch_list <- .ckpt_load("step2_analyzed.rds")
      already_analyzed <- !is.null(switch_list)
    }
    if (is.null(switch_list) && .ckpt_exists("step1_imported.rds")) {
      switch_list <- .ckpt_load("step1_imported.rds")
    }
  }

  if (is.null(switch_list)) {
    message("Building SwitchList from raw data...")

    if (isoform_obj$type == "tximport") {
      count_matrix <- isoform_obj$txi$counts
    } else {
      count_matrix <- isoform_obj$counts
    }

    sample_col <- if ("Sample" %in% colnames(isoform_obj$meta)) {
      "Sample"
    } else {
      "sample_id"
    }

    meta_sub <- isoform_obj$meta[
      isoform_obj$meta[[condition]] %in% c(base, level),
      ,
      drop = FALSE
    ]

    if (nrow(meta_sub) == 0) {
      stop("No samples found for the comparison groups: ", base, " and ", level)
    }

    keep_samples <- intersect(rownames(meta_sub), colnames(count_matrix))

    if (length(keep_samples) == 0) {
      stop(
        "No matching sample IDs between metadata and count matrix after filtering."
      )
    }

    meta_sub <- meta_sub[keep_samples, , drop = FALSE]
    count_matrix <- count_matrix[, keep_samples, drop = FALSE]

    sample_vector <- meta_sub[[sample_col]]
    colnames(count_matrix) <- sample_vector

    group_counts <- table(meta_sub[[condition]])

    if (any(group_counts < 2)) {
      warning(
        "One or both groups have fewer than 2 samples. DTU/switch analysis may be unreliable."
      )
    }

    if (!is.null(custom_transcript_id_map)) {
      if (
        is.character(custom_transcript_id_map) &&
          file.exists(custom_transcript_id_map)
      ) {
        map_df <- data.table::fread(
          custom_transcript_id_map,
          header = TRUE,
          data.table = FALSE
        )
      } else if (is.data.frame(custom_transcript_id_map)) {
        map_df <- custom_transcript_id_map
      } else {
        stop(
          "custom_transcript_id_map must be a data.frame or a path to a TSV/CSV file."
        )
      }

      required_cols <- c("count_id", "fasta_id")

      if (!all(required_cols %in% colnames(map_df))) {
        stop(
          "custom_transcript_id_map must contain columns: ",
          paste(required_cols, collapse = ", ")
        )
      }
      map_df$count_id <- as.character(map_df$count_id)
      map_df$fasta_id <- as.character(map_df$fasta_id)
      if (
        anyNA(map_df$count_id) ||
          anyNA(map_df$fasta_id) ||
          any(!nzchar(trimws(map_df$count_id))) ||
          any(!nzchar(trimws(map_df$fasta_id)))
      ) {
        stop(
          "custom_transcript_id_map must not contain missing or empty IDs.",
          call. = FALSE
        )
      }
      clean_count_ids <- clean_transcript_id(map_df$count_id)
      if (anyDuplicated(clean_count_ids)) {
        stop(
          "custom_transcript_id_map has count IDs that collide after ",
          "transcript-ID normalization.",
          call. = FALSE
        )
      }
      clean_fasta_ids <- clean_transcript_id(map_df$fasta_id)
      if (anyDuplicated(clean_fasta_ids)) {
        stop(
          "custom_transcript_id_map has FASTA IDs that collide after ",
          "transcript-ID normalization.",
          call. = FALSE
        )
      }

      old_rownames <- rownames(count_matrix)
      new_rownames <- clean_fasta_ids[match(
        clean_transcript_id(old_rownames),
        clean_count_ids
      )]
      mapped <- !is.na(new_rownames) & nzchar(new_rownames)

      if (any(!mapped)) {
        warning(
          "Some count IDs were not found in the mapping file; they will be kept unchanged."
        )
        new_rownames[!mapped] <- clean_transcript_id(old_rownames[!mapped])
      }

      rownames(count_matrix) <- new_rownames
      message(
        "  Remapped ",
        sum(mapped),
        " transcript IDs using custom map."
      )
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
        "  Kept ",
        nrow(count_matrix),
        " / ",
        length(clean_rownames),
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
    if (!is.null(custom_transcript_id_map)) {
      tx2gene$tx_clean <- .remap_switch_tx2gene(tx2gene$tx_id, map_df)
    }

    keep_tx <- rownames(count_matrix) %in% tx2gene$tx_clean
    count_matrix <- count_matrix[keep_tx, , drop = FALSE]
    tx2gene <- tx2gene[match(rownames(count_matrix), tx2gene$tx_clean), ]

    gene_counts <- table(tx2gene$gene_id)
    genes_multi <- names(gene_counts[gene_counts > 1])
    keep_genes <- tx2gene$gene_id %in% genes_multi

    count_matrix <- count_matrix[keep_genes, , drop = FALSE]
    tx2gene <- tx2gene[keep_genes, ]

    message(
      "  Kept ",
      nrow(count_matrix),
      " transcripts from ",
      length(unique(tx2gene$gene_id)),
      " genes after removing single-transcript genes."
    )

    isoform_count_matrix <- round(count_matrix)

    design_matrix <- data.frame(
      sampleID = sample_vector,
      condition = factor(meta_sub[[condition]], levels = c(base, level)),
      stringsAsFactors = FALSE
    )

    rownames(design_matrix) <- design_matrix$sampleID

    if (!"package:dplyr" %in% search()) {
      attachNamespace("dplyr")
    }

    .count_fasta_entries <- function(path) {
      tryCatch(length(Biostrings::fasta.seqlengths(path)), error = function(e) {
        NA_integer_
      })
    }

    .count_gtf_transcripts <- function(path) {
      tryCatch(
        {
          con <- file(path, open = "rt")
          on.exit(close(con), add = TRUE)

          n <- 0L

          repeat {
            lines <- readLines(con, n = 200000L, warn = FALSE)
            if (length(lines) == 0L) {
              break
            }
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
          "transcript sets (",
          n_fasta,
          " fasta sequences vs. ",
          n_gtf,
          " GTF transcripts; ",
          "count matrix after filtering has ",
          n_count,
          " transcripts). importRdata() will ",
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
    if (!dir.exists(clean_ref_dir)) {
      dir.create(clean_ref_dir, recursive = TRUE)
    }

    fasta_file_clean <- file.path(
      clean_ref_dir,
      "isoform_fasta_ids_cleaned.fasta"
    )
    gff_file_clean <- file.path(clean_ref_dir, "isoform_gff_ids_cleaned.gtf")

    message(
      "  Pre-cleaning transcript IDs in FASTA/GTF so importRdata() gets exact matches..."
    )

    seqs <- Biostrings::readDNAStringSet(fasta_file)
    names(seqs) <- clean_id(names(seqs))

    dup_fa <- unique(names(seqs)[duplicated(names(seqs))])

    if (length(dup_fa) > 0) {
      stop(
        "Cleaning transcript IDs in isoform_fasta produced ",
        length(dup_fa),
        " duplicate ID(s) -- e.g. ",
        paste(utils::head(dup_fa, 5), collapse = ", "),
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
        if (length(lines) == 0L) {
          break
        }

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
            tx_ids_seen <- c(
              tx_ids_seen,
              sub(id_re, "\\1", regmatches(lines[tx_row], m2))
            )
          }
        }

        writeLines(lines, con_out)
      }

      if (!any_id_seen) {
        stop(
          "No `transcript_id \"...\"` attributes found in ",
          path_in,
          " -- expected standard GTF2-style attributes."
        )
      }

      dup <- unique(tx_ids_seen[duplicated(tx_ids_seen)])

      if (length(dup) > 0) {
        stop(
          "Cleaning transcript IDs in isoform_gff produced ",
          length(dup),
          " duplicate transcript ID(s) among `transcript` feature rows -- e.g. ",
          paste(utils::head(dup, 5), collapse = ", "),
          ". The GTF contains near-duplicate transcript records that only differ in the ",
          "part clean_id() strips."
        )
      }

      invisible(path_out)
    }

    gtf_id_cleaner <- if (is.null(custom_transcript_id_map)) clean_id else
      function(ids) .remap_switch_tx2gene(ids, map_df)
    .rewrite_gtf_transcript_ids(gff_file, gff_file_clean, gtf_id_cleaner)

    importrdata_extra_args <- list(estimateDifferentialGeneRange = FALSE)

    if (
      "detectAndCorrectUnwantedEffects" %in%
        names(formals(IsoformSwitchAnalyzeR::importRdata))
    ) {
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

    if (
      !is.null(switch_list$ntSequence) && length(switch_list$ntSequence) > 0
    ) {
      genome_object <- NULL
    } else if (
      !is.null(bsgenome_name) && requireNamespace(bsgenome_name, quietly = TRUE)
    ) {
      genome_object <- getExportedValue(bsgenome_name, bsgenome_name)
      message("Using BSgenome: ", bsgenome_name)
    } else if (
      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
    ) {
      genome_object <- getExportedValue(
        "BSgenome.Hsapiens.UCSC.hg38",
        "BSgenome.Hsapiens.UCSC.hg38"
      )
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
            no_switches <- grepl(
              "no genes were considered switching",
              conditionMessage(e),
              ignore.case = TRUE
            )

            if (no_switches) {
              message("\n=== No isoform switches detected ===")
              message(
                "isoformSwitchAnalysisCombined() found zero genes meeting the switching ",
                "cutoffs for '",
                level,
                "' vs '",
                base,
                "'. Isoform switch analysis will be skipped, ",
                "but DTE and DTU results are unaffected."
              )
            } else {
              message("\n=== isoformSwitchAnalysisCombined() failed ===")
              message("Error: ", conditionMessage(e))
              message(
                "Isoform switch analysis will be skipped; DTE and DTU results are unaffected."
              )
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
          stop(
            "DRIMSeq results must contain columns: ",
            paste(required_cols, collapse = ", ")
          )
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
            message(
              "Error in isoformSwitchAnalysisCombined with DRIMSeq results: ",
              e$message
            )
            NULL
          }
        )
      },
      satuRn = {
        message("Running satuRn DTU test...")

        if (!requireNamespace("satuRn", quietly = TRUE)) {
          stop(
            "Package 'satuRn' is required for test_engine = 'satuRn'. Install with BiocManager::install('satuRn')."
          )
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
          message(
            "satuRn test did not return a switchList. Skipping further analysis."
          )
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
          "  Could not annotate alternative splicing events: ",
          e$message,
          "\n  Splicing-type summary/enrichment plots will be skipped; all other results are unaffected."
        )
        switch_list
      }
    )

    .ckpt_save(switch_list, "step2_analyzed.rds")
  }

  if (isTRUE(run_predictors)) {
    predictors_loaded <- FALSE
    if (.ckpt_exists("step3_predictors.rds")) {
      message("Loading external predictor results from step-3 checkpoint.")
      cached_switch_list <- .ckpt_load("step3_predictors.rds")
      if (!is.null(cached_switch_list)) {
        switch_list <- cached_switch_list
        predictors_loaded <- TRUE
      }
    }
    if (!predictors_loaded) {
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
    refresh_loaded <- FALSE
    if (.ckpt_exists("step3_5_refreshed.rds")) {
      message("Loading refreshed switch consequences from step-3.5 checkpoint.")
      cached_switch_list <- .ckpt_load("step3_5_refreshed.rds")
      if (!is.null(cached_switch_list)) {
        switch_list <- cached_switch_list
        refresh_loaded <- TRUE
      }
    }
    if (!refresh_loaded) {
      message(
        "Refreshing switch consequence analysis and plots with predictor annotations..."
      )

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
      plot_refresh_dir <- file.path(
        out_dir,
        "DTU_DTE_report",
        "plots",
        "by_consequence"
      )

      switch_list <- tryCatch(
        {
          feat_cols <- colnames(switch_list$isoformFeatures)
          consequences_available <- c(
            "intron_retention",
            "ORF_seq_similarity",
            "NMD_status"
          )

          if (
            any(grepl(
              "coding_potential|coding_prob",
              feat_cols,
              ignore.case = TRUE
            ))
          ) {
            consequences_available <- c(
              consequences_available,
              "coding_potential"
            )
          }

          if (any(grepl("signal_peptide", feat_cols, ignore.case = TRUE))) {
            consequences_available <- c(
              consequences_available,
              "signal_peptide_identified"
            )
          }

          if (any(grepl("^domain", feat_cols, ignore.case = TRUE))) {
            consequences_available <- c(
              consequences_available,
              "domains_identified",
              "domain_isotype"
            )
          }

          if (!is.null(switch_list$topologyAnalysis)) {
            consequences_available <- c(
              consequences_available,
              "isoform_topology"
            )
          }

          sl <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
            switch_list,
            consequencesToAnalyze = consequences_available,
            alpha = 0.05,
            dIFcutoff = 0.1
          )

          if (!dir.exists(plot_refresh_dir)) {
            dir.create(plot_refresh_dir, recursive = TRUE)
          }

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

          if (isTRUE(plot_topology)) {
            run_switch_plots()
          } else {
            suppressMessages(run_switch_plots())
          }

          message(
            "  Refreshed switch plots (now including predictor annotations) saved to: ",
            plot_refresh_dir
          )

          sl
        },
        error = function(e) {
          message(
            "  Could not refresh switch consequences/plots with predictor data: ",
            e$message,
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
    saveRDS(
      switch_checkpoint_signature,
      file.path(save_dir, "switch_list.signature.rds")
    )
    message("Saved final SwitchList to ", save_dir)
  }

  message("Isoform switch analysis completed. Results saved in: ", out_dir)

  return(switch_list)
}
