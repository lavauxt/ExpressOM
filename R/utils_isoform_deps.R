# utils_isoform_deps.R -- locating CPAT/Pfam data files and diagnosing the WSL/conda
# execution environment (debug_wsl() etc.)

#' Find the Pfam-A.hmm database in common paths inside the execution environment
#' @keywords internal
.find_pfam_db <- function(wsl_distro = "Ubuntu-22.04",
                          use_wsl = TRUE,
                          conda_sh = NULL,
                          conda_env = "isoform_tools",
                          log_dir = NULL) {
  result <- .wsl_exec_script(
    bash_body = c(
      'if [ -n "$PFAM_DB" ] && [ -f "$PFAM_DB" ]; then echo "$PFAM_DB"; exit 0; fi',
      'for _p in "$HOME/pfam_db/Pfam-A.hmm" "$HOME/databases/Pfam-A.hmm" "$HOME/.local/share/pfam/Pfam-A.hmm" "/opt/databases/Pfam-A.hmm" "/usr/local/share/pfam/Pfam-A.hmm"; do',
      '  if [ -f "$_p" ]; then echo "$_p"; break; fi',
      'done'
    ),
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    conda_sh = conda_sh,
    conda_env = conda_env,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = log_dir
  )

  result <- trimws(result[nzchar(trimws(result))])

  if (length(result) > 0) result[1] else NULL
}

#' Locate a CPAT data file (hexamer table or logit model) for a given organism
#'
#' CPAT's prebuilt hexamer/logit files are not bundled with the
#' IsoformSwitchAnalyzeR R package (that package only ships a small
#' `cpat_results.txt` example, not the actual model files); they are
#' downloaded by `install_isoform_databases()` into `$HOME/.cpat_data` (or
#' `$CPAT_DATA`). Both `.find_cpat_hexamer()` and `.find_cpat_logit_model()`
#' search the same set of locations.
#' @keywords internal
.find_cpat_data_file <- function(fname,
                                 wsl_distro = "Ubuntu-22.04",
                                 use_wsl = TRUE,
                                 conda_sh = NULL,
                                 conda_env = "isoform_tools",
                                 log_dir = NULL) {
  isa_local <- system.file("extdata", fname, package = "IsoformSwitchAnalyzeR")

  if (nzchar(isa_local) && file.exists(isa_local)) {
    if (.Platform$OS.type == "windows" && use_wsl) {
      return(.to_wsl_path(isa_local, wsl_distro))
    }
    return(isa_local)
  }

  cpat_data_env <- Sys.getenv("CPAT_DATA", "")

  if (nzchar(cpat_data_env)) {
    candidate <- file.path(cpat_data_env, fname)

    if (file.exists(candidate)) {
      if (.Platform$OS.type == "windows" && use_wsl) {
        return(.to_wsl_path(candidate, wsl_distro))
      }
      return(candidate)
    }
  }

  result <- .wsl_exec_script(
    bash_body = sprintf(
      'for _p in "$HOME/.cpat_data/%s" "${CPAT_DATA:-.}/%s"; do if [ -f "$_p" ]; then echo "$_p"; break; fi; done',
      fname,
      fname
    ),
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    conda_sh = conda_sh,
    conda_env = conda_env,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = log_dir
  )

  result <- trimws(result[nzchar(trimws(result))])

  if (length(result) > 0) result[1] else NULL
}

#' Locate the CPAT logit model RData file for a given organism
#' @keywords internal
.find_cpat_logit_model <- function(organism = "Human",
                                   wsl_distro = "Ubuntu-22.04",
                                   use_wsl = TRUE,
                                   conda_sh = NULL,
                                   conda_env = "isoform_tools",
                                   log_dir = NULL) {
  .find_cpat_data_file(
    paste0(organism, "_logitModel.RData"),
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    conda_sh = conda_sh,
    conda_env = conda_env,
    log_dir = log_dir
  )
}

#' Locate the CPAT hexamer frequency table for a given organism
#' @keywords internal
.find_cpat_hexamer <- function(organism = "Human",
                               wsl_distro = "Ubuntu-22.04",
                               use_wsl = TRUE,
                               conda_sh = NULL,
                               conda_env = "isoform_tools",
                               log_dir = NULL) {
  .find_cpat_data_file(
    paste0(organism, "_Hexamer.tsv"),
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    conda_sh = conda_sh,
    conda_env = conda_env,
    log_dir = log_dir
  )
}

#' Check if WSL is available and configured
#' @param distro WSL distribution name
#' @return Logical indicating success
#' @export
check_wsl <- function(distro = "Ubuntu-22.04") {
  if (.Platform$OS.type != "windows") {
    return(FALSE)
  }

  wsl_bin <- unname(Sys.which("wsl"))
  if (!nzchar(wsl_bin)) {
    return(FALSE)
  }

  res <- tryCatch(
    suppressWarnings(
      system2(
        "wsl",
        c("-d", distro, "--", "echo", "OK"),
        stdout = TRUE,
        stderr = FALSE
      )
    ),
    error = function(e) character(0)
  )

  length(res) > 0 && trimws(res[[1]]) == "OK"
}

#' Debug the external-predictor execution environment
#' @export
debug_wsl <- function(distro = "Ubuntu-22.04",
                      out_dir = NULL,
                      log_dir = NULL,
                      conda_env = "isoform_tools",
                      verbose = TRUE,
                      use_wsl = NULL) {

  is_windows <- .Platform$OS.type == "windows"
  if (is.null(use_wsl)) use_wsl <- is_windows

  via_wsl <- is_windows && isTRUE(use_wsl)

  if (verbose) {
    where <- if (via_wsl) paste0("WSL (distro: ", distro, ")") else "the native execution environment"
    message("Checking external predictor tools in ", where, "...")
  }

  results <- list(
    timestamp = Sys.time(),
    platform = if (via_wsl) "windows+wsl" else .Platform$OS.type,
    distro = if (via_wsl) distro else NA_character_,
    wsl_available = FALSE,
    conda_available = FALSE,
    conda_env_exists = FALSE,
    tools = list(),
    errors = character()
  )

  effective_log_dir <- if (!is.null(log_dir)) {
    log_dir
  } else if (!is.null(out_dir)) {
    file.path(out_dir, "Log")
  } else {
    NULL
  }

  if (!is.null(effective_log_dir) && !dir.exists(effective_log_dir)) {
    dir.create(effective_log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  .finish <- function(results) {
    if (!is.null(effective_log_dir)) {
      log_file <- file.path(effective_log_dir, "wsl_debug.json")
      jsonlite::write_json(results, log_file, pretty = TRUE, auto_unbox = TRUE)

      if (verbose) {
        message("  \u2192 Log written to ", log_file)
      }
    } else if (verbose) {
      message(
        "  (no out_dir/log_dir supplied to debug_wsl() -- wsl_debug.json not written; ",
        "pass one of those to get a persistent record of this check)"
      )
    }

    if (verbose && length(results$errors) > 0) {
      message("Predictor environment check found ", length(results$errors), " issue(s):")
      for (e in results$errors) message("  - ", e)
    }

    invisible(results)
  }

  if (via_wsl && !nzchar(Sys.which("wsl"))) {
    results$errors <- c(results$errors, "'wsl' executable not found on PATH")
    if (verbose) message("  \u2717 wsl executable not found")
    return(.finish(results))
  }

  probe <- .wsl_exec_script(
    bash_body = "echo OK",
    wsl_distro = distro,
    use_wsl = via_wsl,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = effective_log_dir
  )

  if (length(probe) == 0 || !any(grepl("OK", probe))) {
    results$errors <- c(
      results$errors,
      if (via_wsl) {
        "WSL not responding or distro not found"
      } else {
        "Could not execute a bash script natively (is 'bash' on PATH?)"
      }
    )

    if (verbose) message("  \u2717 execution environment not responding")

    return(.finish(results))
  }

  results$wsl_available <- TRUE
  if (verbose) message("  \u2713 execution environment reachable")

  conda_check <- .wsl_exec_script(
    bash_body = "command -v conda || echo 'not found'",
    wsl_distro = distro,
    use_wsl = via_wsl,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = effective_log_dir
  )

  if (any(grepl("conda", conda_check)) && !any(grepl("^not found$", trimws(conda_check)))) {
    results$conda_available <- TRUE
    if (verbose) message("  \u2713 conda found")

    env_check <- .wsl_exec_script(
      bash_body = sprintf(
        "conda env list 2>/dev/null | grep -q '%s' && echo 'exists'",
        conda_env
      ),
      wsl_distro = distro,
      use_wsl = via_wsl,
      intern = TRUE,
      ignore_stderr = TRUE,
      log_dir = effective_log_dir
    )

    if (length(env_check) > 0 && any(grepl("exists", env_check))) {
      results$conda_env_exists <- TRUE
      if (verbose) message("  \u2713 conda environment '", conda_env, "' exists")
    } else {
      results$errors <- c(results$errors, paste0("conda environment '", conda_env, "' missing"))
      if (verbose) message("  \u2717 conda environment '", conda_env, "' not found")
    }
  } else {
    results$errors <- c(results$errors, "conda not found")
    if (verbose) message("  \u2717 conda not found (tools may still be available directly on PATH)")
  }

  conda_sh_guess <- NULL

  if (results$conda_env_exists) {
    conda_sh_guess <- .find_conda_sh(
      wsl_distro = distro,
      use_wsl = via_wsl,
      log_dir = effective_log_dir
    )
  }

  tools_list <- c("cpat", "signalp6", "hmmscan", "interproscan.sh")

  .check_tool <- function(tool) {
    found <- .wsl_tool_exists(
      tool,
      wsl_distro = distro,
      use_wsl = via_wsl,
      conda_sh = NULL,
      conda_env = conda_env,
      log_dir = effective_log_dir
    )

    if (!found && !is.null(conda_sh_guess)) {
      found <- .wsl_tool_exists(
        tool,
        wsl_distro = distro,
        use_wsl = via_wsl,
        conda_sh = conda_sh_guess,
        conda_env = conda_env,
        log_dir = effective_log_dir
      )
    }

    found
  }

  for (tool in tools_list) {
    found <- .check_tool(tool)
    results$tools[[tool]] <- found

    if (verbose) {
      message("   ", if (found) "\u2713" else "\u2717", "  ", tool)
    }
  }

  if (isTRUE(results$tools[["cpat"]])) {
    results$tools[["run_cpat.py"]] <- NA

    if (verbose) {
      message("   -  run_cpat.py (skipped, not needed: 'cpat' already found)")
    }
  } else {
    found <- .check_tool("run_cpat.py")
    results$tools[["run_cpat.py"]] <- found

    if (verbose) {
      message("   ", if (found) "\u2713" else "\u2717", "  run_cpat.py")
    }
  }

  tool_alternatives <- list(
    "CPAT (cpat / run_cpat.py)" = c("cpat", "run_cpat.py"),
    "SignalP (signalp6)" = c("signalp6"),
    "hmmscan" = c("hmmscan")
  )

  for (label in names(tool_alternatives)) {
    alts <- tool_alternatives[[label]]

    if (!any(unlist(results$tools[alts]))) {
      results$errors <- c(results$errors, paste0("Missing tool: ", label))
    }
  }

  pfam_path <- .find_pfam_db(
    wsl_distro = distro,
    use_wsl = via_wsl,
    conda_sh = conda_sh_guess,
    conda_env = conda_env,
    log_dir = effective_log_dir
  )

  results$pfam_db_found <- !is.null(pfam_path)

  if (verbose) {
    message("  ", if (results$pfam_db_found) "\u2713" else "\u2717", " Pfam-A.hmm")
  }

  results$pfam_db_indexed <- FALSE

  if (results$pfam_db_found) {
    idx_check <- .wsl_exec_script(
      bash_body = sprintf(
        '[ -f %s ] && [ -f %s ] && [ -f %s ] && [ -f %s ] && echo INDEXED || echo MISSING',
        .dq(paste0(pfam_path, ".h3f")),
        .dq(paste0(pfam_path, ".h3i")),
        .dq(paste0(pfam_path, ".h3m")),
        .dq(paste0(pfam_path, ".h3p"))
      ),
      wsl_distro = distro,
      use_wsl = via_wsl,
      conda_sh = conda_sh_guess,
      conda_env = conda_env,
      intern = TRUE,
      ignore_stderr = TRUE,
      log_dir = effective_log_dir
    )

    results$pfam_db_indexed <- any(grepl("INDEXED", idx_check))

    if (verbose) {
      message(
        "   ",
        if (results$pfam_db_indexed) "\u2713" else "\u2717",
        " Pfam-A.hmm hmmpress-indexed",
        if (!results$pfam_db_indexed) {
          " (found the database but it is NOT indexed -- hmmscan will fail; run `hmmpress <path>` inside the isoform_tools env, or re-run install_isoform_databases())"
        } else {
          ""
        }
      )
    }

    if (!results$pfam_db_indexed) {
      results$errors <- c(
        results$errors,
        paste0("Pfam-A.hmm found at ", pfam_path, " but not hmmpress-indexed")
      )
    }
  }

  cpat_logit_human <- .find_cpat_logit_model(
    "Human",
    distro,
    via_wsl,
    conda_sh_guess,
    conda_env,
    log_dir = effective_log_dir
  )

  cpat_logit_mouse <- .find_cpat_logit_model(
    "Mouse",
    distro,
    via_wsl,
    conda_sh_guess,
    conda_env,
    log_dir = effective_log_dir
  )

  results$cpat_models_found <- !is.null(cpat_logit_human) && !is.null(cpat_logit_mouse)

  if (verbose) {
    message("   ", if (results$cpat_models_found) "\u2713" else "\u2717", " CPAT logit models")
  }

  .finish(results)
}

