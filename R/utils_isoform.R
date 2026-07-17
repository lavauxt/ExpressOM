# ==============================================================================
# utils_isoform.R - WSL/conda helpers + install/database utilities
# ==============================================================================
#
# Cross-platform notes
# ---------------------------------------------------------------------------
# External predictor tools (CPAT, SignalP, hmmscan/InterProScan) are launched
# through a single chokepoint, .wsl_exec_script(). On Windows with
# use_wsl = TRUE, commands are routed through `wsl -d <distro> bash <script>`.
# On every other platform (native Linux, macOS, or R already running inside
# WSL) commands are simply run with the local `bash`. Both branches are
# wrapped defensively so a missing `wsl` or `bash` executable produces a
# clean, informative skip instead of a hard R error that would abort the
# whole pipeline.
# ---------------------------------------------------------------------------

# Small in-session cache to avoid repeatedly shelling out to `find` for the
# same conda.sh / distro combination during a single predictor run.
# @keywords internal
.expressom_cache <- new.env(parent = emptyenv())

# ---- Shared persistent-environment-variable file -----------------------------
#
.ISOFORM_ENV_FILE <- "$HOME/.isoform_tools_env.sh"

#' Double-quote a shell argument while preserving `$VAR`-style expansion
#'
#' Unlike `shQuote(x, type = "sh")` (which single-quotes and therefore
#' freezes any literal `$HOME`/`$VAR` inside `x`), this wraps `x` in double
#' quotes and escapes only the characters that are special inside double
#' quotes (`\\`, `"`, `` ` ``), so a value like `"$HOME/pfam_db"` still
#' expands correctly when the generated script runs.
#' @keywords internal
.dq <- function(x) {
  sprintf('"%s"', gsub('([\\\\"`])', '\\\\\\1', x))
}

#' Persist an environment variable for later predictor runs
#'
#' Appends (or replaces) an `export VAR="value"` line in the dedicated
#' env file that `.wsl_exec_script()` sources on every invocation. This is
#' the supported way for install helpers to make a path available to CPAT /
#' SignalP / hmmscan at prediction time -- do not write to ~/.bashrc (see
#' note above).
#'
#' @param var        Environment variable name (e.g. "CPAT_DATA").
#' @param value      Value to assign (will be shell-quoted).
#' @param wsl_distro WSL distribution name.
#' @param use_wsl    Route through WSL? (only effective on Windows.)
#' @return Invisibly, TRUE/FALSE indicating whether the write succeeded.
#' @keywords internal
.wsl_write_env_var <- function(var, value, wsl_distro = "Ubuntu", use_wsl = TRUE) {

  export_line <- sprintf("export %s=%s", var, .dq(value))
  body <- c(
    sprintf('touch "%s"', .ISOFORM_ENV_FILE),
    sprintf('grep -v "^export %s=" "%s" > "%s.tmp" 2>/dev/null || true', var, .ISOFORM_ENV_FILE, .ISOFORM_ENV_FILE),
    sprintf('echo %s >> "%s.tmp"', shQuote(export_line, type = "sh"), .ISOFORM_ENV_FILE),
    sprintf('mv "%s.tmp" "%s"', .ISOFORM_ENV_FILE, .ISOFORM_ENV_FILE)
  )
  status <- .wsl_exec_script(body, wsl_distro = wsl_distro, use_wsl = use_wsl,
                              intern = FALSE, ignore_stderr = TRUE, log_dir = NULL)
  ok <- isTRUE(status == 0L)
  if (ok) message("  -> Persisted ", var, " for future predictor runs (", .ISOFORM_ENV_FILE, ")")
  else message("  ! Could not persist ", var, " to ", .ISOFORM_ENV_FILE,
               " (non-fatal; export it manually inside the execution environment if predictors can't find it)")
  invisible(ok)
}

#' Best-effort discovery of a conda profile script (conda.sh)
#'
#' Searches the usual conda/mamba/miniforge install locations inside the
#' execution environment (WSL distro or native shell). Shared by debug_wsl()
#' and install_isoform_databases() so both use one, tested search order.
#'
#' @param wsl_distro WSL distribution name.
#' @param use_wsl    Route through WSL? (only effective on Windows.)
#' @return Path to conda.sh inside the execution environment, or NULL.
#' @keywords internal
.find_conda_sh <- function(wsl_distro = "Ubuntu", use_wsl = TRUE) {
  csh <- .wsl_exec_script(
    paste(
      'for c in "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh"',
      '"$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh"',
      '"$HOME/mambaforge/etc/profile.d/conda.sh" "$HOME/miniforge3/etc/profile.d/conda.sh";',
      'do [ -f "$c" ] && echo "$c" && break; done'
    ),
    wsl_distro = wsl_distro, use_wsl = use_wsl, intern = TRUE, ignore_stderr = TRUE,
    log_dir = NULL
  )
  csh <- trimws(csh[nzchar(trimws(csh))])
  if (length(csh) > 0) csh[1] else NULL
}

# ---- Internal WSL + conda execution helpers ----------------------------------

#' Convert a Windows file path to a WSL-compatible Unix path
#'
#' Uses wslpath -u when available, otherwise performs a manual
#' C:/foo -> /mnt/c/foo substitution. On non-Windows platforms
#' the path is returned unchanged.
#' @param win_path A Windows-style path string.
#' @param distro   WSL distribution name (default "Ubuntu").
#' @return A Unix-style path string suitable for use inside WSL.
#' @keywords internal
.to_wsl_path <- function(win_path, distro = "Ubuntu") {
  if (.Platform$OS.type != "windows") return(win_path)
  if (is.null(win_path) || !nzchar(trimws(as.character(win_path)))) return(win_path)

  p <- normalizePath(win_path, winslash = "/", mustWork = FALSE)

  r <- tryCatch(
    suppressWarnings(
      system2("wsl", c("-d", distro, "wslpath", "-u", p),
              stdout = TRUE, stderr = FALSE, timeout = 15)
    ),
    error = function(e) character(0)
  )
  r <- trimws(r[nzchar(trimws(r))])
  if (length(r) > 0) return(r[1])

  # Manual fallback: C:/foo/bar -> /mnt/c/foo/bar
  if (grepl("^[A-Za-z]:/", p))
    return(paste0("/mnt/", tolower(substr(p, 1, 1)), substring(p, 3)))
  p
}

#' Write a WSL command execution to log (for debugging)
#'
#' @param cmd The command string executed
#' @param exit_code Integer exit code
#' @param stdout Character vector of stdout (optional)
#' @param stderr Character vector of stderr (optional)
#' @param log_dir Directory where log files are stored
#' @keywords internal
.log_wsl_command <- function(cmd, exit_code, stdout = NULL, stderr = NULL, log_dir) {
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir, "wsl_commands.log")
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- c(
    paste0("[", timestamp, "]"),
    paste0("  CMD: ", cmd),
    paste0("  EXIT: ", exit_code)
  )
  if (!is.null(stdout) && length(stdout) > 0)
    entry <- c(entry, paste0("  STDOUT: ", paste(stdout, collapse = "\n    ")))
  if (!is.null(stderr) && length(stderr) > 0)
    entry <- c(entry, paste0("  STDERR: ", paste(stderr, collapse = "\n    ")))
  write(entry, file = log_file, append = TRUE)
}

#' Write bash commands to a temp script and execute via WSL or natively
#'
#' Avoids shell-escaping issues by writing commands to a temp .sh file and
#' running it with `wsl -d <distro> bash <script>` (Windows) or
#' `bash <script>` (Linux/macOS). Conda activation lines are prepended
#' automatically when conda_sh is provided.
#'
#' @param bash_body    Character vector of bash commands.
#' @param wsl_distro   WSL distribution name.
#' @param use_wsl      Route through WSL? (only effective on Windows.)
#' @param conda_sh     Full path to conda profile script inside execution env,
#'                     or NULL to skip conda activation.
#' @param conda_env    Conda environment to activate.
#' @param intern       If TRUE return stdout as character vector.
#' @param ignore_stderr Suppress stderr in R console when TRUE.
#' @param log_dir      Optional directory to write command log (if NULL, no log)
#' @return Integer exit code (0 = success) or character vector when intern=TRUE.
#' @keywords internal
.wsl_exec_script <- function(bash_body,
                              wsl_distro    = "Ubuntu",
                              use_wsl       = TRUE,
                              conda_sh      = NULL,
                              conda_env     = "isoform_tools",
                              intern        = FALSE,
                              ignore_stderr = TRUE,
                              log_dir       = NULL) {

  env_file_source <- sprintf('[ -f %s ] && . %s 2>/dev/null || true',
                             .ISOFORM_ENV_FILE, .ISOFORM_ENV_FILE)

  activate <- if (!is.null(conda_sh) && nzchar(conda_sh)) {
    c(env_file_source,
      sprintf('. "%s" 2>/dev/null || true', conda_sh),
      sprintf("conda activate %s 2>/dev/null || true", conda_env))
  } else env_file_source

  script <- c("#!/bin/bash", "set -e", activate, bash_body)
  tmp    <- tempfile(fileext = ".sh")
  on.exit(unlink(tmp), add = TRUE)

  con <- file(tmp, "wb")
  writeLines(script, con, sep = "\n")
  close(con)

  via_wsl <- (.Platform$OS.type == "windows") && isTRUE(use_wsl)

  run_status      <- 127L
  out             <- character(0)
  missing_exe_msg <- NULL

  run_capture <- function(cmd, args) {
    tryCatch({
      o <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
      list(out = o, status = attr(o, "status") %||% 0L)
    }, error = function(e) list(out = character(0), status = 127L,
                                 error = conditionMessage(e)))
  }

  if (via_wsl) {
    if (!nzchar(Sys.which("wsl"))) {
      missing_exe_msg <- paste0(
        "'wsl' executable not found on PATH. Install WSL (wsl --install) or ",
        "set use_wsl = FALSE if the required tools are available natively."
      )
    } else {
      wsl_sh <- .to_wsl_path(tmp, wsl_distro)
      args   <- c("-d", wsl_distro, "bash", wsl_sh)
      res <- run_capture("wsl", args)
      out <- res$out; run_status <- res$status
      if (!is.null(res$error)) missing_exe_msg <- res$error
      if (!ignore_stderr && length(out) > 0) message(paste(out, collapse = "\n"))
    }
  } else {
    bash_bin <- Sys.which("bash")
    if (!nzchar(bash_bin)) {
      missing_exe_msg <- paste0(
        "'bash' executable not found on PATH. External predictor tools ",
        "(CPAT / SignalP / Pfam) require a bash shell -- on native Windows ",
        "without WSL/Git-Bash these steps will be skipped."
      )
    } else {
      res <- run_capture(bash_bin, tmp)
      out <- res$out; run_status <- res$status
      if (!is.null(res$error)) missing_exe_msg <- res$error
      if (!ignore_stderr && length(out) > 0) message(paste(out, collapse = "\n"))
    }
  }

  if (!is.null(missing_exe_msg)) {
    message("WSL/predictor environment error: ", missing_exe_msg)
    warning(missing_exe_msg, call. = FALSE)
    run_status <- 127L
  }

  if (!is.null(log_dir)) {
    .log_wsl_command(paste(bash_body, collapse = "; "),
                     exit_code = run_status,
                     stdout = out,
                     stderr = missing_exe_msg,
                     log_dir = log_dir)
  }

  if (intern) return(out) else return(run_status)
}

#' Check whether a command-line tool is accessible in the execution environment
#' @keywords internal
.wsl_tool_exists <- function(tool_name,
                              wsl_distro = "Ubuntu",
                              use_wsl    = TRUE,
                              conda_sh   = NULL,
                              conda_env  = "isoform_tools") {
  status <- .wsl_exec_script(
    bash_body     = sprintf("command -v %s >/dev/null 2>&1", shQuote(tool_name, type = "sh")),
    wsl_distro    = wsl_distro, use_wsl = use_wsl,
    conda_sh      = conda_sh,  conda_env = conda_env,
    intern        = FALSE,     ignore_stderr = TRUE,
    log_dir       = NULL
  )
  isTRUE(status == 0L)
}

#' Find the Pfam-A.hmm database in common paths inside the execution environment
#' @keywords internal
.find_pfam_db <- function(wsl_distro = "Ubuntu",
                           use_wsl    = TRUE,
                           conda_sh   = NULL,
                           conda_env  = "isoform_tools") {
  result <- .wsl_exec_script(
    bash_body = c(
      'if [ -n "$PFAM_DB" ] && [ -f "$PFAM_DB" ]; then echo "$PFAM_DB"; exit 0; fi',
      'for _p in "$HOME/pfam_db/Pfam-A.hmm" "$HOME/databases/Pfam-A.hmm" "$HOME/.local/share/pfam/Pfam-A.hmm" "/opt/databases/Pfam-A.hmm" "/usr/local/share/pfam/Pfam-A.hmm"; do [ -f "$_p" ] && echo "$_p" && break; done'
    ),
    wsl_distro = wsl_distro, use_wsl = use_wsl,
    conda_sh   = conda_sh, conda_env = conda_env,
    intern = TRUE, ignore_stderr = TRUE,
    log_dir = NULL
  )
  result <- trimws(result[nzchar(trimws(result))])
  if (length(result) > 0) result[1] else NULL
}

#' Locate the CPAT logit model RData file for a given organism
#'
#' Search order: (1) IsoformSwitchAnalyzeR extdata, (2) CPAT_DATA env var,
#' (3) ~/.cpat_data inside the execution environment.
#' Returns a path valid inside the execution environment.
#' @keywords internal
.find_cpat_logit_model <- function(organism  = "Human",
                                    wsl_distro = "Ubuntu",
                                    use_wsl    = TRUE,
                                    conda_sh   = NULL,
                                    conda_env  = "isoform_tools") {
  fname <- paste0(organism, "_logitModel.RData")

  # 1. IsoformSwitchAnalyzeR ships the models in its extdata
  isa_local <- system.file("extdata", fname, package = "IsoformSwitchAnalyzeR")
  if (nzchar(isa_local) && file.exists(isa_local)) {
    if (.Platform$OS.type == "windows" && use_wsl)
      return(.to_wsl_path(isa_local, wsl_distro))
    return(isa_local)
  }

  # 2. CPAT_DATA environment variable (Windows R side)
  cpat_data_env <- Sys.getenv("CPAT_DATA", "")
  if (nzchar(cpat_data_env)) {
    candidate <- file.path(cpat_data_env, fname)
    if (file.exists(candidate)) {
      if (.Platform$OS.type == "windows" && use_wsl)
        return(.to_wsl_path(candidate, wsl_distro))
      return(candidate)
    }
  }

  # 3. Search inside execution environment (~/.cpat_data or $CPAT_DATA)
  result <- .wsl_exec_script(
    bash_body = sprintf(
      'for _p in "$HOME/.cpat_data/%s" "${CPAT_DATA:-.}/%s"; do [ -f "$_p" ] && echo "$_p" && break; done',
      fname, fname
    ),
    wsl_distro = wsl_distro, use_wsl = use_wsl,
    conda_sh   = conda_sh, conda_env = conda_env,
    intern = TRUE, ignore_stderr = TRUE,
    log_dir = NULL
  )
  result <- trimws(result[nzchar(trimws(result))])
  if (length(result) > 0) result[1] else NULL
}


#' Check if WSL is available and configured
#' @param distro WSL distribution name
#' @return Logical indicating success
#' @export
check_wsl <- function(distro = "Ubuntu") {
  if (.Platform$OS.type != "windows") return(FALSE)
  res <- system(paste('wsl -d', distro, '-- echo "OK"'), intern = TRUE, ignore.stderr = TRUE)
  return(length(res) > 0 && trimws(res[[1]]) == "OK")
}

#' Debug the external-predictor execution environment (CPAT / SignalP / Pfam)
#'
#' Checks that the shell used to run external predictor tools is reachable,
#' whether conda is available, whether the `isoform_tools` conda environment
#' exists, whether each required tool is on PATH, and whether the Pfam / CPAT
#' databases can be located. Writes a detailed log file to
#' `out_dir/Log/wsl_debug.json`.
#'
#' This function is fully cross-platform:
#' \itemize{
#'   \item On Windows with \code{use_wsl = TRUE} (the default on Windows),
#'     every check is routed through \code{wsl -d <distro> bash ...}.
#'   \item On native Linux, macOS, or when R itself is already running inside
#'     WSL, every check runs against the local \code{bash} / PATH directly --
#'     no WSL distribution is required or contacted.
#' }
#'
#' @param distro WSL distribution name (only used when routing through WSL)
#' @param out_dir Output directory for logs (if NULL, returns only the list).
#'   Logs are written to \code{out_dir/Log/wsl_debug.json} unless
#'   \code{log_dir} is given, in which case that path is used verbatim.
#' @param log_dir Optional explicit log directory, overriding the
#'   \code{out_dir/Log} default. Lets callers that already manage their own
#'   unified log tree (e.g. the full \code{expressom()} pipeline) point this
#'   at a single canonical location instead of scattering a second copy.
#' @param conda_env Name of the conda environment to check (default "isoform_tools")
#' @param verbose Print progress messages
#' @param use_wsl Logical: route checks through WSL. Defaults to `TRUE` on
#'   Windows and is ignored (treated as `FALSE`) on every other platform.
#' @return Invisibly, a list with check results
#' @export
debug_wsl <- function(distro = "Ubuntu", out_dir = NULL, log_dir = NULL,
                      conda_env = "isoform_tools", verbose = TRUE,
                      use_wsl = NULL) {

  is_windows <- .Platform$OS.type == "windows"
  if (is.null(use_wsl)) use_wsl <- is_windows
  via_wsl <- is_windows && isTRUE(use_wsl)

  if (verbose) {
    where <- if (via_wsl) paste0("WSL (distro: ", distro, ")") else "the native execution environment"
    message("Checking external predictor tools in ", where, "...")
  }

  results <- list(
    timestamp        = Sys.time(),
    platform         = if (via_wsl) "windows+wsl" else .Platform$OS.type,
    distro           = if (via_wsl) distro else NA_character_,
    wsl_available    = FALSE,
    conda_available  = FALSE,
    conda_env_exists = FALSE,
    tools            = list(),
    errors           = character()
  )

  # 1. Confirm the wsl executable exists (Windows+WSL path only), then confirm
  #    the execution shell (WSL distro, or local bash) actually responds.
  if (via_wsl && !nzchar(Sys.which("wsl"))) {
    results$errors <- c(results$errors, "'wsl' executable not found on PATH")
    if (verbose) message("  \u2717 wsl executable not found")
    return(invisible(results))
  }

  probe <- .wsl_exec_script("echo OK", wsl_distro = distro, use_wsl = via_wsl,
                            intern = TRUE, ignore_stderr = TRUE)
  if (length(probe) == 0 || !any(grepl("OK", probe))) {
    results$errors <- c(results$errors,
                        if (via_wsl) "WSL not responding or distro not found"
                        else "Could not execute a bash script natively (is 'bash' on PATH?)")
    if (verbose) message("  \u2717 execution environment not responding")
    return(invisible(results))
  }
  results$wsl_available <- TRUE
  if (verbose) message("  \u2713 execution environment reachable")

  # 2. Conda installation and environment
  conda_check <- .wsl_exec_script("command -v conda || echo 'not found'",
                                  wsl_distro = distro, use_wsl = via_wsl,
                                  intern = TRUE, ignore_stderr = TRUE)
  if (any(grepl("conda", conda_check)) && !any(grepl("^not found$", trimws(conda_check)))) {
    results$conda_available <- TRUE
    if (verbose) message("  \u2713 conda found")

    env_check <- .wsl_exec_script(
      sprintf("conda env list 2>/dev/null | grep -q '%s' && echo 'exists'", conda_env),
      wsl_distro = distro, use_wsl = via_wsl, intern = TRUE, ignore_stderr = TRUE
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
    conda_sh_guess <- .find_conda_sh(wsl_distro = distro, use_wsl = via_wsl)
  }

  # 3. Required tools (CPAT, SignalP, hmmscan, interproscan.sh)
  tools_list <- c("cpat", "run_cpat.py", "signalp6", "signalp", "hmmscan", "interproscan.sh")
  for (tool in tools_list) {
    found <- .wsl_tool_exists(tool, wsl_distro = distro, use_wsl = via_wsl,
                              conda_sh = NULL, conda_env = conda_env)
    if (!found && !is.null(conda_sh_guess)) {
      found <- .wsl_tool_exists(tool, wsl_distro = distro, use_wsl = via_wsl,
                                conda_sh = conda_sh_guess, conda_env = conda_env)
    }
    results$tools[[tool]] <- found
    if (verbose) message("  ", if (found) "\u2713" else "\u2717", " ", tool)
    if (!found && tool %in% c("cpat", "run_cpat.py", "signalp6", "signalp", "hmmscan")) {
      results$errors <- c(results$errors, paste0("Missing tool: ", tool))
    }
  }

  # 4. Additional databases (Pfam, CPAT models)
  pfam_path <- .find_pfam_db(wsl_distro = distro, use_wsl = via_wsl,
                             conda_sh = conda_sh_guess, conda_env = conda_env)
  results$pfam_db_found <- !is.null(pfam_path)
  if (verbose) message("  ", if (results$pfam_db_found) "\u2713" else "\u2717", " Pfam-A.hmm")

  results$pfam_db_indexed <- FALSE
  if (results$pfam_db_found) {
    idx_check <- .wsl_exec_script(
      sprintf('[ -f "%s.h3f" ] && [ -f "%s.h3i" ] && [ -f "%s.h3m" ] && [ -f "%s.h3p" ] && echo INDEXED || echo MISSING',
              pfam_path, pfam_path, pfam_path, pfam_path),
      wsl_distro = distro, use_wsl = via_wsl,
      conda_sh = conda_sh_guess, conda_env = conda_env,
      intern = TRUE, ignore_stderr = TRUE
    )
    results$pfam_db_indexed <- any(grepl("INDEXED", idx_check))
    if (verbose) message("  ", if (results$pfam_db_indexed) "\u2713" else "\u2717",
                         " Pfam-A.hmm hmmpress-indexed",
                         if (!results$pfam_db_indexed)
                           " (found the database but it is NOT indexed -- hmmscan will fail; run `hmmpress <path>` inside the isoform_tools env, or re-run install_isoform_databases())"
                         else "")
    if (!results$pfam_db_indexed)
      results$errors <- c(results$errors, paste0("Pfam-A.hmm found at ", pfam_path, " but not hmmpress-indexed"))
  }

  cpat_logit_human <- .find_cpat_logit_model("Human", distro, via_wsl, conda_sh_guess, conda_env)
  cpat_logit_mouse <- .find_cpat_logit_model("Mouse", distro, via_wsl, conda_sh_guess, conda_env)
  results$cpat_models_found <- !is.null(cpat_logit_human) && !is.null(cpat_logit_mouse)
  if (verbose) message("  ", if (results$cpat_models_found) "\u2713" else "\u2717", " CPAT logit models")

  effective_log_dir <- if (!is.null(log_dir)) log_dir
                        else if (!is.null(out_dir)) file.path(out_dir, "Log")
                        else NULL
  if (!is.null(effective_log_dir)) {
    if (!dir.exists(effective_log_dir)) dir.create(effective_log_dir, recursive = TRUE)
    log_file <- file.path(effective_log_dir, "wsl_debug.json")
    jsonlite::write_json(results, log_file, pretty = TRUE, auto_unbox = TRUE)
    if (verbose) message("  \u2192 Log written to ", log_file)
  }
  if (verbose && length(results$errors) > 0) {
    message("Predictor environment check found ", length(results$errors), " issue(s):")
    for (e in results$errors) message("  - ", e)
  }

  invisible(results)
}

#' Install SignalP from a Windows directory into WSL
#'
#' @param windows_signalp_dir Path to SignalP directory on Windows (e.g., "C:/signalp-5.0")
#' @param distro WSL distribution name (default "Ubuntu")
#' @param install_path Destination path inside WSL (default "/usr/local/signalp")
#' @return Logical indicating success
#' @export
install_signalp_from_windows <- function(windows_signalp_dir,
                                          distro       = "Ubuntu",
                                          install_path = "/usr/local/signalp") {
  if (!check_wsl(distro)) stop("WSL with distro ", distro, " not available.")
  if (!dir.exists(windows_signalp_dir))
    stop("Windows directory does not exist: ", windows_signalp_dir)

  wsl_windows_path <- tryCatch({
    system2("wsl", c("-d", distro, "wslpath", "-u", shQuote(windows_signalp_dir)), stdout = TRUE)
  }, error = function(e) {
    path <- gsub("\\\\", "/", windows_signalp_dir)
    if (grepl("^[A-Za-z]:", path)) {
      drive <- tolower(substr(path, 1, 1))
      path  <- sub("^[A-Za-z]:", paste0("/mnt/", drive), path)
    }
    path
  })
  wsl_windows_path <- trimws(wsl_windows_path[nzchar(trimws(wsl_windows_path))])[1]

  message("Copying SignalP from Windows (", wsl_windows_path, ") to WSL (", install_path, ")...")
  mkdir_status <- .wsl_exec_script(sprintf("sudo mkdir -p %s", .dq(install_path)),
                                   wsl_distro = distro, use_wsl = TRUE)
  if (!isTRUE(mkdir_status == 0L)) {
    message("Failed to create ", install_path, " (exit code ", mkdir_status, "). Check sudo permissions.")
    return(FALSE)
  }

  cp_status <- .wsl_exec_script(
    sprintf("sudo cp -r %s/. %s/", .dq(wsl_windows_path), .dq(install_path)),
    wsl_distro = distro, use_wsl = TRUE
  )
  if (!isTRUE(cp_status == 0L)) {
    message("Failed to copy SignalP files (exit code ", cp_status, "). Check permissions and path.")
    return(FALSE)
  }

  bin_path <- file.path(install_path, "bin", "signalp")
  chmod_status <- .wsl_exec_script(sprintf("sudo chmod +x %s", .dq(bin_path)),
                                   wsl_distro = distro, use_wsl = TRUE)
  ln_status <- .wsl_exec_script(
    sprintf("sudo ln -sf %s /usr/local/bin/signalp", .dq(bin_path)),
    wsl_distro = distro, use_wsl = TRUE
  )
  if (!isTRUE(chmod_status == 0L) || !isTRUE(ln_status == 0L)) {
    message("  ! Warning: chmod/symlink step reported a non-zero exit code (chmod=", chmod_status,
            ", ln=", ln_status, "). signalp may not be directly callable as `signalp`; check ", bin_path)
  }

  data_dir <- file.path(install_path, "data")
  data_dir_exists <- .wsl_exec_script(sprintf("[ -d %s ]", .dq(data_dir)),
                                      wsl_distro = distro, use_wsl = TRUE)
  if (isTRUE(data_dir_exists == 0L)) {
    .wsl_write_env_var("SIGNALP_DIR", data_dir, wsl_distro = distro, use_wsl = TRUE)
  } else {
    message("  ! No 'data' subdirectory found under ", install_path,
            " -- SIGNALP_DIR was not set. If this SignalP version stores models elsewhere, ",
            "set the appropriate env var manually via .wsl_write_env_var().")
  }

  message("SignalP installed to ", install_path, " and linked to /usr/local/bin/signalp")
  TRUE
}

#' Install required external tools inside WSL using mamba (default) or apt/pip
#'
#' @param distro WSL distribution name (default "Ubuntu")
#' @param use_mamba Logical: use mamba/conda (default TRUE)
#' @param install_databases Logical: also install CPAT, Pfam databases
#' @param windows_signalp_dir Optional Windows path to a SignalP installation
#' @param log_dir Optional directory to write a WSL command audit trail
#' @export
install_wsl_isoform_tools <- function(distro              = "Ubuntu",
                                       use_mamba           = TRUE,
                                       install_databases   = TRUE,
                                       windows_signalp_dir = NULL,
                                       log_dir              = NULL) {
  if (!check_wsl(distro)) stop("WSL with distro ", distro, " not available.")

  .run <- function(cmd, conda = FALSE) {
    body <- if (conda)
      sprintf('bash -c "source %s 2>/dev/null && conda activate isoform_tools && %s"',
              .dq(conda_sh_path %||% "$HOME/mambaforge/etc/profile.d/conda.sh"), cmd)
    else cmd
    status <- .wsl_exec_script(body, wsl_distro = distro, use_wsl = TRUE, log_dir = log_dir)
    status
  }
  .run_intern <- function(cmd) {
    .wsl_exec_script(cmd, wsl_distro = distro, use_wsl = TRUE, intern = TRUE,
                     ignore_stderr = TRUE, log_dir = log_dir)
  }

  if (use_mamba) {
    message("Installing tools using mamba/conda in WSL (default)...")

    mamba_check <- .run_intern("command -v mamba || true")
    mamba_check <- trimws(mamba_check[nzchar(trimws(mamba_check))])
    if (length(mamba_check) == 0) {
      message("mamba not found. Installing mambaforge...")
      bootstrap_cmds <- c(
        "wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O /tmp/Mambaforge.sh",
        "bash /tmp/Mambaforge.sh -b -p $HOME/mambaforge"
      )
      bootstrap_ok <- TRUE
      for (cmd in bootstrap_cmds) {
        status <- .run(cmd)
        if (!isTRUE(status == 0L)) {
          bootstrap_ok <- FALSE
          message("  \u2717 FAILED: ", cmd, " (exit code ", status, ")")
        } else {
          message("  \u2713 OK: ", cmd)
        }
      }
      if (!bootstrap_ok) {
        message("mambaforge bootstrap failed -- aborting mamba-based install. ",
                "Re-run with use_mamba = FALSE for the apt/pip path, or install mambaforge manually.")
        return(invisible(FALSE))
      }
    } else {
      message("mamba already installed (", mamba_check[1], ").")
    }

    conda_sh_path <- .find_conda_sh(wsl_distro = distro, use_wsl = TRUE)
    if (is.null(conda_sh_path)) {
      message("  ! Could not locate conda.sh after mamba install/detection -- falling back to ",
              "$HOME/mambaforge/etc/profile.d/conda.sh, which may not exist for this install. ",
              "If subsequent steps report activation failures, locate conda.sh manually and check ",
              ".find_conda_sh()'s search paths.")
    } else {
      message("  Using conda.sh: ", conda_sh_path)
    }

    env_check <- .run_intern(sprintf(
      'bash -c "source %s 2>/dev/null && conda env list | grep isoform_tools || true"',
      .dq(conda_sh_path %||% "$HOME/mambaforge/etc/profile.d/conda.sh")
    ))
    env_check <- trimws(env_check[nzchar(trimws(env_check))])
    if (length(env_check) == 0) {
      message("Creating conda environment 'isoform_tools'...")
      create_status <- .run(sprintf(
        'bash -c "source %s && conda activate base && mamba create -y -n isoform_tools python=3.9"',
        .dq(conda_sh_path %||% "$HOME/mambaforge/etc/profile.d/conda.sh")
      ))
      if (!isTRUE(create_status == 0L)) {
        message("  \u2717 Failed to create conda environment 'isoform_tools' (exit code ",
                create_status, "). Aborting -- fix this before tool installs can proceed.")
        return(invisible(FALSE))
      }
      message("  \u2713 Environment 'isoform_tools' created.")
    } else {
      message("Environment 'isoform_tools' already exists.")
    }

    install_cmds <- c("mamba install -y -c bioconda cpat",
                      "mamba install -y -c bioconda hmmer")

    for (icmd in install_cmds) {
      status <- .run(icmd, conda = TRUE)
      if (isTRUE(status == 0L)) {
        message("  \u2713 ", icmd)
      } else {
        message("  \u2717 FAILED: ", icmd, " (exit code ", status,
                "). Try running this manually inside the 'isoform_tools' conda env to see the full error.")
      }
    }

    if (!is.null(windows_signalp_dir)) {
      message("Installing SignalP from Windows directory...")
      install_signalp_from_windows(windows_signalp_dir, distro)
    } else {
      message("SignalP was not installed (no conda/apt package exists for it -- academic license). ",
              "Call install_signalp_from_windows(windows_signalp_dir, distro) once you have downloaded ",
              "a licensed copy from https://services.healthtech.dtu.dk/services/SignalP-6.0/ to enable SignalP predictions.")
    }

    message("Tool install step complete for conda environment 'isoform_tools'. Run debug_wsl() to verify.")

  } else {
    message("Installing tools using apt/pip (legacy method)...")
    apt_ok <- TRUE
    for (cmd in c("sudo apt update", "sudo apt install -y python3-pip hmmer", "pip3 install cpat")) {
      status <- .run(cmd)
      if (isTRUE(status == 0L)) message("  \u2713 ", cmd)
      else { apt_ok <- FALSE; message("  \u2717 FAILED: ", cmd, " (exit code ", status, ")") }
    }
    if (!apt_ok) message("  ! One or more base package installs failed -- CPAT/hmmscan may not work until resolved.")

    message("Resolving latest InterProScan release...")
    iprscan_url <- .run_intern(paste(
      "curl -fsSL https://api.github.com/repos/ebi-pf-team/interproscan/releases/latest",
      "| grep -o '\"browser_download_url\": *\"[^\"]*-64-bit.tar.gz\"'",
      "| head -1 | sed 's/.*\"\\(https[^\"]*\\)\"/\\1/'"
    ))
    iprscan_url <- trimws(iprscan_url[nzchar(trimws(iprscan_url))])
    if (length(iprscan_url) == 0) {
      message("  \u2717 Could not resolve a current InterProScan download URL (network access to ",
              "api.github.com from inside WSL is required). Skipping InterProScan -- Pfam annotation ",
              "will fall back to hmmscan (already installed above), which is fully supported by this pipeline. ",
              "To install InterProScan manually later, see https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html")
    } else {
      message("  Found: ", iprscan_url[1])
      dl_status <- .run(sprintf("wget -q %s -O /tmp/interproscan.tar.gz", shQuote(iprscan_url[1], type = "sh")))
      if (!isTRUE(dl_status == 0L)) {
        message("  \u2717 Download failed (exit code ", dl_status, "). Skipping InterProScan (hmmscan fallback still available).")
      } else {
        extract_status <- .run("mkdir -p $HOME/interproscan && tar -xzf /tmp/interproscan.tar.gz -C $HOME/interproscan --strip-components=1")
        if (!isTRUE(extract_status == 0L)) {
          message("  \u2717 Extraction failed (exit code ", extract_status, "). Skipping InterProScan.")
        } else {
          java_ok <- .run("command -v java >/dev/null 2>&1")
          if (!isTRUE(java_ok == 0L)) {
            message("  \u2717 Java not found -- InterProScan requires Java 11+. Install it (e.g. `sudo apt install -y default-jre`) ",
                    "then re-run `python3 $HOME/interproscan/initial_setup.py` manually.")
          } else {
            setup_status <- .run("cd $HOME/interproscan && python3 initial_setup.py")
            if (isTRUE(setup_status == 0L)) {
              link_status <- .run("sudo ln -sf $HOME/interproscan/interproscan.sh /usr/local/bin/interproscan.sh")
              if (isTRUE(link_status == 0L)) {
                message("  \u2713 InterProScan installed and linked to /usr/local/bin/interproscan.sh")
              } else {
                message("  ! InterProScan set up but symlink step failed (exit code ", link_status, ")")
              }
            } else {
              message("  \u2717 initial_setup.py failed (exit code ", setup_status, "). InterProScan install incomplete; hmmscan fallback still available.")
            }
          }
        }
      }
    }

    if (!is.null(windows_signalp_dir)) {
      message("Installing SignalP from Windows directory...")
      install_signalp_from_windows(windows_signalp_dir, distro)
    } else {
      message("SignalP was not installed (no apt package exists for it -- academic license). ",
              "Call install_signalp_from_windows(windows_signalp_dir, distro) once you have downloaded ",
              "a licensed copy from https://services.healthtech.dtu.dk/services/SignalP-6.0/ to enable SignalP predictions.")
    }

    message("Tool install step complete via apt/pip.")
  }

  if (install_databases)
    install_isoform_databases(distro = distro, use_wsl = TRUE, log_dir = log_dir)

  message("Installation complete. Run debug_wsl(distro = ", shQuote(distro, type='sh'),
          ", use_wsl = TRUE) to verify every tool and database is now detected.")
  invisible(TRUE)
}

#' Install required databases for CPAT, Pfam, and SignalP
#'
#' Downloads CPAT hexamer tables & logit models (from GitHub) and Pfam-A.hmm
#' (from EBI FTP), then `hmmpress`-indexes the Pfam database so `hmmscan` can
#' actually use it. SignalP models require a separate license.
#'
#' @param distro        WSL distribution name
#' @param use_wsl       Run inside WSL (default: TRUE on Windows, FALSE elsewhere)
#' @param cpat_data_dir Where to store CPAT data (auto-detected if NULL)
#' @param pfam_db_dir   Where to store Pfam database (default: ~/pfam_db)
#' @param conda_env     Conda environment to activate for tool checks (e.g.
#'   `hmmpress`), matching the environment install_wsl_isoform_tools() creates
#'   (default "isoform_tools")
#' @param log_dir       Optional directory to write a WSL command audit trail
#'   (via .log_wsl_command()); if NULL, only console messages are produced.
#' @export
install_isoform_databases <- function(distro        = "Ubuntu",
                                       use_wsl       = (.Platform$OS.type == "windows"),
                                       cpat_data_dir = NULL,
                                       pfam_db_dir   = NULL,
                                       conda_env     = "isoform_tools",
                                       log_dir       = NULL) {
  via_wsl <- use_wsl && .Platform$OS.type == "windows"


  conda_sh <- .find_conda_sh(wsl_distro = distro, use_wsl = via_wsl)
  if (is.null(conda_sh))
    message("  (no conda.sh found -- checks below will use the bare PATH; ",
            "if you installed tools via mamba/conda, make sure the shell ",
            "that runs this function can see them, or ignore if you used the apt/pip path)")

  .run <- function(body, intern = FALSE) {
    res <- .wsl_exec_script(body, wsl_distro = distro, use_wsl = via_wsl,
                             conda_sh = conda_sh, conda_env = conda_env,
                             intern = intern, ignore_stderr = TRUE, log_dir = NULL)
    if (!is.null(log_dir)) {
      .log_wsl_command(if (length(body) > 1) paste(body, collapse = "; ") else body,
                       exit_code = if (intern) (attr(res, "status") %||% 0L) else res,
                       stdout = if (intern) res else NULL, log_dir = log_dir)
    }
    res
  }

  # ---- CPAT databases ----
  message("Installing CPAT hexamer and logit models...")
  find_cpat <- .run("command -v cpat || command -v run_cpat.py || true", intern = TRUE)
  find_cpat <- trimws(find_cpat[nzchar(trimws(find_cpat))])
  if (length(find_cpat) > 0) {
    cpat_base     <- dirname(dirname(find_cpat[1]))
    cpat_data_dir <- file.path(cpat_base, "data")
    message("  CPAT found at: ", find_cpat[1], ". Installing data to: ", cpat_data_dir)
  } else {
    if (is.null(cpat_data_dir)) cpat_data_dir <- "$HOME/.cpat_data"
    message("  CPAT not found on PATH/conda env (this is OK if you haven't installed it yet). ",
            "Installing data to: ", cpat_data_dir)
  }
  mkdir_status <- .run(sprintf("mkdir -p %s", .dq(cpat_data_dir)))
  if (!isTRUE(mkdir_status == 0L))
    message("  ! Could not create ", cpat_data_dir, " (exit code ", mkdir_status, ") -- check permissions.")

  cpat_urls <- c(
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Human_Hexamer.tsv/download",
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Mouse_Hexamer.tsv/download",
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Human_logitModel.RData/download",
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Mouse_logitModel.RData/download"
  )
  cpat_ok <- TRUE
  for (url in cpat_urls) {
    fname  <- basename(sub("/download$", "", url))
    dest   <- paste0(cpat_data_dir, "/", fname)
    status <- .run(sprintf('wget -q -L --max-redirect=20 --user-agent="Mozilla/5.0" -O %s %s',
                          .dq(dest), .dq(url)))
    if (isTRUE(status == 0L)) {
      message("  \u2713 Downloaded: ", fname)
    } else {
      cpat_ok <- FALSE
      message("  \u2717 FAILED to download ", fname, " (exit code ", status,
              "). Check network access to sourceforge.net inside the execution environment, ",
              "or download it manually from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/ ",
              "and place it in ", cpat_data_dir)
    }
  }

  .wsl_write_env_var("CPAT_DATA", cpat_data_dir, wsl_distro = distro, use_wsl = via_wsl)
  message(if (cpat_ok) "CPAT data installed successfully. CPAT_DATA -> " else
          "CPAT data installed with errors (see above). CPAT_DATA -> ", cpat_data_dir)

  # ---- Pfam database ----
  message("Installing Pfam-A.hmm...")
  if (is.null(pfam_db_dir)) pfam_db_dir <- "$HOME/pfam_db"
  .run(sprintf("mkdir -p %s", .dq(pfam_db_dir)))

  pfam_url    <- "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
  pfam_hmm    <- paste0(pfam_db_dir, "/Pfam-A.hmm")
  pfam_gz     <- paste0(pfam_hmm, ".gz")
  dl_status   <- .run(sprintf("wget -q -O %s %s", .dq(pfam_gz), .dq(pfam_url)))
  if (!isTRUE(dl_status == 0L)) {
    message("  \u2717 FAILED to download Pfam-A.hmm.gz (exit code ", dl_status,
            "). Check network access to ftp.ebi.ac.uk inside the execution environment. Aborting Pfam install.")
  } else {
    message("  \u2713 Downloaded Pfam-A.hmm.gz")
    gz_status <- .run(sprintf("gunzip -f %s", .dq(pfam_gz)))
    if (!isTRUE(gz_status == 0L)) {
      message("  \u2717 gunzip failed (exit code ", gz_status, ") -- Pfam-A.hmm.gz may be corrupt or incomplete.")
    } else {
      message("  \u2713 Extracted Pfam-A.hmm")
      hmmpress_found <- .run("command -v hmmpress >/dev/null 2>&1")
      if (isTRUE(hmmpress_found == 0L)) {
        press_status <- .run(sprintf("hmmpress -f %s", .dq(pfam_hmm)))
        if (isTRUE(press_status == 0L)) {
          message("  \u2713 hmmpress indexed ", pfam_hmm, " -- hmmscan is ready to use.")
        } else {
          message("  \u2717 hmmpress FAILED (exit code ", press_status, ") on ", pfam_hmm,
                  " -- hmmscan will not work against this database until this is resolved.")
        }
      } else {
        message("  \u2717 hmmpress not found (looked ",
                if (!is.null(conda_sh)) paste0("inside conda env '", conda_env, "' and ") else "",
                "on PATH). Pfam database NOT indexed -- hmmscan will fail against it. ",
                "Install hmmer (`mamba install -c bioconda hmmer` inside '", conda_env,
                "', or install_wsl_isoform_tools()) and re-run install_isoform_databases().")
      }
    }
  }
  .wsl_write_env_var("PFAM_DB", pfam_hmm, wsl_distro = distro, use_wsl = via_wsl)
  message("Pfam database step complete. Location: ", pfam_db_dir, " (PFAM_DB -> ", pfam_hmm, ")")

  # ---- SignalP ----
  message("SignalP models require a license and cannot be automatically installed.")
  message("Use install_signalp_from_windows() to copy a local SignalP distribution into WSL.")

  message("\nAll database installation steps attempted. Review any \u2717 messages above ",
          "-- run debug_wsl() to re-check the environment once you've resolved them.")
  invisible(NULL)
}

#' Load previously saved isoform analysis step checkpoints
#'
#' Loads any available RDS files from a save directory. Recognises both the
#' legacy three-file layout (dte_results.rds, dtu_results.rds, switch_list.rds)
#' and the new per-step checkpoint files (step1_imported.rds, etc.).
#'
#' @param save_dir Directory containing saved RDS files.
#' @return Named list: isoform_import, dte_results, dtu_results, switch_list,
#'   switch_step1, switch_step2, switch_step3. Missing slots are NULL.
#' @export
load_isoform_results <- function(save_dir) {
  slots <- list(
    list(slot = "isoform_import", file = "isoform_import.rds",  label = "Isoform import"),
    list(slot = "dte_results",    file = "dte_results.rds",      label = "DTE results"),
    list(slot = "dtu_results",    file = "dtu_results.rds",      label = "DTU results"),
    list(slot = "switch_list",    file = "switch_list.rds",      label = "Final SwitchList"),
    list(slot = "switch_step1",   file = "step1_imported.rds",   label = "SwitchList step-1 (imported)"),
    list(slot = "switch_step2",   file = "step2_analyzed.rds",   label = "SwitchList step-2 (analyzed)"),
    list(slot = "switch_step3",   file = "step3_predictors.rds", label = "SwitchList step-3 (predictors)"),
    list(slot = "switch_step3_5", file = "step3_5_refreshed.rds", label = "SwitchList step-3.5 (consequences refreshed)"),
    list(slot = "dexseq_results", file = "dexseq_results.rds",    label = "DEXSeq DTU results")
  )

  res <- lapply(slots, function(s) {
    p <- file.path(save_dir, s$file)
    if (file.exists(p)) {
      message("Loaded ", s$label, " from ", p)
      readRDS(p)
    } else NULL
  })
  names(res) <- vapply(slots, `[[`, character(1), "slot")

  if (all(vapply(res, is.null, logical(1))))
    warning("No saved RDS checkpoint files found in: ", save_dir)

  res
}