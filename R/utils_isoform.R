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

  activate <- if (!is.null(conda_sh) && nzchar(conda_sh)) {
    c(sprintf('. "%s" 2>/dev/null || true', conda_sh),
      sprintf("conda activate %s 2>/dev/null || true", conda_env))
  } else character(0)

  script <- c("#!/bin/bash", "set -e", activate, bash_body)
  tmp    <- tempfile(fileext = ".sh")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(script, tmp, useBytes = TRUE)

  # Route through WSL only on Windows AND when explicitly requested. Every
  # other case (native Linux, macOS, or R already running inside WSL itself)
  # executes the script with the local `bash` -- no WSL concept applies there.
  via_wsl <- (.Platform$OS.type == "windows") && isTRUE(use_wsl)

  run_status      <- 127L
  out             <- character(0)
  missing_exe_msg <- NULL

  if (via_wsl) {
    if (!nzchar(Sys.which("wsl"))) {
      missing_exe_msg <- paste0(
        "'wsl' executable not found on PATH. Install WSL (wsl --install) or ",
        "set use_wsl = FALSE if the required tools are available natively."
      )
    } else {
      wsl_sh <- .to_wsl_path(tmp, wsl_distro)
      args   <- c("-d", wsl_distro, "bash", wsl_sh)
      res <- tryCatch({
        if (intern) {
          o <- suppressWarnings(system2("wsl", args,
                                        stdout = TRUE,
                                        stderr = if (ignore_stderr) FALSE else TRUE))
          list(out = o, status = attr(o, "status") %||% 0L)
        } else {
          s <- system2("wsl", args,
                      stdout = if (ignore_stderr) FALSE else "",
                      stderr = if (ignore_stderr) FALSE else "")
          list(out = character(0), status = s)
        }
      }, error = function(e) list(out = character(0), status = 127L,
                                   error = conditionMessage(e)))
      out <- res$out; run_status <- res$status
      if (!is.null(res$error)) missing_exe_msg <- res$error
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
      res <- tryCatch({
        if (intern) {
          o <- suppressWarnings(system2(bash_bin, tmp,
                                        stdout = TRUE,
                                        stderr = if (ignore_stderr) FALSE else TRUE))
          list(out = o, status = attr(o, "status") %||% 0L)
        } else {
          s <- system2(bash_bin, tmp,
                      stdout = if (ignore_stderr) FALSE else "",
                      stderr = if (ignore_stderr) FALSE else "")
          list(out = character(0), status = s)
        }
      }, error = function(e) list(out = character(0), status = 127L,
                                   error = conditionMessage(e)))
      out <- res$out; run_status <- res$status
      if (!is.null(res$error)) missing_exe_msg <- res$error
    }
  }

  if (!is.null(missing_exe_msg)) {
    warning(missing_exe_msg, call. = FALSE)
    run_status <- 127L
  }

  if (!is.null(log_dir)) {
    .log_wsl_command(paste(bash_body, collapse = "; "),
                     exit_code = run_status,
                     stdout = if (intern) out else NULL,
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
    bash_body = paste(
      'for _p in',
      '"$HOME/pfam_db/Pfam-A.hmm"',
      '"$HOME/databases/Pfam-A.hmm"',
      '"$HOME/.local/share/pfam/Pfam-A.hmm"',
      '"/opt/databases/Pfam-A.hmm"',
      '"/usr/local/share/pfam/Pfam-A.hmm"',
      '; do [ -f "$_p" ] && echo "$_p" && break; done'
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

# ---- Public utilities --------------------------------------------------------

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
#' @param out_dir Output directory for logs (if NULL, returns only the list)
#' @param conda_env Name of the conda environment to check (default "isoform_tools")
#' @param verbose Print progress messages
#' @param use_wsl Logical: route checks through WSL. Defaults to `TRUE` on
#'   Windows and is ignored (treated as `FALSE`) on every other platform.
#' @return Invisibly, a list with check results
#' @export
debug_wsl <- function(distro = "Ubuntu", out_dir = NULL,
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

  # Best-effort conda.sh discovery, reused for tool checks below so tools
  # that only exist inside the conda env (and not on the bare PATH) are
  # still detected correctly.
  conda_sh_guess <- NULL
  if (results$conda_env_exists) {
    csh <- .wsl_exec_script(
      paste(
        'for c in "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh"',
        '"$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh"',
        '"$HOME/mambaforge/etc/profile.d/conda.sh" "$HOME/miniforge3/etc/profile.d/conda.sh";',
        'do [ -f "$c" ] && echo "$c" && break; done'
      ),
      wsl_distro = distro, use_wsl = via_wsl, intern = TRUE, ignore_stderr = TRUE
    )
    csh <- trimws(csh[nzchar(trimws(csh))])
    if (length(csh) > 0) conda_sh_guess <- csh[1]
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

  cpat_logit_human <- .find_cpat_logit_model("Human", distro, via_wsl, conda_sh_guess, conda_env)
  cpat_logit_mouse <- .find_cpat_logit_model("Mouse", distro, via_wsl, conda_sh_guess, conda_env)
  results$cpat_models_found <- !is.null(cpat_logit_human) && !is.null(cpat_logit_mouse)
  if (verbose) message("  ", if (results$cpat_models_found) "\u2713" else "\u2717", " CPAT logit models")

  # 5. Write log file if out_dir provided
  if (!is.null(out_dir)) {
    log_dir <- file.path(out_dir, "Log")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, "wsl_debug.json")
    jsonlite::write_json(results, log_file, pretty = TRUE, auto_unbox = TRUE)
    if (verbose) message("  \u2192 Log written to ", log_file)
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
    system(paste('wsl -d', distro, 'wslpath -u', shQuote(windows_signalp_dir)), intern = TRUE)
  }, error = function(e) {
    path <- gsub("\\\\", "/", windows_signalp_dir)
    if (grepl("^[A-Za-z]:", path)) {
      drive <- tolower(substr(path, 1, 1))
      path  <- sub("^[A-Za-z]:", paste0("/mnt/", drive), path)
    }
    path
  })

  message("Copying SignalP from Windows to WSL...")
  system(paste('wsl -d', distro, '-- sudo mkdir -p', install_path), wait = TRUE)
  status <- system(paste('wsl -d', distro, '-- sudo cp -r', wsl_windows_path, install_path),
                   wait = TRUE)
  if (status != 0) {
    message("Failed to copy SignalP files. Check permissions and path.")
    return(FALSE)
  }

  system(paste('wsl -d', distro, '-- sudo chmod +x',
               file.path(install_path, "bin", "signalp")), wait = TRUE)
  system(paste('wsl -d', distro, '-- sudo ln -sf',
               file.path(install_path, "bin", "signalp"), '/usr/local/bin/signalp'), wait = TRUE)

  data_dir <- file.path(install_path, "data")
  if (dir.exists(data_dir))
    system(paste('wsl -d', distro, '-- echo "export SIGNALP_DIR=', data_dir, '" >> ~/.bashrc'),
           wait = TRUE)

  message("SignalP installed to ", install_path, " and linked to /usr/local/bin/signalp")
  TRUE
}

#' Install required external tools inside WSL using mamba (default) or apt/pip
#'
#' @param distro WSL distribution name (default "Ubuntu")
#' @param use_mamba Logical: use mamba/conda (default TRUE)
#' @param install_databases Logical: also install CPAT, Pfam databases
#' @param windows_signalp_dir Optional Windows path to a SignalP installation
#' @export
install_wsl_isoform_tools <- function(distro              = "Ubuntu",
                                       use_mamba           = TRUE,
                                       install_databases   = TRUE,
                                       windows_signalp_dir = NULL) {
  if (!check_wsl(distro)) stop("WSL with distro ", distro, " not available.")

  if (use_mamba) {
    message("Installing tools using mamba/conda in WSL (default)...")

    mamba_check <- system(paste('wsl -d', distro, '-- which mamba'),
                          intern = TRUE, ignore.stderr = TRUE)
    if (length(mamba_check) == 0 || mamba_check == "") {
      message("mamba not found. Installing mambaforge...")
      for (cmd in c(
        "wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O /tmp/Mambaforge.sh",
        "bash /tmp/Mambaforge.sh -b -p $HOME/mambaforge",
        'eval "$($HOME/mambaforge/bin/conda shell.bash hook)"',
        "conda init",
        "mamba init"
      )) system(paste('wsl -d', distro, '--', cmd), wait = TRUE)
    } else {
      message("mamba already installed.")
    }

    env_check <- system(
      paste('wsl -d', distro,
            '-- bash -c "source $HOME/mambaforge/etc/profile.d/conda.sh && conda env list | grep isoform_tools"'),
      intern = TRUE, ignore.stderr = TRUE
    )
    if (length(env_check) == 0 || !grepl("isoform_tools", env_check[1])) {
      message("Creating conda environment 'isoform_tools'...")
      system(paste(
        'wsl -d', distro,
        '-- bash -c "source $HOME/mambaforge/etc/profile.d/conda.sh && conda activate base && mamba create -y -n isoform_tools python=3.9"'
      ), wait = TRUE)
    }

    install_cmds <- c("mamba install -y -c bioconda cpat",
                      "mamba install -y -c bioconda hmmer")
    if (is.null(windows_signalp_dir))
      install_cmds <- c(install_cmds, "mamba install -y -c bioconda signalp")
    else
      message("Skipping conda SignalP; will install from Windows directory.")

    for (icmd in install_cmds) {
      full_cmd <- paste(
        'wsl -d', distro,
        '-- bash -c "source $HOME/mambaforge/etc/profile.d/conda.sh && conda activate isoform_tools &&',
        icmd, '"'
      )
      status <- system(full_cmd, wait = TRUE)
      if (status != 0)
        message("Warning: ", icmd, " failed. Install manually inside the environment.")
    }

    if (!is.null(windows_signalp_dir))
      install_signalp_from_windows(windows_signalp_dir, distro)

    message("Tools installed in conda environment 'isoform_tools'.")

  } else {
    message("Installing tools using apt/pip (legacy method)...")
    cmds <- c(
      "sudo apt update",
      "sudo apt install -y python3-pip hmmer",
      "pip3 install cpat",
      "wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/ -O interproscan.tar.gz",
      "tar -xzf interproscan.tar.gz",
      "sudo ln -s $PWD/interproscan-*/interproscan.sh /usr/local/bin/"
    )
    if (is.null(windows_signalp_dir))
      cmds <- c(cmds, "sudo apt install -y signalp")
    else
      message("SignalP will be installed from Windows directory.")

    for (cmd in cmds) system(paste('wsl -d', distro, '--', cmd), wait = TRUE)
    if (!is.null(windows_signalp_dir))
      install_signalp_from_windows(windows_signalp_dir, distro)

    message("Tools installed via apt/pip.")
  }

  if (install_databases)
    install_isoform_databases(distro = distro, use_wsl = TRUE)

  message("Installation complete.")
}

#' Install required databases for CPAT, Pfam, and SignalP
#'
#' Downloads CPAT hexamer tables & logit models (from GitHub) and Pfam-A.hmm
#' (from EBI FTP). SignalP models require a separate license.
#'
#' @param distro        WSL distribution name
#' @param use_wsl       Run inside WSL (default FALSE)
#' @param cpat_data_dir Where to store CPAT data (auto-detected if NULL)
#' @param pfam_db_dir   Where to store Pfam database (default: ~/pfam_db)
#' @export
install_isoform_databases <- function(distro        = "Ubuntu",
                                       use_wsl       = FALSE,
                                       cpat_data_dir = NULL,
                                       pfam_db_dir   = NULL) {
  cmd_prefix <- if (use_wsl && .Platform$OS.type == "windows")
    sprintf("wsl -d %s --", distro)
  else ""

  # ---- CPAT databases ----
  message("Installing CPAT hexamer and logit models...")
  find_cpat <- system(paste(cmd_prefix, "which run_cpat.py"),
                      intern = TRUE, ignore.stderr = TRUE)
  if (length(find_cpat) > 0 && find_cpat[1] != "") {
    cpat_base     <- dirname(dirname(find_cpat[1]))
    cpat_data_dir <- file.path(cpat_base, "data")
    message("CPAT found at: ", cpat_base, ". Installing data to: ", cpat_data_dir)
  } else {
    if (is.null(cpat_data_dir)) cpat_data_dir <- file.path(Sys.getenv("HOME"), ".cpat_data")
    message("CPAT not in PATH. Installing data to: ", cpat_data_dir)
  }
  system(paste(cmd_prefix, "mkdir -p", cpat_data_dir), wait = TRUE)

  for (url in c(
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Human_Hexamer.tsv",
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Mouse_Hexamer.tsv",
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Human_logitModel.RData",
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Mouse_logitModel.RData"
  )) {
    dest <- file.path(cpat_data_dir, basename(url))
    system(sprintf('%s wget -O %s %s', cmd_prefix, dest, url), wait = TRUE)
    message("  Downloaded: ", basename(url))
  }

  if (nchar(cmd_prefix) == 0)
    Sys.setenv(CPAT_DATA = cpat_data_dir)
  else
    system(paste(cmd_prefix, 'echo "export CPAT_DATA=', cpat_data_dir, '" >> ~/.bashrc'), wait = TRUE)
  message("CPAT data installed. CPAT_DATA -> ", cpat_data_dir)

  # ---- Pfam database ----
  message("Installing Pfam-A.hmm...")
  if (is.null(pfam_db_dir)) pfam_db_dir <- file.path(Sys.getenv("HOME"), "pfam_db")
  system(paste(cmd_prefix, "mkdir -p", pfam_db_dir), wait = TRUE)

  pfam_url <- "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
  pfam_gz  <- file.path(pfam_db_dir, "Pfam-A.hmm.gz")
  system(sprintf('%s wget -O %s %s', cmd_prefix, pfam_gz, pfam_url), wait = TRUE)
  system(paste(cmd_prefix, "gunzip -f", pfam_gz), wait = TRUE)

  if (system(paste(cmd_prefix, "which hmmpress"),
             ignore.stdout = TRUE, ignore.stderr = TRUE) == 0)
    system(paste(cmd_prefix, "hmmpress", file.path(pfam_db_dir, "Pfam-A.hmm")), wait = TRUE)
  else
    message("hmmpress not found. Pfam database not indexed. Install hmmer to use hmmscan.")
  message("Pfam database installed at: ", pfam_db_dir)

  # ---- SignalP ----
  message("SignalP models require a license and cannot be automatically installed.")
  message("Use install_signalp_from_windows() to copy a local SignalP distribution into WSL.")

  message("\nAll databases installed successfully.")
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
    list(slot = "switch_step3",   file = "step3_predictors.rds", label = "SwitchList step-3 (predictors)")
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