# Small in-session cache to avoid repeatedly shelling out for the same
# conda.sh / distro combination during a single predictor run.
#' @keywords internal
.expressom_cache <- new.env(parent = emptyenv())

# Shared persistent-environment-variable file.
.ISOFORM_ENV_FILE <- "$HOME/.isoform_tools_env.sh"

#' Double-quote a shell argument while preserving `$VAR`-style expansion
#'
#' Unlike shQuote(x, type = "sh"), which single-quotes and freezes `$HOME`,
#' this wraps x in double quotes and escapes only characters special inside
#' double quotes, so `$HOME/path` still expands in the generated shell script.
#' @keywords internal
.dq <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub('"', '\\\"', x, fixed = TRUE)
  x <- gsub('`', '\\`', x, fixed = TRUE)

  sprintf('"%s"', x)
}

#' Resolve where WSL command logs are written when a caller doesn't specify
#' log_dir explicitly.
#'
#' Overridable globally via options(expressom.wsl_log_dir = "/some/path") so
#' a whole pipeline run can point every WSL command at one place; otherwise
#' falls back to a fixed, discoverable folder so command logging is never
#' silently dropped just because a caller left log_dir at its default.
#' @keywords internal
.default_wsl_log_dir <- function() {
  getOption("expressom.wsl_log_dir", file.path(getwd(), "expressom_wsl_logs"))
}

#' Persist an environment variable for later predictor runs
#'
#' Appends or replaces an `export VAR="value"` line in the dedicated env file
#' sourced by `.wsl_exec_script()` on every invocation.
#' @keywords internal
.wsl_write_env_var <- function(var, value, wsl_distro = "Ubuntu", use_wsl = TRUE, log_dir = NULL) {
  export_line <- sprintf("export %s=%s", var, .dq(value))
  tmp_env <- paste0(.ISOFORM_ENV_FILE, ".tmp")

  body <- c(
    sprintf("touch %s", .dq(.ISOFORM_ENV_FILE)),
    sprintf(
      "(grep -v %s %s || true) > %s",
      .dq(paste0("^export ", var, "=")),
      .dq(.ISOFORM_ENV_FILE),
      .dq(tmp_env)
    ),
    paste0(
      "printf '%s\\n' ",
      shQuote(export_line, type = "sh"),
      " >> ",
      .dq(tmp_env)
    ),
    sprintf("mv %s %s", .dq(tmp_env), .dq(.ISOFORM_ENV_FILE))
  )

  status <- .wsl_exec_script(
    bash_body = body,
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    intern = FALSE,
    ignore_stderr = TRUE,
    log_dir = log_dir
  )

  ok <- isTRUE(status == 0L)

  if (ok) {
    message("  -> Persisted ", var, " for future predictor runs (", .ISOFORM_ENV_FILE, ")")
  } else {
    message(
      "  ! Could not persist ", var, " to ", .ISOFORM_ENV_FILE,
      " (non-fatal; export it manually inside the execution environment if predictors can't find it)"
    )
  }

  invisible(ok)
}

#' Best-effort discovery of a conda profile script (conda.sh)
#' @keywords internal
.find_conda_sh <- function(wsl_distro = "Ubuntu", use_wsl = TRUE, log_dir = NULL) {
  body <- c(
    'base="$(conda info --base 2>/dev/null || true)"',
    'if [ -n "$base" ] && [ -f "$base/etc/profile.d/conda.sh" ]; then',
    '  echo "$base/etc/profile.d/conda.sh"',
    '  exit 0',
    'fi',
    '',
    'for c in \\',
    '  "$HOME/miniconda3/etc/profile.d/conda.sh" \\',
    '  "$HOME/anaconda3/etc/profile.d/conda.sh" \\',
    '  "$HOME/mambaforge/etc/profile.d/conda.sh" \\',
    '  "$HOME/miniforge3/etc/profile.d/conda.sh";',
    'do',
    '  if [ -f "$c" ]; then',
    '    echo "$c"',
    '    exit 0',
    '  fi',
    'done',
    '',
    'exit 0'
  )

  csh <- .wsl_exec_script(
    bash_body = body,
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = log_dir
  )

  csh <- trimws(csh[nzchar(trimws(csh))])

  if (length(csh) > 0) csh[1] else NULL
}

#' Convert a Windows file path to a WSL-compatible Unix path
#' @keywords internal
.to_wsl_path <- function(win_path, distro = "Ubuntu") {
  if (.Platform$OS.type != "windows") return(win_path)
  if (is.null(win_path) || !nzchar(trimws(as.character(win_path)))) return(win_path)

  p <- normalizePath(win_path, winslash = "/", mustWork = FALSE)

  r <- tryCatch(
    suppressWarnings(
      system2(
        "wsl",
        c("-d", distro, "wslpath", "-u", p),
        stdout = TRUE,
        stderr = FALSE,
        timeout = 15
      )
    ),
    error = function(e) character(0)
  )

  r <- trimws(r[nzchar(trimws(r))])
  if (length(r) > 0) return(r[1])

  if (grepl("^[A-Za-z]:/", p)) {
    return(paste0("/mnt/", tolower(substr(p, 1, 1)), substring(p, 3)))
  }

  p
}

#' Write a WSL command execution to log
#' @keywords internal
.log_wsl_command <- function(cmd, exit_code, stdout = NULL, stderr = NULL, log_dir) {
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  log_file <- file.path(log_dir, "wsl_commands.log")
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  entry <- c(
    paste0("[", timestamp, "]"),
    paste0("  CMD: ", cmd),
    paste0("  EXIT: ", exit_code)
  )

  if (!is.null(stdout) && length(stdout) > 0) {
    entry <- c(entry, paste0("  STDOUT: ", paste(stdout, collapse = "\n     ")))
  }

  if (!is.null(stderr) && length(stderr) > 0) {
    entry <- c(entry, paste0("  STDERR: ", paste(stderr, collapse = "\n     ")))
  }

  write(entry, file = log_file, append = TRUE)
}

#' Write bash commands to a temp script and execute via WSL or natively
#'
#' Every call is written to wsl_commands.log as it happens -- one entry per
#' invocation, appended immediately after that command finishes -- regardless
#' of what the caller passes for log_dir (see .default_wsl_log_dir()). Set
#' verbose = FALSE to suppress the live console trace for callers (like
#' run_tool() in mod_isoform.R) that already print their own equivalent
#' summary; the file log is unaffected by verbose either way.
#'
#' `ignore_stderr` only controls whether a non-zero exit raises an R
#' warning() -- it no longer hides the command's actual output. Previously,
#' because ignore_stderr defaulted to TRUE and every "does this exist" probe
#' hardcoded it to TRUE too, a genuine failure (as opposed to a routine
#' "not found") printed nothing but an exit code: the one piece of text that
#' would explain *why* was thrown away.
#' @keywords internal
.wsl_exec_script <- function(bash_body,
                             wsl_distro = "Ubuntu",
                             use_wsl = TRUE,
                             conda_sh = NULL,
                             conda_env = "isoform_tools",
                             intern = FALSE,
                             ignore_stderr = TRUE,
                             log_dir = NULL,
                             verbose = TRUE) {

  effective_log_dir <- if (!is.null(log_dir)) log_dir else .default_wsl_log_dir()
  cmd_preview <- paste(bash_body, collapse = "; ")

  if (isTRUE(verbose)) {
    message(
      "  [wsl] $ ",
      if (nchar(cmd_preview) > 200) paste0(substr(cmd_preview, 1, 200), " ...") else cmd_preview
    )
  }

  env_file_source <- sprintf(
    '[ -f %s ] && . %s 2>/dev/null || true',
    .dq(.ISOFORM_ENV_FILE),
    .dq(.ISOFORM_ENV_FILE)
  )

  activate <- env_file_source

  if (!is.null(conda_sh) && nzchar(conda_sh)) {
    activate <- c(
      activate,
      sprintf('. %s 2>/dev/null || true', .dq(conda_sh)),
      sprintf("conda activate %s 2>/dev/null || true", shQuote(conda_env, type = "sh"))
    )
  }

  script <- c("#!/bin/bash", "set -e", activate, bash_body)

  tmp <- tempfile(fileext = ".sh")
  on.exit(unlink(tmp), add = TRUE)

  con <- file(tmp, "wb")
  writeLines(script, con, sep = "\n")
  close(con)

  via_wsl <- (.Platform$OS.type == "windows") && isTRUE(use_wsl)

  run_status <- 127L
  out <- character(0)
  missing_exe_msg <- NULL

  run_capture <- function(cmd, args) {
    tryCatch(
      {
        o <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
        st <- attr(o, "status")
        if (is.null(st)) st <- 0L
        list(out = o, status = as.integer(st))
      },
      error = function(e) {
        list(out = character(0), status = 127L, error = conditionMessage(e))
      }
    )
  }

  if (via_wsl) {
    wsl_bin <- unname(Sys.which("wsl"))

    if (!nzchar(wsl_bin)) {
      missing_exe_msg <- paste0(
        "'wsl' executable not found on PATH. Install WSL (`wsl --install -d ",
        wsl_distro,
        "`) or set use_wsl = FALSE if the required tools are available natively."
      )
    } else {
      wsl_sh <- .to_wsl_path(tmp, wsl_distro)
      args <- c("-d", wsl_distro, "bash", wsl_sh)

      res <- run_capture(wsl_bin, args)
      out <- res$out
      run_status <- res$status

      if (!is.null(res$error)) {
        missing_exe_msg <- res$error
      }
    }
  } else {
    bash_bin <- unname(Sys.which("bash"))

    if (!nzchar(bash_bin)) {
      missing_exe_msg <- paste0(
        "'bash' executable not found on PATH. External predictor tools ",
        "(CPAT / SignalP / Pfam) require a bash shell. On native Windows ",
        "without WSL/Git-Bash these steps will be skipped."
      )
    } else {
      res <- run_capture(bash_bin, tmp)
      out <- res$out
      run_status <- res$status

      if (!is.null(res$error)) {
        missing_exe_msg <- res$error
      }
    }
  }

  if (!is.null(missing_exe_msg)) {
    message("WSL/predictor environment error: ", missing_exe_msg)
    warning(missing_exe_msg, call. = FALSE)
    run_status <- 127L
  } else if (run_status != 0L) {
    if (!isTRUE(ignore_stderr)) {
      warning("WSL/shell command exited with status ", run_status, call. = FALSE)
    }

    if (isTRUE(verbose) && length(out) > 0) {
      message("  [wsl] exit ", run_status, " -- output:")
      message(paste(utils::head(out, 30), collapse = "\n"))

      if (length(out) > 30) {
        message(
          "  ... (", length(out) - 30, " more lines in ",
          file.path(effective_log_dir, "wsl_commands.log"), ")"
        )
      }
    }
  } else if (isTRUE(verbose)) {
    message("  [wsl] exit 0")
  }

  .log_wsl_command(
    cmd = cmd_preview,
    exit_code = run_status,
    stdout = out,
    stderr = missing_exe_msg,
    log_dir = effective_log_dir
  )

  if (intern) {
    attr(out, "status") <- run_status
    if (!is.null(missing_exe_msg)) {
      attr(out, "error") <- missing_exe_msg
    }
    return(out)
  }

  run_status
}

#' Check whether a command-line tool is accessible in the execution environment
#' @keywords internal
.wsl_tool_exists <- function(tool_name,
                             wsl_distro = "Ubuntu",
                             use_wsl = TRUE,
                             conda_sh = NULL,
                             conda_env = "isoform_tools",
                             log_dir = NULL) {
  status <- .wsl_exec_script(
    bash_body = sprintf("command -v %s >/dev/null 2>&1", shQuote(tool_name, type = "sh")),
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    conda_sh = conda_sh,
    conda_env = conda_env,
    intern = FALSE,
    ignore_stderr = TRUE,
    log_dir = log_dir
  )

  isTRUE(status == 0L)
}

#' Find the Pfam-A.hmm database in common paths inside the execution environment
#' @keywords internal
.find_pfam_db <- function(wsl_distro = "Ubuntu",
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

#' Locate the CPAT logit model RData file for a given organism
#' @keywords internal
.find_cpat_logit_model <- function(organism = "Human",
                                   wsl_distro = "Ubuntu",
                                   use_wsl = TRUE,
                                   conda_sh = NULL,
                                   conda_env = "isoform_tools",
                                   log_dir = NULL) {
  fname <- paste0(organism, "_logitModel.RData")

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

#' Check if WSL is available and configured
#' @param distro WSL distribution name
#' @return Logical indicating success
#' @export
check_wsl <- function(distro = "Ubuntu") {
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
debug_wsl <- function(distro = "Ubuntu",
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

  tools_list <- c("cpat", "run_cpat.py", "signalp6", "signalp", "hmmscan", "interproscan.sh")

  for (tool in tools_list) {
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

    results$tools[[tool]] <- found

    if (verbose) {
      message("   ", if (found) "\u2713" else "\u2717", "  ", tool)
    }
  }

  tool_alternatives <- list(
    "CPAT (cpat / run_cpat.py)" = c("cpat", "run_cpat.py"),
    "SignalP (signalp6 / signalp)" = c("signalp6", "signalp"),
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

#' Install SignalP from a Windows directory into WSL
#' @export
install_signalp_from_windows <- function(windows_signalp_dir,
                                         distro = "Ubuntu",
                                         install_path = "/usr/local/signalp",
                                         log_dir = NULL) {

  effective_log_dir <- if (!is.null(log_dir)) log_dir else .default_wsl_log_dir()
  message("Logging every command to: ", file.path(effective_log_dir, "wsl_commands.log"))

  if (!check_wsl(distro)) {
    stop("WSL with distro ", distro, " not available.")
  }

  if (!dir.exists(windows_signalp_dir)) {
    stop("Windows directory does not exist: ", windows_signalp_dir)
  }

  wsl_windows_path <- tryCatch(
    {
      r <- suppressWarnings(
        system2(
          "wsl",
          c("-d", distro, "wslpath", "-u", windows_signalp_dir),
          stdout = TRUE,
          stderr = FALSE
        )
      )
      trimws(r[nzchar(trimws(r))])
    },
    error = function(e) character(0)
  )

  if (length(wsl_windows_path) == 0) {
    path <- gsub("\\", "/", windows_signalp_dir)

    if (grepl("^[A-Za-z]:", path)) {
      drive <- tolower(substr(path, 1, 1))
      path <- sub("^[A-Za-z]:", paste0("/mnt/", drive), path)
    }

    wsl_windows_path <- path
  } else {
    wsl_windows_path <- wsl_windows_path[1]
  }

  message("Copying SignalP from Windows (", wsl_windows_path, ") to WSL (", install_path, ")...")

  mkdir_status <- .wsl_exec_script(
    sprintf("sudo mkdir -p %s", .dq(install_path)),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (!isTRUE(mkdir_status == 0L)) {
    message("Failed to create ", install_path, " (exit code ", mkdir_status, "). Check sudo permissions.")
    return(FALSE)
  }

  cp_status <- .wsl_exec_script(
    sprintf("sudo cp -r %s/. %s/", .dq(wsl_windows_path), .dq(install_path)),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (!isTRUE(cp_status == 0L)) {
    message("Failed to copy SignalP files (exit code ", cp_status, "). Check permissions and path.")
    return(FALSE)
  }

  bin_path <- file.path(install_path, "bin", "signalp")

  chmod_status <- .wsl_exec_script(
    sprintf("sudo chmod +x %s", .dq(bin_path)),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  ln_status <- .wsl_exec_script(
    sprintf("sudo ln -sf %s /usr/local/bin/signalp", .dq(bin_path)),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (!isTRUE(chmod_status == 0L) || !isTRUE(ln_status == 0L)) {
    message(
      "  ! Warning: chmod/symlink step reported a non-zero exit code (chmod=", chmod_status,
      ", ln=", ln_status, "). signalp may not be directly callable as `signalp`; check ", bin_path
    )
  }

  data_dir <- file.path(install_path, "data")

  data_dir_exists <- .wsl_exec_script(
    sprintf("[ -d %s ]", .dq(data_dir)),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (isTRUE(data_dir_exists == 0L)) {
    .wsl_write_env_var("SIGNALP_DIR", data_dir, wsl_distro = distro, use_wsl = TRUE, log_dir = effective_log_dir)
  } else {
    message(
      "  ! No 'data' subdirectory found under ", install_path,
      " -- SIGNALP_DIR was not set. If this SignalP version stores models elsewhere, ",
      "set the appropriate env var manually via .wsl_write_env_var()."
    )
  }

  message("SignalP installed to ", install_path, " and linked to /usr/local/bin/signalp")

  TRUE
}

#' Install required external tools inside WSL using mamba or apt/pip
#' @export
install_wsl_isoform_tools <- function(distro = "Ubuntu",
                                      use_mamba = TRUE,
                                      install_databases = TRUE,
                                      windows_signalp_dir = NULL,
                                      log_dir = NULL) {

  if (!check_wsl(distro)) {
    stop("WSL with distro ", distro, " not available.")
  }

  log_dir <- if (!is.null(log_dir)) log_dir else .default_wsl_log_dir()
  message(
    "Logging every command in real time to: ", file.path(log_dir, "wsl_commands.log"),
    " (tail -f that file from another terminal to watch progress as it happens)"
  )

  conda_sh_path <- NULL

  .run <- function(cmd, conda = FALSE) {
    body <- if (conda) {
      sprintf(
        'bash -c "source %s 2>/dev/null && conda activate isoform_tools && %s"',
        .dq(conda_sh_path %||% "$HOME/mambaforge/etc/profile.d/conda.sh"),
        cmd
      )
    } else {
      cmd
    }

    .wsl_exec_script(body, wsl_distro = distro, use_wsl = TRUE, log_dir = log_dir)
  }

  .run_intern <- function(cmd) {
    .wsl_exec_script(
      cmd,
      wsl_distro = distro,
      use_wsl = TRUE,
      intern = TRUE,
      ignore_stderr = TRUE,
      log_dir = log_dir
    )
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
        message(
          "mambaforge bootstrap failed -- aborting mamba-based install. ",
          "Re-run with use_mamba = FALSE for the apt/pip path, or install mambaforge manually."
        )
        return(invisible(FALSE))
      }
    } else {
      message("mamba already installed (", mamba_check[1], ").")
    }

    conda_sh_path <- .find_conda_sh(wsl_distro = distro, use_wsl = TRUE)

    if (is.null(conda_sh_path)) {
      message(
        "  ! Could not locate conda.sh after mamba install/detection -- falling back to ",
        "$HOME/mambaforge/etc/profile.d/conda.sh, which may not exist for this install. ",
        "If subsequent steps report activation failures, locate conda.sh manually and check ",
        ".find_conda_sh()'s search paths."
      )
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
        message(
          "  \u2717 Failed to create conda environment 'isoform_tools' (exit code ",
          create_status, "). Aborting -- fix this before tool installs can proceed."
        )
        return(invisible(FALSE))
      }

      message("  \u2713 Environment 'isoform_tools' created.")
    } else {
      message("Environment 'isoform_tools' already exists.")
    }

    install_cmds <- c(
      "mamba install -y -c bioconda cpat",
      "mamba install -y -c bioconda hmmer"
    )

    for (icmd in install_cmds) {
      status <- .run(icmd, conda = TRUE)

      if (isTRUE(status == 0L)) {
        message("  \u2713 ", icmd)
      } else {
        message(
          "  \u2717 FAILED: ", icmd, " (exit code ", status,
          "). Try running this manually inside the 'isoform_tools' conda env to see the full error."
        )
      }
    }

    if (!is.null(windows_signalp_dir)) {
      message("Installing SignalP from Windows directory...")
      install_signalp_from_windows(windows_signalp_dir, distro, log_dir = log_dir)
    } else {
      message(
        "SignalP was not installed (no conda/apt package exists for it -- academic license). ",
        "Call install_signalp_from_windows(windows_signalp_dir, distro) once you have downloaded ",
        "a licensed copy from https://services.healthtech.dtu.dk/services/SignalP-6.0/ to enable SignalP predictions."
      )
    }

    message("Tool install step complete for conda environment 'isoform_tools'. Run debug_wsl() to verify.")
  } else {
    message("Installing tools using apt/pip (legacy method)...")

    apt_ok <- TRUE

    for (cmd in c("sudo apt update", "sudo apt install -y python3-pip hmmer", "pip3 install cpat")) {
      status <- .run(cmd)

      if (isTRUE(status == 0L)) {
        message("  \u2713 ", cmd)
      } else {
        apt_ok <- FALSE
        message("  \u2717 FAILED: ", cmd, " (exit code ", status, ")")
      }
    }

    if (!apt_ok) {
      message("  ! One or more base package installs failed -- CPAT/hmmscan may not work until resolved.")
    }

    message("Resolving latest InterProScan release...")

    iprscan_url <- .run_intern(paste(
      "curl -fsSL https://api.github.com/repos/ebi-pf-team/interproscan/releases/latest",
      "| grep -o '\"browser_download_url\": *\"[^\"]*-64-bit.tar.gz\"'",
      "| head -1 | sed 's/.*\"\\(https[^\"]*\\)\"/\\1/'"
    ))

    iprscan_url <- trimws(iprscan_url[nzchar(trimws(iprscan_url))])

    if (length(iprscan_url) == 0) {
      message(
        "  \u2717 Could not resolve a current InterProScan download URL (network access to ",
        "api.github.com from inside WSL is required). Skipping InterProScan -- Pfam annotation ",
        "will fall back to hmmscan (already installed above), which is fully supported by this pipeline. ",
        "To install InterProScan manually later, see https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html"
      )
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
            message(
              "  \u2717 Java not found -- InterProScan requires Java 11+. Install it (e.g. `sudo apt install -y default-jre`) ",
              "then re-run `python3 $HOME/interproscan/initial_setup.py` manually."
            )
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
              message(
                "  \u2717 initial_setup.py failed (exit code ", setup_status,
                "). InterProScan install incomplete; hmmscan fallback still available."
              )
            }
          }
        }
      }
    }

    if (!is.null(windows_signalp_dir)) {
      message("Installing SignalP from Windows directory...")
      install_signalp_from_windows(windows_signalp_dir, distro, log_dir = log_dir)
    } else {
      message(
        "SignalP was not installed (no apt package exists for it -- academic license). ",
        "Call install_signalp_from_windows(windows_signalp_dir, distro) once you have downloaded ",
        "a licensed copy from https://services.healthtech.dtu.dk/services/SignalP-6.0/ to enable SignalP predictions."
      )
    }

    message("Tool install step complete via apt/pip.")
  }

  if (install_databases) {
    install_isoform_databases(distro = distro, use_wsl = TRUE, log_dir = log_dir)
  }

  message(
    "Installation complete. Run debug_wsl(distro = ", shQuote(distro, type = "sh"),
    ", use_wsl = TRUE) to verify every tool and database is now detected."
  )

  invisible(TRUE)
}

#' Install required databases for CPAT, Pfam, and SignalP
#' @export
install_isoform_databases <- function(distro = "Ubuntu",
                                      use_wsl = (.Platform$OS.type == "windows"),
                                      cpat_data_dir = NULL,
                                      pfam_db_dir = NULL,
                                      conda_env = "isoform_tools",
                                      log_dir = NULL) {

  via_wsl <- use_wsl && .Platform$OS.type == "windows"
  log_dir <- if (!is.null(log_dir)) log_dir else .default_wsl_log_dir()
  message(
    "Logging every command in real time to: ", file.path(log_dir, "wsl_commands.log"),
    " (tail -f that file from another terminal to watch progress as it happens)"
  )

  conda_sh <- .find_conda_sh(wsl_distro = distro, use_wsl = via_wsl, log_dir = log_dir)

  if (is.null(conda_sh)) {
    message(
      "  (no conda.sh found -- checks below will use the bare PATH; ",
      "if you installed tools via mamba/conda, make sure the shell ",
      "that runs this function can see them, or ignore if you used the apt/pip path)"
    )
  }

  .run <- function(body, intern = FALSE) {
    .wsl_exec_script(
      body,
      wsl_distro = distro,
      use_wsl = via_wsl,
      conda_sh = conda_sh,
      conda_env = conda_env,
      intern = intern,
      ignore_stderr = TRUE,
      log_dir = log_dir
    )
  }

  message("Installing CPAT hexamer and logit models...")

  find_cpat <- .run("command -v cpat || command -v run_cpat.py || true", intern = TRUE)
  find_cpat <- trimws(find_cpat[nzchar(trimws(find_cpat))])

  if (length(find_cpat) > 0) {
    cpat_base <- dirname(dirname(find_cpat[1]))
    cpat_data_dir <- file.path(cpat_base, "data")
    message("  CPAT found at: ", find_cpat[1], ". Installing data to: ", cpat_data_dir)
  } else {
    if (is.null(cpat_data_dir)) cpat_data_dir <- "$HOME/.cpat_data"
    message(
      "  CPAT not found on PATH/conda env (this is OK if you haven't installed it yet). ",
      "Installing data to: ", cpat_data_dir
    )
  }

  mkdir_status <- .run(sprintf("mkdir -p %s", .dq(cpat_data_dir)))

  if (!isTRUE(mkdir_status == 0L)) {
    message("  ! Could not create ", cpat_data_dir, " (exit code ", mkdir_status, ") -- check permissions.")
  }

  cpat_urls <- c(
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Human_Hexamer.tsv/download",
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Mouse_Hexamer.tsv/download",
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Human_logitModel.RData/download",
    "https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Mouse_logitModel.RData/download"
  )

  cpat_ok <- TRUE

  for (url in cpat_urls) {
    fname <- basename(sub("/download$", "", url))
    dest <- paste0(cpat_data_dir, "/", fname)

    status <- .run(sprintf(
      'wget -q --max-redirect=20 --user-agent="Mozilla/5.0" -O %s %s',
      .dq(dest),
      .dq(url)
    ))

    if (isTRUE(status == 0L)) {
      message("  \u2713 Downloaded: ", fname)
    } else {
      cpat_ok <- FALSE
      message(
        "  \u2717 FAILED to download ", fname, " (exit code ", status,
        "). Check network access to sourceforge.net inside the execution environment, ",
        "or download it manually from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/ ",
        "and place it in ", cpat_data_dir
      )
    }
  }

  .wsl_write_env_var("CPAT_DATA", cpat_data_dir, wsl_distro = distro, use_wsl = via_wsl, log_dir = log_dir)

  message(
    if (cpat_ok) "CPAT data installed successfully. CPAT_DATA -> " else
      "CPAT data installed with errors (see above). CPAT_DATA -> ",
    cpat_data_dir
  )

  message("Installing Pfam-A.hmm...")

  if (is.null(pfam_db_dir)) pfam_db_dir <- "$HOME/pfam_db"

  .run(sprintf("mkdir -p %s", .dq(pfam_db_dir)))

  pfam_url <- "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
  pfam_hmm <- paste0(pfam_db_dir, "/Pfam-A.hmm")
  pfam_gz <- paste0(pfam_hmm, ".gz")

  dl_status <- .run(sprintf("wget -q -O %s %s", .dq(pfam_gz), .dq(pfam_url)))

  if (!isTRUE(dl_status == 0L)) {
    message(
      "  \u2717 FAILED to download Pfam-A.hmm.gz (exit code ", dl_status,
      "). Check network access to ftp.ebi.ac.uk inside the execution environment. Aborting Pfam install."
    )
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
          message(
            "  \u2717 hmmpress FAILED (exit code ", press_status, ") on ", pfam_hmm,
            " -- hmmscan will not work against this database until this is resolved."
          )
        }
      } else {
        message(
          "  \u2717 hmmpress not found (looked ",
          if (!is.null(conda_sh)) paste0("inside conda env '", conda_env, "' and ") else "",
          "on PATH). Pfam database NOT indexed -- hmmscan will fail against it. ",
          "Install hmmer (`mamba install -c bioconda hmmer` inside '", conda_env,
          "', or install_wsl_isoform_tools()) and re-run install_isoform_databases()."
        )
      }
    }
  }

  .wsl_write_env_var("PFAM_DB", pfam_hmm, wsl_distro = distro, use_wsl = via_wsl, log_dir = log_dir)

  message("Pfam database step complete. Location: ", pfam_db_dir, " (PFAM_DB -> ", pfam_hmm, ")")

  message("SignalP models require a license and cannot be automatically installed.")
  message("Use install_signalp_from_windows() to copy a local SignalP distribution into WSL.")
  message(
    "\nAll database installation steps attempted. Review any \u2717 messages above ",
    "-- run debug_wsl() to re-check the environment once you've resolved them."
  )

  invisible(NULL)
}

#' Load previously saved isoform analysis step checkpoints
#' @export
load_isoform_results <- function(save_dir) {
  slots <- list(
    list(slot = "isoform_import", file = "isoform_import.rds", label = "Isoform import"),
    list(slot = "dte_results", file = "dte_results.rds", label = "DTE results"),
    list(slot = "dtu_results", file = "dtu_results.rds", label = "DTU results"),
    list(slot = "switch_list", file = "switch_list.rds", label = "Final SwitchList"),
    list(slot = "switch_step1", file = "step1_imported.rds", label = "SwitchList step-1 (imported)"),
    list(slot = "switch_step2", file = "step2_analyzed.rds", label = "SwitchList step-2 (analyzed)"),
    list(slot = "switch_step3", file = "step3_predictors.rds", label = "SwitchList step-3 (predictors)"),
    list(slot = "switch_step3_5", file = "step3_5_refreshed.rds", label = "SwitchList step-3.5 (consequences refreshed)"),
    list(slot = "dexseq_results", file = "dexseq_results.rds", label = "DEXSeq DTU results")
  )

  res <- lapply(slots, function(s) {
    p <- file.path(save_dir, s$file)

    if (file.exists(p)) {
      message("Loaded ", s$label, " from ", p)
      readRDS(p)
    } else {
      NULL
    }
  })

  names(res) <- vapply(slots, `[[`, character(1), "slot")

  if (all(vapply(res, is.null, logical(1)))) {
    warning("No saved RDS checkpoint files found in: ", save_dir)
  }

  res
}