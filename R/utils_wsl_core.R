# utils_wsl_core.R -- low-level WSL/shell command execution primitives
# shared by every external-tool integration (CPAT, SignalP, DeepTMHMM, Pfam, etc.)

# Small in-session cache to avoid repeatedly shelling out for the same
# conda.sh / distro combination during a single predictor run.
#' @keywords internal
.expressom_cache <- new.env(parent = emptyenv())

# Shared persistent-environment-variable file.
.ISOFORM_ENV_FILE <- "$HOME/.isoform_tools_env.sh"

#' Quote a literal shell argument, preventing variable and command expansion
#' @keywords internal
.dq <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub('"', '\\\"', x, fixed = TRUE)
  x <- gsub('`', '\\`', x, fixed = TRUE)
  x <- gsub("$", "\\$", x, fixed = TRUE)
  sprintf('"%s"', x)
}

.dq_home <- function(x) {
  x <- as.character(x)
  if (length(x) != 1L || !(identical(x, "$HOME") || startsWith(x, "$HOME/"))) {
    return(.dq(x))
  }
  suffix <- substring(x, 6L)
  if (startsWith(suffix, "/")) suffix <- substring(suffix, 2L)
  paste0('"${HOME}"', if (nzchar(suffix)) paste0("/", .dq(suffix)) else "")
}

#' Resolve where WSL command logs are written when a caller doesn't specify
#' log_dir explicitly.
#'
#' Priority:
#'   1. `options(expressom.wsl_log_dir)` -- set by `expressom()` /
#'      `run_expressom_targets()` to `<out_dir>/Log/WSL`, so every WSL command
#'      a pipeline run issues lands under the run's `Log/` folder next to
#'      SessionInfo.txt, Warnings.txt, the DGE/isoform logs, etc.
#'   2. `<getwd()>/expressom_wsl_logs` -- fallback when the installer
#'      functions are called outside any pipeline run (e.g. a one-off
#'      `install_wsl_isoform_tools()`).
#' @keywords internal
.default_wsl_log_dir <- function() {
  opt <- getOption("expressom.wsl_log_dir")
  if (!is.null(opt) && length(opt) == 1L && nzchar(opt)) {
    if (!dir.exists(opt)) {
      dir.create(opt, recursive = TRUE, showWarnings = FALSE)
    }
    return(opt)
  }
  d <- file.path(getwd(), "expressom_wsl_logs")
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  d
}

#' Persist an environment variable for later predictor runs
#' @keywords internal
.wsl_write_env_var <- function(
  var,
  value,
  wsl_distro = "Ubuntu",
  use_wsl = TRUE,
  log_dir = NULL
) {
  export_value <- if (
    length(value) == 1L &&
      !is.na(value) &&
      startsWith(as.character(value), "$HOME/")
  ) {
    .dq_home(value)
  } else {
    .dq(value)
  }
  export_line <- sprintf("export %s=%s", var, export_value)
  tmp_env <- paste0(.ISOFORM_ENV_FILE, ".tmp")

  body <- c(
    sprintf("touch %s", .dq_home(.ISOFORM_ENV_FILE)),
    sprintf(
      "(grep -v %s %s || true) > %s",
      .dq(paste0("^export ", var, "=")),
      .dq_home(.ISOFORM_ENV_FILE),
      .dq_home(tmp_env)
    ),
    paste0(
      "printf '%s\\n' ",
      shQuote(export_line, type = "sh"),
      " >> ",
      .dq_home(tmp_env)
    ),
    sprintf("mv %s %s", .dq_home(tmp_env), .dq_home(.ISOFORM_ENV_FILE))
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
    message("  -> Persisted ", var, " for future predictor runs (",
            .ISOFORM_ENV_FILE, ")")
  } else {
    message(
      "  ! Could not persist ", var, " to ", .ISOFORM_ENV_FILE,
      " (non-fatal; export it manually inside the execution environment if ",
      "predictors can't find it)"
    )
  }
  invisible(ok)
}

#' Best-effort discovery of a conda profile script (conda.sh)
#' @keywords internal
.find_conda_sh <- function(
  wsl_distro = "Ubuntu",
  use_wsl = TRUE,
  log_dir = NULL
) {
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
      system2("wsl", c("-d", distro, "wslpath", "-u", p),
              stdout = TRUE, stderr = FALSE, timeout = 15)
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
.log_wsl_command <- function(
  cmd,
  exit_code,
  stdout = NULL,
  stderr = NULL,
  log_dir
) {
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
#' @keywords internal
.wsl_exec_script <- function(
  bash_body,
  wsl_distro = "Ubuntu",
  use_wsl = TRUE,
  conda_sh = NULL,
  conda_env = "isoform_tools",
  intern = FALSE,
  ignore_stderr = TRUE,
  log_dir = NULL,
  verbose = TRUE
) {
  effective_log_dir <- if (!is.null(log_dir)) log_dir else .default_wsl_log_dir()
  cmd_preview <- paste(bash_body, collapse = "; ")

  if (isTRUE(verbose)) {
    message("  [wsl] $ ",
      if (nchar(cmd_preview) > 200) paste0(substr(cmd_preview, 1, 200), " ...")
      else cmd_preview)
  }

  env_file_source <- sprintf(
    '[ -f %s ] && . %s 2>/dev/null || true',
    .dq_home(.ISOFORM_ENV_FILE), .dq_home(.ISOFORM_ENV_FILE)
  )
  activate <- env_file_source

  if (!is.null(conda_sh) && nzchar(conda_sh)) {
    activate <- c(
      activate,
      sprintf('. %s 2>/dev/null || true', .dq(conda_sh)),
      sprintf("conda activate %s 2>/dev/null || true",
              shQuote(conda_env, type = "sh"))
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
      if (!is.null(res$error)) missing_exe_msg <- res$error
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
      if (!is.null(res$error)) missing_exe_msg <- res$error
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
        message("  ... (", length(out) - 30, " more lines in ",
                file.path(effective_log_dir, "wsl_commands.log"), ")")
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
    if (!is.null(missing_exe_msg)) attr(out, "error") <- missing_exe_msg
    return(out)
  }
  run_status
}

#' Check whether a command-line tool is accessible in the execution environment
#' @keywords internal
.wsl_tool_exists <- function(
  tool_name,
  wsl_distro = "Ubuntu",
  use_wsl = TRUE,
  conda_sh = NULL,
  conda_env = "isoform_tools",
  log_dir = NULL
) {
  status <- .wsl_exec_script(
    bash_body = sprintf("command -v %s >/dev/null 2>&1",
                        shQuote(tool_name, type = "sh")),
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

#' Fingerprint of the external predictor environment
#'
#' @param log_dir Directory for the WSL command log; `NULL` falls back to
#'   `getOption("expressom.wsl_log_dir")` or `<cwd>/expressom_wsl_logs`. The
#'   pipeline always passes `<out_dir>/Log/WSL`.
#' @keywords internal
.predictor_environment_fingerprint <- function(
  wsl_distro = "Ubuntu",
  use_wsl    = (.Platform$OS.type == "windows"),
  conda_env  = "isoform_tools",
  log_dir    = NULL
) {
  if (is.null(wsl_distro) || !nzchar(wsl_distro)) wsl_distro <- "Ubuntu"
  if (is.null(use_wsl)) use_wsl <- (.Platform$OS.type == "windows")
  if (is.null(conda_env) || !nzchar(conda_env)) conda_env <- "isoform_tools"

  conda_sh <- tryCatch(
    .find_conda_sh(wsl_distro = wsl_distro, use_wsl = use_wsl, log_dir = log_dir),
    error = function(e) NULL
  )

  script <- c(
    'for t in cpat run_cpat.py signalp6 hmmscan interproscan.sh biolib; do',
    '  command -v "$t" >/dev/null 2>&1 && echo "TOOL:$t"',
    'done',
    '[ -n "$PFAM_DB" ] && [ -f "$PFAM_DB" ] && echo "PFAM:found"',
    'for _p in "$HOME/pfam_db/Pfam-A.hmm" "$HOME/databases/Pfam-A.hmm" \\',
    '          "$HOME/.local/share/pfam/Pfam-A.hmm" "/opt/databases/Pfam-A.hmm"; do',
    '  [ -f "$_p" ] && { echo "PFAM:found"; break; }',
    'done',
    'for _p in "$HOME/.cpat_data/Human_logitModel.RData" \\',
    '          "${CPAT_DATA:-.}/Human_logitModel.RData"; do',
    '  [ -f "$_p" ] && { echo "CPAT_HUMAN:found"; break; }',
    'done',
    'for _p in "$HOME/.cpat_data/Mouse_logitModel.RData" \\',
    '          "${CPAT_DATA:-.}/Mouse_logitModel.RData"; do',
    '  [ -f "$_p" ] && { echo "CPAT_MOUSE:found"; break; }',
    'done'
  )

  out <- tryCatch(
    .wsl_exec_script(
      bash_body     = script,
      wsl_distro    = wsl_distro,
      use_wsl       = use_wsl,
      conda_sh      = conda_sh,
      conda_env     = conda_env,
      intern        = TRUE,
      ignore_stderr = TRUE,
      log_dir       = log_dir,
      verbose       = FALSE
    ),
    error = function(e) character(0)
  )

  tools <- sort(unique(sub("^TOOL:", "", out[grepl("^TOOL:", out)])))

  .object_md5(list(
    tools      = tools,
    pfam_db    = any(grepl("^PFAM:", out)),
    cpat_human = any(grepl("^CPAT_HUMAN:", out)),
    cpat_mouse = any(grepl("^CPAT_MOUSE:", out))
  ))
}