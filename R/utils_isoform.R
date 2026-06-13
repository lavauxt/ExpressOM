# ==============================================================================
# utils_isoform.R - WSL/conda helpers + install/database utilities
# ==============================================================================

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
#' @return Integer exit code (0 = success) or character vector when intern=TRUE.
#' @keywords internal
.wsl_exec_script <- function(bash_body,
                              wsl_distro    = "Ubuntu",
                              use_wsl       = TRUE,
                              conda_sh      = NULL,
                              conda_env     = "isoform_tools",
                              intern        = FALSE,
                              ignore_stderr = TRUE) {

  activate <- if (!is.null(conda_sh) && nzchar(conda_sh)) {
    c(sprintf('. "%s" 2>/dev/null || true', conda_sh),
      sprintf("conda activate %s 2>/dev/null || true", conda_env))
  } else character(0)

  script <- c("#!/bin/bash", "set -e", activate, bash_body)
  tmp    <- tempfile(fileext = ".sh")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(script, tmp)

  via_wsl <- (.Platform$OS.type == "windows") && use_wsl

  if (via_wsl) {
    wsl_sh <- .to_wsl_path(tmp, wsl_distro)
    args   <- c("-d", wsl_distro, "bash", wsl_sh)
    if (intern)
      system2("wsl", args, stdout = TRUE, stderr = FALSE)
    else
      as.integer(system2("wsl", args,
                          stdout = if (ignore_stderr) FALSE else "",
                          stderr = if (ignore_stderr) FALSE else ""))
  } else {
    if (intern)
      system2("bash", tmp, stdout = TRUE, stderr = FALSE)
    else
      as.integer(system2("bash", tmp,
                          stdout = if (ignore_stderr) FALSE else "",
                          stderr = if (ignore_stderr) FALSE else ""))
  }
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
    intern        = FALSE,     ignore_stderr = TRUE
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
    intern = TRUE, ignore_stderr = TRUE
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
    intern = TRUE, ignore_stderr = TRUE
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
