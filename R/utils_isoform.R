#' Check if WSL is available and configured
#'
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
install_signalp_from_windows <- function(windows_signalp_dir, distro = "Ubuntu", install_path = "/usr/local/signalp") {
  if (!check_wsl(distro)) {
    stop("WSL with distro ", distro, " not available.")
  }
  if (!dir.exists(windows_signalp_dir)) {
    stop("Windows directory does not exist: ", windows_signalp_dir)
  }
  
  # Convert Windows path to WSL-compatible path using wslpath if available
  wsl_windows_path <- tryCatch({
    system(paste('wsl -d', distro, 'wslpath -u', shQuote(windows_signalp_dir)), intern = TRUE)
  }, error = function(e) {
    path <- gsub("\\\\", "/", windows_signalp_dir)
    if (grepl("^[A-Za-z]:", path)) {
      drive <- tolower(substr(path, 1, 1))
      path <- sub("^[A-Za-z]:", paste0("/mnt/", drive), path)
    }
    path
  })
  
  message("Copying SignalP from Windows to WSL...")
  # Create destination directory in WSL
  system(paste('wsl -d', distro, '-- sudo mkdir -p', install_path), wait = TRUE)
  # Copy files recursively (preserving permissions)
  cp_cmd <- paste('wsl -d', distro, '-- sudo cp -r', wsl_windows_path, install_path)
  status <- system(cp_cmd, wait = TRUE)
  if (status != 0) {
    message("Failed to copy SignalP files. Check permissions and path.")
    return(FALSE)
  }
  
  # Set executable permissions on signalp binary
  system(paste('wsl -d', distro, '-- sudo chmod +x', file.path(install_path, "bin", "signalp")), wait = TRUE)
  
  # Create symlink in /usr/local/bin
  system(paste('wsl -d', distro, '-- sudo ln -sf', file.path(install_path, "bin", "signalp"), '/usr/local/bin/signalp'), wait = TRUE)
  
  # Set environment variable for SignalP data (if needed)
  data_dir <- file.path(install_path, "data")
  if (dir.exists(data_dir)) {
    system(paste('wsl -d', distro, '-- echo "export SIGNALP_DIR=', data_dir, '" >> ~/.bashrc'), wait = TRUE)
  }
  
  message("SignalP installed successfully to ", install_path, " and linked to /usr/local/bin/signalp")
  return(TRUE)
}

#' Install required external tools inside WSL using mamba (default) or apt/pip
#'
#' @param distro WSL distribution name (default "Ubuntu")
#' @param use_mamba Logical: if TRUE (default), use mamba/conda for installation; otherwise fall back to apt/pip.
#' @param install_databases Logical: if TRUE (default), also install CPAT, Pfam, and SignalP databases.
#' @param windows_signalp_dir Optional Windows path to SignalP installation (if NULL, SignalP will be skipped or installed via mamba)
#' @export
install_wsl_isoform_tools <- function(distro = "Ubuntu", use_mamba = TRUE, install_databases = TRUE, windows_signalp_dir = NULL) {
  if (!check_wsl(distro)) stop("WSL with distro ", distro, " not available.")
  
  if (use_mamba) {
    message("Installing tools using mamba/conda in WSL (default)...")
    
    # 1. Install mambaforge if not already present
    mamba_check <- system(paste('wsl -d', distro, '-- which mamba'), intern = TRUE, ignore.stderr = TRUE)
    if (length(mamba_check) == 0 || mamba_check == "") {
      message("mamba not found. Installing mambaforge...")
      cmds_install_mamba <- c(
        "wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O /tmp/Mambaforge.sh",
        "bash /tmp/Mambaforge.sh -b -p $HOME/mambaforge",
        "eval \"$($HOME/mambaforge/bin/conda shell.bash hook)\"",
        "conda init",
        "mamba init"
      )
      for (cmd in cmds_install_mamba) {
        system(paste('wsl -d', distro, '--', cmd), wait = TRUE)
      }
    } else {
      message("mamba already installed.")
    }
    
    # 2. Create a conda environment for isoform tools (if not already exists)
    env_check <- system(paste('wsl -d', distro, '-- bash -c "source $HOME/mambaforge/etc/profile.d/conda.sh && conda env list | grep isoform_tools"'), 
                        intern = TRUE, ignore.stderr = TRUE)
    if (length(env_check) == 0 || !grepl("isoform_tools", env_check[1])) {
      message("Creating conda environment 'isoform_tools'...")
      env_cmd <- paste('wsl -d', distro, '-- bash -c "source $HOME/mambaforge/etc/profile.d/conda.sh && conda activate base && mamba create -y -n isoform_tools python=3.9"')
      system(env_cmd, wait = TRUE)
    }
    
    # 3. Install CPAT, hmmer, and (optionally) signalp via conda if no Windows dir provided
    install_cmds <- c(
      "mamba install -y -c bioconda cpat",
      "mamba install -y -c bioconda hmmer"
    )
    # Only try to install signalp via conda if no Windows directory is provided
    if (is.null(windows_signalp_dir)) {
      install_cmds <- c(install_cmds, "mamba install -y -c bioconda signalp")
    } else {
      message("Skipping conda installation of SignalP; will install from Windows directory.")
    }
    
    for (icmd in install_cmds) {
      full_cmd <- paste('wsl -d', distro, '-- bash -c "source $HOME/mambaforge/etc/profile.d/conda.sh && conda activate isoform_tools &&', icmd, '"')
      status <- system(full_cmd, wait = TRUE)
      if (status != 0) {
        message("Warning: ", icmd, " failed. You may need to install manually inside the environment.")
      }
    }
    
    # 4. If Windows SignalP directory provided, install from there
    if (!is.null(windows_signalp_dir)) {
      install_signalp_from_windows(windows_signalp_dir, distro)
    }
    
    # 5. Install InterProScan separately (manual)
    message("Installing InterProScan (manual download)...")
    interpro_cmds <- c(
      "cd /tmp",
      "wget -nc ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/ -O interproscan.tar.gz || true",
      "tar -xzf interproscan.tar.gz -C $HOME/",
      "sudo ln -sf $HOME/interproscan-*/interproscan.sh /usr/local/bin/"
    )
    for (icmd in interpro_cmds) {
      system(paste('wsl -d', distro, '--', icmd), wait = TRUE)
    }
    
    message("Tools installed in conda environment 'isoform_tools'.")
    
  } else {
    # Fallback to original apt/pip method
    message("Installing tools using apt/pip (legacy method)...")
    cmds <- c(
      "sudo apt update",
      "sudo apt install -y python3-pip hmmer",
      "pip3 install cpat",
      "wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/ -O interproscan.tar.gz",
      "tar -xzf interproscan.tar.gz",
      "sudo ln -s $PWD/interproscan-*/interproscan.sh /usr/local/bin/"
    )
    # If SignalP from Windows is provided, handle separately; otherwise try apt
    if (is.null(windows_signalp_dir)) {
      cmds <- c(cmds, "sudo apt install -y signalp")
    } else {
      message("SignalP will be installed from Windows directory.")
    }
    
    for (cmd in cmds) {
      system(paste('wsl -d', distro, '--', cmd), wait = TRUE)
    }
    
    if (!is.null(windows_signalp_dir)) {
      install_signalp_from_windows(windows_signalp_dir, distro)
    }
    
    message("Tools installed via apt/pip.")
  }
  
  # Optionally install databases
  if (install_databases) {
    install_isoform_databases(distro = distro, use_wsl = TRUE)
  }
  
  message("Installation complete.")
}

#' Install required databases for CPAT, Pfam, and SignalP (where possible)
#'
#' This function downloads:
#' - CPAT hexamer tables and logit models (from CPAT GitHub repository)
#' - Pfam HMM database (Pfam-A.hmm) for hmmscan
#' - SignalP models are NOT automatically installed (requires license)
#'
#' @param distro WSL distribution name (if using WSL)
#' @param use_wsl Logical: whether to run inside WSL (default FALSE)
#' @param cpat_data_dir Directory to store CPAT data (default: auto-detected CPAT data folder)
#' @param pfam_db_dir Directory to store Pfam database (default: ~/pfam_db)
#' @export
install_isoform_databases <- function(distro = "Ubuntu", use_wsl = FALSE,
                                      cpat_data_dir = NULL,
                                      pfam_db_dir = NULL) {
  
  # Determine command prefix for WSL
  if (use_wsl && .Platform$OS.type == "windows") {
    cmd_prefix <- sprintf('wsl -d %s --', distro)
  } else {
    cmd_prefix <- ""
  }
  
  # ---- CPAT databases ----
  message("Installing CPAT hexamer and logit models...")
  # Find CPAT installation path
  find_cpat <- system(paste(cmd_prefix, "which run_cpat.py"), intern = TRUE, ignore.stderr = TRUE)
  if (length(find_cpat) > 0 && find_cpat[1] != "") {
    cpat_base <- dirname(dirname(find_cpat[1]))
    cpat_data_dir <- file.path(cpat_base, "data")
    message("CPAT found at: ", cpat_base, ". Installing data into: ", cpat_data_dir)
  } else {
    if (is.null(cpat_data_dir)) {
      cpat_data_dir <- file.path(Sys.getenv("HOME"), ".cpat_data")
    }
    message("CPAT not found in PATH. Installing data to: ", cpat_data_dir)
  }
  
  # Create directory
  system(paste(cmd_prefix, "mkdir -p", cpat_data_dir), wait = TRUE)
  
  # Download hexamer tables and logit models from CPAT GitHub
  urls <- c(
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Human_Hexamer.tsv",
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Mouse_Hexamer.tsv",
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Human_logitModel.RData",
    "https://raw.githubusercontent.com/maqinli/CPAT_data/master/Mouse_logitModel.RData"
  )
  for (url in urls) {
    dest <- file.path(cpat_data_dir, basename(url))
    download_cmd <- sprintf('%s wget -O %s %s', cmd_prefix, dest, url)
    system(download_cmd, wait = TRUE)
    message("  Downloaded: ", basename(url))
  }
  
  # Set environment variable for CPAT
  if (nchar(cmd_prefix) == 0) {
    Sys.setenv(CPAT_DATA = cpat_data_dir)
  } else {
    system(paste(cmd_prefix, 'echo "export CPAT_DATA=', cpat_data_dir, '" >> ~/.bashrc'), wait = TRUE)
  }
  message("CPAT data installed. Set CPAT_DATA environment variable to: ", cpat_data_dir)
  
  # ---- Pfam database ----
  message("Installing Pfam HMM database (Pfam-A.hmm)...")
  if (is.null(pfam_db_dir)) {
    pfam_db_dir <- file.path(Sys.getenv("HOME"), "pfam_db")
  }
  system(paste(cmd_prefix, "mkdir -p", pfam_db_dir), wait = TRUE)
  pfam_url <- "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
  pfam_gz <- file.path(pfam_db_dir, "Pfam-A.hmm.gz")
  download_cmd <- sprintf('%s wget -O %s %s', cmd_prefix, pfam_gz, pfam_url)
  system(download_cmd, wait = TRUE)
  # Uncompress
  system(paste(cmd_prefix, "gunzip -f", pfam_gz), wait = TRUE)
  # Index the HMM database (check hmmpress availability)
  if (system(paste(cmd_prefix, "which hmmpress"), ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
    system(paste(cmd_prefix, "hmmpress", file.path(pfam_db_dir, "Pfam-A.hmm")), wait = TRUE)
  } else {
    message("hmmpress not found. Pfam database not indexed. Install hmmer to use hmmscan.")
  }
  message("Pfam database installed at: ", pfam_db_dir)
  
  # ---- SignalP databases ----
  message("SignalP models require a license and are not automatically installed.")
  message("If you have a SignalP license, place the model files in the SignalP installation directory.")
  message("Alternatively, use install_signalp_from_windows() to copy a pre-downloaded SignalP distribution.")
  
  message("\nAll databases installed successfully.")
  invisible(NULL)
}

#' Load previously saved isoform analysis results
#'
#' @param save_dir Directory containing saved RDS files
#' @return A list with components: dte_results, dtu_results, switch_list (if available)
#' @export
load_isoform_results <- function(save_dir) {
  res <- list()
  
  dte_file <- file.path(save_dir, "dte_results.rds")
  if (file.exists(dte_file)) {
    res$dte_results <- readRDS(dte_file)
    message("Loaded DTE results from ", dte_file)
  } else {
    warning("DTE results not found at ", dte_file)
  }
  
  dtu_file <- file.path(save_dir, "dtu_results.rds")
  if (file.exists(dtu_file)) {
    res$dtu_results <- readRDS(dtu_file)
    message("Loaded DTU results from ", dtu_file)
  } else {
    warning("DTU results not found at ", dtu_file)
  }
  
  switch_file <- file.path(save_dir, "switch_list.rds")
  if (file.exists(switch_file)) {
    res$switch_list <- readRDS(switch_file)
    message("Loaded SwitchList from ", switch_file)
  } else {
    warning("SwitchList not found at ", switch_file)
  }
  
  return(res)
}