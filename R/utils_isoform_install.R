# utils_isoform_install.R -- one-time setup helpers for installing SignalP/CPAT/Pfam
# tools and reference databases into the WSL isoform_tools conda environment

#' Install SignalP from a Windows directory into WSL
#' @export
install_signalp_from_windows <- function(windows_signalp_dir,
                                         distro = "Ubuntu-22.04",
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

#' Install SignalP 6.0 from a `signalp-6.0*.tar.gz` distribution archive
#'
#' Unlike `install_signalp_from_windows()` (which copies an *already installed*
#' SignalP directory, e.g. from a prior Linux/native install, and expects a
#' `bin/signalp` executable + `data/` directory), this installs directly from
#' the `.tar.gz` package DTU distributes: it is a Python package (`pip install
#' signalp-6-package/`) that provides the `signalp6` command, plus a separate
#' `models/` folder of weight files that must be copied into the installed
#' package's `model_weights/` directory (see the archive's own
#' `signalp6_fast/signalp-6-package/README.md`).
#'
#' @param tarball_path Path to the `signalp-6.0*.tar.gz` archive on Windows
#'   (or already inside WSL). Defaults to the copy bundled at
#'   `inst/extdata/signalp-6.0i.fast.tar.gz`, resolved via
#'   `system.file("extdata", ..., package = "ExpressOM")` with a source-checkout
#'   fallback -- pass an explicit path to use a different archive (e.g. a
#'   `.fast.tar.gz`, `.slow.tar.gz`, or a newer release).
#' @param distro WSL distribution name.
#' @param conda_env Conda environment to install the `signalp` Python package
#'   into (created via `install_wsl_isoform_tools()` / `mamba create`).
#' @param log_dir Directory for `wsl_commands.log`; defaults per
#'   `.default_wsl_log_dir()`.
#' @export
install_signalp_from_tarball <- function(tarball_path = NULL,
                                         distro = "Ubuntu-22.04",
                                         conda_env = "isoform_tools",
                                         log_dir = NULL) {

  if (is.null(tarball_path)) {
    tarball_path <- system.file("extdata", "signalp-6.0i.fast.tar.gz", package = "ExpressOM")

    if (!nzchar(tarball_path) || !file.exists(tarball_path)) {
      tarball_path <- "inst/extdata/signalp-6.0i.fast.tar.gz"
    }
  }

  if (!file.exists(tarball_path)) {
    stop(
      "SignalP archive not found: ", tarball_path,
      ". Pass tarball_path = ... pointing at your downloaded ",
      "signalp-6.0*.tar.gz (from https://services.healthtech.dtu.dk/services/SignalP-6.0/)."
    )
  }

  effective_log_dir <- if (!is.null(log_dir)) log_dir else .default_wsl_log_dir()
  message("Logging every command to: ", file.path(effective_log_dir, "wsl_commands.log"))

  if (!check_wsl(distro)) {
    stop("WSL with distro ", distro, " not available.")
  }

  tarball_path <- normalizePath(tarball_path, winslash = "/", mustWork = TRUE)
  wsl_tarball_path <- .to_wsl_path(tarball_path, distro)

  message("Installing SignalP 6.0 from archive: ", tarball_path)

  conda_sh <- .find_conda_sh(wsl_distro = distro, use_wsl = TRUE, log_dir = effective_log_dir)

  if (is.null(conda_sh)) {
    stop(
      "Could not locate conda.sh in WSL -- install conda/mamba first ",
      "(see install_wsl_isoform_tools()) before installing SignalP."
    )
  }

  env_check <- .wsl_exec_script(
    bash_body = sprintf(
      '. %s 2>/dev/null && conda env list | grep -q %s',
      .dq(conda_sh),
      shQuote(conda_env, type = "sh")
    ),
    wsl_distro = distro,
    use_wsl = TRUE,
    intern = FALSE,
    ignore_stderr = TRUE,
    log_dir = effective_log_dir
  )

  if (!isTRUE(env_check == 0L)) {
    stop(
      "Conda environment '", conda_env, "' not found. Run install_wsl_isoform_tools() ",
      "(or create it manually) before installing SignalP into it."
    )
  }

  extract_dir <- "/tmp/expressom_signalp6_install"

  message("  Extracting archive inside WSL to ", extract_dir, "...")

  extract_status <- .wsl_exec_script(
    bash_body = c(
      sprintf("rm -rf %s", .dq(extract_dir)),
      sprintf("mkdir -p %s", .dq(extract_dir)),
      sprintf("tar -xzf %s -C %s", .dq(wsl_tarball_path), .dq(extract_dir))
    ),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (!isTRUE(extract_status == 0L)) {
    stop("Failed to extract ", tarball_path, " inside WSL (exit code ", extract_status, ").")
  }

  pkg_dir_probe <- .wsl_exec_script(
    bash_body = sprintf(
      'find %s -maxdepth 3 -type d -name "signalp-6-package" | head -1',
      .dq(extract_dir)
    ),
    wsl_distro = distro,
    use_wsl = TRUE,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = effective_log_dir
  )

  pkg_dir <- trimws(pkg_dir_probe[nzchar(trimws(pkg_dir_probe))])

  if (length(pkg_dir) == 0) {
    stop(
      "Could not find a 'signalp-6-package' directory inside the extracted archive. ",
      "The archive layout may differ from the expected 'signalp6_fast/signalp-6-package/' ",
      "structure -- inspect ", extract_dir, " inside WSL manually."
    )
  }

  pkg_dir <- pkg_dir[1]
  message("  Found package directory: ", pkg_dir)
  message("  Installing signalp Python package into conda env '", conda_env, "'...")

  pip_status <- .wsl_exec_script(
    bash_body = sprintf(
      '. %s && conda activate %s && pip install %s',
      .dq(conda_sh),
      shQuote(conda_env, type = "sh"),
      .dq(pkg_dir)
    ),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (!isTRUE(pip_status == 0L)) {
    stop(
      "pip install of the signalp package failed (exit code ", pip_status,
      "). See wsl_commands.log for the full pip output."
    )
  }

  message("  Copying model weight files into the installed package...")

  signalp_dir_probe <- .wsl_exec_script(
    bash_body = sprintf(
      '. %s && conda activate %s && python3 -c "import signalp, os; print(os.path.dirname(signalp.__file__))"',
      .dq(conda_sh),
      shQuote(conda_env, type = "sh")
    ),
    wsl_distro = distro,
    use_wsl = TRUE,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = effective_log_dir
  )

  signalp_dir <- trimws(signalp_dir_probe[nzchar(trimws(signalp_dir_probe))])

  if (length(signalp_dir) == 0) {
    stop(
      "pip install appeared to succeed, but could not locate the installed 'signalp' ",
      "Python package to copy model weights into. Check wsl_commands.log."
    )
  }

  signalp_dir <- signalp_dir[1]
  model_weights_dir <- file.path(signalp_dir, "model_weights")

  copy_status <- .wsl_exec_script(
    bash_body = sprintf(
      "mkdir -p %s && cp -r %s/models/* %s/",
      .dq(model_weights_dir),
      .dq(pkg_dir),
      .dq(model_weights_dir)
    ),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (!isTRUE(copy_status == 0L)) {
    message(
      "  ! Failed to copy model weight files into ", model_weights_dir, " (exit code ", copy_status,
      "). signalp6 will likely fail at runtime with a missing-weights error; copy them manually, ",
      "e.g. from ", pkg_dir, "/models/ inside WSL."
    )
  } else {
    message("  Model weights copied to: ", model_weights_dir)
  }

  verify_status <- .wsl_exec_script(
    bash_body = sprintf(
      '. %s && conda activate %s && command -v signalp6 >/dev/null 2>&1',
      .dq(conda_sh),
      shQuote(conda_env, type = "sh")
    ),
    wsl_distro = distro,
    use_wsl = TRUE,
    log_dir = effective_log_dir
  )

  if (isTRUE(verify_status == 0L)) {
    message(
      "SignalP 6.0 installed successfully. `signalp6` is available inside conda env '",
      conda_env, "'. Run debug_wsl() to confirm it is now detected."
    )
  } else {
    message(
      "  ! 'signalp6' command was not found on PATH inside conda env '", conda_env,
      "' after installation. The pip install reported success, but check wsl_commands.log ",
      "for details -- you may need to re-activate the conda env or check for a console-script ",
      "entry-point issue."
    )
  }

  .wsl_exec_script(
    sprintf("rm -rf %s", .dq(extract_dir)),
    wsl_distro = distro,
    use_wsl = TRUE,
    ignore_stderr = TRUE,
    log_dir = effective_log_dir
  )

  invisible(isTRUE(verify_status == 0L))
}

#' Install required external tools inside WSL using mamba or apt/pip
#' @export
install_wsl_isoform_tools <- function(distro = "Ubuntu-22.04",
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
        "If you already have a signalp-6.0*.tar.gz distribution archive, call ",
        "install_signalp_from_tarball(tarball_path, distro) to install it directly; ",
        "if you have an already-installed SignalP directory (e.g. from a prior native/Linux ",
        "install), call install_signalp_from_windows(windows_signalp_dir, distro) instead. ",
        "Download a licensed copy from https://services.healthtech.dtu.dk/services/SignalP-6.0/ ",
        "if you don't have one yet."
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
        "If you already have a signalp-6.0*.tar.gz distribution archive, call ",
        "install_signalp_from_tarball(tarball_path, distro) to install it directly; ",
        "if you have an already-installed SignalP directory (e.g. from a prior native/Linux ",
        "install), call install_signalp_from_windows(windows_signalp_dir, distro) instead. ",
        "Download a licensed copy from https://services.healthtech.dtu.dk/services/SignalP-6.0/ ",
        "if you don't have one yet."
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
install_isoform_databases <- function(distro = "Ubuntu-22.04",
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
