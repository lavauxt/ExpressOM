#' Safely create a directory
#' @keywords internal
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

#' Safely save a ggplot or base R plot to PDF
#'
#' Uses non-standard evaluation so that `expr` is captured unevaluated and only
#' executed *after* the PDF device is opened — fixing the bug where the previous
#' `force(expr)` approach evaluated the plot expression before the device was
#' ready, sending output to whatever device happened to be current at the time.
#' @keywords internal
safe_pdf <- function(path, expr, width = 10, height = 8) {
  expr_sub   <- substitute(expr)        # capture unevaluated
  caller_env <- parent.frame()          # capture caller's environment
  tryCatch({
    pdf(path, width = width, height = height)
    eval(expr_sub, envir = caller_env)  # evaluate AFTER pdf() is open
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("Warning: Failed to generate plot at: ", path, "\n  Error: ", e$message)
  })
}

#' Safely run an expression, returning NULL on error with an optional message
#' @keywords internal
safe_run <- function(expr, label = "") {
  tryCatch(expr, error = function(e) {
    if (nchar(label) > 0) message("Warning: ", label, " failed. Skipping. Error: ", e$message)
    NULL
  })
}

#' Load an OrgDb annotation package by name (e.g. "org.Hs.eg.db")
#' @keywords internal
.load_org_db <- function(org_db_name) {
  if (!requireNamespace(org_db_name, quietly = TRUE)) {
    stop("Package '", org_db_name, "' is required. Please install it.")
  }
  getExportedValue(org_db_name, org_db_name)
}

#' Create and Bundle Homemade Ensembl Database
#'
#' @param species String: "human" or "mouse"
#' @param release Ensembl release version (e.g., "107")
#' @export
create_homemade_db <- function(species = "human", release = "107") {

  spec_prefix <- if(tolower(species) == "human") "Hsapiens" else "Mmusculus"
  pkg_name <- paste0("EnsDb.", spec_prefix, ".v", release)
  tar_name <- paste0(pkg_name, ".tar.gz")

  tmp_dir <- file.path(tempdir(), paste0("build_", pkg_name))
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir, recursive = TRUE)

  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  if (tolower(species) == "human") {
    org_folder     <- "homo_sapiens"
    org_scientific <- "Homo_sapiens"
    genome_ver <- if (as.numeric(release) <= 75) "GRCh37" else "GRCh38"
  } else {
    org_folder     <- "mus_musculus"
    org_scientific <- "Mus_musculus"
    genome_ver <- if (as.numeric(release) <= 102) "GRCm38" else "GRCm39"
  }

  url <- sprintf(
    "https://ftp.ensembl.org/pub/release-%s/gtf/%s/%s.%s.%s.gtf.gz",
    release, org_folder, org_scientific, genome_ver, release
  )
  gtf_path <- file.path(tmp_dir, basename(url))

  message("--- Step 1: Downloading GTF ---")
  message("Downloading from: ", url)
  options(timeout = 900)
  download.file(url, destfile = gtf_path, mode = "wb")

  message("--- Step 2: Generating SQLite Database ---")
  db_file <- ensembldb::ensDbFromGtf(
    gtf = gtf_path,
    organism = org_scientific,
    genomeVersion = genome_ver,
    version = release,
    path = tmp_dir
  )

  message("--- Step 3: Creating R Package Wrapper ---")
  ensembldb::makeEnsembldbPackage(
    ensdb = db_file,
    version = "0.0.1",
    maintainer = "User <user@work.com>",
    author = "ExpressOM Builder",
    destDir = tmp_dir,
    license = "Artistic-2.0"
  )

  message("--- Step 4: Compressing ---")
  if (!dir.exists(file.path(tmp_dir, pkg_name))) {
    stop("Expected package folder '", pkg_name, "' not found in temp directory.")
  }
  withr::with_dir(tmp_dir, {
    utils::tar(tar_name, files = pkg_name, compression = "gzip")
  })
  target_dir <- "inst/extdata"
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  file.copy(file.path(tmp_dir, tar_name), file.path(target_dir, tar_name), overwrite = TRUE)
  message("SUCCESS: Database bundled at ", file.path(target_dir, tar_name))
}

#' Install Bundled Ensembl Database
#'
#' @param pkg_name Optional string to match a specific database (e.g., "Mmusculus" or "Hsapiens")
#' @export
install_internal_db <- function(pkg_name = NULL) {
  ext_path <- system.file("extdata", package = "ExpressOM")
  if (ext_path == "") ext_path <- "inst/extdata"

  tar_files <- list.files(ext_path, pattern = "\\.tar\\.gz$", full.names = TRUE)

  if (length(tar_files) == 0) {
    stop("No .tar.gz database found in inst/extdata. Run create_homemade_db() first.")
  }

  if (!is.null(pkg_name)) {
    matched_files <- tar_files[grepl(pkg_name, basename(tar_files), ignore.case = TRUE)]
    if (length(matched_files) == 0) {
      stop("No database matching '", pkg_name, "' found. Available databases:\n",
           paste(basename(tar_files), collapse = "\n"))
    }
    db_path <- matched_files[1]
  } else {
    db_path <- tar_files[1]
    if (length(tar_files) > 1) {
      warning("Multiple databases found. Defaulting to the first one: ", basename(db_path),
              "\nUse install_internal_db(pkg_name = '...') to specify.")
    }
  }

  message("Installing bundled database: ", basename(db_path))
  remotes::install_local(db_path, upgrade = "never", build = FALSE, force = TRUE)
}

#' Extract organism specific databases and properties
#'
#' @param edb Ensembl database object
#' @export
get_organism_info <- function(edb) {
  detected_org <- ensembldb::organism(edb)

  tf_dbs <- c("ChEA_2022", "TRRUST_Transcription_Factors_2019", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

  if (grepl("Homo sapiens", detected_org, ignore.case = TRUE)) {
    return(list(
      name = detected_org,
      org_db = "org.Hs.eg.db",
      kegg_code = "hsa",
      tf_db = tf_dbs,
      msig_org = "Homo sapiens",
      msig_cat = "H",
      msig_db = "HS"
    ))
  } else if (grepl("Mus musculus", detected_org, ignore.case = TRUE)) {
    return(list(
      name = detected_org,
      org_db = "org.Mm.eg.db",
      kegg_code = "mmu",
      tf_db = tf_dbs,
      msig_org = "Mus musculus",
      msig_cat = "MH",
      msig_db = "MM"
    ))
  } else {
    stop(paste("Organism not supported:", detected_org))
  }
}

# Internal SPIA combination helpers (single authoritative copy; mod_plots.R re-uses these)
#' @keywords internal
.combfunc <- function(p1, p2, method = "fisher") {
  if (method == "fisher") {
    pchisq(-2 * (log(p1) + log(p2)), df = 4, lower.tail = FALSE)
  } else {
    pnorm((qnorm(p1) + qnorm(p2)) / sqrt(2))
  }
}

#' @keywords internal
.getP2 <- function(p, method = "fisher") {
  if (method == "fisher") {
    exp(-qchisq(p, df = 4, lower.tail = FALSE) / 2)
  } else {
    pnorm(sqrt(2) * qnorm(p), lower.tail = FALSE)
  }
}

#' Download Ensembl Reference FASTA and GTF
#'
#' @param ensembl_package_name Name of the Ensembl package (e.g., "EnsDb.Hsapiens.v107")
#' @param out_dir Directory where the reference files should be saved
#' @return A list containing the local paths: `gtf`, `cdna_fasta`, `ncrna_fasta`.
#' @export
download_ensembl_refs <- function(ensembl_package_name, out_dir = "./reference") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  matches    <- regexec("^EnsDb\\.(Hsapiens|Mmusculus)\\.v([0-9]+)$", ensembl_package_name)
  match_parts <- regmatches(ensembl_package_name, matches)[[1]]
  if (length(match_parts) != 3) stop("Could not parse ensembl_package_name. Expected format like 'EnsDb.Hsapiens.v107'")

  species <- if (match_parts[2] == "Hsapiens") "human" else "mouse"
  release <- match_parts[3]

  if (species == "human") {
    org_folder     <- "homo_sapiens"
    org_scientific <- "Homo_sapiens"
    genome_ver     <- if (as.numeric(release) <= 75) "GRCh37" else "GRCh38"
  } else {
    org_folder     <- "mus_musculus"
    org_scientific <- "Mus_musculus"
    genome_ver     <- if (as.numeric(release) <= 102) "GRCm38" else "GRCm39"
  }

  # Use https:// — http:// redirects can silently fail on some systems
  base_url  <- "https://ftp.ensembl.org/pub/release-%s"
  gtf_url   <- sprintf(paste0(base_url, "/gtf/%s/%s.%s.%s.gtf.gz"),   release, org_folder, org_scientific, genome_ver, release)
  cdna_url  <- sprintf(paste0(base_url, "/fasta/%s/cdna/%s.%s.cdna.all.fa.gz"), release, org_folder, org_scientific, genome_ver)
  ncrna_url <- sprintf(paste0(base_url, "/fasta/%s/ncrna/%s.%s.ncrna.fa.gz"),   release, org_folder, org_scientific, genome_ver)

  gtf_dest   <- file.path(out_dir, basename(gtf_url))
  cdna_dest  <- file.path(out_dir, basename(cdna_url))
  ncrna_dest <- file.path(out_dir, basename(ncrna_url))

  old_timeout <- getOption("timeout")
  options(timeout = max(1800, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  if (!file.exists(gtf_dest)) {
    message("Downloading GTF from: ", gtf_url)
    download.file(gtf_url, destfile = gtf_dest, mode = "wb")
  } else message("GTF already exists at: ", gtf_dest)

  if (!file.exists(cdna_dest)) {
    message("Downloading cDNA FASTA from: ", cdna_url)
    download.file(cdna_url, destfile = cdna_dest, mode = "wb")
  } else message("cDNA FASTA already exists at: ", cdna_dest)

  if (!file.exists(ncrna_dest)) {
    message("Downloading ncRNA FASTA from: ", ncrna_url)
    download.file(ncrna_url, destfile = ncrna_dest, mode = "wb")
  } else message("ncRNA FASTA already exists at: ", ncrna_dest)

  message("Reference downloads complete.")
  return(list(gtf = gtf_dest, cdna_fasta = cdna_dest, ncrna_fasta = ncrna_dest))
}

#' ExpressOM Pre-Flight Environment Validation
#'
#' @param run_isoform Logical; checks isoform dependencies if TRUE.
#' @param run_functional Logical; checks enrichment dependencies if TRUE.
#' @return Throws an error if dependencies are missing, returns TRUE silently otherwise.
#' @export
validate_environment <- function(run_isoform = TRUE, run_functional = TRUE) {
  message("Checking ExpressOM environment readiness...")

  core_pkgs <- c("DESeq2", "tximport", "dplyr", "ggplot2", "pheatmap")
  missing_core <- core_pkgs[!sapply(core_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_core) > 0) {
    stop(paste("Critical core packages are missing. Please install:", paste(missing_core, collapse = ", ")))
  }

  if (run_functional) {
    func_pkgs <- c("clusterProfiler", "SPIA", "fgsea", "ReactomePA", "DOSE")
    missing_func <- func_pkgs[!sapply(func_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_func) > 0) {
      stop(paste("Functional module requested, but packages are missing:", paste(missing_func, collapse = ", ")))
    }
  }

  if (run_isoform) {
    iso_pkgs <- c("DRIMSeq", "IsoformSwitchAnalyzeR", "Biostrings")
    missing_iso <- iso_pkgs[!sapply(iso_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_iso) > 0) {
      stop(paste("Isoform module requested, but packages are missing:", paste(missing_iso, collapse = ", ")))
    }
  }

  message("Environment check passed successfully!")
  return(TRUE)
}