#' Safely create a directory
#' @keywords internal
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

#' Safely save a ggplot or base R plot to PDF
#' @keywords internal
safe_pdf <- function(path, expr, width = 10, height = 8) {
  tryCatch({
    pdf(path, width = width, height = height)
    force(expr)
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


#' Create and Bundle Homemade Ensembl Database
#' @param species String: "human" or "mouse"
#' @param release Ensembl release version (e.g., "107")
#' @export
create_homemade_db <- function(species = "human", release = "107") {
  
  # 1. Explicitly define the Package Name first
  # Example: EnsDb.Hsapiens.v107
  spec_prefix <- if(tolower(species) == "human") "Hsapiens" else "Mmusculus"
  pkg_name <- paste0("EnsDb.", spec_prefix, ".v", release)
  tar_name <- paste0(pkg_name, ".tar.gz")
  
  # 2. Setup Sandbox Directory
  tmp_dir <- file.path(tempdir(), paste0("build_", pkg_name))
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir, recursive = TRUE)
  
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  
  # 3. Setup Variables
  if (tolower(species) == "human") {
    org_folder <- "homo_sapiens"
    org_scientific <- "Homo_sapiens"
    genome_ver <- "GRCh38"
  } else {
    org_folder <- "mus_musculus"
    org_scientific <- "Mus_musculus"
    genome_ver <- "GRCm39"
  }

  url <- sprintf(
    "https://ftp.ensembl.org/pub/release-%s/gtf/%s/%s.%s.%s.gtf.gz",
    release, org_folder, org_scientific, genome_ver, release
  )
  gtf_path <- file.path(tmp_dir, basename(url))
  
  # 4. Download
  message("--- Step 1: Downloading GTF ---")
  options(timeout = 900) 
  download.file(url, destfile = gtf_path, mode = "wb")

  # 5. Generate SQLite DB
  message("--- Step 2: Generating SQLite Database ---")
  db_file <- ensembldb::ensDbFromGtf(
    gtf = gtf_path,
    organism = org_scientific,
    genomeVersion = genome_ver,
    version = release,
    path = tmp_dir
  )
  
  # 6. Create Package Wrapper
  message("--- Step 3: Creating R Package Wrapper ---")
  # We don't save the return value to pkg_path anymore to avoid the "TRUE" bug
  ensembldb::makeEnsembldbPackage(
    ensdb = db_file,
    version = "0.0.1",
    maintainer = "User <user@work.com>",
    author = "ExpressOM Builder",
    destDir = tmp_dir,
    license = "Artistic-2.0"
  )
  
  # 7. Compress
  message("--- Step 4: Compressing ---")
  old_wd <- getwd()
  # Use on.exit to ensure we always return to the original working directory
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  
  # Ensure the folder actually exists before tarring
  if (!dir.exists(pkg_name)) {
    stop("Expected package folder '", pkg_name, "' not found in temp directory.")
  }
  
  utils::tar(tar_name, files = pkg_name, compression = "gzip")
  
  # 8. Move to inst/extdata
  setwd(old_wd) # Move back to project root to find inst/
  target_dir <- "inst/extdata"
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  
  file.copy(file.path(tmp_dir, tar_name), file.path(target_dir, tar_name), overwrite = TRUE)
  
  message("SUCCESS: Database bundled at ", file.path(target_dir, tar_name))
}

#' Install Bundled Ensembl Database
#' @param pkg_name Optional string to match a specific database (e.g., "Mmusculus" or "Hsapiens")
#' @export
install_internal_db <- function(pkg_name = NULL) {
  # Look for any .tar.gz in the bundled data folder
  ext_path <- system.file("extdata", package = "ExpressOM")
  if (ext_path == "") ext_path <- "inst/extdata"
  
  tar_files <- list.files(ext_path, pattern = "\\.tar\\.gz$", full.names = TRUE)
  
  if (length(tar_files) == 0) {
    stop("No .tar.gz database found in inst/extdata. Run create_homemade_db() first.")
  }
  
  # Filter by pkg_name if provided
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
#' @param edb Ensembl database object
#' @export
get_organism_info <- function(edb) {
  detected_org <- ensembldb::organism(edb)
  
  # Define robust fallback lists for TFs in case one is temporarily missing from the EnrichR API
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