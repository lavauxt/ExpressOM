`%||%` <- function(x, y) if (is.null(x)) y else x

#' Safely strip Ensembl-style version suffixes from an identifier vector
#'
#' Unlike \code{sub("\\..*$", "", x)}, which truncates at the FIRST dot and
#' therefore destroys any non-Ensembl identifier that legitimately contains a
#' dot (e.g. StringTie loci like "MSTRG.116" or transcripts like "MSTRG.6.1",
#' both of which would otherwise collapse to the bare string "MSTRG"), this
#' only strips a trailing ".<digits>" version suffix when the identifier
#' actually looks like an Ensembl gene/transcript ID (ENSG/ENST, or their
#' species-prefixed variants like ENSMUSG/ENSMUST). Every other identifier
#' (MSTRG IDs, gene symbols, custom IDs) is returned unchanged.
#'
#' @param x Character vector of identifiers
#' @return Character vector, same length as \code{x}
#' @keywords internal
#' @export
strip_ensembl_version <- function(x) {
  is_versioned_ensembl <- grepl("^ENS[A-Z]*[GT][0-9]+\\.[0-9]+$", x)
  x[is_versioned_ensembl] <- sub("\\.[0-9]+$", "", x[is_versioned_ensembl])
  x
}

#' Clean a transcript/gene ID for exact-match purposes
#'
#' \code{strip_ensembl_version()}'s pattern is anchored (\code{^...$}), so it
#' only recognizes a version suffix on a BARE id and is a silent no-op on any
#' compound header where extra fields remain attached after the version
#' number -- e.g. a GENCODE-style pipe-delimited FASTA header
#' ("ENST00000456328.2|ENSG00000223972.5|..."), or a SQANTI3-rescued
#' reference mixing reference-matched and novel transcript IDs. Bar- and
#' space-delimited description text has to be truncated FIRST, then
#' \code{strip_ensembl_version()} runs on what's left. Doing it in the
#' opposite order silently fails to strip the version at all (see the
#' identical fix and write-up in \code{run_isoform_switch()}'s local
#' \code{clean_id()}, mod_isoform.R). This is the single source of truth for
#' that ordering -- used anywhere a transcript ID coming from a
#' quantification file (kallisto/salmon target_id, tx2gene, custom ID maps)
#' needs to be normalized before matching.
#'
#' @param x Character vector of identifiers
#' @return Character vector, same length as \code{x}
#' @keywords internal
#' @export
clean_transcript_id <- function(x) {
  x <- sub("\\|.*$", "", x)
  x <- sub(" .*$",   "", x)
  strip_ensembl_version(x)
}

#' Fill missing Entrez IDs in a gene_map using clusterProfiler::bitr
#'
#' Shared by both the DGE (\code{import_counts()}) and isoform
#' (\code{import_transcript_counts()}) import paths. Previously this was
#' defined as two separately-maintained, drifting copies (one in mod_dge.R,
#' one in mod_isoform.R marked "copied from mod_dge.R") -- consolidated here
#' as the single source of truth.
#'
#' @param gene_map Data frame with at least \code{id_col}, \code{symbol_col},
#'   and \code{entrezid} columns
#' @param org_obj Loaded OrgDb object (e.g. org.Hs.eg.db), or NULL to skip
#' @param id_col Name of the ID column to try mapping from Ensembl (default "ensembl")
#' @param symbol_col Name of the gene symbol column (default "symbol")
#' @return \code{gene_map} with as many \code{entrezid} NAs filled in as possible
#' @keywords internal
.fill_entrez_with_bitr <- function(gene_map, org_obj, id_col = "ensembl", symbol_col = "symbol") {
  if (is.null(org_obj)) return(gene_map)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    message("  clusterProfiler not installed; skipping advanced Entrez mapping.")
    return(gene_map)
  }

  # Identify rows with missing Entrez
  idx_na <- is.na(gene_map$entrezid) | gene_map$entrezid == ""
  if (!any(idx_na)) return(gene_map)

  message("  Attempting to fill missing Entrez IDs using clusterProfiler::bitr...")

  # 1) Try mapping by Ensembl ID (if IDs look like Ensembl)
  ens_ids <- gene_map[[id_col]][idx_na]
  ens_like <- grepl("^ENS", ens_ids)
  if (any(ens_like)) {
    ens_to_map <- unique(ens_ids[ens_like])
    map_df <- tryCatch({
      clusterProfiler::bitr(ens_to_map, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org_obj)
    }, error = function(e) NULL)
    if (!is.null(map_df) && nrow(map_df) > 0) {
      # Map back to gene_map
      for (i in which(idx_na)) {
        if (gene_map[[id_col]][i] %in% map_df$ENSEMBL) {
          gene_map$entrezid[i] <- map_df$ENTREZID[map_df$ENSEMBL == gene_map[[id_col]][i]][1]
        }
      }
      message("    Mapped ", nrow(map_df), " Ensembl IDs to Entrez.")
    }
  }

  # 2) For remaining NAs, try mapping by gene symbol (if symbol is available)
  idx_na2 <- is.na(gene_map$entrezid) | gene_map$entrezid == ""
  if (any(idx_na2)) {
    syms <- gene_map[[symbol_col]][idx_na2]
    # Exclude symbols that are NA, empty, or identical to the ID (likely custom IDs)
    syms <- syms[!is.na(syms) & syms != "" & syms != gene_map[[id_col]][idx_na2]]
    syms <- unique(syms)
    if (length(syms) > 0) {
      map_df <- tryCatch({
        clusterProfiler::bitr(syms, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_obj)
      }, error = function(e) NULL)
      if (!is.null(map_df) && nrow(map_df) > 0) {
        for (i in which(idx_na2)) {
          sym_i <- gene_map[[symbol_col]][i]
          if (sym_i %in% map_df$SYMBOL) {
            gene_map$entrezid[i] <- map_df$ENTREZID[map_df$SYMBOL == sym_i][1]
          }
        }
        message("    Mapped ", nrow(map_df), " symbols to Entrez.")
      }
    }
  }

  # Return updated gene_map
  gene_map
}

#' Apply remove_sample / subset_sample filters to an imported sample table
#'
#' Single source of truth for the "drop explicitly excluded samples, then
#' apply an optional filter expression" logic used by both import_counts()
#' (gene-level, mod_dge.R) and import_transcript_counts() (isoform-level,
#' mod_isoform.R). These were previously two independently-maintained copies
#' that had already drifted apart: the isoform copy silently dropped both the
#' progress messages and the tryCatch() around subset_sample, so a typo'd
#' subset_sample expression there raised R's raw eval() error instead of the
#' friendlier "Failed to evaluate subset_sample condition" message the DGE
#' path gives. Neither copy guarded against subset_sample matching zero rows
#' (only remove_sample did); that gap is closed here for both callers.
#'
#' @param sample_df Data frame read from the sample table (before rownames are set)
#' @param sample_col Name of the column holding sample IDs ("Sample" or "sample_id")
#' @param remove_sample Optional character vector of sample IDs to drop
#' @param subset_sample Optional string, evaluated as an expression in the
#'   data frame's environment (e.g. "cell_type == 'T_cells'")
#' @return The filtered data frame
#' @keywords internal
.apply_sample_filters <- function(sample_df, sample_col, remove_sample = NULL, subset_sample = NULL) {
  if (!is.null(remove_sample)) {
    message("   -> Excluding requested samples: ", paste(remove_sample, collapse = ", "))
    keep_indices <- !(sample_df[[sample_col]] %in% remove_sample)
    sample_df    <- sample_df[keep_indices, , drop = FALSE]
    if (nrow(sample_df) == 0) stop("The remove_sample constraint removed all available samples from your metadata!")
  }

  if (!is.null(subset_sample)) {
    message("   -> Applying subset condition: ", subset_sample)
    sample_df <- tryCatch({
      filter_expr    <- rlang::parse_expr(subset_sample)
      subset_indices <- eval(filter_expr, envir = sample_df)
      sample_df[subset_indices, , drop = FALSE]
    }, error = function(e) {
      stop("Failed to evaluate subset_sample condition. Error: ", e$message)
    })
    if (nrow(sample_df) == 0) stop("The subset_sample condition matched zero samples.")
  }

  sample_df
}

#' Per-row z-score matrix, guarding against division by zero
#'
#' Single source of truth for the "z = (x - row_mean) / row_sd, with a
#' zero-variance row held at z = 0 instead of NaN" computation shared by
#' plot_sample_zscore(), plot_geneset_zscore_avg(), and
#' plot_gene_zscore_individual() in mod_plots.R (previously three separately
#' hand-copied instances of the same four lines).
#'
#' @param expr_mat Numeric matrix, genes/features as rows, samples as columns
#' @return Matrix of the same dimensions, row-wise z-scored
#' @keywords internal
.zscore_matrix <- function(expr_mat) {
  row_means             <- rowMeans(expr_mat)
  row_sds               <- apply(expr_mat, 1, stats::sd)
  row_sds[row_sds == 0] <- 1
  (expr_mat - row_means) / row_sds
}

#' Run an expression while muffling one specific, known-noisy lifecycle warning
#'
#' IsoformSwitchAnalyzeR's own internal code (there is no dplyr::across() call
#' anywhere in ExpressOM itself) uses the deprecated
#' \code{dplyr::filter(dplyr::across())} pattern without supplying
#' \code{.cols}, which since dplyr 1.1.0 emits one lifecycle-deprecation
#' warning per row it's evaluated on. On a genome-wide isoform/exon
#' annotation, a single \code{importRdata()} or
#' \code{analyzeAlternativeSplicing()} call can therefore raise hundreds of
#' thousands of warnings -- collecting, formatting, and (interactively)
#' printing that many condition objects is not free, and it buries any other,
#' genuinely actionable warning in noise. This muffles ONLY warnings matching
#' that one message pattern; every other warning (including ExpressOM's own)
#' still propagates normally, so it's safe to wrap broadly.
#'
#' @param expr Expression to evaluate
#' @keywords internal
.muffle_across_deprecation <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("across\\(\\).*\\.cols", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

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

#' Locate a bundled template/Rmd file under inst/rmd/
#'
#' Thin wrapper around \code{system.file()} with a development-mode fallback,
#' so templates resolve correctly both for an installed package and under
#' \code{devtools::load_all()} (which normally shims \code{system.file()} to
#' find \code{inst/} already, but this keeps things working even if invoked
#' in a context where that shim isn't active -- the same class of bug noted
#' for the SPARKS report templates: "not bundled post-install").
#'
#' @param filename Filename under \code{inst/rmd/} (e.g. \code{"dte_dtu_report.Rmd"})
#' @return Absolute path to the template file
#' @keywords internal
.expressom_rmd_path <- function(filename) {
  p <- system.file("rmd", filename, package = "ExpressOM")
  if (!nzchar(p)) {
    p <- file.path("inst", "rmd", filename)
  }
  if (!file.exists(p)) {
    stop("Could not locate bundled template 'inst/rmd/", filename, "'. ",
         "If running from a development checkout, make sure the working ",
         "directory is the package root, or reinstall the package so ",
         "inst/rmd/ is bundled.")
  }
  p
}

#' Render a {{PLACEHOLDER}}-style template to a temp file with values substituted
#'
#' Used for small Rmd fragments (e.g. regionReport's \code{customCode}
#' argument) that need simple string substitution rather than a full
#' \code{rmarkdown::render(params = ...)} pass.
#'
#' @param template_file Filename under \code{inst/rmd/}
#' @param values Named list/character vector; each name's \code{{{NAME}}}
#'   token in the template is replaced with its (character-coerced) value
#' @param fileext File extension for the returned temp file (default ".Rmd")
#' @return Path to a temp file containing the substituted content. Caller is
#'   responsible for \code{unlink()}-ing it when done.
#' @keywords internal
.render_placeholder_template <- function(template_file, values, fileext = ".Rmd") {
  txt <- readLines(.expressom_rmd_path(template_file), warn = FALSE)
  txt <- paste(txt, collapse = "\n")
  for (nm in names(values)) {
    txt <- gsub(paste0("{{", nm, "}}"), as.character(values[[nm]]), txt, fixed = TRUE)
  }
  out <- tempfile(fileext = fileext)
  writeLines(txt, out)
  out
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
create_homemade_db <- function(species = "human", release = "107",
                               maintainer = "User <user@example.com>",
                               author = "ExpressOM Builder") {

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
    maintainer = maintainer,
    author = author,
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
  } else if (grepl("Rattus norvegicus", detected_org, ignore.case = TRUE)) {
    return(list(
      name      = detected_org,
      org_db    = "org.Rn.eg.db",
      kegg_code = "rno",
      tf_db     = tf_dbs,
      msig_org  = "Rattus norvegicus",
      msig_cat  = "C2",
      msig_db   = "RN"
    ))
  } else {
    stop(paste("Organism not supported:", detected_org))
  }
}

# Internal SPIA combination helpers (single authoritative copy; mod_plots.R re-uses these)
#' @keywords internal
.combfunc <- function(p1, p2, method = "fisher") {
  if (method == "fisher") {
    p1 <- pmax(p1, .Machine$double.xmin)
    p2 <- pmax(p2, .Machine$double.xmin)
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

  old_timeout <- getOption("timeout", 60L)
  options(timeout = max(1800L, old_timeout))
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
      warning(paste("Functional module requested, but packages are missing:", paste(missing_func, collapse = ", "),
                    "\nInstall with: BiocManager::install(c(", paste0('"', missing_func, '"', collapse=", "), "))"))
    }
  }

  if (run_isoform) {
    iso_pkgs <- c("DRIMSeq", "IsoformSwitchAnalyzeR", "Biostrings")
    missing_iso <- iso_pkgs[!sapply(iso_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_iso) > 0) {
      warning(paste("Isoform module requested, but packages are missing:", paste(missing_iso, collapse = ", "),
                    "\nInstall with: BiocManager::install(c(", paste0('"', missing_iso, '"', collapse=", "), "))"))
    }
  }

  message("Environment check passed successfully!")
  return(TRUE)
}