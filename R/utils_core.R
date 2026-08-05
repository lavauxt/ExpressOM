`%||%` <- function(x, y) if (is.null(x)) y else x

#' Safely strip Ensembl-style version suffixes from an identifier vector
#' @keywords internal
#' @export
strip_ensembl_version <- function(x) {
  x <- as.character(x)

  is_versioned_ensembl <- grepl("^ENS[A-Z]*[GT][0-9]+\\.[0-9]+$", x)

  x[is_versioned_ensembl] <- sub(
    "\\.[0-9]+$",
    "",
    x[is_versioned_ensembl]
  )

  x
}

#' Clean a transcript/gene ID for exact-match purposes
#'
#' Handles:
#'   ENST00000456328.2
#'   ENST00000456328.2|ENSG00000223972.5|DDX11L1|...
#'   ENST00000456328.2 some description
#'   ENSG00000223972.5
#'
#' @keywords internal
#' @export
clean_transcript_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  # Remove everything after the first pipe.
  # Example:
  #   ENST00000456328.2|ENSG00000223972.5|DDX11L1|...
  # becomes:
  #   ENST00000456328.2
  x <- sub("\\|.*$", "", x)

  # Remove everything after the first whitespace.
  # Example:
  #   ENST00000456328.2 some description
  # becomes:
  #   ENST00000456328.2
  x <- sub("\\s+.*$", "", x)

  # Strip Ensembl version suffix.
  # Example:
  #   ENST00000456328.2 -> ENST00000456328
  x <- strip_ensembl_version(x)

  x
}
#' Build a display-safe gene label, falling back to the ID when no symbol is known
#' @keywords internal
#' @export
.coalesce_gene_label <- function(symbol, id) {
  ifelse(is.na(symbol) | symbol == "", id, symbol)
}

#' Build a "SYMBOL (id)" label without a redundant repeat when there's no symbol
#' @keywords internal
#' @export
.format_gene_tx_label <- function(gene, id) {
  ifelse(gene == id, id, paste0(gene, " (", id, ")"))
}

#' Fill missing Entrez IDs in a gene_map using clusterProfiler::bitr
#' @keywords internal
.fill_entrez_with_bitr <- function(gene_map, org_obj, id_col = "ensembl", symbol_col = "symbol") {
  if (is.null(org_obj)) return(gene_map)

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    message("  clusterProfiler not installed; skipping advanced Entrez mapping.")
    return(gene_map)
  }

  idx_na <- is.na(gene_map$entrezid) | gene_map$entrezid == ""
  if (!any(idx_na)) return(gene_map)

  message("  Attempting to fill missing Entrez IDs using clusterProfiler::bitr...")

  ens_ids <- gene_map[[id_col]][idx_na]
  ens_like <- grepl("^ENS", ens_ids)

  if (any(ens_like)) {
    ens_to_map <- unique(ens_ids[ens_like])

    map_df <- tryCatch(
      {
        clusterProfiler::bitr(
          ens_to_map,
          fromType = "ENSEMBL",
          toType = "ENTREZID",
          OrgDb = org_obj
        )
      },
      error = function(e) NULL
    )

    if (!is.null(map_df) && nrow(map_df) > 0) {
      for (i in which(idx_na)) {
        if (gene_map[[id_col]][i] %in% map_df$ENSEMBL) {
          gene_map$entrezid[i] <- map_df$ENTREZID[map_df$ENSEMBL == gene_map[[id_col]][i]][1]
        }
      }
      message("    Mapped ", nrow(map_df), " Ensembl IDs to Entrez.")
    }
  }

  idx_na2 <- is.na(gene_map$entrezid) | gene_map$entrezid == ""

  if (any(idx_na2)) {
    syms <- gene_map[[symbol_col]][idx_na2]
    syms <- syms[!is.na(syms) & syms != "" & syms != gene_map[[id_col]][idx_na2]]
    syms <- unique(syms)

    if (length(syms) > 0) {
      map_df <- tryCatch(
        {
          clusterProfiler::bitr(
            syms,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org_obj
          )
        },
        error = function(e) NULL
      )

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

  gene_map
}

#' Apply remove_sample / subset_sample filters to an imported sample table
#' @keywords internal
.apply_sample_filters <- function(sample_df, sample_col, remove_sample = NULL, subset_sample = NULL) {
  if (!is.null(remove_sample)) {
    message("   -> Excluding requested samples: ", paste(remove_sample, collapse = ", "))
    keep_indices <- !(sample_df[[sample_col]] %in% remove_sample)
    sample_df <- sample_df[keep_indices, , drop = FALSE]

    if (nrow(sample_df) == 0) {
      stop("The remove_sample constraint removed all available samples from your metadata!")
    }
  }

  if (!is.null(subset_sample)) {
    message("   -> Applying subset condition: ", subset_sample)

    sample_df <- tryCatch(
      {
        filter_expr <- rlang::parse_expr(subset_sample)
        subset_indices <- eval(filter_expr, envir = sample_df)
        sample_df[subset_indices, , drop = FALSE]
      },
      error = function(e) {
        stop("Failed to evaluate subset_sample condition. Error: ", e$message)
      }
    )

    if (nrow(sample_df) == 0) {
      stop("The subset_sample condition matched zero samples.")
    }
  }

  sample_df
}

#' Per-row z-score matrix, guarding against division by zero
#' @keywords internal
.zscore_matrix <- function(expr_mat) {
  row_means <- rowMeans(expr_mat)
  row_sds <- apply(expr_mat, 1, stats::sd)
  row_sds[row_sds == 0] <- 1
  (expr_mat - row_means) / row_sds
}

#' Run an expression while muffling one specific, known-noisy lifecycle warning
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

#' Safely relevel a condition factor to a specified base/reference level
#'
#' Centralizes the guard that used to live inline (and inconsistently) in
#' both create_dds_object() and run_dte(): skips releveling -- with a
#' message that says *why* -- instead of either (a) silently doing nothing
#' while claiming "not provided" when a value actually was supplied but
#' just didn't apply to this call, or (b) crashing with relevel()'s opaque
#' "'ref' must be an existing level" error when `base` doesn't match any
#' actual level of the factor (e.g. mismatched default level/base values
#' for this dataset's condition labels).
#'
#' @param x a factor (or character vector, which will be coerced) to relevel
#' @param base the level that should become the reference; NULL/"" skips releveling
#' @param label a short name for the factor/column, used in messages only
#' @return the releveled factor, or `x` unchanged (as a factor) if releveling was skipped
#' @keywords internal
.safe_relevel_condition <- function(x, base, label = "condition") {
  if (!is.factor(x)) x <- as.factor(x)

  if (is.null(base) || length(base) == 0 || !nzchar(as.character(base))) {
    message(
      "Note: no base level supplied for '", label,
      "'; skipping releveling (using existing factor level order: ",
      paste(levels(x), collapse = ", "), ")."
    )
    return(x)
  }

  if (!(as.character(base) %in% levels(x))) {
    message(
      "Note: base level '", base, "' not found among the levels of '", label,
      "' (found: ", paste(levels(x), collapse = ", "), "); skipping releveling."
    )
    return(x)
  }

  relevel(x, ref = as.character(base))
}

#' Safely create a directory
#' @keywords internal
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

#' Safely save a ggplot or base R plot to PDF
#' @keywords internal
safe_pdf <- function(path, expr, width = 10, height = 8) {
  expr_sub <- substitute(expr)
  caller_env <- parent.frame()

  tryCatch(
    {
      grDevices::cairo_pdf(path, width = width, height = height)
      eval(expr_sub, envir = caller_env)
      dev.off()
    },
    error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("Warning: Failed to generate plot at: ", path, "\n  Error: ", e$message)
    }
  )
}

#' Safely run an expression, returning NULL on error with an optional message
#' @keywords internal
safe_run <- function(expr, label = "") {
  tryCatch(
    expr,
    error = function(e) {
      if (nchar(label) > 0) {
        message("Warning: ", label, " failed. Skipping. Error: ", e$message)
      }
      NULL
    }
  )
}

#' Locate a bundled template/Rmd file under inst/rmd/
#'
#' More robust than the original version:
#' - checks installed package first
#' - checks common source-checkout locations
#' - walks up from the working directory
#' - honors EXPRESSOM_RMD_DIR
#' @keywords internal
.expressom_rmd_path <- function(filename) {
  env_dir <- Sys.getenv("EXPRESSOM_RMD_DIR", "")

  candidates <- c(
    if (nzchar(env_dir)) file.path(env_dir, filename),
    system.file("rmd", filename, package = "ExpressOM"),
    file.path("inst", "rmd", filename),
    file.path("rmd", filename),
    file.path("..", "inst", "rmd", filename),
    file.path("..", "..", "inst", "rmd", filename),
    file.path("..", "..", "..", "inst", "rmd", filename)
  )

  candidates <- candidates[nzchar(candidates)]

  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/"))
    }
  }

  dir <- normalizePath(getwd(), winslash = "/")

  for (i in 0:6) {
    p1 <- file.path(dir, "inst", "rmd", filename)
    if (file.exists(p1)) return(normalizePath(p1, winslash = "/"))

    p2 <- file.path(dir, "rmd", filename)
    if (file.exists(p2)) return(normalizePath(p2, winslash = "/"))

    dir <- dirname(dir)
  }

  stop(
    "Could not locate bundled template '", filename, "'.\n",
    "Run from the package root, reinstall the package, or set:\n",
    "  Sys.setenv(EXPRESSOM_RMD_DIR = '/path/to/inst/rmd')"
  )
}

#' Render a {{PLACEHOLDER}}-style template to a temp file with values substituted
#' @keywords internal
.render_placeholder_template <- function(template_file, values, fileext = ".Rmd") {
  template_path <- tryCatch(
    .expressom_rmd_path(template_file),
    error = function(e) NULL
  )

  if (is.null(template_path)) {
    warning(
      "Could not locate template '", template_file, "'. ",
      "Skipping custom template step.",
      call. = FALSE
    )
    return(NULL)
  }

  txt <- readLines(template_path, warn = FALSE)
  txt <- paste(txt, collapse = "\n")

  for (nm in names(values)) {
    txt <- gsub(
      paste0("{{", nm, "}}"),
      as.character(values[[nm]]),
      txt,
      fixed = TRUE
    )
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
#' @export
create_homemade_db <- function(species = "human",
                               release = "107",
                               maintainer = "User [user@example.com](mailto:user@example.com)",
                               author = "ExpressOM Builder") {
  spec_prefix <- if (tolower(species) == "human") "Hsapiens" else "Mmusculus"
  pkg_name <- paste0("EnsDb.", spec_prefix, ".v", release)
  tar_name <- paste0(pkg_name, ".tar.gz")

  tmp_dir <- file.path(tempdir(), paste0("build_", pkg_name))
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  if (tolower(species) == "human") {
    org_folder <- "homo_sapiens"
    org_scientific <- "Homo_sapiens"
    genome_ver <- if (as.numeric(release) <= 75) "GRCh37" else "GRCh38"
  } else {
    org_folder <- "mus_musculus"
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

  file.copy(
    file.path(tmp_dir, tar_name),
    file.path(target_dir, tar_name),
    overwrite = TRUE
  )

  message("SUCCESS: Database bundled at ", file.path(target_dir, tar_name))
}

#' Install Bundled Ensembl Database
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
      stop(
        "No database matching '", pkg_name, "' found. Available databases:\n",
        paste(basename(tar_files), collapse = "\n")
      )
    }

    db_path <- matched_files[1]
  } else {
    db_path <- tar_files[1]

    if (length(tar_files) > 1) {
      warning(
        "Multiple databases found. Defaulting to the first one: ", basename(db_path),
        "\nUse install_internal_db(pkg_name = '...') to specify."
      )
    }
  }

  message("Installing bundled database: ", basename(db_path))

  remotes::install_local(
    db_path,
    upgrade = "never",
    build = FALSE,
    force = TRUE
  )
}

#' Extract organism specific databases and properties
#' @export
get_organism_info <- function(edb) {
  detected_org <- ensembldb::organism(edb)

  tf_dbs <- c(
    "ChEA_2022",
    "TRRUST_Transcription_Factors_2019",
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
  )

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
      name = detected_org,
      org_db = "org.Rn.eg.db",
      kegg_code = "rno",
      tf_db = tf_dbs,
      msig_org = "Rattus norvegicus",
      msig_cat = "C2",
      msig_db = "RN"
    ))
  } else {
    stop(paste("Organism not supported:", detected_org))
  }
}

#' Internal SPIA combination helpers
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
#' @export
download_ensembl_refs <- function(ensembl_package_name, out_dir = "./reference") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  matches <- regexec("^EnsDb\\.(Hsapiens|Mmusculus)\\.v([0-9]+)$", ensembl_package_name)
  match_parts <- regmatches(ensembl_package_name, matches)[[1]]

  if (length(match_parts) != 3) {
    stop("Could not parse ensembl_package_name. Expected format like 'EnsDb.Hsapiens.v107'")
  }

  species <- if (match_parts[2] == "Hsapiens") "human" else "mouse"
  release <- match_parts[3]

  if (species == "human") {
    org_folder <- "homo_sapiens"
    org_scientific <- "Homo_sapiens"
    genome_ver <- if (as.numeric(release) <= 75) "GRCh37" else "GRCh38"
  } else {
    org_folder <- "mus_musculus"
    org_scientific <- "Mus_musculus"
    genome_ver <- if (as.numeric(release) <= 102) "GRCm38" else "GRCm39"
  }

  base_url <- "https://ftp.ensembl.org/pub/release-%s"

  gtf_url <- sprintf(
    paste0(base_url, "/gtf/%s/%s.%s.%s.gtf.gz"),
    release, org_folder, org_scientific, genome_ver, release
  )

  cdna_url <- sprintf(
    paste0(base_url, "/fasta/%s/cdna/%s.%s.cdna.all.fa.gz"),
    release, org_folder, org_scientific, genome_ver
  )

  ncrna_url <- sprintf(
    paste0(base_url, "/fasta/%s/ncrna/%s.%s.ncrna.fa.gz"),
    release, org_folder, org_scientific, genome_ver
  )

  gtf_dest <- file.path(out_dir, basename(gtf_url))
  cdna_dest <- file.path(out_dir, basename(cdna_url))
  ncrna_dest <- file.path(out_dir, basename(ncrna_url))

  old_timeout <- getOption("timeout", 60L)
  options(timeout = max(1800L, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  if (!file.exists(gtf_dest)) {
    message("Downloading GTF from: ", gtf_url)
    download.file(gtf_url, destfile = gtf_dest, mode = "wb")
  } else {
    message("GTF already exists at: ", gtf_dest)
  }

  if (!file.exists(cdna_dest)) {
    message("Downloading cDNA FASTA from: ", cdna_url)
    download.file(cdna_url, destfile = cdna_dest, mode = "wb")
  } else {
    message("cDNA FASTA already exists at: ", cdna_dest)
  }

  if (!file.exists(ncrna_dest)) {
    message("Downloading ncRNA FASTA from: ", ncrna_url)
    download.file(ncrna_url, destfile = ncrna_dest, mode = "wb")
  } else {
    message("ncRNA FASTA already exists at: ", ncrna_dest)
  }

  message("Reference downloads complete.")

  return(list(
    gtf = gtf_dest,
    cdna_fasta = cdna_dest,
    ncrna_fasta = ncrna_dest
  ))
}

#' ExpressOM Pre-Flight Environment Validation
#' @export
validate_environment <- function(run_isoform = TRUE, run_functional = TRUE) {
  message("Checking ExpressOM environment readiness...")

  core_pkgs <- c("DESeq2", "tximport", "dplyr", "ggplot2", "pheatmap")
  missing_core <- core_pkgs[!sapply(core_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_core) > 0) {
    stop(paste(
      "Critical core packages are missing. Please install:",
      paste(missing_core, collapse = ", ")
    ))
  }

  if (run_functional) {
    func_pkgs <- c("clusterProfiler", "SPIA", "fgsea", "ReactomePA", "DOSE")
    missing_func <- func_pkgs[!sapply(func_pkgs, requireNamespace, quietly = TRUE)]

    if (length(missing_func) > 0) {
      warning(paste(
        "Functional module requested, but packages are missing:",
        paste(missing_func, collapse = ", "),
        "\nInstall with: BiocManager::install(c(",
        paste0("'", missing_func, "'", collapse = ", "),
        "))"
      ))
    }
  }

  if (run_isoform) {
    iso_pkgs <- c("DRIMSeq", "IsoformSwitchAnalyzeR", "Biostrings")
    missing_iso <- iso_pkgs[!sapply(iso_pkgs, requireNamespace, quietly = TRUE)]

    if (length(missing_iso) > 0) {
      warning(paste(
        "Isoform module requested, but packages are missing:",
        paste(missing_iso, collapse = ", "),
        "\nInstall with: BiocManager::install(c(",
        paste0("'", missing_iso, "'", collapse = ", "),
        "))"
      ))
    }
  }

  if (!requireNamespace("pdftools", quietly = TRUE)) {
    warning(
      "Package 'pdftools' is not installed. ",
      "HTML reports (RegionReport, DTU/DTE report) will not be able to ",
      "embed PDF plots as PNG, causing missing figures. ",
      "Install with: install.packages('pdftools')"
    )
  }

  message("Environment check passed successfully!")

  return(TRUE)
}

#' Cap plot dimensions to avoid device errors
#' @keywords internal
.cap_plot_dims <- function(width, height, max_dim = 30, min_dim = 3) {
  width <- as.numeric(width)
  height <- as.numeric(height)

  if (length(width) == 0 || is.na(width)) width <- 8
  if (length(height) == 0 || is.na(height)) height <- 6

  width <- max(min_dim, min(width, max_dim))
  height <- max(min_dim, min(height, max_dim))

  list(width = width, height = height)
}

#' Safely save a ggplot to file, capping runaway width/height and defaulting
#' to a Cairo-backed PDF device so text renders correctly
#'
#' ggsave()'s default "pdf" device only *references* the base-14 PDF fonts
#' by name instead of embedding glyphs; as soon as a plot uses a character
#' outside those fonts' basic set, later PDF -> PNG rasterization (e.g. for
#' report embedding via convert_pdf_to_png()) fails with "No display font
#' for 'Symbol'" / "'ArialUnicode'" from Poppler/Ghostscript, and the plot
#' silently never appears in the report. cairo_pdf() embeds real glyphs via
#' Cairo/fontconfig, which avoids that class of error entirely.
#' @keywords internal
.safe_ggsave <- function(filename,
                         plot,
                         width,
                         height,
                         device = grDevices::cairo_pdf,
                         ...) {

  if (is.null(plot)) {
    message("   -> Skipping plot save, plot object is NULL: ", basename(filename))
    return(invisible(FALSE))
  }

  dims <- .cap_plot_dims(width, height, max_dim = 30, min_dim = 4)

  tryCatch(
    {
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        device = device,
        width = dims$width,
        height = dims$height,
        ...
      )

      invisible(TRUE)
    },
    error = function(e) {
      message(
        "   -> Failed to save plot: ", basename(filename),
        "\n      Error: ", conditionMessage(e)
      )
      invisible(FALSE)
    }
  )
}

#' Standardize a gene annotation table so downstream code can rely on
#' columns named "ensembl", "symbol", and "entrezid".
#'
#' @keywords internal
#' @export
.standardize_gene_map <- function(gene_map) {
  if (is.null(gene_map)) {
    return(
      data.frame(
        ensembl = character(0),
        symbol = character(0),
        entrezid = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  gene_map <- as.data.frame(gene_map, stringsAsFactors = FALSE)

  if (ncol(gene_map) == 0) {
    stop("gene_map is empty.")
  }

  names(gene_map) <- make.unique(trimws(names(gene_map)))
  if (!"ensembl" %in% names(gene_map)) {
    id_candidates <- c(
      "gene_id",
      "ensembl_id",
      "ensembl_gene_id",
      "ensembl_gene",
      "gene",
      "id",
      "GeneID",
      "gene_id.1"
    )

    id_col <- intersect(id_candidates, names(gene_map))

    if (length(id_col) == 0) {
      stop(
        "gene_map does not contain an Ensembl/gene ID column.\n",
        "Expected one of: ensembl, gene_id, ensembl_id, ensembl_gene_id.\n",
        "Found columns: ",
        paste(names(gene_map), collapse = ", ")
      )
    }

    names(gene_map)[names(gene_map) == id_col[1]] <- "ensembl"
  }

  extra_ensembl_cols <- grep("^ensembl\\.", names(gene_map), value = TRUE)

  if (length(extra_ensembl_cols) > 0) {
    gene_map <- gene_map[, setdiff(names(gene_map), extra_ensembl_cols), drop = FALSE]
  }

  gene_map$ensembl <- strip_ensembl_version(as.character(gene_map$ensembl))

  if (!"symbol" %in% names(gene_map)) {
    symbol_candidates <- c(
      "gene_name",
      "symbol",
      "external_gene_name",
      "gene_symbol",
      "SYMBOL",
      "GeneSymbol",
      "name"
    )

    symbol_col <- intersect(symbol_candidates, names(gene_map))

    if (length(symbol_col) > 0) {
      gene_map$symbol <- as.character(gene_map[[symbol_col[1]]])
    } else {
      gene_map$symbol <- gene_map$ensembl
    }
  }

  gene_map$symbol <- as.character(gene_map$symbol)

  gene_map$symbol[is.na(gene_map$symbol) | gene_map$symbol == ""] <-
    gene_map$ensembl[is.na(gene_map$symbol) | gene_map$symbol == ""]

  if (!"entrezid" %in% names(gene_map)) {
    entrez_candidates <- c(
      "entrezid",
      "ENTREZID",
      "entrez_id",
      "entrezgene",
      "entrez_gene"
    )

    entrez_col <- intersect(entrez_candidates, names(gene_map))

    if (length(entrez_col) > 0) {
      gene_map$entrezid <- as.character(gene_map[[entrez_col[1]]])
    } else {
      gene_map$entrezid <- NA_character_
    }
  }

  gene_map$entrezid <- as.character(gene_map$entrezid)

  gene_map <- gene_map[!is.na(gene_map$ensembl) & gene_map$ensembl != "", , drop = FALSE]
  gene_map <- gene_map[!duplicated(gene_map$ensembl), , drop = FALSE]

  gene_map
}