# mod_isoform_utils.R -- small shared utilities for isoform report/plot generation
#

#' Convert a PDF plot to a PNG alongside it
#' @keywords internal
convert_pdf_to_png <- function(pdf_file, dpi = 200) {
  if (!file.exists(pdf_file)) {
    warning("PDF not found: ", pdf_file)
    return(NULL)
  }

  sz <- file.info(pdf_file)$size

  if (is.na(sz) || sz < 2000) {
    warning("PDF is empty or too small (<2KB), skipping conversion: ", pdf_file)
    return(NULL)
  }

  png_file <- sub("\\.pdf$", ".png", pdf_file)

  if (requireNamespace("pdftools", quietly = TRUE)) {
    produced <- tryCatch(
      {
        pdftools::pdf_convert(
          pdf = pdf_file,
          filenames = png_file,
          dpi = dpi,
          verbose = FALSE
        )
      },
      error = function(e) {
        warning("pdftools failed to convert ", pdf_file, " to PNG: ", conditionMessage(e))
        character(0)
      }
    )

    existing <- produced[file.exists(produced)]

    if (length(existing) > 0) {
      if (length(existing) > 1) {
        message(
          "   -> ", basename(pdf_file), " has ", length(existing),
          " pages; using page 1 (", basename(existing[1]), ") for report embedding."
        )
      }

      return(normalizePath(existing[1], winslash = "/", mustWork = FALSE))
    }
  }

  if (requireNamespace("magick", quietly = TRUE)) {
    ok <- tryCatch(
      {
        img <- magick::image_read(
          pdf_file,
          density = sprintf("%dx%d", dpi, dpi)
        )

        magick::image_write(
          image = img,
          path = png_file,
          format = "png"
        )

        file.exists(png_file)
      },
      error = function(e) {
        warning("magick failed to convert ", pdf_file, " to PNG: ", conditionMessage(e))
        FALSE
      }
    )

    if (isTRUE(ok) && file.exists(png_file)) {
      return(normalizePath(png_file, winslash = "/", mustWork = FALSE))
    }
  }

  warning(
    "Could not convert ", pdf_file, " to PNG. ",
    "Install 'pdftools' (recommended) or 'magick' plus Ghostscript. ",
    "The PDF still exists and can be opened directly."
  )

  NULL
}

# mod_isoform_import.R -- loading transcript-level counts for isoform analysis
#

#' Import transcript-level counts for isoform analysis
#' @export
import_transcript_counts <- function(
  data_dir,
  sample_table,
  ensembl_package_name,
  count_type = "salmon",
  matrix_file = NULL,
  subset_sample = NULL,
  remove_sample = NULL,
  custom_tx2gene = NULL,
  custom_gene_map = NULL
) {
  custom_tx2gene <- .validate_optional_file(custom_tx2gene, "custom_tx2gene")
  custom_gene_map <- .validate_optional_file(custom_gene_map, "custom_gene_map")

  if (!file.exists(sample_table)) {
    stop("Sample table not found: ", sample_table)
  }

  sample_df <- data.table::fread(
    sample_table,
    header = TRUE,
    data.table = FALSE
  )
  sample_col <- .sample_id_column(sample_df)

  sample_df <- .apply_sample_filters(
    sample_df,
    sample_col,
    remove_sample,
    subset_sample
  )
  rownames(sample_df) <- sample_df[[sample_col]]

  edb <- getExportedValue(ensembl_package_name, ensembl_package_name)

  if (!is.null(custom_tx2gene) && file.exists(custom_tx2gene)) {
    message("Using custom tx2gene file: ", custom_tx2gene)

    tx2gene <- data.table::fread(
      custom_tx2gene,
      header = TRUE,
      data.table = FALSE
    )

    if (!all(c("tx_id", "gene_id") %in% colnames(tx2gene))) {
      stop("Custom tx2gene must contain columns 'tx_id' and 'gene_id'")
    }

    tx2gene <- tx2gene[, c("tx_id", "gene_id")]
  } else {
    tx2gene <- ensembldb::transcripts(
      edb,
      columns = c("tx_id", "gene_id"),
      return.type = "DataFrame"
    )

    tx2gene <- as.data.frame(tx2gene)
  }

  tx2gene <- .normalize_tx2gene(tx2gene)

  org_info <- get_organism_info(edb)
  org_db <- org_info$org_db
  org_obj <- if (requireNamespace(org_db, quietly = TRUE)) {
    .load_org_db(org_db)
  } else {
    NULL
  }

  if (!is.null(custom_gene_map) && file.exists(custom_gene_map)) {
    message("Using custom gene annotation file: ", custom_gene_map)

    gene_map <- data.table::fread(
      custom_gene_map,
      header = TRUE,
      data.table = FALSE
    )

    if (
      !("gene_id" %in% colnames(gene_map)) && "ensembl" %in% colnames(gene_map)
    ) {
      colnames(gene_map)[colnames(gene_map) == "ensembl"] <- "gene_id"
    }

    if (!all(c("gene_id", "symbol") %in% colnames(gene_map))) {
      stop(
        "Custom gene map must contain columns 'gene_id' (or 'ensembl') and 'symbol'"
      )
    }

    gene_map$gene_id <- strip_ensembl_version(gene_map$gene_id)
    colnames(gene_map)[colnames(gene_map) == "gene_id"] <- "ensembl"

    if (!"entrezid" %in% colnames(gene_map)) {
      gene_map$entrezid <- NA_character_
    }

    gene_map <- gene_map[, c("ensembl", "symbol", "entrezid")]
    gene_map <- gene_map[!duplicated(gene_map$ensembl), ]

    gene_map$symbol[is.na(gene_map$symbol) | gene_map$symbol == ""] <-
      gene_map$ensembl[is.na(gene_map$symbol) | gene_map$symbol == ""]

    if (!is.null(org_obj)) {
      gene_map <- .fill_entrez_with_bitr(
        gene_map,
        org_obj,
        id_col = "ensembl",
        symbol_col = "symbol"
      )
    }

    entrez_present <- sum(!is.na(gene_map$entrezid) & gene_map$entrezid != "")
    message(
      "  Gene map loaded: ",
      nrow(gene_map),
      " genes, ",
      entrez_present,
      " with Entrez IDs"
    )
  } else {
    gene_map <- ensembldb::genes(
      edb,
      columns = c("gene_id", "gene_name"),
      return.type = "DataFrame"
    )

    gene_map <- as.data.frame(gene_map)
    colnames(gene_map) <- c("ensembl", "symbol")
    gene_map$ensembl <- strip_ensembl_version(gene_map$ensembl)

    if (!is.null(org_obj)) {
      mapped_entrez <- suppressMessages(
        AnnotationDbi::mapIds(
          org_obj,
          keys = gene_map$ensembl,
          column = "ENTREZID",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      )

      gene_map$entrezid <- as.character(mapped_entrez)
    } else {
      gene_map$entrezid <- NA_character_
    }
  }

  if (count_type != "matrix") {
    count_file_name <- switch(
      count_type,
      "salmon" = "quant.sf",
      "kallisto" = "abundance.tsv",
      # NOTE: this must be the *isoform*-level RSEM file, not
      # "quant.genes.results" (the gene-level file used by import_counts()
      # in mod_dge.R). tximport::tximport(type = "rsem") auto-detects
      # gene-level input by checking whether the filename contains "genes"
      # and, if so, force-switches to its gene-level reader regardless of
      # txOut = TRUE below -- so pointing this at the genes.results file
      # silently returns gene-level counts mislabeled as transcript-level
      # data (no error; see the mapping check right after tximport() for
      # the safety net this relies on if that ever regresses).
      "rsem" = "quant.isoforms.results",
      "stringtie" = "t_data.ctab",
      stop("Unsupported count_type for tximport: ", count_type)
    )

    file_list <- .resolve_quantification_files(
      data_dir,
      sample_df[[sample_col]],
      count_type,
      count_file_name
    )

    txi <- tximport::tximport(
      file_list,
      type = count_type,
      txOut = TRUE,
      countsFromAbundance = "lengthScaledTPM"
    )

    rownames(txi$counts) <- clean_transcript_id(rownames(txi$counts))
    if (anyDuplicated(rownames(txi$counts))) {
      stop(
        "Imported feature IDs are duplicated after transcript-ID normalization. ",
        "Check the quantification IDs for version/header variants that collapse ",
        "to the same transcript ID.",
        call. = FALSE
      )
    }

    # Sanity check: rownames(txi$counts) should be transcript IDs that
    # overlap tx2gene$tx_id (already version-stripped above). If a
    # count_type's quantification file turns out to actually be gene-level
    # (as silently happens for "rsem" if a genes.results file is passed
    # in), or count_type/data_dir otherwise point at the wrong files, the
    # ID namespaces won't match and every downstream annotation join fails
    # silently (DTE) or with a hard-to-read error much later in the
    # pipeline (DTU/switch analysis). Catch it here, once, with a message
    # that points straight at the cause.
    .validate_feature_overlap(
      rownames(txi$counts),
      tx2gene,
      paste0("imported '", count_type, "'")
    )

    meta <- sample_df[colnames(txi$counts), , drop = FALSE]

    return(list(
      txi = txi,
      meta = meta,
      tx2gene = tx2gene,
      gene_map = gene_map,
      type = "tximport"
    ))
  } else {
    if (is.null(matrix_file)) {
      stop("matrix_file required for count_type='matrix'")
    }

    counts_df <- data.table::fread(matrix_file, data.table = FALSE)
    rownames(counts_df) <- counts_df[, 1]
    counts_df <- counts_df[, -1, drop = FALSE]

    valid_samples <- .match_matrix_samples(
      colnames(counts_df),
      rownames(sample_df)
    )
    count_mat <- as.matrix(counts_df[, valid_samples, drop = FALSE])
    suppressWarnings(mode(count_mat) <- "numeric")
    count_mat <- .validate_raw_count_matrix(count_mat)

    rownames(count_mat) <- clean_transcript_id(rownames(count_mat))
    if (anyDuplicated(rownames(count_mat))) {
      stop(
        "Matrix feature IDs are duplicated after transcript-ID normalization. ",
        "Check for version or header variants that collapse to the same ID.",
        call. = FALSE
      )
    }
    .validate_feature_overlap(rownames(count_mat), tx2gene, "matrix")
    meta <- sample_df[valid_samples, , drop = FALSE]

    return(list(
      counts = count_mat,
      meta = meta,
      tx2gene = tx2gene,
      gene_map = gene_map,
      type = "matrix"
    ))
  }
}

# utils_isoform_deps.R -- locating CPAT/Pfam data files and diagnosing the WSL/conda
# execution environment (debug_wsl() etc.)

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

#' Cheap sanity check that a local file is plausible downloaded data and not
#' an HTML error/mirror-selection page saved in its place (see the wget
#' validation in install_isoform_databases() for how that happens). Reads
#' raw bytes rather than text so it's safe to call on the binary
#' *_logitModel.RData files too.
#' @keywords internal
.data_file_looks_valid <- function(path, min_size = 100L) {
  if (!file.exists(path)) return(FALSE)

  sz <- file.info(path)$size
  if (is.na(sz) || sz < min_size) return(FALSE)

  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  head_raw <- tryCatch(readBin(con, "raw", n = 512L), error = function(e) raw(0))
  if (length(head_raw) == 0) return(FALSE)

  hay <- as.integer(head_raw)
  n <- length(hay)

  contains_bytes <- function(needle) {
    nb <- as.integer(charToRaw(needle))
    m <- length(nb)
    if (n < m) return(FALSE)
    for (i in seq_len(n - m + 1)) {
      if (identical(hay[i:(i + m - 1)], nb)) return(TRUE)
    }
    FALSE
  }

  !(contains_bytes("<!DOCTYPE") || contains_bytes("<!doctype") ||
    contains_bytes("<html") || contains_bytes("<HTML"))
}

#' Locate a CPAT data file (hexamer table or logit model) for a given organism
#'
#' CPAT's prebuilt hexamer/logit files are not bundled with the
#' IsoformSwitchAnalyzeR R package (that package only ships a small
#' `cpat_results.txt` example, not the actual model files); they are
#' downloaded by `install_isoform_databases()` into `$HOME/.cpat_data` (or
#' `$CPAT_DATA`). ExpressOM's own `inst/extdata` is checked first, so a file
#' placed there manually (e.g. as a workaround when the automated download
#' can't reach its source) takes precedence over a downloaded copy. Both
#' `.find_cpat_hexamer()` and `.find_cpat_logit_model()` search the same set
#' of locations.
#' @keywords internal
.find_cpat_data_file <- function(fname,
                                 wsl_distro = "Ubuntu",
                                 use_wsl = TRUE,
                                 conda_sh = NULL,
                                 conda_env = "isoform_tools",
                                 log_dir = NULL) {
  expressom_local <- system.file("extdata", fname, package = "ExpressOM")

  if (nzchar(expressom_local) && .data_file_looks_valid(expressom_local)) {
    if (.Platform$OS.type == "windows" && use_wsl) {
      return(.to_wsl_path(expressom_local, wsl_distro))
    }
    return(expressom_local)
  }

  isa_local <- system.file("extdata", fname, package = "IsoformSwitchAnalyzeR")

  if (nzchar(isa_local) && .data_file_looks_valid(isa_local)) {
    if (.Platform$OS.type == "windows" && use_wsl) {
      return(.to_wsl_path(isa_local, wsl_distro))
    }
    return(isa_local)
  }

  cpat_data_env <- Sys.getenv("CPAT_DATA", "")

  if (nzchar(cpat_data_env)) {
    candidate <- file.path(cpat_data_env, fname)

    if (.data_file_looks_valid(candidate)) {
      if (.Platform$OS.type == "windows" && use_wsl) {
        return(.to_wsl_path(candidate, wsl_distro))
      }
      return(candidate)
    }
  }

  # A candidate can exist and still be garbage -- e.g. a previous
  # install_isoform_databases() run that saved an HTML error/mirror page in
  # place of the real file (see the wget validation added there). Checking
  # existence alone was letting that silently pass through to CPAT, which
  # then failed deep inside cpmodule instead of here, where it's actually
  # diagnosable. Skip anything too small or HTML-looking and keep searching.
  result <- .wsl_exec_script(
    bash_body = c(
      sprintf('for _p in "$HOME/.cpat_data/%s" "${CPAT_DATA:-.}/%s"; do', fname, fname),
      '  [ -f "$_p" ] || continue',
      '  sz=$(wc -c < "$_p" 2>/dev/null || echo 0)',
      '  [ "$sz" -ge 100 ] || continue',
      '  head -c 512 "$_p" | grep -qi "<!doctype html\\|<html[ >]" && continue',
      '  echo "$_p"; break',
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
  .find_cpat_data_file(
    paste0(organism, "_logitModel.RData"),
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    conda_sh = conda_sh,
    conda_env = conda_env,
    log_dir = log_dir
  )
}

#' Locate the CPAT hexamer frequency table for a given organism
#' @keywords internal
.find_cpat_hexamer <- function(organism = "Human",
                               wsl_distro = "Ubuntu",
                               use_wsl = TRUE,
                               conda_sh = NULL,
                               conda_env = "isoform_tools",
                               log_dir = NULL) {
  .find_cpat_data_file(
    paste0(organism, "_Hexamer.tsv"),
    wsl_distro = wsl_distro,
    use_wsl = use_wsl,
    conda_sh = conda_sh,
    conda_env = conda_env,
    log_dir = log_dir
  )
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
#'
#' Diagnoses whether the WSL/conda environment that hosts the external
#' predictor tools (CPAT, SignalP 6, hmmscan, InterProScan, biolib/
#' DeepTMHMM) is reachable, whether `isoform_tools` exists, and which
#' tools and databases are actually visible.
#'
#' Earlier versions of this function had three false-negative bugs that
#' hid a working environment:
#'   * `command -v conda` was run without first sourcing `conda.sh`.
#'     Conda installs registered with `conda init bash` only inject
#'     themselves into PATH from an *interactive login* shell, but the
#'     shell used by `.wsl_exec_script()` is non-login, non-interactive.
#'     Result: `conda` was reported missing even when it existed, which
#'     in turn left `conda_env_exists = FALSE` and made every tool check
#'     fall back to the bare PATH.
#'   * Every tool probe used `conda_sh = NULL` first and only retried on
#'     the env-derived conda.sh -- but that value was computed *after*
#'     the conda check, which never succeeded (bug above). The tools were
#'     therefore never probed inside `isoform_tools`.
#'   * The CPAT-model check required *both* the Human and Mouse logit
#'     models. A single-species run is complete with one model; the
#'     Human-only case was still being reported as a failure.
#'
#' All three are fixed below: `conda.sh` is located once, up front, and
#' passed to every probe; `conda` is checked after sourcing it; the CPAT
#' model status is reported per species and treated as OK when at least
#' one model is present. `biolib` (DeepTMHMM) is now included in the
#' tool inventory so its status appears alongside the other predictors.
#'
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
    .expressom_log_dir(out_dir, create = TRUE)
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

  # FIX #1: pre-source conda.sh so `command -v conda` sees an install that
  # was registered with `conda init bash`. This also gives us a reusable
  # `pre_conda_sh` value to pass to every tool / database probe below,
  # breaking the circular dependency that made every check fall back to
  # the bare PATH.
  pre_conda_sh <- .find_conda_sh(
    wsl_distro = distro,
    use_wsl = via_wsl,
    log_dir = effective_log_dir
  )

  conda_check <- .wsl_exec_script(
    bash_body = c(
      if (!is.null(pre_conda_sh)) {
        sprintf('. %s 2>/dev/null || true', .dq(pre_conda_sh))
      },
      "command -v conda || echo 'not found'"
    ),
    wsl_distro = distro,
    use_wsl = via_wsl,
    intern = TRUE,
    ignore_stderr = TRUE,
    log_dir = effective_log_dir
  )

  conda_found <- any(grepl("conda", conda_check)) &&
    !any(grepl("^not found$", trimws(conda_check)))

  if (conda_found) {
    results$conda_available <- TRUE
    if (verbose) message("  \u2713 conda found")

    env_check <- .wsl_exec_script(
      bash_body = c(
        if (!is.null(pre_conda_sh)) {
          sprintf('. %s 2>/dev/null || true', .dq(pre_conda_sh))
        },
        sprintf(
          "conda env list 2>/dev/null | grep -q '%s' && echo 'exists'",
          conda_env
        )
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

  # Reuse the pre-computed conda.sh for everything below; keep the old name
  # so downstream code paths that referred to it still work.
  conda_sh_guess <- pre_conda_sh

  # FIX #2: try the bare PATH first (catches apt-installed / /usr/local/bin
  # tools), then the conda env. Previously the conda-env attempt only
  # happened when conda_env_exists was TRUE -- which required the (failing)
  # conda check. Now conda_sh_guess is always available.
  tools_list <- c("cpat", "signalp6", "hmmscan", "interproscan.sh", "biolib")

  .check_tool <- function(tool) {
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

    found
  }

  for (tool in tools_list) {
    found <- .check_tool(tool)
    results$tools[[tool]] <- found

    if (verbose) {
      message("   ", if (found) "\u2713" else "\u2717", "  ", tool)
    }
  }

  if (isTRUE(results$tools[["cpat"]])) {
    results$tools[["run_cpat.py"]] <- NA

    if (verbose) {
      message("   -  run_cpat.py (skipped, not needed: 'cpat' already found)")
    }
  } else {
    found <- .check_tool("run_cpat.py")
    results$tools[["run_cpat.py"]] <- found

    if (verbose) {
      message("   ", if (found) "\u2713" else "\u2717", "  run_cpat.py")
    }
  }

  tool_alternatives <- list(
    "CPAT (cpat / run_cpat.py)" = c("cpat", "run_cpat.py"),
    "SignalP (signalp6)" = c("signalp6"),
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

  # FIX #3: CPAT logit models are per-species; a human-only run does not
  # need the mouse model (and vice versa). Report each species and treat
  # the check as passing when at least one model is present.
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

  results$cpat_human_model <- !is.null(cpat_logit_human)
  results$cpat_mouse_model <- !is.null(cpat_logit_mouse)
  results$cpat_models_found <-
    results$cpat_human_model || results$cpat_mouse_model

  if (verbose) {
    message(
      "   ",
      if (results$cpat_human_model) "\u2713" else "\u2717",
      " CPAT Human logit model",
      if (is.null(cpat_logit_human)) " (needed only for human data)" else ""
    )
    message(
      "   ",
      if (results$cpat_mouse_model) "\u2713" else "\u2717",
      " CPAT Mouse logit model",
      if (is.null(cpat_logit_mouse)) " (needed only for mouse data)" else ""
    )
  }

  if (!results$cpat_models_found) {
    results$errors <- c(
      results$errors,
      "No CPAT logit models found (need at least one species)"
    )
  }

  .finish(results)
}