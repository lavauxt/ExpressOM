`%||%` <- function(x, y) if (is.null(x)) y else x

.parse_ensembl_package_name <- function(ensembl_package_name) {
  if (length(ensembl_package_name) != 1L || is.na(ensembl_package_name)) {
    stop("`ensembl_package_name` must be one package name.", call. = FALSE)
  }
  matches <- regexec(
    "^EnsDb\\.(Hsapiens|Mmusculus)\\.v([0-9]+)$",
    ensembl_package_name
  )
  parts <- regmatches(ensembl_package_name, matches)[[1]]
  if (length(parts) != 3L) {
    stop(
      "Could not parse ensembl_package_name '",
      ensembl_package_name,
      "'. Expected format like 'EnsDb.Hsapiens.v107'.",
      call. = FALSE
    )
  }
  list(
    species = if (parts[[2]] == "Hsapiens") "human" else "mouse",
    release = parts[[3]]
  )
}

.resolve_ensembl_metadata <- function(species, release) {
  species <- tolower(as.character(species))
  release_text <- as.character(release)
  release_num <- suppressWarnings(as.numeric(release_text))
  if (
    length(species) != 1L || is.na(species) || !species %in% c("human", "mouse")
  ) {
    stop("`species` must be either 'human' or 'mouse'.", call. = FALSE)
  }
  if (
    length(release_text) != 1L ||
      !grepl("^[0-9]+$", release_text) ||
      !is.finite(release_num) ||
      release_num < 1 ||
      release_num != floor(release_num)
  ) {
    stop("`release` must be a positive integer Ensembl release.", call. = FALSE)
  }
  release <- sprintf("%.0f", release_num)
  if (species == "human") {
    list(
      species = species,
      release = as.character(release),
      package_prefix = "Hsapiens",
      org_folder = "homo_sapiens",
      org_scientific = "Homo_sapiens",
      genome_version = if (release_num <= 75) "GRCh37" else "GRCh38"
    )
  } else {
    list(
      species = species,
      release = as.character(release),
      package_prefix = "Mmusculus",
      org_folder = "mus_musculus",
      org_scientific = "Mus_musculus",
      genome_version = if (release_num <= 102) "GRCm38" else "GRCm39"
    )
  }
}

.sample_id_column <- function(sample_df) {
  columns <- colnames(sample_df)
  if ("Sample" %in% columns) {
    return("Sample")
  }
  if ("sample_id" %in% columns) {
    return("sample_id")
  }
  stop(
    "Sample table must contain a 'Sample' or 'sample_id' column.",
    call. = FALSE
  )
}

.resolve_quantification_files <- function(
  data_dir,
  sample_ids,
  count_type,
  count_file_name
) {
  files <- vapply(
    sample_ids,
    function(sid) {
      nested <- file.path(
        data_dir,
        sid,
        paste0(sid, ".", count_type, ".quant"),
        count_file_name
      )
      direct <- file.path(data_dir, sid, count_file_name)
      if (file.exists(nested)) nested else direct
    },
    character(1),
    USE.NAMES = FALSE
  )
  names(files) <- sample_ids

  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop(
      "Missing quantification files for samples: ",
      paste(names(missing), collapse = ", "),
      "\nExpected file like: ",
      count_file_name,
      call. = FALSE
    )
  }
  files
}

.match_matrix_samples <- function(matrix_samples, metadata_samples) {
  matched <- intersect(matrix_samples, metadata_samples)
  if (length(matched) == 0L) {
    stop(
      "No matching sample names between matrix and sample table.",
      call. = FALSE
    )
  }
  matched
}

.validate_nbest <- function(nbest, cap = 1000L) {
  if (
    !is.numeric(nbest) ||
      length(nbest) != 1L ||
      is.na(nbest) ||
      !is.finite(nbest) ||
      nbest < 1 ||
      nbest != floor(nbest)
  ) {
    stop("`nBest` must be a single positive integer.", call. = FALSE)
  }
  min(nbest, cap)
}

.local_markdown_path <- function(path) {
  path <- gsub("\\\\", "/", path)
  path <- gsub("%", "%25", path, fixed = TRUE)
  path <- gsub(" ", "%20", path, fixed = TRUE)
  path <- gsub("(", "%28", path, fixed = TRUE)
  path <- gsub(")", "%29", path, fixed = TRUE)
  path <- gsub("#", "%23", path, fixed = TRUE)
  path <- gsub("?", "%3F", path, fixed = TRUE)
  path
}

.validate_isoform_design <- function(design, condition, metadata) {
  design <- if (inherits(design, "formula")) {
    design
  } else {
    stats::as.formula(design)
  }
  if (length(design) != 2L) {
    stop(
      "Isoform design must be a one-sided additive formula, e.g. ~ batch + condition.",
      call. = FALSE
    )
  }
  labels <- attr(stats::terms(design), "term.labels")
  if (any(grepl(":|\\*|\\^|/", labels))) {
    stop(
      "Isoform design currently supports additive main effects only; interactions are not supported.",
      call. = FALSE
    )
  }
  variables <- all.vars(design)
  missing <- setdiff(variables, colnames(metadata))
  if (length(missing)) {
    stop(
      "Isoform design variable(s) missing from sample metadata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!condition %in% variables) {
    stop(
      "Isoform design must include the comparison condition column '",
      condition,
      "'.",
      call. = FALSE
    )
  }
  for (variable in variables) {
    if (is.character(metadata[[variable]])) {
      metadata[[variable]] <- factor(metadata[[variable]])
    }
  }
  list(formula = design, metadata = metadata, variables = variables)
}

.set_contrast_reference <- function(metadata, condition, base, level) {
  values <- as.character(metadata[[condition]])
  if (!all(c(base, level) %in% values)) {
    stop(
      "Both `base` and `level` must be observed in metadata column '",
      condition,
      "'.",
      call. = FALSE
    )
  }
  observed <- unique(values)
  metadata[[condition]] <- factor(
    values,
    levels = c(base, setdiff(observed, base))
  )
  metadata
}

.design_contrast_coef <- function(design, metadata, condition, level, base) {
  terms <- stats::terms(design)
  if (any(grepl(":|\\*|\\^|/", attr(terms, "term.labels")))) {
    stop(
      "DRIMSeq covariate designs currently require additive main effects.",
      call. = FALSE
    )
  }
  metadata <- .set_contrast_reference(metadata, condition, base, level)
  mm <- stats::model.matrix(design, data = metadata)
  expected <- paste0(condition, level)
  if (!expected %in% colnames(mm)) {
    stop(
      "Could not resolve contrast coefficient '",
      expected,
      "' from design matrix columns: ",
      paste(colnames(mm), collapse = ", "),
      call. = FALSE
    )
  }
  expected
}

.dexseq_usage_designs <- function(design, condition) {
  labels <- attr(stats::terms(design), "term.labels")
  covariates <- setdiff(labels, condition)
  feature_terms <- paste0("exon:`", c(covariates, condition), "`")
  reduced_terms <- if (length(covariates)) {
    paste0("exon:`", covariates, "`")
  } else {
    character(0)
  }
  list(
    full = stats::as.formula(paste(
      "~ sample + exon +",
      paste(feature_terms, collapse = " + ")
    )),
    reduced = stats::as.formula(paste(
      "~ sample + exon",
      if (length(reduced_terms)) {
        paste("+", paste(reduced_terms, collapse = " + "))
      } else {
        ""
      }
    ))
  )
}

.de_direction_colors <- function() {
  c(up = "red2", down = "royalblue", ns = "grey70")
}

.de_direction_label <- function(
  log2_fold_change,
  significant,
  test_type = "Wald"
) {
  if (length(test_type) == 1L) {
    test_type <- rep(test_type, length(log2_fold_change))
  }
  significant[!is.na(test_type) & test_type == "LRT"] <- FALSE
  significant <- !is.na(significant) & significant
  direction <- rep("Not significant", length(log2_fold_change))
  direction[
    significant & !is.na(log2_fold_change) & log2_fold_change > 0
  ] <- "Upregulated"
  direction[
    significant & !is.na(log2_fold_change) & log2_fold_change < 0
  ] <- "Downregulated"
  factor(
    direction,
    levels = c("Not significant", "Downregulated", "Upregulated")
  )
}

.de_direction_scale <- function() {
  cols <- .de_direction_colors()
  stats::setNames(
    cols[c("ns", "down", "up")],
    c("Not significant", "Downregulated", "Upregulated")
  )
}

.select_internal_db_archive <- function(tar_files, pkg_name = NULL) {
  if (length(tar_files) == 0L) {
    stop(
      "No .tar.gz database found in inst/extdata. Run create_homemade_db() first.",
      call. = FALSE
    )
  }
  if (is.null(pkg_name)) {
    return(tar_files[[1]])
  }
  expected <- paste0(pkg_name, ".tar.gz")
  matches <- tar_files[tolower(basename(tar_files)) == tolower(expected)]
  if (length(matches) != 1L) {
    stop(
      "No unique exact database archive matching '",
      pkg_name,
      "' found. Available databases:\n",
      paste(basename(tar_files), collapse = "\n"),
      call. = FALSE
    )
  }
  matches[[1]]
}

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
.fill_entrez_with_bitr <- function(
  gene_map,
  org_obj,
  id_col = "ensembl",
  symbol_col = "symbol"
) {
  if (is.null(org_obj)) {
    return(gene_map)
  }

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    message(
      "  clusterProfiler not installed; skipping advanced Entrez mapping."
    )
    return(gene_map)
  }

  idx_na <- is.na(gene_map$entrezid) | gene_map$entrezid == ""
  if (!any(idx_na)) {
    return(gene_map)
  }

  message(
    "  Attempting to fill missing Entrez IDs using clusterProfiler::bitr..."
  )

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
          gene_map$entrezid[i] <- map_df$ENTREZID[
            map_df$ENSEMBL == gene_map[[id_col]][i]
          ][1]
        }
      }
      message("    Mapped ", nrow(map_df), " Ensembl IDs to Entrez.")
    }
  }

  idx_na2 <- is.na(gene_map$entrezid) | gene_map$entrezid == ""

  if (any(idx_na2)) {
    syms <- gene_map[[symbol_col]][idx_na2]
    syms <- syms[
      !is.na(syms) & syms != "" & syms != gene_map[[id_col]][idx_na2]
    ]
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
.apply_sample_filters <- function(
  sample_df,
  sample_col,
  remove_sample = NULL,
  subset_sample = NULL
) {
  if (!is.null(remove_sample)) {
    message(
      "   -> Excluding requested samples: ",
      paste(remove_sample, collapse = ", ")
    )
    keep_indices <- !(sample_df[[sample_col]] %in% remove_sample)
    sample_df <- sample_df[keep_indices, , drop = FALSE]

    if (nrow(sample_df) == 0) {
      stop(
        "The remove_sample constraint removed all available samples from your metadata!"
      )
    }
  }

  if (!is.null(subset_sample)) {
    message("   -> Applying subset condition: ", subset_sample)

    filter_expr <- tryCatch(
      rlang::parse_expr(subset_sample),
      error = function(e) stop("Invalid subset_sample expression: ", e$message)
    )
    if (!.is_safe_subset_expr(filter_expr, names(sample_df))) {
      stop(
        "subset_sample supports only column names, literal values, c(), ",
        "comparisons, %in%, and logical operators (&, |, !).",
        call. = FALSE
      )
    }
    eval_env <- list2env(as.list(sample_df), parent = baseenv())
    subset_indices <- tryCatch(
      eval(filter_expr, envir = eval_env),
      error = function(e) {
        stop("Failed to evaluate subset_sample condition. Error: ", e$message)
      }
    )
    if (
      !is.logical(subset_indices) ||
        length(subset_indices) != nrow(sample_df) ||
        anyNA(subset_indices)
    ) {
      stop(
        "subset_sample must evaluate to a non-missing logical vector matching the sample table.",
        call. = FALSE
      )
    }
    sample_df <- sample_df[subset_indices, , drop = FALSE]

    if (nrow(sample_df) == 0) {
      stop("The subset_sample condition matched zero samples.")
    }
  }

  sample_df
}

.is_safe_subset_expr <- function(expr, columns) {
  if (is.symbol(expr)) {
    return(as.character(expr) %in% c(columns, "TRUE", "FALSE", "NA"))
  }
  if (!is.call(expr)) {
    return(is.atomic(expr) && length(expr) > 0L)
  }

  op <- as.character(expr[[1L]])
  args <- as.list(expr)[-1L]
  allowed <- c("==", "!=", ">", ">=", "<", "<=", "%in%", "&", "|", "!")
  if (identical(op, "c")) {
    return(all(vapply(
      args,
      function(arg) is.atomic(arg) && length(arg) > 0L,
      logical(1)
    )))
  }
  op %in%
    allowed &&
    length(args) >= 1L &&
    all(vapply(args, .is_safe_subset_expr, logical(1), columns = columns))
}

.validate_tx2gene <- function(tx2gene) {
  required <- c("tx_id", "gene_id")
  if (!all(required %in% names(tx2gene))) {
    stop("tx2gene must contain tx_id and gene_id columns.", call. = FALSE)
  }
  tx2gene <- tx2gene[, required, drop = FALSE]
  tx2gene$tx_id <- as.character(tx2gene$tx_id)
  tx2gene$gene_id <- as.character(tx2gene$gene_id)
  if (
    anyNA(tx2gene$tx_id) ||
      anyNA(tx2gene$gene_id) ||
      any(!nzchar(tx2gene$tx_id)) ||
      any(!nzchar(tx2gene$gene_id))
  ) {
    stop(
      "tx2gene contains missing or empty transcript/gene IDs.",
      call. = FALSE
    )
  }
  tx2gene <- unique(tx2gene)
  duplicate_tx <- unique(tx2gene$tx_id[duplicated(tx2gene$tx_id)])
  if (length(duplicate_tx)) {
    stop(
      "Transcript IDs in tx2gene map to multiple genes after ID normalization; e.g. ",
      paste(utils::head(duplicate_tx, 5L), collapse = ", "),
      call. = FALSE
    )
  }
  tx2gene
}

.object_md5 <- function(object) {
  path <- tempfile("expressom_signature_")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 2)
  unname(tools::md5sum(path))
}

.file_fingerprint <- function(paths) {
  paths <- unique(as.character(paths[!is.na(paths) & nzchar(paths)]))
  paths <- paths[file.exists(paths) & !dir.exists(paths)]
  if (!length(paths)) {
    return(character())
  }
  normalized <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  stats <- file.info(paths)
  data.frame(
    path = normalized,
    size = stats$size,
    mtime = as.numeric(stats$mtime),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

.checkpoint_save <- function(object, path, signature) {
  saveRDS(list(signature = signature, object = object), path)
  invisible(object)
}

.checkpoint_load <- function(path, signature) {
  checkpoint <- tryCatch(readRDS(path), error = function(e) NULL)
  if (
    is.null(checkpoint) ||
      !is.list(checkpoint) ||
      !identical(checkpoint$signature, signature) ||
      !"object" %in% names(checkpoint)
  ) {
    message("Ignoring incompatible or legacy checkpoint: ", basename(path))
    return(NULL)
  }
  checkpoint$object
}

.warn_isoform_switch_covariates <- function(model, main_condition) {
  covariates <- setdiff(all.vars(stats::as.formula(model)), main_condition)
  if (length(covariates)) {
    warning(
      "IsoformSwitchAnalyzeR switch tests use condition-only models; ",
      "the supplied design covariates (",
      paste(covariates, collapse = ", "),
      ") are not adjusted for in switch analysis.",
      call. = FALSE
    )
  }
  invisible(covariates)
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
  if (!is.factor(x)) {
    x <- as.factor(x)
  }

  if (is.null(base) || length(base) == 0 || !nzchar(as.character(base))) {
    message(
      "Note: no base level supplied for '",
      label,
      "'; skipping releveling (using existing factor level order: ",
      paste(levels(x), collapse = ", "),
      ")."
    )
    return(x)
  }

  if (!(as.character(base) %in% levels(x))) {
    message(
      "Note: base level '",
      base,
      "' not found among the levels of '",
      label,
      "' (found: ",
      paste(levels(x), collapse = ", "),
      "); skipping releveling."
    )
    return(x)
  }

  relevel(x, ref = as.character(base))
}

#' Identify which term of a design formula carries the level/base contrast
#'
#' Returns the design variable whose column in `meta` contains BOTH `level`
#' and `base` as observed values. This replaces the previous
#' `tail(all.vars(as.formula(model)), 1)` idiom, which silently picks the
#' LAST term in the formula -- fine for `~ condition`, wrong for any model
#' that puts a covariate after the contrasted factor (e.g.
#' `~ condition + donor`, which would try to relevel `donor` to `base` and
#' build a DESeq2 contrast on the wrong factor). Falls back to the last
#' term when nothing matches (e.g. EDA-only calls where level/base are
#' NULL, or when `meta` isn't available at call time).
#'
#' @param design_formula Formula object (not a string)
#' @param meta data.frame / DataFrame of sample metadata (may be NULL)
#' @param level Foreground level of the contrast (may be NULL)
#' @param base Reference level of the contrast (may be NULL)
#' @return Character scalar naming the design variable to relevel/contrast on,
#'   or NULL if the formula has no terms at all.
#' @keywords internal
.resolve_main_condition <- function(
  design_formula,
  meta = NULL,
  level = NULL,
  base = NULL
) {
  design_vars <- all.vars(design_formula)
  if (length(design_vars) == 0) {
    return(NULL)
  }

  if (!is.null(meta) && !is.null(level) && !is.null(base)) {
    for (var in rev(design_vars)) {
      if (var %in% colnames(meta)) {
        vals <- as.character(meta[[var]])
        if (all(c(level, base) %in% vals)) return(var)
      }
    }
  }

  tail(design_vars, 1)
}

#' Safely create a directory
#' @keywords internal
safe_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  invisible(path)
}

#' Whether this R build has Cairo support compiled in
#'
#' cairo_pdf()/cairo_ps() are optional at R's build time (see
#' ?grDevices::cairo) -- R's own docs say packages should check
#' capabilities("cairo") before relying on them. This is common enough to
#' hit in practice on native-Windows R (Cairo not compiled in, or a cairo
#' DLL issue) that a hard dependency on cairo_pdf() would silently produce
#' zero output files on some machines instead of merely mis-rendering a
#' few special characters. Checked once and cached: capabilities() doesn't
#' change during a session, and this is consulted from many call sites.
#' @keywords internal
.cairo_pdf_available <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- isTRUE(capabilities("cairo"))
    }
    cached
  }
})

#' The best available vector PDF device function for this R session
#'
#' cairo_pdf() when available (embeds real glyphs, avoiding "No display
#' font for 'Symbol'/'ArialUnicode'" when the PDF is later rasterized to
#' PNG for report embedding); the base pdf() device otherwise, so plots
#' still get saved -- with the original font-substitution risk, but that
#' is strictly better than not being saved at all.
#' @keywords internal
.pdf_device <- function() {
  if (.cairo_pdf_available()) grDevices::cairo_pdf else grDevices::pdf
}

#' Safely save a ggplot or base R plot to PDF
#' @keywords internal
safe_pdf <- function(path, expr, width = 10, height = 8) {
  expr_sub <- substitute(expr)
  caller_env <- parent.frame()
  dev_fun <- .pdf_device()

  tryCatch(
    {
      dev_fun(path, width = width, height = height)
      eval(expr_sub, envir = caller_env)
      dev.off()
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      message(
        "Warning: Failed to generate plot at: ",
        path,
        "\n  Error: ",
        e$message
      )
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
    if (file.exists(p1)) {
      return(normalizePath(p1, winslash = "/"))
    }

    p2 <- file.path(dir, "rmd", filename)
    if (file.exists(p2)) {
      return(normalizePath(p2, winslash = "/"))
    }

    dir <- dirname(dir)
  }

  stop(
    "Could not locate bundled template '",
    filename,
    "'.\n",
    "Run from the package root, reinstall the package, or set:\n",
    "  Sys.setenv(EXPRESSOM_RMD_DIR = '/path/to/inst/rmd')"
  )
}

#' Render a {{PLACEHOLDER}}-style template to a temp file with values substituted
#' @keywords internal
.render_placeholder_template <- function(
  template_file,
  values,
  fileext = ".Rmd"
) {
  template_path <- tryCatch(
    .expressom_rmd_path(template_file),
    error = function(e) NULL
  )

  if (is.null(template_path)) {
    warning(
      "Could not locate template '",
      template_file,
      "'. ",
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
#' @param species Either "human" or "mouse".
#' @param release Positive integer Ensembl release.
#' @param maintainer Package maintainer string.
#' @param author Package author string.
#' @param output_dir Directory receiving the generated source archive.
#' @return The normalized path to the generated source archive.
#' @export
create_homemade_db <- function(
  species = "human",
  release = "107",
  maintainer = "User <user@example.com>",
  author = "ExpressOM Builder",
  output_dir = "inst/extdata"
) {
  metadata <- .resolve_ensembl_metadata(species, release)
  release <- metadata$release
  pkg_name <- paste0("EnsDb.", metadata$package_prefix, ".v", release)
  tar_name <- paste0(pkg_name, ".tar.gz")

  tmp_dir <- tempfile(
    pattern = paste0("build_", pkg_name, "_"),
    tmpdir = tempdir()
  )
  if (!dir.create(tmp_dir, recursive = TRUE)) {
    stop("Could not create temporary build directory: ", tmp_dir, call. = FALSE)
  }
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  url <- sprintf(
    "https://ftp.ensembl.org/pub/release-%s/gtf/%s/%s.%s.%s.gtf.gz",
    release,
    metadata$org_folder,
    metadata$org_scientific,
    metadata$genome_version,
    release
  )

  gtf_path <- file.path(tmp_dir, basename(url))

  message("--- Step 1: Downloading GTF ---")
  message("Downloading from: ", url)

  old_timeout <- getOption("timeout")
  options(timeout = max(900L, old_timeout %||% 60L))
  on.exit(options(timeout = old_timeout), add = TRUE)
  download.file(url, destfile = gtf_path, mode = "wb")

  message("--- Step 2: Generating SQLite Database ---")

  db_file <- ensembldb::ensDbFromGtf(
    gtf = gtf_path,
    organism = metadata$org_scientific,
    genomeVersion = metadata$genome_version,
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
    stop(
      "Expected package folder '",
      pkg_name,
      "' not found in temp directory."
    )
  }

  withr::with_dir(tmp_dir, {
    utils::tar(tar_name, files = pkg_name, compression = "gzip")
  })
  tar_path <- file.path(tmp_dir, tar_name)
  if (!file.exists(tar_path) || file.info(tar_path)$size <= 0) {
    stop("Failed to create database archive: ", tar_path, call. = FALSE)
  }
  archive_members <- utils::untar(tar_path, list = TRUE)
  if (!any(archive_members == paste0(pkg_name, "/DESCRIPTION"))) {
    stop(
      "Database archive is missing the package DESCRIPTION: ",
      tar_path,
      call. = FALSE
    )
  }

  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop(
      "Could not create archive output directory: ",
      output_dir,
      call. = FALSE
    )
  }
  archive_path <- file.path(output_dir, tar_name)

  copied <- file.copy(
    tar_path,
    archive_path,
    overwrite = TRUE
  )
  if (!isTRUE(copied) || !file.exists(archive_path)) {
    stop(
      "Failed to write Ensembl database archive to: ",
      archive_path,
      call. = FALSE
    )
  }

  message("SUCCESS: Database bundled at ", archive_path)
  normalizePath(archive_path, winslash = "/", mustWork = TRUE)
}

#' Install Bundled Ensembl Database
#' @param pkg_name Optional exact EnsDb package name to install.
#' @param archive_path Optional source archive path returned by
#'   `create_homemade_db()`.
#' @export
install_internal_db <- function(pkg_name = NULL, archive_path = NULL) {
  if (!is.null(archive_path)) {
    if (length(archive_path) != 1L || !file.exists(archive_path)) {
      stop(
        "`archive_path` must name an existing database archive.",
        call. = FALSE
      )
    }
    tar_files <- normalizePath(archive_path, winslash = "/", mustWork = TRUE)
    if (!grepl("\\.tar\\.gz$", tar_files, ignore.case = TRUE)) {
      stop(
        "`archive_path` must point to a .tar.gz source archive.",
        call. = FALSE
      )
    }
    if (
      !is.null(pkg_name) &&
        tolower(basename(tar_files)) != tolower(paste0(pkg_name, ".tar.gz"))
    ) {
      stop(
        "Archive '",
        basename(tar_files),
        "' does not match package '",
        pkg_name,
        "'.",
        call. = FALSE
      )
    }
    members <- tryCatch(
      utils::untar(tar_files, list = TRUE),
      error = function(e) character()
    )
    if (!any(grepl("/DESCRIPTION$", members))) {
      stop(
        "Archive does not contain an R package DESCRIPTION file: ",
        tar_files,
        call. = FALSE
      )
    }
  } else {
    ext_path <- system.file("extdata", package = "ExpressOM")
    if (ext_path == "") {
      ext_path <- "inst/extdata"
    }
    tar_files <- list.files(
      ext_path,
      pattern = "\\.tar\\.gz$",
      full.names = TRUE
    )

    local_path <- "inst/extdata"
    local_files <- if (dir.exists(local_path)) {
      list.files(local_path, pattern = "\\.tar\\.gz$", full.names = TRUE)
    } else {
      character()
    }
    installed_match <- !is.null(pkg_name) &&
      any(
        tolower(basename(tar_files)) == tolower(paste0(pkg_name, ".tar.gz"))
      )
    if (!installed_match && (length(tar_files) == 0L || !is.null(pkg_name))) {
      tar_files <- local_files
    }
  }
  db_path <- .select_internal_db_archive(tar_files, pkg_name)
  if (is.null(pkg_name) && is.null(archive_path) && length(tar_files) > 1L) {
    warning(
      "Multiple databases found. Defaulting to the first one: ",
      basename(db_path),
      "\nUse install_internal_db(pkg_name = '...') to specify."
    )
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
#'
#' @note `msig_cat`/`msig_db` below are per-organism defaults kept for
#'   backward compatibility with any external caller reading them directly;
#'   nothing in this package uses them anymore. The actual, per-collection
#'   MSigDB resolution (human vs. mouse-native "M*" collections, subcategory
#'   fallback, etc.) lives in `.resolve_msigdbr_collection()` in
#'   mod_functional.R, which is the source of truth `run_fgsea_analysis()`
#'   and `run_local_enrichment()` both call through.
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
    sqrt(exp(-qchisq(p, df = 4, lower.tail = FALSE) / 2))
  } else {
    pnorm(qnorm(p) * sqrt(2) / 2)
  }
}

#' Download Ensembl Reference FASTA and GTF
#' @export
download_ensembl_refs <- function(
  ensembl_package_name,
  out_dir = "./reference"
) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  parsed <- .parse_ensembl_package_name(ensembl_package_name)
  metadata <- .resolve_ensembl_metadata(parsed$species, parsed$release)
  release <- metadata$release

  base_url <- "https://ftp.ensembl.org/pub/release-%s"

  gtf_url <- sprintf(
    paste0(base_url, "/gtf/%s/%s.%s.%s.gtf.gz"),
    release,
    metadata$org_folder,
    metadata$org_scientific,
    metadata$genome_version,
    release
  )

  cdna_url <- sprintf(
    paste0(base_url, "/fasta/%s/cdna/%s.%s.cdna.all.fa.gz"),
    release,
    metadata$org_folder,
    metadata$org_scientific,
    metadata$genome_version
  )

  ncrna_url <- sprintf(
    paste0(base_url, "/fasta/%s/ncrna/%s.%s.ncrna.fa.gz"),
    release,
    metadata$org_folder,
    metadata$org_scientific,
    metadata$genome_version
  )

  gtf_dest <- file.path(out_dir, basename(gtf_url))
  cdna_dest <- file.path(out_dir, basename(cdna_url))
  ncrna_dest <- file.path(out_dir, basename(ncrna_url))

  old_timeout <- getOption("timeout")
  options(timeout = max(1800L, old_timeout %||% 60L))
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
  missing_core <- core_pkgs[
    !sapply(core_pkgs, requireNamespace, quietly = TRUE)
  ]

  if (length(missing_core) > 0) {
    stop(paste(
      "Critical core packages are missing. Please install:",
      paste(missing_core, collapse = ", ")
    ))
  }

  if (run_functional) {
    func_pkgs <- c("clusterProfiler", "SPIA", "fgsea", "ReactomePA", "DOSE")
    missing_func <- func_pkgs[
      !sapply(func_pkgs, requireNamespace, quietly = TRUE)
    ]

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

  if (length(width) == 0 || is.na(width)) {
    width <- 8
  }
  if (length(height) == 0 || is.na(height)) {
    height <- 6
  }

  width <- max(min_dim, min(width, max_dim))
  height <- max(min_dim, min(height, max_dim))

  list(width = width, height = height)
}

#' Safely save a ggplot to file, capping runaway width/height and defaulting
#' to the best available PDF device so text renders correctly
#'
#' ggsave()'s default "pdf" device only *references* the base-14 PDF fonts
#' by name instead of embedding glyphs; as soon as a plot uses a character
#' outside those fonts' basic set, later PDF -> PNG rasterization (e.g. for
#' report embedding via convert_pdf_to_png()) fails with "No display font
#' for 'Symbol'" / "'ArialUnicode'" from Poppler/Ghostscript, and the plot
#' silently never appears in the report. cairo_pdf() embeds real glyphs via
#' Cairo/fontconfig, which avoids that class of error entirely -- but it
#' requires Cairo support compiled into R (see .pdf_device()), so the
#' default here falls back to plain pdf() rather than hard-failing when
#' Cairo isn't available.
#' @keywords internal
.safe_ggsave <- function(
  filename,
  plot,
  width,
  height,
  device = .pdf_device(),
  ...
) {
  if (is.null(plot)) {
    message(
      "   -> Skipping plot save, plot object is NULL: ",
      basename(filename)
    )
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
        "   -> Failed to save plot: ",
        basename(filename),
        "\n      Error: ",
        conditionMessage(e)
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
    gene_map <- gene_map[,
      setdiff(names(gene_map), extra_ensembl_cols),
      drop = FALSE
    ]
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

  gene_map <- gene_map[
    !is.na(gene_map$ensembl) & gene_map$ensembl != "",
    ,
    drop = FALSE
  ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl), , drop = FALSE]

  gene_map
}
