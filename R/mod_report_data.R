# mod_report_data.R -- the "result contract" shared by the unified HTML report
# and the Shiny app.

# ---- small internal helpers -------------------------------------------------

.rep_require <- function(pkg, why = "to build the interactive report") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required ", why, ". ",
         "Install it with install.packages('", pkg, "').",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rep_read_json <- function(path) {
  if (is.null(path) || length(path) != 1L || !file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE),
           error = function(e) NULL)
}

#' Read a delimited result table (TSV/CSV/TXT) tolerantly
#' @keywords internal
.rep_read_table <- function(path, max_rows = Inf) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  if (isTRUE(file.info(path)$size == 0)) return(NULL)
  out <- tryCatch(
    data.table::fread(
      path,
      nrows = if (is.finite(max_rows)) as.numeric(max_rows) else Inf,
      data.table = FALSE, check.names = FALSE, fill = TRUE,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0L || ncol(out) == 0L) return(NULL)
  out
}

#' Find a table by stem name, trying .tsv then .txt then .csv.
#' @keywords internal
.rep_find_table <- function(dir, stem, exts = c("tsv", "txt", "csv")) {
  if (is.null(dir) || is.null(stem) || !nzchar(stem)) return(NULL)
  for (ext in exts) {
    p <- file.path(dir, paste0(stem, ".", ext))
    if (file.exists(p)) return(p)
  }
  NULL
}

.rep_rel <- function(paths, root) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", root), "/?"), "", paths)
}

.rep_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit)) hit[[1L]] else NA_character_
}

.rep_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

.rep_named_chr <- function(x) {
  if (is.null(x)) return(NULL)
  x <- unlist(x)
  if (is.null(names(x))) NULL else x
}

.rep_json_safe <- function(x, depth = 0L) {
  if (is.null(x) || depth > 3L) return(NULL)
  if (is.function(x) || is.environment(x) || is.data.frame(x) || isS4(x)) return(NULL)
  if (is.list(x)) {
    out <- lapply(x, .rep_json_safe, depth = depth + 1L)
    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0L) return(NULL)
    return(out)
  }
  if (is.atomic(x) && length(x) <= 200L) {
    nm <- names(x)
    if (!is.null(nm) && length(x) > 0L && all(nzchar(nm))) return(as.list(x))
    return(unname(x))
  }
  NULL
}

#' Write the run-configuration snapshot read back by the report and the app
#' @keywords internal
.write_run_config <- function(out_dir, params) {
  log_dir <- .expressom_log_dir(out_dir, create = TRUE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  safe <- .rep_json_safe(params)
  if (is.null(safe)) safe <- list()
  safe$written_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  safe$expressom_version <- tryCatch(
    as.character(utils::packageVersion("ExpressOM")),
    error = function(e) NA_character_
  )
  path <- file.path(log_dir, "expressom_run.json")
  tryCatch(
    jsonlite::write_json(safe, path, pretty = TRUE, auto_unbox = TRUE,
                         null = "null", na = "null"),
    error = function(e)
      warning("Could not write run configuration: ", conditionMessage(e))
  )
  invisible(path)
}

# ---- comparison discovery ----------------------------------------------------

#' List the comparisons available in an ExpressOM output directory
#'
#' @param out_dir Output directory of [expressom()].
#' @return Character vector of comparison names (`<level>_vs_<base>`).
#' @export
list_expressom_comparisons <- function(out_dir) {
  de_dir <- file.path(out_dir, "DE_raw_results")
  files <- if (dir.exists(de_dir)) {
    list.files(de_dir, pattern = "^DEgenes_raw_.*\\.(txt|tsv)$",
               full.names = TRUE)
  } else {
    character(0)
  }
  log_files <- unlist(lapply(.expressom_log_dirs(out_dir), function(d) {
    list.files(file.path(d, "DGE"), pattern = "^DGE_params_.*\\.json$",
               full.names = TRUE)
  }), use.names = FALSE)

  strip_ext <- function(basename) sub("\\.(txt|tsv)$", "", basename)
  nm <- c(
    strip_ext(sub("^DEgenes_raw_", "", basename(files))),
    sub("\\.json$", "", sub("^DGE_params_", "", basename(log_files)))
  )
  mt <- c(file.mtime(files), file.mtime(log_files))
  if (length(nm) == 0L) return(character(0))

  # a comparison is "available" only if its result table exists
  has_table <- nm %in% strip_ext(sub("^DEgenes_raw_", "", basename(files)))
  nm <- nm[has_table]
  mt <- mt[has_table]
  unique(nm[order(mt, decreasing = TRUE)])
}

# ---- enrichment inventory ---------------------------------------------------

.rep_enrichment_inventory <- function(out_dir) {
  empty <- data.frame(
    path = character(0), source = character(0), method = character(0),
    collection = character(0), kind = character(0), stringsAsFactors = FALSE
  )
  rel <- list.files(out_dir, pattern = "\\.(csv|tsv)$", recursive = TRUE,
                    ignore.case = TRUE)
  if (length(rel) == 0L) return(empty)

  roots <- c("ORA", "GSEA", "SPIA", "Transcription_Factors", "DecoupleR")
  rows <- lapply(rel, function(r) {
    parts <- strsplit(r, "/", fixed = TRUE)[[1L]]
    ri <- which(parts %in% roots)
    if (length(ri) == 0L) return(NULL)
    ri <- ri[[1L]]
    file <- parts[[length(parts)]]
    if (grepl("_significant\\.csv$", file, ignore.case = TRUE)) return(NULL)

    # ---- FIX: exclude the per-TF regulon artefacts -------------------------
    # The files written by .run_decouple_tf_regulon_plots() under
    # <resource>/plots/ (regulon_<TF>.csv, leading_edge_<TF>.csv) and the
    # top-level per-TF summary (tf_regulon_summary.csv) are NOT enrichment
    # tables. They are loaded separately by .rep_load_decouple_regulons()
    # and rendered under "TF regulon plots". Leaving them in the inventory
    # produced ghost entries under "Functional analysis -> COLLECTRI" (one
    # per significant TF) that showed as "no recognised term / adjusted
    # p-value columns".
    if (grepl("[/\\\\]plots[/\\\\]", r)) return(NULL)
    if (identical(file, "tf_regulon_summary.csv")) return(NULL)

    src <- if (any(grepl("^DTE_Enrichment$", parts))) "DTE"
           else if (any(grepl("^DTU_Enrichment$", parts))) "DTU"
           else "DGE"

    root <- parts[[ri]]
    sub_parts <- if (ri < length(parts) - 1L) {
      parts[(ri + 1L):(length(parts) - 1L)]
    } else character(0)

    if (root == "ORA") {
      kind <- "ora"
      method <- if (length(sub_parts)) paste(sub_parts[[1L]], "ORA") else "ORA"
      collection <- .rep_strip_comp(file)
    } else if (root == "GSEA") {
      if (length(sub_parts) && sub_parts[[1L]] == "FGSEA") {
        kind <- "gsea"
        method <- "FGSEA"
        collection <- if (length(sub_parts) > 1L) sub_parts[[2L]]
                      else .rep_strip_comp(file)
      } else {
        kind <- "gsea"
        method <- if (length(sub_parts)) paste(sub_parts[[1L]], "GSEA")
                  else "GSEA"
        collection <- .rep_strip_comp(file)
      }
    } else if (root == "SPIA") {
      kind <- "spia"; method <- "SPIA"; collection <- .rep_strip_comp(file)
    } else if (root == "DecoupleR") {
      kind <- "decouple"; method <- "decoupleR"
      collection <- if (length(sub_parts)) sub_parts[[1L]]
                    else .rep_strip_comp(file)
    } else {
      kind <- "tf"; method <- "Transcription factors"
      collection <- .rep_strip_comp(file)
    }

    data.frame(
      path = r, source = src, method = method,
      collection = collection, kind = kind,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(empty)
  inv <- do.call(rbind, rows)
  rownames(inv) <- NULL
  inv[order(inv$source, inv$method, inv$collection), , drop = FALSE]
}

.rep_strip_comp <- function(file) {
  x <- sub("\\.(csv|tsv)$", "", file, ignore.case = TRUE)
  sub("_[^_]+_vs_[^_]+$", "", x)
}

#' Standardise an enrichment table for plotting
#' @keywords internal
.rep_standardise_enrichment <- function(df, kind = "ora") {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  term_col <- .rep_first_col(df, c("Description", "pathway", "Name", "Term", "ID",
                                    "source"))
  padj_col <- .rep_first_col(df, c("p.adjust", "padj", "pGFdr", "Adjusted.P.value",
                                    "qvalue", "qvalues", "FDR", "p_value"))
  if (is.na(term_col) || is.na(padj_col)) return(NULL)
  size_col <- .rep_first_col(df, c("Count", "size", "setSize", "NDE", "Overlap"))
  score_col <- .rep_first_col(df, c("NES", "enrichmentScore", "tA",
                                     "Combined.Score", "score"))

  size <- if (is.na(size_col)) rep(NA_real_, nrow(df))
          else if (size_col == "Overlap") .rep_num(sub("/.*$", "", as.character(df[[size_col]])))
          else .rep_num(df[[size_col]])

  out <- data.frame(
    term  = as.character(df[[term_col]]),
    padj  = .rep_num(df[[padj_col]]),
    size  = size,
    score = if (is.na(score_col)) rep(NA_real_, nrow(df)) else .rep_num(df[[score_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$padj) & !is.na(out$term), , drop = FALSE]
  if (nrow(out) == 0L) return(NULL)
  out[order(out$padj), , drop = FALSE]
}

# ---- the loader ---------------------------------------------------------------

#' Load the results of an ExpressOM run into one typed list
#' @export
load_expressom_results <- function(
  out_dir,
  comparison = NULL,
  max_expr_genes = 200L,
  max_table_rows = Inf
) {
  if (!is.character(out_dir) || length(out_dir) != 1L || !dir.exists(out_dir)) {
    stop("`out_dir` must be an existing ExpressOM output directory.", call. = FALSE)
  }
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  notes <- character(0)

  comparisons <- list_expressom_comparisons(out_dir)
  if (is.null(comparison)) {
    comparison <- if (length(comparisons)) comparisons[[1L]] else NA_character_
  } else if (!comparison %in% comparisons) {
    stop(
      "Comparison '", comparison, "' not found in '", out_dir, "'. Available: ",
      if (length(comparisons)) paste(comparisons, collapse = ", ") else "none",
      call. = FALSE
    )
  }
  has_cmp <- !is.na(comparison)

  # ---- configuration ----------------------------------------------------------
  run_cfg <- .rep_read_json(.expressom_log_file(out_dir, "expressom_run.json"))
  dge_log <- if (has_cmp) {
    .rep_read_json(.expressom_log_file(
      out_dir, "DGE", paste0("DGE_params_", comparison, ".json")
    ))
  } else NULL

  lvl <- run_cfg$level %||% dge_log$contrast$level
  bse <- run_cfg$base %||% dge_log$contrast$base
  if ((is.null(lvl) || is.null(bse)) && has_cmp &&
      grepl("_vs_", comparison, fixed = TRUE)) {
    pieces <- strsplit(comparison, "_vs_", fixed = TRUE)[[1L]]
    lvl <- lvl %||% pieces[[1L]]
    bse <- bse %||% pieces[[2L]]
  }

  vcols <- .de_direction_colors()
  if (!is.null(.rep_named_chr(run_cfg$volcano_colors))) {
    given <- .rep_named_chr(run_cfg$volcano_colors)
    vcols[intersect(names(given), names(vcols))] <- given[intersect(names(given), names(vcols))]
  }
  pcols <- c(level = "red2", base = "royalblue")
  if (!is.null(.rep_named_chr(run_cfg$pca_colors))) {
    given <- .rep_named_chr(run_cfg$pca_colors)
    pcols[intersect(names(given), names(pcols))] <- given[intersect(names(given), names(pcols))]
  }

  config <- list(
    level = lvl, base = bse,
    model = run_cfg$model %||% dge_log$design_formula,
    test = run_cfg$test %||% dge_log$test_type %||% "Wald",
    reduced = run_cfg$reduced %||% {
      r <- dge_log$reduced_formula
      if (identical(r, "none")) NULL else r
    },
    main_condition = dge_log$contrast$factor %||% run_cfg$group_col,
    shrink_method = run_cfg$shrink_method %||% dge_log$shrinkage_method,
    padj_cutoff = .rep_num(run_cfg$padj_cutoff %||% dge_log$padj_cutoff %||% 0.05),
    lfc_cutoff = .rep_num(run_cfg$deg_lfc_cutoff %||% dge_log$lfc_cutoff %||% 1),
    deg_padj_cutoffs = run_cfg$deg_padj_cutoffs,
    count_type = run_cfg$count_type,
    volcano_colors = vcols, pca_colors = pcols,
    run = run_cfg
  )
  if (length(config$padj_cutoff) != 1L || is.na(config$padj_cutoff)) {
    config$padj_cutoff <- 0.05
  }
  if (length(config$lfc_cutoff) != 1L || is.na(config$lfc_cutoff)) {
    config$lfc_cutoff <- 1
  }

  # ---- differential expression ------------------------------------------------
  dge <- NULL
  if (has_cmp) {
    de_dir <- file.path(out_dir, "DE_raw_results")
    results <- .rep_read_table(
      .rep_find_table(de_dir, paste0("DEgenes_raw_", comparison))
    )
    if (!is.null(results)) {
      for (cn in intersect(
        c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "stat",
          "contrast_log2FoldChange", "contrast_pvalue", "contrast_padj", "contrast_stat"),
        colnames(results)
      )) {
        results[[cn]] <- .rep_num(results[[cn]])
      }
      if (!"gene" %in% colnames(results)) {
        results$gene <- if ("ensembl" %in% colnames(results)) {
          as.character(results$ensembl)
        } else {
          as.character(seq_len(nrow(results)))
        }
      }
      if ("test_type" %in% colnames(results) && !is.null(results$test_type[[1L]])) {
        config$test <- as.character(results$test_type[[1L]])
      }

      deg_counts <- .rep_read_table(
        .rep_find_table(de_dir, paste0("DEG_counts_", comparison))
      )
      sig_lists <- list.files(
        de_dir,
        pattern = paste0("^DEgenes_pval_.*_",
                         gsub("([.])", "\\\\\\1", comparison),
                         "\\.(txt|tsv)$")
      )

      dge <- list(
        results = results,
        deg_counts = deg_counts,
        sig_lists = sig_lists,
        expression = .rep_load_expression(de_dir, comparison, results, max_expr_genes)
      )
    } else {
      notes <- c(notes, paste0("DE result table for '", comparison, "' could not be read."))
    }
  }

  # ---- PCA --------------------------------------------------------------------
  pca <- .rep_load_pca(out_dir)

  if (!is.null(dge$expression) && length(pca)) {
    grp <- .rep_sample_groups(pca[[1L]], config$main_condition)
    if (!is.null(grp)) {
      dge$expression$groups <- grp[match(dge$expression$samples, names(grp))]
    }
  }

  # ---- functional analysis ------------------------------------------------------
  inv <- .rep_enrichment_inventory(out_dir)
  functional <- list(inventory = inv, tables = list())
  if (nrow(inv)) {
    read_one <- function(p) .rep_read_table(file.path(out_dir, p), max_rows = max_table_rows)
    tbl_list <- if (.Platform$OS.type != "windows" &&
                    parallel::detectCores() > 1L &&
                    length(inv$path) > 1L) {
      if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
      parallel::mclapply(inv$path, read_one,
                         mc.cores = min(parallel::detectCores() - 1L, 8L))
    } else {
      lapply(inv$path, read_one)
    }
    for (i in seq_along(tbl_list)) {
      if (!is.null(tbl_list[[i]])) functional$tables[[inv$path[[i]]]] <- tbl_list[[i]]
    }
    inv <- inv[inv$path %in% names(functional$tables), , drop = FALSE]
    functional$inventory <- inv
  }

  # ---- decoupleR ----------------------------------------------------------------
  decouple <- .rep_load_decouple(out_dir, max_table_rows)
  decouple_regulons <- .rep_load_decouple_regulons(out_dir, max_table_rows)

  # ---- isoforms -----------------------------------------------------------------
  isoform <- .rep_load_isoform(out_dir, max_table_rows)

  # ---- figures, logs, files -----------------------------------------------------
  figures <- .rep_list_figures(out_dir)
  logs    <- .rep_load_logs(out_dir, dge_log)
  files   <- .rep_file_inventory(out_dir)

  structure(
    list(
      out_dir = out_dir, comparison = comparison, comparisons = comparisons,
      config = config, dge = dge, pca = pca,
      functional = functional,
      decouple = decouple,
      decouple_regulons = decouple_regulons,
      isoform = isoform,
      figures = figures, logs = logs, files = files, notes = notes
    ),
    class = "expressom_results"
  )
}

.rep_load_expression <- function(de_dir, comparison, results, max_genes,
                                max_qc_genes = 1500L, max_distribution_genes = 5000L) {
  path <- .rep_find_table(de_dir, paste0("DESeq_Normalized_Counts_", comparison))
  if (is.null(path) || max_genes < 2L) return(NULL)
  nc <- .rep_read_table(path)
  if (is.null(nc) || !"gene" %in% colnames(nc)) return(NULL)
  annot <- c("ensembl", "symbol", "entrezid", "gene", "ensembl_clean",
             "gene_biotype", "seqnames", "start", "end", "width", "strand",
             "description")
  sample_cols <- setdiff(colnames(nc), annot)
  sample_cols <- sample_cols[vapply(nc[sample_cols], is.numeric, logical(1))]
  if (length(sample_cols) < 2L) return(NULL)

  # Keep all usable rows for sample-level QC. The heatmap remains focused on
  # the most significant DE genes, while clustering uses the most variable
  # genes across the full normalised matrix (not only DE hits).
  nc <- nc[!is.na(nc$gene) & nzchar(as.character(nc$gene)), , drop = FALSE]
  nc <- nc[!duplicated(as.character(nc$gene)), , drop = FALSE]
  full_mat <- as.matrix(nc[, sample_cols, drop = FALSE])
  storage.mode(full_mat) <- "double"
  rownames(full_mat) <- make.unique(as.character(nc$gene))
  full_mat[!is.finite(full_mat)] <- NA_real_
  log_full <- log2(full_mat + 1)

  row_var <- apply(log_full, 1L, stats::var, na.rm = TRUE)
  row_var[!is.finite(row_var)] <- -Inf
  n_qc <- min(nrow(full_mat), max(2L, as.integer(max_qc_genes)))
  qc_idx <- order(row_var, decreasing = TRUE)[seq_len(n_qc)]
  qc_mat <- full_mat[qc_idx, , drop = FALSE]

  # A deterministic, expression-range-spanning subset keeps the HTML compact
  # while preserving broad per-sample distribution diagnostics.
  n_dist <- min(nrow(full_mat), max(2L, as.integer(max_distribution_genes)))
  if (nrow(full_mat) > n_dist) {
    ord_mean <- order(rowMeans(full_mat, na.rm = TRUE), decreasing = TRUE)
    dist_idx <- unique(round(seq(1, length(ord_mean), length.out = n_dist)))
    dist_idx <- ord_mean[dist_idx]
  } else {
    dist_idx <- seq_len(nrow(full_mat))
  }
  distribution_mat <- full_mat[dist_idx, , drop = FALSE]

  ord <- results[!is.na(results$padj), , drop = FALSE]
  ord <- ord[order(ord$padj, -abs(ord$log2FoldChange)), , drop = FALSE]
  top <- unique(utils::head(ord$gene, max_genes))
  top_idx <- match(top, rownames(full_mat), nomatch = 0L)
  top_idx <- top_idx[top_idx > 0L]
  mat <- full_mat[top_idx, , drop = FALSE]
  if (nrow(mat) < 2L) return(NULL)

  list(
    matrix = mat, samples = sample_cols, groups = NULL,
    qc_matrix = qc_mat,
    distribution_matrix = distribution_mat,
    qc_gene_count = nrow(qc_mat),
    distribution_gene_count = nrow(distribution_mat)
  )
}

.rep_load_pca <- function(out_dir) {
  score_files <- list.files(out_dir, pattern = "_sample_scores\\.tsv$",
                            recursive = TRUE, full.names = TRUE)
  out <- list()
  for (f in score_files) {
    sc <- .rep_read_table(f)
    if (is.null(sc)) next
    pc_cols <- grep("^PC[0-9]+$", colnames(sc), value = TRUE)
    if (length(pc_cols) < 2L) next
    stub <- sub("_sample_scores\\.tsv$", "", basename(f))
    var_tbl <- .rep_read_table(file.path(dirname(f),
                                          paste0(stub, "_variance_explained.tsv")))
    variance <- NULL
    if (!is.null(var_tbl) && all(c("PC", "PercentVariance") %in% colnames(var_tbl))) {
      variance <- stats::setNames(.rep_num(var_tbl$PercentVariance), var_tbl$PC)
    }
    out[[stub]] <- list(
      name = stub, scores = sc, pcs = pc_cols,
      variance = variance, path = .rep_rel(f, out_dir)
    )
  }
  if (length(out) > 1L) {
    pri <- ifelse(grepl("ComparisonOnly", names(out)), 1L,
           ifelse(grepl("BatchCorrected", names(out)), 3L, 2L))
    out <- out[order(pri, names(out))]
  }
  out
}

.rep_sample_groups <- function(pca_entry, main_condition = NULL) {
  sc <- pca_entry$scores
  col <- .rep_pca_group_col(sc, main_condition)
  if (is.na(col)) return(NULL)
  stats::setNames(as.character(sc[[col]]), as.character(sc$sample_label))
}

.rep_pca_group_col <- function(scores, main_condition = NULL) {
  if (!is.null(main_condition) && main_condition %in% colnames(scores)) {
    return(main_condition)
  }
  skip <- c("sample_label", "plot_label", "sizeFactor", "replaceable")
  cand <- setdiff(colnames(scores), c(skip, grep("^PC[0-9]+$", colnames(scores), value = TRUE)))
  cand <- cand[vapply(scores[cand], function(v) {
    !is.numeric(v) && length(unique(v)) > 1L && length(unique(v)) < nrow(scores)
  }, logical(1))]
  if (length(cand)) cand[[1L]] else NA_character_
}

.rep_load_decouple <- function(out_dir, max_rows = Inf) {
  dec_dir <- file.path(out_dir, "DecoupleR")
  if (!dir.exists(dec_dir)) return(NULL)
  files <- list.files(dec_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  # Per-TF regulon / leading-edge CSVs live under <resource>/plots/ and are
  # loaded by .rep_load_decouple_regulons(); keep them out of the activity
  # table index.
  files <- files[!grepl("[/\\\\]plots[/\\\\]", files)]
  files <- files[!grepl("^tf_regulon_summary\\.csv$", basename(files))]
  if (!length(files)) return(NULL)
  out <- list()
  for (f in files) {
    tb <- .rep_read_table(f, max_rows = max_rows)
    if (is.null(tb) || !nrow(tb)) next
    resource <- basename(dirname(f))
    rel <- .rep_rel(f, out_dir)
    out[[rel]] <- list(resource = resource, table = tb, path = rel)
  }
  if (!length(out)) NULL else out
}

#' Load the per-TF regulon plots data written by run_decouple_analysis()
#'
#' Reads `DecoupleR/<resource>/plots/regulon_<TF>.csv`,
#' `DecoupleR/<resource>/plots/leading_edge_<TF>.csv`, and the top-level
#' `DecoupleR/tf_regulon_summary.csv`. Returns a list with `summary` and
#' `by_resource` (one entry per resource with `regulons` and `leading_edges`
#' named lists).
#' @keywords internal
.rep_load_decouple_regulons <- function(out_dir, max_rows = Inf) {
  dec_dir <- file.path(out_dir, "DecoupleR")
  if (!dir.exists(dec_dir)) return(NULL)

  summary_file <- list.files(
    dec_dir, pattern = "^tf_regulon_summary\\.csv$",
    recursive = TRUE, full.names = TRUE
  )
  summary <- if (length(summary_file)) {
    .rep_read_table(summary_file[1], max_rows = max_rows)
  } else NULL

  resource_dirs <- list.dirs(dec_dir, recursive = FALSE)
  resource_dirs <- resource_dirs[file.exists(file.path(resource_dirs, "plots"))]

  by_resource <- list()
  for (rd in resource_dirs) {
    resource  <- basename(rd)
    plots_dir <- file.path(rd, "plots")
    reg_files <- list.files(plots_dir, pattern = "^regulon_.*\\.csv$", full.names = TRUE)
    le_files  <- list.files(plots_dir, pattern = "^leading_edge_.*\\.csv$", full.names = TRUE)

    regulons <- list()
    leading_edges <- list()
    for (f in reg_files) {
      tf <- sub("^regulon_", "", sub("\\.csv$", "", basename(f)))
      tb <- .rep_read_table(f, max_rows = max_rows)
      if (!is.null(tb) && nrow(tb)) regulons[[tf]] <- tb
    }
    for (f in le_files) {
      tf <- sub("^leading_edge_", "", sub("\\.csv$", "", basename(f)))
      tb <- .rep_read_table(f, max_rows = max_rows)
      if (!is.null(tb) && nrow(tb)) leading_edges[[tf]] <- tb
    }
    if (length(regulons) || length(leading_edges)) {
      by_resource[[resource]] <- list(
        resource      = resource,
        regulons      = regulons,
        leading_edges = leading_edges
      )
    }
  }

  if (is.null(summary) && !length(by_resource)) return(NULL)
  list(summary = summary, by_resource = by_resource)
}

.rep_load_isoform <- function(out_dir, max_rows) {
  iso_dir <- file.path(out_dir, "IsoformSwitch")
  rep_dir <- file.path(iso_dir, "DTU_DTE_report")
  if (!dir.exists(iso_dir)) return(NULL)

  dte <- .rep_read_table(file.path(rep_dir, "DTE_results_annotated.csv"))
  dtu <- .rep_read_table(file.path(rep_dir, "DTU_results_annotated.csv"))
  dex <- .rep_read_table(file.path(rep_dir, "DEXSeq_results_annotated.csv"))
  switch_features <- .rep_read_table(file.path(rep_dir, "switch_features.csv"))

  methods <- NULL
  mp <- file.path(rep_dir, "analysis_methods.txt")
  if (file.exists(mp)) {
    methods <- tryCatch(readLines(mp, warn = FALSE), error = function(e) NULL)
  }

  tx_pca <- NULL
  tp <- list.files(file.path(iso_dir), pattern = "^PCA_transcripts_.*\\.tsv$",
                   recursive = TRUE, full.names = TRUE)
  if (length(tp)) tx_pca <- .rep_read_table(tp[[1L]], max_rows = max_rows)

  plots <- list.files(file.path(rep_dir, "plots"),
                      pattern = "\\.(png|jpg|jpeg)$",
                      recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  skipped <- .rep_read_table(file.path(rep_dir, "skipped_plots.csv"))

  if (is.null(dte) && is.null(dtu) && is.null(dex) &&
      is.null(switch_features) && length(plots) == 0L) {
    return(NULL)
  }
  list(
    dte = dte, dtu = dtu, dexseq = dex,
    switch_features = switch_features,
    methods = methods, transcript_pca = tx_pca,
    plots = plots, skipped_plots = skipped
  )
}

.rep_list_figures <- function(out_dir) {
  f <- list.files(out_dir, pattern = "\\.(pdf|png)$",
                  recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f <- f[!grepl("/IsoformSwitch/DTU_DTE_report/plots/", f, fixed = TRUE)]
  f <- f[!grepl("/Save_rdata/|/\\.", f)]
  # exclude the report HTML itself if it was generated here
  f <- f[!grepl("^ExpressOM_report.*\\.html$", basename(f))]
  if (length(f) == 0L) {
    return(data.frame(
      path = character(0), name = character(0), group = character(0),
      type = character(0), size = numeric(0), stringsAsFactors = FALSE
    ))
  }
  info <- file.info(f)
  rel <- .rep_rel(f, out_dir)
  data.frame(
    path = f, name = basename(f),
    group = vapply(strsplit(rel, "/", fixed = TRUE), function(p) {
      if (length(p) > 1L) p[[1L]] else "."
    }, character(1)),
    type = tolower(tools::file_ext(f)),
    size = info$size, stringsAsFactors = FALSE
  )
}

.rep_load_logs <- function(out_dir, dge_log) {
  log_dirs <- .expressom_log_dirs(out_dir)
  if (!length(log_dirs)) log_dirs <- .expressom_log_dir(out_dir)
  read_txt <- function(name) {
    p <- .expressom_log_file(out_dir, name)
    if (!file.exists(p)) return(NULL)
    tryCatch(readLines(p, warn = FALSE), error = function(e) NULL)
  }
  collect <- function(subdir) {
    paths <- unlist(lapply(log_dirs, function(d) {
      path <- file.path(d, subdir)
      if (dir.exists(path)) list.files(path, full.names = TRUE,
                                       recursive = TRUE) else character(0)
    }), use.names = FALSE)
    unique(paths[file.info(paths)$isdir %in% FALSE])
  }
  list(
    dge_params = dge_log,
    session_info = read_txt("SessionInfo.txt"),
    warnings = read_txt("Warnings.txt"),
    isoform_logs = collect("Isoform"),
    wsl_logs = collect("WSL"),
    omnipath_logs = collect("omnipathr"),
    targets_logs = collect("targets")
  )
}

.rep_file_inventory <- function(out_dir) {
  f <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (length(f) == 0L) {
    return(data.frame(
      path = character(0), category = character(0), size_kb = numeric(0),
      modified = character(0), stringsAsFactors = FALSE
    ))
  }
  info <- file.info(f)
  rel <- .rep_rel(f, out_dir)
  data.frame(
    path = rel,
    category = vapply(strsplit(rel, "/", fixed = TRUE), function(p) {
      if (length(p) > 1L) p[[1L]] else "."
    }, character(1)),
    size_kb = round(info$size / 1024, 1),
    modified = format(info$mtime, "%Y-%m-%d %H:%M"),
    stringsAsFactors = FALSE
  )
}

#' @export
print.expressom_results <- function(x, ...) {
  cat("<expressom_results>\n")
  cat("  out_dir    :", x$out_dir, "\n")
  cat("  comparison :",
      if (is.na(x$comparison)) "none (EDA / isoform-only run)" else x$comparison,
      if (length(x$comparisons) > 1L)
        paste0("  [", length(x$comparisons), " available]") else "", "\n")
  cat("  DGE        :",
      if (is.null(x$dge)) "-" else paste(nrow(x$dge$results), "genes"), "\n")
  cat("  PCA        :",
      if (length(x$pca)) paste(length(x$pca), "plot(s)") else "-", "\n")
  cat("  functional :",
      if (nrow(x$functional$inventory)) paste(nrow(x$functional$inventory), "table(s)")
      else "-", "\n")
  cat("  decoupleR  :",
      if (!is.null(x$decouple)) paste(length(x$decouple), "table(s)") else "-", "\n")
  if (!is.null(x$decouple_regulons)) {
    n_reg <- sum(vapply(
      x$decouple_regulons$by_resource,
      function(r) length(r$regulons),
      integer(1)
    ))
    cat("  decoupleR regulons  :", n_reg, "TF(s) with regulon plots\n")
  }
  cat("  isoform    :",
      if (is.null(x$isoform)) "-" else paste(
        c(if (!is.null(x$isoform$dte)) "DTE",
          if (!is.null(x$isoform$dtu)) "DTU",
          if (!is.null(x$isoform$dexseq)) "DEXSeq",
          if (!is.null(x$isoform$switch_features)) "switch_features"),
        collapse = ", "), "\n")
  invisible(x)
}

#' Where are the megabytes of an ExpressOM report?
#' @export
audit_expressom_report <- function(path) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  idx <- grepl('<script type="application/json" data-for="', lines, fixed = TRUE)
  w <- data.frame(
    widget  = sub('.*data-for="([^"]+)".*', '\\1', lines[idx]),
    json_kb = round(nchar(lines[idx]) / 1024),
    stringsAsFactors = FALSE
  )
  w <- w[order(-w$json_kb), , drop = FALSE]
  cat(sprintf("file: %.1f MB | widgets: %d | widget JSON: %.1f MB\n",
              file.size(path) / 1e6, nrow(w), sum(w$json_kb) / 1024))
  if (nrow(w)) print(utils::head(w, 15), row.names = FALSE)
  invisible(w)
}