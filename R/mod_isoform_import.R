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
    tx2gene$tx_id <- strip_ensembl_version(tx2gene$tx_id)
    tx2gene$gene_id <- strip_ensembl_version(tx2gene$gene_id)
  } else {
    tx2gene <- ensembldb::transcripts(
      edb,
      columns = c("tx_id", "gene_id"),
      return.type = "DataFrame"
    )

    tx2gene <- as.data.frame(tx2gene)
    tx2gene$tx_id <- strip_ensembl_version(tx2gene$tx_id)
    tx2gene$gene_id <- strip_ensembl_version(tx2gene$gene_id)
  }

  tx2gene <- .validate_tx2gene(tx2gene)

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

    # Sanity check: rownames(txi$counts) should be transcript IDs that
    # overlap tx2gene$tx_id (already version-stripped above). If a
    # count_type's quantification file turns out to actually be gene-level
    # (as silently happens for "rsem" if a genes.results file is passed
    # in), or count_type/data_dir otherwise point at the wrong files, the
    # ID namespaces won't match and every downstream annotation join fails
    # silently (DTE) or with a hard-to-read error much later in the
    # pipeline (DTU/switch analysis). Catch it here, once, with a message
    # that points straight at the cause.
    n_total_tx <- length(rownames(txi$counts))
    n_mapped_tx <- sum(rownames(txi$counts) %in% tx2gene$tx_id)

    if (n_total_tx > 0 && n_mapped_tx / n_total_tx < 0.5) {
      stop(
        "Only ",
        n_mapped_tx,
        " / ",
        n_total_tx,
        " imported feature IDs match transcript IDs in tx2gene. This usually means ",
        "count_type = '",
        count_type,
        "' resolved to a gene-level quantification file ",
        "instead of an isoform-level one (RSEM in particular has separate ",
        "genes.results / isoforms.results outputs -- import_transcript_counts() needs ",
        "the isoform-level file), or data_dir/count_type point at the wrong files.\n",
        "Example imported feature IDs:\n  ",
        paste(utils::head(rownames(txi$counts), 5), collapse = "\n  "),
        "\n\nExample tx2gene tx_id values:\n  ",
        paste(utils::head(tx2gene$tx_id, 5), collapse = "\n  "),
        call. = FALSE
      )
    }

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
