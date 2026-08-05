# mod_isoform_import.R -- loading transcript-level counts for isoform analysis
#

#' Import transcript-level counts for isoform analysis
#' @export
import_transcript_counts <- function(data_dir,
                                     sample_table,
                                     ensembl_package_name,
                                     count_type = "salmon",
                                     matrix_file = NULL,
                                     subset_sample = NULL,
                                     remove_sample = NULL,
                                     custom_tx2gene = NULL,
                                     custom_gene_map = NULL) {

  if (!file.exists(sample_table)) stop("Sample table not found: ", sample_table)

  sample_df <- data.table::fread(sample_table, header = TRUE, data.table = FALSE)
  sample_col <- if ("Sample" %in% colnames(sample_df)) "Sample" else "sample_id"

  sample_df <- .apply_sample_filters(sample_df, sample_col, remove_sample, subset_sample)
  rownames(sample_df) <- sample_df[[sample_col]]

  edb <- getExportedValue(ensembl_package_name, ensembl_package_name)

  if (!is.null(custom_tx2gene) && file.exists(custom_tx2gene)) {
    message("Using custom tx2gene file: ", custom_tx2gene)

    tx2gene <- data.table::fread(custom_tx2gene, header = TRUE, data.table = FALSE)

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

  org_info <- get_organism_info(edb)
  org_db <- org_info$org_db
  org_obj <- if (requireNamespace(org_db, quietly = TRUE)) .load_org_db(org_db) else NULL

  if (!is.null(custom_gene_map) && file.exists(custom_gene_map)) {
    message("Using custom gene annotation file: ", custom_gene_map)

    gene_map <- data.table::fread(custom_gene_map, header = TRUE, data.table = FALSE)

    if (!("gene_id" %in% colnames(gene_map)) && "ensembl" %in% colnames(gene_map)) {
      colnames(gene_map)[colnames(gene_map) == "ensembl"] <- "gene_id"
    }

    if (!all(c("gene_id", "symbol") %in% colnames(gene_map))) {
      stop("Custom gene map must contain columns 'gene_id' (or 'ensembl') and 'symbol'")
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
    message("  Gene map loaded: ", nrow(gene_map), " genes, ", entrez_present, " with Entrez IDs")
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
      "rsem" = "quant.genes.results",
      "stringtie" = "t_data.ctab",
      stop("Unsupported count_type for tximport: ", count_type)
    )

    file_list <- sapply(sample_df[[sample_col]], function(sid) {
      p_nested <- file.path(data_dir, sid, paste0(sid, ".", count_type, ".quant"), count_file_name)
      p_direct <- file.path(data_dir, sid, count_file_name)

      if (file.exists(p_nested)) return(p_nested)
      return(p_direct)
    })

    names(file_list) <- sample_df[[sample_col]]

    txi <- tximport::tximport(
      file_list,
      type = count_type,
      txOut = TRUE,
      countsFromAbundance = "lengthScaledTPM"
    )

    rownames(txi$counts) <- clean_transcript_id(rownames(txi$counts))
    meta <- sample_df[colnames(txi$counts), , drop = FALSE]

    return(list(
      txi = txi,
      meta = meta,
      tx2gene = tx2gene,
      gene_map = gene_map,
      type = "tximport"
    ))
  } else {
    if (is.null(matrix_file)) stop("matrix_file required for count_type='matrix'")

    counts_df <- data.table::fread(matrix_file, data.table = FALSE)
    rownames(counts_df) <- counts_df[, 1]
    counts_df <- counts_df[, -1, drop = FALSE]

    valid_samples <- intersect(colnames(counts_df), rownames(sample_df))
    count_mat <- as.matrix(counts_df[, valid_samples, drop = FALSE])
    mode(count_mat) <- "numeric"
    count_mat[is.na(count_mat)] <- 0

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

