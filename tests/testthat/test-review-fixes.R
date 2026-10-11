test_that("LRT contrast statistics come from a directional Wald fit", {
  skip_if_not_installed("DESeq2")

  set.seed(918)
  n <- 12L
  condition <- factor(
    rep(c("base", "level"), each = n / 2),
    levels = c("base", "level")
  )
  counts <- matrix(rnbinom(80L * n, mu = 100, size = 10), nrow = 80L)
  counts[1:20, condition == "level"] <- counts[1:20, condition == "level"] * 4L
  counts[21:40, condition == "level"] <- pmax(
    1L,
    counts[21:40, condition == "level"] %/% 4L
  )
  colnames(counts) <- paste0("s", seq_len(n))
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(condition, row.names = colnames(counts)),
    design = ~condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
  lrt <- DESeq2::nbinomLRT(dds, reduced = ~1, quiet = TRUE)

  result <- .lrt_wald_contrast_results(
    lrt,
    c("condition", "level", "base"),
    alpha = 0.05
  )$results
  omnibus <- DESeq2::results(lrt)

  expect_true(any(is.finite(result$stat)))
  expect_true(any(result$stat < 0, na.rm = TRUE))
  expect_false(isTRUE(all.equal(result$stat, omnibus$stat)))
})

test_that("isoform checkpoint signatures gate checkpoint reuse", {
  path <- tempfile(fileext = ".rds")
  .checkpoint_save(list(value = 1), path, signature = "first")
  expect_identical(.checkpoint_load(path, signature = "first"), list(value = 1))
  expect_null(.checkpoint_load(path, signature = "different"))

  saveRDS(list(value = "legacy"), path)
  expect_null(.checkpoint_load(path, signature = "first"))
})

test_that("resume loader unwraps signed checkpoints and keeps legacy files", {
  save_dir <- withr::local_tempdir()
  signed <- list(value = "signed")
  legacy <- list(value = "legacy")
  .checkpoint_save(signed, file.path(save_dir, "isoform_import.rds"), "sig")
  saveRDS(legacy, file.path(save_dir, "dte_results.rds"))

  loaded <- load_isoform_results(save_dir)

  expect_identical(loaded$isoform_import, signed)
  expect_identical(loaded$dte_results, legacy)
  expect_null(loaded$dtu_results)
})

test_that("resume loader can require a matching signature", {
  save_dir <- withr::local_tempdir()
  .checkpoint_save(
    list(value = "matching"),
    file.path(save_dir, "dte_results.rds"),
    "current"
  )
  saveRDS(
    list(value = "legacy"),
    file.path(save_dir, "dtu_results.rds")
  )

  matching <- load_isoform_results(
    save_dir,
    signature = "current",
    require_signature = TRUE
  )
  expect_identical(matching$dte_results, list(value = "matching"))
  expect_null(matching$dtu_results)

  stale <- expect_warning(
    load_isoform_results(
      save_dir,
      signature = "stale",
      require_signature = TRUE
    ),
    "No saved RDS checkpoint files found"
  )
  expect_null(stale$dte_results)
})

test_that("resume loader preserves legacy behavior unless signatures are required", {
  save_dir <- withr::local_tempdir()
  saveRDS(list(value = "legacy"), file.path(save_dir, "dte_results.rds"))

  loaded <- load_isoform_results(save_dir)
  expect_identical(loaded$dte_results, list(value = "legacy"))
})

test_that("gene symbols resolve all mapped transcripts without recycling", {
  gene_map <- data.frame(
    ensembl = c("gene1", "gene2", "gene3"),
    symbol = c("DUP", "DUP", "OTHER")
  )
  tx2gene <- data.frame(
    tx_id = paste0("tx", 1:4),
    gene_id = c("gene1", "gene2", "gene3", "gene3")
  )

  expect_identical(
    .gene_transcript_ids(gene_map, tx2gene, "DUP"),
    c("tx1", "tx2")
  )
})

test_that("custom switch count IDs map to matching FASTA IDs", {
  custom_map <- data.frame(
    count_id = c("ENST0001.1", "ENST0002.1"),
    fasta_id = c("fa_tx1", "fa_tx2")
  )

  expect_identical(
    .remap_switch_tx2gene(c("ENST0001.1", "ENST0003.1"), custom_map),
    c("fa_tx1", "ENST0003")
  )
  expect_identical(
    .remap_switch_tx2gene(c("fa_tx1", "fa_tx2"), custom_map),
    c("fa_tx1", "fa_tx2")
  )
})

test_that("feature overlap validation rejects mismatched matrix namespaces", {
  tx2gene <- data.frame(tx_id = c("tx1", "tx2"), gene_id = c("g1", "g1"))
  expect_error(
    .validate_feature_overlap(c("gene1", "gene2", "tx1"), tx2gene, "matrix"),
    "matrix feature IDs match transcript IDs"
  )
  expect_invisible(.validate_feature_overlap(c("tx1", "tx2"), tx2gene, "matrix"))
})

test_that("empty matrix imports bypass feature-overlap proportion checks", {
  tx2gene <- data.frame(tx_id = "tx1", gene_id = "g1")
  expect_invisible(.validate_feature_overlap(character(), tx2gene, "matrix"))
})

test_that("overlap validation ignores missing IDs in its mapped count", {
  tx2gene <- data.frame(tx_id = "tx1", gene_id = "g1")
  expect_error(
    .validate_feature_overlap(c(NA_character_, "unmatched"), tx2gene, "matrix"),
    "Only 0 / 2"
  )
})

test_that("rank vectors exclude non-finite metrics and invalid identifiers", {
  expect_identical(
    .finite_ranked_vector(c(2, NA, Inf, -1), c("a", "b", "c", "d")),
    c(a = 2, d = -1)
  )
})

test_that("main contrast resolution uses filtered metadata", {
  meta <- data.frame(
    Sample = paste0("s", 1:4),
    condition = c("A", "A", "B", "B"),
    batch = c("A", "B", "A", "A")
  )
  filtered <- .apply_sample_filters(
    meta,
    "Sample",
    subset_sample = "Sample != 's2'"
  )

  expect_identical(
    .resolve_main_condition(
      stats::as.formula("~ condition + batch"),
      filtered,
      level = "B",
      base = "A"
    ),
    "condition"
  )
})

test_that("MSigDB arguments support old and new collection formals", {
  old_api <- .msigdbr_collection_args(
    c("species", "category", "subcategory"),
    "Homo sapiens",
    "C5",
    "GO:BP"
  )
  new_api <- .msigdbr_collection_args(
    c("species", "collection", "subcollection", "db_species"),
    "Homo sapiens",
    "C5",
    "GO:BP"
  )

  expect_identical(old_api$category, "C5")
  expect_identical(old_api$subcategory, "GO:BP")
  expect_identical(new_api$collection, "C5")
  expect_identical(new_api$subcollection, "GO:BP")
  expect_identical(new_api$db_species, "HS")
  expect_identical(
    .msigdbr_collection_args(
      c("species", "collection", "db_species"),
      "Mus musculus",
      "M2"
    )$db_species,
    "MM"
  )
  expect_identical(
    .finite_ranked_vector(numeric(), character()),
    stats::setNames(numeric(), character())
  )
})

test_that("PCA comparison keeps non-syntactic metadata names", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")

  counts <- matrix(
    c(2, 4, 6, 8, 3, 5, 7, 9, 4, 6, 8, 10),
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:4))
  )
  col_data <- S4Vectors::DataFrame(
    condition = factor(c("A", "A", "B", "B")),
    check.names = FALSE
  )
  names(col_data) <- "cell type"
  rownames(col_data) <- colnames(counts)
  vsd <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )
  out_dir <- withr::local_tempdir()

  expect_no_error(
    plot_comparison_pca(
      vsd,
      condition_col = "cell type",
      level = "B",
      base = "A",
      plot_dir = out_dir,
      ntop = NULL
    )
  )
})

test_that("isoform switch covariates are explicitly reported", {
  expect_warning(
    .warn_isoform_switch_covariates("~ batch + condition", "condition"),
    "batch.*not adjusted"
  )
  expect_no_warning(
    .warn_isoform_switch_covariates("~ condition", "condition")
  )
})

test_that("tx2gene mappings must be unique after transcript normalization", {
  expect_identical(
    .validate_tx2gene(data.frame(
      tx_id = c("tx1", "tx1"),
      gene_id = c("gene1", "gene1")
    )),
    data.frame(tx_id = "tx1", gene_id = "gene1")
  )
  expect_error(
    .validate_tx2gene(data.frame(
      tx_id = c("tx1", "tx1"),
      gene_id = c("gene1", "gene2")
    )),
    "multiple genes"
  )
  expect_error(
    .validate_tx2gene(data.frame(tx_id = NA_character_, gene_id = "gene1")),
    "missing or empty"
  )
  expect_error(
    .normalize_tx2gene(data.frame(
      tx_id = c("tx1.1", "tx1.2"),
      gene_id = c("gene1.1", "gene2.1")
    )),
    "multiple genes"
  )
})

test_that("reference resolution preserves a caller-supplied reference", {
  refs <- list(
    gtf = "downloaded.gtf",
    cdna_fasta = "downloaded.cdna.fa",
    ncrna_fasta = "downloaded.ncrna.fa"
  )
  resolved <- .resolve_isoform_references("custom.fa", NULL, refs)
  expect_identical(resolved$fasta, "custom.fa")
  expect_identical(resolved$gff, "downloaded.gtf")
  resolved <- .resolve_isoform_references(NULL, "custom.gtf", refs)
  expect_identical(
    resolved$fasta,
    c("downloaded.cdna.fa", "downloaded.ncrna.fa")
  )
  expect_identical(resolved$gff, "custom.gtf")
})

test_that("switch result fingerprints change when supplied results change", {
  dtu <- data.frame(feature_id = "tx1", adj_pvalue = 0.1)
  first <- .switch_results_fingerprint(NULL, dtu)
  second <- .switch_results_fingerprint(
    NULL,
    transform(dtu, adj_pvalue = 0.01)
  )
  expect_false(identical(first, second))
})

test_that("invalid custom annotation paths are rejected", {
  expect_error(
    .validate_optional_file("missing-tx2gene.tsv", "custom_tx2gene"),
    "custom_tx2gene.*existing file"
  )
  expect_error(
    .validate_optional_file(tempdir(), "custom_gene_map"),
    "custom_gene_map.*existing file"
  )
})

test_that("GMT inputs distinguish collections from missing file paths", {
  expect_identical(.validate_gmt_input(NULL), "H")
  expect_identical(.validate_gmt_input("C5:GO:BP"), "C5:GO:BP")
  gmt <- tempfile(fileext = ".gmt")
  writeLines("pathway\tdescription\tGENE1", gmt)
  on.exit(unlink(gmt), add = TRUE)
  expect_identical(.validate_gmt_input(gmt), gmt)
  expect_error(.validate_gmt_input("missing-pathways.gmt"), "does not exist")
  expect_error(.validate_gmt_input(tempdir()), "does not exist")
})

test_that("ORA Entrez universe retains genes without fold changes", {
  results <- data.frame(
    gene = c("g1", "g2", "g3"),
    entrezid = c("1", "2", NA),
    log2FoldChange = c(1, NA, 2)
  )
  expect_identical(
    .mapped_entrez_universe(results)$entrezid,
    c("1", "2")
  )
})

test_that("native installer defaults expand HOME but WSL defaults do not", {
  native <- .native_or_wsl_default_path("$HOME/.cpat_data", FALSE)
  expect_identical(native, path.expand("~/.cpat_data"))
  expect_identical(
    .native_or_wsl_default_path("$HOME/.cpat_data", TRUE),
    "$HOME/.cpat_data"
  )
})

test_that("invalid cached reference gzip files are replaced", {
  dest <- tempfile(fileext = ".gz")
  writeBin(charToRaw("truncated"), dest)
  on.exit(unlink(dest), add = TRUE)
  local_mocked_bindings(
    download.file = function(url, destfile, mode) {
      con <- gzfile(destfile, open = "wb")
      writeBin(as.raw(rep(1L, 200L)), con)
      close(con)
      0L
    },
    .package = "utils"
  )
  .ensure_valid_gzip_reference(
    "https://example.test/reference.gz",
    dest,
    "test"
  )
  expect_true(.is_valid_gzip_file(dest))
})

test_that("Reactome organism codes include rat", {
  expect_identical(.reactome_organism("hsa"), "human")
  expect_identical(.reactome_organism("mmu"), "mouse")
  expect_identical(.reactome_organism("rno"), "rat")
  expect_error(.reactome_organism("unknown"), "Unsupported Reactome")
})

test_that("sample subset expressions are constrained and validated", {
  samples <- data.frame(
    sample_id = c("s1", "s2", "s3"),
    group = c("A", "B", "A"),
    batch = c(1L, 1L, 2L)
  )

  expect_identical(
    .apply_sample_filters(samples, "sample_id", subset_sample = "group == 'A'"),
    samples[c(1, 3), , drop = FALSE]
  )
  expect_identical(
    .apply_sample_filters(
      samples,
      "sample_id",
      subset_sample = "batch %in% c(1, 2)"
    ),
    samples
  )
  expect_error(
    .apply_sample_filters(
      samples,
      "sample_id",
      subset_sample = "system('whoami')"
    ),
    "supports only"
  )
  expect_error(
    .apply_sample_filters(
      samples,
      "sample_id",
      subset_sample = "group == 'A' | NA"
    ),
    "non-missing logical"
  )
})
