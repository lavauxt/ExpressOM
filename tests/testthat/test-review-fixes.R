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
