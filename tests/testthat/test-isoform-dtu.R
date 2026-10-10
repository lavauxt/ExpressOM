test_that("DEXSeq feature p-values are adjusted across the combined results", {
  results <- data.frame(
    pvalue = c(0.01, 0.02, 0.04, NA_real_),
    padj = c(0.01, 0.02, 0.04, NA_real_)
  )

  adjusted <- .adjust_dexseq_feature_pvalues(results)

  expect_equal(adjusted$padj, stats::p.adjust(results$pvalue, method = "BH"))
  expect_equal(adjusted$padj[1:3], c(0.03, 0.03, 0.04))
})

test_that("DRIMSeq count preparation converts abundance-derived counts to integers", {
  counts <- matrix(
    c(1.2, 2.7, 0.4, 4.5),
    nrow = 2,
    dimnames = list(c("tx1", "tx2"), c("s1", "s2"))
  )

  observed <- .integerize_dtu_counts(counts)

  expect_type(observed, "integer")
  expect_identical(observed, matrix(
    c(1L, 3L, 0L, 4L),
    nrow = 2,
    dimnames = dimnames(counts)
  ))
})

test_that("DTU p-values are normalized for reports from either DRIMSeq spelling", {
  from_pvalue <- .normalize_dtu_pvalues(data.frame(pvalue = c("0.01", "0.2")))
  from_p_value <- .normalize_dtu_pvalues(data.frame(p_value = c("0.01", "0.2")))

  expect_identical(from_pvalue$pvalue, c(0.01, 0.2))
  expect_identical(from_p_value$pvalue, c(0.01, 0.2))
  expect_error(.normalize_dtu_pvalues(data.frame(other = 1)), "No p-value column")
})
