test_that("DEXSeq feature p-values are adjusted across the combined results", {
  results <- data.frame(
    pvalue = c(0.01, 0.02, 0.04, NA_real_),
    padj = c(0.01, 0.02, 0.04, NA_real_)
  )

  adjusted <- .adjust_dexseq_feature_pvalues(results)

  expect_equal(adjusted$padj, stats::p.adjust(results$pvalue, method = "BH"))
  expect_equal(adjusted$padj[1:3], c(0.03, 0.03, 0.04))
})
