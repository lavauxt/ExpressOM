test_that("RegionReport nBest is validated and capped", {
  expect_equal(.validate_nbest(20000), 1000)
  expect_equal(.validate_nbest(25), 25)
  expect_error(.validate_nbest(0), "positive integer")
  expect_error(.validate_nbest(1.5), "positive integer")
  expect_error(.validate_nbest(NA_real_), "positive integer")
  expect_error(.validate_nbest("many"), "positive integer")
})

test_that("report paths are escaped for Markdown destinations", {
  path <- .local_markdown_path("C:\\results with space\\plot (A)#1?.png")
  expect_identical(path, "C:/results%20with%20space/plot%20%28A%29%231%3F.png")
})

test_that("isoform designs retain covariates and isolate the condition term", {
  metadata <- data.frame(
    condition = c("control", "control", "treated", "treated"),
    batch = c("a", "b", "a", "b"),
    row.names = paste0("s", 1:4)
  )
  info <- .validate_isoform_design(~ batch + condition, "condition", metadata)
  metadata <- .set_contrast_reference(info$metadata, "condition", "control", "treated")
  expect_identical(
    .design_contrast_coef(info$formula, metadata, "condition", "treated", "control"),
    "conditiontreated"
  )

  dex <- .dexseq_usage_designs(info$formula, "condition")
  expect_true(all(c("exon", "condition", "batch") %in% all.vars(dex$full)))
  expect_false("condition" %in% all.vars(dex$reduced))
  expect_true("batch" %in% all.vars(dex$reduced))
})

test_that("isoform design validation rejects interactions and missing variables", {
  metadata <- data.frame(condition = c("control", "treated"))
  expect_error(
    .validate_isoform_design(~ batch + condition, "condition", metadata),
    "missing from sample metadata"
  )
  expect_error(
    .validate_isoform_design(~ batch * condition, "condition", metadata),
    "interactions are not supported"
  )
})

test_that("LRT p-values are not converted into contrast-direction labels", {
  expect_identical(
    as.character(.de_direction_label(c(2, -2), c(TRUE, TRUE), test_type = "LRT")),
    c("Not significant", "Not significant")
  )
  counts <- .deg_count_table(
    data.frame(padj = c(0.001, 0.001), log2FoldChange = c(2, -2)),
    0.05,
    test_type = "LRT"
  )
  expect_true(all(is.na(counts$n_up)))
  expect_true(all(is.na(counts$n_down)))
})
