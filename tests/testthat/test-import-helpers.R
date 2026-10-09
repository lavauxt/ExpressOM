test_that("sample ID column resolution accepts supported names", {
  expect_identical(.sample_id_column(data.frame(Sample = "s1")), "Sample")
  expect_identical(.sample_id_column(data.frame(sample_id = "s1")), "sample_id")
  expect_error(
    .sample_id_column(data.frame(id = "s1")),
    "'Sample' or 'sample_id'"
  )
})

test_that("quantification files resolve direct and nested layouts", {
  data_dir <- withr::local_tempdir()
  dir.create(file.path(data_dir, "s1"))
  dir.create(file.path(data_dir, "s2", "s2.salmon.quant"), recursive = TRUE)
  file.create(file.path(data_dir, "s1", "quant.sf"))
  file.create(file.path(data_dir, "s2", "s2.salmon.quant", "quant.sf"))

  files <- .resolve_quantification_files(
    data_dir,
    c("s1", "s2"),
    "salmon",
    "quant.sf"
  )

  expect_identical(names(files), c("s1", "s2"))
  expect_identical(
    unname(files),
    c(
      file.path(data_dir, "s1", "quant.sf"),
      file.path(data_dir, "s2", "s2.salmon.quant", "quant.sf")
    )
  )
})

test_that("missing quantification files identify affected samples", {
  expect_error(
    .resolve_quantification_files(
      withr::local_tempdir(),
      c("s1", "s2"),
      "salmon",
      "quant.sf"
    ),
    "s1, s2.*quant.sf"
  )
})

test_that("matrix sample matching rejects an empty intersection", {
  expect_identical(
    .match_matrix_samples(c("s2", "s1"), c("s1", "s2")),
    c("s2", "s1")
  )
  expect_error(
    .match_matrix_samples(c("matrix_sample"), c("metadata_sample")),
    "No matching sample names"
  )
})
