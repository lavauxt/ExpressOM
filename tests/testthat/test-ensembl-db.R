test_that("Ensembl package names and metadata resolve consistently", {
  expect_identical(
    .parse_ensembl_package_name("EnsDb.Hsapiens.v75"),
    list(species = "human", release = "75")
  )
  expect_identical(
    .parse_ensembl_package_name("EnsDb.Mmusculus.v103"),
    list(species = "mouse", release = "103")
  )

  expect_identical(
    .resolve_ensembl_metadata("human", 75)$genome_version,
    "GRCh37"
  )
  expect_identical(
    .resolve_ensembl_metadata("human", 76)$genome_version,
    "GRCh38"
  )
  expect_identical(
    .resolve_ensembl_metadata("mouse", 102)$genome_version,
    "GRCm38"
  )
  expect_identical(
    .resolve_ensembl_metadata("mouse", 103)$genome_version,
    "GRCm39"
  )

  expect_error(.parse_ensembl_package_name("EnsDb.Rn.v107"), "Expected format")
  expect_error(.resolve_ensembl_metadata("rat", 107), "human.*mouse")
  expect_error(.resolve_ensembl_metadata("human", 0), "positive integer")
  expect_error(.resolve_ensembl_metadata("human", "v107"), "positive integer")
  expect_error(.resolve_ensembl_metadata("human", "107.5"), "positive integer")
})

test_that("database archive selection requires exact package names", {
  paths <- c(
    "EnsDb.Hsapiens.v107.tar.gz",
    "EnsDb.Hsapiens.v1070.tar.gz"
  )

  expect_identical(
    .select_internal_db_archive(paths, "EnsDb.Hsapiens.v107"),
    paths[[1]]
  )
  expect_error(
    .select_internal_db_archive(paths, "EnsDb.Hsapiens.v10"),
    "No unique exact database"
  )
  expect_error(.select_internal_db_archive(character()), "No .tar.gz database")
})

test_that("internal database installation can target an explicit archive", {
  archive <- tempfile(fileext = ".tar.gz")
  file.create(archive)
  on.exit(unlink(archive), add = TRUE)
  expect_error(
    install_internal_db("EnsDb.Hsapiens.v107", archive_path = archive),
    "does not match package"
  )
})

test_that("database timeout is restored after download errors", {
  old_timeout <- getOption("timeout")
  options(timeout = 37L)
  on.exit(options(timeout = old_timeout), add = TRUE)
  local_mocked_bindings(
    download.file = function(...) stop("synthetic download failure"),
    .package = "utils"
  )
  expect_error(
    create_homemade_db(output_dir = withr::local_tempdir()),
    "synthetic download failure"
  )
  expect_identical(getOption("timeout"), 37L)
})

test_that("database maintainer default remains template-safe", {
  expect_identical(
    formals(create_homemade_db)$maintainer,
    "User <user@example.com>"
  )
})
