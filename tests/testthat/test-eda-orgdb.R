test_that("EDA identifier mapping falls back when the optional OrgDb is absent", {
  local_mocked_bindings(
    get_organism_info = function(edb) list(org_db = "org.NotInstalled.eg.db"),
    .load_org_db = function(org_db_name) stop("database unavailable"),
    .package = "ExpressOM"
  )

  ids <- c("ENSG000001", "ENSG000002")
  expect_identical(
    .map_eda_symbols(ids, edb = NULL),
    stats::setNames(ids, ids)
  )
})
