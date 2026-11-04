# --- code walkers ---------------------------------------------------------------

# names accessed as `p$x` / `ctx$p$x` (parameters) and `ctx$x` (derived values)
accessed <- function(fn) {
  params <- character()
  derived <- character()
  is_name <- function(e, nm) is.name(e) && identical(as.character(e), nm)
  walk <- function(e) {
    if (!is.call(e)) {
      return(invisible())
    }
    if (is_name(e[[1L]], "$") && length(e) == 3L) {
      lhs <- e[[2L]]
      key <- as.character(e[[3L]])
      if (is_name(lhs, "p") || (is.call(lhs) && is_name(lhs[[1L]], "$") && is_name(lhs[[2L]], "ctx") && is_name(lhs[[3L]], "p"))) {
        params <<- c(params, key)
      } else if (is_name(lhs, "ctx") && key != "p") {
        derived <<- c(derived, key)
      }
    }
    for (i in seq_along(e)) {
      tryCatch(walk(e[[i]]), error = function(err) NULL)
    }
  }
  walk(body(fn))
  list(params = unique(params), derived = unique(derived))
}

stage_fns <- function() {
  c(import = ".stage_import", eda = ".stage_eda", dge = ".stage_dge",
    dge_plots = ".stage_dge_plots", functional = ".stage_functional",
    fgsea = ".stage_fgsea", isoform = ".stage_isoform", report = ".stage_report",
    rdata = ".stage_save_rdata")
}

# parameters read by helpers that receive the whole `p` / `ctx`
helper_reads <- list(
  import = character(),
  eda = "ensembl_package_name",
  dge = character(),
  dge_plots = "ensembl_package_name",
  functional = "ensembl_package_name",
  fgsea = "ensembl_package_name",
  isoform = c("ensembl_package_name", "count_type", "data_dir"),
  report = character(),
  rdata = "ensembl_package_name"
)

test_that("every parameter group names real arguments of expressom()", {
  groups <- .stage_param_groups()
  known <- names(formals(expressom))
  for (g in names(groups)) {
    expect_true(all(groups[[g]] %in% known), info = paste(g, ":", paste(setdiff(groups[[g]], known), collapse = ", ")))
  }
})

test_that("each stage reads only the parameters of its group (so invalidation is exact)", {
  groups <- .stage_param_groups()
  derived_ok <- .stage_derived_values()
  for (st in names(stage_fns())) {
    fn <- get(stage_fns()[[st]])
    acc <- accessed(fn)
    allowed <- c(groups[[st]], helper_reads[[st]])
    extra <- setdiff(acc$params, allowed)
    expect_identical(extra, character(0), info = paste0("stage '", st, "' reads parameters outside its group: ", paste(extra, collapse = ", ")))
    extra_d <- setdiff(acc$derived, derived_ok[[st]])
    expect_identical(extra_d, character(0), info = paste0("stage '", st, "' reads derived values outside its context: ", paste(extra_d, collapse = ", ")))
  }
})

test_that(".stage_ctx keeps only what a stage needs", {
  ctx <- list(
    p = expressom_default_params(),
    main_condition = "condition", model_eda = "~ condition", comp_name = "avsb"
  )
  imp <- .stage_ctx(ctx, "import")
  expect_setequal(names(imp$p), intersect(.stage_param_groups()$import, names(ctx$p)))
  expect_null(imp$main_condition)
  expect_null(imp$comp_name)
  expect_identical(.stage_ctx(ctx, "eda")$main_condition, "condition")
  expect_identical(.stage_ctx(ctx, "fgsea")$comp_name, "avsb")
  expect_null(.stage_ctx(ctx, "fgsea")$main_condition)
})

test_that("stage hashes change only for the stages that use the changed parameter", {
  p <- expressom_default_params()
  base <- .stage_hashes(p)
  changed <- function(...) {
    q <- p
    for (nm in names(list(...))) q[nm] <- list(list(...)[[nm]])
    names(base)[.stage_hashes(q) != base]
  }
  expect_identical(changed(volcano_colors = c(up = "orange", down = "navy", ns = "grey")), "dge_plots")
  expect_identical(changed(report_settings = list(title = "x")), "report")
  expect_setequal(changed(gsea_metric = "log2FoldChange"), c("functional", "fgsea", "isoform"))
  expect_setequal(
    changed(padj_cutoff = 0.2),
    c("dge", "dge_plots", "functional", "fgsea", "isoform")
  )
  # the count files are not touched by anything but the import / isoform inputs
  expect_false("import" %in% changed(padj_cutoff = 0.2, level = "x", model = "~ y"))
  expect_true("import" %in% changed(count_type = "kallisto"))
  expect_identical(changed(), character(0))
})

test_that(".out_manifest lists files relative to the output directory and tracks changes", {
  d <- withr::local_tempdir()
  expect_equal(nrow(.out_manifest(d, "Plots")), 0L)
  dir.create(file.path(d, "Plots"))
  writeLines("a", file.path(d, "Plots", "a.txt"))
  m1 <- .out_manifest(d, "Plots")
  expect_identical(m1$path, "Plots/a.txt")
  expect_named(m1, c("path", "size", "mtime"))
  writeLines("a longer content", file.path(d, "Plots", "a.txt"))
  writeLines("b", file.path(d, "Plots", "b.txt"))
  m2 <- .out_manifest(d, "Plots")
  expect_identical(m2$path, c("Plots/a.txt", "Plots/b.txt"))
  expect_false(identical(m1, m2))
  # sub-directories that do not exist are ignored
  expect_equal(nrow(.out_manifest(d, c("Plots", "nope"))), 2L)
})

test_that(".quant_input_files finds the quantification file of each tool", {
  d <- withr::local_tempdir()
  for (s in c("s1", "s2")) {
    dir.create(file.path(d, s))
    for (f in c("quant.sf", "abundance.tsv", "quant.genes.results", "quant.isoforms.results", "t_data.ctab")) {
      writeLines("x", file.path(d, s, f))
    }
  }
  n <- function(type, ...) length(.quant_input_files(list(count_type = type, data_dir = d), ...))
  expect_equal(n("salmon"), 2L)
  expect_equal(n("kallisto"), 2L)
  expect_equal(n("stringtie"), 2L)
  expect_equal(n("matrix"), 0L)
  expect_true(all(grepl("isoforms", .quant_input_files(list(count_type = "rsem", data_dir = d)))))
  expect_true(all(grepl("genes", .quant_input_files(list(count_type = "rsem", data_dir = d), gene_level = TRUE))))
  expect_length(.quant_input_files(list(count_type = "salmon", data_dir = file.path(d, "missing"))), 0L)
})

test_that("report settings from expressom(report_settings = ) override the defaults", {
  expect_identical(.as_report_settings(NULL), expressom_report_settings())
  s <- .as_report_settings(list(title = "Mine", top_terms = 7))
  expect_identical(s$title, "Mine")
  expect_equal(s$top_terms, 7)
  expect_identical(s$max_points, expressom_report_settings()$max_points)
})

test_that(".expressom_prepare rejects invalid cut-offs before doing any work", {
  p <- expressom_default_params()
  p$deg_lfc_cutoff <- -1
  expect_error(.expressom_prepare(p), "deg_lfc_cutoff")
  p$deg_lfc_cutoff <- c(1, 2)
  expect_error(.expressom_prepare(p), "deg_lfc_cutoff")
})

test_that(".stage_report builds the report of the comparison that was run", {
  skip_on_cran()
  skip_if_no_report_stack()
  skip_if_no_pandoc()
  out <- withr::local_tempdir()
  make_run_dir(out)
  ctx <- list(p = list(out_dir = out, run_dge = TRUE, level = "treated", base = "control", report_settings = NULL))
  path <- suppressMessages(.stage_report(ctx))
  expect_true(file.exists(path))
  expect_match(basename(path), "treated_vs_control")
  # a comparison that has no results falls back to the most recent one
  ctx$p$level <- "ghost"
  path2 <- suppressMessages(.stage_report(ctx))
  expect_true(file.exists(path2))
})

test_that("validate_environment accepts the report stack when it is installed", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("tximport")
  skip_if_not_installed("pheatmap")
  skip_if_no_report_stack()
  skip_if_no_pandoc()
  expect_no_warning(suppressMessages(
    validate_environment(run_isoform = FALSE, run_functional = FALSE, run_report = TRUE)
  ))
})
