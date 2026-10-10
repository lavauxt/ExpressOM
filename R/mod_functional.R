#' @keywords internal
.write_enrich_csv <- function(obj, path) {
  df <- as.data.frame(obj)
  if (nrow(df) > 0) {
    write.csv(df, path, row.names = FALSE)
  }
  invisible(df)
}

#' @keywords internal
.run_go_ontology <- function(
  ont,
  sigOE_genes,
  allOE_genes,
  org_db,
  go_pvalue_cutoff,
  go_qvalue_cutoff
) {
  message("   -> Running GO ", ont, "...")
  ego <- safe_run(
    suppressMessages(clusterProfiler::enrichGO(
      gene = sigOE_genes,
      universe = allOE_genes,
      keyType = "SYMBOL",
      OrgDb = .load_org_db(org_db),
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = go_pvalue_cutoff,
      qvalueCutoff = go_qvalue_cutoff,
      readable = FALSE
    )),
    label = paste("GO", ont)
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("      No significant results for GO ", ont)
    return(NULL)
  }
  ego
}

#' Write CSV/plot outputs for a completed GO ORA result
#' @keywords internal
.write_go_ontology_outputs <- function(
  ego,
  ont,
  OE_foldchanges,
  top_genes,
  dir_go,
  comp_name
) {
  .write_enrich_csv(
    ego,
    file.path(dir_go, paste0("GO_", ont, "_", comp_name, ".csv"))
  )

  ego_wrapped <- ego
  if (nrow(ego@result) > 5) {
    ego_wrapped@result$Description <- stringr::str_wrap(
      ego_wrapped@result$Description,
      width = 50
    )
    safe_pdf(
      file.path(dir_go, paste0("GO_", ont, "_Dotplot_", comp_name, ".pdf")),
      width = 14,
      height = 14,
      expr = print(enrichplot::dotplot(
        ego_wrapped,
        showCategory = top_genes,
        label_format = 50
      ))
    )

    safe_pdf(
      file.path(dir_go, paste0("GO_", ont, "_Cnetplot_", comp_name, ".pdf")),
      width = 14,
      height = 14,
      expr = print(
        enrichplot::cnetplot(
          ego_wrapped,
          showCategory = top_genes,
          foldChange = OE_foldchanges
        ) +
          ggplot2::scale_color_gradient2(
            low = .de_direction_colors()[["down"]],
            mid = "white",
            high = .de_direction_colors()[["up"]],
            midpoint = 0,
            name = "log2 fold change"
          )
      )
    )
  }
  ego_wrapped
}

#' @keywords internal
.hdo_schema_ok <- function(path) {
  tryCatch(
    {
      con <- DBI::dbConnect(RSQLite::SQLite(), path)
      on.exit(DBI::dbDisconnect(con))
      "gene2allont" %in% DBI::dbListTables(con)
    },
    error = function(e) FALSE
  )
}

#' @keywords internal
.ensure_hdo_sqlite <- function() {
  cache_dir <- tryCatch(
    {
      if (requireNamespace("rappdirs", quietly = TRUE)) {
        rappdirs::user_cache_dir("GOSemSim")
      } else {
        if (Sys.info()["sysname"] == "Windows") {
          file.path(Sys.getenv("LOCALAPPDATA"), "GOSemSim")
        } else {
          file.path("~", ".cache", "GOSemSim")
        }
      }
    },
    error = function(e) file.path(tempdir(), "GOSemSim")
  )

  safe_dir(cache_dir)
  hdo <- file.path(cache_dir, "HDO.sqlite")

  if (file.exists(hdo) && !.hdo_schema_ok(hdo)) {
    message("HDO.sqlite outdated; removing.")
    file.remove(hdo)
  }
  if (file.exists(hdo)) {
    return(hdo)
  }

  message("Downloading HDO.sqlite.gz...")
  safe_dir(cache_dir)
  gz_tmp <- tempfile(fileext = ".sqlite.gz")
  old_to <- getOption("timeout")
  options(timeout = 300)
  on.exit(options(timeout = old_to), add = TRUE)

  urls <- c(
    "https://yulab-smu.top/DOSE/HDO.sqlite.gz",
    "https://raw.githubusercontent.com/YuLab-SMU/DOSE/refs/heads/gh-pages/HDO.sqlite.gz"
  )
  for (url in urls) {
    ok <- tryCatch(
      {
        download.file(url, gz_tmp, mode = "wb", quiet = TRUE)
        TRUE
      },
      error = function(e) FALSE
    )
    if (!ok) {
      next
    }
    tryCatch(
      {
        con_gz <- gzcon(file(gz_tmp, "rb"))
        tryCatch(
          writeBin(readBin(con_gz, "raw", n = 200e6), hdo),
          error = function(e) message("Decompression failed: ", e$message),
          finally = try(close(con_gz), silent = TRUE)
        )
      },
      error = function(e) message("Could not open gz file: ", e$message)
    )
    if (.hdo_schema_ok(hdo)) {
      return(hdo)
    }
    if (file.exists(hdo)) file.remove(hdo)
  }
  NULL
}

#' @keywords internal
.run_disease_ontology <- function(
  sig_entrez_ids,
  universe_entrez,
  go_pvalue_cutoff,
  go_qvalue_cutoff
) {
  hdo <- .ensure_hdo_sqlite()
  if (is.null(hdo)) {
    message("Skipping Disease Ontology (no HDO.sqlite).")
    return(invisible(NULL))
  }

  do_res <- safe_run(
    DOSE::enrichDO(
      gene = sig_entrez_ids,
      universe = universe_entrez,
      pvalueCutoff = go_pvalue_cutoff,
      qvalueCutoff = go_qvalue_cutoff
    ),
    label = "DOSE enrichDO"
  )
  if (is.null(do_res) || nrow(as.data.frame(do_res)) == 0) {
    return(invisible(NULL))
  }
  do_res
}

#' Write CSV/plot outputs for a completed Disease Ontology result
#' @keywords internal
.write_disease_ontology_outputs <- function(
  do_res,
  top_genes,
  dir_ora,
  comp_name
) {
  dir_dose <- safe_dir(file.path(dir_ora, "DOSE"))
  .write_enrich_csv(
    do_res,
    file.path(dir_dose, paste0("DOSE_Enrichment_", comp_name, ".csv"))
  )
  if (nrow(do_res@result) > 5) {
    safe_pdf(
      file.path(dir_dose, paste0("DOSE_Dotplot_", comp_name, ".pdf")),
      width = 14,
      height = 14,
      expr = print(enrichplot::dotplot(
        do_res,
        showCategory = top_genes,
        label_format = 50
      ))
    )
  }
  invisible(do_res)
}

#' Run KEGG Pathview for GSEA results
#' @keywords internal
.run_kegg_pathview <- function(
  gseaKEGG,
  expression_vector,
  kegg_code,
  top_genes,
  dir_kegg
) {
  gseaKEGG_results <- as.data.frame(gseaKEGG)
  if (nrow(gseaKEGG_results) == 0) {
    return(invisible(NULL))
  }

  kegg_ids <- as.character(gseaKEGG_results$ID)
  if (length(kegg_ids) > top_genes) {
    kegg_ids <- kegg_ids[seq_len(top_genes)]
  }
  kegg_ids <- ifelse(
    grepl("^[a-zA-Z]", kegg_ids),
    kegg_ids,
    paste0(kegg_code, kegg_ids)
  )

  overview_map_ids <- c(
    "01100",
    "01110",
    "01120",
    "01200",
    "01210",
    "01212",
    "01230",
    "01232",
    "01240",
    "01250"
  )
  overview_pattern <- paste0(
    "^",
    kegg_code,
    "(",
    paste(overview_map_ids, collapse = "|"),
    ")$"
  )
  is_overview <- grepl(overview_pattern, kegg_ids)
  if (any(is_overview)) {
    message(
      "   -> Skipping ",
      sum(is_overview),
      " KEGG global/overview map(s) ",
      "(known to crash pathview's per-node layout): ",
      paste(kegg_ids[is_overview], collapse = ", ")
    )
    kegg_ids <- kegg_ids[!is_overview]
  }
  if (length(kegg_ids) == 0) {
    return(invisible(NULL))
  }

  gene_data <- expression_vector[!is.na(expression_vector)]
  if (length(gene_data) == 0) {
    message("   -> No valid gene expression data for KEGG pathview.")
    return(invisible(NULL))
  }

  if (!"package:pathview" %in% search()) {
    try(attachNamespace("pathview"), silent = TRUE)
  }

  succeeded <- character(0)
  failed <- character(0)

  for (pid in kegg_ids) {
    devs_before <- grDevices::dev.list()
    result <- safe_run(
      {
        withr::with_dir(dir_kegg, {
          pv_out <- pathview::pathview(
            gene.data = gene_data,
            pathway.id = pid,
            species = kegg_code,
            gene.idtype = "entrez",
            limit = list(gene = 2, cpd = 1),
            low = list(gene = .de_direction_colors()[["down"]], cpd = "blue"),
            mid = list(gene = "white", cpd = "gray"),
            high = list(gene = .de_direction_colors()[["up"]], cpd = "yellow"),
            kegg.dir = "."
          )

          gd <- tryCatch(pv_out$plot.data.gene, error = function(e) NULL)
          if (!is.null(gd) && nrow(gd) > 0) {
            n_matched <- sum(!is.na(gd$mol.data))
            message(
              "   -> ",
              pid,
              ": ",
              n_matched,
              "/",
              nrow(gd),
              " pathway genes had expression data to color."
            )
            if (n_matched == 0) {
              message(
                "      (0 matched -- check that gene_data's Entrez IDs actually cover ",
                "this pathway's genes; a custom/SQANTI3-derived gene_map with many novel ",
                "loci can leave canonical pathway genes without a valid Entrez ID.)"
              )
            }
          }
          TRUE
        })
      },
      label = paste("KEGG pathview", pid)
    )

    withr::with_dir(dir_kegg, {
      for (ext in c(".xml", ".png")) {
        f <- paste0(pid, ext)
        if (file.exists(f)) file.remove(f)
      }
    })

    devs_after <- grDevices::dev.list()
    stray <- setdiff(devs_after, devs_before)
    for (d in stray) {
      grDevices::dev.off(d)
    }

    if (isTRUE(result)) {
      succeeded <- c(succeeded, pid)
    } else {
      failed <- c(failed, pid)
    }
  }

  message(
    "   -> KEGG pathview: ",
    length(succeeded),
    "/",
    length(kegg_ids),
    " pathway map(s) rendered successfully",
    if (length(failed) > 0) {
      paste0(
        "; failed (skipped, others unaffected): ",
        paste(failed, collapse = ", ")
      )
    } else {
      ""
    },
    "."
  )
  invisible(list(succeeded = succeeded, failed = failed))
}

#' Run GSEA for Gene Ontology (BP, MF, CC)
#' @keywords internal
.run_go_gsea <- function(
  gene_list,
  org_db,
  ont,
  pvalue_cutoff,
  out_dir,
  comp_name
) {
  message("   -> Running GO GSEA for ", ont, "...")
  gseago <- safe_run(
    clusterProfiler::gseGO(
      geneList = gene_list,
      OrgDb = .load_org_db(org_db),
      ont = ont,
      keyType = "ENTREZID",
      pvalueCutoff = pvalue_cutoff,
      verbose = FALSE
    ),
    label = paste("GO GSEA", ont)
  )
  if (is.null(gseago) || nrow(gseago@result) == 0) {
    message("      No significant results for GO GSEA ", ont)
    return(NULL)
  }
  dir_go_gsea <- safe_dir(file.path(out_dir, "GSEA", "GO"))
  .write_enrich_csv(
    gseago,
    file.path(dir_go_gsea, paste0("GO_GSEA_", ont, "_", comp_name, ".csv"))
  )
  if (nrow(gseago@result) > 5) {
    safe_pdf(
      file.path(
        dir_go_gsea,
        paste0("GO_GSEA_", ont, "_Dotplot_", comp_name, ".pdf")
      ),
      width = 14,
      height = 14,
      expr = print(enrichplot::dotplot(
        gseago,
        showCategory = 20,
        label_format = 50
      ))
    )
  }
  return(gseago)
}

#' Resolve an MSigDB collection code for the detected organism
#'
#' MSigDB category codes are human-centric ("H", "C1"..."C9"). The native
#' Mouse MSigDB adds a parallel "M" series ("MH", "M1", "M2", "M3", "M5",
#' "M7", "M8" -- there is no native M4/M6/M9). Previously a "C" code was
#' sent to msigdbr as-is regardless of organism, so a mouse run defaulting
#' to e.g. "C2"/"C5"/"C8" would fetch those human-native collections and
#' rely on msigdbr's ortholog mapping via `species =` instead of using the
#' curated mouse-native collection. This maps a requested code to the
#' correct one for `msig_org`: mouse-native when one exists and the
#' organism is mouse, the human/C-style code otherwise (still resolved via
#' `species =` ortholog mapping for organisms/categories with no native
#' collection, e.g. rat, or C4/C6/C9 for mouse).
#'
#' Only ever called with a code already validated against the
#' `valid_collections` whitelist below, so an "M"-prefixed match is always
#' one of the six that actually exist (M1/M2/M3/M5/M7/M8).
#'
#' @param input_cat MSigDB category code, e.g. "C2", "H", "M2", "MH", "HALLMARK"
#' @param msig_org Organism string as returned by get_organism_info()$msig_org
#' @keywords internal
.resolve_msigdbr_collection <- function(input_cat, msig_org) {
  mouse_native <- c("MH", "M1", "M2", "M3", "M5", "M7", "M8")
  is_mouse <- identical(msig_org, "Mus musculus")

  if (input_cat %in% c("H", "MH", "HALLMARK")) {
    return(if (is_mouse) "MH" else "H")
  }

  if (grepl("^M[0-9]", input_cat)) {
    return(if (is_mouse) input_cat else paste0("C", sub("^M", "", input_cat)))
  }

  # Remaining case: a "C"-style code. Only promote it to the mouse-native
  # equivalent if that equivalent actually exists (C4/C6/C9 have none, and
  # fall through to ortholog-mapped "C" for mouse same as any other organism).
  candidate <- paste0("M", sub("^C", "", input_cat))
  if (is_mouse && candidate %in% mouse_native) candidate else input_cat
}

#' Filter a results table down to the ORA "significant gene" input list
#'
#' @description Used identically by GO ORA, Reactome ORA, and Disease
#'   Ontology ORA so all three draw from the same significance definition.
#'   Two knobs are deliberately kept separate from `padj_cutoff` (which
#'   defines "differentially expressed" everywhere else in the pipeline --
#'   GSEA ranking, exported results, plots): ORA is a hypergeometric test on
#'   a fixed gene list and is underpowered by construction when that list is
#'   short, so it benefits from its own, looser threshold.
#'
#' @param tbl A results table (or an Entrez-mapped subset of one) containing
#'   `id_col`, `padj`, `pvalue`, and `log2FoldChange`.
#' @param id_col Name of the identifier column to filter/return by ("gene"
#'   or "entrezid").
#' @param ora_padj_cutoff Adjusted p-value cutoff for the primary filter.
#' @param ora_min_genes Minimum gene count required at `ora_padj_cutoff`
#'   before falling back to a raw `pvalue < 0.05` list. (Previously this
#'   fallback only fired when the padj-filtered list was completely empty --
#'   `nrow(sigOE) == 0` -- so e.g. a 3- or 4-gene list, unlikely to give any
#'   hypergeometric test real power, was never backfilled.)
#' @param ora_lfc_cutoff Optional additional `abs(log2FoldChange) >=` floor,
#'   applied after the p-value filter. `NULL` (default) applies none.
#' @param label Used only in progress messages.
#' @return The subset of `tbl` passing the filter (a data frame, not just an
#'   ID vector, since callers also need other columns e.g. `log2FoldChange`).
#' @keywords internal
.filter_ora_genes <- function(
  tbl,
  id_col,
  ora_padj_cutoff,
  ora_min_genes,
  ora_lfc_cutoff = NULL,
  label = "ORA"
) {
  ids <- as.character(tbl[[id_col]])
  keep <- !is.na(ids) &
    ids != "" &
    !is.na(tbl$padj) &
    tbl$padj < ora_padj_cutoff
  n_padj <- sum(keep)

  if (n_padj < ora_min_genes && "pvalue" %in% colnames(tbl)) {
    message(
      "   -> ",
      label,
      ": only ",
      n_padj,
      " gene(s) at padj < ",
      ora_padj_cutoff,
      " (below ora_min_genes = ",
      ora_min_genes,
      "). Falling back to raw pvalue < 0.05."
    )
    keep <- !is.na(ids) & ids != "" & !is.na(tbl$pvalue) & tbl$pvalue < 0.05
  }

  if (!is.null(ora_lfc_cutoff)) {
    keep <- keep &
      !is.na(tbl$log2FoldChange) &
      abs(tbl$log2FoldChange) >= ora_lfc_cutoff
  }

  out <- tbl[keep, , drop = FALSE]
  message("   -> ", label, " significant genes: ", nrow(out))
  out
}

.reactome_organism <- function(kegg_code) {
  switch(
    kegg_code,
    hsa = "human",
    mmu = "mouse",
    rno = "rat",
    stop("Unsupported Reactome organism code: ", kegg_code, call. = FALSE)
  )
}

#' Run Functional Analysis with Directional Stat Management
#'
#' @export
#' @param res_tbl Results table
#' @param sig_res Significant list
#' @param edb Ensembl database object
#' @param out_dir Output dir
#' @param level Target condition
#' @param base Base condition
#' @param top_genes Limit for GO dotplots
#' @param padj_cutoff Adjusted p-value significance cutoff
#' @param go_pvalue_cutoff GO ORA raw p-value cutoff (default 0.05)
#' @param go_qvalue_cutoff GO ORA q-value cutoff (default 0.2)
#' @param gsea_metric Metric to rank genes for GSEA ("stat", "signed_pval", or "log2FoldChange")
#' @param test_type The upstream test design used: "Wald" or "LRT" (Required for correct stat ranking)
#' @param ora_padj_cutoff Adjusted p-value cutoff for the gene list fed INTO
#'   ORA (GO/Reactome/DO), kept separate from `padj_cutoff` -- see
#'   `.filter_ora_genes()`. Default 0.05 (looser than `padj_cutoff`'s 0.01).
#' @param ora_min_genes Minimum gene count at `ora_padj_cutoff` before
#'   falling back to raw `pvalue < 0.05`; see `.filter_ora_genes()`.
#'   Default 10.
#' @param ora_lfc_cutoff Optional `abs(log2FoldChange) >=` floor applied
#'   uniformly to GO/Reactome/DO ORA input. `NULL` (default) applies none.
#'   (Previously Reactome ORA alone applied an undocumented `abs(LFC) >= 1`
#'   floor, making it stricter than GO/DO ORA for no stated reason.)
#' @param run_gsea Logical: also run the Reactome/KEGG/GO GSEA blocks (and
#'   SPIA, if `run_spia = TRUE`), all of which need a signed per-gene
#'   ranking statistic. Set FALSE for result types with no such statistic
#'   -- e.g. DTU, where DRIMSeq's feature-level output is a directionless
#'   likelihood-ratio stat with no computed proportion-difference/
#'   coefficient column. ORA (GO/Reactome/DO) is unaffected either way.
#'   Default TRUE.
#' @return List with functional results
run_functional_analysis <- function(
  res_tbl,
  sig_res,
  edb,
  out_dir,
  level,
  base,
  top_genes,
  padj_cutoff = 0.01,
  go_pvalue_cutoff = 0.05,
  go_qvalue_cutoff = 0.2,
  gsea_metric = "stat",
  test_type = "Wald",
  ora_padj_cutoff = 0.05,
  ora_min_genes = 10,
  ora_lfc_cutoff = NULL,
  run_gsea = TRUE,
  run_spia = FALSE
) {
  comp_name <- paste0(level, "_vs_", base)
  org_info <- get_organism_info(edb)
  message("Detected organism: ", org_info$name)
  org_db <- org_info$org_db
  org_obj <- .load_org_db(org_db)
  kegg_code <- org_info$kegg_code
  tf_db <- org_info$tf_db

  gseaKEGG <- NULL
  gseaReac <- NULL
  spia_result <- NULL

  if (getOption("ExpressOM.verbose", FALSE)) {
    message("=== DEBUG res_tbl ===")
    message(
      "  Dimensions : ",
      nrow(res_tbl),
      " rows x ",
      ncol(res_tbl),
      " cols"
    )
    message("  Columns    : ", paste(colnames(res_tbl), collapse = ", "))
    for (col in c(
      "gene",
      "log2FoldChange",
      "padj",
      "pvalue",
      "entrezid",
      "stat"
    )) {
      if (col %in% colnames(res_tbl)) {
        vals <- res_tbl[[col]]
        if (is.numeric(vals)) {
          finite_vals <- vals[is.finite(vals)]
          rng <- if (length(finite_vals) > 0) {
            paste(round(range(finite_vals), 4), collapse = " to ")
          } else {
            "no finite values"
          }
          message(
            "  [",
            col,
            "] numeric | NA=",
            sum(is.na(vals)),
            " | range=",
            rng,
            " | head=",
            paste(round(head(finite_vals, 5), 4), collapse = ", ")
          )
        } else {
          non_na <- vals[!is.na(vals) & vals != ""]
          message(
            "  [",
            col,
            "] character | NA/empty=",
            sum(is.na(vals) | vals == ""),
            " | valid=",
            length(non_na),
            " | head=",
            paste(head(non_na, 5), collapse = ", ")
          )
        }
      } else {
        message("  [", col, "] MISSING")
      }
    }
    message("=== END DEBUG ===")
  }

  required_pkgs <- c(
    "ReactomePA",
    "DOSE",
    "enrichplot",
    "enrichR",
    "clusterProfiler",
    "msigdbr"
  )
  if (isTRUE(run_gsea)) {
    required_pkgs <- c(required_pkgs, "pathview")
  }
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "The required package '",
        pkg,
        "' is missing. Please install it to proceed with functional analysis."
      )
    }
  }

  dir_ora <- file.path(out_dir, "ORA")
  dir_gsea <- safe_dir(file.path(out_dir, "GSEA"))

  # DTU has no signed effect-size column at all (DRIMSeq's feature-level
  # output is a directionless likelihood-ratio stat, with no computed
  # proportion-difference/coefficient column). Capture has_lfc BEFORE
  # dummy-filling so res_entrez below (needed for ORA, always runs) doesn't
  # get every row wrongly zeroed out by an `!is.na()` check against a
  # column that was never really there -- and dummy-fill so every
  # downstream `sigOE[..., c("gene","log2FoldChange")]`-style column
  # selection doesn't hard-error with "undefined columns selected".
  has_lfc <- "log2FoldChange" %in% colnames(res_tbl)
  if (!has_lfc) {
    message(
      "   -> No 'log2FoldChange' column in res_tbl (expected for DTU). ",
      "LFC-dependent extras (cnetplot coloring, KEGG pathview, SPIA, ",
      "gsea_metric = \"log2FoldChange\") will be empty/skipped; ORA is unaffected."
    )
    res_tbl$log2FoldChange <- NA_real_
  }

  allOE_genes <- as.character(res_tbl$gene[!is.na(res_tbl$gene)])
  sigOE <- .filter_ora_genes(
    res_tbl,
    "gene",
    ora_padj_cutoff,
    ora_min_genes,
    ora_lfc_cutoff,
    label = "GO ORA"
  )
  sigOE_genes <- unique(as.character(sigOE$gene[
    !is.na(sigOE$gene) & sigOE$gene != ""
  ]))
  if (length(sigOE_genes) == 0) {
    message(
      "No significant genes found for ORA; continuing with ranked GSEA if possible."
    )
  }

  oe_fc_data <- sigOE[
    !is.na(sigOE$log2FoldChange),
    c("gene", "log2FoldChange"),
    drop = FALSE
  ]
  if (nrow(oe_fc_data) > 0) {
    oe_fc <- stats::aggregate(
      log2FoldChange ~ gene,
      data = oe_fc_data,
      FUN = function(x) x[which.max(abs(x))][1]
    )
    OE_foldchanges <- purrr::set_names(oe_fc$log2FoldChange, oe_fc$gene)
    OE_foldchanges <- pmin(pmax(OE_foldchanges, -2), 2)
  } else {
    OE_foldchanges <- numeric(0)
  }

  if (length(sigOE_genes) > 0 && length(tf_db) > 0) {
    message("Running Transcription Factor (TF) enrichment...")
    websiteLive <- getOption("enrichR.live", NA)
    if (is.na(websiteLive)) {
      websiteLive <- tryCatch(
        {
          enrichR::listEnrichrDbs()
          TRUE
        },
        error = function(e) FALSE
      )
      options(enrichR.live = websiteLive)
    }

    if (isTRUE(websiteLive)) {
      message("   Using EnrichR API...")
      enrichr_results <- list()
      for (db in tf_db) {
        res <- safe_run(
          {
            r <- enrichR::enrichr(sigOE_genes, databases = db)
            if (length(r) > 0 && !is.null(r[[1]]) && nrow(r[[1]]) > 0) {
              df <- r[[1]]
              df$Database <- db
              df
            } else {
              NULL
            }
          },
          label = paste("EnrichR TF", db)
        )
        if (!is.null(res)) enrichr_results[[db]] <- res
      }

      if (length(enrichr_results) > 0) {
        df_all <- dplyr::bind_rows(enrichr_results)
        dir_tf <- safe_dir(file.path(out_dir, "Transcription_Factors"))
        write.csv(
          df_all,
          file.path(dir_tf, paste0("Enrichr_TF_", comp_name, ".csv")),
          row.names = FALSE
        )
        message(
          "      EnrichR TF: ",
          length(enrichr_results),
          "/",
          length(tf_db),
          " DBs returned results."
        )
      } else {
        message("      EnrichR TF: No DBs returned valid results.")
      }
    } else {
      msig_org <- org_info$msig_org
      tf_collection <- .resolve_msigdbr_collection("C3", msig_org)
      message(
        "      EnrichR API unreachable; falling back to local TF enrichment using MSigDB ",
        tf_collection,
        " (TFT)."
      )
      message("      Using MSigDB ", tf_collection, " (TFT) for ", msig_org)

      local_res <- safe_run(
        run_local_enrichment(
          gene_list = sigOE_genes,
          universe = allOE_genes,
          organism = msig_org,
          collection = tf_collection,
          subcategory = "TFT",
          pvalue_cutoff = go_pvalue_cutoff,
          qvalue_cutoff = go_qvalue_cutoff
        ),
        label = "Local TF enrichment"
      )

      if (!is.null(local_res) && nrow(as.data.frame(local_res)) > 0) {
        dir_tf <- safe_dir(file.path(out_dir, "Transcription_Factors"))
        write.csv(
          as.data.frame(local_res),
          file.path(
            dir_tf,
            paste0(
              "Local_TF_Enrichment_",
              tf_collection,
              "_",
              comp_name,
              ".csv"
            )
          ),
          row.names = FALSE
        )
        message(
          "      Local TF enrichment completed with ",
          nrow(as.data.frame(local_res)),
          " terms."
        )
        if (nrow(local_res@result) > 5) {
          safe_pdf(
            file.path(dir_tf, paste0("Local_TF_Dotplot_", comp_name, ".pdf")),
            width = 12,
            height = 10,
            expr = print(enrichplot::dotplot(
              local_res,
              showCategory = 20,
              label_format = 50
            ))
          )
        }
      } else {
        message("      No significant TF enrichment found locally.")
      }
    }
  } else if (length(sigOE_genes) == 0) {
    message("Skipping TF enrichment: no significant genes available.")
  } else {
    message("No TF databases specified; skipping TF enrichment.")
  }

  message("Running GO ORA...")
  ego_list <- list()
  has_go_results <- FALSE
  if (length(sigOE_genes) > 0) {
    for (ont in c("BP", "MF", "CC")) {
      ego <- .run_go_ontology(
        ont,
        sigOE_genes,
        allOE_genes,
        org_db,
        go_pvalue_cutoff,
        go_qvalue_cutoff
      )
      if (!is.null(ego)) {
        dir_ora <- safe_dir(dir_ora)
        dir_go <- safe_dir(file.path(dir_ora, "GO"))
        ego_list[[ont]] <- .write_go_ontology_outputs(
          ego,
          ont,
          OE_foldchanges,
          top_genes,
          dir_go,
          comp_name
        )
        has_go_results <- TRUE
      }
    }
    ego_list <- purrr::compact(ego_list)
  } else {
    message("Skipping GO ORA: no significant genes available.")
  }

  if (!"entrezid" %in% colnames(res_tbl)) {
    res_tbl$entrezid <- NA_character_
  }
  missing_entrez <- is.na(res_tbl$entrezid) | res_tbl$entrezid == ""
  if (any(missing_entrez)) {
    message(
      "   -> ",
      sum(missing_entrez),
      "/",
      nrow(res_tbl),
      " genes missing Entrez IDs — mapping gene symbols via AnnotationDbi..."
    )
    if (!is.null(org_obj)) {
      mapped <- suppressMessages(
        AnnotationDbi::mapIds(
          org_obj,
          keys = as.character(res_tbl$gene[missing_entrez]),
          column = "ENTREZID",
          keytype = "SYMBOL",
          multiVals = "first"
        )
      )
      res_tbl$entrezid[missing_entrez] <- mapped[as.character(res_tbl$gene[
        missing_entrez
      ])]
    } else {
      message(
        "   -> WARNING: Could not load org_db object for on-the-fly mapping."
      )
    }
  }

  lfc_ok <- if (has_lfc) !is.na(res_tbl$log2FoldChange) else TRUE
  res_entrez <- res_tbl[
    !is.na(res_tbl$entrezid) & res_tbl$entrezid != "" & lfc_ok,
  ]
  res_entrez <- res_entrez[!duplicated(res_entrez$entrezid), ]

  total_genes <- nrow(res_tbl)
  mapped_genes <- nrow(res_entrez)
  message(
    "   -> Entrez mapping: ",
    mapped_genes,
    "/",
    total_genes,
    " genes have valid Entrez IDs"
  )
  if (mapped_genes == 0) {
    message(
      "   -> WARNING: No Entrez IDs found. Check that org_db (",
      org_db,
      ") is installed and gene symbols match."
    )
  }

  if (isTRUE(run_gsea)) {
    message("   -> Generating ranked list using metric: ", gsea_metric)
    if (gsea_metric == "stat") {
      if (!"stat" %in% colnames(res_entrez)) {
        stop(
          "Column 'stat' not found in results table. Verify that you injected res_unshrunken$stat back into your results."
        )
      }
      if (
        toupper(test_type) == "LRT" &&
          !"contrast_stat" %in% colnames(res_entrez)
      ) {
        stop(
          "LRT GSEA requires a separate contrast-specific statistic; the omnibus LRT statistic is unsigned.",
          call. = FALSE
        )
      }

      if (toupper(test_type) == "LRT") {
        if ("contrast_stat" %in% colnames(res_entrez)) {
          message(
            "   -> LRT detected: ranking by separate contrast-specific Wald statistics."
          )
          metric_vals <- res_entrez$contrast_stat
        } else {
          stop(
            "LRT enrichment ranking requires contrast-specific statistics; omnibus LRT statistics are not signed.",
            call. = FALSE
          )
        }
      } else {
        message(
          "   -> Wald design detected. Using native directional Wald z-scores."
        )
        metric_vals <- res_entrez$stat
      }
    } else if (gsea_metric == "signed_pval") {
      pvalue_col <- if (
        toupper(test_type) == "LRT" &&
          "contrast_pvalue" %in% colnames(res_entrez)
      ) {
        "contrast_pvalue"
      } else {
        "pvalue"
      }
      if (
        toupper(test_type) == "LRT" &&
          !"contrast_pvalue" %in% colnames(res_entrez)
      ) {
        stop(
          "LRT signed-p-value ranking requires contrast-specific p-values.",
          call. = FALSE
        )
      }
      safe_pvals <- ifelse(
        res_entrez[[pvalue_col]] == 0,
        .Machine$double.xmin,
        res_entrez[[pvalue_col]]
      )
      metric_vals <- sign(res_entrez$log2FoldChange) * -log10(safe_pvals)
    } else {
      metric_vals <- res_entrez$log2FoldChange
    }
    set.seed(123456)
    metric_vals <- metric_vals + runif(nrow(res_entrez), -1e-9, 1e-9)
    gsea_list <- sort(
      purrr::set_names(metric_vals, as.character(res_entrez$entrezid)),
      decreasing = TRUE
    )
  } else {
    message(
      "   -> Skipping GSEA-family ranking (run_gsea = FALSE): no valid signed statistic for this result type (e.g. DTU)."
    )
    gsea_list <- NULL
  }

  # GO ORA, Reactome ORA, and Disease Ontology ORA all draw from the same
  # significant-gene definition via .filter_ora_genes(). Previously Reactome
  # alone applied an extra, hardcoded abs(log2FC) >= 1 floor on top of the
  # padj/pvalue filter (reac_lfc_cutoff <- 1.0), making it strictly harder to
  # get a Reactome hit than a GO or DO hit for no stated reason. Set
  # ora_lfc_cutoff to reintroduce an LFC floor -- uniformly, across all
  # three -- if desired.
  sig_entrez_df <- .filter_ora_genes(
    res_entrez,
    "entrezid",
    ora_padj_cutoff,
    ora_min_genes,
    ora_lfc_cutoff,
    label = "Reactome/DO ORA"
  )
  sig_entrez_ids <- as.character(sig_entrez_df$entrezid)
  sig_entrez_reac <- sig_entrez_ids

  message("Running Reactome ORA on curated list...")
  reac_org <- .reactome_organism(kegg_code)
  if (length(sig_entrez_reac) > 0) {
    x <- safe_run(
      suppressMessages(ReactomePA::enrichPathway(
        gene = sig_entrez_reac,
        universe = as.character(res_entrez$entrezid),
        organism = reac_org,
        pvalueCutoff = go_pvalue_cutoff,
        qvalueCutoff = go_qvalue_cutoff,
        readable = TRUE
      )),
      label = "Reactome ORA"
    )
    if (!is.null(x) && nrow(as.data.frame(x)) > 0) {
      dir_ora <- safe_dir(dir_ora)
      dir_reac_ora <- safe_dir(file.path(dir_ora, "Reactome"))
      .write_enrich_csv(
        x,
        file.path(dir_reac_ora, paste0("Reactome_ORA_", comp_name, ".csv"))
      )
      if (nrow(x@result) > 5) {
        safe_pdf(
          file.path(
            dir_reac_ora,
            paste0("Reactome_ORA_Dotplot_", comp_name, ".pdf")
          ),
          width = 14,
          height = 14,
          expr = print(enrichplot::dotplot(
            x,
            showCategory = top_genes,
            label_format = 50,
            color = "pvalue"
          ))
        )

        safe_pdf(
          file.path(
            dir_reac_ora,
            paste0("Reactome_ORA_Cnetplot_", comp_name, ".pdf")
          ),
          width = 12,
          height = 12,
          expr = print(enrichplot::cnetplot(
            x,
            foldChange = OE_foldchanges,
            showCategory = top_genes
          ))
        )
      }
    } else {
      message("      No significant results for Reactome ORA")
    }
  }

  message("Running Disease Ontology ORA...")
  if (kegg_code == "hsa" && length(sig_entrez_ids) > 0) {
    do_res <- .run_disease_ontology(
      sig_entrez_ids,
      universe_entrez = as.character(res_entrez$entrezid),
      go_pvalue_cutoff,
      go_qvalue_cutoff
    )
    if (!is.null(do_res)) {
      dir_ora <- safe_dir(dir_ora)
      .write_disease_ontology_outputs(do_res, top_genes, dir_ora, comp_name)
    }
  } else {
    message("      Skipping Disease Ontology ORA (not human or no sig genes)")
  }

  if (isTRUE(run_gsea)) {
    message("Running Reactome GSEA on complete ranked genome background...")
    set.seed(123456)
    gseaReac <- safe_run(
      suppressMessages(ReactomePA::gsePathway(
        geneList = gsea_list,
        organism = reac_org,
        minGSSize = 10,
        maxGSSize = 300,
        pvalueCutoff = padj_cutoff,
        pAdjustMethod = "BH",
        verbose = FALSE
      )),
      label = "Reactome GSEA"
    )

    if (!is.null(gseaReac) && nrow(as.data.frame(gseaReac)) > 0) {
      dir_reac_gsea <- safe_dir(file.path(dir_gsea, "Reactome"))

      if (nrow(gseaReac@result) > 5) {
        safe_pdf(
          file.path(
            dir_reac_gsea,
            paste0("Reactome_GSEA_Ridgeplot_", comp_name, ".pdf")
          ),
          width = 14,
          height = 14,
          expr = {
            p <- enrichplot::ridgeplot(
              gseaReac,
              showCategory = top_genes,
              label_format = 50
            )
            print(p)
          }
        )
      }

      gseaReac <- suppressMessages(clusterProfiler::setReadable(
        gseaReac,
        OrgDb = org_obj,
        keyType = "ENTREZID"
      ))
      .write_enrich_csv(
        gseaReac,
        file.path(dir_reac_gsea, paste0("Reactome_GSEA_", comp_name, ".csv"))
      )

      if (nrow(gseaReac@result) > 0) {
        safe_pdf(
          file.path(
            dir_reac_gsea,
            paste0("Reactome_GSEA_Dotplot_", comp_name, ".pdf")
          ),
          width = 14,
          height = 14,
          expr = {
            p <- enrichplot::dotplot(
              gseaReac,
              showCategory = top_genes,
              label_format = 50,
              color = "p.adjust"
            )
            print(p)
          }
        )

        safe_pdf(
          file.path(
            dir_reac_gsea,
            paste0("Reactome_GSEA_Gseaplot_", comp_name, ".pdf")
          ),
          width = 12,
          height = 8,
          expr = print(enrichplot::gseaplot2(
            gseaReac,
            geneSetID = 1:min(3, nrow(gseaReac@result))
          ))
        )
      }
    } else {
      message("      No significant results for Reactome GSEA")
    }

    message("Running KEGG GSEA...")
    set.seed(123456)
    gseaKEGG <- safe_run(
      clusterProfiler::gseKEGG(
        geneList = gsea_list,
        organism = kegg_code,
        minGSSize = 5,
        pvalueCutoff = padj_cutoff,
        verbose = FALSE
      ),
      label = "KEGG GSEA"
    )

    if (!is.null(gseaKEGG) && nrow(as.data.frame(gseaKEGG)) > 0) {
      disease_pattern <- paste0("^", kegg_code, "05")
      gseaKEGG <- dplyr::filter(gseaKEGG, !grepl(disease_pattern, ID))

      if (nrow(as.data.frame(gseaKEGG)) == 0) {
        message("      No non-disease KEGG pathways remained after filtering.")
        gseaKEGG <- NULL
      } else {
        message(
          "      Filtered out KEGG Disease pathways. Remaining: ",
          nrow(as.data.frame(gseaKEGG))
        )
      }
    }

    if (!is.null(gseaKEGG) && nrow(as.data.frame(gseaKEGG)) > 0) {
      dir_kegg <- safe_dir(file.path(dir_gsea, "KEGG"))

      if (nrow(gseaKEGG@result) > 5) {
        safe_pdf(
          file.path(
            dir_kegg,
            paste0("KEGG_GSEA_Ridgeplot_", comp_name, ".pdf")
          ),
          width = 14,
          height = 14,
          expr = print(enrichplot::ridgeplot(
            gseaKEGG,
            showCategory = top_genes,
            label_format = 50
          ))
        )
      }

      gseaKEGG <- suppressMessages(clusterProfiler::setReadable(
        gseaKEGG,
        OrgDb = org_obj,
        keyType = "ENTREZID"
      ))
      .write_enrich_csv(
        gseaKEGG,
        file.path(dir_kegg, paste0("KEGG_GSEA_", comp_name, ".csv"))
      )

      if (nrow(gseaKEGG@result) > 5) {
        safe_pdf(
          file.path(dir_kegg, paste0("KEGG_GSEA_Dotplot_", comp_name, ".pdf")),
          width = 14,
          height = 14,
          expr = print(enrichplot::dotplot(
            gseaKEGG,
            showCategory = top_genes,
            label_format = 50
          ))
        )
      }

      message("Generating KEGG pathway maps...")
      all_lfc_vector <- purrr::set_names(
        res_entrez$log2FoldChange,
        as.character(res_entrez$entrezid)
      )
      .run_kegg_pathview(
        gseaKEGG,
        all_lfc_vector,
        kegg_code,
        top_genes,
        dir_kegg
      )
    } else {
      message("      No significant KEGG pathways found in GSEA.")
    }

    message("Running GO GSEA...")
    for (ont in c("BP", "MF", "CC")) {
      .run_go_gsea(gsea_list, org_db, ont, padj_cutoff, out_dir, comp_name)
    }
  } else {
    message(
      "Skipping Reactome/KEGG/GO GSEA (run_gsea = FALSE): no valid signed statistic for this result type (e.g. DTU)."
    )
  }

  if (isTRUE(run_spia) && isTRUE(run_gsea)) {
    message("Running SPIA analysis...")
    if (length(sig_entrez_ids) > 0) {
      spia_de <- purrr::set_names(
        as.numeric(res_entrez$log2FoldChange),
        as.character(res_entrez$entrezid)
      )
      spia_de <- spia_de[names(spia_de) %in% sig_entrez_ids]
      spia_de <- spia_de[!is.na(spia_de) & !duplicated(names(spia_de))]

      message(
        "   -> SPIA: ",
        length(spia_de),
        " DE genes, ",
        length(gsea_list),
        " background genes"
      )
      spia_result <- safe_run(
        SPIA::spia(
          de = spia_de,
          all = as.character(res_entrez$entrezid),
          organism = kegg_code,
          plots = FALSE
        ),
        label = "SPIA"
      )

      if (!is.null(spia_result) && nrow(spia_result) > 0) {
        dir_spia <- safe_dir(file.path(out_dir, "SPIA"))
        write.csv(
          spia_result,
          file.path(dir_spia, paste0("SPIA_Results_", comp_name, ".csv")),
          row.names = FALSE
        )
        message("   -> SPIA: ", nrow(spia_result), " pathways saved to CSV.")

        safe_pdf(
          file.path(dir_spia, paste0("SPIA_Evidence_", comp_name, ".pdf")),
          width = 8,
          height = 8,
          expr = plotP_fork(spia_result, threshold = padj_cutoff)
        )
      } else if (!is.null(spia_result) && nrow(spia_result) == 0) {
        message("      SPIA returned an empty result (no pathways perturbed).")
      } else {
        message("      SPIA failed — check Log/Warnings.txt for details.")
      }
    } else {
      message("      Skipping SPIA: no significant Entrez genes available.")
    }
  } else if (isTRUE(run_spia)) {
    message(
      "Skipping SPIA: run_gsea = FALSE means no directional statistic is available for this result type (e.g. DTU)."
    )
  } else {
    message("Skipping SPIA (run_spia = FALSE).")
  }

  message("Functional analysis complete.")
  list(
    ego_list = ego_list,
    spia_result = spia_result,
    gseaKEGG = gseaKEGG,
    gseaReac = gseaReac
  )
}

#' Run fgsea analysis using DE results and a GMT file (or multiple GMT files)
#'
#' @description Performs Gene Set Enrichment Analysis using both `fgsea` and
#'   `clusterProfiler` using external GMT files (e.g., from MSigDB) or
#'   downloading via msigdbr.
#'
#' @param res_tbl The full results table from `export_significant_results()`
#'   containing 'gene' and 'log2FoldChange'.
#' @param gmt_file Path to a local `.gmt` file, a character vector/list of
#'   multiple `.gmt` files, or MSigDB category names (e.g., "H", "C2") to
#'   download via msigdbr. Category names are resolved per-organism (see
#'   `edb`): for mouse, "C" codes are automatically promoted to their
#'   native Mouse MSigDB "M" equivalent where one exists (e.g. "C2" ->
#'   "M2", "H" -> "MH") via `.resolve_msigdbr_collection()`, rather than
#'   being sent to msigdbr as human-native collections regardless of
#'   species.
#' @param edb Ensembl Database (used to detect species for msigdbr).
#' @param out_dir Output directory for fgsea results.
#' @param comp_name Comparison name for file naming.
#' @param padj_cutoff Adjusted p-value cutoff for filtering significant pathways.
#' @param gsea_metric Metric to rank genes for GSEA ("stat", "signed_pval", or
#'   "log2FoldChange"). Same semantics as `run_functional_analysis()`'s
#'   `gsea_metric` -- kept as a separate argument here (rather than silently
#'   inherited) because this function ranks by gene SYMBOL, not Entrez ID.
#' @param test_type The upstream test design used: "Wald" or "LRT". Only
#'   affects ranking when `gsea_metric = "stat"` (see `run_functional_analysis()`).
#'
#' @export
run_fgsea_analysis <- function(
  res_tbl,
  gmt_file = c("C2", "C5", "C8"),
  edb,
  out_dir,
  comp_name,
  padj_cutoff = 0.01,
  gsea_metric = "stat",
  test_type = "Wald"
) {
  gmt_list <- if (is.null(gmt_file)) list(NULL) else as.list(gmt_file)

  # .safe_ggsave() is defined once, package-level, in utils_core.R.

  # Built once, outside the gmt_list loop, since the ranking does not depend
  # on which pathway collection is being tested (previously this identical
  # computation -- including the fixed jitter seed -- was silently repeated
  # once per collection, e.g. 10x with the default `gmt_file`).
  message(
    "-> Preparing ranked gene list for FGSEA (metric: ",
    gsea_metric,
    ")..."
  )

  res2 <- res_tbl |>
    dplyr::select(
      gene,
      log2FoldChange,
      dplyr::any_of(c("stat", "pvalue", "contrast_stat", "contrast_pvalue"))
    ) |>
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FoldChange)) |>
    dplyr::distinct()

  if (gsea_metric == "stat") {
    if (!"stat" %in% colnames(res2)) {
      stop(
        "Column 'stat' not found in results table. Verify that you injected ",
        "res_unshrunken$stat back into your results, or call run_fgsea_analysis() ",
        "with gsea_metric = \"log2FoldChange\" or \"signed_pval\" instead."
      )
    }
    # One row per gene symbol, keeping whichever duplicate has the most
    # extreme |stat| (mirrors the OE_foldchanges convention used elsewhere
    # in this file), so the log2FoldChange used for the LRT sign transform
    # below always comes from the same row as the stat it signs.
    res2 <- res2 |>
      dplyr::group_by(gene) |>
      dplyr::slice_max(abs(stat), n = 1, with_ties = FALSE) |>
      dplyr::ungroup()

    if (toupper(test_type) == "LRT") {
      if (!"contrast_stat" %in% colnames(res2)) {
        stop(
          "LRT enrichment ranking requires contrast-specific statistics; omnibus LRT statistics are not signed.",
          call. = FALSE
        )
      }
      message(
        "   -> LRT detected: ranking by separate contrast-specific Wald statistics."
      )
      metric_vals <- res2$contrast_stat
    } else {
      metric_vals <- res2$stat
    }
  } else if (gsea_metric == "signed_pval") {
    if (!"pvalue" %in% colnames(res2)) {
      stop(
        "Column 'pvalue' not found in results table. Cannot use gsea_metric = \"signed_pval\"."
      )
    }
    if (
      toupper(test_type) == "LRT" &&
        !"contrast_pvalue" %in% colnames(res2)
    ) {
      stop(
        "LRT signed-p-value ranking requires contrast-specific p-values.",
        call. = FALSE
      )
    }
    res2 <- res2 |>
      dplyr::group_by(gene) |>
      dplyr::slice_min(pvalue, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()

    pvalue_col <- if (
      toupper(test_type) == "LRT" &&
        "contrast_pvalue" %in% colnames(res2)
    ) {
      "contrast_pvalue"
    } else {
      "pvalue"
    }
    safe_pvals <- ifelse(
      res2[[pvalue_col]] == 0,
      .Machine$double.xmin,
      res2[[pvalue_col]]
    )
    metric_vals <- sign(res2$log2FoldChange) * -log10(safe_pvals)
  } else {
    # Unchanged from the original implementation: mean log2FoldChange
    # across duplicate gene symbols.
    res2 <- res2 |>
      dplyr::group_by(gene) |>
      dplyr::summarize(
        log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
        .groups = "drop"
      )
    metric_vals <- res2$log2FoldChange
  }

  ranks <- purrr::set_names(metric_vals, as.character(res2$gene))

  if (length(ranks) == 0) {
    message("-> Skipping FGSEA: no valid ranked genes available.")
    return(invisible(NULL))
  }

  set.seed(123456)
  ranks <- ranks + stats::runif(length(ranks), min = -1e-6, max = 1e-6)
  ranks <- sort(ranks, decreasing = TRUE)

  for (gmt_item in gmt_list) {
    if (is.null(gmt_item) || !file.exists(gmt_item)) {
      if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop(
          "Package 'msigdbr' is required to download pathways. Please install it."
        )
      }

      org_info <- get_organism_info(edb)
      msig_org <- org_info$msig_org

      valid_collections <- c(
        "H",
        "C1",
        "C2",
        "C3",
        "C4",
        "C5",
        "C6",
        "C7",
        "C8",
        "C9",
        "MH",
        "M1",
        "M2",
        "M3",
        "M5",
        "M7",
        "M8",
        "HALLMARK"
      )

      raw_input <- toupper(gmt_item)
      has_subcat <- grepl(":", raw_input)
      input_cat <- if (has_subcat) sub(":.*$", "", raw_input) else raw_input
      input_subcat <- if (has_subcat) sub("^[^:]+:", "", raw_input) else NULL

      if (!is.null(gmt_item) && input_cat %in% valid_collections) {
        msig_cat <- .resolve_msigdbr_collection(input_cat, msig_org)

        # Named after the resolved collection (msig_cat), not the raw input,
        # so mouse runs get e.g. "msigdbr_M2..." output folders/files rather
        # than a "C2" name that would misrepresent what's actually inside.
        gmt_name <- paste0(
          "msigdbr_",
          msig_cat,
          if (!is.null(input_subcat)) {
            paste0("_", gsub(":", "_", input_subcat))
          } else {
            ""
          }
        )

        message(
          "-> Fetching MSigDB category [",
          input_cat,
          "]",
          if (!is.null(input_subcat)) {
            paste0(" subcategory [", input_subcat, "]")
          } else {
            ""
          },
          " (mapped to ",
          msig_cat,
          ") via msigdbr for ",
          msig_org,
          "..."
        )
      } else {
        message(
          "Provided GMT file '",
          gmt_item,
          "' not recognized as MSigDB category. Falling back to Hallmark..."
        )

        msig_cat <- .resolve_msigdbr_collection("H", msig_org)
        input_subcat <- NULL
        gmt_name <- "hallmark_msigdbr"
      }

      .msigdbr_fetch <- function(species, cat, subcat = NULL) {
        fmls <- names(formals(msigdbr::msigdbr))

        args <- list(species = species)

        if ("collection" %in% fmls) {
          args$collection <- cat
        } else {
          args$category <- cat
        }

        if (!is.null(subcat)) {
          if ("subcollection" %in% fmls) {
            args$subcollection <- subcat
          } else {
            args$subcategory <- subcat
          }
        }

        # Mouse-native collections (MH/M1/M2/M3/M5/M7/M8) live in msigdbr's
        # "MM" database division. msigdbr defaults db_species to "HS", under
        # which those collection codes don't exist and the call errors out
        # (e.g. collection = "M2" against the default "HS" division) --
        # this is what was actually causing FGSEA to fail outright on mouse
        # runs rather than just falling back to Hallmark. The collection
        # code itself tells us which division it belongs to regardless of
        # the requested output species (e.g. cat = "C4" for a mouse run
        # correctly still queries "HS", since C4 has no native mouse
        # collection -- see .resolve_msigdbr_collection()).
        if ("db_species" %in% fmls) {
          args$db_species <- if (grepl("^M[H0-9]", cat)) "MM" else "HS"
        }

        do.call(msigdbr::msigdbr, args) |>
          dplyr::select(gs_name, gene_symbol)
      }

      m_t2g <- tryCatch(
        .msigdbr_fetch(msig_org, msig_cat, input_subcat),
        error = function(e) {
          hallmark_cat <- .resolve_msigdbr_collection("H", msig_org)
          message(
            "   -> msigdbr collection '",
            msig_cat,
            if (!is.null(input_subcat)) {
              paste0(" / subcategory '", input_subcat, "'")
            } else {
              ""
            },
            "' not available, falling back to '",
            hallmark_cat,
            "'"
          )

          tryCatch(
            .msigdbr_fetch(msig_org, hallmark_cat),
            error = function(e2) {
              stop(
                "msigdbr failed for both '",
                msig_cat,
                "' and '",
                hallmark_cat,
                "': ",
                e2$message
              )
            }
          )
        }
      )

      if (!is.null(input_subcat) && nrow(m_t2g) == 0) {
        message(
          "   -> msigdbr subcategory '",
          input_subcat,
          "' returned 0 gene sets for '",
          msig_cat,
          "' / ",
          msig_org,
          "; using the full '",
          msig_cat,
          "' collection instead."
        )

        m_t2g <- tryCatch(
          .msigdbr_fetch(msig_org, msig_cat),
          error = function(e) m_t2g
        )
      }

      if (nrow(m_t2g) == 0) {
        message(
          "-> Skipping [",
          gmt_name,
          "]: msigdbr returned 0 gene sets."
        )
        next
      }

      pathways.GSEA <- m_t2g
      pathways.fgsea <- split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)
    } else {
      gmt_name <- tools::file_path_sans_ext(basename(gmt_item))

      message("-> Loading pathways for fgsea/GSEA: ", gmt_name)

      pathways.GSEA <- clusterProfiler::read.gmt(gmt_item)
      pathways.fgsea <- fgsea::gmtPathways(gmt_item)
    }

    if (length(pathways.fgsea) == 0) {
      message("-> Skipping [", gmt_name, "]: no pathways available.")
      next
    }

    fgsea_out <- file.path(out_dir, "GSEA", "FGSEA", gmt_name)

    if (!dir.exists(fgsea_out)) {
      dir.create(fgsea_out, recursive = TRUE)
    }

    message("-> Running clusterProfiler GSEA for [", gmt_name, "]...")

    set.seed(123456)

    gsea_results <- tryCatch(
      {
        suppressWarnings(
          suppressMessages(
            clusterProfiler::GSEA(
              ranks,
              TERM2GENE = pathways.GSEA,
              pvalueCutoff = padj_cutoff
            )
          )
        )
      },
      error = function(e) {
        message("   -> clusterProfiler GSEA failed: ", conditionMessage(e))
        NULL
      }
    )

    message("-> Running fgseaMultilevel for [", gmt_name, "]...")

    set.seed(123456)

    fgseaRes <- tryCatch(
      fgsea::fgseaMultilevel(
        pathways = pathways.fgsea,
        stats = ranks,
        # fgseaMultilevel() defaults to minSize=1/maxSize=length(stats)-1 --
        # i.e. unbounded -- unlike every other enrichment call in this file
        # (GO ORA/GSEA and the clusterProfiler::GSEA() call just above default
        # to 10/500; Reactome GSEA uses 10/300). Left unbounded, fgsea tests
        # and reports tiny/huge gene sets that ORA would never consider,
        # which on its own inflates the apparent hit count relative to ORA.
        minSize = 10,
        maxSize = 500
      ),
      error = function(e) {
        message("   -> fgseaMultilevel failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(fgseaRes)) {
      message(
        "-> Skipping [",
        gmt_name,
        "]: fgseaMultilevel returned no results."
      )
      next
    }

    message("-> Saving results table for [", gmt_name, "]...")

    fgseaResTidy <- tibble::as_tibble(fgseaRes) |>
      dplyr::arrange(dplyr::desc(NES))

    utils::write.csv(
      fgseaResTidy |> dplyr::select(-leadingEdge),
      file.path(fgsea_out, paste0("FGSEA_Results_", gmt_name, ".csv")),
      row.names = FALSE
    )

    # unlike enrichGO()/enricher() (which only ever return terms that already
    # pass pvalueCutoff/qvalueCutoff), fgseaMultilevel() has no significance
    # filter -- the CSV above lists every pathway tested. Export a
    # significant-only companion CSV so it's directly comparable to the ORA
    # CSVs in ORA/GO, ORA/Reactome, etc. (row count vs row count) instead of
    # "everything tested" vs "already-significant-only".
    fgseaResTidy_sig <- fgseaResTidy |> dplyr::filter(padj < padj_cutoff)
    utils::write.csv(
      fgseaResTidy_sig |> dplyr::select(-leadingEdge),
      file.path(
        fgsea_out,
        paste0("FGSEA_Results_", gmt_name, "_significant.csv")
      ),
      row.names = FALSE
    )
    message(
      "   -> ",
      nrow(fgseaResTidy_sig),
      "/",
      nrow(fgseaResTidy),
      " pathways significant at padj < ",
      padj_cutoff
    )

    if (
      requireNamespace("DT", quietly = TRUE) &&
        requireNamespace("htmlwidgets", quietly = TRUE)
    ) {
      tryCatch(
        {
          datatable_object <- fgseaResTidy |>
            dplyr::select(-leadingEdge, -ES) |>
            dplyr::arrange(padj) |>
            DT::datatable()

          html_name <- paste0("GSEA_Table_", gmt_name, ".html")

          withr::with_dir(
            fgsea_out,
            {
              htmlwidgets::saveWidget(
                datatable_object,
                file = html_name,
                selfcontained = TRUE
              )
            }
          )
        },
        error = function(e) {
          message(
            "   -> Skipping interactive HTML table for [",
            gmt_name,
            "]: ",
            conditionMessage(e)
          )
        }
      )
    } else {
      message(
        "Skipping HTML table generation: 'DT' or 'htmlwidgets' is not installed."
      )
    }

    message("-> Plotting NES barplot for [", gmt_name, "]...")

    fgseaResTidy_filtered <- fgseaResTidy_sig

    if (nrow(fgseaResTidy_filtered) > 0) {
      max_barplot_pathways <- getOption(
        "ExpressOM.max_gsea_barplot_pathways",
        80L
      )
      max_barplot_pathways <- as.integer(max_barplot_pathways)

      if (is.na(max_barplot_pathways) || max_barplot_pathways < 5) {
        max_barplot_pathways <- 80L
      }

      total_sig <- nrow(fgseaResTidy_filtered)

      fgseaResTidy_filtered <- fgseaResTidy_filtered[
        order(
          fgseaResTidy_filtered$padj,
          -abs(fgseaResTidy_filtered$NES)
        ),
      ]

      barplot_subtitle <- NULL

      if (total_sig > max_barplot_pathways) {
        message(
          "   -> Limiting NES barplot to top ",
          max_barplot_pathways,
          " of ",
          total_sig,
          " significant pathways.",
          " Set options(ExpressOM.max_gsea_barplot_pathways = N) to change."
        )

        fgseaResTidy_filtered <- utils::head(
          fgseaResTidy_filtered,
          max_barplot_pathways
        )

        barplot_subtitle <- sprintf(
          "Showing top %d of %d significant pathways (ordered by padj)",
          max_barplot_pathways,
          total_sig
        )
      }

      fgseaResTidy_filtered$pathway <- gsub(
        "_",
        " ",
        fgseaResTidy_filtered$pathway
      )

      fgseaResTidy_filtered$pathway <- stringr::str_wrap(
        fgseaResTidy_filtered$pathway,
        width = 60
      )

      fgseaResTidy_filtered$direction <- ifelse(
        fgseaResTidy_filtered$NES > 0,
        "Up",
        "Down"
      )

      p_bar <- ggplot2::ggplot(
        fgseaResTidy_filtered,
        ggplot2::aes(reorder(pathway, NES), NES)
      ) +
        ggplot2::geom_col(
          ggplot2::aes(fill = direction),
          width = 0.6
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          x = "Pathway",
          y = "Normalized Enrichment Score",
          title = paste("GSEA NES:", gmt_name),
          subtitle = barplot_subtitle
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 9, color = "grey30"),
          axis.title.x = ggplot2::element_text(face = "bold"),
          axis.title.y = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(face = "bold"),
          axis.text.y = ggplot2::element_text(face = "bold", size = 8),
          plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")
        ) +
        ggplot2::scale_fill_manual(
          values = c(
            "Up" = .de_direction_colors()[["up"]],
            "Down" = .de_direction_colors()[["down"]]
          ),
          drop = FALSE
        ) +
        ggplot2::guides(fill = "none")

      num_pathways <- nrow(fgseaResTidy_filtered)

      ## Never allow height to exceed 30 inches.
      height_adjustment <- max(8, min(30, num_pathways * 0.3 + 2))

      .safe_ggsave(
        filename = file.path(
          fgsea_out,
          paste0("Barplot_", gmt_name, ".pdf")
        ),
        plot = p_bar,
        width = 12,
        height = height_adjustment
      )
    } else {
      message(
        "   -> No significant pathways at padj < ",
        padj_cutoff,
        " for NES barplot."
      )
    }

    message(
      "-> Generating individual pathway gseaplots for [",
      gmt_name,
      "]..."
    )

    max_gsea_plots <- getOption("ExpressOM.max_gsea_plots", 20L)
    max_gsea_plots <- as.integer(max_gsea_plots)

    if (is.na(max_gsea_plots) || max_gsea_plots < 0) {
      max_gsea_plots <- 20L
    }

    if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
      gene_set_ids <- gsea_results@result$ID
      valid_idx <- which(!is.na(gene_set_ids))

      if (max_gsea_plots == 0) {
        message(
          "   -> Individual gseaplots disabled by options(ExpressOM.max_gsea_plots = 0)."
        )
        valid_idx <- integer(0)
      } else if (length(valid_idx) > max_gsea_plots) {
        message(
          "   -> Limiting individual gseaplots to top ",
          max_gsea_plots,
          " of ",
          length(valid_idx),
          " (set options(ExpressOM.max_gsea_plots = N) to change)"
        )

        valid_idx <- valid_idx[seq_len(max_gsea_plots)]
      }

      for (i in valid_idx) {
        pid <- gene_set_ids[i]
        clean_pid <- gsub("[^A-Za-z0-9_-]", "_", pid)

        safe_run(
          {
            pathway_name <- gsub(
              "_",
              " ",
              gsea_results@result$Description[i]
            )

            p <- enrichplot::gseaplot2(
              gsea_results,
              geneSetID = i,
              title = as.character(pathway_name)
            )

            .safe_ggsave(
              filename = file.path(
                fgsea_out,
                paste0("GSEA_plot_", gmt_name, "_", clean_pid, ".pdf")
              ),
              plot = p,
              width = 8,
              height = 6
            )

            .safe_ggsave(
              filename = file.path(
                fgsea_out,
                paste0("GSEA_plot_", gmt_name, "_", clean_pid, ".tiff")
              ),
              plot = p,
              device = "tiff",
              width = 8,
              height = 6,
              compression = "lzw",
              dpi = 600
            )
          },
          label = paste("gseaplot", pid)
        )
      }
    } else {
      message(
        "   -> Skipping individual gseaplots: clusterProfiler GSEA returned no significant results."
      )
    }

    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
  }

  message("-> FGSEA processing complete across all specified gene sets.")
}


#' Run local MSigDB-based enrichment analysis (offline fallback for EnrichR)
#'
#' Single source of truth for "gene set overrepresentation via msigdbr +
#' clusterProfiler::enricher()". Used both as a standalone export and as the
#' offline TF-enrichment fallback inside \code{run_functional_analysis()}
#' when the EnrichR API is unreachable.
#'
#' @param gene_list Vector of significant gene symbols (e.g. DEGs)
#' @param universe Vector of ALL genes expressed in the experiment (background)
#' @param organism "Homo sapiens" or "Mus musculus" (passed to \code{msigdbr::msigdbr(species=)})
#' @param collection MSigDB collection, e.g. "H" (Hallmark), "C3" (motif/TFT targets), "C5" (GO)
#' @param subcategory Optional MSigDB subcategory within \code{collection}
#'   (e.g. "TFT" within "C3"). If fetching with the subcategory fails or
#'   returns no gene sets, silently falls back to the full collection.
#' @param pvalue_cutoff,qvalue_cutoff Passed to \code{clusterProfiler::enricher()}
#' @return A clusterProfiler \code{enrichResult} object, or \code{NULL} if no
#'   gene sets were available for the requested collection/organism
#' @export
run_local_enrichment <- function(
  gene_list,
  universe,
  organism = "Homo sapiens",
  collection = "C3",
  subcategory = NULL,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2
) {
  for (pkg in c("clusterProfiler", "msigdbr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '",
        pkg,
        "' is required. Install with: BiocManager::install('",
        pkg,
        "')"
      )
    }
  }

  .fetch <- function(subcat) {
    fmls <- names(formals(msigdbr::msigdbr))
    args <- list(species = organism, collection = collection)
    if (!is.null(subcat)) {
      args$subcategory <- subcat
    }
    # Same db_species requirement as .msigdbr_fetch() in run_fgsea_analysis()
    # -- mouse-native collections (M1/M2/M3/M5/M7/M8/MH) only exist under
    # msigdbr's "MM" database division, not the default "HS".
    if ("db_species" %in% fmls) {
      args$db_species <- if (grepl("^M[H0-9]", collection)) "MM" else "HS"
    }
    do.call(msigdbr::msigdbr, args)
  }

  m_df <- tryCatch(.fetch(subcategory), error = function(e) NULL)
  if ((is.null(m_df) || nrow(m_df) == 0) && !is.null(subcategory)) {
    message(
      "  '",
      collection,
      "' subcategory '",
      subcategory,
      "' not available; using the full '",
      collection,
      "' collection."
    )
    m_df <- tryCatch(.fetch(NULL), error = function(e) NULL)
  }
  if (is.null(m_df) || nrow(m_df) == 0) {
    message(
      "  No gene sets found in MSigDB '",
      collection,
      "' for organism ",
      organism
    )
    return(NULL)
  }

  term2gene <- data.frame(
    term = m_df$gs_name,
    gene = m_df$gene_symbol,
    stringsAsFactors = FALSE
  )

  clusterProfiler::enricher(
    gene = gene_list,
    universe = universe,
    TERM2GENE = term2gene,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )
}

#' Adapt a DTE or DTU results table to the DGE-style res_tbl shape expected
#' by `run_functional_analysis()` / `run_fgsea_analysis()`
#'
#' @description Both DTE's `res_df` and DTU's `dtu_results` already carry
#'   the columns those two functions need (`gene`, `padj`, `pvalue`,
#'   `entrezid`, and for DTE only `log2FoldChange`/`stat`) -- multiple
#'   transcripts per gene are already handled by the existing group-by-gene
#'   dedup logic in both functions. Two adjustments are still needed:
#'   \itemize{
#'     \item Both tables' `gene` column falls back to the transcript/feature
#'       ID via `.coalesce_gene_label()` when no symbol was resolved. Left
#'       as-is, those rows would pad ORA's background/universe with IDs
#'       that can never match a pathway member -- use `gene_symbol` and
#'       drop unresolved rows instead.
#'     \item DTU's adjusted p-value column is named `adj_pvalue`, and its
#'       raw p-value column is whichever of `pvalue`/`p_value` DRIMSeq's
#'       `results()` produced -- normalize both to `padj`/`pvalue`.
#'   }
#' @keywords internal
.adapt_isoform_res_tbl <- function(res_tbl, analysis_type) {
  res_tbl <- as.data.frame(res_tbl)

  if ("gene_symbol" %in% colnames(res_tbl)) {
    res_tbl$gene <- res_tbl$gene_symbol
  }
  res_tbl <- res_tbl[!is.na(res_tbl$gene) & res_tbl$gene != "", , drop = FALSE]

  if (identical(analysis_type, "DTU")) {
    if (!"padj" %in% colnames(res_tbl) && "adj_pvalue" %in% colnames(res_tbl)) {
      res_tbl$padj <- res_tbl$adj_pvalue
    }
    if (!"pvalue" %in% colnames(res_tbl)) {
      pcol <- intersect(c("pvalue", "p_value"), colnames(res_tbl))[1]
      if (!is.na(pcol)) res_tbl$pvalue <- res_tbl[[pcol]]
    }
  }

  res_tbl
}

#' Run functional enrichment on isoform-level (DTE/DTU) results
#'
#' @description Thin adapter around `run_functional_analysis()` /
#'   `run_fgsea_analysis()` -- the exact same ORA/FGSEA engines and
#'   GO/Reactome/DO/MSigDB databases used for the DGE-level analysis -- fed
#'   with an isoform-level results table.
#'
#'   DTE (DESeq2, transcript-level) carries the same columns as a DGE-level
#'   `res_tbl` (`log2FoldChange`, `stat`, `padj`, `pvalue`, `entrezid`) and
#'   gets the full ORA + FGSEA treatment. `run_dte()` always runs its
#'   DESeq2 test with `test = "Wald"` regardless of whatever test the
#'   DGE-level pipeline used, so ranking here always uses the Wald branch.
#'
#'   DTU (DRIMSeq, transcript-usage) has no such statistic: DRIMSeq's
#'   feature-level `lr` is a directionless likelihood-ratio value, and no
#'   proportion-difference/coefficient column is computed upstream. FGSEA
#'   is skipped entirely for DTU, and `run_functional_analysis()` is called
#'   with `run_gsea = FALSE` so its Reactome/KEGG/GO GSEA (and SPIA) blocks
#'   are skipped too -- only GO/Reactome/DO ORA run.
#'
#' @param res_tbl DTE's `res_df`, or DTU's `dtu_results` data frame (i.e.
#'   `dtu_res$dtu_results` -- unwrap `run_dtu()`'s `list(dtu_results = ...)`
#'   return value before calling this).
#' @param analysis_type "DTE" or "DTU".
#' @param gmt_file MSigDB collections for FGSEA. Pass the *same* value used
#'   for the DGE-level `run_fgsea_analysis()` call so isoform- and gene-level
#'   FGSEA cover identical databases. Ignored when `analysis_type = "DTU"`.
#' @inheritParams run_functional_analysis
#' @return List with elements `functional` (see `run_functional_analysis()`)
#'   and `fgsea` (see `run_fgsea_analysis()`; always `NULL` for DTU).
#' @export
run_isoform_functional_analysis <- function(
  res_tbl,
  edb,
  out_dir,
  analysis_type = c("DTE", "DTU"),
  level,
  base,
  top_genes = 30,
  padj_cutoff = 0.01,
  go_pvalue_cutoff = 0.05,
  go_qvalue_cutoff = 0.2,
  gsea_metric = "stat",
  ora_padj_cutoff = 0.05,
  ora_min_genes = 10,
  ora_lfc_cutoff = NULL,
  gmt_file = c("C2", "C5", "C8"),
  run_spia = FALSE
) {
  analysis_type <- match.arg(analysis_type)
  comp_name <- paste0(level, "_vs_", base)
  run_gsea <- identical(analysis_type, "DTE")

  if (identical(analysis_type, "DTU") && !is.null(ora_lfc_cutoff)) {
    message(
      "   -> ora_lfc_cutoff ignored for DTU: DRIMSeq DTU results have no log2FoldChange column."
    )
    ora_lfc_cutoff <- NULL
  }

  res_tbl <- .adapt_isoform_res_tbl(res_tbl, analysis_type)
  out_dir_analysis <- safe_dir(file.path(
    out_dir,
    paste0(analysis_type, "_Enrichment")
  ))

  message(
    "\n=== Functional enrichment on ",
    analysis_type,
    " results (",
    comp_name,
    ") ==="
  )

  functional_res <- safe_run(
    run_functional_analysis(
      res_tbl = res_tbl,
      sig_res = NULL,
      edb = edb,
      out_dir = out_dir_analysis,
      level = level,
      base = base,
      top_genes = top_genes,
      padj_cutoff = padj_cutoff,
      go_pvalue_cutoff = go_pvalue_cutoff,
      go_qvalue_cutoff = go_qvalue_cutoff,
      gsea_metric = gsea_metric,
      test_type = "Wald", # run_dte() always tests with test = "Wald"; irrelevant for DTU (run_gsea = FALSE)
      ora_padj_cutoff = ora_padj_cutoff,
      ora_min_genes = ora_min_genes,
      ora_lfc_cutoff = ora_lfc_cutoff,
      run_gsea = run_gsea,
      run_spia = run_spia
    ),
    label = paste0(analysis_type, " ORA", if (run_gsea) "/GSEA" else "")
  )

  fgsea_res <- NULL
  if (run_gsea) {
    fgsea_res <- safe_run(
      run_fgsea_analysis(
        res_tbl = res_tbl,
        gmt_file = gmt_file,
        edb = edb,
        out_dir = out_dir_analysis,
        comp_name = comp_name,
        padj_cutoff = padj_cutoff,
        gsea_metric = gsea_metric,
        test_type = "Wald"
      ),
      label = paste0(analysis_type, " FGSEA")
    )
  } else {
    message(
      "   -> Skipping FGSEA for ",
      analysis_type,
      ": no signed ranking statistic is available from DRIMSeq DTU output."
    )
  }

  list(functional = functional_res, fgsea = fgsea_res)
}
