# ---- internal helpers --------------------------------------------------------

#' @keywords internal
.write_enrich_csv <- function(obj, path) {
  df <- as.data.frame(obj)
  if (nrow(df) > 0) write.csv(df, path, row.names = FALSE)
  invisible(df)
}

#' @keywords internal
.run_go_ontology <- function(ont, sigOE_genes, allOE_genes, org_db,
                             go_pvalue_cutoff, go_qvalue_cutoff,
                             OE_foldchanges, top_genes,
                             dir_go, comp_name) {
  message("   -> Running GO ", ont, "...")
  ego <- safe_run(
    suppressMessages(clusterProfiler::enrichGO(
      gene          = sigOE_genes,
      universe      = allOE_genes,
      keyType       = "SYMBOL",
      OrgDb         = org_db,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = go_pvalue_cutoff,
      qvalueCutoff  = go_qvalue_cutoff,
      readable      = FALSE
    )),
    label = paste("GO", ont)
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("      No significant results for GO ", ont)
    return(NULL)
  }

  .write_enrich_csv(ego, file.path(dir_go, paste0("GO_", ont, "_", comp_name, ".csv")))
  if (nrow(ego@result) > 5) {
    safe_pdf(file.path(dir_go, paste0("GO_", ont, "_Dotplot_", comp_name, ".pdf")),
             width = 10, height = 12,
             expr  = print(enrichplot::dotplot(ego, showCategory = top_genes)))
    safe_pdf(file.path(dir_go, paste0("GO_", ont, "_Cnetplot_", comp_name, ".pdf")),
             width = 12, height = 12,
             expr  = print(enrichplot::cnetplot(ego, showCategory = 10, foldChange = OE_foldchanges)))
  }
  return(ego)
}

#' @keywords internal
.hdo_schema_ok <- function(path) {
  tryCatch({
    con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
    "gene2allont" %in% DBI::dbListTables(con)
  }, error = function(e) FALSE)
}

#' @keywords internal
.ensure_hdo_sqlite <- function() {
  cache_dir <- tryCatch(yulab.utils::user_dir("GOSemSim"),
                        error = function(e) file.path(Sys.getenv("LOCALAPPDATA"), "GOSemSim"))
  hdo <- file.path(cache_dir, "HDO.sqlite")

  if (file.exists(hdo) && !.hdo_schema_ok(hdo)) { message("HDO.sqlite outdated; removing."); file.remove(hdo) }
  if (file.exists(hdo)) return(hdo)

  message("Downloading HDO.sqlite.gz...")
  safe_dir(cache_dir)
  gz_tmp <- tempfile(fileext = ".sqlite.gz")
  old_to <- getOption("timeout"); options(timeout = 300); on.exit(options(timeout = old_to), add = TRUE)

  urls <- c("https://yulab-smu.top/DOSE/HDO.sqlite.gz",
            "https://raw.githubusercontent.com/YuLab-SMU/DOSE/refs/heads/gh-pages/HDO.sqlite.gz")
  for (url in urls) {
    ok <- tryCatch({ download.file(url, gz_tmp, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (!ok) next
    tryCatch({
      con_gz <- gzcon(file(gz_tmp, "rb"))
      writeBin(readBin(con_gz, "raw", n = 200e6), hdo)
      close(con_gz)
    }, error = function(e) message("Decompression failed: ", e$message))
    if (.hdo_schema_ok(hdo)) return(hdo)
    if (file.exists(hdo)) file.remove(hdo)
  }
  NULL
}

#' @keywords internal
.run_disease_ontology <- function(sig_entrez, go_pvalue_cutoff, go_qvalue_cutoff, top_genes,
                                  dir_ora, comp_name) {
  hdo <- .ensure_hdo_sqlite()
  if (is.null(hdo)) { message("Skipping Disease Ontology (no HDO.sqlite)."); return(invisible(NULL)) }

  do_res <- safe_run(DOSE::enrichDO(gene = names(sig_entrez), 
                                    pvalueCutoff = go_pvalue_cutoff,
                                    qvalueCutoff = go_qvalue_cutoff),
                     label = "DOSE enrichDO")
  if (is.null(do_res) || nrow(as.data.frame(do_res)) == 0) return(invisible(NULL))

  dir_dose <- safe_dir(file.path(dir_ora, "DOSE"))
  .write_enrich_csv(do_res, file.path(dir_dose, paste0("DOSE_Enrichment_", comp_name, ".csv")))
  if (nrow(do_res@result) > 5) {
    safe_pdf(file.path(dir_dose, paste0("DOSE_Dotplot_", comp_name, ".pdf")),
             width = 10, height = 8,
             expr  = print(enrichplot::dotplot(do_res, showCategory = top_genes)))
  }
}

#' @keywords internal
.run_kegg_pathview <- function(gseaKEGG, foldchanges, kegg_code, top_genes, dir_kegg) {
  gseaKEGG_results <- as.data.frame(gseaKEGG)
  if (nrow(gseaKEGG_results) == 0) return(invisible(NULL))
  
  kegg_ids <- as.character(gseaKEGG_results$ID)
  if (length(kegg_ids) > top_genes) kegg_ids <- kegg_ids[seq_len(top_genes)]
  kegg_ids <- ifelse(grepl("^[a-zA-Z]", kegg_ids), kegg_ids, paste0(kegg_code, kegg_ids))

  # Force-load pathview internal datasets to avoid "object 'bods' not found"
  suppressMessages({
    utils::data("bods", package = "pathview", envir = environment())
    utils::data("gene.idtype.list", package = "pathview", envir = environment())
    if (!"package:pathview" %in% search()) try(attachNamespace("pathview"), silent = TRUE)
  })

  # Define the robust plotting function just like your snippet
  safe_get_kegg_plots <- purrr::safely(function(pid) {
    pathview::pathview(
      gene.data  = foldchanges, 
      pathway.id = pid, 
      species    = kegg_code,
      limit      = list(gene = 2, cpd = 1), # Enforce limit to fix the color scale issue
      kegg.dir   = "."
    )
    
    # Optional cleanup of leftover .xml/.png files normally dumped by pathview
    for (ext in c(".xml", ".png")) {
      f <- paste0(pid, ext)
      if (file.exists(f)) file.remove(f)
    }
  })

  # Set WD to the KEGG output directory so pathview downloads/saves there
  old_wd <- setwd(dir_kegg)
  on.exit(setwd(old_wd), add = TRUE)
  
  # Run the safe function over all KEGG IDs
  results <- purrr::map(kegg_ids, safe_get_kegg_plots)
  
  # Check for and report errors
  errors <- purrr::compact(purrr::map(results, "error"))
  if (length(errors) > 0) {
    message("   -> Warning: ", length(errors), " KEGG pathview plots failed to generate.")
  }
}



# ---- main --------------------------------------------------------------------

#' Run Functional Analysis
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
#' @return List with functional results
run_functional_analysis <- function(res_tbl, sig_res, edb, out_dir,
                                    level, base, top_genes, padj_cutoff,
                                    go_pvalue_cutoff = 0.05,
                                    go_qvalue_cutoff = 0.2) {
  comp_name <- paste0(level, "_vs_", base)
  org_info  <- get_organism_info(edb)
  message("Detected organism: ", org_info$name)
  org_db    <- org_info$org_db
  kegg_code <- org_info$kegg_code
  tf_db     <- org_info$tf_db

  # ---- DEBUG: inspect res_tbl ------------------------------------------------
  message("=== DEBUG res_tbl ===")
  message("  Dimensions : ", nrow(res_tbl), " rows x ", ncol(res_tbl), " cols")
  message("  Columns    : ", paste(colnames(res_tbl), collapse = ", "))
  for (col in c("gene", "log2FoldChange", "padj", "pvalue", "entrezid")) {
    if (col %in% colnames(res_tbl)) {
      vals <- res_tbl[[col]]
      if (is.numeric(vals)) {
        finite_vals <- vals[is.finite(vals)]
        rng <- if (length(finite_vals) > 0) paste(round(range(finite_vals), 4), collapse = " to ") else "no finite values"
        message("  [", col, "] numeric | NA=", sum(is.na(vals)),
                " | range=", rng,
                " | head=", paste(round(head(finite_vals, 5), 4), collapse=", "))
      } else {
        non_na <- vals[!is.na(vals) & vals != ""]
        message("  [", col, "] character | NA/empty=", sum(is.na(vals) | vals == ""),
                " | valid=", length(non_na),
                " | head=", paste(head(non_na, 5), collapse=", "))
      }
    } else {
      message("  [", col, "] MISSING")
    }
  }
  message("=== END DEBUG ===")
  # ----


  for (pkg in c(org_db, "ReactomePA", "DOSE", "pathview", "enrichplot", "enrichR", "clusterProfiler")) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE)
  }

  dir_ora  <- safe_dir(file.path(out_dir, "ORA"))
  dir_gsea <- safe_dir(file.path(out_dir, "GSEA"))

  # ---- gene sets
  allOE_genes <- as.character(res_tbl$gene[!is.na(res_tbl$gene)])
  sigOE       <- dplyr::filter(res_tbl, padj < padj_cutoff)
  if (nrow(sigOE) == 0 && "pvalue" %in% colnames(res_tbl))
    sigOE <- dplyr::filter(res_tbl, pvalue < 0.05)
  sigOE_genes <- unique(as.character(sigOE$gene[!is.na(sigOE$gene) & sigOE$gene != ""]))
  if (length(sigOE_genes) == 0) {
    message("No significant genes found. Skipping functional analysis."); return(NULL)
  }

  OE_foldchanges <- purrr::set_names(sigOE$log2FoldChange, sigOE$gene)
  OE_foldchanges <- pmin(pmax(OE_foldchanges, -2), 2)

  # ==============================================================================
  # PART 1: ORA (Fisher / Hypergeometric)
  # ==============================================================================

 # ---- 1A. ORA EnrichR TF
  message("Running EnrichR TF ORA...")
  if (length(tf_db) > 0) {
    
    suppressPackageStartupMessages(library(enrichR))
    websiteLive <- getOption("enrichR.live")
    
    if (!isTRUE(websiteLive)) {
      message("      Skipping EnrichR: enrichR.live is FALSE (no internet or API down)")
    } else {
      enrichr_results <- list()
      
      for (db in tf_db) {
        res <- safe_run({
          r <- enrichr(sigOE_genes, databases = db)
          if (length(r) > 0 && !is.null(r[[1]]) && nrow(r[[1]]) > 0) {
            df <- r[[1]]; df$Database <- db; df
          } else NULL
        }, label = paste("EnrichR TF", db))
        if (!is.null(res)) enrichr_results[[db]] <- res
      }
      
      if (length(enrichr_results) > 0) {
        df_all <- dplyr::bind_rows(enrichr_results)
        dir_tf <- safe_dir(file.path(out_dir, "Transcription_Factors"))
        write.csv(df_all, file.path(dir_tf, paste0("Enrichr_TF_", comp_name, ".csv")), row.names = FALSE)
        message("      EnrichR TF: ", length(enrichr_results), "/", length(tf_db), " DBs returned results.")
      } else {
        message("      EnrichR TF: No DBs returned valid results.")
      }
    }
  } else {
    message("      No TF databases specified for EnrichR ORA.")
  }

  # ---- 1B. ORA GO (BP, MF, CC)
  message("Running GO ORA...")
  dir_go <- safe_dir(file.path(dir_ora, "GO"))
  ego_list <- list()
  for (ont in c("BP", "MF", "CC")) {
    ego_list[[ont]] <- .run_go_ontology(ont, sigOE_genes, allOE_genes, org_db,
                                        go_pvalue_cutoff, go_qvalue_cutoff,
                                        OE_foldchanges, top_genes,
                                        dir_go, comp_name)
  }
  ego_list <- purrr::compact(ego_list)

# On-the-fly Entrez mapping fallback if column is missing or all-NA
  if (!"entrezid" %in% colnames(res_tbl) || all(is.na(res_tbl$entrezid)) || all(res_tbl$entrezid == "")) {
    message("   -> entrezid column empty — mapping gene symbols via AnnotationDbi...")
    org_obj <- tryCatch(getExportedValue(org_db, org_db), error = function(e) NULL)
    if (!is.null(org_obj)) {
      mapped <- suppressMessages(
        AnnotationDbi::mapIds(org_obj,
                              keys      = as.character(res_tbl$gene),
                              column    = "ENTREZID",
                              keytype   = "SYMBOL",
                              multiVals = "first")
      )
      res_tbl$entrezid <- mapped[as.character(res_tbl$gene)]
    } else {
      message("   -> WARNING: Could not load org_db object for on-the-fly mapping.")
    }
  }

  res_entrez <- res_tbl[!is.na(res_tbl$entrezid) & res_tbl$entrezid != "" & !is.na(res_tbl$log2FoldChange), ]
  res_entrez <- res_entrez[!duplicated(res_entrez$entrezid), ]
  
  # Diagnostic: report mapping rate
  total_genes  <- nrow(res_tbl)
  mapped_genes <- nrow(res_entrez)
  message("   -> Entrez mapping: ", mapped_genes, "/", total_genes, " genes have valid Entrez IDs")
  if (mapped_genes == 0) {
    message("   -> WARNING: No Entrez IDs found. Check that org_db (", org_db, ") is installed and gene symbols match.")
  }
  
  raw_lfc     <- res_entrez$log2FoldChange + runif(nrow(res_entrez), -1e-6, 1e-6)
  foldchanges <- sort(purrr::set_names(raw_lfc, as.character(res_entrez$entrezid)),
                      decreasing = TRUE)
  sig_entrez  <- foldchanges[names(foldchanges) %in%
                              as.character(res_entrez$entrezid[!is.na(res_entrez$padj) & res_entrez$padj < padj_cutoff])]
  message("   -> sig_entrez: ", length(sig_entrez), " DE genes with Entrez IDs (padj < ", padj_cutoff, ")")

  # ---- 1C. ORA Reactome
  message("Running Reactome ORA...")
  reac_org <- ifelse(kegg_code == "hsa", "human", "mouse")
  if (length(sig_entrez) > 0) {
    x <- safe_run(
      suppressMessages(ReactomePA::enrichPathway(
        gene         = names(sig_entrez), 
        organism     = reac_org, 
        pvalueCutoff = go_pvalue_cutoff, 
        qvalueCutoff = go_qvalue_cutoff, 
        readable     = TRUE)),
      label = "Reactome ORA"
    )
    if (!is.null(x) && nrow(as.data.frame(x)) > 0) {
      dir_reac_ora <- safe_dir(file.path(dir_ora, "Reactome"))
      .write_enrich_csv(x, file.path(dir_reac_ora, paste0("Reactome_ORA_", comp_name, ".csv")))
      if (nrow(x@result) > 5) {
        safe_pdf(file.path(dir_reac_ora, paste0("Reactome_ORA_Dotplot_", comp_name, ".pdf")),
                 width = 10, height = 10,
                 expr  = print(enrichplot::dotplot(x, showCategory = top_genes, label_format = 40, color = "pvalue")))
        safe_pdf(file.path(dir_reac_ora, paste0("Reactome_ORA_Cnetplot_", comp_name, ".pdf")),
                 width = 12, height = 12,
                 expr  = print(enrichplot::cnetplot(x, foldChange = foldchanges, showCategory = 10, colorEdge = TRUE)))
      }
    } else {
      message("      No significant results for Reactome ORA")
    }
  }

  # ---- 1D. ORA Disease Ontology (Human only)
  message("Running Disease Ontology ORA...")
  if (kegg_code == "hsa" && length(sig_entrez) > 0) {
    .run_disease_ontology(sig_entrez, go_pvalue_cutoff, go_qvalue_cutoff, top_genes, dir_ora, comp_name)
  } else {
    message("      Skipping Disease Ontology ORA (not human or no sig genes)")
  }

  # ==============================================================================
  # PART 2: GSEA (Ranked Lists)
  # ==============================================================================

  # ---- 2A. GSEA Reactome
  message("Running Reactome GSEA...")
  gseaReac <- safe_run(
    suppressMessages(ReactomePA::gsePathway(
      geneList     = foldchanges,
      organism     = reac_org,
      minGSSize    = 10,
      pvalueCutoff = padj_cutoff,
      verbose      = FALSE
    )),
    label = "Reactome GSEA"
  )
  if (!is.null(gseaReac) && nrow(as.data.frame(gseaReac)) > 0) {
    gseaReac <- suppressMessages(clusterProfiler::setReadable(gseaReac, OrgDb = org_db, keyType = "ENTREZID"))
    dir_reac_gsea <- safe_dir(file.path(dir_gsea, "Reactome"))
    .write_enrich_csv(gseaReac,
      file.path(dir_reac_gsea, paste0("Reactome_GSEA_", comp_name, ".csv")))
    
    if (nrow(gseaReac@result) > 5) {
      safe_pdf(file.path(dir_reac_gsea, paste0("Reactome_GSEA_Dotplot_", comp_name, ".pdf")),
               width = 10, height = 10,
               expr  = print(enrichplot::dotplot(gseaReac, showCategory = top_genes, label_format = 40, color = "pvalue")))
      safe_pdf(file.path(dir_reac_gsea, paste0("Reactome_GSEA_Ridgeplot_", comp_name, ".pdf")),
               width = 10, height = 10,
               expr  = print(enrichplot::ridgeplot(gseaReac, showCategory = 10, label_format = 40, color = "pvalue")))
      safe_pdf(file.path(dir_reac_gsea, paste0("Reactome_GSEA_EnrichmentMap_", comp_name, ".pdf")),
               width = 12, height = 8,
               expr  = print(enrichplot::gseaplot2(gseaReac, geneSetID = 1:min(3, nrow(gseaReac)))))
    }
  } else {
    message("      No significant results for Reactome GSEA")
  }

  # ---- 2B. GSEA KEGG + Pathview
  message("Running KEGG GSEA...")
  gseaKEGG <- safe_run(
    clusterProfiler::gseKEGG(geneList = foldchanges, organism = kegg_code,
                             minGSSize = 5, pvalueCutoff = padj_cutoff, verbose = FALSE),
    label = "KEGG GSEA"
  )
  if (!is.null(gseaKEGG) && nrow(as.data.frame(gseaKEGG)) > 0) {
    gseaKEGG <- suppressMessages(clusterProfiler::setReadable(gseaKEGG, OrgDb = org_db, keyType = "ENTREZID"))
    dir_kegg <- safe_dir(file.path(dir_gsea, "KEGG"))
    .write_enrich_csv(gseaKEGG, file.path(dir_kegg, paste0("KEGG_GSEA_", comp_name, ".csv")))
    
    if (nrow(gseaKEGG@result) > 5) {
      safe_pdf(file.path(dir_kegg, paste0("KEGG_GSEA_Dotplot_", comp_name, ".pdf")),
               width = 10, height = 10,
               expr  = print(enrichplot::dotplot(gseaKEGG, showCategory = top_genes)))
      safe_pdf(file.path(dir_kegg, paste0("KEGG_GSEA_Ridgeplot_", comp_name, ".pdf")),
               width = 10, height = 10,
               expr  = print(enrichplot::ridgeplot(gseaKEGG, showCategory = 10)))
    }
    
    message("Generating KEGG pathway maps...")
    .run_kegg_pathview(gseaKEGG, foldchanges, kegg_code, top_genes, dir_kegg)
  } else {
    message("      No significant KEGG pathways found in GSEA.")
  }

  # ==============================================================================
  # PART 3: TOPOLOGY SCORING (SPIA)
  # ==============================================================================

  # ---- 3A. SPIA
  message("Running SPIA analysis...")
  if (length(sig_entrez) == 0 && "pvalue" %in% colnames(res_entrez)) {
    sig_entrez <- foldchanges[names(foldchanges) %in%
                                as.character(res_entrez$entrezid[res_entrez$pvalue < 0.05])]
  }
  spia_result <- NULL
  if (length(sig_entrez) > 0) {
    message("   -> SPIA: ", length(sig_entrez), " DE genes, ", length(foldchanges), " background genes")
    spia_result <- safe_run(
      SPIA::spia(de = sig_entrez, all = as.character(res_entrez$entrezid),
                 organism = kegg_code, plots = FALSE),
      label = "SPIA"
    )
    
    if (!is.null(spia_result) && nrow(spia_result) > 0) {
      dir_spia <- safe_dir(file.path(out_dir, "SPIA"))
      
      # Save CSV first, independently of plotting
      write.csv(spia_result, file.path(dir_spia, paste0("SPIA_Results_", comp_name, ".csv")),
                row.names = FALSE)
      message("   -> SPIA: ", nrow(spia_result), " pathways saved to CSV.")
      
      # Plot separately so a plot failure does not block the CSV
      safe_pdf(file.path(dir_spia, paste0("SPIA_Evidence_", comp_name, ".pdf")),
               width = 8, height = 8,
               expr = plotP_fork(spia_result, threshold = padj_cutoff))
               
    } else if (!is.null(spia_result) && nrow(spia_result) == 0) {
      message("      SPIA returned an empty result (no pathways perturbed).")
    } else {
      message("      SPIA failed — check Log/Warnings.txt for details.")
    }
  } else {
    message("      Skipping SPIA: no significant Entrez genes available.")
  }

  message("Functional analysis complete.")
  list(ego_list = ego_list, spia_result = spia_result, gseaKEGG = gseaKEGG, gseaReac = gseaReac)
}