# ==============================================================================
# mod_functional.R - Consolidated functional analysis (ORA, GSEA, SPIA)
# ==============================================================================

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
      OrgDb         = .load_org_db(org_db),
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
  
  ego_wrapped <- ego
  
  if (nrow(ego@result) > 5) {
    ego_wrapped@result$Description <- stringr::str_wrap(ego_wrapped@result$Description, width = 50)
    safe_pdf(file.path(dir_go, paste0("GO_", ont, "_Dotplot_", comp_name, ".pdf")),
             width = 14, height = 14,
             expr  = print(enrichplot::dotplot(ego_wrapped, showCategory = top_genes, label_format = 50)))
             
    safe_pdf(file.path(dir_go, paste0("GO_", ont, "_Cnetplot_", comp_name, ".pdf")),
             width = 14, height = 14,
             expr  = print(enrichplot::cnetplot(ego_wrapped, showCategory = top_genes, foldChange = OE_foldchanges)))
  }
  return(ego_wrapped)
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
  # Cross-platform cache directory
  cache_dir <- tryCatch({
    if (requireNamespace("rappdirs", quietly = TRUE)) {
      rappdirs::user_cache_dir("GOSemSim")
    } else {
      if (Sys.info()["sysname"] == "Windows") {
        file.path(Sys.getenv("LOCALAPPDATA"), "GOSemSim")
      } else {
        file.path("~", ".cache", "GOSemSim")
      }
    }
  }, error = function(e) file.path(tempdir(), "GOSemSim"))
  
  safe_dir(cache_dir)
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
      tryCatch(
        writeBin(readBin(con_gz, "raw", n = 200e6), hdo),
        error   = function(e) message("Decompression failed: ", e$message),
        finally = try(close(con_gz), silent = TRUE)
      )
    }, error = function(e) message("Could not open gz file: ", e$message))
    if (.hdo_schema_ok(hdo)) return(hdo)
    if (file.exists(hdo)) file.remove(hdo)
  }
  NULL
}

#' @keywords internal
.run_disease_ontology <- function(sig_entrez_ids, universe_entrez, go_pvalue_cutoff, go_qvalue_cutoff, top_genes,
                                  dir_ora, comp_name) {
  hdo <- .ensure_hdo_sqlite()
  if (is.null(hdo)) { message("Skipping Disease Ontology (no HDO.sqlite)."); return(invisible(NULL)) }

  do_res <- safe_run(DOSE::enrichDO(gene = sig_entrez_ids, 
                                    universe = universe_entrez,
                                    pvalueCutoff = go_pvalue_cutoff,
                                    qvalueCutoff = go_qvalue_cutoff),
                     label = "DOSE enrichDO")
  if (is.null(do_res) || nrow(as.data.frame(do_res)) == 0) return(invisible(NULL))

  dir_dose <- safe_dir(file.path(dir_ora, "DOSE"))
  .write_enrich_csv(do_res, file.path(dir_dose, paste0("DOSE_Enrichment_", comp_name, ".csv")))
  if (nrow(do_res@result) > 5) {
    safe_pdf(file.path(dir_dose, paste0("DOSE_Dotplot_", comp_name, ".pdf")),
             width = 14, height = 14,
             expr  = print(enrichplot::dotplot(do_res, showCategory = top_genes, label_format = 50)))
  }
}

#' @keywords internal
.run_kegg_pathview <- function(gseaKEGG, expression_vector, kegg_code, top_genes, dir_kegg) {
  gseaKEGG_results <- as.data.frame(gseaKEGG)
  if (nrow(gseaKEGG_results) == 0) return(invisible(NULL))
  
  kegg_ids <- as.character(gseaKEGG_results$ID)
  if (length(kegg_ids) > top_genes) kegg_ids <- kegg_ids[seq_len(top_genes)]
  kegg_ids <- ifelse(grepl("^[a-zA-Z]", kegg_ids), kegg_ids, paste0(kegg_code, kegg_ids))

  suppressMessages({
    utils::data("bods", package = "pathview", envir = environment())
    utils::data("gene.idtype.list", package = "pathview", envir = environment())
    if (!"package:pathview" %in% search()) try(attachNamespace("pathview"), silent = TRUE)
  })

  safe_get_kegg_plots <- purrr::safely(function(pid) {
    withr::with_dir(dir_kegg, {
      pathview::pathview(
        gene.data  = expression_vector, 
        pathway.id = pid, 
        species    = kegg_code,
        limit      = list(gene = 2, cpd = 1),
        kegg.dir   = "."
      )
      for (ext in c(".xml", ".png")) {
        f <- paste0(pid, ext)
        if (file.exists(f)) file.remove(f)
      }
    })
  })
  
  results <- purrr::map(kegg_ids, safe_get_kegg_plots)
  errors <- purrr::compact(purrr::map(results, "error"))
  if (length(errors) > 0) {
    message("   -> Warning: ", length(errors), " KEGG pathview plots failed to generate.")
  }
}

#' Run GSEA for Gene Ontology (BP, MF, CC)
#' @keywords internal
.run_go_gsea <- function(gene_list, org_db, ont, pvalue_cutoff, out_dir, comp_name) {
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
  .write_enrich_csv(gseago, file.path(dir_go_gsea, paste0("GO_GSEA_", ont, "_", comp_name, ".csv")))
  if (nrow(gseago@result) > 5) {
    safe_pdf(file.path(dir_go_gsea, paste0("GO_GSEA_", ont, "_Dotplot_", comp_name, ".pdf")),
             width = 14, height = 14,
             expr = print(enrichplot::dotplot(gseago, showCategory = 20, label_format = 50)))
  }
  return(gseago)
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
#' @return List with functional results
run_functional_analysis <- function(res_tbl, sig_res, edb, out_dir,
                                    level, base, top_genes, padj_cutoff = 0.01,
                                    go_pvalue_cutoff = 0.05,
                                    go_qvalue_cutoff = 0.2,
                                    gsea_metric = "stat",
                                    test_type = "Wald") { 
  
  comp_name <- paste0(level, "_vs_", base)
  org_info  <- get_organism_info(edb)
  message("Detected organism: ", org_info$name)
  org_db    <- org_info$org_db
  org_obj   <- .load_org_db(org_db)
  kegg_code <- org_info$kegg_code
  tf_db     <- org_info$tf_db

  gseaKEGG <- NULL
  gseaReac <- NULL
  spia_result <- NULL

  # Optional debug output (set option("ExpressOM.verbose" = TRUE) to enable)
  if (getOption("ExpressOM.verbose", FALSE)) {
    message("=== DEBUG res_tbl ===")
    message("  Dimensions : ", nrow(res_tbl), " rows x ", ncol(res_tbl), " cols")
    message("  Columns    : ", paste(colnames(res_tbl), collapse = ", "))
    for (col in c("gene", "log2FoldChange", "padj", "pvalue", "entrezid", "stat")) {
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
  }

  # Ensure required packages are available
  required_pkgs <- c("ReactomePA", "DOSE", "pathview", "enrichplot", "enrichR", "clusterProfiler", "msigdbr")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("The required package '", pkg, "' is missing. Please install it to proceed with functional analysis.")
    }
  }

  dir_ora  <- safe_dir(file.path(out_dir, "ORA"))
  dir_gsea <- safe_dir(file.path(out_dir, "GSEA"))

  allOE_genes <- as.character(res_tbl$gene[!is.na(res_tbl$gene)])
  sigOE       <- dplyr::filter(res_tbl, .data$padj < padj_cutoff)
  if (nrow(sigOE) == 0 && "pvalue" %in% colnames(res_tbl))
    sigOE <- dplyr::filter(res_tbl, .data$pvalue < 0.05)
  sigOE_genes <- unique(as.character(sigOE$gene[!is.na(sigOE$gene) & sigOE$gene != ""]))
  if (length(sigOE_genes) == 0) {
    message("No significant genes found. Skipping functional analysis."); return(NULL)
  }

  oe_fc_data <- sigOE[!is.na(sigOE$log2FoldChange), c("gene", "log2FoldChange"), drop = FALSE]
  if (nrow(oe_fc_data) > 0) {
    oe_fc <- stats::aggregate(
      log2FoldChange ~ gene, data = oe_fc_data,
      FUN = function(x) x[which.max(abs(x))][1]
    )
    OE_foldchanges <- purrr::set_names(oe_fc$log2FoldChange, oe_fc$gene)
    OE_foldchanges <- pmin(pmax(OE_foldchanges, -2), 2)
  } else {
    OE_foldchanges <- numeric(0)
  }

  # ==============================================================================
  # PART 1: ORA (Fisher / Hypergeometric)
  # ==============================================================================

  # ---- 1A. Transcription Factor (TF) enrichment ----
  #    Try EnrichR online first; if unreachable, fall back to local MSigDB C3 (TFT)
  message("Running Transcription Factor (TF) enrichment...")
  if (length(tf_db) > 0) {
    websiteLive <- getOption("enrichR.live", NA)
    if (is.na(websiteLive)) {
      websiteLive <- tryCatch({
        enrichR::listEnrichrDbs()
        TRUE
      }, error = function(e) FALSE)
      options(enrichR.live = websiteLive)
    }

    if (isTRUE(websiteLive)) {
      message("   Using EnrichR API...")
      enrichr_results <- list()
      for (db in tf_db) {
        res <- safe_run({
          r <- enrichR::enrichr(sigOE_genes, databases = db)
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
    } else {
      message("      EnrichR API unreachable; falling back to local TF enrichment using MSigDB C3 (TFT).")
      # Local enrichment using MSigDB C3 (TFT subset)
      msig_org <- org_info$msig_org
      message("      Using MSigDB C3 (TFT) for ", msig_org)
      
      m_t2g <- tryCatch({
        msigdbr::msigdbr(species = msig_org, collection = "C3", subcategory = "TFT") %>%
          dplyr::select(gs_name, gene_symbol)
      }, error = function(e) {
        message("      C3 TFT subcategory not available; using full C3 collection.")
        msigdbr::msigdbr(species = msig_org, collection = "C3") %>%
          dplyr::select(gs_name, gene_symbol)
      })
      
      if (nrow(m_t2g) > 0) {
        term2gene <- data.frame(term = m_t2g$gs_name, gene = m_t2g$gene_symbol, stringsAsFactors = FALSE)
        local_res <- tryCatch({
          clusterProfiler::enricher(
            gene = sigOE_genes,
            universe = allOE_genes,
            TERM2GENE = term2gene,
            pvalueCutoff = go_pvalue_cutoff,
            qvalueCutoff = go_qvalue_cutoff
          )
        }, error = function(e) {
          message("      Local TF enrichment failed: ", e$message)
          NULL
        })
        
        if (!is.null(local_res) && nrow(as.data.frame(local_res)) > 0) {
          dir_tf <- safe_dir(file.path(out_dir, "Transcription_Factors"))
          write.csv(as.data.frame(local_res), file.path(dir_tf, paste0("Local_TF_Enrichment_C3_", comp_name, ".csv")), row.names = FALSE)
          message("      Local TF enrichment completed with ", nrow(as.data.frame(local_res)), " terms.")
          if (nrow(local_res@result) > 5) {
            safe_pdf(file.path(dir_tf, paste0("Local_TF_Dotplot_", comp_name, ".pdf")),
                     width = 12, height = 10,
                     expr = print(enrichplot::dotplot(local_res, showCategory = 20, label_format = 50)))
          }
        } else {
          message("      No significant TF enrichment found locally.")
        }
      } else {
        message("      No gene sets found in MSigDB C3 for organism ", msig_org)
      }
    }
  } else {
    message("      No TF databases specified; skipping TF enrichment.")
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

  if (!"entrezid" %in% colnames(res_tbl) || all(is.na(res_tbl$entrezid)) || all(res_tbl$entrezid == "")) {
    message("   -> entrezid column empty — mapping gene symbols via AnnotationDbi...")
    if (!is.null(org_obj)) {
      mapped <- suppressMessages(
        AnnotationDbi::mapIds(org_obj,
                              keys       = as.character(res_tbl$gene),
                              column     = "ENTREZID",
                              keytype    = "SYMBOL",
                              multiVals  = "first")
      )
      res_tbl$entrezid <- mapped[as.character(res_tbl$gene)]
    } else {
      message("   -> WARNING: Could not load org_db object for on-the-fly mapping.")
    }
  }

  res_entrez <- res_tbl[!is.na(res_tbl$entrezid) & res_tbl$entrezid != "" & !is.na(res_tbl$log2FoldChange), ]
  res_entrez <- res_entrez[!duplicated(res_entrez$entrezid), ]
  
  total_genes  <- nrow(res_tbl)
  mapped_genes <- nrow(res_entrez)
  message("   -> Entrez mapping: ", mapped_genes, "/", total_genes, " genes have valid Entrez IDs")
  if (mapped_genes == 0) {
    message("   -> WARNING: No Entrez IDs found. Check that org_db (", org_db, ") is installed and gene symbols match.")
  }

  message("   -> Generating ranked list using metric: ", gsea_metric)
  if (gsea_metric == "stat") {
    if (!"stat" %in% colnames(res_entrez)) {
      stop("Column 'stat' not found in results table. Verify that you injected res_unshrunken$stat back into your results.")
    }
    
    if (toupper(test_type) == "LRT") {
      message("   -> [MANAGEMENT] LRT design detected. Transforming Chi-Square metrics into directional stats via sign(log2FoldChange).")
      metric_vals <- sign(res_entrez$log2FoldChange) * res_entrez$stat
    } else {
      message("   -> Wald design detected. Using native directional Wald z-scores.")
      metric_vals <- res_entrez$stat
    }
    
  } else if (gsea_metric == "signed_pval") {
    safe_pvals  <- ifelse(res_entrez$pvalue == 0, .Machine$double.xmin, res_entrez$pvalue)
    metric_vals <- sign(res_entrez$log2FoldChange) * -log10(safe_pvals)
  } else {
    metric_vals <- res_entrez$log2FoldChange
  }
  set.seed(123456)
  metric_vals <- metric_vals + runif(nrow(res_entrez), -1e-9, 1e-9)
  gsea_list <- sort(purrr::set_names(metric_vals, as.character(res_entrez$entrezid)), 
                    decreasing = TRUE)

  sig_entrez_ids <- as.character(res_entrez$entrezid[!is.na(res_entrez$padj) & res_entrez$padj < padj_cutoff])
  if (length(sig_entrez_ids) == 0 && "pvalue" %in% colnames(res_entrez)) {
    sig_entrez_ids <- as.character(res_entrez$entrezid[!is.na(res_entrez$pvalue) & res_entrez$pvalue < 0.05])
  }
  message("   -> Significant ORA genes with Entrez IDs: ", length(sig_entrez_ids))

  # ==============================================================================
  # PART 1C.1: Curated List specifically for Reactome ORA
  # ==============================================================================
  message(" -> Generating Curated List specifically for Reactome ORA...")
  reac_lfc_cutoff <- 1.0 
  
  res_reac_curated <- res_entrez[as.character(res_entrez$entrezid) %in% sig_entrez_ids & 
                                   abs(res_entrez$log2FoldChange) >= reac_lfc_cutoff, ]
  message(sprintf("   -> Reactome Filter (%s): Candidate background & abs(LFC) >= %s", test_type, reac_lfc_cutoff))
  
  sig_entrez_reac <- as.character(res_reac_curated$entrezid)
  message("   -> Reactome curated ORA genes available: ", length(sig_entrez_reac))

  # ---- 1D. ORA Reactome
  message("Running Reactome ORA on curated list...")
  reac_org <- ifelse(kegg_code == "hsa", "human", "mouse")
  if (length(sig_entrez_reac) > 0) {
    x <- safe_run(
      suppressMessages(ReactomePA::enrichPathway(
        gene         = sig_entrez_reac,          
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
                 width = 14, height = 14,
                 expr  = print(enrichplot::dotplot(x, showCategory = top_genes, label_format = 50, color = "pvalue")))
        
        safe_pdf(file.path(dir_reac_ora, paste0("Reactome_ORA_Cnetplot_", comp_name, ".pdf")),
                 width = 12, height = 12,
                 expr  = print(enrichplot::cnetplot(x, foldChange = OE_foldchanges, showCategory = top_genes)))
      }
    } else {
      message("      No significant results for Reactome ORA")
    }
  }

  # ---- 1E. ORA Disease Ontology (Human only)
  message("Running Disease Ontology ORA...")
  if (kegg_code == "hsa" && length(sig_entrez_ids) > 0) {
    .run_disease_ontology(sig_entrez_ids, universe_entrez = as.character(res_entrez$entrezid), 
                          go_pvalue_cutoff, go_qvalue_cutoff, top_genes, dir_ora, comp_name)
  } else {
    message("      Skipping Disease Ontology ORA (not human or no sig genes)")
  }

  # ==============================================================================
  # PART 2: GSEA (Properly tracking across genome background vector)
  # ==============================================================================

  # ---- 2A. GSEA Reactome
  message("Running Reactome GSEA on complete ranked genome background...")
  set.seed(123456)
  gseaReac <- safe_run(
    suppressMessages(ReactomePA::gsePathway(
      geneList      = gsea_list,          
      organism      = reac_org,
      minGSSize     = 10,                                                                            
      maxGSSize     = 300,          
      pvalueCutoff  = padj_cutoff,
      pAdjustMethod = "BH", 
      verbose       = FALSE
    )),
    label = "Reactome GSEA"
  )

  if (!is.null(gseaReac) && nrow(as.data.frame(gseaReac)) > 0) {
    dir_reac_gsea <- safe_dir(file.path(dir_gsea, "Reactome"))
    
    if (nrow(gseaReac@result) > 5) {
      safe_pdf(file.path(dir_reac_gsea, paste0("Reactome_GSEA_Ridgeplot_", comp_name, ".pdf")),
               width = 14, height = 14,
               expr = {
                 p <- enrichplot::ridgeplot(gseaReac, showCategory = top_genes, label_format = 50)
                 print(p)
               })
    }
    
    gseaReac <- suppressMessages(clusterProfiler::setReadable(gseaReac, OrgDb = org_obj, keyType = "ENTREZID"))
    .write_enrich_csv(gseaReac, file.path(dir_reac_gsea, paste0("Reactome_GSEA_", comp_name, ".csv")))
    
    if (nrow(gseaReac@result) > 0) {
      safe_pdf(file.path(dir_reac_gsea, paste0("Reactome_GSEA_Dotplot_", comp_name, ".pdf")),
               width = 14, height = 14,
               expr = {
                 p <- enrichplot::dotplot(gseaReac, showCategory = top_genes, label_format = 50, color = "p.adjust")
                 print(p)
               })
  
      safe_pdf(file.path(dir_reac_gsea, paste0("Reactome_GSEA_Gseaplot_", comp_name, ".pdf")),
               width = 12, height = 8,
               expr  = print(enrichplot::gseaplot2(gseaReac, geneSetID = 1:min(3, nrow(gseaReac@result)))))
    }
  } else {
    message("      No significant results for Reactome GSEA")
  }

  # ---- 2B. GSEA KEGG + Pathview
  message("Running KEGG GSEA...")
  set.seed(123456)
  gseaKEGG <- safe_run(
    clusterProfiler::gseKEGG(geneList = gsea_list,
                             organism = kegg_code,
                             minGSSize = 5, pvalueCutoff = padj_cutoff, verbose = FALSE),
    label = "KEGG GSEA"
  )

  if (!is.null(gseaKEGG) && nrow(as.data.frame(gseaKEGG)) > 0) {
    disease_pattern <- paste0("^", kegg_code, "05")
    gseaKEGG@result <- gseaKEGG@result[!grepl(disease_pattern, gseaKEGG@result$ID), ]
    
    if (nrow(gseaKEGG@result) == 0) {
      message("      No non-disease KEGG pathways remained after filtering.")
      gseaKEGG <- NULL
    } else {
      message("      Filtered out KEGG Disease pathways. Remaining: ", nrow(gseaKEGG@result))
    }
  }

  if (!is.null(gseaKEGG) && nrow(as.data.frame(gseaKEGG)) > 0) {
    dir_kegg <- safe_dir(file.path(dir_gsea, "KEGG"))

    if (nrow(gseaKEGG@result) > 5) {
      safe_pdf(file.path(dir_kegg, paste0("KEGG_GSEA_Ridgeplot_", comp_name, ".pdf")),
               width = 14, height = 14,
               expr  = print(enrichplot::ridgeplot(gseaKEGG, showCategory = top_genes, label_format = 50)))
    }
    
    gseaKEGG <- suppressMessages(clusterProfiler::setReadable(gseaKEGG, OrgDb = org_obj, keyType = "ENTREZID"))
    .write_enrich_csv(gseaKEGG, file.path(dir_kegg, paste0("KEGG_GSEA_", comp_name, ".csv")))
    
    if (nrow(gseaKEGG@result) > 5) {
      safe_pdf(file.path(dir_kegg, paste0("KEGG_GSEA_Dotplot_", comp_name, ".pdf")),
               width = 14, height = 14,
               expr  = print(enrichplot::dotplot(gseaKEGG, showCategory = top_genes, label_format = 50)))
    }
    
    message("Generating KEGG pathway maps...")
    all_lfc_vector <- purrr::set_names(res_entrez$log2FoldChange, as.character(res_entrez$entrezid))
    .run_kegg_pathview(gseaKEGG, all_lfc_vector, kegg_code, top_genes, dir_kegg)
  } else {
    message("      No significant KEGG pathways found in GSEA.")
  }

  # ---- 2C. GSEA GO (BP, MF, CC)
  message("Running GO GSEA...")
  for (ont in c("BP", "MF", "CC")) {
    .run_go_gsea(gsea_list, org_db, ont, padj_cutoff, out_dir, comp_name)
  }

  # ==============================================================================
  # PART 3: TOPOLOGY SCORING (SPIA)
  # ==============================================================================

  message("Running SPIA analysis...")
  if (length(sig_entrez_ids) > 0) {
    # SPIA requires a named numeric vector of L2FC for DE genes; deduplicate by entrezid
    spia_de <- purrr::set_names(as.numeric(res_entrez$log2FoldChange), as.character(res_entrez$entrezid))
    spia_de <- spia_de[names(spia_de) %in% sig_entrez_ids]
    spia_de <- spia_de[!is.na(spia_de) & !duplicated(names(spia_de))]
    
    message("   -> SPIA: ", length(spia_de), " DE genes, ", length(gsea_list), " background genes")
    spia_result <- safe_run(
      SPIA::spia(de = spia_de, all = as.character(res_entrez$entrezid),
                 organism = kegg_code, plots = FALSE),
      label = "SPIA"
    )
    
    if (!is.null(spia_result) && nrow(spia_result) > 0) {
      dir_spia <- safe_dir(file.path(out_dir, "SPIA"))
      write.csv(spia_result, file.path(dir_spia, paste0("SPIA_Results_", comp_name, ".csv")), row.names = FALSE)
      message("   -> SPIA: ", nrow(spia_result), " pathways saved to CSV.")
      
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

# ------------------------------------------------------------------------------
# FGSEA Analysis (from fgsea.R)
# ------------------------------------------------------------------------------

#' Run fgsea analysis using DE results and a GMT file (or multiple GMT files)
#' 
#' @description Performs Gene Set Enrichment Analysis using both `fgsea` and `clusterProfiler` using external GMT files (e.g., from MSigDB) or downloading via msigdbr.
#' @param res_tbl The full results table from `export_significant_results()` containing 'gene' and 'log2FoldChange'.
#' @param gmt_file Path to a local `.gmt` file, a character vector/list of multiple `.gmt` files, or MSigDB category names (e.g., "H", "C2") to download via msigdbr.
#' @param edb Ensembl Database (used to detect species for msigdbr).
#' @param out_dir Output directory for fgsea results.
#' @param comp_name Comparison name (e.g. "Treated_vs_Control") for file naming.
#' @param padj_cutoff Adjusted p-value cutoff for filtering significant pathways.
#' @export
run_fgsea_analysis <- function(res_tbl, gmt_file = c("C2", "C5", "C8"), edb, out_dir, comp_name, padj_cutoff = 0.01) {
  
  gmt_list <- if (is.null(gmt_file)) list(NULL) else as.list(gmt_file)

  for (gmt_item in gmt_list) {
    
    if (is.null(gmt_item) || !file.exists(gmt_item)) {
      if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("Package 'msigdbr' is required to download pathways. Please install it.")
      }
      
      org_info <- get_organism_info(edb)
      msig_org <- org_info$msig_org
      
      valid_collections <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
                             "MH", "M1", "M2", "M3", "M5", "M7", "M8", "HALLMARK")
      
      input_cat <- toupper(gmt_item)
      
      if (!is.null(gmt_item) && input_cat %in% valid_collections) {    
              if (input_cat %in% c("MH", "HALLMARK")) {
                msig_cat <- "H"
              } else if (grepl("^M[0-9]", input_cat)) {
                msig_cat <- paste0("C", gsub("M", "", input_cat))
              } else {
                msig_cat <- input_cat
              }
        gmt_name <- paste0("msigdbr_", input_cat)
        message("-> Fetching MSigDB category [", input_cat, "] (mapped to ", msig_cat, ") via msigdbr for ", msig_org, "...")
        } else {
          message("Provided GMT file '", gmt_item, "' not recognized as MSigDB category. Falling back to Hallmark...")
          msig_cat <- "H"
          gmt_name <- "hallmark_msigdbr"
      }
      
      .msigdbr_fetch <- function(species, cat) {
        fmls <- names(formals(msigdbr::msigdbr))
        if ("collection" %in% fmls) {
          msigdbr::msigdbr(species = species, collection = cat) %>%
            dplyr::select(gs_name, gene_symbol)
        } else {
          msigdbr::msigdbr(species = species, category = cat) %>%
            dplyr::select(gs_name, gene_symbol)
        }
      }

      m_t2g <- tryCatch(
        .msigdbr_fetch(msig_org, msig_cat),
        error = function(e) {
          message("   -> msigdbr collection '", msig_cat, "' not available, falling back to 'H'")
          tryCatch(
            .msigdbr_fetch(msig_org, "H"),
            error = function(e2) {
              stop("msigdbr failed for both '", msig_cat, "' and 'H': ", e2$message)
            }
          )
        }
      )
      
      pathways.GSEA <- m_t2g
      pathways.fgsea <- split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)
      
    } else {
      gmt_name <- tools::file_path_sans_ext(basename(gmt_item))
      message('-> Loading pathways for fgsea/GSEA: ', gmt_name)
      pathways.GSEA <- clusterProfiler::read.gmt(gmt_item)
      pathways.fgsea <- fgsea::gmtPathways(gmt_item)
    }
     
    fgsea_out <- file.path(out_dir, "GSEA", "FGSEA", gmt_name)
    if(!dir.exists(fgsea_out)) dir.create(fgsea_out, recursive = TRUE)

    message('-> Preparing ranked gene list for [', gmt_name, ']...')
    res2 <- res_tbl |>
      dplyr::select(gene, log2FoldChange) |>
      stats::na.omit() |>
      dplyr::distinct() |>
      dplyr::group_by(gene) |>
      dplyr::summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = 'drop')
    
    ranks <- tibble::deframe(res2)
    
    set.seed(123456)
    ranks <- ranks + runif(length(ranks), min = -1e-6, max = 1e-6)
    ranks <- sort(ranks, decreasing = TRUE)
    
    message('-> Running clusterProfiler GSEA for [', gmt_name, ']...')
    set.seed(123456)
    gsea_results <- tryCatch({
      suppressWarnings(suppressMessages(
        clusterProfiler::GSEA(ranks, TERM2GENE = pathways.GSEA, pvalueCutoff = padj_cutoff)
      ))
    }, error = function(e) NULL)
    
    message('-> Running fgseaMultilevel for [', gmt_name, ']...')
    set.seed(123456)
    fgseaRes <- fgsea::fgseaMultilevel(pathways = pathways.fgsea, stats = ranks)
    
    message('-> Saving results table for [', gmt_name, ']...')
    fgseaResTidy <- tibble::as_tibble(fgseaRes) |> dplyr::arrange(dplyr::desc(NES))
    utils::write.csv(
      fgseaResTidy %>% dplyr::select(-leadingEdge),
      file.path(fgsea_out, paste0("FGSEA_Results_", gmt_name, ".csv")),
      row.names = FALSE
    )
    
    if (requireNamespace("DT", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
      datatable_object <- fgseaResTidy |>
        dplyr::select(-leadingEdge, -ES) |>
        dplyr::arrange(padj) |>
        DT::datatable()
      html_name <- paste0("GSEA_Table_", gmt_name, ".html")
      withr::with_dir(fgsea_out, {
        htmlwidgets::saveWidget(datatable_object, file = html_name, selfcontained = TRUE)
      })
      
    } else {
      message("Skipping HTML table generation: 'DT' or 'htmlwidgets' is not installed.")
    }
    
    message('-> Plotting NES Barplot for [', gmt_name, ']...')
    fgseaResTidy_filtered <- fgseaResTidy |> dplyr::filter(padj < padj_cutoff)
    
    if (nrow(fgseaResTidy_filtered) > 0) {
      fgseaResTidy_filtered$pathway <- gsub("_", " ", fgseaResTidy_filtered$pathway)
      fgseaResTidy_filtered$pathway <- stringr::str_wrap(fgseaResTidy_filtered$pathway, width = 60)
      
      fgseaResTidy_filtered$direction <- ifelse(fgseaResTidy_filtered$NES > 0, "Up", "Down")

      p_bar <- ggplot2::ggplot(fgseaResTidy_filtered, ggplot2::aes(reorder(pathway, NES), NES)) +
        ggplot2::geom_col(ggplot2::aes(fill = direction), width = 0.6) +
        ggplot2::coord_flip() +  
        ggplot2::labs(x = "Pathway", y = "Normalized Enrichment Score", title = paste("GSEA NES:", gmt_name)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          axis.title.x = ggplot2::element_text(face = "bold"),
          axis.title.y = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(face = "bold"),
          axis.text.y = ggplot2::element_text(face = "bold", size = 8),
          plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")
        ) +
        ggplot2::scale_fill_manual(values = c("Up" = "#1E90FF", "Down" = "#FF6347"), drop = FALSE) +
        ggplot2::guides(fill = "none")
      
      num_pathways <- nrow(fgseaResTidy_filtered)
      height_adjustment <- max(8, num_pathways * 0.3)
      ggplot2::ggsave(file.path(fgsea_out, paste0("Barplot_", gmt_name, ".pdf")), plot = p_bar, width = 12, height = height_adjustment, device = "pdf")
    }
    
    message('-> Generating individual pathway gseaplots for [', gmt_name, ']...')
    max_gsea_plots <- getOption("ExpressOM.max_gsea_plots", 20L)
    if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
      gene_set_ids <- gsea_results@result$ID
      valid_idx <- which(!is.na(gene_set_ids))
      if (length(valid_idx) > max_gsea_plots) {
        message("   -> Limiting individual gseaplots to top ", max_gsea_plots,
                " of ", length(valid_idx),
                " (set options(ExpressOM.max_gsea_plots = N) to change)")
        valid_idx <- valid_idx[seq_len(max_gsea_plots)]
      }
      for (i in valid_idx) {
        pid <- gene_set_ids[i]
        clean_pid <- gsub("[^A-Za-z0-9_-]", "_", pid)
        
        safe_run({
          pathway_name <- gsub("_", " ", gsea_results@result$Description[i])
          p <- enrichplot::gseaplot2(gsea_results, geneSetID = i, title = as.character(pathway_name))
          ggplot2::ggsave(filename = file.path(fgsea_out, paste0("GSEA_plot_", gmt_name, "_", clean_pid, ".pdf")),
                          plot = p, device = "pdf", width = 8, height = 6)
          ggplot2::ggsave(filename = file.path(fgsea_out, paste0("GSEA_plot_", gmt_name, "_", clean_pid, ".tiff")),
                          plot = p, device = "tiff", width = 8, height = 6, compression = "lzw", dpi = 600)
        }, label = paste("gseaplot", pid))
      }
    } else {
      message("   -> Skipping individual gseaplots: clusterProfiler GSEA returned no significant results.")
    }
    
    while(grDevices::dev.cur() > 1) grDevices::dev.off()
  }
  message('-> FGSEA processing complete across all specified gene sets.')
}

# ------------------------------------------------------------------------------
# Local enrichment analysis (replaces EnrichR for offline use)
# ------------------------------------------------------------------------------

#' Run local enrichment analysis (replaces EnrichR)
#' @param gene_list Vector of Significant Gene Symbols (e.g., DEGs)
#' @param universe Vector of ALL genes expressed in the experiment (Background)
#' @param organism "Homo sapiens" or "Mus musculus"
#' @param category "H" (Hallmark), "C3" (TFT/ChIP-seq targets), or "C5" (GO)
#' @return A clusterProfiler result object
run_local_enrichment <- function(gene_list, universe, organism = "Homo sapiens", category = "C3") {
  for (pkg in c("clusterProfiler", "msigdbr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')")
  }
  org_pkg <- ifelse(grepl("Homo", organism), "org.Hs.eg.db", "org.Mm.eg.db")
  if (!requireNamespace(org_pkg, quietly = TRUE)) stop("Please install ", org_pkg)
  
  gene_entrez <- clusterProfiler::bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_pkg, drop=FALSE)
  univ_entrez <- clusterProfiler::bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_pkg, drop=FALSE)
  
  msigdbr_df <- msigdbr::msigdbr(species = organism, category = category)
  term2gene <- msigdbr_df[, c("gs_name", "entrez_gene")]
  
  res <- clusterProfiler::enricher(
    gene = na.omit(gene_entrez$ENTREZID),
    universe = na.omit(univ_entrez$ENTREZID),
    TERM2GENE = term2gene,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  return(res)
}