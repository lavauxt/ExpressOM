#' Run fgsea analysis using DE results and a GMT file
#' 
#' @description Performs Gene Set Enrichment Analysis using both `fgsea` and `clusterProfiler` using external GMT files (e.g., from MSigDB).
#' @param res_tbl The full results table from `export_significant_results()` containing 'gene' and 'log2FoldChange'.
#' @param gmt_file Path to a local `.gmt` file, or NULL to download Hallmark via msigdbr.
#' @param edb Ensembl Database (used to detect species for msigdbr).
#' @param out_dir Output directory for fgsea results.
#' @param comp_name Comparison name (e.g. "Treated_vs_Control") for file naming.
#' @param padj_cutoff Adjusted p-value cutoff for filtering significant pathways.
#' @export
run_fgsea_analysis <- function(res_tbl, gmt_file = NULL, edb, out_dir, comp_name, padj_cutoff = 0.01) {
  
  fgsea_out <- file.path(out_dir, "GSEA", "FGSEA")
  if(!dir.exists(fgsea_out)) dir.create(fgsea_out, recursive = TRUE)

  if (is.null(gmt_file) || !file.exists(gmt_file)) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Package 'msigdbr' is required to download pathways. Please install it.")
    }
    message("No valid local GMT provided. Downloading Hallmark pathways via msigdbr...")
    
    org_info <- get_organism_info(edb)
    msig_org <- org_info$msig_org
    msig_cat <- if (!is.null(org_info$msig_cat)) org_info$msig_cat else "H"
    
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
    
    gmt_name <- "hallmark_msigdbr"
    pathways.GSEA <- m_t2g
    pathways.fgsea <- split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)
    
  } else {
    gmt_name <- tools::file_path_sans_ext(basename(gmt_file))
    message('-> Loading pathways for fgsea/GSEA: ', gmt_name)
    pathways.GSEA <- clusterProfiler::read.gmt(gmt_file)
    pathways.fgsea <- fgsea::gmtPathways(gmt_file)
  }
   
  message('-> Preparing ranked gene list...')
  res2 <- res_tbl %>%
    dplyr::select(gene, log2FoldChange) %>%
    na.omit() %>%
    dplyr::distinct() %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = 'drop')
  
  ranks <- tibble::deframe(res2)
  
  set.seed(123456)
  ranks <- ranks + runif(length(ranks), min = -1e-6, max = 1e-6)
  ranks <- sort(ranks, decreasing = TRUE)
  
  message('-> Running clusterProfiler GSEA...')
  set.seed(123456)
  gsea_results <- tryCatch({
    suppressWarnings(suppressMessages(
      clusterProfiler::GSEA(ranks, TERM2GENE = pathways.GSEA, pvalueCutoff = padj_cutoff)
    ))
  }, error = function(e) NULL)
  
  message('-> Running fgseaMultilevel...')
  set.seed(123456)
  fgseaRes <- fgsea::fgseaMultilevel(pathways = pathways.fgsea, stats = ranks)
  
  message('-> Saving HTML table...')
  fgseaResTidy <- tibble::as_tibble(fgseaRes) %>% dplyr::arrange(dplyr::desc(NES))
  
  if (requireNamespace("DT", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    datatable_object <- fgseaResTidy %>%
      dplyr::select(-leadingEdge, -ES) %>%
      dplyr::arrange(padj) %>%
      DT::datatable()

    old_wd <- getwd()
    setwd(fgsea_out)
    html_name <- paste0("GSEA_Table_", gmt_name, ".html")
    htmlwidgets::saveWidget(datatable_object, file = html_name, selfcontained = TRUE)
    setwd(old_wd)
  } else {
    message("Skipping HTML table generation: 'DT' or 'htmlwidgets' is not installed.")
  }
  
  message('-> Plotting NES Barplot...')
  fgseaResTidy_filtered <- fgseaResTidy %>% dplyr::filter(padj < padj_cutoff)
  
  if (nrow(fgseaResTidy_filtered) > 0) {
    fgseaResTidy_filtered$pathway <- gsub("_", " ", fgseaResTidy_filtered$pathway)
    fgseaResTidy_filtered$pathway <- stringr::str_wrap(fgseaResTidy_filtered$pathway, width = 60)
    
    p_bar <- ggplot2::ggplot(fgseaResTidy_filtered, ggplot2::aes(reorder(pathway, NES), NES)) +
      ggplot2::geom_col(ggplot2::aes(fill = NES > 0), width = 0.6) +
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
      ggplot2::scale_fill_manual(values = c("TRUE" = "#1E90FF", "FALSE" = "#FF6347")) +
      ggplot2::guides(fill = "none")
    
    num_pathways <- nrow(fgseaResTidy_filtered)
    width_adjustment <- max(8, num_pathways * 0.25)  
    ggplot2::ggsave(file.path(fgsea_out, paste0("Barplot_", gmt_name, ".pdf")), plot = p_bar, width = width_adjustment, height = 8, device = "pdf")
  }
  
  message('-> Generating individual pathway gseaplots...')
  if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
    gene_set_ids <- gsea_results@result$ID
    # Remove any NA IDs before iterating
    valid_idx <- which(!is.na(gene_set_ids))
    
    for (i in valid_idx) {
      pid <- gene_set_ids[i]
      safe_run({
        pathway_name <- gsub("_", " ", gsea_results@result$Description[i])
        # Use integer index i, NOT the string pid — gseaplot2 requires integer geneSetID
        p <- enrichplot::gseaplot2(gsea_results, geneSetID = i, title = as.character(pathway_name))
        ggplot2::ggsave(filename = file.path(fgsea_out, paste0("GSEA_plot_", gmt_name, "_", pid, ".pdf")),
                        plot = p, device = "pdf", width = 8, height = 6)
        ggplot2::ggsave(filename = file.path(fgsea_out, paste0("GSEA_plot_", gmt_name, "_", pid, ".tiff")),
                        plot = p, device = "tiff", width = 8, height = 6, compression = "lzw", dpi = 600)
      }, label = paste("gseaplot", pid))
    }
  } else {
    message("   -> Skipping individual gseaplots: clusterProfiler GSEA returned no significant results.")
  }
  
  while(dev.cur() > 1) dev.off()
  message('-> FGSEA processing complete.')
}