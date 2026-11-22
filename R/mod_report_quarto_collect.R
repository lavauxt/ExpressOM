# mod_report_quarto_collect.R -- everything the Quarto dashboard displays is
# gathered here, once, as (a) tables for the binary columnar container and
# (b) small JSON "view specs". Heavy statistics (sample correlation, MDS,
# dendrogram, boxplot summaries, densities) are computed in R, so the browser
# only draws.

.xq_arr <- function(x) I(as.vector(x))

# ---- helpers ---------------------------------------------------------------------

.xq_rbind_fill <- function(dfs) {
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (!length(dfs)) return(NULL)
  if (length(dfs) == 1L) return(as.data.frame(dfs[[1L]], stringsAsFactors = FALSE))
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  out <- lapply(cols, function(cn) {
    parts <- lapply(dfs, function(d) {
      if (cn %in% names(d)) d[[cn]] else rep(NA, nrow(d))
    })
    all_num <- all(vapply(parts, function(p) is.numeric(p) || all(is.na(p)), logical(1)))
    if (all_num) {
      unlist(lapply(parts, as.numeric), use.names = FALSE)
    } else {
      unlist(lapply(parts, as.character), use.names = FALSE)
    }
  })
  names(out) <- cols
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Canonical enrichment columns: `collection`, `term`, `padj`, `size`, `score`
#' (the source columns they come from are dropped; everything else is kept).
#' @keywords internal
.xq_enrich_std <- function(tb, collection) {
  tb <- as.data.frame(tb, stringsAsFactors = FALSE, check.names = FALSE)
  term_c <- .rep_first_col(tb, c("Description", "pathway", "Name", "Term", "ID", "source"))
  padj_c <- .rep_first_col(tb, c("p.adjust", "padj", "pGFdr", "Adjusted.P.value",
                                  "qvalue", "qvalues", "FDR", "p_value"))
  if (is.na(term_c) || is.na(padj_c)) {
    return(cbind(data.frame(collection = rep(collection, nrow(tb)), stringsAsFactors = FALSE), tb))
  }
  size_c <- .rep_first_col(tb, c("Count", "size", "setSize", "NDE", "Overlap"))
  score_c <- .rep_first_col(tb, c("NES", "enrichmentScore", "tA", "Combined.Score", "score"))
  size <- if (is.na(size_c)) rep(NA_real_, nrow(tb))
          else if (size_c == "Overlap") .rep_num(sub("/.*$", "", as.character(tb[[size_c]])))
          else .rep_num(tb[[size_c]])
  head_df <- data.frame(
    collection = rep(collection, nrow(tb)),
    term = as.character(tb[[term_c]]),
    padj = .rep_num(tb[[padj_c]]),
    size = size,
    score = if (is.na(score_c)) rep(NA_real_, nrow(tb)) else .rep_num(tb[[score_c]]),
    stringsAsFactors = FALSE
  )
  used <- stats::na.omit(c(term_c, padj_c, size_c, score_c))
  rest <- tb[setdiff(names(tb), used)]
  clash <- names(rest) %in% names(head_df)
  if (any(clash)) names(rest)[clash] <- paste0(names(rest)[clash], ".orig")
  cbind(head_df, rest)
}

.xq_n_ranked <- function(res) {
  de <- res$dge$results
  if (is.null(de)) return(0L)
  is_lrt <- identical(res$config$test, "LRT") && "contrast_stat" %in% names(de)
  sc <- if (is_lrt) "contrast_stat" else "stat"
  if (!sc %in% names(de) || !"gene" %in% names(de)) return(0L)
  keep <- !is.na(de[[sc]]) & !is.na(de$gene) & de$gene != ""
  length(unique(de$gene[keep]))
}

.xq_group_list <- function(groups, samples) {
  if (is.null(groups) || length(groups) != length(samples)) return(NULL)
  g <- as.character(groups)
  names(g) <- samples
  g <- g[!is.na(g)]
  if (!length(g)) NULL else as.list(g)
}

# ---- sample-level QC computed in R ------------------------------------------------------

#' @keywords internal
.xq_box_stats <- function(m) {
  lapply(seq_len(ncol(m)), function(j) {
    v <- m[, j]
    v <- v[is.finite(v)]
    if (length(v) < 5L) return(NULL)
    b <- grDevices::boxplot.stats(v)$stats
    list(label = colnames(m)[j], lo = round(b[1L], 4), q1 = round(b[2L], 4),
         med = round(b[3L], 4), q3 = round(b[4L], 4), hi = round(b[5L], 4))
  })
}

#' Correlation / distance / MDS / dendrogram / box / RLE / density summaries
#' @keywords internal
.xq_qc_payload <- function(ex, groups = NULL) {
  out <- list(similarity = NULL, expr = NULL, mv = NULL)
  sim_mat <- .rep_matrix_for_qc(ex$qc_matrix %||% ex$matrix)
  dist_mat <- .rep_matrix_for_qc(ex$distribution_matrix %||% ex$qc_matrix %||% ex$matrix)

  if (!is.null(sim_mat)) {
    x <- log2(sim_mat + 1)
    sim <- list(groups = groups)
    cm <- suppressWarnings(stats::cor(x, method = "pearson", use = "pairwise.complete.obs"))
    cm[!is.finite(cm)] <- 0
    diag(cm) <- 1
    if (nrow(cm) >= 3L) {
      ord <- tryCatch(stats::hclust(stats::as.dist(1 - cm))$order, error = function(e) seq_len(nrow(cm)))
      cm <- cm[ord, ord, drop = FALSE]
    }
    sim$cor <- list(labels = .xq_arr(rownames(cm)), values = .xq_arr(round(t(cm), 4)))

    d <- as.matrix(stats::dist(t(x), method = "euclidean"))
    if (nrow(d) >= 3L) {
      ord <- tryCatch(stats::hclust(stats::as.dist(d))$order, error = function(e) seq_len(nrow(d)))
      d <- d[ord, ord, drop = FALSE]
    }
    sim$dist <- list(labels = .xq_arr(rownames(d)), values = .xq_arr(round(t(d), 3)))

    if (ncol(sim_mat) >= 3L) {
      dd <- stats::dist(t(x))
      co <- tryCatch(stats::cmdscale(dd, k = 2L), error = function(e) NULL)
      if (!is.null(co) && ncol(co) >= 2L) {
        sim$mds <- lapply(seq_len(nrow(co)), function(i) {
          list(label = rownames(co)[i], x = round(co[i, 1L], 4), y = round(co[i, 2L], 4))
        })
      }
      hc <- tryCatch(stats::hclust(dd, method = "complete"), error = function(e) NULL)
      if (!is.null(hc)) {
        n <- length(hc$order)
        xpos <- rep(NA_real_, 2L * n - 1L)
        ypos <- rep(0, 2L * n - 1L)
        for (i in seq_len(n)) xpos[i] <- match(i, hc$order)
        segs <- vector("list", 3L * (n - 1L))
        for (k in seq_len(n - 1L)) {
          child_pos <- function(v) if (v < 0L) -v else n + v
          a <- child_pos(hc$merge[k, 1L]); b <- child_pos(hc$merge[k, 2L])
          parent <- n + k
          xpos[parent] <- mean(c(xpos[a], xpos[b]))
          ypos[parent] <- hc$height[k]
          h <- hc$height[k]
          segs[[3L * k - 2L]] <- list(x1 = xpos[a], y1 = ypos[a], x2 = xpos[a], y2 = h)
          segs[[3L * k - 1L]] <- list(x1 = xpos[b], y1 = ypos[b], x2 = xpos[b], y2 = h)
          segs[[3L * k]]      <- list(x1 = xpos[a], y1 = h, x2 = xpos[b], y2 = h)
        }
        sim$dendro <- list(
          segs = lapply(segs, function(s) lapply(s, function(v) round(v, 4))),
          leaves = .xq_arr(colnames(sim_mat)[hc$order])
        )
      }
    }
    out$similarity <- sim
  }

  if (!is.null(dist_mat)) {
    lx <- log2(dist_mat + 1)
    expr <- list(groups = groups)
    box <- .xq_box_stats(lx)
    expr$box <- box[!vapply(box, is.null, logical(1))]
    med <- apply(lx, 1L, stats::median, na.rm = TRUE)
    rle <- .xq_box_stats(sweep(lx, 1L, med, "-"))
    expr$rle <- rle[!vapply(rle, is.null, logical(1))]
    expr$dens <- lapply(colnames(lx), function(s) {
      v <- lx[, s]
      v <- v[is.finite(v)]
      if (length(v) < 5L) return(NULL)
      de <- stats::density(v, n = 128L)
      list(name = s, x = .xq_arr(round(de$x, 3)), y = .xq_arr(round(de$y, 5)))
    })
    expr$dens <- expr$dens[!vapply(expr$dens, is.null, logical(1))]
    if (!length(expr$box)) expr$box <- NULL
    if (!length(expr$rle)) expr$rle <- NULL
    if (!length(expr$dens)) expr$dens <- NULL
    out$expr <- expr
  }

  mv_mat <- .rep_matrix_for_qc(ex$qc_matrix %||% ex$matrix)
  if (!is.null(mv_mat)) {
    means <- rowMeans(mv_mat, na.rm = TRUE)
    vars <- apply(mv_mat, 1L, stats::var, na.rm = TRUE)
    keep <- is.finite(means) & is.finite(vars)
    if (any(keep)) {
      out$mv <- data.frame(
        gene = rownames(mv_mat)[keep],
        lmean = log10(means[keep] + 1), lvar = log10(vars[keep] + 1),
        stringsAsFactors = FALSE
      )
    }
  }
  out
}

# ---- the collector ----------------------------------------------------------------------

#' Gather every table and view spec of the dashboard
#'
#' @return list(tables = key -> list(df, no_p), views = id -> spec,
#'   cfg = list, groups = ..., enrich = ..., regulon = ..., pca = ...)
#' @keywords internal
.xq_collect <- function(res, st) {
  ctx <- .xq_help_ctx(res, st)
  max_rows <- st$table_rows %||% Inf
  tabs <- list()
  views <- list()
  put_table <- function(key, df, no_p = character(0)) {
    if (is.null(df) || !is.data.frame(df) || !nrow(df) || !ncol(df)) return(FALSE)
    tabs[[key]] <<- list(df = df, no_p = no_p)
    TRUE
  }
  put_view <- function(id, spec) { views[[id]] <<- spec; invisible(id) }

  cfg <- res$config
  vcol <- .rep_hex(cfg$volcano_colors)
  pcol <- .rep_hex(cfg$pca_colors)
  out_cfg <- list(
    title = st$title %||% "ExpressOM report",
    comparison = if (is.na(res$comparison)) NULL else res$comparison,
    level = cfg$level, base = cfg$base, test = ctx$test,
    padj = ctx$padj, lfc = ctx$lfc,
    colors = list(up = unname(vcol[["up"]]), down = unname(vcol[["down"]]), ns = unname(vcol[["ns"]])),
    pcaColors = list(level = unname(pcol[["level"]]), base = unname(pcol[["base"]])),
    duckdb_url = st$duckdb_url
  )
  out <- list(cfg = out_cfg, ctx = ctx)

  # ---- differential expression -----------------------------------------------------
  de <- res$dge$results
  groups <- NULL
  if (!is.null(de)) {
    de <- de[order(de$padj, na.last = TRUE), , drop = FALSE]
    lead <- intersect(c("gene", "symbol", "ensembl"), names(de))
    de <- de[c(lead, setdiff(names(de), lead))]
    if (is.finite(max_rows) && nrow(de) > max_rows) de <- de[seq_len(max_rows), , drop = FALSE]
    put_table("de", de)
    show_cols <- intersect(c("gene", "symbol", "ensembl", "baseMean", "log2FoldChange", "lfcSE",
                             "pvalue", "padj", "contrast_log2FoldChange", "contrast_padj",
                             "gene_biotype"), names(de))
    put_view("volcano", list(
      view = "deExplorer", table = "de", idCol = "gene", lfcCol = "log2FoldChange", pCol = "padj",
      meanCol = if ("baseMean" %in% names(de)) "baseMean" else NULL,
      padj = ctx$padj, lfc = ctx$lfc, lrt = ctx$lrt, topN = 10L,
      columns = .xq_arr(show_cols),
      tipCols = .xq_arr(intersect(c("symbol", "ensembl", "gene_biotype"), names(de))),
      title = paste0("Volcano plot: ", cfg$level %||% "", " vs ", cfg$base %||% ""),
      filename = paste0("volcano_", res$comparison), minHeight = 620
    ))
    put_view("diag", list(
      view = "pDiag", table = "de",
      hists = list(list(col = "pvalue", title = "Raw p-value distribution")),
      qq = list(col = "pvalue", title = "p-value Q-Q plot"),
      ranked = list(idCol = "gene", lfcCol = "log2FoldChange", n = 25L), minHeight = 400
    ))
    put_view("de_table", list(
      view = "table", table = "de", sort = "padj", sqlName = "de",
      filename = paste0("DEgenes_", res$comparison),
      examples = .xq_arr(c("padj < 0.01 AND abs(log2FoldChange) > 1",
                           "baseMean > 100 AND padj < 0.05",
                           "gene LIKE 'HOX%'")),
      sqlExample = "SELECT sign(log2FoldChange) AS direction, count(*) AS n FROM de WHERE padj < 0.05 GROUP BY 1"
    ))
    if (!is.null(res$dge$deg_counts)) {
      put_table("deg_counts", res$dge$deg_counts)
      put_view("deg_counts", list(view = "degCounts", table = "deg_counts", lrt = ctx$lrt))
      put_view("deg_counts_table", list(view = "table", table = "deg_counts", filter = FALSE, sql = FALSE,
                                        search = FALSE, filename = "DEG_counts"))
    }

    ex <- res$dge$expression
    if (!is.null(ex)) {
      groups <- .xq_group_list(ex$groups, ex$samples)
      keep_n <- min(nrow(ex$matrix), 100L)
      hm <- as.data.frame(ex$matrix[seq_len(keep_n), , drop = FALSE], stringsAsFactors = FALSE,
                          check.names = FALSE)
      hm <- cbind(data.frame(gene = rownames(ex$matrix)[seq_len(keep_n)], stringsAsFactors = FALSE), hm)
      if (put_table("heat", hm, no_p = ex$samples)) {
        put_view("heat", list(
          view = "heatmap", table = "heat", idCol = "gene", samples = .xq_arr(ex$samples),
          groups = groups, groupName = cfg$main_condition %||% "group",
          n = min(as.integer(st$heatmap_genes %||% 20L), keep_n), nMax = keep_n
        ))
      }
      qc <- .xq_qc_payload(ex, groups)
      if (!is.null(qc$similarity)) put_view("similarity", c(list(view = "similarity"), qc$similarity))
      if (!is.null(qc$expr)) {
        spec <- c(list(view = "exprQC"), qc$expr)
        if (!is.null(qc$mv) && put_table("mv", qc$mv, no_p = names(qc$mv))) spec$mv <- "mv"
        put_view("expr_qc", spec)
      }
    }
  }
  out$groups <- groups

  # ---- PCA -------------------------------------------------------------------------
  pca_ids <- character(0)
  for (i in seq_along(res$pca)) {
    entry <- res$pca[[i]]
    key <- paste0("pca_", i)
    sc <- entry$scores
    if (!put_table(key, sc, no_p = names(sc))) next
    gcol <- .rep_pca_group_col(sc, cfg$main_condition)
    pcs <- utils::head(entry$pcs, 10L)
    put_view(key, list(
      view = "pca", table = key, pcs = .xq_arr(pcs),
      variance = if (!is.null(entry$variance)) as.list(entry$variance[intersect(names(entry$variance), pcs)]) else NULL,
      groupCol = if (is.na(gcol)) NULL else gcol,
      labelCol = if ("plot_label" %in% names(sc)) "plot_label" else "sample_label",
      title = .xp_pca_title(entry$name), filename = entry$name
    ))
    pca_ids <- c(pca_ids, key)
  }
  out$pca <- pca_ids

  # ---- functional enrichment (one merged table per source + method) -----------------
  enrich <- list()
  inv <- res$functional$inventory
  if (nrow(inv)) {
    n_rank <- .xq_n_ranked(res)
    grp <- paste(inv$source, inv$method, sep = "\u001f")
    k <- 0L
    for (g in unique(grp)) {
      sub <- inv[grp == g, , drop = FALSE]
      parts <- lapply(seq_len(nrow(sub)), function(r) {
        tb <- res$functional$tables[[sub$path[[r]]]]
        if (is.null(tb) || !nrow(tb)) return(NULL)
        .xq_enrich_std(tb, sub$collection[[r]])
      })
      merged <- .xq_rbind_fill(parts)
      if (is.null(merged) || !nrow(merged)) next
      if (is.finite(max_rows) && nrow(merged) > max_rows) merged <- merged[seq_len(max_rows), , drop = FALSE]
      k <- k + 1L
      key <- paste0("enr_", k)
      put_table(key, merged)
      kind <- sub$kind[[1L]]
      id <- key
      put_view(id, list(
        view = "enrichment", table = key, kind = kind, method = sub$method[[1L]],
        source = sub$source[[1L]], collections = .xq_arr(unique(sub$collection)),
        padj = ctx$padj, topN = as.integer(st$top_terms %||% 15L),
        rank = if (identical(kind, "gsea") && !is.null(de)) list(table = "de", geneCol = "gene", statCols = .xq_arr(c("contrast_stat", "stat"))) else NULL,
        filename = paste0("enrichment_", gsub("[^A-Za-z0-9]+", "_", sub$method[[1L]])),
        minHeight = 560
      ))
      enrich[[length(enrich) + 1L]] <- list(id = id, source = sub$source[[1L]], method = sub$method[[1L]], kind = kind)
    }
    out$n_ranked <- n_rank
  }
  out$enrich <- enrich

  # ---- decoupleR activity + regulons ------------------------------------------------
  dec <- list()
  for (i in seq_along(res$decouple)) {
    entry <- res$decouple[[i]]
    key <- paste0("dec_", i)
    if (!put_table(key, entry$table)) next
    put_view(key, list(view = "decouple", table = key, resource = entry$resource,
                       padj = ctx$padj, topN = as.integer(st$top_terms %||% 15L),
                       title = paste0("TF activity - ", entry$resource),
                       filename = paste0("decoupler_", entry$resource), minHeight = 480))
    dec[[length(dec) + 1L]] <- list(id = key, resource = entry$resource, path = entry$path)
  }
  # The network-status table (live / offline copy / failed per resource)
  # is shown as a plain table card; it is not part of res$decouple because it
  # has no activity columns.
  if (!is.null(res$decouple_status) &&
      put_table("dec_status", res$decouple_status, no_p = names(res$decouple_status))) {
    put_view("dec_status", list(view = "table", table = "dec_status",
                                filename = "decoupleR_network_status", sql = FALSE))
    out$dec_status <- TRUE
  }
  out$decouple <- dec

  reg <- list()
  rr <- res$decouple_regulons$by_resource
  for (rn in names(rr)) {
    entry <- rr[[rn]]
    regs <- lapply(names(entry$regulons), function(tf) {
      d <- entry$regulons[[tf]]
      if (is.null(d) || !nrow(d)) return(NULL)
      cbind(data.frame(tf = rep(tf, nrow(d)), stringsAsFactors = FALSE), as.data.frame(d, stringsAsFactors = FALSE))
    })
    rt <- .xq_rbind_fill(regs)
    need <- c("tf", "target", "log2FoldChange", "padj", "Regulation")
    if (is.null(rt) || !all(need %in% names(rt))) next
    les <- lapply(names(entry$leading_edges), function(tf) {
      d <- entry$leading_edges[[tf]]
      if (is.null(d) || !nrow(d)) return(NULL)
      cbind(data.frame(tf = rep(tf, nrow(d)), stringsAsFactors = FALSE), as.data.frame(d, stringsAsFactors = FALSE))
    })
    lt <- .xq_rbind_fill(les)
    key <- paste0("reg_", length(reg) + 1L)
    put_table(key, rt)
    lkey <- NULL
    if (!is.null(lt) && all(c("tf", "rank") %in% names(lt)) && put_table(paste0("le_", length(reg) + 1L), lt, no_p = names(lt))) {
      lkey <- paste0("le_", length(reg) + 1L)
    }
    put_view(key, list(view = "regulon", reg = key, le = lkey, tfs = .xq_arr(unique(rt$tf)),
                       nRanked = .xq_n_ranked(res), padj = ctx$padj, lfc = ctx$lfc, minHeight = 700))
    reg[[length(reg) + 1L]] <- list(id = key, resource = rn, n_tf = length(unique(rt$tf)))
  }
  out$regulon <- reg
  if (!is.null(res$decouple_regulons$summary) &&
      put_table("reg_summary", res$decouple_regulons$summary)) {
    put_view("reg_summary", list(view = "table", table = "reg_summary", filename = "tf_regulon_summary", sql = FALSE))
    out$reg_summary <- TRUE
  }

  # ---- isoform analysis -------------------------------------------------------------
  iso <- res$isoform
  if (!is.null(iso)) {
    ctab <- function(key, df, sort = NULL) {
      if (is.null(df)) return(FALSE)
      if (is.finite(max_rows) && nrow(df) > max_rows) df <- df[seq_len(max_rows), , drop = FALSE]
      put_table(key, df)
    }
    if (ctab("sw", iso$switch_features)) put_view("sw", list(view = "switch", table = "sw", padj = ctx$padj, dif = 0.1, minHeight = 700))
    if (ctab("dte", iso$dte)) {
      put_view("dte", list(view = "dte", table = "dte", padj = ctx$padj, lfc = ctx$lfc, minHeight = 900))
      put_view("dte_table", list(view = "table", table = "dte", filename = "DTE_results", sqlName = "dte",
                                 sort = if ("padj" %in% names(iso$dte)) "padj" else NULL))
    }
    if (ctab("dtu", iso$dtu)) {
      put_view("dtu", list(view = "dtu", table = "dtu", padj = ctx$padj, minHeight = 900))
      pc <- .rep_first_col(iso$dtu, c("adj_pvalue", "padj", "FDR", "qvalue", "p.adjust"))
      put_view("dtu_table", list(view = "table", table = "dtu", filename = "DTU_results", sqlName = "dtu",
                                 sort = if (is.na(pc)) NULL else pc))
    }
    if (ctab("dex", iso$dexseq)) {
      pc <- .rep_first_col(iso$dexseq, c("padj", "adj_pvalue", "FDR", "qvalue", "pvalue"))
      put_view("dex", list(view = "table", table = "dex", filename = "DEXSeq_results", sqlName = "dex",
                           sort = if (is.na(pc)) NULL else pc))
    }
    if (ctab("txpca", iso$transcript_pca)) {
      put_view("txpca", list(view = "table", table = "txpca", filename = "PCA_transcripts", sqlName = "txpca"))
    }
  }

  # ---- run info ---------------------------------------------------------------------
  if (put_table("files", res$files)) {
    put_view("files", list(view = "table", table = "files", sort = "path", pageSize = 25L, sqlName = "files",
                           filename = "file_inventory",
                           examples = .xq_arr(c("category = 'ORA'", "size_kb > 1000"))))
  }

  # free SQL over the main tables
  sqlt <- list()
  for (nm in c("de", "dte", "dtu", "dex", "sw", "files")) if (!is.null(tabs[[nm]])) sqlt[[nm]] <- nm
  for (e in enrich) sqlt[[e$id]] <- e$id
  if (length(sqlt)) {
    put_view("sql", list(
      view = "sql", tables = sqlt,
      intro = paste0("Available tables: ", paste(names(sqlt), collapse = ", "), "."),
      sql = if (!is.null(tabs$de)) "SELECT gene, baseMean, log2FoldChange, padj\nFROM de\nWHERE padj < 0.05\nORDER BY abs(log2FoldChange) DESC\nLIMIT 50" else NULL
    ))
  }

  out$tables <- tabs
  out$views <- views
  out
}
