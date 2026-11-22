# mod_report_quarto_help.R -- the "How to read this" note shown at the top of
# every tab / card of the Quarto dashboard.
#
# Each note is static HTML (a collapsible <details class="xq-help">), so it is
# readable even before any JavaScript runs, and it is written *for the person
# looking at the plot*: what is on each axis, what a "good" or "worrying"
# picture looks like, and what the colours / thresholds mean. Values that
# depend on the run (contrast, cut-offs, test) are filled in from `ctx`.

.xq_e <- function(x) htmltools::htmlEscape(as.character(x))

#' Context used to word the notes
#' @keywords internal
.xq_help_ctx <- function(res, st) {
  cfg <- res$config
  co <- .xp_cutoffs(res, st)
  list(
    level = cfg$level, base = cfg$base, padj = co$padj, lfc = co$lfc,
    test = cfg$test %||% "Wald", lrt = identical(cfg$test, "LRT"),
    group = cfg$main_condition
  )
}

#' Build one collapsible note
#' @param bullets Character vector of list items (inline HTML allowed; escape
#'   dynamic values with `.xq_e()`).
#' @param intro Optional lead sentence.
#' @param title Summary line.
#' @keywords internal
.xq_help <- function(bullets, intro = NULL, title = "How to read this", open = TRUE) {
  paste0(
    "<details class=\"xq-help\"", if (open) " open" else "", ">",
    "<summary><span class=\"xq-help-ic\" aria-hidden=\"true\">i</span>", .xq_e(title), "</summary>",
    "<div class=\"xq-help-body\">",
    if (!is.null(intro)) paste0("<p>", intro, "</p>") else "",
    if (length(bullets)) paste0("<ul>", paste0("<li>", bullets, "</li>", collapse = ""), "</ul>") else "",
    "</div></details>"
  )
}

# direction wording shared by several notes
.xq_dir_text <- function(ctx) {
  if (!is.null(ctx$level) && !is.null(ctx$base)) {
    paste0("a <b>positive</b> log2 fold change means higher expression in <b>", .xq_e(ctx$level),
           "</b> than in <b>", .xq_e(ctx$base), "</b>, a <b>negative</b> one the opposite")
  } else {
    "a positive log2 fold change means higher expression in the tested group than in the reference group"
  }
}

#' The help-text dictionary
#'
#' @param key Which card the note belongs to.
#' @param ctx From `.xq_help_ctx()`.
#' @param ... Extra values some notes use (`name`, `kind`, `resource`, `tf`).
#' @return One HTML string.
#' @keywords internal
.xq_help_html <- function(key, ctx, ...) {
  a <- list(...)
  padj <- .xq_e(ctx$padj)
  lfc <- .xq_e(ctx$lfc)
  dir <- .xq_dir_text(ctx)
  switch(
    key,

    # ---------------------------------------------------------------- overview
    overview = .xq_help(c(
      "The tiles summarise the main contrast: how many genes could be tested and how many pass the adjusted p-value and fold-change cut-offs (<b>up</b> = higher in the tested group, <b>down</b> = lower).",
      "Counts depend on the cut-offs, so read them together with the <em>DEGs by cut-off</em> card and the volcano plot rather than as a single \"number of genes\".",
      "Every page is a set of tabs. Each tab starts with a note like this one; use <b>Hide hints</b> (bottom right) once you know the layout.",
      "Tabs load only when you open them, tables filter with plain SQL conditions, and every chart can be exported (PNG/SVG, or CSV for tables)."
    ), title = "How to use this report"),

    deg_counts = .xq_help(c(
      "Each group of bars is one adjusted p-value cut-off (smaller = stricter); bar height is the number of genes that pass it together with the |log2 fold change| threshold shown in the tooltip.",
      if (ctx$lrt) "The LRT tests <em>any</em> change across the design, so it has no direction: one bar per cut-off." else "Red bars count up-regulated genes, blue bars down-regulated ones.",
      "A steep drop between neighbouring cut-offs means many genes sit just below the usual 0.05 line; a stable count means the call is robust to the threshold."
    ), title = "How to read this chart"),

    contents = .xq_help(c(
      "A tick means the corresponding result files were found in the run folder and are shown in this report; a dash means the step was not run (or produced no output).",
      "Everything here is computed by the ExpressOM pipeline in R; this dashboard only displays, filters and re-plots those results. The numbers in the tables are the ones in the CSV files of the run."
    ), title = "What is in this report"),

    # ------------------------------------------------------------------- QC
    pca = {
      nm <- a$name %||% ""
      kind <- if (grepl("BatchCorrected", nm)) "batch-corrected expression: batch effects have been removed, so any remaining separation should be biological"
              else if (grepl("ComparisonOnly", nm)) "only the samples of this contrast, so the axes capture the variation that matters for the comparison"
              else "all samples, before any batch correction (exploratory)"
      .xq_help(c(
        paste0("Each dot is a <b>sample</b>; the plot is built from ", kind, "."),
        "Samples that are <b>close</b> have similar expression profiles. Replicates of one group should sit together, and the groups should separate along PC1 or PC2 if the condition is the main source of variation.",
        "The axis labels give the share of the total variance carried by each component: a large PC1 (&gt; 30-40 %) is usually the strongest effect. Use the X / Y selectors to look at other components.",
        "A sample far from its own group, or groups that split along a technical variable (batch, library size, date) instead of the condition, are warning signs to check before trusting the DE results.",
        "The bar chart is a scree plot: components with a clearly higher bar than the rest are the ones worth interpreting."
      ))
    },

    similarity = .xq_help(c(
      "All four views use the most variable genes (log2 of normalised counts + 1) and ask the same question: <b>which samples look alike?</b>",
      "<b>Correlation heat map</b> - dark = highly correlated. Replicates of the same group should form dark blocks (typically r &gt; 0.9); a sample that is paler than its group-mates is an outlier candidate.",
      "<b>Distance heat map</b> - the mirror image (small distance = similar). Samples are ordered by clustering, so blocks along the diagonal are groups.",
      "<b>MDS</b> - a 2-D map of the distances, comparable to a PCA: distance on the plot approximates distance in expression.",
      "<b>Dendrogram</b> - samples that join low are most similar. Group members should merge with each other before merging with other groups; a sample that joins its group only near the top is a possible outlier or mislabel.",
      "This is exploratory: use it with the metadata (batches, preparation dates) before excluding anything."
    )),

    expr_qc = .xq_help(c(
      "Values are log2(normalised count + 1) over the same set of genes for every sample, so profiles can be compared directly.",
      "<b>Box plot / density</b> - after normalisation, medians, spread and curve shapes should look alike. One sample shifted sideways or with an extra bump (e.g. a spike of zeros or a rRNA / globin peak) points to a library or normalisation problem.",
      "<b>RLE</b> (relative log expression: each gene minus its median across samples) - boxes should be centred on 0 with similar height. A box that is off zero or much taller than the others indicates a sample-wide bias.",
      "<b>Mean-variance</b> - variability grows with expression; for count data on a log scale the cloud should be smooth. Isolated strong outliers are genes with unstable counts, not a sample problem."
    )),

    # ------------------------------------------------------------------ DGE
    volcano = .xq_help(c(
      paste0("<b>Volcano</b>: x = log2 fold change (", dir, "); y = -log10 adjusted p-value, so <b>higher = more significant</b>. The most interesting genes are in the top-left and top-right corners."),
      if (ctx$lrt) paste0("The dashed horizontal line is the adjusted p-value cut-off (", padj, "). The LRT tests the whole model, so the fold change shown is only the contrast estimate and colours simply mark the significant genes.")
      else paste0("Dashed lines are the cut-offs (adj. p &lt; ", padj, ", |log2FC| &ge; ", lfc, "). <span style=\"color:#c94a4a\"><b>Red</b></span> = significantly up, <span style=\"color:#2b5c8f\"><b>blue</b></span> = significantly down, grey = not significant."),
      "<b>MA plot</b> (switch with the View buttons): x = mean expression, y = log2 fold change. Large fold changes at very low expression are usually noise; trustworthy hits stand out at moderate or high expression.",
      "Move the cut-off controls to see how many genes pass; the tiles and the table below update instantly. Hover a point for its values, drag to pan, scroll to zoom, <b>Lasso</b> (or Alt-drag) to list a group of genes, and type gene names in <em>Find</em> to ring them.",
      "Fold changes are estimates (shrunk when the run used LFC shrinkage); adjusted p-values control the false-discovery rate across all tested genes."
    )),

    diagnostics = .xq_help(c(
      "<b>p-value histogram</b> - with no effect the bars are flat (dashed line). A flat floor with a <b>spike on the left</b> is the healthy picture: mostly null genes plus real signal. A U shape, a pile-up on the right, or a hump in the middle suggests a problem with the model, filtering or hidden batch effects.",
      "<b>Q-Q plot</b> - observed vs expected p-values on a -log10 scale. Points on the diagonal follow the null; a tail that lifts off at the top right is real signal. Points that leave the diagonal early indicate inflation (confounding or an under-fitted model).",
      "<b>Largest effects</b> - the genes with the greatest |log2 fold change|. They can be low-count genes, so always check their adjusted p-value in the table before drawing conclusions."
    )),

    heatmap = .xq_help(c(
      "Rows are the most significant genes, columns are samples. Colour is a <b>row z-score</b> of log2(normalised count + 1): red = above that gene's average across samples, blue = below, white = average.",
      "z-scores show the <b>pattern</b> of each gene, not its absolute level, so a gene with low counts can look as strong as a highly expressed one.",
      "Rows and columns are clustered. Samples should group by condition (see the coloured bar on top); genes forming a red-blue block share one pattern, e.g. up in one group and down in the other.",
      "Choose how many genes to display and whether columns are clustered or grouped by condition with the controls above the plot."
    )),

    de_table = .xq_help(c(
      "One row per gene, sorted by adjusted p-value. <b>baseMean</b> = mean normalised count, <b>log2FoldChange</b> = effect size, <b>lfcSE</b> its standard error, <b>pvalue</b> the raw test, <b>padj</b> the Benjamini-Hochberg adjusted value.",
      "An empty (NA) <b>padj</b> means the gene was removed by independent filtering (too few counts) or flagged as an outlier - it was not tested, which is different from \"not significant\".",
      "Click a header to sort. Type a SQL condition to filter, for example <code>padj &lt; 0.01 AND abs(log2FoldChange) &gt; 1</code> or <code>gene LIKE 'HOX%'</code>; the CSV button exports exactly the rows you filtered."
    )),

    gene_lists = .xq_help(c(
      "These are the gene lists the pipeline wrote for each adjusted p-value cut-off (one gene per line). Use them as input for other tools or for the enrichment analyses.",
      "Links are relative to the run folder: they work when this HTML file stays next to the result folders."
    )),

    sql_de = .xq_help(c(
      "Full SQL (DuckDB-Wasm) over the results table <code>de</code>: GROUP BY, joins, window functions, quantiles.",
      "The engine loads the first time you press <b>Run SQL</b> and needs internet access (or a local copy set with <code>duckdb_url</code>). The WHERE filters elsewhere in the report never need it."
    ), title = "About the SQL explorer"),

    # ----------------------------------------------------------- enrichment
    enrichment = {
      kind <- a$kind %||% "ora"
      method <- a$method %||% ""
      common <- "Use <b>Collection</b> to switch between gene-set databases, <b>Terms shown</b> to change how many are plotted and the table at the bottom (SQL filter, CSV export) for all rows."
      switch(
        kind,
        gsea = .xq_help(c(
          paste0("<b>GSEA / FGSEA</b> (", .xq_e(method), ") ranks <em>all</em> genes by their test statistic and asks whether the members of a gene set sit towards one end of that list, without any cut-off."),
          paste0("<b>NES</b> (normalised enrichment score): <span style=\"color:#c94a4a\"><b>positive</b></span> = the set is concentrated among genes up in ", .xq_e(ctx$level %||% "the tested group"), "; <span style=\"color:#2b5c8f\"><b>negative</b></span> = among genes that are down. Bigger |NES| = stronger enrichment. The bar chart shows the terms with the smallest adjusted p-values."),
          "<b>Running score</b> (below the plots): the curve climbs each time a gene of the set is met going down the ranked list and falls otherwise. A single early peak means the set is driven by the top-ranked genes; ticks mark the set members. Only the leading-edge genes reported by the tool are drawn.",
          "<b>Gene overlap / map</b>: terms that share many genes (high Jaccard) are redundant and tell the same biological story; tight clusters in the network are one theme.",
          common
        )),
        spia = .xq_help(c(
          "<b>SPIA</b> combines two pieces of evidence for a signalling pathway: whether many pathway genes are differentially expressed and whether the <em>perturbation</em> propagated through the pathway topology is larger than expected.",
          "The adjusted value shown is the combined global p-value (pGFdr). <b>tA</b> is the total perturbation accumulation: positive = pathway predicted <b>activated</b>, negative = <b>inhibited</b> (see also the Status column).",
          common
        )),
        tf = .xq_help(c(
          "Transcription factors whose known targets are over-represented among the differentially expressed genes (or enriched at one end of the ranked list).",
          "A small adjusted p-value means the targets of that factor change more than expected by chance; it suggests, but does not prove, that the factor's activity changed. Check the direction of the targets in the regulon explorer when available.",
          common
        )),
        .xq_help(c(
          paste0("<b>Over-representation analysis</b> (", .xq_e(method), "): takes the list of differentially expressed genes and tests whether a pathway / GO term contains more of them than expected by chance."),
          "<b>Ranked bars</b>: -log10 adjusted p-value, so a longer bar = stronger enrichment. <b>Bubble plot</b>: x = number of DE genes in the term, bubble size = size of the term, colour = significance.",
          "Results depend on the DE gene list (cut-offs) and on the background gene set used by the tool; very large generic terms are easy to enrich, specific small terms are more informative.",
          "<b>Gene overlap / map</b>: terms that share many genes (high Jaccard) are redundant and tell the same biological story; tight clusters in the network are one theme.",
          common
        ))
      )
    },

    decouple = {
      resource <- a$resource %||% ""
      .xq_help(c(
        paste0("decoupleR / <b>", .xq_e(resource), "</b>: each bar is one transcription factor (or one pathway, for PROGENy) and its length is the inferred <b>activity score</b>, a t-like statistic computed from the signed gene-level statistics of the source's target genes (its regulon)."),
        "A <b>positive</b> score (red) means the source behaves as more active than chance: its activated targets go up and/or its repressed targets go down. A <b>negative</b> score (blue) is the opposite, a net repressive footprint.",
        paste0("By default only sources with p &lt; ", padj, " are shown, ordered by |score|; untick the switch to see the strongest ones regardless of significance. The table lists every tested source."),
        "Activity is inferred from the targets' expression, not measured directly; it is more reliable when many targets of the source are present in the data."
      ))
    },

    regulon = .xq_help(c(
      "Pick a transcription factor (or source) to see its regulon, i.e. its known direct targets, in two views.",
      paste0("<b>Volcano of targets</b>: only the target genes are plotted. <span style=\"color:#1f78b4\"><b>Activated (expected)</b></span> = the target went up and the factor is an activator (mor &gt; 0); <span style=\"color:#e31a1c\"><b>Repressed (expected)</b></span> = it went down and the factor is a repressor (mor &lt; 0); orange = significant change opposite to the network's expectation; grey = not significant (adj. p &lt; ", padj, ", |log2FC| &ge; ", lfc, "). A majority of expected-direction points agrees with a non-zero activity score."),
      "<b>Leading edge</b>: all genes are ranked by their test statistic (left = most up). The curve is a running enrichment score that rises at every target; red ticks mark where the targets sit. A peak on the left means targets are concentrated among up-regulated genes (activating role), a trough on the right among down-regulated genes (repressing role).",
      "The tables list every matched target and the leading-edge genes in rank order."
    ), title = "How to read this explorer"),

    # ------------------------------------------------------------ isoforms
    switch = .xq_help(c(
      "An <b>isoform switch</b> is a change in which transcript of a gene is mainly used. <b>dIF</b> (x axis) is the change in the isoform's share of its gene's expression between the conditions: from -1 to +1, with a positive value meaning the isoform is used more in the tested group.",
      "y = -log10 of the switch q-value (higher = more significant). Points beyond both the q and |dIF| thresholds (dashed lines) are called <b>switches</b> (red). A large dIF with a weak q-value is an unreliable effect; a tiny dIF with a strong q-value is statistically clear but biologically small.",
      "The table lists the significant switches; click a point legend entry to hide a class and use the controls to tighten or relax the call."
    )),

    dte = .xq_help(c(
      paste0("<b>Differential transcript expression</b>: the same test as for genes, but every point is one <b>transcript</b>. Read the volcano and MA plot exactly as in the gene-level tab (", dir, ")."),
      "<b>Abundance vs significance</b>: highly abundant transcripts can reach very small p-values with modest changes; low-abundance ones need big changes. A cloud of significant points only at high abundance is expected.",
      "<b>Top transcripts</b> are ranked by adjusted p-value (red = up, blue = down). The p-value histograms should be flat with a spike near 0 when there is signal.",
      "A gene can have DTE in one isoform without a change in the gene total, so compare with the gene-level results and with the DTU tab."
    )),

    dtu = .xq_help(c(
      "<b>Differential transcript usage</b> asks whether the <em>proportions</em> of a gene's transcripts change between conditions, even if total gene expression does not. The test reports a <b>likelihood ratio (LR)</b>; this table has no log2 fold change, so there is no volcano.",
      "<b>LR vs adjusted p-value</b>: larger LR = stronger evidence; points above the dashed line pass the cut-off. <b>LR vs degrees of freedom</b>: genes with more transcripts have more degrees of freedom and larger LR values by construction.",
      "<b>Top features</b> are the strongest DTU calls; the <b>gene summary</b> counts how many tested features a gene has (a descriptive count, not a gene-level FDR).",
      "The p-value histograms should be flat with a left spike when there is signal. Confirm interesting genes in the Switch overview or the isoform plots."
    )),

    dexseq = .xq_help(c(
      "DEXSeq treats each transcript as an exon-like bin and tests whether its share of the gene's counts differs between conditions (differential usage).",
      "Look at the <b>adjusted p-value</b> to select features and at the log2 fold-change column (when present) for the direction of the usage change. Several features of one gene can be significant together when usage moves from one isoform to another."
    )),

    tx_pca = .xq_help(c(
      "The most variable transcripts used for the transcript-level PCA. Their loadings tell you which transcripts drive the separation of samples seen in the transcript PCA.",
      "Sort by the loading columns to find the transcripts that contribute most to a component."
    )),

    iso_plots = .xq_help(c(
      "Static plots produced by the isoform analysis (volcanoes, switch plots, transcript structures and usage). Click a picture to enlarge it.",
      "They show the same results as the interactive tabs, drawn by the original tools, so they keep tool-specific annotations such as domains or coding potential."
    )),

    iso_methods = .xq_help(c(
      "Statistical notes recorded by the isoform analysis: which test, which adjustment and which filters were used. Read them before comparing DTE, DTU and DEXSeq results, because the units of the effect sizes differ."
    ), title = "About these notes"),

    # ------------------------------------------------------------- figures
    figures = .xq_help(c(
      "Static figures written by the pipeline (PDF pages are shown as low-resolution previews, PNG files as they are). Click a picture to enlarge it; the originals are in the run folder.",
      "These are the publication-style versions of the interactive plots, useful for slides and papers."
    )),

    # ------------------------------------------------------------ run info
    params = .xq_help(c(
      "How the analysis was run: the contrast, the DESeq2 design formula, the test (Wald or LRT), LFC shrinkage and the cut-offs used for the DE gene lists. These values decide how to read every other tab.",
      "The <b>design formula</b> lists the factors that were modelled; effects of listed covariates (batch, sex...) have been accounted for in the fold changes. With an LRT the <b>reduced formula</b> is what the full model was compared against."
    )),
    logs = .xq_help(c(
      "Messages recorded while the pipeline ran. <b>Warnings</b> are worth a quick look (for example genes with extreme counts or factors that were collapsed); <b>session information</b> lists package versions for reproducibility.",
      "Log files of the optional steps (isoform, WSL predictors, OmniPath, targets) are listed as links."
    )),
    files = .xq_help(c(
      "Every file the run produced, with its size and modification time. Click a path to open it (relative to this report), or filter by folder with a SQL condition such as <code>category = 'ORA'</code>."
    )),

    stop("No help text for key: ", key)
  )
}
