<p align="center">
  <img src="assets/logo.png" alt="ExpressOM logo" width="350"/>
</p>

<h1 align="center">ExpressOM</h1>

<p align="center">
  Bulk RNA-seq analysis toolkit
</p>

<p align="center">
  <img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="GPL-3 License"/>
</p>

**ExpressOM** is a fully modularized R package for bulk RNA-seq analysis. It wraps `tximport`, `DESeq2`, `clusterProfiler`, `SPIA`, `fgsea`, and `Enrichr` into a streamlined, automated workflow. 

It effortlessly handles counts import, exploratory data analysis (EDA), differential expression analysis, comprehensive visualizations, and functional enrichment (GO, Reactome, Disease Ontology, KEGG, Transcription Factors, SPIA, and GSEA). Species currently supported: **Human** and **Mouse**.

---

## Installation

You can install the package directly from your local source directory using `devtools` or `remotes`:

```R
# Install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install directly from GitHub
devtools::install_github("lavauxt/ExpressOM")
```

### 1. Generate and Install the Custom Ensembl Database
`ExpressOM` relies on a custom-built Ensembl database package to annotate transcripts and map identifiers securely. You must build and install it locally before processing your data.

```R
library(ExpressOM)

ExpressOM::run_bulk_pipeline(
  ensembl_package_name = "EnsDb.Hsapiens.v107", 
  data_dir = "./raw_data",
  out_dir = "./results"
)
```

### 1. Automated Ensembl Annotation Database Handling

The pipeline fully automates annotation package installation. If the package assigned to the `ensembl_package_name` parameter (e.g., `"EnsDb.Hsapiens.v107"`) is missing from the local R environment, the runtime automatically:
1. Parses out the target species (`human` or `mouse`) and Ensembl release version (`107`).
2. Fetches appropriate remote `.gtf` files and generates a custom backend mapping table (`create_homemade_db`).
3. Bundles and local-installs the structured R package wrapper securely on the fly.

> 💡 **Tip:** Use https://www.gencodegenes.org/ to check for proper version of Ensembl fur Human and Mouse

Every run starts with a pre-flight dependency check (`validate_environment()`): missing core packages (`DESeq2`, `tximport`, `dplyr`, `ggplot2`, `pheatmap`) stop the run immediately with a clear install message, while missing functional-analysis or isoform-analysis packages (only relevant if `run_dge`/`run_isoform` are enabled) produce an early warning instead of failing hours into a run.

---

## Data Preparation

The pipeline expects a specific directory structure for your quantification tools (`salmon`, `kallisto`, `rsem`, `stringtie`, etc.) or raw count matrices (`matrix` mode).

**Sample Table (`sample_table.csv` / `.tsv`) Format:**
Your sample table must contain a column identifying each sample (sample_id or Sample), along with factors used in your model (e.g., genotype, cell type, batch). This format is directly compatible with DESeq2, which uses it to perform differential expression analysis (DEG) in the pipeline.

| sample_id | genotype | cell_type  | batch |
|-----------|----------|------------|-------|
| Sample1   | WT       | T_cell     | 1     |
| Sample2   | KO       | T_cell     | 1     |
| Sample3   | WT       | B_cell     | 2     |

- **`sample_id`**: Unique identifier for each sample. This should match the filenames of your count or expression data.  
- **Factors (`genotype`, `cell_type`, `batch`)**: Variables to include in your DESeq2 model for differential expression analysis. These are used to build the design formula (e.g., `~ batch + cell_type + genotype`).  
- **File type**: Accepts **CSV** or **TSV**.  

> 💡 **Tip:** DESeq2 expects this table to be complete, with no missing sample IDs, and all factor columns properly formatted (as `factor` in R).

---

## Usage

### Example 1: Matrix Input (e.g. standard Read Counts)
Used when importing pre-quantified text files (like `raw_counts.csv`).

```R
library(ExpressOM)

run_bulk_pipeline(
  ensembl_package_name = "EnsDb.Mmusculus.v114",        # Ensure you use the exact database name
  count_type = "matrix",
  matrix_file = "./data/matrix_raw_counts.csv",
  sample_table = "./data/sample_table.tsv",
  out_dir = "./result_matrix",
  model = "~ cell_type + genotype + cell_type:genotype", # Design formula DESeq2 style
  level = "KO",                                          # Active/Treatment group
  base = "WT",                                           # Reference/Baseline group
  padj_cutoff = 0.05
)
```

### Example 2: Transcriptomic Input (`salmon`, `kallisto`, `rsem`)
Provide the base `data_dir` containing subfolders matching your `sample_id` column. `tximport` will automatically look for the expected file name (e.g. `quant.sf` for `salmon`, `abundance.tsv` for `kallisto`).

```R
run_bulk_pipeline(
  ensembl_package_name = "EnsDb.Hsapiens.v107", # Ensure you use the exact database name
  count_type = "salmon",                        # Options: "salmon", "kallisto", "rsem", "stringtie", etc.
  data_dir = "./data/counts",  
  out_dir = "./results_classic",
  sample_table = "./data/sample_table.csv",
  model = "~ condition",                        # Design formula DESeq2 style
  level = "Treated",                            # Active/Treatment group
  base = "Control"                              # Reference/Baseline group
)
```

## Multi-GMT Pathway Profiling (FGSEA)

The gmt_file parameter accepts a string vector or a standard list object mapping multiple independent gene matrices simultaneously. The workflow loops through each matrix file and saves separate outputs under dedicated subfolders named after the file

```R
# Profile multiple collections in one execution pass
ExpressOM::run_bulk_pipeline(
  gmt_file = c(
    "./gene_sets/h.all.v2023.1.Hs.symbols.gmt",
    "./gene_sets/c2.cp.kegg.v2023.1.Hs.symbols.gmt"
  ),
  out_dir = "./results"
)
```

# Resulting Workspace Tree
```txt
results/
└── GSEA/
    └── FGSEA/
        ├── h.all.v2023.1.Hs.symbols/
        │   ├── Barplot_h.all.v2023.1.Hs.symbols.pdf
        │   ├── GSEA_Table_h.all.v2023.1.Hs.symbols.html
        │   └── GSEA_plot_...pdf
        └── c2.cp.kegg.v2023.1.Hs.symbols/
            ├── Barplot_c2.cp.kegg.v2023.1.Hs.symbols.pdf
            ├── GSEA_Table_c2.cp.kegg.v2023.1.Hs.symbols.html
            └── GSEA_plot_...pdf
```
## Pipeline Outputs

The `run_bulk_pipeline` function generates a heavily organized output structure in your designated `out_dir`:

* **`DE_raw_results/`**: TSV tables of raw/filtered differential expression results and normalized counts.
* **`Plots/`**: PCA plots, Sample Correlation Heatmaps, MA plots, Volcano plots, and top DE gene boxplots.
* **`ORA/` & `GSEA/`**: Extensive targets and plots for GO mapping, Reactome, Disease Ontology, and KEGG generic pathways (Dotplots, Ridgeplots, Pathway Graphs).
* **`SPIA/`**: Signaling Pathway Impact Analysis graphs and Evidence CSVs.
* **`Transcription_Factors/`**: Output mapping active transcription factor perturbations based on dynamically mapped Enrichr libraries.
* **`RegionReport/`**: Auto-generated interactive `RegionReport` HTML documents. Figures embedded in this report (and in the isoform `DTU_DTE_report/report.html`) are automatically converted from PDF to PNG for display, since browsers can't render a PDF through an `<img>` tag; install the `pdftools` package if you want these figures to appear inline (otherwise the underlying `.pdf` plot files are still written to `Plots/`, just not embedded in the HTML report).
* **`Save_rdata/`**: The complete populated R environment serialized as an `.RData` file.
* **`Log/`**: A single, unified log tree for the whole run — always check here first if something looks like it didn't run:
  * `Log/SessionInfo.txt`, `Log/Warnings.txt` — captured at the very end of the run, regardless of what was enabled.
  * `Log/DGE/DGE_params_<comparison>.json` — one file per DGE comparison (design formula, contrast, filters, up/down counts).
  * `Log/Isoform/wsl_debug.json` — output of `debug_wsl()`: which predictor tools were found, whether the Pfam database is present *and* `hmmpress`-indexed, whether the `isoform_tools` conda env exists, etc. Re-run `debug_wsl()` any time to refresh this without re-running the whole pipeline.
  * `Log/Isoform/wsl_commands.log` — a timestamped audit trail of every command run inside WSL/bash for CPAT, SignalP, hmmscan, and the install helpers, with exit codes and captured output. This is the first place to look if a predictor "silently" didn't produce results.

## Windows / WSL Troubleshooting

If CPAT, SignalP, or Pfam/hmmscan results are missing after a run with `run_predictors = TRUE`, check `Log/Isoform/wsl_debug.json` and `Log/Isoform/wsl_commands.log` first — they will usually show exactly which tool or database was missing or failed, and why. Common causes that are now checked and reported explicitly:

* **Pfam database found but not indexed.** `hmmscan` requires `Pfam-A.hmm` to be `hmmpress`-indexed (companion `.h3f/.h3i/.h3m/.h3p` files). `debug_wsl()` now checks this specifically and `install_isoform_databases()` performs the indexing itself (activating the `isoform_tools` conda environment first, since `hmmer` lives inside it rather than on the bare WSL PATH).
* **Custom `CPAT_DATA` / `SIGNALP_DIR` paths not taking effect.** These are now persisted to a dedicated file (`~/.isoform_tools_env.sh`) that every predictor invocation sources, instead of `~/.bashrc` (which is never read by the non-interactive scripts this package runs).
* **`use_wsl` not defaulting the way you'd expect.** All `use_wsl` parameters across the package now default to `TRUE` on Windows and `FALSE` elsewhere; you only need to pass it explicitly to override that.

## Functional Analysis Methods

| Section            | Method                      | Type                |
|--------------------|-----------------------------|---------------------|
| GO BP/MF/CC        | `clusterProfiler::enrichGO` | ORA (Hypergeom.)    |
| Reactome           | `ReactomePA::enrichPathway` | ORA                 |
| Disease Ontology   | `DOSE::enrichDO`            | ORA (Human Only)    |
| Reactome & KEGG    | `gsePathway` / `gseKEGG`    | GSEA (Ranked lists) |
| MSigDB (Hallmark)  | `fgsea::fgseaMultilevel`    | GSEA (Native MH/H)  |
| SPIA               | `SPIA::spia`                | Topology Modeling   |
| EnrichR TF         | `enrichR::enrichr`          | Fisher's Exact      |

## Advanced Subgroup Analyses (Subsampling Data)

The pipeline includes a built-in `subset_sample` parameter. This allows you to evaluate any standard R logical condition as a string to instantly filter your sample table. Any samples that do not match your subsetting criteria are dropped dynamically along with their corresponding columns in the raw count matrix, completely bypassing the need to create temporary metadata files.

This is highly effective when you want to loop over distinct tissue types, cell cohorts, or time-points to look at treatments independently.

### Fictional Example: Target Analysis Across Brain Regions

Imagine you have a complex study examining a **Drug Treatment** (`Treated` vs `Vehicle`) across multiple distinct **Brain Regions** (`Cortex`, `Striatum`, `Hippocampus`, `Cerebellum`), and you want to execute 4 entirely separate differential expression analyses.

#### Sample Table (`brain_sample_table.csv`):
| sample_id | treatment | brain_region | batch |
|-----------|-----------|--------------|-------|
|ctx_v_1    | Vehicle   | Cortex       | B1    |
|ctx_t_1    | Treated   | Cortex       | B1    |
|str_v_1    | Vehicle   | Striatum     | B1    |
|str_t_1    | Treated   | Striatum     | B2    |
|hip_v_1    | Vehicle   | Hippocampus  | B2    |
|hip_t_1    | Treated   | Hippocampus  | B2    |

#### The Automated Loop Command:

By utilizing the `subset_sample` option within a standard R `for` loop, you can process every brain region independently in a single execution block:

```R
library(ExpressOM)

# Define the subgroups you want to analyze independently
target_regions <- c("Cortex", "Striatum", "Hippocampus", "Cerebellum")

# Run all independent sub-analyses sequentially
for (region in target_regions) {
  
  message(paste("--- Running Pipeline for Region:", region, "---"))
  
  run_bulk_pipeline(
    ensembl_package_name = "EnsDb.Hsapiens.v107",
    count_type           = "matrix",
    matrix_file          = "./data/brain_raw_counts.csv",
    sample_table         = "./data/brain_sample_table.csv",
    
    # Organizes outputs into region-specific target folders
    out_dir              = paste0("./results/analysis_", region), 
    
    # Models the treatment effect within that region specifically
    model                = "~ batch + treatment", 
    level                = "Treated", 
    base                 = "Vehicle",
    
    # Dynamically subsets the metadata columns per iteration
    subset_sample        = paste0("brain_region == '", region, "'")
  )
}
```

# Resulting Workspace Tree
```txt
The workflow isolates the data on the fly and segments your outputs cleanly:
results/
├── analysis_Cortex/
│   ├── DE_raw_results/  # Only Cortex Treated vs Vehicle metrics
│   ├── Plots/           # Cortex specific PCA, Volcanoplots, and Boxplots
│   └── ORA/             # Functional pathways enriched specifically in Cortex
├── analysis_Striatum/
│   └── ...
├── analysis_Hippocampus/
│   └── ...
└── analysis_Cerebellum/
    └── ...
```
### Excluding Outlier Samples on the Fly

If exploratory data analysis (EDA) reveals an outlier sample (e.g., due to low sequencing depth or batch contamination), you can blacklist it directly in the pipeline command using `remove_sample`. 

This safely omits the sample from the metadata matrix and cross-filters your counts automatically:

```R
library(ExpressOM)

run_bulk_pipeline(
  ensembl_package_name = "EnsDb.Mmusculus.v112",
  count_type           = "matrix",
  matrix_file          = "./data/raw_counts.csv",
  sample_table         = "./data/sample_table.tsv",
  out_dir              = "./results_cleaned",
  model                = "~ genotype",
  level                = "KO",
  base                 = "WT",
  
  # Clean up data without modifying files on disk:
  remove_sample        = c("Sample1", "Sample2") 
)
```

## Isoform-Level Analysis: External Predictors (CPAT, SignalP, Pfam)

When `run_isoform = TRUE` and `run_predictors = TRUE`, `run_isoform_switch()`
can enrich the `IsoformSwitchAnalyzeR` result with coding-potential (CPAT),
signal-peptide (SignalP), and protein-domain (Pfam via InterProScan or
`hmmscan`) predictions. **This now works the same way on Windows and on
Linux/macOS:**

* **Windows**: the tools run inside a WSL Ubuntu distribution. Set
  `use_wsl = TRUE` (the default on Windows) and `wsl_distro` to the name of
  your installed distribution.
* **Linux / macOS** (including an R session already running *inside* WSL):
  the tools run natively against `PATH` / a conda environment named
  `isoform_tools`. No `use_wsl` / `wsl_distro` configuration is needed.

```R
# Check tool/database availability first (works on any platform):
debug_wsl(use_wsl = FALSE)          # native Linux / macOS
debug_wsl(distro = "Ubuntu-22.04")  # Windows, routes through WSL

# Install tools + reference databases:
install_wsl_isoform_tools(distro = "Ubuntu-22.04")   # Windows / WSL
install_isoform_databases(use_wsl = FALSE)           # Linux / macOS (tools via mamba/conda separately)

switch_list <- run_isoform_switch(
  isoform_obj    = isoform_import,
  condition      = "condition",
  level          = "Treated",
  base           = "Control",
  fasta_file     = "reference/transcripts.fa",
  gff_file       = "reference/annotation.gtf",
  out_dir        = "results/isoform",
  run_predictors = TRUE,
  use_wsl        = TRUE,             # ignored on Linux/macOS; defaults to TRUE on Windows anyway
  wsl_distro     = "Ubuntu-22.04",   # ignored on Linux/macOS
  predictor_cpu  = NULL,             # NULL auto-detects available CPU cores
  log_dir        = NULL              # NULL -> out_dir/Log; override to point at a shared log tree
)
```

See `vignette("isoform-predictors", package = "ExpressOM")` for a full
walkthrough, including troubleshooting and performance notes.

> 💡 **Custom `fasta_file`/`gff_file` compatibility:** `run_isoform_switch()`
> normalizes transcript IDs (stripping Ensembl-style version suffixes,
> GENCODE-style `|`-delimited header fields, and space-separated
> descriptions) consistently across the count matrix, the FASTA, and the
> annotation before they're compared, and asks `importRdata()` for the same
> tolerance (`ignoreAfterBar`/`ignoreAfterSpace`/`ignoreAfterPeriod`). If you
> still see a "no transcripts match" / near-zero-overlap error with a custom
> reference, it usually means the FASTA/GTF and the quantification were
> genuinely built from different transcriptome versions rather than just a
> naming-convention mismatch.

Independently of `run_predictors`, `run_isoform_switch()` also runs
`IsoformSwitchAnalyzeR::analyzeAlternativeSplicing()` on every call (exon
skipping, mutually exclusive exons, intron retention, alternative 5'/3'
splice sites, alternative TSS/TES), and the report generator adds four
genome-wide summary figures alongside the existing dIF-vs-q-value overview
plot: `ConsequenceSummary.pdf`, `ConsequenceEnrichment.pdf`,
`SplicingSummary.pdf`, and `SplicingEnrichment.pdf` — matching the equivalent
figures in the `IsoformSwitchAnalyzeR` vignette. Each is best-effort: on a
small or single-direction dataset one or more may legitimately be skipped
(reported in the console and in `Log/Isoform/`), without affecting anything
else in the run.

## Isoform-Level Analysis: Enhanced Visualization (DEXSeq, Switch Plots, Sashimi, Exon Usage)

Beyond DTE/DTU tables and the transcript-proportion barplots, the isoform
pipeline can generate a richer set of graphics for interpreting *how* a
gene's isoform usage differs between conditions. These are produced
automatically as part of the DTE/DTU HTML report whenever `run_isoform = TRUE`,
and are also available as standalone functions for ad-hoc use.

| Plot | Function | What it shows |
|------|----------|----------------|
| **Isoform switch summary** | `plot_isoform_switch_summary()` | Composite figure per gene: isoform/transcript structure (with ORF, coding potential, protein domains, signal peptides if `run_predictors = TRUE`) alongside gene expression, isoform expression, and isoform usage (IF). Wraps `IsoformSwitchAnalyzeR::switchPlotTopSwitches()`/`switchPlot()`. |
| **DEXSeq transcript usage** | `plot_dexseq_gene()` | Classic DEXSeq-style plot comparing fitted expression (or usage, via `splicing = TRUE`) of every transcript in a gene across conditions, with significant transcripts highlighted. Requires `run_dexseq = TRUE`. |
| **Sashimi-style junction usage** | `plot_isoform_sashimi()` | Mirrored splice-junction diagram: exon structure on a shared genomic axis, with arcs above/below representing isoform-fraction-weighted junction usage in `level` vs `base`. |
| **Exon-bin usage comparison** | `plot_exon_usage_comparison()` | Bar chart comparing isoform-expression-weighted signal across non-overlapping exon bins (DEXSeq-style flattened annotation) between conditions — a coverage-style comparison built from transcript-level quantification. |

> 💡 **Note on the sashimi-style plot:** `plot_isoform_sashimi()` and
> `plot_exon_usage_comparison()` are built entirely from transcript-level
> quantification (the isoform-fraction values IsoformSwitchAnalyzeR already
> computes), **not** from aligned-read/junction coverage in a BAM file. They
> give a fast, alignment-free approximation of splicing differences; for
> base-pair-resolution coverage you would need a BAM-based tool (e.g.
> ggsashimi, ggbio, Gviz) run separately.

### Turning on the DEXSeq engine and enhanced plots

```R
switch_list <- run_isoform_switch(
  isoform_obj    = isoform_import,
  condition      = "condition",
  level          = "Treated",
  base           = "Control",
  fasta_file     = "reference/transcripts.fa",
  gff_file       = "reference/annotation.gtf",
  out_dir        = "results/isoform"
)

# Complementary DEXSeq-based DTU test (transcripts treated as DEXSeq exonic bins,
# following the Soneson/Love/Robinson "Swimming downstream" workflow)
dexseq_res <- run_dexseq_dtu(isoform_import, condition = "condition",
                              level = "Treated", base = "Control")

# Standalone plots for a gene of interest
plot_isoform_switch_summary(switch_list, plot_dir = "results/isoform/plots",
                             genes_of_interest = c("TP53", "BCL2"))
plot_isoform_sashimi(switch_list, gene = "TP53", level = "Treated", base = "Control",
                      plot_dir = "results/isoform/plots")
plot_exon_usage_comparison(switch_list, gene = "TP53", level = "Treated", base = "Control",
                            plot_dir = "results/isoform/plots")
plot_dexseq_gene(dexseq_res$dxr_list, gene_id = "ENSG00000141510",
                  plot_dir = "results/isoform/plots", gene_symbol = "TP53")
```

From the full pipeline, the same behaviour is reached with two extra arguments
to `run_bulk_pipeline()`:

```R
run_bulk_pipeline(
  ensembl_package_name = "EnsDb.Hsapiens.v107",
  count_type           = "salmon",
  data_dir             = "./data/counts",
  out_dir              = "./results",
  sample_table         = "./data/sample_table.csv",
  model                = "~ condition",
  level                = "Treated",
  base                 = "Control",
  run_isoform          = TRUE,
  run_dexseq           = TRUE,                       # enable the complementary DEXSeq DTU engine
  isoform_report_genes = c("TP53", "BCL2"),          # gets every plot in the table above
  isoform_plot_top_n   = 15                          # auto-plot the top 15 isoform switches genome-wide
)
```

### Where the files land

```txt
results/
└── IsoformSwitch/
    ├── DTU_DTE_report/
    │   ├── report.html                       # now includes DEXSeq / Switch Overview / gene sections
    │   ├── DTE_results_annotated.csv
    │   ├── DTU_results_annotated.csv
    │   ├── DEXSeq_results_annotated.csv       # only if run_dexseq = TRUE
    │   └── plots/
    │       ├── DTE_volcano.pdf, DTE_MA.pdf, DTU_pvalue_hist.pdf, ...
    │       ├── IsoformSwitch_overview.pdf     # dIF vs. switch q-value, genome-wide
    │       ├── ConsequenceSummary.pdf         # genome-wide functional consequence bar chart (NEW)
    │       ├── ConsequenceEnrichment.pdf      # gain-vs-loss enrichment per consequence (NEW)
    │       ├── SplicingSummary.pdf            # genome-wide alt. splicing event-type bar chart (NEW)
    │       ├── SplicingEnrichment.pdf         # gain-vs-loss enrichment per splice type (NEW)
    │       ├── top_switches/                  # auto-selected top isoform_plot_top_n switches
    │       │   └── *.pdf
    │       ├── proportions_<GENE>.pdf
    │       ├── SwitchPlot_<GENE>.pdf
    │       ├── DEXSeq_<GENE>.pdf              # only if run_dexseq = TRUE
    │       ├── Sashimi_<GENE>.pdf
    │       └── ExonUsage_<GENE>.pdf
    ├── plots/switch_plots_with_predictors/    # refreshed switch plots incl. domains/coding potential
    │   └── ...                                # (only if run_predictors = TRUE)
    ├── dexseq_results.rds                     # only if run_dexseq = TRUE
    └── switch_list.rds
```

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).