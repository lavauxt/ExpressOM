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

* **`DE_Results/`**: TSV tables of raw/filtered differential expression results and normalized counts.
* **`Plots/`**: PCA plots, Sample Correlation Heatmaps, MA plots, Volcano plots, and top DE gene boxplots.
* **`ORA/` & `GSEA/`**: Extensive targets and plots for GO mapping, Reactome, Disease Ontology, and KEGG generic pathways (Dotplots, Ridgeplots, Pathway Graphs).
* **`SPIA/`**: Signaling Pathway Impact Analysis graphs and Evidence CSVs.
* **`Transcription_Factors/`**: Output mapping active transcription factor perturbations based on dynamically mapped Enrichr libraries.
* **`RegionReport/`**: Auto-generated interactive `RegionReport` HTML documents.
* **`Save_rdata/`**: The complete populated R environment serialized as an `.RData` file.
* **`Log/`**: Final model constraints, environment data, session info, and processing warnings/errors log.

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
│   ├── DE_Results/      # Only Cortex Treated vs Vehicle metrics
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

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).