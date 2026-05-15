# ExpressOM

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

# 1. Build the database versions you need 
# This downloads the GTF mapping, builds an SQLite DB, and compresses it into a local tarball
create_homemade_db(species = "human", release = "107")
create_homemade_db(species = "mouse", release = "114")

# 2. Install the newly generated database packages explicitly
#    (Matches strings dynamically so version collisions are avoided)
install_internal_db("Hsapiens")
install_internal_db("Mmusculus")
```

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
  ensembl_package_name = "EnsDb.Mmusculus.v114", # Ensure you use the exact database name
  count_type = "matrix",
  matrix_file = "./data/HELIOS/HELIOS_raw_counts.csv",
  sample_table = "./data/HELIOS/sample_table.tsv",
  out_dir = "./resultsHELIOS",
  model = "~ cell_type + genotype + cell_type:genotype",
  level = "KO",
  base = "WT",
  padj_cutoff = 0.05
)
```

### Example 2: Transcriptomic Input (`salmon`, `kallisto`, `rsem`)
Provide the base `data_dir` containing subfolders matching your `sample_id` column. `tximport` will automatically look for the expected file name (e.g. `quant.sf` for `salmon`, `abundance.tsv` for `kallisto`).

```R
run_bulk_pipeline(
  ensembl_package_name = "EnsDb.Hsapiens.v107",
  count_type = "salmon",            # Options: "salmon", "kallisto", "rsem", "stringtie", etc.
  data_dir = "./data/counts",  
  out_dir = "./results_classic",
  sample_table = "./data/sample_table.csv",
  model = "~ condition",            # Design formula
  level = "Treated",                # Active/Treatment group
  base = "Control",                 # Reference/Baseline group
  padj_cutoff = 0.05,
  test = "Wald"                     # Default statistical test
)
```

---

## Pipeline Outputs

The `run_bulk_pipeline` function generates a heavily organized output structure in your designated `out_dir`:

* **`DE_Results/`**: TSV tables of raw/filtered differential expression results and normalized counts.
* **`Plots/`**: PCA plots, Sample Correlation Heatmaps, MA plots, Volcano plots, and top DE gene boxplots.
* **`ORA/` & `GSEA/`**: Extensive targets and plots for GO mapping, Reactome, Disease Ontology, and KEGG generic pathways (Dotplots, Ridgeplots, Pathway Graphs).
* **`SPIA/`**: Signaling Pathway Impact Analysis graphs and Evidence CSVs.
* **`Transcription_Factors/`**: Output mapping active transcription factor perturbations based on dynamically mapped Enrichr libraries.
* **`Reports/`**: Auto-generated interactive `RegionReport` HTML documents and RMarkdown standard summary reports.
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