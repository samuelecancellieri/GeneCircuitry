---
layout: default
title: Quick Start
nav_order: 2
parent: Getting Started
---

# Quick Start

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

## Run the complete pipeline

The simplest way to run GeneCircuitry is through the entry-point wrapper:

```bash
# With example data (PBMC 3k auto-downloaded)
python run_complete_analysis.py

# With your own data
python run_complete_analysis.py \
    --input data/my_cells.h5ad \
    --output results/my_run \
    --name "My Experiment"

# Skip optional analyses
python run_complete_analysis.py --skip-celloracle  # no GRN inference
python run_complete_analysis.py --skip-hotspot      # no module identification

# See all options
python run_complete_analysis.py --help
```

---

## Run specific pipeline steps

```bash
# Run only QC + normalization + clustering (skip GRN and Hotspot)
python -m genecircuitry.pipeline \
    --input data/my_cells.h5ad \
    --output results/ \
    --steps load preprocessing clustering

# Resume from checkpoints — already-completed steps are skipped automatically
python -m genecircuitry.pipeline \
    --input data/my_cells.h5ad \
    --output results/ \
    --steps celloracle hotspot
```

**Available step names:** `load`, `preprocessing`, `stratification`, `clustering`, `atac_peaks`, `celloracle`, `hotspot`, `grn_analysis`, `report`, `summary`

---

## Stratified (per-cluster) analysis

Analyse each cell type independently in parallel:

```bash
python -m genecircuitry.pipeline \
    --input data/my_cells.h5ad \
    --output results/ \
    --cluster-key-stratification cell_type \
    --parallel \
    --n-jobs 4
```

Each stratification gets its own subdirectory under `results/stratified_analysis/<ClusterName>/`.

---

## Python API

```python
import scanpy as sc
from genecircuitry import config, set_random_seed, set_scanpy_settings
from genecircuitry.preprocessing import perform_qc, perform_normalization

# 1. Setup
set_random_seed(42)
set_scanpy_settings()

# 2. Load data
adata = sc.read_h5ad("data/my_cells.h5ad")

# 3. QC (all thresholds come from config)
adata = perform_qc(adata)

# Override thresholds per-call (doesn't alter global config)
adata = perform_qc(adata, min_genes=300, pct_counts_mt_max=15.0)

# 4. Normalization
adata = perform_normalization(adata)

# 5. GRN preprocessing + CellOracle (requires celloracle installed)
from genecircuitry.celloracle_processing import (
    perform_grn_pre_processing,
    create_oracle_object,
    run_PCA, run_KNN, run_links,
)
adata_grn = perform_grn_pre_processing(adata, cluster_key="leiden")
oracle = create_oracle_object(adata_grn, cluster_column_name="leiden",
                              embedding_name="X_umap")
oracle, n_comps = run_PCA(oracle)
oracle = run_KNN(oracle, n_comps=n_comps)
links = run_links(oracle, cluster_column_name="leiden")
```

---

## Expected output structure

```
results/
├── preprocessed_adata.h5ad        # QC + normalized AnnData
├── clustered_adata.h5ad           # After dim-reduction & clustering
├── report.html                    # Interactive analysis report
├── report.pdf                     # PDF report (requires weasyprint)
├── analysis_summary.txt
├── logs/
│   ├── pipeline.log               # All steps with timestamps
│   ├── error.log                  # Errors with full tracebacks
│   └── *.checkpoint               # Auto-resume markers
├── figures/
│   ├── qc/                        # QC violin/scatter plots
│   ├── grn/                       # GRN network plots, rank plots
│   └── hotspot/                   # Module heatmaps
├── celloracle/
│   ├── grn_merged_scores.csv
│   └── grn_filtered_links.pkl
├── hotspot/
│   ├── autocorrelation_results.csv
│   └── gene_modules.csv
└── stratified_analysis/
    └── <ClusterName>/             # One folder per stratification
        ├── report.html
        ├── figures/
        ├── celloracle/
        └── hotspot/
```

---

## Checkpointing

GeneCircuitry writes `.checkpoint` files to `logs/` after each step. If the pipeline is interrupted, re-running the same command resumes from where it left off. The checkpoint is invalidated if the input file or key parameters change.
