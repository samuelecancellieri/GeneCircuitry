---
layout: default
title: Pipeline Overview
nav_order: 1
parent: User Guide
---

# Pipeline Overview

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

GeneCircuitry is orchestrated by the `PipelineController` class (`genecircuitry/pipeline/controller.py`). It handles step execution, checkpointing, state management, logging, and parallel stratified runs from a single entry point.

---

## Running the pipeline

### Command-line (recommended)

```bash
# Full pipeline with example data (PBMC 3k auto-downloaded)
python run_complete_analysis.py

# Full pipeline with your own data
python run_complete_analysis.py \
    --input data/my_cells.h5ad \
    --output results/run1 \
    --name "My Experiment"

# Skip optional analyses
python run_complete_analysis.py --skip-celloracle    # no GRN inference
python run_complete_analysis.py --skip-hotspot        # no gene modules
python run_complete_analysis.py --skip-celloracle --skip-hotspot   # preprocessing only

# Run as a module (alternative entry point)
python -m genecircuitry.pipeline \
    --input data/my_cells.h5ad \
    --output results/ \
    --steps load preprocessing clustering
```

### Python API

```python
from genecircuitry.pipeline.controller import PipelineController
from argparse import Namespace
from datetime import datetime

args = Namespace(
    input="data/my_cells.h5ad",
    output="results/",
    n_jobs=8,
    seed=42,
    cluster_key_stratification=None,
    parallel=False,
    steps=None,          # None = all steps
)

controller = PipelineController(args, start_time=datetime.now())
controller.run_complete_pipeline()
```

---

## Available step names

Use these with `--steps` to run only selected stages:

| Step name        | Description                             |
| ---------------- | --------------------------------------- |
| `load`           | Read the input `.h5ad` / `.h5` file     |
| `preprocessing`  | QC filtering + normalization            |
| `stratification` | Split AnnData by cluster key (optional) |
| `clustering`     | PCA → UMAP → Leiden clustering          |
| `atac_peaks`     | ATAC-seq motif scanning (optional)      |
| `celloracle`     | CellOracle GRN inference                |
| `hotspot`        | Hotspot gene module detection           |
| `grn_analysis`   | NetworkX deep GRN analysis              |
| `report`         | HTML/PDF report generation              |
| `summary`        | Text summary of all outputs             |

```bash
# Example: only preprocessing and clustering
python -m genecircuitry.pipeline \
    --input data.h5ad --output results/ \
    --steps load preprocessing clustering

# Example: only the analysis steps (assumes preprocessing done)
python -m genecircuitry.pipeline \
    --input results/clustered_adata.h5ad --output results/ \
    --steps celloracle hotspot grn_analysis report
```

---

## Step dependencies

Steps must be run in dependency order:

```
load
  └── preprocessing
        └── stratification (optional — required for stratified runs)
              └── clustering
                    ├── atac_peaks (optional — feeds into celloracle)
                    ├── celloracle
                    └── hotspot
                          └── grn_analysis
                                └── report
                                      └── summary
```

If checkpoint files exist for earlier steps, they are automatically skipped — see [Checkpoints & Resume](../checkpoints/).

---

## All CLI flags

| Flag                               | Short | Description                            | Default                |
| ---------------------------------- | ----- | -------------------------------------- | ---------------------- |
| `--input PATH`                     | `-i`  | Input `.h5ad` or `.h5` file            | PBMC 3k example        |
| `--output DIR`                     | `-o`  | Output directory                       | `output/`              |
| `--name STR`                       |       | Run name (appears in report title)     |                        |
| `--steps STEP…`                    |       | Only run these steps                   | all steps              |
| `--skip-celloracle`                |       | Skip CellOracle GRN inference          |                        |
| `--skip-hotspot`                   |       | Skip Hotspot module detection          |                        |
| `--cluster-key-stratification KEY` |       | `adata.obs` column for stratified runs | disabled               |
| `--parallel`                       |       | Stratified runs in parallel            | sequential             |
| `--n-jobs N`                       |       | Number of parallel workers             | `8`                    |
| `--seed N`                         |       | Global random seed                     | `42`                   |
| `--min-genes N`                    |       | QC: minimum genes per cell             | `config.QC_MIN_GENES`  |
| `--min-counts N`                   |       | QC: minimum counts per cell            | `config.QC_MIN_COUNTS` |

---

## PipelineController API

### Step methods

Each pipeline stage has a dedicated `run_step_*()` method. Call them programmatically to build custom workflows:

```python
controller = PipelineController(args, datetime.now())

controller.run_step_load()                       # → sets controller.adata
controller.run_step_preprocessing()              # → sets controller.adata_preprocessed
controller.run_step_stratification()             # → sets controller.adata_list
controller.run_step_clustering(adata)            # → returns AnnData with embeddings
controller.run_step_celloracle(adata)            # → runs CellOracle, saves results
controller.run_step_hotspot(adata)               # → runs Hotspot, saves results
controller.run_step_grn_analysis(csv_path)       # → NetworkX visualizations
```

### Execution modes

```python
# Run the complete pipeline (respects checkpoints)
controller.run_complete_pipeline(steps=None, parallel=False)

# Run stratified pipeline sequentially (one cluster at a time)
controller.run_stratified_pipeline_sequential()

# Run stratified pipeline in parallel
controller.run_stratified_pipeline_parallel(n_jobs=4)

# Process a single stratification (used internally)
controller.process_single_stratification(adata, cluster_name="TypeA")
```

### Controller state

After steps run, results are accessible as attributes:

```python
controller.adata                  # loaded AnnData
controller.adata_preprocessed     # after QC + normalization
controller.adata_list             # list of per-cluster AnnData (after stratification)
controller.adata_stratification_list  # corresponding cluster names
```

---

## Adding a new step

To add a custom step to the pipeline:

```python
# In genecircuitry/pipeline/controller.py — add to PipelineController class
def run_step_my_analysis(self, adata, log_dir=None):
    """Execute my analysis step."""
    from genecircuitry.pipeline import log_step, log_error
    log_step("Controller.MyAnalysis", "STARTED", {"n_cells": adata.n_obs})
    try:
        result = my_analysis_function(adata)
        log_step("Controller.MyAnalysis", "COMPLETED")
        return result
    except Exception as e:
        log_error("Controller.MyAnalysis", e)
        raise
```

Then add `"my_analysis"` to the step routing dict in `run_complete_pipeline()`.

---

## Output structure

```
results/
├── preprocessed_adata.h5ad        # QC + normalized AnnData
├── clustered_adata.h5ad           # After dimensionality reduction & clustering
├── report.html                    # Interactive HTML report
├── report.pdf                     # PDF report (requires weasyprint)
├── analysis_summary.txt
├── logs/
│   ├── pipeline.log               # All steps with timestamps and metrics
│   ├── error.log                  # Errors with full Python tracebacks
│   └── *.checkpoint               # Auto-resume markers (one per step)
├── figures/
│   ├── qc/                        # QC violin and scatter plots
│   ├── grn/                       # GRN network plots, rank plots, heatmaps
│   └── hotspot/                   # Module heatmaps, enrichment plots
├── celloracle/
│   ├── grn_merged_scores.csv
│   └── grn_filtered_links.pkl
├── hotspot/
│   ├── autocorrelation_results.csv
│   └── gene_modules.csv
└── stratified_analysis/
    └── <ClusterName>/             # One directory per stratification
        ├── report.html
        ├── figures/
        ├── celloracle/
        └── hotspot/
```

---

## Parallel execution

When using `--parallel`, each stratification (cell type) is processed by a separate worker process. Workers are spawned after the shared preprocessing/clustering stage completes.

```bash
python -m genecircuitry.pipeline \
    --input data.h5ad \
    --output results/ \
    --cluster-key-stratification cell_type \
    --parallel \
    --n-jobs 4
```

| Aspect   | Parallel                     | Sequential                     |
| -------- | ---------------------------- | ------------------------------ |
| Speed    | Faster (near-linear speedup) | Slower                         |
| Memory   | High (`n_jobs × data size`)  | Low                            |
| Use when | Many cell types, ample RAM   | Few types or limited resources |

{: .warning }

> Each parallel worker loads a full copy of the data. For large datasets, reduce `--n-jobs` or run sequentially.

---

## Troubleshooting

**Out-of-memory errors with `--parallel`**: reduce `--n-jobs 2` or drop the `--parallel` flag.

**Step fails but checkpoint exists**: delete the checkpoint file and re-run:

```bash
rm results/logs/clustering.checkpoint
python -m genecircuitry.pipeline --input data.h5ad --output results/
```

**`ModuleNotFoundError: celloracle`**: install the optional dep or add `--skip-celloracle`.

**`ModuleNotFoundError: hotspot`**: install the optional dep or add `--skip-hotspot`.
