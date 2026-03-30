---
layout: default
title: Stratified Analysis
nav_order: 7
parent: User Guide
---

# Stratified Analysis

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

**Stratified analysis** runs the GeneCircuitry downstream pipeline (CellOracle + Hotspot + GRN analysis + Reporting) independently for each group in a categorical `adata.obs` column — typically a cell type or condition label. This allows you to compare GRN structure across cell populations.

Each stratification gets its own output directory, checkpoint files, and HTML report. Stratifications can run **sequentially or in parallel**.

---

## When to use stratified analysis

Use stratified analysis when you want to:

- Infer separate GRNs for each cell type (e.g. endothelial vs epithelial)
- Compare Hotspot gene modules across conditions (treated vs untreated)
- Identify cell-type-specific transcription factor usage
- Run reproducible per-cluster analyses as part of an atlas project

---

## CLI usage

```bash
# Sequential stratified analysis (one cluster at a time)
python -m genecircuitry.pipeline \
    --input data.h5ad \
    --output results/ \
    --cluster-key-stratification cell_type

# Parallel stratified analysis (all clusters at once)
python -m genecircuitry.pipeline \
    --input data.h5ad \
    --output results/ \
    --cluster-key-stratification cell_type \
    --parallel \
    --n-jobs 4
```

The `--cluster-key-stratification` value must be a column in `adata.obs`. Use `adata.obs.columns.tolist()` to list all available columns.

---

## Python API

```python
from genecircuitry.pipeline.controller import PipelineController
from argparse import Namespace
from datetime import datetime

args = Namespace(
    input="data.h5ad",
    output="results/",
    cluster_key_stratification="cell_type",  # obs column
    parallel=True,
    n_jobs=4,
    seed=42,
)

controller = PipelineController(args, start_time=datetime.now())

# Run preprocessing and stratification first
controller.run_step_load()
controller.run_step_preprocessing()
controller.run_step_stratification()

print(f"Stratifications: {controller.adata_stratification_list}")

# Sequential
controller.run_stratified_pipeline_sequential()

# Or parallel
controller.run_stratified_pipeline_parallel(n_jobs=4)
```

### Inspecting stratification results

```python
# After run_step_stratification(), the controller holds:
controller.adata_list                   # list of per-cluster AnnData objects
controller.adata_stratification_list    # matching cluster names

# Process only specific clusters
for adata_sub, name in zip(controller.adata_list, controller.adata_stratification_list):
    if name in ["T_cells", "NK_cells"]:
        controller.process_single_stratification(adata_sub, name)
```

---

## Output structure

Each cluster gets its own subdirectory under `<output>/stratified_analysis/`:

```
results/
└── stratified_analysis/
    ├── T_cells/
    │   ├── report.html
    │   ├── clustered_adata.h5ad
    │   ├── logs/
    │   │   ├── pipeline.log
    │   │   ├── error.log
    │   │   └── *.checkpoint
    │   ├── figures/
    │   │   ├── grn/
    │   │   └── hotspot/
    │   ├── celloracle/
    │   │   ├── grn_merged_scores.csv
    │   │   └── grn_filtered_links.pkl
    │   └── hotspot/
    │       ├── autocorrelation_results.csv
    │       └── gene_modules.csv
    ├── NK_cells/
    │   └── ...
    └── B_cells/
        └── ...
```

Checkpoint files are **per-stratification** — if one cluster fails, re-running the pipeline resumes from where it left off for that cluster only.

---

## Parallel vs sequential

| Aspect                | Parallel (`--parallel`)                | Sequential (default)                        |
| --------------------- | -------------------------------------- | ------------------------------------------- |
| **Speed**             | Near-linear speedup with `--n-jobs`    | One cluster finishes before the next starts |
| **Memory**            | `n_jobs × per-cluster data size`       | `1 × per-cluster data size`                 |
| **Failure isolation** | A failed cluster does not block others | Failure in one cluster stops the sequence   |
| **Log output**        | Interleaved (prefixed by cluster name) | Sequential, easy to follow                  |
| **Best for**          | Many clusters, ample RAM               | Few clusters or memory-constrained machines |

{: .warning }

> Each parallel worker receives a full copy of the AnnData for its cluster. For large datasets, start with `--n-jobs 2` and monitor memory usage before scaling up.

---

## Configuration

Stratified analysis respects all global config parameters. Per-cluster overrides can be applied before calling `process_single_stratification()`:

```python
from genecircuitry import config

# Apply stricter QC for a small cluster
config.update_config(QC_MIN_CELLS=5)
controller.process_single_stratification(adata_sub, "Rare_population")

# Reset to default
config.update_config(QC_MIN_CELLS=3)
```

---

## Troubleshooting

**Out-of-memory with `--parallel`**: reduce `--n-jobs` or run without `--parallel`.

**Some clusters produce no GRN links**: the cluster may be too small for CellOracle regression. Check cluster sizes; cells under ~100 are unlikely to yield reliable networks.

**Resume after partial failure**: delete the failed cluster's checkpoint files and re-run. Successful clusters are skipped:

```bash
rm -rf results/stratified_analysis/Rare_population/logs/*.checkpoint
python -m genecircuitry.pipeline --input data.h5ad --output results/ \
    --cluster-key-stratification cell_type
```
