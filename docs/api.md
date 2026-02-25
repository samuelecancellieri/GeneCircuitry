---
layout: default
title: API Reference
nav_order: 6
has_children: false
---

# API Reference

{: .no_toc }

Key public functions by module. All parameters with `None` defaults fall back to the corresponding `config.*` value.

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

## `genecircuitry.preprocessing`

### `perform_qc(adata, ...)`

Perform quality control: calculate QC metrics, filter cells/genes, and optionally generate plots.

```python
from genecircuitry.preprocessing import perform_qc

adata = perform_qc(
    adata,
    min_genes=None,          # int  — config.QC_MIN_GENES (100)
    min_counts=None,         # int  — config.QC_MIN_COUNTS (200)
    max_counts=None,         # int  — config.QC_MAX_COUNTS (None)
    pct_counts_mt_max=None,  # float — config.QC_PCT_MT_MAX (20.0)
    min_cells=None,          # int  — config.QC_MIN_CELLS (3)
    plot=True,               # bool — generate violin/scatter plots
    figsize=None,            # tuple — config.PLOT_FIGSIZE_MEDIUM
    save_plots=None,         # str  — path prefix for saving plots
)
```

**Adds to `adata`:** `obs['n_genes_by_counts']`, `obs['total_counts']`, `obs['pct_counts_mt']`, `var['mt']`, `var['ribo']`, `var['hb']`.

---

### `perform_normalization(adata, ...)`

Normalize library size, log1p-transform, select HVGs, scale.

```python
from genecircuitry.preprocessing import perform_normalization

adata = perform_normalization(
    adata,
    target_sum=None,     # float — config.NORMALIZE_TARGET_SUM (1e4)
    n_top_genes=None,    # int   — config.HVGS_N_TOP_GENES (2000)
)
```

Stores raw counts in `adata.layers['raw_counts']` before normalization.

---

### `perform_dimensionality_reduction_clustering(adata, ...)`

Run PCA → neighbors → UMAP → Leiden clustering.

```python
from genecircuitry.preprocessing import perform_dimensionality_reduction_clustering

adata = perform_dimensionality_reduction_clustering(
    adata,
    n_comps=None,       # int — config.PCA_N_COMPS (50)
    n_neighbors=None,   # int — config.NEIGHBORS_N_NEIGHBORS (15)
    n_pcs=None,         # int — config.NEIGHBORS_N_PCS (40)
    resolution=None,    # float — config.LEIDEN_RESOLUTION (0.5)
)
```

Adds `adata.obsm['X_pca']`, `adata.obsm['X_umap']`, `adata.obs['leiden']`.

---

## `genecircuitry.celloracle_processing`

{: .note }

> Requires `pip install -e ".[grn]"`. Returns `None` at module level if CellOracle is not installed.

### `perform_grn_pre_processing(adata, cluster_key, ...)`

Downsample, reselect HVGs, compute diffusion maps — prepares AnnData for CellOracle.

```python
from genecircuitry.celloracle_processing import perform_grn_pre_processing

adata_grn = perform_grn_pre_processing(
    adata,
    cluster_key="leiden",   # required
    cell_downsample=20000,
    top_genes=None,         # config.HVGS_N_TOP_GENES
    n_neighbors=None,       # config.NEIGHBORS_N_NEIGHBORS
    n_pcs=None,             # config.NEIGHBORS_N_PCS
)
```

### `create_oracle_object(adata, cluster_column_name, embedding_name, ...)`

Create a CellOracle `Oracle` object with base GRN loaded.

```python
oracle = create_oracle_object(
    adata,
    cluster_column_name="leiden",
    embedding_name="X_umap",
    species="human",             # "human" | "mouse"
    TG_to_TF_dictionary=None,    # path to ATAC PKL, or None
    raw_count_layer=None,        # layer name for raw counts
)
```

### `run_PCA(oracle)` → `(oracle, n_comps)`

Perform PCA and auto-select number of components (elbow method, max 50).

### `run_KNN(oracle, n_comps)` → `oracle`

KNN imputation with `k = 2.5% × n_cells`.

### `run_links(oracle, cluster_column_name, p_cutoff=0.001)` → `links`

Infer per-cluster GRN links using ridge regression; filter by p-value; compute network scores.

---

## `genecircuitry.hotspot_processing`

{: .note }

> Requires `pip install -e ".[hotspot]"`.

### `create_hotspot_object(adata, ...)`

Create a `Hotspot` object from AnnData.

```python
from genecircuitry.hotspot_processing import create_hotspot_object

hs = create_hotspot_object(
    adata,
    top_genes=None,         # config.HOTSPOT_TOP_GENES (500)
    layer_key="raw_count",
    model="danb",
    embedding_key="X_pca",
    normalization_key="n_counts",
)
```

### `run_hotspot_analysis(hotspot_obj, adata, cluster_key)` → `hotspot_obj`

Run autocorrelation test, create modules, compute local correlations.

---

## `genecircuitry.reporting`

### `generate_report(output_dir, ...)`

Generate HTML and/or PDF reports.

```python
from genecircuitry.reporting import generate_report

outputs = generate_report(
    output_dir="results/",
    title="My Analysis",
    subtitle="",
    adata=adata,
    celloracle_result=(oracle, links),  # or None
    hotspot_result=hs,                  # or None
    log_file="results/logs/pipeline.log",
    embed_images=True,    # False = smaller HTML via relative paths
    formats=["html", "pdf"],
)
# outputs = {"html": "results/report.html", "pdf": "results/report.pdf"}
```

---

## `genecircuitry.pipeline`

### `log_step(step_name, status, details=None)`

```python
from genecircuitry.pipeline import log_step, log_error

log_step("MyStep", "STARTED", {"n_cells": adata.n_obs})
log_error("MyStep", exception)
```

### `PipelineController(args, start_time)`

Central orchestrator class. See [Architecture](../architecture/) for the full step method list.

---

## `genecircuitry.config`

```python
from genecircuitry import config

config.update_config(QC_MIN_GENES=300)  # runtime update
config.print_config()                   # print all parameters
cfg = config.get_config()               # returns dict
```
