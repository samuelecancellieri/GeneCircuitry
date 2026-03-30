---
layout: default
title: Preprocessing & QC
nav_order: 2
parent: User Guide
---

# Preprocessing & QC

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

The preprocessing stage transforms raw single-cell count data into a clean, normalized, clustered `AnnData` object ready for GRN inference and module detection. All functions live in `genecircuitry/preprocessing.py` and use `genecircuitry/config.py` for default parameter values.

---

## Input requirements

Your input file must be a valid `AnnData` object (`.h5ad` or `.h5`):

- `adata.X` — raw or normalized count matrix (cells × genes)
- `adata.obs_names` — unique cell barcodes
- `adata.var_names` — gene symbols (required for CellOracle base GRN; must include `MT-` prefix for human mitochondrial genes)

If you supply raw counts, GeneCircuitry stores them in `adata.layers['raw_count']` before any normalization.

---

## Step 1 — Quality control (`perform_qc`)

Calculates per-cell QC metrics, filters low-quality cells and lowly expressed genes, and optionally generates violin/scatter plots.

```python
from genecircuitry.preprocessing import perform_qc

adata = perform_qc(
    adata,
    min_genes=None,          # int  — minimum genes per cell (config.QC_MIN_GENES = 100)
    min_counts=None,         # int  — minimum UMI counts per cell (config.QC_MIN_COUNTS = 200)
    max_counts=None,         # int  — maximum counts, guards against doublets (None = no limit)
    pct_counts_mt_max=None,  # float — max % mitochondrial reads (config.QC_PCT_MT_MAX = 20.0)
    min_cells=None,          # int  — minimum cells expressing a gene (config.QC_MIN_CELLS = 3)
    plot=True,               # bool — generate QC violin and scatter plots
    figsize=None,            # tuple — figure size (config.PLOT_FIGSIZE_LARGE)
    save_plots=None,         # str  — path prefix for saving plots
)
```

**What it adds to `adata`:**

| Location    | Key                 | Description                        |
| ----------- | ------------------- | ---------------------------------- |
| `adata.obs` | `n_genes_by_counts` | Genes detected per cell            |
| `adata.obs` | `total_counts`      | Total UMI counts per cell          |
| `adata.obs` | `pct_counts_mt`     | Percentage of mitochondrial counts |
| `adata.var` | `mt`                | Boolean flag: mitochondrial gene   |
| `adata.var` | `ribo`              | Boolean flag: ribosomal gene       |
| `adata.var` | `hb`                | Boolean flag: haemoglobin gene     |

**Config-driven defaults** — override them per-call or globally:

```python
# Per-call override (does not change global config)
adata = perform_qc(adata, min_genes=300, pct_counts_mt_max=15.0)

# Global override (affects all subsequent calls)
from genecircuitry import config
config.update_config(QC_MIN_GENES=300, QC_PCT_MT_MAX=15.0)
adata = perform_qc(adata)
```

### Recommended QC thresholds

| Metric              | Typical range | Notes                               |
| ------------------- | ------------- | ----------------------------------- |
| `min_genes`         | 200–500       | Fewer genes → likely empty droplets |
| `min_counts`        | 500–1 000     | Low counts → poor capture           |
| `max_counts`        | 20 000–50 000 | Very high counts → likely doublets  |
| `pct_counts_mt_max` | 10–20 %       | High % → dying or stressed cells    |

Inspect the generated QC plots before finalising thresholds — distributions vary by tissue and protocol.

### Mitochondrial gene detection

By default, mitochondrial genes are detected by the `MT-` prefix (human). For mouse data:

```python
import scanpy as sc
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True,
    var_type_prefixes={"mt": "mt-"}    # lowercase for mouse
)
```

---

## Step 2 — Normalization (`perform_normalization`)

Library-size normalizes counts, log1p-transforms, identifies highly variable genes (HVGs), and scales the data.

```python
from genecircuitry.preprocessing import perform_normalization

adata = perform_normalization(
    adata,
    target_sum=None,    # float — library-size target (config.NORMALIZE_TARGET_SUM = 1e4)
    n_top_genes=None,   # int   — number of HVGs to keep (config.HVGS_N_TOP_GENES = 2000)
)
```

**What happens internally:**

1. **Store raw counts** — `adata.layers['raw_count'] = adata.X.copy()`
2. **Library-size normalization** — `sc.pp.normalize_total(adata, target_sum=1e4)`
3. **Log1p transform** — `sc.pp.log1p(adata)`
4. **HVG selection** — `sc.pp.highly_variable_genes(adata, n_top_genes=2000)`
5. **Scale** — `sc.pp.scale(adata)`

{: .note }

> Raw counts are stored in `adata.layers['raw_count']` and are required by CellOracle. Do not skip normalization before running CellOracle.

---

## Step 3 — Dimensionality reduction & clustering

Runs PCA, builds the KNN graph, computes UMAP, and performs Leiden clustering.

```python
from genecircuitry.preprocessing import perform_dimensionality_reduction_clustering

adata = perform_dimensionality_reduction_clustering(
    adata,
    n_comps=None,       # int   — PCA components (config.PCA_N_COMPS = 50)
    n_neighbors=None,   # int   — KNN neighbors (config.NEIGHBORS_N_NEIGHBORS = 15)
    n_pcs=None,         # int   — PCs for neighbor graph (config.NEIGHBORS_N_PCS = 40)
    resolution=None,    # float — Leiden resolution (config.LEIDEN_RESOLUTION = 0.5)
)
```

**What it adds:**

| Location     | Key      | Description           |
| ------------ | -------- | --------------------- |
| `adata.obsm` | `X_pca`  | PCA embedding         |
| `adata.obsm` | `X_umap` | UMAP embedding        |
| `adata.obs`  | `leiden` | Leiden cluster labels |

**Tuning clustering resolution:**

- Lower resolution (e.g. 0.2) → fewer, larger clusters
- Higher resolution (e.g. 1.0) → more, smaller clusters

---

## GRN preprocessing (`perform_grn_pre_processing`)

This optional step prepares a processed `AnnData` specifically for CellOracle by downsampling cells, reselecting HVGs optimized for TF analysis, and computing diffusion maps.

```python
from genecircuitry.celloracle_processing import perform_grn_pre_processing

adata_grn = perform_grn_pre_processing(
    adata,
    cluster_key="leiden",   # required: adata.obs column with cluster labels
    cell_downsample=20000,  # max cells (downsample if larger)
    top_genes=None,         # HVGs to select (config.HVGS_N_TOP_GENES)
    n_neighbors=None,       # config.NEIGHBORS_N_NEIGHBORS
    n_pcs=None,             # config.NEIGHBORS_N_PCS
)
```

See [GRN Inference](../celloracle/) for the full CellOracle workflow that follows.

---

## Complete preprocessing example

```python
import scanpy as sc
from genecircuitry import config, set_random_seed
from genecircuitry.preprocessing import (
    perform_qc,
    perform_normalization,
    perform_dimensionality_reduction_clustering,
)

# 1. Reproducibility
set_random_seed(42)

# 2. Load data
adata = sc.read_h5ad("data/my_cells.h5ad")
print(f"Raw shape: {adata.shape}")  # (n_cells, n_genes)

# 3. QC filtering
adata = perform_qc(
    adata,
    min_genes=300,
    min_counts=1000,
    pct_counts_mt_max=15.0,
    plot=True,
    save_plots="results/figures/qc/qc_metrics",
)
print(f"After QC: {adata.shape}")

# 4. Normalization
adata = perform_normalization(adata)

# 5. Dimensionality reduction and clustering
adata = perform_dimensionality_reduction_clustering(
    adata,
    resolution=0.6,
)
print(f"Clusters: {adata.obs['leiden'].nunique()}")

# 6. Save
adata.write("results/clustered_adata.h5ad")
```

---

## Configuration reference

| Config parameter        | Default | Description                       |
| ----------------------- | ------- | --------------------------------- |
| `QC_MIN_GENES`          | `100`   | Minimum genes per cell            |
| `QC_MIN_COUNTS`         | `200`   | Minimum total counts per cell     |
| `QC_MAX_COUNTS`         | `None`  | Maximum counts per cell           |
| `QC_PCT_MT_MAX`         | `20.0`  | Max mitochondrial %               |
| `QC_MIN_CELLS`          | `3`     | Min cells expressing a gene       |
| `NORMALIZE_TARGET_SUM`  | `1e4`   | Library-size normalization target |
| `HVGS_N_TOP_GENES`      | `2000`  | Number of highly variable genes   |
| `PCA_N_COMPS`           | `50`    | PCA components                    |
| `NEIGHBORS_N_NEIGHBORS` | `15`    | KNN graph k                       |
| `NEIGHBORS_N_PCS`       | `40`    | PCs used for neighbor graph       |
| `LEIDEN_RESOLUTION`     | `0.5`   | Leiden clustering resolution      |

See [Configuration](../configuration/) for the full parameter table.
