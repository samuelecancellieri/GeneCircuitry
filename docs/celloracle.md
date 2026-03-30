---
layout: default
title: GRN Inference (CellOracle)
nav_order: 3
parent: User Guide
---

# GRN Inference with CellOracle

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

GeneCircuitry uses [CellOracle](https://celloracle.org/) to infer **per-cluster gene regulatory networks (GRNs)**: directed graphs where transcription factors (TFs) regulate target genes. The output is a `Links` object containing TF–target pairs, regression coefficients, and significance values for every cluster in your dataset.

Functions live in `genecircuitry/celloracle_processing.py`.

{: .note }

> CellOracle is an optional dependency. Install it with `pip install -e ".[grn]"` before using this module. The pipeline skips GRN inference gracefully when CellOracle is not installed.

---

## Prerequisites

Your `AnnData` object must have:

| Location                    | Key                                        | Required       |
| --------------------------- | ------------------------------------------ | -------------- |
| `adata.X`                   | Normalized, log1p-transformed counts       | ✅             |
| `adata.layers['raw_count']` | Raw counts before normalization            | ✅ recommended |
| `adata.obs`                 | Cluster labels (e.g. `'leiden'`)           | ✅             |
| `adata.obsm`                | UMAP or PCA embedding                      | ✅             |
| `adata.var_names`           | Human gene symbols (for base GRN matching) | ✅             |

Use [`perform_normalization()`](../preprocessing/#step-2--normalization-perform_normalization) to produce the correct `layers['raw_count']` automatically.

---

## Step 0 — GRN preprocessing (`perform_grn_pre_processing`)

Before creating the Oracle object, prepare a CellOracle-optimized `AnnData`: downsample large datasets, reselect HVGs, and compute diffusion maps.

```python
from genecircuitry.celloracle_processing import perform_grn_pre_processing

adata_grn = perform_grn_pre_processing(
    adata,
    cluster_key="leiden",       # adata.obs column with clusters
    cell_downsample=20000,      # max cells — downsamples if larger
    top_genes=None,             # HVGs to use (config.HVGS_N_TOP_GENES = 2000)
    n_neighbors=None,           # config.NEIGHBORS_N_NEIGHBORS = 15
    n_pcs=None,                 # config.NEIGHBORS_N_PCS = 40
)
```

---

## Step 1 — Create Oracle object (`create_oracle_object`)

Initialises a CellOracle `Oracle` with the human promoter base GRN and your expression data.

```python
from genecircuitry.celloracle_processing import create_oracle_object

oracle = create_oracle_object(
    adata_grn,
    cluster_column_name="leiden",    # cluster column name
    embedding_name="X_umap",         # embedding key in adata.obsm
    species="human",                 # "human" | "mouse"
    TG_to_TF_dictionary=None,        # path to ATAC PKL file (see ATAC Peaks doc), or None
    raw_count_layer="raw_count",     # layer with raw counts (None = use adata.X)
)
```

**What happens internally:**

1. Loads the human (or mouse) promoter base GRN from CellOracle
2. Imports expression data from `adata.layers['raw_count']` (or `adata.X`)
3. Adds TF information from the base GRN
4. Optionally overlays custom TF–target relationships from an ATAC-derived motif matrix

{: .note }

> If you have ATAC-seq data, pass the path to the PKL file produced by [`process_atac_peaks()`](../atac-peaks/) as `TG_to_TF_dictionary` to use chromatin-informed base GRNs instead of promoter-only models.

---

## Step 2 — PCA on Oracle (`run_PCA`)

Fits PCA on the Oracle expression matrix and automatically selects the optimal number of components using the elbow method (maximum 50).

```python
from genecircuitry.celloracle_processing import run_PCA

oracle, n_comps = run_PCA(oracle)
# n_comps: auto-selected number of components (int, capped at 50)
```

A variance-explained plot is saved to `config.FIGURES_DIR_GRN/pca_variance.png`.

---

## Step 3 — KNN imputation (`run_KNN`)

Smooths gene expression via KNN-based imputation to reduce technical noise before network inference.

```python
from genecircuitry.celloracle_processing import run_KNN

oracle = run_KNN(
    oracle,
    n_comps=n_comps,   # from run_PCA()
)
```

**How k is chosen:** automatically set to 2.5 % of total cells (`k = int(0.025 × n_cells)`). Uses balanced KNN with parallel processing controlled by `config.GRN_N_JOBS`.

---

## Step 4 — Link inference (`run_links`)

Infers per-cluster TF–target links using ridge regression (α = 10), filters by p-value, and computes network centrality scores.

```python
from genecircuitry.celloracle_processing import run_links

links = run_links(
    oracle,
    cluster_column_name="leiden",
    p_cutoff=0.001,     # p-value cutoff — lower = stricter (fewer but more significant links)
)
```

**Output:** a `co.Links` object containing:

- TF–target pairs and regression coefficients
- Per-link p-values and significance flags
- Network topology metrics (in-degree, out-degree, betweenness centrality)

The degree distribution plot is saved to `config.FIGURES_DIR_GRN/grn_degree_distribution.png`.

---

## Complete workflow example

```python
import scanpy as sc
from genecircuitry import config, set_random_seed
from genecircuitry.celloracle_processing import (
    perform_grn_pre_processing,
    create_oracle_object,
    run_PCA,
    run_KNN,
    run_links,
)

set_random_seed(42)

# Load clustered data
adata = sc.read_h5ad("results/clustered_adata.h5ad")

# GRN-specific preprocessing
adata_grn = perform_grn_pre_processing(adata, cluster_key="leiden")

# Create Oracle
oracle = create_oracle_object(
    adata_grn,
    cluster_column_name="leiden",
    embedding_name="X_umap",
    raw_count_layer="raw_count",
)

# PCA (auto-selects n_comps)
oracle, n_comps = run_PCA(oracle)

# KNN imputation
oracle = run_KNN(oracle, n_comps=n_comps)

# Infer GRN links
links = run_links(oracle, cluster_column_name="leiden", p_cutoff=0.001)

# Save results
oracle.to_hdf5(f"{config.OUTPUT_DIR}/celloracle/oracle_object.celloracle.oracle")
links.to_hdf5(f"{config.OUTPUT_DIR}/celloracle/grn_links.celloracle.links")

print("GRN inference complete.")
```

---

## Output files

| File                                         | Description                        |
| -------------------------------------------- | ---------------------------------- |
| `celloracle/oracle_object.celloracle.oracle` | Serialised Oracle object (HDF5)    |
| `celloracle/grn_links.celloracle.links`      | Serialised Links object (HDF5)     |
| `celloracle/grn_merged_scores.csv`           | Centrality scores for all clusters |
| `celloracle/grn_filtered_links.pkl`          | Filtered links as DataFrame        |
| `figures/grn/pca_variance.png`               | Variance explained by PCA          |
| `figures/grn/grn_degree_distribution.png`    | Node degree distribution           |

---

## Using ATAC-seq data

If you have matched ATAC-seq peaks, use them to build a chromatin-informed base GRN:

```python
from genecircuitry.atac_peaks_processing import process_atac_peaks

# Step 1: convert BED peaks → TF motif matrix
tf_info_path = process_atac_peaks(
    bed_path="data/atac_peaks.bed",
    species="human",
    output_dir="results/atac/",
)

# Step 2: pass the PKL path to create_oracle_object
oracle = create_oracle_object(
    adata_grn,
    cluster_column_name="leiden",
    embedding_name="X_umap",
    TG_to_TF_dictionary=tf_info_path,    # ← ATAC-derived motifs
)
```

See [ATAC Peaks Processing](../atac-peaks/) for full documentation.

---

## Configuration reference

| Config parameter           | Default              | Description                                         |
| -------------------------- | -------------------- | --------------------------------------------------- |
| `GRN_N_JOBS`               | `8`                  | Parallel jobs for KNN imputation and link inference |
| `GRN_MIN_TARGETS`          | `5`                  | Minimum target genes per TF (post-filter)           |
| `GRN_CONFIDENCE_THRESHOLD` | `0.5`                | Link confidence filter threshold                    |
| `FIGURES_DIR_GRN`          | `output/figures/grn` | Directory for GRN plots                             |
| `OUTPUT_DIR`               | `output`             | Root directory for all outputs                      |

---

## Troubleshooting

**`ImportError: celloracle`** — install with `pip install -e ".[grn]"` or skip with `--skip-celloracle`.

**`KeyError: 'raw_count'`** — your AnnData is missing the raw count layer. Either run `perform_normalization()` (which stores raw counts automatically) or set `raw_count_layer=None` to use `adata.X` directly.

**Memory errors during KNN** — reduce `config.GRN_N_JOBS` or use `perform_grn_pre_processing()` to downsample cells.

**No links inferred for a cluster** — the cluster may be too small. CellOracle requires a minimum number of cells per cluster for reliable regression. Check cluster sizes with `adata.obs['leiden'].value_counts()`.
