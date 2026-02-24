---
layout: default
title: Configuration
nav_order: 5
---

# Configuration Reference

{: .no_toc }

All pipeline parameters live in `genecircuitry/config.py`. Import with `from genecircuitry import config`.

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

## Quality Control

| Parameter       | Default | Description                            |
| --------------- | ------- | -------------------------------------- |
| `QC_MIN_GENES`  | `100`   | Minimum genes expressed per cell       |
| `QC_MIN_COUNTS` | `200`   | Minimum total UMI counts per cell      |
| `QC_MAX_COUNTS` | `None`  | Maximum total counts (None = no limit) |
| `QC_PCT_MT_MAX` | `20.0`  | Maximum % mitochondrial counts         |
| `QC_MIN_CELLS`  | `3`     | Minimum cells expressing a gene        |

---

## Preprocessing

| Parameter              | Default    | Description                                  |
| ---------------------- | ---------- | -------------------------------------------- |
| `NORMALIZE_TARGET_SUM` | `1e4`      | Library-size normalization target (CPM-like) |
| `HVGS_N_TOP_GENES`     | `2000`     | Number of highly variable genes              |
| `HVGS_MIN_MEAN`        | `0.0125`   | Min mean expression for HVG selection        |
| `HVGS_MAX_MEAN`        | `3`        | Max mean expression for HVG selection        |
| `HVGS_MIN_DISP`        | `0.5`      | Min dispersion for HVG selection             |
| `PCA_N_COMPS`          | `50`       | Number of PCA components                     |
| `PCA_SVDSOLVE`         | `'arpack'` | SVD solver for PCA                           |

---

## Clustering / Neighbors

| Parameter               | Default       | Description                       |
| ----------------------- | ------------- | --------------------------------- |
| `NEIGHBORS_N_NEIGHBORS` | `15`          | `k` for KNN graph                 |
| `NEIGHBORS_N_PCS`       | `40`          | PCs used for neighbor computation |
| `NEIGHBORS_METHOD`      | `'umap'`      | Neighbor computation method       |
| `NEIGHBORS_METRIC`      | `'euclidean'` | Distance metric                   |
| `LEIDEN_RESOLUTION`     | `0.5`         | Leiden clustering resolution      |
| `LOUVAIN_RESOLUTION`    | `0.5`         | Louvain clustering resolution     |

---

## UMAP

| Parameter           | Default | Description                        |
| ------------------- | ------- | ---------------------------------- |
| `UMAP_MIN_DIST`     | `0.5`   | Minimum distance in UMAP embedding |
| `UMAP_SPREAD`       | `1.0`   | Spread parameter                   |
| `UMAP_N_COMPONENTS` | `2`     | Number of UMAP dimensions          |

---

## CellOracle / GRN

| Parameter                  | Default | Description                       |
| -------------------------- | ------- | --------------------------------- |
| `GRN_N_JOBS`               | `8`     | Parallel jobs for GRN computation |
| `GRN_MIN_TARGETS`          | `5`     | Minimum target genes per TF       |
| `GRN_CONFIDENCE_THRESHOLD` | `0.5`   | Link confidence filter            |

---

## Hotspot

| Parameter                      | Default | Description                      |
| ------------------------------ | ------- | -------------------------------- |
| `HOTSPOT_N_JOBS`               | `8`     | Parallel jobs                    |
| `HOTSPOT_N_NEIGHBORS`          | `30`    | Neighbors for spatial graph      |
| `HOTSPOT_FDR_THRESHOLD`        | `0.05`  | FDR cutoff for significance      |
| `HOTSPOT_MIN_GENES_PER_MODULE` | `10`    | Minimum genes per module         |
| `HOTSPOT_CORE_ONLY`            | `True`  | Use only core genes per module   |
| `HOTSPOT_TOP_GENES`            | `500`   | Top autocorrelated genes to keep |

---

## Plotting

| Parameter                    | Default     | Description           |
| ---------------------------- | ----------- | --------------------- |
| `PLOT_DPI`                   | `200`       | Screen DPI            |
| `SAVE_DPI`                   | `600`       | File save DPI         |
| `PLOT_FORMAT`                | `'png'`     | Output format         |
| `PLOT_FIGSIZE_SMALL`         | `(6, 4)`    | Small figure          |
| `PLOT_FIGSIZE_MEDIUM`        | `(10, 7)`   | Medium figure         |
| `PLOT_FIGSIZE_LARGE`         | `(20, 15)`  | Large figure          |
| `PLOT_FIGSIZE_WIDE`          | `(20, 8)`   | Wide landscape figure |
| `PLOT_FIGSIZE_SQUARED`       | `(6, 6)`    | Square figure         |
| `PLOT_FIGSIZE_SQUARED_LARGE` | `(10, 10)`  | Large square figure   |
| `PLOT_COLOR_PALETTE`         | `'viridis'` | Default colormap      |

---

## File I/O

| Parameter             | Default                    | Description              |
| --------------------- | -------------------------- | ------------------------ |
| `OUTPUT_DIR`          | `'output'`                 | Default output directory |
| `FIGURES_DIR`         | `'output/figures'`         | Figures root             |
| `FIGURES_DIR_QC`      | `'output/figures/qc'`      | QC figures               |
| `FIGURES_DIR_GRN`     | `'output/figures/grn'`     | GRN figures              |
| `FIGURES_DIR_HOTSPOT` | `'output/figures/hotspot'` | Hotspot figures          |
| `LOG_PATH`            | `'output/logs'`            | Log directory            |

---

## Runtime update

```python
from genecircuitry import config

# Update one or more parameters globally
config.update_config(QC_MIN_GENES=300, LEIDEN_RESOLUTION=0.8)

# Inspect current config
config.print_config()

# Get config as dict
cfg = config.get_config()
```

---

## Adding new parameters

1. Add the constant to `genecircuitry/config.py` with a docstring.
2. Add it to the `get_config()` return dict in the same file.
3. Add a test in `tests/test_config.py`.
