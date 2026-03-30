---
layout: default
title: Plotting System
nav_order: 5
parent: User Guide
---

# Plotting System

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

All canonical plots are generated through the `genecircuitry/plotting/` subpackage. This replaces the legacy inline plotting code that was scattered across `preprocessing.py`, `grn_deep_analysis.py`, and `hotspot_processing.py`.

```
genecircuitry/plotting/
├── __init__.py          # Re-exports all public plot functions
├── qc_plots.py          # QC violin and scatter plots
├── grn_plots.py         # GRN network graphs, centrality heatmaps, rank plots
├── hotspot_plots.py     # Hotspot local correlation heatmaps, module plots
└── utils.py             # Shared helpers: save_plot(), plot_exists()
```

---

## Shared utilities (`plotting/utils.py`)

### `save_plot(fig, filepath, plot_type, metadata, skip_existing)`

Saves a matplotlib figure to disk, creates parent directories automatically, logs the save event, and returns `True` on success.

```python
from genecircuitry.plotting.utils import save_plot

saved = save_plot(
    fig=fig,
    filepath="results/figures/grn/my_network.png",
    plot_type="grn",          # label for logging
    metadata={"cluster": "0"},
    skip_existing=True,       # skip if file already exists
)
```

### `plot_exists(filepath, skip_existing)`

Returns `True` if the file already exists and `skip_existing=True`. Used internally in every plot function to avoid re-generating expensive plots.

```python
from genecircuitry.plotting.utils import plot_exists

if not plot_exists("results/figures/grn/my_plot.png", skip_existing=True):
    # generate and save
```

---

## QC plots (`plotting/qc_plots.py`)

### `plot_qc_violin_pre_filter(adata, save_name, figsize, skip_existing)`

Three-panel violin plot showing **genes per cell**, **total counts**, and **mitochondrial %** _before_ QC filtering.

```python
from genecircuitry.plotting.qc_plots import plot_qc_violin_pre_filter

plot_qc_violin_pre_filter(
    adata,
    save_name="run1",          # appended to filename: violin_pre_filter_run1.png
    figsize=None,              # config.PLOT_FIGSIZE_LARGE if None
    skip_existing=True,
)
# Saved to: config.FIGURES_DIR_QC/violin_pre_filter_run1.png
```

### `plot_qc_violin_post_filter(adata, save_name, figsize, skip_existing)`

Same three-panel violin plot _after_ QC filtering — for before/after comparison.

```python
from genecircuitry.plotting.qc_plots import plot_qc_violin_post_filter

plot_qc_violin_post_filter(adata, save_name="run1")
# Saved to: config.FIGURES_DIR_QC/violin_post_filter_run1.png
```

### `plot_qc_scatter_pre_filter(adata, save_name, figsize, skip_existing)`

Scatter plots: **counts vs genes**, **counts vs MT%**, **genes vs MT%** — useful for spotting doublets and dying cells.

```python
from genecircuitry.plotting.qc_plots import plot_qc_scatter_pre_filter

plot_qc_scatter_pre_filter(adata, save_name="run1")
# Saved to: config.FIGURES_DIR_QC/scatter_pre_filter_run1.png
```

---

## GRN plots (`plotting/grn_plots.py`)

### `generate_all_grn_plots(links, oracle, output_dir, stratification)`

Convenience function that generates the full suite of GRN plots for a given `Links` object:

- Network graphs (one per cluster)
- Centrality score heatmaps
- TF rank plots
- Cross-cluster score comparisons

```python
from genecircuitry.plotting.grn_plots import generate_all_grn_plots

generate_all_grn_plots(
    links=links,
    oracle=oracle,
    output_dir="results/figures/grn/",
    stratification="all_cells",    # label for filenames
)
```

### `plot_network_graph_single(graph, cluster, stratification, ...)`

Plots a single cluster's regulatory network as a directed graph using NetworkX.

```python
from genecircuitry.plotting.grn_plots import plot_network_graph_single
import networkx as nx

g = nx.DiGraph()  # your GRN graph
plot_network_graph_single(
    graph=g,
    cluster="0",
    stratification="TypeA",
    top_n_nodes=30,        # show top N hubs
    figsize=None,          # config.PLOT_FIGSIZE_SQUARED_LARGE
)
```

### `plot_enriched_tf_network(links, cluster, ...)`

Network graph where TF nodes are annotated with their top enriched pathway term (requires `gseapy`).

### `plot_tf_shared_target_network(links, cluster_a, cluster_b, ...)`

Bipartite network showing TFs and targets shared between two clusters — useful for comparing regulatory states.

---

## Hotspot plots (`plotting/hotspot_plots.py`)

### `plot_hotspot_local_correlations(hotspot_obj, skip_existing)`

Heatmap of pairwise local correlation z-scores between significant genes, with genes ordered by module.

```python
from genecircuitry.plotting.hotspot_plots import plot_hotspot_local_correlations

plot_hotspot_local_correlations(hs, skip_existing=True)
# Saved to: config.FIGURES_DIR_HOTSPOT/hotspot_local_correlations.png
```

### `plot_module_scores_violin(hotspot_obj, adata, cluster_key, figsize, skip_existing)`

Violin plots of per-module scores split by cell cluster — shows which clusters are enriched for each gene module.

```python
from genecircuitry.plotting.hotspot_plots import plot_module_scores_violin

plot_module_scores_violin(
    hotspot_obj=hs,
    adata=adata,
    cluster_key="leiden",
    skip_existing=True,
)
# Saved to: config.FIGURES_DIR_HOTSPOT/module_scores_violin.png
```

---

## Figure sizing and DPI

All plot functions accept a `figsize` parameter; when `None`, they fall back to a config preset:

| Config key                   | Default    | Use                                   |
| ---------------------------- | ---------- | ------------------------------------- |
| `PLOT_FIGSIZE_SMALL`         | `(6, 4)`   | Inline/thumbnail plots                |
| `PLOT_FIGSIZE_MEDIUM`        | `(10, 7)`  | Standard single plots                 |
| `PLOT_FIGSIZE_LARGE`         | `(20, 15)` | QC panels with many subplots          |
| `PLOT_FIGSIZE_WIDE`          | `(20, 8)`  | Landscape / wide figures              |
| `PLOT_FIGSIZE_SQUARED`       | `(6, 6)`   | Square (network graphs)               |
| `PLOT_FIGSIZE_SQUARED_LARGE` | `(10, 10)` | Large square (network graphs)         |
| `PLOT_DPI`                   | `200`      | Screen rendering DPI                  |
| `SAVE_DPI`                   | `600`      | File save DPI for publication quality |
| `PLOT_FORMAT`                | `'png'`    | Output format                         |

Override globally:

```python
from genecircuitry import config
config.update_config(SAVE_DPI=300, PLOT_FORMAT="svg")
```

---

## Skipping existing plots

All plot functions accept `skip_existing=True` (default). When a file already exists at the target path, the function returns `False` without regenerating the plot. Set `skip_existing=False` to force regeneration:

```python
plot_qc_violin_pre_filter(adata, skip_existing=False)  # always regenerate
```

This is integrated with the checkpoint system — completed plot steps are not re-run on pipeline resume.

---

## Adding new plots

New plot functions must live in the appropriate `genecircuitry/plotting/` module (not in `preprocessing.py` or other analysis modules). Follow this template:

```python
# In genecircuitry/plotting/qc_plots.py (or grn_plots.py / hotspot_plots.py)

def plot_my_new_figure(adata, save_name="default", figsize=None, skip_existing=True):
    from .. import config
    from .utils import save_plot, plot_exists

    filepath = f"{config.FIGURES_DIR_QC}/my_figure_{save_name}.png"
    if plot_exists(filepath, skip_existing):
        return False

    if figsize is None:
        figsize = config.PLOT_FIGSIZE_MEDIUM

    fig, ax = plt.subplots(figsize=figsize)
    # ... plotting code ...

    return save_plot(fig, filepath, plot_type="qc", metadata={"save_name": save_name})
```

Then export it from `genecircuitry/plotting/__init__.py`.
