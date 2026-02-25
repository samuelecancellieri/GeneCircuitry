---
layout: default
title: Architecture
nav_order: 4
---

# Architecture

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

## Key components

| File / Directory                         | Purpose                                                                       |
| ---------------------------------------- | ----------------------------------------------------------------------------- |
| `genecircuitry/config.py`                | **Single source of truth** — all parameters live here                         |
| `genecircuitry/pipeline/controller.py`   | `PipelineController` class — central orchestrator (~1900 lines)               |
| `run_complete_analysis.py`               | CLI entry point — thin wrapper around `genecircuitry.pipeline.main`           |
| `genecircuitry/preprocessing.py`         | Scanpy wrappers: QC, normalization, dimensionality reduction, clustering      |
| `genecircuitry/celloracle_processing.py` | CellOracle GRN inference (optional dep)                                       |
| `genecircuitry/hotspot_processing.py`    | Hotspot gene module analysis (optional dep)                                   |
| `genecircuitry/grn_deep_analysis.py`     | NetworkX network visualization                                                |
| `genecircuitry/plotting/`                | Canonical plot generation (`qc_plots.py`, `grn_plots.py`, `hotspot_plots.py`) |
| `genecircuitry/reporting/`               | HTML / PDF report generation                                                  |
| `genecircuitry/atac_peaks_processing.py` | ATAC-seq peak motif scanning (optional)                                       |
| `genecircuitry/enrichment_analysis.py`   | Gene set enrichment via gseapy (optional)                                     |

---

## Data flow

```
AnnData (.h5ad)
    │
    ▼
[load_data()]
    │
    ▼
[preprocessing_pipeline()]         genecircuitry/preprocessing.py
  ├─ perform_qc()
  └─ perform_normalization()
    │
    ▼ (checkpoint: preprocessed_adata.h5ad)
    │
    ├─► [stratification_pipeline()]  — optional, splits by cell type
    │
    ▼
[dimensionality_reduction_clustering()]
  └─ perform_dimensionality_reduction_clustering()
    │
    ▼ (checkpoint: clustered_adata.h5ad)
    │
    ├──────────────────────────┐
    ▼                          ▼
[celloracle_pipeline()]    [hotspot_pipeline()]
  ├─ perform_grn_pre_processing   create_hotspot_object()
  ├─ create_oracle_object()       run_hotspot_analysis()
  ├─ run_PCA()
  ├─ run_KNN()
  └─ run_links()
    │                          │
    ▼ (checkpoint: oracle + links)
    │
    ▼
[grn_deep_analysis_pipeline()]   genecircuitry/grn_deep_analysis.py
  └─ process, plot, compare networks
    │
    ▼
[generate_report()]              genecircuitry/reporting/
  └─ HTML + PDF with embedded figures
```

---

## PipelineController

`PipelineController` (in `genecircuitry/pipeline/controller.py`) is the central class. It:

- Holds pipeline state (`self.adata`, `self.adata_preprocessed`, `self.stratification_results`)
- Exposes one `run_step_*()` method per pipeline step
- Handles checkpointing, logging, and error recovery per step
- Supports stratified (per-cluster) runs via `run_stratified_pipeline_sequential()`

**Adding a new step:**

```python
# In genecircuitry/pipeline/controller.py — add to PipelineController class
def run_step_my_analysis(self, adata, log_dir=None):
    """Execute my custom analysis step."""
    log_step("Controller.MyAnalysis", "STARTED")
    try:
        result = my_analysis_function(adata)
        log_step("Controller.MyAnalysis", "COMPLETED")
        return result
    except Exception as e:
        log_error("Controller.MyAnalysis", e)
        raise
```

Then register the step name in `run_complete_pipeline()`.

---

## Config as single source of truth

All numeric thresholds, figure sizes, and directory paths are constants in `genecircuitry/config.py`.  
Functions use them as defaults, allowing per-call overrides without mutating global state:

```python
def perform_qc(adata, min_genes=None, ...):
    if min_genes is None:
        min_genes = config.QC_MIN_GENES   # falls back to config
```

To change a constant globally at runtime: `config.update_config(QC_MIN_GENES=500)`.

---

## Optional dependencies

CellOracle and Hotspot are wrapped in `try/except ImportError` in `genecircuitry/__init__.py`:

```python
try:
    from . import celloracle_processing
except ImportError:
    celloracle_processing = None
```

Downstream code checks `if genecircuitry.celloracle_processing is not None` before use. New optional-dep modules must follow the same pattern.

---

## Plotting subpackage

The canonical location for all plots is `genecircuitry/plotting/`:

| Module             | Contents                                                                                    |
| ------------------ | ------------------------------------------------------------------------------------------- |
| `qc_plots.py`      | `plot_qc_violin_pre_filter()`, `plot_qc_scatter_pre_filter()`, post-filter variants         |
| `grn_plots.py`     | `generate_all_grn_plots()`, `plot_enriched_tf_network()`, `plot_tf_shared_target_network()` |
| `hotspot_plots.py` | Hotspot module heatmaps and per-module gene plots                                           |
| `utils.py`         | `plot_exists()`, `save_plot()` — shared helpers                                             |

Legacy inline plotting in `grn_deep_analysis.py` and `hotspot_processing.py` is kept for backward compatibility but delegates to the `plotting/` subpackage.

---

## Logging

Two log files are created in `<output>/logs/`:

| File           | Contents                                                 |
| -------------- | -------------------------------------------------------- |
| `pipeline.log` | Every step with start/complete timestamp and key metrics |
| `error.log`    | Errors with full Python tracebacks                       |

```python
from genecircuitry.pipeline import log_step, log_error

log_step("MyStep", "STARTED", {"n_cells": adata.n_obs})
log_step("MyStep", "COMPLETED", {"result_count": len(results)})
log_error("MyStep", exception)
```
