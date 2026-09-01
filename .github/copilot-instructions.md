# GeneCircuitry AI Coding Instructions

**GeneCircuitry** is a Python package for transcriptional regulatory network (TRN) analysis using single-cell data. It integrates Scanpy, CellOracle, and Hotspot into a modular, checkpoint-enabled pipeline.

## Architecture

### Data Flow

`AnnData (.h5ad)` → **Preprocessing** (QC/normalization) → **Clustering** → **CellOracle** (GRN inference) + **Hotspot** (gene modules) → **Deep Analysis** (visualization) → **Reporting** (HTML/PDF)

### Key Components

| File                                     | Purpose                                                         |
| ---------------------------------------- | --------------------------------------------------------------- |
| `genecircuitry/config.py`                | **SINGLE SOURCE OF TRUTH** - all parameters                     |
| `genecircuitry/pipeline/controller.py`   | `PipelineController` class - central orchestrator (~1900 lines) |
| `run_complete_analysis.py`               | Entry point wrapper (calls `genecircuitry.pipeline.main`)       |
| `genecircuitry/preprocessing.py`         | Scanpy wrappers (QC, normalize, cluster)                        |
| `genecircuitry/celloracle_processing.py` | GRN inference with CellOracle (optional dep)                    |
| `genecircuitry/hotspot_processing.py`    | Spatial autocorrelation modules (optional dep)                  |
| `genecircuitry/grn_deep_analysis.py`     | Network visualization (NetworkX)                                |
| `genecircuitry/plotting/`                | Canonical plot generation (`qc_plots.py`, `grn_plots.py`, etc.) |
| `genecircuitry/reporting/`               | HTML/PDF report generation                                      |

## Critical Patterns

### 1. Never Hardcode Parameters

```python
# ✅ CORRECT - use config
from genecircuitry import config
sc.pp.filter_cells(adata, min_genes=config.QC_MIN_GENES)
plt.figure(figsize=config.PLOT_FIGSIZE_MEDIUM)
plt.savefig(f"{config.FIGURES_DIR_GRN}/plot.png", dpi=config.SAVE_DPI)

# ❌ WRONG - hardcoded values
sc.pp.filter_cells(adata, min_genes=200)
```

### 2. Extend PipelineController for New Steps

```python
# In genecircuitry/pipeline/controller.py - add new method to PipelineController:
def run_step_new_analysis(self, adata, log_dir=None):
    """Execute new analysis step."""
    log_step("Controller.NewAnalysis", "STARTED")
    try:
        result = my_analysis_function(adata)
        log_step("Controller.NewAnalysis", "COMPLETED")
        return result
    except Exception as e:
        log_error("Controller.NewAnalysis", e)
        raise
```

**NEVER** create standalone scripts - always integrate into the controller.

### 3. Logging Pattern

```python
# In PipelineController methods only — log_step lives in controller.py:
from genecircuitry.pipeline import log_step
log_step("MyStep", "STARTED", {"n_cells": adata.n_obs})
log_step("MyStep", "COMPLETED", {"result_count": len(results)})

# In all other submodules — use logging_utils (safe before setup_logging):
from genecircuitry.logging_utils import log_error, log_warning
log_error("MyModule", exception)   # Logs to error.log with traceback
log_warning("MyModule", "message")  # Logs to pipeline.log
```

### 4. Function Signatures with Config Defaults

```python
def my_function(
    adata: AnnData,
    param: Optional[int] = None,  # Allow override
) -> AnnData:
    if param is None:
        param = config.MY_PARAM  # Fall back to config
```

### 5. Optional Dependency Pattern

`celloracle_processing` and `hotspot_processing` use optional dependencies. `genecircuitry/__init__.py` wraps them in `try/except`; access via `genecircuitry.celloracle_processing` (may be `None` if not installed). New optional-dep modules must follow the same pattern.

Install optional groups: `pip install genecircuitry[grn]` (CellOracle), `pip install genecircuitry[hotspot]`, `pip install genecircuitry[all]`.

### 6. Plotting — Use `genecircuitry/plotting/` Subpackage

Canonical plotting lives in `genecircuitry/plotting/` (`qc_plots.py`, `grn_plots.py`, `hotspot_plots.py`). Legacy inline plotting in `grn_deep_analysis.py` and `hotspot_processing.py` is deprecated. Always add new plots to the `plotting/` subpackage.

## Developer Workflows

### Run Pipeline

```bash
# Full run with example data
python run_complete_analysis.py

# With custom data
python run_complete_analysis.py --input data.h5ad --output results

# Specific steps only (checkpoints auto-resume)
python -m genecircuitry.pipeline --steps load preprocessing clustering

# Stratified parallel analysis (by cell type)
python -m genecircuitry.pipeline --cluster-key-stratification celltype --parallel --n-jobs 4
```

**Available step names:** `load`, `preprocessing`, `stratification`, `clustering`, `celloracle`, `hotspot`, `grn_analysis`, `summary`

### Development Commands (preferred: `pixi`)

```bash
pixi run -e dev test        # pytest tests/ -v
pixi run -e dev lint        # flake8 genecircuitry/
pixi run -e dev format      # black genecircuitry/
pixi run -e dev typecheck   # mypy genecircuitry/
```

Direct fallback: `pytest tests/ -v`

**Python version**: `>=3.9, <3.11` only.

When adding new config parameters, add corresponding tests to `tests/test_config.py`.

### Output Structure

```
output/
├── preprocessed_adata.h5ad
├── report.html                  # Interactive HTML report
├── report.pdf                   # PDF report (requires weasyprint)
├── logs/{pipeline.log, error.log, *.checkpoint}
├── figures/{qc/, grn/, hotspot/}
├── celloracle/{grn_merged_scores.csv, grn_filtered_links.pkl}
├── hotspot/{autocorrelation_results.csv, gene_modules.csv}
└── stratified_analysis/<ClusterName>/  # Per-cluster results
```

## Reporting Module

```python
from genecircuitry.reporting import generate_report

outputs = generate_report(
    output_dir="results/",
    title="My Analysis Report",
    adata=adata,
    celloracle_result=(oracle, links),
    hotspot_result=hs,
    log_file="results/logs/pipeline.log",
    embed_images=True,   # False = use relative paths (smaller file)
    formats=["html", "pdf"],
)
```

## AnnData Conventions

- **Metrics**: `.obs` (per-cell), `.var` (per-gene)
- **Embeddings**: `.obsm['X_pca']`, `.obsm['X_umap']`
- **Raw counts**: Store in `.layers['raw_count']` before normalization
- **Import**: `import scanpy as sc`, access as `adata`

## Adding New Config Parameters

1. Add to `genecircuitry/config.py` with docstring
2. Add to `get_config()` return dict in same file
3. Add test in `tests/test_config.py`

## Checkpoint System

Checkpoints skip completed steps on re-runs. Stored in `output/logs/<step_name>.checkpoint` (JSON). Key functions in `controller.py`:

- `compute_input_hash(input_path, **params)` — MD5 of file path + mtime + params
- `write_checkpoint(log_dir, step_name, input_hash, **metadata)` — saves checkpoint
- `check_checkpoint(log_dir, step_name, input_hash)` — returns `True` to skip step

See [docs/checkpoints.md](../docs/checkpoints.md) for full details.

## Where to Make Changes

| Task                                | Location                                                                  |
| ----------------------------------- | ------------------------------------------------------------------------- |
| Add/change pipeline parameter       | `genecircuitry/config.py` + `get_config()` + `tests/test_config.py`       |
| Add pipeline step                   | `PipelineController` in `genecircuitry/pipeline/controller.py`            |
| Add a plot                          | `genecircuitry/plotting/<module>.py` + export in `plotting/__init__.py`   |
| Add a report section                | `genecircuitry/reporting/sections.py` + export in `reporting/__init__.py` |
| Log from any submodule              | `from genecircuitry.logging_utils import log_error, log_warning`          |
| Log step progress (controller only) | `from genecircuitry.pipeline import log_step`                             |

## Notes & Known Dependencies (see [docs/implementation.md](../docs/implementation.md))

- `pyproject.toml` optional-dep groups (`enrichment`, `atac`) are in optional extras (`gseapy>=1.0.0`, `genomepy`, `gimmemotifs`).
- PDF generation requires `weasyprint` which depends on system-level `pango` — use `pixi` env or install separately.
