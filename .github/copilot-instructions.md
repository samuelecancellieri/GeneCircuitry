# GeneCircuitry AI Coding Instructions

**GeneCircuitry** is a Python package for transcriptional regulatory network (TRN) analysis using single-cell data. It integrates Scanpy, CellOracle, and Hotspot into a modular, checkpoint-enabled pipeline.

## Architecture

### Data Flow

`AnnData (.h5ad)` → **Preprocessing** (QC/normalization) → **Clustering** → **CellOracle** (GRN inference) + **Hotspot** (gene modules) → **Deep Analysis** (visualization) → **Reporting** (HTML/PDF)

### Key Components

| File                               | Purpose                                                         |
| ---------------------------------- | --------------------------------------------------------------- |
| `genecircuitry/config.py`                | **SINGLE SOURCE OF TRUTH** - all parameters                     |
| `genecircuitry/pipeline/controller.py`   | `PipelineController` class - central orchestrator (~1900 lines) |
| `run_complete_analysis.py`         | Entry point wrapper (calls `genecircuitry.pipeline.main`)             |
| `genecircuitry/preprocessing.py`         | Scanpy wrappers (QC, normalize, cluster)                        |
| `genecircuitry/celloracle_processing.py` | GRN inference with CellOracle (optional dep)                    |
| `genecircuitry/hotspot_processing.py`    | Spatial autocorrelation modules (optional dep)                  |
| `genecircuitry/grn_deep_analysis.py`     | Network visualization (NetworkX, Marsilea)                      |
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
from genecircuitry.pipeline import log_step, log_error

log_step("MyStep", "STARTED", {"n_cells": adata.n_obs})
log_step("MyStep", "COMPLETED", {"result_count": len(results)})
log_error("MyStep", exception)  # Logs to error.log with traceback
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

### Testing

```bash
pytest tests/                    # Run all tests
pytest tests/test_config.py     # Config tests specifically
```

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

## Known Remaining Issues (see `RESTRUCTURING_PLAN.md`)

- `pyproject.toml` optional-dep groups (`enrichment`, `atac`) are not pinned — `gseapy>=1.0.0`, `genomepy`, `gimmemotifs` are in optional extras only, not core deps
