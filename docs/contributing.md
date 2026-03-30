---
layout: default
title: Contributing
nav_order: 7
---

# Contributing

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

## Setup for development

```bash
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd GeneCircuitry
python -m venv venv && source venv/bin/activate
pip install -e ".[dev,grn,hotspot]"
```

---

## Running tests

```bash
pytest tests/                      # all tests
pytest tests/test_config.py       # config-specific
pytest -v --tb=short              # verbose with short tracebacks
```

When adding new config parameters, add a corresponding test to `tests/test_config.py`.

---

## Code conventions

### Parameters must use config

```python
# ✅ correct
from genecircuitry import config

def my_function(adata, threshold=None):
    if threshold is None:
        threshold = config.MY_THRESHOLD
```

### New pipeline steps go in PipelineController

```python
# In genecircuitry/pipeline/controller.py
def run_step_my_analysis(self, adata, log_dir=None):
    log_step("Controller.MyAnalysis", "STARTED")
    try:
        result = my_analysis_function(adata)
        log_step("Controller.MyAnalysis", "COMPLETED")
        return result
    except Exception as e:
        log_error("Controller.MyAnalysis", e)
        raise
```

Then add `"my_analysis"` to the `steps` list in `run_complete_pipeline()`.

### New plots go in `genecircuitry/plotting/`

Do **not** add plotting code to `preprocessing.py`, `grn_deep_analysis.py`, or `hotspot_processing.py`. Create or extend the relevant file in `genecircuitry/plotting/`.

### Optional dependencies

Wrap new optional-dep modules in `genecircuitry/__init__.py`:

```python
try:
    from . import my_new_module
except ImportError:
    my_new_module = None
```

---

## Adding a new config parameter

1. Add the constant to `genecircuitry/config.py` with a docstring.
2. Add it to the `get_config()` return dict in the same file.
3. Add `assert "MY_PARAM" in config` to `tests/test_config.py`.

---

## AnnData conventions

| Location               | Usage                                                  |
| ---------------------- | ------------------------------------------------------ |
| `.obs`                 | Per-cell metrics (QC values, cluster labels)           |
| `.var`                 | Per-gene flags (`mt`, `ribo`, `hb`, `highly_variable`) |
| `.obsm['X_pca']`       | PCA embedding                                          |
| `.obsm['X_umap']`      | UMAP embedding                                         |
| `.layers['raw_count']` | Raw counts (stored before normalization)               |

---

## Pull request checklist

- [ ] No hardcoded numeric values — all thresholds reference `config.*`
- [ ] New config parameters have tests in `tests/test_config.py`
- [ ] New plots added to `genecircuitry/plotting/`, not inline in processing modules
- [ ] New pipeline steps integrated into `PipelineController`, not standalone scripts
- [ ] Optional dependencies wrapped in `try/except ImportError`
- [ ] Docstrings include example imports from the correct module
