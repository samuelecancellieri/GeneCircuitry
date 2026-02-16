# TRNspot Repository Restructuring Plan

> **Generated**: 2026-02-13  
> **Status**: IMPLEMENTED — all restructuring changes have been applied. Files to remove are renamed with `_TODELETE` suffix. See summary below.

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Files Marked TO REMOVE](#2-files-marked-to-remove)
3. [Files to Move](#3-files-to-move)
4. [Bugs Found](#4-bugs-found)
5. [Structural Issues & Fixes](#5-structural-issues--fixes)
6. [Proposed Clean Directory Layout](#6-proposed-clean-directory-layout)
7. [Action Checklist](#7-action-checklist)

---

## 1. Executive Summary

The repo has accumulated technical debt from an abandoned subpackage refactoring, broken example scripts, duplicated plotting code, stale development notes at the root, and a redundant `setup.py`. This plan categorizes every issue and proposes a clean target structure.

**Key findings:**

- **6 ghost subpackages** (`core/`, `preprocessing/`, `computation/`, `computation/celloracle/`, `computation/hotspot/`, `pipeline/`) contain only `__pycache__/` — source files were deleted but bytecache remains
- **5 broken example/test files** import functions that don't exist (`plot_qc_violin`, `plot_qc_scatter`, `perform_grn_pre_processing` from wrong module)
- **Redundant build config**: both `setup.py` and `pyproject.toml`
- **Duplicate plotting code** in `grn_deep_analysis.py` / `hotspot_processing.py` vs. dedicated `plotting/` subpackage
- **Dev notes cluttering root**: 3 implementation docs belong in `docs/`
- **Mislabeled CHANGELOG.md**: contains preprocessing config notes, not version history
- **1 dev-only verification script** (`verify_controller.py`) at root
- **1 lab-specific test script** (`run_pipe_regev_test.sh`) at root
- **Eager imports of optional deps** in `__init__.py` (CellOracle, Hotspot)

---

## 2. Files Marked TO REMOVE

Each file below has been annotated with a `# TO REMOVE` (or `<!-- TO REMOVE -->` for markdown/toml) comment at the top explaining why.

### 2.1 Ghost Subpackages (empty — only `__pycache__/`)

| Path                                  | Reason                                                                                                                           |
| ------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `trnspot/core/`                       | Abandoned refactoring — only `__pycache__/`, no source files                                                                     |
| `trnspot/preprocessing/` (subpackage) | Abandoned refactoring — only `__pycache__/`, no source files. Not to be confused with `trnspot/preprocessing.py` which is active |
| `trnspot/computation/`                | Abandoned refactoring — only `__pycache__/`, no source files                                                                     |
| `trnspot/pipeline/`                   | Abandoned refactoring — only `__pycache__/`, no source files                                                                     |

**Cleanup command:** `find trnspot/ -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null; rmdir trnspot/core trnspot/preprocessing trnspot/computation/celloracle trnspot/computation/hotspot trnspot/computation trnspot/pipeline 2>/dev/null`

### 2.2 Redundant / Stale Files

| File                        | Reason                                                                                                                                                                  | Marked |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ |
| `setup.py`                  | Redundant — fully superseded by `pyproject.toml`                                                                                                                        | ✅     |
| `verify_controller.py`      | One-time dev verification script, not needed post-implementation                                                                                                        | ✅     |
| `CHANGELOG.md`              | **Mislabeled** — contains preprocessing config integration notes, not a real changelog. Content belongs in `docs/PREPROCESSING_CONFIG_UPDATE.md` (which already exists) | ✅     |
| `trnspot.egg-info/`         | Build artifact — should be in `.gitignore`, not committed                                                                                                               | ✅     |
| `tests/test_preprocessing/` | Empty directory — no files                                                                                                                                              | ✅     |
| `examples/__pycache__/`     | Build artifact — should be in `.gitignore`                                                                                                                              | ✅     |

### 2.3 Broken Example Scripts

| File                                  | Broken Import                                                                                                                      | Marked |
| ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ------ |
| `examples/quick_demo.py`              | `from trnspot.preprocessing import plot_qc_violin, plot_qc_scatter` — these functions don't exist in `preprocessing.py`            | ✅     |
| `examples/example_qc.py`              | Same broken imports as above                                                                                                       | ✅     |
| `examples/test_config_integration.py` | `from trnspot.preprocessing import perform_grn_pre_processing` — function is in `celloracle_processing.py`, not `preprocessing.py` | ✅     |

### 2.4 Lab-Specific Scripts

| File                     | Reason                                                                                                      | Marked |
| ------------------------ | ----------------------------------------------------------------------------------------------------------- | ------ |
| `run_pipe_regev_test.sh` | Hardcoded lab paths — not portable, should be a personal script or in a `scripts/` folder with `.gitignore` | ✅     |

---

## 3. Files to Move (not remove)

| File                            | Current Location | Proposed Location                       | Reason                                                                                                                                    |
| ------------------------------- | ---------------- | --------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `CONTROLLER_IMPLEMENTATION.md`  | Root             | `docs/dev/CONTROLLER_IMPLEMENTATION.md` | Dev notes — cluttering root                                                                                                               |
| `LOGGING_SYSTEM.md`             | Root             | `docs/dev/LOGGING_SYSTEM.md`            | Dev notes — cluttering root                                                                                                               |
| `PARALLEL_PROCESSING_FIX.md`    | Root             | `docs/dev/PARALLEL_PROCESSING_FIX.md`   | Dev notes — cluttering root                                                                                                               |
| `examples/complete_pipeline.py` | `examples/`      | `trnspot/pipeline/controller.py`        | This is the **core pipeline engine** (1748 lines), not an "example". It's the central orchestrator imported by `run_complete_analysis.py` |

---

## 4. Bugs Found

These are not about removal but need fixing:

### 4.1 Missing Config Constant

- **File**: `trnspot/grn_deep_analysis.py` (referenced on line ~79)
- **Bug**: `config.PLOT_FIGSIZE_WIDE` is used but never defined in `config.py`
- **Fix**: Add `PLOT_FIGSIZE_WIDE = (20, 8)` to `config.py` or replace with `PLOT_FIGSIZE_LARGE`

### 4.2 Broken `perform_qc()` Signature

- **File**: `trnspot/preprocessing.py`
- **Bug**: Docstring documents `plot: bool` parameter but function signature uses `save_plots: str`. Tests/examples calling `plot=False` will crash with `TypeError`
- **Fix**: Either add `plot` parameter or update docstring and all callers

### 4.3 Broken Test File

- **File**: `tests/test_preprocessing.py`
- **Bug**: Imports `plot_qc_violin`, `plot_qc_scatter` from `trnspot.preprocessing` (don't exist). Calls `perform_qc(adata, plot=False)` (parameter doesn't exist)
- **Fix**: Update imports to use `trnspot.plotting.qc_plots` and fix `perform_qc` calls

### 4.4 Broken Example Imports (3 files)

- `examples/celloracle_workflow.py` — imports `perform_grn_pre_processing` from `trnspot.preprocessing` (it's in `trnspot.celloracle_processing`)
- `examples/preprocess_workflow.py` — same broken import
- `examples/test_config_integration.py` — same broken import

### 4.5 Eager Optional Dependency Imports

- **File**: `trnspot/__init__.py`
- **Bug**: `from . import celloracle_processing` and `from . import hotspot_processing` run at import time, so `import trnspot` fails if CellOracle/Hotspot aren't installed
- **Fix**: Use lazy imports or wrap in try/except

### 4.6 Incomplete Dependencies

- **File**: `pyproject.toml` and `requirements.txt`
- **Missing from both**: `networkx`, `marsilea`, `gseapy`
- **Missing from `pyproject.toml`**: `matplotlib`, `seaborn`, `scipy`, `adjustText`
- **Fix**: Sync `pyproject.toml` dependencies with `requirements.txt` and add missing ones. Consider `[project.optional-dependencies]` for `celloracle`, `hotspot`, `gseapy`, `genomepy`, `gimmemotifs`

### 4.7 Missing LICENSE file

- **File**: `MANIFEST.in` references `include LICENSE` but no LICENSE file exists
- **Fix**: Add MIT LICENSE file (as declared in `pyproject.toml`)

---

## 5. Structural Issues & Fixes

### 5.1 Duplicate Plotting Code

**Problem**: Plotting functions exist in **two places**:

1. Legacy inline code in `grn_deep_analysis.py` and `hotspot_processing.py`
2. Dedicated `trnspot/plotting/` subpackage (`grn_plots.py`, `hotspot_plots.py`)

**Resolution**: After confirming `plotting/` subpackage covers all functionality, deprecate the inline plotting functions in the legacy modules. Add `# DEPRECATED: Use trnspot.plotting.grn_plots instead` comments.

### 5.2 Monolithic `complete_pipeline.py` (1748 lines)

**Problem**: The entire pipeline engine lives in `examples/complete_pipeline.py` — a misleading location for the core orchestrator.

**Resolution**: Move to `trnspot/pipeline/controller.py` and split into:

- `controller.py` — `PipelineController` class
- `steps.py` — Individual pipeline step functions
- `cli.py` — CLI argument parser and `main()`
- `logging.py` — `log_step()`, `log_error()`, checkpoint system

### 5.3 `trnspot.egg-info/` Committed to Repo

**Problem**: Build artifacts shouldn't be version-controlled.

**Resolution**: Add to `.gitignore` and remove from repo.

---

## 6. Proposed Clean Directory Layout

```
TRNspot_simple/
├── .gitignore                      # Add: *.egg-info/, __pycache__/, *.pyc
├── LICENSE                         # NEW — MIT license (currently missing)
├── README.md                       # Keep as-is
├── pyproject.toml                  # Keep — FIX dependencies
├── requirements.txt                # Keep — sync with pyproject.toml
├── requirements-dev.txt            # Keep as-is
├── run_complete_analysis.py        # Keep — entry point
│
├── docs/
│   ├── CONFIG.md
│   ├── COMPLETE_PIPELINE_GUIDE.md
│   ├── CONTROLLER_QUICK_REF.md
│   ├── PACKAGE_STRUCTURE.md
│   ├── PIPELINE_CONTROLLER_GUIDE.md
│   ├── QC_FUNCTIONS.md
│   ├── CELLORACLE_MODULE_SUMMARY.md
│   ├── CELLORACLE_PROCESSING.md
│   ├── CELLORACLE_QUICK_REF.md
│   ├── PREPROCESSING_CONFIG_UPDATE.md
│   └── dev/                         # NEW — move dev notes here
│       ├── CONTROLLER_IMPLEMENTATION.md
│       ├── LOGGING_SYSTEM.md
│       └── PARALLEL_PROCESSING_FIX.md
│
├── examples/
│   ├── config_example.py            # Keep
│   ├── controller_usage_examples.py # Keep
│   ├── celloracle_workflow.py       # Keep — FIX imports
│   ├── hotspot_workflow.py          # Keep
│   └── preprocess_workflow.py       # Keep — FIX imports
│
├── tests/
│   ├── __init__.py
│   ├── test_config.py               # Keep
│   ├── test_celloracle.py           # Keep
│   └── test_preprocessing.py        # Keep — FIX imports
│
├── trnspot/
│   ├── __init__.py                  # FIX — lazy imports for optional deps
│   ├── config.py                    # FIX — add PLOT_FIGSIZE_WIDE
│   ├── preprocessing.py             # FIX — plot parameter
│   ├── celloracle_processing.py     # Keep
│   ├── hotspot_processing.py        # Keep — deprecate inline plots
│   ├── grn_deep_analysis.py         # Keep — deprecate inline plots
│   ├── atac_peaks_processing.py     # Keep
│   ├── enrichment_analysis.py       # Keep
│   ├── plotting/                    # Keep — canonical plotting location
│   │   ├── __init__.py
│   │   ├── utils.py
│   │   ├── qc_plots.py
│   │   ├── grn_plots.py
│   │   └── hotspot_plots.py
│   ├── reporting/                   # Keep as-is
│   │   ├── __init__.py
│   │   ├── generator.py
│   │   └── sections.py
│   └── pipeline/                    # NEW — move complete_pipeline.py here
│       ├── __init__.py
│       └── controller.py            # Was examples/complete_pipeline.py
│
├── data/                            # Keep — example data
│   └── pbmc3k_raw.h5ad
│
└── results/                         # Keep — gitignored output
    └── ...
```

### What's Removed in Clean Layout

- `setup.py`
- `verify_controller.py`
- `CHANGELOG.md` (mislabeled)
- `CONTROLLER_IMPLEMENTATION.md` → moved to `docs/dev/`
- `LOGGING_SYSTEM.md` → moved to `docs/dev/`
- `PARALLEL_PROCESSING_FIX.md` → moved to `docs/dev/`
- `run_pipe_regev_test.sh`
- `trnspot.egg-info/`
- `trnspot/core/` (ghost)
- `trnspot/preprocessing/` subpackage (ghost)
- `trnspot/computation/` (ghost)
- `trnspot/pipeline/` (ghost — will be recreated with real source)
- `tests/test_preprocessing/` (empty)
- `examples/quick_demo.py` (broken)
- `examples/example_qc.py` (broken)
- `examples/test_config_integration.py` (broken)
- `examples/__pycache__/`
- All `__pycache__/` directories

---

## 7. Action Checklist

- [ ] **Phase 1: Cleanup** — Remove files marked `TO REMOVE`
  - [ ] Delete ghost subpackages (`trnspot/core/`, `trnspot/preprocessing/`, `trnspot/computation/`, `trnspot/pipeline/`)
  - [ ] Delete `setup.py`, `verify_controller.py`, `CHANGELOG.md`
  - [ ] Delete broken examples (`quick_demo.py`, `example_qc.py`, `test_config_integration.py`)
  - [ ] Delete `run_pipe_regev_test.sh`
  - [ ] Delete `trnspot.egg-info/`, `tests/test_preprocessing/`
  - [ ] Clean all `__pycache__/` directories

- [ ] **Phase 2: Move**
  - [ ] Move `CONTROLLER_IMPLEMENTATION.md` → `docs/dev/`
  - [ ] Move `LOGGING_SYSTEM.md` → `docs/dev/`
  - [ ] Move `PARALLEL_PROCESSING_FIX.md` → `docs/dev/`
  - [ ] Move `examples/complete_pipeline.py` → `trnspot/pipeline/controller.py`
  - [ ] Update `run_complete_analysis.py` import path

- [ ] **Phase 3: Fix Bugs**
  - [ ] Add `PLOT_FIGSIZE_WIDE` to `config.py`
  - [ ] Fix `perform_qc()` signature (`plot` parameter)
  - [ ] Fix broken imports in `celloracle_workflow.py`, `preprocess_workflow.py`
  - [ ] Fix `test_preprocessing.py` imports and calls
  - [ ] Make `__init__.py` imports lazy for optional deps
  - [ ] Complete dependencies in `pyproject.toml` and `requirements.txt`
  - [ ] Add `LICENSE` file
  - [ ] Fix `MANIFEST.in`

- [ ] **Phase 4: .gitignore**
  - [ ] Add `*.egg-info/`, `__pycache__/`, `*.pyc`, `results/`, `regev_test_run/`
