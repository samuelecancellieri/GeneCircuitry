# Changelog

All notable changes to GeneCircuitry are documented here.

## [0.1.7] - 2026-03-31

### Fixed

- Corrected hotspot install hint in `controller.py`: `pip install hotspot-sc` → `pip install hotspotsc`

### Dependencies

- Pinned `pandas>=1.0.3,<=1.5.3` (was `>=1.3.0`) for compatibility with CellOracle
- Pinned `matplotlib>=3.6.3,<3.7` for stable rendering
- Added `velocyto>=0.17.17` explicitly to the `grn` optional-dependency group
- Constrained `setuptools>=62.1.0,<81.0.0` in build requirements (setuptools 81+ has breaking changes)
- Added `Cython` to build requirements (needed by velocyto's `setup.py`)

### Infrastructure

- **Docker**: Added `gcc` installation; pinned `setuptools<81` and `cython` in build layer; improved layer caching
- **CI – pixi-test**: New GitHub Actions workflow that installs via pixi, runs the test suite, and caches the pixi environment
- **CI – conda-build**: New workflow that validates the conda recipe on every push/PR
- **CI – docker**: Refactored workflow; removed redundant `docker-image.yml`
- Dropped macOS (`osx-arm64`, `osx-64`) from pixi platform list — Linux-only until optional deps support macOS
- Added `pango` conda dependency (required by WeasyPrint for PDF export)
- Added `no-build-isolation = ["velocyto"]` to `[pypi-options]` so velocyto builds against the conda env's numpy/cython
- Moved per-feature dev tasks (`test`, `lint`, `format`, `typecheck`) under `[feature.dev.tasks]` in `pixi.toml`

## [0.1.6] - 2026-02-13

### Fixed

- Fixed velocyto dependency in Dockerfile
- Corrected `hotspotsc` package name in pixi/pyproject configs

### Infrastructure

- Moved velocyto to `pip install` step inside pixi environment
- Updated pixi.toml for compatibility with velocyto build requirements
- Updated pandas dependency constraints across all config files

## [0.1.5] - 2025-11-20

### Added

- Stratified parallel analysis support (`--cluster-key-stratification`, `--parallel`, `--n-jobs`)
- Centralized plotting subpackage (`genecircuitry/plotting/`)
- HTML and PDF report generation (`genecircuitry/reporting/`)
- Checkpoint-enabled pipeline resume

### Changed

- Migrated from `setup.py` to `pyproject.toml`
- All parameters now sourced from `genecircuitry/config.py`

## [0.1.0] – [0.1.4]

Initial development releases establishing the core pipeline, CellOracle integration, Hotspot integration, and GRN deep analysis.
