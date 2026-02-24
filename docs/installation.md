---
layout: default
title: Installation
nav_order: 2
---

# Installation

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

## Requirements

- Python ≥ 3.9
- Core dependencies installed automatically via pip: `numpy`, `pandas`, `scanpy`, `anndata`, `matplotlib`, `seaborn`, `scipy`, `adjustText`, `networkx`, `marsilea`

---

## From source (recommended)

```bash
git clone https://github.com/samuelecancellieri/GeneCircuitry.git
cd GeneCircuitry

# Create a virtual environment
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate

# Install core package
pip install -e .

# Install with GRN inference support (CellOracle)
pip install -e ".[grn]"

# Install with Hotspot module support
pip install -e ".[hotspot]"

# Install everything
pip install -e ".[grn,hotspot,enrichment,atac]"

# Development dependencies (testing, linting)
pip install -e ".[dev]"
```

---

## Optional dependency groups

| Extra        | Packages                            | Description                  |
| ------------ | ----------------------------------- | ---------------------------- |
| `grn`        | `celloracle`                        | GRN inference via CellOracle |
| `hotspot`    | `hotspot-sc`                        | Gene module identification   |
| `enrichment` | `gseapy>=1.0.0`                     | Gene set enrichment analysis |
| `atac`       | `genomepy`, `gimmemotifs`           | ATAC-seq peak motif analysis |
| `dev`        | `pytest`, `black`, `flake8`, `mypy` | Development tools            |

{: .note }

> **CellOracle and Hotspot are optional.** `import genecircuitry` will succeed even when they are not installed. The relevant processing modules (`genecircuitry.celloracle_processing`, `genecircuitry.hotspot_processing`) will be `None` at runtime and gracefully skipped by the pipeline.

---

## Verify installation

```python
import genecircuitry
print(genecircuitry.__version__)           # 0.1.0

# Check optional deps
print(genecircuitry.celloracle_processing) # None if not installed
print(genecircuitry.hotspot_processing)    # None if not installed
```

---

## Example data

The repository ships with `data/pbmc3k_raw.h5ad` (PBMC 3k dataset) and `data/paul15/paul15.h5` (Paul et al. 2015, mouse hematopoiesis) for testing. The pipeline also falls back to `sc.datasets.pbmc3k()` if no input file is provided.
