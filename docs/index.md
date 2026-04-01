---
layout: home
title: Home
nav_order: 1
---

# GeneCircuitry

**Transcriptional Regulatory Network analysis for single-cell RNA-seq data**

GeneCircuitry is a Python package that integrates [Scanpy](https://scanpy.readthedocs.io/), [CellOracle](https://celloracle.org/), and [Hotspot](https://hotspot.readthedocs.io/) into a single **modular, checkpoint-enabled pipeline** for gene regulatory network (GRN) inference from single-cell data. Run the full workflow with one command, or execute individual steps selectively and resume from where you left off.

---

## What GeneCircuitry does

| Step                                  | Tool                       | Output                                       |
| ------------------------------------- | -------------------------- | -------------------------------------------- |
| Quality control & normalization       | Scanpy                     | Filtered, normalized `AnnData`               |
| Dimensionality reduction & clustering | Scanpy (PCA, UMAP, Leiden) | Cell embeddings + cluster labels             |
| GRN inference                         | CellOracle                 | Per-cluster transcription factor networks    |
| Gene module identification            | Hotspot                    | Spatially autocorrelated gene modules        |
| Network visualization                 | NetworkX                   | Centrality plots, network graphs, rank plots |
| Reporting                             | Built-in HTML/PDF engine   | Interactive analysis report                  |
| Stratified analysis                   | `PipelineController`       | Per-cell-type parallel runs                  |

---

## Quick install

**Pixi (recommended)** — installs all conda and pip dependencies in one step:

```bash
curl -fsSL https://pixi.sh/install.sh | bash   # install pixi once
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd genecircuitry
pixi install
pixi run run        # launch the pipeline
```

**pip / venv:**

```bash
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd genecircuitry
pip install -e ".[grn,hotspot]"
```

**Docker:**

```bash
docker run --rm \
    -v /path/to/your/data:/data \
    -v /path/to/output:/output \
    zanathos/genecircuitry:latest \
    --input /data/your_data.h5ad --output /output
```

**Singularity / Apptainer** (HPC clusters):

```bash
singularity pull genecircuitry.sif docker://zanathos/genecircuitry:latest
singularity exec \
    --bind /path/to/your/data:/data \
    --bind /path/to/output:/output \
    genecircuitry.sif genecircuitry --input /data/your_data.h5ad --output /output
```

See [Getting Started → Installation](installation/) for all five options, including Conda, version pins, and HPC-specific tips.

---

## Run your first analysis

```bash
# Run the full pipeline on the bundled PBMC 3k demo dataset
python run_complete_analysis.py

# Run on your own data
python run_complete_analysis.py --input my_data.h5ad --output results/

# Run only specific steps
python -m genecircuitry.pipeline --input my_data.h5ad --output results/ \
    --steps load preprocessing clustering
```

See [Getting Started → Quick Start](quickstart/) for the full walkthrough.

---

## Documentation

### Getting Started

- [Installation](installation/) — pip, conda, pixi, Docker
- [Quick Start](quickstart/) — CLI and Python API walkthrough

### User Guide

- [Pipeline Overview](pipeline/) — step names, CLI flags, parallel execution, architecture
- [Preprocessing & QC](preprocessing/) — filtering, normalization, clustering
- [GRN Inference](celloracle/) — CellOracle Oracle objects and link inference
- [Gene Modules](hotspot/) — Hotspot autocorrelation and module detection
- [Plotting System](plotting/) — `genecircuitry/plotting/` subpackage
- [Reporting](reporting/) — HTML/PDF report generation
- [Stratified Analysis](stratification/) — per-cluster parallel runs
- [Checkpoints & Resume](checkpoints/) — auto-resume and checkpoint management

### Reference

- [Configuration](configuration/) — all `config.py` parameters
- [API Reference](api/) — function signatures across all modules
- [ATAC Peaks Processing](atac-peaks/) — BED → TF motif matrix for CellOracle

### Other

- [Architecture](architecture/) — codebase design, data flow, design patterns
- [Contributing](contributing/) — code conventions, PR checklist

---

## Data flow

```
AnnData (.h5ad)
    │
    ▼
[Preprocessing & QC]            genecircuitry/preprocessing.py
  ├─ perform_qc()
  ├─ perform_normalization()
  └─ perform_dimensionality_reduction_clustering()
    │
    ▼  (checkpoint: preprocessed_adata.h5ad)
    │
    ├──────────────────────────────────────────────────┐
    ▼                                                  │
[Optional: Stratification]       split by cell type    │
  └─ per-cluster AnnData objects                       │
    │                                                  │
    ├────────────────────────┐                         │
    ▼                        ▼                         │
[CellOracle GRN]        [Hotspot Modules]              │
  create_oracle_object()  create_hotspot_object()      │
  run_PCA()               run_hotspot_analysis()       │
  run_KNN()                                            │
  run_links()                                          │
    │                        │                         │
    └────────────┬───────────┘                         │
                 ▼                                     │
         [GRN Deep Analysis]   genecircuitry/grn_deep_analysis.py
         [Reporting]           genecircuitry/reporting/ → HTML + PDF
                 ▲
                 └─────────────────────────────────────┘
```
