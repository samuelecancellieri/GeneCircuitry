---
layout: home
title: GeneCircuitry
nav_order: 1
---

# GeneCircuitry

**Transcriptional Regulatory Network analysis for single-cell data**

GeneCircuitry is a Python package that integrates [Scanpy](https://scanpy.readthedocs.io/), [CellOracle](https://celloracle.org/), and [Hotspot](https://hotspot.readthedocs.io/) into a modular, checkpoint-enabled pipeline for gene regulatory network (GRN) inference from single-cell RNA-seq data.

---

## What GeneCircuitry does

| Step                                  | Tool                       | Output                                    |
| ------------------------------------- | -------------------------- | ----------------------------------------- |
| Quality control & normalization       | Scanpy                     | Filtered, normalized `AnnData`            |
| Dimensionality reduction & clustering | Scanpy (PCA, UMAP, Leiden) | Cell embeddings + cluster labels          |
| GRN inference                         | CellOracle                 | Per-cluster transcription factor networks |
| Gene module identification            | Hotspot                    | Spatially autocorrelated gene modules     |
| Visualization                         | NetworkX, Marsilea         | Network plots, heatmaps, scatter plots    |
| Reporting                             | Custom HTML/PDF engine     | Interactive analysis report               |

---

## Quick navigation

- [Installation](installation/) — how to install GeneCircuitry and its dependencies
- [Quick Start](quickstart/) — run your first analysis in minutes
- [Architecture](architecture/) — understand the codebase design
- [API Reference](api/) — module-level function documentation
- [Configuration](configuration/) — all parameters in `genecircuitry/config.py`
- [Contributing](contributing/) — how to extend GeneCircuitry

---

## Data flow

```
AnnData (.h5ad)
    │
    ▼
Preprocessing    ← genecircuitry/preprocessing.py
(QC, normalize, HVG, PCA, clustering)
    │
    ├──────────────────────┐
    ▼                      ▼
CellOracle GRN         Hotspot modules
(celloracle_processing) (hotspot_processing)
    │                      │
    └──────────┬───────────┘
               ▼
          Deep Analysis    ← genecircuitry/grn_deep_analysis.py
          (NetworkX, Marsilea plots)
               │
               ▼
          HTML / PDF Report ← genecircuitry/reporting/
```

---

## Quick install

```bash
pip install -e ".[grn,hotspot]"   # with CellOracle + Hotspot
pip install -e "."                # core only
```

See [Installation](installation/) for full instructions.
