---
layout: default
title: User Guide
nav_order: 3
has_children: true
---

# User Guide

In-depth documentation for each stage of the GeneCircuitry pipeline, from QC through reporting.

| Page                                         | What you'll find                                                                                     |
| -------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| [Pipeline Overview](../pipeline/)            | `PipelineController` architecture, step names, CLI flags, parallel execution, checkpoint integration |
| [Preprocessing & QC](../preprocessing/)      | Cell/gene filtering, normalization, HVG selection, dimensionality reduction, clustering              |
| [GRN Inference (CellOracle)](../celloracle/) | Oracle object creation, PCA, KNN imputation, link inference                                          |
| [Gene Modules (Hotspot)](../hotspot/)        | Autocorrelation testing, module detection, enrichment-annotated heatmaps                             |
| [Plotting System](../plotting/)              | `genecircuitry/plotting/` subpackage — QC, GRN, and Hotspot canonical plot functions                 |
| [Reporting](../reporting/)                   | HTML and PDF report generation with `generate_report()`                                              |
| [Stratified Analysis](../stratification/)    | Per-cluster parallel analysis, output layout, configuration                                          |
| [Checkpoints & Resume](../checkpoints/)      | How `.checkpoint` files work, re-running steps, clearing state                                       |

---

## Pipeline at a glance

```
AnnData (.h5ad)
    │
    ▼
[Preprocessing & QC]          QC filtering · normalization · HVG selection
    │
    ▼
[Dimensionality Reduction]     PCA · UMAP · Leiden clustering
    │
    ├──────────────────────────┬──────────────────────────┐
    ▼                          ▼                          │
[CellOracle GRN]          [Hotspot Modules]               │
 per-cluster TF networks   autocorrelated gene modules     │
    │                          │                          │
    └──────────────────────────┘                          │
    │                                                     │
    ▼                                                     │
[GRN Deep Analysis]           NetworkX plots, scoring     │
    │                                                     │
    ▼                                                     │
[Reporting]                   HTML · PDF report           │
                                                          │
    ◄─────────────────────────────────────────────────────┘
    Optional: stratify by cell type, run each in parallel
```

All stages are optional and can be run selectively using `--steps` or `--skip-*` flags.
