---
layout: default
title: Getting Started
nav_order: 2
has_children: true
---

# Getting Started

Everything you need to install GeneCircuitry and run your first analysis.

| Page                             | What you'll find                                               |
| -------------------------------- | -------------------------------------------------------------- |
| [Installation](../installation/) | Four install methods, optional dependency groups, verification |
| [Quick Start](../quickstart/)    | CLI and Python API walkthroughs, expected output structure     |

---

## What is GeneCircuitry?

GeneCircuitry is a Python package for **transcriptional regulatory network (TRN) analysis** of single-cell RNA-seq data. It combines three established tools into a single checkpointed pipeline:

| Tool                                       | Role                                                                 |
| ------------------------------------------ | -------------------------------------------------------------------- |
| [Scanpy](https://scanpy.readthedocs.io/)   | Quality control, normalization, dimensionality reduction, clustering |
| [CellOracle](https://celloracle.org/)      | Gene regulatory network (GRN) inference per cluster                  |
| [Hotspot](https://hotspot.readthedocs.io/) | Identification of spatially autocorrelated gene modules              |

The pipeline is built around the `PipelineController` class, which manages checkpoints, logging, parallelism, and per-cluster stratification so you can run analyses reproducibly and resume interrupted jobs automatically.

---

## Minimum requirements

- Python ≥ 3.9
- An AnnData `.h5ad` file (or use the bundled PBMC 3k example)

CellOracle and Hotspot are **optional** — the core pipeline runs without them.
