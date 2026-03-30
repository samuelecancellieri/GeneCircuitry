---
layout: default
title: Checkpoints & Resume
nav_order: 8
parent: User Guide
---

# Checkpoints & Resume

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

GeneCircuitry automatically saves **checkpoint files** after each pipeline step completes. On the next run with the same `--output` directory, completed steps are detected and skipped — letting you resume interrupted analyses without recomputing expensive steps.

---

## How checkpoints work

After a step finishes successfully, the `PipelineController` writes a `.checkpoint` file to `<output>/logs/`:

```
results/logs/
├── pipeline.log
├── error.log
├── load.checkpoint
├── preprocessing.checkpoint
├── clustering.checkpoint
├── celloracle.checkpoint
└── hotspot.checkpoint
```

On the next run, before executing a step, the controller checks whether the corresponding `.checkpoint` file exists. If it does, the step is logged as skipped and execution moves to the next step.

For [stratified analyses](../stratification/), each cluster has its own checkpoint namespace:

```
results/stratified_analysis/<ClusterName>/logs/
├── pipeline.log
├── clustering.checkpoint
├── celloracle.checkpoint
└── hotspot.checkpoint
```

---

## Resume an interrupted run

Simply re-run the same command — no flags needed:

```bash
# First run: interrupted during celloracle step
python -m genecircuitry.pipeline --input data.h5ad --output results/

# Second run: automatically skips load, preprocessing, clustering
python -m genecircuitry.pipeline --input data.h5ad --output results/
# → Pipeline resumes from celloracle
```

For stratified analyses, per-cluster checkpoints work independently:

```bash
# One cluster failed — re-run the stratified pipeline
python -m genecircuitry.pipeline \
    --input data.h5ad --output results/ \
    --cluster-key-stratification cell_type
# → Only the failed/incomplete clusters are re-run
```

---

## Running specific steps

Use `--steps` to run only selected steps, ignoring checkpoints:

```bash
# Explicitly run only celloracle and report (even if checkpointed)
python -m genecircuitry.pipeline \
    --input results/clustered_adata.h5ad \
    --output results/ \
    --steps celloracle report
```

{: .note }

> When `--steps` is provided, only the listed steps are evaluated for checkpoints. Steps not in the list are skipped entirely (regardless of checkpoint state).

---

## Re-running a specific step

To force a step to re-run, delete its checkpoint file:

```bash
# Re-run clustering with updated resolution
rm results/logs/clustering.checkpoint
python -m genecircuitry.pipeline --input data.h5ad --output results/
```

Or delete all checkpoints to restart from scratch:

```bash
rm results/logs/*.checkpoint
python -m genecircuitry.pipeline --input data.h5ad --output results/
```

---

## Checkpoint files

Checkpoint files are plain text files containing a JSON payload with metadata about the completed step (timestamp, input hash, key metrics). They are not large and can be inspected:

```bash
cat results/logs/preprocessing.checkpoint
# {"step": "preprocessing", "timestamp": "2026-03-30T14:22:11", "n_cells": 4521, "n_genes": 18302}
```

{: .warning }

> Do not move checkpoint files to a different `--output` directory — they are path-relative. If you copy results to a new location, delete the checkpoint files first to avoid stale state.

---

## Using a new output directory

To start a fresh run (regardless of any previous checkpoints) while preserving old results:

```bash
python -m genecircuitry.pipeline \
    --input data.h5ad \
    --output results/run_v2/    # new directory → no checkpoints exist → full run
```

---

## Programmatic checkpoint control

The `PipelineController` exposes checkpoint helpers used internally:

```python
from genecircuitry.pipeline.controller import PipelineController

controller = PipelineController(args, start_time)

# Check if a step is checkpointed
if controller.is_step_checkpointed("celloracle"):
    print("CellOracle already done")

# Mark a step as done (write checkpoint)
controller.write_checkpoint("my_custom_step", metadata={"n_links": 1234})

# Clear a checkpoint
controller.clear_checkpoint("celloracle")
```
