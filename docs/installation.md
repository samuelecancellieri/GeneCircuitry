---
layout: default
title: Installation
nav_order: 1
parent: Getting Started
---

# Installation

{: .no_toc }

## Table of contents

{: .no_toc .text-delta }

1. TOC
   {:toc}

---

## Requirements

genecircuitry requires **Python ≥ 3.9, < 3.11**. Most dependencies are available
on [conda-forge](https://conda-forge.org/) and [bioconda](https://bioconda.github.io/).
Two optional analysis engines — [CellOracle](https://github.com/morris-lab/CellOracle)
and [hotspotsc](https://github.com/YosefLab/Hotspot) — are **only available via pip**
and must be installed as a separate step after the conda environment is set up.

---

## Option 1 — Pixi (recommended)

[Pixi](https://prefix.dev/docs/pixi/overview) manages conda and pip dependencies
together in a single reproducible environment. It is the easiest and cleanest way
to get a fully working installation.

```bash
# Install pixi (one-time, see https://prefix.dev/docs/pixi/installation)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone the repository
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd genecircuitry

# Create the environment and install all dependencies (conda + pip) in one step
pixi install

# Run the pipeline inside the pixi environment
pixi run run

# Run with the bundled test dataset
pixi run genecircuitry

# Or drop into an interactive shell
pixi shell
```

{: .note }

> **Developer environment** (adds pytest, black, flake8, mypy):
>
> ```bash
> pixi install -e dev
> pixi run -e dev test
> ```

---

## Option 2 — Conda

Install genecircuitry and its conda-available dependencies from
[bioconda](https://bioconda.github.io/) and [conda-forge](https://conda-forge.org/),
then install the pip-only dependencies manually.

```bash
# 1. Create a fresh environment (Python 3.9 is recommended)
conda create -n genecircuitry python=3.9
conda activate genecircuitry

# 2. Install genecircuitry and all conda-available dependencies
conda install -c bioconda -c conda-forge genecircuitry

# 3. Install the pip-only optional analysis engines
#    (CellOracle for GRN inference, hotspotsc for gene modules)
pip install celloracle==0.18.0 hotspotsc==1.1.3
```

{: .note }

> Skip step 3 if you only need preprocessing/QC and do not require GRN inference
> or gene module analysis.

---

## Option 3 — pip / venv

```bash
# Clone the repository
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd genecircuitry

# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install the package with all optional dependencies
pip install -e ".[grn,hotspot]"

# Or install core only (no CellOracle / hotspotsc)
pip install -e .

# Install with development dependencies
pip install -e ".[dev]"
```

---

## Option 4 — Docker

A pre-built image is available that ships all dependencies (including CellOracle
and hotspotsc) and works out of the box:

```bash
# Pull and run (bind-mount your data and output directories)
docker run --rm \
    -v /path/to/your/data:/data \
    -v /path/to/output:/output \
    zanathos/genecircuitry:latest \
    --input /data/your_data.h5ad --output /output

# Check available options
docker run --rm zanathos/genecircuitry:latest --help
```

Build the image locally from source:

```bash
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd genecircuitry
docker build -t genecircuitry .
docker run --rm genecircuitry --help
```

---

## Option 5 — Singularity / Apptainer

On HPC clusters where Docker is not available, use [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) or its drop-in replacement [Apptainer](https://apptainer.org/) to convert the pre-built Docker image into a `.sif` container:

```bash
# Pull the Docker image and convert to a Singularity image file
singularity pull genecircuitry.sif docker://zanathos/genecircuitry:latest

# Equivalent with Apptainer
apptainer pull genecircuitry.sif docker://zanathos/genecircuitry:latest
```

Run the analysis, binding your data and output directories into the container:

```bash
singularity exec \
    --bind /path/to/your/data:/data \
    --bind /path/to/output:/output \
    genecircuitry.sif \
    genecircuitry --input /data/your_data.h5ad --output /output

# Check available options
singularity exec genecircuitry.sif genecircuitry --help
```

For GPU nodes or when wrapping environment variables, use `--nv` and `--env`:

```bash
singularity exec \
    --nv \
    --bind /scratch/data:/data \
    --bind /scratch/results:/output \
    genecircuitry.sif \
    genecircuitry --input /data/cells.h5ad --output /output --n-jobs 16
```

Build the `.sif` locally from the Dockerfile (requires Docker and root/`fakeroot`):

```bash
# Build Docker image first, then export to Singularity
docker build -t genecircuitry .
docker save genecircuitry | singularity build genecircuitry.sif docker-archive:/dev/stdin
```

{: .note }

> Singularity/Apptainer containers run as the calling user by default, so output files
> are always written with your own permissions — no `chown` step needed.

---

## Optional dependency groups

| Extra        | Packages                            | Description                  |
| ------------ | ----------------------------------- | ---------------------------- |
| `grn`        | `celloracle==0.18.0`                | GRN inference via CellOracle |
| `hotspot`    | `hotspotsc==1.1.3`                  | Gene module identification   |
| `enrichment` | `gseapy>=1.0.0`                     | Gene set enrichment analysis |
| `atac`       | `genomepy`, `gimmemotifs`           | ATAC-seq peak motif analysis |
| `dev`        | `pytest`, `black`, `flake8`, `mypy` | Development tools            |

{: .note }

> **CellOracle and Hotspot are optional.** `import genecircuitry` will succeed even when
> they are not installed. The relevant modules (`genecircuitry.celloracle_processing`,
> `genecircuitry.hotspot_processing`) will be `None` at runtime and gracefully skipped
> by the pipeline.

---

## Verify installation

```python
import genecircuitry
print(genecircuitry.__version__)           # e.g. 0.1.7

# Check optional deps
print(genecircuitry.celloracle_processing) # None if not installed
print(genecircuitry.hotspot_processing)    # None if not installed
```

---

## Example data

The repository ships with `data/pbmc3k_raw.h5ad` (PBMC 3k dataset) and
`data/paul15/paul15.h5` (Paul et al. 2015, mouse hematopoiesis) for testing.
The pipeline also falls back to `sc.datasets.pbmc3k()` if no input file is provided.
