# GeneCircuitry

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9 | 3.10](https://img.shields.io/badge/python-3.9%20%7C%203.10-blue.svg)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-brightgreen.svg)](https://samuelecancellieri.github.io/GeneCircuitry/)

**GeneCircuitry** is a modular Python framework for transcriptional regulatory network (TRN) inference and spatial/phenotypic gene module discovery from single-cell data. It seamlessly integrates **[Scanpy](https://scanpy.readthedocs.io/)**, **[CellOracle](https://celloracle.org/)**, and **[Hotspot](https://hotspot.readthedocs.io/)** into a checkpoint-enabled, parallelized pipeline.

---

## Key Features

- ⚡ **Modular & Checkpointed Execution**: Run full workflows or selective steps (`--steps load preprocessing clustering`) with automatic checkpoint resumption.
- 🚀 **Parallel Stratified Analysis**: Process multiple cell types or experimental conditions simultaneously across CPU workers with linear speedup.
- 🔀 **Multi-Key Grouping**: Native support for compound keys (e.g., `--cluster-key-stratification cell_type,condition` or `["cell_type", "condition"]`) with filesystem-safe sanitization.
- 📝 **Dual-Stream Logging**: Automated structured step tracking in `pipeline.log` and contextual stack traces in `error.log`.
- 📊 **Publication-Ready Reports**: Generates interactive HTML and PDF summaries with embedded quality metrics, GRN graphs, and module heatmaps.

---

## Installation

GeneCircuitry requires **Python >=3.9, <3.11**.

### Option 1: Pixi (Recommended)
[Pixi](https://prefix.dev/) installs both Conda and PyPI dependencies into a reproducible environment in a single command:

```bash
# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and install
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd genecircuitry
pixi install

# Run the pipeline or enter interactive shell
pixi run run
pixi shell
```

### Option 2: pip / venv
```bash
git clone https://github.com/samuelecancellieri/genecircuitry.git
cd genecircuitry
python -m venv venv && source venv/bin/activate

# Install with optional analysis engines (CellOracle + Hotspot)
pip install -e ".[grn,hotspot]"
```

### Option 3: Conda
```bash
conda create -n genecircuitry python=3.9 -y && conda activate genecircuitry
conda install -c bioconda -c conda-forge genecircuitry
pip install celloracle==0.18.0 hotspotsc==1.1.3
```

### Option 4: Docker
```bash
docker run --rm \
    -v /path/to/data:/data \
    -v /path/to/output:/output \
    zanathos/genecircuitry:latest \
    --input /data/input.h5ad --output /output
```

---

## Quick Start

### 1. Command-Line Interface (CLI)

```bash
# Run full analysis on demo PBMC3k dataset
python run_complete_analysis.py

# Run on your dataset
python run_complete_analysis.py --input data.h5ad --output results/

# Run selective steps only (auto-resumes from previous checkpoints)
python -m genecircuitry.pipeline \
    --input data.h5ad --output results/ \
    --steps load preprocessing clustering

# Multi-key parallel stratified analysis across 4 cores
python -m genecircuitry.pipeline \
    --input data.h5ad --output results/ \
    --cluster-key-stratification cell_type,condition \
    --parallel --n-jobs 4

# Target specific composite subgroups
python -m genecircuitry.pipeline \
    --input data.h5ad --output results/ \
    --cluster-key-stratification cell_type,condition \
    --clusters CD4_Control,CD8_Treated
```

**Available modular steps**: `load`, `preprocessing`, `stratification`, `clustering`, `celloracle`, `hotspot`, `grn_analysis`, `summary`.

---

### 2. Python API

```python
from datetime import datetime
from argparse import Namespace
from genecircuitry import set_random_seed, config
from genecircuitry.pipeline.controller import PipelineController

# Set reproducibility seed
set_random_seed(42)

# Configure parameters
config.update_config(QC_MIN_GENES=300, LEIDEN_RESOLUTION=0.8)

# Initialize controller
args = Namespace(
    input="data.h5ad",
    output="results/",
    cluster_key="cell_type",
    cluster_key_stratification=None,
    force_dim_reduction=False,
    steps=None,
    parallel=False,
    n_jobs=4,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_COUNTS_MT_MAX,
    min_cells=config.QC_MIN_CELLS,
    downsample_cells=config.DOWNSAMPLE_CELLS,
    cell_downsample=config.GRN_CELL_DOWNSAMPLE,
    alpha=config.GRN_ALPHA,
    p_cutoff=config.GRN_P_CUTOFF,
    hotspot_model=config.HOTSPOT_MODEL,
    hotspot_top_genes=config.HOTSPOT_TOP_GENES,
    report_title="Analysis Report",
    report_formats=["html", "pdf"],
)

controller = PipelineController(args, datetime.now())
controller.run_complete_pipeline()
```

---

## Comparative & Aggregation Analysis

GeneCircuitry features an integrated comparative analysis engine that bridges **Hotspot co-expression modules** and **CellOracle regulatory networks** across cell types (clusters) and experimental conditions (stratifications):

- 🧬 **Gene-Overlap Module Alignment**: Resolves arbitrary sequential module IDs (`Module 1`, `Module 2` across runs) by computing pairwise Jaccard similarity matrices ($J = \frac{|A \cap B|}{|A \cup B|}$) to align truly homologous biological programs.
- 🎯 **TF-to-Module Regulatory Concordance**: Cross-references module expression activity with TF target coverage to identify which **Master Driver TF** actively regulates each module in each cell state.
- 🌊 **4-Stage Biological Sankey Flow**: Visualizes end-to-end gene provenance: *Hotspot Significance $\rightarrow$ Co-expression Module $\rightarrow$ GRN Master Regulator $\rightarrow$ Cell-State Activation $\rightarrow$ Enriched Pathway*.
- 🔄 **Differential TF Rewiring**: Identifies whether transcription factors regulate a static set of core targets or dynamically rewire their downstream networks across conditions.
- 📊 **Integrated 6-Panel Dashboard**: Synthesizes module activities, TF centralities, provenance breakdowns, and top regulatory circuits into a single publication-ready figure.

For full mathematical formulations and interpretation guides, see the [Comparative Analysis Guide](docs/comparative-analysis.md).

---

## Output Structure

Each run produces an organized directory structure:

```
output/
├── preprocessed_adata.h5ad             # Filtered, normalized, and clustered AnnData
├── report.html                         # Interactive HTML report
├── report.pdf                          # Compiled PDF summary
├── logs/
│   ├── pipeline.log                    # Execution trace and step metrics
│   ├── error.log                       # Stack traces for failures
│   └── *.checkpoint                    # JSON state hashes for smart resume
├── figures/
│   ├── qc/                             # QC violin and scatter plots
│   ├── grn/                            # Regulatory network graphs & centrality plots
│   ├── hotspot/                        # Module heatmaps and expression violins
│   └── comparative/                    # Concordance bubbles, alignment heatmaps, Sankey flow
├── comparative/
│   ├── module_activity_matrix.csv      # Mean module expression scores
│   ├── cross_stratification_module_jaccard.csv # Pairwise module Jaccard similarity
│   ├── cross_stratification_module_alignment.csv # Homologous module matching table
│   ├── tf_to_module_matrix.csv         # TF-to-module target count matrix
│   ├── tf_specificity_summary.csv      # Master vs. group-specific TF classification
│   ├── differential_tf_targets.csv     # Conserved vs. rewired TF targets
│   ├── gene_selection_provenance.csv   # 4-stage gene lineage tracing table
│   └── tf_module_concordance.csv       # TF driver and module concordance metrics
├── celloracle/
│   ├── oracle_object.celloracle.oracle # CellOracle object
│   ├── grn_links.celloracle.links      # Inferred regulatory links
│   └── grn_merged_scores.csv           # Edge weights and scores
├── hotspot/
│   ├── autocorrelation_results.csv     # Autocorrelation p-values & z-scores
│   └── gene_modules.csv                # Gene module cluster assignments
└── stratified_analysis/                # Subgroup outputs (when stratified)
    └── <Subgroup_Name>/
```

---

## Documentation

Full documentation is hosted at **[samuelecancellieri.github.io/GeneCircuitry](https://samuelecancellieri.github.io/GeneCircuitry/)**:

- 📖 **[Getting Started](https://samuelecancellieri.github.io/GeneCircuitry/installation/)** — Installation methods and troubleshooting
- ⚡ **[Quick Start](https://samuelecancellieri.github.io/GeneCircuitry/quickstart/)** — Step-by-step tutorial
- 🛠️ **[Pipeline Overview](https://samuelecancellieri.github.io/GeneCircuitry/pipeline/)** — Modular stages and CLI flags
- ⚙️ **[Configuration](https://samuelecancellieri.github.io/GeneCircuitry/configuration/)** — Full parameters reference
- 🔬 **[GRN Inference](https://samuelecancellieri.github.io/GeneCircuitry/celloracle/)** — CellOracle workflows
- 🧩 **[Gene Modules](https://samuelecancellieri.github.io/GeneCircuitry/hotspot/)** — Hotspot workflows
- 📊 **[Comparative Analysis](https://samuelecancellieri.github.io/GeneCircuitry/comparative-analysis/)** — Cross-cluster and cross-stratification integration
- 🎨 **[Plotting System](https://samuelecancellieri.github.io/GeneCircuitry/plotting/)** — Complete visualization catalog
- 🏗️ **[Implementation & Architecture](https://samuelecancellieri.github.io/GeneCircuitry/implementation/)** — Internal design, multiprocessing, logging, and APIs

---

## Development & Testing

```bash
# Install developer environment
pixi install -e dev

# Run test suite
pixi run -e dev test

# Linting and formatting
pixi run -e dev lint
pixi run -e dev format
pixi run -e dev typecheck
```

---

## License

This project is licensed under the [MIT License](LICENSE).

## Author

**Samuele Cancellieri** ([GitHub](https://github.com/samuelecancellieri))
