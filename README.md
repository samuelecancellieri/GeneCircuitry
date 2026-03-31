# genecircuitry

A Python package for transcriptional regulatory network analysis.

## Installation

genecircuitry requires **Python >=3.9, <3.11**. Most dependencies are available
on [conda-forge](https://conda-forge.org/) and [bioconda](https://bioconda.github.io/).
Two optional analysis engines — [CellOracle](https://github.com/morris-lab/CellOracle)
and [hotspotsc](https://github.com/YosefLab/Hotspot) — are **only available via pip**
and must be installed as a separate step after the conda environment is set up.

---

### Option 1 — Pixi (recommended)

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

# Or drop into an interactive shell
pixi shell
```

> **Developer environment** (adds pytest, black, flake8, mypy):
>
> ```bash
> pixi install -e dev
> pixi run -e dev test
> ```

---

### Option 2 — Conda

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

> Skip step 3 if you only need preprocessing/QC and do not require GRN inference
> or gene module analysis.

---

### Option 3 — pip / venv

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

### Option 4 — Docker

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

## Quick Start

### Complete Analysis Pipeline

Run the full analysis pipeline from preprocessing through CellOracle and Hotspot:

```bash
# Run with example dataset (default)
python run_complete_analysis.py

# Or use the script directly
python examples/complete_pipeline.py

# Run with your own data
python run_complete_analysis.py --input your_data.h5ad --output results

# Skip specific analyses
python run_complete_analysis.py --skip-celloracle  # Skip GRN inference
python run_complete_analysis.py --skip-hotspot      # Skip module identification

# Custom parameters
python run_complete_analysis.py --seed 123 --n-jobs 16 --min-genes 300

# See all options
python run_complete_analysis.py --help
```

### Modular Execution & Parallel Processing

NEW: The pipeline now supports modular execution and parallel processing:

```bash
# Run only specific steps
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps load preprocessing clustering

# Stratified analysis in parallel (multiple cell types/clusters)
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 4

# Resume from checkpoints
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps celloracle hotspot  # Skips preprocessing if checkpoint exists
```

**Available step names:** `load`, `preprocessing`, `stratification`, `clustering`, `celloracle`, `hotspot`, `grn_analysis`, `summary`

**Parallel benefits:**

- Process multiple stratifications simultaneously
- Linear speedup with number of workers
- Automatic checkpoint integration
- See [Controller Guide](docs/CONTROLLER_QUICK_REF.md) for details

The complete pipeline includes:

1. **Data Loading** - Load h5ad/h5 files or use example dataset
2. **Quality Control** - Cell and gene filtering with QC metrics
3. **Preprocessing** - Normalization, HVG selection, PCA, clustering
4. **CellOracle GRN Inference** - Gene regulatory network prediction
5. **Hotspot Module Analysis** - Spatially autocorrelated gene modules
6. **Summary Report** - Comprehensive analysis summary with output files

Output structure:

```
output/
├── preprocessed_adata.h5ad        # Preprocessed dataset
├── analysis_summary.txt           # Analysis report
├── celloracle/
│   ├── oracle_object.celloracle.oracle
│   └── grn_links.celloracle.links
└── hotspot/
    ├── autocorrelation_results.csv
    ├── significant_genes.csv
    └── gene_modules.csv
figures/
├── qc/                            # QC plots
└── grn_analysis/                  # GRN visualizations
```

## Usage

### Configuration and Reproducibility

Set random seed and configure default parameters:

```python
from genecircuitry import set_random_seed, config

# Set random seed for reproducibility
set_random_seed(42)

# View all configuration parameters
config.print_config()

# Update specific parameters
config.update_config(
    QC_MIN_GENES=300,
    QC_MIN_COUNTS=1000,
    PLOT_DPI=600
)
```

### Quality Control

Perform comprehensive quality control on single-cell RNA-seq data:

```python
import scanpy as sc
from genecircuitry import config
from genecircuitry.preprocessing import perform_qc, plot_qc_violin, plot_qc_scatter

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Perform QC - uses config defaults automatically
adata_qc = perform_qc(adata)
# Equivalent to: min_genes=200, min_counts=500, pct_counts_mt_max=20.0, min_cells=10

# Or override specific parameters
adata_qc = perform_qc(
    adata,
    min_genes=300,      # Override
    min_counts=1000     # Override
    # Other params use config defaults
)
```

### Complete Workflow Example

```python
import scanpy as sc
from genecircuitry import set_random_seed, config
from genecircuitry.preprocessing import perform_qc, perform_grn_pre_processing

# 1. Set up reproducibility
set_random_seed(42)

# 2. Optionally customize config
config.update_config(QC_MIN_GENES=300, QC_MIN_COUNTS=1000)

# 3. Load and process data
adata = sc.read_h5ad('your_data.h5ad')
adata = perform_qc(adata)  # Uses config defaults

# 4. Normalize
sc.pp.normalize_total(adata, target_sum=config.NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata)

# 5. GRN preprocessing
adata = perform_grn_pre_processing(adata)  # Uses config defaults
```

### CellOracle Integration

Perform gene regulatory network inference using CellOracle:

```python
from genecircuitry.celloracle_processing import (
    create_oracle_object,
    run_PCA,
    run_KNN,
    run_links
)

# Note: Requires CellOracle installation
# pip install celloracle

# Create Oracle object
oracle = create_oracle_object(
    adata=adata,
    cluster_column_name='leiden',
    embedding_name='X_umap',
    raw_count_layer='raw_counts'
)

# Perform PCA and KNN imputation
oracle = run_PCA(oracle)
run_KNN(oracle, n_comps=50)

# Infer regulatory links
links = run_links(
    oracle,
    cluster_column_name='leiden',
    p_cutoff=0.001
)

# Save results
oracle.to_hdf5('oracle_object.celloracle.oracle')
links.to_hdf5('grn_links.celloracle.links')
```

### Running Examples

```bash
# Activate environment
source venv/bin/activate

# Run QC example
python examples/example_qc.py

# Run configuration example
python examples/config_example.py

# Run CellOracle workflow (requires celloracle)
python examples/celloracle_workflow.py

# Run quick demo
python examples/quick_demo.py

# Test config integration
python examples/test_config_integration.py
```

## Features

- **Configuration Management**: Centralized configuration for reproducibility
  - Global random seed setting
  - Default parameters for all analysis steps
  - Easy parameter updates
  - Configuration profiles for different analysis types
- **Quality Control**: Comprehensive QC with multiple visualization options
  - Cell filtering based on gene count, total counts, and mitochondrial percentage
  - Automated QC metrics calculation
  - Before/after filtering comparison plots
  - Violin and scatter plots for detailed inspection
- **Data Preprocessing**: Complete preprocessing pipeline
  - Normalization and scaling for single-cell RNA-seq data
  - Highly variable genes selection
  - PCA and dimensionality reduction
  - Neighborhood graph construction
- **CellOracle Integration**: Gene regulatory network inference
  - Oracle object creation with raw or normalized counts
  - Automated PCA component selection
  - KNN imputation for noise reduction
  - Regulatory link inference with statistical filtering
  - Network visualization and quality metrics
- **Gene Regulatory Network Analysis**: Network construction and analysis tools
  - TF-target gene relationship inference
  - Network topology analysis
  - Cluster-specific GRN construction

## Documentation

- **[Configuration Guide](docs/CONFIG.md)** - Complete configuration documentation
- **[QC Functions](docs/QC_FUNCTIONS.md)** - Quality control functions guide
- **[CellOracle Processing](docs/CELLORACLE_PROCESSING.md)** - CellOracle integration guide
- **[Package Structure](docs/PACKAGE_STRUCTURE.md)** - Package organization overview
- **[Preprocessing Updates](docs/PREPROCESSING_CONFIG_UPDATE.md)** - Config integration details

## Development

### Running Tests

```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_config.py -v
pytest tests/test_preprocessing.py -v
pytest tests/test_celloracle.py -v

# Run with coverage
pytest tests/ --cov=genecircuitry --cov-report=html
```

### Code Quality

```bash
# Format code
black genecircuitry/

# Check linting
flake8 genecircuitry/

# Type checking
mypy genecircuitry/
```

### Testing Notes

- CellOracle tests use mocking when CellOracle is not installed
- Some tests are skipped if CellOracle is not available
- Use `pytest -v` for verbose output
- Tests cover all major functionality with both unit and integration tests

## License

MIT License

## Authors

Samuele Cancellieri
