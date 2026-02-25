# GeneCircuitry Package Structure

## Overview

GeneCircuitry is a Python package for transcriptional regulatory network analysis with built-in quality control and preprocessing capabilities.

## Directory Structure

```
GeneCircuitry_simple/
├── venv/                         # Virtual environment
├── genecircuitry/                      # Main package
│   ├── __init__.py              # Package initialization
│   ├── config.py                # Configuration module (NEW!)
│   ├── preprocessing.py         # Data preprocessing & QC
│   ├── grn_analysis.py          # GRN analysis
│   └── celloracle_processing.py # CellOracle integration
├── tests/                        # Test suite
│   ├── __init__.py
│   ├── test_config.py           # Config tests (NEW!)
│   └── test_preprocessing.py    # Preprocessing tests
├── examples/                     # Example scripts
│   ├── config_example.py        # Config usage example (NEW!)
│   ├── example_qc.py            # QC workflow example
│   └── quick_demo.py            # Quick demo
├── docs/                         # Documentation
│   ├── CONFIG.md                # Config documentation (NEW!)
│   └── QC_FUNCTIONS.md          # QC documentation
├── setup.py                      # Package setup
├── pyproject.toml               # Modern Python config
├── requirements.txt             # Dependencies
├── requirements-dev.txt         # Dev dependencies
├── README.md                    # Main documentation
├── .gitignore                   # Git ignore rules
└── MANIFEST.in                  # Package manifest
```

## Key Components

### 1. Configuration Module (`config.py`)

**Purpose**: Centralized configuration management

**Features**:

- Random seed management for reproducibility
- Default parameters for QC, preprocessing, plotting
- Easy parameter updates
- Configuration profiles

**Key Functions**:

- `set_random_seed(seed)`: Set global random seed
- `get_config()`: Get all config as dictionary
- `print_config()`: Display all parameters
- `update_config(**kwargs)`: Update specific parameters

**Parameters Categories**:

- Random seed (RANDOM_SEED=42)
- QC parameters (min_genes, min_counts, pct_mt_max)
- Plotting settings (DPI, figsize, colors)
- Preprocessing (normalization, HVG selection, PCA)
- Neighbors & clustering (k-NN, resolution)
- UMAP (min_dist, spread)
- GRN inference (jobs, thresholds)
- File I/O (output directories)
- Advanced (jobs, memory, logging)

### 2. Preprocessing Module (`preprocessing.py`)

**Purpose**: Data quality control and preprocessing

**Key Functions**:

- `perform_qc()`: Complete QC pipeline with filtering
- `plot_qc_violin()`: Violin plots for QC metrics
- `plot_qc_scatter()`: Scatter plots showing metric relationships

**Features**:

- Automatic mitochondrial gene detection
- Cell/gene filtering
- Before/after comparison plots
- Multiple visualization options

### 3. Example Scripts

- `config_example.py`: Demonstrates configuration usage
- `example_qc.py`: Complete QC workflow
- `quick_demo.py`: Quick interactive demo

### 4. Documentation

- `CONFIG.md`: Comprehensive config documentation
- `QC_FUNCTIONS.md`: QC functions documentation
- `README.md`: Package overview and quick start

### 5. Tests

- `test_config.py`: Tests for configuration module
- `test_preprocessing.py`: Tests for QC functions

## Usage Workflow

### Basic Workflow

```python
# 1. Import and configure
from genecircuitry import set_random_seed, config
from genecircuitry.preprocessing import perform_qc
import scanpy as sc

# 2. Set reproducibility
set_random_seed(42)

# 3. Load data
adata = sc.read_h5ad('data.h5ad')

# 4. Perform QC
adata_qc = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_MT_MAX,
    plot=True
)

# 5. Continue analysis...
```

### Custom Configuration Workflow

```python
from genecircuitry import config

# Update for specific analysis
config.update_config(
    RANDOM_SEED=123,
    QC_MIN_GENES=500,
    QC_MIN_COUNTS=1000,
    QC_PCT_MT_MAX=15.0
)

# Use updated config
adata_qc = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_MT_MAX
)
```

## Development Commands

```bash
# Activate environment
source venv/bin/activate

# Install package in development mode
pip install -e ".[dev]"

# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_config.py -v

# Run with coverage
pytest tests/ --cov=genecircuitry

# Format code
black genecircuitry/

# Check linting
flake8 genecircuitry/

# Type checking
mypy genecircuitry/
```

## Next Steps

1. **Add More Preprocessing Functions**:

   - Normalization
   - Highly variable genes selection
   - PCA/dimensionality reduction
   - Batch correction

2. **Implement GRN Analysis**:

   - Network inference
   - TF-target identification
   - Network visualization

3. **CellOracle Integration**:

   - Oracle construction
   - Perturbation simulation
   - Vector field analysis

4. **Enhanced Documentation**:

   - API reference
   - Tutorials
   - Best practices guide

5. **Additional Features**:
   - Pipeline orchestration
   - Report generation
   - Interactive visualizations

## Package Design Principles

1. **Reproducibility**: Random seed management
2. **Configurability**: Centralized parameters
3. **Modularity**: Separate concerns
4. **Documentation**: Comprehensive docs
5. **Testing**: Unit tests for all modules
6. **Usability**: Simple, intuitive API
7. **Extensibility**: Easy to add new features
