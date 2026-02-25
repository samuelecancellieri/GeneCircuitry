# GeneCircuitry Complete Pipeline Usage Guide

## Overview

The complete pipeline script (`run_complete_analysis.py` / `examples/complete_pipeline.py`) provides a single command to run the entire GeneCircuitry analysis workflow from raw data through GRN inference and module identification.

## Quick Start

### Basic Usage

```bash
# Run with example dataset
python run_complete_analysis.py

# Run with your own data
python run_complete_analysis.py --input your_data.h5ad
```

### Common Use Cases

#### 1. Full Analysis with Custom Data

```bash
python run_complete_analysis.py \
    --input data/my_experiment.h5ad \
    --output results/experiment1 \
    --figures plots/experiment1 \
    --seed 42 \
    --n-jobs 16
```

#### 2. Skip Specific Analyses

```bash
# Only preprocessing and clustering (skip GRN and modules)
python run_complete_analysis.py --skip-celloracle --skip-hotspot

# Only GRN analysis (skip module identification)
python run_complete_analysis.py --skip-hotspot

# Only module analysis (skip GRN)
python run_complete_analysis.py --skip-celloracle
```

#### 3. Custom QC Parameters

```bash
python run_complete_analysis.py \
    --min-genes 300 \
    --min-counts 1000 \
    --input your_data.h5ad
```

#### 4. Pre-processed Data

```bash
# If QC already performed
python run_complete_analysis.py \
    --input preprocessed_data.h5ad \
    --skip-qc
```

## Command-Line Options

### Input/Output

- `--input PATH` or `-i PATH`: Input data file (.h5ad or .h5)
  - If not provided, uses Paul et al. 2015 example dataset
- `--output DIR` or `-o DIR`: Output directory (default: `output`)
- `--figures DIR` or `-f DIR`: Figures directory (default: `figures`)

### Analysis Control

- `--skip-qc`: Skip quality control step (use if already performed)
- `--skip-celloracle`: Skip CellOracle GRN inference
- `--skip-hotspot`: Skip Hotspot module identification

### Configuration

- `--seed N`: Random seed for reproducibility (default: 42)
- `--n-jobs N`: Number of parallel jobs (default: 8)
- `--min-genes N`: Minimum genes per cell for QC (default: 200)
- `--min-counts N`: Minimum counts per cell for QC (default: 500)

### Help

- `--help` or `-h`: Show help message and exit

## Pipeline Steps

### Step 0: Configuration Setup

- Sets random seed for reproducibility
- Configures analysis parameters
- Creates output directories

### Step 1: Data Loading

- Loads data from .h5ad or .h5 file
- Falls back to example dataset if no input provided
- Reports dataset dimensions

### Step 2: Quality Control and Preprocessing

- Performs cell and gene filtering
- Calculates QC metrics (MT%, counts, genes)
- Normalizes and log-transforms data
- Stores raw counts for downstream analysis

### Step 3: Dimensionality Reduction and Clustering

- Selects highly variable genes
- Performs PCA
- Computes diffusion map and PAGA
- Generates UMAP embedding
- Performs Leiden clustering
- Saves preprocessed data

### Step 4: CellOracle GRN Inference (Optional)

- Creates Oracle object with base GRN
- Runs PCA on Oracle
- Performs KNN imputation
- Infers gene regulatory links
- Saves Oracle object and links

### Step 5: Hotspot Module Identification (Optional)

- Creates Hotspot object
- Computes gene autocorrelations
- Identifies spatially correlated genes
- Discovers gene modules
- Saves results in multiple formats

### Step 6: Summary Report

- Generates comprehensive analysis summary
- Lists all output files
- Provides analysis statistics
- Suggests next steps

## Output Files

After running the pipeline, you'll find:

### Main Output Directory

```
output/
├── preprocessed_adata.h5ad        # Preprocessed AnnData object
├── analysis_summary.txt           # Text summary of analysis
├── celloracle/                    # CellOracle results
│   ├── oracle_object.celloracle.oracle
│   └── grn_links.celloracle.links
└── hotspot/                       # Hotspot results
    ├── hotspot_result.pkl         # Module scores (pickle)
    ├── autocorrelation_results.csv
    ├── significant_genes.csv
    └── gene_modules.csv
```

### Figures Directory

```
figures/
├── qc/                            # QC plots
│   ├── qc_before_filtering.png
│   └── qc_after_filtering.png
├── grn_analysis/                  # GRN visualizations
│   └── grn_link_count.png
└── hotspot_local_correlations.png
```

## Error Handling

The pipeline gracefully handles:

- **Missing optional dependencies**: CellOracle and Hotspot are optional
- **Data format issues**: Checks file format before loading
- **Invalid parameters**: Validates inputs before processing
- **Runtime errors**: Catches and reports errors without stopping entire pipeline

## Tips and Best Practices

### 1. Start Small

```bash
# Test with example dataset first
python run_complete_analysis.py --skip-celloracle --skip-hotspot
```

### 2. Monitor Resource Usage

```bash
# Adjust parallel jobs based on available cores
python run_complete_analysis.py --n-jobs 8  # For 8-core machine
```

### 3. Modular Execution

```bash
# Run preprocessing first, then add analyses later
python run_complete_analysis.py --skip-celloracle --skip-hotspot
# Later, load preprocessed data for specific analyses
```

### 4. Check Dependencies

Before running CellOracle or Hotspot:

```bash
# Install CellOracle
pip install celloracle

# Install Hotspot
pip install hotspot-sc
```

### 5. Review Summary

Always check `output/analysis_summary.txt` for:

- Dataset statistics
- Analysis completion status
- Output file locations
- Suggested next steps

## Troubleshooting

### "ModuleNotFoundError: No module named 'celloracle'"

Solution: Either install CellOracle (`pip install celloracle`) or skip with `--skip-celloracle`

### "ModuleNotFoundError: No module named 'hotspot'"

Solution: Either install Hotspot (`pip install hotspot-sc`) or skip with `--skip-hotspot`

### "Memory Error"

Solution: Reduce dataset size or number of parallel jobs:

```bash
python run_complete_analysis.py --n-jobs 4
```

### Pipeline Fails Mid-Run

- Check `output/analysis_summary.txt` to see which steps completed
- Use skip flags to avoid re-running completed steps
- Review error messages for specific issues

## Integration with Existing Workflows

### Loading Results in Python

```python
import scanpy as sc
import pickle

# Load preprocessed data
adata = sc.read_h5ad('output/preprocessed_adata.h5ad')

# Load CellOracle results
from celloracle import Oracle
oracle = Oracle()
oracle.load_hdf5('output/celloracle/oracle_object.celloracle.oracle')

# Load Hotspot results
with open('output/hotspot/hotspot_result.pkl', 'rb') as f:
    module_scores = pickle.load(f)
```

### Continuing Analysis

```python
# Load and continue with additional analysis
adata = sc.read_h5ad('output/preprocessed_adata.h5ad')

# Differential expression
sc.tl.rank_genes_groups(adata, 'leiden')

# Additional visualizations
sc.pl.umap(adata, color='leiden')
```

## Advanced Usage

### Custom Configuration

Create a Python script to run with custom settings:

```python
from genecircuitry import config
config.update_config(
    HOTSPOT_TOP_GENES=5000,
    HOTSPOT_FDR_THRESHOLD=0.01,
    GRN_N_JOBS=16
)

# Then run pipeline functions...
```

### Batch Processing

```bash
#!/bin/bash
for file in data/*.h5ad; do
    name=$(basename "$file" .h5ad)
    python run_complete_analysis.py \
        --input "$file" \
        --output "results/$name" \
        --figures "plots/$name"
done
```

## See Also

- Individual workflow examples:

  - `examples/celloracle_workflow.py` - CellOracle-specific workflow
  - `examples/hotspot_workflow.py` - Hotspot-specific workflow
  - `examples/preprocess_workflow.py` - Preprocessing only

- Documentation:
  - `docs/CONFIG.md` - Configuration options
  - `docs/CELLORACLE_PROCESSING.md` - CellOracle documentation
  - `README.md` - Package overview
