# Configuration Module

## Overview

The `config.py` module provides centralized configuration management for GeneCircuitry, including default parameters for quality control, preprocessing, plotting, and reproducibility settings.

## Key Features

- **Reproducibility**: Set and manage random seeds across all libraries
- **Default Parameters**: Centralized defaults for all analysis steps
- **Easy Customization**: Simple functions to update configuration
- **Type Safety**: All parameters are documented with their expected types

## Quick Start

```python
import genecircuitry
from genecircuitry import config

# Print all configuration parameters
config.print_config()

# Set random seed for reproducibility
genecircuitry.set_random_seed(42)

# Access configuration values
min_genes = config.QC_MIN_GENES
print(f"Default min genes: {min_genes}")

# Update configuration
config.update_config(RANDOM_SEED=123, QC_MIN_GENES=300)
```

## Configuration Categories

### 1. Random Seed

**Purpose**: Ensure reproducibility across analyses

```python
from genecircuitry import set_random_seed

# Set seed for numpy, random, and scanpy
set_random_seed(42)
```

**Parameters**:

- `RANDOM_SEED` (int, default=42): Global random seed

### 2. Quality Control Parameters

Default thresholds for cell and gene filtering:

```python
from genecircuitry import config

# QC parameters
QC_MIN_GENES = 200          # Min genes per cell
QC_MIN_COUNTS = 500         # Min counts per cell
QC_MAX_COUNTS = None        # Max counts per cell
QC_PCT_MT_MAX = 20.0        # Max mitochondrial %
QC_MIN_CELLS = 1            # Min cells per gene
```

**Usage**:

```python
from genecircuitry.preprocessing import perform_qc

adata_qc = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_MT_MAX
)
```

### 3. Plotting Configuration

Default settings for figure generation:

```python
PLOT_DPI = 300                      # Resolution
PLOT_FORMAT = 'png'                  # File format
PLOT_FIGSIZE_SMALL = (6, 4)         # Small figures
PLOT_FIGSIZE_MEDIUM = (12, 4)       # Medium figures
PLOT_FIGSIZE_LARGE = (15, 10)       # Large figures
PLOT_COLOR_PALETTE = 'viridis'      # Color scheme
```

### 4. Preprocessing Parameters

Normalization and feature selection defaults:

```python
NORMALIZE_TARGET_SUM = 1e4      # Target sum (10,000)
HVGS_N_TOP_GENES = 2000         # Number of HVGs
HVGS_MIN_MEAN = 0.0125          # Min mean for HVG
HVGS_MAX_MEAN = 3               # Max mean for HVG
HVGS_MIN_DISP = 0.5             # Min dispersion
PCA_N_COMPS = 50                # PCA components
```

### 5. Neighborhood Graph Parameters

Settings for computing cell-cell relationships:

```python
NEIGHBORS_N_NEIGHBORS = 15      # Number of neighbors
NEIGHBORS_N_PCS = 40            # PCs to use
NEIGHBORS_METHOD = 'umap'       # Method (umap/gauss)
NEIGHBORS_METRIC = 'euclidean'  # Distance metric
```

### 6. Clustering Parameters

Default resolution parameters:

```python
LEIDEN_RESOLUTION = 1.0         # Leiden resolution
LOUVAIN_RESOLUTION = 1.0        # Louvain resolution
```

### 7. UMAP Parameters

Dimensionality reduction settings:

```python
UMAP_MIN_DIST = 0.5             # Min distance
UMAP_SPREAD = 1.0               # Spread
UMAP_N_COMPONENTS = 2           # Dimensions
```

### 8. Gene Regulatory Network Parameters

GRN inference settings:

```python
GRN_N_JOBS = -1                      # Parallel jobs
GRN_MIN_TARGETS = 10                 # Min target genes
GRN_CONFIDENCE_THRESHOLD = 0.5       # Min confidence
```

### 9. File I/O Configuration

Directory paths for outputs:

```python
OUTPUT_DIR = 'output'           # Main output directory
CACHE_DIR = '.cache'            # Cache directory
FIGURES_DIR = 'figures'         # Figures directory
```

### 10. Advanced Settings

Performance and logging:

```python
N_JOBS = -1                     # Parallel jobs
CHUNK_SIZE = 1000               # Batch size
LOW_MEMORY = False              # Memory mode
LOG_LEVEL = 'INFO'              # Logging level
VERBOSE = True                  # Verbose output
```

## Functions

### `set_random_seed(seed)`

Set random seed across all libraries for reproducibility.

```python
from genecircuitry import set_random_seed

set_random_seed(42)
# Random seed set to: 42
```

### `get_config()`

Get all configuration parameters as a dictionary.

```python
from genecircuitry import get_config

config_dict = get_config()
print(config_dict['RANDOM_SEED'])  # 42
```

### `print_config()`

Print all configuration parameters in a formatted table.

```python
from genecircuitry import print_config

print_config()
# ============================================================
# GeneCircuitry Configuration
# ============================================================
# RANDOM_SEED                    = 42
# QC_MIN_GENES                   = 200
# ...
```

### `update_config(**kwargs)`

Update one or more configuration parameters.

```python
from genecircuitry import config

config.update_config(
    RANDOM_SEED=123,
    QC_MIN_GENES=300,
    PLOT_DPI=600
)
# Updated RANDOM_SEED = 123
# Updated QC_MIN_GENES = 300
# Updated PLOT_DPI = 600
```

## Usage Patterns

### Pattern 1: Use Default Configuration

```python
from genecircuitry import config
from genecircuitry.preprocessing import perform_qc

# Use defaults directly
adata_qc = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_MT_MAX
)
```

### Pattern 2: Customize for Specific Analysis

```python
from genecircuitry import config
from genecircuitry.preprocessing import perform_qc

# Update config for stricter QC
config.update_config(
    QC_MIN_GENES=500,
    QC_MIN_COUNTS=1000,
    QC_PCT_MT_MAX=10.0
)

# Apply customized thresholds
adata_qc = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_MT_MAX
)
```

### Pattern 3: Reproducible Analysis

```python
from genecircuitry import set_random_seed, config

# Set seed at start of analysis
set_random_seed(42)

# All subsequent random operations will be reproducible
# ... your analysis code ...
```

### Pattern 4: Configuration Profiles

```python
from genecircuitry import config

def set_strict_qc():
    """Profile for strict QC"""
    config.update_config(
        QC_MIN_GENES=500,
        QC_MIN_COUNTS=1000,
        QC_PCT_MT_MAX=10.0
    )

def set_lenient_qc():
    """Profile for lenient QC"""
    config.update_config(
        QC_MIN_GENES=100,
        QC_MIN_COUNTS=200,
        QC_PCT_MT_MAX=30.0
    )

# Use profiles
set_strict_qc()
# ... run analysis ...
```

## Best Practices

1. **Set Random Seed Early**: Call `set_random_seed()` at the beginning of your script
2. **Document Changes**: When updating config, document why you changed defaults
3. **Use Config for Defaults**: Reference config values rather than hardcoding numbers
4. **Create Profiles**: For common analysis types, create configuration profiles
5. **Version Control**: Include configuration parameters in analysis reports

## Integration with Analysis Pipeline

```python
#!/usr/bin/env python
"""
Complete analysis pipeline with configuration
"""

from genecircuitry import set_random_seed, config
from genecircuitry.preprocessing import perform_qc
import scanpy as sc

# 1. Set reproducibility
set_random_seed(42)
config.print_config()

# 2. Load data
adata = sc.read_h5ad('data.h5ad')

# 3. QC with config defaults
adata = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_MT_MAX,
    plot=True,
    save_plots=f"{config.FIGURES_DIR}/qc.png"
)

# 4. Further analysis...
```

## Notes

- Configuration is loaded and random seed is set automatically on import
- Changes to configuration are global and affect all subsequent operations
- For parallel execution, consider the `N_JOBS` parameter
- The config module uses Python's module-level variables for simplicity
