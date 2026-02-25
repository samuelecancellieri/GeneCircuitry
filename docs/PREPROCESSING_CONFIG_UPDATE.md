# Preprocessing Module - Config Integration

## Overview

The preprocessing module has been updated to use configuration defaults from `genecircuitry.config`. This provides consistent parameter management and improves reproducibility.

## Updated Functions

### `perform_qc()`

**Before:**

```python
def perform_qc(
    adata: AnnData,
    min_genes: int = 200,
    min_counts: int = 500,
    max_counts: Optional[int] = None,
    pct_counts_mt_max: float = 20.0,
    ...
) -> AnnData:
```

**After:**

```python
def perform_qc(
    adata: AnnData,
    min_genes: Optional[int] = None,  # Uses config.QC_MIN_GENES
    min_counts: Optional[int] = None,  # Uses config.QC_MIN_COUNTS
    max_counts: Optional[int] = None,  # Uses config.QC_MAX_COUNTS
    pct_counts_mt_max: Optional[float] = None,  # Uses config.QC_PCT_MT_MAX
    min_cells: Optional[int] = None,  # Uses config.QC_MIN_CELLS (NEW!)
    figsize: Optional[Tuple[int, int]] = None,  # Uses config.PLOT_FIGSIZE_LARGE
    ...
) -> AnnData:
```

**Config Parameters Used:**

- `config.QC_MIN_GENES` (default: 200)
- `config.QC_MIN_COUNTS` (default: 500)
- `config.QC_MAX_COUNTS` (default: None)
- `config.QC_PCT_MT_MAX` (default: 20.0)
- `config.QC_MIN_CELLS` (default: 10)
- `config.PLOT_FIGSIZE_LARGE` (default: (15, 10))

### `perform_grn_pre_processing()`

**Before:**

```python
def perform_grn_pre_processing(
    adata: AnnData,
    cluster_key: Optional[str] = None,
    cell_downsample: int = 20000,
    top_genes: int = 3000,
) -> AnnData:
```

**After:**

```python
def perform_grn_pre_processing(
    adata: AnnData,
    cluster_key: Optional[str] = None,
    cell_downsample: int = 20000,
    top_genes: Optional[int] = None,  # Uses config.HVGS_N_TOP_GENES
    n_neighbors: Optional[int] = None,  # Uses config.NEIGHBORS_N_NEIGHBORS (NEW!)
    n_pcs: Optional[int] = None,  # Uses config.NEIGHBORS_N_PCS (NEW!)
    svd_solver: Optional[str] = None,  # Uses config.PCA_SVDSOLVE (NEW!)
) -> AnnData:
```

**Config Parameters Used:**

- `config.HVGS_N_TOP_GENES` (default: 2000)
- `config.NEIGHBORS_N_NEIGHBORS` (default: 15)
- `config.NEIGHBORS_N_PCS` (default: 40)
- `config.PCA_SVDSOLVE` (default: 'arpack')

**Additional Changes:**

- Added return statement to return preprocessed AnnData
- Added print statements for better logging
- Uses config parameters throughout the function

## Usage Examples

### Example 1: Use All Defaults

```python
from genecircuitry.preprocessing import perform_qc
import scanpy as sc

adata = sc.read_h5ad('data.h5ad')

# Uses all config defaults
adata_qc = perform_qc(adata)
# min_genes=200, min_counts=500, pct_counts_mt_max=20.0, min_cells=10
```

### Example 2: Override Specific Parameters

```python
# Override only specific parameters, others use config defaults
adata_qc = perform_qc(
    adata,
    min_genes=300,  # Override
    min_counts=1000  # Override
    # pct_counts_mt_max uses config default (20.0)
    # min_cells uses config default (10)
)
```

### Example 3: Update Config Globally

```python
from genecircuitry import config

# Update config for all subsequent calls
config.update_config(
    QC_MIN_GENES=500,
    QC_MIN_COUNTS=2000,
    QC_PCT_MT_MAX=15.0
)

# Now uses updated defaults
adata_qc = perform_qc(adata)
# min_genes=500, min_counts=2000, pct_counts_mt_max=15.0
```

### Example 4: GRN Preprocessing with Config

```python
from genecircuitry.preprocessing import perform_grn_pre_processing

# Uses all config defaults
adata_grn = perform_grn_pre_processing(adata_qc)
# top_genes=2000, n_neighbors=15, n_pcs=40, svd_solver='arpack'

# Override specific parameters
adata_grn = perform_grn_pre_processing(
    adata_qc,
    top_genes=5000,  # Override
    n_pcs=50  # Override
    # n_neighbors uses config default (15)
    # svd_solver uses config default ('arpack')
)
```

### Example 5: Configuration Profiles

```python
from genecircuitry import config

def set_strict_qc_profile():
    """Profile for strict quality control"""
    config.update_config(
        QC_MIN_GENES=500,
        QC_MIN_COUNTS=2000,
        QC_PCT_MT_MAX=10.0,
        QC_MIN_CELLS=20
    )

def set_lenient_qc_profile():
    """Profile for lenient quality control"""
    config.update_config(
        QC_MIN_GENES=100,
        QC_MIN_COUNTS=200,
        QC_PCT_MT_MAX=30.0,
        QC_MIN_CELLS=3
    )

# Apply a profile
set_strict_qc_profile()
adata_strict = perform_qc(adata)

# Switch to different profile
set_lenient_qc_profile()
adata_lenient = perform_qc(adata)
```

## Benefits

### 1. Consistency

All functions use the same default values from a central location.

### 2. Reproducibility

Configuration can be saved and shared to ensure identical analysis parameters.

```python
import json
from genecircuitry import config

# Save configuration
with open('analysis_config.json', 'w') as f:
    json.dump(config.get_config(), f, indent=2, default=str)

# Load and apply configuration
with open('analysis_config.json', 'r') as f:
    cfg = json.load(f)
    config.update_config(**{k: v for k, v in cfg.items() if k in dir(config)})
```

### 3. Less Code Duplication

No need to repeatedly specify the same parameters.

### 4. Easy Customization

Change defaults globally without modifying function calls.

### 5. Better Documentation

Config parameters are documented in one place.

## Migration Guide

### Old Code

```python
from genecircuitry.preprocessing import perform_qc

# Had to specify all parameters
adata_qc = perform_qc(
    adata,
    min_genes=200,
    min_counts=500,
    pct_counts_mt_max=20.0
)
```

### New Code - Option 1: Use Defaults

```python
from genecircuitry.preprocessing import perform_qc

# Simply use defaults
adata_qc = perform_qc(adata)
```

### New Code - Option 2: Update Config

```python
from genecircuitry import config
from genecircuitry.preprocessing import perform_qc

# Set defaults once
config.update_config(
    QC_MIN_GENES=200,
    QC_MIN_COUNTS=500,
    QC_PCT_MT_MAX=20.0
)

# Use throughout analysis
adata_qc = perform_qc(adata)
```

### New Code - Option 3: Override Specific

```python
from genecircuitry.preprocessing import perform_qc

# Override only what you need
adata_qc = perform_qc(
    adata,
    min_genes=300  # Override, others use config defaults
)
```

## Breaking Changes

**None** - The changes are backward compatible. If you were passing explicit parameter values, your code will continue to work exactly as before.

The only difference is that if you were relying on the hardcoded defaults and now want different values, you can either:

1. Pass them explicitly (as before)
2. Update the config (new option)

## Additional Improvements

1. **Added `min_cells` parameter** to `perform_qc()` - now configurable via `config.QC_MIN_CELLS`
2. **Added return statement** to `perform_grn_pre_processing()` - function now properly returns the processed AnnData
3. **Added logging** - Better print statements throughout preprocessing functions
4. **Consistent parameter names** - All parameters follow config naming conventions

## Testing

Updated tests include:

- `test_perform_qc_with_config_defaults()` - Tests using config defaults
- Updated existing tests to include `min_cells` parameter
- All tests pass with new parameter handling

## See Also

- `docs/CONFIG.md` - Complete configuration documentation
- `examples/preprocessing_with_config.py` - Comprehensive examples
- `genecircuitry/config.py` - Configuration module source
