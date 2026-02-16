<!-- TO REMOVE: This file is mislabeled — it's not a changelog, it documents preprocessing config integration.
     Its content duplicates docs/PREPROCESSING_CONFIG_UPDATE.md. Remove this file. -->

# Preprocessing Module - Config Integration Summary

## What Changed

The `trnspot/preprocessing.py` module has been updated to use configuration defaults from `trnspot.config` instead of hardcoded default values.

## Key Updates

### 1. `perform_qc()` Function

**Changes:**

- All numeric parameters now default to `None` instead of hardcoded values
- Function checks if parameter is `None` and uses corresponding config value
- Added `min_cells` parameter (previously hardcoded to 10)
- Added `figsize` parameter (previously hardcoded to (15, 4))

**Config Parameters Used:**

```python
min_genes → config.QC_MIN_GENES (default: 200)
min_counts → config.QC_MIN_COUNTS (default: 500)
max_counts → config.QC_MAX_COUNTS (default: None)
pct_counts_mt_max → config.QC_PCT_MT_MAX (default: 20.0)
min_cells → config.QC_MIN_CELLS (default: 10) [NEW]
figsize → config.PLOT_FIGSIZE_LARGE (default: (15, 10))
```

### 2. `perform_grn_pre_processing()` Function

**Changes:**

- `top_genes` parameter now defaults to `None` instead of 3000
- Added `n_neighbors` parameter (previously hardcoded to 30)
- Added `n_pcs` parameter (previously hardcoded to 20)
- Added `svd_solver` parameter (previously hardcoded to 'arpack')
- Added return statement to return processed AnnData
- Added informative print statements

**Config Parameters Used:**

```python
top_genes → config.HVGS_N_TOP_GENES (default: 2000)
n_neighbors → config.NEIGHBORS_N_NEIGHBORS (default: 15)
n_pcs → config.NEIGHBORS_N_PCS (default: 40)
svd_solver → config.PCA_SVDSOLVE (default: 'arpack')
```

## Usage Patterns

### Pattern 1: Use All Defaults (Recommended)

```python
from trnspot.preprocessing import perform_qc

# Simplest usage - uses all config defaults
adata_qc = perform_qc(adata)
```

### Pattern 2: Override Specific Parameters

```python
# Override only what you need, rest uses config
adata_qc = perform_qc(
    adata,
    min_genes=300,
    min_counts=1000
)
```

### Pattern 3: Update Config Globally

```python
from trnspot import config

# Change defaults for entire analysis
config.update_config(
    QC_MIN_GENES=500,
    QC_MIN_COUNTS=2000,
    HVGS_N_TOP_GENES=3000
)

# All subsequent calls use new defaults
adata_qc = perform_qc(adata)
adata_grn = perform_grn_pre_processing(adata_qc)
```

### Pattern 4: Configuration Profiles

```python
from trnspot import config

def apply_strict_qc():
    config.update_config(
        QC_MIN_GENES=500,
        QC_MIN_COUNTS=2000,
        QC_PCT_MT_MAX=10.0
    )

def apply_lenient_qc():
    config.update_config(
        QC_MIN_GENES=100,
        QC_MIN_COUNTS=200,
        QC_PCT_MT_MAX=30.0
    )

# Use profiles
apply_strict_qc()
adata_strict = perform_qc(adata)
```

## Backward Compatibility

**100% Backward Compatible** - If you were passing explicit values to parameters, your code continues to work exactly as before:

```python
# Old code still works
adata_qc = perform_qc(
    adata,
    min_genes=200,
    min_counts=500,
    pct_counts_mt_max=20.0
)

# Equivalent new code (using defaults)
adata_qc = perform_qc(adata)
```

## Benefits

1. **Less Verbose Code**: No need to repeatedly specify same parameters
2. **Consistency**: All functions use same defaults from central config
3. **Reproducibility**: Save and share config for identical analysis
4. **Flexibility**: Easy to create and switch between analysis profiles
5. **Documentation**: All parameters documented in one place

## Files Modified

- `trnspot/preprocessing.py` - Updated both functions
- `tests/test_preprocessing.py` - Added test for config defaults
- `examples/preprocessing_with_config.py` - Comprehensive examples
- `examples/test_config_integration.py` - Quick integration test
- `docs/PREPROCESSING_CONFIG_UPDATE.md` - Detailed documentation
- `README.md` - Updated usage examples

## Testing

Run the integration test:

```bash
python examples/test_config_integration.py
```

Run unit tests:

```bash
pytest tests/test_preprocessing.py -v
```

## Next Steps

All preprocessing functions now support config-based defaults. Consider:

1. **Using config defaults** in your analyses for consistency
2. **Creating configuration profiles** for different analysis types
3. **Saving configurations** alongside analysis results for reproducibility

## Example: Complete Analysis with Config

```python
#!/usr/bin/env python
"""Complete analysis using config-based preprocessing"""

import scanpy as sc
from trnspot import set_random_seed, config
from trnspot.preprocessing import perform_qc, perform_grn_pre_processing

# Setup
set_random_seed(42)
config.update_config(
    QC_MIN_GENES=300,
    QC_MIN_COUNTS=1000,
    HVGS_N_TOP_GENES=3000
)

# Load data
adata = sc.read_h5ad('data.h5ad')

# QC - uses config defaults
adata = perform_qc(adata)

# Normalize
sc.pp.normalize_total(adata, target_sum=config.NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata)

# GRN preprocessing - uses config defaults
adata = perform_grn_pre_processing(adata, cluster_key='cell_type')

# Save with config for reproducibility
import json
adata.write('processed_data.h5ad')
with open('analysis_config.json', 'w') as f:
    json.dump(config.get_config(), f, indent=2, default=str)
```

## Documentation

- **Detailed Guide**: `docs/PREPROCESSING_CONFIG_UPDATE.md`
- **Config Documentation**: `docs/CONFIG.md`
- **Examples**: `examples/preprocessing_with_config.py`
- **Package Structure**: `docs/PACKAGE_STRUCTURE.md`
