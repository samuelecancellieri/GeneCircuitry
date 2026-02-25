# CellOracle Processing Module - Summary

## Overview

The `celloracle_processing.py` module provides a complete interface for gene regulatory network (GRN) inference using CellOracle. This document summarizes the module's functionality, documentation, and testing.

## Module Contents

### Functions

1. **`create_oracle_object()`**

   - Creates a CellOracle Oracle object from AnnData
   - Supports raw or normalized counts
   - Loads CellOracle's base GRN
   - Optionally integrates custom TF-target dictionaries

2. **`run_PCA()`**

   - Performs PCA on Oracle object
   - Automatically selects optimal number of components
   - Visualizes explained variance

3. **`run_KNN()`**

   - Performs KNN imputation for noise reduction
   - Auto-calculates optimal k (2.5% of cells)
   - Uses config for parallel processing

4. **`run_links()`**
   - Infers gene regulatory network
   - Filters links by statistical significance
   - Generates network visualizations
   - Calculates network quality metrics

### Configuration Integration

Uses the following config parameters:

- `config.GRN_N_JOBS` - Parallel processing (-1 = all cores)
- `config.FIGURES_DIR` - Figure output directory
- `config.OUTPUT_DIR` - Results output directory

## Documentation

### Main Documentation

**File:** `docs/CELLORACLE_PROCESSING.md`

**Contents:**

- Detailed function descriptions with all parameters
- Complete workflow examples
- Data requirements and preprocessing guidelines
- Tips and best practices
- Troubleshooting guide
- Configuration information

**Sections:**

1. Installation requirements
2. Function reference (4 functions)
3. Complete workflow example
4. Configuration parameters
5. Data requirements
6. Tips and best practices
7. Troubleshooting
8. References

### Example Scripts

**File:** `examples/celloracle_workflow.py`

**Features:**

- Complete end-to-end workflow
- 7-step process with detailed output
- Error handling and validation
- Result saving and configuration export
- Summary statistics and next steps

**Steps:**

1. Setup and configuration
2. Load and preprocess data
3. Create Oracle object
4. Perform PCA
5. KNN imputation
6. Infer regulatory links
7. Save all results

## Testing

### Test File

**File:** `tests/test_celloracle.py`

**Test Coverage:**

- Oracle object creation (raw/normalized counts)
- Oracle object with custom TF dictionary
- PCA computation
- KNN imputation
- Regulatory link inference
- Different p-value cutoffs
- Configuration integration
- Automatic k calculation
- Mock tests for environments without CellOracle

**Test Types:**

1. **Real tests** - Run when CellOracle is installed
2. **Mock tests** - Always run, test logic without CellOracle
3. **Integration tests** - Test config parameter usage
4. **Signature tests** - Verify function interfaces

**Total Tests:** 15+ tests covering all functionality

### Running Tests

```bash
# Run all CellOracle tests
pytest tests/test_celloracle.py -v

# Run only mock tests (no CellOracle needed)
pytest tests/test_celloracle.py -v -k "mock"

# Skip tests requiring CellOracle
pytest tests/test_celloracle.py -v -k "not skipif"
```

## Usage Patterns

### Pattern 1: Basic Workflow

```python
from genecircuitry.celloracle_processing import *

oracle = create_oracle_object(adata, 'leiden', 'X_umap')
oracle = run_PCA(oracle)
run_KNN(oracle, n_comps=50)
links = run_links(oracle, 'leiden', p_cutoff=0.001)
```

### Pattern 2: With Custom Dictionary

```python
oracle = create_oracle_object(
    adata,
    cluster_column_name='cell_type',
    embedding_name='X_umap',
    TG_to_TF_dictionary='path/to/dict.pkl'
)
```

### Pattern 3: Different P-value Cutoffs

```python
# Strict filtering
links_strict = run_links(oracle, 'leiden', p_cutoff=0.0001)

# Standard filtering
links_standard = run_links(oracle, 'leiden', p_cutoff=0.001)

# Lenient filtering
links_lenient = run_links(oracle, 'leiden', p_cutoff=0.01)
```

### Pattern 4: With Config Control

```python
from genecircuitry import config

config.update_config(
    GRN_N_JOBS=8,
    FIGURES_DIR='my_figures',
    OUTPUT_DIR='my_output'
)

# Functions automatically use updated config
oracle = create_oracle_object(adata, 'leiden', 'X_umap')
# ... continue workflow
```

## Data Requirements

### Minimum Requirements

1. AnnData object with gene expression data
2. Cluster annotations in `adata.obs`
3. Dimensionality reduction in `adata.obsm`
4. Gene names as symbols (for base GRN matching)

### Recommended Setup

```python
# Save raw counts
adata.layers['raw_counts'] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Clustering
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

# UMAP
sc.tl.umap(adata)

# Ready for CellOracle
```

## Key Features

### Strengths

1. **Easy to Use** - Simple function calls
2. **Automatic Selection** - PCA components and k automatically chosen
3. **Config Integration** - Uses centralized configuration
4. **Well Documented** - Comprehensive docs and examples
5. **Fully Tested** - 15+ tests with mocking support
6. **Error Handling** - Graceful failures with informative messages
7. **Visualization** - Automatic plot generation
8. **Flexible** - Supports raw or normalized counts

### Limitations

1. **Requires CellOracle** - Must install separately
2. **Memory Intensive** - Large datasets may need HPC
3. **Human Data** - Base GRN optimized for human genes
4. **Computation Time** - Can be slow for large datasets

## Output Files

### Generated Files

1. **Oracle Object** - `.celloracle.oracle` file

   - Contains processed expression data
   - PCA results
   - KNN imputation results

2. **Links Object** - `.celloracle.links` file

   - Inferred regulatory links
   - Statistical significance values
   - Network topology information

3. **Figures** - PNG plots

   - Explained variance (from PCA)
   - Degree distribution (from links)

4. **Configuration** - JSON file
   - Analysis parameters
   - For reproducibility

## Best Practices

### Before Starting

1. Ensure proper data preprocessing
2. Save raw counts before normalization
3. Perform quality control
4. Generate meaningful clusters
5. Check gene name format

### During Analysis

1. Monitor explained variance plot
2. Verify automatic k selection
3. Check network statistics
4. Examine degree distribution
5. Validate with known biology

### After Analysis

1. Save all objects (oracle, links, config)
2. Document parameters used
3. Perform validation analyses
4. Compare across conditions
5. Share configuration for reproducibility

## Common Issues and Solutions

### Issue: Gene Names Not Recognized

**Solution:** Convert to standard gene symbols

```python
# Example conversion
adata.var_names = [name.upper() for name in adata.var_names]
```

### Issue: Memory Error

**Solution:** Downsample cells or use HPC

```python
# Downsample
sc.pp.subsample(adata, n_obs=10000)
```

### Issue: Few Significant Links

**Solution:** Try more lenient p_cutoff

```python
links = run_links(oracle, 'leiden', p_cutoff=0.01)
```

### Issue: Slow Performance

**Solution:** Optimize config parameters

```python
config.update_config(GRN_N_JOBS=8)  # Use fewer cores
```

## Integration with GeneCircuitry

### Complete Pipeline

```python
# 1. Config and QC
from genecircuitry import set_random_seed, config
from genecircuitry.preprocessing import perform_qc

set_random_seed(42)
adata = perform_qc(adata)

# 2. Standard preprocessing
sc.pp.normalize_total(adata, target_sum=config.NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata)

# 3. CellOracle GRN inference
from genecircuitry.celloracle_processing import *

oracle = create_oracle_object(adata, 'leiden', 'X_umap')
oracle = run_PCA(oracle)
run_KNN(oracle, n_comps=50)
links = run_links(oracle, 'leiden')

# 4. Save everything
import json
oracle.to_hdf5(f"{config.OUTPUT_DIR}/oracle.celloracle.oracle")
links.to_hdf5(f"{config.OUTPUT_DIR}/links.celloracle.links")
with open(f"{config.OUTPUT_DIR}/config.json", 'w') as f:
    json.dump(config.get_config(), f, indent=2, default=str)
```

## Future Enhancements

Potential improvements:

1. Support for mouse genes
2. Custom base GRN integration
3. Batch processing utilities
4. Network comparison tools
5. Automated validation metrics
6. Interactive visualizations
7. Perturbation simulation wrappers
8. Multi-condition analysis helpers

## Related Documentation

- Main README: `README.md`
- Config guide: `docs/CONFIG.md`
- QC functions: `docs/QC_FUNCTIONS.md`
- Package structure: `docs/PACKAGE_STRUCTURE.md`
- Preprocessing updates: `docs/PREPROCESSING_CONFIG_UPDATE.md`

## References

1. CellOracle: https://morris-lab.github.io/CellOracle.documentation/
2. Kamimoto et al., Nature, 2023 - Original CellOracle paper
3. GeneCircuitry documentation: `docs/` directory

## Summary Statistics

- **Functions:** 4 main functions
- **Documentation:** 350+ lines
- **Tests:** 15+ test cases
- **Examples:** 1 complete workflow script
- **Config Parameters:** 3 integrated parameters
- **Lines of Code:** ~150 in module

## Status

✅ **Complete** - Module is fully functional, documented, and tested.

**Ready for:**

- Production use
- Further development
- Integration into larger workflows
- Publication-ready analyses
