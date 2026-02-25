# CellOracle Processing - Quick Reference

## Installation

```bash
pip install celloracle
```

## Quick Start

```python
from genecircuitry.celloracle_processing import *

# 1. Create Oracle
oracle = create_oracle_object(adata, 'leiden', 'X_umap')

# 2. Run PCA
oracle = run_PCA(oracle)

# 3. KNN imputation
run_KNN(oracle, n_comps=50)

# 4. Infer links
links = run_links(oracle, 'leiden', p_cutoff=0.001)

# 5. Save
oracle.to_hdf5('oracle.celloracle.oracle')
links.to_hdf5('links.celloracle.links')
```

## Function Reference

### create_oracle_object()

```python
oracle = create_oracle_object(
    adata=adata,                    # AnnData object
    cluster_column_name='leiden',   # Cluster key in adata.obs
    embedding_name='X_umap',        # Embedding in adata.obsm
    TG_to_TF_dictionary=None,       # Optional custom dict
    raw_count_layer='raw_counts'    # Layer name or None
)
```

### run_PCA()

```python
oracle = run_PCA(oracle)  # Auto-selects n_comps
```

### run_KNN()

```python
run_KNN(
    oracle,
    n_comps=50  # Number of PCs to use
)
```

### run_links()

```python
links = run_links(
    oracle,
    cluster_column_name='leiden',
    p_cutoff=0.001  # Significance threshold
)
```

## Data Preparation Checklist

- [ ] Save raw counts: `adata.layers['raw_counts'] = adata.X.copy()`
- [ ] Normalize: `sc.pp.normalize_total()` + `sc.pp.log1p()`
- [ ] Cluster: `sc.tl.leiden()` or `sc.tl.louvain()`
- [ ] UMAP: `sc.tl.umap()`
- [ ] Gene symbols: Ensure genes are in symbol format

## Common P-value Cutoffs

- `0.0001` - Very strict (high confidence only)
- `0.001` - Standard (recommended)
- `0.01` - Lenient (more links)

## Config Parameters

```python
from genecircuitry import config

config.update_config(
    GRN_N_JOBS=-1,           # Parallel jobs
    FIGURES_DIR='figures',   # Output directory
    OUTPUT_DIR='output'      # Save directory
)
```

## Troubleshooting

| Issue                 | Solution                     |
| --------------------- | ---------------------------- |
| Gene names not found  | Convert to gene symbols      |
| Memory error          | Downsample cells or use HPC  |
| Few significant links | Increase p_cutoff to 0.01    |
| Slow performance      | Reduce n_comps or GRN_N_JOBS |

## Output Files

- `oracle.celloracle.oracle` - Oracle object
- `links.celloracle.links` - Regulatory links
- `grn_degree_distribution.png` - Network plot

## Full Documentation

- Comprehensive guide: `docs/CELLORACLE_PROCESSING.md`
- Module summary: `docs/CELLORACLE_MODULE_SUMMARY.md`
- Example workflow: `examples/celloracle_workflow.py`
- Tests: `tests/test_celloracle.py`
