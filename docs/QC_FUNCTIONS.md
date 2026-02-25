# Quality Control Functions

## Overview

The `preprocessing.py` module provides comprehensive quality control (QC) functions for single-cell RNA-seq data analysis. These functions help identify and filter low-quality cells based on standard QC metrics.

## Functions

### `perform_qc()`

Main function to perform quality control on AnnData objects.

**Parameters:**

- `adata` (AnnData): Input annotated data matrix
- `min_genes` (int, default=200): Minimum number of genes per cell
- `min_counts` (int, default=500): Minimum total counts per cell
- `max_counts` (int, optional): Maximum total counts per cell (helps identify doublets)
- `pct_counts_mt_max` (float, default=20.0): Maximum percentage of mitochondrial reads
- `plot` (bool, default=True): Whether to generate QC plots
- `figsize` (tuple, default=(15, 4)): Figure size for plots
- `save_plots` (str, optional): Path to save QC plots

**Returns:**

- Filtered AnnData object with QC metrics in `.obs` and `.var`

**QC Metrics Calculated:**

- `n_genes_by_counts`: Number of genes detected per cell
- `total_counts`: Total UMI counts per cell
- `pct_counts_mt`: Percentage of mitochondrial reads
- `mt`: Boolean flag for mitochondrial genes

**Example:**

```python
import scanpy as sc
from genecircuitry.preprocessing import perform_qc

# Load data
adata = sc.read_h5ad('data.h5ad')

# Perform QC
adata_qc = perform_qc(
    adata,
    min_genes=300,
    min_counts=1000,
    max_counts=50000,
    pct_counts_mt_max=15.0,
    plot=True,
    save_plots='qc_results.png'
)
```

### `plot_qc_violin()`

Generate violin plots for QC metrics, optionally grouped by a categorical variable.

**Parameters:**

- `adata` (AnnData): AnnData object with QC metrics calculated
- `groupby` (str, optional): Column in `adata.obs` to group by (e.g., 'sample', 'batch')
- `figsize` (tuple, default=(12, 4)): Figure size
- `save_path` (str, optional): Path to save the figure

**Example:**

```python
# Violin plots for all cells
plot_qc_violin(adata_qc)

# Violin plots grouped by sample
plot_qc_violin(adata_qc, groupby='sample', save_path='violin_by_sample.png')
```

### `plot_qc_scatter()`

Generate scatter plots to visualize relationships between QC metrics.

**Parameters:**

- `adata` (AnnData): AnnData object with QC metrics calculated
- `color_by` (str, optional): Metric to use for coloring points (default='pct_counts_mt')
- `figsize` (tuple, default=(12, 4)): Figure size
- `save_path` (str, optional): Path to save the figure

**Plots Generated:**

1. Total counts vs Number of genes
2. Total counts vs Mitochondrial percentage
3. Number of genes vs Mitochondrial percentage

**Example:**

```python
# Scatter plots colored by MT percentage
plot_qc_scatter(adata_qc, color_by='pct_counts_mt')

# Scatter plots colored by batch
plot_qc_scatter(adata_qc, color_by='batch', save_path='scatter_qc.png')
```

## Typical Workflow

```python
import scanpy as sc
from genecircuitry.preprocessing import perform_qc, plot_qc_violin, plot_qc_scatter

# 1. Load data
adata = sc.read_h5ad('raw_data.h5ad')
print(f"Initial shape: {adata.shape}")

# 2. Perform QC with filtering
adata_qc = perform_qc(
    adata,
    min_genes=200,
    min_counts=500,
    max_counts=30000,
    pct_counts_mt_max=20.0,
    plot=True,
    save_plots='qc_histograms.png'
)

# 3. Generate additional QC visualizations
plot_qc_violin(
    adata_qc,
    groupby='sample',
    save_path='qc_violin.png'
)

plot_qc_scatter(
    adata_qc,
    color_by='pct_counts_mt',
    save_path='qc_scatter.png'
)

# 4. Save filtered data
adata_qc.write('filtered_data.h5ad')
```

## Interpreting QC Plots

### Histogram Plots (from `perform_qc`)

- **Red bars**: Cells before filtering
- **Blue bars**: Cells after filtering
- Shows distribution of:
  - Genes per cell
  - Counts per cell
  - Mitochondrial percentage

### Violin Plots

- Shows distribution of QC metrics across all cells or by group
- Wider sections indicate more cells with those values
- Horizontal lines show median and mean

### Scatter Plots

- Reveals relationships between metrics
- Can help identify:
  - Doublets (high counts + high genes)
  - Low-quality cells (low counts + high MT%)
  - Batch effects (different patterns by group)

## Quality Control Guidelines

### Typical Thresholds

- **min_genes**: 200-500 (cells with fewer genes are likely empty droplets)
- **min_counts**: 500-1000 (low counts indicate poor capture)
- **max_counts**: 20000-50000 (very high counts may indicate doublets)
- **pct_counts_mt_max**: 10-20% (high MT% indicates dying/stressed cells)

### Adjusting Thresholds

- Visualize distributions before setting thresholds
- Consider tissue type (brain cells typically have lower MT%)
- Account for experimental protocol (10x, Smart-seq2, etc.)
- Balance between removing low-quality cells and retaining biological variation

## Notes

- The function creates a copy of the input data, so the original is not modified
- Mitochondrial genes are detected by the "MT-" prefix (human genes)
- For mouse data, use "mt-" prefix or modify the code accordingly
- All QC metrics are stored in `adata.obs` for downstream analysis
