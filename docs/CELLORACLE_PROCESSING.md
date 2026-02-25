# CellOracle Processing Module

## Overview

The `celloracle_processing.py` module provides integration with CellOracle for gene regulatory network (GRN) analysis and cell fate simulation. It includes functions for creating Oracle objects, performing dimensionality reduction, KNN imputation, and inferring regulatory links.

## Installation Requirements

CellOracle must be installed separately:

```bash
pip install celloracle
```

## Functions

### `create_oracle_object()`

Create an Oracle object for CellOracle analysis.

**Signature:**

```python
def create_oracle_object(
    adata: AnnData,
    cluster_column_name: str,
    embedding_name: str,
    TG_to_TF_dictionary: Optional[str] = None,
    raw_count_layer: Optional[str] = "raw_counts",
) -> co.Oracle
```

**Parameters:**

- `adata` (AnnData): Annotated data object containing gene expression data
- `cluster_column_name` (str): Name of the column in `adata.obs` containing cluster information (e.g., 'leiden', 'louvain', 'cell_type')
- `embedding_name` (str): Name of the embedding to use for analysis (e.g., 'X_umap', 'X_pca')
- `TG_to_TF_dictionary` (str, optional): Path to pickle file containing dictionary mapping target genes to transcription factors. If None, uses CellOracle's base GRN only
- `raw_count_layer` (str, optional): Name of the layer in `adata.layers` containing raw count data. If None, uses normalized counts from `adata.X`

**Returns:**

- `oracle` (co.Oracle): Configured Oracle object ready for analysis

**Description:**

This function initializes a CellOracle Oracle object with the provided single-cell data. It:

1. Loads the human promoter base GRN from CellOracle
2. Imports the AnnData object (either raw or normalized counts)
3. Adds TF information from the base GRN
4. Optionally adds custom TF-target gene relationships

**Example:**

```python
import scanpy as sc
from genecircuitry.celloracle_processing import create_oracle_object

# Load preprocessed data
adata = sc.read_h5ad('preprocessed_data.h5ad')

# Create Oracle object with raw counts
oracle = create_oracle_object(
    adata=adata,
    cluster_column_name='leiden',
    embedding_name='X_umap',
    raw_count_layer='raw_counts'
)

# Or with normalized counts
oracle = create_oracle_object(
    adata=adata,
    cluster_column_name='cell_type',
    embedding_name='X_umap',
    raw_count_layer=None  # Uses adata.X
)

# With custom TF-target dictionary
oracle = create_oracle_object(
    adata=adata,
    cluster_column_name='leiden',
    embedding_name='X_umap',
    TG_to_TF_dictionary='path/to/TG_TF_dict.pkl'
)
```

### `run_PCA()`

Perform PCA on the Oracle object and automatically select optimal number of components.

**Signature:**

```python
def run_PCA(oracle: co.Oracle) -> co.Oracle
```

**Parameters:**

- `oracle` (co.Oracle): Oracle object to perform PCA on

**Returns:**

- `oracle_cc` (co.Oracle): Copy of Oracle object with PCA performed

**Description:**

This function:

1. Creates a copy of the input Oracle object
2. Performs PCA on the gene expression data
3. Automatically selects the optimal number of principal components based on the cumulative explained variance
4. Visualizes the explained variance and selected cutoff
5. Caps the maximum number of components at 50

The automatic selection finds the "elbow" point where the rate of change in cumulative variance levels off.

**Example:**

```python
from genecircuitry.celloracle_processing import create_oracle_object, run_PCA

# Create and process Oracle object
oracle = create_oracle_object(adata, 'leiden', 'X_umap')
oracle_pca = run_PCA(oracle)

# The plot shows cumulative variance and selected n_comps
# Continue with downstream analysis using oracle_pca
```

**Notes:**

- Generates a matplotlib plot showing explained variance
- Automatically determines optimal number of components (max 50)
- Returns a copy, leaving the original object unchanged

### `run_KNN()`

Perform KNN imputation on the Oracle object to smooth gene expression data.

**Signature:**

```python
def run_KNN(
    oracle: co.Oracle,
    n_comps: int = 50
) -> None
```

**Parameters:**

- `oracle` (co.Oracle): Oracle object to perform KNN imputation on
- `n_comps` (int, default=50): Number of PCA components to use for KNN. Should match the value from `run_PCA()`

**Returns:**

- None (modifies Oracle object in place)

**Description:**

This function performs KNN-based imputation to smooth gene expression data. It:

1. Creates a copy of the Oracle object
2. Automatically calculates optimal k (number of neighbors) as 2.5% of total cells
3. Performs balanced KNN imputation with optimized sight and max length parameters
4. Uses parallel processing based on `config.GRN_N_JOBS`

The imputation helps reduce noise and improve the quality of GRN inference.

**Example:**

```python
from genecircuitry.celloracle_processing import (
    create_oracle_object,
    run_PCA,
    run_KNN
)

# Complete preprocessing pipeline
oracle = create_oracle_object(adata, 'leiden', 'X_umap')
oracle = run_PCA(oracle)

# Get number of components from PCA
n_comps = oracle.pca.n_components_

# Run KNN imputation
run_KNN(oracle, n_comps=n_comps)

# Oracle is now ready for link inference
```

**Notes:**

- K is automatically calculated as 2.5% of total cells
- Uses balanced KNN for better performance
- Parallel processing controlled by `config.GRN_N_JOBS`
- Modifies the Oracle object in place

### `run_links()`

Infer gene regulatory network links and filter by significance.

**Signature:**

```python
def run_links(
    oracle: co.Oracle,
    cluster_column_name: str,
    p_cutoff: float = 0.001
) -> co.Links
```

**Parameters:**

- `oracle` (co.Oracle): Oracle object with PCA and KNN imputation completed
- `cluster_column_name` (str): Name of the column in oracle.adata.obs containing cluster information
- `p_cutoff` (float, default=0.001): P-value cutoff for filtering regulatory links. Lower values = stricter filtering

**Returns:**

- `links` (co.Links): Links object containing the inferred GRN

**Description:**

This function:

1. Infers gene regulatory links using ridge regression (alpha=10)
2. Filters links by p-value significance
3. Plots the degree distribution of the network
4. Calculates network scores

The resulting Links object contains:

- TF-target gene relationships
- Regulatory coefficients
- Statistical significance values
- Network topology metrics

**Example:**

```python
from genecircuitry.celloracle_processing import (
    create_oracle_object,
    run_PCA,
    run_KNN,
    run_links
)
from genecircuitry import config

# Complete GRN inference pipeline
oracle = create_oracle_object(adata, 'leiden', 'X_umap')
oracle = run_PCA(oracle)
run_KNN(oracle, n_comps=50)

# Infer regulatory links
links = run_links(
    oracle,
    cluster_column_name='leiden',
    p_cutoff=0.001
)

# Links object contains the GRN
print(f"Network has {len(links.links_dict)} regulatory links")

# Save the network
links.to_hdf5(f"{config.OUTPUT_DIR}/grn_links.celloracle.links")
```

**Output:**

- Saves degree distribution plot to `config.FIGURES_DIR/grn_degree_distribution.png`
- Prints network statistics and scores

**Notes:**

- Uses ridge regression with alpha=10 for regularization
- Stricter p_cutoff (lower value) = fewer but more significant links
- Network scores help assess GRN quality
- Verbose output shows detailed progress

## Complete Workflow Example

```python
import scanpy as sc
from genecircuitry import config, set_random_seed
from genecircuitry.celloracle_processing import (
    create_oracle_object,
    run_PCA,
    run_KNN,
    run_links
)

# Set up
set_random_seed(42)

# Load preprocessed single-cell data
adata = sc.read_h5ad('preprocessed_data.h5ad')

# Ensure you have:
# - Clustering results in adata.obs (e.g., 'leiden')
# - Dimensionality reduction in adata.obsm (e.g., 'X_umap')
# - Raw counts in adata.layers (optional but recommended)

# Step 1: Create Oracle object
print("Creating Oracle object...")
oracle = create_oracle_object(
    adata=adata,
    cluster_column_name='leiden',
    embedding_name='X_umap',
    raw_count_layer='raw_counts'  # or None for normalized
)

# Step 2: Perform PCA
print("Running PCA...")
oracle = run_PCA(oracle)

# Step 3: KNN imputation
print("Running KNN imputation...")
run_KNN(oracle, n_comps=50)

# Step 4: Infer regulatory links
print("Inferring gene regulatory network...")
links = run_links(
    oracle,
    cluster_column_name='leiden',
    p_cutoff=0.001
)

# Step 5: Save results
print("Saving results...")
oracle.to_hdf5(f"{config.OUTPUT_DIR}/oracle_object.celloracle.oracle")
links.to_hdf5(f"{config.OUTPUT_DIR}/grn_links.celloracle.links")

print("CellOracle analysis complete!")
```

## Configuration Parameters

The module uses the following config parameters:

```python
config.GRN_N_JOBS        # Parallel jobs for KNN (-1 = all cores)
config.FIGURES_DIR       # Directory for saving plots
config.OUTPUT_DIR        # Directory for saving results
```

Update these as needed:

```python
from genecircuitry import config

config.update_config(
    GRN_N_JOBS=8,
    FIGURES_DIR='my_figures',
    OUTPUT_DIR='my_results'
)
```

## Data Requirements

### Input AnnData Requirements

Your AnnData object should have:

1. **Gene expression data**:

   - `adata.X`: Normalized counts (log-transformed)
   - `adata.layers['raw_counts']`: Raw counts (recommended)

2. **Cell metadata** (`adata.obs`):

   - Cluster assignments (e.g., 'leiden', 'louvain', 'cell_type')
   - Any other cell annotations

3. **Dimensionality reduction** (`adata.obsm`):

   - UMAP coordinates ('X_umap')
   - Or PCA coordinates ('X_pca')

4. **Gene metadata** (`adata.var`):
   - Gene names as index
   - Preferably gene symbols (required for CellOracle base GRN)

### Example Preprocessing

```python
import scanpy as sc
from genecircuitry.preprocessing import perform_qc, perform_grn_pre_processing

# Load and preprocess
adata = sc.read_h5ad('raw_data.h5ad')

# Store raw counts before normalization
adata.layers['raw_counts'] = adata.X.copy()

# QC and filtering
adata = perform_qc(adata)

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Feature selection and PCA
adata = perform_grn_pre_processing(adata)

# Clustering
sc.tl.leiden(adata, resolution=1.0)

# UMAP
sc.tl.umap(adata)

# Now ready for CellOracle
```

## Tips and Best Practices

### 1. Choose Appropriate Clusters

- Use biologically meaningful clusters (cell types, states)
- Avoid over-clustering (too many small clusters)
- Ensure clusters have sufficient cells (>50 recommended)

### 2. PCA Components

- Default automatic selection works well in most cases
- Can manually override if needed
- More components = more detail but slower computation

### 3. KNN Parameters

- Auto-calculated k (2.5% of cells) is usually optimal
- For small datasets (<1000 cells), may need manual adjustment
- Larger k = more smoothing but may lose biological signal

### 4. P-value Cutoff

- Default 0.001 provides good balance
- Stricter (0.0001) for high-confidence links only
- More lenient (0.01) to include more potential regulators
- Consider multiple cutoffs and compare results

### 5. Computational Considerations

- CellOracle is memory-intensive
- For large datasets (>50k cells), consider:
  - Downsampling cells
  - Running on HPC cluster
  - Increasing available RAM
- Use `config.GRN_N_JOBS` to control parallel processing

### 6. Quality Control

- Examine degree distribution plot
- Check network scores
- Validate predicted links with known biology
- Compare across different p-value cutoffs

## Troubleshooting

### Error: "Gene names not found in base GRN"

- Ensure genes are human gene symbols
- Convert Ensembl IDs to symbols if needed
- Check gene name formatting (case-sensitive)

### Error: "Not enough memory"

- Reduce number of cells
- Close other applications
- Use a machine with more RAM
- Reduce n_comps

### Warning: "Few significant links found"

- Try more lenient p_cutoff
- Check data quality (noise, depth)
- Ensure proper normalization
- Verify cluster quality

### Slow Performance

- Reduce `n_comps`
- Set `config.GRN_N_JOBS` appropriately
- Consider downsampling cells
- Use balanced KNN (already default)

## References

1. CellOracle documentation: https://morris-lab.github.io/CellOracle.documentation/
2. Original paper: Kamimoto et al., Nature, 2023
3. Base GRN sources: CellOracle human promoter database

## See Also

- `genecircuitry.preprocessing` - Data preprocessing functions
- `genecircuitry.config` - Configuration management
- `docs/CONFIG.md` - Configuration documentation
- `examples/` - Example workflows
