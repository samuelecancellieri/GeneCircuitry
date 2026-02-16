# TO REMOVE: Broken script — imports plot_qc_violin/plot_qc_scatter from trnspot.preprocessing
# which don't exist. Also uses nonexistent plot=True parameter. Remove this file.
"""
Quick demo of TRNspot QC functions
Run this in an interactive Python session or Jupyter notebook
"""

import scanpy as sc
from trnspot.preprocessing import perform_qc, plot_qc_violin, plot_qc_scatter

# Load example dataset (PBMC 3k cells from 10x Genomics)
print("Loading example dataset...")
adata = sc.datasets.pbmc3k()
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Basic QC with default parameters
print("\n" + "=" * 60)
print("Running QC with default parameters...")
print("=" * 60)
adata_qc = perform_qc(adata, plot=True)

# More stringent filtering
print("\n" + "=" * 60)
print("Running QC with stringent parameters...")
print("=" * 60)
adata_strict = perform_qc(
    adata,
    min_genes=300,
    min_counts=1000,
    max_counts=25000,
    pct_counts_mt_max=10.0,
    plot=True,
)

# Generate additional plots
print("\n" + "=" * 60)
print("Generating additional QC visualizations...")
print("=" * 60)
plot_qc_violin(adata_qc)
plot_qc_scatter(adata_qc, color_by="pct_counts_mt")

print("\nQC analysis complete!")
