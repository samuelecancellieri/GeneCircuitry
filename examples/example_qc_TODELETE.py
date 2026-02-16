# TO REMOVE: Broken script — imports plot_qc_violin/plot_qc_scatter from trnspot.preprocessing
# which don't exist. Also uses nonexistent plot=True parameter. Remove this file.
"""
Example script demonstrating QC functionality in TRNspot
"""

import scanpy as sc
from trnspot.preprocessing import perform_qc, plot_qc_violin, plot_qc_scatter

# Load example data (you can replace this with your own data)
# Example: adata = sc.read_h5ad('your_data.h5ad')
# For demonstration, we'll use a built-in dataset
adata = sc.datasets.pbmc3k()

print("=" * 60)
print("TRNspot QC Example")
print("=" * 60)

# Perform QC with default parameters
adata_qc = perform_qc(
    adata,
    min_genes=200,
    min_counts=500,
    pct_counts_mt_max=20.0,
    plot=True,
    save_plots="qc_histograms.png",
)

# Generate additional QC plots
print("\n" + "=" * 60)
print("Generating violin plots...")
print("=" * 60)
plot_qc_violin(adata_qc, figsize=(12, 4), save_path="qc_violin.png")

print("\n" + "=" * 60)
print("Generating scatter plots...")
print("=" * 60)
plot_qc_scatter(
    adata_qc, color_by="pct_counts_mt", figsize=(12, 4), save_path="qc_scatter.png"
)

print("\n" + "=" * 60)
print("QC complete!")
print("=" * 60)
print(f"Final data: {adata_qc.n_obs} cells × {adata_qc.n_vars} genes")

# Optionally save the filtered data
# adata_qc.write('adata_qc_filtered.h5ad')
