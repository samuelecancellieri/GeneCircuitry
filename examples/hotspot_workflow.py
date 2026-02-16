"""
Complete Hotspot workflow example for gene module identification
This workflow demonstrates the full process of identifying co-expressed gene modules
using Hotspot, which identifies genes with spatially autocorrelated expression patterns.
"""

import os
import scanpy as sc
import pandas as pd
from trnspot import set_random_seed, set_scanpy_settings, config

from trnspot.preprocessing import (
    perform_qc,
    perform_normalization,
)

from trnspot.hotspot_processing import (
    create_hotspot_object,
    run_hotspot_analysis,
)

print("=" * 70)
print("Hotspot Gene Module Identification Workflow")
print("=" * 70)

# ============================================================================
# Setup and Configuration
# ============================================================================
print("\n[Step 1/6] Setting up configuration...")
set_random_seed(42)
set_scanpy_settings()

# Configure for Hotspot analysis
config.update_config(
    HOTSPOT_N_JOBS=8,  # Use 8 cores for parallel processing
    HOTSPOT_N_NEIGHBORS=30,  # Number of neighbors for KNN graph
    HOTSPOT_FDR_THRESHOLD=0.05,  # FDR threshold for significant genes
    HOTSPOT_MIN_GENES_PER_MODULE=30,  # Minimum genes per module
    HOTSPOT_TOP_GENES=3000,  # Number of highly variable genes to use
    FIGURES_DIR="figures",
    OUTPUT_DIR="output",
)

print(f"  Hotspot parallel jobs: {config.HOTSPOT_N_JOBS}")
print(f"  KNN neighbors: {config.HOTSPOT_N_NEIGHBORS}")
print(f"  FDR threshold: {config.HOTSPOT_FDR_THRESHOLD}")
print(f"  Min genes per module: {config.HOTSPOT_MIN_GENES_PER_MODULE}")
print(f"  Top genes to analyze: {config.HOTSPOT_TOP_GENES}")
print(f"  Figures directory: {config.FIGURES_DIR}")
print(f"  Output directory: {config.OUTPUT_DIR}")

# Create output directories
os.makedirs(config.FIGURES_DIR, exist_ok=True)
os.makedirs(f"{config.OUTPUT_DIR}/hotspot", exist_ok=True)

# ============================================================================
# Load and Preprocess Data
# ============================================================================
print("\n[Step 2/6] Loading and preprocessing data...")

# Load example data (Paul et al. 2015 hematopoietic cells)
adata = sc.datasets.paul15()
print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Perform QC
print("\n  Performing quality control...")
adata = perform_qc(
    adata,
    min_genes=200,
    min_counts=500,
    save_plots=None,  # Skip plots for workflow
)
print(f"  After QC: {adata.n_obs} cells × {adata.n_vars} genes")

# Normalize and log-transform
print("\n  Normalizing data...")
adata = perform_normalization(adata)

# ============================================================================
# Dimensionality Reduction and Clustering
# ============================================================================
print("\n[Step 3/6] Performing dimensionality reduction and clustering...")

# Identify highly variable genes
print("  Identifying highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)

# PCA
print("  Running PCA...")
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
print(f"    Computed {adata.obsm['X_pca'].shape[1]} principal components")

# Compute neighbors graph
print("  Computing neighborhood graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

# UMAP embedding
print("  Computing UMAP embedding...")
sc.tl.umap(adata)

# Leiden clustering
print("  Performing Leiden clustering...")
sc.tl.leiden(adata, resolution=1.0)
n_clusters = len(adata.obs["leiden"].unique())
print(f"    Identified {n_clusters} clusters")

# ============================================================================
# Create Hotspot Object
# ============================================================================
print("\n[Step 4/6] Creating Hotspot object...")

try:
    hotspot_obj = create_hotspot_object(
        adata,
        top_genes=500,
        layer_key="raw_counts",
        model="danb",  # Depth-adjusted negative binomial model
        embedding_key="X_pca",
        normalization_key="total_counts",
    )
    print("  ✓ Hotspot object created successfully")
    print(f"    Model: danb (depth-adjusted negative binomial)")
    print(f"    Embedding: PCA")
    print(f"    Analyzing top 500 genes")

except Exception as e:
    print(f"  ✗ Error creating Hotspot object: {e}")
    print("\nTroubleshooting tips:")
    print("  - Ensure hotspot package is installed: pip install hotspot")
    print("  - Check that raw counts are available in adata.layers['raw_counts']")
    print("  - Verify that n_counts is in adata.obs")
    exit(1)

# ============================================================================
# Run Hotspot Analysis
# ============================================================================
print("\n[Step 5/6] Running Hotspot analysis...")
print("  This may take several minutes depending on dataset size...")

hotspot_obj = run_hotspot_analysis(hotspot_obj)

# try:
#     # Run the complete Hotspot analysis pipeline
#     hotspot_obj = run_hotspot_analysis(hotspot_obj)

#     print("  ✓ Hotspot analysis completed successfully")

#     Get results summary
#     autocorr_results = hotspot_obj.results
#     significant_genes = autocorr_results[
#         autocorr_results.FDR < config.HOTSPOT_FDR_THRESHOLD
#     ]

#     print(f"\n  Analysis Summary:")
#     print(f"    Total genes analyzed: {len(autocorr_results)}")
#     print(
#         f"    Significant genes (FDR < {config.HOTSPOT_FDR_THRESHOLD}): {len(significant_genes)}"
#     )

#     if hasattr(hotspot_obj, "modules"):
#         n_modules = len(hotspot_obj.modules.unique())
#         print(f"    Gene modules identified: {n_modules}")

#         # Module size distribution
#         module_sizes = hotspot_obj.modules.value_counts()
#         print(f"\n  Module Size Distribution:")
#         for module, size in module_sizes.head(10).items():
#             print(f"    Module {module}: {size} genes")

#     # Save local correlations plot
#     print(
#         f"\n  Local correlations plot saved to: {config.FIGURES_DIR}/hotspot_local_correlations.png"
#     )

# except Exception as e:
#     print(f"  ✗ Error during Hotspot analysis: {e}")
#     print("\nTroubleshooting tips:")
#     print("  - Check memory usage (Hotspot can be memory-intensive)")
#     print("  - Try reducing HOTSPOT_TOP_GENES if out of memory")
#     print("  - Ensure sufficient compute resources for parallel jobs")
#     exit(1)

# ============================================================================
# Save Results
# ============================================================================
print("\n[Step 6/6] Saving results...")

try:
    # Save module scores (already saved by run_hotspot_analysis)
    print(f"  ✓ Module scores saved to: {config.OUTPUT_DIR}/hotspot/hotspot_result.pkl")

    # Save autocorrelation results to CSV
    results_df = hotspot_obj.results
    results_path = f"{config.OUTPUT_DIR}/hotspot/autocorrelation_results.csv"
    results_df.to_csv(results_path)
    print(f"  ✓ Autocorrelation results saved to: {results_path}")

    # Save significant genes
    significant_genes = results_df[results_df.FDR < config.HOTSPOT_FDR_THRESHOLD]
    significant_path = f"{config.OUTPUT_DIR}/hotspot/significant_genes.csv"
    significant_genes.to_csv(significant_path)
    print(f"  ✓ Significant genes saved to: {significant_path}")

    # Save module assignments if available
    if hasattr(hotspot_obj, "modules"):
        modules_path = f"{config.OUTPUT_DIR}/hotspot/gene_modules.csv"
        hotspot_obj.modules.to_csv(modules_path)
        print(f"  ✓ Gene modules saved to: {modules_path}")

    # Add module scores to AnnData object
    if hasattr(hotspot_obj, "module_scores"):
        for module in hotspot_obj.module_scores.columns:
            adata.obs[f"module_{module}"] = hotspot_obj.module_scores[module].values

        # Save updated AnnData
        adata_path = f"{config.OUTPUT_DIR}/hotspot/adata_with_modules.h5ad"
        adata.write(adata_path)
        print(f"  ✓ AnnData with module scores saved to: {adata_path}")

except Exception as e:
    print(f"  ⚠ Warning: Could not save some results: {e}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Hotspot Workflow Completed Successfully!")
print("=" * 70)

print("\nOutput files:")
print(f"  1. Module scores (pickle): {config.OUTPUT_DIR}/hotspot/hotspot_result.pkl")
print(
    f"  2. Autocorrelation results (CSV): {config.OUTPUT_DIR}/hotspot/autocorrelation_results.csv"
)
print(
    f"  3. Significant genes (CSV): {config.OUTPUT_DIR}/hotspot/significant_genes.csv"
)
if hasattr(hotspot_obj, "modules"):
    print(f"  4. Gene modules (CSV): {config.OUTPUT_DIR}/hotspot/gene_modules.csv")
    print(
        f"  5. AnnData with modules (h5ad): {config.OUTPUT_DIR}/hotspot/adata_with_modules.h5ad"
    )
print(
    f"  6. Local correlations plot (PNG): {config.FIGURES_DIR}/hotspot_local_correlations.png"
)

print("\nNext steps:")
print("  - Examine significant genes and their autocorrelation scores")
print("  - Explore identified gene modules and their biological functions")
print("  - Visualize module scores on UMAP/t-SNE embeddings")
print("  - Perform gene ontology enrichment analysis on modules")
print("  - Compare modules across different cell types or conditions")

print("\nTo load results in a new session:")
print("  import pickle")
print("  with open('output/hotspot/hotspot_result.pkl', 'rb') as f:")
print("      module_scores = pickle.load(f)")
