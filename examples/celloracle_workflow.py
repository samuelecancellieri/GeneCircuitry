"""
Complete CellOracle workflow example for GRN inference
"""

import scanpy as sc
from genecircuitry import set_random_seed, set_scanpy_settings, config

from genecircuitry.preprocessing import (
    perform_qc,
    perform_normalization,
)
from genecircuitry.celloracle_processing import (
    perform_grn_pre_processing,
)

from genecircuitry.celloracle_processing import (
    create_oracle_object,
    run_PCA,
    run_KNN,
    run_links,
)

print("=" * 70)
print("CellOracle GRN Inference Workflow")
print("=" * 70)

# ============================================================================
# Setup and Configuration
# ============================================================================
print("\n[Step 1/7] Setting up configuration...")
set_random_seed(42)
set_scanpy_settings()

# Configure for GRN analysis
config.update_config(
    GRN_N_JOBS=10, FIGURES_DIR="figures", OUTPUT_DIR="output"  # Use 10 cores
)

print(f"  GRN parallel jobs: {config.GRN_N_JOBS}")
print(f"  Figures directory: {config.FIGURES_DIR}")
print(f"  Output directory: {config.OUTPUT_DIR}")

# ============================================================================
# Load and Preprocess Data
# ============================================================================
print("\n[Step 2/7] Loading and preprocessing data...")

# Load example data
adata = sc.datasets.paul15()
print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Perform QC
print("  Performing quality control...")
adata = perform_qc(adata, save_plots="qc_grn.png")
adata = perform_normalization(adata)
adata = perform_grn_pre_processing(adata, cluster_key="paul15_clusters")

print(f"  Data prepared for CellOracle: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  Clusters: {adata.obs['paul15_clusters'].nunique()} clusters")

# ============================================================================
# Create Oracle Object
# ============================================================================
print("\n[Step 3/7] Creating CellOracle Oracle object...")

try:
    oracle = create_oracle_object(
        adata=adata,
        cluster_column_name="paul15_clusters",
        embedding_name="X_draw_graph_fa",
        species="mouse",
        raw_count_layer="raw_counts",
        TG_to_TF_dictionary=None,  # Use CellOracle's base GRN
    )
    print("  ✓ Oracle object created successfully")
    print(f"  Cells: {oracle.adata.n_obs}")
    print(f"  Genes: {oracle.adata.n_vars}")
except Exception as e:
    print(f"  ✗ Error creating Oracle object: {e}")
    print("  Note: CellOracle must be installed for this step")
    print("  Install with: pip install celloracle")
    exit(1)

# ============================================================================
# Perform PCA on Oracle Object
# ============================================================================
print("\n[Step 4/7] Performing PCA on Oracle object...")

try:
    oracle, n_comps = run_PCA(oracle)
    print(f"  ✓ PCA complete")
    print(f"  Selected {n_comps} principal components")
except Exception as e:
    print(f"  ✗ Error during PCA: {e}")
    exit(1)

# ============================================================================
# KNN Imputation
# ============================================================================
print("\n[Step 5/7] Performing KNN imputation...")

try:
    oracle = run_KNN(oracle, n_comps=n_comps)
    print("  ✓ KNN imputation complete")
    print(f"  Used {n_comps} PCA components")
except Exception as e:
    print(f"  ✗ Error during KNN imputation: {e}")
    exit(1)

# ============================================================================
# Infer Regulatory Links
# ============================================================================
print("\n[Step 6/7] Inferring gene regulatory network...")

try:
    links = run_links(oracle, cluster_column_name="paul15_clusters", p_cutoff=0.001)
    print("  ✓ GRN inference complete")
    print(f"  P-value cutoff: 0.001")

    # Get network statistics
    if hasattr(links, "links_dict"):
        n_links = sum(len(targets) for targets in links.links_dict.values())
        n_tfs = len(links.links_dict)
        print(f"  Network contains:")
        print(f"    - {n_tfs} transcription factors")
        print(f"    - {n_links} regulatory links")
except Exception as e:
    print(f"  ✗ Error during link inference: {e}")
    exit(1)

# ============================================================================
# Save Results
# ============================================================================
print("\n[Step 7/7] Saving results...")

import os

os.makedirs(config.OUTPUT_DIR, exist_ok=True)
os.makedirs(config.FIGURES_DIR, exist_ok=True)

try:
    # Save Oracle object
    oracle_path = f"{config.OUTPUT_DIR}/oracle_object.celloracle.oracle"
    oracle.to_hdf5(oracle_path)
    print(f"  ✓ Saved Oracle object to: {oracle_path}")

    # Save Links object
    links_path = f"{config.OUTPUT_DIR}/grn_links.celloracle.links"
    links.to_hdf5(links_path)
    print(f"  ✓ Saved GRN links to: {links_path}")

    # Save preprocessed AnnData
    adata_path = f"{config.OUTPUT_DIR}/preprocessed_adata.h5ad"
    adata.write(adata_path)
    print(f"  ✓ Saved AnnData to: {adata_path}")

    # Save configuration
    import json

    config_path = f"{config.OUTPUT_DIR}/analysis_config.json"
    with open(config_path, "w") as f:
        json.dump(config.get_config(), f, indent=2, default=str)
    print(f"  ✓ Saved configuration to: {config_path}")

except Exception as e:
    print(f"  ✗ Error saving results: {e}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("CellOracle Analysis Complete!")
print("=" * 70)
print("\nWorkflow Summary:")
print(f"  1. ✓ Configuration and setup")
print(f"  2. ✓ Data preprocessing ({adata.n_obs} cells, {adata.n_vars} genes)")
print(f"  3. ✓ Oracle object created")
print(f"  4. ✓ PCA performed ({n_comps} components)")
print(f"  5. ✓ KNN imputation completed")
print(f"  6. ✓ GRN inferred (p < 0.001)")
print(f"  7. ✓ Results saved")

print("\nOutput Files:")
print(f"  - {oracle_path}")
print(f"  - {links_path}")
print(f"  - {adata_path}")
print(f"  - {config_path}")
print(f"  - {config.FIGURES_DIR}/grn_degree_distribution/degree_distribution.png")

print("\nNext Steps:")
print("  1. Examine the degree distribution plot")
print("  2. Validate predicted regulatory links")
print("  3. Perform perturbation simulations")
print("  4. Visualize gene regulatory networks")
print("  5. Compare across cell types/clusters")

print("\nFor more information:")
print("  - Documentation: docs/CELLORACLE_PROCESSING.md")
print("  - CellOracle docs: https://morris-lab.github.io/CellOracle.documentation/")
print("=" * 70)
