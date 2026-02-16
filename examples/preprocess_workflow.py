"""
Example demonstrating preprocessing with config defaults
"""

import scanpy as sc
from trnspot import set_random_seed, set_scanpy_settings, config
from trnspot.preprocessing import (
    perform_normalization,
    perform_qc,
)
from trnspot.celloracle_processing import (
    perform_grn_pre_processing,
)

# ============================================================================
# Setup
# ============================================================================
print("=" * 70)
print("Preprocessing with Config Defaults Example")
print("=" * 70)

# Set random seed
set_random_seed(42)
set_scanpy_settings()

# Load example data
print("\nLoading data...")
adata = sc.datasets.paul15()
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# ============================================================================
# Example 1: Using Default Config Values
# ============================================================================
print("\n" + "=" * 70)
print("Example 1: Using Default Config Values")
print("=" * 70)

# Perform QC with all defaults from config
adata_qc = perform_qc(adata, save_plots="qc1.png")

print(f"\nConfig defaults used:")
print(f"  min_genes: {config.QC_MIN_GENES}")
print(f"  min_counts: {config.QC_MIN_COUNTS}")
print(f"  pct_counts_mt_max: {config.QC_PCT_MT_MAX}")
print(f"  min_cells: {config.QC_MIN_CELLS}")

# ============================================================================
# Example 2: Override Specific Parameters
# ============================================================================
print("\n" + "=" * 70)
print("Example 2: Override Specific Parameters")
print("=" * 70)

# Override only specific parameters
adata_qc2 = perform_qc(
    adata,
    min_genes=300,  # Override this
    min_counts=1000,  # Override this
    # pct_counts_mt_max uses default from config
    # min_cells uses default from config
    save_plots="qc2.png",
)

# ============================================================================
# Example 3: Update Config Globally
# ============================================================================
print("\n" + "=" * 70)
print("Example 3: Update Config Globally")
print("=" * 70)

# Update config for stricter QC
config.update_config(QC_MIN_GENES=500, QC_MIN_COUNTS=2000, QC_PCT_MT_MAX=15.0)

print("\nUpdated config:")
print(f"  min_genes: {config.QC_MIN_GENES}")
print(f"  min_counts: {config.QC_MIN_COUNTS}")
print(f"  pct_counts_mt_max: {config.QC_PCT_MT_MAX}")

# Now perform_qc will use these updated defaults
adata_qc3 = perform_qc(adata, save_plots="qc3.png")

# ============================================================================
# Example 4: GRN Preprocessing with Config Defaults
# ============================================================================
print("\n" + "=" * 70)
print("Example 4: GRN Preprocessing with Config Defaults")
print("=" * 70)

# Reset config to defaults
config.update_config(QC_MIN_GENES=200, QC_MIN_COUNTS=500, QC_PCT_MT_MAX=20.0)

# First do basic QC
adata_qc4 = perform_qc(adata, save_plots="qc4.png")

# Normalize
adata_qc4 = perform_normalization(adata_qc4)

# Then do GRN preprocessing with config defaults
print("\nPerforming GRN preprocessing with config defaults...")
print(f"  top_genes: {config.HVGS_N_TOP_GENES}")
print(f"  n_neighbors: {config.NEIGHBORS_N_NEIGHBORS}")
print(f"  n_pcs: {config.NEIGHBORS_N_PCS}")
print(f"  svd_solver: {config.PCA_SVDSOLVE}")

adata_grn = perform_grn_pre_processing(
    adata_qc4, cluster_key="paul15_clusters"  # Will use all cells
)

# ============================================================================
# Example 5: Complete Config-Based Workflow
# ============================================================================
print("\n" + "=" * 70)
print("Example 5: Complete Config-Based Workflow")
print("=" * 70)

# Set up a specific analysis profile
print("\nSetting up 'strict_qc' profile...")
config.update_config(
    RANDOM_SEED=42,
    QC_MIN_GENES=300,
    QC_MIN_COUNTS=1000,
    QC_PCT_MT_MAX=10.0,
    QC_MIN_CELLS=5,
    HVGS_N_TOP_GENES=3000,
    NEIGHBORS_N_NEIGHBORS=20,
    NEIGHBORS_N_PCS=50,
)

# Print the profile
print("\nStrict QC Profile:")
print(f"  QC_MIN_GENES: {config.QC_MIN_GENES}")
print(f"  QC_MIN_COUNTS: {config.QC_MIN_COUNTS}")
print(f"  QC_PCT_MT_MAX: {config.QC_PCT_MT_MAX}")
print(f"  HVGS_N_TOP_GENES: {config.HVGS_N_TOP_GENES}")

# Run complete workflow with this profile
print("\nRunning complete workflow...")
adata_workflow = perform_qc(adata, save_plots="qc5.png")
print(f"After QC: {adata_workflow.n_obs} cells × {adata_workflow.n_vars} genes")

# Normalize
adata_workflow = perform_normalization(adata_workflow)

# GRN preprocessing
adata_workflow = perform_grn_pre_processing(
    adata_workflow, cluster_key="paul15_clusters"
)
print(
    f"After GRN preprocessing: {adata_workflow.n_obs} cells × {adata_workflow.n_vars} genes"
)

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(
    """
Key Takeaways:
1. All functions use config defaults when parameters are not specified
2. You can override individual parameters while using defaults for others
3. Update config globally to change defaults for all subsequent calls
4. Create configuration profiles for different analysis types
5. Config-based workflow ensures consistency and reproducibility

Benefits:
- Consistent parameters across analysis
- Easy to share and reproduce analysis settings
- Less code duplication
- Centralized parameter management
"""
)

print("=" * 70)
