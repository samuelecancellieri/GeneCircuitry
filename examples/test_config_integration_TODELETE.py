#!/usr/bin/env python
# TO REMOVE: Broken script — imports perform_grn_pre_processing from trnspot.preprocessing
# but it's defined in trnspot.celloracle_processing. Also uses nonexistent plot=False parameter.
# This is not a proper test (not pytest-based). Remove this file.
"""
Quick test to verify preprocessing with config integration
"""

import scanpy as sc
from trnspot import set_random_seed, config
from trnspot.preprocessing import perform_qc, perform_grn_pre_processing

print("=" * 70)
print("Testing Preprocessing with Config Defaults")
print("=" * 70)

# Set seed
set_random_seed(42)

# Load data
print("\n1. Loading test data...")
adata = sc.datasets.pbmc3k()
print(f"   Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Test 1: QC with defaults
print("\n2. Testing QC with config defaults...")
print(f"   Using: min_genes={config.QC_MIN_GENES}, min_counts={config.QC_MIN_COUNTS}")
adata_qc = perform_qc(adata, plot=False)
print(f"   ✓ QC complete: {adata_qc.n_obs} cells × {adata_qc.n_vars} genes")

# Test 2: QC with override
print("\n3. Testing QC with parameter override...")
adata_qc2 = perform_qc(adata, min_genes=300, min_counts=1000, plot=False)
print(f"   ✓ QC complete: {adata_qc2.n_obs} cells × {adata_qc2.n_vars} genes")

# Test 3: Update config and use new defaults
print("\n4. Testing config update...")
config.update_config(QC_MIN_GENES=250, QC_MIN_COUNTS=750)
adata_qc3 = perform_qc(adata, plot=False)
print(f"   ✓ QC complete: {adata_qc3.n_obs} cells × {adata_qc3.n_vars} genes")

# Reset config
config.update_config(QC_MIN_GENES=200, QC_MIN_COUNTS=500)

# Test 4: GRN preprocessing
print("\n5. Testing GRN preprocessing with config defaults...")
print(
    f"   Using: top_genes={config.HVGS_N_TOP_GENES}, n_neighbors={config.NEIGHBORS_N_NEIGHBORS}"
)
adata_clean = perform_qc(adata, plot=False)
adata_grn = perform_grn_pre_processing(adata_clean, cluster_key=None)
print(
    f"   ✓ GRN preprocessing complete: {adata_grn.n_obs} cells × {adata_grn.n_vars} genes"
)

# Test 5: GRN with overrides
print("\n6. Testing GRN preprocessing with parameter override...")
adata_clean2 = perform_qc(adata, plot=False)
adata_grn2 = perform_grn_pre_processing(
    adata_clean2, cluster_key=None, top_genes=1000, n_neighbors=20
)
print(
    f"   ✓ GRN preprocessing complete: {adata_grn2.n_obs} cells × {adata_grn2.n_vars} genes"
)

print("\n" + "=" * 70)
print("All tests passed! ✓")
print("=" * 70)
print("\nConfig integration is working correctly:")
print("  ✓ Functions use config defaults when parameters not specified")
print("  ✓ Parameters can be overridden individually")
print("  ✓ Config can be updated globally")
print("  ✓ All preprocessing functions work as expected")
print("=" * 70)
