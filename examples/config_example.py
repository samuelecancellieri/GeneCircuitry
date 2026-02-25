"""
Example demonstrating configuration usage in GeneCircuitry. This script shows how to access, update, and utilize configuration settings
"""

import genecircuitry
from genecircuitry import config

# Print current configuration
print("=" * 60)
print("Default Configuration")
print("=" * 60)
config.print_config()

# Access individual config values
print("\n" + "=" * 60)
print("Accessing Individual Config Values")
print("=" * 60)
print(f"Random seed: {config.RANDOM_SEED}")
print(f"Min genes for QC: {config.QC_MIN_GENES}")
print(f"Min counts for QC: {config.QC_MIN_COUNTS}")
print(f"Max MT percentage: {config.QC_PCT_MT_MAX}")
print(f"Plot DPI: {config.PLOT_DPI}")

# Update configuration
print("\n" + "=" * 60)
print("Updating Configuration")
print("=" * 60)
config.update_config(RANDOM_SEED=123, QC_MIN_GENES=300, QC_MIN_COUNTS=1000)

# Set random seed
print("\n" + "=" * 60)
print("Setting Random Seed")
print("=" * 60)
genecircuitry.set_random_seed(456)

# Get config as dictionary
print("\n" + "=" * 60)
print("Get Config as Dictionary")
print("=" * 60)
cfg = genecircuitry.get_config()
print(f"Config type: {type(cfg)}")
print(f"Number of config parameters: {len(cfg)}")
print(f"\nSample parameters:")
for key in list(cfg.keys())[:5]:
    print(f"  {key}: {cfg[key]}")

# Using config in preprocessing
print("\n" + "=" * 60)
print("Using Config with Preprocessing")
print("=" * 60)
print(
    """
Example usage:

from genecircuitry import config
from genecircuitry.preprocessing import perform_qc

# Use default config values
adata_qc = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    min_counts=config.QC_MIN_COUNTS,
    pct_counts_mt_max=config.QC_PCT_MT_MAX
)

# Or customize for specific analysis
config.update_config(
    QC_MIN_GENES=500,
    QC_PCT_MT_MAX=15.0
)
adata_qc = perform_qc(
    adata,
    min_genes=config.QC_MIN_GENES,
    pct_counts_mt_max=config.QC_PCT_MT_MAX
)
"""
)

print("Configuration example complete!")
