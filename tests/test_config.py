"""
Tests for configuration module
"""

import pytest
import numpy as np
from trnspot import config


def test_random_seed():
    """Test random seed setting"""
    config.set_random_seed(42)

    # Generate random numbers
    rand1 = np.random.rand(10)

    # Reset seed
    config.set_random_seed(42)
    rand2 = np.random.rand(10)

    # Should be identical
    assert np.allclose(rand1, rand2)


def test_get_config():
    """Test getting configuration as dictionary"""
    cfg = config.get_config()

    # Check it's a dictionary
    assert isinstance(cfg, dict)

    # Check essential keys exist
    assert "RANDOM_SEED" in cfg
    assert "QC_MIN_GENES" in cfg
    assert "QC_MIN_COUNTS" in cfg
    assert "PLOT_DPI" in cfg

    # Check values
    assert isinstance(cfg["RANDOM_SEED"], int)
    assert isinstance(cfg["QC_MIN_GENES"], int)
    assert isinstance(cfg["PLOT_DPI"], int)


def test_update_config():
    """Test updating configuration"""
    # Store original values
    original_seed = config.RANDOM_SEED
    original_min_genes = config.QC_MIN_GENES

    # Update config
    config.update_config(RANDOM_SEED=999, QC_MIN_GENES=500)

    # Check updated values
    assert config.RANDOM_SEED == 999
    assert config.QC_MIN_GENES == 500

    # Restore original values
    config.update_config(RANDOM_SEED=original_seed, QC_MIN_GENES=original_min_genes)


def test_config_values():
    """Test that config values are reasonable"""
    # QC parameters
    assert config.QC_MIN_GENES > 0
    assert config.QC_MIN_COUNTS > 0
    assert 0 <= config.QC_PCT_MT_MAX <= 100

    # Plotting parameters
    assert config.PLOT_DPI > 0
    assert config.PLOT_FORMAT in ["png", "pdf", "svg", "jpg"]

    # Preprocessing parameters
    assert config.NORMALIZE_TARGET_SUM > 0
    assert config.HVGS_N_TOP_GENES > 0
    assert config.PCA_N_COMPS > 0

    # Neighbor parameters
    assert config.NEIGHBORS_N_NEIGHBORS > 0
    assert config.NEIGHBORS_N_PCS > 0


def test_print_config(capsys):
    """Test printing configuration"""
    config.print_config()

    # Capture printed output
    captured = capsys.readouterr()

    # Check that output contains expected text
    assert "TRNspot Configuration" in captured.out
    assert "RANDOM_SEED" in captured.out
    assert "QC_MIN_GENES" in captured.out


def test_config_types():
    """Test that config values have correct types"""
    cfg = config.get_config()

    # Integer values
    int_params = [
        "RANDOM_SEED",
        "QC_MIN_GENES",
        "QC_MIN_COUNTS",
        "QC_MIN_CELLS",
        "PLOT_DPI",
        "HVGS_N_TOP_GENES",
        "PCA_N_COMPS",
        "NEIGHBORS_N_NEIGHBORS",
        "NEIGHBORS_N_PCS",
    ]
    for param in int_params:
        if cfg[param] is not None:
            assert isinstance(cfg[param], (int, float)), f"{param} should be numeric"

    # Float values
    float_params = [
        "QC_PCT_MT_MAX",
        "HVGS_MIN_MEAN",
        "HVGS_MAX_MEAN",
        "HVGS_MIN_DISP",
        "LEIDEN_RESOLUTION",
        "LOUVAIN_RESOLUTION",
        "UMAP_MIN_DIST",
        "UMAP_SPREAD",
    ]
    for param in float_params:
        assert isinstance(cfg[param], (int, float)), f"{param} should be numeric"

    # String values
    string_params = [
        "PLOT_FORMAT",
        "NEIGHBORS_METHOD",
        "NEIGHBORS_METRIC",
        "OUTPUT_DIR",
        "CACHE_DIR",
        "FIGURES_DIR",
        "LOG_LEVEL",
    ]
    for param in string_params:
        assert isinstance(cfg[param], str), f"{param} should be string"

    # Boolean values
    bool_params = ["VERBOSE", "LOW_MEMORY"]
    for param in bool_params:
        assert isinstance(cfg[param], bool), f"{param} should be boolean"

    # List values
    list_params = ["ENRICHMENT_GENE_SETS"]
    for param in list_params:
        assert isinstance(cfg[param], list), f"{param} should be list"


def test_invalid_config_update():
    """Test updating with invalid parameter name"""
    # This should print a warning but not raise an error
    config.update_config(INVALID_PARAM=123)

    # Check that invalid parameter wasn't added
    cfg = config.get_config()
    assert "INVALID_PARAM" not in cfg


def test_enrichment_gene_sets_default():
    """Test that ENRICHMENT_GENE_SETS has expected default"""
    assert isinstance(config.ENRICHMENT_GENE_SETS, list)
    assert len(config.ENRICHMENT_GENE_SETS) >= 1
    assert "MSigDB_Hallmark_2020" in config.ENRICHMENT_GENE_SETS


def test_enrichment_gene_sets_in_get_config():
    """Test that ENRICHMENT_GENE_SETS is present in get_config()"""
    cfg = config.get_config()
    assert "ENRICHMENT_GENE_SETS" in cfg
    assert cfg["ENRICHMENT_GENE_SETS"] == config.ENRICHMENT_GENE_SETS


def test_enrichment_gene_sets_update():
    """Test updating ENRICHMENT_GENE_SETS via update_config"""
    original = list(config.ENRICHMENT_GENE_SETS)

    new_sets = ["MSigDB_Hallmark_2020", "KEGG_2021_Human"]
    config.update_config(ENRICHMENT_GENE_SETS=new_sets)
    assert config.ENRICHMENT_GENE_SETS == new_sets

    # Restore
    config.update_config(ENRICHMENT_GENE_SETS=original)
