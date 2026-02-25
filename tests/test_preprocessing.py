"""
Tests for preprocessing module
"""

import pytest
import numpy as np
import scanpy as sc
from anndata import AnnData
from genecircuitry.preprocessing import perform_qc


@pytest.fixture
def sample_adata():
    """Create a sample AnnData object for testing"""
    np.random.seed(42)
    n_cells = 100
    n_genes = 50

    # Create count matrix
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))

    # Create gene names (include some MT genes)
    var_names = [f"Gene_{i}" for i in range(n_genes - 5)]
    var_names.extend([f"MT-Gene{i}" for i in range(5)])

    # Create AnnData object
    adata = AnnData(X=X)
    adata.var_names = var_names
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    return adata


def test_perform_qc_basic(sample_adata):
    """Test basic QC functionality"""
    adata_qc = perform_qc(
        sample_adata,
        min_genes=5,
        min_counts=10,
        min_cells=1,
        pct_counts_mt_max=50.0,
        plot=False,
    )

    # Check that QC metrics were calculated
    assert "n_genes_by_counts" in adata_qc.obs.columns
    assert "total_counts" in adata_qc.obs.columns
    assert "pct_counts_mt" in adata_qc.obs.columns

    # Check that filtering was applied
    assert adata_qc.n_obs <= sample_adata.n_obs
    assert adata_qc.n_vars <= sample_adata.n_vars

    # Check that MT genes were identified
    assert "mt" in adata_qc.var.columns


def test_perform_qc_with_config_defaults(sample_adata):
    """Test QC with config defaults"""
    from genecircuitry import config

    # Update config
    config.update_config(
        QC_MIN_GENES=5, QC_MIN_COUNTS=10, QC_MIN_CELLS=1, QC_PCT_MT_MAX=50.0
    )

    # Use defaults from config
    adata_qc = perform_qc(sample_adata, plot=False)

    # Check that filtering was applied
    assert adata_qc.n_obs <= sample_adata.n_obs
    assert adata_qc.n_vars <= sample_adata.n_vars

    # Reset config
    config.update_config(
        QC_MIN_GENES=200, QC_MIN_COUNTS=500, QC_MIN_CELLS=10, QC_PCT_MT_MAX=20.0
    )


def test_perform_qc_strict_filtering(sample_adata):
    """Test QC with strict filtering parameters"""
    adata_qc = perform_qc(
        sample_adata,
        min_genes=20,
        min_counts=100,
        min_cells=1,
        max_counts=1000,
        pct_counts_mt_max=10.0,
        plot=False,
    )

    # Check filtering conditions
    assert all(adata_qc.obs["n_genes_by_counts"] >= 20)
    assert all(adata_qc.obs["total_counts"] >= 100)
    assert all(adata_qc.obs["total_counts"] < 1000)
    assert all(adata_qc.obs["pct_counts_mt"] < 10.0)


def test_plot_qc_violin(sample_adata):
    """Test violin plot generation via plotting subpackage"""
    from genecircuitry.plotting.qc_plots import plot_qc_violin_pre_filter

    # First perform QC to get metrics
    adata_qc = perform_qc(sample_adata, plot=False)

    # Test violin plot
    plot_qc_violin_pre_filter(adata_qc, skip_existing=False)


def test_plot_qc_scatter(sample_adata):
    """Test scatter plot generation via plotting subpackage"""
    from genecircuitry.plotting.qc_plots import plot_qc_scatter_pre_filter

    # First perform QC to get metrics
    adata_qc = perform_qc(sample_adata, plot=False)

    # Test scatter plot
    plot_qc_scatter_pre_filter(adata_qc, skip_existing=False)
