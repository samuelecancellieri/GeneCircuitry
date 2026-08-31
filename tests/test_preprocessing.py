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


def test_perform_dimensionality_reduction_basic(sample_adata):
    """Test standard dimensionality reduction and clustering pipeline"""
    from genecircuitry.preprocessing import (
        perform_qc,
        perform_normalization,
        perform_dimensionality_reduction_clustering,
    )

    adata = perform_qc(sample_adata, min_genes=1, min_counts=1, min_cells=1, plot=False)
    adata.X = (adata.X.astype(np.float32) + 1.0) * 50.0
    adata = perform_normalization(adata)
    adata = perform_dimensionality_reduction_clustering(adata, cluster_key="leiden")

    assert "highly_variable" in adata.var.columns
    assert "X_pca" in adata.obsm
    assert "neighbors" in adata.uns
    assert "X_umap" in adata.obsm
    assert "leiden" in adata.obs.columns


def test_perform_dimensionality_reduction_skip_existing(sample_adata):
    """Test that existing embeddings are preserved when force=False"""
    from genecircuitry.preprocessing import (
        perform_qc,
        perform_normalization,
        perform_dimensionality_reduction_clustering,
    )

    adata = perform_qc(sample_adata, min_genes=1, min_counts=1, min_cells=1, plot=False)
    adata.X = (adata.X.astype(np.float32) + 1.0) * 50.0
    adata = perform_normalization(adata)

    # Inject dummy PCA embedding
    dummy_pca = np.ones((adata.n_obs, 5), dtype=np.float32) * 42.0
    adata.obsm["X_pca"] = dummy_pca.copy()

    # Run without force
    adata_out = perform_dimensionality_reduction_clustering(adata, force=False)

    # PCA should remain unchanged (dummy)
    np.testing.assert_array_equal(adata_out.obsm["X_pca"], dummy_pca)


def test_perform_dimensionality_reduction_force(sample_adata):
    """Test that existing embeddings are overwritten when force=True"""
    from genecircuitry.preprocessing import (
        perform_qc,
        perform_normalization,
        perform_dimensionality_reduction_clustering,
    )

    adata = perform_qc(sample_adata, min_genes=1, min_counts=1, min_cells=1, plot=False)
    adata.X = (adata.X.astype(np.float32) + 1.0) * 50.0
    adata = perform_normalization(adata)

    # Inject dummy PCA and UMAP embeddings
    dummy_pca = np.ones((adata.n_obs, 5), dtype=np.float32) * 42.0
    dummy_umap = np.ones((adata.n_obs, 2), dtype=np.float32) * 99.0
    adata.obsm["X_pca"] = dummy_pca.copy()
    adata.obsm["X_umap"] = dummy_umap.copy()

    # Run with force=True
    adata_out = perform_dimensionality_reduction_clustering(adata, force=True)

    # Embeddings should be recomputed and NOT equal to dummy
    assert not np.allclose(adata_out.obsm["X_pca"], dummy_pca)
    assert not np.allclose(adata_out.obsm["X_umap"], dummy_umap)
    assert "leiden" in adata_out.obs.columns


def test_perform_dimensionality_reduction_force_alias(sample_adata):
    """Test that force_dim_reduction kwarg works as an alias"""
    from genecircuitry.preprocessing import (
        perform_qc,
        perform_normalization,
        perform_dimensionality_reduction_clustering,
    )

    adata = perform_qc(sample_adata, min_genes=1, min_counts=1, min_cells=1, plot=False)
    adata.X = (adata.X.astype(np.float32) + 1.0) * 50.0
    adata = perform_normalization(adata)

    dummy_pca = np.ones((adata.n_obs, 5), dtype=np.float32) * 42.0
    adata.obsm["X_pca"] = dummy_pca.copy()

    adata_out = perform_dimensionality_reduction_clustering(
        adata, force_dim_reduction=True
    )
    assert not np.allclose(adata_out.obsm["X_pca"], dummy_pca)

