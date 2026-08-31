"""
Tests for PipelineController dimensionality reduction force flag and stratification re-clustering.
"""

import os
import pytest
import numpy as np
from types import SimpleNamespace
from datetime import datetime
from anndata import AnnData

from genecircuitry.preprocessing import perform_qc, perform_normalization
from genecircuitry.pipeline.controller import (
    PipelineController,
    dimensionality_reduction_clustering,
    stratification_pipeline,
)


def make_default_args(**kwargs):
    """Helper to create complete SimpleNamespace for PipelineController tests"""
    default_kwargs = dict(
        output="output",
        name="test_run",
        species="human",
        cluster_key="leiden",
        clusters="all",
        cluster_key_stratification=None,
        embedding_grn="X_draw_graph_fa",
        embedding_hotspot="X_umap",
        normalization_key="n_counts",
        raw_count_layer="raw_counts",
        tf_dictionary=None,
        atac_peaks=None,
        no_base_grn=False,
        min_genes=1,
        min_counts=1,
        seed=42,
        n_jobs=1,
        skip_qc=False,
        skip_celloracle=True,
        skip_hotspot=True,
        use_hvgs=False,
        debug=False,
        steps=None,
        force_dim_reduction=False,
    )
    default_kwargs.update(kwargs)
    return SimpleNamespace(**default_kwargs)


@pytest.fixture
def sample_adata_for_clustering():
    """Create a sample AnnData object for testing controller clustering"""
    np.random.seed(42)
    n_cells = 60
    n_genes = 30

    X = (
        np.random.negative_binomial(5, 0.3, (n_cells, n_genes)).astype(np.float32) + 1.0
    ) * 50.0
    var_names = [f"Gene_{i}" for i in range(n_genes)]
    obs_names = [f"Cell_{i}" for i in range(n_cells)]

    adata = AnnData(X=X)
    adata.var_names = var_names
    adata.obs_names = obs_names

    # Add stratification label: half 'TypeA', half 'TypeB'
    adata.obs["cell_type"] = ["TypeA"] * 30 + ["TypeB"] * 30

    adata = perform_qc(adata, min_genes=1, min_counts=1, min_cells=1, plot=False)
    adata = perform_normalization(adata)
    return adata


def test_dimensionality_reduction_clustering_force(tmp_path, sample_adata_for_clustering):
    """Test dimensionality_reduction_clustering with force=True bypassing checkpoints and pre-existing embeddings"""
    adata = sample_adata_for_clustering.copy()
    dummy_pca = np.ones((adata.n_obs, 5), dtype=np.float32) * 55.0
    adata.obsm["X_pca"] = dummy_pca.copy()

    # With force=False, dummy PCA is kept
    adata_noforce = dimensionality_reduction_clustering(adata.copy(), force=False)
    np.testing.assert_array_equal(adata_noforce.obsm["X_pca"], dummy_pca)

    # With force=True, PCA is recalculated
    adata_force = dimensionality_reduction_clustering(adata.copy(), force=True)
    assert not np.allclose(adata_force.obsm["X_pca"], dummy_pca)


def test_pipeline_controller_force_flag(tmp_path, sample_adata_for_clustering):
    """Test that PipelineController respects force_dim_reduction flag"""
    adata = sample_adata_for_clustering.copy()
    dummy_pca = np.ones((adata.n_obs, 5), dtype=np.float32) * 77.0
    adata.obsm["X_pca"] = dummy_pca.copy()

    output_dir = str(tmp_path / "output_test")
    os.makedirs(output_dir, exist_ok=True)

    args = make_default_args(
        output=output_dir,
        name="test_run",
        cluster_key="leiden",
        cluster_key_stratification=None,
        clusters="all",
        force_dim_reduction=True,
    )

    controller = PipelineController(args, datetime.now())
    controller.adata_preprocessed = adata

    adata_clustered = controller.run_step_clustering()
    assert not np.allclose(adata_clustered.obsm["X_pca"], dummy_pca)


def test_stratification_forces_dim_reduction(tmp_path, sample_adata_for_clustering):
    """Test that stratified analysis forces re-run of dimensionality reduction for each stratification adata"""
    adata = sample_adata_for_clustering.copy()

    # Pre-compute global embeddings on full dataset
    dummy_pca = np.ones((adata.n_obs, 5), dtype=np.float32) * 88.0
    dummy_umap = np.ones((adata.n_obs, 2), dtype=np.float32) * 99.0
    adata.obsm["X_pca"] = dummy_pca.copy()
    adata.obsm["X_umap"] = dummy_umap.copy()

    output_dir = str(tmp_path / "output_strat_test")
    os.makedirs(output_dir, exist_ok=True)

    args = make_default_args(
        output=output_dir,
        name="strat_test",
        cluster_key="leiden",
        cluster_key_stratification="cell_type",
        clusters="all",
        force_dim_reduction=False,
        skip_celloracle=True,
        skip_hotspot=True,
    )

    controller = PipelineController(args, datetime.now())
    controller.adata_preprocessed = adata

    # Run stratification step
    adata_list, name_list = controller.run_step_stratification()
    assert len(adata_list) == 2

    # Process first stratification
    strat_dir = controller.process_single_stratification(adata_list[0], name_list[0])
    assert strat_dir is not None

    # Check stratification results stored in controller
    assert len(controller.stratification_results) == 1
    strat_result = controller.stratification_results[0]
    strat_adata = strat_result["adata"]

    # The stratified adata's PCA and UMAP should have been recomputed, NOT equal to dummy
    assert not np.allclose(strat_adata.obsm["X_pca"], np.ones((strat_adata.n_obs, 5)) * 88.0)
    assert not np.allclose(strat_adata.obsm["X_umap"], np.ones((strat_adata.n_obs, 2)) * 99.0)

