"""Unit and integration tests for multi-modal stratification (scRNA and scATAC)."""

import os
from datetime import datetime
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from genecircuitry.pipeline.controller import (
    PipelineController,
    stratification_pipeline,
)


@pytest.fixture
def sample_multimodal_data():
    """Create paired scRNA and scATAC AnnData objects with partial cell overlap."""
    # RNA has cells c1, c2, c3, c4, c5
    # ATAC has cells c2, c3, c4, c5, c6
    # Common cells: c2, c3, c4, c5
    rna_cells = ["c1", "c2", "c3", "c4", "c5"]
    rna_genes = ["geneA", "geneB", "geneC"]
    rna_X = np.ones((5, 3), dtype=np.float32)
    rna_obs = pd.DataFrame(
        {
            "cell_type": ["TypeA", "TypeA", "TypeB", "TypeB", "TypeA"],
            "batch": ["b1", "b1", "b2", "b2", "b1"],
        },
        index=rna_cells,
    )
    rna_adata = AnnData(X=rna_X, obs=rna_obs)

    atac_cells = ["c2", "c3", "c4", "c5", "c6"]
    atac_peaks = ["chr1_100_200", "chr1_300_400", "chr2_500_600"]
    atac_X = np.ones((5, 3), dtype=np.float32)
    atac_obs = pd.DataFrame(index=atac_cells)
    atac_adata = AnnData(X=atac_X, obs=atac_obs)

    return rna_adata, atac_adata


@pytest.fixture
def sample_unpaired_data():
    """Create unpaired scRNA and scATAC AnnData objects with zero cell overlap."""
    rna_cells = ["rna_1", "rna_2", "rna_3", "rna_4"]
    rna_X = np.ones((4, 3), dtype=np.float32)
    rna_obs = pd.DataFrame(
        {"cell_type": ["TypeA", "TypeA", "TypeB", "TypeB"]},
        index=rna_cells,
    )
    rna_adata = AnnData(X=rna_X, obs=rna_obs)

    atac_cells = ["atac_1", "atac_2", "atac_3"]
    atac_X = np.ones((3, 3), dtype=np.float32)
    atac_obs = pd.DataFrame(index=atac_cells)
    atac_adata = AnnData(X=atac_X, obs=atac_obs)

    return rna_adata, atac_adata


class TestMultiModalStratificationPipeline:
    """Tests for stratification_pipeline multi-modal behavior."""

    def test_stratification_without_scatac_backward_compatible(self, sample_multimodal_data):
        """Calling stratification_pipeline without scatac should return exactly 2 items."""
        rna_adata, _ = sample_multimodal_data
        result = stratification_pipeline(
            rna_adata,
            cluster_key_stratification="cell_type",
            clusters="all",
        )
        assert isinstance(result, tuple)
        assert len(result) == 2
        adata_list, names = result
        assert len(adata_list) == 2
        assert set(names) == {"TypeA", "TypeB"}

    def test_stratification_with_scatac_partial_overlap(self, sample_multimodal_data, capsys):
        """When scATAC is provided with partial overlap, filter to common cells and stratify both."""
        rna_adata, atac_adata = sample_multimodal_data
        result = stratification_pipeline(
            rna_adata,
            cluster_key_stratification="cell_type",
            clusters="all",
            scatac=atac_adata,
        )
        assert isinstance(result, tuple)
        assert len(result) == 3
        adata_list, names, scatac_list = result

        assert len(adata_list) == 2
        assert len(scatac_list) == 2
        assert set(names) == {"TypeA", "TypeB"}

        # Original RNA had c1 (dropped) and ATAC had c6 (dropped)
        # Common cells: c2, c3, c4, c5
        # TypeA: c2, c5 (c1 dropped) -> 2 cells
        # TypeB: c3, c4 -> 2 cells
        for rna_sub, atac_sub in zip(adata_list, scatac_list):
            assert list(rna_sub.obs_names) == list(atac_sub.obs_names)
            assert "c1" not in rna_sub.obs_names
            assert "c6" not in atac_sub.obs_names

        captured = capsys.readouterr().out
        assert "Multi-modal cell intersection" in captured or "common cells" in captured

    def test_stratification_with_scatac_zero_overlap_fallback(self, sample_unpaired_data, capsys):
        """When 0 common cells exist, do not stratify scATAC; log message and return None for scatac_list."""
        rna_adata, atac_adata = sample_unpaired_data
        result = stratification_pipeline(
            rna_adata,
            cluster_key_stratification="cell_type",
            clusters="all",
            scatac=atac_adata,
        )
        assert isinstance(result, tuple)
        assert len(result) == 3
        adata_list, names, scatac_list = result

        assert len(adata_list) == 2
        assert set(names) == {"TypeA", "TypeB"}
        assert scatac_list is None

        # RNA should retain all its cells because no paired filtering was applied
        assert sum(sub.n_obs for sub in adata_list) == rna_adata.n_obs

        captured = capsys.readouterr().out
        assert "0 common cell barcodes found" in captured or "not from a multiome assay" in captured

    def test_stratification_multi_key_with_scatac(self, sample_multimodal_data):
        """Multi-key compound stratification correctly stratifies both RNA and scATAC."""
        rna_adata, atac_adata = sample_multimodal_data
        result = stratification_pipeline(
            rna_adata,
            cluster_key_stratification=["cell_type", "batch"],
            clusters="all",
            scatac=atac_adata,
        )
        assert isinstance(result, tuple)
        assert len(result) == 3
        adata_list, names, scatac_list = result

        assert len(adata_list) > 0
        assert len(scatac_list) == len(adata_list)
        for rna_sub, atac_sub in zip(adata_list, scatac_list):
            assert list(rna_sub.obs_names) == list(atac_sub.obs_names)
            assert len(rna_sub.obs_names) > 0


class TestControllerMultiModalIntegration:
    """Tests for PipelineController integration with multi-modal stratification."""

    @patch("genecircuitry.pipeline.controller.dimensionality_reduction_clustering", side_effect=lambda ad, **kwargs: ad)
    @patch("genecircuitry.atac_peaks_processing.process_scatac_data")
    def test_controller_stratified_scatac_runs_per_stratification(
        self, mock_process_scatac, mock_dim_red, tmp_path, sample_multimodal_data
    ):
        """In paired multiome mode, global CIRCE is skipped and per-stratification CIRCE is run."""
        mock_process_scatac.return_value = "/mock/dir/enriched_scatac_peaks_dict.pkl"

        rna_adata, atac_adata = sample_multimodal_data
        # Pre-compute dummy PCA/UMAP
        rna_adata.obsm["X_pca"] = np.ones((rna_adata.n_obs, 2))
        rna_adata.obsm["X_umap"] = np.ones((rna_adata.n_obs, 2))

        output_dir = str(tmp_path / "output_test")
        os.makedirs(output_dir, exist_ok=True)

        args = MagicMock()
        args.output = output_dir
        args.input = None
        args.name = "test_run"
        args.species = "human"
        args.cluster_key = "cell_type"
        args.cluster_key_stratification = "cell_type"
        args.clusters = "all"
        args.scatac_data = atac_adata
        args.atac_peaks = None
        args.tf_dictionary = None
        args.no_base_grn = False
        args.keep_promoter_grn = False
        args.skip_celloracle = True
        args.skip_hotspot = True
        args.use_hvgs = False
        args.raw_layer_atac = None
        args.scatac_coaccess_threshold = "0.95q"
        args.scatac_metacells = False
        args.n_jobs = 1

        controller = PipelineController(args, datetime.now())
        controller.adata_preprocessed = rna_adata.copy()

        # 1. Stratification step
        adata_list, name_list = controller.run_step_stratification()
        assert len(adata_list) == 2
        assert len(controller.scatac_list) == 2

        # 2. ATAC peaks step (global) should be skipped because stratified scATAC is active
        global_result = controller.run_step_atac_peaks()
        assert global_result is None
        mock_process_scatac.assert_not_called()

        # 3. Process single stratification
        strat_dir = controller.process_single_stratification(
            adata_list[0], name_list[0], scatac_cluster=controller.scatac_list[0]
        )
        assert strat_dir is not None
        # process_scatac_data should be called for this stratification
        mock_process_scatac.assert_called_once()
        _, call_kwargs = mock_process_scatac.call_args
        assert call_kwargs["scatac_data"].n_obs == adata_list[0].n_obs
        assert os.path.join(output_dir, "stratified_analysis", name_list[0]) in call_kwargs["output_dir"]

    @patch("genecircuitry.atac_peaks_processing.process_scatac_data")
    def test_controller_unpaired_scatac_runs_global_circe(
        self, mock_process_scatac, tmp_path, sample_unpaired_data
    ):
        """In unpaired mode (0 common cells), global CIRCE runs once and supplies shared base GRN."""
        mock_process_scatac.return_value = "/mock/dir/enriched_scatac_peaks_dict.pkl"

        rna_adata, atac_adata = sample_unpaired_data
        rna_adata.obsm["X_pca"] = np.ones((rna_adata.n_obs, 2))
        rna_adata.obsm["X_umap"] = np.ones((rna_adata.n_obs, 2))

        output_dir = str(tmp_path / "output_test_unpaired")
        os.makedirs(output_dir, exist_ok=True)

        args = MagicMock()
        args.output = output_dir
        args.input = None
        args.name = "test_run_unpaired"
        args.species = "human"
        args.cluster_key = "cell_type"
        args.cluster_key_stratification = "cell_type"
        args.clusters = "all"
        args.scatac_data = atac_adata
        args.atac_peaks = None
        args.tf_dictionary = None
        args.no_base_grn = False
        args.keep_promoter_grn = False
        args.skip_celloracle = True
        args.skip_hotspot = True
        args.use_hvgs = False
        args.raw_layer_atac = None
        args.scatac_coaccess_threshold = "0.95q"
        args.scatac_metacells = False
        args.n_jobs = 1

        controller = PipelineController(args, datetime.now())
        controller.adata_preprocessed = rna_adata.copy()

        # 1. Stratification step
        adata_list, name_list = controller.run_step_stratification()
        assert len(adata_list) == 2
        assert len(controller.scatac_list) == 0  # Unpaired -> not stratified

        # 2. ATAC peaks step (global) should NOT be skipped because scatac_list is empty
        global_result = controller.run_step_atac_peaks()
        assert global_result == "/mock/dir/enriched_scatac_peaks_dict.pkl"
        assert controller.atac_peaks_pkl == "/mock/dir/enriched_scatac_peaks_dict.pkl"
        mock_process_scatac.assert_called_once()
