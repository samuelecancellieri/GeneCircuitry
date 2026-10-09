"""Unit tests for scATAC-seq CIRCE co-accessibility integration and CellOracle pipeline."""

import os
import pickle
import sys
from datetime import datetime
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

# Ensure optional dependencies can be safely mocked if not installed in the test environment
mock_celloracle = MagicMock()
mock_ma = MagicMock()
mock_celloracle.motif_analysis = mock_ma
sys.modules.setdefault("celloracle", mock_celloracle)
sys.modules.setdefault("celloracle.motif_analysis", mock_ma)

mock_gimme = MagicMock()
mock_gimme_motif = MagicMock()
mock_gimme.motif = mock_gimme_motif
sys.modules.setdefault("gimmemotifs", mock_gimme)
sys.modules.setdefault("gimmemotifs.motif", mock_gimme_motif)

mock_circe = MagicMock()
sys.modules.setdefault("circe", mock_circe)

from genecircuitry.atac_peaks_processing import (
    _load_scatac_data,
    _parse_threshold,
    process_scatac_data,
)
from genecircuitry.pipeline.controller import (
    PipelineController,
    check_checkpoint,
    compute_input_hash,
    create_parser,
)


@pytest.fixture
def dummy_scatac_adata():
    """Create a dummy scATAC AnnData object with peak coordinates."""
    n_cells = 20
    peaks = ["chr1:1000-2000", "chr1:3000-4000", "chr2:5000-6000"]
    X = np.random.randint(0, 5, size=(n_cells, len(peaks)))
    adata = AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=peaks),
    )
    return adata


class TestLoadScAtacData:
    """Test loading scATAC data from various formats."""

    def test_in_memory_anndata(self, dummy_scatac_adata):
        """Test passing an in-memory AnnData object returns a copy."""
        loaded = _load_scatac_data(dummy_scatac_adata)
        assert isinstance(loaded, AnnData)
        assert list(loaded.var_names) == list(dummy_scatac_adata.var_names)
        assert loaded is not dummy_scatac_adata

    def test_h5ad_file(self, dummy_scatac_adata, tmp_path):
        """Test loading .h5ad file."""
        h5ad_path = tmp_path / "scatac.h5ad"
        dummy_scatac_adata.write_h5ad(str(h5ad_path))

        loaded = _load_scatac_data(str(h5ad_path))
        assert isinstance(loaded, AnnData)
        assert len(loaded.var_names) == len(dummy_scatac_adata.var_names)

    def test_h5mu_file_with_atac_modality(self, dummy_scatac_adata, tmp_path):
        """Test loading .h5mu extracts the 'atac' modality."""
        h5mu_path = tmp_path / "multiome.h5mu"
        h5mu_path.touch()

        mock_mudata = MagicMock()
        mock_mudata.mod = {"atac": dummy_scatac_adata, "rna": MagicMock()}

        with patch.dict("sys.modules", {"muon": MagicMock()}):
            import muon as mu

            mu.read_h5mu.return_value = mock_mudata
            loaded = _load_scatac_data(str(h5mu_path))
            assert isinstance(loaded, AnnData)
            assert list(loaded.var_names) == list(dummy_scatac_adata.var_names)

    def test_h5mu_file_missing_atac_modality_raises(self, tmp_path):
        """Test .h5mu without 'atac' modality raises ValueError when multiple modalities exist."""
        h5mu_path = tmp_path / "multiome.h5mu"
        h5mu_path.touch()

        mock_mudata = MagicMock()
        mock_mudata.mod = {"rna": MagicMock(), "protein": MagicMock()}

        with patch.dict("sys.modules", {"muon": MagicMock()}):
            import muon as mu

            mu.read_h5mu.return_value = mock_mudata
            with pytest.raises(ValueError, match="no 'atac' modality found"):
                _load_scatac_data(str(h5mu_path))

    def test_file_not_found(self):
        """Test non-existent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            _load_scatac_data("non_existent_file.h5ad")

    def test_invalid_type(self):
        """Test invalid input type raises TypeError."""
        with pytest.raises(TypeError, match="Unsupported type"):
            _load_scatac_data(12345)

    def test_auto_detect_raw_layer_counts(self, dummy_scatac_adata):
        """Test auto-detecting and extracting raw counts from layers['counts']."""
        raw_counts = np.array([[10, 20, 30]] * dummy_scatac_adata.n_obs)
        dummy_scatac_adata.layers["counts"] = raw_counts.copy()
        dummy_scatac_adata.X = np.log1p(raw_counts)

        loaded = _load_scatac_data(dummy_scatac_adata)
        np.testing.assert_array_equal(loaded.X, raw_counts)

    def test_auto_detect_raw_layer_raw_counts(self, dummy_scatac_adata):
        """Test auto-detecting and extracting raw counts from layers['raw_counts']."""
        raw_counts = np.array([[15, 25, 35]] * dummy_scatac_adata.n_obs)
        dummy_scatac_adata.layers["raw_counts"] = raw_counts.copy()
        dummy_scatac_adata.X = np.log1p(raw_counts)

        loaded = _load_scatac_data(dummy_scatac_adata)
        np.testing.assert_array_equal(loaded.X, raw_counts)

    def test_auto_detect_dot_raw(self, dummy_scatac_adata):
        """Test auto-detecting and extracting raw counts from adata.raw."""
        raw_counts = np.array([[100, 200, 300]] * dummy_scatac_adata.n_obs)
        dummy_scatac_adata.X = raw_counts.copy()
        dummy_scatac_adata.raw = dummy_scatac_adata
        dummy_scatac_adata.X = np.log1p(raw_counts)

        loaded = _load_scatac_data(dummy_scatac_adata)
        np.testing.assert_array_equal(loaded.X, raw_counts)

    def test_explicit_raw_layer_atac(self, dummy_scatac_adata):
        """Test extracting raw counts from explicitly specified raw_layer_atac."""
        custom_raw = np.array([[50, 60, 70]] * dummy_scatac_adata.n_obs)
        dummy_scatac_adata.layers["my_raw_layer"] = custom_raw.copy()
        dummy_scatac_adata.X = np.zeros_like(custom_raw)

        loaded = _load_scatac_data(dummy_scatac_adata, raw_layer_atac="my_raw_layer")
        np.testing.assert_array_equal(loaded.X, custom_raw)

    def test_explicit_raw_layer_raw_keyword(self, dummy_scatac_adata):
        """Test specifying raw_layer_atac='raw' explicitly extracts from adata.raw."""
        raw_counts = np.array([[77, 88, 99]] * dummy_scatac_adata.n_obs)
        dummy_scatac_adata.X = raw_counts.copy()
        dummy_scatac_adata.raw = dummy_scatac_adata
        dummy_scatac_adata.X = np.zeros_like(raw_counts)

        loaded = _load_scatac_data(dummy_scatac_adata, raw_layer_atac="raw")
        np.testing.assert_array_equal(loaded.X, raw_counts)

    def test_explicit_raw_layer_not_found_raises(self, dummy_scatac_adata):
        """Test specifying non-existent raw layer raises ValueError."""
        with pytest.raises(ValueError, match="Specified raw layer 'non_existent' not found"):
            _load_scatac_data(dummy_scatac_adata, raw_layer_atac="non_existent")

    def test_fallback_to_X_when_no_raw_present(self, dummy_scatac_adata):
        """Test keeping adata.X when neither raw layer nor .raw is found."""
        orig_X = dummy_scatac_adata.X.copy()
        loaded = _load_scatac_data(dummy_scatac_adata)
        np.testing.assert_array_equal(loaded.X, orig_X)


class TestParseThreshold:
    """Test co-accessibility threshold parsing."""

    def test_quantile_syntax(self):
        """Test '0.95q' and '95%' syntax computes quantile."""
        series = pd.Series([0.1, 0.2, 0.5, 0.8, 0.9, 1.0])
        val_q = _parse_threshold(series, "0.95q")
        val_pct = _parse_threshold(series, "95%")
        expected = float(series.quantile(0.95))
        assert val_q == pytest.approx(expected)
        assert val_pct == pytest.approx(expected)

    def test_absolute_score(self):
        """Test numeric score threshold."""
        series = pd.Series([0.1, 0.2, 0.5, 0.8])
        assert _parse_threshold(series, 0.8) == 0.8
        assert _parse_threshold(series, "0.8") == 0.8

    def test_invalid_spec_raises(self):
        """Test invalid threshold specification raises ValueError."""
        series = pd.Series([0.1, 0.2])
        with pytest.raises(ValueError, match="Invalid co-accessibility threshold"):
            _parse_threshold(series, "invalid_threshold")


class TestProcessScAtacData:
    """Test process_scatac_data end-to-end with CIRCE and CellOracle mocking."""

    @patch("celloracle.motif_analysis.TFinfo")
    @patch("celloracle.motif_analysis.check_peak_format")
    @patch("celloracle.motif_analysis.integrate_tss_peak_with_cicero")
    @patch("celloracle.motif_analysis.get_tss_info")
    @patch("genecircuitry.atac_peaks_processing._ensure_genome_installed")
    def test_process_scatac_data_creates_files_and_checkpoint(
        self,
        mock_ensure_genome,
        mock_get_tss,
        mock_integrate,
        mock_check_format,
        mock_tfinfo_cls,
        dummy_scatac_adata,
        tmp_path,
    ):
        """Test process_scatac_data performs full workflow and writes checkpoint."""
        mock_ci = MagicMock()
        mock_ci.extract_atac_links.return_value = pd.DataFrame(
            {
                "Peak1": ["chr1_1000_2000", "chr1_3000_4000"],
                "Peak2": ["chr1_3000_4000", "chr2_5000_6000"],
                "score": [0.85, 0.92],
            }
        )

        mock_get_tss.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000"], "gene_short_name": ["GeneA"]}
        )
        mock_integrate.return_value = pd.DataFrame(
            {
                "peak_id": ["chr1_1000_2000", "chr1_3000_4000"],
                "gene_short_name": ["GeneA", "GeneB"],
                "coaccess": [0.85, 0.92],
            }
        )
        mock_check_format.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000", "chr1_3000_4000"], "gene_short_name": ["GeneA", "GeneB"]}
        )

        mock_tfi = MagicMock()
        mock_tfinfo_cls.return_value = mock_tfi
        dummy_df = pd.DataFrame({"TF1": [1, 0], "TF2": [0, 1]})
        dummy_dict = {"GeneA": ["TF1"], "GeneB": ["TF2"]}
        mock_tfi.to_dataframe.return_value = dummy_df
        mock_tfi.to_dictionary.return_value = dummy_dict

        output_dir = tmp_path / "output"
        log_dir = tmp_path / "output" / "logs"

        with patch.dict("sys.modules", {"circe": mock_ci}):
            dict_pkl_path = process_scatac_data(
                scatac_data=dummy_scatac_adata,
                species="human",
                output_dir=str(output_dir),
                coaccess_threshold="0.95q",
                compute_metacells=False,
                log_dir=str(log_dir),
            )

        assert dict_pkl_path.endswith("enriched_scatac_peaks_dict.pkl")
        assert os.path.exists(dict_pkl_path)
        df_pkl_path = os.path.join(str(output_dir), "celloracle", "enriched_scatac_peaks_df.pkl")
        assert os.path.exists(df_pkl_path)

        # Check CIRCE link file
        links_csv = os.path.join(str(output_dir), "celloracle", "scatac_coaccessibility_links.csv")
        assert os.path.exists(links_csv)

        # Checkpoint file exists
        checkpoint_file = log_dir / "scatac_data.checkpoint"
        assert checkpoint_file.exists()

        # Metacells was not called because compute_metacells=False
        mock_ci.metacells.compute_metacells.assert_not_called()

    @patch("celloracle.motif_analysis.TFinfo")
    @patch("celloracle.motif_analysis.check_peak_format")
    @patch("celloracle.motif_analysis.integrate_tss_peak_with_cicero")
    @patch("celloracle.motif_analysis.get_tss_info")
    @patch("genecircuitry.atac_peaks_processing._ensure_genome_installed")
    def test_process_scatac_data_with_metacells(
        self,
        mock_ensure_genome,
        mock_get_tss,
        mock_integrate,
        mock_check_format,
        mock_tfinfo_cls,
        dummy_scatac_adata,
        tmp_path,
    ):
        """Test compute_metacells=True invokes CIRCE metacell aggregation."""
        mock_ci = MagicMock()
        mock_ci.metacells.compute_metacells.return_value = dummy_scatac_adata
        mock_ci.extract_atac_links.return_value = pd.DataFrame(
            {"Peak1": ["chr1_1000_2000"], "Peak2": ["chr1_3000_4000"], "score": [0.9]}
        )
        mock_integrate.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000"], "gene_short_name": ["GeneA"], "coaccess": [0.9]}
        )
        mock_check_format.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000"], "gene_short_name": ["GeneA"]}
        )

        mock_tfi = MagicMock()
        mock_tfinfo_cls.return_value = mock_tfi
        mock_tfi.to_dataframe.return_value = pd.DataFrame({"TF1": [1]})
        mock_tfi.to_dictionary.return_value = {"GeneA": ["TF1"]}

        output_dir = tmp_path / "output_mc"

        with patch.dict("sys.modules", {"circe": mock_ci}):
            process_scatac_data(
                scatac_data=dummy_scatac_adata,
                species="human",
                output_dir=str(output_dir),
                compute_metacells=True,
            )

        mock_ci.metacells.compute_metacells.assert_called_once()

    @patch("genecircuitry.atac_peaks_processing._ensure_genome_installed")
    def test_process_scatac_data_checkpoint_skip(
        self,
        mock_ensure_genome,
        dummy_scatac_adata,
        tmp_path,
    ):
        """Test that existing checkpoint and files skip processing."""
        output_dir = tmp_path / "output"
        celloracle_dir = output_dir / "celloracle"
        celloracle_dir.mkdir(parents=True)
        log_dir = output_dir / "logs"
        log_dir.mkdir(parents=True)

        dict_path = celloracle_dir / "enriched_scatac_peaks_dict.pkl"
        df_path = celloracle_dir / "enriched_scatac_peaks_df.pkl"

        dummy_dict = {"GeneA": ["TF1"]}
        dummy_df = pd.DataFrame({"TF1": [1]})
        with open(dict_path, "wb") as f:
            pickle.dump(dummy_dict, f)
        dummy_df.to_pickle(df_path)

        input_repr = str(dummy_scatac_adata.shape)
        step_hash = compute_input_hash(
            input_repr,
            species="human",
            fpr=0.02,
            threshold=10,
            coaccess_threshold="0.95q",
            compute_metacells=False,
        )
        from genecircuitry.pipeline.controller import write_checkpoint

        write_checkpoint(str(log_dir), "scatac_data", step_hash, pkl_path=str(dict_path))

        # Call process_scatac_data without CIRCE mock — should return early from checkpoint
        result_pkl = process_scatac_data(
            scatac_data=dummy_scatac_adata,
            species="human",
            output_dir=str(output_dir),
            log_dir=str(log_dir),
            coaccess_threshold="0.95q",
            compute_metacells=False,
            force=False,
        )

        assert result_pkl == str(dict_path)


class TestControllerScAtacIntegration:
    """Test CLI argument parsing and Controller execution with scATAC."""

    def test_cli_parser_has_scatac_args(self):
        """Test that CLI parser registers all scATAC arguments."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "--scatac-data",
                "test.h5ad",
                "--raw-layer-atac",
                "counts",
                "--scatac-coaccess-threshold",
                "0.90q",
                "--scatac-metacells",
                "--keep-promoter-grn",
            ]
        )
        assert args.scatac_data == "test.h5ad"
        assert args.raw_layer_atac == "counts"
        assert args.scatac_coaccess_threshold == "0.90q"
        assert args.scatac_metacells is True
        assert args.keep_promoter_grn is True

    @patch("genecircuitry.atac_peaks_processing.process_scatac_data")
    def test_controller_runs_scatac_step(self, mock_process_scatac, tmp_path):
        """Test Controller executes process_scatac_data and sets no_base_grn=True by default."""
        mock_process_scatac.return_value = "/path/to/enriched_scatac_peaks_dict.pkl"

        parser = create_parser()
        args = parser.parse_args([
            "--scatac-data", "data.h5ad",
            "--raw-layer-atac", "counts",
            "--output", str(tmp_path),
            "--scatac-coaccess-threshold", "0.95q",
        ])

        controller = PipelineController(args, datetime.now())
        controller.log_dir = str(tmp_path / "logs")

        dict_path = controller.run_step_atac_peaks()

        assert dict_path == "/path/to/enriched_scatac_peaks_dict.pkl"
        assert controller.atac_peaks_pkl == "/path/to/enriched_scatac_peaks_dict.pkl"
        # Since keep_promoter_grn is False, no_base_grn should be set to True
        assert controller.args.no_base_grn is True
        mock_process_scatac.assert_called_once()
        _, call_kwargs = mock_process_scatac.call_args
        assert call_kwargs["raw_layer_atac"] == "counts"

    @patch("genecircuitry.atac_peaks_processing.process_scatac_data")
    def test_controller_scatac_with_keep_promoter_grn(self, mock_process_scatac, tmp_path):
        """Test that --keep-promoter-grn preserves no_base_grn=False."""
        mock_process_scatac.return_value = "/path/to/enriched_scatac_peaks_dict.pkl"

        parser = create_parser()
        args = parser.parse_args([
            "--scatac-data", "data.h5ad",
            "--output", str(tmp_path),
            "--keep-promoter-grn",
        ])

        controller = PipelineController(args, datetime.now())
        controller.log_dir = str(tmp_path / "logs")

        controller.run_step_atac_peaks()
        # With keep_promoter_grn=True, no_base_grn remains False
        assert controller.args.no_base_grn is False
