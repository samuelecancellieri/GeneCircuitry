"""Unit tests for ATAC peaks processing and checkpointing."""

import os
import pickle
from datetime import datetime
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from genecircuitry.atac_peaks_processing import process_atac_peaks
from genecircuitry.pipeline.controller import (
    PipelineController,
    check_checkpoint,
    compute_input_hash,
    write_checkpoint,
)


@pytest.fixture
def dummy_bed_file(tmp_path):
    """Create a dummy BED file for testing."""
    bed_content = (
        "chr1\t1000\t2000\tpeak_1\n"
        "chr1\t3000\t4000\tpeak_2\n"
        "chr2\t5000\t6000\tpeak_3\n"
    )
    bed_file = tmp_path / "test_peaks.bed"
    bed_file.write_text(bed_content)
    return str(bed_file)


class TestAtacPeaksCheckpoint:
    """Test checkpointing in ATAC peaks processing."""

    @patch("celloracle.motif_analysis.TFinfo")
    @patch("celloracle.motif_analysis.check_peak_format")
    @patch("genecircuitry.atac_peaks_processing._annotate_bed_peaks")
    @patch("genecircuitry.atac_peaks_processing._ensure_genome_installed")
    def test_process_atac_peaks_creates_files_and_checkpoint(
        self,
        mock_ensure_genome,
        mock_annotate,
        mock_check_format,
        mock_tfinfo_cls,
        dummy_bed_file,
        tmp_path,
    ):
        """Test process_atac_peaks creates df and dict pkl files and saves checkpoint."""
        mock_annotate.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000", "chr1_3000_4000"], "gene_short_name": ["GeneA", "GeneB"]}
        )
        mock_check_format.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000", "chr1_3000_4000"], "gene_short_name": ["GeneA", "GeneB"]}
        )

        mock_tfi_inst = MagicMock()
        mock_tfinfo_cls.return_value = mock_tfi_inst
        dummy_df = pd.DataFrame({"TF1": [1, 0], "TF2": [0, 1]})
        dummy_dict = {"GeneA": ["TF1"], "GeneB": ["TF2"]}
        mock_tfi_inst.to_dataframe.return_value = dummy_df
        mock_tfi_inst.to_dictionary.return_value = dummy_dict

        output_dir = tmp_path / "output"
        log_dir = tmp_path / "output" / "logs"

        pkl_path = process_atac_peaks(
            bed_path=dummy_bed_file,
            species="human",
            output_dir=str(output_dir),
            log_dir=str(log_dir),
        )

        assert pkl_path.endswith("enriched_atac_peaks_dict.pkl")
        assert os.path.exists(pkl_path)
        df_pkl_path = os.path.join(str(output_dir), "celloracle", "enriched_atac_peaks_df.pkl")
        assert os.path.exists(df_pkl_path)

        # Checkpoint file should exist in log_dir
        checkpoint_file = log_dir / "atac_peaks.checkpoint"
        assert checkpoint_file.exists()

        # Checkpoint should be valid
        step_hash = compute_input_hash(dummy_bed_file, species="human", fpr=0.02, threshold=10)
        assert check_checkpoint(str(log_dir), "atac_peaks", step_hash) is True

    @patch("celloracle.motif_analysis.TFinfo")
    @patch("genecircuitry.atac_peaks_processing._annotate_bed_peaks")
    @patch("genecircuitry.atac_peaks_processing._ensure_genome_installed")
    def test_process_atac_peaks_skips_when_checkpoint_exists(
        self,
        mock_ensure_genome,
        mock_annotate,
        mock_tfinfo_cls,
        dummy_bed_file,
        tmp_path,
    ):
        """Test process_atac_peaks skips motif analysis when checkpoint and output files exist."""
        output_dir = tmp_path / "output"
        celloracle_dir = output_dir / "celloracle"
        celloracle_dir.mkdir(parents=True)
        log_dir = output_dir / "logs"
        log_dir.mkdir(parents=True)

        dict_path = celloracle_dir / "enriched_atac_peaks_dict.pkl"
        df_path = celloracle_dir / "enriched_atac_peaks_df.pkl"

        dummy_dict = {"GeneA": ["TF1"]}
        dummy_df = pd.DataFrame({"TF1": [1]})
        with open(dict_path, "wb") as f:
            pickle.dump(dummy_dict, f)
        dummy_df.to_pickle(df_path)

        step_hash = compute_input_hash(dummy_bed_file, species="human", fpr=0.02, threshold=10)
        write_checkpoint(str(log_dir), "atac_peaks", step_hash, pkl_path=str(dict_path))

        # Call process_atac_peaks
        result_pkl = process_atac_peaks(
            bed_path=dummy_bed_file,
            species="human",
            output_dir=str(output_dir),
            log_dir=str(log_dir),
            force=False,
        )

        # Should return existing path without invoking annotation or motif scanning
        assert result_pkl == str(dict_path)
        mock_annotate.assert_not_called()
        mock_tfinfo_cls.assert_not_called()

    @patch("celloracle.motif_analysis.TFinfo")
    @patch("celloracle.motif_analysis.check_peak_format")
    @patch("genecircuitry.atac_peaks_processing._annotate_bed_peaks")
    @patch("genecircuitry.atac_peaks_processing._ensure_genome_installed")
    def test_process_atac_peaks_force_reruns(
        self,
        mock_ensure_genome,
        mock_annotate,
        mock_check_format,
        mock_tfinfo_cls,
        dummy_bed_file,
        tmp_path,
    ):
        """Test process_atac_peaks re-runs when force=True even if files exist."""
        output_dir = tmp_path / "output"
        celloracle_dir = output_dir / "celloracle"
        celloracle_dir.mkdir(parents=True)
        log_dir = output_dir / "logs"
        log_dir.mkdir(parents=True)

        dict_path = celloracle_dir / "enriched_atac_peaks_dict.pkl"
        df_path = celloracle_dir / "enriched_atac_peaks_df.pkl"
        with open(dict_path, "wb") as f:
            pickle.dump({"GeneA": ["TF1"]}, f)
        pd.DataFrame({"TF1": [1]}).to_pickle(df_path)

        mock_annotate.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000"], "gene_short_name": ["GeneA"]}
        )
        mock_check_format.return_value = pd.DataFrame(
            {"peak_id": ["chr1_1000_2000"], "gene_short_name": ["GeneA"]}
        )
        mock_tfi_inst = MagicMock()
        mock_tfinfo_cls.return_value = mock_tfi_inst
        mock_tfi_inst.to_dataframe.return_value = pd.DataFrame({"TF2": [1]})
        mock_tfi_inst.to_dictionary.return_value = {"GeneA": ["TF2"]}

        result_pkl = process_atac_peaks(
            bed_path=dummy_bed_file,
            species="human",
            output_dir=str(output_dir),
            log_dir=str(log_dir),
            force=True,
        )

        assert result_pkl == str(dict_path)
        mock_annotate.assert_called_once()
        mock_tfi_inst.scan.assert_called_once()

    @patch("genecircuitry.atac_peaks_processing.process_atac_peaks")
    def test_pipeline_controller_run_step_atac_peaks_checkpoint(
        self,
        mock_process,
        dummy_bed_file,
        tmp_path,
    ):
        """Test PipelineController.run_step_atac_peaks uses checkpoint."""
        output_dir = tmp_path / "output_controller"
        celloracle_dir = output_dir / "celloracle"
        celloracle_dir.mkdir(parents=True)
        log_dir = output_dir / "logs"
        log_dir.mkdir(parents=True)

        dict_path = celloracle_dir / "enriched_atac_peaks_dict.pkl"
        with open(dict_path, "wb") as f:
            pickle.dump({"GeneA": ["TF1"]}, f)

        step_hash = compute_input_hash(dummy_bed_file, species="human", fpr=0.02, threshold=10)
        write_checkpoint(str(log_dir), "atac_peaks", step_hash, pkl_path=str(dict_path))

        args = SimpleNamespace(
            output=str(output_dir),
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
            atac_peaks=dummy_bed_file,
            no_base_grn=False,
            min_genes=1,
            min_counts=1,
            seed=42,
            n_jobs=1,
            skip_qc=False,
            skip_celloracle=False,
            skip_hotspot=False,
            use_hvgs=False,
            debug=False,
            steps=None,
            force_dim_reduction=False,
        )

        controller = PipelineController(args, datetime.now())
        res = controller.run_step_atac_peaks()

        assert res == str(dict_path)
        assert controller.atac_peaks_pkl == str(dict_path)
        # process_atac_peaks was NOT called because checkpoint hit
        mock_process.assert_not_called()
