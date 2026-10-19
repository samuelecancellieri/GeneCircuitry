"""Tests for CellOracle in-silico perturbation integration in PipelineController."""

import os
import pytest
from unittest.mock import MagicMock, patch
import pandas as pd
from genecircuitry.pipeline.controller import create_parser, PipelineController


def test_create_parser_perturbation_args():
    """Test that create_parser has perturbation CLI arguments with correct defaults."""
    parser = create_parser()

    # Test defaults
    args = parser.parse_args([])
    assert hasattr(args, "run_perturbation")
    assert args.run_perturbation is False
    assert hasattr(args, "compare_clusters")
    assert args.compare_clusters is None
    assert hasattr(args, "top_n_perturbed_tfs")
    assert args.top_n_perturbed_tfs == 3
    assert hasattr(args, "perturbation_mode")
    assert args.perturbation_mode == "directional"
    assert hasattr(args, "n_propagation")
    assert args.n_propagation == 3
    assert hasattr(args, "pseudotime_key")
    assert args.pseudotime_key == "pseudotime"

    # Test custom arguments
    custom_args = parser.parse_args([
        "--run-perturbation",
        "--compare-clusters", "MEP,GMP;Ery,Gran",
        "--top-n-perturbed-tfs", "5",
        "--perturbation-mode", "ko",
        "--n-propagation", "4",
        "--pseudotime-key", "dpt_pseudotime",
    ])
    assert custom_args.run_perturbation is True
    assert custom_args.compare_clusters == "MEP,GMP;Ery,Gran"
    assert custom_args.top_n_perturbed_tfs == 5
    assert custom_args.perturbation_mode == "ko"
    assert custom_args.n_propagation == 4
    assert custom_args.pseudotime_key == "dpt_pseudotime"


def test_controller_parse_compare_clusters():
    """Test parsing string cluster pairs in PipelineController."""
    controller = PipelineController.__new__(PipelineController)
    assert controller.parse_compare_clusters(None) is None
    assert controller.parse_compare_clusters("") is None
    assert controller.parse_compare_clusters("MEP,GMP") == [("MEP", "GMP")]
    assert controller.parse_compare_clusters("MEP,GMP;Ery,Gran") == [("MEP", "GMP"), ("Ery", "Gran")]
    assert controller.parse_compare_clusters(" MEP , GMP ; Ery , Gran ") == [("MEP", "GMP"), ("Ery", "Gran")]


def test_controller_run_step_perturbation_calls_analysis(tmp_path):
    """Test that controller.run_step_perturbation wires inputs and calls run_perturbation_analysis."""
    out_dir = str(tmp_path / "output")
    os.makedirs(out_dir, exist_ok=True)

    args = MagicMock()
    args.output = out_dir
    args.cluster_key = "celltype"
    args.run_perturbation = True
    args.compare_clusters = "MEP,GMP"
    args.top_n_perturbed_tfs = 2
    args.perturbation_mode = "directional"
    args.n_propagation = 3
    args.pseudotime_key = "pseudotime"
    args.skip_celloracle = False

    controller = PipelineController.__new__(PipelineController)
    controller.args = args
    controller.log_dir = str(tmp_path / "logs")
    os.makedirs(controller.log_dir, exist_ok=True)

    mock_oracle = MagicMock()
    mock_links = MagicMock()
    mock_links.merged_score = pd.DataFrame({
        "cluster": ["MEP", "GMP"],
        "gene": ["Gata1", "Spi1"],
        "degree_centrality": [0.8, 0.2],
        "betweenness_centrality": [0.7, 0.1],
        "eigenvector_centrality": [0.9, 0.3],
    })
    mock_links.filtered_links = {
        "MEP": pd.DataFrame({"source": ["Gata1"], "target": ["Klf1"], "coef_abs": [1.0]}),
        "GMP": pd.DataFrame({"source": ["Spi1"], "target": ["Cebpe"], "coef_abs": [1.0]}),
    }

    with patch("genecircuitry.perturbation.run_perturbation_analysis") as mock_run:
        mock_run.return_value = {"comparisons": {}, "summary_df": pd.DataFrame()}

        res = controller.run_step_perturbation(
            oracle=mock_oracle,
            links=mock_links,
            output_dir=out_dir,
        )

        assert mock_run.called
        kwargs = mock_run.call_args[1]
        assert kwargs["oracle"] is mock_oracle
        assert kwargs["links"] is mock_links
        assert kwargs["compare_clusters"] == [("MEP", "GMP")]
        assert kwargs["top_n"] == 2
        assert kwargs["perturbation_mode"] == "directional"
        assert kwargs["cluster_key"] == "celltype"
        assert kwargs["output_dir"] == out_dir
        assert res is not None
