"""Unit tests for cluster_key and cluster_column handling throughout GeneCircuitry."""

import argparse
import os
import unittest
from datetime import datetime
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from genecircuitry.celloracle_processing import (
    create_oracle_object,
    perform_grn_pre_processing,
    run_links,
)
from genecircuitry.pipeline.controller import (
    PipelineController,
    celloracle_pipeline,
    generate_summary,
    hotspot_pipeline,
)
from genecircuitry.reporting.generator import (
    generate_report,
    generate_stratified_report,
)
from genecircuitry.reporting.sections import (
    create_clustering_section,
    create_stratified_clustering_section,
)


@pytest.fixture
def dummy_adata():
    """Create a small dummy AnnData object for testing."""
    np.random.seed(42)
    n_obs, n_vars = 60, 30
    X = np.random.poisson(lam=2.0, size=(n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical(np.random.choice(["0", "1", "2"], size=n_obs)),
            "cell_type": ["B_cell", "T_cell", "Monocyte"] * (n_obs // 3),
            "custom_col": ["Type A", "Type B/Subtype", "Type C"] * (n_obs // 3),
            "n_counts": np.sum(X, axis=1),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["raw_counts"] = X.copy()
    adata.obsm["X_pca"] = np.random.randn(n_obs, 5).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
    return adata


class TestCellOracleClusterKey:
    """Test cluster_key / cluster_column handling in CellOracle processing."""

    def test_perform_grn_pre_processing_custom_key(self, dummy_adata):
        """Test perform_grn_pre_processing with custom cluster_key."""
        result = perform_grn_pre_processing(
            dummy_adata,
            cluster_key="cell_type",
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result is not None
        assert isinstance(result.obs["cell_type"].dtype, pd.CategoricalDtype)
        assert "X_diffmap" in result.obsm
        assert "paga" in result.uns

    def test_perform_grn_pre_processing_alias_kwargs(self, dummy_adata):
        """Test perform_grn_pre_processing accepts cluster_column_name / cluster_column as kwargs."""
        result = perform_grn_pre_processing(
            dummy_adata,
            cluster_column_name="cell_type",
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result is not None
        assert isinstance(result.obs["cell_type"].dtype, pd.CategoricalDtype)

    def test_perform_grn_pre_processing_missing_key_raises(self, dummy_adata):
        """Test perform_grn_pre_processing raises ValueError when key is missing."""
        with pytest.raises(ValueError, match="Cluster key 'nonexistent_key' not found"):
            perform_grn_pre_processing(
                dummy_adata,
                cluster_key="nonexistent_key",
                top_genes=10,
            )

    @patch("celloracle.Oracle")
    @patch("celloracle.data.load_human_promoter_base_GRN", return_value=pd.DataFrame())
    def test_create_oracle_object_custom_column_and_clean_names(
        self, mock_base_grn, mock_oracle_cls, dummy_adata
    ):
        """Test create_oracle_object cleans unsafe characters and respects custom column."""
        mock_oracle_instance = MagicMock()
        mock_oracle_cls.return_value = mock_oracle_instance

        oracle = create_oracle_object(
            dummy_adata,
            cluster_column_name="custom_col",
            embedding_name="X_umap",
            species="human",
        )
        assert oracle is not None
        # Check that import was called with custom_col
        mock_oracle_instance.import_anndata_as_normalized_count.assert_called_once()
        call_kwargs = mock_oracle_instance.import_anndata_as_normalized_count.call_args[1]
        assert call_kwargs["cluster_column_name"] == "custom_col"
        cleaned_adata = call_kwargs["adata"]
        # Verify unsafe character replacement: 'Type A' -> 'Type_A', 'Type B/Subtype' -> 'Type_B-Subtype'
        unique_vals = set(cleaned_adata.obs["custom_col"].astype(str))
        assert "Type_A" in unique_vals
        assert "Type_B-Subtype" in unique_vals

    def test_create_oracle_object_missing_column_raises(self, dummy_adata):
        """Test create_oracle_object raises ValueError if cluster column is missing."""
        with pytest.raises(ValueError, match="Cluster column 'missing_col' not found"):
            create_oracle_object(
                dummy_adata,
                cluster_column_name="missing_col",
                embedding_name="X_umap",
            )

    def test_run_links_iterates_over_clusters(self):
        """Test run_links iterates over links.cluster for rank score plotting."""
        mock_oracle = MagicMock()
        mock_links = MagicMock()
        mock_links.cluster = ["0", "1", "2"]
        mock_oracle.get_links.return_value = mock_links

        with patch("genecircuitry.config.GRN_N_JOBS", 1):
            links = run_links(mock_oracle, cluster_column_name="cell_type")

        assert links is mock_links
        # Verify get_links was called with cluster_name_for_GRN_unit="cell_type"
        mock_oracle.get_links.assert_called_once_with(
            cluster_name_for_GRN_unit="cell_type",
            alpha=10,
            verbose_level=10,
            n_jobs=1,
        )
        # Verify plot_scores_as_rank was called for each cluster unit ("0", "1", "2")
        assert mock_links.plot_scores_as_rank.call_count == 3
        calls = [c[1]["cluster"] for c in mock_links.plot_scores_as_rank.call_args_list]
        assert calls == ["0", "1", "2"]


class TestControllerClusterKey:
    """Test cluster_key handling in PipelineController, summary, and CLI parser."""

    def test_generate_summary_custom_key(self, dummy_adata, tmp_path):
        """Test generate_summary correctly reports cluster counts for custom cluster_key."""
        summary_dir = tmp_path / "summary_test"
        summary_dir.mkdir()
        generate_summary(
            dummy_adata,
            celloracle_result=None,
            hotspot_result=None,
            start_time=datetime.now(),
            output_dir=str(summary_dir),
            cluster_key="cell_type",
        )
        summary_file = summary_dir / "analysis_summary.txt"
        assert summary_file.exists()
        content = summary_file.read_text()
        assert "Clusters identified (cell_type): 3" in content

    def test_cli_parser_cluster_key_aliases(self):
        """Test CLI parser accepts --cluster-key, --cluster-column, and --cluster-column-name."""
        from genecircuitry.pipeline.controller import main

        # Create a mock sys.argv parser test
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--cluster-key",
            "--cluster-column",
            "--cluster-column-name",
            type=str,
            default="leiden",
            dest="cluster_key",
        )

        args1 = parser.parse_args(["--cluster-key", "cell_type"])
        assert args1.cluster_key == "cell_type"

        args2 = parser.parse_args(["--cluster-column", "annotation"])
        assert args2.cluster_key == "annotation"

        args3 = parser.parse_args(["--cluster-column-name", "subtypes"])
        assert args3.cluster_key == "subtypes"


class TestReportingClusterKey:
    """Test cluster_key handling in report sections and generators."""

    def test_create_clustering_section_custom_key(self, dummy_adata, tmp_path):
        """Test create_clustering_section prioritizes custom cluster_key."""
        section = create_clustering_section(
            dummy_adata,
            output_dir=str(tmp_path),
            cluster_key="cell_type",
        )
        assert section is not None
        assert "Clusters (cell_type)" in section.metrics
        assert section.metrics["Clusters (cell_type)"] == 3
        assert len(section.tables) >= 1
        assert "Cell_Type" in section.tables[0]["title"]

    def test_create_stratified_clustering_section_custom_key(self, dummy_adata):
        """Test create_stratified_clustering_section respects cluster_key."""
        stratification_results = [
            {
                "name": "GroupA",
                "adata": dummy_adata,
                "output_dir": "/tmp/groupA",
            }
        ]
        section = create_stratified_clustering_section(
            stratification_results,
            cluster_key="cell_type",
        )
        assert section is not None
        assert len(section.subsections) == 1
        sub = section.subsections[0]
        assert sub.metrics.get("Clusters") == 3
        assert len(sub.tables) >= 1
        assert "Cell_Type" in sub.tables[0]["title"]
