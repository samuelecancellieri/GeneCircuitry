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
from genecircuitry.preprocessing import (
    parse_cluster_keys,
    sanitize_identifier,
    resolve_cluster_key_name,
    resolve_cluster_key,
    ensure_categorical_obs,
)
from genecircuitry.pipeline.controller import (
    PipelineController,
    celloracle_pipeline,
    dimensionality_reduction_clustering,
    generate_summary,
    hotspot_pipeline,
    stratification_pipeline,
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

    def test_perform_grn_pre_processing_cell_downsample_explicit(self, dummy_adata):
        """Test perform_grn_pre_processing with explicit cell_downsample."""
        result = perform_grn_pre_processing(
            dummy_adata,
            cluster_key="cell_type",
            cell_downsample=20,
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result is not None
        assert result.n_obs == 20

    def test_perform_grn_pre_processing_cell_downsample_from_config(self, dummy_adata):
        """Test perform_grn_pre_processing falls back to config.GRN_CELL_DOWNSAMPLE."""
        from genecircuitry import config

        orig = config.GRN_CELL_DOWNSAMPLE
        try:
            config.update_config(GRN_CELL_DOWNSAMPLE=25)
            result = perform_grn_pre_processing(
                dummy_adata,
                cluster_key="cell_type",
                cell_downsample=None,
                top_genes=10,
                n_neighbors=5,
                n_pcs=5,
            )
            assert result is not None
            assert result.n_obs == 25
        finally:
            config.update_config(GRN_CELL_DOWNSAMPLE=orig)

    def test_perform_grn_pre_processing_cell_downsample_disabled(self, dummy_adata):
        """Test perform_grn_pre_processing when downsampling is disabled (0 or None in config)."""
        result = perform_grn_pre_processing(
            dummy_adata,
            cluster_key="cell_type",
            cell_downsample=0,
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result is not None
        assert result.n_obs == dummy_adata.n_obs

    def test_perform_grn_pre_processing_cell_downsample_larger_than_data(self, dummy_adata):
        """Test perform_grn_pre_processing when downsample target exceeds cell count."""
        result = perform_grn_pre_processing(
            dummy_adata,
            cluster_key="cell_type",
            cell_downsample=1000,
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result is not None
        assert result.n_obs == dummy_adata.n_obs

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

    def test_cli_cell_downsample_flags(self):
        """Test CLI argument parsing for --cell-downsample and --grn-cell-downsample."""
        from genecircuitry import config

        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--cell-downsample",
            "--grn-cell-downsample",
            type=int,
            default=config.GRN_CELL_DOWNSAMPLE,
            dest="cell_downsample",
        )

        args_default = parser.parse_args([])
        assert args_default.cell_downsample == 30000
        args1 = parser.parse_args(["--cell-downsample", "50000"])
        assert args1.cell_downsample == 50000

        args2 = parser.parse_args(["--grn-cell-downsample", "100000"])
        assert args2.cell_downsample == 100000



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


class TestMultiKeyParsingAndHelpers:
    """Test parse_cluster_keys, sanitize_identifier, resolve_cluster_key_name, resolve_cluster_key, ensure_categorical_obs."""

    def test_parse_cluster_keys(self):
        assert parse_cluster_keys(None) == []
        assert parse_cluster_keys("cell_type") == ["cell_type"]
        assert parse_cluster_keys("cell_type,condition") == ["cell_type", "condition"]
        assert parse_cluster_keys(" cell_type , condition , batch ") == [
            "cell_type",
            "condition",
            "batch",
        ]
        assert parse_cluster_keys("k1, k2, k1") == ["k1", "k2"]
        assert parse_cluster_keys(["cell_type", "condition"]) == ["cell_type", "condition"]
        assert parse_cluster_keys(("cell_type", "condition")) == ["cell_type", "condition"]
        assert parse_cluster_keys({"b_key", "a_key"}) == ["a_key", "b_key"]
        assert parse_cluster_keys(np.array(["cell_type", "condition"])) == [
            "cell_type",
            "condition",
        ]
        assert parse_cluster_keys(pd.Index(["cell_type", "condition"])) == [
            "cell_type",
            "condition",
        ]
        # Dict keys view and generator expressions
        assert parse_cluster_keys({"cell_type": 1, "condition": 2}.keys()) == [
            "cell_type",
            "condition",
        ]
        assert parse_cluster_keys(k for k in ["cell_type", "condition"]) == [
            "cell_type",
            "condition",
        ]
        assert parse_cluster_keys(k for k in ["cell_type, condition", "batch"]) == [
            "cell_type",
            "condition",
            "batch",
        ]
        # Mixed type sets and frozensets
        assert parse_cluster_keys({"key1", 2}) == ["2", "key1"]
        assert parse_cluster_keys(frozenset(["b", "a"])) == ["a", "b"]
        assert parse_cluster_keys(frozenset([1, "a"])) == ["1", "a"]
        assert parse_cluster_keys(123) == ["123"]

    def test_sanitize_identifier(self):
        assert sanitize_identifier("Type A") == "Type_A"
        assert sanitize_identifier("Type B/Subtype") == "Type_B-Subtype"
        assert sanitize_identifier("val?*<>|:") == "val______"
        assert sanitize_identifier(123) == "123"

    def test_resolve_cluster_key_name(self):
        assert resolve_cluster_key_name(None) == "leiden"
        assert resolve_cluster_key_name("cell_type") == "cell_type"
        assert resolve_cluster_key_name("cell_type,condition") == "cell_type_condition"
        assert (
            resolve_cluster_key_name(["cell_type", "condition"])
            == "cell_type_condition"
        )
        assert resolve_cluster_key_name({"b_key", "a_key"}) == "a_key_b_key"

    def test_resolve_cluster_key(self, dummy_adata):
        # Single key existing
        adata_out, col = resolve_cluster_key(dummy_adata, "cell_type")
        assert col == "cell_type"
        assert isinstance(adata_out.obs["cell_type"].dtype, pd.CategoricalDtype)

        # Single key missing
        with pytest.raises(ValueError, match="Cluster key 'missing' not found"):
            resolve_cluster_key(dummy_adata, "missing")

        # Multi-key string
        adata_out, col = resolve_cluster_key(dummy_adata, "cell_type,custom_col")
        assert col == "cell_type_custom_col"
        assert col in adata_out.obs.columns
        assert isinstance(adata_out.obs[col].dtype, pd.CategoricalDtype)
        # Check values
        assert "B_cell_Type_A" in adata_out.obs[col].values

        # Multi-key sequence
        adata_out, col = resolve_cluster_key(dummy_adata, ["cell_type", "custom_col"])
        assert col == "cell_type_custom_col"

        # Multi-key missing one column
        with pytest.raises(ValueError, match="Cluster key 'missing' not found"):
            resolve_cluster_key(dummy_adata, "cell_type,missing")

        # Multi-key missing multiple columns
        with pytest.raises(ValueError, match="Cluster keys.*not found"):
            resolve_cluster_key(dummy_adata, "missing1,missing2")

        # View adata safety
        adata_view = dummy_adata[0:10]
        adata_out, col = resolve_cluster_key(adata_view, "cell_type,custom_col")
        assert not adata_out.is_view
        assert col in adata_out.obs.columns

    def test_ensure_categorical_obs_string_and_sequences(self, dummy_adata):
        dummy_adata.obs["str_col1"] = ["a", "b"] * (dummy_adata.n_obs // 2)
        dummy_adata.obs["str_col2"] = ["x", "y"] * (dummy_adata.n_obs // 2)
        assert dummy_adata.obs["str_col1"].dtype == "object"
        assert dummy_adata.obs["str_col2"].dtype == "object"

        # String single column
        adata = ensure_categorical_obs(dummy_adata, columns="str_col1")
        assert isinstance(adata.obs["str_col1"].dtype, pd.CategoricalDtype)

        # String comma-separated columns
        adata = ensure_categorical_obs(dummy_adata, columns="str_col1,str_col2")
        assert isinstance(adata.obs["str_col2"].dtype, pd.CategoricalDtype)

        # Dict keys view and generator inputs
        dummy_adata.obs["str_col3"] = ["m", "n"] * (dummy_adata.n_obs // 2)
        dummy_adata.obs["str_col4"] = ["p", "q"] * (dummy_adata.n_obs // 2)
        adata = ensure_categorical_obs(
            dummy_adata, columns={"str_col3": 1, "str_col4": 2}.keys()
        )
        assert isinstance(adata.obs["str_col3"].dtype, pd.CategoricalDtype)
        assert isinstance(adata.obs["str_col4"].dtype, pd.CategoricalDtype)

        dummy_adata.obs["str_col5"] = ["u", "v"] * (dummy_adata.n_obs // 2)
        adata = ensure_categorical_obs(
            dummy_adata, columns=(k for k in ["str_col5"])
        )
        assert isinstance(adata.obs["str_col5"].dtype, pd.CategoricalDtype)

        # AnnData view input
        view = dummy_adata[0:10]
        assert view.is_view
        adata_from_view = ensure_categorical_obs(view, columns="str_col1")
        assert not adata_from_view.is_view
        assert isinstance(adata_from_view.obs["str_col1"].dtype, pd.CategoricalDtype)


class TestMultiKeyStratification:
    """Test multi-key stratification in stratification_pipeline."""

    def test_stratification_pipeline_multi_key_string(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,condition",
        )
        assert len(adata_list) > 0
        assert len(adata_list) == len(names)
        assert "B_cell_Ctrl" in names

    def test_stratification_pipeline_multi_key_list(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification=["cell_type", "condition"],
        )
        assert len(adata_list) > 0
        assert len(adata_list) == len(names)

    def test_stratification_pipeline_multi_key_set(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification={"cell_type", "condition"},
        )
        assert len(adata_list) > 0
        assert len(adata_list) == len(names)

    def test_stratification_pipeline_multi_key_generator_and_dict_keys(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification=(k for k in ["cell_type", "condition"]),
            clusters=(c for c in ["B_cell_Ctrl"]),
        )
        assert len(adata_list) == 1
        assert names == ["B_cell_Ctrl"]

        adata_list2, names2 = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification={"cell_type": 1, "condition": 2}.keys(),
            clusters={"B_cell_Ctrl": 1, "T_cell_Treat": 2}.keys(),
        )
        assert len(adata_list2) == 2
        assert set(names2) == {"B_cell_Ctrl", "T_cell_Treat"}

    def test_stratification_pipeline_three_keys(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        dummy_adata.obs["batch"] = ["b1"] * dummy_adata.n_obs
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification=["cell_type", "condition", "batch"],
        )
        assert len(adata_list) > 0
        for name in names:
            assert name.endswith("_b1")

    def test_stratification_pipeline_special_chars_sanitization(self, dummy_adata):
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,custom_col",
        )
        # custom_col has 'Type A', 'Type B/Subtype', 'Type C'
        # Check that sanitized names contain Type_A, Type_B-Subtype
        assert any("Type_A" in n for n in names)
        assert any("Type_B-Subtype" in n for n in names)
        for n in names:
            assert " " not in n
            assert "/" not in n

    def test_stratification_pipeline_clusters_filter_multi_key(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,condition",
            clusters="B_cell_Ctrl,T_cell_Treat",
        )
        assert len(adata_list) == 2
        assert set(names) == {"B_cell_Ctrl", "T_cell_Treat"}

    def test_stratification_pipeline_clusters_filter_list_and_set(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,condition",
            clusters=["B_cell_Ctrl", "T_cell_Treat"],
        )
        assert len(adata_list) == 2
        assert set(names) == {"B_cell_Ctrl", "T_cell_Treat"}

        adata_list_set, names_set = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,condition",
            clusters={"B_cell_Ctrl"},
        )
        assert len(adata_list_set) == 1
        assert names_set == ["B_cell_Ctrl"]

    def test_stratification_pipeline_filter_by_individual_key_cluster(self, dummy_adata):
        """Test filtering multi-key combinations using an individual cluster from one key."""
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        
        # Select only 'B_cell' -> should get both 'B_cell_Ctrl' and 'B_cell_Treat'
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,condition",
            clusters="B_cell",
        )
        assert len(adata_list) == 2
        assert set(names) == {"B_cell_Ctrl", "B_cell_Treat"}
        for ad_sub in adata_list:
            assert set(ad_sub.obs["cell_type"].unique()) == {"B_cell"}

        # Select only 'Treat' from condition -> should get B_cell_Treat, T_cell_Treat, Monocyte_Treat
        adata_list2, names2 = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification=["cell_type", "condition"],
            clusters=["Treat"],
        )
        assert len(adata_list2) == 3
        assert set(names2) == {"B_cell_Treat", "T_cell_Treat", "Monocyte_Treat"}

    def test_stratification_pipeline_filter_by_key_value_syntax(self, dummy_adata):
        """Test filtering multi-key combinations using key:val or key=val syntax."""
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)

        adata_list1, names1 = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,condition",
            clusters="cell_type:B_cell",
        )
        assert len(adata_list1) == 2
        assert set(names1) == {"B_cell_Ctrl", "B_cell_Treat"}

        adata_list2, names2 = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type,condition",
            clusters="condition=Ctrl",
        )
        assert len(adata_list2) == 3
        assert set(names2) == {"B_cell_Ctrl", "T_cell_Ctrl", "Monocyte_Ctrl"}

    def test_stratification_pipeline_unobserved_combinations_skipped(self, dummy_adata):
        # Create non-overlapping groups
        dummy_adata.obs["group1"] = ["A"] * 30 + ["B"] * 30
        dummy_adata.obs["group2"] = ["X"] * 30 + ["Y"] * 30
        # A only occurs with X, B only occurs with Y (A_Y and B_X are unobserved)
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="group1,group2",
        )
        assert len(adata_list) == 2
        assert set(names) == {"A_X", "B_Y"}

    def test_stratification_pipeline_single_key_backward_compatible(self, dummy_adata):
        adata_list, names = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="cell_type",
            clusters="B_cell",
        )
        assert len(adata_list) == 1
        assert names == ["B_cell"]

    def test_stratification_pipeline_view_safety(self, dummy_adata):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        view = dummy_adata[0:20]
        assert view.is_view
        adata_list, names = stratification_pipeline(
            view,
            cluster_key_stratification="cell_type,condition",
        )
        assert len(adata_list) > 0
        for ad_sub in adata_list:
            assert not ad_sub.is_view

    def test_stratification_pipeline_numeric_clusters_filter(self, dummy_adata):
        dummy_adata.obs["num_cat"] = [0, 1] * (dummy_adata.n_obs // 2)
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)

        # Single key numeric filter
        adata_list1, names1 = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="num_cat",
            clusters="0",
        )
        assert len(adata_list1) == 1
        assert names1 == ["0"]

        adata_list2, names2 = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification="num_cat",
            clusters=[0],
        )
        assert len(adata_list2) == 1
        assert names2 == ["0"]

        # Multi-key numeric + string filter
        adata_list3, names3 = stratification_pipeline(
            dummy_adata,
            cluster_key_stratification=["num_cat", "condition"],
            clusters="0_Ctrl",
        )
        assert len(adata_list3) == 1
        assert names3 == ["0_Ctrl"]


class TestMultiKeyCellOracle:
    """Test multi-key clustering in CellOracle processing."""

    def test_perform_grn_pre_processing_multi_key_string(self, dummy_adata):
        result = perform_grn_pre_processing(
            dummy_adata,
            cluster_key="cell_type,custom_col",
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result is not None
        assert "cell_type_custom_col" in result.obs.columns
        assert isinstance(
            result.obs["cell_type_custom_col"].dtype, pd.CategoricalDtype
        )
        assert "paga" in result.uns

    def test_perform_grn_pre_processing_multi_key_list(self, dummy_adata):
        result = perform_grn_pre_processing(
            dummy_adata,
            cluster_key=["cell_type", "custom_col"],
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result is not None
        assert "cell_type_custom_col" in result.obs.columns

    def test_perform_grn_pre_processing_multi_key_generator_and_dict_keys(self, dummy_adata):
        result_gen = perform_grn_pre_processing(
            dummy_adata,
            cluster_key=(k for k in ["cell_type", "custom_col"]),
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result_gen is not None
        assert "cell_type_custom_col" in result_gen.obs.columns

        result_dict = perform_grn_pre_processing(
            dummy_adata,
            cluster_key={"cell_type": 1, "custom_col": 2}.keys(),
            top_genes=10,
            n_neighbors=5,
            n_pcs=5,
        )
        assert result_dict is not None
        assert "cell_type_custom_col" in result_dict.obs.columns

    @patch("celloracle.Oracle")
    @patch("celloracle.data.load_human_promoter_base_GRN", return_value=pd.DataFrame())
    def test_create_oracle_object_multi_key(
        self, mock_base_grn, mock_oracle_cls, dummy_adata
    ):
        mock_oracle_instance = MagicMock()
        mock_oracle_cls.return_value = mock_oracle_instance

        oracle = create_oracle_object(
            dummy_adata,
            cluster_column_name="cell_type,custom_col",
            embedding_name="X_umap",
            species="human",
        )
        assert oracle is not None
        mock_oracle_instance.import_anndata_as_normalized_count.assert_called_once()
        call_kwargs = mock_oracle_instance.import_anndata_as_normalized_count.call_args[1]
        assert call_kwargs["cluster_column_name"] == "cell_type_custom_col"
        cleaned_adata = call_kwargs["adata"]
        assert "cell_type_custom_col" in cleaned_adata.obs.columns

    def test_run_links_multi_key(self):
        mock_oracle = MagicMock()
        mock_links = MagicMock()
        mock_links.cluster = ["A_X", "B_Y"]
        mock_oracle.get_links.return_value = mock_links

        with patch("genecircuitry.config.GRN_N_JOBS", 1):
            links = run_links(mock_oracle, cluster_column_name="group1,group2")

        assert links is mock_links
        mock_oracle.get_links.assert_called_once_with(
            cluster_name_for_GRN_unit="group1_group2",
            alpha=10,
            verbose_level=10,
            n_jobs=1,
        )


class TestMultiKeyPipelineControllerAndReporting:
    """Test multi-key clustering and stratification in controller and reporting."""

    def test_controller_multi_key_clustering(self, dummy_adata, tmp_path):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        output_dir = str(tmp_path / "controller_multikey")
        os.makedirs(output_dir, exist_ok=True)

        args = argparse.Namespace(
            output=output_dir,
            name="test_multikey",
            cluster_key="cell_type,condition",
            cluster_key_stratification=None,
            clusters="all",
            force_dim_reduction=True,
        )
        controller = PipelineController(args, datetime.now())
        controller.adata_preprocessed = dummy_adata

        adata_clustered = controller.run_step_clustering()
        assert "cell_type_condition" in adata_clustered.obs.columns
        assert isinstance(
            adata_clustered.obs["cell_type_condition"].dtype, pd.CategoricalDtype
        )

    def test_controller_multi_key_stratification(self, dummy_adata, tmp_path):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        output_dir = str(tmp_path / "controller_strat_multikey")
        os.makedirs(output_dir, exist_ok=True)

        args = argparse.Namespace(
            output=output_dir,
            name="test_strat_multikey",
            cluster_key="leiden",
            cluster_key_stratification="cell_type,condition",
            clusters="B_cell_Ctrl",
            force_dim_reduction=False,
            skip_celloracle=True,
            skip_hotspot=True,
        )
        controller = PipelineController(args, datetime.now())
        controller.adata_preprocessed = dummy_adata

        adata_list, names = controller.run_step_stratification()
        assert len(adata_list) == 1
        assert names == ["B_cell_Ctrl"]

    def test_controller_multi_key_stratification_single_cluster_filter(self, dummy_adata, tmp_path):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        output_dir = str(tmp_path / "controller_strat_single_cluster_filter")
        os.makedirs(output_dir, exist_ok=True)

        args = argparse.Namespace(
            output=output_dir,
            name="test_strat_filter",
            cluster_key="leiden",
            cluster_key_stratification="cell_type,condition",
            clusters="B_cell",
            force_dim_reduction=False,
            skip_celloracle=True,
            skip_hotspot=True,
        )
        controller = PipelineController(args, datetime.now())
        controller.adata_preprocessed = dummy_adata

        adata_list, names = controller.run_step_stratification()
        assert len(adata_list) == 2
        assert set(names) == {"B_cell_Ctrl", "B_cell_Treat"}

    def test_reporting_multi_key_clustering_section(self, dummy_adata, tmp_path):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        from genecircuitry.preprocessing import resolve_cluster_key

        adata_out, col = resolve_cluster_key(dummy_adata, "cell_type,condition")

        section = create_clustering_section(
            adata_out,
            output_dir=str(tmp_path),
            cluster_key="cell_type,condition",
        )
        assert section is not None
        assert "Clusters (cell_type_condition)" in section.metrics

    def test_generate_summary_multi_key(self, dummy_adata, tmp_path):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        from genecircuitry.preprocessing import resolve_cluster_key

        adata_out, col = resolve_cluster_key(dummy_adata, "cell_type,condition")

        summary_dir = tmp_path / "summary_multikey"
        summary_dir.mkdir()
        generate_summary(
            adata_out,
            celloracle_result=None,
            hotspot_result=None,
            start_time=datetime.now(),
            output_dir=str(summary_dir),
            cluster_key="cell_type,condition",
        )
        summary_file = summary_dir / "analysis_summary.txt"
        assert summary_file.exists()
        content = summary_file.read_text()
        assert "Clusters identified (cell_type_condition):" in content

    def test_cli_multikey_parsing(self):
        """Test CLI parser accepts multi-key strings for cluster-key and stratification."""
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--cluster-key",
            "--cluster-column",
            "--cluster-column-name",
            type=str,
            default="leiden",
            dest="cluster_key",
        )
        parser.add_argument(
            "--cluster-key-stratification",
            type=str,
            default=None,
        )
        parser.add_argument(
            "--clusters",
            type=str,
            default="all",
        )

        args = parser.parse_args([
            "--cluster-key", "cell_type,condition",
            "--cluster-key-stratification", "group1,group2",
            "--clusters", "val1_val2,val3_val4",
        ])
        assert args.cluster_key == "cell_type,condition"
        assert args.cluster_key_stratification == "group1,group2"
        assert args.clusters == "val1_val2,val3_val4"

        # Verify parsing through helper
        assert parse_cluster_keys(args.cluster_key) == ["cell_type", "condition"]
        assert parse_cluster_keys(args.cluster_key_stratification) == ["group1", "group2"]

    def test_generate_stratified_report_multi_key(self, dummy_adata, tmp_path):
        dummy_adata.obs["condition"] = ["Ctrl", "Treat"] * (dummy_adata.n_obs // 2)
        from genecircuitry.preprocessing import resolve_cluster_key

        adata_out, col = resolve_cluster_key(dummy_adata, "cell_type,condition")
        report_dir = tmp_path / "strat_report_multikey"
        report_dir.mkdir()

        strat_results = [
            {
                "name": "B_cell_Ctrl",
                "adata": adata_out,
                "output_dir": str(report_dir / "B_cell_Ctrl"),
            }
        ]
        outputs = generate_stratified_report(
            output_dir=str(report_dir),
            title="Test MultiKey Strat Report",
            adata_preprocessed=dummy_adata,
            stratification_results=strat_results,
            formats=["html"],
            embed_images=False,
            cluster_key="cell_type,condition",
        )
        assert "html" in outputs
        assert os.path.exists(outputs["html"])


