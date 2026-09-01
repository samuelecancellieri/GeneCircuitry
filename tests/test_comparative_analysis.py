"""
Unit Tests for Comparative Analysis Engine
==========================================

Tests data processing functions in genecircuitry.comparative_analysis.
"""

import os
import shutil
import tempfile
import numpy as np
import pandas as pd
from anndata import AnnData
import pytest

from genecircuitry.comparative_analysis import (
    compute_module_activity_matrix,
    compute_module_pathway_enrichments,
    compute_tf_centrality_matrix,
    compute_tf_to_module_mapping,
    compute_differential_tf_targets,
    compute_module_gene_overlap_matrix,
    compute_module_tf_integration,
    compute_gene_selection_provenance,
    compute_cross_cluster_regulatory_summary,
    run_comparative_analysis,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    d = tempfile.mkdtemp(prefix="test_comparative_analysis_")
    yield d
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def mock_adata():
    """Create synthetic AnnData with cluster and module scores."""
    n_cells = 50
    n_genes = 20
    X = np.random.randn(n_cells, n_genes)
    obs = pd.DataFrame({
        "leiden": ["0"] * 25 + ["1"] * 25,
        "Module_1": np.random.uniform(0.5, 2.0, n_cells),
        "Module_2": np.random.uniform(-1.0, 0.5, n_cells),
    })
    var = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
    return AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def mock_score_df():
    """Create synthetic CellOracle centrality scores DataFrame."""
    rows = []
    # Cluster 0
    rows.append({"gene": "TF_GLOBAL", "cluster": "0", "degree_centrality_all": 0.9, "betweenness_centrality": 0.8})
    rows.append({"gene": "TF_CLU0", "cluster": "0", "degree_centrality_all": 0.85, "betweenness_centrality": 0.7})
    rows.append({"gene": "TF_CLU1", "cluster": "0", "degree_centrality_all": 0.1, "betweenness_centrality": 0.05})
    # Cluster 1
    rows.append({"gene": "TF_GLOBAL", "cluster": "1", "degree_centrality_all": 0.88, "betweenness_centrality": 0.75})
    rows.append({"gene": "TF_CLU0", "cluster": "1", "degree_centrality_all": 0.05, "betweenness_centrality": 0.02})
    rows.append({"gene": "TF_CLU1", "cluster": "1", "degree_centrality_all": 0.92, "betweenness_centrality": 0.85})
    return pd.DataFrame(rows)


@pytest.fixture
def mock_links_df():
    """Create synthetic GRN regulatory links DataFrame."""
    rows = []
    # TF_GLOBAL regulates Module 1 genes in both clusters
    for g in ["Gene_0", "Gene_1", "Gene_2"]:
        rows.append({"source": "TF_GLOBAL", "target": g, "cluster": "0", "coef_abs": 0.8})
        rows.append({"source": "TF_GLOBAL", "target": g, "cluster": "1", "coef_abs": 0.75})
    # TF_GLOBAL regulates Gene_3 only in cluster 0
    rows.append({"source": "TF_GLOBAL", "target": "Gene_3", "cluster": "0", "coef_abs": 0.5})
    # TF_CLU1 regulates Module 2 genes in cluster 1
    for g in ["Gene_5", "Gene_6", "Gene_7"]:
        rows.append({"source": "TF_CLU1", "target": g, "cluster": "1", "coef_abs": 0.9})
    return pd.DataFrame(rows)


@pytest.fixture
def mock_modules_dict():
    """Create synthetic module-to-genes mapping."""
    return {
        "1": ["Gene_0", "Gene_1", "Gene_2", "Gene_3", "Gene_4"],
        "2": ["Gene_5", "Gene_6", "Gene_7", "Gene_8", "Gene_9"],
    }


class MockHotspotObj:
    """Mock Hotspot object with modules and results attributes."""
    def __init__(self):
        # modules: gene -> module_id mapping
        self.modules = pd.Series({
            "Gene_0": 1, "Gene_1": 1, "Gene_2": 1, "Gene_3": 1, "Gene_4": 1,
            "Gene_5": 2, "Gene_6": 2, "Gene_7": 2, "Gene_8": 2, "Gene_9": 2,
            "Gene_10": -1, "Gene_11": -1,
        })
        # results: autocorrelation results with FDR
        self.results = pd.DataFrame({
            "Z": np.random.randn(12),
            "Pval": np.random.uniform(0, 0.1, 12),
            "FDR": [0.001, 0.002, 0.003, 0.01, 0.02,
                    0.005, 0.006, 0.008, 0.015, 0.025,
                    0.04, 0.06],
        }, index=[f"Gene_{i}" for i in range(12)])
        self.module_scores = None


@pytest.fixture
def mock_hotspot_obj():
    """Create a mock Hotspot object with modules and results."""
    return MockHotspotObj()


def test_compute_module_activity_matrix(mock_adata):
    """Test module activity calculation across clusters."""
    activity_df = compute_module_activity_matrix(mock_adata, cluster_key="leiden")
    assert not activity_df.empty
    assert "Cluster 0" in activity_df.columns
    assert "Cluster 1" in activity_df.columns
    assert "Module 1" in activity_df.index
    assert "Module 2" in activity_df.index


def test_compute_module_pathway_enrichments(mock_modules_dict):
    """Test module pathway enrichment fallback / calculation."""
    enr_df = compute_module_pathway_enrichments(mock_modules_dict, top_n_terms=2)
    assert isinstance(enr_df, pd.DataFrame)
    assert not enr_df.empty
    assert "module" in enr_df.columns
    assert "term" in enr_df.columns


def test_compute_tf_centrality_matrix(mock_score_df):
    """Test TF centrality matrix and Global vs Group-Specific classification."""
    pivot_df, summary_df = compute_tf_centrality_matrix(mock_score_df, score="degree_centrality_all")
    assert not pivot_df.empty
    assert "0" in pivot_df.columns
    assert "1" in pivot_df.columns
    assert "TF_GLOBAL" in pivot_df.index
    assert not summary_df.empty
    assert "classification" in summary_df.columns

    # TF_GLOBAL should be classified as Global Master
    global_row = summary_df[summary_df["gene"] == "TF_GLOBAL"]
    assert len(global_row) == 1
    assert global_row.iloc[0]["classification"] == "Global Master"


def test_compute_tf_to_module_mapping(mock_links_df, mock_modules_dict):
    """Test mapping TFs to co-expression modules."""
    matrix_df, summary_df = compute_tf_to_module_mapping(mock_links_df, mock_modules_dict)
    assert not matrix_df.empty
    assert "Module 1" in matrix_df.columns
    assert "TF_GLOBAL" in matrix_df.index
    assert matrix_df.loc["TF_GLOBAL", "Module 1"] >= 3
    assert not summary_df.empty
    assert "tf" in summary_df.columns
    assert "module" in summary_df.columns


def test_compute_differential_tf_targets(mock_links_df):
    """Test differential TF target conservation calculation."""
    diff_df = compute_differential_tf_targets(mock_links_df)
    assert isinstance(diff_df, pd.DataFrame)
    assert not diff_df.empty
    assert "TF_GLOBAL" in diff_df["tf"].values
    tf_row = diff_df[diff_df["tf"] == "TF_GLOBAL"].iloc[0]
    assert tf_row["shared_targets_count"] >= 3
    assert tf_row["specific_targets_count"] >= 1


def test_run_comparative_analysis(temp_dir, mock_adata, mock_score_df, mock_links_df, mock_modules_dict):
    """Test end-to-end comparative analysis and table file generation."""
    # Add modules Series to adata
    results = run_comparative_analysis(
        adata=mock_adata,
        score_df=mock_score_df,
        links_df=mock_links_df,
        output_dir=temp_dir,
        cluster_key="leiden",
        save_tables=True,
    )

    assert isinstance(results, dict)
    assert "module_activity" in results
    assert "tf_centrality" in results
    assert "tf_to_module_matrix" in results
    assert "differential_tf_targets" in results

    # Verify CSV files exist on disk
    comp_dir = os.path.join(temp_dir, "comparative")
    assert os.path.exists(os.path.join(comp_dir, "module_activity_matrix.csv"))
    assert os.path.exists(os.path.join(comp_dir, "tf_centrality_matrix.csv"))
    assert os.path.exists(os.path.join(comp_dir, "differential_tf_targets.csv"))


def test_stratification_results_list_format(temp_dir, mock_adata):
    """Test compute_module_activity_matrix and run_comparative_analysis with list-of-dicts stratification_results."""
    adata1 = mock_adata[:25].copy()
    adata2 = mock_adata[25:].copy()

    # Pass as list of dicts (the format controller.stratification_results uses)
    strat_list = [
        {"name": "StratA", "adata": adata1, "output_dir": os.path.join(temp_dir, "StratA")},
        {"name": "StratB", "adata": adata2, "output_dir": os.path.join(temp_dir, "StratB")},
    ]

    act_df = compute_module_activity_matrix(stratification_results=strat_list)
    assert isinstance(act_df, pd.DataFrame)
    assert not act_df.empty
    assert "StratA" in act_df.columns
    assert "StratB" in act_df.columns

    # Also test run_comparative_analysis
    res = run_comparative_analysis(
        stratification_results=strat_list,
        output_dir=temp_dir,
        save_tables=True,
    )
    assert isinstance(res, dict)
    assert "module_activity" in res
    assert not res["module_activity"].empty


def test_compute_module_gene_overlap_matrix(mock_hotspot_obj, mock_links_df):
    """Test module gene overlap and Jaccard similarity computation."""
    coverage_df, jaccard_df = compute_module_gene_overlap_matrix(
        hotspot_obj=mock_hotspot_obj,
        links_df=mock_links_df,
    )
    assert isinstance(coverage_df, pd.DataFrame)
    assert isinstance(jaccard_df, pd.DataFrame)
    assert not coverage_df.empty
    assert "Module 1" in coverage_df.index or "Module 2" in coverage_df.index
    # Coverage values should be between 0 and 1
    assert coverage_df.min().min() >= 0.0
    assert coverage_df.max().max() <= 1.0


def test_compute_module_gene_overlap_matrix_empty():
    """Test module gene overlap with no data returns empty DataFrames."""
    coverage_df, jaccard_df = compute_module_gene_overlap_matrix()
    assert coverage_df.empty
    assert jaccard_df.empty


def test_compute_module_tf_integration(mock_hotspot_obj, mock_links_df, mock_score_df):
    """Test module-TF integration cross-reference."""
    integration_df = compute_module_tf_integration(
        links_df=mock_links_df,
        hotspot_obj=mock_hotspot_obj,
        score_df=mock_score_df,
    )
    assert isinstance(integration_df, pd.DataFrame)
    assert not integration_df.empty
    assert "cluster" in integration_df.columns
    assert "module" in integration_df.columns
    assert "tf" in integration_df.columns
    assert "n_targets_in_module" in integration_df.columns
    assert "coverage_pct" in integration_df.columns
    # All coverage percentages should be between 0 and 1
    assert integration_df["coverage_pct"].min() >= 0.0
    assert integration_df["coverage_pct"].max() <= 1.0


def test_compute_gene_selection_provenance(mock_hotspot_obj, mock_links_df):
    """Test gene selection provenance tracing."""
    provenance_df = compute_gene_selection_provenance(
        hotspot_obj=mock_hotspot_obj,
        links_df=mock_links_df,
    )
    assert isinstance(provenance_df, pd.DataFrame)
    assert not provenance_df.empty
    assert "gene" in provenance_df.columns
    assert "hotspot_fdr" in provenance_df.columns
    assert "hotspot_module" in provenance_df.columns
    assert "stage" in provenance_df.columns
    assert "is_tf" in provenance_df.columns
    assert "is_target" in provenance_df.columns
    # All significant genes should appear (FDR < 0.05)
    assert len(provenance_df) >= 10  # 11 genes have FDR < 0.05 in our mock
    # Stages should be valid
    valid_stages = {"TF & Target", "TF", "Target", "Module Member", "Significant Only"}
    assert set(provenance_df["stage"].unique()).issubset(valid_stages)


def test_compute_gene_selection_provenance_empty():
    """Test gene provenance with no hotspot returns empty DataFrame."""
    provenance_df = compute_gene_selection_provenance()
    assert provenance_df.empty


def test_compute_cross_cluster_regulatory_summary(mock_hotspot_obj, mock_links_df, mock_score_df):
    """Test cross-cluster regulatory summary computation."""
    summary_df = compute_cross_cluster_regulatory_summary(
        links_df=mock_links_df,
        score_df=mock_score_df,
        hotspot_obj=mock_hotspot_obj,
    )
    assert isinstance(summary_df, pd.DataFrame)
    assert not summary_df.empty
    assert "cluster" in summary_df.columns
    assert "n_regulatory_edges" in summary_df.columns
    assert "n_unique_targets" in summary_df.columns
    assert "top_tfs" in summary_df.columns
    assert "n_active_tfs" in summary_df.columns
    # Should have one row per cluster
    assert len(summary_df) == 2  # clusters "0" and "1"


def test_run_comparative_analysis_with_new_functions(
    temp_dir, mock_adata, mock_score_df, mock_links_df, mock_hotspot_obj
):
    """Test end-to-end comparative analysis includes new result keys."""
    results = run_comparative_analysis(
        adata=mock_adata,
        score_df=mock_score_df,
        links_df=mock_links_df,
        hotspot_obj=mock_hotspot_obj,
        output_dir=temp_dir,
        cluster_key="leiden",
        save_tables=True,
    )
    assert isinstance(results, dict)
    # Check new result keys exist
    assert "module_coverage" in results
    assert "module_jaccard" in results
    assert "module_tf_integration" in results
    assert "gene_provenance" in results
    assert "regulatory_summary" in results
    # Check new CSV files were written
    comp_dir = os.path.join(temp_dir, "comparative")
    assert os.path.exists(os.path.join(comp_dir, "module_gene_coverage.csv"))
    assert os.path.exists(os.path.join(comp_dir, "module_tf_integration.csv"))
    assert os.path.exists(os.path.join(comp_dir, "gene_selection_provenance.csv"))
    assert os.path.exists(os.path.join(comp_dir, "cross_cluster_regulatory_summary.csv"))
