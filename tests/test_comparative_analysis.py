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


