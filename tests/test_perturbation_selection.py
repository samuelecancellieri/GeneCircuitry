"""
Unit Tests for TF Selection in Perturbation Engine
===================================================

Tests genecircuitry.perturbation.tf_selection functions for selecting top candidate TFs
based on differential centrality and Hotspot module overlap.
"""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from genecircuitry.perturbation.tf_selection import (
    select_top_perturbed_tfs,
    compute_cluster_pair_centrality_diff,
    compute_tf_hotspot_enrichment,
)


@pytest.fixture
def sample_centrality_score_df():
    """Create a sample score DataFrame with centralities across two clusters."""
    data = [
        {"gene": "GATA1", "cluster": "Erythroid", "degree_centrality_all": 0.8, "betweenness_centrality": 0.6, "eigenvector_centrality": 0.75},
        {"gene": "GATA1", "cluster": "Myeloid", "degree_centrality_all": 0.1, "betweenness_centrality": 0.05, "eigenvector_centrality": 0.1},
        {"gene": "SPI1", "cluster": "Erythroid", "degree_centrality_all": 0.15, "betweenness_centrality": 0.1, "eigenvector_centrality": 0.1},
        {"gene": "SPI1", "cluster": "Myeloid", "degree_centrality_all": 0.85, "betweenness_centrality": 0.7, "eigenvector_centrality": 0.8},
        {"gene": "KLF1", "cluster": "Erythroid", "degree_centrality_all": 0.7, "betweenness_centrality": 0.5, "eigenvector_centrality": 0.65},
        {"gene": "KLF1", "cluster": "Myeloid", "degree_centrality_all": 0.2, "betweenness_centrality": 0.1, "eigenvector_centrality": 0.15},
        {"gene": "CEBPA", "cluster": "Erythroid", "degree_centrality_all": 0.1, "betweenness_centrality": 0.05, "eigenvector_centrality": 0.05},
        {"gene": "CEBPA", "cluster": "Myeloid", "degree_centrality_all": 0.65, "betweenness_centrality": 0.45, "eigenvector_centrality": 0.6},
        {"gene": "STAT3", "cluster": "Erythroid", "degree_centrality_all": 0.4, "betweenness_centrality": 0.3, "eigenvector_centrality": 0.4},
        {"gene": "STAT3", "cluster": "Myeloid", "degree_centrality_all": 0.4, "betweenness_centrality": 0.3, "eigenvector_centrality": 0.4},
    ]
    return pd.DataFrame(data)


@pytest.fixture
def sample_links_df():
    """Create sample links DataFrame connecting TFs to targets."""
    data = [
        # GATA1 targets in Erythroid
        {"source": "GATA1", "target": "HBB", "cluster": "Erythroid", "coef_abs": 0.8},
        {"source": "GATA1", "target": "HBA1", "cluster": "Erythroid", "coef_abs": 0.75},
        {"source": "GATA1", "target": "ALAS2", "cluster": "Erythroid", "coef_abs": 0.65},
        {"source": "GATA1", "target": "KLF1", "cluster": "Erythroid", "coef_abs": 0.6},
        # KLF1 targets in Erythroid
        {"source": "KLF1", "target": "HBB", "cluster": "Erythroid", "coef_abs": 0.7},
        {"source": "KLF1", "target": "HBA1", "cluster": "Erythroid", "coef_abs": 0.65},
        # SPI1 targets in Myeloid
        {"source": "SPI1", "target": "MPO", "cluster": "Myeloid", "coef_abs": 0.85},
        {"source": "SPI1", "target": "LYZ", "cluster": "Myeloid", "coef_abs": 0.75},
        {"source": "SPI1", "target": "CSF1R", "cluster": "Myeloid", "coef_abs": 0.6},
        # CEBPA targets in Myeloid
        {"source": "CEBPA", "target": "MPO", "cluster": "Myeloid", "coef_abs": 0.7},
        {"source": "CEBPA", "target": "NE", "cluster": "Myeloid", "coef_abs": 0.65},
    ]
    return pd.DataFrame(data)


@pytest.fixture
def sample_module_genes():
    """Hotspot module gene sets."""
    return {
        "Module_1": {"GATA1", "HBB", "HBA1", "ALAS2", "KLF1"},  # Erythroid module
        "Module_2": {"SPI1", "CEBPA", "MPO", "LYZ", "CSF1R", "NE"},  # Myeloid module
    }


@pytest.fixture
def sample_adata():
    """Create minimal AnnData with cluster annotations and expression data."""
    genes = ["GATA1", "SPI1", "KLF1", "CEBPA", "STAT3", "HBB", "HBA1", "MPO", "LYZ"]
    n_ery = 20
    n_mye = 20
    n_cells = n_ery + n_mye
    
    X = np.zeros((n_cells, len(genes)))
    # GATA1, KLF1 high in Erythroid
    X[:n_ery, 0] = 3.5  # GATA1
    X[:n_ery, 2] = 2.8  # KLF1
    X[:n_ery, 5] = 4.0  # HBB
    X[:n_ery, 6] = 3.8  # HBA1
    
    # SPI1, CEBPA high in Myeloid
    X[n_ery:, 1] = 3.2  # SPI1
    X[n_ery:, 3] = 2.5  # CEBPA
    X[n_ery:, 7] = 3.9  # MPO
    X[n_ery:, 8] = 3.6  # LYZ
    
    obs = pd.DataFrame(
        {"cell_type": ["Erythroid"] * n_ery + ["Myeloid"] * n_mye},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=genes)
    return AnnData(X=X, obs=obs, var=var)


def test_compute_cluster_pair_centrality_diff(sample_centrality_score_df):
    """Test differential centrality calculation between two clusters."""
    diff_df = compute_cluster_pair_centrality_diff(
        sample_centrality_score_df,
        cluster_a="Erythroid",
        cluster_b="Myeloid",
        cluster_key="cluster",
    )
    assert not diff_df.empty
    assert "delta_degree" in diff_df.columns
    assert "delta_betweenness" in diff_df.columns
    assert "delta_eigenvector" in diff_df.columns
    assert "composite_centrality_diff" in diff_df.columns

    # GATA1 should have strongly positive delta (Erythroid > Myeloid)
    gata1_row = diff_df[diff_df["gene"] == "GATA1"].iloc[0]
    assert gata1_row["delta_degree"] == pytest.approx(0.7, rel=1e-2)
    assert gata1_row["delta_eigenvector"] == pytest.approx(0.65, rel=1e-2)

    # SPI1 should have strongly negative delta (Erythroid < Myeloid)
    spi1_row = diff_df[diff_df["gene"] == "SPI1"].iloc[0]
    assert spi1_row["delta_degree"] == pytest.approx(-0.7, rel=1e-2)
    assert spi1_row["delta_eigenvector"] == pytest.approx(-0.7, rel=1e-2)


def test_compute_tf_hotspot_enrichment(sample_links_df, sample_module_genes):
    """Test TF Hotspot module membership and downstream target enrichment calculation."""
    enrichment_df = compute_tf_hotspot_enrichment(
        links_df=sample_links_df,
        module_genes=sample_module_genes,
        cluster="Erythroid",
    )
    assert not enrichment_df.empty
    assert "tf" in enrichment_df.columns
    assert "tf_in_modules" in enrichment_df.columns
    assert "top_target_module" in enrichment_df.columns
    assert "target_overlap_ratio" in enrichment_df.columns

    # GATA1 in Erythroid has targets in Module_1
    gata1_row = enrichment_df[enrichment_df["tf"] == "GATA1"].iloc[0]
    assert "Module_1" in gata1_row["tf_in_modules"]
    assert gata1_row["top_target_module"] == "Module_1"
    assert gata1_row["target_overlap_ratio"] > 0.5


def test_select_top_perturbed_tfs(
    sample_centrality_score_df, sample_links_df, sample_module_genes, sample_adata
):
    """Test end-to-end TF selection across compared clusters."""
    selection_df = select_top_perturbed_tfs(
        score_df=sample_centrality_score_df,
        links_df=sample_links_df,
        module_genes=sample_module_genes,
        adata=sample_adata,
        compare_clusters=[("Erythroid", "Myeloid")],
        cluster_key="cell_type",
        top_n_tfs=2,
        mode="directional",
    )
    assert not selection_df.empty
    assert len(selection_df) <= 4  # Up to 2 top TFs per direction
    assert "tf" in selection_df.columns
    assert "comparison" in selection_df.columns
    assert "cluster_high" in selection_df.columns
    assert "perturbation_type" in selection_df.columns
    assert "selection_reason" in selection_df.columns
    assert "composite_score" in selection_df.columns

    # Top drivers for Erythroid should include GATA1
    ery_drivers = selection_df[selection_df["cluster_high"] == "Erythroid"]["tf"].tolist()
    assert "GATA1" in ery_drivers

    # Top drivers for Myeloid should include SPI1
    mye_drivers = selection_df[selection_df["cluster_high"] == "Myeloid"]["tf"].tolist()
    assert "SPI1" in mye_drivers

    # Check that selection rationale string is informative
    gata1_reason = selection_df[selection_df["tf"] == "GATA1"]["selection_reason"].iloc[0]
    assert "centrality" in gata1_reason.lower()
    assert "module" in gata1_reason.lower()


def test_select_top_perturbed_tfs_auto_pairs(sample_centrality_score_df, sample_links_df):
    """Test automatic pair selection when compare_clusters is None."""
    selection_df = select_top_perturbed_tfs(
        score_df=sample_centrality_score_df,
        links_df=sample_links_df,
        compare_clusters=None,
        cluster_key="cluster",
        top_n_tfs=1,
    )
    assert not selection_df.empty
    assert "comparison" in selection_df.columns
