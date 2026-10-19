"""
Unit Tests for Perturbation Plotting Module
============================================

Tests figure generation in genecircuitry.plotting.perturbation_plots.
"""

import os
import shutil
import tempfile
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from genecircuitry import config
from genecircuitry.plotting.perturbation_plots import (
    plot_tf_selection_rationale,
    plot_perturbation_vector_field,
    plot_perturbation_developmental_flow,
    plot_cluster_shift_magnitude,
    plot_target_gene_shifts,
    plot_comparative_perturbation_dashboard,
    generate_all_perturbation_plots,
)


class MockOracleForPlot:
    """Mock Oracle object for plotting tests."""

    def __init__(self):
        n_cells = 30
        genes = ["GATA1", "SPI1", "HBB", "MPO"]
        X = np.random.uniform(0.1, 4.0, (n_cells, len(genes)))
        coords = np.random.randn(n_cells, 2)
        obs = pd.DataFrame({
            "cell_type": ["Erythroid"] * 15 + ["Myeloid"] * 15,
            "pseudotime": np.linspace(0.0, 1.0, n_cells),
        }, index=[f"cell_{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=genes)
        self.adata = AnnData(X=X, obs=obs, var=var)
        self.adata.obsm["X_draw_graph_fa"] = coords
        self.cluster_column_name = "cell_type"
        self.embedding_name = "X_draw_graph_fa"

        n_steps = 20
        self.flow = {
            "grid_x": np.linspace(-3, 3, n_steps),
            "grid_y": np.linspace(-3, 3, n_steps),
            "flow_x": np.ones((n_steps, n_steps)) * 0.1,
            "flow_y": np.ones((n_steps, n_steps)) * 0.1,
        }
        self.delta_embedding = np.ones((n_cells, 2)) * 0.2


class MockDevModuleForPlot:
    """Mock Oracle_development_module for developmental layout plotting."""

    def __init__(self, oracle):
        self.oracle = oracle
        self.inner_product = np.random.normal(loc=0.5, scale=0.2, size=oracle.adata.n_obs)
        self.inner_product_df = pd.DataFrame({
            "cell_id": oracle.adata.obs_names,
            "inner_product": self.inner_product,
            "cluster": oracle.adata.obs[oracle.cluster_column_name],
        })


@pytest.fixture
def temp_pert_fig_dir():
    """Create a temporary output directory and update config paths."""
    d = tempfile.mkdtemp(prefix="test_pert_plots_")
    old_output = config.OUTPUT_DIR
    old_figures = config.FIGURES_DIR

    pert_fig_dir = os.path.join(d, "figures", "perturbation")
    os.makedirs(pert_fig_dir, exist_ok=True)

    config.update_config(
        OUTPUT_DIR=d,
        FIGURES_DIR=os.path.join(d, "figures"),
        SAVE_PDF=True,
    )

    yield d, pert_fig_dir

    config.update_config(
        OUTPUT_DIR=old_output,
        FIGURES_DIR=old_figures,
    )
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def sample_selection_df():
    """Sample selection table for plotting."""
    return pd.DataFrame([
        {
            "tf": "GATA1",
            "comparison": "Erythroid vs Myeloid",
            "cluster_high": "Erythroid",
            "cluster_low": "Myeloid",
            "perturbation_type": "KO",
            "delta_degree": 0.7,
            "delta_betweenness": 0.55,
            "delta_eigenvector": 0.65,
            "tf_hotspot_module": "Module_1",
            "top_target_module": "Module_1",
            "target_module_overlap_ratio": 0.8,
            "sample_target_genes": "HBB, HBA1, ALAS2",
            "log2fc": 3.5,
            "composite_score": 0.92,
            "selection_reason": "Primary driver of Erythroid (+0.70 deg, +0.65 eig, 80% targets in Module_1).",
        },
        {
            "tf": "KLF1",
            "comparison": "Erythroid vs Myeloid",
            "cluster_high": "Erythroid",
            "cluster_low": "Myeloid",
            "perturbation_type": "KO",
            "delta_degree": 0.5,
            "delta_betweenness": 0.4,
            "delta_eigenvector": 0.5,
            "tf_hotspot_module": "Module_1",
            "top_target_module": "Module_1",
            "target_module_overlap_ratio": 0.65,
            "sample_target_genes": "HBB, HBA1",
            "log2fc": 2.8,
            "composite_score": 0.78,
            "selection_reason": "Secondary driver of Erythroid.",
        },
    ])


@pytest.fixture
def sample_shift_df():
    """Sample displacement score table."""
    return pd.DataFrame({
        "cell_id": [f"cell_{i}" for i in range(20)],
        "cluster_source": ["Erythroid"] * 20,
        "cluster_target": ["Myeloid"] * 20,
        "displacement_score": np.random.normal(loc=0.4, scale=0.15, size=20),
        "shift_magnitude": np.random.uniform(0.1, 0.5, size=20),
    })


@pytest.fixture
def sample_delta_x_df():
    """Sample gene shift Delta X table."""
    return pd.DataFrame([
        {"tf": "GATA1", "target": "HBB", "cluster": "Erythroid", "mean_delta_x": -1.8, "abs_delta_x": 1.8, "direction": "downregulated"},
        {"tf": "GATA1", "target": "HBA1", "cluster": "Erythroid", "mean_delta_x": -1.5, "abs_delta_x": 1.5, "direction": "downregulated"},
        {"tf": "GATA1", "target": "ALAS2", "cluster": "Erythroid", "mean_delta_x": -1.2, "abs_delta_x": 1.2, "direction": "downregulated"},
        {"tf": "GATA1", "target": "MPO", "cluster": "Erythroid", "mean_delta_x": 0.4, "abs_delta_x": 0.4, "direction": "upregulated"},
    ])


def test_plot_tf_selection_rationale(sample_selection_df, temp_pert_fig_dir):
    """Test TF selection rationale panel generation."""
    _, pert_dir = temp_pert_fig_dir
    save_path = os.path.join(pert_dir, "tf_selection_rationale.png")

    fig = plot_tf_selection_rationale(
        sample_selection_df,
        comparison="Erythroid vs Myeloid",
        save_path=save_path,
    )
    assert fig is not None
    assert os.path.exists(save_path)


def test_plot_perturbation_vector_field(temp_pert_fig_dir):
    """Test 2D vector field grid arrow plot generation."""
    _, pert_dir = temp_pert_fig_dir
    oracle = MockOracleForPlot()
    save_path = os.path.join(pert_dir, "GATA1_KO_vector_field.png")

    fig = plot_perturbation_vector_field(
        oracle,
        tf="GATA1",
        mode="KO",
        cluster_key="cell_type",
        save_path=save_path,
    )
    assert fig is not None
    assert os.path.exists(save_path)


def test_plot_cluster_shift_magnitude(sample_shift_df, temp_pert_fig_dir):
    """Test cluster shift magnitude and transition score plot."""
    _, pert_dir = temp_pert_fig_dir
    save_path = os.path.join(pert_dir, "cluster_shift_magnitude.png")

    fig = plot_cluster_shift_magnitude(
        sample_shift_df,
        comparison="Erythroid vs Myeloid",
        save_path=save_path,
    )
    assert fig is not None
    assert os.path.exists(save_path)


def test_plot_target_gene_shifts(sample_delta_x_df, temp_pert_fig_dir):
    """Test downstream target expression shifts Delta X plot."""
    _, pert_dir = temp_pert_fig_dir
    save_path = os.path.join(pert_dir, "GATA1_target_shifts.png")

    fig = plot_target_gene_shifts(
        sample_delta_x_df,
        tf="GATA1",
        save_path=save_path,
    )
    assert fig is not None
    assert os.path.exists(save_path)


def test_plot_perturbation_developmental_flow(temp_pert_fig_dir):
    """Test developmental flow and inner product multi-panel layout."""
    _, pert_dir = temp_pert_fig_dir
    oracle = MockOracleForPlot()
    dev = MockDevModuleForPlot(oracle)
    save_path = os.path.join(pert_dir, "GATA1_developmental_flow.png")

    fig = plot_perturbation_developmental_flow(
        dev_module=dev,
        oracle=oracle,
        tf="GATA1",
        save_path=save_path,
    )
    assert fig is not None
    assert os.path.exists(save_path)


def test_plot_comparative_perturbation_dashboard(temp_pert_fig_dir):
    """Test multi-TF comparative summary dashboard."""
    _, pert_dir = temp_pert_fig_dir
    oracle = MockOracleForPlot()
    save_path = os.path.join(pert_dir, "comparative_perturbation_dashboard.png")

    results_dict = {
        "GATA1": {"oracle": oracle, "displacement": 0.45, "direction": "Erythroid -> Myeloid"},
        "KLF1": {"oracle": oracle, "displacement": 0.32, "direction": "Erythroid -> Myeloid"},
    }

    fig = plot_comparative_perturbation_dashboard(
        results_dict,
        comparison="Erythroid vs Myeloid",
        save_path=save_path,
    )
    assert fig is not None
    assert os.path.exists(save_path)
