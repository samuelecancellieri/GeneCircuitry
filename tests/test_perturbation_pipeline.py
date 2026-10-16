"""
Unit Tests for Perturbation Pipeline Orchestration
===================================================

Tests genecircuitry.perturbation.pipeline.run_perturbation_analysis.
"""

import os
import shutil
import tempfile
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from genecircuitry import config
from genecircuitry.perturbation.pipeline import run_perturbation_analysis


class MockOracleForPipeline:
    """Mock Oracle object for pipeline execution."""

    def __init__(self, adata):
        self.adata = adata
        self.cluster_column_name = "cell_type"
        self.embedding_name = "X_draw_graph_fa"
        self.is_fitted = False
        self.perturb_condition = None
        self.flow = {
            "grid_x": np.linspace(-3, 3, 20),
            "grid_y": np.linspace(-3, 3, 20),
            "flow_x": np.ones((20, 20)) * 0.1,
            "flow_y": np.ones((20, 20)) * 0.1,
        }
        self.delta_embedding = np.ones((adata.n_obs, 2)) * 0.2

    def copy(self):
        new_obj = MockOracleForPipeline(self.adata.copy())
        new_obj.cluster_column_name = self.cluster_column_name
        new_obj.embedding_name = self.embedding_name
        new_obj.is_fitted = self.is_fitted
        return new_obj

    def get_cluster_specific_TFdict_from_Links(self, links_object, ignore_warning=False):
        pass

    def fit_GRN_for_simulation(self, **kwargs):
        self.is_fitted = True

    def simulate_shift(self, perturb_condition=None, **kwargs):
        self.perturb_condition = perturb_condition
        tf = list(perturb_condition.keys())[0]
        val = list(perturb_condition.values())[0]
        sim_data = self.adata.layers["imputed_count"].copy()
        if tf in self.adata.var_names:
            idx = self.adata.var_names.get_loc(tf)
            sim_data[:, idx] = val
        self.adata.layers["simulated_count"] = sim_data

    def estimate_transition_prob(self, **kwargs):
        pass

    def calculate_embedding_shift(self, **kwargs):
        pass

    def calculate_p_mass(self, **kwargs):
        pass

    def calculate_grid_arrows(self, **kwargs):
        pass

    def get_markov_simulation_cell_transition_table(self, **kwargs):
        return pd.DataFrame([[0.8, 0.2], [0.1, 0.9]], index=["Erythroid", "Myeloid"], columns=["Erythroid", "Myeloid"])


@pytest.fixture
def pipeline_temp_dir():
    """Create temp directory for pipeline output."""
    d = tempfile.mkdtemp(prefix="test_pert_pipeline_")
    old_output = config.OUTPUT_DIR
    config.update_config(OUTPUT_DIR=d)
    yield d
    config.update_config(OUTPUT_DIR=old_output)
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def pipeline_data():
    """Prepare mock data for end-to-end pipeline run."""
    n_cells = 40
    genes = ["GATA1", "SPI1", "HBB", "MPO"]
    X = np.random.uniform(0.1, 4.0, (n_cells, len(genes)))
    coords = np.random.randn(n_cells, 2)
    obs = pd.DataFrame({
        "cell_type": ["Erythroid"] * 20 + ["Myeloid"] * 20,
        "pseudotime": np.linspace(0.0, 1.0, n_cells),
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=genes)
    adata = AnnData(X=X, obs=obs, var=var)
    adata.layers["imputed_count"] = X.copy()
    adata.obsm["X_draw_graph_fa"] = coords

    oracle = MockOracleForPipeline(adata)
    links = object()

    score_df = pd.DataFrame([
        {"gene": "GATA1", "cluster": "Erythroid", "degree_centrality_all": 0.8, "betweenness_centrality": 0.6, "eigenvector_centrality": 0.75},
        {"gene": "GATA1", "cluster": "Myeloid", "degree_centrality_all": 0.1, "betweenness_centrality": 0.05, "eigenvector_centrality": 0.1},
        {"gene": "SPI1", "cluster": "Erythroid", "degree_centrality_all": 0.15, "betweenness_centrality": 0.1, "eigenvector_centrality": 0.1},
        {"gene": "SPI1", "cluster": "Myeloid", "degree_centrality_all": 0.85, "betweenness_centrality": 0.7, "eigenvector_centrality": 0.8},
    ])

    links_df = pd.DataFrame([
        {"source": "GATA1", "target": "HBB", "cluster": "Erythroid", "coef_abs": 0.8},
        {"source": "SPI1", "target": "MPO", "cluster": "Myeloid", "coef_abs": 0.85},
    ])

    module_genes = {
        "Module_1": {"GATA1", "HBB"},
        "Module_2": {"SPI1", "MPO"},
    }

    return oracle, links, score_df, links_df, module_genes, adata


def test_run_perturbation_analysis(pipeline_temp_dir, pipeline_data):
    """Test full perturbation analysis workflow."""
    oracle, links, score_df, links_df, module_genes, adata = pipeline_data

    res = run_perturbation_analysis(
        oracle=oracle,
        links=links,
        score_df=score_df,
        links_df=links_df,
        module_genes=module_genes,
        adata=adata,
        compare_clusters=[("Erythroid", "Myeloid")],
        cluster_key="cell_type",
        top_n_tfs=1,
        output_dir=pipeline_temp_dir,
        save_figures=True,
    )

    assert res is not None
    assert "selection_df" in res
    assert "simulation_summary_df" in res
    assert "generated_figures" in res

    # Verify generated output files
    pert_data_dir = os.path.join(pipeline_temp_dir, "perturbation")
    assert os.path.exists(os.path.join(pert_data_dir, "tf_selection_rationale.csv"))
    assert os.path.exists(os.path.join(pert_data_dir, "tf_selection_rationale.md"))
    assert os.path.exists(os.path.join(pert_data_dir, "simulation_summary.csv"))
    assert os.path.exists(os.path.join(pert_data_dir, "target_gene_shifts.csv"))

    pert_fig_dir = os.path.join(pipeline_temp_dir, "figures", "perturbation")
    assert os.path.exists(pert_fig_dir)
    assert len(res["generated_figures"]) > 0
