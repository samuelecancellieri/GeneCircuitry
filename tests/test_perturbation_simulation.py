"""
Unit Tests for In-Silico Perturbation Simulation
=================================================

Tests genecircuitry.perturbation.simulation and development_analysis functions.
"""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from genecircuitry.perturbation.simulation import (
    setup_oracle_for_simulation,
    simulate_single_tf_perturbation,
    extract_target_gene_shifts,
)
from genecircuitry.perturbation.development_analysis import (
    compute_cluster_displacement_score,
    compute_perturbation_cluster_transitions,
)


class MockOracle:
    """Mock CellOracle Oracle object for fast unit testing."""

    def __init__(self, adata):
        self.adata = adata
        self.cluster_column_name = "cell_type"
        self.embedding_name = "X_draw_graph_fa"
        self.is_fitted = False
        self.perturb_condition = None
        self.flow = None
        self.delta_embedding = None
        self.total_velocity_mass = None

    def copy(self):
        new_obj = MockOracle(self.adata.copy())
        new_obj.cluster_column_name = self.cluster_column_name
        new_obj.embedding_name = self.embedding_name
        new_obj.is_fitted = self.is_fitted
        new_obj.perturb_condition = self.perturb_condition
        return new_obj

    def get_cluster_specific_TFdict_from_Links(self, links_object, ignore_warning=False):
        self.has_tf_dict = True

    def fit_GRN_for_simulation(self, GRN_unit="cluster", alpha=1, use_cluster_specific_TFdict=False, verbose_level=1):
        self.is_fitted = True

    def simulate_shift(self, perturb_condition=None, GRN_unit=None, n_propagation=3, **kwargs):
        self.perturb_condition = perturb_condition
        # Create simulated_count in layers: reduced expression for perturbed gene
        tf = list(perturb_condition.keys())[0]
        val = list(perturb_condition.values())[0]
        sim_data = self.adata.layers["imputed_count"].copy()
        if tf in self.adata.var_names:
            idx = self.adata.var_names.get_loc(tf)
            sim_data[:, idx] = val
        self.adata.layers["simulated_count"] = sim_data

    def estimate_transition_prob(self, n_neighbors=200, **kwargs):
        pass

    def calculate_embedding_shift(self, sigma_corr=0.05):
        # Synthetic 2D shift: move cells towards a direction
        n_cells = self.adata.n_obs
        self.delta_embedding = np.ones((n_cells, 2)) * 0.5

    def calculate_p_mass(self, smooth=0.8, n_grid=40, n_neighbors=200):
        self.total_velocity_mass = np.ones((40, 40))

    def calculate_grid_arrows(self, smooth=0.8, steps=(40, 40), n_neighbors=200):
        n_steps = steps[0]
        grid_x = np.linspace(-5, 5, n_steps)
        grid_y = np.linspace(-5, 5, n_steps)
        flow_x = np.ones((n_steps, n_steps)) * 0.1
        flow_y = np.ones((n_steps, n_steps)) * 0.1
        self.flow = {
            "grid_x": grid_x,
            "grid_y": grid_y,
            "flow_x": flow_x,
            "flow_y": flow_y,
        }

    def get_markov_simulation_cell_transition_table(self, cluster_column_name=None, end=-1, return_df=True):
        clusters = sorted(list(self.adata.obs[self.cluster_column_name].unique()))
        df = pd.DataFrame(
            [[0.7, 0.3], [0.2, 0.8]],
            index=clusters,
            columns=clusters,
        )
        return df


@pytest.fixture
def mock_oracle_adata():
    """Create test AnnData with embeddings and clusters."""
    n_cells = 40
    genes = ["GATA1", "SPI1", "HBB", "MPO"]
    X = np.random.uniform(0.1, 4.0, (n_cells, len(genes)))
    
    # 2D embedding
    obsm_coords = np.zeros((n_cells, 2))
    obsm_coords[:20] = np.random.normal(loc=[-2.0, -2.0], scale=0.5, size=(20, 2))  # Cluster 0
    obsm_coords[20:] = np.random.normal(loc=[2.0, 2.0], scale=0.5, size=(20, 2))   # Cluster 1

    obs = pd.DataFrame({
        "cell_type": ["Erythroid"] * 20 + ["Myeloid"] * 20,
        "pseudotime": np.linspace(0.0, 1.0, n_cells),
    }, index=[f"cell_{i}" for i in range(n_cells)])
    
    var = pd.DataFrame(index=genes)
    adata = AnnData(X=X, obs=obs, var=var)
    adata.layers["imputed_count"] = X.copy()
    adata.obsm["X_draw_graph_fa"] = obsm_coords
    return adata


def test_setup_oracle_for_simulation(mock_oracle_adata):
    """Test setting up oracle object with fitted GRNs."""
    oracle = MockOracle(mock_oracle_adata)
    mock_links = object()
    oracle = setup_oracle_for_simulation(oracle, mock_links, alpha=10)
    assert oracle.is_fitted is True


def test_simulate_single_tf_perturbation(mock_oracle_adata):
    """Test simulating a single TF knockout perturbation."""
    oracle = MockOracle(mock_oracle_adata)
    oracle.fit_GRN_for_simulation()

    res_oracle = simulate_single_tf_perturbation(
        oracle,
        tf="GATA1",
        condition_val=0.0,
        n_propagation=3,
        cluster_key="cell_type",
    )
    assert res_oracle is not None
    assert res_oracle.perturb_condition == {"GATA1": 0.0}
    assert res_oracle.delta_embedding is not None
    assert res_oracle.flow is not None


def test_extract_target_gene_shifts(mock_oracle_adata):
    """Test extracting downstream target shifts Delta X."""
    oracle = MockOracle(mock_oracle_adata)
    oracle.simulate_shift(perturb_condition={"GATA1": 0.0})

    links_df = pd.DataFrame([
        {"source": "GATA1", "target": "HBB", "cluster": "Erythroid"},
        {"source": "GATA1", "target": "MPO", "cluster": "Erythroid"},
    ])
    
    shift_df = extract_target_gene_shifts(
        oracle,
        tf="GATA1",
        links_df=links_df,
        cluster="Erythroid",
    )
    assert not shift_df.empty
    assert "target" in shift_df.columns
    assert "mean_delta_x" in shift_df.columns
    assert "mean_wt_expr" in shift_df.columns
    assert "mean_sim_expr" in shift_df.columns


def test_compute_cluster_displacement_score(mock_oracle_adata):
    """Test calculating projection score of cell shifts towards target cluster centroid."""
    oracle = MockOracle(mock_oracle_adata)
    oracle.calculate_embedding_shift()

    scores_df = compute_cluster_displacement_score(
        oracle,
        cluster_source="Erythroid",
        cluster_target="Myeloid",
        cluster_key="cell_type",
    )
    assert not scores_df.empty
    assert "cell_id" in scores_df.columns
    assert "displacement_score" in scores_df.columns
    assert "shift_magnitude" in scores_df.columns
    # All cells in source cluster should have a score
    assert len(scores_df) == 20


def test_compute_perturbation_cluster_transitions(mock_oracle_adata):
    """Test Markov transition table extraction for perturbation."""
    oracle = MockOracle(mock_oracle_adata)
    trans_df = compute_perturbation_cluster_transitions(
        oracle,
        cluster_key="cell_type",
    )
    assert not trans_df.empty
    assert "Erythroid" in trans_df.index
    assert "Myeloid" in trans_df.columns
