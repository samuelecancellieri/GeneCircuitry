"""
Unit Tests for Comparative Pipeline Controller Integration
===========================================================

Tests PipelineController execution of comparative analysis and report generation.
"""

import os
import shutil
import tempfile
from argparse import Namespace
from datetime import datetime
import numpy as np
import pandas as pd
from anndata import AnnData
import pytest

from genecircuitry import config
from genecircuitry.pipeline.controller import PipelineController


@pytest.fixture
def temp_run_dir():
    """Create a temporary run directory for controller execution."""
    d = tempfile.mkdtemp(prefix="test_controller_comp_")
    old_output = config.OUTPUT_DIR
    old_figures = config.FIGURES_DIR
    old_comp_figs = config.FIGURES_DIR_COMPARATIVE

    comp_dir = os.path.join(d, "comparative")
    fig_comp_dir = os.path.join(d, "figures", "comparative")
    os.makedirs(comp_dir, exist_ok=True)
    os.makedirs(fig_comp_dir, exist_ok=True)

    config.update_config(
        OUTPUT_DIR=d,
        FIGURES_DIR=os.path.join(d, "figures"),
        FIGURES_DIR_COMPARATIVE=fig_comp_dir,
        SAVE_PDF=True,
    )

    yield d

    config.update_config(
        OUTPUT_DIR=old_output,
        FIGURES_DIR=old_figures,
        FIGURES_DIR_COMPARATIVE=old_comp_figs,
    )
    shutil.rmtree(d, ignore_errors=True)


def test_controller_comparative_step_non_stratified(temp_run_dir):
    """Test controller running comparative step on a single dataset with clusters."""
    # Synthetic AnnData
    n_cells = 40
    n_genes = 15
    X = np.random.randn(n_cells, n_genes)
    obs = pd.DataFrame({
        "leiden": ["0"] * 20 + ["1"] * 20,
        "Module_1": np.random.uniform(0.5, 2.0, n_cells),
        "Module_2": np.random.uniform(-1.0, 0.5, n_cells),
    })
    var = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
    adata = AnnData(X=X, obs=obs, var=var)

    # Synthetic scores and links
    co_dir = os.path.join(temp_run_dir, "celloracle")
    os.makedirs(co_dir, exist_ok=True)

    scores_df = pd.DataFrame([
        {"gene": "TF1", "cluster": "0", "degree_centrality_all": 0.8},
        {"gene": "TF2", "cluster": "0", "degree_centrality_all": 0.2},
        {"gene": "TF1", "cluster": "1", "degree_centrality_all": 0.85},
        {"gene": "TF2", "cluster": "1", "degree_centrality_all": 0.9},
    ])
    scores_df.to_csv(os.path.join(co_dir, "grn_merged_scores.csv"), index=False)

    links_df = pd.DataFrame([
        {"source": "TF1", "target": "Gene_1", "cluster": "0"},
        {"source": "TF1", "target": "Gene_1", "cluster": "1"},
        {"source": "TF2", "target": "Gene_2", "cluster": "1"},
    ])
    import pickle
    with open(os.path.join(co_dir, "grn_filtered_links.pkl"), "wb") as f:
        pickle.dump(links_df, f)

    args = Namespace(
        input="dummy.h5ad",
        output=temp_run_dir,
        cluster_key="leiden",
        cluster_key_stratification=None,
        name="test_experiment",
        skip_hotspot=False,
        skip_celloracle=False,
        atac_peaks=False,
    )

    controller = PipelineController(args=args, start_time=datetime.now())
    comp_results = controller.run_step_comparative_analysis(
        adata=adata,
        output_dir=temp_run_dir,
    )

    assert comp_results is not None
    assert "module_activity" in comp_results
    assert not comp_results["module_activity"].empty

    # Check generated files
    assert os.path.exists(os.path.join(temp_run_dir, "comparative", "module_activity_matrix.csv"))
    assert os.path.exists(os.path.join(temp_run_dir, "figures", "comparative", "comparative_summary_dashboard.png"))
    assert os.path.exists(os.path.join(temp_run_dir, "figures", "comparative", "comparative_summary_dashboard.pdf"))

