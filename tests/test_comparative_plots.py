"""
Unit Tests for Comparative Plotting Module
==========================================

Tests plot generation and PDF + PNG file creation in genecircuitry.plotting.comparative_plots.
"""

import os
import shutil
import tempfile
import numpy as np
import pandas as pd
import pytest

from genecircuitry import config
from genecircuitry.plotting.comparative_plots import (
    plot_comparative_module_activity,
    plot_comparative_pathway_enrichment,
    plot_tf_module_regulatory_matrix,
    plot_comparative_tf_centrality,
    plot_differential_tf_targets,
    plot_comparative_summary_dashboard,
    generate_all_comparative_plots,
)


@pytest.fixture
def temp_plot_dir():
    """Create a temporary output directory and update config paths."""
    d = tempfile.mkdtemp(prefix="test_comparative_plots_")
    old_output = config.OUTPUT_DIR
    old_figures = config.FIGURES_DIR
    old_comp_figs = config.FIGURES_DIR_COMPARATIVE

    comp_fig_dir = os.path.join(d, "figures", "comparative")
    os.makedirs(comp_fig_dir, exist_ok=True)

    config.update_config(
        OUTPUT_DIR=d,
        FIGURES_DIR=os.path.join(d, "figures"),
        FIGURES_DIR_COMPARATIVE=comp_fig_dir,
        SAVE_PDF=True,
    )

    yield d, comp_fig_dir

    config.update_config(
        OUTPUT_DIR=old_output,
        FIGURES_DIR=old_figures,
        FIGURES_DIR_COMPARATIVE=old_comp_figs,
    )
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def sample_comparative_data():
    """Create synthetic comparative datasets for testing."""
    activity_df = pd.DataFrame(
        [[1.5, -0.5], [-0.2, 1.8]],
        index=["Module 1", "Module 2"],
        columns=["Cluster 0", "Cluster 1"],
    )

    enrichment_df = pd.DataFrame([
        {"module": "Module 1", "term": "Immune Response", "adjusted_p_value": 1e-5},
        {"module": "Module 2", "term": "Cell Cycle", "adjusted_p_value": 1e-4},
    ])

    tf_mod_matrix = pd.DataFrame(
        [[12, 1], [0, 15]],
        index=["STAT1", "MYC"],
        columns=["Module 1", "Module 2"],
    )

    tf_pivot_df = pd.DataFrame(
        [[0.8, 0.75], [0.1, 0.9]],
        index=["STAT1", "MYC"],
        columns=["Cluster 0", "Cluster 1"],
    )

    diff_targets_df = pd.DataFrame([
        {"tf": "STAT1", "shared_targets_count": 8, "specific_targets_count": 4},
        {"tf": "MYC", "shared_targets_count": 10, "specific_targets_count": 5},
    ])

    return {
        "module_activity": activity_df,
        "module_enrichment": enrichment_df,
        "tf_to_module_matrix": tf_mod_matrix,
        "tf_centrality": tf_pivot_df,
        "differential_tf_targets": diff_targets_df,
    }


def test_plot_comparative_module_activity(temp_plot_dir, sample_comparative_data):
    """Test module activity heatmap generates both PNG and PDF."""
    _, comp_fig_dir = temp_plot_dir
    success = plot_comparative_module_activity(
        sample_comparative_data["module_activity"],
        save_name="test",
        skip_existing=False,
    )
    assert success is True
    png_file = os.path.join(comp_fig_dir, "comparative_module_activity_test.png")
    pdf_file = os.path.join(comp_fig_dir, "comparative_module_activity_test.pdf")
    assert os.path.exists(png_file)
    assert os.path.exists(pdf_file)


def test_plot_comparative_pathway_enrichment(temp_plot_dir, sample_comparative_data):
    """Test pathway enrichment plot generates both PNG and PDF."""
    _, comp_fig_dir = temp_plot_dir
    success = plot_comparative_pathway_enrichment(
        sample_comparative_data["module_enrichment"],
        save_name="test",
        skip_existing=False,
    )
    assert success is True
    png_file = os.path.join(comp_fig_dir, "comparative_pathway_enrichment_test.png")
    pdf_file = os.path.join(comp_fig_dir, "comparative_pathway_enrichment_test.pdf")
    assert os.path.exists(png_file)
    assert os.path.exists(pdf_file)


def test_plot_tf_module_regulatory_matrix(temp_plot_dir, sample_comparative_data):
    """Test TF-to-Module matrix plot generates both PNG and PDF."""
    _, comp_fig_dir = temp_plot_dir
    success = plot_tf_module_regulatory_matrix(
        sample_comparative_data["tf_to_module_matrix"],
        save_name="test",
        skip_existing=False,
    )
    assert success is True
    png_file = os.path.join(comp_fig_dir, "tf_module_regulatory_matrix_test.png")
    pdf_file = os.path.join(comp_fig_dir, "tf_module_regulatory_matrix_test.pdf")
    assert os.path.exists(png_file)
    assert os.path.exists(pdf_file)


def test_plot_comparative_tf_centrality(temp_plot_dir, sample_comparative_data):
    """Test TF centrality heatmap generates both PNG and PDF."""
    _, comp_fig_dir = temp_plot_dir
    success = plot_comparative_tf_centrality(
        sample_comparative_data["tf_centrality"],
        save_name="test",
        skip_existing=False,
    )
    assert success is True
    png_file = os.path.join(comp_fig_dir, "comparative_tf_centrality_test.png")
    pdf_file = os.path.join(comp_fig_dir, "comparative_tf_centrality_test.pdf")
    assert os.path.exists(png_file)
    assert os.path.exists(pdf_file)


def test_plot_differential_tf_targets(temp_plot_dir, sample_comparative_data):
    """Test differential TF target plot generates both PNG and PDF."""
    _, comp_fig_dir = temp_plot_dir
    success = plot_differential_tf_targets(
        sample_comparative_data["differential_tf_targets"],
        save_name="test",
        skip_existing=False,
    )
    assert success is True
    png_file = os.path.join(comp_fig_dir, "differential_tf_targets_test.png")
    pdf_file = os.path.join(comp_fig_dir, "differential_tf_targets_test.pdf")
    assert os.path.exists(png_file)
    assert os.path.exists(pdf_file)


def test_plot_comparative_summary_dashboard(temp_plot_dir, sample_comparative_data):
    """Test summary dashboard generates both PNG and PDF."""
    _, comp_fig_dir = temp_plot_dir
    success = plot_comparative_summary_dashboard(
        sample_comparative_data,
        save_name="test",
        skip_existing=False,
    )
    assert success is True
    png_file = os.path.join(comp_fig_dir, "comparative_summary_dashboard_test.png")
    pdf_file = os.path.join(comp_fig_dir, "comparative_summary_dashboard_test.pdf")
    assert os.path.exists(png_file)
    assert os.path.exists(pdf_file)


def test_generate_all_comparative_plots(temp_plot_dir, sample_comparative_data):
    """Test batch generation of all comparative plots."""
    results = generate_all_comparative_plots(
        sample_comparative_data,
        save_name="all",
        skip_existing=False,
    )
    assert isinstance(results, dict)
    assert results.get("module_activity") is True
    assert results.get("pathway_enrichment") is True
    assert results.get("tf_module_matrix") is True
    assert results.get("tf_centrality") is True
    assert results.get("differential_targets") is True
    assert results.get("summary_dashboard") is True

