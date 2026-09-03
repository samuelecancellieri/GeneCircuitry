"""
Tests for parallel plotting refactoring in genecircuitry.plotting.
Verifies serial (n_jobs=1) and parallel (n_jobs=2) execution, task execution,
logger atomicity, and GRN/QC/Hotspot plot generation with synthetic data.
"""

import os
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock
from anndata import AnnData

from genecircuitry import config
from genecircuitry.plotting.utils import (
    _resolve_n_jobs,
    run_parallel_tasks,
    PlotLogger,
    get_plot_logger,
)
from genecircuitry.plotting.grn_plots import (
    plot_heatmap_scores,
    plot_scatter_scores,
    plot_difference_cluster_scores,
    plot_compare_cluster_scores,
    plot_network_graph,
    plot_enriched_tf_network,
    plot_tf_shared_target_network,
    generate_all_grn_plots,
)
from genecircuitry.plotting.qc_plots import (
    generate_all_qc_plots,
)
from genecircuitry.plotting.hotspot_plots import (
    plot_hotspot_local_correlations,
    plot_hotspot_annotation,
    plot_module_scores_violin,
    generate_all_hotspot_plots,
)


def _dummy_task_worker(task):
    """Top-level worker for testing run_parallel_tasks."""
    return task["val"] * 2


def _failing_task_worker(task):
    """Top-level worker that raises exception on odd values."""
    if task["val"] % 2 != 0:
        raise ValueError(f"Odd value error: {task['val']}")
    return task["val"] * 10


class TestResolveNJobs(unittest.TestCase):
    def test_resolve_none_uses_config(self):
        n = _resolve_n_jobs(None)
        self.assertIsInstance(n, int)
        self.assertGreater(n, 0)

    def test_resolve_negative_one(self):
        expected = os.cpu_count() or 1
        self.assertEqual(_resolve_n_jobs(-1), expected)

    def test_resolve_explicit(self):
        self.assertEqual(_resolve_n_jobs(4), 4)

    def test_resolve_zero_or_negative(self):
        self.assertEqual(_resolve_n_jobs(0), 1)
        self.assertEqual(_resolve_n_jobs(-2), 1)


class TestRunParallelTasks(unittest.TestCase):
    def test_empty_tasks(self):
        res = run_parallel_tasks(_dummy_task_worker, [], n_jobs=2)
        self.assertEqual(res, [])

    def test_single_task(self):
        res = run_parallel_tasks(_dummy_task_worker, [{"val": 5}], n_jobs=2)
        self.assertEqual(res, [10])

    def test_serial_execution(self):
        tasks = [{"val": 1}, {"val": 2}, {"val": 3}]
        res = run_parallel_tasks(_dummy_task_worker, tasks, n_jobs=1)
        self.assertEqual(res, [2, 4, 6])

    def test_parallel_execution(self):
        tasks = [{"val": 1}, {"val": 2}, {"val": 3}, {"val": 4}]
        res = run_parallel_tasks(_dummy_task_worker, tasks, n_jobs=2)
        self.assertEqual(res, [2, 4, 6, 8])

    def test_failing_task_returns_none(self):
        tasks = [{"val": 2}, {"val": 3}, {"val": 4}]
        res = run_parallel_tasks(_failing_task_worker, tasks, n_jobs=2)
        self.assertEqual(res, [20, None, 40])


class TestPlotLoggerBatchAndAtomic(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_register_plots_batch_and_save(self):
        logger = PlotLogger(output_dir=self.temp_dir.name, log_file="test_registry.json")
        p1 = os.path.join(self.temp_dir.name, "fig1.png")
        p2 = os.path.join(self.temp_dir.name, "fig2.png")
        p3 = os.path.join(self.temp_dir.name, "fig3.png")
        batch = [
            (p1, "qc", {"name": "violin1"}),
            (p2, "grn", {"name": "scatter1"}),
            (p3, "hotspot", {"name": "heatmap1"}),
        ]
        logger.register_plots_batch(batch)
        self.assertEqual(logger.get_plot_count(), 3)
        logger.save()
        self.assertTrue(os.path.exists(logger.log_file))

        # Re-read with new logger instance to test persistence
        logger2 = PlotLogger(output_dir=self.temp_dir.name, log_file="test_registry.json")
        self.assertEqual(logger2.get_plot_count(), 3)
        self.assertTrue(logger2.is_registered(p1))


class TestParallelGRNPlots(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.orig_figures_grn = config.FIGURES_DIR_GRN
        self.orig_output_dir = config.OUTPUT_DIR
        config.FIGURES_DIR_GRN = os.path.join(self.temp_dir.name, "figures", "grn")
        config.OUTPUT_DIR = os.path.join(self.temp_dir.name, "output")
        os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)
        os.makedirs(f"{config.OUTPUT_DIR}/grn_deep_analysis", exist_ok=True)

        # Synthetic scores DataFrame
        clusters = ["Cluster_A", "Cluster_B", "Cluster_C"]
        genes = [f"Gene_{i}" for i in range(25)]
        rows = []
        for c in clusters:
            for idx, g in enumerate(genes):
                # Ensure Gene_0 and Gene_1 have high scores (above 90th percentile)
                base = 100.0 if idx < 3 else 1.0
                rows.append({
                    "gene": g,
                    "cluster": c,
                    "stratification": "Combined",
                    "degree_all": base * (np.random.rand() + 1),
                    "degree_centrality_all": base * (np.random.rand() + 1),
                    "degree_centrality_in": base * (np.random.rand() + 1),
                    "degree_centrality_out": base * (np.random.rand() + 1),
                    "betweenness_centrality": base * (np.random.rand() + 1),
                    "eigenvector_centrality": base * (np.random.rand() + 1),
                    "degree_out": int(base),
                })
        self.score_df = pd.DataFrame(rows)
        self.score_df = self.score_df.set_index("gene", drop=False)

        # Synthetic links DataFrame with sources from top-scoring genes
        link_rows = []
        for c in clusters:
            for i in range(15):
                link_rows.append({
                    "source": f"Gene_{i % 3}",
                    "target": f"Gene_{(i + 1) % 20}",
                    "cluster": c,
                    "coef_abs": np.random.rand() + 0.1,
                })
        self.links_df = pd.DataFrame(link_rows)

    def tearDown(self):
        config.FIGURES_DIR_GRN = self.orig_figures_grn
        config.OUTPUT_DIR = self.orig_output_dir
        self.temp_dir.cleanup()

    def test_plot_heatmap_scores_serial_and_parallel(self):
        # Serial (n_jobs=1)
        count_serial = plot_heatmap_scores(
            self.score_df,
            scores=["degree_centrality_all"],
            skip_existing=False,
            n_jobs=1,
        )
        self.assertGreater(count_serial, 0)

        # Parallel (n_jobs=2) with skip_existing=False
        count_parallel = plot_heatmap_scores(
            self.score_df,
            scores=["degree_centrality_all"],
            skip_existing=False,
            n_jobs=2,
        )
        self.assertEqual(count_serial, count_parallel)

    def test_plot_scatter_scores_parallel(self):
        count = plot_scatter_scores(
            self.score_df,
            scores_list=["eigenvector_centrality"],
            skip_existing=False,
            n_jobs=2,
        )
        self.assertGreater(count, 0)

    def test_plot_difference_cluster_scores_parallel(self):
        count = plot_difference_cluster_scores(
            self.score_df,
            scores=["betweenness_centrality"],
            skip_existing=False,
            n_jobs=2,
        )
        self.assertGreater(count, 0)

    def test_plot_compare_cluster_scores_parallel(self):
        count = plot_compare_cluster_scores(
            self.score_df,
            scores=["degree_all"],
            skip_existing=False,
            n_jobs=2,
        )
        self.assertGreater(count, 0)

    def test_plot_network_graph_parallel(self):
        count = plot_network_graph(
            self.score_df,
            self.links_df,
            scores=["degree_all"],
            skip_existing=False,
            n_jobs=2,
        )
        self.assertGreater(count, 0)

    def test_generate_all_grn_plots_parallel(self):
        results = generate_all_grn_plots(
            self.score_df,
            self.links_df,
            scores=["degree_centrality_all"],
            gene_sets=[],
            skip_existing=False,
            n_jobs=2,
        )
        self.assertIn("scatter", results)
        self.assertIn("difference", results)
        self.assertIn("comparison", results)
        self.assertIn("heatmap", results)
        self.assertIn("network", results)
        self.assertGreater(sum(results.values()), 0)


class TestParallelQCPlots(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.orig_figures_qc = config.FIGURES_DIR_QC
        config.FIGURES_DIR_QC = os.path.join(self.temp_dir.name, "figures", "qc")
        os.makedirs(config.FIGURES_DIR_QC, exist_ok=True)

        n_cells = 50
        n_genes = 100
        X = np.random.poisson(1, size=(n_cells, n_genes)).astype(np.float32)
        obs = pd.DataFrame({
            "total_counts": np.random.randint(100, 1000, size=n_cells),
            "n_genes_by_counts": np.random.randint(20, 90, size=n_cells),
            "pct_counts_mt": np.random.uniform(0.1, 5.0, size=n_cells),
        })
        self.adata_pre = AnnData(X=X, obs=obs)
        self.adata_post = AnnData(X=X[:40], obs=obs.iloc[:40].copy())

    def tearDown(self):
        config.FIGURES_DIR_QC = self.orig_figures_qc
        self.temp_dir.cleanup()

    def test_generate_all_qc_plots_serial_and_parallel(self):
        res_serial = generate_all_qc_plots(
            self.adata_pre,
            self.adata_post,
            save_name="test_serial",
            skip_existing=False,
            n_jobs=1,
        )
        self.assertEqual(len(res_serial), 4)
        self.assertTrue(all(res_serial.values()))

        res_parallel = generate_all_qc_plots(
            self.adata_pre,
            self.adata_post,
            save_name="test_parallel",
            skip_existing=False,
            n_jobs=2,
        )
        self.assertEqual(len(res_parallel), 4)
        self.assertTrue(all(res_parallel.values()))


class TestParallelHotspotPlots(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.orig_figures_hotspot = config.FIGURES_DIR_HOTSPOT
        self.orig_output_dir = config.OUTPUT_DIR
        config.FIGURES_DIR_HOTSPOT = os.path.join(self.temp_dir.name, "figures", "hotspot")
        config.OUTPUT_DIR = os.path.join(self.temp_dir.name, "output")
        os.makedirs(config.FIGURES_DIR_HOTSPOT, exist_ok=True)
        os.makedirs(f"{config.OUTPUT_DIR}/hotspot", exist_ok=True)

        n_cells = 60
        cell_ids = [f"Cell_{i}" for i in range(n_cells)]
        clusters = ["Cluster_1", "Cluster_2", "Cluster_3"] * (n_cells // 3)

        obs = pd.DataFrame({"leiden": clusters}, index=cell_ids)
        self.adata = AnnData(X=np.zeros((n_cells, 10)), obs=obs)

        # Mock Hotspot object
        self.hotspot_obj = MagicMock()
        genes = [f"Gene_{i}" for i in range(5)]
        self.hotspot_obj.modules = pd.Series(
            [1, 1, 2, 2, -1],
            index=genes,
        )
        self.hotspot_obj.local_correlation_z = pd.DataFrame(
            np.array([
                [8.0, 3.5, -1.2, -2.0, 0.1],
                [3.5, 8.0, -1.0, -2.5, 0.0],
                [-1.2, -1.0, 8.0, 4.2, -0.3],
                [-2.0, -2.5, 4.2, 8.0, -0.1],
                [0.1, 0.0, -0.3, -0.1, 8.0],
            ]),
            index=genes,
            columns=genes,
        )
        self.hotspot_obj.linkage = None
        self.hotspot_obj.module_scores = pd.DataFrame(
            {
                1: np.random.randn(n_cells),
                2: np.random.randn(n_cells),
            },
            index=cell_ids,
        )

    def tearDown(self):
        config.FIGURES_DIR_HOTSPOT = self.orig_figures_hotspot
        config.OUTPUT_DIR = self.orig_output_dir
        self.temp_dir.cleanup()

    def test_plot_hotspot_local_correlations(self):
        # Generate plot
        res = plot_hotspot_local_correlations(
            self.hotspot_obj,
            skip_existing=False,
        )
        self.assertTrue(res)
        png_path = os.path.join(
            config.FIGURES_DIR_HOTSPOT, "hotspot_local_correlations.png"
        )
        self.assertTrue(os.path.exists(png_path))

        # Test skip existing
        res_skip = plot_hotspot_local_correlations(
            self.hotspot_obj,
            skip_existing=True,
        )
        self.assertFalse(res_skip)

    def test_plot_hotspot_annotation(self):
        res = plot_hotspot_annotation(
            self.hotspot_obj,
            gene_sets=[],
            skip_existing=False,
        )
        self.assertTrue(res)
        png_path = os.path.join(
            config.FIGURES_DIR_HOTSPOT,
            "hotspot_local_correlation_heatmap_with_annotations.png",
        )
        self.assertTrue(os.path.exists(png_path))

        # Test skip existing
        res_skip = plot_hotspot_annotation(
            self.hotspot_obj,
            gene_sets=[],
            skip_existing=True,
        )
        self.assertFalse(res_skip)

    def test_generate_all_hotspot_plots(self):
        results = generate_all_hotspot_plots(
            self.hotspot_obj,
            adata=self.adata,
            cluster_key="leiden",
            gene_sets=[],
            skip_existing=False,
        )
        self.assertIn("local_correlations", results)
        self.assertIn("annotated_heatmap", results)
        self.assertTrue(results["local_correlations"])
        self.assertTrue(results["annotated_heatmap"])

    def test_plot_module_scores_violin_parallel(self):
        results = plot_module_scores_violin(
            self.hotspot_obj,
            self.adata,
            cluster_key="leiden",
            gene_sets=[],
            skip_existing=False,
            n_jobs=2,
        )
        self.assertIn("per_cluster", results)
        self.assertIn("all_clusters", results)
        self.assertIn("horizontal", results)
        self.assertTrue(all(results.values()))


if __name__ == "__main__":
    unittest.main()
