"""
Tests that plotting functions in grn_deep_analysis and hotspot_processing
correctly delegate to the canonical genecircuitry.plotting subpackage.
"""

import importlib
import inspect
import unittest
try:
    import pytest
except ImportError:
    pytest = None


def test_grn_deep_analysis_no_inline_plt():
    """grn_deep_analysis must not contain any direct plt.savefig or plt.figure calls."""
    import genecircuitry.grn_deep_analysis as mod

    source = inspect.getsource(mod)
    assert "plt.savefig" not in source, "grn_deep_analysis still has plt.savefig"
    assert "plt.figure" not in source, "grn_deep_analysis still has plt.figure"
    assert "plt.show" not in source, "grn_deep_analysis still has plt.show"


def test_hotspot_processing_no_inline_plt():
    """hotspot_processing must not contain any direct plt.savefig or plt.figure calls."""
    try:
        import genecircuitry.hotspot_processing as mod
    except ImportError:
        if pytest is not None:
            pytest.skip("hotspot optional dependency not installed")
        return

    source = inspect.getsource(mod)
    assert "plt.savefig" not in source, "hotspot_processing still has plt.savefig"
    assert "plt.show" not in source, "hotspot_processing still has plt.show"


def test_grn_deep_analysis_delegates_plot_network_graph():
    """plot_network_graph in grn_deep_analysis should delegate to grn_plots."""
    import genecircuitry.grn_deep_analysis as legacy
    from genecircuitry.plotting import grn_plots as canonical

    # The legacy wrapper should import the canonical implementation
    source = inspect.getsource(legacy.plot_network_graph)
    assert "grn_plots" in source or "_impl" in source


def test_grn_deep_analysis_delegates_plot_scatter_scores():
    """plot_scatter_scores in grn_deep_analysis should delegate to grn_plots."""
    import genecircuitry.grn_deep_analysis as legacy

    source = inspect.getsource(legacy.plot_scatter_scores)
    assert "grn_plots" in source or "_impl" in source


def test_grn_deep_analysis_delegates_plot_heatmap_scores():
    """plot_heatmap_scores in grn_deep_analysis should delegate to grn_plots."""
    import genecircuitry.grn_deep_analysis as legacy

    source = inspect.getsource(legacy.plot_heatmap_scores)
    assert "grn_plots" in source or "_impl" in source


def test_grn_deep_analysis_delegates_plot_difference_cluster_scores():
    """plot_difference_cluster_scores should delegate to grn_plots."""
    import genecircuitry.grn_deep_analysis as legacy

    source = inspect.getsource(legacy.plot_difference_cluster_scores)
    assert "grn_plots" in source or "_impl" in source


def test_grn_deep_analysis_delegates_plot_compare_cluster_scores():
    """plot_compare_cluster_scores should delegate to grn_plots."""
    import genecircuitry.grn_deep_analysis as legacy

    source = inspect.getsource(legacy.plot_compare_cluster_scores)
    assert "grn_plots" in source or "_impl" in source


def test_hotspot_processing_delegates_plot_hotspot_annotation():
    """plot_hotspot_annotation in hotspot_processing should delegate to hotspot_plots."""
    try:
        import genecircuitry.hotspot_processing as mod
    except ImportError:
        if pytest is not None:
            pytest.skip("hotspot optional dependency not installed")
        return

    source = inspect.getsource(mod.plot_hotspot_annotation)
    assert "hotspot_plots" in source or "_impl" in source


def test_hotspot_processing_delegates_plot_module_scores_violin():
    """plot_module_scores_violin in hotspot_processing should delegate to hotspot_plots."""
    try:
        import genecircuitry.hotspot_processing as mod
    except ImportError:
        if pytest is not None:
            pytest.skip("hotspot optional dependency not installed")
        return

    source = inspect.getsource(mod.plot_module_scores_violin)
    assert "hotspot_plots" in source or "_impl" in source


def test_grn_deep_analysis_no_deprecated_plot_exists():
    """grn_deep_analysis should not define a local _plot_exists function."""
    import genecircuitry.grn_deep_analysis as mod

    assert not hasattr(mod, "_plot_exists"), (
        "_plot_exists deprecated wrapper still present in grn_deep_analysis"
    )


def test_grn_deep_analysis_no_plot_score_comparison_2d():
    """plot_score_comparison_2D (legacy, had plt.show) should be removed."""
    import genecircuitry.grn_deep_analysis as mod

    assert not hasattr(mod, "plot_score_comparison_2D"), (
        "Legacy plot_score_comparison_2D still exported from grn_deep_analysis"
    )


def test_canonical_grn_plots_has_no_plt_show():
    """Canonical grn_plots must not call plt.show()."""
    from genecircuitry.plotting import grn_plots

    source = inspect.getsource(grn_plots)
    assert "plt.show" not in source, "grn_plots.py still calls plt.show()"


def test_canonical_hotspot_plots_uses_save_plot():
    """Canonical hotspot_plots should use save_plot() not plt.savefig()."""
    from genecircuitry.plotting import hotspot_plots

    source = inspect.getsource(hotspot_plots)
    assert "save_plot" in source, "hotspot_plots.py does not use save_plot()"
    assert "plt.savefig" not in source, "hotspot_plots.py still calls plt.savefig()"


def test_grn_deep_analysis_imports_are_clean():
    """grn_deep_analysis should not import matplotlib, seaborn, or networkx at module level."""
    import genecircuitry.grn_deep_analysis as mod

    # These heavy plotting libs should not be imported at module level
    assert not hasattr(mod, "plt"), "matplotlib.pyplot imported at module level"
    assert not hasattr(mod, "sns"), "seaborn imported at module level"
    assert not hasattr(mod, "nx"), "networkx imported at module level"


class TestPlottingDelegation(unittest.TestCase):
    def test_grn_deep_analysis_no_inline_plt(self):
        test_grn_deep_analysis_no_inline_plt()

    def test_hotspot_processing_no_inline_plt(self):
        test_hotspot_processing_no_inline_plt()

    def test_grn_deep_analysis_delegates_plot_network_graph(self):
        test_grn_deep_analysis_delegates_plot_network_graph()

    def test_grn_deep_analysis_delegates_plot_scatter_scores(self):
        test_grn_deep_analysis_delegates_plot_scatter_scores()

    def test_grn_deep_analysis_delegates_plot_heatmap_scores(self):
        test_grn_deep_analysis_delegates_plot_heatmap_scores()

    def test_grn_deep_analysis_delegates_plot_difference_cluster_scores(self):
        test_grn_deep_analysis_delegates_plot_difference_cluster_scores()

    def test_grn_deep_analysis_delegates_plot_compare_cluster_scores(self):
        test_grn_deep_analysis_delegates_plot_compare_cluster_scores()

    def test_hotspot_processing_delegates_plot_hotspot_annotation(self):
        test_hotspot_processing_delegates_plot_hotspot_annotation()

    def test_hotspot_processing_delegates_plot_module_scores_violin(self):
        test_hotspot_processing_delegates_plot_module_scores_violin()

    def test_grn_deep_analysis_no_deprecated_plot_exists(self):
        test_grn_deep_analysis_no_deprecated_plot_exists()

    def test_grn_deep_analysis_no_plot_score_comparison_2d(self):
        test_grn_deep_analysis_no_plot_score_comparison_2d()

    def test_canonical_grn_plots_has_no_plt_show(self):
        test_canonical_grn_plots_has_no_plt_show()

    def test_canonical_hotspot_plots_uses_save_plot(self):
        test_canonical_hotspot_plots_uses_save_plot()

    def test_grn_deep_analysis_imports_are_clean(self):
        test_grn_deep_analysis_imports_are_clean()


if __name__ == "__main__":
    unittest.main()
