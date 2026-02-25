"""
GeneCircuitry Plotting Module
=======================

Centralized plotting functions for GeneCircuitry analysis pipeline.
Organized by analysis type: QC, GRN, and Hotspot.

Each plotting module provides:
- Plot generation functions
- Plot existence checking (to avoid overwrites)
- Plot logging for tracking generated figures
"""

from .utils import (
    PlotLogger,
    get_plot_logger,
    plot_exists,
    save_plot,
    get_plot_registry,
)

from .qc_plots import (
    plot_qc_violin_pre_filter,
    plot_qc_violin_post_filter,
    plot_qc_scatter_pre_filter,
    plot_qc_scatter_post_filter,
    generate_all_qc_plots,
)

from .grn_plots import (
    plot_network_graph,
    plot_enriched_tf_network,
    plot_tf_shared_target_network,
    plot_heatmap_scores,
    plot_scatter_scores,
    plot_difference_cluster_scores,
    plot_compare_cluster_scores,
    generate_all_grn_plots,
)

from .hotspot_plots import (
    plot_hotspot_local_correlations,
    plot_hotspot_annotation,
    plot_module_scores_violin,
    generate_all_hotspot_plots,
)

__all__ = [
    # Utilities
    "PlotLogger",
    "get_plot_logger",
    "plot_exists",
    "save_plot",
    "get_plot_registry",
    # QC plots
    "plot_qc_violin_pre_filter",
    "plot_qc_violin_post_filter",
    "plot_qc_scatter_pre_filter",
    "plot_qc_scatter_post_filter",
    "generate_all_qc_plots",
    # GRN plots
    "plot_network_graph",
    "plot_enriched_tf_network",
    "plot_tf_shared_target_network",
    "plot_heatmap_scores",
    "plot_scatter_scores",
    "plot_difference_cluster_scores",
    "plot_compare_cluster_scores",
    "generate_all_grn_plots",
    # Hotspot plots
    "plot_hotspot_local_correlations",
    "plot_hotspot_annotation",
    "plot_module_scores_violin",
    "generate_all_hotspot_plots",
]
