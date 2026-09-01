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
    run_parallel_tasks,
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

from .comparative_plots import (
    plot_comparative_module_activity,
    plot_comparative_pathway_enrichment,
    plot_tf_module_regulatory_matrix,
    plot_comparative_tf_centrality,
    plot_differential_tf_targets,
    plot_comparative_summary_dashboard,
    plot_module_overlap_heatmap,
    plot_module_tf_regulatory_network,
    plot_gene_selection_sankey,
    plot_cross_cluster_regulatory_comparison,
    plot_tf_module_concordance,
    plot_cross_stratification_module_overlap,
    plot_integrated_regulatory_dashboard,
    generate_all_comparative_plots,
)

__all__ = [
    # Utilities
    "PlotLogger",
    "get_plot_logger",
    "plot_exists",
    "save_plot",
    "get_plot_registry",
    "run_parallel_tasks",
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
    # Comparative plots
    "plot_comparative_module_activity",
    "plot_comparative_pathway_enrichment",
    "plot_tf_module_regulatory_matrix",
    "plot_comparative_tf_centrality",
    "plot_differential_tf_targets",
    "plot_comparative_summary_dashboard",
    "plot_module_overlap_heatmap",
    "plot_module_tf_regulatory_network",
    "plot_gene_selection_sankey",
    "plot_cross_cluster_regulatory_comparison",
    "plot_tf_module_concordance",
    "plot_cross_stratification_module_overlap",
    "plot_integrated_regulatory_dashboard",
    "generate_all_comparative_plots",
]

