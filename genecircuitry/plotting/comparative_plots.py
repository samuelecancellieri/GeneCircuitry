"""
Comparative Plotting Module for GeneCircuitry
=============================================

Visualizations for cross-cluster and cross-stratification comparative analysis:
- Comparative module activity heatmaps
- Functional pathway enrichment summaries
- TF-to-Module regulatory connection heatmaps
- TF centrality & specificity comparative heatmaps
- Differential TF-Target Gene (TF-TG) conservation plots
- Integrated multi-panel summary dashboard
"""

import os
from typing import Optional, Dict, Any, List, Tuple, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .. import config
from .utils import save_plot, plot_exists

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


def _clean_save_suffix(save_name: str) -> str:
    """Format save name suffix cleanly."""
    return f"_{save_name}" if save_name and save_name != "default" else ""


def plot_comparative_module_activity(
    activity_df: pd.DataFrame,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Plot clustered heatmap of co-expression module activity scores across clusters/stratifications.

    Parameters
    ----------
    activity_df : pd.DataFrame
        Module activity DataFrame (rows = modules, columns = groups).
    save_name : str, default="default"
        Suffix for saved figure filename.
    figsize : tuple, optional
        Figure size (width, height).
    skip_existing : bool, default=True
        Skip if file already exists.

    Returns
    -------
    bool
        True if generated, False if skipped.
    """
    suffix = _clean_save_suffix(save_name)
    filepath = (
        f"{config.FIGURES_DIR_COMPARATIVE}/comparative_module_activity{suffix}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    if activity_df is None or activity_df.empty:
        return False

    if figsize is None:
        n_rows, n_cols = activity_df.shape
        w = max(8, min(24, 4 + n_cols * 1.2))
        h = max(6, min(20, 3 + n_rows * 0.8))
        figsize = (w, h)

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        activity_df,
        cmap=config.PLOT_COLOR_PALETTE,
        annot=False,
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={"label": "Mean Module Score"},
        ax=ax,
    )

    ax.set_title("Co-expression Module Activity Across Groups", fontsize=13, pad=12)
    ax.set_xlabel("Clusters / Stratifications", fontsize=11)
    ax.set_ylabel("Gene Module", fontsize=11)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="comparative",
        metadata={"plot_name": "module_activity", "n_modules": len(activity_df)},
        skip_existing=False,
    )


def plot_comparative_pathway_enrichment(
    enrichment_df: pd.DataFrame,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
    max_terms_per_module: int = 3,
) -> bool:
    """
    Plot summary dotplot / barplot of top biological pathways enriched in co-expression modules.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Table of module enrichment terms.
    save_name : str, default="default"
        Suffix for saved filename.
    figsize : tuple, optional
        Figure size.
    skip_existing : bool, default=True
        Skip if file exists.
    max_terms_per_module : int, default=3
        Max terms to show per module.

    Returns
    -------
    bool
        True if generated, False if skipped.
    """
    suffix = _clean_save_suffix(save_name)
    filepath = (
        f"{config.FIGURES_DIR_COMPARATIVE}/comparative_pathway_enrichment{suffix}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    if enrichment_df is None or enrichment_df.empty:
        return False

    # Filter top terms per module
    plot_data = (
        enrichment_df.groupby("module", observed=False)
        .head(max_terms_per_module)
        .copy()
    )

    if plot_data.empty:
        return False

    # Calculate -log10 p-value
    plot_data["neg_log_p"] = -np.log10(
        np.clip(plot_data["adjusted_p_value"], 1e-30, 1.0)
    )
    plot_data["label"] = plot_data["module"] + ": " + plot_data["term"]

    if figsize is None:
        h = max(6, min(24, 2 + len(plot_data) * 0.4))
        figsize = (10, h)

    fig, ax = plt.subplots(figsize=figsize)

    bars = ax.barh(
        range(len(plot_data)),
        plot_data["neg_log_p"],
        color=sns.color_palette(config.PLOT_COLOR_PALETTE, len(plot_data)),
        edgecolor="black",
        linewidth=0.5,
    )

    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels(plot_data["label"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
    ax.set_title(
        "Top Biological Pathways Enriched in Co-expression Modules", fontsize=13, pad=12
    )
    ax.axvline(
        x=-np.log10(0.05), color="red", linestyle="--", linewidth=1, label="p = 0.05"
    )
    ax.legend(loc="lower right")

    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="comparative",
        metadata={"plot_name": "pathway_enrichment", "n_terms": len(plot_data)},
        skip_existing=False,
    )


def plot_tf_module_regulatory_matrix(
    tf_module_matrix: pd.DataFrame,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Plot heatmap of TF-to-Module regulatory connections based on target gene overlaps.

    Parameters
    ----------
    tf_module_matrix : pd.DataFrame
        Matrix of TFs (rows) by Modules (columns) target counts.
    save_name : str, default="default"
        Suffix for saved filename.
    figsize : tuple, optional
        Figure size.
    skip_existing : bool, default=True
        Skip if file exists.

    Returns
    -------
    bool
        True if generated, False if skipped.
    """
    suffix = _clean_save_suffix(save_name)
    filepath = (
        f"{config.FIGURES_DIR_COMPARATIVE}/tf_module_regulatory_matrix{suffix}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    if tf_module_matrix is None or tf_module_matrix.empty:
        return False

    if figsize is None:
        n_rows, n_cols = tf_module_matrix.shape
        w = max(8, min(24, 4 + n_cols * 1.2))
        h = max(6, min(22, 3 + n_rows * 0.6))
        figsize = (w, h)

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        tf_module_matrix,
        cmap="YlOrRd",
        annot=False,
        fmt="d",
        linewidths=0.5,
        cbar_kws={"label": "Number of Target Genes"},
        ax=ax,
    )

    ax.set_title(
        "TF-to-Module Regulatory Mapping (Target Overlap)", fontsize=13, pad=12
    )
    ax.set_xlabel("Co-expression Module", fontsize=11)
    ax.set_ylabel("Transcription Factor (TF)", fontsize=11)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="comparative",
        metadata={"plot_name": "tf_to_module", "n_tfs": len(tf_module_matrix)},
        skip_existing=False,
    )


def plot_comparative_tf_centrality(
    tf_pivot_df: pd.DataFrame,
    tf_summary_df: Optional[pd.DataFrame] = None,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Plot comparative heatmap of top active TFs across clusters and stratifications.

    Parameters
    ----------
    tf_pivot_df : pd.DataFrame
        Matrix of TFs (rows) by Groups (columns) centrality scores.
    tf_summary_df : pd.DataFrame, optional
        Summary with classifications (Global Master vs Group-Specific).
    save_name : str, default="default"
        Suffix for saved filename.
    figsize : tuple, optional
        Figure size.
    skip_existing : bool, default=True
        Skip if file exists.

    Returns
    -------
    bool
        True if generated, False if skipped.
    """
    suffix = _clean_save_suffix(save_name)
    filepath = f"{config.FIGURES_DIR_COMPARATIVE}/comparative_tf_centrality{suffix}.png"

    if plot_exists(filepath, skip_existing):
        return False

    if tf_pivot_df is None or tf_pivot_df.empty:
        return False

    # Sort TFs by mean centrality
    sorted_df = tf_pivot_df.loc[
        tf_pivot_df.mean(axis=1).sort_values(ascending=False).index
    ]

    if figsize is None:
        n_rows, n_cols = sorted_df.shape
        w = max(8, min(24, 4 + n_cols * 1.2))
        h = max(7, min(24, 3 + n_rows * 0.55))
        figsize = (w, h)

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        sorted_df,
        cmap="Blues",
        annot=False,
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={"label": "Centrality Score"},
        ax=ax,
    )

    ax.set_title("Transcription Factor Centrality Across Groups", fontsize=13, pad=12)
    ax.set_xlabel("Clusters / Stratifications", fontsize=11)
    ax.set_ylabel("Transcription Factor", fontsize=11)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="comparative",
        metadata={"plot_name": "tf_centrality", "n_tfs": len(sorted_df)},
        skip_existing=False,
    )


def plot_differential_tf_targets(
    diff_targets_df: pd.DataFrame,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
    top_n: int = 15,
) -> bool:
    """
    Plot stacked bar chart of shared vs group-specific target genes for top TFs.

    Parameters
    ----------
    diff_targets_df : pd.DataFrame
        Table of differential TF target metrics.
    save_name : str, default="default"
        Suffix for saved filename.
    figsize : tuple, optional
        Figure size.
    skip_existing : bool, default=True
        Skip if file exists.
    top_n : int, default=15
        Number of top TFs to display.

    Returns
    -------
    bool
        True if generated, False if skipped.
    """
    suffix = _clean_save_suffix(save_name)
    filepath = f"{config.FIGURES_DIR_COMPARATIVE}/differential_tf_targets{suffix}.png"

    if plot_exists(filepath, skip_existing):
        return False

    if diff_targets_df is None or diff_targets_df.empty:
        return False

    plot_df = diff_targets_df.head(top_n).copy()

    if figsize is None:
        h = max(6, min(20, 2 + len(plot_df) * 0.5))
        figsize = (10, h)

    fig, ax = plt.subplots(figsize=figsize)

    y_pos = range(len(plot_df))
    ax.barh(
        y_pos,
        plot_df["shared_targets_count"],
        label="Conserved Targets (Shared)",
        color="#2b83ba",
        edgecolor="black",
        linewidth=0.5,
    )
    ax.barh(
        y_pos,
        plot_df["specific_targets_count"],
        left=plot_df["shared_targets_count"],
        label="Group-Specific Targets (Rewired)",
        color="#fdae61",
        edgecolor="black",
        linewidth=0.5,
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df["tf"], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Number of Target Genes", fontsize=11)
    ax.set_title(
        "TF Target Gene Conservation vs. Rewiring Across Groups", fontsize=13, pad=12
    )
    ax.legend(loc="lower right")

    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="comparative",
        metadata={"plot_name": "differential_targets", "n_tfs": len(plot_df)},
        skip_existing=False,
    )


def plot_comparative_summary_dashboard(
    comparative_results: Dict[str, Any],
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Plot an integrated multi-panel publication-ready summary dashboard bringing together
    (A) Module Activity, (B) TF Centrality, (C) TF-to-Module Mapping, and (D) Differential Targets.

    Parameters
    ----------
    comparative_results : dict
        Results dictionary returned by `run_comparative_analysis()`.
    save_name : str, default="default"
        Suffix for saved filename.
    figsize : tuple, optional
        Figure size.
    skip_existing : bool, default=True
        Skip if file exists.

    Returns
    -------
    bool
        True if generated, False if skipped.
    """
    suffix = _clean_save_suffix(save_name)
    filepath = (
        f"{config.FIGURES_DIR_COMPARATIVE}/comparative_summary_dashboard{suffix}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    if not comparative_results:
        return False

    activity_df = comparative_results.get("module_activity")
    tf_pivot_df = comparative_results.get("tf_centrality")
    tf_mod_matrix = comparative_results.get("tf_to_module_matrix")
    enrichment_df = comparative_results.get("module_enrichment")

    if figsize is None:
        figsize = (20, 15)

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Panel A: Module Activity
    if activity_df is not None and not activity_df.empty:
        sns.heatmap(
            activity_df,
            cmap=config.PLOT_COLOR_PALETTE,
            annot=False,
            fmt=".2f",
            linewidths=0.5,
            cbar_kws={"label": "Mean Score"},
            ax=axes[0, 0],
        )
        axes[0, 0].set_title(
            "A. Co-expression Module Activity", fontsize=12, fontweight="bold"
        )
        axes[0, 0].tick_params(axis="x", rotation=30)
    else:
        axes[0, 0].text(0.5, 0.5, "No Module Activity Data", ha="center", va="center")
        axes[0, 0].set_title(
            "A. Co-expression Module Activity", fontsize=12, fontweight="bold"
        )

    # Panel B: TF Centrality
    if tf_pivot_df is not None and not tf_pivot_df.empty:
        top_tfs_sub = tf_pivot_df.head(15)
        sns.heatmap(
            top_tfs_sub,
            cmap="Blues",
            annot=False,
            fmt=".2f",
            linewidths=0.5,
            cbar_kws={"label": "Centrality"},
            ax=axes[0, 1],
        )
        axes[0, 1].set_title(
            "B. Transcription Factor Centrality", fontsize=12, fontweight="bold"
        )
        axes[0, 1].tick_params(axis="x", rotation=30)
    else:
        axes[0, 1].text(0.5, 0.5, "No TF Centrality Data", ha="center", va="center")
        axes[0, 1].set_title(
            "B. Transcription Factor Centrality", fontsize=12, fontweight="bold"
        )

    # Panel C: TF-to-Module Regulation
    if tf_mod_matrix is not None and not tf_mod_matrix.empty:
        sns.heatmap(
            tf_mod_matrix.head(15),
            cmap="YlOrRd",
            annot=False,
            fmt="d",
            linewidths=0.5,
            cbar_kws={"label": "Target Count"},
            ax=axes[1, 0],
        )
        axes[1, 0].set_title(
            "C. TF-to-Module Regulation", fontsize=12, fontweight="bold"
        )
        axes[1, 0].tick_params(axis="x", rotation=30)
    else:
        axes[1, 0].text(0.5, 0.5, "No TF-to-Module Data", ha="center", va="center")
        axes[1, 0].set_title(
            "C. TF-to-Module Regulation", fontsize=12, fontweight="bold"
        )

    # Panel D: Top Functional Pathways
    if enrichment_df is not None and not enrichment_df.empty:
        top_enr = enrichment_df.head(12).copy()
        top_enr["neg_log_p"] = -np.log10(
            np.clip(top_enr["adjusted_p_value"], 1e-30, 1.0)
        )
        top_enr["label"] = top_enr["module"] + ": " + top_enr["term"]
        y_pos = range(len(top_enr))
        axes[1, 1].barh(
            y_pos,
            top_enr["neg_log_p"],
            color=sns.color_palette("viridis", len(top_enr)),
            edgecolor="black",
            linewidth=0.5,
        )
        axes[1, 1].set_yticks(y_pos)
        axes[1, 1].set_yticklabels(top_enr["label"], fontsize=8)
        axes[1, 1].invert_yaxis()
        axes[1, 1].set_xlabel("-log10(Adjusted P-value)", fontsize=10)
        axes[1, 1].set_title(
            "D. Top Enriched Module Pathways", fontsize=12, fontweight="bold"
        )
    else:
        axes[1, 1].text(
            0.5, 0.5, "No Pathway Enrichment Data", ha="center", va="center"
        )
        axes[1, 1].set_title(
            "D. Top Enriched Module Pathways", fontsize=12, fontweight="bold"
        )

    fig.suptitle(
        "GeneCircuitry Comparative Analysis Dashboard",
        fontsize=16,
        fontweight="bold",
        y=0.99,
    )
    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="comparative",
        metadata={"plot_name": "summary_dashboard"},
        skip_existing=False,
    )


def plot_module_overlap_heatmap(
    coverage_df: pd.DataFrame,
    jaccard_df: pd.DataFrame,
    save_name: str = 'default',
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """Plot module overlap heatmaps showing coverage and Jaccard similarity."""
    suffix = _clean_save_suffix(save_name)
    filepath = f'{config.FIGURES_DIR_COMPARATIVE}/module_overlap_heatmap{suffix}.png'
    
    if plot_exists(filepath, skip_existing):
        return False
        
    cov_empty = coverage_df is None or coverage_df.empty
    jac_empty = jaccard_df is None or jaccard_df.empty
    
    if cov_empty and jac_empty:
        return False
        
    n_panels = 2 if not cov_empty and not jac_empty else 1
    if figsize is None:
        figsize = (16, 8) if n_panels == 2 else (8, 8)
            
    fig, axes = plt.subplots(1, n_panels, figsize=figsize)
    if n_panels == 1:
        axes = [axes]
        
    idx = 0
    if not cov_empty:
        sns.heatmap(coverage_df, cmap='YlGn', annot=False, fmt='.2f', ax=axes[idx])
        axes[idx].set_title('Module Gene GRN Coverage per Cluster')
        idx += 1
        
    if not jac_empty:
        sns.heatmap(jaccard_df, cmap='PuBu', annot=False, fmt='.2f', ax=axes[idx])
        axes[idx].set_title('Module Gene Overlap Between Clusters (Jaccard)')
        
    plt.tight_layout()
    return save_plot(
        fig=fig, filepath=filepath, plot_type="comparative",
        metadata={"plot_name": "module_overlap_heatmap"}, skip_existing=False
    )


def plot_module_tf_regulatory_network(
    integration_df: pd.DataFrame,
    save_name: str = 'default',
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
    top_n_tfs: int = 15,
    top_n_edges: int = 50,
) -> bool:
    """Bipartite network graph for TFs and Modules."""
    if not HAS_NETWORKX:
        return False
        
    suffix = _clean_save_suffix(save_name)
    filepath = f'{config.FIGURES_DIR_COMPARATIVE}/module_tf_regulatory_network{suffix}.png'
    
    if plot_exists(filepath, skip_existing):
        return False
        
    if integration_df is None or integration_df.empty:
        return False
        
    tf_totals = integration_df.groupby('tf', observed=False)['n_targets_in_module'].sum()
    top_tfs = tf_totals.nlargest(top_n_tfs).index.tolist()
    
    df_filtered = integration_df[integration_df['tf'].isin(top_tfs)].copy()
    df_filtered = df_filtered.nlargest(top_n_edges, 'n_targets_in_module')
    
    if df_filtered.empty:
        return False
        
    G = nx.Graph()
    for _, row in df_filtered.iterrows():
        tf_node = f"TF: {row['tf']}"
        mod_node = f"Mod: {row['module']}"
        if not G.has_node(tf_node):
            G.add_node(tf_node, bipartite=0, type='tf', size=tf_totals.get(row['tf'], 0))
        if not G.has_node(mod_node):
            G.add_node(mod_node, bipartite=1, type='module', size=row['total_module_genes'])
        G.add_edge(tf_node, mod_node, weight=row['n_targets_in_module'], cluster=row['cluster'])
        
    if figsize is None:
        figsize = (12, 10)
        
    fig, ax = plt.subplots(figsize=figsize)
    tf_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'tf']
    mod_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'module']
    
    pos = {}
    for i, node in enumerate(sorted(tf_nodes)):
        pos[node] = (0, i / max(1, len(tf_nodes)-1))
    for i, node in enumerate(sorted(mod_nodes)):
        pos[node] = (1, i / max(1, len(mod_nodes)-1))
        
    node_sizes = [d['size'] * 10 for n, d in G.nodes(data=True)]
    node_colors = ['#4C72B0' if d['type'] == 'tf' else '#DD8452' for n, d in G.nodes(data=True)]
    
    clusters = df_filtered['cluster'].unique()
    cluster_cmap = plt.get_cmap('tab10')
    cluster_color_map = {c: cluster_cmap(i % 10) for i, c in enumerate(clusters)}
    
    edge_colors = [cluster_color_map[G[u][v]['cluster']] for u, v in G.edges()]
    edge_widths = [G[u][v]['weight'] / 5.0 for u, v in G.edges()]
    
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, ax=ax)
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color=edge_colors, alpha=0.6, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=9, ax=ax)
    
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color=c, lw=2, label=str(cl)) for cl, c in cluster_color_map.items()]
    ax.legend(handles=legend_elements, title='Clusters', loc='best')
    ax.set_title('TF → Module Regulatory Network')
    ax.axis('off')
    
    plt.tight_layout()
    return save_plot(
        fig=fig, filepath=filepath, plot_type="comparative",
        metadata={"plot_name": "module_tf_network"}, skip_existing=False
    )


def _draw_alluvial_band(ax, x0, y0_bot, y0_top, x1, y1_bot, y1_top, color, alpha=0.4):
    """Draw a curved band between two stages."""
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    mid = (x0 + x1) / 2
    verts = [
        (x0, y0_bot), (mid, y0_bot), (mid, y1_bot), (x1, y1_bot),
        (x1, y1_top), (mid, y1_top), (mid, y0_top), (x0, y0_top), (x0, y0_bot)
    ]
    codes = [
        Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY
    ]
    patch = PathPatch(Path(verts, codes), facecolor=color, alpha=alpha, edgecolor='none', linewidth=0)
    ax.add_patch(patch)


def plot_gene_selection_sankey(
    provenance_df: pd.DataFrame,
    save_name: str = 'default',
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """Plot alluvial/flow diagram showing gene journey through pipeline stages."""
    suffix = _clean_save_suffix(save_name)
    filepath = f'{config.FIGURES_DIR_COMPARATIVE}/gene_selection_sankey{suffix}.png'
    
    if plot_exists(filepath, skip_existing):
        return False
        
    if provenance_df is None or provenance_df.empty:
        return False
        
    if figsize is None:
        figsize = (14, 8)
        
    fig, ax = plt.subplots(figsize=figsize)
    stages = ['Hotspot Significant', 'Module Assignment', 'GRN Role', 'Pathway']
    
    def get_groups(col):
        if col == 'all':
            return ['All Genes'] * len(provenance_df)
        val_counts = provenance_df[col].value_counts()
        if len(val_counts) > 10:
            top_cats = set(val_counts.head(8).index)
            return provenance_df[col].apply(lambda x: x if x in top_cats else 'Other').tolist()
        return provenance_df[col].tolist()
        
    provenance_df = provenance_df.copy()
    provenance_df['_stage0'] = get_groups('all')
    provenance_df['_stage1'] = get_groups('hotspot_module')
    provenance_df['_stage2'] = get_groups('stage')
    provenance_df['_stage3'] = get_groups('pathway')
    
    cols = ['_stage0', '_stage1', '_stage2', '_stage3']
    node_positions = {}
    colors = {}
    spacing = len(provenance_df) * 0.05
    
    for i, col in enumerate(cols):
        cats = provenance_df[col].value_counts()
        y = 0
        cat_colors = sns.color_palette('husl', len(cats))
        for j, (cat, count) in enumerate(cats.items()):
            y_top = y + count
            node_positions[(i, cat)] = (y, y_top)
            colors[(i, cat)] = cat_colors[j]
            ax.add_patch(plt.Rectangle((i - 0.05, y), 0.1, count, facecolor=colors[(i, cat)], edgecolor='black'))
            ax.text(i, y + count/2, f"{cat}\n(n={count})", ha='center', va='center', fontsize=8, color='black')
            y = y_top + spacing
            
    for i in range(len(cols) - 1):
        col1 = cols[i]
        col2 = cols[i+1]
        flows = provenance_df.groupby([col1, col2], observed=False).size()
        y_out = {cat: node_positions[(i, cat)][0] for cat in provenance_df[col1].unique()}
        y_in = {cat: node_positions[(i+1, cat)][0] for cat in provenance_df[col2].unique()}
        
        for (c1, c2), count in flows.items():
            if count == 0: continue
            y0_bot = y_out[c1]
            y0_top = y0_bot + count
            y_out[c1] = y0_top
            
            y1_bot = y_in[c2]
            y1_top = y1_bot + count
            y_in[c2] = y1_top
            
            _draw_alluvial_band(ax, i + 0.05, y0_bot, y0_top, i + 1 - 0.05, y1_bot, y1_top, color=colors[(i, c1)], alpha=0.4)
            
    ax.set_xticks(range(4))
    ax.set_xticklabels(stages, fontweight='bold')
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    ax.set_title('Gene Selection Provenance: Hotspot → Modules → GRN → Pathways')
    plt.tight_layout()
    return save_plot(
        fig=fig, filepath=filepath, plot_type="comparative",
        metadata={"plot_name": "gene_selection_sankey"}, skip_existing=False
    )


def plot_cross_cluster_regulatory_comparison(
    reg_summary_df: pd.DataFrame,
    save_name: str = 'default',
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """Multi-panel comparison showing regulatory statistics per cluster."""
    suffix = _clean_save_suffix(save_name)
    filepath = f'{config.FIGURES_DIR_COMPARATIVE}/cross_cluster_regulatory_comparison{suffix}.png'
    
    if plot_exists(filepath, skip_existing):
        return False
        
    if reg_summary_df is None or reg_summary_df.empty:
        return False
        
    if figsize is None:
        figsize = (16, 12)
        
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    clusters = reg_summary_df['cluster'].tolist()
    x = np.arange(len(clusters))
    width = 0.35
    
    axes[0, 0].bar(x - width/2, reg_summary_df['n_regulatory_edges'], width, label='Edges')
    axes[0, 0].bar(x + width/2, reg_summary_df['n_unique_targets'], width, label='Unique Targets')
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(clusters)
    axes[0, 0].legend()
    axes[0, 0].set_title('Regulatory Edges and Targets')
    
    axes[0, 1].bar(x - width/2, reg_summary_df['n_active_tfs'], width, label='Active TFs')
    axes[0, 1].bar(x + width/2, reg_summary_df['n_active_modules'], width, label='Active Modules')
    axes[0, 1].set_xticks(x)
    axes[0, 1].set_xticklabels(clusters)
    axes[0, 1].legend()
    axes[0, 1].set_title('Active TFs and Modules')
    
    axes[1, 0].axis('off')
    cell_text = []
    for tf_val in reg_summary_df['top_tfs']:
        if isinstance(tf_val, (list, tuple)):
            cell_text.append([", ".join(map(str, tf_val[:5]))])
        elif isinstance(tf_val, str) and tf_val.startswith('['):
            try:
                tfs = eval(tf_val)
                cell_text.append([", ".join(map(str, tfs[:5]))])
            except:
                cell_text.append([str(tf_val)[:50]])
        else:
            cell_text.append([str(tf_val)[:50]])
    axes[1, 0].table(cellText=cell_text, rowLabels=clusters, colLabels=['Top TFs'], loc='center')
    axes[1, 0].set_title('Top TFs per Cluster')
    
    axes[1, 1].axis('off')
    cell_text_edges = []
    for edge_val in reg_summary_df['top_edges']:
        if isinstance(edge_val, (list, tuple)):
            cell_text_edges.append([", ".join(map(str, edge_val[:3]))])
        elif isinstance(edge_val, str) and edge_val.startswith('['):
            try:
                edges = eval(edge_val)
                cell_text_edges.append([", ".join(map(str, edges[:3]))])
            except:
                cell_text_edges.append([str(edge_val)[:50]])
        else:
            cell_text_edges.append([str(edge_val)[:50]])
    axes[1, 1].table(cellText=cell_text_edges, rowLabels=clusters, colLabels=['Top Edges'], loc='center')
    axes[1, 1].set_title('Top Edges per Cluster')
    
    fig.suptitle('Cross-Cluster Regulatory Comparison')
    plt.tight_layout()
    return save_plot(
        fig=fig, filepath=filepath, plot_type="comparative",
        metadata={"plot_name": "cross_cluster_regulatory_comparison"}, skip_existing=False
    )


def plot_integrated_regulatory_dashboard(
    comparative_results: Dict[str, Any],
    save_name: str = 'default',
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """Larger, richer dashboard replacing the existing minimal 4-panel one. 3x2 = 6 panels."""
    suffix = _clean_save_suffix(save_name)
    filepath = f'{config.FIGURES_DIR_COMPARATIVE}/integrated_regulatory_dashboard{suffix}.png'
    
    if plot_exists(filepath, skip_existing):
        return False
        
    if not comparative_results:
        return False
        
    if figsize is None:
        figsize = (22, 24)
        
    fig, axes = plt.subplots(3, 2, figsize=figsize)
    
    def set_no_data(ax, title):
        ax.text(0.5, 0.5, "No Data", ha="center", va="center")
        ax.set_title(title, fontsize=12, fontweight="bold")
        
    activity_df = comparative_results.get('module_activity')
    if activity_df is not None and not activity_df.empty:
        sns.heatmap(activity_df, cmap=config.PLOT_COLOR_PALETTE, annot=False, fmt=".2f", ax=axes[0, 0])
        axes[0, 0].set_title("Module Activity", fontsize=12, fontweight="bold")
    else:
        set_no_data(axes[0, 0], "Module Activity")
        
    coverage_df = comparative_results.get('module_coverage')
    if coverage_df is not None and not coverage_df.empty:
        sns.heatmap(coverage_df, cmap='YlGn', annot=False, fmt=".2f", ax=axes[0, 1])
        axes[0, 1].set_title("Module Coverage", fontsize=12, fontweight="bold")
    else:
        set_no_data(axes[0, 1], "Module Coverage")
        
    tf_centrality = comparative_results.get('tf_centrality')
    if tf_centrality is not None and not tf_centrality.empty:
        sns.heatmap(tf_centrality.head(10), cmap='Blues', annot=False, fmt=".2f", ax=axes[1, 0])
        axes[1, 0].set_title("TF Centrality (Top 10)", fontsize=12, fontweight="bold")
    else:
        set_no_data(axes[1, 0], "TF Centrality (Top 10)")
        
    provenance_df = comparative_results.get('gene_provenance')
    if provenance_df is not None and not provenance_df.empty and 'stage' in provenance_df:
        stage_counts = provenance_df['stage'].value_counts()
        stage_counts.plot(kind='barh', ax=axes[1, 1], color=sns.color_palette('husl', len(stage_counts)))
        axes[1, 1].set_title("Gene Provenance Summary", fontsize=12, fontweight="bold")
        axes[1, 1].set_xlabel("Gene Count")
    else:
        set_no_data(axes[1, 1], "Gene Provenance Summary")
        
    reg_summary_df = comparative_results.get('regulatory_summary')
    if reg_summary_df is not None and not reg_summary_df.empty:
        x = np.arange(len(reg_summary_df))
        width = 0.35
        axes[2, 0].bar(x - width/2, reg_summary_df['n_regulatory_edges'], width, label='Edges')
        axes[2, 0].bar(x + width/2, reg_summary_df['n_unique_targets'], width, label='Targets')
        axes[2, 0].set_xticks(x)
        axes[2, 0].set_xticklabels(reg_summary_df['cluster'].tolist())
        axes[2, 0].legend()
        axes[2, 0].set_title("Regulatory Summary", fontsize=12, fontweight="bold")
    else:
        set_no_data(axes[2, 0], "Regulatory Summary")
        
    integration_df = comparative_results.get('module_tf_integration')
    if integration_df is not None and not integration_df.empty:
        top_mods = integration_df['module'].value_counts().head(3).index
        mod_df = integration_df[integration_df['module'].isin(top_mods)]
        if not mod_df.empty:
            sns.barplot(data=mod_df, x='module', y='n_targets_in_module', hue='tf', ax=axes[2, 1])
            axes[2, 1].set_title("Module-TF Integration", fontsize=12, fontweight="bold")
            axes[2, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            set_no_data(axes[2, 1], "Module-TF Integration")
    else:
        set_no_data(axes[2, 1], "Module-TF Integration")
        
    fig.suptitle('GeneCircuitry Integrated Regulatory Dashboard', fontsize=16, fontweight='bold')
    plt.tight_layout()
    return save_plot(
        fig=fig, filepath=filepath, plot_type="comparative",
        metadata={"plot_name": "integrated_regulatory_dashboard"}, skip_existing=False
    )


def generate_all_comparative_plots(
    comparative_results: Dict[str, Any],
    save_name: str = "default",
    skip_existing: bool = True,
) -> Dict[str, bool]:
    """
    Generate all comparative visualizations and return status dictionary.

    Parameters
    ----------
    comparative_results : dict
        Dictionary of comparative analysis results from `run_comparative_analysis()`.
    save_name : str, default="default"
        Suffix for saved filenames.
    skip_existing : bool, default=True
        Whether to skip already existing plots.

    Returns
    -------
    dict
        Mapping of plot names to generation status (True = created, False = skipped/failed).
    """
    results: Dict[str, bool] = {}

    print("Generating Comparative visualizations...")

    # 1. Module Activity Heatmap
    activity_df = comparative_results.get("module_activity")
    if activity_df is not None and not activity_df.empty:
        results["module_activity"] = plot_comparative_module_activity(
            activity_df, save_name=save_name, skip_existing=skip_existing
        )

    # 2. Pathway Enrichment Summary
    enrichment_df = comparative_results.get("module_enrichment")
    if enrichment_df is not None and not enrichment_df.empty:
        results["pathway_enrichment"] = plot_comparative_pathway_enrichment(
            enrichment_df, save_name=save_name, skip_existing=skip_existing
        )

    # 3. TF-to-Module Regulatory Matrix
    tf_mod_matrix = comparative_results.get("tf_to_module_matrix")
    if tf_mod_matrix is not None and not tf_mod_matrix.empty:
        results["tf_module_matrix"] = plot_tf_module_regulatory_matrix(
            tf_mod_matrix, save_name=save_name, skip_existing=skip_existing
        )

    # 4. TF Centrality Matrix
    tf_pivot_df = comparative_results.get("tf_centrality")
    tf_summary_df = comparative_results.get("tf_summary")
    if tf_pivot_df is not None and not tf_pivot_df.empty:
        results["tf_centrality"] = plot_comparative_tf_centrality(
            tf_pivot_df,
            tf_summary_df=tf_summary_df,
            save_name=save_name,
            skip_existing=skip_existing,
        )

    # 5. Differential TF Targets
    diff_targets_df = comparative_results.get("differential_tf_targets")
    if diff_targets_df is not None and not diff_targets_df.empty:
        results["differential_targets"] = plot_differential_tf_targets(
            diff_targets_df, save_name=save_name, skip_existing=skip_existing
        )

    # 6. Integrated Summary Dashboard
    results["summary_dashboard"] = plot_comparative_summary_dashboard(
        comparative_results, save_name=save_name, skip_existing=skip_existing
    )
    
    # 7. Module Overlap Heatmap
    coverage_df = comparative_results.get('module_coverage')
    jaccard_df = comparative_results.get('module_jaccard')
    if (coverage_df is not None and not coverage_df.empty) or (jaccard_df is not None and not jaccard_df.empty):
        results['module_overlap'] = plot_module_overlap_heatmap(
            coverage_df=coverage_df if coverage_df is not None else pd.DataFrame(),
            jaccard_df=jaccard_df if jaccard_df is not None else pd.DataFrame(),
            save_name=save_name, skip_existing=skip_existing,
        )

    # 8. Module-TF Regulatory Network
    integration_df = comparative_results.get('module_tf_integration')
    if integration_df is not None and not integration_df.empty:
        results['module_tf_network'] = plot_module_tf_regulatory_network(
            integration_df, save_name=save_name, skip_existing=skip_existing,
        )

    # 9. Gene Selection Sankey
    provenance_df = comparative_results.get('gene_provenance')
    if provenance_df is not None and not provenance_df.empty:
        results['gene_sankey'] = plot_gene_selection_sankey(
            provenance_df, save_name=save_name, skip_existing=skip_existing,
        )

    # 10. Cross-Cluster Regulatory Comparison
    reg_summary_df = comparative_results.get('regulatory_summary')
    if reg_summary_df is not None and not reg_summary_df.empty:
        results['regulatory_comparison'] = plot_cross_cluster_regulatory_comparison(
            reg_summary_df, save_name=save_name, skip_existing=skip_existing,
        )

    # 11. Integrated Regulatory Dashboard
    results['integrated_dashboard'] = plot_integrated_regulatory_dashboard(
        comparative_results, save_name=save_name, skip_existing=skip_existing,
    )

    generated = sum(1 for status in results.values() if status)
    print(
        f"  Comparative plots: {generated} generated, {len(results) - generated} skipped/unchanged"
    )

    return results
