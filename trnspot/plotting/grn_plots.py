"""
GRN Plotting Module for TRNspot
===============================

Functions for generating Gene Regulatory Network visualizations including:
- Network graphs
- Heatmaps of centrality scores
- Scatter plots comparing scores across clusters
- Rank plots showing score differences
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from adjustText import adjust_text
from itertools import combinations
from typing import Optional, List, Dict, Tuple
from matplotlib.lines import Line2D

from .. import config
from .utils import save_plot, plot_exists


def _scale_to_01(values: pd.Series) -> pd.Series:
    """
    Min-max scale a Series to [0, 1].

    Parameters
    ----------
    values : pd.Series
        Numeric values to scale.

    Returns
    -------
    pd.Series
        Scaled values in [0, 1].
    """
    vmin, vmax = values.min(), values.max()
    if vmax == vmin:
        return pd.Series(0.5, index=values.index)
    return (values - vmin) / (vmax - vmin)


def _enrich_tf_targets(
    tf: str,
    targets: List[str],
    gene_sets: list = config.ENRICHMENT_GENE_SETS,
    pval_cutoff: float = 0.05,
) -> str:
    """
    Run ORA enrichment on targets of a TF and return the top term.

    Parameters
    ----------
    tf : str
        Transcription factor name.
    targets : list of str
        List of target gene names.
    gene_sets : list of str, optional
        Gene set libraries for enrichment. Uses defaults from enrichment module.
    pval_cutoff : float
        Adjusted p-value cutoff for significance.

    Returns
    -------
    str
        Top enriched term, or empty string if no significant result.
    """
    try:
        from .. import enrichment_analysis as ea
    except ImportError as e:
        from ..logging_utils import log_warning

        log_warning(
            "GRNPlotting.Enrichment",
            f"enrichment_analysis not available ({type(e).__name__}): {e}",
        )
        return ""

    if len(targets) < 2:
        return ""

    kwargs = {"gene_list": targets, "pval_cutoff": pval_cutoff}
    if gene_sets is not None:
        kwargs["gene_sets"] = gene_sets

    try:
        enr = ea.gseapy_ora_enrichment_analysis(**kwargs)
        if enr.results is not None and not enr.results.empty:
            if "Combined_Score" in enr.results.columns:
                top_term = enr.results.nlargest(1, "Combined_Score")["Term"].iloc[0]
            else:
                top_term = enr.results.nsmallest(1, "Adjusted P-value")["Term"].iloc[0]
            # Clean up term name
            top_term = top_term.replace("HALLMARK_", "").replace("_", " ").title()
            if len(top_term) > 50:
                top_term = top_term[:47] + "..."
            return top_term
    except Exception as e:
        from ..logging_utils import log_warning

        log_warning(
            "GRNPlotting.Enrichment",
            f"ORA enrichment failed for {len(targets)} targets ({type(e).__name__}): {e}",
        )
    return ""


# =============================================================================
# Network Graph Plots
# =============================================================================


def plot_network_graph_single(
    graph: nx.DiGraph,
    cluster: str,
    stratification: str,
    score: str,
    node_size_factor: float = 1000.0,
    edge_width_factor: float = 1.0,
    skip_existing: bool = True,
) -> bool:
    """
    Plot a network graph for a single score.

    Parameters
    ----------
    graph : nx.DiGraph
        The directed graph to plot.
    cluster : str
        Cluster name for labeling.
    stratification : str
        Stratification name for labeling.
    score : str
        Score being visualized.
    node_size_factor : float
        Factor to scale node sizes.
    edge_width_factor : float
        Factor to scale edge widths.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = (
        f"{config.FIGURES_DIR_GRN}/grn_network_{score}_{stratification}_{cluster}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    fig = plt.figure(figsize=config.PLOT_FIGSIZE_LARGE)

    pos = nx.spring_layout(graph, seed=42)

    # Node sizes based on degree
    degrees = dict(graph.degree())
    node_sizes = [degrees[node] * node_size_factor for node in graph.nodes()]

    # Edge widths based on weight
    edge_weights = nx.get_edge_attributes(graph, "weight")
    edge_widths = [
        edge_weights[edge] * edge_width_factor if edge in edge_weights else 1.0
        for edge in graph.edges()
    ]

    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=node_sizes,
        node_color="skyblue",
        alpha=0.7,
        edgecolors="black",
    )
    nx.draw_networkx_edges(
        graph,
        pos,
        width=edge_widths,
        alpha=0.5,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=10,
    )
    nx.draw_networkx_labels(graph, pos, font_size=8)

    plt.title(f"GRN - Cluster: {cluster}, Stratification: {stratification}")
    plt.axis("off")
    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="grn",
        metadata={
            "plot_name": "network_graph",
            "score": score,
            "cluster": cluster,
            "stratification": stratification,
        },
        skip_existing=False,
    )


def plot_enriched_tf_network(
    score_df: pd.DataFrame,
    links_df: pd.DataFrame,
    score: str = "eigenvector_centrality",
    percentile_threshold: float = 90,
    top_n_links: int = 50,
    gene_sets: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot an enriched TF-target network graph.

    Filters TFs by removing those with degree_out == 0, selects the top 10th
    percentile (above the 90th percentile) of TFs per cluster based on the
    chosen score, retrieves their targets, runs ORA enrichment on the target
    gene lists, and draws a directed network with:
    - TF nodes sized by the 0-1 scaled score
    - TF nodes coloured by cluster
    - Target nodes in light grey
    - Enrichment labels annotated on TF nodes

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame with GRN centrality scores.  Must contain columns:
        ``gene``, ``cluster``, and one or more score columns (e.g.
        ``degree_out``, ``eigenvector_centrality``).
    links_df : pd.DataFrame
        Regulatory links with at least ``source``, ``target``, ``cluster``,
        and ``coef_abs`` columns.
    score : str
        Score column to rank TFs by. Default ``eigenvector_centrality``.
    percentile_threshold : float
        Percentile above which TFs are kept (e.g., 90 keeps top 10%).
    top_n_links : int
        Maximum number of links per cluster to include (by ``coef_abs``).
    gene_sets : list of str, optional
        Gene-set libraries for enrichment. Defaults to ``config.ENRICHMENT_GENE_SETS``.
    skip_existing : bool
        If True, skip if the output file already exists.

    Returns
    -------
    int
        Number of plots generated.
    """
    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)

    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    score_df = score_df.copy()
    if "gene" not in score_df.columns:
        score_df = score_df.reset_index()
        if "index" in score_df.columns:
            score_df = score_df.rename(columns={"index": "gene"})

    score_df["cluster"] = score_df["cluster"].astype(str)
    links_df = links_df.copy()
    links_df["cluster"] = links_df["cluster"].astype(str)

    # Get stratification name for filepath
    stratification = "combined"
    if "stratification" in score_df.columns:
        strats = score_df["stratification"].unique()
        if len(strats) == 1:
            stratification = str(strats[0])

    plots_generated = 0

    for cluster in score_df["cluster"].unique():
        cluster_str = str(cluster)
        filepath = (
            f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
            f"grn_enriched_network_{score}_{stratification}_{cluster_str}.png"
        )
        if plot_exists(filepath, skip_existing):
            continue

        cluster_scores = score_df[score_df["cluster"] == cluster_str].copy()

        if cluster_scores.empty:
            continue

        # --- Step 1: Filter out TFs with degree_out == 0 -----------------------
        if "degree_out" in cluster_scores.columns:
            cluster_scores = cluster_scores[cluster_scores["degree_out"] > 0]
        else:
            # Compute degree_out from links_df if column not present
            out_counts = (
                links_df[links_df["cluster"] == cluster_str]
                .groupby("source")
                .size()
                .rename("_degree_out")
            )
            cluster_scores = cluster_scores.merge(
                out_counts, left_on="gene", right_index=True, how="left"
            )
            cluster_scores["_degree_out"] = cluster_scores["_degree_out"].fillna(0)
            cluster_scores = cluster_scores[cluster_scores["_degree_out"] > 0]
            cluster_scores = cluster_scores.drop(columns=["_degree_out"])

        if cluster_scores.empty or score not in cluster_scores.columns:
            continue

        valid_scores = cluster_scores[score].dropna()
        if valid_scores.empty:
            continue

        # --- Step 2: Select top 10th percentile TFs ----------------------------
        pct_val = np.percentile(valid_scores, percentile_threshold)
        top_tfs_df = cluster_scores[cluster_scores[score] >= pct_val].copy()

        if top_tfs_df.empty:
            continue

        top_tf_names = top_tfs_df["gene"].tolist()

        # --- Step 3: Get links for top TFs and build directed graph ------------
        cluster_links = links_df[
            (links_df["cluster"] == cluster_str)
            & (links_df["source"].isin(top_tf_names))
        ].copy()

        if cluster_links.empty:
            continue

        cluster_links = cluster_links.nlargest(top_n_links, "coef_abs")

        graph = nx.from_pandas_edgelist(
            cluster_links,
            source="source",
            target="target",
            edge_attr="coef_abs",
            create_using=nx.DiGraph(),
        )

        if graph.number_of_nodes() == 0:
            continue

        # --- Step 4: Scale score to [0, 1] for node sizing --------------------
        tf_score_map = top_tfs_df.set_index("gene")[score]
        # Check if score values are already in [0, 1]
        if tf_score_map.min() < 0 or tf_score_map.max() > 1:
            tf_score_scaled = _scale_to_01(tf_score_map)
        else:
            tf_score_scaled = tf_score_map

        # --- Step 5: Enrich targets per TF ------------------------------------
        tf_enrichment_labels: Dict[str, str] = {}
        for tf in top_tf_names:
            targets = cluster_links[cluster_links["source"] == tf]["target"].tolist()
            label = _enrich_tf_targets(tf, targets, gene_sets=gene_sets)
            tf_enrichment_labels[tf] = label

        # --- Step 6: Draw network ---------------------------------------------
        fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_LARGE)

        pos = nx.spring_layout(graph, seed=config.RANDOM_SEED, k=1.5)

        # Separate TF and target nodes
        tf_nodes = [n for n in graph.nodes() if n in top_tf_names]
        target_nodes = [n for n in graph.nodes() if n not in top_tf_names]

        # TF node sizes from scaled score (min 300, max 3000)
        tf_sizes = []
        for n in tf_nodes:
            s = tf_score_scaled.get(n, 0.5)
            tf_sizes.append(300 + s * 2700)

        # Draw target nodes (small, grey)
        if target_nodes:
            nx.draw_networkx_nodes(
                graph,
                pos,
                nodelist=target_nodes,
                node_size=80,
                node_color="#d9d9d9",
                alpha=0.6,
                edgecolors="grey",
                linewidths=0.5,
                ax=ax,
            )

        # Draw TF nodes (coloured, sized by score)
        nx.draw_networkx_nodes(
            graph,
            pos,
            nodelist=tf_nodes,
            node_size=tf_sizes,
            node_color=[tf_score_scaled.get(n, 0.5) for n in tf_nodes],
            cmap=plt.cm.YlOrRd,
            vmin=0,
            vmax=1,
            alpha=0.9,
            edgecolors="black",
            linewidths=1.2,
            ax=ax,
        )

        # Colour bar for score
        sm = plt.cm.ScalarMappable(
            cmap=plt.cm.YlOrRd, norm=plt.Normalize(vmin=0, vmax=1)
        )
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label(f"{score} (scaled 0-1)", fontsize=10)

        # Draw edges
        edge_weights = nx.get_edge_attributes(graph, "coef_abs")
        if edge_weights:
            max_w = max(edge_weights.values()) if edge_weights else 1.0
            widths = [
                0.5 + 2.5 * (edge_weights.get(e, 0) / max_w) for e in graph.edges()
            ]
        else:
            widths = [1.0] * graph.number_of_edges()

        nx.draw_networkx_edges(
            graph,
            pos,
            width=widths,
            alpha=0.4,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=12,
            edge_color="#888888",
            ax=ax,
        )

        # Labels for TF nodes (bold, with enrichment annotation)
        tf_labels = {}
        for tf in tf_nodes:
            enr = tf_enrichment_labels.get(tf, "")
            if enr:
                tf_labels[tf] = f"{tf}\n({enr})"
            else:
                tf_labels[tf] = tf

        nx.draw_networkx_labels(
            graph,
            pos,
            labels=tf_labels,
            font_size=7,
            font_weight="bold",
            ax=ax,
        )

        # Labels for target nodes (smaller)
        if target_nodes:
            target_labels = {n: n for n in target_nodes}
            nx.draw_networkx_labels(
                graph,
                pos,
                labels=target_labels,
                font_size=5,
                font_color="#555555",
                ax=ax,
            )

        # Legend
        legend_elements = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="#d9d9d9",
                markersize=8,
                markeredgecolor="grey",
                label="Target gene",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="#e6550d",
                markersize=12,
                markeredgecolor="black",
                label=f"TF (top {100 - percentile_threshold:.0f}%)",
            ),
        ]
        ax.legend(
            handles=legend_elements,
            loc="upper left",
            fontsize=9,
            framealpha=0.8,
        )

        ax.set_title(
            f"Enriched TF Network — {cluster_str}\n"
            f"Score: {score} | Stratification: {stratification}",
            fontsize=13,
        )
        ax.axis("off")
        fig.tight_layout()

        if save_plot(
            fig=fig,
            filepath=filepath,
            plot_type="grn",
            metadata={
                "plot_name": "enriched_tf_network",
                "score": score,
                "cluster": cluster_str,
                "stratification": stratification,
                "percentile_threshold": percentile_threshold,
                "n_tfs": len(tf_nodes),
                "n_targets": len(target_nodes),
                "enrichment_labels": tf_enrichment_labels,
            },
            skip_existing=False,
        ):
            plots_generated += 1

    return plots_generated


def plot_network_graph(
    score_df: pd.DataFrame,
    links_df: pd.DataFrame,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot network graphs for multiple scores.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    links_df : pd.DataFrame
        DataFrame containing links.
    scores : list, optional
        List of scores to plot. If None, uses default centrality scores.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_all",
            "degree_centrality_all",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    plots_generated = 0
    score_df = score_df.copy()
    score_df["cluster"] = score_df["cluster"].astype(str)

    for score in scores:
        top_genes_by_cluster = pd.DataFrame()
        filtered_links_df = pd.DataFrame()

        for cluster in score_df["cluster"].unique():
            cluster_str = str(cluster)
            cluster_scores = score_df.query(f"cluster == '{cluster_str}'").copy()

            if cluster_scores.empty or cluster_scores[score].dropna().empty:
                continue

            percentile = np.percentile(cluster_scores[score].dropna(), 90)
            cluster_scores = cluster_scores.query(f"{score} > {percentile}")

            if cluster_scores.empty:
                continue

            cluster_scores["gene"] = cluster_scores.index
            top_genes_by_cluster = pd.concat(
                [top_genes_by_cluster, cluster_scores], axis=0
            )

            filtered_links_cluster = links_df[links_df["cluster"] == cluster_str]
            filtered_links_cluster = filtered_links_cluster[
                filtered_links_cluster["source"].isin(
                    top_genes_by_cluster.query("cluster==@cluster_str")["gene"]
                )
            ]
            filtered_links_cluster["cluster"] = cluster_str
            filtered_links_cluster = filtered_links_cluster.nlargest(20, "coef_abs")
            filtered_links_df = pd.concat(
                [filtered_links_df, filtered_links_cluster], axis=0
            )

        if top_genes_by_cluster.empty or filtered_links_df.empty:
            continue

        graph = nx.from_pandas_edgelist(
            filtered_links_df,
            source="source",
            target="target",
            create_using=nx.Graph(),
        )

        if graph.number_of_nodes() == 0:
            continue

        components = list(nx.connected_components(graph))
        if not components:
            continue

        # Get stratification name
        stratification = "combined"
        if "stratification" in score_df.columns:
            strats = score_df["stratification"].unique()
            if len(strats) == 1:
                stratification = str(strats[0])

        # Plot for each cluster
        for cluster in score_df["cluster"].unique():
            if plot_network_graph_single(
                graph=graph,
                cluster=str(cluster),
                stratification=stratification,
                score=score,
                skip_existing=skip_existing,
            ):
                plots_generated += 1

    return plots_generated


# =============================================================================
# Heatmap Plots
# =============================================================================


def plot_heatmap_single_score(
    heatmap_data: pd.DataFrame,
    cluster1: str,
    cluster2: str,
    score: str,
    top_n_genes: int = 5,
    skip_existing: bool = True,
) -> bool:
    """
    Plot a heatmap for a single score comparing two clusters.

    Parameters
    ----------
    heatmap_data : pd.DataFrame
        DataFrame with melted score data.
    cluster1 : str
        First cluster name.
    cluster2 : str
        Second cluster name.
    score : str
        Score column to plot.
    top_n_genes : int
        Number of top genes to display per stratification.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = (
        f"{config.FIGURES_DIR_GRN}/grn_heatmap_{score}_difference_"
        f"top{top_n_genes}_{cluster1}_vs_{cluster2}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    heatmap_data_score = heatmap_data.query("score == @score").copy()

    if heatmap_data_score.empty:
        return False

    # Transform values (negative for cluster2)
    heatmap_transformed = heatmap_data_score.copy()
    heatmap_transformed.loc[heatmap_transformed["cluster"] == cluster2, "value"] *= -1

    # Aggregate by gene and stratification
    heatmap_final = (
        heatmap_transformed.groupby(["gene", "stratification", "score"])["value"]
        .sum()
        .reset_index()
    )

    # Create pivot table
    pivot_df = heatmap_final.pivot(
        index="gene", columns="stratification", values="value"
    ).fillna(0)

    # Select top genes per stratification
    selected = set()
    for strat in pivot_df.columns:
        selected.update(pivot_df[strat].abs().nlargest(top_n_genes).index)

    pivot_top = pivot_df.loc[list(selected)].reindex(
        pivot_df.loc[list(selected)]
        .abs()
        .sum(axis=1)
        .sort_values(ascending=False)
        .index
    )

    # Create figure
    n_genes = len(pivot_top.index)
    n_strats = len(pivot_top.columns)
    fig_width = max(8, n_strats * 0.8)
    fig_height = max(6, n_genes * 0.3)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        pivot_top,
        annot=False,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        cbar_kws={"label": f"{score} (difference)"},
        ax=ax,
    )

    ax.set_ylabel("Gene")
    ax.set_xlabel("Stratification")
    plt.xticks(rotation=90)

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="grn",
        metadata={
            "plot_name": "heatmap",
            "score": score,
            "cluster1": cluster1,
            "cluster2": cluster2,
            "top_n_genes": top_n_genes,
        },
        skip_existing=False,
    )


def plot_heatmap_scores(
    scores_df: pd.DataFrame,
    top_n_genes: int = 10,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot heatmaps for multiple scores across cluster combinations.

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing scores.
    top_n_genes : int
        Number of top genes to display.
    scores : list, optional
        List of scores to plot. If None, uses default centrality scores.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = scores_df["cluster"].unique()
    cluster_combinations = list(combinations(clusters, 2))

    heatmap_data_melt = scores_df.melt(
        id_vars=["gene", "cluster", "stratification"],
        value_vars=scores,
        var_name="score",
        value_name="value",
    )

    plots_generated = 0

    for cluster1, cluster2 in cluster_combinations:
        for score in scores:
            filtered_data = heatmap_data_melt[
                (heatmap_data_melt["cluster"].isin([cluster1, cluster2]))
                & (heatmap_data_melt["score"] == score)
            ]

            if plot_heatmap_single_score(
                heatmap_data=filtered_data,
                cluster1=str(cluster1),
                cluster2=str(cluster2),
                score=score,
                top_n_genes=top_n_genes,
                skip_existing=skip_existing,
            ):
                plots_generated += 1

    return plots_generated


# =============================================================================
# Scatter Plots
# =============================================================================


def plot_score_comparison_2d(
    score_df: pd.DataFrame,
    value: str,
    cluster1: str,
    cluster2: str,
    percentile: float = 99,
    dot_color: str = "#0096FF",
    edge_color: Optional[str] = "black",
    fillna_with_zero: bool = True,
) -> Tuple[np.ndarray, plt.Figure, List]:
    """
    Create a 2D scatter plot comparing scores between two clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame with scores.
    value : str
        Score column to compare.
    cluster1 : str
        First cluster (x-axis).
    cluster2 : str
        Second cluster (y-axis).
    percentile : float
        Percentile threshold for annotation.
    dot_color : str
        Color for scatter dots.
    edge_color : str, optional
        Edge color for dots.
    fillna_with_zero : bool
        If True, fill NaN with 0. Otherwise, fill with mean.

    Returns
    -------
    tuple
        (genes_of_interest, figure, text_annotations)
    """
    res = score_df[score_df.cluster.isin([cluster1, cluster2])][[value, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=value, columns="cluster", index="gene")

    if fillna_with_zero:
        piv = piv.fillna(0)
    else:
        piv = piv.fillna(piv.mean(axis=0))

    goi1 = piv[piv[cluster1] > np.percentile(piv[cluster1].values, percentile)].index
    goi2 = piv[piv[cluster2] > np.percentile(piv[cluster2].values, percentile)].index
    gois = np.union1d(goi1, goi2)

    x, y = piv[cluster1], piv[cluster2]

    fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_SQUARED)
    sns.scatterplot(
        x=x, y=y, markers="o", s=50, edgecolor=edge_color, color=dot_color, ax=ax
    )
    ax.set_title(f"{value}")

    # Add diagonal line
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),
        np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    ax.plot(lims, lims, "r--", alpha=0.75, zorder=0, linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Add annotations
    texts = []
    for goi in gois:
        x_val, y_val = piv.loc[goi, cluster1], piv.loc[goi, cluster2]
        texts.append(ax.text(x_val, y_val, goi))

    adjust_text(
        texts=texts,
        arrowprops=dict(arrowstyle="->", color="r", lw=0.5),
        ax=ax,
        time_lim=2,
    )

    return gois, fig, texts


def plot_scatter_scores(
    score_df: pd.DataFrame,
    scores_list: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot scatter plots comparing scores between cluster pairs.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    scores_list : list, optional
        List of scores to plot. If None, uses default centrality scores.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores_list is None:
        scores_list = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = sorted(set(score_df["cluster"].tolist()))
    cluster_combinations = list(combinations(clusters, 2))

    # Ensure output directory exists
    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    plots_generated = 0

    for cluster1, cluster2 in cluster_combinations:
        for score in scores_list:
            cluster1_clean = str(cluster1).replace(" ", "_").strip()
            cluster2_clean = str(cluster2).replace(" ", "_").strip()
            filepath = (
                f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
                f"grn_scatter_{score}_{cluster1_clean}_vs_{cluster2_clean}.png"
            )

            if plot_exists(filepath, skip_existing):
                continue

            genes, fig, texts = plot_score_comparison_2d(
                score_df=score_df,
                value=score,
                cluster1=cluster1,
                cluster2=cluster2,
                percentile=99,
                edge_color="black",
                dot_color="#0096FF",
                fillna_with_zero=True,
            )

            fig.axes[0].grid(False)
            sns.despine(ax=fig.axes[0])

            if save_plot(
                fig=fig,
                filepath=filepath,
                plot_type="grn",
                metadata={
                    "plot_name": "scatter",
                    "score": score,
                    "cluster1": cluster1,
                    "cluster2": cluster2,
                },
                skip_existing=False,
            ):
                plots_generated += 1

    return plots_generated


# =============================================================================
# Difference and Comparison Plots
# =============================================================================


def plot_difference_cluster_scores(
    score_df: pd.DataFrame,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot rank plots showing score differences between clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    scores : list, optional
        List of scores to plot.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = sorted(set(score_df["cluster"].tolist()))
    stratifications = sorted(set(score_df["stratification"].tolist()))
    stratification_name = (
        stratifications[0] if len(stratifications) == 1 else "combined"
    )

    if len(clusters) < 2:
        print("  Not enough clusters to compare")
        return 0

    cluster_combinations = list(combinations(clusters, 2))
    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    plots_generated = 0

    for cluster_pair in cluster_combinations:
        cluster1, cluster2 = cluster_pair

        cluster_data_1 = score_df[score_df["cluster"] == cluster1]
        cluster_data_2 = score_df[score_df["cluster"] == cluster2]

        tf_overlap = set(cluster_data_1.index.tolist()) & set(
            cluster_data_2.index.tolist()
        )
        common_tfs = list(tf_overlap)

        if not common_tfs:
            continue

        for score in scores:
            filepath = (
                f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
                f"top_genes_difference_{cluster1}-{cluster2}_{score}_rankplot.png"
            )

            if plot_exists(filepath, skip_existing):
                continue

            diff_df = pd.DataFrame(
                {
                    cluster1: cluster_data_1.loc[common_tfs, score],
                    cluster2: cluster_data_2.loc[common_tfs, score],
                }
            )
            diff_df["Difference"] = diff_df[cluster1] - diff_df[cluster2]
            diff_df_sorted = diff_df.sort_values("Difference", ascending=False)

            fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_SQUARED)

            top_x = 5
            top_positive = diff_df_sorted.head(top_x)
            top_negative = diff_df_sorted.tail(top_x)

            # Color scheme
            colors = ["gray"] * len(diff_df_sorted)
            for i in range(top_x):
                colors[i] = "#0096FF"
            for i in range(len(diff_df_sorted) - top_x, len(diff_df_sorted)):
                colors[i] = "#fb3310"

            ranks = range(len(diff_df_sorted))
            ax.scatter(ranks, diff_df_sorted["Difference"], c=colors, alpha=0.7, s=20)

            # Annotations
            texts = []
            for i, (gene, row) in enumerate(top_positive.iterrows()):
                texts.append(
                    ax.annotate(
                        gene,
                        (i, row["Difference"]),
                        fontsize=8,
                        bbox=dict(
                            boxstyle="round,pad=0.2", facecolor="lightgreen", alpha=0.7
                        ),
                    )
                )

            for i, (gene, row) in enumerate(top_negative.iterrows()):
                rank_pos = len(diff_df_sorted) - top_x + i
                texts.append(
                    ax.annotate(
                        gene,
                        (rank_pos, row["Difference"]),
                        fontsize=8,
                        bbox=dict(
                            boxstyle="round,pad=0.2", facecolor="lightcoral", alpha=0.7
                        ),
                    )
                )

            adjust_text(
                texts, arrowprops=dict(arrowstyle="->", color="black", lw=0.5), ax=ax
            )

            ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
            ax.set_xlabel("Rank")
            ax.set_ylabel("Difference")
            ax.set_title(
                f"Rank Plot: {score}\n{cluster1} vs {cluster2} in {stratification_name}"
            )
            sns.despine(ax=ax)
            ax.grid(False)

            if save_plot(
                fig=fig,
                filepath=filepath,
                plot_type="grn",
                metadata={
                    "plot_name": "difference_rankplot",
                    "score": score,
                    "cluster1": cluster1,
                    "cluster2": cluster2,
                },
                skip_existing=False,
            ):
                plots_generated += 1

    return plots_generated


def plot_compare_cluster_scores(
    score_df: pd.DataFrame,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot bar plots comparing top genes across clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    scores : list, optional
        List of scores to plot.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = sorted(set(score_df["cluster"].tolist()))
    stratifications = sorted(set(score_df["stratification"].tolist()))
    stratification_name = (
        stratifications[0] if len(stratifications) == 1 else "combined"
    )

    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    plots_generated = 0

    for score in scores:
        filepath = (
            f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
            f"grn_barplot_{score}_{stratification_name}_top10.png"
        )

        if plot_exists(filepath, skip_existing):
            continue

        n_clusters = len(clusters)
        n_cols = int(np.ceil(np.sqrt(n_clusters)))
        n_rows = int(np.ceil(n_clusters / n_cols))

        fig, axes = plt.subplots(
            nrows=n_rows,
            ncols=n_cols,
            figsize=(7 * n_cols, 7 * n_rows),
            squeeze=False,
        )
        axes = axes.flatten()

        for idx, cluster in enumerate(clusters):
            cluster_data = score_df[score_df["cluster"] == cluster]
            top10_genes = cluster_data.nlargest(10, score)

            ax = axes[idx]
            sns.scatterplot(
                data=top10_genes,
                x=score,
                y="gene",
                ax=ax,
                s=100,
                edgecolor="black",
                palette="viridis",
            )
            sns.despine(ax=ax)
            ax.grid(False)
            ax.set_title(f"{cluster}")
            ax.set_xlabel(f"{score}")
            ax.set_ylabel("Gene")

        # Hide unused axes
        for idx in range(len(clusters), len(axes)):
            axes[idx].set_visible(False)

        fig.subplots_adjust(wspace=0.2, top=0.9)
        fig.suptitle(
            f"Top 10 Genes by {score} Across Clusters\n{stratification_name}",
            fontsize=16,
        )

        if save_plot(
            fig=fig,
            filepath=filepath,
            plot_type="grn",
            metadata={
                "plot_name": "barplot_comparison",
                "score": score,
                "stratification": stratification_name,
            },
            skip_existing=False,
        ):
            plots_generated += 1

    return plots_generated


# =============================================================================
# Main Entry Point
# =============================================================================


def generate_all_grn_plots(
    score_df: pd.DataFrame,
    links_df: Optional[pd.DataFrame] = None,
    skip_existing: bool = True,
    scores: Optional[List[str]] = None,
    gene_sets: Optional[List[str]] = None,
) -> Dict[str, int]:
    """
    Generate all GRN plots from score data.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing GRN scores.
    links_df : pd.DataFrame, optional
        DataFrame containing GRN links (for network plots).
    skip_existing : bool
        If True, skip existing plots.
    scores : list, optional
        List of scores to plot.
    gene_sets : list of str, optional
        Gene-set libraries for enrichment. Defaults to ``config.ENRICHMENT_GENE_SETS``.

    Returns
    -------
    dict
        Dictionary mapping plot types to counts generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    # Ensure gene column exists
    if "gene" not in score_df.columns:
        score_df = score_df.reset_index()
        if "index" in score_df.columns:
            score_df = score_df.rename(columns={"index": "gene"})

    results = {}

    print("Generating GRN plots...")

    # Scatter plots
    results["scatter"] = plot_scatter_scores(
        score_df, scores_list=scores, skip_existing=skip_existing
    )

    # Difference plots
    results["difference"] = plot_difference_cluster_scores(
        score_df, scores=scores, skip_existing=skip_existing
    )

    # Comparison plots
    results["comparison"] = plot_compare_cluster_scores(
        score_df, scores=scores, skip_existing=skip_existing
    )

    # Heatmap plots
    if "stratification" in score_df.columns:
        results["heatmap"] = plot_heatmap_scores(
            score_df, scores=scores, skip_existing=skip_existing
        )

    # Network plots (if links provided)
    if links_df is not None:
        results["network"] = plot_network_graph(
            score_df, links_df, scores=scores, skip_existing=skip_existing
        )

        # Enriched TF network plots
        enriched_count = 0
        for s in scores:
            enriched_count += plot_enriched_tf_network(
                score_df,
                links_df,
                score=s,
                gene_sets=gene_sets,
                skip_existing=skip_existing,
            )
        results["enriched_network"] = enriched_count

        # TF shared-target enrichment network plots
        shared_count = 0
        for s in scores:
            shared_count += plot_tf_shared_target_network(
                score_df,
                links_df,
                score=s,
                gene_sets=gene_sets,
                skip_existing=skip_existing,
            )
        results["tf_shared_target_network"] = shared_count

    total = sum(results.values())
    print(f"  GRN plots: {total} total generated")
    for plot_type, count in results.items():
        if count > 0:
            print(f"    - {plot_type}: {count}")

    return results


# =============================================================================
# TF Shared-Target Enrichment Network
# =============================================================================


def plot_tf_shared_target_network(
    score_df: pd.DataFrame,
    links_df: pd.DataFrame,
    score: str = "eigenvector_centrality",
    top_n_tfs: int = 20,
    gene_sets: Optional[List[str]] = None,
    min_targets: int = 5,
    skip_existing: bool = True,
) -> int:
    """
    Plot a TF interaction network coloured by target-gene enrichment.

    For each cluster the function:

    1. Selects the top *top_n_tfs* TFs ranked by *score*,
       filtering out TFs with ``degree_out == 0``.
    2. Enriches the targets of every selected TF (ORA, MSigDB
       Hallmark by default) and assigns each TF the top enriched
       term.
    3. Builds a directed graph with **only TF-TF edges** (edges
       where both source and target belong to the top-TF set),
       capturing shared regulatory interactions.
    4. Draws the network using a Kamada-Kawai layout that groups
       TFs with the same enrichment label closer together.
       Node colour = enrichment term, node size = scaled score.
    5. Adds a side-panel table of TF centrality values.

    Parameters
    ----------
    score_df : pd.DataFrame
        GRN centrality scores.  Must contain ``gene``, ``cluster``,
        and one or more score columns.
    links_df : pd.DataFrame
        Regulatory links with ``source``, ``target``, ``cluster``,
        ``coef_abs``.
    score : str
        Column in *score_df* to rank TFs by.
    top_n_tfs : int
        Number of top TFs per cluster.
    gene_sets : list of str, optional
        Gene-set libraries for enrichment.  Defaults to
        ``config.ENRICHMENT_GENE_SETS``.
    min_targets : int
        Minimum number of targets required to run enrichment for
        a TF.
    skip_existing : bool
        Skip plots that already exist on disk.

    Returns
    -------
    int
        Number of plots generated.
    """
    from matplotlib.patches import Patch

    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)

    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    # --- prepare DataFrames --------------------------------------------------
    score_df = score_df.copy()
    if "gene" not in score_df.columns:
        score_df = score_df.reset_index()
        if "index" in score_df.columns:
            score_df = score_df.rename(columns={"index": "gene"})
    score_df["cluster"] = score_df["cluster"].astype(str)

    links_df = links_df.copy()
    links_df["cluster"] = links_df["cluster"].astype(str)

    stratification = "combined"
    if "stratification" in score_df.columns:
        strats = score_df["stratification"].unique()
        if len(strats) == 1:
            stratification = str(strats[0])

    # Persistent colour map across clusters
    palette = sns.color_palette("Set2", 20)
    term_color_map: Dict[str, tuple] = {}
    palette_idx = 0

    plots_generated = 0

    for cluster in score_df["cluster"].unique():
        cluster_str = str(cluster)
        cluster_clean = cluster_str.replace(" ", "_").replace("/", "-")
        filepath = (
            f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
            f"network_tf_enrichment_"
            f"{score}_{stratification}_{cluster_clean}.png"
        )
        if plot_exists(filepath, skip_existing):
            continue

        # --- filter cluster scores -------------------------------------------
        df_cluster = score_df[score_df["cluster"] == cluster_str].copy()
        if df_cluster.empty or score not in df_cluster.columns:
            continue

        # Remove TFs with degree_out == 0
        if "degree_out" in df_cluster.columns:
            df_cluster = df_cluster[df_cluster["degree_out"] > 0]
        else:
            out_counts = (
                links_df[links_df["cluster"] == cluster_str]
                .groupby("source")
                .size()
                .rename("_deg_out")
            )
            df_cluster = df_cluster.merge(
                out_counts,
                left_on="gene",
                right_index=True,
                how="left",
            )
            df_cluster["_deg_out"] = df_cluster["_deg_out"].fillna(0)
            df_cluster = df_cluster[df_cluster["_deg_out"] > 0]
            df_cluster = df_cluster.drop(columns=["_deg_out"])

        if df_cluster.empty:
            continue

        # Top N TFs by score
        df_top = df_cluster.nlargest(top_n_tfs, score)
        top_tfs = df_top["gene"].tolist()
        centrality_map = dict(zip(df_top["gene"], df_top[score]))

        if not top_tfs:
            continue

        # --- enrichment per TF -----------------------------------------------
        links_cluster = links_df[links_df["cluster"] == cluster_str]
        tf_enrichment_map: Dict[str, str] = {}

        for tf in top_tfs:
            targets = links_cluster[links_cluster["source"] == tf]["target"].tolist()
            if len(targets) >= min_targets:
                label = _enrich_tf_targets(tf, targets, gene_sets=gene_sets)
                if label:
                    tf_enrichment_map[tf] = label
                else:
                    tf_enrichment_map[tf] = "No Significant Enrichment"
            else:
                tf_enrichment_map[tf] = f"Too few targets (<{min_targets})"

        # --- build TF-TF graph -----------------------------------------------
        G = nx.DiGraph()
        for tf in top_tfs:
            G.add_node(
                tf,
                enrichment=tf_enrichment_map.get(tf, "Not Analyzed"),
            )

        # Edges only between top TFs (shared-target interactions)
        edges = links_cluster[
            links_cluster["source"].isin(top_tfs)
            & links_cluster["target"].isin(top_tfs)
        ]
        for _, row in edges.iterrows():
            G.add_edge(
                row["source"],
                row["target"],
                weight=row["coef_abs"],
            )

        if G.number_of_nodes() == 0:
            continue

        # --- layout (Kamada-Kawai with grouping) -----------------------------
        nodes = list(G.nodes())
        intra_dist = 0.2
        inter_dist = 0.5
        dist: Dict[str, Dict[str, float]] = {}
        for n1 in nodes:
            dist[n1] = {}
            t1 = G.nodes[n1].get("enrichment", "")
            for n2 in nodes:
                if n1 == n2:
                    dist[n1][n2] = 0.0
                else:
                    t2 = G.nodes[n2].get("enrichment", "")
                    dist[n1][n2] = intra_dist if t1 == t2 else inter_dist

        pos = nx.kamada_kawai_layout(G, dist=dist, scale=1.0)

        # --- colours ---------------------------------------------------------
        unique_terms = sorted(set(tf_enrichment_map.values()))
        no_result_labels = {
            "No Significant Enrichment",
            f"Too few targets (<{min_targets})",
            "Enrichment Error",
            "Not Analyzed",
        }
        for term in unique_terms:
            if term not in term_color_map:
                if term in no_result_labels:
                    term_color_map[term] = (0.8, 0.8, 0.8, 1.0)
                else:
                    term_color_map[term] = palette[palette_idx % len(palette)]
                    palette_idx += 1

        node_colors = [term_color_map[G.nodes[n]["enrichment"]] for n in G.nodes()]

        # --- node sizes (scaled centrality) ----------------------------------
        c_vals = [centrality_map.get(n, 0) for n in G.nodes()]
        c_min, c_max = (min(c_vals), max(c_vals)) if c_vals else (0, 1)
        if c_max > c_min:
            node_sizes = [300 + (c - c_min) / (c_max - c_min) * 1200 for c in c_vals]
        else:
            node_sizes = [600] * len(G.nodes())

        # --- draw figure -----------------------------------------------------
        fig = plt.figure(figsize=(14, 11))
        gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])
        ax_net = fig.add_subplot(gs[0])
        ax_table = fig.add_subplot(gs[1])
        ax_table.axis("off")

        nx.draw_networkx_nodes(
            G,
            pos,
            ax=ax_net,
            node_color=node_colors,
            node_size=node_sizes,
            alpha=0.9,
            edgecolors="white",
        )
        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax_net,
            edge_color="gray",
            alpha=0.4,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=12,
        )
        nx.draw_networkx_labels(
            G,
            pos,
            ax=ax_net,
            font_size=10,
            font_weight="bold",
            font_color="black",
        )

        # Legend
        legend_elements = [
            Patch(
                facecolor=term_color_map[t],
                label=t,
                edgecolor="white",
            )
            for t in unique_terms
        ]
        ax_net.legend(
            handles=legend_elements,
            title="Top Enrichment (Targets)",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0.0,
            fontsize=9,
        )
        ax_net.set_title(
            f"TF Network by Target Enrichment\n"
            f"Cluster: {cluster_str} "
            f"(Top {len(top_tfs)} TFs)\n"
            f"Node size ∝ {score}",
            fontsize=13,
            fontweight="bold",
        )
        ax_net.axis("off")

        # Centrality table
        sorted_tfs = sorted(
            [(tf, centrality_map.get(tf, 0)) for tf in G.nodes()],
            key=lambda x: x[1],
            reverse=True,
        )
        table_data = [[tf, f"{val:.4f}"] for tf, val in sorted_tfs]
        table = ax_table.table(
            cellText=table_data,
            colLabels=["TF", score],
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.2)
        ax_table.set_title(
            score.replace("_", " ").title(),
            fontsize=11,
            fontweight="bold",
        )

        fig.tight_layout()

        saved = save_plot(
            fig=fig,
            filepath=filepath,
            plot_type="grn",
            metadata={
                "plot_name": "tf_shared_target_network",
                "score": score,
                "cluster": cluster_str,
                "stratification": stratification,
                "top_n_tfs": len(top_tfs),
                "enrichment_labels": tf_enrichment_map,
            },
            skip_existing=False,
        )
        if saved:
            plots_generated += 1

        # Save enrichment mapping CSV
        csv_path = (
            f"{config.OUTPUT_DIR}/grn_deep_analysis/"
            f"tf_enrichment_mapping_"
            f"{score}_{stratification}_{cluster_clean}.csv"
        )
        df_mapping = pd.DataFrame(
            list(tf_enrichment_map.items()),
            columns=["TF", "Top_Enrichment_Term"],
        )
        df_mapping["Cluster"] = cluster_str
        df_mapping["Score"] = score
        df_mapping["Centrality"] = df_mapping["TF"].map(centrality_map)
        df_mapping.to_csv(csv_path, index=False)

    return plots_generated
