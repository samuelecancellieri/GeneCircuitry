"""
Hotspot Plotting Module for TRNspot
===================================

Functions for generating Hotspot gene module visualizations including:
- Local correlation heatmaps
- Module annotation plots with enrichment
- Module score violin plots per cluster
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Any, Tuple
from matplotlib.patches import Patch
from anndata import AnnData

from .. import config
from ..logging_utils import log_error, log_warning
from .utils import save_plot, plot_exists

# Import hotspot only when needed to avoid hard dependency
try:
    import hotspot as hs
except ImportError:
    hs = None


def plot_hotspot_local_correlations(
    hotspot_obj,
    skip_existing: bool = True,
) -> bool:
    """
    Plot the local correlation heatmap from Hotspot analysis.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with computed local correlations.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = f"{config.FIGURES_DIR_HOTSPOT}/hotspot_local_correlations.png"

    if plot_exists(filepath, skip_existing):
        return False

    plt.close("all")
    hotspot_obj.plot_local_correlations()
    fig = plt.gcf()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="hotspot",
        metadata={"plot_name": "local_correlations"},
        skip_existing=False,
    )


def _get_module_enrichment_labels(
    hotspot_obj,
    gene_sets: List[str] = ["MSigDB_Hallmark_2020"],
    max_label_length: int = 30,
) -> Dict[int, str]:
    """
    Get enrichment-based labels for each module.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with modules.
    gene_sets : list
        Gene sets for enrichment analysis.
    max_label_length : int
        Maximum length of label text.

    Returns
    -------
    dict
        Mapping from module number to enrichment label.
    """
    # Import enrichment analysis
    try:
        from .. import enrichment_analysis as ea
    except ImportError as e:
        log_warning(
            "HotspotPlotting.ModuleLabels",
            f"enrichment_analysis not available ({type(e).__name__}): {e}",
        )
        return {m: f"Module {m}" for m in hotspot_obj.modules.unique() if m != -1}

    module_labels = {}

    # Try to load existing enrichment results
    enrichment_file = (
        f"{config.OUTPUT_DIR}/hotspot/hotspot_module_enrichment_results.csv"
    )

    if os.path.exists(enrichment_file):
        try:
            df_enrichment = pd.read_csv(enrichment_file)
            for module in hotspot_obj.modules.unique():
                if module == -1:
                    continue
                module_df = df_enrichment[df_enrichment["module"] == module]
                if not module_df.empty:
                    if "Combined_Score" in module_df.columns:
                        top_term = module_df.nlargest(1, "Combined_Score")["Term"].iloc[
                            0
                        ]
                    else:
                        top_term = module_df.nsmallest(1, "Adjusted_P-value")[
                            "Term"
                        ].iloc[0]
                    # Clean up term name
                    top_term = (
                        top_term.replace("HALLMARK_", "").replace("_", " ").title()
                    )
                    if len(top_term) > max_label_length:
                        top_term = top_term[: max_label_length - 3] + "..."
                    module_labels[module] = f"M{module}: {top_term}"
                else:
                    module_labels[module] = f"Module {module}"
            return module_labels
        except Exception as e:
            log_error("HotspotPlotting.LoadEnrichmentFile", e)

    # If no file exists, compute enrichment on the fly
    for module in hotspot_obj.modules.unique():
        if module == -1:
            continue
        genes = hotspot_obj.modules[hotspot_obj.modules == module].index.tolist()
        try:
            enr_result = ea.gseapy_ora_enrichment_analysis(genes, gene_sets=gene_sets)
            if enr_result.results is not None and not enr_result.results.empty:
                if "Combined_Score" in enr_result.results.columns:
                    top_term = enr_result.results.nlargest(1, "Combined_Score")[
                        "Term"
                    ].iloc[0]
                else:
                    top_term = enr_result.results.nsmallest(1, "Adjusted P-value")[
                        "Term"
                    ].iloc[0]
                top_term = top_term.replace("HALLMARK_", "").replace("_", " ").title()
                if len(top_term) > max_label_length:
                    top_term = top_term[: max_label_length - 3] + "..."
                module_labels[module] = f"M{module}: {top_term}"
            else:
                module_labels[module] = f"Module {module}"
        except Exception as e:
            log_warning(
                f"HotspotPlotting.ModuleEnrichment(module={module})",
                f"Enrichment failed ({type(e).__name__}): {e}",
            )
            module_labels[module] = f"Module {module}"

    return module_labels


def plot_hotspot_annotation(
    hotspot_obj,
    gene_sets: List[str] = ["MSigDB_Hallmark_2020"],
    top_n_annotations: int = 1,
    skip_existing: bool = True,
) -> bool:
    """
    Plot Hotspot local correlation heatmap with enrichment annotations.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with analysis results.
    gene_sets : list
        Gene sets for enrichment analysis.
    top_n_annotations : int
        Number of top annotations per module.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = (
        f"{config.FIGURES_DIR_HOTSPOT}/"
        "hotspot_local_correlation_heatmap_with_annotations.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    print("  Generating annotated local correlation heatmap...")

    # Import enrichment analysis
    try:
        from .. import enrichment_analysis as ea
    except ImportError as e:
        log_warning(
            "HotspotPlotting.Annotation",
            f"enrichment_analysis not available ({type(e).__name__}): {e}",
        )
        print("  Warning: Enrichment analysis not available")
        return False

    # Perform enrichment analysis for each module
    df_enrichment = pd.DataFrame()
    for module in hotspot_obj.modules.unique():
        if module == -1:
            continue
        genes = hotspot_obj.modules[hotspot_obj.modules == module].index.tolist()
        try:
            enr_result = ea.gseapy_ora_enrichment_analysis(genes, gene_sets=gene_sets)
            if enr_result.results is not None and not enr_result.results.empty:
                df_module_enrichment = enr_result.results.copy()
                df_module_enrichment.columns = [
                    x.replace(" ", "_") for x in df_module_enrichment.columns
                ]
                df_module_enrichment["module"] = module
                df_enrichment = pd.concat([df_enrichment, df_module_enrichment])
        except Exception as e:
            log_error(f"HotspotPlotting.Enrichment(module={module})", e)
            print(
                f"  Warning: Enrichment failed for module {module} "
                f"({type(e).__name__}): {e}"
            )
            continue

    # Save enrichment results
    if not df_enrichment.empty:
        df_enrichment.to_csv(
            f"{config.OUTPUT_DIR}/hotspot/hotspot_module_enrichment_results.csv",
            index=False,
        )

    # Create module colors using a large, distinct color palette
    # Combine multiple palettes to ensure enough unique colors
    n_modules_total = len([m for m in hotspot_obj.modules.unique() if m != -1])
    if n_modules_total <= 10:
        colors = sns.color_palette("tab10", n_colors=10)
    elif n_modules_total <= 20:
        colors = sns.color_palette("tab20", n_colors=20)
    else:
        # For many modules, use husl which generates evenly spaced hues
        colors = sns.color_palette("husl", n_colors=max(n_modules_total, 30))

    module_colors = {
        i: colors[(i - 1) % len(colors)] for i in hotspot_obj.modules.unique()
    }
    module_colors[-1] = "#ffffff"

    row_colors1 = pd.Series(
        [module_colors[i] for i in hotspot_obj.modules],
        index=hotspot_obj.local_correlation_z.index,
    )

    # Get top annotations for each module
    module_annotations = {}
    for module in hotspot_obj.modules.unique():
        if module == -1:
            module_annotations[module] = ""
            continue
        if not df_enrichment.empty and module in df_enrichment["module"].values:
            module_df = df_enrichment[df_enrichment["module"] == module]
            if not module_df.empty:
                if "Combined_Score" in module_df.columns:
                    top_terms = module_df.nlargest(top_n_annotations, "Combined_Score")[
                        "Term"
                    ].tolist()
                else:
                    top_terms = module_df.nsmallest(
                        top_n_annotations, "Adjusted_P-value"
                    )["Term"].tolist()
                top_terms = [
                    t.replace("HALLMARK_", "").replace("_", " ").title()
                    for t in top_terms
                ]
                module_annotations[module] = "; ".join(top_terms[:top_n_annotations])
            else:
                module_annotations[module] = "No enrichment"
        else:
            module_annotations[module] = "No enrichment"

    # Create the clustermap
    g = sns.clustermap(
        hotspot_obj.local_correlation_z,
        row_linkage=hotspot_obj.linkage,
        col_linkage=hotspot_obj.linkage,
        row_colors=row_colors1,
        cmap="RdBu_r",
        vmin=-8,
        vmax=8,
        xticklabels=False,
        yticklabels=False,
        rasterized=True,
        figsize=(12, 12),
    )

    # Remove dendrograms
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    # Create legend
    legend_elements = []
    for module, color in sorted(module_colors.items()):
        if module == -1:
            continue
        annotation = module_annotations.get(module, "")
        if annotation and annotation != "No enrichment":
            if len(annotation) > 50:
                annotation = annotation[:47] + "..."
            label = f"M{module}: {annotation}"
        else:
            label = f"Module {module}"
        legend_elements.append(Patch(facecolor=color, edgecolor="k", label=label))

    g.ax_heatmap.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        title="Modules (Enrichment)",
        frameon=False,
        fontsize=8,
    )

    fig = g.fig

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="hotspot",
        metadata={
            "plot_name": "annotated_heatmap",
            "gene_sets": gene_sets,
            "top_n_annotations": top_n_annotations,
        },
        skip_existing=False,
    )


def plot_module_scores_violin(
    hotspot_obj,
    adata: AnnData,
    cluster_key: str = "leiden",
    gene_sets: List[str] = ["MSigDB_Hallmark_2020"],
    figsize_per_cluster: Tuple[int, int] = (16, 8),
    skip_existing: bool = True,
) -> Dict[str, bool]:
    """
    Plot violin plots of module scores for each cluster.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with computed module scores.
    adata : AnnData
        AnnData object with cluster annotations.
    cluster_key : str
        Column name in adata.obs containing cluster assignments.
    gene_sets : list
        Gene sets for enrichment annotation labels.
    figsize_per_cluster : tuple
        Figure size for each cluster subplot.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    dict
        Dictionary mapping plot names to generation status.
    """
    results = {}

    # Get module scores
    module_scores = hotspot_obj.module_scores

    if module_scores is None or module_scores.empty:
        log_warning(
            "HotspotPlotting.ViolinPlot", "No module scores available for violin plots"
        )
        print("  Warning: No module scores available for violin plots")
        return results

    # Check cluster key
    if cluster_key not in adata.obs.columns:
        log_warning(
            "HotspotPlotting.ViolinPlot",
            f"Cluster key '{cluster_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}",
        )
        print(f"  Warning: Cluster key '{cluster_key}' not found in adata.obs")
        return results

    # Align cell indices
    common_cells = module_scores.index.intersection(adata.obs.index)
    if len(common_cells) == 0:
        log_warning(
            "HotspotPlotting.ViolinPlot",
            f"No common cells between module scores ({len(module_scores)}) and adata ({adata.n_obs})",
        )
        print("  Warning: No common cells between module scores and adata")
        return results

    # Get module enrichment annotations
    module_labels = _get_module_enrichment_labels(hotspot_obj, gene_sets)

    # Create combined dataframe
    scores_df = module_scores.loc[common_cells].copy()
    scores_df["cluster"] = adata.obs.loc[common_cells, cluster_key].values

    # Get unique modules (excluding -1)
    modules = [col for col in module_scores.columns if col != -1]
    if not modules:
        log_warning(
            "HotspotPlotting.ViolinPlot", "No valid modules found (all modules are -1)"
        )
        print("  Warning: No valid modules found for violin plots")
        return results

    # Create module to label mapping
    module_to_label = {m: module_labels.get(m, f"Module {m}") for m in modules}

    # Rename columns to use enrichment labels
    rename_dict = {m: module_to_label[m] for m in modules}
    scores_df_labeled = scores_df.rename(columns=rename_dict)
    labeled_modules = [module_to_label[m] for m in modules]

    # Melt dataframe
    scores_melted = scores_df_labeled.melt(
        id_vars=["cluster"],
        value_vars=labeled_modules,
        var_name="Module",
        value_name="Score",
    )

    clusters = sorted(scores_df["cluster"].unique())
    n_clusters = len(clusters)
    n_modules = len(modules)

    # Plot 1: Faceted plot by cluster with shared x-axis
    filepath1 = (
        f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_per_cluster.png"
    )

    if not plot_exists(filepath1, skip_existing):
        # Calculate optimal layout - prefer more rows for better x-axis visibility
        n_cols = min(2, n_clusters)  # Max 2 columns for wider violins
        n_rows = (n_clusters + n_cols - 1) // n_cols

        # Larger figure size for better visibility
        fig_width = max(16, n_modules * 1.8) * n_cols / 2
        fig_height = 6 * n_rows

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(fig_width, fig_height),
            squeeze=False,
            sharex=True,  # Share x-axis for consistency
            sharey=False,  # Share y-axis for comparison
        )
        axes = axes.flatten()

        # Use a distinct color palette for modules
        module_palette = sns.color_palette("husl", n_colors=n_modules)

        # Get global y-axis limits for consistency
        y_min = scores_melted["Score"].min()
        y_max = scores_melted["Score"].max()
        y_margin = (y_max - y_min) * 0.1

        for idx, cluster in enumerate(clusters):
            ax = axes[idx]
            cluster_data = scores_melted[scores_melted["cluster"] == cluster]

            if not cluster_data.empty:
                sns.violinplot(
                    data=cluster_data,
                    x="Module",
                    y="Score",
                    palette=module_palette,
                    ax=ax,
                    inner="box",
                    cut=0,
                    linewidth=1.5,
                    width=0.85,  # Wider violins
                    saturation=0.9,
                )

                # Style the plot
                ax.set_title(
                    f"Cluster: {cluster}",
                    fontsize=14,
                    fontweight="bold",
                    pad=10,
                )
                ax.set_ylabel("Module Score", fontsize=11)
                ax.set_ylim(y_min - y_margin, y_max + y_margin)

                # Only show x-label on bottom row
                if idx >= (n_rows - 1) * n_cols:
                    ax.set_xlabel("Module", fontsize=11)
                    ax.tick_params(axis="x", rotation=45, labelsize=9)
                    for label in ax.get_xticklabels():
                        label.set_ha("right")
                else:
                    ax.set_xlabel("")

                # Add grid for readability
                ax.yaxis.grid(True, linestyle="--", alpha=0.3)
                ax.set_axisbelow(True)

                # Style spines
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
            else:
                ax.set_visible(False)

        # Hide unused axes
        for idx in range(n_clusters, len(axes)):
            axes[idx].set_visible(False)

        plt.suptitle(
            "Hotspot Module Scores per Cluster",
            fontsize=16,
            fontweight="bold",
            y=1.01,
        )
        plt.tight_layout()

        results["per_cluster"] = save_plot(
            fig=fig,
            filepath=filepath1,
            plot_type="hotspot",
            metadata={
                "plot_name": "module_scores_violin_per_cluster",
                "cluster_key": cluster_key,
                "n_clusters": n_clusters,
            },
            skip_existing=False,
        )
    else:
        results["per_cluster"] = False

    # Plot 2: All clusters combined - horizontal layout with larger violins
    filepath2 = (
        f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_all_clusters.png"
    )

    if not plot_exists(filepath2, skip_existing):
        # Make figure wider to accommodate all modules and clusters
        fig_width = max(30, n_modules * 3.5)
        fig_height = 10

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # Use a color palette that distinguishes clusters well
        cluster_palette = sns.color_palette("Set2", n_colors=n_clusters)

        sns.violinplot(
            data=scores_melted,
            x="Module",
            y="Score",
            hue="cluster",
            palette=cluster_palette,
            ax=ax,
            inner="box",
            cut=0,
            linewidth=1.2,
            width=0.9,  # Wider violins
            saturation=0.85,
        )

        ax.set_title(
            "Hotspot Module Scores by Cluster",
            fontsize=16,
            fontweight="bold",
            pad=15,
        )
        ax.set_xlabel("Module", fontsize=13)
        ax.set_ylabel("Module Score", fontsize=13)
        ax.tick_params(axis="x", rotation=45, labelsize=10)
        ax.tick_params(axis="y", labelsize=10)
        for label in ax.get_xticklabels():
            label.set_ha("right")

        # Move legend outside and make it more visible
        ax.legend(
            title="Cluster",
            title_fontsize=11,
            fontsize=10,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=True,
            fancybox=True,
            shadow=True,
        )

        # Add grid for readability
        ax.yaxis.grid(True, linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        # Style spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        results["all_clusters"] = save_plot(
            fig=fig,
            filepath=filepath2,
            plot_type="hotspot",
            metadata={
                "plot_name": "module_scores_violin_all_clusters",
                "cluster_key": cluster_key,
            },
            skip_existing=False,
        )
    else:
        results["all_clusters"] = False

    # Plot 3: Horizontal violins for better label readability
    filepath3 = (
        f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_horizontal.png"
    )

    if not plot_exists(filepath3, skip_existing):
        # Horizontal layout - modules on y-axis, better for long labels
        fig_height = max(10, n_modules * 1.5)
        fig_width = max(14, n_clusters * 2.5)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        cluster_palette = sns.color_palette("Set2", n_colors=n_clusters)

        sns.violinplot(
            data=scores_melted,
            y="Module",  # Modules on y-axis
            x="Score",  # Score on x-axis
            hue="cluster",
            palette=cluster_palette,
            ax=ax,
            inner="box",
            cut=0,
            linewidth=1.2,
            width=0.85,
            saturation=0.85,
            orient="h",  # Horizontal orientation
        )

        ax.set_title(
            "Hotspot Module Scores by Cluster (Horizontal)",
            fontsize=16,
            fontweight="bold",
            pad=15,
        )
        ax.set_xlabel("Module Score", fontsize=13)
        ax.set_ylabel("Module", fontsize=13)
        ax.tick_params(axis="both", labelsize=10)

        # Move legend outside
        ax.legend(
            title="Cluster",
            title_fontsize=11,
            fontsize=10,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=True,
            fancybox=True,
            shadow=True,
        )

        # Add grid for readability
        ax.xaxis.grid(True, linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        # Style spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        results["horizontal"] = save_plot(
            fig=fig,
            filepath=filepath3,
            plot_type="hotspot",
            metadata={
                "plot_name": "module_scores_violin_horizontal",
                "cluster_key": cluster_key,
            },
            skip_existing=False,
        )
    else:
        results["horizontal"] = False

    return results


def generate_all_hotspot_plots(
    hotspot_obj,
    adata: Optional[AnnData] = None,
    cluster_key: str = "leiden",
    gene_sets: List[str] = ["MSigDB_Hallmark_2020"],
    skip_existing: bool = True,
) -> Dict[str, Any]:
    """
    Generate all Hotspot plots.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with analysis results.
    adata : AnnData, optional
        AnnData object with cluster annotations (for violin plots).
    cluster_key : str
        Column name for cluster assignments.
    gene_sets : list
        Gene sets for enrichment analysis.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    dict
        Dictionary mapping plot types to generation status.
    """
    results = {}

    print("Generating Hotspot plots...")

    # Local correlations
    results["local_correlations"] = plot_hotspot_local_correlations(
        hotspot_obj, skip_existing=skip_existing
    )

    # Annotated heatmap
    results["annotated_heatmap"] = plot_hotspot_annotation(
        hotspot_obj, gene_sets=gene_sets, skip_existing=skip_existing
    )

    # Violin plots (if adata provided)
    if adata is not None:
        violin_results = plot_module_scores_violin(
            hotspot_obj,
            adata,
            cluster_key=cluster_key,
            gene_sets=gene_sets,
            skip_existing=skip_existing,
        )
        results.update(violin_results)

    generated = sum(1 for v in results.values() if v is True)
    skipped = sum(1 for v in results.values() if v is False)
    print(f"  Hotspot plots: {generated} generated, {skipped} skipped")

    return results
