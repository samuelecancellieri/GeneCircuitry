"""
Hotspot Plotting Module for GeneCircuitry
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
from typing import Optional, List, Dict, Any, Tuple, Union, Sequence
from matplotlib.patches import Patch
from anndata import AnnData

from .. import config
from ..logging_utils import log_error, log_warning
from .utils import save_plot, plot_exists, run_parallel_tasks

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


def _enrich_single_module_worker(task: Dict[str, Any]) -> Tuple[int, Optional[pd.DataFrame]]:
    """Worker function for single module ORA enrichment."""
    module = task["module"]
    genes = task["genes"]
    gene_sets = task.get("gene_sets")

    try:
        from .. import enrichment_analysis as ea

        enr_result = ea.gseapy_ora_enrichment_analysis(genes, gene_sets=gene_sets)
        if enr_result.results is not None and not enr_result.results.empty:
            df_module_enrichment = enr_result.results.copy()
            df_module_enrichment.columns = [
                x.replace(" ", "_") for x in df_module_enrichment.columns
            ]
            df_module_enrichment["module"] = module
            return (module, df_module_enrichment)
    except Exception as e:
        log_warning(
            f"HotspotPlotting.ModuleEnrichment(module={module})",
            f"Enrichment failed ({type(e).__name__}): {e}",
        )

    return (module, None)


def _get_module_enrichment_labels(
    hotspot_obj,
    gene_sets: list = config.ENRICHMENT_GENE_SETS,
    max_label_length: int = 30,
    n_jobs: Optional[int] = None,
) -> Dict[int, str]:
    """
    Get enrichment-based labels for each module.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with modules.
    gene_sets : list
        Gene sets for enrichment analysis. Defaults to ``config.ENRICHMENT_GENE_SETS``.
    max_label_length : int
        Maximum length of label text.
    n_jobs : int, optional
        Number of parallel worker processes. Default uses config.N_JOBS.

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

    # Try to load existing enrichment results
    enrichment_file = (
        f"{config.OUTPUT_DIR}/hotspot/hotspot_module_enrichment_results.csv"
    )

    if os.path.exists(enrichment_file):
        try:
            df_enrichment = pd.read_csv(enrichment_file)
            module_labels = {}
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

    # If no file exists, compute enrichment on the fly in parallel
    unique_modules = [m for m in hotspot_obj.modules.unique() if m != -1]
    tasks = [
        {
            "module": module,
            "genes": hotspot_obj.modules[hotspot_obj.modules == module].index.tolist(),
            "gene_sets": gene_sets,
        }
        for module in unique_modules
    ]

    enrichment_results = run_parallel_tasks(
        _enrich_single_module_worker,
        tasks,
        n_jobs=n_jobs,
        desc="module_enrichment_labels",
    )

    module_labels = {}
    for module, df_res in enrichment_results:
        if df_res is not None and not df_res.empty:
            if "Combined_Score" in df_res.columns:
                top_term = df_res.nlargest(1, "Combined_Score")["Term"].iloc[0]
            else:
                adj_col = "Adjusted_P-value" if "Adjusted_P-value" in df_res.columns else "Adjusted P-value"
                top_term = df_res.nsmallest(1, adj_col)["Term"].iloc[0]
            top_term = top_term.replace("HALLMARK_", "").replace("_", " ").title()
            if len(top_term) > max_label_length:
                top_term = top_term[: max_label_length - 3] + "..."
            module_labels[module] = f"M{module}: {top_term}"
        else:
            module_labels[module] = f"Module {module}"

    return module_labels


def plot_hotspot_annotation(
    hotspot_obj,
    gene_sets: List[str] = None,
    top_n_annotations: int = 1,
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
) -> bool:
    """
    Plot Hotspot local correlation heatmap with enrichment annotations in parallel.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with analysis results.
    gene_sets : list
        Gene sets for enrichment analysis. Defaults to ``config.ENRICHMENT_GENE_SETS``.
    top_n_annotations : int
        Number of top annotations per module.
    skip_existing : bool
        If True, skip if file already exists.
    n_jobs : int, optional
        Number of parallel worker processes. Default uses config.N_JOBS.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)

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

    # Perform enrichment analysis for each module in parallel
    unique_modules = [m for m in hotspot_obj.modules.unique() if m != -1]
    tasks = [
        {
            "module": module,
            "genes": hotspot_obj.modules[hotspot_obj.modules == module].index.tolist(),
            "gene_sets": gene_sets,
        }
        for module in unique_modules
    ]

    enrichment_results = run_parallel_tasks(
        _enrich_single_module_worker,
        tasks,
        n_jobs=n_jobs,
        desc="hotspot_module_enrichment",
    )

    df_list = [df for _, df in enrichment_results if df is not None and not df.empty]
    df_enrichment = pd.concat(df_list, ignore_index=True) if df_list else pd.DataFrame()

    # Save enrichment results
    if not df_enrichment.empty:
        df_enrichment.to_csv(
            f"{config.OUTPUT_DIR}/hotspot/hotspot_module_enrichment_results.csv",
            index=False,
        )

    # Create module colors using a large, distinct color palette
    n_modules_total = len(unique_modules)
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
                    adj_col = (
                        "Adjusted_P-value"
                        if "Adjusted_P-value" in module_df.columns
                        else "Adjusted P-value"
                    )
                    top_terms = module_df.nsmallest(
                        top_n_annotations, adj_col
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


def _plot_hotspot_violin_single_worker(task: Dict[str, Any]) -> Tuple[str, bool]:
    """Worker function for a single hotspot violin plot variant."""
    variant = task["variant"]
    scores_melted = task["scores_melted"]
    clusters = task["clusters"]
    modules = task["modules"]
    cluster_key = task["cluster_key"]
    skip_existing = task.get("skip_existing", True)
    n_clusters = len(clusters)
    n_modules = len(modules)

    if variant == "per_cluster":
        filepath = (
            f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_per_cluster.png"
        )
        if plot_exists(filepath, skip_existing):
            return ("per_cluster", False)

        n_cols = min(2, n_clusters)
        n_rows = (n_clusters + n_cols - 1) // n_cols
        fig_width = max(16, n_modules * 1.8) * n_cols / 2
        fig_height = 6 * n_rows

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(fig_width, fig_height),
            squeeze=False,
            sharex=True,
            sharey=False,
        )
        axes = axes.flatten()

        module_palette = sns.color_palette("husl", n_colors=n_modules)
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
                    hue="Module",
                    palette=module_palette,
                    ax=ax,
                    inner="box",
                    cut=0,
                    linewidth=1.5,
                    width=0.85,
                    saturation=0.9,
                )
                ax.set_title(
                    f"Cluster: {cluster}",
                    fontsize=14,
                    fontweight="bold",
                    pad=10,
                )
                ax.set_ylabel("Module Score", fontsize=11)
                ax.set_ylim(y_min - y_margin, y_max + y_margin)

                if idx >= (n_rows - 1) * n_cols:
                    ax.set_xlabel("Module", fontsize=11)
                    ax.tick_params(axis="x", rotation=45, labelsize=9)
                    for label in ax.get_xticklabels():
                        label.set_ha("right")
                else:
                    ax.set_xlabel("")

                ax.yaxis.grid(True, linestyle="--", alpha=0.3)
                ax.set_axisbelow(True)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)

        for ax in axes:
            if not ax.has_data():
                ax.set_visible(False)

        plt.suptitle(
            "Hotspot Module Scores per Cluster",
            fontsize=16,
            fontweight="bold",
            y=1.01,
        )
        plt.tight_layout()

        saved = save_plot(
            fig=fig,
            filepath=filepath,
            plot_type="hotspot",
            metadata={
                "plot_name": "module_scores_violin_per_cluster",
                "cluster_key": cluster_key,
                "n_clusters": n_clusters,
            },
            skip_existing=False,
        )
        return ("per_cluster", saved)

    elif variant == "all_clusters":
        filepath = (
            f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_all_clusters.png"
        )
        if plot_exists(filepath, skip_existing):
            return ("all_clusters", False)

        fig_width = max(30, n_modules * 3.5)
        fig_height = 10

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
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
            width=0.9,
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

        ax.yaxis.grid(True, linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if not ax.has_data():
            ax.set_visible(False)

        plt.tight_layout()

        saved = save_plot(
            fig=fig,
            filepath=filepath,
            plot_type="hotspot",
            metadata={
                "plot_name": "module_scores_violin_all_clusters",
                "cluster_key": cluster_key,
            },
            skip_existing=False,
        )
        return ("all_clusters", saved)

    elif variant == "horizontal":
        filepath = (
            f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_horizontal.png"
        )
        if plot_exists(filepath, skip_existing):
            return ("horizontal", False)

        fig_height = max(10, n_modules * 1.5)
        fig_width = max(14, n_clusters * 2.5)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        cluster_palette = sns.color_palette("Set2", n_colors=n_clusters)

        sns.violinplot(
            data=scores_melted,
            y="Module",
            x="Score",
            hue="cluster",
            palette=cluster_palette,
            ax=ax,
            inner="box",
            cut=0,
            linewidth=1.2,
            width=0.85,
            saturation=0.85,
            orient="h",
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

        ax.xaxis.grid(True, linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if not ax.has_data():
            ax.set_visible(False)

        plt.tight_layout()

        saved = save_plot(
            fig=fig,
            filepath=filepath,
            plot_type="hotspot",
            metadata={
                "plot_name": "module_scores_violin_horizontal",
                "cluster_key": cluster_key,
            },
            skip_existing=False,
        )
        return ("horizontal", saved)

    return (variant, False)


def plot_module_scores_violin(
    hotspot_obj,
    adata: AnnData,
    cluster_key: Union[str, Sequence[str]] = "leiden",
    gene_sets: List[str] = None,
    figsize_per_cluster: Tuple[int, int] = (16, 8),
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
) -> Dict[str, bool]:
    """
    Plot violin plots of module scores for each cluster in parallel.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with computed module scores.
    adata : AnnData
        AnnData object with cluster annotations.
    cluster_key : str or sequence of str
        Column name(s) in adata.obs containing cluster assignments.
    gene_sets : list
        Gene sets for enrichment annotation labels. Defaults to ``config.ENRICHMENT_GENE_SETS``.
    figsize_per_cluster : tuple
        Figure size for each cluster subplot.
    skip_existing : bool
        If True, skip existing plots.
    n_jobs : int, optional
        Number of parallel worker processes. Default uses config.N_JOBS.

    Returns
    -------
    dict
        Dictionary mapping plot names to generation status.
    """
    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)

    results = {}

    # Get module scores
    module_scores = hotspot_obj.module_scores

    if module_scores is None or module_scores.empty:
        log_warning(
            "HotspotPlotting.ViolinPlot", "No module scores available for violin plots"
        )
        print("  Warning: No module scores available for violin plots")
        return results

    # Check and resolve cluster key if multi-key
    if not isinstance(cluster_key, str) or "," in cluster_key or cluster_key not in adata.obs.columns:
        try:
            from ..preprocessing import resolve_cluster_key
            adata, cluster_key = resolve_cluster_key(adata, cluster_key)
        except Exception:
            from ..preprocessing import resolve_cluster_key_name
            cluster_key = resolve_cluster_key_name(cluster_key)

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
    module_labels = _get_module_enrichment_labels(
        hotspot_obj, gene_sets, n_jobs=n_jobs
    )

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

    tasks = [
        {
            "variant": variant,
            "scores_melted": scores_melted,
            "clusters": clusters,
            "modules": modules,
            "cluster_key": cluster_key,
            "skip_existing": skip_existing,
        }
        for variant in ["per_cluster", "all_clusters", "horizontal"]
    ]

    results_list = run_parallel_tasks(
        _plot_hotspot_violin_single_worker,
        tasks,
        n_jobs=n_jobs,
        desc="hotspot_violins",
    )

    return dict(results_list)


def generate_all_hotspot_plots(
    hotspot_obj,
    adata: Optional[AnnData] = None,
    cluster_key: Union[str, Sequence[str]] = "leiden",
    gene_sets: List[str] = None,
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Generate all Hotspot plots in parallel.

    Parameters
    ----------
    hotspot_obj : hotspot.Hotspot
        Hotspot object with analysis results.
    adata : AnnData, optional
        AnnData object with cluster annotations (for violin plots).
    cluster_key : str or sequence of str
        Column name(s) for cluster assignments.
    gene_sets : list
        Gene sets for enrichment analysis. Defaults to ``config.ENRICHMENT_GENE_SETS``.
    skip_existing : bool
        If True, skip existing plots.
    n_jobs : int, optional
        Number of parallel worker processes. Default uses config.N_JOBS.

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
        hotspot_obj, gene_sets=gene_sets, skip_existing=skip_existing, n_jobs=n_jobs
    )

    # Violin plots (if adata provided)
    if adata is not None:
        violin_results = plot_module_scores_violin(
            hotspot_obj,
            adata,
            cluster_key=cluster_key,
            gene_sets=gene_sets,
            skip_existing=skip_existing,
            n_jobs=n_jobs,
        )
        results.update(violin_results)

    generated = sum(1 for v in results.values() if v is True)
    skipped = sum(1 for v in results.values() if v is False)
    print(f"  Hotspot plots: {generated} generated, {skipped} skipped")

    return results
