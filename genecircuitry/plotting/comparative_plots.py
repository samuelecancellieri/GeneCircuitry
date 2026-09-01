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
    filepath = f"{config.FIGURES_DIR_COMPARATIVE}/comparative_module_activity{suffix}.png"

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
        annot=True,
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
    filepath = f"{config.FIGURES_DIR_COMPARATIVE}/comparative_pathway_enrichment{suffix}.png"

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
    plot_data["neg_log_p"] = -np.log10(np.clip(plot_data["adjusted_p_value"], 1e-30, 1.0))
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
    ax.set_title("Top Biological Pathways Enriched in Co-expression Modules", fontsize=13, pad=12)
    ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", linewidth=1, label="p = 0.05")
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
    filepath = f"{config.FIGURES_DIR_COMPARATIVE}/tf_module_regulatory_matrix{suffix}.png"

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
        annot=True,
        fmt="d",
        linewidths=0.5,
        cbar_kws={"label": "Number of Target Genes"},
        ax=ax,
    )

    ax.set_title("TF-to-Module Regulatory Mapping (Target Overlap)", fontsize=13, pad=12)
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
    sorted_df = tf_pivot_df.loc[tf_pivot_df.mean(axis=1).sort_values(ascending=False).index]

    if figsize is None:
        n_rows, n_cols = sorted_df.shape
        w = max(8, min(24, 4 + n_cols * 1.2))
        h = max(7, min(24, 3 + n_rows * 0.55))
        figsize = (w, h)

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        sorted_df,
        cmap="Blues",
        annot=True,
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
    ax.set_title("TF Target Gene Conservation vs. Rewiring Across Groups", fontsize=13, pad=12)
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
    filepath = f"{config.FIGURES_DIR_COMPARATIVE}/comparative_summary_dashboard{suffix}.png"

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
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            cbar_kws={"label": "Mean Score"},
            ax=axes[0, 0],
        )
        axes[0, 0].set_title("A. Co-expression Module Activity", fontsize=12, fontweight="bold")
        axes[0, 0].tick_params(axis="x", rotation=30)
    else:
        axes[0, 0].text(0.5, 0.5, "No Module Activity Data", ha="center", va="center")
        axes[0, 0].set_title("A. Co-expression Module Activity", fontsize=12, fontweight="bold")

    # Panel B: TF Centrality
    if tf_pivot_df is not None and not tf_pivot_df.empty:
        top_tfs_sub = tf_pivot_df.head(15)
        sns.heatmap(
            top_tfs_sub,
            cmap="Blues",
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            cbar_kws={"label": "Centrality"},
            ax=axes[0, 1],
        )
        axes[0, 1].set_title("B. Transcription Factor Centrality", fontsize=12, fontweight="bold")
        axes[0, 1].tick_params(axis="x", rotation=30)
    else:
        axes[0, 1].text(0.5, 0.5, "No TF Centrality Data", ha="center", va="center")
        axes[0, 1].set_title("B. Transcription Factor Centrality", fontsize=12, fontweight="bold")

    # Panel C: TF-to-Module Regulation
    if tf_mod_matrix is not None and not tf_mod_matrix.empty:
        sns.heatmap(
            tf_mod_matrix.head(15),
            cmap="YlOrRd",
            annot=True,
            fmt="d",
            linewidths=0.5,
            cbar_kws={"label": "Target Count"},
            ax=axes[1, 0],
        )
        axes[1, 0].set_title("C. TF-to-Module Regulation", fontsize=12, fontweight="bold")
        axes[1, 0].tick_params(axis="x", rotation=30)
    else:
        axes[1, 0].text(0.5, 0.5, "No TF-to-Module Data", ha="center", va="center")
        axes[1, 0].set_title("C. TF-to-Module Regulation", fontsize=12, fontweight="bold")

    # Panel D: Top Functional Pathways
    if enrichment_df is not None and not enrichment_df.empty:
        top_enr = enrichment_df.head(12).copy()
        top_enr["neg_log_p"] = -np.log10(np.clip(top_enr["adjusted_p_value"], 1e-30, 1.0))
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
        axes[1, 1].set_title("D. Top Enriched Module Pathways", fontsize=12, fontweight="bold")
    else:
        axes[1, 1].text(0.5, 0.5, "No Pathway Enrichment Data", ha="center", va="center")
        axes[1, 1].set_title("D. Top Enriched Module Pathways", fontsize=12, fontweight="bold")

    fig.suptitle("GeneCircuitry Comparative Analysis Dashboard", fontsize=16, fontweight="bold", y=0.99)
    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="comparative",
        metadata={"plot_name": "summary_dashboard"},
        skip_existing=False,
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
            tf_pivot_df, tf_summary_df=tf_summary_df, save_name=save_name, skip_existing=skip_existing
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

    generated = sum(1 for status in results.values() if status)
    print(f"  Comparative plots: {generated} generated, {len(results) - generated} skipped/unchanged")

    return results
