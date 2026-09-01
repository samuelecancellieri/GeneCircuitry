"""
QC Plotting Module for GeneCircuitry
==============================

Functions for generating quality control visualizations including:
- Violin plots of QC metrics (pre/post filtering)
- Scatter plots of counts vs genes
- Distribution plots
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Tuple, Dict, Any
from anndata import AnnData

from .. import config
from .utils import save_plot, plot_exists, run_parallel_tasks


def plot_qc_violin_pre_filter(
    adata: AnnData,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Create violin plot for pre-filtering QC metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC metrics computed.
    save_name : str
        Suffix for the saved filename.
    figsize : tuple, optional
        Figure size (width, height). If None, uses config.PLOT_FIGSIZE_LARGE.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = f"{config.FIGURES_DIR_QC}/violin_pre_filter_{save_name}.png"

    if plot_exists(filepath, skip_existing):
        return False

    if figsize is None:
        figsize = config.PLOT_FIGSIZE_LARGE

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    pastel = sns.color_palette(getattr(config, "PLOT_CATEGORICAL_PALETTE", "pastel"))

    sns.violinplot(data=adata.obs, y="n_genes_by_counts", ax=axes[0], inner="box", color=pastel[0])
    axes[0].set_ylabel("Number of genes")
    axes[0].set_title("Genes per cell (Pre-filter)")

    sns.violinplot(data=adata.obs, y="total_counts", ax=axes[1], inner="box", color=pastel[1 % len(pastel)])
    axes[1].set_ylabel("Total counts")
    axes[1].set_title("Total counts per cell (Pre-filter)")

    sns.violinplot(data=adata.obs, y="pct_counts_mt", ax=axes[2], inner="box", color=pastel[2 % len(pastel)])
    axes[2].set_ylabel("% Mitochondrial counts")
    axes[2].set_title("Mitochondrial percentage (Pre-filter)")

    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="qc",
        metadata={"step": "pre_filter", "plot_name": "violin"},
        skip_existing=False,  # Already checked above
    )


def plot_qc_violin_post_filter(
    adata: AnnData,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Create violin plot for post-filtering QC metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC metrics computed (after filtering).
    save_name : str
        Suffix for the saved filename.
    figsize : tuple, optional
        Figure size (width, height). If None, uses config.PLOT_FIGSIZE_LARGE.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = f"{config.FIGURES_DIR_QC}/violin_post_filter_{save_name}.png"

    if plot_exists(filepath, skip_existing):
        return False

    if figsize is None:
        figsize = config.PLOT_FIGSIZE_LARGE

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    pastel = sns.color_palette(getattr(config, "PLOT_CATEGORICAL_PALETTE", "pastel"))

    sns.violinplot(data=adata.obs, y="n_genes_by_counts", ax=axes[0], inner="box", color=pastel[0])
    axes[0].set_ylabel("Number of genes")
    axes[0].set_title("Genes per cell (Post-filter)")

    sns.violinplot(data=adata.obs, y="total_counts", ax=axes[1], inner="box", color=pastel[1 % len(pastel)])
    axes[1].set_ylabel("Total counts")
    axes[1].set_title("Total counts per cell (Post-filter)")

    sns.violinplot(data=adata.obs, y="pct_counts_mt", ax=axes[2], inner="box", color=pastel[2 % len(pastel)])
    axes[2].set_ylabel("% Mitochondrial counts")
    axes[2].set_title("Mitochondrial percentage (Post-filter)")

    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="qc",
        metadata={"step": "post_filter", "plot_name": "violin"},
        skip_existing=False,
    )


def plot_qc_scatter_pre_filter(
    adata: AnnData,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Create scatter plot of counts vs genes (pre-filtering).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC metrics computed.
    save_name : str
        Suffix for the saved filename.
    figsize : tuple, optional
        Figure size (width, height). If None, uses config.PLOT_FIGSIZE_MEDIUM.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = f"{config.FIGURES_DIR_QC}/scatter_pre_filter_{save_name}.png"

    if plot_exists(filepath, skip_existing):
        return False

    if figsize is None:
        figsize = config.PLOT_FIGSIZE_MEDIUM

    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        cmap=config.PLOT_COLOR_PALETTE,
        alpha=0.7,
        s=5,
    )

    ax.set_xlabel("Total counts")
    ax.set_ylabel("Number of genes")
    ax.set_title("Genes vs Total Counts (Pre-filter)")

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("% Mitochondrial counts")

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="qc",
        metadata={"step": "pre_filter", "plot_name": "scatter"},
        skip_existing=False,
    )


def plot_qc_scatter_post_filter(
    adata: AnnData,
    save_name: str = "default",
    figsize: Optional[Tuple[int, int]] = None,
    skip_existing: bool = True,
) -> bool:
    """
    Create scatter plot of counts vs genes (post-filtering).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC metrics computed (after filtering).
    save_name : str
        Suffix for the saved filename.
    figsize : tuple, optional
        Figure size (width, height). If None, uses config.PLOT_FIGSIZE_MEDIUM.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = f"{config.FIGURES_DIR_QC}/scatter_post_filter_{save_name}.png"

    if plot_exists(filepath, skip_existing):
        return False

    if figsize is None:
        figsize = config.PLOT_FIGSIZE_MEDIUM

    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        cmap=config.PLOT_COLOR_PALETTE,
        alpha=0.7,
        s=5,
    )

    ax.set_xlabel("Total counts")
    ax.set_ylabel("Number of genes")
    ax.set_title("Genes vs Total Counts (Post-filter)")

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("% Mitochondrial counts")

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="qc",
        metadata={"step": "post_filter", "plot_name": "scatter"},
        skip_existing=False,
    )


def _plot_qc_single_worker(task: Dict[str, Any]) -> Tuple[str, bool]:
    """Worker function for a single QC plot."""
    plot_type = task["plot_type"]
    adata = task["adata"]
    save_name = task["save_name"]
    skip_existing = task.get("skip_existing", True)
    figsize = task.get("figsize")

    if plot_type == "violin_pre_filter":
        res = plot_qc_violin_pre_filter(
            adata, save_name=save_name, figsize=figsize, skip_existing=skip_existing
        )
    elif plot_type == "violin_post_filter":
        res = plot_qc_violin_post_filter(
            adata, save_name=save_name, figsize=figsize, skip_existing=skip_existing
        )
    elif plot_type == "scatter_pre_filter":
        res = plot_qc_scatter_pre_filter(
            adata, save_name=save_name, figsize=figsize, skip_existing=skip_existing
        )
    elif plot_type == "scatter_post_filter":
        res = plot_qc_scatter_post_filter(
            adata, save_name=save_name, figsize=figsize, skip_existing=skip_existing
        )
    else:
        res = False

    return (plot_type, res)


def generate_all_qc_plots(
    adata_pre: AnnData,
    adata_post: AnnData,
    save_name: str = "default",
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
) -> Dict[str, bool]:
    """
    Generate all QC plots for pre and post-filtering data in parallel.

    Parameters
    ----------
    adata_pre : AnnData
        Data before filtering (with QC metrics computed).
    adata_post : AnnData
        Data after filtering.
    save_name : str
        Suffix for saved filenames.
    skip_existing : bool
        If True, skip existing plots.
    n_jobs : int, optional
        Number of parallel worker processes. Default uses config.N_JOBS.

    Returns
    -------
    dict
        Dictionary mapping plot names to generation status (True/False).
    """
    print("Generating QC plots...")

    tasks = [
        {
            "plot_type": "violin_pre_filter",
            "adata": adata_pre,
            "save_name": save_name,
            "skip_existing": skip_existing,
        },
        {
            "plot_type": "violin_post_filter",
            "adata": adata_post,
            "save_name": save_name,
            "skip_existing": skip_existing,
        },
        {
            "plot_type": "scatter_pre_filter",
            "adata": adata_pre,
            "save_name": save_name,
            "skip_existing": skip_existing,
        },
        {
            "plot_type": "scatter_post_filter",
            "adata": adata_post,
            "save_name": save_name,
            "skip_existing": skip_existing,
        },
    ]

    results_list = run_parallel_tasks(
        _plot_qc_single_worker,
        tasks,
        n_jobs=n_jobs,
        desc="qc_plots",
    )

    results = {plot_type: status for plot_type, status in results_list}

    generated = sum(1 for status in results.values() if status)
    skipped = len(results) - generated
    print(f"  QC plots: {generated} generated, {skipped} skipped")

    return results
