"""
GRN Deep Analysis module for GeneCircuitry

This module provides functions for deep analysis of Gene Regulatory Networks (GRN).
Plotting functions have been consolidated in genecircuitry.plotting.grn_plots;
the names below are kept for backward compatibility and delegate to that module.
"""

import os
import pickle
import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Any

from . import config
from .logging_utils import log_error, log_warning


# ---------------------------------------------------------------------------
# Data helpers (non-plotting)
# ---------------------------------------------------------------------------


def process_single_links_file(
    links_file: str,
) -> pd.DataFrame:
    """
    Read a links pickle file and return a DataFrame.

    Args:
        links_file (str): Path to the links file.
    Returns:
        pd.DataFrame: DataFrame constructed from the links file.
    """
    links_df = pd.DataFrame()
    with open(links_file, "rb") as f:
        links_pickle = pickle.load(f)

    for key in links_pickle.keys():
        df_tmp = pd.DataFrame(links_pickle[key])
        df_tmp["cluster"] = key
        links_df = pd.concat([links_df, df_tmp], axis=0)

    return links_df


def process_single_score_file(score_file: str) -> pd.DataFrame:
    """
    Process a single score CSV file and return a DataFrame.

    Args:
        score_file (str): Path to the score CSV file.

    Returns:
        pd.DataFrame: Processed DataFrame with scores.
    """
    score_df = pd.read_csv(score_file, index_col=0)
    score_df.index.name = "gene"
    stratification_name = score_file.split("/")[-3]
    score_df["stratification"] = stratification_name

    return score_df


def merge_scores(tracked_file: str) -> pd.DataFrame:
    """
    Merge multiple score CSV files into a single DataFrame.
    This function reads multiple CSV files containing scores, adds a stratification
    column based on the parent directory name of each file, and concatenates them
    into a single DataFrame.
    Args:
        tracked_file (str): A file path to a text file containing paths to CSV files with scores.
    Returns:
        pd.DataFrame: A concatenated DataFrame containing all scores from the input
            files, with an additional 'stratification' column indicating the source
            stratification (extracted from the parent directory name of each file).
    Note:
        The stratification name is extracted from the second-to-last component of
        the file path (i.e., file.split("/")[-2]).
        Files are concatenated along axis=0 (rows).
    """
    with open(tracked_file) as f:
        file_list = f.readlines()
    file_list = [file.strip() for file in file_list]

    merged_scores = list()
    for file in file_list:
        if file.count("grn_merged_scores.csv") == 0:
            continue
        print(f"Processing file: {file}")
        score_tmp_df = pd.read_csv(file)
        score_tmp_df.rename(columns={"Unnamed: 0": "gene"}, inplace=True)
        stratification_name = file.split("/")[-3]
        score_tmp_df["stratification"] = stratification_name
        merged_scores.append(score_tmp_df)

    # concatenate all DataFrames
    merged_scores_df = pd.concat(merged_scores, axis=0)
    merged_scores_df.to_csv(
        f"{config.OUTPUT_DIR}/celloracle/total_merged_scores.csv", index=False
    )
    print(
        f"\n✓ Merged GRN scores saved to: {config.OUTPUT_DIR}/celloracle/total_merged_scores.csv"
    )

    return merged_scores_df


# ---------------------------------------------------------------------------
# Plotting helpers — delegating to genecircuitry.plotting.grn_plots
# ---------------------------------------------------------------------------


def plot_network_graph(
    score_df: pd.DataFrame,
    links_df: pd.DataFrame,
    scores=None,
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
):
    """
    Plot network graphs for multiple scores.

    Delegates to genecircuitry.plotting.grn_plots.plot_network_graph.
    """
    from .plotting.grn_plots import plot_network_graph as _impl

    return _impl(
        score_df=score_df,
        links_df=links_df,
        scores=scores,
        skip_existing=skip_existing,
        n_jobs=n_jobs,
    )


def plot_heatmap_single_score(
    heatmap_data: pd.DataFrame,
    cluster1: str,
    cluster2: str,
    score: str,
    top_n_genes: int = 5,
    skip_existing: bool = True,
):
    """
    Plot a heatmap for a single score and stratification.

    Delegates to genecircuitry.plotting.grn_plots.plot_heatmap_single_score.
    """
    from .plotting.grn_plots import plot_heatmap_single_score as _impl

    return _impl(
        heatmap_data=heatmap_data,
        cluster1=cluster1,
        cluster2=cluster2,
        score=score,
        top_n_genes=top_n_genes,
        skip_existing=skip_existing,
    )


def plot_heatmap_scores(
    scores_df: pd.DataFrame,
    top_n_genes: int = 10,
    scores=None,
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
):
    """
    Plot heatmaps for multiple scores across cluster combinations.

    Delegates to genecircuitry.plotting.grn_plots.plot_heatmap_scores.
    """
    from .plotting.grn_plots import plot_heatmap_scores as _impl

    return _impl(
        scores_df=scores_df,
        top_n_genes=top_n_genes,
        scores=scores,
        skip_existing=skip_existing,
        n_jobs=n_jobs,
    )


def plot_scatter_scores(
    score_df: pd.DataFrame,
    scores_list=None,
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
):
    """
    Plot scatter plots comparing scores between cluster pairs.

    Delegates to genecircuitry.plotting.grn_plots.plot_scatter_scores.
    """
    from .plotting.grn_plots import plot_scatter_scores as _impl

    return _impl(
        score_df=score_df,
        scores_list=scores_list,
        skip_existing=skip_existing,
        n_jobs=n_jobs,
    )


def plot_difference_cluster_scores(
    score_df: pd.DataFrame,
    scores=None,
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
):
    """
    Plot rank plots showing score differences between clusters.

    Delegates to genecircuitry.plotting.grn_plots.plot_difference_cluster_scores.
    """
    from .plotting.grn_plots import plot_difference_cluster_scores as _impl

    return _impl(
        score_df=score_df,
        scores=scores,
        skip_existing=skip_existing,
        n_jobs=n_jobs,
    )


def plot_compare_cluster_scores(
    score_df: pd.DataFrame,
    scores=None,
    skip_existing: bool = True,
    n_jobs: Optional[int] = None,
):
    """
    Compare scores across clusters and plot differences.

    Delegates to genecircuitry.plotting.grn_plots.plot_compare_cluster_scores.
    """
    from .plotting.grn_plots import plot_compare_cluster_scores as _impl

    return _impl(
        score_df=score_df,
        scores=scores,
        skip_existing=skip_existing,
        n_jobs=n_jobs,
    )
