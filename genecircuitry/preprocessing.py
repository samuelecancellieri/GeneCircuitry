"""
Data preprocessing module for genecircuitry
"""

import os
import re
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from typing import Optional, Tuple, List, Union, Any
from collections.abc import Iterable, Sequence
from anndata import AnnData

from . import config


def perform_qc(
    adata: AnnData,
    min_genes: Optional[int] = None,
    min_counts: Optional[int] = None,
    max_counts: Optional[int] = None,
    pct_counts_mt_max: Optional[float] = None,
    min_cells: Optional[int] = None,
    plot: bool = True,
    figsize: Optional[Tuple[int, int]] = None,
    save_plots: Optional[str] = None,
) -> AnnData:
    """
    Perform quality control on AnnData object.

    This function calculates QC metrics, filters cells based on specified thresholds,
    and generates QC visualization plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations and genes as variables.
    min_genes : int, optional
        Minimum number of genes expressed required for a cell to pass filtering.
        If None, uses config.QC_MIN_GENES (default: 200).
    min_counts : int, optional
        Minimum number of counts required for a cell to pass filtering.
        If None, uses config.QC_MIN_COUNTS (default: 500).
    max_counts : int, optional
        Maximum number of counts allowed for a cell to pass filtering.
        If None, uses config.QC_MAX_COUNTS (default: None, no upper limit).
    pct_counts_mt_max : float, optional
        Maximum percentage of mitochondrial counts allowed for a cell.
        If None, uses config.QC_PCT_MT_MAX (default: 20.0).
    min_cells : int, optional
        Minimum number of cells expressing a gene for the gene to be kept.
        If None, uses config.QC_MIN_CELLS (default: 10).
    plot : bool, default=True
        Whether to generate QC plots.
    figsize : tuple, optional
        Figure size for QC plots (width, height).
        If None, uses config.PLOT_FIGSIZE_LARGE (default: (15, 10)).
    save_plots : str, optional
        Path to save the QC plots. If None, plots are not saved.

    Returns
    -------
    AnnData
        Filtered AnnData object with QC metrics stored in .obs and .var.

    Examples
    --------
    >>> import scanpy as sc
    >>> from genecircuitry.preprocessing import perform_qc
    >>> from genecircuitry import config
    >>>
    >>> # Use default config values
    >>> adata = sc.read_h5ad('data.h5ad')
    >>> adata_qc = perform_qc(adata)
    >>>
    >>> # Override specific parameters
    >>> adata_qc = perform_qc(adata, min_genes=300, min_counts=1000)
    >>>
    >>> # Use updated config
    >>> config.update_config(QC_MIN_GENES=500)
    >>> adata_qc = perform_qc(adata)
    """

    # Use config defaults if not specified
    if min_genes is None:
        min_genes = config.QC_MIN_GENES
    if min_counts is None:
        min_counts = config.QC_MIN_COUNTS
    if max_counts is None:
        max_counts = config.QC_MAX_COUNTS
    if pct_counts_mt_max is None:
        pct_counts_mt_max = config.QC_PCT_MT_MAX
    if min_cells is None:
        min_cells = config.QC_MIN_CELLS
    if figsize is None:
        figsize = config.PLOT_FIGSIZE_MEDIUM

    # Make a copy to avoid modifying the original
    adata_cc = adata.copy()

    # Store initial cell and gene counts
    n_cells_initial = adata_cc.n_obs
    n_genes_initial = adata_cc.n_vars

    print(f"Initial data shape: {n_cells_initial} cells × {n_genes_initial} genes")

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata_cc.var["mt"] = adata_cc.var_names.str.startswith("MT-")
    # ribosomal genes
    adata_cc.var["ribo"] = adata_cc.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata_cc.var["hb"] = adata_cc.var_names.str.contains("^HB[^(P)]")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata_cc,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=True,
        inplace=True,
    )

    # Create violin plot for pre-filtering metrics
    if plot:
        fig, axes = plt.subplots(1, 3, figsize=config.PLOT_FIGSIZE_LARGE)

        pastel = sns.color_palette(getattr(config, "PLOT_CATEGORICAL_PALETTE", "pastel"))
        sns.violinplot(data=adata_cc.obs, y="n_genes_by_counts", ax=axes[0], inner="box", color=pastel[0])
        axes[0].set_ylabel("Number of genes")
        axes[0].set_title("Genes per cell")

        sns.violinplot(data=adata_cc.obs, y="total_counts", ax=axes[1], inner="box", color=pastel[1 % len(pastel)])
        axes[1].set_ylabel("Total counts")
        axes[1].set_title("Total counts per cell")

        sns.violinplot(data=adata_cc.obs, y="pct_counts_mt", ax=axes[2], inner="box", color=pastel[2 % len(pastel)])
        axes[2].set_ylabel("% Mitochondrial counts")
        axes[2].set_title("Mitochondrial percentage")

        if save_plots:
            fig.savefig(
                f"{config.FIGURES_DIR_QC}/violin_pre_filter_{save_plots}.png",
                dpi=config.SAVE_DPI,
                bbox_inches="tight",
            )

        # Create scatter plot for pre-filtering metrics
        fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_MEDIUM)
        scatter = ax.scatter(
            adata_cc.obs["total_counts"],
            adata_cc.obs["n_genes_by_counts"],
            c=adata_cc.obs["pct_counts_mt"],
            cmap=config.PLOT_COLOR_PALETTE,
            alpha=0.7,
            s=5,
        )
        ax.set_xlabel("Total counts")
        ax.set_ylabel("Number of genes")
        ax.set_title("Genes vs Total Counts colored by % Mitochondrial")
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("% Mitochondrial counts")
        if save_plots:
            fig.savefig(
                f"{config.FIGURES_DIR_QC}/scatter_pre_filter_{save_plots}.png",
                dpi=config.SAVE_DPI,
                bbox_inches="tight",
            )

    # Apply filtering
    print("\nApplying filters...")
    print(f"  - Minimum genes per cell: {min_genes}")
    print(f"  - Minimum counts per cell: {min_counts}")
    if max_counts:
        print(f"  - Maximum counts per cell: {max_counts}")
    print(f"  - Maximum mitochondrial percentage: {pct_counts_mt_max}%")

    # Filter cells
    sc.pp.filter_cells(adata_cc, min_genes=min_genes)
    sc.pp.filter_cells(adata_cc, min_counts=min_counts)

    if max_counts:
        sc.pp.filter_cells(adata_cc, max_counts=max_counts)

    adata_cc = adata_cc[adata_cc.obs["pct_counts_mt"] < pct_counts_mt_max, :].copy()

    # Filter genes (keep genes expressed in at least min_cells to preserve rare cell populations)
    sc.pp.filter_genes(adata_cc, min_cells=min_cells)

    # Report filtering results
    n_cells_filtered = adata_cc.n_obs
    n_genes_filtered = adata_cc.n_vars
    cells_removed = n_cells_initial - n_cells_filtered
    genes_removed = n_genes_initial - n_genes_filtered

    print(f"\nFiltering results:")
    print(
        f"  - Cells removed: {cells_removed} ({cells_removed/n_cells_initial*100:.2f}%)"
    )
    print(
        f"  - Genes removed: {genes_removed} ({genes_removed/n_genes_initial*100:.2f}%)"
    )
    print(f"  - Final shape: {n_cells_filtered} cells × {n_genes_filtered} genes")

    # Create seaborn violin plot for post-filtering metrics
    if plot:
        fig, axes = plt.subplots(1, 3, figsize=config.PLOT_FIGSIZE_LARGE)
        pastel = sns.color_palette(getattr(config, "PLOT_CATEGORICAL_PALETTE", "pastel"))

        sns.violinplot(
            data=adata_cc.obs, y="n_genes_by_counts", ax=axes[0], inner="box", color=pastel[0]
        )
        axes[0].set_ylabel("Number of genes")
        axes[0].set_title("Genes per cell")

        sns.violinplot(data=adata_cc.obs, y="total_counts", ax=axes[1], inner="box", color=pastel[1 % len(pastel)])
        axes[1].set_ylabel("Total counts")
        axes[1].set_title("Total counts per cell")

        sns.violinplot(data=adata_cc.obs, y="pct_counts_mt", ax=axes[2], inner="box", color=pastel[2 % len(pastel)])
        axes[2].set_ylabel("% Mitochondrial counts")
        axes[2].set_title("Mitochondrial percentage")

        if save_plots:
            fig.savefig(
                f"{config.FIGURES_DIR_QC}/violin_post_filter_{save_plots}.png",
                dpi=config.SAVE_DPI,
                bbox_inches="tight",
            )

        # scatter of post-filtering metrics
        fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_MEDIUM)
        scatter = ax.scatter(
            adata_cc.obs["total_counts"],
            adata_cc.obs["n_genes_by_counts"],
            c=adata_cc.obs["pct_counts_mt"],
            cmap=config.PLOT_COLOR_PALETTE,
            alpha=0.7,
        )
        ax.set_xlabel("Total counts")
        ax.set_ylabel("Number of genes")
        ax.set_title("Genes vs Total Counts colored by % Mitochondrial")
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("% Mitochondrial counts")
        if save_plots:
            fig.savefig(
                f"{config.FIGURES_DIR_QC}/scatter_post_filter_{save_plots}.png",
                dpi=config.SAVE_DPI,
                bbox_inches="tight",
            )

        plt.close("all")

    return adata_cc


def perform_normalization(adata: AnnData) -> AnnData:
    """
    Perform normalization on AnnData object.

    This function normalizes the data to a target sum and applies log transformation.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations and genes as variables.

    Returns
    -------
    AnnData
        Normalized AnnData object.

    Examples
    --------
    >>> import scanpy as sc
    >>> from genecircuitry.preprocessing import perform_normalization
    >>> from genecircuitry import config
    >>>
    >>> adata = sc.read_h5ad('data.h5ad')
    >>> adata_normalized = perform_normalization(adata)
    """

    print("\nPerforming normalization...")
    # Make a copy to avoid modifying the original
    adata_cc = adata.copy()

    # Check if data is already normalized (heuristic: max count > 100 indicates raw counts)
    if adata_cc.X.max() > 100:
        print(
            "Warning: Data appears to be unnormalized (max count > 100). Proceeding with normalization."
        )
    else:
        print(
            "Data appears to be already normalized (max count <= 100). Skipping normalization."
        )
        return adata_cc

    # Saving raw count data
    if "raw_counts" not in adata_cc.layers:
        adata_cc.layers["raw_counts"] = adata_cc.X.copy()
        print("Stored raw counts in layer 'raw_counts'")

    # Normalize to target sum
    sc.pp.normalize_total(adata_cc, target_sum=config.NORMALIZE_TARGET_SUM)
    print(f"Normalized to target sum: {config.NORMALIZE_TARGET_SUM}")

    # Log transform
    sc.pp.log1p(adata_cc)
    print("Applied log1p transformation")

    return adata_cc


def perform_dimensionality_reduction_clustering(
    adata: AnnData,
    cluster_key: str = "leiden",
    skip_dimensionality_reduction: bool = False,
    force: bool = False,
    **kwargs,
) -> AnnData:
    """
    Perform dimensionality reduction and clustering on AnnData object.

    This function identifies highly variable genes, computes PCA, constructs a
    neighborhood graph, computes UMAP embedding, and performs Leiden clustering.
    If `skip_dimensionality_reduction` is True, all steps are skipped.
    If `force` is False, each step checks whether the corresponding result is already
    present in `adata` and skips that individual step if already completed.
    If `force` is True, all dimensionality reduction and clustering steps are re-run.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations and genes as variables.
    cluster_key : str, default="leiden"
        Key to store cluster labels in `adata.obs`. Also used to check if clustering
        is already present.
    skip_dimensionality_reduction : bool, default=False
        Whether to skip all dimensionality reduction and clustering steps.
    force : bool, default=False
        Whether to force re-computation of dimensionality reduction and clustering
        steps even if results already exist in `adata`.

    Returns
    -------
    AnnData
        Annotated data matrix with PCA, UMAP, and clustering results.
    """
    force = (
        force
        or kwargs.get("force_dim_reduction", False)
        or kwargs.get("force_dimensionality_reduction", False)
    )
    adata_cc = adata.copy()
    if hasattr(adata_cc.X, "dtype") and not np.issubdtype(
        adata_cc.X.dtype, np.floating
    ):
        adata_cc.X = adata_cc.X.astype(np.float32)

    if skip_dimensionality_reduction:
        print("\nSkipping dimensionality reduction and clustering (as requested)...")
        resolved_cluster_key = resolve_cluster_key_name(cluster_key)
        parsed = parse_cluster_keys(cluster_key)
        if len(parsed) > 1 and all(k in adata_cc.obs.columns for k in parsed):
            adata_cc, resolved_cluster_key = resolve_cluster_key(adata_cc, cluster_key)
        adata_cc = ensure_categorical_obs(
            adata_cc, columns=[resolved_cluster_key, cluster_key]
        )
        return adata_cc

    print("\nPerforming dimensionality reduction and clustering...")

    # Step 1: Highly Variable Genes (HVG)
    if not force and "highly_variable" in adata_cc.var.columns:
        n_hvgs = int(adata_cc.var["highly_variable"].sum())
        print(
            f"  - Highly variable genes already identified ({n_hvgs} HVGs found in .var['highly_variable']). Skipping HVG selection."
        )
    else:
        sc.pp.highly_variable_genes(
            adata_cc,
            n_top_genes=config.HVGS_N_TOP_GENES,
            subset=False,
        )
        print(f"Identified top {config.HVGS_N_TOP_GENES} highly variable genes")

    # Step 2: PCA
    if not force and "X_pca" in adata_cc.obsm:
        n_comps = adata_cc.obsm["X_pca"].shape[1]
        print(
            f"  - PCA already computed ({n_comps} components found in .obsm['X_pca']). Skipping PCA."
        )
    else:
        n_comps = config.PCA_N_COMPS
        if n_comps is not None and config.PCA_SVD_SOLVER == "arpack":
            n_comps = min(n_comps, min(adata_cc.shape) - 1)
        sc.pp.pca(adata_cc, n_comps=n_comps, svd_solver=config.PCA_SVD_SOLVER)
        print(f"Computed PCA with {n_comps} components")

    # Step 3: Neighborhood graph
    if not force and (
        "neighbors" in adata_cc.uns
        or ("connectivities" in adata_cc.obsp and "distances" in adata_cc.obsp)
    ):
        print(
            "  - Neighborhood graph already constructed (found in .uns['neighbors'] / .obsp). Skipping neighbors."
        )
    else:
        n_pcs = config.NEIGHBORS_N_PCS
        if "X_pca" in adata_cc.obsm and n_pcs is not None:
            n_pcs = min(n_pcs, adata_cc.obsm["X_pca"].shape[1])
        sc.pp.neighbors(
            adata_cc,
            metric=config.NEIGHBORS_METRIC,
            method=config.NEIGHBORS_METHOD,
            n_neighbors=config.NEIGHBORS_N_NEIGHBORS,
            n_pcs=n_pcs,
        )
        print(
            f"Constructed neighborhood graph with {config.NEIGHBORS_N_NEIGHBORS} neighbors"
        )

    # Step 4: UMAP
    if not force and "X_umap" in adata_cc.obsm:
        print(
            "  - UMAP embedding already computed (found in .obsm['X_umap']). Skipping UMAP."
        )
    else:
        sc.tl.umap(
            adata_cc,
            min_dist=config.UMAP_MIN_DIST,
            spread=config.UMAP_SPREAD,
            n_components=config.UMAP_N_COMPONENTS,
        )
        print("Computed UMAP embedding")

    # Resolve cluster_key if multi-key or custom
    resolved_cluster_key = resolve_cluster_key_name(cluster_key)
    parsed = parse_cluster_keys(cluster_key)
    if len(parsed) > 1 and all(k in adata_cc.obs.columns for k in parsed):
        adata_cc, resolved_cluster_key = resolve_cluster_key(adata_cc, cluster_key)

    # Step 5: Clustering
    if not force and resolved_cluster_key in adata_cc.obs.columns:
        n_clusters = len(adata_cc.obs[resolved_cluster_key].unique())
        print(
            f"  - Clustering '{resolved_cluster_key}' already exists ({n_clusters} clusters found in .obs['{resolved_cluster_key}']). Skipping clustering."
        )
    else:
        sc.tl.leiden(
            adata_cc,
            key_added="leiden",
            resolution=config.LEIDEN_RESOLUTION,
            flavor="igraph",
            n_iterations=2,
        )
        print(
            f"Performed Leiden clustering with resolution {config.LEIDEN_RESOLUTION} (key: 'leiden')"
        )

    # Convert categorical columns for stratification compatibility
    adata_cc = ensure_categorical_obs(
        adata_cc, columns=[resolved_cluster_key, cluster_key]
    )

    return adata_cc


def ensure_categorical_obs(
    adata: AnnData,
    columns: Optional[list] = None,
) -> AnnData:
    """
    Convert object/string columns in adata.obs to pandas Categorical type.

    This ensures consistent behavior during stratification, avoiding conflicts
    between string comparisons and numeric values.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    columns : list, optional
        Specific columns to convert. If None, converts all object/string columns
        and common clustering columns (leiden, louvain, cell_type, etc.).

    Returns
    -------
    AnnData
        AnnData with categorical columns in .obs.

    Examples
    --------
    >>> adata = ensure_categorical_obs(adata)
    >>> adata = ensure_categorical_obs(adata, columns=['cell_type', 'batch'])
    """
    if adata.is_view:
        adata = adata.copy()

    # Common stratification/clustering columns to always convert if present
    default_categorical_cols = [
        "leiden",
        "louvain",
        "cell_type",
        "celltype",
        "cluster",
        "clusters",
        "batch",
        "sample",
        "condition",
    ]

    if columns is None:
        # Auto-detect: object dtype columns + known categorical columns
        columns_to_convert = []

        # Add object/string dtype columns
        for col in adata.obs.columns:
            if (
                adata.obs[col].dtype == "object"
                or adata.obs[col].dtype.name == "string"
            ):
                columns_to_convert.append(col)

        # Add default categorical columns if they exist
        for col in default_categorical_cols:
            if col in adata.obs.columns and col not in columns_to_convert:
                columns_to_convert.append(col)
    else:
        columns_to_convert = []
        for k in parse_cluster_keys(columns):
            if k in adata.obs.columns and k not in columns_to_convert:
                columns_to_convert.append(k)

    converted = []
    for col in columns_to_convert:
        if not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            adata.obs[col] = adata.obs[col].astype("category")
            converted.append(col)

    if converted:
        print(f"Converted to categorical: {', '.join(converted)}")

    return adata


def parse_cluster_keys(keys: Any) -> List[str]:
    """
    Parse cluster key(s) from None, str, or sequence/iterable into a list of unique strings.

    Supports comma-separated strings (e.g. 'key1,key2' or 'key1, key2')
    and sequences/iterables of keys (e.g. ['key1', 'key2'], ('k1', 'k2'), {'k1', 'k2'},
    pd.Index, np.ndarray, dict_keys, generators, etc.).
    """
    if keys is None:
        return []
    if isinstance(keys, (str, bytes)):
        keys_str = keys.decode() if isinstance(keys, bytes) else keys
        parsed = []
        for k in keys_str.split(","):
            k_clean = k.strip()
            if k_clean and k_clean not in parsed:
                parsed.append(k_clean)
        return parsed
    if isinstance(keys, (Iterable, Sequence, np.ndarray, pd.Index, pd.Series)):
        parsed = []
        # For deterministic behavior if a set/frozenset is passed, sort by string representation
        items = (
            sorted(list(keys), key=lambda x: str(x))
            if isinstance(keys, (set, frozenset))
            else list(keys)
        )
        for item in items:
            if isinstance(item, (str, bytes)):
                item_str = item.decode() if isinstance(item, bytes) else item
                for sub in item_str.split(","):
                    sub_clean = sub.strip()
                    if sub_clean and sub_clean not in parsed:
                        parsed.append(sub_clean)
            elif item is not None:
                item_clean = str(item).strip()
                if item_clean and item_clean not in parsed:
                    parsed.append(item_clean)
        return parsed
    item_str = str(keys).strip()
    return [item_str] if item_str else []


def sanitize_identifier(val: Any) -> str:
    """
    Sanitize a value for filesystem-safe identifier and composite cluster naming.

    Replaces spaces with underscores, slashes with hyphens, and removes/replaces
    unsafe filename characters.
    """
    s = str(val).strip()
    s = s.replace(" ", "_").replace("/", "-")
    s = re.sub(r'[\\:*?"<>|]', "_", s)
    return s


def resolve_cluster_key_name(cluster_key: Any) -> str:
    """
    Return the resolved cluster key column name (e.g. 'key1_key2' for multi-key).
    """
    keys = parse_cluster_keys(cluster_key)
    if not keys:
        return "leiden"
    if len(keys) == 1:
        return keys[0]
    return "_".join(keys)


def resolve_cluster_key(
    adata: AnnData,
    cluster_key: Any,
    key_term: str = "Cluster key",
) -> Tuple[AnnData, str]:
    """
    Resolve cluster_key to a single column name in adata.obs.

    If cluster_key contains multiple keys (comma-separated string or sequence),
    constructs/ensures a composite categorical column in adata.obs representing
    the combined grouping (e.g. 'key1_key2') and returns (adata, composite_key_name).
    If cluster_key is a single key, ensures it is categorical and returns (adata, key).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cluster_key : str or sequence of str
        Single key, comma-separated keys, or sequence of keys.
    key_term : str, default="Cluster key"
        Descriptive term for error messages (e.g. 'Cluster key' or 'Cluster column').

    Returns
    -------
    Tuple[AnnData, str]
        (adata, resolved_cluster_key_name)
    """
    keys = parse_cluster_keys(cluster_key)
    if not keys:
        raise ValueError(f"{key_term} is required.")

    if len(keys) == 1:
        single_key = keys[0]
        if single_key not in adata.obs.columns:
            raise ValueError(
                f"{key_term} '{single_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        if adata.is_view:
            adata = adata.copy()
        if not isinstance(adata.obs[single_key].dtype, pd.CategoricalDtype):
            adata.obs[single_key] = adata.obs[single_key].astype("category")
            print(f"  Converted '{single_key}' to categorical")
        return adata, single_key

    # Multi-key case:
    missing = [k for k in keys if k not in adata.obs.columns]
    if missing:
        if len(missing) == 1:
            raise ValueError(
                f"{key_term} '{missing[0]}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        else:
            raise ValueError(
                f"{key_term}s {missing} not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

    if adata.is_view:
        adata = adata.copy()

    composite_col = "_".join(keys)

    # Ensure all constituent keys are categorical
    for k in keys:
        if not isinstance(adata.obs[k].dtype, pd.CategoricalDtype):
            adata.obs[k] = adata.obs[k].astype("category")

    # Build composite series from sanitized values
    combined = None
    for k in keys:
        col_str = adata.obs[k].astype(str).apply(sanitize_identifier)
        if combined is None:
            combined = col_str
        else:
            combined = combined + "_" + col_str

    adata.obs[composite_col] = combined.astype("category")
    print(
        f"  Constructed composite {key_term.lower()} column '{composite_col}' "
        f"with {len(adata.obs[composite_col].cat.categories)} categories"
    )
    return adata, composite_col
