"""
CellOracle processing module for GeneCircuitry
"""

import scanpy as sc
import celloracle as co
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from typing import Optional, Tuple, Union, Sequence, List
from anndata import AnnData
import pickle

from . import config
from .preprocessing import (
    resolve_cluster_key,
    resolve_cluster_key_name,
    parse_cluster_keys,
    sanitize_identifier,
)


def load_celloracle_results(
    oracle_path: str,
    links_path: str,
) -> Tuple[co.Oracle, co.Links]:
    """
    Load CellOracle analysis results from specified files.

    Parameters:
        oracle_path (str): Path to the Oracle object file.
        links_path (str): Path to the Links object file.

    Returns:
        Tuple[co.Oracle, co.Links]: Loaded Oracle and Links objects.

    """

    # Load Oracle object
    oracle = co.load_hdf5(oracle_path)
    print(f"\n✓ Oracle object loaded from: {oracle_path}")

    # Load Links object
    links = co.load_hdf5(links_path)
    print(f"✓ GRN links loaded from: {links_path}")

    return oracle, links


def save_celloracle_results(
    oracle: co.Oracle,
    links: co.Links,
):
    """
    Save CellOracle analysis results to specified output directory.

    Parameters:
        oracle (co.Oracle): An instance of the Oracle class containing analysis results.
        links (co.Links): An instance of the Links class containing GRN links.

    """

    # Save Oracle object
    oracle_path = f"{config.OUTPUT_DIR}/celloracle/oracle_object.celloracle.oracle"
    oracle.to_hdf5(oracle_path)
    print(f"\n✓ Oracle object saved to: {oracle_path}")

    # Save links
    links_path = f"{config.OUTPUT_DIR}/celloracle/oracle_object.celloracle.links"
    links.to_hdf5(links_path)
    print(f"✓ GRN links saved to: {links_path}")

    merged_score = links.merged_score
    merged_score_path = f"{config.OUTPUT_DIR}/celloracle/grn_merged_scores.csv"
    merged_score.to_csv(merged_score_path)
    print(f"✓ Merged scores saved to: {merged_score_path}")

    links_filtered = links.filtered_links
    links_filtered_path = f"{config.OUTPUT_DIR}/celloracle/grn_filtered_links.pkl"
    with open(links_filtered_path, "wb") as f:
        pickle.dump(links_filtered, f)
    print(f"✓ Filtered links saved to: {links_filtered_path}")


# Move content from grn_celloracle_processing.py here
def create_oracle_object(
    adata: AnnData,
    cluster_column_name: str = "leiden",
    embedding_name: str = "X_draw_graph_fa",
    species: str = "human",
    TG_to_TF_dictionary: Optional[str] = None,
    raw_count_layer: Optional[str] = None,
    no_base_grn: bool = False,
    **kwargs,
):
    """
    Creates an Oracle object for CellOracle analysis.

    Parameters:
        adata (AnnData): Annotated data object containing gene expression data.
        TG_to_TF_dictionary (str): Path to a pickle file containing either:
            - A dictionary mapping target genes to transcription factors, or
            - An enriched ATAC peaks DataFrame (from process_atac_peaks).
            The file is loaded and added via oracle.addTFinfo_dictionary().
        no_base_grn (bool): If True, skip loading the species promoter base GRN.
            In this case a TG_to_TF_dictionary (or ATAC peaks pkl) must be
            provided to supply TF information; otherwise the Oracle object
            will have no TF data at all. Default: False.
        cluster_column_name (str): Name of the column in `adata.obs`
            that contains cluster information.
        embedding_name (str): Name of the embedding to be used.
        raw_count_layer (str, optional): Name of the layer in
            `adata.layers` that contains raw count data.
            Defaults to None.

    Returns:
        oracle (Oracle): An instance of the Oracle class.

    """
    cluster_column_name = kwargs.get(
        "cluster_key", kwargs.get("cluster_column", cluster_column_name)
    )
    if cluster_column_name is None:
        cluster_column_name = "leiden"

    # Create a copy of the AnnData object to avoid modifying the original
    adata_cc = adata.copy()

    # Resolve cluster column name and construct composite if multi-key
    adata_cc, cluster_column_name = resolve_cluster_key(
        adata_cc, cluster_column_name, key_term="Cluster column"
    )

    # Replace unsafe characters in cluster names
    adata_cc.obs[cluster_column_name] = (
        adata_cc.obs[cluster_column_name]
        .astype(str)
        .apply(sanitize_identifier)
        .astype("category")
    )

    # Create Oracle object
    oracle = co.Oracle()

    if raw_count_layer:
        print("Using raw counts layer for Oracle object creation.")
        # use raw counts to build the oracle object
        adata_cc.X = adata_cc.layers[raw_count_layer].copy()
        oracle.import_anndata_as_raw_count(
            adata=adata_cc,
            cluster_column_name=cluster_column_name,
            embedding_name=embedding_name,
        )
    else:
        print("Using normalized counts for Oracle object creation.")
        # use normalized counts to build the oracle object
        oracle.import_anndata_as_normalized_count(
            adata=adata_cc,
            cluster_column_name=cluster_column_name,
            embedding_name=embedding_name,
        )

    # Load base GRN based on species (skip when --no-base-grn is set)
    base_GRN = None
    if no_base_grn:
        print(
            "⊘ Skipping base GRN loading (--no-base-grn). "
            "TF information must be supplied via TG_to_TF_dictionary or ATAC peaks."
        )
    elif species == "human":
        base_GRN = co.data.load_human_promoter_base_GRN()
    elif species == "mouse":
        base_GRN = co.data.load_mouse_promoter_base_GRN()
    else:
        print("species is not human or mouse; no base GRN is loaded")

    # Import base GRN (skipped when no_base_grn=True)
    if base_GRN is not None:
        print(f"Importing base GRN for species: {species} from CellOracle data module")
        oracle.import_TF_data(TF_info_matrix=base_GRN)
    elif no_base_grn and TG_to_TF_dictionary is not None:
        print(
            f"Loading TG to TF dictionary from: {TG_to_TF_dictionary} and "
            "importing as TF data since --no-base-grn is set"
        )
        TG_to_TF_dictionary_open = pickle.load(
            open(
                TG_to_TF_dictionary,
                "rb",
            )
        )
        oracle.import_TF_data(TFdict=TG_to_TF_dictionary_open)
    elif no_base_grn and TG_to_TF_dictionary is None:
        print(
            "⚠ Warning: --no-base-grn is set but no TG_to_TF_dictionary "
            "was provided. The Oracle object has no TF data."
        )

    if TG_to_TF_dictionary is not None and base_GRN is not None:
        print(
            f"Loading TG to TF dictionary from: {TG_to_TF_dictionary} to enhance base GRN"
        )
        # Load the TG to TF dictionary
        TG_to_TF_dictionary = pickle.load(
            open(
                TG_to_TF_dictionary,
                "rb",
            )
        )
        # Add the TG to TF dictionary to the oracle object
        oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

    return oracle


def run_PCA(oracle: co.Oracle):
    """
    Perform Principal Component Analysis (PCA) on a CellOracle object and determine optimal number of components.
    This function creates a copy of the input Oracle object, performs PCA analysis, and automatically
    selects the optimal number of principal components based on the explained variance ratio. It visualizes
    the cumulative explained variance and identifies the elbow point to determine where the rate of
    variance explanation significantly decreases.
    Args:
        oracle (co.Oracle): A CellOracle Oracle object containing gene expression data
                           to be analyzed with PCA.
    Returns:
        co.Oracle: A copy of the input Oracle object with PCA analysis performed.
                   The object will have PCA results stored and the optimal number of
                   components determined (capped at maximum 50 components).
    Notes:
        - The function automatically determines the optimal number of components by finding
          the elbow point in the cumulative explained variance curve
        - A plot showing cumulative explained variance is displayed with a vertical line
          indicating the selected number of components
        - The number of components is limited to a maximum of 50
        - The original Oracle object is preserved; a copy is returned with PCA results
    """

    # Copy the oracle object to preserve the original data
    oracle_cc = oracle.copy()

    # Perform PCA
    oracle_cc.perform_PCA()

    # Select important PCs
    fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_SMALL)
    ax.plot(np.cumsum(oracle_cc.pca.explained_variance_ratio_)[:200])
    n_comps = np.where(
        np.diff(np.diff(np.cumsum(oracle_cc.pca.explained_variance_ratio_)) > 0.002)
    )[0][0]
    ax.axvline(n_comps, c="k")
    fig.savefig(config.FIGURES_DIR_GRN + "/pca_explained_variance.png")
    plt.close("all")

    n_comps = min(n_comps, 50)

    return oracle_cc, n_comps


def run_KNN(oracle: co.Oracle, n_comps: int = 50):
    """
    Perform K-Nearest Neighbors (KNN) imputation on a CellOracle object.

    This function creates a copy of the input Oracle object and applies KNN imputation
    to smooth gene expression data based on cellular similarity. The number of neighbors
    is automatically calculated as 2.5% of the total cell count.

    Parameters
    ----------
    oracle : co.Oracle
        A CellOracle Oracle object containing single-cell data and analysis results.
    n_comps : int, optional
        Number of principal components to use for KNN calculation. Default is 50.

    Returns
    -------
    None
        The function modifies the oracle_cc object in-place but does not return it.

    Notes
    -----
    - The function automatically calculates k as 2.5% of the total number of cells
    - The b_sight parameter is set to k * 8
    - The b_maxl parameter is set to k * 4
    - The imputation uses balanced KNN with parallel processing
    - The number of parallel jobs is controlled by config.GRN_N_JOBS

    Examples
    --------
    >>> import celloracle as co
    >>> oracle = co.Oracle()
    >>> # Load and prepare your data
    >>> run_KNN(oracle, n_comps=30)
    cell number is :5000
    Auto-selected k is :125
    """

    oracle_cc = oracle.copy()

    n_cell = oracle_cc.adata.shape[0]
    print(f"cell number is :{n_cell}")
    k = int(0.025 * n_cell)
    print(f"Auto-selected k is :{k}")

    oracle_cc.knn_imputation(
        n_pca_dims=n_comps,
        k=k,
        balanced=True,
        b_sight=k * 8,
        b_maxl=k * 4,
        n_jobs=config.GRN_N_JOBS,
    )

    return oracle_cc


def run_links(
    oracle: co.Oracle,
    cluster_column_name: str = "leiden",
    p_cutoff: float = 0.001,
    **kwargs,
):
    """
    Calculate and filter gene regulatory network (GRN) links using CellOracle.
    This function computes regulatory links between transcription factors and target genes
    for each cluster in the dataset, filters them based on statistical significance, and
    generates diagnostic plots and network scores.
    Parameters
    ----------
    oracle : co.Oracle
        A CellOracle Oracle object containing the preprocessed single-cell data and
        base GRN information.
    cluster_column_name : str
        Name of the column in the oracle object that contains cluster assignments
        for cells. Used to compute cluster-specific GRN units.
    p_cutoff : float, optional
        P-value threshold for filtering significant regulatory links. Links with
        p-values above this cutoff will be removed. Default is 0.001.
    Returns
    -------
    links : co.Links
        A CellOracle Links object containing the filtered gene regulatory network
        with methods for further analysis and visualization.
    Notes
    -----
    - The function uses an alpha value of 10 for ridge regression regularization
    - Filtering is based on the absolute value of coefficients ('coef_abs')
    - A degree distribution plot is automatically saved to the figures directory
    - Network scores are computed to assess the quality of the inferred GRN
    """
    cluster_column_name = kwargs.get(
        "cluster_key", kwargs.get("cluster_column", cluster_column_name)
    )
    cluster_column_name = resolve_cluster_key_name(cluster_column_name)
    print(f"Calculating GRN links for cluster column: {cluster_column_name}")

    # Calculate GRN links
    links = oracle.get_links(
        cluster_name_for_GRN_unit=cluster_column_name,
        alpha=10,
        verbose_level=10,
        n_jobs=config.GRN_N_JOBS,
    )
    # Filter links based on p-value cutoff
    links.filter_links(p=p_cutoff, weight="coef_abs")
    # Calculate network scores like centrality, etc
    links.get_network_score()
    # Plot some stats over the network
    links.plot_degree_distributions(
        plot_model=True, save=config.FIGURES_DIR_GRN + "/grn_degree_distribution/"
    )

    # Plot rank scores for each cluster unit
    if hasattr(links, "cluster") and links.cluster:
        for c in links.cluster:
            try:
                links.plot_scores_as_rank(
                    cluster=str(c),
                    n_gene=10,
                    save=config.FIGURES_DIR_GRN + "/grn_ranks/",
                )
            except Exception as e:
                print(f"  ⚠ Could not plot rank score for cluster '{c}': {e}")

    return links


def load_hotspot_genes(hotspot_genes_path: str) -> list:
    """
    Load significant autocorrelated genes identified by Hotspot.

    Parameters
    ----------
    hotspot_genes_path : str
        Path to the ``significant_genes.csv`` file produced by Hotspot
        (typically ``<output_dir>/hotspot/significant_genes.csv``).

    Returns
    -------
    list
        Gene names of Hotspot-significant autocorrelated genes.

    Examples
    --------
    >>> from genecircuitry.celloracle_processing import load_hotspot_genes
    >>> genes = load_hotspot_genes("output/hotspot/significant_genes.csv")
    >>> print(f"Loaded {len(genes)} Hotspot genes")
    """
    df = pd.read_csv(hotspot_genes_path, index_col=0)
    genes = df.index.tolist()
    print(f"  Loaded {len(genes)} Hotspot significant genes from: {hotspot_genes_path}")
    return genes


def perform_grn_pre_processing(
    adata: AnnData,
    cluster_key: Union[str, Sequence[str]] = "leiden",
    cell_downsample: int = 20000,
    cell_downsample: Optional[int] = None,
    top_genes: Optional[int] = None,
    gene_list: Optional[list] = None,
    n_neighbors: Optional[int] = None,
    n_pcs: Optional[int] = None,
    svd_solver: Optional[str] = None,
    **kwargs,
) -> AnnData:
    """
    Perform preprocessing for GRN analysis on AnnData object.

    This function identifies highly variable genes (or uses a custom gene list),
    computes PCA and diffusion maps, in preparation for gene regulatory network
    analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations and genes as variables.
    cluster_key : str or sequence of str, default="leiden"
        Key(s) in adata.obs to use for clustering over GRN. Can be a single column,
        a comma-separated string, or a sequence/list of columns.
    cell_downsample : int, default=20000
        Number of cells to downsample to (default: 20000).
    cell_downsample : int, optional
        Number of cells to downsample to.
        If None, uses config.GRN_CELL_DOWNSAMPLE (default: 20000).
        Set to None, 0, or negative to disable downsampling.
    top_genes : int, optional
        Number of highly variable genes to select.
        If None, uses config.HVGS_N_TOP_GENES (default: 2000).
        Ignored when ``gene_list`` is provided.
    gene_list : list, optional
        Custom list of genes to use instead of HVG selection.
        When provided, overrides ``top_genes``.  Useful for passing Hotspot
        autocorrelated genes via :func:`load_hotspot_genes`.
    n_neighbors : int, optional
        Number of neighbors for KNN graph.
        If None, uses config.NEIGHBORS_N_NEIGHBORS (default: 15).
    n_pcs : int, optional
        Number of principal components to use.
        If None, uses config.NEIGHBORS_N_PCS (default: 40).
    svd_solver : str, optional
        SVD solver for PCA.
        If None, uses config.PCA_SVD_SOLVER (default: 'arpack').

    Returns
    -------
    adata : AnnData
        Preprocessed AnnData object ready for GRN analysis.

    Examples
    --------
    >>> import scanpy as sc
    >>> from genecircuitry.celloracle_processing import perform_grn_pre_processing
    >>> from genecircuitry import config
    >>>
    >>> # Use default config values (HVG selection)
    >>> adata = sc.read_h5ad('data.h5ad')
    >>> adata_preprocessed = perform_grn_pre_processing(adata, cluster_key='louvain')
    >>>
    >>> # Override specific parameters
    >>> adata_preprocessed = perform_grn_pre_processing(adata, top_genes=5000, n_pcs=50)
    >>>
    >>> # Use Hotspot genes instead of HVGs
    >>> from genecircuitry.celloracle_processing import load_hotspot_genes
    >>> hotspot_genes = load_hotspot_genes("output/hotspot/significant_genes.csv")
    >>> adata_preprocessed = perform_grn_pre_processing(
    ...     adata, cluster_key='louvain', gene_list=hotspot_genes
    ... )
    """
    cluster_key = kwargs.get(
        "cluster_column_name", kwargs.get("cluster_column", cluster_key)
    )
    if cluster_key is None:
        cluster_key = "leiden"

    # Use config defaults if not specified
    if cell_downsample is None:
        cell_downsample = config.GRN_CELL_DOWNSAMPLE
    if top_genes is None:
        top_genes = config.HVGS_N_TOP_GENES
    if n_neighbors is None:
        n_neighbors = config.NEIGHBORS_N_NEIGHBORS
    if n_pcs is None:
        n_pcs = config.NEIGHBORS_N_PCS
    if svd_solver is None:
        svd_solver = config.PCA_SVD_SOLVER

    # Make a copy to avoid modifying the original
    adata_cc = adata.copy()

    # Resolve cluster_key (handles single key, comma-separated keys, or sequence of keys)
    adata_cc, cluster_key = resolve_cluster_key(
        adata_cc, cluster_key, key_term="Cluster key"
    )

    if gene_list is not None:
        # Use provided gene list (e.g., Hotspot autocorrelated genes)
        valid_genes = [g for g in gene_list if g in adata_cc.var_names]
        if len(valid_genes) == 0:
            raise ValueError(
                "None of the genes in gene_list are present in adata.var_names. "
                "Check that gene names match."
            )
        if len(valid_genes) < len(gene_list):
            print(
                f"  Warning: {len(gene_list) - len(valid_genes)} gene(s) from "
                f"gene_list not found in data and were excluded."
            )
        adata_cc = adata_cc[:, valid_genes]
        print(f"  Subsetting to {len(valid_genes)} provided genes (gene_list)")
    else:
        # Standard HVG selection
        sc.pp.highly_variable_genes(adata_cc, n_top_genes=top_genes, subset=True)
        print(f"  Selected {top_genes} highly variable genes")

    # PCA
    sc.tl.pca(adata_cc, svd_solver=svd_solver)
    print(f"  Computed PCA with {svd_solver} solver")

    # Diffusion map
    sc.pp.neighbors(adata_cc, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.diffmap(adata_cc)
    print(f"  Computed diffusion map with {n_neighbors} neighbors and {n_pcs} PCs")

    # Calculate neighbors again based on diffusion map
    sc.pp.neighbors(adata_cc, n_neighbors=n_neighbors, use_rep="X_diffmap")
    print(f"  Recomputed neighbors using diffusion map representation")

    # Clustering
    sc.tl.paga(adata_cc, groups=cluster_key)
    sc.pl.paga(adata_cc, show=None, save=None)
    print(f"  Computed PAGA for cluster key: {cluster_key}")

    sc.tl.draw_graph(adata_cc, init_pos="paga")
    sc.pl.draw_graph(adata_cc, color=cluster_key, show=None, save=None)
    print(f"  Computed draw_graph for cluster key: {cluster_key}")

    if adata_cc.n_obs > cell_downsample:
        if cell_downsample is not None and cell_downsample > 0 and adata_cc.n_obs > cell_downsample:
            print(f"  Downsampling to {cell_downsample} cells for GRN analysis")
            sc.pp.subsample(adata_cc, n_obs=cell_downsample)

    return adata_cc
