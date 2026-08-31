"""
Hotspot processing module for GeneCircuitry

This module provides functions for Hotspot analysis of gene expression data.
Plotting functions have been consolidated in genecircuitry.plotting.hotspot_plots;
the names below are kept for backward compatibility and delegate to that module.
"""

from typing import Optional, Union, Sequence
import os
import scanpy as sc
import hotspot as hs
from scipy.sparse import csc_matrix
from anndata import AnnData
import pickle
import pandas as pd

from . import config
from .logging_utils import log_error, log_warning

from . import enrichment_analysis as ea


# ---------------------------------------------------------------------------
# Analysis helpers (non-plotting)
# ---------------------------------------------------------------------------


def save_hotspot_results(
    hotspot_obj: hs.Hotspot,
):
    """
    Save Hotspot analysis results to specified output directory.

    Parameters:
        hotspot_obj (hs.Hotspot): An instance of the Hotspot class containing analysis results.

    """
    # Get results summary
    autocorr_results = hotspot_obj.results

    significant_genes = autocorr_results[
        autocorr_results.FDR < config.HOTSPOT_FDR_THRESHOLD
    ]

    module_scores = hotspot_obj.module_scores
    module_scores.index.name = "cell_id"

    gene_modules = hotspot_obj.modules

    # Save additional results
    results_path = f"{config.OUTPUT_DIR}/hotspot/autocorrelation_results.csv"
    autocorr_results.to_csv(results_path)
    print(f"\n✓ Autocorrelation results saved to: {results_path}")

    significant_path = f"{config.OUTPUT_DIR}/hotspot/significant_genes.csv"
    significant_genes.to_csv(significant_path)
    print(f"✓ Significant genes saved to: {significant_path}")

    module_scores_path = f"{config.OUTPUT_DIR}/hotspot/hotspot_module_scores.csv"
    module_scores.to_csv(module_scores_path)
    print(f"✓ Module scores saved to: {module_scores_path}")

    modules_path = f"{config.OUTPUT_DIR}/hotspot/gene_modules.csv"
    gene_modules.to_csv(modules_path)
    print(f"✓ Gene modules saved to: {modules_path}")

    hotspot_object_path = f"{config.OUTPUT_DIR}/hotspot/hotspot_object.pkl"
    with open(hotspot_object_path, "wb") as f:
        pickle.dump(hotspot_obj, f)
    print(f"✓ Hotspot object saved to: {hotspot_object_path}")


def create_hotspot_object(
    adata: AnnData,
    top_genes: int = config.HOTSPOT_TOP_GENES,
    layer_key: str = "raw_counts",
    model: str = "danb",
    embedding_key: str = "X_pca",
    normalization_key: str = "total_counts",
):
    """
    Creates a Hotspot object for spatial gene module analysis.

    Parameters:
        adata (AnnData): Annotated data object containing gene expression data.
        top_genes (int): Number of top highly variable genes to select. If None, use all genes. Default is 3000.
        layer_key (str): Name of the layer in `adata.layers` that contains the expression data to be used.
        model (str): Statistical model to use for Hotspot analysis. Default is "danb (depth-adjusted negative binomial model)".
        embedding_key (str): Name of the embedding in `adata.obsm` to be used for spatial analysis. Default is "X_pca".
        normalization_key (str): Key in `adata.obs` for normalization. Default is "total_counts".

    Returns:
        hotspot_obj (hs.Hotspot): An instance of the Hotspot class.

    """

    # Create a copy of the AnnData object to avoid modifying the original
    adata_cc = adata.copy()

    if top_genes:
        print(f"Selecting top {top_genes} highly variable genes for Hotspot analysis.")
        sc.pp.highly_variable_genes(adata_cc, n_top_genes=top_genes, subset=True)
    else:
        print("Using all genes for Hotspot analysis.")

    # create a csv matrix from the specified layer or from adata.X
    if layer_key:
        print(f"Using layer {layer_key} for Hotspot analysis.")
        adata_cc.layers[f"csc_{layer_key}"] = csc_matrix(adata_cc.layers[layer_key])
    else:
        print("Using adata.X for Hotspot analysis.")
        adata_cc.layers["csc_X"] = csc_matrix(adata_cc.X)

    # Create Hotspot object
    hotspot_obj = hs.Hotspot(
        adata_cc,
        layer_key=f"csc_{layer_key}" if layer_key else "csc_X",
        model=model,
        latent_obsm_key=embedding_key,
        umi_counts_obs_key=normalization_key,
    )

    return hotspot_obj


def _get_module_enrichment_labels(
    hotspot_obj: hs.Hotspot,
    gene_sets: list = None,
    max_label_length: int = 30,
) -> dict:
    """
    Get enrichment-based labels for each module.

    Parameters:
        hotspot_obj: Hotspot object with modules.
        gene_sets: Gene sets for enrichment analysis.
        max_label_length: Maximum length of label text.

    Returns:
        dict: Mapping from module number to enrichment label.
    """
    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)
    module_labels = {}

    # First try to load existing enrichment results
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
            log_error("Hotspot.LoadEnrichmentFile", e)
            print(
                f"  Warning: Could not load enrichment file "
                f"({type(e).__name__}): {e}"
            )

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
                f"Hotspot.ModuleEnrichment(module={module})",
                f"Enrichment failed ({type(e).__name__}): {e}",
            )
            module_labels[module] = f"Module {module}"

    return module_labels


def run_hotspot_analysis(
    hotspot_obj,
    adata: Optional[AnnData] = None,
    cluster_key: Union[str, Sequence[str]] = "leiden",
    n_jobs: Optional[int] = None,
):
    """
    Run Hotspot analysis on the given Hotspot object.

    Parameters:
        hotspot_obj: An instance of the Hotspot class.
        adata: Optional AnnData object with cluster annotations for violin plots.
        cluster_key: Column name(s) in adata.obs containing cluster assignments.
        n_jobs: Number of parallel worker processes.

    Returns:
        hotspot_obj: The updated Hotspot object with analysis results.
    """
    from .plotting.hotspot_plots import generate_all_hotspot_plots
    from .preprocessing import resolve_cluster_key, resolve_cluster_key_name

    if adata is not None and cluster_key:
        try:
            adata, cluster_key = resolve_cluster_key(adata, cluster_key)
        except Exception:
            cluster_key = resolve_cluster_key_name(cluster_key)

    # Create KNN graph
    hotspot_obj.create_knn_graph(
        weighted_graph=False, n_neighbors=config.HOTSPOT_N_NEIGHBORS
    )
    print("  KNN graph created successfully")

    # Compute autocorrelations
    hs_results = hotspot_obj.compute_autocorrelations(jobs=config.HOTSPOT_N_JOBS)
    print("  Autocorrelations computed successfully")

    # Identify significant genes
    hs_genes = hs_results.loc[
        hs_results.FDR < config.HOTSPOT_FDR_THRESHOLD
    ].index  # Select genes
    local_correlations = hotspot_obj.compute_local_correlations(
        hs_genes, jobs=config.HOTSPOT_N_JOBS
    )
    print(
        f"  Identified {len(hs_genes)} significant genes (FDR < {config.HOTSPOT_FDR_THRESHOLD})"
    )

    modules = hotspot_obj.create_modules(
        min_gene_threshold=config.HOTSPOT_MIN_GENES_PER_MODULE,
        core_only=config.HOTSPOT_CORE_ONLY,
        fdr_threshold=config.HOTSPOT_FDR_THRESHOLD,
    )
    print(f"  Identified {len(modules.unique())} gene modules")

    module_scores = hotspot_obj.calculate_module_scores()
    print("  Module scores calculated")

    # save hotspot results before plotting to ensure results are saved even if plotting fails
    save_hotspot_results(hotspot_obj)

    # Generate all Hotspot plots using the canonical plotting module
    print("\n  Generating Hotspot visualizations...")
    plot_results = generate_all_hotspot_plots(
        hotspot_obj,
        adata=adata,
        cluster_key=cluster_key,
        skip_existing=False,  # Always regenerate for fresh analysis
        n_jobs=n_jobs,
    )

    # Report which plots were generated
    generated = [k for k, v in plot_results.items() if v is True]
    if generated:
        print(f"  Generated plots: {', '.join(generated)}")

    return hotspot_obj


# ---------------------------------------------------------------------------
# Plotting helpers — delegating to genecircuitry.plotting.hotspot_plots
# ---------------------------------------------------------------------------


def plot_hotspot_annotation(
    hs_obj: hs.Hotspot,
    gene_sets: list = None,
    top_n_annotations: int = 1,
    n_jobs: Optional[int] = None,
):
    """
    Plot Hotspot gene module annotations with enrichment analysis results.

    Delegates to genecircuitry.plotting.hotspot_plots.plot_hotspot_annotation.
    """
    from .plotting.hotspot_plots import plot_hotspot_annotation as _impl

    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)
    return _impl(
        hotspot_obj=hs_obj,
        gene_sets=gene_sets,
        top_n_annotations=top_n_annotations,
        n_jobs=n_jobs,
    )


def plot_module_scores_violin(
    hotspot_obj: hs.Hotspot,
    adata: AnnData,
    cluster_key: str = "leiden",
    figsize_per_cluster: tuple = (12, 6),
    gene_sets: list = None,
    n_jobs: Optional[int] = None,
):
    """
    Plot violin plots of module scores for each cluster.

    Delegates to genecircuitry.plotting.hotspot_plots.plot_module_scores_violin.
    """
    from .plotting.hotspot_plots import plot_module_scores_violin as _impl

    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)
    return _impl(
        hotspot_obj=hotspot_obj,
        adata=adata,
        cluster_key=cluster_key,
        figsize_per_cluster=figsize_per_cluster,
        gene_sets=gene_sets,
        n_jobs=n_jobs,
    )
