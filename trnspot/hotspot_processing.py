"""
Hotspot processing module for TRNspot

This module provides functions for Hotspot analysis of gene expression data.

NOTE: Plotting functions have been moved to trnspot.plotting.hotspot_plots.
The functions here are kept for backward compatibility but delegate to the new module.
"""

from typing import Optional
import os
import scanpy as sc
import hotspot as hs
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
from anndata import AnnData
import pickle
import pandas as pd
import seaborn as sns
import marsilea as ma

from . import config

from . import enrichment_analysis as ea
from matplotlib.patches import Patch

# Import from new plotting module for backward compatibility
from .plotting import plot_exists as _plot_exists
from .plotting import save_plot
from .plotting import (
    plot_module_scores_violin as _new_plot_module_scores_violin,
    generate_all_hotspot_plots,
)


def plot_hotspot_annotation(
    hs_obj: hs.Hotspot,
    gene_sets: list = ["MSigDB_Hallmark_2020"],
    top_n_annotations: int = 1,
):
    """
    Plot Hotspot gene module annotations with enrichment analysis results.

    Parameters:
        hs_obj (hs.Hotspot): An instance of the Hotspot class containing analysis results.
        gene_sets (list): Gene sets to use for enrichment analysis.
            Default is ["MSigDB_Hallmark_2020"]. Other options include:
            - "Reactome_Pathways_2024"
            - "KEGG_2021_Human"
            - "GO_Biological_Process_2023"
            - "GO_Molecular_Function_2023"
        top_n_annotations (int): Number of top annotations to display per module.

    """

    print(
        f"Generating Hotspot local correlation heatmap with {gene_sets} annotations..."
    )

    # Perform enrichment analysis for each module
    df_enrichment = pd.DataFrame()
    for module in hs_obj.modules.unique():
        if module == -1:  # Skip unassigned genes
            continue
        genes = hs_obj.modules[hs_obj.modules == module].index.tolist()
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
            print(f"  Warning: Enrichment analysis failed for module {module}: {e}")
            continue

    # Save enrichment results
    if not df_enrichment.empty:
        df_enrichment.to_csv(
            f"{config.OUTPUT_DIR}/hotspot/hotspot_module_enrichment_results.csv",
            index=False,
        )

    # Create module colors using a large, distinct color palette
    n_modules_total = len([m for m in hs_obj.modules.unique() if m != -1])
    if n_modules_total <= 10:
        colors = sns.color_palette("tab10", n_colors=10)
    elif n_modules_total <= 20:
        colors = sns.color_palette("tab20", n_colors=20)
    else:
        # For many modules, use husl which generates evenly spaced hues
        colors = sns.color_palette("husl", n_colors=max(n_modules_total, 30))

    module_colors = {i: colors[(i - 1) % len(colors)] for i in hs_obj.modules.unique()}
    module_colors[-1] = "#ffffff"

    row_colors1 = pd.Series(
        [module_colors[i] for i in hs_obj.modules],
        index=hs_obj.local_correlation_z.index,
    )

    # Get top annotations for each module from enrichment results
    module_annotations = {}
    for module in hs_obj.modules.unique():
        if module == -1:  # Skip unassigned genes
            module_annotations[module] = ""
            continue
        if not df_enrichment.empty and module in df_enrichment["module"].values:
            module_df = df_enrichment[df_enrichment["module"] == module]
            if not module_df.empty:
                # Sort by Combined_Score or Adjusted_P-value
                if "Combined_Score" in module_df.columns:
                    top_terms = module_df.nlargest(top_n_annotations, "Combined_Score")[
                        "Term"
                    ].tolist()
                else:
                    top_terms = module_df.nsmallest(
                        top_n_annotations, "Adjusted_P-value"
                    )["Term"].tolist()
                # Clean up term names (remove prefix like "HALLMARK_")
                top_terms = [
                    t.replace("HALLMARK_", "").replace("_", " ").title()
                    for t in top_terms
                ]
                module_annotations[module] = "; ".join(top_terms[:top_n_annotations])
            else:
                module_annotations[module] = "No enrichment"
        else:
            module_annotations[module] = "No enrichment"

    # Create annotation labels for genes based on their module's enrichment
    gene_annotations = pd.Series(
        [module_annotations.get(i, "") for i in hs_obj.modules],
        index=hs_obj.local_correlation_z.index,
        name="Top_Annotations",
    )

    # Create row colors DataFrame
    row_colors = pd.DataFrame(
        {
            "Modules": row_colors1,
        }
    )

    # Create the clustermap
    g = sns.clustermap(
        hs_obj.local_correlation_z,
        row_linkage=hs_obj.linkage,
        col_linkage=hs_obj.linkage,
        row_colors=row_colors["Modules"],
        cmap="RdBu_r",
        vmin=-8,
        vmax=8,
        xticklabels=False,
        yticklabels=False,
        rasterized=True,
        figsize=(12, 12),
    )

    # Remove axis dendrogram
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    # Create legend elements with module annotations
    legend_elements = []
    for module, color in sorted(module_colors.items()):
        if module == -1:
            continue
        annotation = module_annotations.get(module, "")
        if annotation and annotation != "No enrichment":
            # Truncate long annotations
            if len(annotation) > 50:
                annotation = annotation[:47] + "..."
            label = f"M{module}: {annotation}"
        else:
            label = f"Module {module}"
        legend_elements.append(Patch(facecolor=color, edgecolor="k", label=label))

    # Add the legend to the figure
    g.ax_heatmap.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        title="Modules (Hallmark Enrichment)",
        frameon=False,
        fontsize=8,
    )

    # Save the figure
    plt.savefig(
        f"{config.FIGURES_DIR_HOTSPOT}/hotspot_local_correlation_heatmap_with_annotations.png",
        dpi=config.SAVE_DPI,
        bbox_inches="tight",
    )
    plt.close()
    print(
        f"  ✓ Saved annotated heatmap to: {config.FIGURES_DIR_HOTSPOT}/hotspot_local_correlation_heatmap_with_annotations.png"
    )


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
    gene_sets: list = ["MSigDB_Hallmark_2020"],
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
            print(f"  Warning: Could not load enrichment file: {e}")

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
        except Exception:
            module_labels[module] = f"Module {module}"

    return module_labels


def run_hotspot_analysis(
    hotspot_obj, adata: Optional[AnnData] = None, cluster_key: str = "leiden"
):
    """
    Run Hotspot analysis on the given Hotspot object.

    Parameters:
        hotspot_obj: An instance of the Hotspot class.
        adata: Optional AnnData object with cluster annotations for violin plots.
        cluster_key: Column name in adata.obs containing cluster assignments.

    Returns:
        hotspot_obj: The updated Hotspot object with analysis results.
    """

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

    # Generate all Hotspot plots using the new plotting module
    print("\n  Generating Hotspot visualizations...")
    plot_results = generate_all_hotspot_plots(
        hotspot_obj,
        adata=adata,
        cluster_key=cluster_key,
        skip_existing=False,  # Always regenerate for fresh analysis
    )

    # Report which plots were generated
    generated = [k for k, v in plot_results.items() if v is True]
    if generated:
        print(f"  Generated plots: {', '.join(generated)}")

    save_hotspot_results(hotspot_obj)

    return hotspot_obj


def plot_module_scores_violin(
    hotspot_obj: hs.Hotspot,
    adata: AnnData,
    cluster_key: str = "leiden",
    figsize_per_cluster: tuple = (12, 6),
    gene_sets: list = ["MSigDB_Hallmark_2020"],
):
    """
    Plot violin plots of module scores for each cluster.

    For each cluster, creates a violin plot with modules on x-axis
    and module scores on y-axis. X-axis labels show enrichment annotations.

    Parameters:
        hotspot_obj: Hotspot object with computed module scores.
        adata: AnnData object with cluster annotations.
        cluster_key: Column name in adata.obs containing cluster assignments.
        figsize_per_cluster: Figure size for each cluster plot.
        gene_sets: Gene sets for enrichment annotation labels.
    """
    print("  Generating module score violin plots per cluster...")

    # Get module scores
    module_scores = hotspot_obj.module_scores

    if module_scores is None or module_scores.empty:
        print("  Warning: No module scores available for violin plots")
        return

    # Get clusters from adata
    if cluster_key not in adata.obs.columns:
        print(f"  Warning: Cluster key '{cluster_key}' not found in adata.obs")
        return

    # Align cell indices
    common_cells = module_scores.index.intersection(adata.obs.index)
    if len(common_cells) == 0:
        print("  Warning: No common cells between module scores and adata")
        return

    # Get module enrichment annotations
    module_labels = _get_module_enrichment_labels(hotspot_obj, gene_sets)

    # Create combined dataframe
    scores_df = module_scores.loc[common_cells].copy()
    scores_df["cluster"] = adata.obs.loc[common_cells, cluster_key].values

    # Get unique modules (excluding -1 if present)
    modules = [col for col in module_scores.columns if col != -1]
    if not modules:
        print("  Warning: No valid modules found for violin plots")
        return

    # Create a mapping from module number to label
    module_to_label = {m: module_labels.get(m, f"Module {m}") for m in modules}

    # Rename columns in scores_df to use enrichment labels
    rename_dict = {m: module_to_label[m] for m in modules}
    scores_df_labeled = scores_df.rename(columns=rename_dict)
    labeled_modules = [module_to_label[m] for m in modules]

    # Melt dataframe for seaborn
    scores_melted = scores_df_labeled.melt(
        id_vars=["cluster"],
        value_vars=labeled_modules,
        var_name="Module",
        value_name="Score",
    )

    # Get unique clusters
    clusters = scores_df["cluster"].unique()

    # Plot 1: All clusters combined in one plot with facets
    n_clusters = len(clusters)
    n_cols = min(3, n_clusters)
    n_rows = (n_clusters + n_cols - 1) // n_cols

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(
            figsize_per_cluster[0] * n_cols / 2,
            figsize_per_cluster[1] * n_rows / 2,
        ),
        squeeze=False,
    )
    axes = axes.flatten()

    # Color palette for modules
    module_palette = sns.color_palette("husl", n_colors=len(modules))

    for idx, cluster in enumerate(sorted(clusters)):
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
            )
            ax.set_title(f"Cluster: {cluster}", fontsize=12, fontweight="bold")
            ax.set_xlabel("Module (Enrichment)", fontsize=10)
            ax.set_ylabel("Module Score", fontsize=10)
            ax.tick_params(axis="x", rotation=60, labelsize=8)
            # Ensure x-axis labels are fully visible
            for label in ax.get_xticklabels():
                label.set_ha("right")
        else:
            ax.set_visible(False)

    # Hide unused axes
    for idx in range(len(clusters), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle(
        "Hotspot Module Scores per Cluster", fontsize=14, fontweight="bold", y=1.02
    )
    plt.tight_layout()
    plt.savefig(
        f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_per_cluster.png",
        dpi=config.SAVE_DPI,
        bbox_inches="tight",
    )
    plt.close()
    print(
        f"  ✓ Saved violin plots to: {config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_per_cluster.png"
    )

    # Plot 2: All clusters in one violin plot (modules on x, clusters as hue)
    fig, ax = plt.subplots(figsize=(max(14, len(modules) * 2.5), 8))

    # Create a palette for clusters
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
    )

    ax.set_title("Hotspot Module Scores by Cluster", fontsize=14, fontweight="bold")
    ax.set_xlabel("Module (Enrichment)", fontsize=12)
    ax.set_ylabel("Module Score", fontsize=12)
    ax.tick_params(axis="x", rotation=60, labelsize=9)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.legend(title="Cluster", bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()
    plt.savefig(
        f"{config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_all_clusters.png",
        dpi=config.SAVE_DPI,
        bbox_inches="tight",
    )
    plt.close()
    print(
        f"  ✓ Saved combined violin plot to: {config.FIGURES_DIR_HOTSPOT}/hotspot_module_scores_violin_all_clusters.png"
    )
