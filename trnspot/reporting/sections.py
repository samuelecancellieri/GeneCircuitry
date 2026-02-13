"""
Report Section Builders for TRNspot
===================================

Functions to create individual report sections with detailed information.
"""

import os
import glob
import re
from datetime import datetime
from typing import Optional, Dict, Any, List
from pathlib import Path

from .generator import ReportSection
from .. import config


def _find_figures(
    output_dir: str,
    subdirs: List[str],
    patterns: List[str] = ["*.png", "*.jpg", "*.svg"],
    max_figures: int = 20,
) -> List[str]:
    """Find figure files in specified subdirectories."""
    figures = []
    for subdir in subdirs:
        search_dir = os.path.join(output_dir, subdir)
        if not os.path.exists(search_dir):
            continue
        for pattern in patterns:
            # Search recursively
            figures.extend(
                glob.glob(os.path.join(search_dir, "**", pattern), recursive=True)
            )
    return sorted(set(figures))[:max_figures]


def _find_embedding_figures(
    output_dir: str,
    max_figures: int = 15,
) -> List[str]:
    """
    Find embedding/clustering-specific figures.

    Excludes GRN, hotspot, and other analysis plots.
    """
    embedding_keywords = [
        "umap",
        "pca",
        "cluster",
        "leiden",
        "louvain",
        "neighbors",
        "embedding",
        "tsne",
    ]
    exclude_keywords = [
        "grn",
        "network",
        "centrality",
        "hotspot",
        "module",
        "autocorrelation",
        "heatmap",
        "tf_",
        "regulatory",
    ]

    figures = []
    figures_dir = os.path.join(output_dir, "figures")

    if os.path.exists(figures_dir):
        for pattern in ["*.png", "*.jpg", "*.svg"]:
            # Only search top-level figures directory, not subdirectories
            for f in glob.glob(os.path.join(figures_dir, pattern)):
                fname_lower = os.path.basename(f).lower()
                # Include if matches embedding keywords or is a general plot
                if any(kw in fname_lower for kw in embedding_keywords):
                    # Exclude if matches analysis-specific keywords
                    if not any(ex in fname_lower for ex in exclude_keywords):
                        figures.append(f)

    return sorted(set(figures))[:max_figures]


def _format_number(value) -> str:
    """Format a number for display."""
    if isinstance(value, float):
        if value >= 1000:
            return f"{value:,.0f}"
        elif value >= 1:
            return f"{value:.2f}"
        else:
            return f"{value:.4f}"
    elif isinstance(value, int):
        return f"{value:,}"
    return str(value)


def create_data_summary_section(
    adata,
    output_dir: str,
) -> ReportSection:
    """
    Create Input Data Summary section.

    Shows information about the input dataset including dimensions,
    data types, and available annotations.
    """
    # Basic metrics
    metrics = {
        "Total Cells": adata.n_obs,
        "Total Genes": adata.n_vars,
    }

    # Count layers
    if adata.layers:
        metrics["Data Layers"] = len(adata.layers)

    # Build content
    content_parts = []

    # Dataset dimensions
    content_parts.append(
        """
    <h3>Dataset Overview</h3>
    <p>Summary of the input single-cell RNA-seq dataset.</p>
    """
    )

    # Observation annotations table
    if adata.obs.columns.any():
        obs_rows = []
        for col in adata.obs.columns[:15]:  # Limit to 15 columns
            dtype = str(adata.obs[col].dtype)
            n_unique = adata.obs[col].nunique()
            sample_vals = ", ".join(str(x) for x in adata.obs[col].unique()[:3])
            if len(adata.obs[col].unique()) > 3:
                sample_vals += "..."
            obs_rows.append([col, dtype, n_unique, sample_vals])

        content_parts.append("<h3>Cell Annotations (obs)</h3>")
        content_parts.append(
            f"<p>Found {len(adata.obs.columns)} annotation columns:</p>"
        )

    # Variable annotations
    if adata.var.columns.any():
        var_cols = list(adata.var.columns[:10])
        content_parts.append("<h3>Gene Annotations (var)</h3>")
        content_parts.append(
            f"<p>Found {len(adata.var.columns)} gene annotation columns: <code>{', '.join(var_cols)}</code></p>"
        )

    # Embeddings
    if adata.obsm:
        embeddings = []
        for key in adata.obsm.keys():
            shape = adata.obsm[key].shape
            embeddings.append(f"<code>{key}</code> ({shape[1]} dims)")
        content_parts.append("<h3>Embeddings (obsm)</h3>")
        content_parts.append(f"<p>Available embeddings: {', '.join(embeddings)}</p>")

    # Layers
    if adata.layers:
        layers = [f"<code>{key}</code>" for key in adata.layers.keys()]
        content_parts.append("<h3>Data Layers</h3>")
        content_parts.append(f"<p>Available layers: {', '.join(layers)}</p>")

    # Tables
    tables = []
    if adata.obs.columns.any():
        tables.append(
            {
                "title": "Cell Annotation Columns",
                "headers": ["Column", "Data Type", "Unique Values", "Sample Values"],
                "rows": obs_rows if "obs_rows" in dir() else [],
            }
        )

    return ReportSection(
        title="Input Data Summary",
        section_id="data-summary",
        content="\n".join(content_parts),
        metrics=metrics,
        tables=tables,
    )


def create_settings_section() -> ReportSection:
    """
    Create Configuration Settings section.

    Shows all configuration parameters used in the analysis.
    """
    cfg = config.get_config()

    # Organize settings by category
    categories = {
        "Quality Control": {
            "Min Genes per Cell": cfg.get("QC_MIN_GENES"),
            "Min Counts per Cell": cfg.get("QC_MIN_COUNTS"),
            "Max MT %": cfg.get("QC_PCT_MT_MAX"),
            "Min Cells per Gene": cfg.get("QC_MIN_CELLS"),
        },
        "Preprocessing": {
            "Normalization Target": _format_number(cfg.get("NORMALIZE_TARGET_SUM")),
            "HVG Count": cfg.get("HVGS_N_TOP_GENES"),
            "HVG Min Mean": cfg.get("HVGS_MIN_MEAN"),
            "HVG Max Mean": cfg.get("HVGS_MAX_MEAN"),
        },
        "Dimensionality Reduction": {
            "PCA Components": cfg.get("PCA_N_COMPS"),
            "UMAP Min Distance": cfg.get("UMAP_MIN_DIST"),
            "UMAP Components": cfg.get("UMAP_N_COMPONENTS"),
        },
        "Clustering": {
            "N Neighbors": cfg.get("NEIGHBORS_N_NEIGHBORS"),
            "N PCs Used": cfg.get("NEIGHBORS_N_PCS"),
            "Leiden Resolution": cfg.get("LEIDEN_RESOLUTION"),
        },
        "GRN Analysis": {
            "N Jobs": cfg.get("GRN_N_JOBS"),
            "Min Targets": cfg.get("GRN_MIN_TARGETS"),
            "Confidence Threshold": cfg.get("GRN_CONFIDENCE_THRESHOLD"),
        },
        "Hotspot Analysis": {
            "N Neighbors": cfg.get("HOTSPOT_N_NEIGHBORS", 30),
            "FDR Threshold": cfg.get("HOTSPOT_FDR_THRESHOLD", 0.05),
            "Top Genes": cfg.get("HOTSPOT_TOP_GENES", 3000),
        },
        "Output": {
            "Output Directory": cfg.get("OUTPUT_DIR"),
            "Plot DPI": cfg.get("PLOT_DPI"),
            "Save DPI": cfg.get("SAVE_DPI", 600),
            "Random Seed": cfg.get("RANDOM_SEED"),
        },
    }

    # Build HTML content
    content_parts = ['<div class="settings-grid">']

    for category, settings in categories.items():
        settings_html = ""
        for key, value in settings.items():
            if value is not None:
                settings_html += f"<dt>{key}</dt><dd>{value}</dd>"

        content_parts.append(
            f"""
        <div class="settings-category">
            <h4>{category}</h4>
            <dl>{settings_html}</dl>
        </div>
        """
        )

    content_parts.append("</div>")

    return ReportSection(
        title="Configuration Settings",
        section_id="settings",
        content="\n".join(content_parts),
    )


def create_qc_section(
    adata,
    output_dir: str,
) -> ReportSection:
    """
    Create Quality Control section.

    Shows QC metrics and filtering results.
    """
    metrics = {}

    # Collect QC metrics if available
    if "n_genes_by_counts" in adata.obs.columns:
        metrics["Median Genes/Cell"] = int(adata.obs["n_genes_by_counts"].median())
        metrics["Mean Genes/Cell"] = int(adata.obs["n_genes_by_counts"].mean())

    if "total_counts" in adata.obs.columns:
        metrics["Median Counts/Cell"] = int(adata.obs["total_counts"].median())

    if "pct_counts_mt" in adata.obs.columns:
        metrics["Mean MT %"] = round(adata.obs["pct_counts_mt"].mean(), 2)
        metrics["Max MT %"] = round(adata.obs["pct_counts_mt"].max(), 2)

    metrics["Cells After QC"] = adata.n_obs
    metrics["Genes After QC"] = adata.n_vars

    content = """
    <p>Quality control was performed to filter out low-quality cells and uninformative genes.
    The following metrics were calculated and filtering thresholds applied:</p>

    <h3>QC Metrics Calculated</h3>
    <ul>
        <li><strong>n_genes_by_counts:</strong> Number of genes with at least 1 count per cell</li>
        <li><strong>total_counts:</strong> Total number of counts per cell</li>
        <li><strong>pct_counts_mt:</strong> Percentage of counts from mitochondrial genes</li>
    </ul>

    <h3>Filtering Criteria</h3>
    <ul>
        <li>Removed cells with too few detected genes</li>
        <li>Removed cells with abnormally low/high total counts</li>
        <li>Removed cells with high mitochondrial content (potential dying cells)</li>
        <li>Removed genes expressed in too few cells</li>
    </ul>
    """

    # Find QC figures
    figures = _find_figures(output_dir, ["figures/qc", "qc"])

    return ReportSection(
        title="Quality Control",
        section_id="qc",
        content=content,
        figures=figures,
        metrics=metrics,
    )


def create_preprocessing_section(
    adata,
    output_dir: str,
) -> ReportSection:
    """
    Create Preprocessing section.

    Shows normalization, HVG selection, and scaling steps.
    """
    metrics = {}

    # Check for HVGs
    if "highly_variable" in adata.var.columns:
        n_hvg = int(adata.var["highly_variable"].sum())
        metrics["Highly Variable Genes"] = n_hvg

    # Check for normalization
    if "log1p" in adata.uns:
        metrics["Log Transformed"] = "Yes"

    content = """
    <p>Data preprocessing was performed to normalize and prepare the data for downstream analysis.</p>

    <h3>Preprocessing Steps</h3>
    <ol>
        <li><strong>Normalization:</strong> Library size normalization to account for differences in sequencing depth</li>
        <li><strong>Log transformation:</strong> Log(x + 1) transformation to stabilize variance</li>
        <li><strong>Highly Variable Genes:</strong> Selection of genes with high cell-to-cell variation</li>
        <li><strong>Scaling:</strong> Z-score scaling for PCA computation</li>
    </ol>

    <h3>Technical Details</h3>
    <ul>
        <li>Normalization target sum: 10,000 counts per cell (CPM-like)</li>
        <li>HVG selection method: Seurat v3 (highly_variable_genes with flavor='seurat_v3')</li>
        <li>Raw counts preserved in <code>adata.layers['raw_count']</code> for downstream tools</li>
    </ul>
    """

    return ReportSection(
        title="Preprocessing",
        section_id="preprocessing",
        content=content,
        metrics=metrics,
    )


def create_clustering_section(
    adata,
    output_dir: str,
) -> ReportSection:
    """
    Create Clustering & Visualization section.
    """
    metrics = {}

    # Check for clustering results
    cluster_key = None
    if "leiden" in adata.obs.columns:
        cluster_key = "leiden"
        n_clusters = len(adata.obs["leiden"].unique())
        metrics["Clusters (Leiden)"] = n_clusters
    elif "louvain" in adata.obs.columns:
        cluster_key = "louvain"
        n_clusters = len(adata.obs["louvain"].unique())
        metrics["Clusters (Louvain)"] = n_clusters

    # Check embeddings
    if "X_pca" in adata.obsm:
        metrics["PCA Components"] = adata.obsm["X_pca"].shape[1]

    if "X_umap" in adata.obsm:
        metrics["UMAP Computed"] = "Yes"

    # Cluster sizes table
    tables = []
    if cluster_key:
        cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
        rows = [
            [str(cluster), count, f"{100*count/adata.n_obs:.1f}%"]
            for cluster, count in cluster_counts.items()
        ]
        tables.append(
            {
                "title": f"Cluster Sizes ({cluster_key.title()} Clustering)",
                "headers": ["Cluster", "Cell Count", "Percentage"],
                "rows": rows,
            }
        )

    content = """
    <p>Dimensionality reduction and clustering were performed to identify cell populations.</p>

    <h3>Dimensionality Reduction</h3>
    <ol>
        <li><strong>PCA:</strong> Principal Component Analysis for initial dimensionality reduction</li>
        <li><strong>Neighbor Graph:</strong> K-nearest neighbor graph construction</li>
        <li><strong>UMAP:</strong> Uniform Manifold Approximation and Projection for 2D visualization</li>
    </ol>

    <h3>Clustering</h3>
    <p>Graph-based clustering using the Leiden algorithm was applied to identify cell populations
    with similar transcriptional profiles.</p>
    """

    # Find only clustering/embedding plots, not GRN or other analysis plots
    figures = _find_embedding_figures(output_dir, max_figures=10)

    return ReportSection(
        title="Clustering & Visualization",
        section_id="clustering",
        content=content,
        figures=figures,
        metrics=metrics,
        tables=tables,
    )


def _find_grn_figures_by_score(
    output_dir: str,
    score_type: str,
    max_figures: int = 30,
) -> List[str]:
    """
    Find GRN figures for a specific score type.

    Parameters
    ----------
    output_dir : str
        Output directory
    score_type : str
        Score type to filter by (e.g., 'eigenvector', 'betweenness', 'degree_in')
    max_figures : int
        Maximum figures to return

    Returns
    -------
    list
        List of figure paths matching the score type
    """
    figures = []
    grn_dirs = [
        os.path.join(output_dir, "figures/grn"),
        os.path.join(output_dir, "figures/grn/grn_deep_analysis"),
        os.path.join(output_dir, "grn_deep_analysis"),
    ]

    for grn_dir in grn_dirs:
        if not os.path.exists(grn_dir):
            continue
        for pattern in ["*.png", "*.jpg", "*.svg"]:
            for f in glob.glob(os.path.join(grn_dir, "**", pattern), recursive=True):
                fname_lower = os.path.basename(f).lower()
                # Match score type in filename
                if score_type.lower() in fname_lower:
                    figures.append(f)

    return sorted(set(figures))[:max_figures]


def create_celloracle_section(
    celloracle_result,
    output_dir: str,
) -> ReportSection:
    """
    Create CellOracle GRN Analysis section with centrality score categories only.

    Organizes plots by score type (eigenvector, betweenness, degree centralities)
    and excludes network graphs and other miscellaneous plots.
    """
    metrics = {}

    if celloracle_result is not None:
        oracle, links = celloracle_result

        if hasattr(links, "merged_score") and links.merged_score is not None:
            merged_scores = links.merged_score
            if "cluster" in merged_scores.columns:
                metrics["Clusters Analyzed"] = len(merged_scores["cluster"].unique())

            # Count unique TFs
            if hasattr(merged_scores, "index"):
                metrics["Transcription Factors"] = len(merged_scores.index.unique())

        if hasattr(links, "filtered_links") and links.filtered_links is not None:
            if isinstance(links.filtered_links, dict):
                total_links = sum(
                    len(v)
                    for v in links.filtered_links.values()
                    if hasattr(v, "__len__")
                )
                metrics["Regulatory Links"] = total_links

    content = """
    <p>Gene Regulatory Network (GRN) inference was performed using CellOracle to identify
    transcription factor regulatory relationships from single-cell expression data.</p>

    <h3>Analysis Pipeline</h3>
    <ol>
        <li><strong>Preprocessing:</strong> Diffusion map computation and PAGA graph construction</li>
        <li><strong>Oracle Initialization:</strong> Integration of expression data with base GRN from promoter analysis</li>
        <li><strong>Network Inference:</strong> Ridge regression-based network inference per cluster</li>
        <li><strong>Centrality Analysis:</strong> Calculation of network centrality scores for TF ranking</li>
    </ol>

    <h3>Network Scores Computed</h3>
    <ul>
        <li><strong>Degree centrality:</strong> Number of connections (in/out/all)</li>
        <li><strong>Betweenness centrality:</strong> How often a TF lies on shortest paths</li>
        <li><strong>Eigenvector centrality:</strong> Influence based on connections to influential nodes</li>
    </ul>

    <h3>Output Files</h3>
    <ul>
        <li><code>celloracle/oracle_object.celloracle.oracle</code> - Oracle object</li>
        <li><code>celloracle/oracle_object.celloracle.links</code> - Links object</li>
        <li><code>celloracle/grn_merged_scores.csv</code> - Network centrality scores</li>
        <li><code>celloracle/grn_filtered_links.pkl</code> - Filtered regulatory links</li>
    </ul>
    """

    # Define score categories with centrality scores AND difference plots
    score_categories = [
        ("Eigenvector Centrality", "eigenvector", ["eigenvector"]),
        ("Betweenness Centrality", "betweenness", ["betweenness"]),
        ("Degree Centrality (All)", "degree_centrality_all", ["degree_centrality_all"]),
        (
            "Degree Centrality (In)",
            "degree_centrality_in",
            ["degree_centrality_in", "degree_in"],
        ),
        (
            "Degree Centrality (Out)",
            "degree_centrality_out",
            ["degree_centrality_out", "degree_out"],
        ),
        ("Cluster Comparisons", "comparison", ["compare", "comparison", "scatter"]),
        ("Difference Analysis", "difference", ["difference", "diff_"]),
        ("Network Graphs", "network", ["network", "grn_network"]),
    ]

    # Collect figures for each score category
    subsections = []
    used_figures = set()

    for display_name, section_id, patterns in score_categories:
        category_figures = []
        for pattern in patterns:
            figures_for_pattern = _find_grn_figures_by_score(output_dir, pattern)
            category_figures.extend(figures_for_pattern)

        # Remove duplicates and already used figures
        category_figures = list(
            dict.fromkeys(category_figures)
        )  # preserve order, remove dups
        category_figures = [f for f in category_figures if f not in used_figures]

        if category_figures:
            used_figures.update(category_figures)

            # Create descriptive content based on category
            if "difference" in section_id:
                desc = f"<p>{len(category_figures)} plots showing score differences between clusters for comparative analysis.</p>"
            elif "comparison" in section_id:
                desc = f"<p>{len(category_figures)} scatter plots comparing TF scores across different clusters.</p>"
            elif "network" in section_id:
                desc = f"<p>{len(category_figures)} network visualizations showing regulatory connections.</p>"
            else:
                desc = f"<p>{len(category_figures)} plots showing {display_name.lower()} scores across clusters.</p>"

            subsections.append(
                ReportSection(
                    title=display_name,
                    section_id=f"grn-{section_id.replace('_', '-')}",
                    content=desc,
                    figures=category_figures[:30],  # Increased limit for galleries
                )
            )

    return ReportSection(
        title="CellOracle GRN Analysis",
        section_id="celloracle",
        content=content,
        metrics=metrics,
        subsections=subsections,
    )


def _find_hotspot_figures_by_category(
    output_dir: str,
    category: str,
    max_figures: int = 20,
) -> List[str]:
    """
    Find Hotspot figures for a specific category.
    """
    figures = []
    hotspot_dirs = [
        os.path.join(output_dir, "figures/hotspot"),
        os.path.join(output_dir, "hotspot"),
    ]

    for hs_dir in hotspot_dirs:
        if not os.path.exists(hs_dir):
            continue
        for pattern in ["*.png", "*.jpg", "*.svg"]:
            for f in glob.glob(os.path.join(hs_dir, "**", pattern), recursive=True):
                fname_lower = os.path.basename(f).lower()
                if category.lower() in fname_lower:
                    figures.append(f)

    return sorted(set(figures))[:max_figures]


def create_hotspot_section(
    hotspot_result,
    output_dir: str,
) -> ReportSection:
    """
    Create Hotspot Gene Module section with organized galleries.
    """
    metrics = {}

    if hotspot_result is not None:
        if hasattr(hotspot_result, "results") and hotspot_result.results is not None:
            n_tested = len(hotspot_result.results)
            n_significant = (hotspot_result.results["FDR"] < 0.05).sum()
            metrics["Genes Tested"] = n_tested
            metrics["Significant Genes"] = int(n_significant)

        if hasattr(hotspot_result, "modules") and hotspot_result.modules is not None:
            n_modules = len([m for m in hotspot_result.modules.unique() if m != -1])
            metrics["Gene Modules"] = n_modules

    content = """
    <p>Hotspot analysis was performed to identify groups of genes that show spatially
    coordinated expression patterns across cells.</p>

    <h3>Analysis Steps</h3>
    <ol>
        <li><strong>KNN Graph:</strong> Cell-cell similarity graph construction from PCA space</li>
        <li><strong>Autocorrelation:</strong> Testing genes for local autocorrelation using Geary's C</li>
        <li><strong>Local Correlations:</strong> Computing pairwise gene correlations within local neighborhoods</li>
        <li><strong>Module Detection:</strong> Hierarchical clustering of correlated genes into modules</li>
    </ol>

    <h3>Module Analysis</h3>
    <p>Gene modules represent co-expressed gene programs that may correspond to biological
    processes, cell states, or regulatory programs.</p>

    <h3>Output Files</h3>
    <ul>
        <li><code>hotspot/autocorrelation_results.csv</code> - Gene autocorrelation statistics</li>
        <li><code>hotspot/significant_genes.csv</code> - Significantly autocorrelated genes</li>
        <li><code>hotspot/gene_modules.csv</code> - Gene-to-module assignments</li>
        <li><code>hotspot/hotspot_module_scores.csv</code> - Module scores per cell</li>
    </ul>
    """

    # Define hotspot figure categories with specific patterns for violin plots
    hotspot_categories = [
        ("Local Correlations", "correlation", ["local_correlation"]),
        ("Module Heatmaps", "heatmap", ["heatmap", "annotation"]),
        (
            "Module Scores - Per Cluster",
            "violin-per-cluster",
            ["violin_per_cluster"],
        ),
        (
            "Module Scores - All Clusters",
            "violin-all",
            ["violin_all_clusters"],
        ),
        (
            "Module Scores - Horizontal",
            "violin-horizontal",
            ["violin_horizontal"],
        ),
        ("Other Plots", "other", []),  # Catch-all for remaining figures
    ]

    subsections = []
    used_figures = set()

    for display_name, section_id, patterns in hotspot_categories:
        category_figures = []

        if patterns:
            for pattern in patterns:
                figures_for_pattern = _find_hotspot_figures_by_category(
                    output_dir, pattern
                )
                category_figures.extend(figures_for_pattern)
        else:
            # "Other" category - get all hotspot figures not yet used
            all_hotspot_figures = _find_figures(
                output_dir, ["figures/hotspot", "hotspot"], max_figures=50
            )
            category_figures = [f for f in all_hotspot_figures if f not in used_figures]

        # Remove duplicates and already used figures
        category_figures = list(dict.fromkeys(category_figures))
        category_figures = [f for f in category_figures if f not in used_figures]

        if category_figures:
            used_figures.update(category_figures)

            # Custom descriptions for violin plot categories
            if "violin" in section_id:
                if "per-cluster" in section_id:
                    desc = "<p>Violin plots showing module score distributions for each cluster separately with shared axes for comparison.</p>"
                elif "all" in section_id:
                    desc = "<p>Combined violin plot showing all clusters side-by-side for each module.</p>"
                elif "horizontal" in section_id:
                    desc = "<p>Horizontal violin plot with modules on y-axis for better label readability.</p>"
                else:
                    desc = f"<p>{len(category_figures)} {display_name.lower()} visualizations.</p>"
            else:
                desc = f"<p>{len(category_figures)} {display_name.lower()} visualizations.</p>"

            subsections.append(
                ReportSection(
                    title=display_name,
                    section_id=f"hotspot-{section_id}",
                    content=desc,
                    figures=category_figures[:20],
                )
            )

    # If no subsections were created, fall back to single gallery
    if not subsections:
        figures = _find_figures(
            output_dir, ["figures/hotspot", "hotspot"], max_figures=20
        )
        return ReportSection(
            title="Hotspot Gene Modules",
            section_id="hotspot",
            content=content,
            figures=figures,
            metrics=metrics,
        )

    return ReportSection(
        title="Hotspot Gene Modules",
        section_id="hotspot",
        content=content,
        metrics=metrics,
        subsections=subsections,
    )


def create_operations_log_section(
    log_file: str,
) -> ReportSection:
    """
    Create Operations Log section from pipeline log file.
    """
    log_entries = []

    try:
        with open(log_file, "r") as f:
            for line in f:
                # Parse log entries
                match = re.match(
                    r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}),\d+ - (\w+) - (.+)",
                    line.strip(),
                )
                if match:
                    timestamp, level, message = match.groups()
                    log_entries.append(
                        {
                            "timestamp": timestamp,
                            "level": level,
                            "message": message,
                        }
                    )
    except Exception as e:
        log_entries = [
            {
                "timestamp": "",
                "level": "ERROR",
                "message": f"Could not read log file: {e}",
            }
        ]

    # Build log HTML
    log_html = ['<div class="operations-log">']

    for entry in log_entries[-100:]:  # Last 100 entries
        level = entry.get("level", "INFO")
        if level == "ERROR":
            status_class = "status-error"
        elif "COMPLETED" in entry.get("message", ""):
            status_class = "status-completed"
        elif "SKIPPED" in entry.get("message", ""):
            status_class = "status-skipped"
        else:
            status_class = ""

        log_html.append(
            f"""
        <div class="log-entry">
            <span class="timestamp">{entry.get('timestamp', '')}</span>
            <span class="step">{entry.get('message', '')}</span>
            <span class="status-badge {status_class}">{level}</span>
        </div>
        """
        )

    log_html.append("</div>")

    content = f"""
    <p>Complete log of all operations performed during the pipeline execution.</p>
    <p><strong>Log file:</strong> <code>{log_file}</code></p>
    <p><strong>Showing:</strong> Last {min(100, len(log_entries))} entries</p>

    {"".join(log_html)}
    """

    return ReportSection(
        title="Operations Log",
        section_id="operations-log",
        content=content,
    )


def create_plot_gallery_section(
    output_dir: str,
) -> ReportSection:
    """
    Create a summary section with links to plot locations.

    Note: Individual plots are shown in their respective analysis sections
    (QC, Clustering, GRN, Hotspot) to avoid duplication.
    """
    # Count figures by category without including them
    qc_count = len(_find_figures(output_dir, ["figures/qc", "qc"], max_figures=100))
    grn_count = len(
        _find_figures(
            output_dir,
            ["figures/grn", "figures/grn/grn_deep_analysis", "grn_deep_analysis"],
            max_figures=100,
        )
    )
    hotspot_count = len(
        _find_figures(output_dir, ["figures/hotspot", "hotspot"], max_figures=100)
    )
    embedding_count = len(_find_embedding_figures(output_dir, max_figures=100))

    total_figures = qc_count + grn_count + hotspot_count + embedding_count

    content = f"""
    <p>All plots are shown in their respective analysis sections above.</p>
    
    <h3>Plot Summary</h3>
    <ul>
        <li><strong>Quality Control:</strong> {qc_count} plots 
            (see <a href="#qc">QC Section</a>)</li>
        <li><strong>Clustering & Embedding:</strong> {embedding_count} plots 
            (see <a href="#clustering">Clustering Section</a>)</li>
        <li><strong>GRN Analysis:</strong> {grn_count} plots 
            (see <a href="#celloracle">CellOracle Section</a>)</li>
        <li><strong>Hotspot Modules:</strong> {hotspot_count} plots 
            (see <a href="#hotspot">Hotspot Section</a>)</li>
    </ul>
    
    <h3>Figure Locations</h3>
    <pre><code>{output_dir}/figures/
├── qc/           # Quality control plots
├── grn/          # GRN and network analysis plots
└── hotspot/      # Gene module analysis plots</code></pre>
    
    <p><strong>Total plots generated:</strong> {total_figures}</p>
    """

    return ReportSection(
        title="Plot Summary",
        section_id="plot-gallery",
        content=content,
        metrics={"Total Plots": total_figures},
    )


def create_grn_deep_analysis_section(
    output_dir: str,
    merged_scores=None,
) -> ReportSection:
    """
    Create GRN Deep Analysis section with merged score heatmaps.

    Displays cross-stratification heatmaps generated by plot_heatmap_scores(),
    organized by score type as subsections (tabs).

    Parameters
    ----------
    output_dir : str
        Output directory containing figures/grn/ heatmaps.
    merged_scores : pd.DataFrame, optional
        Merged scores DataFrame from merge_scores().
    """
    import pandas as pd

    metrics = {}

    # Load merged scores CSV if DataFrame not passed
    merged_csv = os.path.join(output_dir, "celloracle", "total_merged_scores.csv")
    if merged_scores is None and os.path.exists(merged_csv):
        try:
            merged_scores = pd.read_csv(merged_csv)
        except Exception:
            merged_scores = None

    if merged_scores is not None and not merged_scores.empty:
        if "stratification" in merged_scores.columns:
            metrics["Stratifications"] = len(merged_scores["stratification"].unique())
        if "cluster" in merged_scores.columns:
            metrics["Clusters"] = len(merged_scores["cluster"].unique())
        if "gene" in merged_scores.columns:
            metrics["Transcription Factors"] = len(merged_scores["gene"].unique())

    # Find all heatmap figures
    grn_dir = os.path.join(output_dir, "figures", "grn")
    heatmap_figures = []
    if os.path.exists(grn_dir):
        for pattern in ["*.png", "*.jpg", "*.svg"]:
            for f in glob.glob(os.path.join(grn_dir, "**", pattern), recursive=True):
                if "grn_heatmap_" in os.path.basename(f).lower():
                    heatmap_figures.append(f)
    heatmap_figures = sorted(set(heatmap_figures))

    content = """
    <p>Cross-stratification GRN analysis comparing transcription factor centrality scores
    across different stratifications. Heatmaps show the difference in centrality scores
    between pairs of clusters, highlighting TFs with divergent regulatory roles.</p>

    <h3>Score Types</h3>
    <ul>
        <li><strong>Degree centrality:</strong> Number of regulatory connections (in/out/all)</li>
        <li><strong>Betweenness centrality:</strong> How often a TF lies on shortest paths between other genes</li>
        <li><strong>Eigenvector centrality:</strong> Influence based on connections to other influential nodes</li>
    </ul>

    <h3>Interpretation</h3>
    <p>Positive values (red) indicate higher centrality in the first cluster;
    negative values (blue) indicate higher centrality in the second cluster.
    Genes are ranked by absolute difference across stratifications.</p>
    """

    if os.path.exists(merged_csv):
        content += f"""
    <h3>Output Files</h3>
    <ul>
        <li><code>celloracle/total_merged_scores.csv</code> - Merged centrality scores across all stratifications</li>
    </ul>
    """

    # Group heatmaps by score type
    score_categories = [
        ("Eigenvector Centrality", "eigenvector", ["eigenvector_centrality"]),
        ("Betweenness Centrality", "betweenness", ["betweenness_centrality"]),
        ("Degree (All)", "degree-all", ["degree_all", "degree_centrality_all"]),
        ("Degree (In)", "degree-in", ["degree_in", "degree_centrality_in"]),
        ("Degree (Out)", "degree-out", ["degree_out", "degree_centrality_out"]),
    ]

    subsections = []
    used_figures = set()

    for display_name, section_id, patterns in score_categories:
        category_figures = []
        for fig in heatmap_figures:
            fname = os.path.basename(fig).lower()
            # Match heatmap files for this score type
            for pat in patterns:
                if pat in fname and fig not in used_figures:
                    category_figures.append(fig)
                    break

        if category_figures:
            used_figures.update(category_figures)
            desc = (
                f"<p>{len(category_figures)} heatmap(s) showing {display_name.lower()} "
                f"differences between cluster pairs across stratifications.</p>"
            )
            subsections.append(
                ReportSection(
                    title=display_name,
                    section_id=f"grn-deep-{section_id}",
                    content=desc,
                    figures=sorted(category_figures)[:30],
                )
            )

    # Catch any remaining heatmaps not matched above
    remaining = [f for f in heatmap_figures if f not in used_figures]
    if remaining:
        subsections.append(
            ReportSection(
                title="Other Heatmaps",
                section_id="grn-deep-other",
                content=f"<p>{len(remaining)} additional heatmap(s).</p>",
                figures=sorted(remaining)[:30],
            )
        )

    metrics["Heatmap Plots"] = len(heatmap_figures)

    return ReportSection(
        title="GRN Deep Analysis (Merged Scores)",
        section_id="grn-deep-analysis",
        content=content,
        metrics=metrics,
        subsections=subsections,
    )


def create_stratification_summary_section(
    stratification_names: List[str],
    output_dir: str,
) -> ReportSection:
    """
    Create Stratified Analysis Summary section.
    """
    n_strats = len(stratification_names)

    # Build stratification list with links
    strat_items = []
    for name in stratification_names:
        strat_dir = os.path.join(output_dir, "stratified_analysis", name)
        if os.path.exists(strat_dir):
            strat_items.append(f"<li><strong>{name}</strong> - Analysis completed</li>")
        else:
            strat_items.append(f"<li><strong>{name}</strong> - Pending</li>")

    content = f"""
    <p>The data was stratified into {n_strats} subsets for independent analysis.
    Each stratification was processed through the full analysis pipeline.</p>

    <h3>Stratifications Analyzed</h3>
    <ul>{"".join(strat_items)}</ul>

    <h3>Output Structure</h3>
    <p>Results for each stratification are saved in:</p>
    <pre><code>{output_dir}/stratified_analysis/&lt;stratification_name&gt;/
├── clustered_adata.h5ad
├── celloracle/
├── hotspot/
├── figures/
└── logs/</code></pre>
    """

    return ReportSection(
        title="Stratified Analysis Summary",
        section_id="stratification",
        content=content,
        metrics={"Stratifications": n_strats},
    )
