"""
Perturbation Plotting Module for GeneCircuitry
==============================================

Provides visualization routines for in-silico TF perturbation analysis:
- TF selection rationale dashboard (differential centralities, Hotspot module overlap, expression)
- 2D embedding vector field plots with quiver / grid arrows
- Developmental gradient flow and inner product (PS) layouts
- Cluster shift magnitude and displacement projection scores
- Downstream target gene expression shifts (Delta X waterfall/heatmap)
- Multi-TF comparative summary dashboard
"""

import os
from typing import Optional, Dict, Any, List, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

from .. import config
from .utils import save_plot


def _get_save_path(filename: str, save_path: Optional[str] = None) -> str:
    """Resolve file path for perturbation figures."""
    if save_path:
        return save_path
    fig_dir = os.path.join(config.OUTPUT_DIR, "figures", "perturbation")
    os.makedirs(fig_dir, exist_ok=True)
    return os.path.join(fig_dir, filename)


def plot_tf_selection_rationale(
    selection_df: pd.DataFrame,
    comparison: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Plot a multi-panel dashboard illustrating why and how candidate TFs were selected,
    showing differential network centrality metrics, Hotspot module overlap, and expression.

    Parameters
    ----------
    selection_df : pd.DataFrame
        Table from select_top_perturbed_tfs.
    comparison : str, optional
        Comparison label (e.g. "ClusterA vs ClusterB").
    save_path : str, optional
        Output file path.

    Returns
    -------
    plt.Figure or None
    """
    if selection_df is None or selection_df.empty:
        return None

    df = selection_df.copy()
    if comparison is not None and "comparison" in df.columns:
        df_comp = df[df["comparison"] == comparison]
        if not df_comp.empty:
            df = df_comp

    tfs = df["tf"].unique().tolist()
    if not tfs:
        return None

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), constrained_layout=True)
    comp_title = comparison or df["comparison"].iloc[0] if "comparison" in df.columns else ""
    fig.suptitle(f"TF Selection Rationale: {comp_title}", fontsize=14, fontweight="bold")

    y_pos = np.arange(len(df))
    labels = [f"{r['tf']} ({r.get('cluster_high', '')})" for _, r in df.iterrows()]

    # Panel 1: Differential Centralities (Degree, Betweenness, Eigenvector)
    ax1 = axes[0]
    bar_width = 0.25
    if "delta_degree" in df.columns:
        ax1.barh(y_pos - bar_width, df["delta_degree"], height=bar_width, label=r"$\Delta$ Degree", color="#1f77b4")
    if "delta_betweenness" in df.columns:
        ax1.barh(y_pos, df["delta_betweenness"], height=bar_width, label=r"$\Delta$ Betweenness", color="#ff7f0e")
    if "delta_eigenvector" in df.columns:
        ax1.barh(y_pos + bar_width, df["delta_eigenvector"], height=bar_width, label=r"$\Delta$ Eigenvector", color="#2ca02c")

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels, fontsize=10, fontweight="medium")
    ax1.set_xlabel("Centrality Delta", fontsize=11)
    ax1.set_title("Network Centrality Differences", fontsize=12)
    ax1.axvline(0, color="gray", linestyle="--", linewidth=0.8)
    ax1.legend(loc="lower right", frameon=True, fontsize=9)
    ax1.grid(axis="x", linestyle=":", alpha=0.6)

    # Panel 2: Module Overlap & Expression Log2FC
    ax2 = axes[1]
    if "target_module_overlap_ratio" in df.columns:
        ax2.barh(y_pos - 0.15, df["target_module_overlap_ratio"] * 100, height=0.3, label="Target Module Overlap (%)", color="#9467bd")
    if "composite_score" in df.columns:
        ax2.barh(y_pos + 0.15, df["composite_score"] * 100, height=0.3, label="Composite Driver Score", color="#8c564b")

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([])
    ax2.set_xlabel("Percentage / Score (%)", fontsize=11)
    ax2.set_title("Module Overlap & Driver Score", fontsize=12)
    ax2.legend(loc="lower right", frameon=True, fontsize=9)
    ax2.grid(axis="x", linestyle=":", alpha=0.6)

    # Panel 3: Explanatory Selection Reasons
    ax3 = axes[2]
    ax3.axis("off")
    ax3.set_title("Selection Rationale & Evidence", fontsize=12, fontweight="medium")

    text_lines = []
    for idx, (_, r) in enumerate(df.iterrows()):
        tf_name = r["tf"]
        reason = r.get("selection_reason", f"Selected as top regulator for {r.get('cluster_high', 'group')}")
        text_lines.append(f"• {tf_name}:\n  {reason}\n")

    summary_text = "\n".join(text_lines)
    ax3.text(
        0.05,
        0.95,
        summary_text,
        transform=ax3.transAxes,
        fontsize=9,
        verticalalignment="top",
        wrap=True,
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#f8f9fa", edgecolor="#ced4da"),
    )

    out_path = _get_save_path("tf_selection_rationale.png", save_path)
    save_plot(fig, out_path, plot_type="perturbation", close_fig=False, skip_existing=False)
    return fig


def plot_perturbation_vector_field(
    oracle: Any,
    tf: str,
    mode: str = "KO",
    cluster_key: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Plot 2D embedding vector field (quiver / grid arrows) representing simulated
    cell identity transitions following TF perturbation.

    Parameters
    ----------
    oracle : co.Oracle
        Oracle object after simulation.
    tf : str
        Perturbed gene name.
    mode : str, default="KO"
        Perturbation condition.
    cluster_key : str, optional
        Cluster column in adata.obs.
    save_path : str, optional
        Output file path.

    Returns
    -------
    plt.Figure or None
    """
    if not hasattr(oracle, "adata"):
        return None

    adata = oracle.adata
    emb_key = getattr(oracle, "embedding_name", "X_draw_graph_fa")
    c_col = cluster_key or getattr(oracle, "cluster_column_name", "cluster")

    if emb_key not in adata.obsm:
        return None

    coords = adata.obsm[emb_key]
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)

    # Plot cell background scatter
    if c_col in adata.obs.columns:
        clusters = adata.obs[c_col].astype(str)
        unique_c = sorted(list(clusters.unique()))
        cmap = plt.cm.get_cmap("tab20", len(unique_c))
        for i, c in enumerate(unique_c):
            c_mask = clusters == c
            ax.scatter(
                coords[c_mask, 0],
                coords[c_mask, 1],
                s=18,
                color=cmap(i),
                label=c,
                alpha=0.65,
                edgecolors="none",
            )
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, markerscale=2)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=15, color="#ced4da", alpha=0.5)

    # Overlay flow arrows if available
    flow = getattr(oracle, "flow", None)
    if flow is not None and "grid_x" in flow and "flow_x" in flow:
        gx = flow["grid_x"]
        gy = flow["grid_y"]
        fx = flow["flow_x"]
        fy = flow["flow_y"]
        # Meshgrid if 1D coords
        if gx.ndim == 1 and gy.ndim == 1:
            X_grid, Y_grid = np.meshgrid(gx, gy)
        else:
            X_grid, Y_grid = gx, gy

        speed = np.sqrt(fx**2 + fy**2)
        ax.quiver(
            X_grid,
            Y_grid,
            fx,
            fy,
            color="#212529",
            alpha=0.85,
            scale=None,
            pivot="mid",
        )
    elif getattr(oracle, "delta_embedding", None) is not None:
        delta = oracle.delta_embedding
        # Subsample for clear visualization
        sub_idx = np.random.choice(coords.shape[0], size=min(300, coords.shape[0]), replace=False)
        ax.quiver(
            coords[sub_idx, 0],
            coords[sub_idx, 1],
            delta[sub_idx, 0],
            delta[sub_idx, 1],
            color="#212529",
            alpha=0.7,
            pivot="mid",
        )

    ax.set_title(f"In-Silico Perturbation: {tf} {mode} Flow", fontsize=13, fontweight="bold")
    ax.set_xlabel(f"{emb_key} 1", fontsize=10)
    ax.set_ylabel(f"{emb_key} 2", fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.4)

    out_path = _get_save_path(f"{tf}_{mode}_vector_field.png", save_path)
    save_plot(fig, out_path, plot_type="perturbation", close_fig=False, skip_existing=False)
    return fig


def plot_cluster_shift_magnitude(
    shift_df: pd.DataFrame,
    comparison: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Plot displacement projection scores and shift magnitudes towards the target cluster.

    Parameters
    ----------
    shift_df : pd.DataFrame
        Table from compute_cluster_displacement_score.
    comparison : str, optional
        Comparison title.
    save_path : str, optional
        Output file path.

    Returns
    -------
    plt.Figure or None
    """
    if shift_df is None or shift_df.empty:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    title = f"Shift Towards Target Cluster: {comparison}" if comparison else "Cluster Shift Magnitude"
    fig.suptitle(title, fontsize=13, fontweight="bold")

    # Panel 1: Displacement projection score
    ax1 = axes[0]
    sns.violinplot(
        data=shift_df,
        y="displacement_score",
        ax=ax1,
        color="#4a90e2",
        inner="quartile",
    )
    ax1.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax1.set_ylabel("Displacement Projection Score", fontsize=11)
    ax1.set_title("Directional Shift Toward Target Centroid", fontsize=11)
    ax1.grid(axis="y", linestyle=":", alpha=0.6)

    # Panel 2: Total velocity magnitude
    ax2 = axes[1]
    sns.boxplot(
        data=shift_df,
        y="shift_magnitude",
        ax=ax2,
        color="#50e3c2",
        width=0.4,
    )
    ax2.set_ylabel("Total Shift Magnitude", fontsize=11)
    ax2.set_title("Overall Transition Magnitude", fontsize=11)
    ax2.grid(axis="y", linestyle=":", alpha=0.6)

    out_path = _get_save_path("cluster_shift_magnitude.png", save_path)
    save_plot(fig, out_path, plot_type="perturbation", close_fig=False, skip_existing=False)
    return fig


def plot_target_gene_shifts(
    delta_x_df: pd.DataFrame,
    tf: Optional[str] = None,
    top_n: int = 15,
    save_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Plot waterfall horizontal barplot of downstream target expression shifts (Delta X).

    Parameters
    ----------
    delta_x_df : pd.DataFrame
        Table from extract_target_gene_shifts.
    tf : str, optional
        Target TF name.
    top_n : int, default=15
        Number of top target genes to display.
    save_path : str, optional
        Output file path.

    Returns
    -------
    plt.Figure or None
    """
    if delta_x_df is None or delta_x_df.empty:
        return None

    df = delta_x_df.sort_values("abs_delta_x", ascending=False).head(top_n).copy()
    df = df.sort_values("mean_delta_x", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(4.5, len(df) * 0.3)), constrained_layout=True)
    tf_name = tf or df["tf"].iloc[0] if "tf" in df.columns else "TF"

    colors = ["#e74c3c" if x > 0 else "#3498db" for x in df["mean_delta_x"]]
    y_pos = np.arange(len(df))

    ax.barh(y_pos, df["mean_delta_x"], color=colors, height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["target"], fontsize=10)
    ax.axvline(0, color="black", linestyle="-", linewidth=0.8)
    ax.set_xlabel(r"Simulated Expression Shift ($\Delta X$)", fontsize=11)
    ax.set_title(f"Downstream Target Shifts upon {tf_name} Perturbation", fontsize=12, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.5)

    out_path = _get_save_path(f"{tf_name}_target_shifts.png", save_path)
    save_plot(fig, out_path, plot_type="perturbation", close_fig=False, skip_existing=False)
    return fig


def plot_perturbation_developmental_flow(
    dev_module: Any,
    oracle: Any,
    tf: str,
    mode: str = "KO",
    save_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Multi-panel layout comparing perturbation simulation against developmental gradients
    following the CellOracle tutorial notebook layout.

    Parameters
    ----------
    dev_module : co.Oracle_development_module
        Development module object after inner product calculation.
    oracle : co.Oracle
        Oracle object.
    tf : str
        Perturbed gene symbol.
    mode : str, default="KO"
        Perturbation condition.
    save_path : str, optional
        Output file path.

    Returns
    -------
    plt.Figure or None
    """
    if not hasattr(oracle, "adata"):
        return None

    adata = oracle.adata
    emb_key = getattr(oracle, "embedding_name", "X_draw_graph_fa")
    c_col = getattr(oracle, "cluster_column_name", "cell_type")

    if emb_key not in adata.obsm:
        return None

    coords = adata.obsm[emb_key]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    fig.suptitle(f"Developmental Flow & Perturbation Score (PS): {tf} {mode}", fontsize=14, fontweight="bold")

    # Panel 1: Cluster Map
    ax1 = axes[0, 0]
    if c_col in adata.obs.columns:
        clusters = adata.obs[c_col].astype(str)
        unique_c = sorted(list(clusters.unique()))
        cmap = plt.cm.get_cmap("tab20", len(unique_c))
        for i, c in enumerate(unique_c):
            mask = clusters == c
            ax1.scatter(coords[mask, 0], coords[mask, 1], s=14, color=cmap(i), label=c, alpha=0.6)
        ax1.legend(loc="upper right", frameon=False, markerscale=1.5, fontsize=8)
    else:
        ax1.scatter(coords[:, 0], coords[:, 1], s=14, color="#ced4da", alpha=0.6)
    ax1.set_title("Cell Clusters", fontsize=11, fontweight="medium")

    # Panel 2: Pseudotime
    ax2 = axes[0, 1]
    pt_col = next((c for c in ["pseudotime", "dpt_pseudotime", "latent_time", "_inferred_pseudotime"] if c in adata.obs.columns), None)
    if pt_col:
        pt_vals = adata.obs[pt_col]
        sc = ax2.scatter(coords[:, 0], coords[:, 1], c=pt_vals, cmap="viridis", s=14, alpha=0.8)
        fig.colorbar(sc, ax=ax2, label="Pseudotime", fraction=0.046, pad=0.04)
        ax2.set_title("Differentiation Pseudotime", fontsize=11, fontweight="medium")
    else:
        ax2.text(0.5, 0.5, "No pseudotime available", ha="center", va="center")
        ax2.set_title("Pseudotime", fontsize=11)

    # Panel 3: Perturbation Simulation Flow
    ax3 = axes[1, 0]
    ax3.scatter(coords[:, 0], coords[:, 1], s=10, color="#ced4da", alpha=0.4)
    flow = getattr(oracle, "flow", None)
    if flow is not None and "grid_x" in flow and "flow_x" in flow:
        gx, gy = flow["grid_x"], flow["grid_y"]
        if gx.ndim == 1 and gy.ndim == 1:
            X_grid, Y_grid = np.meshgrid(gx, gy)
        else:
            X_grid, Y_grid = gx, gy
        ax3.quiver(X_grid, Y_grid, flow["flow_x"], flow["flow_y"], color="#212529", alpha=0.8, pivot="mid")
    ax3.set_title(f"Perturbation Flow ({tf} {mode})", fontsize=11, fontweight="medium")

    # Panel 4: Inner Product (PS) Distribution
    ax4 = axes[1, 1]
    ip_vals = None
    if dev_module is not None and hasattr(dev_module, "inner_product"):
        ip_vals = dev_module.inner_product
    elif dev_module is not None and hasattr(dev_module, "inner_product_df"):
        ip_vals = dev_module.inner_product_df["inner_product"]

    if ip_vals is not None and len(ip_vals) == adata.n_obs:
        if c_col in adata.obs.columns:
            plot_df = pd.DataFrame({"inner_product": ip_vals, "cluster": adata.obs[c_col]})
            sns.boxplot(data=plot_df, x="cluster", y="inner_product", ax=ax4, palette="Set2")
            ax4.set_xticklabels(ax4.get_xticklabels(), rotation=45, ha="right")
        else:
            sns.histplot(ip_vals, ax=ax4, kde=True, color="#3498db")
        ax4.axhline(0, color="red", linestyle="--", linewidth=0.8)
        ax4.set_title("Perturbation Score (PS) Distribution", fontsize=11, fontweight="medium")
        ax4.set_ylabel("Inner Product (PS)")
    else:
        ax4.text(0.5, 0.5, "No inner product data", ha="center", va="center")
        ax4.set_title("Inner Product", fontsize=11)

    out_path = _get_save_path(f"{tf}_{mode}_developmental_flow.png", save_path)
    save_plot(fig, out_path, plot_type="perturbation", close_fig=False, skip_existing=False)
    return fig


def plot_comparative_perturbation_dashboard(
    results_dict: Dict[str, Dict[str, Any]],
    comparison: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Side-by-side comparative dashboard comparing vector fields and displacement
    metrics across all perturbed candidate TFs.

    Parameters
    ----------
    results_dict : dict
        Dictionary mapping TF symbol to result dict containing 'oracle', 'displacement', etc.
    comparison : str, optional
        Comparison title.
    save_path : str, optional
        Output file path.

    Returns
    -------
    plt.Figure or None
    """
    if not results_dict:
        return None

    tfs = list(results_dict.keys())
    n_tfs = len(tfs)
    if n_tfs == 0:
        return None

    cols = min(3, n_tfs)
    rows = (n_tfs + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 5 * rows), constrained_layout=True)
    if rows == 1 and cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    comp_title = comparison or "All Comparisons"
    fig.suptitle(f"In-Silico Perturbation Multi-TF Comparison: {comp_title}", fontsize=14, fontweight="bold")

    for i, tf in enumerate(tfs):
        ax = axes[i]
        res = results_dict[tf]
        oracle = res.get("oracle")

        if oracle is not None and hasattr(oracle, "adata"):
            adata = oracle.adata
            emb_key = getattr(oracle, "embedding_name", "X_draw_graph_fa")
            c_col = getattr(oracle, "cluster_column_name", "cell_type")

            if emb_key in adata.obsm:
                coords = adata.obsm[emb_key]
                if c_col in adata.obs.columns:
                    clusters = adata.obs[c_col].astype(str)
                    unique_c = sorted(list(clusters.unique()))
                    cmap = plt.cm.get_cmap("tab20", len(unique_c))
                    for c_idx, c in enumerate(unique_c):
                        mask = clusters == c
                        ax.scatter(coords[mask, 0], coords[mask, 1], s=8, color=cmap(c_idx), alpha=0.4)
                else:
                    ax.scatter(coords[:, 0], coords[:, 1], s=8, color="#ced4da", alpha=0.4)

                flow = getattr(oracle, "flow", None)
                if flow is not None and "grid_x" in flow and "flow_x" in flow:
                    gx, gy = flow["grid_x"], flow["grid_y"]
                    if gx.ndim == 1 and gy.ndim == 1:
                        X_grid, Y_grid = np.meshgrid(gx, gy)
                    else:
                        X_grid, Y_grid = gx, gy
                    ax.quiver(X_grid, Y_grid, flow["flow_x"], flow["flow_y"], color="#212529", alpha=0.8, pivot="mid")

        score = res.get("displacement", "")
        score_str = f" (Shift Score: {score:.2f})" if isinstance(score, (int, float)) else ""
        ax.set_title(f"{tf} KO{score_str}", fontsize=11, fontweight="bold")
        ax.axis("off")

    # Turn off unused subplots
    for j in range(n_tfs, len(axes)):
        axes[j].axis("off")

    out_path = _get_save_path("comparative_perturbation_dashboard.png", save_path)
    save_plot(fig, out_path, plot_type="perturbation", close_fig=False, skip_existing=False)
    return fig


def generate_all_perturbation_plots(
    selection_df: pd.DataFrame,
    perturbation_results: Dict[str, Any],
    output_dir: Optional[str] = None,
) -> List[str]:
    """
    Generate and save all perturbation comparative figures.

    Parameters
    ----------
    selection_df : pd.DataFrame
        Table from select_top_perturbed_tfs.
    perturbation_results : dict
        Dictionary of perturbation simulation results per TF.
    output_dir : str, optional
        Base output directory.

    Returns
    -------
    list of str
        List of generated image file paths.
    """
    generated_files = []
    fig_dir = os.path.join(output_dir or config.OUTPUT_DIR, "figures", "perturbation")
    os.makedirs(fig_dir, exist_ok=True)

    # 1. Selection Rationale Dashboard
    if selection_df is not None and not selection_df.empty:
        p1 = os.path.join(fig_dir, "tf_selection_rationale.png")
        fig1 = plot_tf_selection_rationale(selection_df, save_path=p1)
        if fig1 is not None and os.path.exists(p1):
            generated_files.append(p1)

    # 2. Per-TF Figures (Vector field, target shifts, developmental flow)
    for tf, res in perturbation_results.items():
        oracle_sim = res.get("oracle_sim") or res.get("oracle")
        mode = res.get("mode", "KO")

        if oracle_sim is not None:
            # Vector field
            p_vec = os.path.join(fig_dir, f"{tf}_{mode}_vector_field.png")
            fig_vec = plot_perturbation_vector_field(oracle_sim, tf=tf, mode=mode, save_path=p_vec)
            if fig_vec is not None and os.path.exists(p_vec):
                generated_files.append(p_vec)

            # Developmental flow
            dev = res.get("dev_module")
            if dev is not None:
                p_dev = os.path.join(fig_dir, f"{tf}_{mode}_developmental_flow.png")
                fig_dev = plot_perturbation_developmental_flow(dev, oracle_sim, tf=tf, mode=mode, save_path=p_dev)
                if fig_dev is not None and os.path.exists(p_dev):
                    generated_files.append(p_dev)

        # Target gene shifts
        delta_x_df = res.get("delta_x_df")
        if delta_x_df is not None and not delta_x_df.empty:
            p_tgt = os.path.join(fig_dir, f"{tf}_target_shifts.png")
            fig_tgt = plot_target_gene_shifts(delta_x_df, tf=tf, save_path=p_tgt)
            if fig_tgt is not None and os.path.exists(p_tgt):
                generated_files.append(p_tgt)

        # Cluster shift magnitude
        shift_df = res.get("shift_df")
        if shift_df is not None and not shift_df.empty:
            p_shift = os.path.join(fig_dir, f"{tf}_cluster_shift_magnitude.png")
            fig_shift = plot_cluster_shift_magnitude(shift_df, comparison=res.get("comparison"), save_path=p_shift)
            if fig_shift is not None and os.path.exists(p_shift):
                generated_files.append(p_shift)

    # 3. Comparative Dashboard
    if perturbation_results:
        p_dash = os.path.join(fig_dir, "comparative_perturbation_dashboard.png")
        fig_dash = plot_comparative_perturbation_dashboard(perturbation_results, save_path=p_dash)
        if fig_dash is not None and os.path.exists(p_dash):
            generated_files.append(p_dash)

    return generated_files
