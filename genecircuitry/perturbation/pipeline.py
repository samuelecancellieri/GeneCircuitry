"""
Perturbation Pipeline Orchestration Module
==========================================

Orchestrates the complete in-silico perturbation analysis workflow:
- TF selection with differential centrality and Hotspot module enrichment
- In-silico signal propagation simulation (KO / OE)
- Downstream target shift analysis
- Cluster transition quantification and developmental inner product scoring
- Saving data tables and generating publication-ready comparative figures
"""

import os
import pickle
from typing import Optional, Dict, Any, List, Tuple, Union
import numpy as np
import pandas as pd
from anndata import AnnData

from .. import config
from ..logging_utils import log_info, log_warning, log_error
from .tf_selection import select_top_perturbed_tfs
from .simulation import (
    setup_oracle_for_simulation,
    simulate_single_tf_perturbation,
    extract_target_gene_shifts,
)
from .development_analysis import (
    compute_cluster_displacement_score,
    compute_perturbation_cluster_transitions,
    compute_developmental_flow_inner_product,
)
from ..plotting.perturbation_plots import generate_all_perturbation_plots


def _save_markdown_selection_report(selection_df: pd.DataFrame, filepath: str) -> None:
    """Save an explanatory markdown report for TF selection."""
    lines = [
        "# Transcription Factor Selection Rationale for In-Silico Perturbation",
        "",
        "This report details the quantitative criteria and biological rationale used to identify",
        "top regulatory drivers for in-silico perturbation between compared clusters.",
        "",
        "| TF | Comparison | Direction | Delta Degree | Delta Eigenvector | Target Module | Module Overlap | Log2FC | Selection Reason |",
        "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    for _, r in selection_df.iterrows():
        tf = r.get("tf", "")
        comp = r.get("comparison", "")
        ptype = r.get("perturbation_type", "KO")
        d_deg = f"{r.get('delta_degree', 0.0):+.2f}"
        d_eig = f"{r.get('delta_eigenvector', 0.0):+.2f}"
        mod = r.get("top_target_module", "None")
        ratio = f"{r.get('target_module_overlap_ratio', 0.0)*100:.1f}%"
        lfc_val = r.get("log2fc")
        lfc_str = f"{lfc_val:+.2f}" if lfc_val is not None and not pd.isna(lfc_val) else "N/A"
        reason = r.get("selection_reason", "")

        lines.append(f"| **{tf}** | {comp} | {ptype} | {d_deg} | {d_eig} | {mod} | {ratio} | {lfc_str} | {reason} |")

    lines.append("")
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def run_perturbation_analysis(
    oracle: Optional[Any] = None,
    links: Optional[Any] = None,
    score_df: Optional[Union[pd.DataFrame, str]] = None,
    links_df: Optional[Union[pd.DataFrame, str]] = None,
    module_genes: Optional[Dict[str, Any]] = None,
    adata: Optional[AnnData] = None,
    compare_clusters: Optional[Union[List[Tuple[str, str]], str]] = None,
    cluster_key: str = "cluster",
    top_n_tfs: int = 3,
    perturbation_mode: str = "directional",
    n_propagation: int = 3,
    pseudotime_key: str = "pseudotime",
    output_dir: Optional[str] = None,
    save_figures: bool = True,
    **kwargs,
) -> Dict[str, Any]:
    """
    Run end-to-end in-silico perturbation analysis.

    Parameters
    ----------
    oracle : co.Oracle or str, optional
        Oracle object or path to saved oracle file.
    links : co.Links or str, optional
        Links object or path to saved links file.
    score_df : pd.DataFrame or str, optional
        Centrality scores DataFrame or path to grn_merged_scores.csv.
    links_df : pd.DataFrame or str, optional
        GRN filtered links DataFrame or path to grn_filtered_links.pkl.
    module_genes : dict, optional
        Hotspot module gene sets.
    adata : AnnData, optional
        Expression AnnData object.
    compare_clusters : list of tuple or str, optional
        Cluster pairs to compare (e.g. [('ClusterA', 'ClusterB')] or 'ClusterA,ClusterB').
    cluster_key : str, default="cluster"
        Cluster column name in adata.obs and score_df.
    top_n_tfs : int, default=3
        Number of top TFs to select per comparison.
    perturbation_mode : str, default="directional"
        Mode: 'directional' (KO for high, OE for low), 'ko', 'oe', or 'both'.
    n_propagation : int, default=3
        Signal propagation steps.
    pseudotime_key : str, default="pseudotime"
        Pseudotime key in adata.obs.
    output_dir : str, optional
        Base output directory (defaults to config.OUTPUT_DIR).
    save_figures : bool, default=True
        Whether to generate and save comparative figures.

    Returns
    -------
    dict
        Dictionary containing selection_df, simulation_summary_df, results_dict,
        and list of generated figure paths.
    """
    base_out = output_dir or config.OUTPUT_DIR
    pert_dir = os.path.join(base_out, "perturbation")
    os.makedirs(pert_dir, exist_ok=True)

    log_info("Perturbation.Pipeline", f"STARTED - output_dir: {pert_dir}")

    # Load score_df if string path
    if isinstance(score_df, str) and os.path.exists(score_df):
        score_df = pd.read_csv(score_df, index_col=0)
        if "gene" not in score_df.columns:
            score_df.index.name = "gene"
            score_df = score_df.reset_index()

    # Load links_df if string path
    if isinstance(links_df, str) and os.path.exists(links_df):
        try:
            with open(links_df, "rb") as f:
                loaded_links = pickle.load(f)
            if isinstance(loaded_links, dict):
                dfs = []
                for c_name, c_df in loaded_links.items():
                    if isinstance(c_df, pd.DataFrame):
                        df_c = c_df.copy()
                        df_c["cluster"] = str(c_name)
                        dfs.append(df_c)
                links_df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
            elif isinstance(loaded_links, pd.DataFrame):
                links_df = loaded_links
        except Exception as e:
            log_warning("Perturbation.Pipeline", f"Could not load links_df: {e}")
            links_df = None

    # Load Oracle / Links if string paths
    if isinstance(oracle, str) and os.path.exists(oracle):
        try:
            import celloracle as co
            oracle = co.load_hdf5(oracle)
        except Exception as e:
            log_error("Perturbation.LoadOracle", e)

    if isinstance(links, str) and os.path.exists(links):
        try:
            import celloracle as co
            links = co.load_hdf5(links)
        except Exception as e:
            log_error("Perturbation.LoadLinks", e)

    if oracle is None:
        log_error("Perturbation.Pipeline", ValueError("Oracle object must be provided for perturbation."))
        return {}

    # Setup Oracle GRNs for simulation
    oracle = setup_oracle_for_simulation(oracle, links)

    # Resolve AnnData
    if adata is None and hasattr(oracle, "adata"):
        adata = oracle.adata

    # Step 1: Select top TFs
    log_info("Perturbation.Pipeline", "Selecting top candidate TFs for perturbation...")
    selection_df = select_top_perturbed_tfs(
        score_df=score_df,
        links_df=links_df,
        module_genes=module_genes,
        adata=adata,
        compare_clusters=compare_clusters,
        cluster_key=cluster_key,
        top_n_tfs=top_n_tfs,
        mode=perturbation_mode,
    )

    if selection_df.empty:
        log_warning("Perturbation.Pipeline", "No candidate TFs selected; terminating perturbation pipeline.")
        return {"selection_df": pd.DataFrame(), "simulation_summary_df": pd.DataFrame(), "generated_figures": []}

    # Save selection tables
    sel_csv = os.path.join(pert_dir, "tf_selection_rationale.csv")
    sel_md = os.path.join(pert_dir, "tf_selection_rationale.md")
    selection_df.to_csv(sel_csv, index=False)
    _save_markdown_selection_report(selection_df, sel_md)
    log_info("Perturbation.Pipeline", f"Saved selection rationale to {sel_csv} and {sel_md}")

    # Step 2: Run in-silico simulations
    results_dict: Dict[str, Any] = {}
    summary_rows: List[Dict[str, Any]] = []
    all_shifts_list: List[pd.DataFrame] = []
    all_displacements_list: List[pd.DataFrame] = []

    for idx, (_, row) in enumerate(selection_df.iterrows()):
        tf = row["tf"]
        comp = row["comparison"]
        c_high = row["cluster_high"]
        c_low = row["cluster_low"]
        ptype = row["perturbation_type"]

        # Determine condition value: 0.0 for KO, 95th percentile for OE
        if ptype == "OE":
            cond_val = 5.0
            if adata is not None and tf in adata.var_names:
                vals = adata[:, tf].X.toarray().flatten() if hasattr(adata[:, tf].X, "toarray") else adata[:, tf].X.flatten()
                cond_val = float(np.percentile(vals, 95))
        else:
            cond_val = 0.0

        log_info("Perturbation.Pipeline", f"[{idx+1}/{len(selection_df)}] Simulating {tf} ({ptype}={cond_val}) for {comp}")

        # Simulate shift
        oracle_sim = simulate_single_tf_perturbation(
            oracle=oracle,
            tf=tf,
            condition_val=cond_val,
            n_propagation=n_propagation,
            cluster_key=cluster_key,
        )

        # Target gene shifts
        delta_x_df = extract_target_gene_shifts(
            oracle=oracle_sim,
            tf=tf,
            links_df=links_df,
            cluster=c_high,
        )
        if not delta_x_df.empty:
            delta_x_df["comparison"] = comp
            all_shifts_list.append(delta_x_df)

        # Cluster displacement score towards target centroid
        disp_df = compute_cluster_displacement_score(
            oracle=oracle_sim,
            cluster_source=c_high,
            cluster_target=c_low,
            cluster_key=cluster_key,
        )
        mean_disp = float(disp_df["displacement_score"].mean()) if not disp_df.empty else 0.0
        if not disp_df.empty:
            disp_df["tf"] = tf
            disp_df["comparison"] = comp
            all_displacements_list.append(disp_df)

        # Developmental flow & inner product
        dev_module, ip_df = compute_developmental_flow_inner_product(
            oracle=oracle_sim,
            tf=tf,
            mode=ptype,
            pseudotime_key=pseudotime_key,
            cluster_key=cluster_key,
        )

        # Markov transition table
        trans_df = compute_perturbation_cluster_transitions(
            oracle=oracle_sim,
            cluster_key=cluster_key,
        )

        results_dict[tf] = {
            "tf": tf,
            "comparison": comp,
            "mode": ptype,
            "condition_val": cond_val,
            "oracle_sim": oracle_sim,
            "oracle": oracle_sim,
            "delta_x_df": delta_x_df,
            "displacement": mean_disp,
            "shift_df": disp_df,
            "dev_module": dev_module,
            "ip_df": ip_df,
            "trans_df": trans_df,
        }

        summary_rows.append({
            "tf": tf,
            "comparison": comp,
            "perturbation_type": ptype,
            "condition_val": cond_val,
            "mean_displacement_score": mean_disp,
            "cluster_high": c_high,
            "cluster_low": c_low,
            "selection_reason": row.get("selection_reason", ""),
        })

    # Save summary tables
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(pert_dir, "simulation_summary.csv")
    summary_df.to_csv(summary_csv, index=False)

    if all_shifts_list:
        all_shifts_df = pd.concat(all_shifts_list, ignore_index=True)
        all_shifts_df.to_csv(os.path.join(pert_dir, "target_gene_shifts.csv"), index=False)

    if all_displacements_list:
        all_disp_df = pd.concat(all_displacements_list, ignore_index=True)
        all_disp_df.to_csv(os.path.join(pert_dir, "cluster_displacement_scores.csv"), index=False)

    # Step 3: Generate comparative figures
    generated_figures = []
    if save_figures:
        log_info("Perturbation.Pipeline", "Generating comparative perturbation figures...")
        generated_figures = generate_all_perturbation_plots(
            selection_df=selection_df,
            perturbation_results=results_dict,
            output_dir=base_out,
        )

    log_info(
        "Perturbation.Pipeline",
        f"COMPLETED - n_simulated_tfs: {len(results_dict)}, n_figures: {len(generated_figures)}",
    )

    return {
        "selection_df": selection_df,
        "simulation_summary_df": summary_df,
        "results_dict": results_dict,
        "generated_figures": generated_figures,
    }
