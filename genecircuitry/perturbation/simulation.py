"""
In-Silico Signal Propagation Simulation Engine
===============================================

Wraps CellOracle's perturbation simulation engine to calculate signal propagation,
cell identity transitions, and embedding shifts following transcription factor perturbation.
"""

from typing import Optional, Dict, Any, List, Union
import numpy as np
import pandas as pd
from anndata import AnnData

from ..logging_utils import log_info, log_warning, log_error


def setup_oracle_for_simulation(
    oracle: Any,
    links: Any,
    alpha: float = 10.0,
    use_cluster_specific_TFdict: bool = True,
) -> Any:
    """
    Fit GRN ridge regression models on Oracle object to prepare for in-silico simulation.

    Parameters
    ----------
    oracle : co.Oracle
        Preprocessed CellOracle Oracle object with imputed counts and embedding.
    links : co.Links
        Inferred and filtered CellOracle Links object.
    alpha : float, default=10.0
        Regularization strength for simulation GRN fitting.
    use_cluster_specific_TFdict : bool, default=True
        Whether to use cluster-specific GRN links.

    Returns
    -------
    oracle : co.Oracle
        Fitted Oracle object ready for simulate_shift calls.
    """
    log_info("Perturbation.Setup", "Preparing Oracle object for simulation")

    if hasattr(oracle, "get_cluster_specific_TFdict_from_Links"):
        try:
            oracle.get_cluster_specific_TFdict_from_Links(
                links_object=links, ignore_warning=True
            )
        except Exception as e:
            log_warning("Perturbation.Setup", f"Could not extract cluster TF dict: {e}")

    if hasattr(oracle, "fit_GRN_for_simulation"):
        oracle.fit_GRN_for_simulation(
            GRN_unit="cluster",
            alpha=alpha,
            use_cluster_specific_TFdict=use_cluster_specific_TFdict,
            verbose_level=0,
        )

    log_info("Perturbation.Setup", "Oracle GRN fitting for simulation complete")
    return oracle


def simulate_single_tf_perturbation(
    oracle: Any,
    tf: str,
    condition_val: float = 0.0,
    n_propagation: int = 3,
    n_neighbors: int = 200,
    smooth: float = 0.8,
    n_grid: int = 40,
    cluster_key: Optional[str] = None,
    clip_delta_X: bool = True,
) -> Any:
    """
    Simulate single-gene perturbation and compute cell identity transition vectors.

    Parameters
    ----------
    oracle : co.Oracle
        Fitted Oracle object.
    tf : str
        Target transcription factor to perturb.
    condition_val : float, default=0.0
        Perturbation expression value (e.g. 0.0 for knockout).
    n_propagation : int, default=3
        Number of signal propagation iterations through the GRN.
    n_neighbors : int, default=200
        Number of neighbors for transition probability calculation.
    smooth : float, default=0.8
        Kernel smoothing parameter for grid arrows.
    n_grid : int, default=40
        Grid resolution for vector field arrow interpolation.
    cluster_key : str, optional
        Cluster column in adata.obs.
    clip_delta_X : bool, default=True
        Clip simulated shifts to observed physiological ranges.

    Returns
    -------
    oracle_sim : co.Oracle
        Copied Oracle object containing perturbation results and flow vector fields.
    """
    log_info(
        "Perturbation.Simulate",
        f"Simulating perturbation for {tf} = {condition_val} (n_prop={n_propagation})",
    )

    oracle_sim = oracle.copy()
    if cluster_key is not None and hasattr(oracle_sim, "cluster_column_name"):
        oracle_sim.cluster_column_name = cluster_key

    # Check if TF is in adata
    if tf not in oracle_sim.adata.var_names:
        log_warning(
            "Perturbation.Simulate",
            f"Gene {tf} is not present in adata.var_names; simulation cannot proceed.",
        )
        return oracle_sim

    perturb_cond = {tf: float(condition_val)}

    try:
        oracle_sim.simulate_shift(
            perturb_condition=perturb_cond,
            n_propagation=n_propagation,
            clip_delta_X=clip_delta_X,
            ignore_warning=True,
        )
    except Exception as e:
        log_error("Perturbation.SimulateShift", e)
        return oracle_sim

    # Transition probabilities and 2D embedding shifts
    try:
        if hasattr(oracle_sim, "estimate_transition_prob"):
            oracle_sim.estimate_transition_prob(
                n_neighbors=n_neighbors, knn_random=True, sampled_fraction=1.0
            )

        if hasattr(oracle_sim, "calculate_embedding_shift"):
            oracle_sim.calculate_embedding_shift(sigma_corr=0.05)

        if hasattr(oracle_sim, "calculate_p_mass"):
            oracle_sim.calculate_p_mass(
                smooth=smooth, n_grid=n_grid, n_neighbors=n_neighbors
            )

        if hasattr(oracle_sim, "calculate_grid_arrows"):
            oracle_sim.calculate_grid_arrows(
                smooth=smooth, steps=(n_grid, n_grid), n_neighbors=n_neighbors
            )

        log_info(
            "Perturbation.Simulate",
            f"Successfully computed shift vectors and grid arrows for {tf}",
        )
    except Exception as e:
        log_warning(
            "Perturbation.SimulateDownstream",
            f"Error computing transition probabilities or grid arrows for {tf}: {e}",
        )

    return oracle_sim


def extract_target_gene_shifts(
    oracle: Any,
    tf: str,
    links_df: Optional[pd.DataFrame] = None,
    target_genes: Optional[List[str]] = None,
    cluster: Optional[str] = None,
) -> pd.DataFrame:
    """
    Extract downstream gene expression shifts (Delta X = simulated - imputed)
    for targets of the perturbed transcription factor.

    Parameters
    ----------
    oracle : co.Oracle
        Oracle object after simulation.
    tf : str
        Perturbed transcription factor name.
    links_df : pd.DataFrame, optional
        GRN links DataFrame to look up targets of the TF.
    target_genes : list of str, optional
        Explicit list of downstream targets to extract.
    cluster : str, optional
        Cluster name to subset cells.

    Returns
    -------
    pd.DataFrame
        Table of target genes with mean WT expression, mean simulated expression, and Delta X.
    """
    cols = [
        "tf",
        "target",
        "cluster",
        "mean_wt_expr",
        "mean_sim_expr",
        "mean_delta_x",
        "abs_delta_x",
        "direction",
    ]
    if not hasattr(oracle, "adata"):
        return pd.DataFrame(columns=cols)

    adata = oracle.adata
    if "simulated_count" not in adata.layers or "imputed_count" not in adata.layers:
        return pd.DataFrame(columns=cols)

    # Determine targets to analyze
    targets_to_check = set()
    if target_genes is not None:
        targets_to_check.update(target_genes)

    if links_df is not None and not links_df.empty:
        src_col = "source" if "source" in links_df.columns else links_df.columns[0]
        tgt_col = "target" if "target" in links_df.columns else links_df.columns[1]
        c_links = links_df[links_df[src_col] == tf]
        if cluster is not None and "cluster" in c_links.columns:
            c_links = c_links[c_links["cluster"].astype(str) == str(cluster)]
        targets_to_check.update(c_links[tgt_col].dropna().unique())

    if not targets_to_check:
        # Fallback: check all genes in adata
        targets_to_check.update([g for g in adata.var_names if g != tf][:50])

    # Subset cells by cluster if specified
    mask = np.ones(adata.n_obs, dtype=bool)
    c_str = str(cluster) if cluster is not None else "All"
    if cluster is not None and hasattr(oracle, "cluster_column_name"):
        c_col = getattr(oracle, "cluster_column_name", "cluster")
        if c_col in adata.obs.columns:
            mask = adata.obs[c_col].astype(str) == c_str

    if mask.sum() == 0:
        mask = np.ones(adata.n_obs, dtype=bool)

    imputed = adata.layers["imputed_count"][mask]
    simulated = adata.layers["simulated_count"][mask]

    # Convert sparse matrices to dense if needed
    if hasattr(imputed, "toarray"):
        imputed = imputed.toarray()
    if hasattr(simulated, "toarray"):
        simulated = simulated.toarray()

    rows = []
    for tgt in sorted(list(targets_to_check)):
        if tgt not in adata.var_names:
            continue
        idx = adata.var_names.get_loc(tgt)
        wt_vals = imputed[:, idx]
        sim_vals = simulated[:, idx]

        mean_wt = float(np.mean(wt_vals))
        mean_sim = float(np.mean(sim_vals))
        delta_x = mean_sim - mean_wt
        direction = "upregulated" if delta_x > 0 else "downregulated"

        rows.append({
            "tf": tf,
            "target": tgt,
            "cluster": c_str,
            "mean_wt_expr": mean_wt,
            "mean_sim_expr": mean_sim,
            "mean_delta_x": delta_x,
            "abs_delta_x": abs(delta_x),
            "direction": direction,
        })

    if not rows:
        return pd.DataFrame(columns=cols)

    df_res = pd.DataFrame(rows)
    return df_res.sort_values("abs_delta_x", ascending=False).reset_index(drop=True)
