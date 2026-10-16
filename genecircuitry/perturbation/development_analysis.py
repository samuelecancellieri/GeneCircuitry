"""
Developmental Trajectory and Cluster Transition Analysis Module
===============================================================

Calculates:
1. Inner products (Perturbation Score - PS) comparing simulation vectors against
   developmental pseudotime gradients (following the official CellOracle tutorial notebook).
2. Cluster-to-cluster centroid displacement scores and Markov transition probabilities.
"""

from typing import Optional, Dict, Any, List, Tuple, Union
import numpy as np
import pandas as pd
from anndata import AnnData

from ..logging_utils import log_info, log_warning, log_error


def compute_cluster_displacement_score(
    oracle: Any,
    cluster_source: str,
    cluster_target: str,
    cluster_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Project cell perturbation displacement vectors onto the unit direction vector
    between cluster centroids to assess whether cells shift toward the target cluster.

    Parameters
    ----------
    oracle : co.Oracle
        Oracle object containing 2D embedding and delta_embedding vectors.
    cluster_source : str
        Source cluster (e.g. cluster where TF is normally highly active).
    cluster_target : str
        Target compared cluster (e.g. cluster towards which shift is hypothesized).
    cluster_key : str, optional
        Column in adata.obs for cluster assignments.

    Returns
    -------
    pd.DataFrame
        Table with per-cell displacement projection scores and shift magnitudes.
    """
    cols = ["cell_id", "cluster_source", "cluster_target", "displacement_score", "shift_magnitude"]
    if not hasattr(oracle, "adata") or not hasattr(oracle, "delta_embedding"):
        return pd.DataFrame(columns=cols)

    adata = oracle.adata
    c_col = cluster_key or getattr(oracle, "cluster_column_name", "cluster")
    emb_key = getattr(oracle, "embedding_name", "X_draw_graph_fa")

    if c_col not in adata.obs.columns or emb_key not in adata.obsm:
        return pd.DataFrame(columns=cols)

    coords = adata.obsm[emb_key]
    delta_emb = oracle.delta_embedding
    if delta_emb is None or coords is None:
        return pd.DataFrame(columns=cols)

    obs_clusters = adata.obs[c_col].astype(str)
    src_mask = obs_clusters == str(cluster_source)
    tgt_mask = obs_clusters == str(cluster_target)

    if src_mask.sum() == 0 or tgt_mask.sum() == 0:
        return pd.DataFrame(columns=cols)

    # Compute centroids
    centroid_src = np.mean(coords[src_mask], axis=0)
    centroid_tgt = np.mean(coords[tgt_mask], axis=0)
    diff_vec = centroid_tgt - centroid_src
    dist = np.linalg.norm(diff_vec)
    if dist < 1e-6:
        unit_vec = np.zeros_like(diff_vec)
    else:
        unit_vec = diff_vec / dist

    src_indices = np.where(src_mask)[0]
    cell_names = adata.obs_names[src_mask]

    scores = []
    magnitudes = []
    for idx in src_indices:
        v = delta_emb[idx]
        proj = float(np.dot(v, unit_vec))
        mag = float(np.linalg.norm(v))
        scores.append(proj)
        magnitudes.append(mag)

    return pd.DataFrame({
        "cell_id": cell_names,
        "cluster_source": str(cluster_source),
        "cluster_target": str(cluster_target),
        "displacement_score": scores,
        "shift_magnitude": magnitudes,
    })


def compute_perturbation_cluster_transitions(
    oracle: Any,
    cluster_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Extract or calculate Markov chain cell transition probability table across clusters.

    Parameters
    ----------
    oracle : co.Oracle
        Oracle object after simulation.
    cluster_key : str, optional
        Cluster column in adata.obs.

    Returns
    -------
    pd.DataFrame
        Cluster transition probability matrix (rows: origin, columns: destination).
    """
    c_col = cluster_key or getattr(oracle, "cluster_column_name", "cluster")

    if hasattr(oracle, "get_markov_simulation_cell_transition_table"):
        try:
            df = oracle.get_markov_simulation_cell_transition_table(
                cluster_column_name=c_col, return_df=True
            )
            if df is not None and isinstance(df, pd.DataFrame):
                return df
        except Exception as e:
            log_warning("Perturbation.Transitions", f"Could not get transition table: {e}")

    # Fallback: estimate from cluster assignments
    if hasattr(oracle, "adata") and c_col in oracle.adata.obs.columns:
        clusters = sorted(list(oracle.adata.obs[c_col].unique()))
        n = len(clusters)
        # Identity fallback
        return pd.DataFrame(np.eye(n), index=clusters, columns=clusters)

    return pd.DataFrame()


def compute_developmental_flow_inner_product(
    oracle: Any,
    tf: str,
    mode: str = "KO",
    pseudotime_key: str = "pseudotime",
    cluster_key: Optional[str] = None,
) -> Tuple[Optional[Any], Optional[pd.DataFrame]]:
    """
    Calculate inner products comparing perturbation simulation vectors against the
    developmental gradient vector field following the CellOracle tutorial notebook.

    Parameters
    ----------
    oracle : co.Oracle
        Fitted Oracle object with simulated shift and grid arrows.
    tf : str
        Perturbed gene name.
    mode : str, default="KO"
        Perturbation condition name.
    pseudotime_key : str, default="pseudotime"
        Column in adata.obs containing pseudotime or differentiation trajectory values.
    cluster_key : str, optional
        Cluster column name in adata.obs.

    Returns
    -------
    Tuple[Optional[Any], Optional[pd.DataFrame]]
        - dev: Oracle_development_module object (if successful, else None).
        - ip_df: DataFrame of inner product values per cell/cluster (or None).
    """
    if not hasattr(oracle, "adata"):
        return None, None

    adata = oracle.adata

    # Check for pseudotime in adata.obs
    pt_col = pseudotime_key
    if pt_col not in adata.obs.columns:
        # Fallback to common trajectory names
        candidates = ["dpt_pseudotime", "pseudotime", "latent_time", "velocity_pseudotime"]
        found = next((c for c in candidates if c in adata.obs.columns), None)
        if found:
            pt_col = found
        elif "X_diffmap" in adata.obsm:
            # Derive 1D pseudotime from DC1
            adata.obs["_inferred_pseudotime"] = adata.obsm["X_diffmap"][:, 1]
            pt_col = "_inferred_pseudotime"
            log_info("Perturbation.DevFlow", "Using first non-trivial diffusion component as pseudotime")
        else:
            log_info(
                "Perturbation.DevFlow",
                "No pseudotime column found in adata.obs; skipping developmental gradient comparison.",
            )
            return None, None

    try:
        from celloracle.applications import Gradient_calculator, Oracle_development_module
    except ImportError:
        log_warning(
            "Perturbation.DevFlow",
            "CellOracle applications module could not be imported; skipping developmental flow.",
        )
        return None, None

    try:
        log_info(
            "Perturbation.DevFlow",
            f"Computing developmental gradient and inner product for {tf}_{mode}",
        )
        gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key=pt_col)
        gradient.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
        gradient.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=200)

        dev = Oracle_development_module()
        dev.load_differentiation_reference_data(gradient_object=gradient)
        dev.load_perturb_simulation_data(oracle_object=oracle, name=f"{tf}_{mode}")

        dev.calculate_inner_product()
        dev.calculate_digitized_ip()

        ip_df = None
        if hasattr(dev, "inner_product_df"):
            ip_df = dev.inner_product_df.copy()
        elif hasattr(dev, "inner_product"):
            ip_df = pd.DataFrame({
                "cell_id": adata.obs_names,
                "inner_product": dev.inner_product,
            })

        log_info("Perturbation.DevFlow", "Developmental inner product calculation complete")
        return dev, ip_df

    except Exception as e:
        log_warning("Perturbation.DevFlow", f"Error in developmental gradient calculation: {e}")
        return None, None
