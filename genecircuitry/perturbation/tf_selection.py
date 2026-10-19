"""
Transcription Factor (TF) Selection Module for In-Silico Perturbation
======================================================================

Selects top candidate TFs for perturbation by evaluating differences in GRN network
centrality (degree, betweenness, eigenvector) between compared clusters and overlaying
Hotspot co-expression module activity and downstream target enrichment.
"""

from typing import Optional, Dict, Any, List, Tuple, Union, Set
import numpy as np
import pandas as pd
from anndata import AnnData
from itertools import combinations

from ..logging_utils import log_info, log_warning


def compute_cluster_pair_centrality_diff(
    score_df: pd.DataFrame,
    cluster_a: str,
    cluster_b: str,
    cluster_key: str = "cluster",
) -> pd.DataFrame:
    """
    Compute differential centrality metrics between two clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        Centrality scores DataFrame from CellOracle (containing degree, betweenness,
        eigenvector centrality).
    cluster_a : str
        First cluster name.
    cluster_b : str
        Second cluster name.
    cluster_key : str, default="cluster"
        Column in score_df indicating cluster assignment.

    Returns
    -------
    pd.DataFrame
        DataFrame with delta centrality metrics for all TFs present in either cluster.
    """
    if score_df is None or score_df.empty:
        return pd.DataFrame()

    df = score_df.copy()
    gene_col = "gene" if "gene" in df.columns else (df.index.name or "gene")
    if gene_col not in df.columns:
        df = df.reset_index()
        gene_col = df.columns[0]

    if cluster_key not in df.columns:
        if "cluster" in df.columns:
            cluster_key = "cluster"
        elif "stratification" in df.columns:
            cluster_key = "stratification"
        else:
            return pd.DataFrame()

    df[cluster_key] = df[cluster_key].astype(str)
    c_a_str = str(cluster_a)
    c_b_str = str(cluster_b)

    df_a = df[df[cluster_key] == c_a_str].set_index(gene_col)
    df_b = df[df[cluster_key] == c_b_str].set_index(gene_col)

    all_tfs = sorted(list(set(df_a.index).union(set(df_b.index))))
    if not all_tfs:
        return pd.DataFrame()

    # Identify metric columns
    deg_col = next((c for c in ["degree_centrality_all", "degree_centrality"] if c in df.columns), None)
    bet_col = next((c for c in ["betweenness_centrality"] if c in df.columns), None)
    eig_col = next((c for c in ["eigenvector_centrality"] if c in df.columns), None)

    rows = []
    for tf in all_tfs:
        row_a = df_a.loc[tf] if tf in df_a.index else None
        row_b = df_b.loc[tf] if tf in df_b.index else None

        deg_a = float(row_a[deg_col]) if row_a is not None and deg_col else 0.0
        deg_b = float(row_b[deg_col]) if row_b is not None and deg_col else 0.0
        delta_deg = deg_a - deg_b

        bet_a = float(row_a[bet_col]) if row_a is not None and bet_col else 0.0
        bet_b = float(row_b[bet_col]) if row_b is not None and bet_col else 0.0
        delta_bet = bet_a - bet_b

        eig_a = float(row_a[eig_col]) if row_a is not None and eig_col else 0.0
        eig_b = float(row_b[eig_col]) if row_b is not None and eig_col else 0.0
        delta_eig = eig_a - eig_b

        # Composite difference: weighted average of deltas
        deltas = [delta_deg, delta_bet, delta_eig]
        composite_diff = float(np.mean(deltas))

        rows.append({
            "gene": tf,
            "cluster_a": c_a_str,
            "cluster_b": c_b_str,
            "degree_a": deg_a,
            "degree_b": deg_b,
            "delta_degree": delta_deg,
            "betweenness_a": bet_a,
            "betweenness_b": bet_b,
            "delta_betweenness": delta_bet,
            "eigenvector_a": eig_a,
            "eigenvector_b": eig_b,
            "delta_eigenvector": delta_eig,
            "composite_centrality_diff": composite_diff,
        })

    return pd.DataFrame(rows)


def compute_tf_hotspot_enrichment(
    links_df: Optional[pd.DataFrame] = None,
    module_genes: Optional[Dict[str, Union[Set[str], List[str]]]] = None,
    cluster: Optional[str] = None,
) -> pd.DataFrame:
    """
    Map TFs and their downstream targets to Hotspot co-expression modules.

    Parameters
    ----------
    links_df : pd.DataFrame, optional
        CellOracle GRN links DataFrame with 'source', 'target', and optionally 'cluster'.
    module_genes : dict, optional
        Dictionary mapping module name to list/set of gene names.
    cluster : str, optional
        Cluster name to filter links by.

    Returns
    -------
    pd.DataFrame
        Table summarizing TF module membership, top target module, and target overlap ratio.
    """
    cols = [
        "tf",
        "tf_in_modules",
        "n_targets",
        "top_target_module",
        "target_overlap_count",
        "target_overlap_ratio",
        "targets_in_module",
    ]
    if links_df is None or links_df.empty or not module_genes:
        return pd.DataFrame(columns=cols)

    df = links_df.copy()
    if cluster is not None and "cluster" in df.columns:
        df = df[df["cluster"].astype(str) == str(cluster)]

    if df.empty:
        return pd.DataFrame(columns=cols)

    src_col = "source" if "source" in df.columns else df.columns[0]
    tgt_col = "target" if "target" in df.columns else df.columns[1]

    # Convert module_genes values to sets
    mod_sets = {str(k): set(v) for k, v in module_genes.items()}

    rows = []
    for tf, tf_group in df.groupby(src_col):
        targets = set(tf_group[tgt_col].dropna().unique())
        n_targets = len(targets)

        # Check which modules contain the TF itself
        tf_mods = [mod for mod, gset in mod_sets.items() if tf in gset]
        tf_in_modules_str = ", ".join(tf_mods) if tf_mods else "None"

        # Check target overlap across modules
        top_mod = "None"
        max_overlap = 0
        top_targets_in_mod: List[str] = []

        if n_targets > 0:
            for mod, gset in mod_sets.items():
                common = targets.intersection(gset)
                if len(common) > max_overlap:
                    max_overlap = len(common)
                    top_mod = mod
                    top_targets_in_mod = list(common)

        overlap_ratio = (max_overlap / n_targets) if n_targets > 0 else 0.0

        rows.append({
            "tf": tf,
            "tf_in_modules": tf_in_modules_str,
            "n_targets": n_targets,
            "top_target_module": top_mod,
            "target_overlap_count": max_overlap,
            "target_overlap_ratio": overlap_ratio,
            "targets_in_module": ", ".join(top_targets_in_mod[:8]),
        })

    return pd.DataFrame(rows)


def select_top_perturbed_tfs(
    score_df: pd.DataFrame,
    links_df: Optional[pd.DataFrame] = None,
    module_genes: Optional[Dict[str, Union[Set[str], List[str]]]] = None,
    adata: Optional[AnnData] = None,
    compare_clusters: Optional[Union[List[Tuple[str, str]], str]] = None,
    cluster_key: str = "cluster",
    top_n_tfs: int = 3,
    mode: str = "directional",
    weight_centrality: float = 0.5,
    weight_module: float = 0.3,
    weight_expression: float = 0.2,
) -> pd.DataFrame:
    """
    Select top candidate TFs for perturbation by integrating differential centrality
    and Hotspot module overlap, producing a selection explanation table.

    Parameters
    ----------
    score_df : pd.DataFrame
        Centrality scores DataFrame.
    links_df : pd.DataFrame, optional
        CellOracle GRN links DataFrame.
    module_genes : dict, optional
        Hotspot module gene sets.
    adata : AnnData, optional
        Expression AnnData to compute log2FC between clusters.
    compare_clusters : list of tuple or str, optional
        Cluster pairs to compare (e.g. [('ClusterA', 'ClusterB')] or 'ClusterA,ClusterB').
        If None, automatically extracts unique cluster combinations.
    cluster_key : str, default="cluster"
        Column in score_df or adata.obs for cluster assignments.
    top_n_tfs : int, default=3
        Number of top candidate TFs to pick per direction.
    mode : str, default="directional"
        Perturbation selection mode: 'directional' (KO for upregulated, OE for downregulated),
        'ko' (KO only), 'oe' (OE only), or 'both'.

    Returns
    -------
    pd.DataFrame
        Table of selected TFs with composite scores and selection rationale reasons.
    """
    if score_df is None or score_df.empty:
        return pd.DataFrame()

    # Resolve cluster_key column in score_df
    actual_cluster_key = cluster_key
    if actual_cluster_key not in score_df.columns:
        if "cluster" in score_df.columns:
            actual_cluster_key = "cluster"
        elif "stratification" in score_df.columns:
            actual_cluster_key = "stratification"

    available_clusters = score_df[actual_cluster_key].astype(str).unique().tolist()
    if len(available_clusters) < 2:
        log_warning("TFSelection", f"Fewer than 2 clusters found in score_df: {available_clusters}")
        return pd.DataFrame()

    # Parse compare_clusters
    pairs: List[Tuple[str, str]] = []
    if compare_clusters is not None:
        if isinstance(compare_clusters, str):
            # Format: "C1,C2;C3,C4" or "C1,C2"
            raw_pairs = [p.strip() for p in compare_clusters.split(";") if p.strip()]
            for rp in raw_pairs:
                parts = [c.strip() for c in rp.split(",") if c.strip()]
                if len(parts) >= 2:
                    pairs.append((parts[0], parts[1]))
        elif isinstance(compare_clusters, (list, tuple)):
            for item in compare_clusters:
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    pairs.append((str(item[0]), str(item[1])))
                elif isinstance(item, str) and "," in item:
                    parts = [c.strip() for c in item.split(",") if c.strip()]
                    if len(parts) >= 2:
                        pairs.append((parts[0], parts[1]))

    if not pairs:
        # Default to all pairwise combinations (limited to first 6 pairs to avoid explosion)
        pairs = list(combinations(available_clusters, 2))[:6]

    selected_rows = []

    for c_a, c_b in pairs:
        if c_a not in available_clusters or c_b not in available_clusters:
            continue

        diff_df = compute_cluster_pair_centrality_diff(
            score_df, cluster_a=c_a, cluster_b=c_b, cluster_key=actual_cluster_key
        )
        if diff_df.empty:
            continue

        # Module enrichment for both clusters
        enrich_a = compute_tf_hotspot_enrichment(links_df, module_genes, cluster=c_a).set_index("tf")
        enrich_b = compute_tf_hotspot_enrichment(links_df, module_genes, cluster=c_b).set_index("tf")

        # Expression log2FC from adata if provided
        expr_means_a = {}
        expr_means_b = {}
        log2fc_dict = {}
        if adata is not None and cluster_key in adata.obs.columns:
            obs_col = adata.obs[cluster_key].astype(str)
            mask_a = obs_col == str(c_a)
            mask_b = obs_col == str(c_b)
            if mask_a.sum() > 0 and mask_b.sum() > 0:
                for gene in diff_df["gene"]:
                    if gene in adata.var_names:
                        val_a = float(np.mean(adata[mask_a, gene].X.toarray() if hasattr(adata[mask_a, gene].X, "toarray") else adata[mask_a, gene].X))
                        val_b = float(np.mean(adata[mask_b, gene].X.toarray() if hasattr(adata[mask_b, gene].X, "toarray") else adata[mask_b, gene].X))
                        expr_means_a[gene] = val_a
                        expr_means_b[gene] = val_b
                        log2fc_dict[gene] = float(np.log2((val_a + 1e-4) / (val_b + 1e-4)))

        # Normalize metrics for scoring
        max_abs_cent = max(diff_df["composite_centrality_diff"].abs().max(), 1e-6)

        # Direction 1: Up in C_A (higher centrality/expression in C_A) -> Drivers of C_A (KO candidates)
        df_up_a = diff_df[diff_df["composite_centrality_diff"] > 0].copy()
        if not df_up_a.empty:
            scores_a = []
            for _, r in df_up_a.iterrows():
                tf = r["gene"]
                norm_cent = r["composite_centrality_diff"] / max_abs_cent
                
                # Module overlap in C_A
                mod_ratio = float(enrich_a.loc[tf, "target_overlap_ratio"]) if tf in enrich_a.index else 0.0
                tf_in_mod = bool(enrich_a.loc[tf, "tf_in_modules"] != "None") if tf in enrich_a.index else False
                mod_bonus = 0.2 if tf_in_mod else 0.0
                mod_score = min(1.0, mod_ratio + mod_bonus)

                # Expression component
                lfc = log2fc_dict.get(tf, 0.0)
                norm_expr = min(1.0, max(-1.0, lfc / 3.0)) if lfc != 0.0 else 0.0

                comp_score = (
                    weight_centrality * norm_cent
                    + weight_module * mod_score
                    + weight_expression * max(0.0, norm_expr)
                )
                scores_a.append(comp_score)

            df_up_a["composite_score"] = scores_a
            top_a = df_up_a.sort_values("composite_score", ascending=False).head(top_n_tfs)

            for _, r in top_a.iterrows():
                tf = r["gene"]
                tf_mod = enrich_a.loc[tf, "tf_in_modules"] if tf in enrich_a.index else "None"
                tgt_mod = enrich_a.loc[tf, "top_target_module"] if tf in enrich_a.index else "None"
                tgt_ratio = enrich_a.loc[tf, "target_overlap_ratio"] if tf in enrich_a.index else 0.0
                tgt_sample = enrich_a.loc[tf, "targets_in_module"] if tf in enrich_a.index else ""
                lfc = log2fc_dict.get(tf, None)

                ptype = "KO" if mode in ["directional", "ko", "both"] else "OE"
                reason_parts = [
                    f"Selected as primary driver of {c_a}:",
                    f"+{r['delta_degree']:.2f} delta degree centrality",
                    f"+{r['delta_eigenvector']:.2f} delta eigenvector centrality",
                ]
                if tf_mod != "None":
                    reason_parts.append(f"member of {tf_mod}")
                if tgt_mod != "None" and tgt_ratio > 0:
                    reason_parts.append(f"{tgt_ratio*100:.1f}% of downstream targets enrich {tgt_mod}")
                if lfc is not None:
                    reason_parts.append(f"log2FC={lfc:+.2f}")

                reason = ", ".join(reason_parts) + "."

                selected_rows.append({
                    "tf": tf,
                    "comparison": f"{c_a} vs {c_b}",
                    "cluster_high": c_a,
                    "cluster_low": c_b,
                    "perturbation_type": ptype,
                    "delta_degree": r["delta_degree"],
                    "delta_betweenness": r["delta_betweenness"],
                    "delta_eigenvector": r["delta_eigenvector"],
                    "composite_centrality_diff": r["composite_centrality_diff"],
                    "tf_hotspot_module": tf_mod,
                    "top_target_module": tgt_mod,
                    "target_module_overlap_ratio": tgt_ratio,
                    "sample_target_genes": tgt_sample,
                    "mean_expression_high": expr_means_a.get(tf, None),
                    "mean_expression_low": expr_means_b.get(tf, None),
                    "log2fc": lfc,
                    "composite_score": r["composite_score"],
                    "selection_reason": reason,
                })

        # Direction 2: Up in C_B (negative delta in diff_df) -> Drivers of C_B
        df_up_b = diff_df[diff_df["composite_centrality_diff"] < 0].copy()
        if not df_up_b.empty:
            scores_b = []
            for _, r in df_up_b.iterrows():
                tf = r["gene"]
                norm_cent = abs(r["composite_centrality_diff"]) / max_abs_cent
                
                # Module overlap in C_B
                mod_ratio = float(enrich_b.loc[tf, "target_overlap_ratio"]) if tf in enrich_b.index else 0.0
                tf_in_mod = bool(enrich_b.loc[tf, "tf_in_modules"] != "None") if tf in enrich_b.index else False
                mod_bonus = 0.2 if tf_in_mod else 0.0
                mod_score = min(1.0, mod_ratio + mod_bonus)

                # Expression component (negative LFC means higher in B)
                lfc = log2fc_dict.get(tf, 0.0)
                norm_expr = min(1.0, max(-1.0, -lfc / 3.0)) if lfc != 0.0 else 0.0

                comp_score = (
                    weight_centrality * norm_cent
                    + weight_module * mod_score
                    + weight_expression * max(0.0, norm_expr)
                )
                scores_b.append(comp_score)

            df_up_b["composite_score"] = scores_b
            top_b = df_up_b.sort_values("composite_score", ascending=False).head(top_n_tfs)

            for _, r in top_b.iterrows():
                tf = r["gene"]
                tf_mod = enrich_b.loc[tf, "tf_in_modules"] if tf in enrich_b.index else "None"
                tgt_mod = enrich_b.loc[tf, "top_target_module"] if tf in enrich_b.index else "None"
                tgt_ratio = enrich_b.loc[tf, "target_overlap_ratio"] if tf in enrich_b.index else 0.0
                tgt_sample = enrich_b.loc[tf, "targets_in_module"] if tf in enrich_b.index else ""
                lfc = log2fc_dict.get(tf, None)

                ptype = "KO" if mode in ["directional", "ko", "both"] else "OE"
                reason_parts = [
                    f"Selected as primary driver of {c_b}:",
                    f"{abs(r['delta_degree']):.2f} higher degree centrality in {c_b}",
                    f"{abs(r['delta_eigenvector']):.2f} higher eigenvector centrality in {c_b}",
                ]
                if tf_mod != "None":
                    reason_parts.append(f"member of {tf_mod}")
                if tgt_mod != "None" and tgt_ratio > 0:
                    reason_parts.append(f"{tgt_ratio*100:.1f}% of downstream targets enrich {tgt_mod}")
                if lfc is not None:
                    reason_parts.append(f"log2FC={-lfc:+.2f} relative to {c_b}")

                reason = ", ".join(reason_parts) + "."

                selected_rows.append({
                    "tf": tf,
                    "comparison": f"{c_a} vs {c_b}",
                    "cluster_high": c_b,
                    "cluster_low": c_a,
                    "perturbation_type": ptype,
                    "delta_degree": -r["delta_degree"],
                    "delta_betweenness": -r["delta_betweenness"],
                    "delta_eigenvector": -r["delta_eigenvector"],
                    "composite_centrality_diff": abs(r["composite_centrality_diff"]),
                    "tf_hotspot_module": tf_mod,
                    "top_target_module": tgt_mod,
                    "target_module_overlap_ratio": tgt_ratio,
                    "sample_target_genes": tgt_sample,
                    "mean_expression_high": expr_means_b.get(tf, None),
                    "mean_expression_low": expr_means_a.get(tf, None),
                    "log2fc": -lfc if lfc is not None else None,
                    "composite_score": r["composite_score"],
                    "selection_reason": reason,
                })

    if not selected_rows:
        return pd.DataFrame()

    return pd.DataFrame(selected_rows)
