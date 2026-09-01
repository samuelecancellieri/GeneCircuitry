"""
Comparative Analysis Module for GeneCircuitry
=============================================

Provides unified cross-cluster and cross-stratification comparative analysis for:
- Co-expression module activity and functional pathway enrichment
- Transcription Factor (TF) centrality, master vs. group-specific classification
- TF-to-Module regulatory connection mapping
- Differential TF-Target Gene (TF-TG) rewiring analysis
"""

import os
from typing import Optional, Dict, Any, List, Tuple, Union
import numpy as np
import pandas as pd
from anndata import AnnData

from . import config
from .logging_utils import log_error, log_warning, log_info


def _normalize_stratification_results(
    stratification_results: Optional[Union[Dict[str, Any], List[Any], Tuple[Any, ...]]],
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Normalize stratification_results whether passed as a dict, list of dicts, or list of tuples.
    Returns a list of (strat_name, data_dict) tuples.
    """
    if not stratification_results:
        return []

    items = []
    if isinstance(stratification_results, dict):
        for k, v in stratification_results.items():
            if isinstance(v, dict):
                items.append((str(k), v))
            else:
                items.append(
                    (
                        str(k),
                        {
                            "adata": getattr(v, "adata", None),
                            "hotspot": getattr(v, "hotspot", None),
                        },
                    )
                )
    elif isinstance(stratification_results, (list, tuple)):
        for item in stratification_results:
            if isinstance(item, dict):
                name = str(
                    item.get("name", item.get("stratification", f"Group_{len(items)}"))
                )
                items.append((name, item))
            elif isinstance(item, tuple) and len(item) == 2:
                items.append(
                    (
                        str(item[0]),
                        item[1] if isinstance(item[1], dict) else {"adata": item[1]},
                    )
                )
            elif hasattr(item, "name"):
                items.append((str(item.name), {"adata": getattr(item, "adata", None)}))

    return items


def compute_module_activity_matrix(
    adata: Optional[AnnData] = None,
    cluster_key: str = "leiden",
    hotspot_obj: Any = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """
    Compute mean co-expression module activity scores per cluster or stratification.

    Parameters
    ----------
    adata : AnnData, optional
        Annotated data matrix containing module scores in `.obs`.
    cluster_key : str, default="leiden"
        Observation column specifying cell clusters.
    hotspot_obj : Hotspot, optional
        Hotspot object with computed module scores.
    stratification_results : dict or list, optional
        Stratified analysis results (dictionary or list of per-stratification dicts).

    Returns
    -------
    pd.DataFrame
        DataFrame where rows are Module IDs and columns are groups (clusters/stratifications).
    """
    activity_dict: Dict[str, pd.Series] = {}

    # Case 1: Stratified results provided
    strat_items = _normalize_stratification_results(stratification_results)
    if strat_items:
        for strat_name, s_data in strat_items:
            s_adata = (
                s_data.get("adata")
                if isinstance(s_data, dict)
                else getattr(s_data, "adata", None)
            )
            s_hs = (
                s_data.get("hotspot_result") or s_data.get("hotspot")
                if isinstance(s_data, dict)
                else getattr(s_data, "hotspot_result", getattr(s_data, "hotspot", None))
            )

            if s_adata is not None:
                # Find module score columns in s_adata.obs
                mod_cols = [
                    c
                    for c in s_adata.obs.columns
                    if str(c).startswith("Module_")
                    or (
                        isinstance(c, (int, np.integer))
                        and f"Module_{c}" in s_adata.obs.columns
                    )
                ]
                if (
                    not mod_cols
                    and s_hs is not None
                    and hasattr(s_hs, "module_scores")
                    and s_hs.module_scores is not None
                ):
                    # Use hotspot module_scores directly
                    hs_scores = s_hs.module_scores
                    for m_col in hs_scores.columns:
                        col_name = (
                            f"Module_{m_col}"
                            if not str(m_col).startswith("Module_")
                            else str(m_col)
                        )
                        activity_dict[f"{strat_name}"] = hs_scores.mean(axis=0)
                    continue

                if mod_cols:
                    mean_scores = s_adata.obs[mod_cols].mean(axis=0)
                    mean_scores.index = [
                        (
                            str(c).replace("Module_", "Module ")
                            if not str(c).startswith("Module ")
                            else str(c)
                        )
                        for c in mean_scores.index
                    ]
                    activity_dict[str(strat_name)] = mean_scores

        if activity_dict:
            df_activity = pd.DataFrame(activity_dict).fillna(0.0)
            df_activity.index.name = "Module"
            return df_activity

    # Case 2: Single dataset with clusters
    if adata is not None:
        # Check if module columns exist in adata.obs
        mod_cols = [c for c in adata.obs.columns if str(c).startswith("Module_")]

        # If not in obs, check if hotspot_obj has them
        if (
            not mod_cols
            and hotspot_obj is not None
            and hasattr(hotspot_obj, "module_scores")
            and hotspot_obj.module_scores is not None
        ):
            hs_scores = hotspot_obj.module_scores.copy()
            hs_cols = []
            for col in hs_scores.columns:
                cname = (
                    f"Module_{col}" if not str(col).startswith("Module_") else str(col)
                )
                adata.obs[cname] = (
                    hs_scores[col].values if len(hs_scores) == len(adata) else np.nan
                )
                hs_cols.append(cname)
            mod_cols = hs_cols

        if mod_cols and cluster_key in adata.obs.columns:
            # Group by cluster and compute mean module scores
            grouped = adata.obs.groupby(cluster_key, observed=False)[mod_cols].mean().T
            grouped.index = [
                str(c).replace("Module_", "Module ") for c in grouped.index
            ]
            grouped.columns = [
                f"Cluster {c}" if not str(c).startswith("Cluster") else str(c)
                for c in grouped.columns
            ]
            grouped.index.name = "Module"
            return grouped
        elif mod_cols:
            mean_all = adata.obs[mod_cols].mean(axis=0).to_frame(name="All Cells")
            mean_all.index = [
                str(c).replace("Module_", "Module ") for c in mean_all.index
            ]
            mean_all.index.name = "Module"
            return mean_all

    # Fallback empty DataFrame
    return pd.DataFrame(columns=["Module"]).set_index("Module")


def compute_module_pathway_enrichments(
    modules_dict_or_df: Union[Dict[Any, List[str]], pd.DataFrame, pd.Series, Any],
    gene_sets: Optional[List[str]] = None,
    top_n_terms: int = 3,
    pval_cutoff: float = 0.05,
) -> pd.DataFrame:
    """
    Compute functional pathway enrichment for co-expression modules.

    Parameters
    ----------
    modules_dict_or_df : dict or pd.Series or pd.DataFrame or Hotspot
        Mapping of genes to module IDs or Hotspot object.
    gene_sets : list of str, optional
        Gene set libraries for ORA. Defaults to config.ENRICHMENT_GENE_SETS.
    top_n_terms : int, default=3
        Number of top enriched terms to retain per module.
    pval_cutoff : float, default=0.05
        Significance cutoff for Adjusted P-value.

    Returns
    -------
    pd.DataFrame
        Table of top enriched pathways per module.
    """
    if gene_sets is None:
        gene_sets = list(config.ENRICHMENT_GENE_SETS)

    # Extract module -> gene list mapping
    module_genes_map: Dict[str, List[str]] = {}

    if (
        hasattr(modules_dict_or_df, "modules")
        and modules_dict_or_df.modules is not None
    ):
        mod_series = modules_dict_or_df.modules
        for m in mod_series.unique():
            if m != -1 and m != "-1":
                module_genes_map[f"Module {m}"] = mod_series[
                    mod_series == m
                ].index.tolist()
    elif isinstance(modules_dict_or_df, pd.Series):
        for m in modules_dict_or_df.unique():
            if m != -1 and m != "-1":
                module_genes_map[f"Module {m}"] = modules_dict_or_df[
                    modules_dict_or_df == m
                ].index.tolist()
    elif isinstance(modules_dict_or_df, pd.DataFrame):
        if "module" in modules_dict_or_df.columns:
            for m in modules_dict_or_df["module"].unique():
                if m != -1 and m != "-1":
                    module_genes_map[f"Module {m}"] = modules_dict_or_df[
                        modules_dict_or_df["module"] == m
                    ].index.tolist()
    elif isinstance(modules_dict_or_df, dict):
        for k, v in modules_dict_or_df.items():
            mod_name = f"Module {k}" if not str(k).startswith("Module") else str(k)
            if isinstance(v, (list, set, tuple, pd.Index)):
                module_genes_map[mod_name] = list(v)

    if not module_genes_map:
        return pd.DataFrame(
            columns=[
                "module",
                "term",
                "adjusted_p_value",
                "combined_score",
                "overlap_genes",
                "gene_count",
            ]
        )

    # Run enrichment for each module
    rows = []
    try:
        from . import enrichment_analysis as ea

        for mod_name, genes in module_genes_map.items():
            if len(genes) < 3:
                continue
            try:
                enr = ea.gseapy_ora_enrichment_analysis(
                    gene_list=genes,
                    gene_sets=gene_sets,
                    pval_cutoff=pval_cutoff,
                )
                if enr.results is not None and not enr.results.empty:
                    df_res = enr.results.copy()
                    # Clean up column names
                    df_res.columns = [c.replace(" ", "_") for c in df_res.columns]
                    sort_col = (
                        "Adjusted_P-value"
                        if "Adjusted_P-value" in df_res.columns
                        else "P-value"
                    )
                    df_res = df_res.sort_values(sort_col).head(top_n_terms)

                    for _, r in df_res.iterrows():
                        term_clean = (
                            str(r.get("Term", ""))
                            .replace("HALLMARK_", "")
                            .replace("_", " ")
                            .title()
                        )
                        rows.append(
                            {
                                "module": mod_name,
                                "term": term_clean,
                                "adjusted_p_value": float(
                                    r.get("Adjusted_P-value", r.get("P-value", 1.0))
                                ),
                                "combined_score": float(r.get("Combined_Score", 0.0)),
                                "overlap_genes": str(r.get("Genes", "")),
                                "gene_count": len(genes),
                            }
                        )
            except Exception as e:
                log_warning(
                    "ComparativeAnalysis.ModuleEnrichment",
                    f"Enrichment for {mod_name} failed ({type(e).__name__}): {e}",
                )
    except Exception as e:
        log_warning(
            "ComparativeAnalysis.ModuleEnrichment",
            f"Enrichment module unavailable ({type(e).__name__}): {e}",
        )

    if not rows:
        # Fallback with gene count summary
        for mod_name, genes in module_genes_map.items():
            rows.append(
                {
                    "module": mod_name,
                    "term": f"{len(genes)} core genes",
                    "adjusted_p_value": 1.0,
                    "combined_score": 0.0,
                    "overlap_genes": ", ".join(genes[:5]),
                    "gene_count": len(genes),
                }
            )

    return pd.DataFrame(rows)


def compute_tf_centrality_matrix(
    score_df: pd.DataFrame,
    score: str = "degree_centrality_all",
    top_n_tfs: int = 15,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute comparative TF centrality across clusters/stratifications and classify
    regulators into Global Master Regulators vs. Group-Specific Drivers.

    Parameters
    ----------
    score_df : pd.DataFrame
        Centrality scores DataFrame from CellOracle analysis.
    score : str, default="degree_centrality_all"
        Centrality score column to rank TFs by.
    top_n_tfs : int, default=15
        Number of top TFs to select per group.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - pivot_df: Matrix of TF (rows) by Group (columns) centrality scores.
        - tf_summary_df: Summary DataFrame with specificity scores and classifications.
    """
    if score_df is None or score_df.empty:
        empty = pd.DataFrame(columns=["gene"]).set_index("gene")
        return empty, pd.DataFrame()

    # Determine group column
    if "stratification" in score_df.columns and "cluster" in score_df.columns:
        score_df = score_df.copy()
        # If stratification has multiple values, combine
        if score_df["stratification"].nunique() > 1:
            score_df["_group"] = (
                score_df["stratification"].astype(str)
                + " - "
                + score_df["cluster"].astype(str)
            )
        else:
            score_df["_group"] = score_df["cluster"].astype(str)
    elif "cluster" in score_df.columns:
        score_df = score_df.copy()
        score_df["_group"] = score_df["cluster"].astype(str)
    elif "stratification" in score_df.columns:
        score_df = score_df.copy()
        score_df["_group"] = score_df["stratification"].astype(str)
    else:
        score_df = score_df.copy()
        score_df["_group"] = "All"

    # Select fallback score if not found
    if score not in score_df.columns:
        num_cols = [
            c
            for c in score_df.columns
            if pd.api.types.is_numeric_dtype(score_df[c]) and c != "_group"
        ]
        score = num_cols[0] if num_cols else score_df.columns[0]

    gene_col = "gene" if "gene" in score_df.columns else score_df.index.name or "index"
    if gene_col not in score_df.columns:
        score_df = score_df.reset_index()
        gene_col = score_df.columns[0]

    # Select top TFs per group
    top_tf_set = set()
    for grp, grp_df in score_df.groupby("_group"):
        top_in_grp = grp_df.nlargest(top_n_tfs, score)[gene_col].tolist()
        top_tf_set.update(top_in_grp)

    # Filter to selected top TFs and pivot
    filtered_df = score_df[score_df[gene_col].isin(top_tf_set)]
    pivot_df = filtered_df.pivot_table(
        index=gene_col,
        columns="_group",
        values=score,
        aggfunc="mean",
    ).fillna(0.0)

    # Compute specificity and classification
    summary_rows = []
    for tf_name, row in pivot_df.iterrows():
        mean_val = row.mean()
        std_val = row.std() if len(row) > 1 else 0.0
        max_val = row.max()
        max_grp = row.idxmax()
        # Specificity index (Coefficient of variation or max-to-mean ratio)
        cv = std_val / (mean_val + 1e-9) if mean_val > 0 else 0.0
        non_zero_ratio = (row > (0.5 * max_val)).sum() / len(row)

        classification = "Global Master" if non_zero_ratio >= 0.6 else "Group-Specific"

        summary_rows.append(
            {
                "gene": tf_name,
                "mean_centrality": mean_val,
                "max_centrality": max_val,
                "top_group": max_grp,
                "specificity_cv": cv,
                "classification": classification,
            }
        )

    tf_summary_df = pd.DataFrame(summary_rows).sort_values(
        "mean_centrality", ascending=False
    )
    pivot_df.index.name = "TF"

    return pivot_df, tf_summary_df


def compute_tf_to_module_mapping(
    links_df: pd.DataFrame,
    modules_dict_or_df: Union[Dict[Any, List[str]], pd.DataFrame, pd.Series, Any],
    top_tfs: Optional[List[str]] = None,
    top_n_tfs: int = 15,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Map transcription factors to the co-expression modules they regulate based on
    TF-Target Gene overlaps.

    Parameters
    ----------
    links_df : pd.DataFrame
        Regulatory links DataFrame (columns: source, target, coef_abs, etc.).
    modules_dict_or_df : dict or pd.Series or pd.DataFrame or Hotspot
        Mapping of genes to module IDs.
    top_tfs : list of str, optional
        List of specific TFs to map. If None, top TFs by link count are used.
    top_n_tfs : int, default=15
        Number of top TFs to evaluate if top_tfs is None.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - matrix_df: Matrix of TF (rows) by Module (columns) target counts.
        - mapping_df: Long-form summary of TF-module regulatory connections.
    """
    if links_df is None or links_df.empty:
        empty = pd.DataFrame(columns=["TF"]).set_index("TF")
        return empty, pd.DataFrame()

    # Extract module -> gene list mapping
    gene_to_module: Dict[str, str] = {}
    if (
        hasattr(modules_dict_or_df, "modules")
        and modules_dict_or_df.modules is not None
    ):
        for g, m in modules_dict_or_df.modules.items():
            if m != -1 and m != "-1":
                gene_to_module[str(g)] = f"Module {m}"
    elif isinstance(modules_dict_or_df, pd.Series):
        for g, m in modules_dict_or_df.items():
            if m != -1 and m != "-1":
                gene_to_module[str(g)] = f"Module {m}"
    elif (
        isinstance(modules_dict_or_df, pd.DataFrame)
        and "module" in modules_dict_or_df.columns
    ):
        for g, row in modules_dict_or_df.iterrows():
            m = row["module"]
            if m != -1 and m != "-1":
                gene_to_module[str(g)] = f"Module {m}"
    elif isinstance(modules_dict_or_df, dict):
        for m, genes in modules_dict_or_df.items():
            mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
            for g in genes:
                gene_to_module[str(g)] = mod_name

    if not gene_to_module:
        empty = pd.DataFrame(columns=["TF"]).set_index("TF")
        return empty, pd.DataFrame()

    # Identify TFs to map
    if top_tfs is None:
        tf_counts = links_df["source"].value_counts()
        top_tfs = tf_counts.head(top_n_tfs).index.tolist()

    # Filter links to selected TFs
    filtered_links = links_df[links_df["source"].isin(top_tfs)].copy()
    filtered_links["module"] = filtered_links["target"].map(gene_to_module)

    # Count targets per TF per module
    mapped_links = filtered_links.dropna(subset=["module"])

    if mapped_links.empty:
        empty = pd.DataFrame(
            index=top_tfs, columns=sorted(list(set(gene_to_module.values())))
        ).fillna(0)
        empty.index.name = "TF"
        return empty, pd.DataFrame()

    matrix_df = mapped_links.pivot_table(
        index="source",
        columns="module",
        values="target",
        aggfunc="nunique",
        fill_value=0,
    )
    matrix_df.index.name = "TF"

    # Compute summary list
    summary_rows = []
    for tf in top_tfs:
        tf_targets = set(filtered_links[filtered_links["source"] == tf]["target"])
        total_targets = len(tf_targets)
        if total_targets == 0:
            continue
        for mod in matrix_df.columns:
            n_in_mod = int(matrix_df.loc[tf, mod]) if tf in matrix_df.index else 0
            if n_in_mod > 0:
                summary_rows.append(
                    {
                        "tf": tf,
                        "module": mod,
                        "target_count": n_in_mod,
                        "total_tf_targets": total_targets,
                        "module_target_ratio": n_in_mod / total_targets,
                    }
                )

    mapping_df = pd.DataFrame(summary_rows).sort_values("target_count", ascending=False)
    return matrix_df, mapping_df


def compute_differential_tf_targets(
    links_df: pd.DataFrame,
    top_tfs: Optional[List[str]] = None,
    top_n_tfs: int = 10,
) -> pd.DataFrame:
    """
    Compare downstream target genes of top TFs across clusters or stratifications to
    identify conserved vs. condition-specific regulatory connections.

    Parameters
    ----------
    links_df : pd.DataFrame
        Regulatory links DataFrame with 'source', 'target', and 'cluster' / 'stratification'.
    top_tfs : list of str, optional
        Specific TFs to compare. If None, top TFs across groups are selected.
    top_n_tfs : int, default=10
        Number of top TFs to evaluate.

    Returns
    -------
    pd.DataFrame
        Summary table comparing target conservation across groups for each TF.
    """
    if links_df is None or links_df.empty:
        return pd.DataFrame(
            columns=[
                "tf",
                "group_count",
                "shared_targets_count",
                "specific_targets_count",
                "conservation_ratio",
            ]
        )

    # Determine group column
    if "stratification" in links_df.columns and "cluster" in links_df.columns:
        if links_df["stratification"].nunique() > 1:
            links_df = links_df.copy()
            links_df["_group"] = (
                links_df["stratification"].astype(str)
                + " - "
                + links_df["cluster"].astype(str)
            )
        else:
            links_df = links_df.copy()
            links_df["_group"] = links_df["cluster"].astype(str)
    elif "cluster" in links_df.columns:
        links_df = links_df.copy()
        links_df["_group"] = links_df["cluster"].astype(str)
    elif "stratification" in links_df.columns:
        links_df = links_df.copy()
        links_df["_group"] = links_df["stratification"].astype(str)
    else:
        return pd.DataFrame(
            columns=[
                "tf",
                "group_count",
                "shared_targets_count",
                "specific_targets_count",
                "conservation_ratio",
            ]
        )

    groups = links_df["_group"].unique()
    if len(groups) < 2:
        return pd.DataFrame(
            columns=[
                "tf",
                "group_count",
                "shared_targets_count",
                "specific_targets_count",
                "conservation_ratio",
            ]
        )

    if top_tfs is None:
        top_tfs = links_df["source"].value_counts().head(top_n_tfs).index.tolist()

    rows = []
    for tf in top_tfs:
        tf_links = links_df[links_df["source"] == tf]
        group_targets: Dict[str, set] = {}
        all_targets: set = set()

        for grp in groups:
            targets = set(tf_links[tf_links["_group"] == grp]["target"])
            if targets:
                group_targets[grp] = targets
                all_targets.update(targets)

        if len(group_targets) >= 2:
            # Intersection across groups
            shared = set.intersection(*group_targets.values())
            # Targets present in only one group
            specific = set()
            for grp, targets in group_targets.items():
                other_targets = set.union(
                    *[t for g, t in group_targets.items() if g != grp]
                )
                specific.update(targets - other_targets)

            rows.append(
                {
                    "tf": tf,
                    "group_count": len(group_targets),
                    "total_targets": len(all_targets),
                    "shared_targets_count": len(shared),
                    "specific_targets_count": len(specific),
                    "conservation_ratio": len(shared) / max(1, len(all_targets)),
                    "sample_shared_targets": ", ".join(list(shared)[:5]),
                }
            )

    return pd.DataFrame(rows).sort_values("conservation_ratio", ascending=False)


def run_comparative_analysis(
    adata: Optional[AnnData] = None,
    score_df: Optional[pd.DataFrame] = None,
    links_df: Optional[pd.DataFrame] = None,
    hotspot_obj: Any = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
    cluster_key: str = "leiden",
    output_dir: Optional[str] = None,
    gene_sets: Optional[List[str]] = None,
    save_tables: bool = True,
) -> Dict[str, Any]:
    """
    Master orchestrator for comparative analysis.

    Executes all cross-group comparative modules, generates processed tables,
    and optionally saves them to `<output_dir>/comparative/`.

    Parameters
    ----------
    adata : AnnData, optional
        Annotated data matrix.
    score_df : pd.DataFrame, optional
        Centrality scores DataFrame.
    links_df : pd.DataFrame, optional
        GRN regulatory links DataFrame.
    hotspot_obj : Hotspot, optional
        Hotspot object.
    stratification_results : dict or list, optional
        Stratified analysis results.
    cluster_key : str, default="leiden"
        Cluster key in `adata.obs`.
    output_dir : str, optional
        Base output directory. Defaults to config.OUTPUT_DIR.
    gene_sets : list of str, optional
        Gene set libraries for enrichment.
    save_tables : bool, default=True
        Whether to save CSV tables to disk.

    Returns
    -------
    dict
        Dictionary of computed comparative data structures.
    """
    if output_dir is None:
        output_dir = config.OUTPUT_DIR

    comp_dir = os.path.join(output_dir, "comparative")
    if save_tables:
        os.makedirs(comp_dir, exist_ok=True)

    log_info(
        "ComparativeAnalysis",
        "Running comparative analysis across clusters/stratifications...",
    )
    print("\nRunning Comparative Analysis...")

    strat_items = _normalize_stratification_results(stratification_results)

    # In stratified mode, auto-gather missing inputs from stratification_results
    if strat_items:
        if hotspot_obj is None:
            for _, s_data in strat_items:
                s_hs = s_data.get("hotspot_result") or s_data.get("hotspot")
                if (
                    s_hs is not None
                    and hasattr(s_hs, "modules")
                    and s_hs.modules is not None
                ):
                    hotspot_obj = s_hs
                    break

        if links_df is None:
            collected_links = []
            for s_name, s_data in strat_items:
                s_out = s_data.get("output_dir")
                if s_out:
                    links_path = os.path.join(
                        s_out, "celloracle", "grn_filtered_links.pkl"
                    )
                    if os.path.exists(links_path):
                        try:
                            import pickle

                            with open(links_path, "rb") as f:
                                df_l = pickle.load(f)
                                if isinstance(df_l, pd.DataFrame):
                                    df_l = df_l.copy()
                                    df_l["stratification"] = str(s_name)
                                    collected_links.append(df_l)
                        except Exception:
                            pass
            if collected_links:
                links_df = pd.concat(collected_links, ignore_index=True)

        if score_df is None:
            collected_scores = []
            for s_name, s_data in strat_items:
                s_out = s_data.get("output_dir")
                if s_out:
                    score_path = os.path.join(
                        s_out, "celloracle", "grn_merged_scores.csv"
                    )
                    if os.path.exists(score_path):
                        try:
                            df_s = pd.read_csv(score_path)
                            df_s["stratification"] = str(s_name)
                            collected_scores.append(df_s)
                        except Exception:
                            pass
            if collected_scores:
                score_df = pd.concat(collected_scores, ignore_index=True)

    results: Dict[str, Any] = {}

    # 1. Module Activity Matrix
    try:
        activity_df = compute_module_activity_matrix(
            adata=adata,
            cluster_key=cluster_key,
            hotspot_obj=hotspot_obj,
            stratification_results=stratification_results,
        )
        results["module_activity"] = activity_df
        if save_tables and not activity_df.empty:
            activity_df.to_csv(os.path.join(comp_dir, "module_activity_matrix.csv"))
            print("  ✓ Computed module activity matrix")
    except Exception as e:
        log_error("ComparativeAnalysis.ModuleActivity", e)
        results["module_activity"] = pd.DataFrame()

    # 2. Module Pathway Enrichment
    try:
        mod_source = (
            hotspot_obj
            if hotspot_obj is not None
            else (adata.obs if adata is not None else None)
        )
        enrichment_df = compute_module_pathway_enrichments(
            modules_dict_or_df=mod_source,
            gene_sets=gene_sets,
        )
        results["module_enrichment"] = enrichment_df
        if save_tables and not enrichment_df.empty:
            enrichment_df.to_csv(
                os.path.join(comp_dir, "module_pathway_enrichment.csv"), index=False
            )
            print("  ✓ Computed module pathway enrichments")
    except Exception as e:
        log_error("ComparativeAnalysis.ModuleEnrichment", e)
        results["module_enrichment"] = pd.DataFrame()

    # 3. TF Centrality Matrix
    try:
        tf_pivot, tf_summary = compute_tf_centrality_matrix(score_df=score_df)
        results["tf_centrality"] = tf_pivot
        results["tf_summary"] = tf_summary
        if save_tables and not tf_pivot.empty:
            tf_pivot.to_csv(os.path.join(comp_dir, "tf_centrality_matrix.csv"))
            tf_summary.to_csv(
                os.path.join(comp_dir, "tf_specificity_summary.csv"), index=False
            )
            print("  ✓ Computed TF centrality & specificity matrix")
    except Exception as e:
        log_error("ComparativeAnalysis.TFCentrality", e)
        results["tf_centrality"] = pd.DataFrame()
        results["tf_summary"] = pd.DataFrame()

    # 4. TF-to-Module Mapping
    try:
        mod_source = (
            hotspot_obj
            if hotspot_obj is not None
            else (adata.obs if adata is not None else None)
        )
        tf_mod_matrix, tf_mod_summary = compute_tf_to_module_mapping(
            links_df=links_df,
            modules_dict_or_df=mod_source,
        )
        results["tf_to_module_matrix"] = tf_mod_matrix
        results["tf_to_module_summary"] = tf_mod_summary
        if save_tables and not tf_mod_matrix.empty:
            tf_mod_matrix.to_csv(os.path.join(comp_dir, "tf_to_module_matrix.csv"))
            tf_mod_summary.to_csv(
                os.path.join(comp_dir, "tf_to_module_mapping.csv"), index=False
            )
            print("  ✓ Computed TF-to-Module regulatory mapping")
    except Exception as e:
        log_error("ComparativeAnalysis.TFToModule", e)
        results["tf_to_module_matrix"] = pd.DataFrame()
        results["tf_to_module_summary"] = pd.DataFrame()

    # 5. Differential TF Targets
    try:
        diff_targets_df = compute_differential_tf_targets(links_df=links_df)
        results["differential_tf_targets"] = diff_targets_df
        if save_tables and not diff_targets_df.empty:
            diff_targets_df.to_csv(
                os.path.join(comp_dir, "differential_tf_targets.csv"), index=False
            )
            print("  ✓ Computed differential TF target rewiring")
    except Exception as e:
        log_error("ComparativeAnalysis.DiffTargets", e)
        results["differential_tf_targets"] = pd.DataFrame()

    print("  Comparative analysis complete.")
    return results
