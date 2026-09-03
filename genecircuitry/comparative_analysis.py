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
from itertools import combinations

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


def _normalize_links_df(
    links_df: Optional[Union[pd.DataFrame, Dict[str, Any], str]] = None,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """Normalize links_df from DataFrame, dict of cluster DataFrames, pickle file, or search paths."""
    if links_df is not None:
        if isinstance(links_df, dict):
            dfs = []
            for c_name, c_df in links_df.items():
                if isinstance(c_df, pd.DataFrame):
                    df_copy = c_df.copy()
                    df_copy["cluster"] = str(c_name)
                    dfs.append(df_copy)
                elif isinstance(c_df, dict):
                    df_copy = pd.DataFrame(c_df)
                    df_copy["cluster"] = str(c_name)
                    dfs.append(df_copy)
            if dfs:
                links_df = pd.concat(dfs, ignore_index=True)
            else:
                links_df = pd.DataFrame()
        elif isinstance(links_df, str) and os.path.exists(links_df):
            try:
                import pickle
                with open(links_df, "rb") as f:
                    data = pickle.load(f)
                return _normalize_links_df(data, output_dir=output_dir, stratification_results=stratification_results)
            except Exception:
                links_df = pd.DataFrame()

    if links_df is None or (isinstance(links_df, pd.DataFrame) and links_df.empty):
        strat_items = _normalize_stratification_results(stratification_results)
        collected = []
        if strat_items:
            for s_name, s_data in strat_items:
                s_out = s_data.get("output_dir") if isinstance(s_data, dict) else getattr(s_data, "output_dir", None)
                if s_out:
                    p = os.path.join(s_out, "celloracle", "grn_filtered_links.pkl")
                    if os.path.exists(p):
                        try:
                            import pickle
                            with open(p, "rb") as f:
                                l_data = pickle.load(f)
                            if isinstance(l_data, dict):
                                for c_name, c_df in l_data.items():
                                    if isinstance(c_df, pd.DataFrame):
                                        df_c = c_df.copy()
                                        df_c["cluster"] = str(c_name)
                                        df_c["stratification"] = str(s_name)
                                        collected.append(df_c)
                            elif isinstance(l_data, pd.DataFrame):
                                df_c = l_data.copy()
                                df_c["stratification"] = str(s_name)
                                collected.append(df_c)
                        except Exception:
                            pass
        if collected:
            links_df = pd.concat(collected, ignore_index=True)

    if (links_df is None or (isinstance(links_df, pd.DataFrame) and links_df.empty)) and output_dir:
        p = os.path.join(output_dir, "celloracle", "grn_filtered_links.pkl")
        if os.path.exists(p):
            try:
                import pickle
                with open(p, "rb") as f:
                    l_data = pickle.load(f)
                return _normalize_links_df(l_data)
            except Exception:
                pass

    if not isinstance(links_df, pd.DataFrame) or links_df.empty:
        return pd.DataFrame()

    df = links_df.copy()
    if "source" not in df.columns or "target" not in df.columns:
        return pd.DataFrame()

    if "cluster" not in df.columns:
        if "stratification" in df.columns:
            df["cluster"] = df["stratification"].astype(str)
        else:
            df["cluster"] = "0"
    else:
        df["cluster"] = df["cluster"].astype(str)

    if "coef_abs" not in df.columns:
        if "coef_mean" in df.columns:
            df["coef_abs"] = df["coef_mean"].abs()
        else:
            df["coef_abs"] = 1.0

    return df


def _normalize_score_df(
    score_df: Optional[Union[pd.DataFrame, str]] = None,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """Normalize score_df from DataFrame, CSV path, or search paths, ensuring 'gene' and 'cluster' columns exist."""
    if score_df is not None and isinstance(score_df, str) and os.path.exists(score_df):
        try:
            score_df = pd.read_csv(score_df)
        except Exception:
            score_df = None

    if score_df is None or (isinstance(score_df, pd.DataFrame) and score_df.empty):
        strat_items = _normalize_stratification_results(stratification_results)
        collected = []
        if strat_items:
            for s_name, s_data in strat_items:
                s_out = s_data.get("output_dir") if isinstance(s_data, dict) else getattr(s_data, "output_dir", None)
                if s_out:
                    p = os.path.join(s_out, "celloracle", "grn_merged_scores.csv")
                    if os.path.exists(p):
                        try:
                            df_s = pd.read_csv(p)
                            df_s["stratification"] = str(s_name)
                            collected.append(df_s)
                        except Exception:
                            pass
        if collected:
            score_df = pd.concat(collected, ignore_index=True)

    if (score_df is None or (isinstance(score_df, pd.DataFrame) and score_df.empty)) and output_dir:
        for fname in ["total_merged_scores.csv", "grn_merged_scores.csv"]:
            p = os.path.join(output_dir, "celloracle", fname)
            if os.path.exists(p):
                try:
                    score_df = pd.read_csv(p)
                    break
                except Exception:
                    pass

    if not isinstance(score_df, pd.DataFrame) or score_df.empty:
        return pd.DataFrame()

    df = score_df.copy()
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "gene"})
    if "gene" not in df.columns:
        if df.index.name == "gene":
            df = df.reset_index()
        elif len(df) > 0 and isinstance(df.index[0], str):
            df = df.reset_index().rename(columns={"index": "gene"})

    if "cluster" not in df.columns:
        if "stratification" in df.columns:
            df["cluster"] = df["stratification"].astype(str)
        else:
            df["cluster"] = "0"
    else:
        df["cluster"] = df["cluster"].astype(str)

    if "degree_centrality_all" not in df.columns:
        for alt_col in ["degree_centrality", "degree", "betweenness_centrality", "eigenvector_centrality"]:
            if alt_col in df.columns:
                df["degree_centrality_all"] = df[alt_col]
                break
        if "degree_centrality_all" not in df.columns:
            df["degree_centrality_all"] = 0.0

    return df


def _extract_module_genes(
    hotspot_obj: Any = None,
    adata: Optional[AnnData] = None,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> Dict[str, set]:
    """Extract module_name -> set_of_genes mapping from Hotspot object, AnnData, disk CSVs, or stratification results."""
    module_genes: Dict[str, set] = {}

    if hotspot_obj is not None:
        if hasattr(hotspot_obj, "modules") and hotspot_obj.modules is not None:
            mod_series = hotspot_obj.modules
            for m in mod_series.unique():
                if m != -1 and m != "-1" and str(m) != "-1":
                    mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                    module_genes[mod_name] = set(mod_series[mod_series == m].index.tolist())
        elif isinstance(hotspot_obj, dict):
            for m, g_list in hotspot_obj.items():
                if m != -1 and m != "-1" and str(m) != "-1":
                    mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                    module_genes[mod_name] = set(g_list) if isinstance(g_list, (list, set, pd.Series)) else set()

    if not module_genes and stratification_results:
        strat_items = _normalize_stratification_results(stratification_results)
        for _, s_data in strat_items:
            s_hs = s_data.get("hotspot_result") or s_data.get("hotspot") if isinstance(s_data, dict) else getattr(s_data, "hotspot_result", getattr(s_data, "hotspot", None))
            if s_hs is not None and hasattr(s_hs, "modules") and s_hs.modules is not None:
                mod_series = s_hs.modules
                for m in mod_series.unique():
                    if m != -1 and m != "-1" and str(m) != "-1":
                        mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                        module_genes.setdefault(mod_name, set()).update(mod_series[mod_series == m].index.tolist())
            s_out = s_data.get("output_dir") if isinstance(s_data, dict) else getattr(s_data, "output_dir", None)
            if s_out:
                csv_path = os.path.join(s_out, "hotspot", "gene_modules.csv")
                if os.path.exists(csv_path):
                    try:
                        m_df = pd.read_csv(csv_path, index_col=0)
                        col = m_df.columns[0] if len(m_df.columns) > 0 else None
                        if col:
                            for g, m in m_df[col].items():
                                if m != -1 and m != "-1" and str(m) != "-1":
                                    mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                                    module_genes.setdefault(mod_name, set()).add(str(g))
                    except Exception:
                        pass

    if not module_genes and output_dir:
        csv_path = os.path.join(output_dir, "hotspot", "gene_modules.csv")
        if os.path.exists(csv_path):
            try:
                m_df = pd.read_csv(csv_path, index_col=0)
                col = m_df.columns[0] if len(m_df.columns) > 0 else None
                if col:
                    for g, m in m_df[col].items():
                        if m != -1 and m != "-1" and str(m) != "-1":
                            mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                            module_genes.setdefault(mod_name, set()).add(str(g))
            except Exception:
                pass

    if not module_genes and adata is not None:
        if "hotspot_modules" in adata.uns:
            hm = adata.uns["hotspot_modules"]
            if isinstance(hm, dict):
                for m, g_list in hm.items():
                    mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                    module_genes[mod_name] = set(g_list)
        elif "hotspot_module" in adata.var:
            for g, m in adata.var["hotspot_module"].items():
                if m != -1 and m != "-1" and str(m) != "-1":
                    mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                    module_genes.setdefault(mod_name, set()).add(str(g))

    return module_genes


def _extract_per_stratification_modules(
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
    hotspot_obj: Any = None,
    adata: Optional[AnnData] = None,
    output_dir: Optional[str] = None,
) -> Dict[str, Dict[str, set]]:
    """
    Extract module_name -> set_of_genes mapping distinctly per stratification or dataset.
    
    Returns
    -------
    dict: {strat_name: {module_name: set(genes)}}
    """
    strat_modules: Dict[str, Dict[str, set]] = {}

    strat_items = _normalize_stratification_results(stratification_results)
    if strat_items:
        for s_name, s_data in strat_items:
            m_dict: Dict[str, set] = {}
            s_hs = s_data.get("hotspot_result") or s_data.get("hotspot") if isinstance(s_data, dict) else getattr(s_data, "hotspot_result", getattr(s_data, "hotspot", None))
            if s_hs is not None and hasattr(s_hs, "modules") and s_hs.modules is not None:
                mod_series = s_hs.modules
                for m in mod_series.unique():
                    if m != -1 and m != "-1" and str(m) != "-1":
                        mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                        m_dict[mod_name] = set(mod_series[mod_series == m].index.tolist())
            
            s_out = s_data.get("output_dir") if isinstance(s_data, dict) else getattr(s_data, "output_dir", None)
            if not m_dict and s_out:
                csv_path = os.path.join(s_out, "hotspot", "gene_modules.csv")
                if os.path.exists(csv_path):
                    try:
                        m_df = pd.read_csv(csv_path, index_col=0)
                        col = m_df.columns[0] if len(m_df.columns) > 0 else None
                        if col:
                            for g, m in m_df[col].items():
                                if m != -1 and m != "-1" and str(m) != "-1":
                                    mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                                    m_dict.setdefault(mod_name, set()).add(str(g))
                    except Exception:
                        pass
            
            s_adata = s_data.get("adata") if isinstance(s_data, dict) else getattr(s_data, "adata", None)
            if not m_dict and s_adata is not None:
                if "hotspot_modules" in s_adata.uns:
                    hm = s_adata.uns["hotspot_modules"]
                    if isinstance(hm, dict):
                        for m, g_list in hm.items():
                            mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                            m_dict[mod_name] = set(g_list)
                elif "hotspot_module" in s_adata.var:
                    for g, m in s_adata.var["hotspot_module"].items():
                        if m != -1 and m != "-1" and str(m) != "-1":
                            mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                            m_dict.setdefault(mod_name, set()).add(str(g))
                            
            if m_dict:
                strat_modules[str(s_name)] = m_dict

    # Single-dataset fallback
    if not strat_modules:
        single_m = _extract_module_genes(hotspot_obj, adata, output_dir=output_dir)
        if single_m:
            strat_modules["Dataset"] = single_m

    return strat_modules


def _extract_autocorr_results(
    hotspot_obj: Any = None,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> Tuple[Optional[pd.DataFrame], Optional[pd.Series]]:
    """Extract results DataFrame (with FDR) and modules Series from Hotspot object or disk CSVs."""
    results_df = None
    modules_series = None

    if hotspot_obj is not None:
        if hasattr(hotspot_obj, "results") and hotspot_obj.results is not None:
            results_df = hotspot_obj.results
        if hasattr(hotspot_obj, "modules") and hotspot_obj.modules is not None:
            modules_series = hotspot_obj.modules

    if results_df is None and stratification_results:
        strat_items = _normalize_stratification_results(stratification_results)
        for _, s_data in strat_items:
            s_out = s_data.get("output_dir") if isinstance(s_data, dict) else getattr(s_data, "output_dir", None)
            if s_out:
                p = os.path.join(s_out, "hotspot", "autocorrelation_results.csv")
                if not os.path.exists(p):
                    p = os.path.join(s_out, "hotspot", "significant_genes.csv")
                if os.path.exists(p):
                    try:
                        results_df = pd.read_csv(p, index_col=0)
                        if "FDR" not in results_df.columns and "fdr" in results_df.columns:
                            results_df["FDR"] = results_df["fdr"]
                    except Exception:
                        pass
                m_path = os.path.join(s_out, "hotspot", "gene_modules.csv")
                if os.path.exists(m_path):
                    try:
                        m_df = pd.read_csv(m_path, index_col=0)
                        modules_series = m_df.iloc[:, 0]
                    except Exception:
                        pass
                if results_df is not None and modules_series is not None:
                    break

    if results_df is None and output_dir:
        p = os.path.join(output_dir, "hotspot", "autocorrelation_results.csv")
        if not os.path.exists(p):
            p = os.path.join(output_dir, "hotspot", "significant_genes.csv")
        if os.path.exists(p):
            try:
                results_df = pd.read_csv(p, index_col=0)
                if "FDR" not in results_df.columns and "fdr" in results_df.columns:
                    results_df["FDR"] = results_df["fdr"]
            except Exception:
                pass
        m_path = os.path.join(output_dir, "hotspot", "gene_modules.csv")
        if os.path.exists(m_path):
            try:
                m_df = pd.read_csv(m_path, index_col=0)
                modules_series = m_df.iloc[:, 0]
            except Exception:
                pass

    return results_df, modules_series


def compute_module_activity_matrix(
    adata: Optional[AnnData] = None,
    cluster_key: str = "leiden",
    hotspot_obj: Any = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
    output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compute mean co-expression module activity scores per cluster or stratification.
    """
    activity_dict: Dict[str, pd.Series] = {}

    strat_items = _normalize_stratification_results(stratification_results)
    if strat_items:
        for strat_name, s_data in strat_items:
            s_adata = s_data.get("adata") if isinstance(s_data, dict) else getattr(s_data, "adata", None)
            s_hs = s_data.get("hotspot_result") or s_data.get("hotspot") if isinstance(s_data, dict) else getattr(s_data, "hotspot_result", getattr(s_data, "hotspot", None))
            s_out = s_data.get("output_dir") if isinstance(s_data, dict) else getattr(s_data, "output_dir", None)

            # Load module scores if needed
            scores_df = None
            if s_hs is not None and hasattr(s_hs, "module_scores") and s_hs.module_scores is not None:
                scores_df = s_hs.module_scores
            elif s_out:
                csv_p = os.path.join(s_out, "hotspot", "hotspot_module_scores.csv")
                if os.path.exists(csv_p):
                    try:
                        scores_df = pd.read_csv(csv_p, index_col=0)
                    except Exception:
                        pass

            if s_adata is not None:
                # Find or attach module columns
                mod_cols = [c for c in s_adata.obs.columns if str(c).startswith("Module_")]
                if not mod_cols and scores_df is not None:
                    for col in scores_df.columns:
                        cname = f"Module_{col}" if not str(col).startswith("Module_") else str(col)
                        try:
                            s_adata.obs[cname] = scores_df[col].reindex(s_adata.obs.index).values
                            mod_cols.append(cname)
                        except Exception:
                            pass

                # If cluster_key exists, compute per-cluster activity
                if mod_cols and cluster_key in s_adata.obs.columns:
                    for c_val in s_adata.obs[cluster_key].unique():
                        c_cells = s_adata.obs[s_adata.obs[cluster_key] == c_val]
                        if len(c_cells) > 0:
                            mean_s = c_cells[mod_cols].mean(axis=0)
                            mean_s.index = [str(c).replace("Module_", "Module ") for c in mean_s.index]
                            prefix = f"{strat_name} - Cl {c_val}" if len(strat_items) > 1 else f"Cluster {c_val}"
                            activity_dict[prefix] = mean_s
                elif mod_cols:
                    mean_s = s_adata.obs[mod_cols].mean(axis=0)
                    mean_s.index = [str(c).replace("Module_", "Module ") for c in mean_s.index]
                    activity_dict[str(strat_name)] = mean_s
            elif scores_df is not None:
                mean_s = scores_df.mean(axis=0)
                mean_s.index = [f"Module {c}" if not str(c).startswith("Module") else str(c) for c in mean_s.index]
                activity_dict[str(strat_name)] = mean_s

        if activity_dict:
            df_activity = pd.DataFrame(activity_dict).fillna(0.0)
            df_activity.index.name = "Module"
            return df_activity

    # Case 2: Single dataset
    if adata is not None:
        mod_cols = [c for c in adata.obs.columns if str(c).startswith("Module_")]
        if not mod_cols:
            scores_df = None
            if hotspot_obj is not None and hasattr(hotspot_obj, "module_scores") and hotspot_obj.module_scores is not None:
                scores_df = hotspot_obj.module_scores
            elif output_dir:
                csv_p = os.path.join(output_dir, "hotspot", "hotspot_module_scores.csv")
                if os.path.exists(csv_p):
                    try:
                        scores_df = pd.read_csv(csv_p, index_col=0)
                    except Exception:
                        pass
            if scores_df is not None:
                for col in scores_df.columns:
                    cname = f"Module_{col}" if not str(col).startswith("Module_") else str(col)
                    try:
                        adata.obs[cname] = scores_df[col].reindex(adata.obs.index).values
                        mod_cols.append(cname)
                    except Exception:
                        pass

        if mod_cols and cluster_key in adata.obs.columns:
            grouped = adata.obs.groupby(cluster_key, observed=False)[mod_cols].mean().T
            grouped.index = [str(c).replace("Module_", "Module ") for c in grouped.index]
            grouped.columns = [f"Cluster {c}" if not str(c).startswith("Cluster") else str(c) for c in grouped.columns]
            grouped.index.name = "Module"
            return grouped
        elif mod_cols:
            mean_all = adata.obs[mod_cols].mean(axis=0).to_frame(name="All Cells")
            mean_all.index = [str(c).replace("Module_", "Module ") for c in mean_all.index]
            mean_all.index.name = "Module"
            return mean_all

    return pd.DataFrame(columns=["Module"]).set_index("Module")


def compute_module_pathway_enrichments(
    modules_dict_or_df: Union[Dict[Any, List[str]], pd.DataFrame, pd.Series, Any],
    gene_sets: Optional[List[str]] = None,
    top_n_terms: int = 3,
    pval_cutoff: float = 0.05,
    online: Optional[bool] = None,
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
    online : bool, optional
        Whether to use the online Enrichr API. Defaults to config.ENRICHMENT_ONLINE.

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
                    online=online,
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
    links_df: Optional[pd.DataFrame] = None,
    modules_dict_or_df: Union[Dict[Any, List[str]], pd.DataFrame, pd.Series, Any] = None,
    top_tfs: Optional[List[str]] = None,
    top_n_tfs: int = 15,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Map transcription factors to the co-expression modules they regulate based on
    TF-Target Gene overlaps.
    """
    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)
    if links_df.empty:
        empty = pd.DataFrame(columns=["TF"]).set_index("TF")
        return empty, pd.DataFrame()

    # Extract module -> gene list mapping
    module_genes = _extract_module_genes(
        hotspot_obj=modules_dict_or_df if not isinstance(modules_dict_or_df, (dict, pd.DataFrame, pd.Series)) else None,
        output_dir=output_dir,
        stratification_results=stratification_results,
    )
    if not module_genes and isinstance(modules_dict_or_df, dict):
        for m, genes in modules_dict_or_df.items():
            mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
            module_genes[mod_name] = set(genes)
    elif not module_genes and isinstance(modules_dict_or_df, pd.Series):
        for g, m in modules_dict_or_df.items():
            if m != -1 and m != "-1" and str(m) != "-1":
                mod_name = f"Module {m}" if not str(m).startswith("Module") else str(m)
                module_genes.setdefault(mod_name, set()).add(str(g))

    if not module_genes:
        empty = pd.DataFrame(columns=["TF"]).set_index("TF")
        return empty, pd.DataFrame()

    gene_to_module = {}
    for mod_name, genes in module_genes.items():
        for g in genes:
            gene_to_module[str(g)] = mod_name

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
            index=top_tfs, columns=sorted(list(module_genes.keys()))
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
                mod_total_genes = len(module_genes.get(mod, set()))
                summary_rows.append(
                    {
                        "tf": tf,
                        "module": mod,
                        "target_count": n_in_mod,
                        "total_tf_targets": total_targets,
                        "module_target_ratio": n_in_mod / total_targets,
                        "module_coverage_pct": (n_in_mod / mod_total_genes) if mod_total_genes > 0 else 0.0,
                    }
                )

    mapping_df = pd.DataFrame(summary_rows).sort_values("target_count", ascending=False)
    return matrix_df, mapping_df


def compute_differential_tf_targets(
    links_df: Optional[pd.DataFrame] = None,
    top_tfs: Optional[List[str]] = None,
    top_n_tfs: int = 10,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """
    Compare downstream target genes of top TFs across clusters or stratifications to
    identify conserved vs. condition-specific regulatory connections.
    """
    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)
    if links_df.empty:
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
        links_df = links_df.copy()
        if links_df["stratification"].nunique() > 1:
            links_df["_group"] = (
                links_df["stratification"].astype(str)
                + " - "
                + links_df["cluster"].astype(str)
            )
        else:
            links_df["_group"] = links_df["cluster"].astype(str)
    elif "cluster" in links_df.columns:
        links_df = links_df.copy()
        links_df["_group"] = links_df["cluster"].astype(str)
    else:
        return pd.DataFrame()

    groups = links_df["_group"].unique()
    if len(groups) < 2:
        return pd.DataFrame()

    if top_tfs is None:
        tf_counts = links_df.groupby("source", observed=False)["_group"].nunique()
        top_tfs = tf_counts.nlargest(top_n_tfs * 2).index.tolist()

    rows = []
    for tf in top_tfs:
        tf_links = links_df[links_df["source"] == tf]
        group_targets = {}
        all_targets: set = set()

        for grp in groups:
            targets = set(tf_links[tf_links["_group"] == grp]["target"])
            if targets:
                group_targets[grp] = targets
                all_targets.update(targets)

        if len(group_targets) >= 2:
            shared = set.intersection(*group_targets.values())
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

    if not rows:
        return pd.DataFrame(columns=["tf", "group_count", "total_targets", "shared_targets_count", "specific_targets_count", "conservation_ratio", "sample_shared_targets"])
    return pd.DataFrame(rows).sort_values("conservation_ratio", ascending=False)


def compute_module_gene_overlap_matrix(
    hotspot_obj: Any = None,
    links_df: Optional[pd.DataFrame] = None,
    adata: Optional[AnnData] = None,
    cluster_key: str = "leiden",
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute module gene GRN coverage and cross-cluster Jaccard similarity."""
    module_genes = _extract_module_genes(hotspot_obj, adata, output_dir=output_dir, stratification_results=stratification_results)
    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)

    if not module_genes or links_df.empty or "cluster" not in links_df.columns:
        empty1 = pd.DataFrame(columns=["Module"]).set_index("Module")
        empty2 = pd.DataFrame(columns=["Module"]).set_index("Module")
        return empty1, empty2

    clusters = sorted(list(links_df["cluster"].unique()))
    cluster_targets = {}
    for c in clusters:
        cluster_targets[c] = set(links_df[links_df["cluster"] == c]["target"].unique())

    # Build coverage_df
    coverage_rows = {}
    for mod_name, genes in module_genes.items():
        row = {}
        for c in clusters:
            intersection = genes.intersection(cluster_targets[c])
            col_label = f"Cluster {c}" if not str(c).startswith("Cluster") else str(c)
            row[col_label] = len(intersection) / len(genes) if len(genes) > 0 else 0.0
        coverage_rows[mod_name] = row
    coverage_df = pd.DataFrame.from_dict(coverage_rows, orient="index").fillna(0.0)
    coverage_df.index.name = "Module"

    # Build jaccard_df
    jaccard_rows = {}
    for mod_name, genes in module_genes.items():
        row = {}
        for c1, c2 in combinations(clusters, 2):
            g1 = genes.intersection(cluster_targets[c1])
            g2 = genes.intersection(cluster_targets[c2])
            intersection = len(g1.intersection(g2))
            union = len(g1.union(g2))
            row[f"{c1} vs {c2}"] = intersection / union if union > 0 else 0.0
        jaccard_rows[mod_name] = row
    jaccard_df = pd.DataFrame.from_dict(jaccard_rows, orient="index").fillna(0.0)
    jaccard_df.index.name = "Module"

    return coverage_df, jaccard_df


def compute_cross_stratification_module_overlap(
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
    hotspot_obj: Any = None,
    adata: Optional[AnnData] = None,
    output_dir: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute pairwise gene overlap (Jaccard similarity and Overlap coefficient) between
    Hotspot modules across different stratifications or runs to align truly similar modules
    by gene content rather than arbitrary numeric names.

    Returns
    -------
    jaccard_df : pd.DataFrame
        Matrix of Jaccard similarity indices between (strat1, mod1) and (strat2, mod2).
    overlap_coef_df : pd.DataFrame
        Matrix of Overlap/Simpson coefficients |A ∩ B| / min(|A|, |B|).
    alignment_summary_df : pd.DataFrame
        Summary table detailing pairwise overlap metrics, shared gene count, sample genes,
        and alignment status (Conserved, Related, Distinct).
    """
    strat_modules = _extract_per_stratification_modules(
        stratification_results=stratification_results,
        hotspot_obj=hotspot_obj,
        adata=adata,
        output_dir=output_dir,
    )

    if not strat_modules or sum(len(m) for m in strat_modules.values()) == 0:
        empty_df = pd.DataFrame()
        return empty_df, empty_df, empty_df

    all_mod_keys: List[Tuple[str, str, set]] = []
    for s_name, m_dict in strat_modules.items():
        for m_name, genes in m_dict.items():
            if genes:
                all_mod_keys.append((s_name, m_name, genes))

    labels = [f"{s}: {m}" if len(strat_modules) > 1 else m for s, m, _ in all_mod_keys]
    n = len(all_mod_keys)

    jaccard_mat = np.zeros((n, n), dtype=float)
    overlap_mat = np.zeros((n, n), dtype=float)
    summary_rows = []

    for i in range(n):
        s1, m1, g1 = all_mod_keys[i]
        for j in range(n):
            s2, m2, g2 = all_mod_keys[j]
            inter = len(g1.intersection(g2))
            union = len(g1.union(g2))
            min_len = min(len(g1), len(g2))

            jacc = inter / union if union > 0 else 0.0
            ov_coef = inter / min_len if min_len > 0 else 0.0

            jaccard_mat[i, j] = jacc
            overlap_mat[i, j] = ov_coef

            if (s1 != s2 and i < j) or (s1 == s2 and i < j and len(strat_modules) == 1):
                shared_genes = sorted(list(g1.intersection(g2)))
                sample_str = ", ".join(shared_genes[:5])
                if len(shared_genes) > 5:
                    sample_str += f" (+{len(shared_genes)-5} more)"

                if jacc >= 0.3 or ov_coef >= 0.5:
                    status = "Conserved Program (High Overlap)"
                elif jacc >= 0.15 or ov_coef >= 0.25:
                    status = "Related Program (Partial Overlap)"
                else:
                    status = "Condition-Specific / Distinct"

                summary_rows.append({
                    "stratification_1": s1,
                    "module_1": m1,
                    "stratification_2": s2,
                    "module_2": m2,
                    "jaccard_similarity": round(jacc, 4),
                    "overlap_coefficient": round(ov_coef, 4),
                    "n_shared_genes": inter,
                    "n_module_1_genes": len(g1),
                    "n_module_2_genes": len(g2),
                    "shared_genes_sample": sample_str,
                    "alignment_status": status,
                })

    jaccard_df = pd.DataFrame(jaccard_mat, index=labels, columns=labels)
    overlap_coef_df = pd.DataFrame(overlap_mat, index=labels, columns=labels)

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows).sort_values(
            ["jaccard_similarity", "overlap_coefficient"], ascending=[False, False]
        )
    else:
        summary_df = pd.DataFrame(columns=[
            "stratification_1", "module_1", "stratification_2", "module_2",
            "jaccard_similarity", "overlap_coefficient", "n_shared_genes",
            "n_module_1_genes", "n_module_2_genes", "shared_genes_sample", "alignment_status"
        ])

    return jaccard_df, overlap_coef_df, summary_df


def compute_module_tf_integration(
    links_df: Optional[pd.DataFrame] = None,
    hotspot_obj: Any = None,
    score_df: Optional[pd.DataFrame] = None,
    adata: Optional[AnnData] = None,
    top_n_tfs: int = 15,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """Cross-reference module genes with TF-TG links and centrality scores per cluster."""
    module_genes = _extract_module_genes(hotspot_obj, adata, output_dir=output_dir, stratification_results=stratification_results)
    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)
    score_df = _normalize_score_df(score_df, output_dir=output_dir, stratification_results=stratification_results)

    cols = ["cluster", "module", "tf", "n_targets_in_module", "total_module_genes", "coverage_pct", "tf_centrality"]
    if not module_genes or links_df.empty or "cluster" not in links_df.columns:
        return pd.DataFrame(columns=cols)

    rows = []
    for cluster in links_df["cluster"].unique():
        c_links = links_df[links_df["cluster"] == cluster]
        top_tfs = c_links["source"].value_counts().head(top_n_tfs).index.tolist()
        for mod_name, genes in module_genes.items():
            for tf in top_tfs:
                tf_targets = set(c_links[c_links["source"] == tf]["target"].unique())
                overlap = genes.intersection(tf_targets)
                if len(overlap) > 0:
                    cent = np.nan
                    if not score_df.empty and "gene" in score_df.columns and "cluster" in score_df.columns:
                        s_df = score_df[(score_df["gene"] == tf) & (score_df["cluster"] == cluster)]
                        if not s_df.empty and "degree_centrality_all" in s_df.columns:
                            cent = s_df["degree_centrality_all"].values[0]
                    rows.append({
                        "cluster": cluster,
                        "module": mod_name,
                        "tf": tf,
                        "n_targets_in_module": len(overlap),
                        "total_module_genes": len(genes),
                        "coverage_pct": len(overlap) / len(genes),
                        "tf_centrality": cent if not np.isnan(cent) else 0.0,
                    })

    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows).sort_values("n_targets_in_module", ascending=False)


def compute_gene_selection_provenance(
    hotspot_obj: Any = None,
    links_df: Optional[pd.DataFrame] = None,
    enrichment_df: Optional[pd.DataFrame] = None,
    adata: Optional[AnnData] = None,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """Trace each gene through the pipeline: Hotspot significance -> Module assignment -> GRN role -> Pathway."""
    cols = ["gene", "hotspot_fdr", "hotspot_module", "stage", "is_tf", "is_target", "n_clusters_active", "regulated_by_tfs", "pathway"]

    results_df, modules_series = _extract_autocorr_results(hotspot_obj, output_dir=output_dir, stratification_results=stratification_results)
    if results_df is None:
        return pd.DataFrame(columns=cols)

    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)

    fdr_col = "FDR" if "FDR" in results_df.columns else ("fdr" if "fdr" in results_df.columns else None)
    if fdr_col:
        sig_genes = results_df[results_df[fdr_col] < 0.05].index.tolist()
    else:
        sig_genes = results_df.index.tolist()

    all_tfs = set()
    all_targets = set()
    gene_clusters: Dict[str, set] = {}
    target_tfs: Dict[str, Dict[str, int]] = {}
    if not links_df.empty:
        all_tfs = set(links_df["source"].unique())
        all_targets = set(links_df["target"].unique())
        if "cluster" in links_df.columns:
            src_clusters = links_df.groupby("source")["cluster"].apply(lambda x: set(x.unique())).to_dict()
            tgt_clusters = links_df.groupby("target")["cluster"].apply(lambda x: set(x.unique())).to_dict()
            for g, cl_set in src_clusters.items():
                gene_clusters.setdefault(g, set()).update(cl_set)
            for g, cl_set in tgt_clusters.items():
                gene_clusters.setdefault(g, set()).update(cl_set)
            tf_target_counts = links_df.groupby(["target", "source"]).size().reset_index(name="count")
            for _, r in tf_target_counts.iterrows():
                target_tfs.setdefault(r["target"], {})[r["source"]] = int(r["count"])

    pathway_map = {}
    if enrichment_df is not None and not enrichment_df.empty and "module" in enrichment_df.columns and "term" in enrichment_df.columns:
        for m in enrichment_df["module"].unique():
            m_df = enrichment_df[enrichment_df["module"] == m]
            if not m_df.empty:
                pathway_map[m] = m_df.iloc[0]["term"]

    rows = []
    for gene in sig_genes:
        fdr = results_df.loc[gene, fdr_col] if fdr_col and gene in results_df.index else 0.01
        if isinstance(fdr, pd.Series):
            fdr = fdr.iloc[0]
        mod = "Unassigned"
        if modules_series is not None and gene in modules_series.index:
            m_val = modules_series.loc[gene]
            if isinstance(m_val, pd.Series):
                m_val = m_val.iloc[0]
            if m_val != -1 and m_val != "-1" and str(m_val) != "-1":
                mod = f"Module {m_val}" if not str(m_val).startswith("Module") else str(m_val)

        is_tf = gene in all_tfs
        is_target = gene in all_targets
        n_clust = len(gene_clusters.get(gene, set()))

        reg_tfs = ""
        if gene in target_tfs:
            sorted_tfs = sorted(target_tfs[gene].items(), key=lambda x: x[1], reverse=True)[:5]
            reg_tfs = ", ".join([t[0] for t in sorted_tfs])

        if is_tf and is_target:
            stage = "TF & Target"
        elif is_tf:
            stage = "TF"
        elif is_target:
            stage = "Target"
        elif mod != "Unassigned":
            stage = "Module Member"
        else:
            stage = "Significant Only"

        pathway = pathway_map.get(mod, "")

        rows.append({
            "gene": gene,
            "hotspot_fdr": fdr,
            "hotspot_module": mod,
            "stage": stage,
            "is_tf": is_tf,
            "is_target": is_target,
            "n_clusters_active": n_clust,
            "regulated_by_tfs": reg_tfs,
            "pathway": pathway,
        })

    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows).sort_values("hotspot_fdr", ascending=True)


def compute_cross_cluster_regulatory_summary(
    links_df: Optional[pd.DataFrame] = None,
    score_df: Optional[pd.DataFrame] = None,
    hotspot_obj: Any = None,
    activity_df: Optional[pd.DataFrame] = None,
    adata: Optional[AnnData] = None,
    top_n: int = 5,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """Per-cluster summary of regulatory activity for side-by-side comparison."""
    cols = ["cluster", "n_active_modules", "top_modules", "n_active_tfs", "top_tfs", "n_regulatory_edges", "n_unique_targets", "top_edges"]
    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)
    score_df = _normalize_score_df(score_df, output_dir=output_dir, stratification_results=stratification_results)
    module_genes = _extract_module_genes(hotspot_obj, adata, output_dir=output_dir, stratification_results=stratification_results)

    if links_df.empty or "cluster" not in links_df.columns:
        return pd.DataFrame(columns=cols)

    clusters = links_df["cluster"].unique()
    rows = []
    for c in clusters:
        c_links = links_df[links_df["cluster"] == c]
        n_edges = len(c_links)
        unique_targets = c_links["target"].nunique()
        unique_tfs = c_links["source"].nunique()

        if not score_df.empty and "gene" in score_df.columns and "cluster" in score_df.columns and "degree_centrality_all" in score_df.columns:
            s_df = score_df[score_df["cluster"] == c].sort_values("degree_centrality_all", ascending=False)
            top_tfs_list = s_df["gene"].head(top_n).tolist()
            if not top_tfs_list:
                top_tfs_list = c_links["source"].value_counts().head(top_n).index.tolist()
        else:
            top_tfs_list = c_links["source"].value_counts().head(top_n).index.tolist()

        active_mods = []
        if module_genes:
            c_targets = set(c_links["target"].unique())
            for mod_name, genes in module_genes.items():
                if len(genes) > 0:
                    overlap = len(genes.intersection(c_targets))
                    if overlap / len(genes) > 0.2:
                        active_mods.append(mod_name)

        if "coef_abs" in c_links.columns:
            top_e_df = c_links.sort_values("coef_abs", ascending=False).head(top_n)
            top_edges_list = [f"{r['source']}->{r['target']}" for _, r in top_e_df.iterrows()]
        else:
            top_edges_list = []

        rows.append({
            "cluster": c,
            "n_active_modules": len(active_mods),
            "top_modules": ", ".join(active_mods) if active_mods else "None",
            "n_active_tfs": unique_tfs,
            "top_tfs": ", ".join(top_tfs_list),
            "n_regulatory_edges": n_edges,
            "n_unique_targets": unique_targets,
            "top_edges": ", ".join(top_edges_list),
        })

    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows)


def compute_tf_module_concordance(
    links_df: Optional[pd.DataFrame] = None,
    score_df: Optional[pd.DataFrame] = None,
    activity_df: Optional[pd.DataFrame] = None,
    hotspot_obj: Any = None,
    adata: Optional[AnnData] = None,
    top_n_tfs_per_module: int = 3,
    output_dir: Optional[str] = None,
    stratification_results: Optional[Union[Dict[str, Any], List[Any]]] = None,
) -> pd.DataFrame:
    """
    Directly correlate and integrate CellOracle TF-TG regulation with Hotspot Module activity
    per cluster to identify the driving regulatory circuits.
    """
    cols = [
        "cluster", "module", "module_activity", "top_driver_tf",
        "n_targets_in_module", "module_coverage_pct", "tf_centrality", "concordance_score"
    ]
    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)
    score_df = _normalize_score_df(score_df, output_dir=output_dir, stratification_results=stratification_results)
    module_genes = _extract_module_genes(hotspot_obj, adata, output_dir=output_dir, stratification_results=stratification_results)

    if links_df.empty or not module_genes:
        return pd.DataFrame(columns=cols)

    clusters = links_df["cluster"].unique()
    rows = []

    for c in clusters:
        c_links = links_df[links_df["cluster"] == c]
        for mod_name, genes in module_genes.items():
            if not genes:
                continue

            # Check module activity in this cluster
            mod_act = 0.0
            if activity_df is not None and not activity_df.empty:
                # Find matching column
                col_match = [col for col in activity_df.columns if str(c) in str(col)]
                row_match = [idx for idx in activity_df.index if str(mod_name) in str(idx)]
                if col_match and row_match:
                    try:
                        mod_act = float(activity_df.loc[row_match[0], col_match[0]])
                    except Exception:
                        mod_act = 0.0

            # Find TFs targeting this module
            mod_links = c_links[c_links["target"].isin(genes)]
            if mod_links.empty:
                continue

            tf_target_counts = mod_links.groupby("source")["target"].nunique().sort_values(ascending=False)
            top_tf = tf_target_counts.index[0]
            n_targets = int(tf_target_counts.iloc[0])
            coverage_pct = n_targets / len(genes)

            tf_cent = 0.0
            if not score_df.empty and "gene" in score_df.columns and "cluster" in score_df.columns:
                s_match = score_df[(score_df["gene"] == top_tf) & (score_df["cluster"] == c)]
                if not s_match.empty and "degree_centrality_all" in s_match.columns:
                    tf_cent = float(s_match["degree_centrality_all"].values[0])

            concordance = (abs(mod_act) + 0.1) * coverage_pct * (1.0 + tf_cent)

            rows.append({
                "cluster": c,
                "module": mod_name,
                "module_activity": mod_act,
                "top_driver_tf": top_tf,
                "n_targets_in_module": n_targets,
                "module_coverage_pct": coverage_pct,
                "tf_centrality": tf_cent,
                "concordance_score": concordance,
            })

    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows).sort_values("concordance_score", ascending=False)


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

    # Normalize inputs
    links_df = _normalize_links_df(links_df, output_dir=output_dir, stratification_results=stratification_results)
    score_df = _normalize_score_df(score_df, output_dir=output_dir, stratification_results=stratification_results)

    results: Dict[str, Any] = {}

    # 1. Module Activity Matrix
    try:
        activity_df = compute_module_activity_matrix(
            adata=adata,
            cluster_key=cluster_key,
            hotspot_obj=hotspot_obj,
            stratification_results=stratification_results,
            output_dir=output_dir,
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
        tf_mod_matrix, tf_mod_summary = compute_tf_to_module_mapping(
            links_df=links_df,
            modules_dict_or_df=hotspot_obj,
            output_dir=output_dir,
            stratification_results=stratification_results,
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
        diff_targets_df = compute_differential_tf_targets(
            links_df=links_df,
            output_dir=output_dir,
            stratification_results=stratification_results,
        )
        results["differential_tf_targets"] = diff_targets_df
        if save_tables and not diff_targets_df.empty:
            diff_targets_df.to_csv(
                os.path.join(comp_dir, "differential_tf_targets.csv"), index=False
            )
            print("  ✓ Computed differential TF target rewiring")
    except Exception as e:
        log_error("ComparativeAnalysis.DiffTargets", e)
        results["differential_tf_targets"] = pd.DataFrame()

    # 6. Module Gene Overlap
    try:
        coverage_df, jaccard_df = compute_module_gene_overlap_matrix(
            hotspot_obj=hotspot_obj,
            links_df=links_df,
            adata=adata,
            cluster_key=cluster_key,
            output_dir=output_dir,
            stratification_results=stratification_results,
        )
        results["module_coverage"] = coverage_df
        results["module_jaccard"] = jaccard_df
        if save_tables and not coverage_df.empty:
            coverage_df.to_csv(os.path.join(comp_dir, "module_gene_coverage.csv"))
            jaccard_df.to_csv(os.path.join(comp_dir, "module_gene_jaccard.csv"))
            print("  ✓ Computed module gene overlap matrix")
    except Exception as e:
        log_error("ComparativeAnalysis.ModuleOverlap", e)
        results["module_coverage"] = pd.DataFrame()
        results["module_jaccard"] = pd.DataFrame()

    # 7. Module-TF Integration
    try:
        integration_df = compute_module_tf_integration(
            links_df=links_df,
            hotspot_obj=hotspot_obj,
            score_df=score_df,
            adata=adata,
            output_dir=output_dir,
            stratification_results=stratification_results,
        )
        results["module_tf_integration"] = integration_df
        if save_tables and not integration_df.empty:
            integration_df.to_csv(os.path.join(comp_dir, "module_tf_integration.csv"), index=False)
            print("  ✓ Computed module-TF regulatory integration")
    except Exception as e:
        log_error("ComparativeAnalysis.ModuleTFIntegration", e)
        results["module_tf_integration"] = pd.DataFrame()

    # 8. Gene Selection Provenance
    try:
        provenance_df = compute_gene_selection_provenance(
            hotspot_obj=hotspot_obj,
            links_df=links_df,
            enrichment_df=results.get("module_enrichment"),
            adata=adata,
            output_dir=output_dir,
            stratification_results=stratification_results,
        )
        results["gene_provenance"] = provenance_df
        if save_tables and not provenance_df.empty:
            provenance_df.to_csv(os.path.join(comp_dir, "gene_selection_provenance.csv"), index=False)
            print("  ✓ Computed gene selection provenance")
    except Exception as e:
        log_error("ComparativeAnalysis.GeneProvenance", e)
        results["gene_provenance"] = pd.DataFrame()

    # 9. Cross-Cluster Regulatory Summary
    try:
        reg_summary_df = compute_cross_cluster_regulatory_summary(
            links_df=links_df,
            score_df=score_df,
            hotspot_obj=hotspot_obj,
            activity_df=results.get("module_activity"),
            adata=adata,
            output_dir=output_dir,
            stratification_results=stratification_results,
        )
        results["regulatory_summary"] = reg_summary_df
        if save_tables and not reg_summary_df.empty:
            reg_summary_df.to_csv(os.path.join(comp_dir, "cross_cluster_regulatory_summary.csv"), index=False)
            print("  ✓ Computed cross-cluster regulatory summary")
    except Exception as e:
        log_error("ComparativeAnalysis.RegulatorySummary", e)
        results["regulatory_summary"] = pd.DataFrame()

    # 10. TF-Module Concordance Integration
    try:
        concordance_df = compute_tf_module_concordance(
            links_df=links_df,
            score_df=score_df,
            activity_df=results.get("module_activity"),
            hotspot_obj=hotspot_obj,
            adata=adata,
            output_dir=output_dir,
            stratification_results=stratification_results,
        )
        results["tf_module_concordance"] = concordance_df
        if save_tables and not concordance_df.empty:
            concordance_df.to_csv(os.path.join(comp_dir, "tf_module_concordance.csv"), index=False)
            print("  ✓ Computed TF-Module regulatory concordance")
    except Exception as e:
        log_error("ComparativeAnalysis.Concordance", e)
        results["tf_module_concordance"] = pd.DataFrame()

    # 11. Cross-Stratification Module Gene Alignment
    try:
        mod_jaccard_df, mod_ov_coef_df, mod_align_summary = compute_cross_stratification_module_overlap(
            stratification_results=stratification_results,
            hotspot_obj=hotspot_obj,
            adata=adata,
            output_dir=output_dir,
        )
        results["cross_strat_module_jaccard"] = mod_jaccard_df
        results["cross_strat_module_overlap"] = mod_ov_coef_df
        results["cross_strat_module_alignment"] = mod_align_summary
        if save_tables and not mod_jaccard_df.empty:
            mod_jaccard_df.to_csv(os.path.join(comp_dir, "cross_stratification_module_jaccard.csv"))
            mod_ov_coef_df.to_csv(os.path.join(comp_dir, "cross_stratification_module_overlap_coefficient.csv"))
            mod_align_summary.to_csv(os.path.join(comp_dir, "cross_stratification_module_alignment.csv"), index=False)
            print("  ✓ Computed cross-stratification module gene alignment")
    except Exception as e:
        log_error("ComparativeAnalysis.CrossStratModuleOverlap", e)
        results["cross_strat_module_jaccard"] = pd.DataFrame()
        results["cross_strat_module_overlap"] = pd.DataFrame()
        results["cross_strat_module_alignment"] = pd.DataFrame()

    print("  Comparative analysis complete.")
    return results
