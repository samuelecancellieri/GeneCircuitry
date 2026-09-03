"""
Local Gene Set Over-Representation Analysis (ORA) Engine
========================================================

Provides offline-first, fast, and resilient functional pathway enrichment
analysis using local gene set libraries (e.g., bundled MSigDB Hallmark,
custom GMT files, Python dictionaries, or locally cached Enrichr libraries).

Features:
- Pure local computation (via hypergeometric / Fisher's exact test) in milliseconds.
- Zero network dependencies during analysis runs.
- Bundled default gene sets (MSigDB Hallmark 2020 for Human and Mouse).
- Local disk caching for additional libraries.
- Support for custom GMT files, dicts, and background gene universes.
"""

from __future__ import annotations

import glob
import os
from pathlib import Path
import time
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import scipy.stats as stats

from . import config
from .logging_utils import log_error, log_info, log_warning

# Standard result column schema
RESULT_COLUMNS = [
    "Gene_set",
    "Term",
    "Overlap",
    "P-value",
    "Adjusted P-value",
    "Odds Ratio",
    "Combined Score",
    "Genes",
]


class EnrichmentResult:
    """
    Container class for gene set enrichment analysis results.

    Attributes
    ----------
    results : pd.DataFrame
        Table of enriched terms with statistics (P-value, Adjusted P-value, Combined Score, etc.).
    res2d : pd.DataFrame
        Alias for ``results`` for compatibility with gseapy Enrichr objects.
    gene_list : list of str
        The query gene list analyzed.
    gene_sets : list or dict or str
        The gene set(s) used for enrichment.
    organism : str
        The target organism/species.
    """

    def __init__(
        self,
        results: Optional[pd.DataFrame] = None,
        gene_list: Optional[List[str]] = None,
        gene_sets: Optional[Union[str, List[str], Dict[str, List[str]]]] = None,
        organism: str = "human",
    ):
        if results is None or results.empty:
            self.results = pd.DataFrame(columns=RESULT_COLUMNS)
        else:
            self.results = results.copy()
            # Ensure standard column names exist
            for col in RESULT_COLUMNS:
                if col not in self.results.columns:
                    # Check for underscore variants
                    alt_col = col.replace(" ", "_")
                    if alt_col in self.results.columns:
                        self.results[col] = self.results[alt_col]

        self.res2d = self.results
        self.gene_list = list(gene_list) if gene_list is not None else []
        self.gene_sets = gene_sets
        self.organism = organism

    @property
    def empty(self) -> bool:
        """Return True if results are empty."""
        return self.results.empty

    def to_dataframe(self) -> pd.DataFrame:
        """Return results as a pandas DataFrame."""
        return self.results.copy()

    def head(self, n: int = 5) -> pd.DataFrame:
        """Return top n enriched terms."""
        return self.results.head(n)

    def __repr__(self) -> str:
        num_terms = len(self.results)
        return f"<EnrichmentResult: {num_terms} enriched terms across {self.gene_sets}>"


def get_bundled_gene_sets_dir() -> str:
    """
    Get the directory containing bundled gene set GMT files.

    Returns
    -------
    str
        Path to genecircuitry/data/gene_sets/
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    bundled_dir = os.path.join(pkg_dir, "data", "gene_sets")
    return bundled_dir


def get_cache_gene_sets_dir(cache_dir: Optional[str] = None) -> str:
    """
    Get or create the directory used for caching downloaded gene set GMT files.

    Parameters
    ----------
    cache_dir : str, optional
        Base cache directory. If None, uses ``config.CACHE_DIR`` (or ~/.cache/genecircuitry).

    Returns
    -------
    str
        Path to gene sets cache directory.
    """
    base = cache_dir or getattr(config, "CACHE_DIR", ".cache")
    gene_sets_cache = os.path.join(base, "gene_sets")
    os.makedirs(gene_sets_cache, exist_ok=True)
    return gene_sets_cache


def read_gmt(path: Union[str, Path]) -> Dict[str, List[str]]:
    """
    Read a gene set database file in GMT format.

    Parameters
    ----------
    path : str or Path
        Path to a .gmt file.

    Returns
    -------
    dict
        Mapping of term name to list of gene symbols.

    Raises
    ------
    FileNotFoundError
        If the GMT file does not exist.
    """
    path_str = str(path)
    if not os.path.exists(path_str):
        raise FileNotFoundError(f"GMT file not found: {path_str}")

    gene_sets = {}
    with open(path_str, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            # parts[1] is description / URL
            genes = [g.strip() for g in parts[2:] if g.strip()]
            if term and genes:
                gene_sets[term] = genes

    return gene_sets


def write_gmt(
    gene_sets: Dict[str, List[str]],
    path: Union[str, Path],
    description: str = "",
) -> str:
    """
    Write a gene set dictionary to a GMT file.

    Parameters
    ----------
    gene_sets : dict
        Mapping from term name to list of gene symbols.
    path : str or Path
        Destination file path.
    description : str, optional
        Description / URL string placed in column 2 of the GMT file.

    Returns
    -------
    str
        The file path written to.
    """
    path_str = str(path)
    os.makedirs(os.path.dirname(os.path.abspath(path_str)), exist_ok=True)

    with open(path_str, "w", encoding="utf-8") as f:
        for term, genes in gene_sets.items():
            desc = description or f"https://genecircuitry.org/geneset/{term}"
            genes_clean = [str(g).strip() for g in genes if str(g).strip()]
            if genes_clean:
                f.write(f"{term}\t{desc}\t" + "\t".join(genes_clean) + "\n")

    return path_str


def list_available_gene_sets(
    include_cached: bool = True,
    cache_dir: Optional[str] = None,
) -> Dict[str, str]:
    """
    List all locally available gene set libraries (bundled and cached).

    Parameters
    ----------
    include_cached : bool, default=True
        Whether to include libraries found in the local cache directory.
    cache_dir : str, optional
        Custom cache directory to inspect.

    Returns
    -------
    dict
        Mapping of library name to absolute file path.
    """
    available = {}

    # 1. Inspect bundled directory
    bundled_dir = get_bundled_gene_sets_dir()
    if os.path.isdir(bundled_dir):
        for gmt_path in glob.glob(os.path.join(bundled_dir, "*.gmt")):
            name = os.path.splitext(os.path.basename(gmt_path))[0]
            available[name] = os.path.abspath(gmt_path)

    # 2. Inspect cache directory
    if include_cached:
        c_dir = get_cache_gene_sets_dir(cache_dir=cache_dir)
        if os.path.isdir(c_dir):
            for gmt_path in glob.glob(os.path.join(c_dir, "*.gmt")):
                name = os.path.splitext(os.path.basename(gmt_path))[0]
                if name not in available:
                    available[name] = os.path.abspath(gmt_path)

    return available


def cache_gene_sets(
    library_name: str,
    organism: str = "human",
    cache_dir: Optional[str] = None,
    force: bool = False,
) -> str:
    """
    Pre-download an Enrichr gene set library and save it locally as a GMT file.

    Parameters
    ----------
    library_name : str
        Name of the library (e.g. 'Reactome_Pathways_2024', 'KEGG_2021_Human').
    organism : str, default='human'
        Target organism (e.g. 'Human', 'Mouse').
    cache_dir : str, optional
        Custom cache directory.
    force : bool, default=False
        If True, re-downloads even if already cached.

    Returns
    -------
    str
        Path to the cached GMT file.

    Raises
    ------
    RuntimeError
        If download fails or library is not available.
    """
    c_dir = get_cache_gene_sets_dir(cache_dir=cache_dir)
    target_path = os.path.join(c_dir, f"{library_name}.gmt")

    if os.path.exists(target_path) and not force:
        return target_path

    try:
        from gseapy.parser import download_library

        org_name = organism.capitalize()
        lib_dict = download_library(library_name, organism=org_name)
        if not lib_dict:
            raise ValueError(
                f"No gene sets found in library '{library_name}' for organism '{org_name}'."
            )

        write_gmt(
            lib_dict, target_path, description=f"Enrichr {library_name} ({org_name})"
        )
        log_info(
            "EnrichmentAnalysis.Cache",
            f"Cached gene set library '{library_name}' to {target_path}",
        )
        return target_path
    except Exception as e:
        log_error("EnrichmentAnalysis.Cache", e)
        raise RuntimeError(
            f"Failed to download and cache gene set library '{library_name}' ({type(e).__name__}): {e}"
        ) from e


def load_gene_set(
    name_or_path_or_dict: Union[str, Path, Dict[str, List[str]]],
    organism: str = "human",
    cache_dir: Optional[str] = None,
) -> Tuple[str, Dict[str, List[str]]]:
    """
    Resolve and load a gene set specification into a dictionary.

    Supports:
    1. Python dict: ``{"Term": ["GENE1", "GENE2", ...]}`` -> used directly.
    2. Local file path: ``"path/to/custom.gmt"`` -> parsed locally.
    3. Bundled library name: ``"MSigDB_Hallmark_2020"`` -> loaded from bundled data.
    4. Cached library name: loaded from disk cache without network.
    5. Online library name: downloaded once, cached to disk, then loaded.

    Parameters
    ----------
    name_or_path_or_dict : str or Path or dict
        Gene set name, file path, or dictionary.
    organism : str, default='human'
        Target organism (e.g. 'human', 'mouse').
    cache_dir : str, optional
        Custom cache directory.

    Returns
    -------
    tuple of (str, dict)
        (library_name, gene_set_dictionary)
    """
    if isinstance(name_or_path_or_dict, dict):
        return ("custom", name_or_path_or_dict)

    spec = str(name_or_path_or_dict).strip()

    # Direct file path
    if os.path.isfile(spec):
        lib_name = os.path.splitext(os.path.basename(spec))[0]
        return (lib_name, read_gmt(spec))

    # Check available bundled and cached libraries
    available = list_available_gene_sets(include_cached=True, cache_dir=cache_dir)

    # Candidate name variations
    org_clean = organism.lower()
    candidates = [
        spec,
        f"{spec}_{organism.capitalize()}",
        f"{spec}_{organism.upper()}",
        f"{spec}_{org_clean}",
    ]

    for cand in candidates:
        if cand in available:
            return (spec, read_gmt(available[cand]))

    # Check bundled dir with case-insensitive matching
    bundled_dir = get_bundled_gene_sets_dir()
    if os.path.isdir(bundled_dir):
        for fname in os.listdir(bundled_dir):
            base, ext = os.path.splitext(fname)
            if ext.lower() == ".gmt":
                if base.lower() in [c.lower() for c in candidates]:
                    full_p = os.path.join(bundled_dir, fname)
                    return (spec, read_gmt(full_p))

    # Check cache dir with case-insensitive matching
    c_dir = get_cache_gene_sets_dir(cache_dir=cache_dir)
    if os.path.isdir(c_dir):
        for fname in os.listdir(c_dir):
            base, ext = os.path.splitext(fname)
            if ext.lower() == ".gmt":
                if base.lower() in [c.lower() for c in candidates]:
                    full_p = os.path.join(c_dir, fname)
                    return (spec, read_gmt(full_p))

    # Not found locally: attempt download once and cache to disk
    try:
        cached_file = cache_gene_sets(spec, organism=organism, cache_dir=cache_dir)
        return (spec, read_gmt(cached_file))
    except Exception as e:
        raise ValueError(
            f"Gene set library '{spec}' could not be found locally or downloaded.\n"
            f"  Locally available libraries: {list(available.keys())}\n"
            f"  To use offline, provide a local .gmt file path or dictionary, or pre-cache with:\n"
            f"    genecircuitry.enrichment_analysis.cache_gene_sets('{spec}')"
        ) from e


def _calculate_ora_local(
    gene_list: List[str],
    gene_set_dict: Dict[str, List[str]],
    background: Union[int, List[str], Set[str], None] = None,
    pval_cutoff: float = 0.05,
    lib_name: str = "custom",
) -> pd.DataFrame:
    """
    Perform local Over-Representation Analysis using hypergeometric / Fisher's exact test.

    Parameters
    ----------
    gene_list : list of str
        Query gene symbols.
    gene_set_dict : dict
        Mapping of term name to list of gene symbols.
    background : int, list, set, optional
        Background universe for hypergeometric test.
    pval_cutoff : float, default=0.05
        Cutoff for Adjusted P-value.
    lib_name : str, default='custom'
        Name of the gene set library.

    Returns
    -------
    pd.DataFrame
        Table of enrichment statistics.
    """
    # Clean query genes (uppercase, stripped, unique)
    genes_query = {str(g).strip().upper() for g in gene_list if g and str(g).strip()}
    if not genes_query or not gene_set_dict:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    # Standardize gene set dictionary to uppercase
    gs_upper = {
        str(term).strip(): [
            str(g).strip().upper() for g in glist if g and str(g).strip()
        ]
        for term, glist in gene_set_dict.items()
    }

    # Determine background universe size N
    if background is None:
        bg_set = set.union(*[set(v) for v in gs_upper.values()]) | genes_query
        N = len(bg_set)
    elif isinstance(background, int):
        N = max(background, len(genes_query))
    elif isinstance(background, (list, set, tuple, pd.Index)):
        bg_set = {
            str(g).strip().upper() for g in background if g and str(g).strip()
        } | genes_query
        N = len(bg_set)
    else:
        N = 20000

    n = len(genes_query)  # number of query genes (draws)
    rows = []

    for term, term_genes in gs_upper.items():
        term_set = set(term_genes)
        K = len(term_set)  # total genes in category
        if K == 0:
            continue

        overlap = genes_query & term_set
        k = len(overlap)  # overlapping genes (successes in sample)
        if k == 0:
            continue

        # Hypergeometric survival function: P(X >= k)
        # scipy.stats.hypergeom.sf(k - 1, M, n, N) where M=total universe (N), n=total category genes (K), N=draws (n)
        pval = float(stats.hypergeom.sf(k - 1, N, K, n))

        # Odds Ratio from 2x2 contingency table:
        # [[k (query & in_set), n - k (query & not_in_set)],
        #  [K - k (not_query & in_set), N - K - (n - k) (not_query & not_in_set)]]
        a = k
        b = n - k
        c = K - k
        d = max(0, N - K - (n - k))
        odds_ratio = float((a * d) / (b * c)) if (b * c) > 0 else float("inf")

        # Combined Score = -ln(P) * Z-score / Odds Ratio
        log_p = -float(np.log(max(pval, 1e-300)))
        combined_score = float(
            log_p * (odds_ratio if np.isfinite(odds_ratio) else 100.0)
        )

        rows.append(
            {
                "Gene_set": lib_name,
                "Term": term,
                "Overlap": f"{k}/{K}",
                "P-value": pval,
                "Odds Ratio": odds_ratio,
                "Combined Score": combined_score,
                "Genes": ";".join(sorted(overlap)),
            }
        )

    if not rows:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    df = pd.DataFrame(rows)

    # Benjamini-Hochberg FDR correction
    df = df.sort_values("P-value").reset_index(drop=True)
    pvals = df["P-value"].values
    m = len(gs_upper)  # Total hypotheses tested
    qvals = np.minimum.accumulate((pvals * m / np.arange(1, len(pvals) + 1))[::-1])[
        ::-1
    ]
    df["Adjusted P-value"] = np.clip(qvals, 0.0, 1.0)

    # Reorder columns to standard schema
    df = df[RESULT_COLUMNS]

    # Filter by pval_cutoff
    df = df[df["Adjusted P-value"] <= pval_cutoff].reset_index(drop=True)

    # Sort by Combined Score descending or Adjusted P-value ascending
    df = df.sort_values(by=["Adjusted P-value", "P-value"]).reset_index(drop=True)

    return df


def _run_enrichr_online(
    gene_list: List[str],
    gene_sets: Union[str, Path, Dict[str, List[str]], List[Any]],
    species: str = "human",
    pval_cutoff: float = 0.05,
    background: Optional[Union[int, List[str], Set[str]]] = None,
    top_n_terms: Optional[int] = None,
    max_retries: int = 3,
    retry_delay: float = 2.0,
) -> pd.DataFrame:
    """
    Execute Over-Representation Analysis using the online Enrichr web API via gseapy.

    Parameters
    ----------
    gene_list : list of str
        Query gene symbols.
    gene_sets : str, Path, dict, or list
        Gene set specification (Enrichr library name(s), GMT file path, or dictionary).
    species : str, default='human'
        Target organism for Enrichr.
    pval_cutoff : float, default=0.05
        Significance cutoff for Adjusted P-value (FDR).
    background : int, list, set, optional
        Background universe (used by gseapy if supported).
    top_n_terms : int, optional
        Maximum number of terms to retain per gene set library.
    max_retries : int, default=3
        Maximum retry attempts for transient API failures or rate limits (HTTP 429).
    retry_delay : float, default=2.0
        Initial backoff delay in seconds.

    Returns
    -------
    pd.DataFrame
        Table of enrichment statistics matching RESULT_COLUMNS.
    """
    import gseapy

    query_genes = [str(g).strip() for g in gene_list if g and str(g).strip()]
    if not query_genes:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    # Format background argument for gseapy.enrichr
    bg_arg: Optional[Union[int, List[str], str]] = None
    if background is not None:
        if isinstance(background, set):
            bg_arg = list(background)
        elif isinstance(background, (int, list, str)):
            bg_arg = background

    # Normalize gene_sets for gseapy
    if isinstance(gene_sets, Path):
        gs_arg: Any = str(gene_sets)
    elif isinstance(gene_sets, list):
        gs_arg = [str(x) if isinstance(x, Path) else x for x in gene_sets]
    else:
        gs_arg = gene_sets

    sleep_time = retry_delay
    enr = None

    for attempt in range(max_retries):
        try:
            enr = gseapy.enrichr(
                gene_list=query_genes,
                gene_sets=gs_arg,
                organism=species,
                background=bg_arg,
                no_plot=True,
                verbose=False,
            )
            break
        except Exception as e:
            err_str = str(e).lower()
            is_retryable = (
                "429" in err_str
                or "timeout" in err_str
                or "connection" in err_str
                or "503" in err_str
                or "502" in err_str
                or "too many requests" in err_str
            )
            if is_retryable and attempt < max_retries - 1:
                log_warning(
                    "EnrichmentAnalysis.OnlineEnrichr",
                    f"Enrichr API request failed ({type(e).__name__}: {e}); "
                    f"retrying in {sleep_time:.1f}s (attempt {attempt + 1}/{max_retries})...",
                )
                time.sleep(sleep_time)
                sleep_time *= 2
            else:
                log_error("EnrichmentAnalysis.OnlineEnrichr", e)
                raise RuntimeError(
                    f"Online Enrichr API analysis failed ({type(e).__name__}): {e}.\n"
                    f"To run offline without network access, use online=False or "
                    f"set config.ENRICHMENT_ONLINE = False."
                ) from e

    if enr is None or getattr(enr, "results", None) is None or enr.results.empty:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    df = enr.results.copy()

    # Standardize column names (Enrichr returns 'Adjusted P-value', 'P-value', etc.)
    for col in RESULT_COLUMNS:
        if col not in df.columns:
            alt_col = col.replace(" ", "_")
            if alt_col in df.columns:
                df[col] = df[alt_col]

    # Filter by pval_cutoff
    if "Adjusted P-value" in df.columns:
        df = df[df["Adjusted P-value"] <= pval_cutoff].reset_index(drop=True)

    # Apply top_n_terms per gene set library if specified
    if top_n_terms is not None and top_n_terms > 0 and not df.empty:
        if "Gene_set" in df.columns:
            filtered_dfs = []
            for _, group_df in df.groupby("Gene_set", sort=False):
                filtered_dfs.append(group_df.head(top_n_terms))
            df = pd.concat(filtered_dfs, ignore_index=True)
        else:
            df = df.head(top_n_terms)

    # Sort by Adjusted P-value ascending, P-value ascending
    sort_cols = [c for c in ["Adjusted P-value", "P-value"] if c in df.columns]
    if sort_cols:
        df = df.sort_values(by=sort_cols).reset_index(drop=True)

    # Reorder columns to ensure RESULT_COLUMNS come first
    ordered_cols = [c for c in RESULT_COLUMNS if c in df.columns] + [
        c for c in df.columns if c not in RESULT_COLUMNS
    ]
    df = df[ordered_cols]

    return df


def run_enrichr_online(
    gene_list: Iterable[str],
    gene_sets: Optional[
        Union[
            str,
            Path,
            Dict[str, List[str]],
            List[Union[str, Dict[str, List[str]]]],
        ]
    ] = None,
    background: Optional[Union[int, List[str], Set[str]]] = None,
    pval_cutoff: float = 0.05,
    species: str = "human",
    top_n_terms: Optional[int] = None,
    max_retries: int = 3,
    retry_delay: float = 2.0,
) -> EnrichmentResult:
    """
    Perform online pathway enrichment analysis using the Enrichr web API via gseapy.

    Parameters
    ----------
    gene_list : iterable of str
        List of query gene symbols.
    gene_sets : str, Path, dict, or list, optional
        Gene set specification. Defaults to ``config.ENRICHMENT_GENE_SETS``.
    background : int, list, set, optional
        Background universe.
    pval_cutoff : float, default=0.05
        Significance cutoff for Adjusted P-value.
    species : str, default='human'
        Target organism.
    top_n_terms : int, optional
        Number of top terms to retain per gene set.
    max_retries : int, default=3
        Maximum retries for network or rate limit errors.
    retry_delay : float, default=2.0
        Initial backoff sleep time in seconds.

    Returns
    -------
    EnrichmentResult
        Container with results DataFrame.
    """
    if gene_sets is None:
        gene_sets = getattr(config, "ENRICHMENT_GENE_SETS", ["MSigDB_Hallmark_2020"])
    if species is None:
        species = getattr(config, "ENRICHMENT_SPECIES", "human")

    query_genes = [str(g).strip() for g in gene_list if g and str(g).strip()]
    df = _run_enrichr_online(
        gene_list=query_genes,
        gene_sets=gene_sets,
        species=species,
        pval_cutoff=pval_cutoff,
        background=background,
        top_n_terms=top_n_terms,
        max_retries=max_retries,
        retry_delay=retry_delay,
    )
    return EnrichmentResult(
        results=df,
        gene_list=query_genes,
        gene_sets=gene_sets,
        organism=species,
    )


def perform_ora_enrichment(
    gene_list: Iterable[str],
    gene_sets: Optional[
        Union[
            str,
            Path,
            Dict[str, List[str]],
            List[Union[str, Dict[str, List[str]]]],
        ]
    ] = None,
    background: Optional[Union[int, List[str], Set[str]]] = None,
    pval_cutoff: float = 0.05,
    species: str = "human",
    top_n_terms: Optional[int] = None,
    cache_dir: Optional[str] = None,
    online: Optional[bool] = None,
) -> EnrichmentResult:
    """
    Perform Over-Representation Analysis (ORA) on a gene list.

    By default, executes statistical hypergeometric enrichment entirely locally
    on CPU with zero remote network requests. Supports bundled gene sets (e.g.
    MSigDB Hallmark), custom GMT files, Python dictionaries, or disk-cached libraries.
    When ``online=True`` (or ``config.ENRICHMENT_ONLINE=True``), performs enrichment
    via the online Enrichr web API.

    Parameters
    ----------
    gene_list : iterable of str
        List of query gene symbols (e.g., Hotspot module genes or TF targets).
    gene_sets : str, Path, dict, or list, optional
        Gene set specification. Defaults to ``config.ENRICHMENT_GENE_SETS``.
        Can be:
        - A library name (e.g., ``"MSigDB_Hallmark_2020"``)
        - A path to a ``.gmt`` file
        - A dictionary: ``{"Term": ["GENE1", "GENE2"]}``
        - A list of any of the above.
    background : int, list, set, optional
        Background gene universe. Defaults to ``config.ENRICHMENT_BACKGROUND`` (or None).
        - If None: Uses the union of genes in the provided gene set library.
        - If int: Total number of genes tested (e.g. 20000).
        - If list/set: Specific background gene symbols (e.g. ``adata.var_names``).
    pval_cutoff : float, default=0.05
        Significance cutoff for Adjusted P-value (FDR).
    species : str, default='human'
        Target organism (e.g. 'human', 'mouse').
    top_n_terms : int, optional
        If specified, retains only the top N terms per gene set library.
    cache_dir : str, optional
        Custom cache directory for gene sets.
    online : bool, optional
        Whether to use the online Enrichr web API instead of local calculation.
        Defaults to ``config.ENRICHMENT_ONLINE`` (False).

    Returns
    -------
    EnrichmentResult
        Result container with ``.results`` DataFrame.

    Examples
    --------
    >>> from genecircuitry.enrichment_analysis import perform_ora_enrichment
    >>> res = perform_ora_enrichment(["TP53", "MDM2", "BAX", "CASP3"])
    >>> print(res.results.head())
    """
    if gene_sets is None:
        gene_sets = getattr(config, "ENRICHMENT_GENE_SETS", ["MSigDB_Hallmark_2020"])

    if background is None:
        background = getattr(config, "ENRICHMENT_BACKGROUND", None)

    if species is None:
        species = getattr(config, "ENRICHMENT_SPECIES", "human")

    if online is None:
        online = bool(getattr(config, "ENRICHMENT_ONLINE", False))

    query_genes = [str(g).strip() for g in gene_list if g and str(g).strip()]
    if not query_genes:
        return EnrichmentResult(
            results=pd.DataFrame(columns=RESULT_COLUMNS),
            gene_list=[],
            gene_sets=gene_sets,
            organism=species,
        )

    # Online Enrichr web API mode
    if online:
        df_online = _run_enrichr_online(
            gene_list=query_genes,
            gene_sets=gene_sets,
            species=species,
            pval_cutoff=pval_cutoff,
            background=background,
            top_n_terms=top_n_terms,
        )
        return EnrichmentResult(
            results=df_online,
            gene_list=query_genes,
            gene_sets=gene_sets,
            organism=species,
        )

    # Local offline ORA mode
    if isinstance(gene_sets, (str, Path, dict)):
        gene_sets_list = [gene_sets]
    else:
        gene_sets_list = list(gene_sets)

    all_dfs = []
    for gs_spec in gene_sets_list:
        try:
            lib_name, gs_dict = load_gene_set(
                gs_spec, organism=species, cache_dir=cache_dir
            )
            df_lib = _calculate_ora_local(
                gene_list=query_genes,
                gene_set_dict=gs_dict,
                background=background,
                pval_cutoff=pval_cutoff,
                lib_name=lib_name,
            )
            if not df_lib.empty:
                if top_n_terms is not None and top_n_terms > 0:
                    df_lib = df_lib.head(top_n_terms)
                all_dfs.append(df_lib)
        except Exception as e:
            log_warning(
                "EnrichmentAnalysis.Library",
                f"Failed to process gene set '{gs_spec}' ({type(e).__name__}): {e}",
            )

    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        # Re-sort across all libraries
        combined_df = combined_df.sort_values(
            by=["Adjusted P-value", "P-value"]
        ).reset_index(drop=True)
    else:
        combined_df = pd.DataFrame(columns=RESULT_COLUMNS)

    return EnrichmentResult(
        results=combined_df,
        gene_list=query_genes,
        gene_sets=gene_sets,
        organism=species,
    )


def gseapy_ora_enrichment_analysis(
    gene_list: list,
    gene_sets: list = config.ENRICHMENT_GENE_SETS,
    pval_cutoff: float = 0.05,
    species: str = "human",
    background: Optional[Union[int, List[str]]] = None,
    online: Optional[bool] = None,
) -> EnrichmentResult:
    """
    Perform ORA enrichment analysis.

    Backward-compatible wrapper for ``perform_ora_enrichment``.
    By default, executes local statistical enrichment without remote network calls.
    When ``online=True`` (or ``config.ENRICHMENT_ONLINE=True``), performs enrichment
    via the online Enrichr web API.

    Parameters
    ----------
    gene_list : list
        List of gene symbols.
    gene_sets : list or str or dict
        Gene set libraries or GMT paths/dicts.
    pval_cutoff : float, default=0.05
        Significance cutoff for Adjusted P-value.
    species : str, default='human'
        Target organism.
    background : int or list, optional
        Background gene universe.
    online : bool, optional
        Whether to run via online Enrichr API. Defaults to config.ENRICHMENT_ONLINE.

    Returns
    -------
    EnrichmentResult
        The result object containing ``.results`` DataFrame.
    """
    return perform_ora_enrichment(
        gene_list=gene_list,
        gene_sets=gene_sets,
        background=background,
        pval_cutoff=pval_cutoff,
        species=species,
        online=online,
    )
