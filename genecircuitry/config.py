"""
Configuration settings for GeneCircuitry package
"""

import numpy as np
import random
import os

# ============================================================================
# Random Seed Configuration
# ============================================================================

RANDOM_SEED = 42
"""Default random seed for reproducibility"""


def set_random_seed(seed: int = RANDOM_SEED):
    """
    Set random seed for reproducibility across all libraries.

    Parameters
    ----------
    seed : int, default=42
        Random seed value

    Examples
    --------
    >>> from genecircuitry.config import set_random_seed
    >>> set_random_seed(123)
    """
    np.random.seed(seed)
    random.seed(seed)

    print(f"Random seed set to: {seed}")


# ============================================================================
# Quality Control Default Parameters
# ============================================================================

QC_MIN_GENES = 100
"""Minimum number of genes expressed per cell"""

QC_MIN_COUNTS = 200
"""Minimum total counts per cell"""

QC_MAX_COUNTS = None
"""Maximum total counts per cell (None = no limit)"""

QC_PCT_MT_MAX = 20.0
"""Maximum percentage of mitochondrial counts"""

QC_MIN_CELLS = 3
"""Minimum number of cells expressing a gene"""


# ============================================================================
# Plotting Configuration
# ============================================================================

PLOT_DPI = 200
"""Default DPI for plot figures"""

SAVE_DPI = 600
"""DPI for saving figures"""

PLOT_FORMAT = "png"
"""Default format for saved figures"""

SAVE_PDF = True
"""Whether to automatically save figures in PDF format alongside PNG"""

PLOT_FIGSIZE_SQUARED = (6, 6)
"""Square figure size"""

PLOT_FIGSIZE_SQUARED_LARGE = (10, 10)
"""Large Square figure size"""

PLOT_FIGSIZE_SMALL = (6, 4)
"""Small figure size"""

PLOT_FIGSIZE_MEDIUM = (10, 7)
"""Medium figure size"""

PLOT_FIGSIZE_LARGE = (20, 15)
"""Large figure size"""

PLOT_FIGSIZE_WIDE = (20, 8)
"""Wide figure size (landscape layout)"""

PLOT_COLOR_PALETTE = "viridis"
"""Default color palette for plots"""


# ============================================================================
# Preprocessing Configuration
# ============================================================================

NORMALIZE_TARGET_SUM = 1e4
"""Target sum for count normalization (e.g., 10,000 for CPM-like normalization)"""

HVGS_N_TOP_GENES = 2000
"""Number of highly variable genes to select"""

HVGS_MIN_MEAN = 0.0125
"""Minimum mean expression for HVG selection"""

HVGS_MAX_MEAN = 3
"""Maximum mean expression for HVG selection"""

HVGS_MIN_DISP = 0.5
"""Minimum dispersion for HVG selection"""

PCA_N_COMPS = 50
"""Number of principal components to compute"""

PCA_SVD_SOLVER = "arpack"
"""SVD solver for PCA"""


# ============================================================================
# Neighborhood Graph Configuration
# ============================================================================

NEIGHBORS_N_NEIGHBORS = 15
"""Number of neighbors for KNN graph"""

NEIGHBORS_N_PCS = 40
"""Number of principal components to use for neighbor computation"""

NEIGHBORS_METHOD = "umap"
"""Method for computing connectivities (umap or gauss)"""

NEIGHBORS_METRIC = "euclidean"
"""Distance metric for neighbor search"""


# ============================================================================
# Clustering Configuration
# ============================================================================

LEIDEN_RESOLUTION = 1.0
"""Resolution parameter for Leiden clustering"""

LOUVAIN_RESOLUTION = 1.0
"""Resolution parameter for Louvain clustering"""

FORCE_DIM_REDUCTION = False
"""Whether to force re-run of dimensionality reduction and clustering"""


# ============================================================================
# UMAP Configuration
# ============================================================================

UMAP_MIN_DIST = 0.5
"""Minimum distance for UMAP"""

UMAP_SPREAD = 1.0
"""Spread parameter for UMAP"""

UMAP_N_COMPONENTS = 2
"""Number of dimensions for UMAP embedding"""


# ============================================================================
# Gene Regulatory Network Configuration
# ============================================================================

GRN_N_JOBS = 8
"""Number of parallel jobs for GRN inference (-1 = all cores)"""

GRN_MIN_TARGETS = 10
"""Minimum number of target genes for a TF to be considered"""

GRN_CONFIDENCE_THRESHOLD = 0.5
"""Minimum confidence score for regulatory interactions"""

GRN_CELL_DOWNSAMPLE = 30000
"""Number of cells to downsample to for GRN analysis (default: 20000, None or <=0 to disable)"""

# ============================================================================
# Hotspot Configuration
# ============================================================================

HOTSPOT_N_JOBS = 8
"""Number of parallel jobs for Hotspot inference (-1 = all cores)"""

HOTSPOT_N_NEIGHBORS = 30
"""Number of neighbors for KNN graph"""

HOTSPOT_FDR_THRESHOLD = 0.05
"""FDR threshold for selecting significant modules"""

HOTSPOT_MIN_GENES_PER_MODULE = 30
"""Minimum number of genes per module"""

HOTSPOT_CORE_ONLY = True
"""Whether to use only core genes in modules"""

HOTSPOT_TOP_GENES = 5000
"""Number of top genes to select for Hotspot analysis"""


# ============================================================================
# Enrichment Analysis Configuration
# ============================================================================

ENRICHMENT_GENE_SETS = ["MSigDB_Hallmark_2020"]
"""Gene set libraries for ORA enrichment analysis (e.g. gseapy enrichr).
Users can extend this list with additional libraries such as:
- "Reactome_Pathways_2024"
- "KEGG_2021_Human"
- "GO_Biological_Process_2023"
- "GO_Molecular_Function_2023"
"""

# ============================================================================
# ATAC Peaks Processing Configuration
# ============================================================================

ATAC_MOTIF_SCAN_FPR = 0.02
"""False positive rate for motif scanning in ATAC peak analysis"""

ATAC_MOTIF_SCORE_THRESHOLD = 10
"""Minimum motif score threshold for filtering enriched motifs"""


# ============================================================================
# File I/O Configuration
# ============================================================================

OUTPUT_DIR = "output"
"""Default output directory for results"""

CACHE_DIR = ".cache"
"""Directory for caching intermediate results"""

FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")
"""Default directory for saving figures"""

FIGURES_DIR_QC = os.path.join(FIGURES_DIR, "qc")
"""Directory for saving QC figures"""

FIGURES_DIR_GRN = os.path.join(FIGURES_DIR, "grn")
"""Directory for saving GRN figures"""

FIGURES_DIR_HOTSPOT = os.path.join(FIGURES_DIR, "hotspot")
"""Directory for saving Hotspot figures"""

FIGURES_DIR_COMPARATIVE = os.path.join(FIGURES_DIR, "comparative")
"""Directory for saving Comparative figures"""


# ============================================================================
# Logging Configuration
# ============================================================================

LOG_LEVEL = "INFO"
"""Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)"""

VERBOSE = True
"""Whether to print verbose output"""

LOG_PATH = os.path.join(OUTPUT_DIR, "logs")
"""Directory for saving log files"""

LOG_FILE = os.path.join(LOG_PATH, "genecircuitry_log.txt")
"""Path to the main log file"""


# ============================================================================
# Advanced Settings
# ============================================================================

N_JOBS = 8
"""Number of parallel jobs for computation (-1 = all cores)"""

CHUNK_SIZE = 1000
"""Chunk size for batch processing"""

LOW_MEMORY = False
"""Use low-memory mode (slower but uses less RAM)"""

# ============================================================================
# Scanpy Settings
# ============================================================================


def set_scanpy_settings():
    """
    Configure scanpy settings with GeneCircuitry defaults.

    This function sets common scanpy settings including figure parameters,
    verbosity, and caching.

    Examples
    --------
    >>> from genecircuitry.config import set_scanpy_settings
    >>> set_scanpy_settings()
    """
    try:
        import scanpy as sc
        import seaborn as sns
        import matplotlib as mpl

        mpl.rcdefaults()
        sns.reset_defaults()
        sns.set_theme(
            context="paper",
            style="white",
            rc={"figure.dpi": PLOT_DPI, "savefig.dpi": SAVE_DPI},
        )

        # Set figure parameters
        sc.settings.set_figure_params(
            dpi=PLOT_DPI,
            dpi_save=SAVE_DPI,
            frameon=False,
            vector_friendly=True,
            format=PLOT_FORMAT,
            facecolor="white",
        )

        # Set verbosity
        sc.settings.verbosity = 3 if VERBOSE else 1

        # Set number of jobs
        sc.settings.n_jobs = N_JOBS

        # Set cache directory
        sc.settings.cachedir = CACHE_DIR

        # Set figure dir
        sc.settings.figdir = FIGURES_DIR

        print("Scanpy settings configured")

    except ImportError:
        print("Warning: scanpy not installed, skipping scanpy configuration")


# ============================================================================
# Configuration Dictionary
# ============================================================================


def get_config():
    """
    Get all configuration parameters as a dictionary.

    Returns
    -------
    dict
        Dictionary containing all configuration parameters

    Examples
    --------
    >>> from genecircuitry.config import get_config
    >>> config = get_config()
    >>> print(config['RANDOM_SEED'])
    42
    """
    return {
        # Random seed
        "RANDOM_SEED": RANDOM_SEED,
        # QC parameters
        "QC_MIN_GENES": QC_MIN_GENES,
        "QC_MIN_COUNTS": QC_MIN_COUNTS,
        "QC_MAX_COUNTS": QC_MAX_COUNTS,
        "QC_PCT_MT_MAX": QC_PCT_MT_MAX,
        "QC_MIN_CELLS": QC_MIN_CELLS,
        # Plotting
        "PLOT_DPI": PLOT_DPI,
        "SAVE_DPI": SAVE_DPI,
        "PLOT_FORMAT": PLOT_FORMAT,
        "SAVE_PDF": SAVE_PDF,
        "PLOT_FIGSIZE_SQUARED": PLOT_FIGSIZE_SQUARED,
        "PLOT_FIGSIZE_SQUARED_LARGE": PLOT_FIGSIZE_SQUARED_LARGE,
        "PLOT_FIGSIZE_SMALL": PLOT_FIGSIZE_SMALL,
        "PLOT_FIGSIZE_MEDIUM": PLOT_FIGSIZE_MEDIUM,
        "PLOT_FIGSIZE_LARGE": PLOT_FIGSIZE_LARGE,
        "PLOT_FIGSIZE_WIDE": PLOT_FIGSIZE_WIDE,
        "PLOT_COLOR_PALETTE": PLOT_COLOR_PALETTE,
        # Preprocessing
        "NORMALIZE_TARGET_SUM": NORMALIZE_TARGET_SUM,
        "HVGS_N_TOP_GENES": HVGS_N_TOP_GENES,
        "HVGS_MIN_MEAN": HVGS_MIN_MEAN,
        "HVGS_MAX_MEAN": HVGS_MAX_MEAN,
        "HVGS_MIN_DISP": HVGS_MIN_DISP,
        "PCA_N_COMPS": PCA_N_COMPS,
        "PCA_SVD_SOLVER": PCA_SVD_SOLVER,
        # Neighbors
        "NEIGHBORS_N_NEIGHBORS": NEIGHBORS_N_NEIGHBORS,
        "NEIGHBORS_N_PCS": NEIGHBORS_N_PCS,
        "NEIGHBORS_METHOD": NEIGHBORS_METHOD,
        "NEIGHBORS_METRIC": NEIGHBORS_METRIC,
        # Clustering
        "LEIDEN_RESOLUTION": LEIDEN_RESOLUTION,
        "LOUVAIN_RESOLUTION": LOUVAIN_RESOLUTION,
        "FORCE_DIM_REDUCTION": FORCE_DIM_REDUCTION,
        # UMAP
        "UMAP_MIN_DIST": UMAP_MIN_DIST,
        "UMAP_SPREAD": UMAP_SPREAD,
        "UMAP_N_COMPONENTS": UMAP_N_COMPONENTS,
        # GRN
        "GRN_N_JOBS": GRN_N_JOBS,
        "GRN_MIN_TARGETS": GRN_MIN_TARGETS,
        "GRN_CONFIDENCE_THRESHOLD": GRN_CONFIDENCE_THRESHOLD,
        "GRN_CELL_DOWNSAMPLE": GRN_CELL_DOWNSAMPLE,
        # Hotspot
        "HOTSPOT_N_JOBS": HOTSPOT_N_JOBS,
        "HOTSPOT_N_NEIGHBORS": HOTSPOT_N_NEIGHBORS,
        "HOTSPOT_FDR_THRESHOLD": HOTSPOT_FDR_THRESHOLD,
        "HOTSPOT_MIN_GENES_PER_MODULE": HOTSPOT_MIN_GENES_PER_MODULE,
        "HOTSPOT_CORE_ONLY": HOTSPOT_CORE_ONLY,
        "HOTSPOT_TOP_GENES": HOTSPOT_TOP_GENES,
        # Enrichment
        "ENRICHMENT_GENE_SETS": ENRICHMENT_GENE_SETS,
        # ATAC Peaks
        "ATAC_MOTIF_SCAN_FPR": ATAC_MOTIF_SCAN_FPR,
        "ATAC_MOTIF_SCORE_THRESHOLD": ATAC_MOTIF_SCORE_THRESHOLD,
        # File I/O
        "OUTPUT_DIR": OUTPUT_DIR,
        "CACHE_DIR": CACHE_DIR,
        "FIGURES_DIR": FIGURES_DIR,
        "FIGURES_DIR_QC": FIGURES_DIR_QC,
        "FIGURES_DIR_GRN": FIGURES_DIR_GRN,
        "FIGURES_DIR_HOTSPOT": FIGURES_DIR_HOTSPOT,
        "FIGURES_DIR_COMPARATIVE": FIGURES_DIR_COMPARATIVE,
        # Logging
        "LOG_LEVEL": LOG_LEVEL,
        "VERBOSE": VERBOSE,
        "LOG_PATH": LOG_PATH,
        "LOG_FILE": LOG_FILE,
        # Advanced
        "N_JOBS": N_JOBS,
        "CHUNK_SIZE": CHUNK_SIZE,
        "LOW_MEMORY": LOW_MEMORY,
    }


def print_config():
    """
    Print all configuration parameters.

    Examples
    --------
    >>> from genecircuitry.config import print_config
    >>> print_config()
    """
    config = get_config()
    print("=" * 60)
    print("GeneCircuitry Configuration")
    print("=" * 60)
    for key, value in config.items():
        print(f"{key:30s} = {value}")
    print("=" * 60)


def update_config(**kwargs):
    """
    Update configuration parameters.

    Parameters
    ----------
    **kwargs
        Configuration parameters to update

    Examples
    --------
    >>> from genecircuitry.config import update_config
    >>> update_config(RANDOM_SEED=123, QC_MIN_GENES=300)
    """
    global_vars = globals()
    for key, value in kwargs.items():
        if key in global_vars:
            global_vars[key] = value
            print(f"Updated {key} = {value}")
        else:
            print(f"Warning: {key} is not a valid configuration parameter")


def load_config_file(path: str) -> dict:
    """
    Load configuration parameters from a YAML or JSON file and apply them.

    Keys in the file must match existing config parameter names (e.g.
    ``QC_MIN_GENES``, ``LEIDEN_RESOLUTION``).  Unknown keys raise a
    ``ValueError`` so typos are caught early.

    CLI arguments always take precedence over values set by this function
    when the pipeline is run via the command line.

    Parameters
    ----------
    path : str
        Path to a ``.yaml`` / ``.yml`` or ``.json`` config file.

    Returns
    -------
    dict
        The parameters that were loaded and applied.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file format is unsupported or contains unknown parameter names.
    ImportError
        If a YAML file is provided but ``pyyaml`` is not installed.

    Examples
    --------
    >>> from genecircuitry.config import load_config_file
    >>> load_config_file("my_run.yaml")
    {'QC_MIN_GENES': 200, 'LEIDEN_RESOLUTION': 0.5}
    """
    import json as _json

    if not os.path.exists(path):
        raise FileNotFoundError(f"Config file not found: {path}")

    ext = os.path.splitext(path)[1].lower()

    if ext in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError:
            raise ImportError(
                "pyyaml is required to load YAML config files. "
                "Install it with: pip install pyyaml"
            )
        with open(path) as f:
            params = yaml.safe_load(f)
    elif ext == ".json":
        with open(path) as f:
            params = _json.load(f)
    else:
        raise ValueError(
            f"Unsupported config file format: {ext!r}. Use .yaml, .yml, or .json"
        )

    if params is None:
        return {}

    if not isinstance(params, dict):
        raise ValueError(
            "Config file must contain a mapping of parameter names to values"
        )

    valid_keys = set(get_config().keys())
    unknown = set(params.keys()) - valid_keys
    if unknown:
        raise ValueError(
            f"Unknown config parameters: {sorted(unknown)}.\n"
            f"Valid parameters are listed in genecircuitry/config.py or via "
            f"genecircuitry.config.get_config()."
        )

    update_config(**params)
    return params
