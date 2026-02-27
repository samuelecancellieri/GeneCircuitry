"""
GeneCircuitry - A package for transcriptional regulatory network analysis
===================================================================

GeneCircuitry is a Python package for transcriptional regulatory network (TRN)
analysis using single-cell data. It integrates Scanpy, CellOracle, and
Hotspot into a modular, checkpoint-enabled pipeline.

Modules
-------
- preprocessing: Scanpy wrappers for QC and normalization
- celloracle_processing: GRN inference with CellOracle (optional dep)
- hotspot_processing: Gene module analysis with Hotspot (optional dep)
- grn_deep_analysis: Network visualization and analysis
- reporting: HTML and PDF report generation
- plotting: Centralized plot generation with logging
- pipeline: Pipeline orchestration (PipelineController)
- config: Central configuration and parameters
"""

__version__ = "0.1.0"
__author__ = "Samuele Cancellieri"

# Import core modules (always available)
from . import preprocessing
from . import grn_deep_analysis
from . import atac_peaks_processing
from . import config
from . import reporting
from . import plotting
from . import logging_utils

# Lazy import for pipeline module (avoids circular import since
# controller.py imports from genecircuitry.config and genecircuitry.preprocessing)
import importlib as _importlib


def __getattr__(name):
    if name == "pipeline":
        return _importlib.import_module(".pipeline", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# Lazy imports for optional dependencies (CellOracle, Hotspot)
# These will only fail when actually accessed, not on `import genecircuitry`
try:
    from . import celloracle_processing
except ImportError as _co_err:
    import warnings as _warnings
    _warnings.warn(
        f"CellOracle optional dependency could not be loaded — GRN inference will be unavailable.\n"
        f"  Missing module: {_co_err.name!r}\n"
        f"  Full error   : {_co_err}\n"
        f"  To fix       : pip install celloracle",
        ImportWarning,
        stacklevel=2,
    )
    celloracle_processing = None  # type: ignore[assignment]

try:
    from . import hotspot_processing
except ImportError as _hs_err:
    import warnings as _warnings
    _warnings.warn(
        f"Hotspot optional dependency could not be loaded — gene-module detection will be unavailable.\n"
        f"  Missing module: {_hs_err.name!r}\n"
        f"  Full error   : {_hs_err}\n"
        f"  To fix       : pip install hotspot-sc",
        ImportWarning,
        stacklevel=2,
    )
    hotspot_processing = None  # type: ignore[assignment]

# Import commonly used functions
from .config import set_random_seed, set_scanpy_settings, get_config, print_config

# Import reporting functions
from .reporting import (
    ReportGenerator,
    generate_report,
    generate_html_report,
    generate_pdf_report,
)

# Import plotting utilities
from .plotting import (
    PlotLogger,
    get_plot_logger,
    plot_exists,
    save_plot,
)

__all__ = [
    # Modules
    "preprocessing",
    "celloracle_processing",
    "hotspot_processing",
    "grn_deep_analysis",
    "config",
    "reporting",
    "plotting",
    "pipeline",
    # Config functions
    "set_random_seed",
    "set_scanpy_settings",
    "get_config",
    "print_config",
    # Reporting functions
    "ReportGenerator",
    "generate_report",
    "generate_html_report",
    "generate_pdf_report",
    # Plotting utilities
    "PlotLogger",
    "get_plot_logger",
    "plot_exists",
    "save_plot",
]
