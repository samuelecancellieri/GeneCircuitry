"""
GeneCircuitry - A package for transcriptional regulatory network analysis
===================================================================

GeneCircuitry is a Python package for transcriptional regulatory network (TRN)
analysis using single-cell data. It integrates Scanpy, CellOracle, and
Hotspot into a modular, checkpoint-enabled pipeline.

Modules
-------
- preprocessing: Scanpy wrappers for QC and normalization
- celloracle_processing: GRN inference with CellOracle
- hotspot_processing: Gene module analysis with Hotspot
- grn_deep_analysis: Network visualization and analysis
- reporting: HTML and PDF report generation
- plotting: Centralized plot generation with logging
- pipeline: Pipeline orchestration (PipelineController)
- config: Central configuration and parameters
"""

from importlib.metadata import (
    version as _pkg_version,
    PackageNotFoundError as _PkgNotFound,
)

try:
    __version__ = _pkg_version("genecircuitry")
except _PkgNotFound:
    __version__ = "0.1.2"  # fallback for editable / source installs
__author__ = "Samuele Cancellieri"

# Import core modules (always available)
from . import preprocessing
from . import grn_deep_analysis
from . import atac_peaks_processing
from . import config
from . import reporting
from . import plotting
from . import logging_utils
from . import celloracle_processing
from . import hotspot_processing

# Lazy import for pipeline module (avoids circular import since
# controller.py imports from genecircuitry.config and genecircuitry.preprocessing)
import importlib as _importlib


def __getattr__(name):
    if name == "pipeline":
        return _importlib.import_module(".pipeline", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


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
