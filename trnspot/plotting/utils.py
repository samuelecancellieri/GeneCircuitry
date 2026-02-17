"""
Plot Utilities for TRNspot
==========================

Provides centralized utilities for plot management including:
- Plot existence checking to avoid overwrites
- Plot logging to track generated figures
- Consistent save functions with configurable DPI
"""

import os
import json
from datetime import datetime
from typing import Optional, Dict, Any, List
from pathlib import Path
import matplotlib.pyplot as plt

from .. import config


class PlotLogger:
    """
    Logger for tracking generated plots.

    Maintains a registry of all plots generated during a pipeline run,
    including metadata like timestamps, file paths, and plot types.

    Parameters
    ----------
    output_dir : str
        Base output directory for the analysis.
    log_file : str, optional
        Name of the log file. Default is "plot_registry.json".

    Examples
    --------
    >>> logger = PlotLogger(output_dir="results/")
    >>> logger.register_plot("figures/qc/violin_pre_filter.png", "qc", {"step": "pre_filter"})
    >>> logger.save()
    """

    def __init__(
        self,
        output_dir: str,
        log_file: str = "plot_registry.json",
    ):
        self.output_dir = output_dir
        self.log_file = os.path.join(output_dir, "logs", log_file)
        self.registry: Dict[str, Dict[str, Any]] = {}
        self._load_existing()

    def _load_existing(self) -> None:
        """Load existing registry if available."""
        if os.path.exists(self.log_file):
            try:
                with open(self.log_file, "r") as f:
                    self.registry = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                import logging

                logging.getLogger("error").error(
                    f"[PlotLogger] Failed to load plot registry from "
                    f"'{self.log_file}' ({type(e).__name__}): {e}",
                    exc_info=True,
                )
                self.registry = {}

    def register_plot(
        self,
        filepath: str,
        plot_type: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Register a generated plot in the registry.

        Parameters
        ----------
        filepath : str
            Path to the generated plot file.
        plot_type : str
            Type of plot (e.g., "qc", "grn", "hotspot").
        metadata : dict, optional
            Additional metadata about the plot.
        """
        abs_path = os.path.abspath(filepath)
        rel_path = (
            os.path.relpath(filepath, self.output_dir) if self.output_dir else filepath
        )

        self.registry[rel_path] = {
            "absolute_path": abs_path,
            "plot_type": plot_type,
            "generated_at": datetime.now().isoformat(),
            "exists": os.path.exists(filepath),
            "metadata": metadata or {},
        }

    def is_registered(self, filepath: str) -> bool:
        """Check if a plot is already registered."""
        rel_path = (
            os.path.relpath(filepath, self.output_dir) if self.output_dir else filepath
        )
        return rel_path in self.registry

    def get_plots_by_type(self, plot_type: str) -> List[str]:
        """Get all registered plots of a specific type."""
        return [
            path
            for path, info in self.registry.items()
            if info.get("plot_type") == plot_type
        ]

    def save(self) -> None:
        """Save the registry to disk."""
        os.makedirs(os.path.dirname(self.log_file), exist_ok=True)
        with open(self.log_file, "w") as f:
            json.dump(self.registry, f, indent=2, default=str)

    def get_summary(self) -> Dict[str, int]:
        """Get summary count of plots by type."""
        summary: Dict[str, int] = {}
        for info in self.registry.values():
            plot_type = info.get("plot_type", "unknown")
            summary[plot_type] = summary.get(plot_type, 0) + 1
        return summary

    def __len__(self) -> int:
        return len(self.registry)

    def __repr__(self) -> str:
        return f"PlotLogger(plots={len(self)}, types={list(self.get_summary().keys())})"


# Global logger instance (lazily initialized)
_global_logger: Optional[PlotLogger] = None


def get_plot_logger(output_dir: Optional[str] = None) -> PlotLogger:
    """
    Get or create the global plot logger instance.

    Parameters
    ----------
    output_dir : str, optional
        Output directory. If None, uses config.OUTPUT_DIR.

    Returns
    -------
    PlotLogger
        The global plot logger instance.
    """
    global _global_logger

    if output_dir is None:
        output_dir = config.OUTPUT_DIR

    if _global_logger is None or _global_logger.output_dir != output_dir:
        _global_logger = PlotLogger(output_dir=output_dir)

    return _global_logger


def plot_exists(
    filepath: str,
    skip_existing: bool = True,
    verbose: bool = True,
) -> bool:
    """
    Check if a plot file already exists.

    Parameters
    ----------
    filepath : str
        Path to the plot file to check.
    skip_existing : bool
        If True, check for existence. If False, always return False
        (allowing overwrite).
    verbose : bool
        If True, print a message when skipping existing plots.

    Returns
    -------
    bool
        True if file exists and skip_existing is True, False otherwise.

    Examples
    --------
    >>> if plot_exists("figures/qc/plot.png"):
    ...     print("Plot already exists, skipping")
    ... else:
    ...     # Generate plot
    ...     pass
    """
    if not skip_existing:
        return False

    exists = os.path.exists(filepath)

    if exists and verbose:
        print(f"  Skipping existing: {os.path.basename(filepath)}")

    return exists


def save_plot(
    fig: plt.Figure,
    filepath: str,
    plot_type: str,
    dpi: Optional[int] = None,
    bbox_inches: str = "tight",
    metadata: Optional[Dict[str, Any]] = None,
    close_fig: bool = True,
    skip_existing: bool = True,
    verbose: bool = True,
) -> bool:
    """
    Save a matplotlib figure with logging and existence checking.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    filepath : str
        Path where the figure should be saved.
    plot_type : str
        Type of plot for logging (e.g., "qc", "grn", "hotspot").
    dpi : int, optional
        Resolution for saving. If None, uses config.SAVE_DPI.
    bbox_inches : str
        Bounding box setting for savefig. Default is "tight".
    metadata : dict, optional
        Additional metadata to log with the plot.
    close_fig : bool
        Whether to close the figure after saving. Default is True.
    skip_existing : bool
        If True, skip saving if file already exists. Default is True.
    verbose : bool
        If True, print status messages. Default is True.

    Returns
    -------
    bool
        True if plot was saved, False if skipped.

    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> ax.plot([1, 2, 3], [1, 2, 3])
    >>> save_plot(fig, "figures/test.png", "test", metadata={"step": "example"})
    """
    # Check if plot already exists
    if plot_exists(filepath, skip_existing=skip_existing, verbose=verbose):
        if close_fig:
            plt.close(fig)
        return False

    # Use config DPI if not specified
    if dpi is None:
        dpi = config.SAVE_DPI

    # Ensure directory exists
    os.makedirs(os.path.dirname(filepath) or ".", exist_ok=True)

    # Save the figure
    fig.savefig(filepath, dpi=dpi, bbox_inches=bbox_inches)

    if verbose:
        print(f"  Saved: {os.path.basename(filepath)}")

    # Register in the logger
    logger = get_plot_logger()
    logger.register_plot(filepath, plot_type, metadata)
    logger.save()

    # Close figure if requested
    if close_fig:
        plt.close(fig)

    return True


def get_plot_registry(output_dir: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
    """
    Get the current plot registry.

    Parameters
    ----------
    output_dir : str, optional
        Output directory. If None, uses config.OUTPUT_DIR.

    Returns
    -------
    dict
        The plot registry dictionary.
    """
    logger = get_plot_logger(output_dir)
    return logger.registry.copy()


def ensure_plot_dirs(output_dir: Optional[str] = None) -> None:
    """
    Ensure all plotting directories exist.

    Parameters
    ----------
    output_dir : str, optional
        Base output directory. If None, uses config.OUTPUT_DIR.
    """
    if output_dir is None:
        output_dir = config.OUTPUT_DIR

    dirs = [
        os.path.join(output_dir, "figures", "qc"),
        os.path.join(output_dir, "figures", "grn"),
        os.path.join(output_dir, "figures", "grn", "grn_deep_analysis"),
        os.path.join(output_dir, "figures", "hotspot"),
        os.path.join(output_dir, "logs"),
    ]

    for d in dirs:
        os.makedirs(d, exist_ok=True)
