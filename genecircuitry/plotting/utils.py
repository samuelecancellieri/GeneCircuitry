"""
Plot Utilities for GeneCircuitry
==========================

Provides centralized utilities for plot management including:
- Plot existence checking to avoid overwrites
- Plot logging to track generated figures
- Consistent save functions with configurable DPI
"""

import os
import json
from datetime import datetime
from typing import Optional, Dict, Any, List, Tuple, Callable, Union
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

    def register_plots_batch(
        self,
        plots: List[Tuple[str, str, Optional[Dict[str, Any]]]],
    ) -> None:
        """
        Register multiple generated plots in the registry at once.

        Parameters
        ----------
        plots : list of tuple
            List of (filepath, plot_type, metadata) tuples.
        """
        for filepath, plot_type, metadata in plots:
            self.register_plot(filepath, plot_type, metadata)

    def save(self) -> None:
        """Save the registry to disk atomically."""
        log_dir = os.path.dirname(self.log_file)
        os.makedirs(log_dir, exist_ok=True)
        # Reload latest from disk if available to merge parallel updates
        if os.path.exists(self.log_file):
            try:
                with open(self.log_file, "r") as f:
                    disk_registry = json.load(f)
                if isinstance(disk_registry, dict):
                    disk_registry.update(self.registry)
                    self.registry = disk_registry
            except Exception:
                pass
        # Write atomically using a temporary file
        temp_file = f"{self.log_file}.tmp.{os.getpid()}"
        try:
            with open(temp_file, "w") as f:
                json.dump(self.registry, f, indent=2, default=str)
            os.replace(temp_file, self.log_file)
        except Exception:
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except OSError:
                    pass
            # Fallback direct write
            with open(self.log_file, "w") as f:
                json.dump(self.registry, f, indent=2, default=str)

    def get_summary(self) -> Dict[str, int]:
        """Get summary count of plots by type."""
        summary: Dict[str, int] = {}
        for info in self.registry.values():
            plot_type = info.get("plot_type", "unknown")
            summary[plot_type] = summary.get(plot_type, 0) + 1
        return summary

    def get_plot_count(self) -> int:
        """Get total number of registered plots."""
        return len(self.registry)

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
    check_pdf: Optional[bool] = None,
) -> bool:
    """
    Check if a plot file already exists (including PDF counterpart if enabled).

    Parameters
    ----------
    filepath : str
        Path to the plot file to check.
    skip_existing : bool
        If True, check for existence. If False, always return False
        (allowing overwrite).
    verbose : bool
        If True, print a message when skipping existing plots.
    check_pdf : bool, optional
        If True, also require the .pdf counterpart to exist. Defaults to config.SAVE_PDF.

    Returns
    -------
    bool
        True if file exists (and PDF counterpart if required) and skip_existing is True, False otherwise.

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

    if check_pdf is None:
        check_pdf = getattr(config, "SAVE_PDF", False)

    exists = os.path.exists(filepath)

    if exists and check_pdf:
        base, ext = os.path.splitext(filepath)
        if ext.lower() != ".pdf":
            pdf_path = f"{base}.pdf"
            if not os.path.exists(pdf_path):
                exists = False

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
    save_pdf: Optional[bool] = None,
) -> bool:
    """
    Save a matplotlib figure with logging, existence checking, and dual PNG/PDF output.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    filepath : str
        Path where the figure should be saved.
    plot_type : str
        Type of plot for logging (e.g., "qc", "grn", "hotspot", "comparative").
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
    save_pdf : bool, optional
        Whether to also save a PDF version. If None, uses config.SAVE_PDF.

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
    if save_pdf is None:
        save_pdf = getattr(config, "SAVE_PDF", True)

    if plot_exists(
        filepath,
        skip_existing=skip_existing,
        verbose=verbose,
        check_pdf=save_pdf,
    ):
        if close_fig:
            plt.close(fig)
        return False

    # Use config DPI if not specified
    if dpi is None:
        dpi = getattr(config, "SAVE_DPI", 600)

    # Ensure directory exists
    os.makedirs(os.path.dirname(filepath) or ".", exist_ok=True)

    base, ext = os.path.splitext(filepath)
    if not ext:
        ext = f".{getattr(config, 'PLOT_FORMAT', 'png')}"
        filepath = f"{base}{ext}"

    # Save primary figure
    fig.savefig(filepath, dpi=dpi, bbox_inches=bbox_inches)
    if verbose:
        print(f"  Saved plot: {filepath}")

    # Register in the logger
    logger = get_plot_logger()
    logger.register_plot(filepath, plot_type, metadata)

    # Save PDF version if enabled
    if save_pdf and ext.lower() != ".pdf":
        pdf_path = f"{base}.pdf"
        fig.savefig(pdf_path, dpi=dpi, bbox_inches=bbox_inches)
        logger.register_plot(pdf_path, plot_type, metadata)
        if verbose:
            print(f"  Saved plot (PDF): {pdf_path}")

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
        os.path.join(output_dir, "figures", "comparative"),
        os.path.join(output_dir, "comparative"),
        os.path.join(output_dir, "logs"),
    ]

    for d in dirs:
        os.makedirs(d, exist_ok=True)


def _resolve_n_jobs(n_jobs: Optional[int] = None) -> int:
    """
    Resolve effective number of parallel workers.

    Parameters
    ----------
    n_jobs : int, optional
        Requested worker count. If None, defaults to config.N_JOBS (or 1).
        If -1, uses os.cpu_count().

    Returns
    -------
    int
        Positive integer worker count >= 1.
    """
    if n_jobs is None:
        n_jobs = getattr(config, "N_JOBS", 1)
    if n_jobs == -1:
        n_jobs = os.cpu_count() or 1
    try:
        n_jobs = int(n_jobs)
    except (ValueError, TypeError):
        n_jobs = 1
    return max(1, n_jobs)


def run_parallel_tasks(
    worker_fn,
    tasks: List[Any],
    n_jobs: Optional[int] = None,
    desc: Optional[str] = None,
) -> List[Any]:
    """
    Execute a worker function across tasks either sequentially or in parallel.

    Parameters
    ----------
    worker_fn : callable
        Module-level picklable worker function accepting a single argument.
    tasks : list
        List of task inputs passed one-by-one to worker_fn.
    n_jobs : int, optional
        Number of worker processes. If None, uses config.N_JOBS.
        If 1 or len(tasks) <= 1, executes sequentially.
    desc : str, optional
        Description of tasks for progress/logging.

    Returns
    -------
    list
        Results collected in the same order as tasks.
    """
    if not tasks:
        return []

    effective_jobs = _resolve_n_jobs(n_jobs)

    if effective_jobs <= 1 or len(tasks) <= 1:
        results = []
        for task in tasks:
            try:
                results.append(worker_fn(task))
            except Exception as e:
                import logging

                logging.getLogger("error").error(
                    f"Task execution failed ({type(e).__name__}): {e}",
                    exc_info=True,
                )
                results.append(None)
        return results

    import concurrent.futures

    max_workers = min(effective_jobs, len(tasks))
    results = [None] * len(tasks)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {
            executor.submit(worker_fn, task): idx
            for idx, task in enumerate(tasks)
        }
        for future in concurrent.futures.as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                import logging

                logging.getLogger("error").error(
                    f"Parallel task {idx} failed ({type(e).__name__}): {e}",
                    exc_info=True,
                )
                results[idx] = None

    return results
