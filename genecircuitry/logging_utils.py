"""
Centralized logging utilities for GeneCircuitry modules.

All GeneCircuitry submodules should import `log_error` and `log_warning` from
here to ensure errors/warnings are recorded in the pipeline's error.log
and pipeline.log files regardless of which module raises them.

The loggers are the same instances created by
``genecircuitry.pipeline.controller.setup_logging`` — they share the
``"pipeline"`` and ``"error"`` logger names so that handlers registered
by the controller are reused automatically.
"""

import logging
import traceback


def get_pipeline_logger():
    """Return the shared pipeline logger (INFO-level step tracking)."""
    return logging.getLogger("pipeline")


def get_error_logger():
    """Return the shared error logger (ERROR-level with tracebacks)."""
    return logging.getLogger("error")


def log_error(context: str, exception: Exception) -> None:
    """
    Log an error with full traceback to both error.log and pipeline.log.

    This is safe to call even if ``setup_logging`` has not been invoked
    yet — in that case the messages are simply discarded (no handlers).

    Parameters
    ----------
    context : str
        Human-readable description of *where* the error occurred,
        e.g. ``"Hotspot.CreateObject"`` or ``"GRNPlotting.scatter"``.
    exception : Exception
        The exception instance.
    """
    error_logger = get_error_logger()
    pipeline_logger = get_pipeline_logger()

    msg = f"[{context}] {type(exception).__name__}: {exception}"

    if error_logger.handlers:
        error_logger.error(msg, exc_info=True)
    if pipeline_logger.handlers:
        pipeline_logger.error(msg)


def log_warning(context: str, message: str) -> None:
    """
    Log a warning to the pipeline log.

    Parameters
    ----------
    context : str
        Module/step context string.
    message : str
        Warning description.
    """
    pipeline_logger = get_pipeline_logger()
    if pipeline_logger.handlers:
        pipeline_logger.warning(f"[{context}] {message}")
