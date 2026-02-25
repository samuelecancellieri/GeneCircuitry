"""
GeneCircuitry Pipeline Module
========================

Central pipeline orchestration for GeneCircuitry analysis workflows.

This module provides the PipelineController class and supporting functions
for running the complete GeneCircuitry analysis pipeline.

Previously located at examples/complete_pipeline.py, this is the core
pipeline engine — not an example.
"""

from .controller import (
    PipelineController,
    setup_directories,
    setup_logging,
    log_step,
    log_error,
    main,
)

__all__ = [
    "PipelineController",
    "setup_directories",
    "setup_logging",
    "log_step",
    "log_error",
    "main",
]
