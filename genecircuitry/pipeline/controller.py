#!/usr/bin/env python
"""
Complete GeneCircuitry Analysis Pipeline
====================================
This script runs the full analysis workflow:
1. Data loading and quality control
2. Normalization and preprocessing
3. Dimensionality reduction and clustering
4. GRN inference with CellOracle
5. Gene module identification with Hotspot

Usage:
    python -m genecircuitry.pipeline --input data/paul15/paul15.h5 --output output

Or with default paths:
    python -m genecircuitry.pipeline
"""

import os
import shutil
import argparse
import sys
from datetime import datetime
from typing import Optional, List, Dict, Any, Union
from collections.abc import Iterable, Sequence
import scanpy as sc
import pandas as pd
import numpy as np
import pickle
import hashlib
import json
import logging
import traceback

# Import genecircuitry modules (use direct submodule imports to avoid circular import)
from genecircuitry.config import set_random_seed, set_scanpy_settings
from genecircuitry import config
from genecircuitry.preprocessing import (
    perform_qc,
    perform_normalization,
    perform_dimensionality_reduction_clustering,
    ensure_categorical_obs,
    resolve_cluster_key,
    resolve_cluster_key_name,
    parse_cluster_keys,
    sanitize_identifier,
)
from genecircuitry.reporting import generate_report, generate_stratified_report


# Global logger instances
pipeline_logger = None
error_logger = None


def setup_logging(output_dir):
    """
    Setup logging system for pipeline execution and errors.

    Creates two log files:
    - pipeline.log: Records all pipeline steps with timestamps
    - error.log: Records all errors with full tracebacks

    Parameters:
    -----------
    output_dir : str
        Directory where log files will be created
    """
    global pipeline_logger, error_logger

    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    # Pipeline logger - tracks all steps
    pipeline_logger = logging.getLogger("pipeline")
    pipeline_logger.setLevel(logging.INFO)
    pipeline_logger.handlers.clear()

    pipeline_handler = logging.FileHandler(
        os.path.join(log_dir, "pipeline.log"), mode="a"
    )
    pipeline_handler.setLevel(logging.INFO)
    pipeline_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    pipeline_handler.setFormatter(pipeline_formatter)
    pipeline_logger.addHandler(pipeline_handler)

    # Error logger - tracks all errors with tracebacks
    error_logger = logging.getLogger("error")
    error_logger.setLevel(logging.ERROR)
    error_logger.handlers.clear()

    error_handler = logging.FileHandler(os.path.join(log_dir, "error.log"), mode="a")
    error_handler.setLevel(logging.ERROR)
    error_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s\n%(exc_info)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    error_handler.setFormatter(error_formatter)
    error_logger.addHandler(error_handler)

    # Log the start of a new session
    pipeline_logger.info("=" * 70)
    pipeline_logger.info("NEW PIPELINE SESSION STARTED")
    pipeline_logger.info("=" * 70)

    return pipeline_logger, error_logger


def log_step(step_name, status="STARTED", details=None):
    """
    Log a pipeline step with timestamp.

    Parameters:
    -----------
    step_name : str
        Name of the pipeline step
    status : str
        Status of the step (STARTED, COMPLETED, FAILED, etc.)
    details : dict, optional
        Additional details to log (e.g., n_obs, n_vars)
    """
    if pipeline_logger is None:
        return

    message = f"[{step_name}] {status}"
    if details:
        detail_str = ", ".join([f"{k}={v}" for k, v in details.items()])
        message += f" - {detail_str}"

    pipeline_logger.info(message)


def log_error(error_context, exception):
    """
    Log an error with full traceback.

    Parameters:
    -----------
    error_context : str
        Description of where/when the error occurred
    exception : Exception
        The exception that was raised
    """
    if error_logger is None:
        return

    error_logger.error(f"ERROR in {error_context}: {str(exception)}", exc_info=True)

    # Also log to pipeline logger
    if pipeline_logger:
        pipeline_logger.error(f"ERROR in {error_context}: {str(exception)}")


def compute_input_hash(input_path, **params):
    """Compute hash of input file and parameters for checkpoint verification."""
    hash_obj = hashlib.md5()

    # Hash input file path
    if input_path:
        hash_obj.update(str(input_path).encode())
        # If file exists, hash its modification time
        if os.path.exists(input_path):
            mtime = os.path.getmtime(input_path)
            hash_obj.update(str(mtime).encode())

    # Hash relevant parameters
    for key, value in sorted(params.items()):
        hash_obj.update(f"{key}:{value}".encode())

    return hash_obj.hexdigest()


def write_checkpoint(log_dir, step_name, input_hash, **metadata):
    """Write checkpoint log file for a completed step."""
    os.makedirs(log_dir, exist_ok=True)
    checkpoint_file = os.path.join(log_dir, f"{step_name}.checkpoint")

    checkpoint_data = {
        "step_name": step_name,
        "input_hash": input_hash,
        "timestamp": datetime.now().isoformat(),
        "metadata": metadata,
    }

    with open(checkpoint_file, "w") as f:
        json.dump(checkpoint_data, f, indent=2)


def check_checkpoint(log_dir, step_name, input_hash):
    """Check if step has already been completed with same input."""
    checkpoint_file = os.path.join(log_dir, f"{step_name}.checkpoint")

    if not os.path.exists(checkpoint_file):
        return False

    try:
        with open(checkpoint_file, "r") as f:
            checkpoint_data = json.load(f)

        # Check if input hash matches
        if checkpoint_data.get("input_hash") == input_hash:
            print(f"  ⏭ Checkpoint found - skipping {step_name}")
            print(f"     (completed at {checkpoint_data.get('timestamp')})")
            return True
    except Exception as e:
        log_error(f"Checkpoint.Read({step_name})", e)
        print(f"  ⚠ Error reading checkpoint for '{step_name}': {e}")

    return False


class PipelineController:
    """Controller for managing GeneCircuitry pipeline execution."""

    def __init__(self, args, start_time):
        """Initialize pipeline controller with arguments."""
        self.args = args
        self.start_time = start_time
        self.log_dir = os.path.join(args.output, "logs")

        # Data containers
        self.adata = None
        self.adata_preprocessed = None
        self.adata_list = []
        self.adata_stratification_list = []
        self.stratification_results = []
        self.atac_peaks_pkl = None

        log_step(
            "PipelineController",
            "INITIALIZED",
            {
                "output_dir": args.output,
                "stratification_key": args.cluster_key_stratification,
            },
        )

    def run_step_load(self):
        """Execute Step 1: Data Loading."""
        log_step("Controller.LoadData", "STARTED")
        try:
            self.adata = load_data(self.args.input)
            log_step("Controller.LoadData", "COMPLETED")
            return self.adata
        except Exception as e:
            log_error("Controller.LoadData", e)
            raise

    def run_step_preprocessing(self):
        """Execute Step 2: Preprocessing."""
        log_step("Controller.Preprocessing", "STARTED")
        try:
            if self.adata is None:
                raise ValueError("Data not loaded. Run step_load first.")
            self.adata_preprocessed = preprocessing_pipeline(
                self.adata,
                self.args.name,
                skip_qc=self.args.skip_qc,
                log_dir=self.log_dir,
            )
            log_step("Controller.Preprocessing", "COMPLETED")
            return self.adata_preprocessed
        except Exception as e:
            log_error("Controller.Preprocessing", e)
            raise

    def run_step_stratification(self):
        """Execute Step 2.5: Stratification."""
        log_step("Controller.Stratification", "STARTED")
        try:
            if self.adata_preprocessed is None:
                raise ValueError("Data not preprocessed. Run step_preprocessing first.")
            self.adata_list, self.adata_stratification_list = stratification_pipeline(
                self.adata_preprocessed,
                self.args.cluster_key_stratification,
                self.args.clusters,
            )
            log_step(
                "Controller.Stratification",
                "COMPLETED",
                {"n_stratifications": len(self.adata_list)},
            )
            return self.adata_list, self.adata_stratification_list
        except Exception as e:
            log_error("Controller.Stratification", e)
            raise

    def run_step_clustering(self, adata=None, log_dir=None, force=None):
        """Execute Step 3: Clustering."""
        log_step("Controller.Clustering", "STARTED")
        try:
            if adata is None:
                adata = self.adata_preprocessed
            if adata is None:
                raise ValueError("No data available for clustering.")
            if log_dir is None:
                log_dir = self.log_dir
            if force is None:
                force = getattr(self.args, "force_dim_reduction", False) or getattr(
                    self.args, "force_dimensionality_reduction", False
                ) or config.FORCE_DIM_REDUCTION

            cluster_key = getattr(self.args, "cluster_key", "leiden")
            parsed = parse_cluster_keys(cluster_key)
            if len(parsed) > 1 and all(k in adata.obs.columns for k in parsed):
                adata, cluster_key = resolve_cluster_key(adata, cluster_key)

            result = dimensionality_reduction_clustering(
                adata,
                cluster_key=cluster_key,
                log_dir=log_dir,
                force=force,
            )
            log_step("Controller.Clustering", "COMPLETED")
            return result
        except Exception as e:
            log_error("Controller.Clustering", e)
            raise

    def run_step_atac_peaks(self, log_dir=None):
        """Execute ATAC Peaks Processing step.

        Processes a BED file with pre-called ATAC peaks through CellOracle
        motif analysis to generate an enriched TF info matrix (PKL). The
        resulting PKL path is stored so CellOracle can use it as a custom
        base GRN instead of the default promoter-based one.
        """
        log_step("Controller.ATACPeaks", "STARTED")
        try:
            if log_dir is None:
                log_dir = self.log_dir

            bed_path = self.args.atac_peaks
            if not bed_path:
                log_step("Controller.ATACPeaks", "SKIPPED", {"reason": "no BED file"})
                return None

            output_dir = getattr(self.args, "output", config.OUTPUT_DIR)

            print(f"\n{'='*70}")
            print("STEP 3.5: ATAC Peaks Processing")
            print(f"{'='*70}")
            print(f"  BED file: {bed_path}")
            print(f"  Species: {self.args.species}")

            # Check checkpoint
            step_hash = compute_input_hash(
                bed_path,
                species=self.args.species,
                fpr=config.ATAC_MOTIF_SCAN_FPR,
                threshold=config.ATAC_MOTIF_SCORE_THRESHOLD,
            )

            dict_pkl_path = os.path.join(
                output_dir, "celloracle", "enriched_atac_peaks_dict.pkl"
            )
            df_pkl_path = os.path.join(
                output_dir, "celloracle", "enriched_atac_peaks_df.pkl"
            )
            if (
                log_dir
                and check_checkpoint(log_dir, "atac_peaks", step_hash)
                and os.path.exists(dict_pkl_path)
            ):
                log_step(
                    "Controller.ATACPeaks",
                    "LOADED_FROM_CHECKPOINT",
                    {"pkl_path": dict_pkl_path, "df_pkl_path": df_pkl_path},
                )
                print(f"  ✓ Loading enriched ATAC peaks from checkpoint: {dict_pkl_path}")
                self.atac_peaks_pkl = dict_pkl_path
                return dict_pkl_path

            from genecircuitry.atac_peaks_processing import process_atac_peaks

            dict_pkl_path = process_atac_peaks(
                bed_path=bed_path,
                species=self.args.species,
                output_dir=output_dir,
                log_dir=log_dir,
            )

            # Save checkpoint
            if log_dir:
                write_checkpoint(
                    log_dir,
                    "atac_peaks",
                    step_hash,
                    bed_path=bed_path,
                    pkl_path=dict_pkl_path,
                    df_pkl_path=df_pkl_path,
                )

            self.atac_peaks_pkl = dict_pkl_path
            log_step(
                "Controller.ATACPeaks",
                "COMPLETED",
                {"pkl_path": dict_pkl_path},
            )
            return dict_pkl_path
        except Exception as e:
            log_error("Controller.ATACPeaks", e)
            raise

    def run_step_celloracle(
        self, adata, log_dir=None, hotspot_genes_path=None, cluster_key=None
    ):
        """Execute Step 4: CellOracle."""
        log_step("Controller.CellOracle", "STARTED")
        try:
            if log_dir is None:
                log_dir = self.log_dir
            if cluster_key is None:
                cluster_key = getattr(self.args, "cluster_key", "leiden")

            parsed = parse_cluster_keys(cluster_key)
            if len(parsed) > 1 and all(k in adata.obs.columns for k in parsed):
                adata, cluster_key = resolve_cluster_key(adata, cluster_key)

            # Determine TF dictionary: ATAC peaks PKL takes
            # priority over --tf-dictionary CLI argument
            tf_dictionary = self.args.tf_dictionary
            if hasattr(self, "atac_peaks_pkl") and self.atac_peaks_pkl:
                tf_dictionary = self.atac_peaks_pkl
                log_step(
                    "Controller.CellOracle",
                    "USING_ATAC_PEAKS_AS_TF_DICT",
                    {"tf_dictionary": tf_dictionary},
                )

            result = celloracle_pipeline(
                adata,
                cluster_key=cluster_key,
                species=self.args.species,
                raw_count_layer=self.args.raw_count_layer,
                embedding_name=self.args.embedding_grn,
                TG_to_TF_dictionary=tf_dictionary,
                skip_celloracle=self.args.skip_celloracle,
                no_base_grn=self.args.no_base_grn,
                hotspot_genes_path=hotspot_genes_path,
                log_dir=log_dir,
            )
            log_step("Controller.CellOracle", "COMPLETED")
            return result
        except Exception as e:
            log_error("Controller.CellOracle", e)
            raise

    def run_step_hotspot(self, adata, log_dir=None, cluster_key=None):
        """Execute Step 5: Hotspot."""
        log_step("Controller.Hotspot", "STARTED")
        try:
            if log_dir is None:
                log_dir = self.log_dir
            if cluster_key is None:
                cluster_key = getattr(self.args, "cluster_key", "leiden")

            parsed = parse_cluster_keys(cluster_key)
            if len(parsed) > 1 and all(k in adata.obs.columns for k in parsed):
                adata, cluster_key = resolve_cluster_key(adata, cluster_key)

            result = hotspot_pipeline(
                adata,
                layer_key=self.args.raw_count_layer,
                embedding_key=self.args.embedding_hotspot,
                normalization_key=self.args.normalization_key,
                cluster_key=cluster_key,
                skip_hotspot=self.args.skip_hotspot,
                log_dir=log_dir,
            )
            log_step("Controller.Hotspot", "COMPLETED")
            return result
        except Exception as e:
            log_error("Controller.Hotspot", e)
            raise

    def run_step_grn_analysis(self, grn_score_path, grn_links_path):
        """Execute Step 6: GRN Deep Analysis."""
        log_step(
            "Controller.GRNAnalysis", "STARTED", {"grn_score_path": grn_score_path}
        )
        try:
            result = grn_deep_analysis_pipeline(grn_score_path, grn_links_path)
            log_step("Controller.GRNAnalysis", "COMPLETED")
            return result
        except Exception as e:
            log_error("Controller.GRNAnalysis", e)
            raise

    def run_step_report(
        self,
        output_dir,
        adata=None,
        celloracle_result=None,
        hotspot_result=None,
        title="GeneCircuitry Analysis Report",
        subtitle="",
    ):
        """Execute Step 7: Generate HTML/PDF Report."""
        log_step("Controller.Report", "STARTED", {"output_dir": output_dir})
        try:
            log_file = os.path.join(output_dir, "logs", "pipeline.log")
            if not os.path.exists(log_file):
                log_file = None

            cluster_key = getattr(self.args, "cluster_key", "leiden")

            outputs = generate_report(
                output_dir=output_dir,
                title=title,
                subtitle=subtitle,
                adata=adata,
                celloracle_result=celloracle_result,
                hotspot_result=hotspot_result,
                log_file=log_file,
                formats=["html", "pdf"],
                cluster_key=cluster_key,
            )

            log_step(
                "Controller.Report",
                "COMPLETED",
                {"html": outputs.get("html"), "pdf": outputs.get("pdf")},
            )
            return outputs
        except Exception as e:
            log_error("Controller.Report", e)
            log_step(
                "Controller.Report",
                "FAILED",
                {
                    "error": str(e),
                    "error_type": type(e).__name__,
                    "output_dir": output_dir,
                },
            )
            # Don't raise - report generation failure shouldn't stop pipeline
            print(f"  ⚠ Report generation failed ({type(e).__name__}): {e}")
            return {}

    def process_single_stratification(self, adata_cluster, stratification_name):
        """Process a single stratified dataset."""
        log_step(
            "Controller.Stratification.Process",
            "STARTED",
            {"stratification": stratification_name, "n_cells": adata_cluster.n_obs},
        )
        try:
            return self._process_single_stratification_impl(
                adata_cluster, stratification_name
            )
        except Exception as e:
            log_error(f"Controller.Stratification.Process({stratification_name})", e)
            log_step(
                "Controller.Stratification.Process",
                "FAILED",
                {
                    "stratification": stratification_name,
                    "error": str(e),
                    "error_type": type(e).__name__,
                },
            )
            print(
                f"\n⚠ Stratification '{stratification_name}' failed "
                f"({type(e).__name__}): {e}"
            )
            print("  Continuing with remaining stratifications...")
            return None

    def _process_single_stratification_impl(self, adata_cluster, stratification_name):
        """Internal implementation for processing a single stratified dataset."""
        # Create stratified folder INSIDE the main output directory
        stratified_output_dir = os.path.join(
            self.args.output, "stratified_analysis", str(stratification_name)
        )
        stratified_figures_dir = os.path.join(stratified_output_dir, "figures")
        stratified_log_dir = os.path.join(stratified_output_dir, "logs")

        print(f"\n{'='*70}")
        print(f"Processing stratified dataset: {stratification_name}")
        print(f"{'='*70}")
        print(f"  Output directory: {stratified_output_dir}")

        # Create directories for this stratification
        setup_directories(stratified_output_dir, stratified_figures_dir)

        # Update config for this stratification
        config.update_config(
            OUTPUT_DIR=stratified_output_dir,
            FIGURES_DIR=stratified_figures_dir,
            FIGURES_DIR_QC=os.path.join(stratified_figures_dir, "qc"),
            FIGURES_DIR_GRN=os.path.join(stratified_figures_dir, "grn"),
            FIGURES_DIR_HOTSPOT=os.path.join(stratified_figures_dir, "hotspot"),
        )

        # Resolve cluster key on adata_cluster if multi-key
        parsed = parse_cluster_keys(self.args.cluster_key)
        if len(parsed) > 1 and all(k in adata_cluster.obs.columns for k in parsed):
            adata_cluster, resolved_cluster_key = resolve_cluster_key(
                adata_cluster, self.args.cluster_key
            )
        else:
            resolved_cluster_key = resolve_cluster_key_name(self.args.cluster_key)

        # Run pipeline steps
        adata_clustered = dimensionality_reduction_clustering(
            adata_cluster,
            cluster_key=resolved_cluster_key,
            log_dir=stratified_log_dir,
            force=True,
        )

        use_hvgs = getattr(self.args, "use_hvgs", False)

        if use_hvgs:
            # Legacy workflow: CellOracle (HVGs) → Hotspot
            celloracle_result = self.run_step_celloracle(
                adata_clustered,
                log_dir=stratified_log_dir,
                cluster_key=resolved_cluster_key,
            )
            hotspot_result = self.run_step_hotspot(
                adata_clustered,
                log_dir=stratified_log_dir,
                cluster_key=resolved_cluster_key,
            )
        else:
            # Default workflow: Hotspot (identify genes) → CellOracle (use genes)
            hotspot_result = self.run_step_hotspot(
                adata_clustered,
                log_dir=stratified_log_dir,
                cluster_key=resolved_cluster_key,
            )
            hotspot_genes_path = None
            if not self.args.skip_hotspot:
                hotspot_genes_path = os.path.join(
                    config.OUTPUT_DIR, "hotspot", "significant_genes.csv"
                )
            celloracle_result = self.run_step_celloracle(
                adata_clustered,
                log_dir=stratified_log_dir,
                hotspot_genes_path=hotspot_genes_path,
                cluster_key=resolved_cluster_key,
            )

        # Generate summary
        generate_summary(
            adata_clustered,
            celloracle_result,
            hotspot_result,
            self.start_time,
            stratified_output_dir,
            cluster_key=resolved_cluster_key,
        )

        # Run GRN deep analysis
        grn_score_file = os.path.join(
            stratified_output_dir, "celloracle", "grn_merged_scores.csv"
        )
        grn_links_file = os.path.join(
            stratified_output_dir, "celloracle", "grn_filtered_links.pkl"
        )
        if os.path.exists(grn_score_file) and os.path.exists(grn_links_file):
            grn_deep_analysis_pipeline(grn_score_file, grn_links_file)

        # Collect results for unified tabbed report
        self.stratification_results.append(
            {
                "name": str(stratification_name),
                "output_dir": stratified_output_dir,
                "adata": adata_clustered,
                "celloracle_result": celloracle_result,
                "hotspot_result": hotspot_result,
            }
        )

        log_step(
            "Controller.Stratification.Process",
            "COMPLETED",
            {"stratification": stratification_name},
        )
        return stratified_output_dir

    def run_stratified_pipeline_sequential(self):
        """Run stratified pipeline sequentially."""
        log_step(
            "Controller.StratifiedPipeline",
            "STARTED",
            {"n_stratifications": len(self.adata_list)},
        )
        results = []
        failed = []
        for adata_cluster, stratification_name in zip(
            self.adata_list, self.adata_stratification_list
        ):
            result = self.process_single_stratification(
                adata_cluster, stratification_name
            )
            if result is not None:
                results.append(result)
            else:
                failed.append(stratification_name)

        if failed:
            log_step(
                "Controller.StratifiedPipeline",
                "COMPLETED_WITH_FAILURES",
                {
                    "succeeded": len(results),
                    "failed": len(failed),
                    "failed_names": ", ".join(failed),
                },
            )
            print(
                f"\n⚠ {len(failed)}/{len(self.adata_list)} stratifications failed: "
                f"{', '.join(failed)}"
            )
        else:
            log_step(
                "Controller.StratifiedPipeline",
                "COMPLETED",
                {"n_succeeded": len(results)},
            )
        return results

    def run_complete_pipeline(self, steps=None):
        """
        Run complete pipeline or specific steps.

        Parameters:
        -----------
        steps : list of str, optional
            Specific steps to run. Options:
            'load', 'preprocessing', 'stratification', 'clustering',
            'atac_peaks', 'celloracle', 'hotspot', 'grn_analysis',
            'report', 'summary'
            If None, runs all steps.
        """
        if steps is None:
            steps = [
                "load",
                "preprocessing",
                "stratification",
                "clustering",
                "atac_peaks",
                "celloracle",
                "hotspot",
                "grn_analysis",
                "report",
                "summary",
            ]

        # Step 1: Load data
        if "load" in steps:
            self.run_step_load()

        # Step 2: Preprocessing
        if "preprocessing" in steps:
            self.run_step_preprocessing()

        # Step 2.5: Stratification
        if "stratification" in steps:
            self.run_step_stratification()

        # Process ATAC peaks (before CellOracle, applies to all modes)
        if "atac_peaks" in steps and self.args.atac_peaks:
            self.run_step_atac_peaks()

        # Process stratified or non-stratified
        if self.args.cluster_key_stratification and self.adata_list:
            # Stratified analysis (sequential)
            self.run_stratified_pipeline_sequential()
        else:
            # Non-stratified analysis
            if "clustering" in steps:
                adata_clustered = self.run_step_clustering()
            else:
                adata_clustered = self.adata_preprocessed

            celloracle_result = None
            hotspot_result = None

            use_hvgs = getattr(self.args, "use_hvgs", False)

            if use_hvgs:
                # Legacy workflow: CellOracle (HVGs) → Hotspot
                if "celloracle" in steps:
                    celloracle_result = self.run_step_celloracle(adata_clustered)

                if "hotspot" in steps:
                    hotspot_result = self.run_step_hotspot(adata_clustered)
            else:
                # Default workflow: Hotspot (identify genes) → CellOracle (use genes)
                if "hotspot" in steps:
                    hotspot_result = self.run_step_hotspot(adata_clustered)

                hotspot_genes_path = None
                if not self.args.skip_hotspot:
                    hotspot_genes_path = os.path.join(
                        config.OUTPUT_DIR, "hotspot", "significant_genes.csv"
                    )

                if "celloracle" in steps:
                    celloracle_result = self.run_step_celloracle(
                        adata_clustered,
                        hotspot_genes_path=hotspot_genes_path,
                    )

            if "summary" in steps:
                generate_summary(
                    adata_clustered,
                    celloracle_result,
                    hotspot_result,
                    self.start_time,
                    self.args.output,
                    cluster_key=self.args.cluster_key,
                )

            if "grn_analysis" in steps and not self.args.skip_celloracle:
                grn_score_file = os.path.join(
                    self.args.output, "celloracle", "grn_merged_scores.csv"
                )
                grn_links_file = os.path.join(
                    self.args.output, "celloracle", "grn_filtered_links.pkl"
                )
                if os.path.exists(grn_score_file) and os.path.exists(grn_links_file):
                    self.run_step_grn_analysis(grn_score_path=grn_score_file, grn_links_path=grn_links_file)

            # Step: Comparative analysis across clusters
            if "comparative" in steps:
                try:
                    self.run_step_comparative_analysis(
                        adata=adata_clustered,
                        hotspot_obj=hotspot_result,
                        output_dir=self.args.output,
                    )
                except Exception as e:
                    log_error("Controller.ComparativeAnalysis", e)
                    print(f"  ⚠ Comparative analysis failed ({type(e).__name__}): {e}")

            if "report" in steps:
                self.run_step_report(
                    output_dir=self.args.output,
                    adata=adata_clustered,
                    celloracle_result=celloracle_result,
                    hotspot_result=hotspot_result,
                    title="GeneCircuitry Analysis Report",
                    subtitle=self.args.name,
                )

        # Final summary
        if "summary" in steps:
            self.print_final_summary()

    def run_step_comparative_analysis(
        self,
        adata=None,
        score_df=None,
        links_df=None,
        hotspot_obj=None,
        stratification_results=None,
        output_dir=None,
    ):
        """
        Execute comparative analysis and generate aggregation plots.

        Parameters
        ----------
        adata : AnnData, optional
            Input AnnData object.
        score_df : pd.DataFrame, optional
            TF centrality score DataFrame.
        links_df : pd.DataFrame, optional
            CellOracle links DataFrame.
        hotspot_obj : Hotspot, optional
            Hotspot analysis result object.
        stratification_results : list, optional
            List of per-stratification result dicts.
        output_dir : str, optional
            Base output directory.

        Returns
        -------
        dict
            Dictionary of computed comparative analysis results.
        """
        from genecircuitry.comparative_analysis import run_comparative_analysis
        from genecircuitry.plotting.comparative_plots import generate_all_comparative_plots

        if output_dir is None:
            output_dir = getattr(self.args, "output", config.OUTPUT_DIR)

        cluster_key = getattr(self.args, "cluster_key", "leiden")

        # Automatically resolve scores and links if not provided
        if score_df is None:
            co_score_path = os.path.join(output_dir, "celloracle", "grn_merged_scores.csv")
            if os.path.exists(co_score_path):
                import pandas as pd
                score_df = pd.read_csv(co_score_path)

        if links_df is None:
            co_links_path = os.path.join(output_dir, "celloracle", "grn_filtered_links.pkl")
            if os.path.exists(co_links_path):
                import pickle
                with open(co_links_path, "rb") as f:
                    links_df = pickle.load(f)

        comp_results = run_comparative_analysis(
            adata=adata,
            score_df=score_df,
            links_df=links_df,
            hotspot_obj=hotspot_obj,
            stratification_results=stratification_results,
            cluster_key=cluster_key,
            output_dir=output_dir,
            save_tables=True,
        )

        generate_all_comparative_plots(
            comparative_results=comp_results,
            save_name="default",
            skip_existing=False,
        )

        return comp_results

    def print_final_summary(self):
        """Print final pipeline summary."""
        print(f"\n{'='*70}")
        print("Pipeline completed successfully! ✓")
        print(f"{'='*70}")

        if self.adata_stratification_list:
            print(
                f"\nProcessed {len(self.adata_stratification_list)} "
                f"stratified datasets:"
            )
            for stratification in self.adata_stratification_list:
                print(f"  - {os.path.join(self.args.output, str(stratification))}/")
        else:
            print(f"\nResults saved to: {self.args.output}/")

        # Track output files
        track_files(self.args.output)

        # Reset config and generate overall analysis
        config.update_config(
            OUTPUT_DIR=self.args.output,
            FIGURES_DIR=os.path.join(self.args.output, "figures"),
            FIGURES_DIR_QC=os.path.join(self.args.output, "figures", "qc"),
            FIGURES_DIR_GRN=os.path.join(self.args.output, "figures", "grn"),
            FIGURES_DIR_HOTSPOT=os.path.join(self.args.output, "figures", "hotspot"),
            FIGURES_DIR_COMPARATIVE=os.path.join(self.args.output, "figures", "comparative"),
        )

        if self.adata_stratification_list and not self.args.skip_celloracle:
            from genecircuitry.grn_deep_analysis import (
                merge_scores,
                plot_heatmap_scores,
            )

            tracked_files_path = os.path.join(self.args.output, "tracked_files.txt")
            if os.path.exists(tracked_files_path):
                try:
                    total_merged_scores = merge_scores(tracked_files_path)
                    plot_heatmap_scores(total_merged_scores)
                    print("  ✓ Generated overall GRN deep analysis heatmap")
                except Exception as e:
                    log_error("Controller.MergeScores", e)
                    log_step(
                        "Controller.MergeScores",
                        "FAILED",
                        {"error": str(e), "error_type": type(e).__name__},
                    )
                    print(
                        f"  ⚠ GRN score merging/heatmap failed "
                        f"({type(e).__name__}): {e}"
                    )
                    total_merged_scores = None
            else:
                total_merged_scores = None
        else:
            total_merged_scores = None

        # Cross-stratification comparative analysis
        if self.adata_stratification_list and self.stratification_results:
            try:
                self.run_step_comparative_analysis(
                    score_df=total_merged_scores,
                    stratification_results=self.stratification_results,
                    output_dir=self.args.output,
                )
            except Exception as e:
                log_error("Controller.ComparativeAnalysis.Stratified", e)
                print(f"  ⚠ Cross-stratification comparative analysis failed ({type(e).__name__}): {e}")

        # Generate unified report with stratification tabs
        if self.adata_stratification_list and self.stratification_results:
            try:
                log_file = os.path.join(self.args.output, "logs", "pipeline.log")
                if not os.path.exists(log_file):
                    log_file = None

                cluster_key = getattr(self.args, "cluster_key", "leiden")

                outputs = generate_stratified_report(
                    output_dir=self.args.output,
                    title="GeneCircuitry Analysis Report",
                    subtitle=(
                        f"{self.args.name} - Stratified Analysis "
                        f"({len(self.stratification_results)} stratifications)"
                    ),
                    adata_preprocessed=self.adata_preprocessed,
                    stratification_results=self.stratification_results,
                    merged_scores=total_merged_scores,
                    log_file=log_file,
                    formats=["html", "pdf"],
                    cluster_key=cluster_key,
                )
                if outputs.get("html"):
                    print(f"  ✓ Generated unified HTML report: {outputs['html']}")
                if outputs.get("pdf"):
                    print(f"  ✓ Generated unified PDF report: {outputs['pdf']}")
            except Exception as e:
                log_error("Controller.StratifiedReport", e)
                log_step(
                    "Controller.StratifiedReport",
                    "FAILED",
                    {"error": str(e), "error_type": type(e).__name__},
                )
                print(
                    f"  ⚠ Unified report generation failed "
                    f"({type(e).__name__}): {e}"
                )


def setup_directories(output_dir, figures_dir, debug=False):
    """Create necessary output directories."""
    directories = [
        output_dir,
        figures_dir,
        f"{output_dir}/logs",
        f"{output_dir}/celloracle",
        f"{output_dir}/hotspot",
        f"{output_dir}/comparative",
        f"{figures_dir}/qc",
        f"{figures_dir}/grn",
        f"{figures_dir}/hotspot",
        f"{figures_dir}/comparative",
    ]

    if debug:
        for directory in directories:
            shutil.rmtree(directory, ignore_errors=True)
            print(f"✓ Removed directory: {directory}")

    for directory in directories:
        os.makedirs(directory, exist_ok=True)

    print(f"✓ Created output directories")


def load_data(input_path):
    """Load data from file or use example dataset."""
    log_step("Data Loading", "STARTED", {"input_path": input_path})

    print(f"\n{'='*70}")
    print("STEP 1: Data Loading")
    print(f"{'='*70}")

    try:
        if input_path and os.path.exists(input_path):
            print(f"Loading data from: {input_path}")
            log_step("Data Loading", "READING", {"file": input_path})

            if input_path.endswith(".h5ad"):
                adata = sc.read_h5ad(input_path)
            elif input_path.endswith(".h5"):
                adata = sc.read_10x_h5(input_path)
            else:
                raise ValueError(f"Unsupported file format: {input_path}")
        else:
            print("No input file specified or file not found.")
            print(
                "Loading example dataset: PBMC 3k (peripheral blood mononuclear cells)"
            )
            log_step("Data Loading", "USING_EXAMPLE_DATASET")
            adata = sc.datasets.pbmc3k()

        print(f"✓ Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        log_step(
            "Data Loading", "COMPLETED", {"n_obs": adata.n_obs, "n_vars": adata.n_vars}
        )

        return adata

    except Exception as e:
        log_error("Data Loading", e)
        raise


def preprocessing_pipeline(adata, name=None, skip_qc=False, log_dir=None):
    """Run preprocessing pipeline."""
    log_step(
        "Preprocessing",
        "STARTED",
        {"n_obs": adata.n_obs, "n_vars": adata.n_vars, "skip_qc": skip_qc},
    )

    print(f"\n{'='*70}")
    print("STEP 2: Quality Control and Preprocessing")
    print(f"{'='*70}")

    try:
        # Check checkpoint for preprocessing
        step_hash = compute_input_hash(
            None,
            n_obs=adata.n_obs,
            n_vars=adata.n_vars,
            min_genes=config.QC_MIN_GENES,
            min_counts=config.QC_MIN_COUNTS,
            pct_mt_max=config.QC_PCT_MT_MAX,
            skip_qc=skip_qc,
        )

        checkpoint_file = os.path.join(config.OUTPUT_DIR, "preprocessed_adata.h5ad")
        if (
            log_dir
            and check_checkpoint(log_dir, "preprocessing", step_hash)
            and os.path.exists(checkpoint_file)
        ):
            print(f"  Loading preprocessed data from: {checkpoint_file}")
            log_step(
                "Preprocessing",
                "LOADED_FROM_CHECKPOINT",
                {"checkpoint_file": checkpoint_file},
            )
            return sc.read_h5ad(checkpoint_file)

        if not skip_qc:
            print("\n[2.1] Performing quality control...")
            log_step("Preprocessing.QC", "STARTED")

            adata = perform_qc(
                adata,
                min_genes=config.QC_MIN_GENES,
                min_counts=config.QC_MIN_COUNTS,
                pct_counts_mt_max=config.QC_PCT_MT_MAX,
                save_plots=name,
            )
            print(f"  After QC: {adata.n_obs} cells × {adata.n_vars} genes")
            log_step(
                "Preprocessing.QC",
                "COMPLETED",
                {"n_obs": adata.n_obs, "n_vars": adata.n_vars},
            )
        else:
            print("\n[2.1] Skipping QC (already performed)")
            log_step("Preprocessing.QC", "SKIPPED")

        # Normalization
        print("\n[2.2] Normalizing data...")
        log_step("Preprocessing.Normalization", "STARTED")

        adata = perform_normalization(adata)
        print("  ✓ Normalization complete")
        log_step("Preprocessing.Normalization", "COMPLETED")

        # Save checkpoint
        if log_dir:
            adata.write(checkpoint_file)
            write_checkpoint(
                log_dir,
                "preprocessing",
                step_hash,
                n_obs=adata.n_obs,
                n_vars=adata.n_vars,
            )
            log_step("Preprocessing", "CHECKPOINT_SAVED", {"file": checkpoint_file})

        log_step(
            "Preprocessing", "COMPLETED", {"n_obs": adata.n_obs, "n_vars": adata.n_vars}
        )

        return adata

    except Exception as e:
        log_error("Preprocessing", e)
        raise


def stratification_pipeline(adata, cluster_key_stratification=None, clusters="all"):
    """Perform stratification based on specified clustering key(s).

    Supports single keys, comma-separated key strings (e.g. 'key1,key2'),
    or sequences/lists of keys (e.g. ['key1', 'key2']).
    """
    keys = parse_cluster_keys(cluster_key_stratification)
    if not keys:
        print("\nNo stratification performed (no valid clustering key provided).")
        return [], []

    missing = [k for k in keys if k not in adata.obs.columns]
    if missing:
        print(
            f"\nNo stratification performed (stratification key(s) {missing} not found in adata.obs)."
        )
        return [], []

    print(f"\n{'='*70}")
    print("STEP 2.5: Stratification by Clusters")
    print(f"{'='*70}")

    if adata.is_view:
        adata = adata.copy()

    # Ensure all stratification columns are categorical to avoid str/numeric conflicts
    for k in keys:
        if not isinstance(adata.obs[k].dtype, pd.CategoricalDtype):
            adata.obs[k] = adata.obs[k].astype("category")
            print(f"  Converted '{k}' to categorical")

    adata_list = list()
    adata_stratification_list = list()

    if clusters is None or clusters == "all" or clusters == ["all"] or clusters == ("all",):
        requested_clusters = None
    elif isinstance(clusters, (str, bytes)):
        c_str = clusters.decode() if isinstance(clusters, bytes) else clusters
        requested_clusters = set(c.strip() for c in c_str.split(",") if c.strip())
    elif isinstance(clusters, (Iterable, Sequence, np.ndarray, pd.Index, pd.Series)):
        requested_clusters = set()
        for item in clusters:
            if isinstance(item, (str, bytes)):
                item_str = item.decode() if isinstance(item, bytes) else item
                for sub in item_str.split(","):
                    sub_clean = sub.strip()
                    if sub_clean:
                        requested_clusters.add(sub_clean)
            elif item is not None:
                item_clean = str(item).strip()
                if item_clean:
                    requested_clusters.add(item_clean)
        if "all" in requested_clusters and len(requested_clusters) == 1:
            requested_clusters = None
    else:
        c_str = str(clusters).strip()
        requested_clusters = {c_str} if c_str and c_str != "all" else None

    if len(keys) == 1:
        key = keys[0]
        categories = adata.obs[key].cat.categories
        if requested_clusters is not None:
            unique_clusters = []
            for c in categories:
                c_str = str(c)
                c_san = sanitize_identifier(c)
                c_sp = c_str.replace(" ", "_")
                matched = False
                for rc in requested_clusters:
                    rc_str = str(rc).strip()
                    if not rc_str:
                        continue
                    rc_san = sanitize_identifier(rc_str)
                    if (
                        rc_str == c_str
                        or rc_san == c_san
                        or rc_str == c_sp
                        or rc_san == c_str
                    ):
                        matched = True
                        break
                    if ":" in rc_str or "=" in rc_str:
                        sep = ":" if ":" in rc_str else "="
                        target_key, target_val = rc_str.split(sep, 1)
                        if target_key.strip() == key:
                            target_val = target_val.strip()
                            if (
                                target_val == c_str
                                or sanitize_identifier(target_val) == c_san
                                or target_val == c_sp
                            ):
                                matched = True
                                break
                if matched:
                    unique_clusters.append(c)
        else:
            unique_clusters = list(categories)

        for cluster in unique_clusters:
            mask = (adata.obs[key] == cluster).values
            if not mask.any():
                continue
            adata_cluster = adata[mask].copy()
            adata_list.append(adata_cluster)
            # Convert to string for folder naming (handles numeric cluster IDs)
            strat_name = str(cluster).replace(" ", "_")
            adata_stratification_list.append(strat_name)
            print(
                f"  ✓ Cluster '{cluster}': {adata_cluster.n_obs} cells × {adata_cluster.n_vars} genes"
            )

        return adata_list, adata_stratification_list
    else:
        # Multi-key stratification: Cartesian combinations of categories
        import itertools

        key_categories = [list(adata.obs[k].cat.categories) for k in keys]
        combinations = itertools.product(*key_categories)

        for comb in combinations:
            # Composite subgroup identifier joined with underscore
            comp_name = "_".join(sanitize_identifier(v) for v in comb)
            raw_name = "_".join(str(v) for v in comb)

            if requested_clusters is not None:
                comb_strs = [str(v) for v in comb]
                comb_sanitized = [sanitize_identifier(v) for v in comb]
                comb_spaced = [str(v).replace(" ", "_") for v in comb]
                key_val_map = {k: v for k, v in zip(keys, comb)}
                key_val_san = {k: sanitize_identifier(v) for k, v in zip(keys, comb)}

                matched = False
                for rc in requested_clusters:
                    rc_str = str(rc).strip()
                    if not rc_str:
                        continue
                    rc_san = sanitize_identifier(rc_str)

                    # 1. Exact or sanitized composite name match
                    if (
                        rc_str == comp_name
                        or rc_str == raw_name
                        or rc_san == comp_name
                        or rc_san == sanitize_identifier(raw_name)
                    ):
                        matched = True
                        break

                    # 2. Key-value format check, e.g. "col1:clu1" or "col1=clu1"
                    if ":" in rc_str or "=" in rc_str:
                        sep = ":" if ":" in rc_str else "="
                        target_key, target_val = rc_str.split(sep, 1)
                        target_key = target_key.strip()
                        target_val = target_val.strip()
                        if target_key in key_val_map:
                            v = key_val_map[target_key]
                            v_san = key_val_san[target_key]
                            target_val_san = sanitize_identifier(target_val)
                            if (
                                target_val == str(v)
                                or target_val_san == v_san
                                or target_val == str(v).replace(" ", "_")
                                or target_val_san == sanitize_identifier(str(v).replace(" ", "_"))
                            ):
                                matched = True
                                break
                    else:
                        # 3. Individual value match across any key in the combination
                        for v_str, v_san_val, v_sp in zip(
                            comb_strs, comb_sanitized, comb_spaced
                        ):
                            if (
                                rc_str == v_str
                                or rc_san == v_san_val
                                or rc_str == v_sp
                                or rc_san == v_str
                            ):
                                matched = True
                                break
                        if matched:
                            break

                if not matched:
                    continue

            # Find cells matching this combination
            mask = np.ones(adata.n_obs, dtype=bool)
            for k, v in zip(keys, comb):
                mask &= (adata.obs[k] == v).values

            if not mask.any():
                continue  # Skip unobserved combinations

            adata_cluster = adata[mask].copy()
            adata_list.append(adata_cluster)
            adata_stratification_list.append(comp_name)
            print(
                f"  ✓ Subgroup '{comp_name}': {adata_cluster.n_obs} cells × {adata_cluster.n_vars} genes"
            )

        return adata_list, adata_stratification_list


def dimensionality_reduction_clustering(
    adata, cluster_key="leiden", log_dir=None, force=False, **kwargs
):
    """Perform dimensionality reduction and clustering."""
    force = (
        force
        or kwargs.get("force_dim_reduction", False)
        or kwargs.get("force_dimensionality_reduction", False)
    )
    parsed = parse_cluster_keys(cluster_key)
    if len(parsed) > 1 and all(k in adata.obs.columns for k in parsed):
        adata, cluster_key = resolve_cluster_key(adata, cluster_key)
    else:
        cluster_key = resolve_cluster_key_name(cluster_key)

    log_step(
        "DimReduction_Clustering",
        "STARTED",
        {
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "cluster_key": cluster_key,
            "force": force,
        },
    )

    print(f"\n{'='*70}")
    print("STEP 3: Dimensionality Reduction and Clustering")
    print(f"{'='*70}")

    try:
        # Check checkpoint for clustering
        step_hash = compute_input_hash(
            None,
            n_obs=adata.n_obs,
            n_vars=adata.n_vars,
            cluster_key=cluster_key,
            top_genes=config.HVGS_N_TOP_GENES,
            n_neighbors=config.NEIGHBORS_N_NEIGHBORS,
            n_pcs=config.NEIGHBORS_N_PCS,
        )

        checkpoint_file = f"{config.OUTPUT_DIR}/clustered_adata.h5ad"
        if (
            not force
            and log_dir
            and check_checkpoint(log_dir, "clustering", step_hash)
            and os.path.exists(checkpoint_file)
        ):
            log_step(
                "DimReduction_Clustering.Checkpoint",
                "LOADED",
                {"checkpoint_file": checkpoint_file},
            )
            print(f"  Loading clustered data from: {checkpoint_file}")
            adata_loaded = sc.read_h5ad(checkpoint_file)
            # Ensure cluster_key is categorical after loading from checkpoint
            adata_loaded = ensure_categorical_obs(adata_loaded, columns=[cluster_key])
            return adata_loaded

        # Perform dimensionality reduction and clustering
        log_step("DimReduction_Clustering.Processing", "STARTED")
        adata = perform_dimensionality_reduction_clustering(
            adata, cluster_key=cluster_key, force=force
        )
        # Ensure cluster_key is categorical (in addition to default columns)
        adata = ensure_categorical_obs(adata, columns=[cluster_key])
        log_step("DimReduction_Clustering.Processing", "COMPLETED")

        # Get cluster count
        n_clusters = len(adata.obs[cluster_key].unique()) if cluster_key in adata.obs.columns else 0
        print(f"✓ Identified {n_clusters} clusters")

        # Save checkpoint
        if log_dir:
            log_step("DimReduction_Clustering.Checkpoint", "SAVING")
            adata.write(checkpoint_file)
            write_checkpoint(
                log_dir,
                "clustering",
                step_hash,
                n_obs=adata.n_obs,
                n_vars=adata.n_vars,
                n_clusters=n_clusters,
            )
            log_step(
                "DimReduction_Clustering.Checkpoint",
                "SAVED",
                {"checkpoint_file": checkpoint_file},
            )

        log_step(
            "DimReduction_Clustering",
            "COMPLETED",
            {"n_obs": adata.n_obs, "n_vars": adata.n_vars, "n_clusters": n_clusters},
        )
        return adata
    except Exception as e:
        log_error("DimReduction_Clustering", e)
        raise


def celloracle_pipeline(
    adata,
    cluster_key="leiden",
    species="human",
    embedding_name="X_draw_graph_fa",
    raw_count_layer="raw_counts",
    TG_to_TF_dictionary=None,
    skip_celloracle=False,
    no_base_grn=False,
    hotspot_genes_path=None,
    log_dir=None,
    **kwargs,
):
    """Run CellOracle GRN inference pipeline.

    Parameters
    ----------
    hotspot_genes_path : str, optional
        Path to ``significant_genes.csv`` produced by Hotspot.  When provided,
        the genes listed in that file are used for GRN preprocessing instead of
        standard highly-variable-gene selection.  This is the default workflow
        (Hotspot-first); pass ``None`` to use HVGs (legacy ``--use-hvgs`` mode).
    """
    cluster_key = kwargs.get(
        "cluster_column_name", kwargs.get("cluster_column", cluster_key)
    )
    parsed = parse_cluster_keys(cluster_key)
    if len(parsed) > 1 and all(k in adata.obs.columns for k in parsed):
        adata, cluster_key = resolve_cluster_key(adata, cluster_key)
    else:
        cluster_key = resolve_cluster_key_name(cluster_key)

    log_step(
        "CellOracle",
        "STARTED",
        {"n_obs": adata.n_obs, "cluster_key": cluster_key, "species": species},
    )

    try:
        print(f"\n{'='*70}")
        print("STEP 4: CellOracle GRN Inference")
        print(f"{'='*70}")

        if skip_celloracle:
            log_step("CellOracle", "SKIPPED", {"reason": "--skip-celloracle"})
            print("⊘ Skipping CellOracle analysis (--skip-celloracle)")
            return None

        # Check checkpoint for CellOracle
        step_hash = compute_input_hash(
            None,
            n_obs=adata.n_obs,
            cluster_key=cluster_key,
            species=species,
            embedding_name=embedding_name,
            hotspot_genes_path=hotspot_genes_path,
        )

        oracle_file = f"{config.OUTPUT_DIR}/celloracle/oracle_object.celloracle.oracle"
        links_file = f"{config.OUTPUT_DIR}/celloracle/oracle_object.celloracle.links"

        if log_dir and check_checkpoint(log_dir, "celloracle", step_hash):
            if os.path.exists(oracle_file) and os.path.exists(links_file):
                try:
                    from genecircuitry.celloracle_processing import (
                        load_celloracle_results,
                    )

                    log_step("CellOracle.Checkpoint", "LOADING")
                    print(f"  Loading CellOracle results from checkpoint...")
                    oracle, links = load_celloracle_results(
                        oracle_path=oracle_file, links_path=links_file
                    )
                    log_step("CellOracle.Checkpoint", "LOADED")
                    return oracle, links
                except Exception as e:
                    log_step("CellOracle.Checkpoint", "FAILED", {"error": str(e)})
                    print(f"  ⚠ Error loading checkpoint: {e}")
                    print(f"  Re-running CellOracle analysis...")

        try:
            from genecircuitry.celloracle_processing import (
                perform_grn_pre_processing,
                create_oracle_object,
                run_PCA,
                run_KNN,
                run_links,
                save_celloracle_results,
                load_hotspot_genes,
            )

            # Load Hotspot genes if path is provided (Hotspot-first workflow)
            gene_list = None
            if hotspot_genes_path is not None:
                if os.path.exists(hotspot_genes_path):
                    log_step(
                        "CellOracle.LoadHotspotGenes",
                        "STARTED",
                        {"path": hotspot_genes_path},
                    )
                    gene_list = load_hotspot_genes(hotspot_genes_path)
                    log_step(
                        "CellOracle.LoadHotspotGenes",
                        "COMPLETED",
                        {"n_genes": len(gene_list)},
                    )
                    print(
                        f"  Using {len(gene_list)} Hotspot autocorrelated genes "
                        f"for GRN preprocessing"
                    )
                else:
                    log_step(
                        "CellOracle.LoadHotspotGenes",
                        "SKIPPED",
                        {
                            "reason": "file not found",
                            "path": hotspot_genes_path,
                        },
                    )
                    print(
                        f"  ⚠ Hotspot genes file not found: {hotspot_genes_path}"
                    )
                    print(
                        "    Falling back to HVG selection for CellOracle preprocessing."
                    )

            log_step("CellOracle.Preprocessing", "STARTED")
            print("\n[4.1] Preprocessing data for CellOracle...")
            adata = perform_grn_pre_processing(
                adata,
                cluster_key=cluster_key,
                gene_list=gene_list,
            )
            print("  ✓ Preprocessing complete")

            log_step("CellOracle.CreateObject", "STARTED")
            print("\n[4.2] Creating Oracle object...")
            oracle = create_oracle_object(
                adata,
                cluster_column_name=cluster_key,
                embedding_name=embedding_name,
                raw_count_layer=raw_count_layer,
                species=species,
                TG_to_TF_dictionary=TG_to_TF_dictionary,
                no_base_grn=no_base_grn,
            )
            print("  ✓ Oracle object created")
            log_step("CellOracle.CreateObject", "COMPLETED")

            log_step("CellOracle.PCA", "STARTED")
            print("\n[4.3] Running PCA on Oracle object...")
            oracle, n_comps = run_PCA(oracle)
            print("  ✓ PCA complete")
            log_step("CellOracle.PCA", "COMPLETED", {"n_comps": n_comps})

            log_step("CellOracle.KNN", "STARTED")
            print("\n[4.4] Running KNN imputation...")
            oracle = run_KNN(oracle, n_comps=n_comps)
            print("  ✓ KNN imputation complete")
            log_step("CellOracle.KNN", "COMPLETED")

            log_step("CellOracle.InferGRN", "STARTED")
            print("\n[4.5] Inferring GRN links...")
            links = run_links(
                oracle,
                cluster_column_name=cluster_key,
                p_cutoff=0.001,
            )
            print("  ✓ GRN inference complete")
            log_step("CellOracle.InferGRN", "COMPLETED")

            log_step("CellOracle.SaveResults", "STARTED")
            print("\n[4.6] Saving CellOracle results...")
            save_celloracle_results(oracle, links)
            print("  ✓ CellOracle results saved")
            log_step("CellOracle.SaveResults", "COMPLETED")

            # Save checkpoint
            if log_dir:
                write_checkpoint(
                    log_dir,
                    "celloracle",
                    step_hash,
                    n_obs=adata.n_obs,
                    cluster_key=cluster_key,
                )

            log_step("CellOracle", "COMPLETED")
            return oracle, links

        except ImportError as ie:
            log_error("CellOracle.Import", ie)
            log_step(
                "CellOracle",
                "SKIPPED",
                {
                    "reason": "CellOracle import failed",
                    "missing_module": getattr(ie, "name", None),
                    "error": str(ie),
                },
            )
            print("\n⚠ CellOracle could not be loaded — GRN inference skipped.")
            print(f"  Missing module : {ie.name!r}")
            print(f"  ImportError    : {ie}")
            print("  To fix         : pip install celloracle")
            print(
                "  Note: if celloracle is already installed, one of its dependencies\n"
                "        may be missing or incompatible. Run:\n"
                "          python -c 'import celloracle'\n"
                "        to reproduce the error outside of this pipeline."
            )
            return None
    except Exception as e:
        log_error("CellOracle", e)
        log_step(
            "CellOracle", "FAILED", {"error": str(e), "error_type": type(e).__name__}
        )
        print(f"\n⚠ CellOracle analysis failed ({type(e).__name__}): {e}")
        print("  Continuing with remaining analysis...")
        return None


def hotspot_pipeline(
    adata,
    layer_key="raw_counts",
    embedding_key="X_pca",
    normalization_key="n_counts",
    cluster_key="leiden",
    skip_hotspot=False,
    log_dir=None,
    **kwargs,
):
    """Run Hotspot gene module identification pipeline."""
    cluster_key = kwargs.get(
        "cluster_column_name", kwargs.get("cluster_column", cluster_key)
    )
    parsed = parse_cluster_keys(cluster_key)
    if len(parsed) > 1 and all(k in adata.obs.columns for k in parsed):
        adata, cluster_key = resolve_cluster_key(adata, cluster_key)
    else:
        cluster_key = resolve_cluster_key_name(cluster_key)

    log_step(
        "Hotspot", "STARTED", {"n_obs": adata.n_obs, "embedding_key": embedding_key}
    )

    try:
        print(f"\n{'='*70}")
        print("STEP 5: Hotspot Gene Module Identification")
        print(f"{'='*70}")

        if skip_hotspot:
            log_step("Hotspot", "SKIPPED", {"reason": "--skip-hotspot"})
            print("⊘ Skipping Hotspot analysis (--skip-hotspot)")
            return None

        # Check checkpoint for Hotspot
        step_hash = compute_input_hash(
            None,
            n_obs=adata.n_obs,
            top_genes=config.HOTSPOT_TOP_GENES,
            embedding_key=embedding_key,
            fdr_threshold=config.HOTSPOT_FDR_THRESHOLD,
        )

        hotspot_file = f"{config.OUTPUT_DIR}/hotspot/hotspot_object.pkl"
        if (
            log_dir
            and check_checkpoint(log_dir, "hotspot", step_hash)
            and os.path.exists(hotspot_file)
        ):
            log_step("Hotspot.Checkpoint", "LOADING")
            print(f"  Loading Hotspot results from checkpoint...")
            with open(hotspot_file, "rb") as f:
                hotspot_obj = pickle.load(f)
            log_step("Hotspot.Checkpoint", "LOADED")
            return hotspot_obj

        try:
            from genecircuitry.hotspot_processing import (
                create_hotspot_object,
                run_hotspot_analysis,
            )

            log_step("Hotspot.CreateObject", "STARTED")
            print("\n[5.1] Creating Hotspot object...")
            hotspot_obj = create_hotspot_object(
                adata,
                top_genes=config.HOTSPOT_TOP_GENES,
                layer_key=layer_key,
                model="danb",
                embedding_key=embedding_key,
                normalization_key=normalization_key,
            )
            print("  ✓ Hotspot object created")
            log_step(
                "Hotspot.CreateObject",
                "COMPLETED",
                {"top_genes": config.HOTSPOT_TOP_GENES},
            )

            log_step("Hotspot.Analysis", "STARTED")
            print("\n[5.2] Running Hotspot analysis...")
            print("  (This may take several minutes...)")
            hotspot_obj = run_hotspot_analysis(
                hotspot_obj, adata=adata, cluster_key=cluster_key
            )
            print("  ✓ Hotspot analysis complete")

            # Get results summary
            autocorr_results = hotspot_obj.results
            significant_genes = autocorr_results[
                autocorr_results.FDR < config.HOTSPOT_FDR_THRESHOLD
            ]

            print(f"\n  Analysis Summary:")
            print(f"    Total genes analyzed: {len(autocorr_results)}")
            print(f"    Significant genes: {len(significant_genes)}")

            n_modules = 0
            if hasattr(hotspot_obj, "modules"):
                n_modules = len(hotspot_obj.modules.unique())
                print(f"    Gene modules identified: {n_modules}")

            log_step(
                "Hotspot.Analysis",
                "COMPLETED",
                {
                    "total_genes": len(autocorr_results),
                    "significant_genes": len(significant_genes),
                    "n_modules": n_modules,
                },
            )

            # Save checkpoint
            if log_dir:
                write_checkpoint(
                    log_dir,
                    "hotspot",
                    step_hash,
                    n_genes=len(autocorr_results),
                    n_significant=len(significant_genes),
                )

            log_step(
                "Hotspot",
                "COMPLETED",
                {
                    "total_genes": len(autocorr_results),
                    "significant_genes": len(significant_genes),
                },
            )
            return hotspot_obj

        except ImportError as ie:
            log_error("Hotspot.Import", ie)
            log_step(
                "Hotspot",
                "SKIPPED",
                {"reason": "Hotspot not installed", "error": str(ie)},
            )
            print("\n⚠ Hotspot not installed. Skipping module identification.")
            print("  To install: pip install hotspotsc")
            return None
    except Exception as e:
        log_error("Hotspot", e)
        log_step("Hotspot", "FAILED", {"error": str(e), "error_type": type(e).__name__})
        print(f"\n⚠ Hotspot analysis failed ({type(e).__name__}): {e}")
        print("  Continuing with remaining analysis...")
        return None


def generate_summary(
    adata,
    celloracle_result,
    hotspot_result,
    start_time,
    output_dir,
    cluster_key="leiden",
    **kwargs,
):
    """Generate analysis summary report."""
    cluster_key = kwargs.get(
        "cluster_column_name", kwargs.get("cluster_column", cluster_key)
    )
    cluster_key = resolve_cluster_key_name(cluster_key)
    log_step("GenerateSummary", "STARTED")

    try:
        print(f"\n{'='*70}")
        print("STEP 6: Generating Analysis Summary")
        print(f"{'='*70}")

        end_time = datetime.now()
        duration = end_time - start_time

        summary = []
        summary.append("=" * 70)
        summary.append("GeneCircuitry Complete Pipeline - Analysis Summary")
        summary.append("=" * 70)
        summary.append(
            f"\nAnalysis completed at: {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
        )
        summary.append(f"Total runtime: {duration}")
        summary.append(f"\n{'='*70}")
        summary.append("Dataset Information")
        summary.append("=" * 70)
        summary.append(f"Final dataset: {adata.n_obs} cells × {adata.n_vars} genes")

        if cluster_key in adata.obs.columns:
            n_clusters = len(adata.obs[cluster_key].unique())
            summary.append(f"Clusters identified ({cluster_key}): {n_clusters}")
        elif "leiden" in adata.obs.columns:
            n_clusters = len(adata.obs["leiden"].unique())
            summary.append(f"Clusters identified (leiden): {n_clusters}")
        elif "louvain" in adata.obs.columns:
            n_clusters = len(adata.obs["louvain"].unique())
            summary.append(f"Clusters identified (louvain): {n_clusters}")

        # CellOracle results
        summary.append(f"\n{'='*70}")
        summary.append("CellOracle GRN Inference")
        summary.append("=" * 70)
        if celloracle_result:
            oracle, links = celloracle_result
            summary.append(f"Status: ✓ Completed")
            summary.append(
                f"Oracle object: {output_dir}/celloracle/oracle_object.celloracle.oracle"
            )
            summary.append(
                f"GRN links: {output_dir}/celloracle/grn_links.celloracle.links"
            )
        else:
            summary.append(f"Status: ⊘ Skipped or failed")

        # Hotspot results
        summary.append(f"\n{'='*70}")
        summary.append("Hotspot Module Identification")
        summary.append("=" * 70)
        if hotspot_result:
            autocorr_results = hotspot_result.results
            significant = autocorr_results[
                autocorr_results.FDR < config.HOTSPOT_FDR_THRESHOLD
            ]
            summary.append(f"Status: ✓ Completed")
            summary.append(f"Genes analyzed: {len(autocorr_results)}")
            summary.append(f"Significant genes: {len(significant)}")
            if hasattr(hotspot_result, "modules"):
                n_modules = len(hotspot_result.modules.unique())
                summary.append(f"Modules identified: {n_modules}")
            summary.append(f"Results: {output_dir}/hotspot/")
        else:
            summary.append(f"Status: ⊘ Skipped or failed")

        # Output files
        summary.append(f"\n{'='*70}")
        summary.append("Output Files")
        summary.append("=" * 70)
        summary.append(f"Preprocessed data: {output_dir}/preprocessed_adata.h5ad")
        summary.append(f"Figures directory: {output_dir}/figures/")
        summary.append(f"Results directory: {output_dir}/")

        # Save summary to file
        summary_text = "\n".join(summary)
        summary_path = f"{output_dir}/analysis_summary.txt"
        with open(summary_path, "w") as f:
            f.write(summary_text)

        print(summary_text)
        print(f"\n✓ Summary saved to: {summary_path}")
        log_step("GenerateSummary", "COMPLETED", {"summary_file": summary_path})

    except Exception as e:
        log_error("GenerateSummary", e)
        log_step(
            "GenerateSummary", "FAILED", {"error": str(e), "error_type": type(e).__name__}
        )
        print(f"⚠ Summary generation failed ({type(e).__name__}): {e}")


def track_files(output_dir):
    """Track output files generated during the pipeline."""
    tracked_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            tracked_files.append(os.path.join(root, file))

    tracked_files_path = os.path.join(output_dir, "tracked_files.txt")
    if os.path.exists(tracked_files_path):
        os.remove(tracked_files_path)
    with open(tracked_files_path, "w") as f:
        for file in tracked_files:
            f.write(f"{file}\n")
    print(f"\n✓ Tracked output files saved to: {tracked_files_path}")

    return tracked_files


def grn_deep_analysis_pipeline(grn_score_path, grn_links_path):
    """Run GRN deep analysis using tracked output files."""
    log_step(
        "GRNDeepAnalysis",
        "STARTED",
        {"score_path": grn_score_path, "links_path": grn_links_path},
    )

    try:
        from genecircuitry.grn_deep_analysis import (
            process_single_score_file,
            process_single_links_file,
            plot_scatter_scores,
            plot_compare_cluster_scores,
            plot_difference_cluster_scores,
            plot_network_graph,
        )
        from genecircuitry.plotting.grn_plots import (
            generate_all_grn_plots,
            plot_enriched_tf_network,
            plot_tf_shared_target_network,
        )

        print(f"\n{'='*70}")
        print("STEP 7: GRN Deep Analysis")
        print(f"{'='*70}")

        os.makedirs(f"{config.OUTPUT_DIR}/grn_deep_analysis", exist_ok=True)
        os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

        score_df = process_single_score_file(grn_score_path)
        print("  ✓ Processed GRN score file")
        links_df = process_single_links_file(grn_links_path)
        print("  ✓ Processed GRN links file")

        # --- Legacy individual plots ---
        plot_compare_cluster_scores(score_df)
        print("  ✓ Generated cluster comparison plots")
        plot_difference_cluster_scores(score_df)
        print("  ✓ Generated cluster difference plots")
        plot_scatter_scores(score_df)
        print("  ✓ Generated scatter plots")
        plot_network_graph(score_df, links_df)
        print("  ✓ Generated GRN network graph")

        # --- New plotting module: enriched TF network ---
        log_step("GRNPlotting", "STARTED")
        try:
            results = generate_all_grn_plots(
                score_df=score_df,
                links_df=links_df,
                skip_existing=True,
            )
            total = sum(results.values())
            log_step(
                "GRNPlotting",
                "COMPLETED",
                {"plots_generated": total, "details": results},
            )
            print(f"  ✓ Generated {total} GRN plots via plotting module")
        except Exception as e:
            log_error("GRNPlotting", e)
            log_step(
                "GRNPlotting",
                "FAILED",
                {"error": str(e), "error_type": type(e).__name__},
            )
            print(f"  ⚠ GRN plotting module failed ({type(e).__name__}): {e}")

        log_step("GRNDeepAnalysis", "COMPLETED")

    except Exception as e:
        log_error("GRNDeepAnalysis", e)
        log_step(
            "GRNDeepAnalysis",
            "FAILED",
            {"error": str(e), "error_type": type(e).__name__},
        )
        print(f"\n⚠ GRN deep analysis failed ({type(e).__name__}): {e}")
        print("  Continuing with remaining pipeline steps...")


def main():
    """Main pipeline execution."""
    parser = argparse.ArgumentParser(
        description="Run complete GeneCircuitry analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with example dataset
  python -m genecircuitry.pipeline
  
  # Run with custom input
  python -m genecircuitry.pipeline --input data/my_data.h5ad
  
  # Skip specific analyses
  python -m genecircuitry.pipeline --skip-celloracle --skip-hotspot
  
  # Run stratified analysis
  python -m genecircuitry.pipeline --cluster-key-stratification celltype
  
  # Run specific steps only
  python -m genecircuitry.pipeline --steps load preprocessing clustering
  
  # Custom configuration
  python -m genecircuitry.pipeline --seed 123 --n-jobs 16
        """,
    )

    # Input/Output arguments
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default=None,
        help="Input data file (.h5ad). If not provided, uses example dataset.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser.add_argument(
        "--name",
        "-n",
        type=str,
        default="test_run",
        help="Name of the analysis run (optional)",
    )
    parser.add_argument(
        "--config",
        "-c",
        type=str,
        default=None,
        help="Path to a YAML or JSON config file with parameter overrides "
        "(e.g. QC_MIN_GENES, LEIDEN_RESOLUTION). "
        "Explicit CLI arguments take precedence over values in this file.",
    )

    # Analysis parameters
    parser.add_argument(
        "--species",
        "-s",
        type=str,
        default="human",
        help="Species for GRN inference (default: human)",
    )
    parser.add_argument(
        "--cluster-key",
        "--cluster-column",
        "--cluster-column-name",
        type=str,
        default="leiden",
        dest="cluster_key",
        help="Clustering column name in adata.obs (default: leiden)",
    )
    parser.add_argument(
        "--clusters",
        type=str,
        default="all",
        help="Specific clusters to analyze (comma-separated, default: all)",
    )
    parser.add_argument(
        "--cluster-key-stratification",
        type=str,
        default=None,
        help="Clustering column name to perform stratification on (default: None)",
    )
    parser.add_argument(
        "--embedding-grn",
        type=str,
        default="X_draw_graph_fa",
        help="Embedding name for GRN inference (default: X_draw_graph_fa)",
    )
    parser.add_argument(
        "--embedding-hotspot",
        type=str,
        default="X_umap",
        help="Embedding name for Hotspot analysis (default: X_umap)",
    )
    parser.add_argument(
        "--normalization-key",
        type=str,
        default="n_counts",
        help="Column name in adata.obs for normalization counts (default: n_counts)",
    )
    parser.add_argument(
        "--raw-count-layer",
        type=str,
        default="raw_counts",
        help="Layer name for raw counts (default: raw_counts)",
    )
    parser.add_argument(
        "--tf-dictionary",
        type=str,
        default=None,
        help="Path to TF to target gene dictionary pickle file (default: None)",
    )
    parser.add_argument(
        "--atac-peaks",
        type=str,
        default=None,
        help="Path to BED file with pre-called ATAC peaks. "
        "When provided, peaks are processed through CellOracle motif analysis "
        "to generate an enriched TF info matrix used as custom base GRN.",
    )
    parser.add_argument(
        "--no-base-grn",
        action="store_true",
        default=False,
        help="Do not use a base GRN (use if already provided)",
    )
    parser.add_argument(
        "--cell-downsample",
        "--grn-cell-downsample",
        type=int,
        default=config.GRN_CELL_DOWNSAMPLE,
        dest="cell_downsample",
        help=f"Number of cells to downsample to for GRN analysis (default: {config.GRN_CELL_DOWNSAMPLE})",
    )

    # Quality control parameters
    parser.add_argument(
        "--min-genes",
        type=int,
        default=config.QC_MIN_GENES,
        help=f"Minimum genes per cell for QC (default: {config.QC_MIN_GENES})",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=config.QC_MIN_COUNTS,
        help=f"Minimum counts per cell for QC (default: {config.QC_MIN_COUNTS})",
    )

    # Computational parameters
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=config.N_JOBS,
        help=f"Number of parallel jobs (default: {config.N_JOBS})",
    )

    # Pipeline control flags
    parser.add_argument(
        "--skip-qc",
        action="store_true",
        default=False,
        help="Skip quality control (use if already performed)",
    )
    parser.add_argument(
        "--force-dim-reduction",
        "--force-dimensionality-reduction",
        action="store_true",
        default=False,
        dest="force_dim_reduction",
        help="Force re-run of dimensionality reduction and clustering even if results already exist",
    )
    parser.add_argument(
        "--skip-celloracle",
        action="store_true",
        default=False,
        help="Skip CellOracle GRN inference",
    )
    parser.add_argument(
        "--skip-hotspot",
        action="store_true",
        default=False,
        help="Skip Hotspot module identification",
    )
    parser.add_argument(
        "--use-hvgs",
        action="store_true",
        default=False,
        help=(
            "Use highly variable genes (HVGs) for CellOracle instead of "
            "Hotspot autocorrelated genes. Reverts to the legacy workflow "
            "(CellOracle first, then Hotspot). By default, Hotspot runs first "
            "and its significant genes are used as input to CellOracle."
        ),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enable debug mode",
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        default=None,
        help="Specific pipeline steps to run (space-separated): "
        "load preprocessing stratification clustering atac_peaks "
        "celloracle hotspot grn_analysis summary",
    )

    args = parser.parse_args()

    # Capture which args were explicitly provided by the user on the command line.
    # We do this by parsing a second time with all defaults set to a sentinel so
    # that only arguments that actually appeared on sys.argv get a non-sentinel
    # value — avoiding the false-negative when a user passes a value that equals
    # the default (e.g. --seed 42 when the default is already 42).
    _sentinel = object()
    _sentinel_parser = argparse.ArgumentParser(add_help=False)
    for action in parser._actions:
        if action.dest == "help" or not action.option_strings:
            continue
        _sentinel_parser.add_argument(
            *action.option_strings,
            dest=action.dest,
            nargs=action.nargs if action.nargs else (None if action.const is None else "?"),
            default=_sentinel,
            const=action.const,
        )
    _sentinel_args, _ = _sentinel_parser.parse_known_args()
    _explicit_args = {
        k for k, v in vars(_sentinel_args).items() if v is not _sentinel
    }

    # Print header
    print("\n" + "=" * 70)
    print("GeneCircuitry Complete Analysis Pipeline")
    print("=" * 70)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    start_time = datetime.now()

    # Setup configuration
    print(f"\n{'='*70}")
    print("STEP 0: Configuration Setup")
    print(f"{'='*70}")

    # Load config file BEFORE directory setup so all subsequent calls see updated values.
    # CLI arguments that were explicitly provided will override the file below.
    if args.config:
        try:
            loaded = config.load_config_file(args.config)
            print(f"  Loaded {len(loaded)} parameter(s) from config file: {args.config}")
        except (FileNotFoundError, ValueError, ImportError) as e:
            print(f"\nError loading config file: {e}")
            return 1

    # Create directories
    setup_directories(
        output_dir=args.output,
        figures_dir=os.path.join(args.output, "figures"),
        debug=args.debug,
    )

    # Resolve seed now (after config file is loaded) so logging and setup both use
    # the same value: explicit CLI flag wins, otherwise fall back to config file value.
    args.seed = args.seed if "seed" in _explicit_args else config.RANDOM_SEED

    # Initialize logging system
    setup_logging(args.output)
    log_step(
        "Pipeline Initialization",
        "STARTED",
        {"output_dir": args.output, "random_seed": args.seed, "n_jobs": args.n_jobs},
    )

    set_random_seed(args.seed)
    set_scanpy_settings()

    sc.settings.logfile = os.path.join(args.output, "logs", "scanpy_log.txt")

    # Update GeneCircuitry configuration.
    # Output/figures dirs always come from CLI.  QC and compute params only
    # override the config file when the user explicitly passed the CLI flag.
    _config_overrides = dict(
        OUTPUT_DIR=args.output,
        FIGURES_DIR=os.path.join(args.output, "figures"),
        FIGURES_DIR_QC=os.path.join(args.output, "figures", "qc"),
        FIGURES_DIR_GRN=os.path.join(args.output, "figures", "grn"),
        FIGURES_DIR_HOTSPOT=os.path.join(args.output, "figures", "hotspot"),
    )
    if "n_jobs" in _explicit_args:
        _config_overrides.update(
            GRN_N_JOBS=args.n_jobs,
            HOTSPOT_N_JOBS=args.n_jobs,
            N_JOBS=args.n_jobs,
        )
    if "min_genes" in _explicit_args:
        _config_overrides["QC_MIN_GENES"] = args.min_genes
    if "min_counts" in _explicit_args:
        _config_overrides["QC_MIN_COUNTS"] = args.min_counts
    if "force_dim_reduction" in _explicit_args:
        _config_overrides["FORCE_DIM_REDUCTION"] = args.force_dim_reduction
    if "cell_downsample" in _explicit_args:
        _config_overrides["GRN_CELL_DOWNSAMPLE"] = args.cell_downsample
    config.update_config(**_config_overrides)

    print(f"Random seed: {args.seed}")
    print(f"Output directory: {args.output}")
    print(f"Figures directory: {os.path.join(args.output, 'figures')}")
    print(f"Parallel jobs: {args.n_jobs}")

    # Print active workflow mode
    if args.use_hvgs:
        print(
            "\nWorkflow mode: LEGACY (--use-hvgs)"
        )
        print(
            "  Order: Load → Preprocess → Cluster → CellOracle (HVGs) → Hotspot"
        )
        log_step("Pipeline", "WORKFLOW_MODE", {"mode": "legacy_hvgs"})
    else:
        print("\nWorkflow mode: DEFAULT (Hotspot-first)")
        print(
            "  Order: Load → Preprocess → Cluster → Hotspot → CellOracle (Hotspot genes)"
        )
        print("  Use --use-hvgs to revert to the legacy CellOracle-first workflow.")
        log_step("Pipeline", "WORKFLOW_MODE", {"mode": "hotspot_first"})

    # Print non-default input arguments
    non_default_args = []
    for arg, value in vars(args).items():
        default_value = parser.get_default(arg)
        if value != default_value:
            non_default_args.append(f"  {arg}: {value}")

    if non_default_args:
        print("\nNon-default arguments:")
        for arg_info in non_default_args:
            print(arg_info)
        log_step(
            "Pipeline Initialization",
            "NON_DEFAULT_ARGS",
            {
                "args": ", ".join(
                    [
                        f"{k}={v}"
                        for k, v in [
                            (arg.split(": ")[0].strip(), arg.split(": ")[1])
                            for arg in non_default_args
                        ]
                    ]
                )
            },
        )
    else:
        print("\nAll arguments using default values")

    log_step("Pipeline Initialization", "COMPLETED")

    try:
        # Create pipeline controller
        log_step("Controller Creation", "STARTED")
        controller = PipelineController(args, start_time)
        log_step("Controller Creation", "COMPLETED")

        # Run pipeline
        log_step(
            "Pipeline Execution",
            "STARTED",
            {"steps": args.steps if args.steps else "all"},
        )

        controller.run_complete_pipeline(steps=args.steps)

        log_step("Pipeline Execution", "COMPLETED")
        log_step(
            "PIPELINE", "SUCCESS", {"total_duration": str(datetime.now() - start_time)}
        )

        return 0

    except Exception as e:
        print(f"\n{'='*70}")
        print(f"Pipeline failed with error: {e}")
        print(f"{'='*70}\n")

        # Log the error
        log_error("Pipeline Execution", e)
        log_step(
            "PIPELINE",
            "FAILED",
            {"error": str(e), "duration": str(datetime.now() - start_time)},
        )

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
