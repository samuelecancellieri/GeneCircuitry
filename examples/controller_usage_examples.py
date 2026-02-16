#!/usr/bin/env python
"""
Examples of using the PipelineController for flexible pipeline execution
=========================================================================

This script demonstrates various ways to use the PipelineController
to run specific parts of the pipeline or customize execution.
"""

import sys
import os
from argparse import Namespace
from datetime import datetime

# Import the controller from its proper package location
from trnspot.pipeline import PipelineController, setup_directories
from trnspot import set_random_seed, set_scanpy_settings, config


def example_1_run_specific_steps():
    """Example 1: Run only specific pipeline steps"""
    print("\n" + "=" * 70)
    print("Example 1: Running only preprocessing and clustering")
    print("=" * 70)

    # Create arguments
    args = Namespace(
        input="data/paul15/paul15.h5",
        output="example_output_1",
        name="partial_run",
        species="human",
        cluster_key="leiden",
        cluster_key_stratification=None,
        embedding_grn="X_draw_graph_fa",
        embedding_hotspot="X_umap",
        raw_count_layer="raw_counts",
        tf_dictionary=None,
        min_genes=200,
        min_counts=500,
        seed=42,
        n_jobs=8,
        skip_qc=False,
        skip_celloracle=False,
        skip_hotspot=False,
        debug=False,
        parallel=False,
        steps=None,
    )

    # Setup
    setup_directories(args.output, os.path.join(args.output, "figures"))
    set_random_seed(args.seed)
    set_scanpy_settings()
    config.update_config(OUTPUT_DIR=args.output)

    # Create controller
    controller = PipelineController(args, datetime.now())

    # Run only specific steps
    controller.run_complete_pipeline(
        steps=["load", "preprocessing", "clustering"], parallel=False
    )

    print("\n✓ Completed preprocessing and clustering only")
    print(f"  Preprocessed data saved to: {args.output}/preprocessed_adata.h5ad")
    print(f"  Clustered data saved to: {args.output}/clustered_adata.h5ad")


def example_2_stratified_parallel():
    """Example 2: Run stratified analysis in parallel"""
    print("\n" + "=" * 70)
    print("Example 2: Stratified analysis with parallel execution")
    print("=" * 70)

    args = Namespace(
        input="data/paul15/paul15.h5",
        output="example_output_2",
        name="parallel_run",
        species="human",
        cluster_key="leiden",
        cluster_key_stratification="paul15_clusters",  # Use existing clustering
        embedding_grn="X_draw_graph_fa",
        embedding_hotspot="X_umap",
        raw_count_layer="raw_counts",
        tf_dictionary=None,
        min_genes=200,
        min_counts=500,
        seed=42,
        n_jobs=4,  # Number of parallel workers
        skip_qc=False,
        skip_celloracle=True,  # Skip to speed up example
        skip_hotspot=False,
        debug=False,
        parallel=True,  # Enable parallel processing
        steps=None,
    )

    # Setup
    setup_directories(args.output, os.path.join(args.output, "figures"))
    set_random_seed(args.seed)
    set_scanpy_settings()
    config.update_config(OUTPUT_DIR=args.output)

    # Create controller
    controller = PipelineController(args, datetime.now())

    # Run with parallel execution for stratified datasets
    controller.run_complete_pipeline(parallel=True)

    print(f"\n✓ Completed stratified analysis with {args.n_jobs} parallel workers")


def example_3_step_by_step_with_inspection():
    """Example 3: Manual step-by-step execution with data inspection"""
    print("\n" + "=" * 70)
    print("Example 3: Step-by-step execution with inspection")
    print("=" * 70)

    args = Namespace(
        input="data/paul15/paul15.h5",
        output="example_output_3",
        name="manual_run",
        species="human",
        cluster_key="leiden",
        cluster_key_stratification=None,
        embedding_grn="X_draw_graph_fa",
        embedding_hotspot="X_umap",
        raw_count_layer="raw_counts",
        tf_dictionary=None,
        min_genes=200,
        min_counts=500,
        seed=42,
        n_jobs=8,
        skip_qc=False,
        skip_celloracle=False,
        skip_hotspot=False,
        debug=False,
        parallel=False,
        steps=None,
    )

    # Setup
    setup_directories(args.output, os.path.join(args.output, "figures"))
    set_random_seed(args.seed)
    set_scanpy_settings()
    config.update_config(OUTPUT_DIR=args.output)

    # Create controller
    controller = PipelineController(args, datetime.now())

    # Step 1: Load data
    print("\n[Manual Step 1] Loading data...")
    adata = controller.run_step_load()
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Step 2: Preprocessing
    print("\n[Manual Step 2] Preprocessing...")
    adata_preprocessed = controller.run_step_preprocessing()
    print(
        f"  After preprocessing: {adata_preprocessed.n_obs} cells × "
        f"{adata_preprocessed.n_vars} genes"
    )

    # Custom inspection
    print("\n[Inspection] Checking data quality...")
    if "n_genes" in adata_preprocessed.obs.columns:
        print(
            f"  Mean genes per cell: " f"{adata_preprocessed.obs['n_genes'].mean():.1f}"
        )

    # Step 3: Clustering
    print("\n[Manual Step 3] Clustering...")
    adata_clustered = controller.run_step_clustering()
    n_clusters = len(adata_clustered.obs[args.cluster_key].unique())
    print(f"  Identified {n_clusters} clusters")

    # Custom analysis between steps
    print("\n[Custom Analysis] Cluster sizes:")
    cluster_sizes = adata_clustered.obs[args.cluster_key].value_counts()
    for cluster, size in cluster_sizes.items():
        print(f"    Cluster {cluster}: {size} cells")

    # Continue with remaining steps only if desired
    print("\n[Manual Step 4] Running CellOracle...")
    celloracle_result = controller.run_step_celloracle(adata_clustered)

    print("\n✓ Completed manual step-by-step execution")


def example_4_resume_from_checkpoint():
    """Example 4: Resume pipeline from checkpoint"""
    print("\n" + "=" * 70)
    print("Example 4: Resuming from checkpoints")
    print("=" * 70)

    args = Namespace(
        input="data/paul15/paul15.h5",
        output="example_output_4",
        name="resume_run",
        species="human",
        cluster_key="leiden",
        cluster_key_stratification=None,
        embedding_grn="X_draw_graph_fa",
        embedding_hotspot="X_umap",
        raw_count_layer="raw_counts",
        tf_dictionary=None,
        min_genes=200,
        min_counts=500,
        seed=42,
        n_jobs=8,
        skip_qc=False,
        skip_celloracle=False,
        skip_hotspot=False,
        debug=False,
        parallel=False,
        steps=None,
    )

    # Setup
    setup_directories(args.output, os.path.join(args.output, "figures"))
    set_random_seed(args.seed)
    set_scanpy_settings()
    config.update_config(OUTPUT_DIR=args.output)

    # Create controller
    controller = PipelineController(args, datetime.now())

    # First run (will create checkpoints)
    print("\n[First Run] Running full pipeline...")
    controller.run_complete_pipeline(steps=["load", "preprocessing", "clustering"])

    # Second run (will skip completed steps due to checkpoints)
    print("\n[Second Run] Running again - should skip completed steps...")
    controller2 = PipelineController(args, datetime.now())
    controller2.run_complete_pipeline(
        steps=["load", "preprocessing", "clustering", "celloracle"]
    )

    print("\n✓ Checkpoint system automatically skipped completed steps!")


def example_5_custom_stratification_processing():
    """Example 5: Custom processing of stratified datasets"""
    print("\n" + "=" * 70)
    print("Example 5: Custom stratification processing")
    print("=" * 70)

    args = Namespace(
        input="data/paul15/paul15.h5",
        output="example_output_5",
        name="custom_strat",
        species="human",
        cluster_key="leiden",
        cluster_key_stratification="paul15_clusters",
        embedding_grn="X_draw_graph_fa",
        embedding_hotspot="X_umap",
        raw_count_layer="raw_counts",
        tf_dictionary=None,
        min_genes=200,
        min_counts=500,
        seed=42,
        n_jobs=8,
        skip_qc=False,
        skip_celloracle=True,
        skip_hotspot=True,
        debug=False,
        parallel=False,
        steps=None,
    )

    # Setup
    setup_directories(args.output, os.path.join(args.output, "figures"))
    set_random_seed(args.seed)
    set_scanpy_settings()
    config.update_config(OUTPUT_DIR=args.output)

    # Create controller
    controller = PipelineController(args, datetime.now())

    # Run preprocessing and stratification
    controller.run_step_load()
    controller.run_step_preprocessing()
    controller.run_step_stratification()

    # Process only specific stratifications
    print(f"\n[Custom] Found {len(controller.adata_list)} stratifications")
    print(f"  Stratifications: {controller.adata_stratification_list}")

    # Process only the first stratification
    if controller.adata_list:
        print(
            f"\n[Custom] Processing only first stratification: "
            f"{controller.adata_stratification_list[0]}"
        )
        result = controller.process_single_stratification(
            controller.adata_list[0], controller.adata_stratification_list[0]
        )
        print(f"  ✓ Processed: {result}")

    print("\n✓ Completed custom stratification processing")


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("PipelineController Usage Examples")
    print("=" * 70)
    print("\nThis script demonstrates various ways to use the controller.")
    print("Uncomment the example you want to run:\n")

    # Uncomment one example at a time to run it:

    # example_1_run_specific_steps()
    # example_2_stratified_parallel()
    # example_3_step_by_step_with_inspection()
    # example_4_resume_from_checkpoint()
    # example_5_custom_stratification_processing()

    print("\nTo run an example, uncomment it in the __main__ block.")
