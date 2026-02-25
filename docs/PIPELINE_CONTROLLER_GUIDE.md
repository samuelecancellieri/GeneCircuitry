# Pipeline Controller Guide

## Overview

The `PipelineController` class provides a flexible, modular way to execute the GeneCircuitry analysis pipeline. It allows you to:

- Run specific pipeline steps independently
- Execute stratified analyses in parallel
- Resume from checkpoints automatically
- Customize execution flow for specific needs
- Inspect data between pipeline steps

## Quick Start

### Basic Usage (Complete Pipeline)

```bash
# Run complete pipeline
python examples/complete_pipeline.py --input data.h5ad --output results

# Run with stratification
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --cluster-key-stratification celltype

# Run stratified analysis in parallel
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 4
```

### Running Specific Steps

```bash
# Run only preprocessing and clustering
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps load preprocessing clustering

# Run all except GRN analysis
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --skip-celloracle

# Run only hotspot analysis (assumes preprocessing done)
python examples/complete_pipeline.py \
    --input preprocessed_data.h5ad \
    --output results \
    --steps hotspot
```

## PipelineController Class

### Architecture

The controller pattern separates pipeline execution logic from the individual analysis functions. This provides:

1. **Modularity**: Each step is independent
2. **Flexibility**: Run any combination of steps
3. **Parallelization**: Stratified analyses can run in parallel
4. **State Management**: Controller manages data flow between steps
5. **Checkpoint Integration**: Automatic resume from previous runs

### Available Methods

#### Pipeline Step Methods

```python
controller = PipelineController(args, start_time)

# Individual step execution
controller.run_step_load()                    # Step 1: Load data
controller.run_step_preprocessing()           # Step 2: Preprocessing
controller.run_step_stratification()          # Step 2.5: Stratification
controller.run_step_clustering(adata)         # Step 3: Clustering
controller.run_step_celloracle(adata)         # Step 4: CellOracle
controller.run_step_hotspot(adata)            # Step 5: Hotspot
controller.run_step_grn_analysis(csv_path)    # Step 6: GRN Analysis
```

#### Execution Methods

```python
# Run complete pipeline
controller.run_complete_pipeline(
    steps=['load', 'preprocessing', 'clustering'],
    parallel=False
)

# Process single stratification
controller.process_single_stratification(adata, "TypeA")

# Stratified execution
controller.run_stratified_pipeline_sequential()  # Sequential
controller.run_stratified_pipeline_parallel(n_jobs=4)  # Parallel
```

## Parallel Execution

### When to Use Parallel Execution

Parallel execution is beneficial when:

- You have multiple stratified datasets (e.g., cell types)
- Each stratification is independent
- You have sufficient CPU cores and memory
- Time is more important than resource usage

### How Parallel Execution Works

```python
# Enable parallel execution with --parallel flag
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 4  # Number of parallel workers
```

**Process Flow:**

1. Load and preprocess data once (sequential)
2. Stratify data by specified key
3. Spawn worker processes (up to `n_jobs`)
4. Each worker processes one stratification independently
5. All results are collected when complete

**Resource Considerations:**

- Each worker needs its own memory for the data
- Memory usage ≈ `n_jobs × dataset_size`
- Recommended: `n_jobs` ≤ number of CPU cores
- Monitor memory usage for large datasets

### Parallel vs Sequential

| Aspect         | Parallel                | Sequential                               |
| -------------- | ----------------------- | ---------------------------------------- |
| **Speed**      | Faster (linear speedup) | Slower                                   |
| **Memory**     | High (n_jobs × data)    | Low (1 × data)                           |
| **CPU Usage**  | High (n_jobs cores)     | Low (1 core)                             |
| **Complexity** | Higher                  | Lower                                    |
| **Best For**   | Many stratifications    | Few stratifications or limited resources |

## Pipeline Steps

### Available Step Names

Use these names with the `--steps` flag:

1. **`load`**: Load input data
2. **`preprocessing`**: QC and normalization
3. **`stratification`**: Split by clustering key
4. **`clustering`**: Dimensionality reduction and clustering
5. **`celloracle`**: GRN inference with CellOracle
6. **`hotspot`**: Gene module identification
7. **`grn_analysis`**: Deep GRN analysis
8. **`summary`**: Generate summary report

### Step Dependencies

```
load
  └── preprocessing
       └── stratification (optional)
            └── clustering
                 ├── celloracle
                 └── hotspot
                      └── grn_analysis
                           └── summary
```

## Advanced Usage Examples

### Example 1: Programmatic Usage

```python
from examples.complete_pipeline import PipelineController
from argparse import Namespace
from datetime import datetime

# Create arguments
args = Namespace(
    input="data.h5ad",
    output="results",
    cluster_key_stratification="celltype",
    parallel=True,
    n_jobs=4,
    # ... other args
)

# Create controller
controller = PipelineController(args, datetime.now())

# Run specific steps
controller.run_complete_pipeline(
    steps=['load', 'preprocessing', 'stratification'],
    parallel=False
)

# Access intermediate results
print(f"Loaded {len(controller.adata_list)} stratifications")
print(f"Names: {controller.adata_stratification_list}")
```

### Example 2: Custom Stratification Processing

```python
# Run preprocessing and stratification
controller.run_step_load()
controller.run_step_preprocessing()
controller.run_step_stratification()

# Process only specific stratifications
for i, (adata, name) in enumerate(zip(
    controller.adata_list,
    controller.adata_stratification_list
)):
    if name in ["TypeA", "TypeB"]:  # Only process specific types
        controller.process_single_stratification(adata, name)
```

### Example 3: Step-by-Step with Inspection

```python
# Load and preprocess
controller.run_step_load()
controller.run_step_preprocessing()

# Inspect preprocessed data
print(f"Cells: {controller.adata_preprocessed.n_obs}")
print(f"Genes: {controller.adata_preprocessed.n_vars}")

# Continue only if quality is acceptable
if controller.adata_preprocessed.n_obs > 1000:
    controller.run_step_clustering()
    controller.run_step_celloracle(controller.adata_preprocessed)
```

### Example 4: Parallel Execution with Custom Workers

```python
# Create controller
controller = PipelineController(args, datetime.now())

# Run preprocessing
controller.run_step_load()
controller.run_step_preprocessing()
controller.run_step_stratification()

# Run parallel with specific number of workers
results = controller.run_stratified_pipeline_parallel(n_jobs=8)

print(f"Processed {len(results)} stratifications in parallel")
```

## Checkpoint System Integration

The controller automatically integrates with the checkpoint system:

```bash
# First run - completes preprocessing and clustering
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps load preprocessing clustering

# Second run - skips completed steps, runs CellOracle
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps load preprocessing clustering celloracle
```

**Checkpoint Behavior:**

- Checkpoints are stored in `{output}/logs/*.checkpoint`
- Each stratification has independent checkpoints
- Checkpoints are validated by input hash
- If inputs change, step is re-executed
- Parallel workers share checkpoint logic

## Performance Optimization

### Memory Management

For large datasets with many stratifications:

```bash
# Run with limited parallelization
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 2  # Limit workers to save memory

# Or run sequentially
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype
    # No --parallel flag
```

### CPU Utilization

```bash
# Maximum parallelization for many stratifications
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 16  # Use all available cores

# Also set per-step parallelization
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 8  # 8 parallel stratifications
    --n-jobs 4  # Note: currently --n-jobs controls both
```

### Disk I/O

Checkpoints reduce redundant computation:

```bash
# First run - saves checkpoints
python examples/complete_pipeline.py --input data.h5ad --output run1

# Resume or extend analysis
python examples/complete_pipeline.py --input data.h5ad --output run1 \
    --steps celloracle hotspot  # Skips preprocessing/clustering
```

## Troubleshooting

### Parallel Execution Issues

**Problem**: Out of memory errors with parallel execution

**Solution**: Reduce `--n-jobs` or disable `--parallel`

```bash
# Reduce workers
python examples/complete_pipeline.py --parallel --n-jobs 2

# Or disable parallel
python examples/complete_pipeline.py  # Sequential by default
```

**Problem**: Parallel processes hang or fail

**Solution**: Check for file locking or shared resource conflicts

```bash
# Run in debug mode
python examples/complete_pipeline.py --parallel --debug

# Or run sequentially to isolate issue
python examples/complete_pipeline.py --cluster-key-stratification celltype
```

### Step Execution Issues

**Problem**: Step fails but checkpoint was created

**Solution**: Delete checkpoint to re-run step

```bash
rm results/logs/preprocessing.checkpoint
python examples/complete_pipeline.py --input data.h5ad --output results
```

**Problem**: Want to re-run specific step with new parameters

**Solution**: Delete checkpoint or use new output directory

```bash
# Option 1: Delete checkpoint
rm results/logs/clustering.checkpoint

# Option 2: New output directory
python examples/complete_pipeline.py --output results_v2
```

## Best Practices

1. **Start Small**: Test with one stratification before running all in parallel
2. **Monitor Resources**: Use `htop` or similar to monitor CPU/memory during parallel runs
3. **Use Checkpoints**: Let checkpoint system save intermediate results
4. **Modular Development**: Use step-by-step execution during development
5. **Parallel for Production**: Use parallel execution for final analysis runs
6. **Document Runs**: Use `--name` flag to identify different runs
7. **Version Outputs**: Use different output directories for experiments

## See Also

- [Complete Pipeline Guide](COMPLETE_PIPELINE_GUIDE.md) - Full pipeline documentation
- [Configuration Guide](CONFIG.md) - Configuration options
- [Examples](../examples/) - Example scripts and workflows
