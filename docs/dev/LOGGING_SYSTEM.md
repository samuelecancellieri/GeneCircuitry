<!-- TO REMOVE (from root): Dev notes — move to docs/dev/LOGGING_SYSTEM.md -->

# GeneCircuitry Logging System

## Overview

A comprehensive logging system has been implemented throughout the GeneCircuitry pipeline to track execution progress, record intermediate steps, and capture errors with full tracebacks.

## Log Files

The logging system creates two separate log files in the `output/logs/` directory:

### 1. `pipeline.log`

- **Purpose**: Records all pipeline steps and operations
- **Content**:
  - Pipeline initialization
  - Step execution (STARTED, COMPLETED, SKIPPED)
  - Intermediate operations (QC, Normalization, Clustering, etc.)
  - Data metrics (n_obs, n_vars, n_clusters, etc.)
  - Checkpoint operations (LOADING, LOADED, SAVING, SAVED)
  - Pipeline completion status (SUCCESS/FAILED)
- **Format**: `YYYY-MM-DD HH:MM:SS - LEVEL - MESSAGE`

### 2. `error.log`

- **Purpose**: Records all errors with full stack traces
- **Content**:
  - Error context (which operation failed)
  - Error message
  - Full Python traceback
- **Format**: `YYYY-MM-DD HH:MM:SS - ERROR - ERROR in [context]: [message]` followed by traceback

## Log Format

```
2024-01-15 14:23:45 - INFO - [Step Name] STATUS - key1=value1, key2=value2
```

Where:

- **Timestamp**: `YYYY-MM-DD HH:MM:SS`
- **Level**: INFO, ERROR
- **Step Name**: Operation identifier (e.g., "Data Loading", "Preprocessing.QC")
- **Status**: STARTED, COMPLETED, LOADING, LOADED, SAVING, SAVED, SKIPPED, FAILED, SUCCESS
- **Details**: Optional key-value pairs with relevant metrics

## Logging Functions

### `setup_logging(output_dir)`

Initializes the logging system. Must be called once at the start of the pipeline.

```python
setup_logging(args.output)
```

Creates:

- `output/logs/pipeline.log` - Main log file
- `output/logs/error.log` - Error log file

### `log_step(step_name, status="STARTED", details=None)`

Logs a pipeline step with optional details.

```python
# Simple step logging
log_step("Data Loading", "STARTED")
log_step("Data Loading", "COMPLETED")

# With details
log_step("Data Loading", "COMPLETED", {
    "n_obs": 1000,
    "n_vars": 2000,
    "input_path": "/path/to/data.h5ad"
})

# Sub-step logging
log_step("Preprocessing.QC", "STARTED")
log_step("Preprocessing.QC", "COMPLETED", {
    "cells_before": 1000,
    "cells_after": 950
})
```

### `log_error(error_context, exception)`

Logs an error with full traceback.

```python
try:
    # ... operation ...
except Exception as e:
    log_error("Data Loading", e)
    raise
```

## Logged Pipeline Steps

### Main Pipeline

- **Pipeline Initialization**: Configuration and setup
- **Data Loading**: Input file reading
- **Preprocessing**: QC, normalization, checkpoint operations
- **Dimensionality Reduction & Clustering**: GRN preprocessing, UMAP, Leiden clustering
- **CellOracle**: Oracle creation, PCA, KNN, GRN inference
- **Hotspot**: Object creation, analysis
- **Generate Summary**: Summary report creation
- **Pipeline Execution**: Overall success/failure with duration

### Controller Steps

- **Controller Initialization**: Controller setup with parameters
- **Controller.LoadData**: Data loading via controller
- **Controller.Preprocessing**: Preprocessing via controller
- **Controller.Stratification**: Data stratification
- **Controller.Clustering**: Clustering via controller
- **Controller.CellOracle**: CellOracle via controller
- **Controller.Hotspot**: Hotspot via controller
- **Controller.GRNAnalysis**: GRN deep analysis via controller

### Intermediate Steps

- **Preprocessing.Checkpoint**: Loading/saving preprocessing checkpoints
- **Preprocessing.QC**: Quality control filtering
- **Preprocessing.Normalization**: Data normalization
- **DimReduction_Clustering.Checkpoint**: Clustering checkpoint operations
- **DimReduction_Clustering.GRNPreprocessing**: GRN-specific preprocessing
- **DimReduction_Clustering.UMAP**: UMAP embedding
- **DimReduction_Clustering.Leiden**: Leiden clustering
- **CellOracle.Checkpoint**: CellOracle checkpoint operations
- **CellOracle.CreateObject**: Oracle object creation
- **CellOracle.PCA**: PCA on Oracle object
- **CellOracle.KNN**: KNN imputation
- **CellOracle.InferGRN**: GRN link inference
- **CellOracle.SaveResults**: Saving CellOracle results
- **Hotspot.Checkpoint**: Hotspot checkpoint operations
- **Hotspot.CreateObject**: Hotspot object creation
- **Hotspot.Analysis**: Hotspot analysis execution

## Example Log Output

### pipeline.log

```
2024-01-15 14:23:30 - INFO - [Pipeline Initialization] STARTED - output_dir=/path/to/output, input=/path/to/data.h5ad
2024-01-15 14:23:35 - INFO - [Pipeline Initialization] COMPLETED
2024-01-15 14:23:35 - INFO - [PipelineController] INITIALIZED - output_dir=/path/to/output, stratification_key=leiden
2024-01-15 14:23:35 - INFO - [Controller Creation] STARTED
2024-01-15 14:23:35 - INFO - [Controller Creation] COMPLETED
2024-01-15 14:23:35 - INFO - [Pipeline Execution] STARTED - mode=complete
2024-01-15 14:23:35 - INFO - [Data Loading] STARTED - input_path=/path/to/data.h5ad
2024-01-15 14:23:40 - INFO - [Data Loading] READING - format=h5ad
2024-01-15 14:23:45 - INFO - [Data Loading] COMPLETED - n_obs=1000, n_vars=2000
2024-01-15 14:23:45 - INFO - [Preprocessing] STARTED - n_obs=1000, n_vars=2000
2024-01-15 14:23:45 - INFO - [Preprocessing.QC] STARTED
2024-01-15 14:24:00 - INFO - [Preprocessing.QC] COMPLETED - cells_before=1000, cells_after=950, genes_before=2000, genes_after=1800
2024-01-15 14:24:00 - INFO - [Preprocessing.Normalization] STARTED
2024-01-15 14:24:15 - INFO - [Preprocessing.Normalization] COMPLETED
2024-01-15 14:24:15 - INFO - [Preprocessing.Checkpoint] SAVED - checkpoint_file=/path/to/output/preprocessed_adata.h5ad
2024-01-15 14:24:15 - INFO - [Preprocessing] COMPLETED - n_obs=950, n_vars=1800
2024-01-15 14:24:15 - INFO - [DimReduction_Clustering] STARTED - n_obs=950, n_vars=1800, cluster_key=leiden
2024-01-15 14:24:15 - INFO - [DimReduction_Clustering.GRNPreprocessing] STARTED
2024-01-15 14:25:00 - INFO - [DimReduction_Clustering.GRNPreprocessing] COMPLETED
2024-01-15 14:25:00 - INFO - [DimReduction_Clustering.UMAP] STARTED
2024-01-15 14:25:30 - INFO - [DimReduction_Clustering.UMAP] COMPLETED
2024-01-15 14:25:30 - INFO - [DimReduction_Clustering.Leiden] STARTED
2024-01-15 14:25:45 - INFO - [DimReduction_Clustering.Leiden] COMPLETED - n_clusters=8
2024-01-15 14:25:45 - INFO - [DimReduction_Clustering] COMPLETED - n_obs=950, n_vars=1800, n_clusters=8
2024-01-15 14:30:00 - INFO - [Pipeline Execution] COMPLETED
2024-01-15 14:30:00 - INFO - [PIPELINE] SUCCESS - total_duration=0:06:30
```

### error.log

```
2024-01-15 14:23:40 - ERROR - ERROR in Data Loading: [Errno 2] No such file or directory: '/path/to/data.h5ad'
Traceback (most recent call last):
  File "/path/to/complete_pipeline.py", line 750, in load_data
    adata = sc.read_h5ad(input_path)
  File "/path/to/scanpy/readwrite.py", line 123, in read_h5ad
    with h5py.File(filename, "r") as f:
  File "/path/to/h5py/_hl/files.py", line 406, in __init__
    fid = make_fid(name, mode, userblock_size, fapl, swmr=swmr)
  File "/path/to/h5py/_hl/files.py", line 173, in make_fid
    fid = h5f.open(name, flags, fapl=fapl)
FileNotFoundError: [Errno 2] No such file or directory: '/path/to/data.h5ad'
```

## Usage Pattern

### Standard Function Pattern

```python
def pipeline_function(adata, **kwargs):
    """Pipeline function with logging."""
    log_step("FunctionName", "STARTED", {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "param1": kwargs.get("param1")
    })

    try:
        # Main operation
        log_step("FunctionName.SubOperation", "STARTED")
        result = perform_operation()
        log_step("FunctionName.SubOperation", "COMPLETED", {
            "result_metric": len(result)
        })

        # Save checkpoint
        if checkpoint:
            log_step("FunctionName.Checkpoint", "SAVING")
            save_checkpoint()
            log_step("FunctionName.Checkpoint", "SAVED", {
                "checkpoint_file": checkpoint_path
            })

        log_step("FunctionName", "COMPLETED", {
            "final_metric": final_value
        })
        return result
    except Exception as e:
        log_error("FunctionName", e)
        raise
```

### Controller Method Pattern

```python
def run_step_operation(self):
    """Controller step with logging."""
    log_step("Controller.Operation", "STARTED")
    try:
        # Validate state
        if self.data is None:
            raise ValueError("Data not available")

        # Execute operation
        result = operation_function(self.data, **params)

        log_step("Controller.Operation", "COMPLETED")
        return result
    except Exception as e:
        log_error("Controller.Operation", e)
        raise
```

## Parallel Processing Notes

In parallel/multiprocess environments (e.g., `_process_stratification_worker`):

- Global logger instances are NOT available in worker processes
- Each stratified analysis creates its own log files in:
  - `output/stratified_analysis/{cluster_name}/logs/pipeline.log`
  - `output/stratified_analysis/{cluster_name}/logs/error.log`
- Worker errors are printed to stderr instead of the main error log

## Benefits

1. **Complete Audit Trail**: Every step is logged with timestamps
2. **Progress Monitoring**: Track pipeline execution in real-time
3. **Error Debugging**: Full tracebacks for error diagnosis
4. **Performance Analysis**: Duration tracking for optimization
5. **Checkpoint Verification**: Log checkpoint load/save operations
6. **Intermediate Step Visibility**: Sub-operations are logged separately
7. **Parallel Execution Tracking**: Each stratification has separate logs

## Integration with Existing System

The logging system integrates seamlessly with:

- **Checkpoint System**: Logs all checkpoint operations
- **Configuration System**: Records config updates
- **Error Handling**: All try-except blocks log errors
- **Controller Pattern**: All controller methods logged
- **Parallel Processing**: Worker processes handle logging independently

## Best Practices

1. **Always initialize logging first**: Call `setup_logging()` before any pipeline operations
2. **Log at function boundaries**: Start and completion of each function
3. **Log intermediate steps**: Major sub-operations within functions
4. **Include relevant metrics**: n_obs, n_vars, file paths, counts, etc.
5. **Wrap in try-except**: Ensure errors are logged before re-raising
6. **Use consistent naming**: Follow the hierarchical naming pattern (Parent.Child)
7. **Log checkpoint operations**: LOADING, LOADED, SAVING, SAVED status

## Future Enhancements

Potential improvements:

- Structured logging (JSON format) for machine parsing
- Log rotation for long-running pipelines
- Performance metrics (memory usage, CPU time)
- Progress bars integrated with logging
- Remote logging support (send logs to monitoring service)
- Log level configuration (DEBUG, INFO, WARNING, ERROR)
- Per-stratification log aggregation in main log file
