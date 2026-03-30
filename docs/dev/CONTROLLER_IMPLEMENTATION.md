---
layout: default
title: Controller Implementation
nav_order: 1
parent: Dev Notes
permalink: /GeneCircuitry/dev/controller-implementation/
---

# Pipeline Controller Implementation

## What Was Implemented

This update transforms the GeneCircuitry pipeline from a monolithic script into a flexible, modular system with the `PipelineController` class.

### Key Features

1. **Modular Execution**
   - Run any combination of pipeline steps
   - Skip already completed steps
   - Resume from intermediate stages

2. **Parallel Processing**
   - Process multiple stratifications simultaneously
   - Configurable worker count
   - Linear speedup for independent datasets

3. **Flexible Control**
   - Command-line interface (`--steps`, `--parallel`)
   - Programmatic API (PipelineController class)
   - Step-by-step execution with inspection

4. **Checkpoint Integration**
   - Automatic resume from previous runs
   - Per-stratification checkpoints
   - Input validation via hashing

## Usage Examples

### Command-Line

```bash
# Complete pipeline
python examples/complete_pipeline.py --input data.h5ad --output results

# Specific steps only
python examples/complete_pipeline.py \
    --steps load preprocessing clustering \
    --input data.h5ad --output results

# Parallel stratified analysis
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 4 \
    --input data.h5ad --output results
```

### Programmatic

```python
from examples.complete_pipeline import PipelineController
from argparse import Namespace
from datetime import datetime

# Setup
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

# Inspect intermediate results
print(f"Stratifications: {controller.adata_stratification_list}")

# Continue with parallel processing
controller.run_stratified_pipeline_parallel(n_jobs=8)
```

## Files Modified/Created

### Modified

- `examples/complete_pipeline.py`
  - Added `PipelineController` class (270+ lines)
  - Added parallel execution support
  - Added `--parallel` and `--steps` CLI flags
  - Refactored main() to use controller

### Created

- `examples/controller_usage_examples.py` (340+ lines)
  - 5 comprehensive usage examples
  - Production-ready code samples

- `docs/pipeline.md`
  - Complete API documentation (merged from PIPELINE_CONTROLLER_GUIDE.md + CONTROLLER_QUICK_REF.md)
  - Performance optimization guide
  - Troubleshooting section

- `verify_controller.py`
  - Automated verification script
  - Checks all components

### Updated

- `README.md`
  - Added modular execution section
  - Added parallel processing examples
  - Updated quick start guide

## Architecture

### PipelineController Class

```python
class PipelineController:
    """Controller for managing GeneCircuitry pipeline execution."""

    # Data management
    adata: AnnData
    adata_preprocessed: AnnData
    adata_list: List[AnnData]
    adata_stratification_list: List[str]

    # Step methods
    run_step_load() -> AnnData
    run_step_preprocessing() -> AnnData
    run_step_stratification() -> Tuple[List, List]
    run_step_clustering(adata) -> AnnData
    run_step_celloracle(adata) -> Tuple
    run_step_hotspot(adata) -> object
    run_step_grn_analysis(path) -> None

    # Execution methods
    run_complete_pipeline(steps=None, parallel=False)
    run_stratified_pipeline_sequential() -> List
    run_stratified_pipeline_parallel(n_jobs) -> List
    process_single_stratification(adata, name) -> str
```

### Parallel Execution Flow

```
1. Load & Preprocess (sequential)
         ↓
2. Stratify by key
         ↓
3. Create worker pool
         ↓
4. Distribute stratifications to workers
         ↓
5. Each worker processes independently:
   - Clustering
   - CellOracle
   - Hotspot
   - Summary
         ↓
6. Collect results
         ↓
7. Generate overall analysis
```

## Performance

### Parallel Speedup

For N stratifications with M workers:

- Sequential time: N × T (where T = time per stratification)
- Parallel time: ⌈N/M⌉ × T
- Speedup: ~M× (ideal case)

### Memory Usage

- Sequential: ~1× dataset size
- Parallel (M workers): ~M× dataset size

### Recommendations

- **Small datasets (<1GB)**: Use parallel with max workers
- **Medium datasets (1-10GB)**: Use parallel with 2-4 workers
- **Large datasets (>10GB)**: Use sequential or 2 workers

## Benefits

1. **Development**
   - Test individual steps quickly
   - Inspect intermediate results
   - Faster iteration cycles

2. **Production**
   - Process multiple cell types in parallel
   - Resume from failures
   - Optimize resource usage

3. **Research**
   - Run ablation studies
   - Compare different parameters
   - Reuse preprocessing results

## Verification

Run the verification script to confirm installation:

```bash
python verify_controller.py
```

Expected output:

```
✓ Class Implementation: PASS
✓ Method Implementation: PASS
✓ CLI Flags: PASS
✓ Required Imports: PASS
✓ Documentation Files: PASS
```

## Next Steps

1. **Try Examples**: Run `examples/controller_usage_examples.py`
2. **Read Guides**: Check `docs/pipeline.md` (Pipeline Overview)
3. **Run Pipeline**: Use new flags with your data
4. **Optimize**: Tune `--n-jobs` for your system

## Support

For questions or issues:

- See documentation in `docs/`
- Check examples in `examples/`
- Review this implementation summary

---

**Implementation Date**: October 2025  
**Version**: 1.0.0  
**Status**: Production Ready ✓
