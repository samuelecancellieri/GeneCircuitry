---
layout: default
title: Parallel Processing Fix
nav_order: 3
parent: Dev Notes
permalink: /GeneCircuitry/dev/parallel-processing-fix/
---

# Parallel Processing Fix — multiprocess Integration

## Problem Solved

**Error:** `AttributeError: Can't pickle local object 'PipelineController.run_stratified_pipeline_parallel.<locals>.worker'`

**Cause:** The nested worker function inside the method couldn't be serialized by Python's pickle module, which is required for multiprocessing.

## Solution

Refactored parallel processing to use:

1. **Module-level worker function** - Fully picklable
2. **multiprocess library** - Better serialization with dill
3. **Dictionary-based work items** - Cleaner data passing

## Changes Made

### 1. Import Strategy

```python
# Try multiprocess first (uses dill for better serialization)
try:
    from multiprocess import Pool
except ImportError:
    # Fallback to standard multiprocessing
    from multiprocessing import Pool
```

### 2. Module-Level Worker Function

```python
def _process_stratification_worker(work_data):
    """
    Worker function for parallel stratification processing.
    Module-level function ensures it can be pickled.
    """
    adata_cluster = work_data['adata_cluster']
    stratification_name = work_data['stratification_name']
    args = work_data['args']
    start_time = work_data['start_time']

    # Process stratification...
    return stratified_output_dir
```

### 3. Updated Parallel Method

```python
def run_stratified_pipeline_parallel(self, n_jobs=None):
    """Run stratified pipeline in parallel using multiprocess."""
    # Prepare work items as dictionaries
    work_items = []
    for adata_cluster, stratification_name in zip(
        self.adata_list, self.adata_stratification_list
    ):
        work_items.append({
            'adata_cluster': adata_cluster,
            'stratification_name': stratification_name,
            'args': self.args,
            'start_time': self.start_time
        })

    # Use module-level worker function
    with Pool(n_jobs) as pool:
        results = pool.map(_process_stratification_worker, work_items)

    return results
```

## Installation

Install the multiprocess library for best results:

```bash
pip install multiprocess
```

Or install all requirements:

```bash
pip install -r requirements.txt
```

## Why multiprocess?

| Feature          | multiprocessing | multiprocess         |
| ---------------- | --------------- | -------------------- |
| Serialization    | pickle          | dill (more powerful) |
| Complex objects  | Limited         | Better support       |
| Lambdas/closures | ❌              | ✅                   |
| Nested functions | ❌              | ✅ (with dill)       |
| AnnData objects  | ✓               | ✓                    |

**Note:** Even with standard multiprocessing, the new module-level worker function works correctly. multiprocess just provides additional robustness.

## Usage

No changes needed from user perspective:

```bash
# Command-line
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 4

# Programmatic
controller = PipelineController(args, start_time)
controller.run_stratified_pipeline_parallel(n_jobs=4)
```

## Technical Details

### Serialization Flow

1. **Main process** prepares work items (list of dictionaries)
2. **Pool.map()** serializes each dictionary
3. **Worker processes** deserialize and execute
4. **Results** are serialized back to main process

### Why Module-Level?

- **Nested functions** exist in local scope, can't be pickled
- **Module-level functions** have global scope, fully picklable
- **Class methods** can be tricky to pickle (need special handling)
- **Top-level functions** are the most reliable for multiprocessing

### Data Structure

Work items are dictionaries for easy serialization:

```python
work_item = {
    'adata_cluster': AnnData(...),      # The dataset
    'stratification_name': 'TypeA',     # Name/identifier
    'args': Namespace(...),             # Controller arguments
    'start_time': datetime(...)         # Pipeline start time
}
```

## Backward Compatibility

✅ **Fully backward compatible**

- Same CLI interface
- Same API methods
- Same output structure
- Falls back to multiprocessing if multiprocess not available

## Testing

```bash
# Verify syntax
python -m py_compile examples/complete_pipeline.py

# Test with data
python examples/complete_pipeline.py \
    --input data/test.h5ad \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 2
```

## Benefits

1. ✅ **No pickle errors** - Module-level worker is fully serializable
2. ✅ **Better serialization** - multiprocess uses dill when available
3. ✅ **More robust** - Handles complex scientific computing objects
4. ✅ **Production ready** - Tested and validated
5. ✅ **User friendly** - Automatic fallback, clear warnings

---

**Status:** ✅ Fixed and Production Ready  
**Date:** October 2025
