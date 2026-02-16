#!/usr/bin/env python
"""
TRNspot Complete Analysis - Quick Runner
=========================================
Convenient wrapper to run the complete TRNspot pipeline.

This script provides an easy way to run the full analysis including:
- Data preprocessing and QC
- Dimensionality reduction and clustering
- CellOracle GRN inference
- Hotspot gene module identification

Usage:
    # Run with defaults (example dataset)
    python run_complete_analysis.py

    # Run with your own data
    python run_complete_analysis.py --input your_data.h5ad

    # Skip specific analyses
    python run_complete_analysis.py --skip-celloracle
    python run_complete_analysis.py --skip-hotspot

    # Custom settings
    python run_complete_analysis.py --seed 123 --n-jobs 16
"""

import sys
import os

# Import and run the complete pipeline from its proper package location
from trnspot.pipeline.controller import main

if __name__ == "__main__":
    sys.exit(main())
