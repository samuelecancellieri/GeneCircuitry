#!/bin/bash
# TO REMOVE: Lab-specific script with hardcoded paths. Not portable. Remove or move to personal scripts.

# RUN SMALL TEST ON REGEV DATASET
python run_complete_analysis.py \
 --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/regev_original_with_merged_metadata_scored.h5ad \
 --output regev_test_run \
 --species human \
 --cluster-key merged_treatment \
 --clusters "CAF" \
 --cluster-key-stratification level_2_annotation \
 --raw-count-layer counts \
 --debug

# python run_complete_analysis.py \
#  --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/analyzed/regev_malignant.h5ad \
#  --output regev_test_run_only_malignant \
#  --species human \
#  --cluster-key merged_response \
#  --raw-count-layer counts \
#  --debug
