#!/bin/bash
# TO REMOVE: Lab-specific script with hardcoded paths. Not portable. Remove or move to personal scripts.

# RUN SMALL TEST ON REGEV DATASET
python run_complete_analysis.py \
 --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/regev_original_with_merged_metadata_scored.h5ad \
 --output results/regev_test_run_update \
 --species human \
 --cluster-key merged_treatment \
 --cluster-key-stratification level_2_annotation \
 --clusters "myCAF,CAF,Acinar,Malignant" \
 --raw-count-layer counts \
 --tf-dictionary /storage/sahuarea/samuele_cancellieri/reference_data/TFs/TG_to_TF_dictionary.pkl \
 --skip-qc \
 --debug

# python run_complete_analysis.py \
#  --input "/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/analyzed/regev_caf-epithelial(malignant-non_malignant).h5ad" \
#  --output results/regev_test_epi-caf_run \
#  --species human \
#  --cluster-key-stratification merged_treatment \
#  --cluster-key level_2_annotation \
#  --raw-count-layer counts \
#  --tf-dictionary /storage/sahuarea/samuele_cancellieri/reference_data/TFs/TG_to_TF_dictionary.pkl \
#  --skip-qc
#  --debug
