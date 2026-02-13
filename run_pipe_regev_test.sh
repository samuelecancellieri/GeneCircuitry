#!/bin/bash

# RUN TEST ON REGEV DATASET
# python run_complete_analysis.py \
#  --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/TRNspot_input.h5ad \
#  --output regev_test_run \
#  --species human \
#  --cluster-key merged_response \
#  --cluster-key-stratification level_2_annotation \
#  --raw-count-layer counts

# RUN SMALL TEST ON REGEV DATASET
python run_complete_analysis.py \
 --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/regev_original_with_merged_metadata_scored.h5ad \
 --output regev_test_run \
 --species human \
 --cluster-key merged_response \
 --clusters "Ductal,Acinar" \
 --cluster-key-stratification level_2_annotation \
 --raw-count-layer counts \
 --atac-peaks /storage/sahuarea/samuele_cancellieri/helsinki_lab_cell-lines/cell_lines_counts/PANC1_ATAC_peaks_blacklisted_top1k.narrowPeak \
 --debug
