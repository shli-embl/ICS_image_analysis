#!/bin/bash
### use example ###
### ./plot_umap.sh data/all_tag_libs.ds.tsv output/output.umap.pdf

./R/umap.R -i $1 -f data/feature_set_S8_GFP_relevant.txt -o $2 -b S8_training_run1,S8_calibration_run1