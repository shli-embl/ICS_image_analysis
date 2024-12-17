#!/bin/bash
### use example ###
### ./predict_dc_tree.sh data/all_tag_libs.ds.txt output/output_prefix

./R/predict.R -i $1 -f data/feature_set_S8_GFP_relevant.txt -o $2 -b S8_training_run1,S8_calibration_run1 -t data/decision_forest.txt
