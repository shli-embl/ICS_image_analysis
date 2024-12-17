#!/bin/bash
### use example ###
### ./train_dc_tree.sh data/all_tag_libs.ds.txt output/output.umap.pdf

./R/decision_tree_training.R -i $1 -f data/feature_set_S8_GFP_relevant.txt -o $2 -b S8_training_run1 -g data/groupings.txt