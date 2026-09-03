#!/bin/bash
### use example ###
### ./group_tag.sh data/all_tag_libs.ds.txt output/output_prefix

./R/cluster_gini_gain.R -i $1 -f data/feature_set_S8_GFP_relevant.txt -o $2_hclust -b S8_training_run1 -l data/localization.txt
./R/cluster_gini_gain_k_means.R -i $1 -f data/feature_set_S8_GFP_relevant.txt -o $2_kmeans -b S8_training_run1 -l data/localization.txt
