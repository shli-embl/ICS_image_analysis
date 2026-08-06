# CODIAC data visualization and gate finding

This workflow is designed to process FCS recordings of single cell image-derived FACS data from BD image-enabled cell sorting (ICS) machine S8. Part 1 involves preprocessing steps of the input data. Part 2 contains data analysis and visualization, and machine learning guided gate finding algorithm used to train models used in the CODIAC study.

## 1. Install conda environment
run 'conda env create -n {your_env_name} -f envs/environment.yaml' to create environment.

## 2. Pre-processing of FCS data from ICS
A snakemake pipeline is provided. 
Dry run: 'snakemake -n'

## 3. Data visualization and gate finding
Individual bash scripts are provided for different analytic steps:
1_plot_umap.sh: Generate UMAP of recorded cells from 50 training tag populations.
2_group_tag.sh: Group tags based on similarity. 
3_train_dc_tree.sh: Generate decision tree models using pre-grouped tag labels. 
4_predict_dc_tree.sh: Predict FACS flows based on defined decision trees/branches. 
5_prune_tree.sh: Remove redundant decision tree/branches to generate final set of gating rules. 
