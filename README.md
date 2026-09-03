# CODIAC image-enabled cell sorting analysis

This repository contains the data-processing, visualization, and gate-finding workflow used for the **CODIAC** analysis of single-cell image-derived flow-cytometry measurements acquired on the BD S8 image-enabled cell sorting (ICS) platform.

The workflow has two main components:

1. **Preprocessing of FCS recordings** — converts FCS files and FlowJo workspaces into cell-level tables, adds gate membership, filters/downsamples cells, and merges samples into a single analysis table.
2. **Visualization and gate finding** — visualizes protein-tag populations, groups tags according to image-derived phenotypes, trains decision-tree classifiers, evaluates a decision forest, and identifies a reduced set of informative sorting gates.

A processed example dataset is included so that the analysis and visualization workflow can be explored without access to the original FCS files.

---

## Repository structure

```text
ICS_image_analysis/
├── README.md
├── License
├── envs/
│   └── environment.yaml
├── p1_preprocess/
│   ├── Snakefile
│   ├── inputs/
│   │   ├── input.yaml
│   │   ├── samples.txt
│   │   ├── feature_set_S8.txt
│   │   ├── feature_set_S8_GFP_relevant.txt
│   │   └── feature_set_S8_imaging_relevant.txt
│   ├── profile/
│   │   └── config.yaml
│   ├── python/
│   │   └── convert_fcs_to_tsv.py
│   └── R/
│       ├── downsampling.R
│       └── merge_table.R
└── p2_visualization/
    ├── 1_plot_umap.sh
    ├── 2_group_tag.sh
    ├── 3_train_dc_tree.sh
    ├── 4_predict_dc_tree.sh
    ├── 5_prune_tree.sh
    ├── R/
    │   ├── umap.R
    │   ├── cluster_gini_gain.R
    │   ├── cluster_gini_gain_k_means.R
    │   ├── decision_tree_training.R
    │   ├── predict.R
    │   └── draw_forest.R
    └── data/
        ├── all_tag_libs.ds.tsv
        ├── localization.txt
        ├── groupings.txt
        ├── decision_forest.txt
        ├── pruned_tree.txt
        ├── tag_group_localization.txt
        └── parameters.S8
```

---

## Installation

The workflow uses Conda to manage Python, R, Bioconductor, and Snakemake dependencies.

```bash
git clone https://github.com/shli-embl/ICS_image_analysis.git
cd ICS_image_analysis

conda env create -n codiac -f envs/environment.yaml
conda activate codiac
```

---

# Quick start: reproduce the analysis using the bundled example data

The simplest way to explore the repository is to start with the processed example dataset:

```text
p2_visualization/data/all_tag_libs.ds.tsv
```

This table contains **62,000 cells from 50 protein-tag populations**. It includes 50,000 cells from `S8_training_run1`, 11,000 cells from `S8_calibration_run1`, and 1,000 cells from `S8_calibration_run2`.

## Step 1. Visualize the tagged cell populations with UMAP

Run either: 
```bash
./1_plot_umap.sh data/all_tag_libs.ds.tsv output/tag_umap.pdf
```
Or
```bash
Rscript R/umap.R \
    -i data/all_tag_libs.ds.tsv \
    -f ../p1_preprocess/inputs/feature_set_S8_GFP_relevant.txt \
    -o output/tag_umap.pdf \
    -b S8_training_run1,S8_calibration_run1
```

The script:

- selects the GFP/NeonGreen-relevant image features defined in `feature_set_S8_GFP_relevant.txt`;
- log-transforms features marked with `log_scale = 1`;
- standardizes the selected features;
- samples 10,000 cells using the fixed random seed used in the study; and
- generates a UMAP containing the protein-tag populations.


---

## Step 2. Group protein tags according to image-derived phenotypes

Two related grouping strategies are provided: hierarchical clustering and k-means clustering.

Run:
```bash
./group_tag.sh data/all_tag_libs.ds.txt output/tag_group
```

This will create grouping information across the 50 training tags using either hierachical clustering or kmeans iteratively at different cluster numbers. Alternatively, one could also run the grouping separately: 

### Hierarchical clustering

```bash
Rscript R/cluster_gini_gain.R \
    -i data/all_tag_libs.ds.tsv \
    -f ../p1_preprocess/inputs/feature_set_S8_GFP_relevant.txt \
    -o output/tag_group_hclust \
    -b S8_training_run1 \
    -l data/localization.txt
```

### k-means clustering

```bash
Rscript R/cluster_gini_gain_k_means.R \
    -i data/all_tag_libs.ds.tsv \
    -f ../p1_preprocess/inputs/feature_set_S8_GFP_relevant.txt \
    -o output/tag_group_kmeans \
    -b S8_training_run1 \
    -l data/localization.txt
```

The scripts first summarize relationships between tagged cell populations using the selected image-derived features and then evaluate cluster solutions containing 2–40 groups. Protein-localization annotations in `data/localization.txt` are used to assess how strongly the inferred groups correspond to broad subcellular localization classes.

Representative outputs include:

```text
output/tag_group_hclust.distmat.txt
output/tag_group_hclust.local.txt
output/tag_group_hclust.lindex.pdf

output/tag_group_kmeans.local.txt
output/tag_group_kmeans.lindex.pdf
```

The grouping assignments used for downstream model training in the CODIAC analysis are provided in:

```text
data/groupings.txt
```

---

## Step 3. Train decision-tree models

Decision trees are trained to distinguish the tag groups defined in `data/groupings.txt` using the GFP/NeonGreen-relevant cytometric and image-derived features.

Run either:
```bash
./train_dc_tree.sh data/all_tag_libs.ds.txt output/decision_tree
```
Or:
```bash
Rscript R/decision_tree_training.R \
    -i data/all_tag_libs.ds.tsv \
    -f ../p1_preprocess/inputs/feature_set_S8_GFP_relevant.txt \
    -o output/decision_tree \
    -b S8_training_run1 \
    -g data/groupings.txt
```

For each grouping scheme, the script:

1. retains cells belonging to the specified training batch;
2. assigns each protein tag to its corresponding group;
3. partitions cells into training and held-out sets;
4. trains a decision tree using repeated five-fold cross-validation; and
5. saves both the trained model and a graphical representation of the tree.

**Outputs**

```text
output/decision_tree.models.RDS
output/decision_tree.dc_tree.<grouping>.pdf
```

The repository also provides the decision-tree rules used in the published analysis in tabular form:

```text
data/decision_forest.txt
```

`decision_forest.txt` contains the image-feature thresholds and branch directions (`>` or `<`) for the decision-tree gates. The current repository provides this rule table as a reference file; conversion of `decision_tree.models.RDS` into `decision_forest.txt` is not currently automated.

---

## Step 4. Apply the decision forest to the training and calibration populations

The precomputed decision-forest rule table can be applied to the processed single-cell data using:

Run either:
```bash
./predict_dc_tree.sh data/all_tag_libs.ds.txt output/decision_forest_prediction
```
Or:
```bash
Rscript R/predict.R \
    -i data/all_tag_libs.ds.tsv \
    -f ../p1_preprocess/inputs/feature_set_S8_GFP_relevant.txt \
    -o output/decision_forest.predict.txt \
    -b S8_training_run1,S8_calibration_run1 \
    -t data/decision_forest.txt
```

For each protein-tag population and batch, the output records the number of cells assigned to each decision-tree branch/gate.


---

## Step 5. Reduce redundancy in the decision forest and generate sorting gates

Use the prediction matrix generated in Step 4 as input to the forest visualization/pruning script:

Run either:
```bash
./prune_tree.sh data/decision_forest.txt output/pruned_tree
```
Or:
```bash
Rscript R/draw_forest.R \
    -i output/decision_forest.predict.txt \
    -o output/decision_tree
```

The script compares the sorting profiles produced by the candidate decision-tree branches and removes highly correlated/redundant dimensions for visualization.

**Outputs**

```text
output/decision_tree.forest.pdf
output/decision_tree.pruned.pdf
output/decision_tree.pruned_withfullgene.pdf
output/decision_tree.pca.full.pdf
output/decision_tree.pca.ensemble.pdf
```

The final reduced gate set used in the CODIAC analysis is provided as:

```text
data/pruned_tree.txt
```

**Note that the pruned tree is a minimal set of gates that contains information to separate training tags, while the original decision forest contains a more rich although redundant informtion to increase the classification power. Practically the sort bin numbers need to be limited to ensure compatibility with the cell sorting capacity of ICS machine. In the article, 10 gates are generated from the automated tree pruning + 2 additional gates were added from the decision forest to increase power. 

---

# Part 1: preprocessing new BD S8 FCS data

Users who want to analyze new S8 recordings can start from FCS files together with their corresponding FlowJo workspace (`.wsp`) files.

## Step 1. Prepare the sample sheet

Edit:

```text
p1_preprocess/inputs/samples.txt
```

The file is tab-separated and contains one sample per row.

| Column | Description |
| --- | --- |
| `id` | Unique sample identifier used for output filenames |
| `protein_id` | Protein/tag identity assigned to the sample |
| `batch` | Experimental batch or recording run |
| `FCS` | Path to the input `.fcs` file |
| `WSP` | Path to the corresponding FlowJo `.wsp` workspace |
| `image_folder` | Root directory containing exported single-cell TIFF images |
| `image_file_prefix` | Prefix used for the exported image filenames |

For example:

```text
id	protein_id	batch	FCS	WSP	image_folder	image_file_prefix
TOMM20_run1	TOMM20	training_run1	/path/to/TOMM20.fcs	/path/to/TOMM20.wsp	/path/to/TOMM20_images	TOMM20
```

The sample identifiers must be unique.

## Step 2. Configure the preprocessing workflow

Edit:

```text
p1_preprocess/inputs/input.yaml
```

At minimum, set the sample sheet, output directory, and feature set:

```yaml
samplesheet: "inputs/samples.txt"
output_dir: "/path/to/codiac_output"
feature_set: "inputs/feature_set_S8.txt"
```

The default preprocessing workflow uses `feature_set_S8.txt` to define the cell-level S8 features retained after FCS conversion.

Before running the current Snakefile, create the expected output directories:

```bash
mkdir -p /path/to/codiac_output/{tables,downsample_tables,data}
```

The R scripts called by the Snakefile are stored without executable permissions in the current repository snapshot. Either make them executable once:

```bash
cd p1_preprocess
chmod +x R/downsampling.R R/merge_table.R
```

or modify the Snakefile to invoke them with `Rscript`.

## Step 3. Convert FCS measurements and FlowJo gates into cell-level tables

For each sample, `python/convert_fcs_to_tsv.py`:

1. reads the FCS event matrix using `readfcs`;
2. loads the matching FlowJo workspace with FlowKit;
3. evaluates the FlowJo gating hierarchy; and
4. appends a Boolean gate-membership column for each gate to the cell-level feature table.

The resulting file is written to:

```text
<output_dir>/tables/<sample_id>.tsv
```

## Step 4. Quality filter and downsample cells

`R/downsampling.R` applies the analysis gate and retains the selected S8 features. In the current workflow, cells are required to pass the gate:

```text
root.Cells.Singlets.SSC.Singlets.FSC.Singlets.RM_Ecc.mNG..Dapi.
```

The gate name must therefore be present in the FlowJo workspace used to generate the cell table, or the script should be edited to match the gate used in a new experiment.

By default, the Snakefile keeps at most **1,000 gated cells per sample** (`n = 1000`). The random seed is fixed to 123 for reproducibility.

The script also reconstructs image paths from the event number, image folder, and image filename prefix.

Outputs include:

```text
<output_dir>/downsample_tables/<sample_id>.ds.tsv
<output_dir>/downsample_tables/<sample_id>.ds.tsv.complete.tsv
```

The `.complete.tsv` file contains the complete cell table with reconstructed image paths, whereas `.ds.tsv` contains the gated/downsampled analysis table.

## Step 5. Merge samples

`R/merge_table.R` combines the downsampled sample tables and appends `prot_id` and `batch` metadata from the sample sheet.

The final preprocessing output is:

```text
<output_dir>/data/all_tag_libs.ds.tsv
```

This table can be used directly as input for the Part 2 visualization and gate-finding workflow.

---

## Running the Snakemake workflow

From the preprocessing directory:

```bash
cd p1_preprocess
```

First perform a dry run:

```bash
snakemake -n -p
```

Run locally:

```bash
snakemake --cores 4
```

A Slurm profile is also provided in `profile/config.yaml`. On a compatible Slurm cluster, the workflow can be submitted with:

```bash
snakemake --profile profile
```

The supplied profile contains cluster-specific resource settings and may need to be adjusted for other computing environments.

---

# Input and reference files

## Feature sets

Three S8 feature lists are provided in `p1_preprocess/inputs/`:

- `feature_set_S8.txt` — broad set of S8 cytometric and image-derived parameters used during preprocessing.
- `feature_set_S8_GFP_relevant.txt` — 16 NeonGreen/GFP-related features used for CODIAC visualization and model training; the `log_scale` column indicates which features are log-transformed for visualization/grouping analyses.
- `feature_set_S8_imaging_relevant.txt` — broader image-relevant feature set.

## Protein localization annotations

```text
p2_visualization/data/localization.txt
```

contains reference subcellular localization annotations for the tagged proteins. These annotations are used to assess the biological coherence of inferred tag clusters.

## Group assignments

```text
p2_visualization/data/groupings.txt
```

contains the tag-to-group assignments evaluated during decision-tree training.

## Decision-tree gate tables

```text
p2_visualization/data/decision_forest.txt
p2_visualization/data/pruned_tree.txt
```

represent, respectively, the full candidate decision forest and the reduced set of gates retained for the CODIAC analysis.

---

# License

This software is distributed under the **GNU General Public License v3.0 (GPL-3.0)**. See `License` for details.

---

# Citation

If you use this workflow, please cite the CODIAC study associated with this repository.

> Citation information will be added upon publication.

---
