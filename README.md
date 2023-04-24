# Medulloblastoma Single Sample Predictor

## Project overview

### Problem

Medulloblastoma (MB) is the most common form of pediatric brain cancer, with over 1000 new cases affecting children and families worldwide each year.
The four main subtypes of MB (WNT, SHH, Group 3, and Group 4) show different prognoses, especially when stratified by age and sex, and may inform treatment options based on the risk of recurrence.
While MB subtype prediction with gene expression data has been well-characterized in large cohort studies, few single sample predictors have been developed.
Single sample predictors make subtype predictions on an individual sample basis and do not require normalization with the model's training data.
This approach facilitates analyzing new patient data from different platforms (e.g., microarray and RNA-seq), sample handling processes, and tumor purity levels.

### Approach

Existing MB single sample predictors use activated gene pathways or gene expression ratios to classify tumors and were built using transcriptome-wide gene sets.
To complement existing methods, our true single sample predictors incorporate a diverse set of publicly-available data from a wide range of studies, including bulk samples generated from microarray and RNA-seq platforms, as well as single-cell RNA-seq.
The models are built using k top-scoring pairs (kTSP) and Random Forest (RF) approaches that identifies gene pairs whose relative expression levels (e.g., Gene A < Gene B) make informative rules for subtype classification.

### Data

We rely on publicly-available data from the following sources:

- Bulk
  - Microarray: [Weishaupt, H. et al. Batch-normalization of cerebellar and medulloblastoma gene expression datasets utilizing empirically defined negative control genes. Bioinformatics 35, 3357–3364 (2019)](https://academic.oup.com/bioinformatics/article/35/18/3357/5305636) (via [refine.bio](refine.bio))
  - RNA-seq:
    - [Luo, Z. et al. Genomic and Transcriptomic Analyses Reveals ZNF124 as a Critical Regulator in Highly Aggressive Medulloblastomas. Front Cell Dev Biol 9, 634056 (2021)](https://www.frontiersin.org/articles/10.3389/fcell.2021.634056/full)
    - [The Open Pediatric Brain Tumor Atlas (OpenPBTA) Project](https://github.com/AlexsLemonade/OpenPBTA-analysis)
    - [St. Jude Cloud](stjude.cloud)
- Single-cell RNA-seq
  - [Hovestadt, V. et al. Resolving medulloblastoma cellular architecture by single-cell genomics. Nature 572, 74–79 (2019)](https://www.nature.com/articles/s41586-019-1434-6)

## Getting started

### Code repository

Our code is available under a BSD-3-Clause license through this repository: [AlexsLemonade/medulloblastoma-classifier](https://github.com/AlexsLemonade/medulloblastoma-classifier).

### Docker image

We recommend using the `envest/mbssp:R-4.2.2` docker image to create an environment with all necessary dependencies.

### Download data

All the data used in this study is publicly available, and most of it can be downloaded using the `00-download_data.sh` script.
RNA-seq expression count data originating from [St. Jude Cloud](stjude.cloud) requires account registration before selecting data to download.
Follow [these directions](download_stjude_cloud_data.md) to download medulloblastoma expression data from St. Jude Cloud.

After moving St. Jude Cloud data into the folder `data/stjudecloud`, run the following command to download the remaining data:

`bash 00-download_data.sh`

This script ends with an md5sum check of all data based on `data/md5_check_sums.tsv`.

### Gather metadata

To collect metadata from all projects, run

`Rscript 01-gather_metadata.R`

The output goes to two metadata files in the `processed_data/` directory: `bulk_metadata.tsv` and `pseudobulk_metadata.tsv`.
The `bulk_` prefix refers to data generated using bulk sampling methods and processed with array or RNA-seq technologies.
The `pseudobulk_` prefix refers to single-cell RNA-seq (scRNA) data we will analyze at the pseudo-bulk and individual cell levels.

### Gather data

The next step is to bind together bulk gene expression data from each data source and create one expression matrix for all samples (`processed_data/bulk_genex.tsv`).
We also create one Single Cell Experiment object for each scRNA sample (`processed_data/pseudobulk_sce/*`) and one pseudobulk gene expression matrix (`processed_data/pseudobulk_genex.tsv`) where each sample's cells have been collapsed into one column at the gene level.

To gather data, run

`Rscript 02-gather_data.R`

## Data exploration

## Bulk experiments

### Baseline

### Targeted gene panels

### Model dynamics

### ESTIMATE tumor purity

## Single-cell experiments

### Pseudobulk subtype prediction

### Individual cell subtype prediction
