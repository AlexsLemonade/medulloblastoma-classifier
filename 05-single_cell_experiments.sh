#!/bin/bash
#
# Jaclyn Taroni
# 2025
#
# Run single-cell predictions and analysis

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set up directories
predict_dir="predict"
analysis_notebooks_dir="analysis_notebooks"
scripts_dir="scripts"
results_dir="results/single_cells_filtered"
mkdir -p $results_dir

# Run prediction on single-cell and pseudobulk data
Rscript "${predict_dir}/predict_pseudobulk.R"

# Render notebook for pseudobulk analysis
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/pseudobulk_analysis.Rmd')"

# Vary the proportion of genes in a model that need to be observed to classify
# an individual cell
for prop_observed in 0 0.05 0.10 0.15 0.20 0.25 0.5; do

  # At each proportion observed value, run kTSP using two different modes of
  # filtering out cells
  for filtering_type in gene rule; do

    Rscript "${predict_dir}/predict_single_cells.R" \
      --output_file "${results_dir}/ktsp_${prop_observed}_${filtering_type}.tsv" \
      --model_type ktsp \
      --prop_observed ${prop_observed} \
      --ktsp_mode ${filtering_type}

  done

  # Random forest models using gene filtering per proportion observed value
  Rscript "${predict_dir}/predict_single_cells.R" \
    --output_file "${results_dir}/rf_${prop_observed}.tsv" \
    --model_type rf \
    --prop_observed ${prop_observed}

done

# Render notebook analyzing results of varying filtering strategies
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/analyze_prop_observed_single_cells.Rmd')"

# Extract UMAP and cluster information from objects
Rscript "${scripts_dir}/extract_single_cell_cluster_umap.R"

# Single-cell visualizations
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/single_cell_viz.Rmd')"

# Set up cell type annotation references
Rscript "${scripts_dir}/set_up_cell_annotation_references.R"
