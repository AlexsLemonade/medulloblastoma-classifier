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
results_dir="results"
filtered_results_dir="${results_dir}/single_cells_filtered"
processed_data_dir="processed_data"
sc_data_dir="${processed_data_dir}/single_cell"
mkdir -p $filtered_results_dir

# Set up file
metaprogram_file="${processed_data_dir}/hovestadt-et-al-group-3-4-metaprogram-genes.tsv"

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
      --output_file "${filtered_results_dir}/ktsp_${prop_observed}_${filtering_type}.tsv" \
      --model_type ktsp \
      --prop_observed ${prop_observed} \
      --ktsp_mode ${filtering_type}

  done

  # Random forest models using gene filtering per proportion observed value
  Rscript "${predict_dir}/predict_single_cells.R" \
    --output_file "${filtered_results_dir}/rf_${prop_observed}.tsv" \
    --model_type rf \
    --prop_observed ${prop_observed}

done

# Render notebook analyzing results of varying filtering strategies
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/analyze_prop_observed_single_cells.Rmd')"

# Extract UMAP and cluster information from objects
Rscript "${scripts_dir}/extract_single_cell_cluster_umap.R"

# Single-cell visualizations
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/single_cell_viz.Rmd')"

# Hovestadt et al. G3/G4 metaprograms
for study in GSE119926 GSE155446; do

  # Handle platform based on study
  if [ ${study} = 'GSE119926' ]; then
	  platform="Smart-seq2"
  else
	  platform="10x"
  fi

  # Generate dataset-specific control gene lists
  Rscript "${scripts_dir}/generate_metaprogram_control_sets.R" \
    --pseudobulk_input_file "${sc_data_dir}/${study}/${study}_pseudobulk_genex.tsv" \
    --metaprogram_file ${metaprogram_file} \
    --platform "${platform}" \
    --output_file "${sc_data_dir}/${study}/${study}_hovestadt-et-al-control-genes.tsv"

  # Calculate program scores
  Rscript "${scripts_dir}/calculate_metaprogram_scores.R" \
    --sce_input_dir "${sc_data_dir}/${study}/sce" \
    --metaprogram_file ${metaprogram_file} \
    --controls_file "${sc_data_dir}/${study}/${study}_hovestadt-et-al-control-genes.tsv" \
    --platform "${platform}" \
    --output_file "${results_dir}/${study}_hovestadt-et-al-group-3-4-metaprogram-scores.tsv"

done
