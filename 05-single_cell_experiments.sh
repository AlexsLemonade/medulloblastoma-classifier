#!/bin/bash
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

# Run prediction on single-cell and pseudobulk data
Rscript "${predict_dir}/predict_pseudobulk_and_single_cells.R"

# Render notebook that inspects the genes included in models
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/inspect_gene_pairs.Rmd')"

# Render notebook for pseudobulk analysis
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/pseudobulk_analysis.Rmd')"

# Render notebook for analyzing individual cell results
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/individual_cells_analysis.Rmd')"