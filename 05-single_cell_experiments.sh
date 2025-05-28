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

# Run prediction on single-cell and pseudobulk data
Rscript "${predict_dir}/predict_pseudobulk.R"

# Render notebook for pseudobulk analysis
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/pseudobulk_analysis.Rmd')"
