#!/bin/bash
#
# Steven Foltz
# 2025
#
# Run targeted gene set experiments

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set up directories
predict_dir="predict"
analysis_notebooks_dir="analysis_notebooks"

# Run targeted gene set model prediction
Rscript "${predict_dir}/predict_targeted_gene_set.R"

# Render notebook for targeted gene set analysis
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/targeted_gene_set.Rmd')"
