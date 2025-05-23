#!/bin/bash
#
# Steven Foltz
# 2025
#
# Run bulk experiments

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set up directories
analysis_notebooks_dir="analysis_notebooks"

# Render notebook for bulk analysis
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/bulk_analysis.Rmd')"

# Render notebook for ESTIMATE analysis
Rscript -e "rmarkdown::render('${analysis_notebooks_dir}/ESTIMATE_tumor_purity.Rmd')"
