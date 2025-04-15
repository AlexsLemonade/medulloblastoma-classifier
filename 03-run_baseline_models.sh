#!/bin/bash
#
# Steven Foltz
# 2025
#
# Run baseline models

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set up directories
predict_dir="predict"

# Run baseline model prediction
Rscript "${predict_dir}/predict_baseline_models.R"
