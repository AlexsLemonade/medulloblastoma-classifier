# Run experiments using pseudobulk and single-cell gene expression data
#
# Chante Bethell, Steven Foltz
# 2023

# set up directories
predict_dir <- here::here("predict")
nb_dir <- here::here("analysis_notebooks")

# script file paths
predict_pseudobulk_and_single_cells_filepath <- file.path(predict_dir, "predict_pseudobulk_and_single_cells.R")
pseudobulk_and_single_cells_analysis_notebook_filepath <- file.path(nb_dir, "predict_pseudobulk_and_single_cells.Rmd")

# Predict pseudobulk and single cells; create plot data
source(predict_pseudobulk_and_single_cells_filepath)

# Predict pseudobulk and single cells
rmarkdown::render(pseudobulk_and_single_cells_analysis_notebook_filepath)
