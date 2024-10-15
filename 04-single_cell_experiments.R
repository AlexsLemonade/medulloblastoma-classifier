# Run experiments using pseudobulk and single-cell gene expression data
#
# Chante Bethell, Steven Foltz
# 2023

# set up directories
nb_dir <- here::here("analysis_notebooks")

# notebook file paths
predict_pseudobulk_and_single_cells_filepath <- file.path(nb_dir, "predict_pseudobulk_and_single_cells.Rmd")

# Predict pseudobulk and single cells
rmarkdown::render(predict_pseudobulk_and_single_cells_filepath)
