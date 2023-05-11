# Run experiments using pseudobulk and single-cell gene expression data
#
# Chante Bethell, Steven Foltz
# 2023

# set up directories
nb_dir <- here::here("analysis_notebooks")

# notebook file paths
test_pseudobulk_and_single_cells_filepath <- file.path(nb_dir, "test_pseudobulk_and_single_cells.Rmd")

# Test pseudobulk and single cells
rmarkdown::render(test_pseudobulk_and_single_cells_filepath)
