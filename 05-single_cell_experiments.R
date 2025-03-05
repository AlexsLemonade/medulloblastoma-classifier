# Run experiments using pseudobulk and single-cell gene expression data
#
# Chante Bethell, Steven Foltz, Jaclyn Taroni
# 2023 - 2025

# set up directories
predict_dir <- here::here("predict")
nb_dir <- here::here("analysis_notebooks")

# script and notebook file paths
predict_pseudobulk_and_single_cells_filepath <- file.path(predict_dir, "predict_pseudobulk_and_single_cells.R")
inspect_genes_filepath <- file.path(nb_dir, "inspect_gene_pairs.Rmd")
pseudobulk_filepath <- file.path(nb_dir, "pseudobulk_analysis.Rmd")
individual_cells_filepath <- file.path(nb_dir, "individual_cells_analysis.Rmd")

# Predict pseudobulk and single cells; create plot data
source(predict_pseudobulk_and_single_cells_filepath)

# Get gene information from kTSP and RF models
rmarkdown::render(inspect_genes_filepath)

# Pseudobulk analysis
rmarkdown::render(pseudobulk_filepath)

# Individual cells analysis
rmarkdown::render(individual_cells_filepath)
