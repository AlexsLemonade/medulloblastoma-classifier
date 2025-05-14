# Run experiments using bulk gene expression data
#
# Steven Foltz
# 2023

# set up directories
predict_dir <- here::here("predict")
nb_dir <- here::here("analysis_notebooks")

# notebook file paths
estimate_tumor_purity_filepath <- file.path(nb_dir, "ESTIMATE_tumor_purity.Rmd")
predict_baseline_models_filepath <- file.path(predict_dir, "predict_baseline_models.R")
baseline_models_filepath <- file.path(nb_dir, "baseline_models.Rmd")

# ESTIMATE tumor purity
#rmarkdown::render(estimate_tumor_purity_filepath)

# Run baseline models
source(predict_baseline_models_filepath)
#rmarkdown::render(baseline_models_filepath, params = list(create_models = TRUE))

