# Train and test baseline models to predict MB subgroups
# Steven Foltz, updated 2025-02

# set some variables that guide what downstream steps are followed
create_models <- TRUE # train new models (if FALSE, reads existing models from file)
overwrite <- TRUE # if create_models is also TRUE, overwrite existing models file
seed <- 44 # set initial seed for run_many_models() and before stochastic plots
n_repeats <- 2 # set number of times to repeat in run_many_models()
n_cores <- 3 # set number of cores to use in run_many_models()
ah_date <- "2022-10-30"

# Rationale ----

#We want to use single sample prediction models to predict tumor subgroup in medulloblastoma (MB) patient samples.
#Our gene expression data comes from multiple studies, and different studies were generated using different platforms (array and RNA-seq).

#Our goals are to train different prediction models (kTSP, RF, MM2S, and LASSO) and to assess model performance on a per-platform and per-subgroup basis.
#Ideally, models would perform well regardless of the platform or subgroup of the sample, but differences in sample size mean some groups may be underrepresented in the data available for training.
#Inverse-weighting may boost performance for smaller groups, and we can test this in kTSP and RF models using parameters for weighted (w) and unweighted (unw) analyses.

#To train and test our models, we use the `run_many_models()` function.
#This function allows us to specify what model types we want to analyze, various model parameters, and segregate samples into training/testing sets via a random seed.
#We train and test each model a number of times (e.g. 10 repeats) using different training/testing sets to understand the range of model performance as measured by Kappa and Balanced Accuracy.
#Training/testing sets segregate at the study level, rather than the sample level.

# Setup ----

# Source code and libraries
source(here::here("utils/color_schemes.R"))
source(here::here("utils/convert_gene_names.R"))
source(here::here("utils/modeling.R"))

# Define directories and input/output file paths
processed_data_dir <- here::here("processed_data")
models_dir <- here::here("models")

bulk_genex_filepath <- file.path(processed_data_dir, "bulk_genex.tsv")
bulk_metadata_filepath <- file.path(processed_data_dir, "bulk_metadata.tsv")

baseline_filepath <- file.path(models_dir, "baseline.rds")

# check that existing baseline model file will not be overwritten if overwrite is FALSE
if (create_models &
    (file.exists(baseline_filepath)) &
    !overwrite) {

  stop("Model output file already exists and overwrite is set to FALSE in analysis_notebooks/baseline.Rmd.")

}

# check that files exist if create_models is FALSE
if (!create_models & !file.exists(baseline_filepath)) {

  stop("Model output file does not exist and create_models is set to FALSE in analysis_notebooks/baseline.Rmd.")

}

# set subgroups analyzed in this notebook (canonical MB subgroups)
mb_subgroups <- c("G3", "G4", "SHH", "WNT")

# set up AnnotationHub
set_up_AnnotationHub(ah_date = ah_date)

# Read in essential data ----
bulk_genex_df <- readr::read_tsv(bulk_genex_filepath,
                                 show_col_types = FALSE) |>
  tibble::column_to_rownames(var = "gene")

bulk_metadata_df <- readr::read_tsv(bulk_metadata_filepath,
                                    show_col_types = FALSE) |>
  dplyr::filter(sample_accession %in% names(bulk_genex_df),
                subgroup %in% mb_subgroups)

bulk_genex_df <- bulk_genex_df |>
  dplyr::select(dplyr::all_of(bulk_metadata_df$sample_accession))

check_input_files(genex_df = bulk_genex_df,
                  metadata_df = bulk_metadata_df)

# Train and test MB subgroup models ----

### Create kTSP, RF, MM2S, LASSO, and medullopackage models

#The following code creates both weighted and unweighted kTSP and RF models.
#The MM2S, LASSO, and medullopackage models use default settings.
#Out of all the repeats, we select one 'official' repeat for visualization.
#The 'official' repeat maximizes the median Kappa across all models within the repeat.
#Supplying the same `initial_seed` to each instance of `run_many_models()` ensures the training/testing sets of each repeat are the same for each model.
#We always include the array study GSE37418 in training because a comparison tool (medulloPackage) also trained on that study.

model_types <- c("ktsp_weighted", "ktsp_unweighted",
                 "rf_weighted", "rf_unweighted",
                 "lasso", "mm2s", "medullopackage")

if (create_models) {

  # kTSP and RF models with ktsp_weighted = TRUE and rf_weighted = TRUE (defaults)
  message("kTSP and RF models with ktsp_weighted = TRUE and rf_weighted = TRUE (defaults) ", Sys.time())
  weighted_kTSP_RF_models_list <- run_many_models(genex_df = bulk_genex_df,
                                                  metadata_df = bulk_metadata_df,
                                                  labels = mb_subgroups,
                                                  model_types = c("ktsp", "rf"),
                                                  array_studies_for_training = "GSE37418",
                                                  initial_seed = seed,
                                                  n_repeats = n_repeats,
                                                  n_cores = n_cores,
                                                  ktsp_featureNo = 1000,
                                                  ktsp_n_rules_min = 5,
                                                  ktsp_n_rules_max = 50,
                                                  ktsp_weighted = TRUE,
                                                  rf_num.trees = 500,
                                                  rf_genes_altogether = 50,
                                                  rf_genes_one_vs_rest = 50,
                                                  rf_gene_repetition = 1,
                                                  rf_rules_altogether = 50,
                                                  rf_rules_one_vs_rest = 50,
                                                  rf_weighted = TRUE)

  # Add _weighted to names of kTSP and RF list elements
  weighted_kTSP_RF_models_list <- weighted_kTSP_RF_models_list |>
    purrr::map(\(x) setNames(x,
                             dplyr::case_match(names(x),
                                               "ktsp" ~ "ktsp_weighted",
                                               "rf" ~ "rf_weighted")))

  # kTSP and RF models with ktsp_weighted = FALSE and rf_weighted = FALSE
  message("kTSP and RF models with ktsp_weighted = FALSE and rf_weighted = FALSE ", Sys.time())
  unweighted_kTSP_RF_models_list <- run_many_models(genex_df = bulk_genex_df,
                                                    metadata_df = bulk_metadata_df,
                                                    labels = mb_subgroups,
                                                    model_types = c("ktsp", "rf"),
                                                    array_studies_for_training = "GSE37418",
                                                    initial_seed = seed,
                                                    n_repeats = n_repeats,
                                                    n_cores = n_cores,
                                                    ktsp_featureNo = 1000,
                                                    ktsp_n_rules_min = 5,
                                                    ktsp_n_rules_max = 50,
                                                    ktsp_weighted = FALSE,
                                                    rf_num.trees = 500,
                                                    rf_genes_altogether = 50,
                                                    rf_genes_one_vs_rest = 50,
                                                    rf_gene_repetition = 1,
                                                    rf_rules_altogether = 50,
                                                    rf_rules_one_vs_rest = 50,
                                                    rf_weighted = FALSE)

  # Add _unweighted to names of kTSP and RF list elements
  unweighted_kTSP_RF_models_list <- unweighted_kTSP_RF_models_list |>
    purrr::map(\(x) setNames(x,
                             dplyr::case_match(names(x),
                                               "ktsp" ~ "ktsp_unweighted",
                                               "rf" ~ "rf_unweighted")))

  # merge weighted and unweighted kTSP and RF model lists
  kTSP_RF_models_list <- purrr::map2(weighted_kTSP_RF_models_list,
                                     unweighted_kTSP_RF_models_list,
                                     c)

  # MM2S, LASSO, and medulloPackage models
  message("MM2S, LASSO, and medulloPackage models ", Sys.time())
  mm2s_lasso_medulloPackage_models_list <- run_many_models(genex_df = bulk_genex_df,
                                                           metadata_df = bulk_metadata_df,
                                                           labels = mb_subgroups,
                                                           model_types = c("mm2s",
                                                                           "lasso",
                                                                           "medullopackage"),
                                                           array_studies_for_training = "GSE37418",
                                                           initial_seed = seed,
                                                           n_repeats = n_repeats,
                                                           n_cores = n_cores)

  # merge kTSP, RF, MM2S, and LASSO model lists
  message("combining all models ", Sys.time())
  baseline_list <- purrr::map2(kTSP_RF_models_list,
                               mm2s_lasso_medulloPackage_models_list,
                               c) |>
    purrr::map(\(x) x[!duplicated(names(x))]) # remove duplicate list items

  # select "official" repeat with greatest median Kappa across models
  official_model <- baseline_list |>
    purrr::map_dbl(\(x) model_types |>
                     purrr::map_dbl(
                       \(mtype) x[[mtype]]$cm$overall[["Kappa"]]) |>
                     median()) |>
    which.max()

  for (model_index in 1:length(baseline_list)) {

    if (model_index == official_model) {

      baseline_list[[model_index]][["official_model"]] <- TRUE

    } else {

      baseline_list[[model_index]][["official_model"]] <- FALSE

    }

  }

  # write models to file
  message("write models to file ", Sys.time())
  readr::write_rds(x = baseline_list,
                   file = baseline_filepath)

} else { # if create_models = FALSE and models already exist

  message("reading models inexplicably ", Sys.time())
  baseline_list <- readr::read_rds(baseline_filepath)

}
