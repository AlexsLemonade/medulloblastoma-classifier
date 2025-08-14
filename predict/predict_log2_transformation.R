# J. Taroni 2025
#
# Test classifiers performance on RNA-seq with and without log2-transformation
# Only predicts on classifier-transformation pairs that are unseen in baseline
# model training to save time
#
# Command line usage:
#     Rscript predict/predict_log2_transformation.R

#### Inputs/Outputs ------------------------------------------------------------

# Directories
processed_data_dir <- here::here("processed_data")
singlecell_data_dir <- here::here(processed_data_dir, "single_cell")
smartseq_data_dir <- here::here(singlecell_data_dir, "GSE119926")
tenx_data_dir <- here::here(singlecell_data_dir, "GSE155446")
models_dir <- here::here("models")
results_dir <- here::here("results")

# Bulk data
bulk_genex_filepath <- file.path(processed_data_dir, "bulk_genex.tsv")
bulk_metadata_filepath <- file.path(processed_data_dir, "bulk_metadata.tsv")

# Pseudobulk data
pseudobulk_metadata_filepath <- file.path(singlecell_data_dir,
                                          "pseudobulk_metadata.tsv")
smartseq_genex_filepath <- file.path(smartseq_data_dir,
                                     "GSE119926_pseudobulk_genex.tsv")
tenx_genex_filepath <- file.path(tenx_data_dir,
                                 "GSE155446_pseudobulk_genex.tsv")

# Models
baseline_filepath <- file.path(models_dir, "baseline.rds")

# Output
results_output_filepath <- file.path(
  results_dir,
  "log2_transformation_prediction_results.rds"
)

#### Functions -----------------------------------------------------------------

# Sourced from files
source(here::here("utils/modeling.R"))
source(here::here("utils/convert_gene_names.R"))

test_rnaseq_prediction <- function(classifier_list,
                                   genex_df,
                                   metadata_df,
                                   model_types = c("ktsp",
                                                   "rf",
                                                   "lasso",
                                                   "mm2s",
                                                   "medullopackage")) {
  # A wrapper around model test functions for convenience
  #
  # For each model_type present, return the results of the corresponding
  # test_*() function in `test_results` and the confusion matrix in `cm`
  #
  # Input:
  #   classifier_list: A list that contains a classifier in the `classifier`
  #                    slot
  #   genex_df: A data frame of gene expression data where the rownames are
  #             Ensembl gene ids -- test functions handle gene id conversion
  #             when necessary
  #   metadata_df: A data frame of metadata that contains the true labels in
  #                `subgroup` and sample ids in `sample_accession`
  #   model_types: What models should be used to predict?
  #
  # Returns: A list where each element corresponds to a model type

  # Labels
  mb_subgroups <- c("G3", "G4", "SHH", "WNT")

  # Initialize list to hold results
  results_list <- list()

  if ("ktsp" %in% model_types) {

    # kTSP (unw)
    results_list$ktsp_unweighted$test_results <- test_ktsp(
      genex_df_test = genex_df,
      metadata_df_test = metadata_df,
      classifier = classifier_list$ktsp_unweighted$classifier,
      labels = mb_subgroups
    )

    # Calculate confusion matrix
    results_list$ktsp_unweighted$cm <- calculate_confusion_matrix(
      predicted_labels = results_list$ktsp_unweighted$test_results$predicted_labels_df$predicted_labels,
      true_labels = metadata_df$subgroup,
      labels = mb_subgroups
    )

  }

  if ("rf" %in% model_types) {

    # RF (w)
    results_list$rf_weighted$test_results <- test_rf(
      genex_df_test = genex_df,
      metadata_df_test = metadata_df,
      classifier = classifier_list$rf_weighted$classifier
    )

    # Calculate confusion matrix
    results_list$rf_weighted$cm <- calculate_confusion_matrix(
      predicted_labels = results_list$rf_weighted$test_results$predicted_labels_df$predicted_labels,
      true_labels = metadata_df$subgroup,
      labels = mb_subgroups
    )

  }

  if ("lasso" %in% model_types) {

    # LASSO
    results_list$lasso$test_results <- test_lasso(
      genex_df_test = genex_df,
      metadata_df_test = metadata_df,
      classifier = classifier_list$lasso$classifier
    )

    # Calculate confusion matrix
    results_list$lasso$cm <- calculate_confusion_matrix(
      predicted_labels = results_list$lasso$test_results$predicted_labels_df$predicted_labels,
      true_labels = metadata_df$subgroup,
      labels = mb_subgroups
    )

  }

  if ("mm2s" %in% model_types) {

    # MM2S
    results_list$mm2s$test_results <- test_mm2s(
      genex_df_test = genex_df,
      metadata_df_test = metadata_df
    )

    # Calculate confusion matrix
    results_list$mm2s$cm <- calculate_confusion_matrix(
      predicted_labels = results_list$mm2s$test_results$predicted_labels_df$predicted_labels,
      true_labels = metadata_df$subgroup,
      labels = c("G3", "G4", "NORMAL", "SHH", "WNT")
    )

  }

  if ("medullopackage" %in% model_types) {

    # medulloPackage
    results_list$medullopackage$test_results <- test_medullo(
      genex_df_test = genex_df,
      metadata_df_test = metadata_df,
      # This is false because we will handle whether or not something is
      # log2 transformed at the input level
      log_transform = FALSE
    )

    # Calculate confusion matrix
    results_list$medullopackage$cm <-calculate_confusion_matrix(
      predicted_labels = results_list$medullopackage$test_results$predicted_labels_df$predicted_labels,
      true_labels = metadata_df$subgroup,
      labels = mb_subgroups
    )

  }

  return(results_list)

}

#### Read in and prepare data --------------------------------------------------

# All bulk gene expression data and metadata
bulk_genex_df <- readr::read_tsv(bulk_genex_filepath,
                                 progress = FALSE,
                                 show_col_types = FALSE)
bulk_metadata_df <- readr::read_tsv(bulk_metadata_filepath,
                                    progress = FALSE,
                                    show_col_types = FALSE)

# All pseudobulk metadata
pseudobulk_metadata_df <- readr::read_tsv(pseudobulk_metadata_filepath,
                                          progress = FALSE,
                                          show_col_types = FALSE)

# Smart-seq2 TPM data -- move gene ids to rownames
smartseq_genex_df <- readr::read_tsv(smartseq_genex_filepath,
                                     progress = FALSE,
                                     show_col_types = FALSE) |>
  tibble::column_to_rownames("gene")

# 10X counts data -- move gene ids to rownames
tenx_genex_df <- readr::read_tsv(tenx_genex_filepath,
                                 progress = FALSE,
                                 show_col_types = FALSE) |>
  tibble::column_to_rownames("gene")

# Read in baseline models
baseline_models <- readr::read_rds(baseline_filepath)

# We need to pick three models to test that each have a different RNA-seq
# dataset held out -- 8 happens to be the official model
models_to_test <- baseline_models[c(6, 8, 10)]
rm(baseline_models)

# Grab the test RNA-seq data from the model object, which has metadata for the
# training and test datasets
heldout_rnaseq_metadata <- models_to_test |>
  purrr::map(\(model)
             model$test_metadata |> dplyr::filter(platform == "RNA-seq")
  )

# Get the heldout RNA-seq as is (i.e., TPM)
heldout_rnaseq <- heldout_rnaseq_metadata |>
  purrr::map(\(metadata_df)
             bulk_genex_df[, c("gene", metadata_df |>
                                 dplyr::pull(sample_accession))] |>
               tibble::column_to_rownames("gene")
  )

# log2 transform all bulk data
heldout_rnaseq_log <- heldout_rnaseq |>
  purrr::map(\(genex_df)
             log2(genex_df + 1)
  )

# log2 transform pseudobulk RNA-seq data
smartseq_log_df <- log2(smartseq_genex_df + 1)
tenx_log_df <- log2(tenx_genex_df + 1)

# Split up metadata for scRNA-seq data
smartseq_metadata_df <- pseudobulk_metadata_df |>
  dplyr::mutate(sample_accession = title) |>
  dplyr::filter(sample_accession %in% colnames(smartseq_genex_df))
tenx_metadata_df <- pseudobulk_metadata_df |>
  dplyr::filter(sample_accession %in% colnames(tenx_genex_df))


#### Prediction ----------------------------------------------------------------

#### Bulk data ####

bulk_tpm_results_list <- list()
bulk_log_results_list <- list()

for (single_repeat in seq_along(models_to_test)) {

  # What dataset will we be testing on?
  single_repeat_dataset <- unique(
    heldout_rnaseq_metadata[[single_repeat]]$study
  )

  # Only medulloPackage has not been tested on straight TPM data
  bulk_tpm_results_list[[single_repeat]] <- test_rnaseq_prediction(
    classifier_list = models_to_test[[single_repeat]],
    genex_df = heldout_rnaseq[[single_repeat]],
    metadata_df = heldout_rnaseq_metadata[[single_repeat]],
    model_types = "medullopackage"
  )
  bulk_tpm_results_list[[single_repeat]]$transformation <- "untransformed"
  bulk_tpm_results_list[[single_repeat]]$dataset <- single_repeat_dataset

  # Other classifiers did not get tested on log2-transformed data
  bulk_log_results_list[[single_repeat]] <- test_rnaseq_prediction(
    classifier_list = models_to_test[[single_repeat]],
    genex_df = heldout_rnaseq_log[[single_repeat]],
    metadata_df = heldout_rnaseq_metadata[[single_repeat]],
    model_types = c("rf", "ktsp", "lasso", "mm2s")
  )
  bulk_log_results_list[[single_repeat]]$transformation <- "log2"
  bulk_log_results_list[[single_repeat]]$dataset <- single_repeat_dataset

}

#### Pseudobulk data ####

# Smart-seq2 pseudobulk log transformed
smartseq_log_result <- test_rnaseq_prediction(
  classifier_list = models_to_test[[2]],  # official model
  genex_df = smartseq_log_df,
  metadata_df = smartseq_metadata_df,
  model_types = c("rf", "ktsp", "mm2s", "medullopackage")
)
smartseq_log_result$transformation <- "log2"
smartseq_log_result$dataset <- "GSE119926"

# Smart-seq2 pseudobulk TPM
smartseq_tpm_result <- test_rnaseq_prediction(
  classifier_list = models_to_test[[2]],  # official model
  genex_df = smartseq_genex_df,
  metadata_df = smartseq_metadata_df,
  model_types = c("rf", "ktsp", "mm2s", "medullopackage")
)
smartseq_tpm_result$transformation <- "untransformed"
smartseq_tpm_result$dataset <- "GSE119926"

# 10x log2 transformed pseudobulk data
tenx_log_result <- test_rnaseq_prediction(
  classifier_list = models_to_test[[2]],  # official model
  genex_df = tenx_log_df,
  metadata_df = tenx_metadata_df,
  model_types = c("rf", "ktsp", "mm2s", "medullopackage")
)
tenx_log_result$transformation <- "log2"
tenx_log_result$dataset <- "GSE155446"

# 10x count pseudobulk data
tenx_count_result <- test_rnaseq_prediction(
  classifier_list = models_to_test[[2]],  # official model
  genex_df = tenx_genex_df,
  metadata_df = tenx_metadata_df,
  model_types = c("rf", "ktsp", "mm2s", "medullopackage")
)
tenx_count_result$transformation <- "untransformed"
tenx_count_result$dataset <- "GSE155446"

#### Compile results and save --------------------------------------------------

# Compile results at the same level of organization
all_results_list <- list(
  bulk_log_results_list[[1]],
  bulk_log_results_list[[2]],
  bulk_log_results_list[[3]],
  bulk_tpm_results_list[[1]],
  bulk_tpm_results_list[[2]],
  bulk_tpm_results_list[[3]],
  smartseq_log_result,
  smartseq_tpm_result,
  tenx_log_result,
  tenx_count_result
)

# Save file
readr::write_rds(all_results_list, file = results_output_filepath)
