# Train and test models to predict MB subgroups using a targeted gene set
# Steven Foltz, updated 2025-06

# set some variables that guide what downstream steps are followed
create_models <- TRUE # train new models (if FALSE, reads existing models from file)
overwrite <- TRUE # if create_models is also TRUE, overwrite existing models file
seed <- 44 # set initial seed for run_many_models()
n_repeats <- 10 # set number of times to repeat in run_many_models()
n_cores <- 10 # set number of cores to use in run_many_models()
ah_date <- "2022-10-30"

# Rationale

# Experiments that generate RNA expression data may use a targeted gene panel as a more cost-effective approach to analyze a subset of genes thought to have greater relevance for disease, especially in clinical settings.
# However, when a model is trained using whole transcriptome data, the genes selected for that model will generally show only partial or no overlap with a targeted gene panel.
# Even though kTSP and RF models can tolerate missing genes, prediction performance on new samples will likely suffer as a result of having fewer gene:gene comparisons to score.
# A model trained using a targeted panel of genes from the start could be applied to sample data that includes more genes, however there may be a decrease in model performance due to starting with a smaller feature space.
# This leads to two key questions:

# 1. Given a targeted gene set of cancer-related genes, how does model performance compare between targeted gene set models and full gene set models?
# 2. Does the identity of the genes in the targeted gene set matter? In other words, can a random gene set of the same size perform just as well as the targeted gene set?

# We will use the [Nanostring nCounter PanCancer IO 360 Panel](https://nanostring.com/products/ncounter-assays-panels/oncology/pancancer-io-360/) targeted gene set published by Nanostring, available for download [here](https://nanostring.com/products/ncounter-assays-panels/panel-selection-tool/) and found in this repository at `processed_data/NS_IO_360_v1.0_Genes.tsv`.

# Setup

# Source code and libraries
source(here::here("utils/convert_gene_names.R"))
source(here::here("utils/modeling.R"))

# Define directories and input/output file paths
processed_data_dir <- here::here("processed_data")
models_dir <- here::here("models")
random_dir <- here::here(models_dir, "targeted_gene_set_random")

bulk_genex_filepath <- file.path(processed_data_dir, "bulk_genex.tsv")
bulk_metadata_filepath <- file.path(processed_data_dir, "bulk_metadata.tsv")

baseline_filepath <- file.path(models_dir, "baseline.rds")
targeted_filepath <- file.path(models_dir, "targeted_gene_set.rds")

# check that existing targeted gene set model file will not be overwritten if overwrite is FALSE
if (create_models &
    (file.exists(targeted_filepath)) &
    !overwrite) {

  stop("Targeted gene set model output file already exists and overwrite is set to FALSE in predict/predict_targeted_gene_set.R.")

}

# check that targeted file exist if create_models is FALSE
if (!create_models & !file.exists(targeted_filepath)) {

  stop("Targeted model output file does not exist and create_models is set to FALSE in predict/predict_targeted_gene_set.R.")

}

# We will need this random gene set function to run inside tryCatch()
# Given a random index (this sets the seed), we pick a random
train_and_write_random_model <- function(random_index) {

  # The targeted gene panel comprises n_targeted_genes_original genes.
  # After converting from SYMBOL to ENSEMBL gene IDs and overlapping with our gene expression data, there are n_targeted_genes_overlap genes available for training using the targeted gene set.
  # We use the same number of genes when selecting genes for the random gene set models.

  # # select random genes from bulk genex df
  set.seed(random_index) # randomly chosen seed
  random_gene_set <- sample(row.names(bulk_genex_df),
                            size = n_targeted_genes_overlap,
                            replace = FALSE)

  # subset bulk_genex_df to random gene set
  random_gene_set_bulk_genex_df <- bulk_genex_df[random_gene_set, ]

  # train and test models using random gene set
  random_kTSP_RF_weighted_models_list <- run_many_models(
    genex_df = random_gene_set_bulk_genex_df,
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
  random_kTSP_RF_weighted_models_list <- random_kTSP_RF_weighted_models_list |>
    purrr::map(\(x) setNames(x,
                             dplyr::case_match(names(x),
                                               "ktsp" ~ "ktsp_weighted",
                                               "rf" ~ "rf_weighted",
                                               .default = names(x))))

  # Random gene sets for unweighted kTSP, RF, and lasso
  random_kTSP_RF_lasso_unweighted_models_list <- run_many_models(
    genex_df = random_gene_set_bulk_genex_df,
    metadata_df = bulk_metadata_df,
    labels = mb_subgroups,
    model_types = c("ktsp", "rf", "lasso"),
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
  random_kTSP_RF_lasso_unweighted_models_list <- random_kTSP_RF_lasso_unweighted_models_list |>
    purrr::map(\(x) setNames(x,
                             dplyr::case_match(names(x),
                                               "ktsp" ~ "ktsp_unweighted",
                                               "rf" ~ "rf_unweighted",
                                               .default = names(x))))

  # combine random lists
  random_list <- purrr::map2(
    random_kTSP_RF_weighted_models_list,
    random_kTSP_RF_lasso_unweighted_models_list,
    c) |>
    purrr::map(\(y) y[!duplicated(names(y))]) # remove duplicate list items

  # write random list to random folder
  random_filepath <- here::here(random_dir,
                                glue::glue("targeted_gene_set.random_",
                                           random_index, ".rds"))

  readr::write_rds(x = random_list,
                   file = random_filepath)

}

# set subgroups analyzed in this notebook (canonical MB subgroups)
mb_subgroups <- c("G3", "G4", "SHH", "WNT")

# Read in essential data
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

# targeted genes (Nanostring IO 360 gene panel)
targeted_genes_filepath <- file.path(processed_data_dir,
                                     "NS_IO_360_v1.0_Genes.tsv")

targeted_genes_df <- readr::read_tsv(file = targeted_genes_filepath,
                                     show_col_types = FALSE)

n_targeted_genes_original <- nrow(targeted_genes_df)

# convert targeted gene set from SYMBOL to ENSEMBL
targeted_genes_df <- targeted_genes_df |>
  convert_gene_names(gene_column_before = "Gene",
                     gene_column_after = "gene",
                     map_from = "SYMBOL",
                     map_to = "ENSEMBL",
                     ah_date = ah_date)

# reduce bulk genex df to overlap with targeted genes
bulk_genex_df_targeted <- bulk_genex_df[row.names(bulk_genex_df) %in% targeted_genes_df$gene, ]

n_targeted_genes_overlap <- nrow(bulk_genex_df_targeted)

# Train and test MB subgroup models using targeted gene set

### Can we use MM2S with this targeted gene set?

# convert overlapping targeted gene set from ENSEMBL to ENTREZID for MM2S
targeted_genes_df_ENTREZID <- bulk_genex_df_targeted |>
  tibble::rownames_to_column(var = "gene") |>
  convert_gene_names(gene_column_before = "gene",
                     gene_column_after = "gene",
                     map_from = "ENSEMBL",
                     map_to = "ENTREZID",
                     ah_date = ah_date)

# MM2S uses pathways defined by genes in HumanGMT$genesets
# Only pathways with between 20-100 genes overlapping with expression data are scored
n_mm2s_common_pathways <- purrr::map(MM2S::HumanGMT$genesets,
                                     \(x) sum(x %in% targeted_genes_df_ENTREZID$gene)) |>
  purrr::map_lgl(\(x) x >= 20 & x <= 100) |>
  sum() # = total number of pathways matching gene overlap criteria for scoring

# MM2S assumes there are at least 24 pathways that get scored
ifelse(n_mm2s_common_pathways >= 24,
       "We can use MM2S!", "We cannot use MM2S.")

# MM2S uses a pathway-based approach, and those pathways are pre-defined in `MM2S::HumanGMT$genesets`.
# To return a GSVA score for a particular pathway, MM2S requires that between 20 and 100 genes overlap between the test data and the set of genes in that pathway.
# MM2S later requires that there are at least 24 such pathways in the data.
# However, with our targeted gene set, only n_mm2s_common_pathways pathways match the `20 <= overlapping genes <= 100` criteria.
# MM2S fails to run because an `NA` value is introduced to `FeatureSelection` and `NorthcottFeatures` inside `MM2S::MM2S.human`, resulting in invalid column selection in `HumanTestGSVA <- HumanTestGSVA[, NorthcottFeatures, drop = FALSE]`.

# medullopackage also fails to run

### Create kTSP, RF, and LASSO models

# The following code creates kTSP (weighted and unweighted), RF (weighted and unweighted), and LASSO models using default settings.
# Supplying the same `initial_seed` to each instance of `run_many_models()` ensures the training/testing sets of each repeat are the same for each model.

if (create_models) {

  # Targeted gene set for weighted kTSP and RF
  targeted_kTSP_RF_weighted_models_list <- run_many_models(
    genex_df = bulk_genex_df_targeted,
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
  targeted_kTSP_RF_weighted_models_list <- targeted_kTSP_RF_weighted_models_list |>
    purrr::map(\(x) setNames(x,
                             dplyr::case_match(names(x),
                                               "ktsp" ~ "ktsp_weighted",
                                               "rf" ~ "rf_weighted",
                                               .default = names(x))))

  # Targeted gene set for unweighted kTSP, RF, and lasso
  targeted_kTSP_RF_lasso_unweighted_models_list <- run_many_models(
    genex_df = bulk_genex_df_targeted,
    metadata_df = bulk_metadata_df,
    labels = mb_subgroups,
    model_types = c("ktsp", "rf", "lasso"),
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
  targeted_kTSP_RF_lasso_unweighted_models_list <- targeted_kTSP_RF_lasso_unweighted_models_list |>
    purrr::map(\(x) setNames(x,
                             dplyr::case_match(names(x),
                                               "ktsp" ~ "ktsp_unweighted",
                                               "rf" ~ "rf_unweighted",
                                               .default = names(x))))

  # combine targeted lists
  targeted_list <- purrr::map2(targeted_kTSP_RF_weighted_models_list,
                               targeted_kTSP_RF_lasso_unweighted_models_list,
                               c) |>
    purrr::map(\(x) x[!duplicated(names(x))]) # remove duplicate list items

  # write models to file
  readr::write_rds(x = targeted_list,
                   file = targeted_filepath)

  # Random gene sets for weighted kTSP and RF

  # run train_and_write_random_model() until it succeeds 1000 times
  # tryCatch() allows the random gene set to fail and then tries the next random index
  # after 1000 successes, exit the while loop
  n_success <- 0
  random_index <- 0
  while(n_success < 1000) {

    random_index <- random_index + 1

    tryCatch(
      {
        print(glue::glue("Trying random index ", random_index))
        train_and_write_random_model(random_index)
        n_success <- n_success + 1
        print(glue::glue("Success! ", random_index))
      },
      error = function(msg) {
        print(glue::glue("Random index ", random_index, " failed"))
      }
    )

  }

}
