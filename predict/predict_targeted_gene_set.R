# Train and test models to predict MB subgroups using a targeted gene set
# Steven Foltz, updated 2025-06

# set some variables that guide what downstream steps are followed
create_models <- TRUE # train new models (if FALSE, reads existing models from file)
overwrite <- TRUE # if create_models is also TRUE, overwrite existing models file
seed <- 6293 # set initial seed for run_many_models()
n_repeats <- 10 # set number of times to repeat in run_many_models()
n_cores <- 3 # set number of cores to use in run_many_models()
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
library(patchwork)

# Define directories and input/output file paths
processed_data_dir <- here::here("processed_data")
models_dir <- here::here("models")

bulk_genex_filepath <- file.path(processed_data_dir, "bulk_genex.tsv")
bulk_metadata_filepath <- file.path(processed_data_dir, "bulk_metadata.tsv")

baseline_filepath <- file.path(models_dir, "baseline.rds")
targeted_filepath <- file.path(models_dir, "targeted_gene_set.rds")
random_filepath <- file.path(models_dir, "random_gene_sets.rds")

# check that existing targeted gene set model file will not be overwritten if overwrite is FALSE
if (create_models &
    (file.exists(targeted_filepath)) &
    !overwrite) {
  
  stop("Targeted gene set model output file already exists and overwrite is set to FALSE in predict/predict_targeted_gene_set.R.")
  
}

# check that existing random gene sets model file will not be overwritten if overwrite is FALSE
if (create_models &
    (file.exists(random_filepath)) &
    !overwrite) {
  
  stop("Random gene sets model output file already exists and overwrite is set to FALSE in predict/predict_targeted_gene_set.R.")
  
}

# check that targeted file exist if create_models is FALSE
if (!create_models & !file.exists(targeted_filepath)) {
  
  stop("Targeted model output file does not exist and create_models is set to FALSE in predict/predict_targeted_gene_set.R.")
  
}

# check that random file exist if create_models is FALSE
if (!create_models & !file.exists(random_filepath)) {
  
  stop("Random model output file does not exist and create_models is set to FALSE in predict/predict_targeted_gene_set.R.")
  
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

# Number of genes in the panel

# The targeted gene panel comprises n_targeted_genes_original genes.
# After converting from SYMBOL to ENSEMBL gene IDs and overlapping with our gene expression data, there are n_targeted_genes_overlap genes available for training using the targeted gene set.
# We use the same number of genes when selecting genes for the random gene set models.

# generate list of 10 random gene sets
set.seed(1034) # randomly chosen seed
random_gene_set_list <- purrr::map(1:10,
                                   \(x) sample(row.names(bulk_genex_df),
                                               size = n_targeted_genes_overlap,
                                               replace = FALSE))

# select random genes from bulk genex df as a list
random_gene_set_bulk_genex_df_list <- random_gene_set_list |>
  purrr::map(\(x) bulk_genex_df[x, ])

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

# The following code creates kTSP (unw), RF (w), and LASSO models using default settings.
# Supplying the same `initial_seed` to each instance of `run_many_models()` ensures the training/testing sets of each repeat are the same for each model.
# Models that are trained and tested on the same sets of samples can be paired to assess relative model performance using each gene set (full, targeted, and random).

model_types <- c("ktsp", "rf", "lasso")

if (create_models) {
  
  # uses the bulk_genex_df reduced to targeted gene set
  targeted_list <- run_many_models(genex_df = bulk_genex_df_targeted,
                                   metadata_df = bulk_metadata_df,
                                   labels = mb_subgroups,
                                   model_types = model_types,
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
                                   rf_weighted = TRUE)
  
  # write models to file
  readr::write_rds(x = targeted_list,
                   file = targeted_filepath)
  
  # uses the bulk_genex_df reduced to random gene sets
  random_list <- random_gene_set_bulk_genex_df_list |>
    purrr::map(\(x) run_many_models(genex_df = x,
                                    metadata_df = bulk_metadata_df,
                                    labels = mb_subgroups,
                                    model_types = model_types,
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
                                    rf_weighted = TRUE))
  
  # write models to file
  readr::write_rds(x = random_list,
                   file = random_filepath)
  
} else { # if create_models = FALSE and models already exist
  
  targeted_list <- readr::read_rds(targeted_filepath)
  random_list <- readr::read_rds(random_filepath)
  
}
