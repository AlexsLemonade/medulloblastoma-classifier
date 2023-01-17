# Run kTSP, Random Forest, MM2S, and LASSO models
# full gene set
# up to 100 rules (revise higher if models include close to that many)
# repeat 10 times

library(tidyverse)
source("utils/modeling.R")

processed_data_dir <- here::here("processed_data")
models_dir <- here::here("models")

genex_df_input_filepath <- file.path(processed_data_dir, "bulk_genex.tsv")
metadata_df_input_filepath <- file.path(processed_data_dir, "bulk_metadata.tsv")

models_list_output_filepath <- file.path(models_dir, "baseline.rds")

# read in data
genex_df <- read_tsv(genex_df_input_filepath) %>%
  column_to_rownames(var = "gene")

metadata_df <- read_tsv(metadata_df_input_filepath) %>%
  filter(sample_accession %in% names(genex_df),
         subgroup %in% c("G3", "G4", "SHH", "WNT"),
         !is_duplicate)

genex_df <- genex_df %>%
  select(all_of(metadata_df$sample_accession))

# run models
models_list <- run_many_models(genex_df,
                               metadata_df,
                               model_types = c("ktsp", "rf", "mm2s", "lasso"),
                               initial_seed = 44,
                               n_repeats = 10,
                               n_cores = 2,
                               n_rules_min = 5,
                               n_rules_max = 100)

# save models
readr::write_rds(models_list,
                 file = models_list_output_filepath)