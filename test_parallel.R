source("utils/modeling.R")

genex_df <- readr::read_tsv("processed_data/bulk_genex.tsv") %>%
  tibble::column_to_rownames(var = "gene")
metadata_df <- readr::read_tsv("processed_data/bulk_metadata.tsv") %>%
  dplyr::filter(sample_accession %in% names(genex_df),
                subgroup %in% c("G3", "G4", "SHH", "WNT"))
genex_df <- genex_df %>%
  dplyr::select(all_of(metadata_df$sample_accession))


run_many_models(genex_df,
                metadata_df,
                model_types = "ktsp", #c("ktsp", "rf", "mm2s", "lasso"),
                initial_seed = 44,
                n_repeats = 1,
                n_cores = 1,
                n_rules_min = 5,
                n_rules_max = NA)