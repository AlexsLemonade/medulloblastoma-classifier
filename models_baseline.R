# Run kTSP, Random Forest, MM2S, and LASSO models
# full gene set
# up to 100 rules (revise higher if models include close to that many)
# repeat 10 times

models_dir <- here::here("models")
models_list_output_filepath <- file.path(models, "baseline.rds")

models_list <- run_many_models(genex_df,
                               metadata_df,
                               model_types = c("ktsp", "rf", "mm2s", "lasso"),
                               initial_seed = 44,
                               n_repeats = 10,
                               n_cores = 2,
                               n_rules_min = 5,
                               n_rules_max = 100)

readr::write_rds(models_list,
                 file = models_list_output_filepath)