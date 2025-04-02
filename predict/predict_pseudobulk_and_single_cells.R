# Use existing prediction models trained on bulk data to predict subgroups
# of pseudobulk and single cells using scRNA-seq data from Hovestadt, et. al
# https://www.nature.com/articles/s41586-019-1434-6
# and Riemondy, et al. https://doi.org/10.1093/neuonc/noab135
#
# Steven Foltz
# 2023

set.seed(8222)

# Load libraries
source(here::here("utils/modeling.R"))
source(here::here("utils/single-cell.R"))
source(here::here("utils/convert_gene_names.R"))

# Directories
models_dir <- here::here("models")
plots_dir <- here::here("plots")
plots_data_dir <- here::here(plots_dir, "data")
processed_data_dir <- here::here("processed_data")
single_cell_data_dir <- here::here(processed_data_dir, "single_cell")
smartseq_data_dir <- here::here(single_cell_data_dir, "GSE119926")
tenx_data_dir <- here::here(single_cell_data_dir, "GSE155446")

# Input files
# pseudobulk data for each sample was generated in 02-gather_data.R
singlecell_metadata_filepath <- here::here(single_cell_data_dir,
                                           "pseudobulk_metadata.tsv")
smartseq_genex_filepath <- here::here(smartseq_data_dir,
                                      "GSE119926_pseudobulk_genex.tsv")
tenx_genex_filepath <- here::here(tenx_data_dir,
                                  "GSE155446_pseudobulk_genex.tsv")

# baseline models are generated prior to this in the notebook analysis_notebooks/baseline_models.Rmd
# there is one baseline model designated as the official model that is used to test single cells here
# "baseline" means that we used all genes and all available training/test data in modeling
baseline_models_filepath <- here::here(models_dir, "baseline.rds")

# Output files
single_cell_plot_data_filepath <- here::here(plots_data_dir, "single_cell_test_predictions.tsv")
pseudobulk_plot_data_filepath <- here::here(plots_data_dir, "pseudobulk_test_predictions.tsv")

# Get file paths of individual SingleCellExperiment RDS objects
sce_files <- c(fs::dir_ls(path = fs::path(smartseq_data_dir,
                                          "sce"),
                          glob = "*_sce.rds"),
               fs::dir_ls(path = fs::path(tenx_data_dir,
                                          "sce"),
                          glob = "*_sce.rds")
)
# Extract sample titles from the file names
sce_sample_titles <- stringr::word(names(sce_files), start = -1, sep = "/") |>
  stringr::str_remove_all("\\_sce.rds")

# Fail if there are any sample title collisions
if (any(duplicated(sce_sample_titles))) {

  stop("Duplicate single-cell sample title detected!")

} else {

  # Name the vector using the sample title
  names(sce_files) <- sce_sample_titles

}

# Read in data
# Smart-Seq2 pseudobulk data
smartseq_genex_df <- readr::read_tsv(smartseq_genex_filepath,
                                       show_col_types = FALSE) |>
  tibble::column_to_rownames(var = "gene")

# 10x pseudobulk data
tenx_genex_df <- readr::read_tsv(tenx_genex_filepath,
                                 show_col_types = FALSE) |>
  tibble::column_to_rownames(var = "gene")

# Metadata
singlecell_metadata_df <- readr::read_tsv(singlecell_metadata_filepath,
                                          show_col_types = FALSE) |>
  dplyr::mutate(sample_accession = title)

# Split up by study to pass to test_*()
smartseq_metadata_df <- singlecell_metadata_df |>
  dplyr::filter(study == "GSE119926")
tenx_metadata_df <- singlecell_metadata_df |>
  dplyr::filter(study == "GSE155446")

# Read in and prepare models
baseline_models <- readr::read_rds(baseline_models_filepath)
official_model <- which(purrr::map_lgl(baseline_models, \(x) x[["official_model"]]))
classifier_list <- list(ktsp = baseline_models[[official_model]]$ktsp_weighted$classifier,
                        rf = baseline_models[[official_model]]$rf_weighted$classifier)

mb_subgroups <- c("G3", "G4", "SHH", "WNT")

# Predict the subgroup of pseudobulk samples

pseudobulk_test_list <- list(
  "GSE119926" = list("ktsp" = test_ktsp(genex_df_test = smartseq_genex_df,
                                        metadata_df_test = smartseq_metadata_df,
                                        classifier = classifier_list[["ktsp"]],
                                        labels = mb_subgroups),
                     "rf" = test_rf(genex_df_test = smartseq_genex_df,
                                    metadata_df_test = smartseq_metadata_df,
                                    classifier = classifier_list[["rf"]]),
                     "medullopackage" = test_medullo(genex_df_test = smartseq_genex_df,
                                                     metadata_df_test = smartseq_metadata_df)),
  "GSE155446" = list("ktsp" = test_ktsp(genex_df_test = tenx_genex_df,
                                        metadata_df_test = tenx_metadata_df,
                                        classifier = classifier_list[["ktsp"]],
                                        labels = mb_subgroups),
                     "rf" = test_rf(genex_df_test = tenx_genex_df,
                                    metadata_df_test = tenx_metadata_df,
                                    classifier = classifier_list[["rf"]]),
                     "medullopackage" = test_medullo(genex_df_test = tenx_genex_df,
                                                     metadata_df_test = tenx_metadata_df))
)

# Prep pseudobulk plot data

pseudobulk_plot_df <- purrr::map2(pseudobulk_test_list,  # each dataset
                                  names(pseudobulk_test_list),  # and series accession
                                  \(singlecell_dataset, study)
                                  purrr::map2(singlecell_dataset, # list of test objects and
                                              names(singlecell_dataset), # their classifier model types
                                              \(pseudobulk_test, model_type) pseudobulk_test$model_output |>
                                                as.data.frame() |>
                                                tibble::rownames_to_column(var = "sample_accession") |>
                                                tibble::as_tibble() |>
                                                dplyr::select(dplyr::any_of(c(
                                                  "sample_accession",
                                                  "G3",
                                                  "G4",
                                                  "SHH",
                                                  "WNT",
                                                  "best.fit"))
                                                ) |>
                                                dplyr::mutate(model_type = model_type,
                                                              study = study))) |>
  purrr::flatten_df() |>
  dplyr::bind_rows() |>
  dplyr::left_join(singlecell_metadata_df |>
                     dplyr::select(sample_accession,
                                   subgroup),
                   by = "sample_accession") |>
  dplyr::mutate(model_type = dplyr::case_match(model_type,
                                               "ktsp" ~ "kTSP (w)",
                                               "rf" ~ "RF (w)",
                                               "medullopackage" ~ "medulloPackage"))

readr::write_tsv(x = pseudobulk_plot_df,
                 file = pseudobulk_plot_data_filepath)

# Predict the subgroup of single cells

model_test_list <- purrr::map(names(classifier_list), # classifier model types
                              \(model_type) purrr::map(singlecell_metadata_df$title,
                                                       \(x) test_single_cells(sample_acc = x,
                                                                              sce_filepath = sce_files[x],
                                                                              metadata_df = singlecell_metadata_df,
                                                                              labels = mb_subgroups,
                                                                              classifier = classifier_list[[model_type]],
                                                                              platform = "scRNA-seq")) |>
                                purrr::set_names(singlecell_metadata_df$title)) |>
  purrr::set_names(names(classifier_list))

# Prep single cell plot data

get_umap_coord_and_cluster <- function(sample_title) {

  sce_filepath <- sce_files[sample_title]

  sce_object <- readr::read_rds(sce_filepath)

  umap_df <- data.frame(sample = sample_title,
                        SingleCellExperiment::reducedDim(sce_object, "UMAP"),
                        cluster = SummarizedExperiment::colData(sce_object)$louvain_10) |>
    dplyr::mutate(cell_index = dplyr::row_number())

  names(umap_df) <- c("sample_accession", "UMAP_1", "UMAP_2", "cluster", "cell_index")
  row.names(umap_df) <- NULL

  return(umap_df)

}

umap_list <- purrr::map(singlecell_metadata_df$title,
                        \(sample_title) get_umap_coord_and_cluster(sample_title)) |>
  dplyr::bind_rows()

single_cell_plot_df <- purrr::map2(model_test_list, # list of test objects and
                                   names(model_test_list), # their classifier model types
                                   \(model_test, model_type) purrr::map(model_test,
                                                                        \(x) dplyr::bind_cols(x$predicted_labels_df,
                                                                                              tibble::as_tibble(x$model_output)) |>
                                                                          tidyr::separate_wider_delim(cols = sample_accession,
                                                                                                      delim = "_",
                                                                                                      names = c("sample_accession",
                                                                                                                "cell_index"))) |>
                                     dplyr::bind_rows() |>
                                     dplyr::mutate(model_type = model_type)) |>
  dplyr::bind_rows() |>
  dplyr::rowwise() |>
  dplyr::mutate(max = max(G3, G4, SHH, WNT),
                total = sum(G3, G4, SHH, WNT),
                confidence = max/total,
                cell_index = as.numeric(cell_index)) |>
  dplyr::ungroup() |>
  dplyr::left_join(umap_list,
                   by = c("sample_accession", "cell_index")) |>
  dplyr::left_join(singlecell_metadata_df |>
                     dplyr::select(sample_accession,
                                   study,
                                   subgroup,
                                   is_PDX),
                   by = "sample_accession") |>
  dplyr::mutate(model_type = dplyr::case_match(model_type,
                                               "ktsp" ~ "kTSP (w)",
                                               "rf" ~ "RF (w)"))

readr::write_tsv(x = single_cell_plot_df,
                 file = single_cell_plot_data_filepath)
