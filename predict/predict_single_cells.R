#! /usr/bin/env Rscript
# J. Taroni 2025

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("-o", "--output_file"),
    type = "character", 
    help = "File path for TSV file that will be output with the results"
  ),
  make_option(
    opt_str = c("-t", "--model_type"),
    type = "character",
    help = "Which model type to use (must be 'rf' or 'ktsp')"
  ),
  make_option(
    opt_str = c("-p", "--prop_observed"),
    type = "character",
    default = 0.1,
    help = "The proportion of genes or rules observed in the classifier that need to be detected to retain a cell for classification (default: %default)"
  ),
  make_option(
    opt_str = c("-k, --ktsp_mode"),
    type = "character",
    default = "rule",
    help = "What method should be used to filter out cells in a kTSP model? 'gene' or 'rule' (default: %default)"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Make certain options lowercase to make checking easier
opt$model_type <- tolower(opt$model_type)
opt$ktsp_mode <- tolower(opt$ktsp_mode)

# Error handling
stopifnot(
  "Unsupported model type" = (opt$model_type %in% c("rf", "ktsp")),
  "Unsupported ktsp_mode argument" = (opt$ktsp_mode %in% c("gene", "rule")) 
)

#### Functions -----------------------------------------------------------------

# Load functions
source(here::here("utils/modeling.R"))
source(here::here("utils/single-cell.R"))
source(here::here("utils/convert_gene_names.R"))

#### Directories and files -----------------------------------------------------

# Create output directory if it doesn't exist
dir.create(dirname(opt$output_file), showWarnings = FALSE, recursive = TRUE)

# Directories
models_dir <- here::here("models")
processed_data_dir <- here::here("processed_data")
single_cell_data_dir <- here::here(processed_data_dir, "single_cell")
smartseq_data_dir <- here::here(single_cell_data_dir, "GSE119926")
tenx_data_dir <- here::here(single_cell_data_dir, "GSE155446")

# Baseline models
baseline_models_filepath <- here::here(models_dir, "baseline.rds")

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

#### Read in data --------------------------------------------------------------

mb_subgroups <- c("G3", "G4", "SHH", "WNT")

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

# Use model type option to determine which classifier to use
if (opt$model_type == "ktsp") {
  
  classifier <- baseline_models[[official_model]]$ktsp_unweighted$classifier
  
  
} else {

  classifier <- baseline_models[[official_model]]$rf_weighted$classifier
  
}

# Remove somewhat large object
rm(baseline_models)


if (model_type == "ktsp") {

  # Extract rules in data frame form
  rules_df <- classifier$classifiers |>  
    purrr::map(\(x)
               data.frame(x$TSP)) |>
    dplyr::bind_rows(.id = "subgroup") |>
    dplyr::group_by(subgroup) |>
    dplyr::mutate(rule_num = stringr::str_c("rule", dplyr::row_number())) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("gene"),
                        names_to = "gene_in_pair",
                        values_to = "gene") |>
    dplyr::ungroup()

  test_objects <- sce_files |>
    purrr::map2()
  
}



genes_in_classifier <- classifier$RF_scheme$genes