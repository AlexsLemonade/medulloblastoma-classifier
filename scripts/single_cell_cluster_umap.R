#! /usr/bin/env Rscript
# J. Taroni 2025

# Load SCE library
library(SingleCellExperiment)

# Set seed
set.seed(8222)

# Source single-cell functions
source(here::here("utils/single-cell.R"))

# Function
get_umap_coord_and_cluster <- function(sample_title) {
  # Extract the UMAP coordinates, cluster label, and -- for 10X data -- barcode
  # from SingleCellExperiment objects that are read in as part of this function
  # Returns a data frame

  # Read in SCE
  sce_filepath <- sce_files[sample_title]
  sce_object <- readr::read_rds(sce_filepath)

  # Extract UMAP coordinates and cluster information
  umap_df <- data.frame(sample = sample_title,
                        reducedDim(sce_object, "UMAP"),
                        cluster = colData(sce_object)$louvain_10) |>
    dplyr::mutate(cell_index = dplyr::row_number())

  # Tend to names
  names(umap_df) <- c("sample_accession",
                      "UMAP_1",
                      "UMAP_2",
                      "cluster",
                      "cell_index")
  row.names(umap_df) <- NULL

  # Grab the cell identifiers in the SCE object
  cell_identifiers <- colnames(counts(sce_object))

  # If there are underscores, the string following the underscore is a cell
  # barcode that we will use to
  if (any(stringr::str_detect(cell_identifiers, "\\_"))) {
    umap_df$barcode <- stringr::word(cell_identifiers, 2, sep = "_")
  }

  return(umap_df)

}

# Directories
plots_dir <- here::here("plots")
plots_data_dir <- here::here(plots_dir, "data")
processed_data_dir <- here::here("processed_data")
single_cell_data_dir <- here::here(processed_data_dir, "single_cell")
smartseq_data_dir <- here::here(single_cell_data_dir, "GSE119926")
tenx_data_dir <- here::here(single_cell_data_dir, "GSE155446")

# Metadata file
singlecell_metadata_filepath <- here::here(single_cell_data_dir,
                                           "pseudobulk_metadata.tsv")

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

# Output file
output_file <- here::here(plots_data_dir, "single_cell_umap_coord_cluster.tsv")

# Read in metadata
singlecell_metadata_df <- readr::read_tsv(singlecell_metadata_filepath,
                                          show_col_types = FALSE) |>
  dplyr::mutate(sample_accession = title)

# Extract UMAP coords, etc. from each SCE object with custom function
umap_df <- purrr::map(singlecell_metadata_df$title,
                        \(s) get_umap_coord_and_cluster(s)) |>
  dplyr::bind_rows()

# Write out data frame
readr::write_tsv(umap_df, file = output_file)
