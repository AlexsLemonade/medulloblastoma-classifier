# J. Taroni
# 2025
#
# Calculate AUCell scores for each set of human cerebellar developmental cell
# types marker genes from Aldinger et al. 2021.
# Marker genes are generated via extract_cbl_dev_marker_genes.R
#

#### Libraries -----------------------------------------------------------------

library(AUCell)
library(SingleCellExperiment)

#### Functions -----------------------------------------------------------------

calculate_aucell_scores <- function(sce,
                                    sample_title) {
  # Given a SingleCellExperiment object, return a data frame of AUCell scores
  # for all genesets from the CBL dev marker genes
  # The first two columns will contain the sample_title and cell_index for
  # downstream plotting applications

  # Run AUCell using the CBL dev gene sets
  auc_results <- AUCell_run(sce, geneSets = cbl_dev_genesets)

  auc_df <- auc_results@assays@data$AUC |>
    # Transpose
    t() |>
    # Convert to data frame
    as.data.frame() |>
    # Add in sample accession for downstream plotting
    dplyr::mutate(sample_accession = sample_title) |>
    # Add cell index for downstream plotting
    tibble::rowid_to_column("cell_index") |>
    # reordering
    dplyr::select(sample_accession, cell_index, dplyr::everything())

  return(auc_df)

}

#### Files and directories -----------------------------------------------------

processed_data_dir <- here::here("processed_data")
single_cell_data_dir <- here::here(processed_data_dir, "single_cell")
smartseq_data_dir <- here::here(single_cell_data_dir, "GSE119926/sce")
tenx_data_dir <- here::here(single_cell_data_dir, "GSE155446/sce")
plots_data_dir <- here::here("plots/data")

# Get file paths of individual SingleCellExperiment RDS objects
sce_files <- c(fs::dir_ls(path = smartseq_data_dir,
                          glob = "*_sce.rds"),
               fs::dir_ls(path = tenx_data_dir,
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

# CBL dev marker genes
markers_genesets_file <- here::here(
  single_cell_data_dir,
  "aldinger_cbl_dev_genesets.rds"
)

# Output file for AUCell scores
aucell_scores_file <- here::here(
  plots_data_dir,
  "individual_cells_cbl_dev_aucell_scores.tsv"
)

#### Read in data --------------------------------------------------------------

# SingleCellExperiment from both datasets
sce_list <- sce_files |>
  purrr::map(readr::read_rds)

# GeneSetCollection
cbl_dev_genesets <- readr::read_rds(markers_genesets_file)

# Data frame of AUCell scores for all genesets from all samples
auc_df <- sce_list |>
  purrr::imap(
    \(sce, st) calculate_aucell_scores(sce, st)
  ) |>
  dplyr::bind_rows()

# Write AUCell scores to file
readr::write_tsv(auc_df, file = aucell_scores_file)

