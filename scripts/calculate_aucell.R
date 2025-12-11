# J. Taroni
# 2025
#
# Calculate AUCell scores for each geneset supplied in an RDS that is a
# GSEABase::GeneSetCollection
# Gene sets are generated/prepared in `scripts/prepare_marker_gene_sets.R`
#
# Usage:
#   Rscript calculate_aucell.R \
#     --geneset_file processed_data/single_cell/joshi_genesets.rds \
#     --output_file plots/data/joshi_geneset_aucell_scores.tsv

#### Libraries -----------------------------------------------------------------

library(AUCell)
library(SingleCellExperiment)
library(optparse)

#### Options and parsing -------------------------------------------------------

option_list <- list(
  make_option(
    opt_str = c("-g", "--geneset_file"),
    type = "character",
    help = "File path for RDS file that contains genesets as "
  ),
  make_option(
    opt_str = c("-o", "--output_file"),
    type = "character",
    help = "File path for TSV file that will be output with the results"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

#### Functions -----------------------------------------------------------------

calculate_aucell_scores <- function(sce,
                                    sample_title) {
  # Given a SingleCellExperiment object, return a data frame of AUCell scores
  # for all genesets in the GSEABase::GeneSetCollection()
  # The first two columns will contain the sample_title and cell_index for
  # downstream plotting applications

  # Run AUCell using the gene sets
  auc_results <- AUCell_run(sce, geneSets = genesets)

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

#### Read in data --------------------------------------------------------------

# SingleCellExperiment from both datasets
sce_list <- sce_files |>
  purrr::map(readr::read_rds)

# GeneSetCollection
genesets <- readr::read_rds(opt$genesets_file)

# Data frame of AUCell scores for all genesets from all samples
auc_df <- sce_list |>
  purrr::imap(
    \(sce, st) calculate_aucell_scores(sce, st)
  ) |>
  dplyr::bind_rows()

# Write AUCell scores to file
readr::write_tsv(auc_df, file = opt$output_file)
