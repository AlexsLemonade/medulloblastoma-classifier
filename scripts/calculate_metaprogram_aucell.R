# J. Taroni 2025
#
# For a set of SingleCellExperiment objects and metaprograms, generate a data
# frame that contains the metaprogram AUCell scores for each program in every
# cell in that set of SingleCellExperiment objects

library(optparse)
library(AUCell)

#### Command line options ------------------------------------------------------

option_list <- list(
  make_option(
    opt_str = c("-i", "--sce_input_dir"),
    type = "character",
    help = "Path to directory that contains RDS of SingleCellExperiment objects"
  ),
  make_option(
    opt_str = c("-m", "--metaprogram_file"),
    type = "character",
    help = "File path for TSV file that contains metaprogram genes",
    default = "processed_data/hovestadt-et-al-group-3-4-metaprogram-genes.tsv"
  ),
  make_option(
    opt_str = c("-o", "--output_file"),
    type = "character",
    help = "File path for TSV file that contains metaprogram scores"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Create output directory if it doesn't exist
dir.create(dirname(opt$output_file), showWarnings = FALSE, recursive = TRUE)

#### Functions -----------------------------------------------------------------

source(here::here("utils/convert_gene_names.R"))

process_gene_expression <- function(sce) {
  # Given a SingleCellExperiment object with Ensembl gene identifiers,
  # return a matrix of the contents of the counts slot with genes converted to
  # gene symbols

  # Extract matrix from counts slot
  gene_exp <- SingleCellExperiment::counts(sce)

  gene_exp <- gene_exp |>
    # Gene identifier conversion
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    convert_gene_names(gene_column_before = "gene",
                       gene_column_after = "gene",
                       map_from = "ENSEMBL",
                       map_to = "SYMBOL") |>
    tibble::column_to_rownames("gene") |>
    # Convert back to matrix
    as.matrix()

  # Return matrix
  return(gene_exp)

}

#### Read and process SCEs -----------------------------------------------------

sce_files <- fs::dir_ls(path = fs::path(opt$sce_input_dir),
                        glob = "*_sce.rds")

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

# Get a list of gene expression matrices with symbols as rownames
gene_exp_list <- sce_files |>
  purrr::map(readr::read_rds) |>
  purrr::map(\(sce) process_gene_expression(sce))

#### Read in and prep gene sets ------------------------------------------------

# Read in tidy data frame of metaprograms
gene_set_df <- readr::read_tsv(opt$metaprogram_file)

# Turn this into a GeneSetCollection we can use with AUCell
gene_set_collection <- unique(gene_set_df$metaprogram) |>
  purrr::map(
    # For each gene set
    \(gene_set_name) {
      gene_set_df |>
        # Subset to the rows in that gene set
        dplyr::filter(metaprogram == gene_set_name) |>
        # Grab the Ensembl gene identifiers
        dplyr::pull(gene) |>
        # Create a GeneSet object
        GSEABase::GeneSet(setName = gene_set_name,
                          geneIdType = GSEABase::SymbolIdentifier())
    }
  ) |>
  # Turn the list of GeneSet objects into a GeneSet collection
  GSEABase::GeneSetCollection()

#### Calculate AUCell scores ---------------------------------------------------

auc_df <- gene_exp_list |>
  # Calculate rankings
  purrr::map(\(gene_exp_mat)
             AUCell_buildRankings(gene_exp_mat,
                                  plotStats = FALSE)) |>
  # Calculate AUC
  purrr::map(\(cell_rankings)
             AUCell_calcAUC(rankings = cell_rankings,
                            geneSets = gene_set_collection,
                            # Make the default explicit -- it is appropriate
                            # for both single-cell datasets used in this
                            # project
                            aucMaxRank = ceiling(nrow(cell_rankings) * 0.05))) |>
  # Extract the AUC matrix and transpose such that columns are gene sets
  purrr::map(\(auc_object)
             t(getAUC(auc_object))) |>
  # Add the row index as a column called cell index for downstream joining
  purrr::map(\(auc_matrix)
             auc_matrix |>
               tibble::as_tibble() |>
               tibble::rowid_to_column("cell_index")) |>
  # Create a data frame where the sample identifier (currently the list name)
  # becomes a column called sample_accession
  dplyr::bind_rows(.id = "sample_accession")

#### Write results -------------------------------------------------------------

readr::write_tsv(auc_df, file = opt$output_file)
