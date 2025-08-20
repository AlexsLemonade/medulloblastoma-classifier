# J. Taroni 2025
#
# For a set of SingleCellExperiment objects and metaprograms, generate a data
# frame that contains the metaprogram scores for each program in every cell
# in that set of SingleCellExperiment objects
#
# If 10X data, the normalized, log transformed data in the logcounts is used to
# calculate the scores
# If Smart-seq2, the log2(TPM + 1) is used to calculate the scores
#
# Use generate_metaprogram_control_sets.R to generate the control file passed
# to this script!
#
# Example usage:
#
#   Rscript scripts/calculate_metaprogram_scores.R
#     --sce_input_dir processed_data/single_cell/GSE155446/sce \
#     --metaprogram_file processed_data/hovestadt-et-al-group-3-4-metaprogram-genes.tsv \
#     --controls_file processed_data/single_cell/GSE155446/GSE155446_hovestadt-et-al-control-genes.tsv \
#     --platform 10X \
#     --output_file results/GSE155446_hovestadt-et-al-group-3-4-metaprogram-scores.tsv
#

library(optparse)

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
    opt_str = c("-c", "--controls_file"),
    type = "character",
    help = "File path for TSV file that contains control gene lists"
  ),
  make_option(
    opt_str = c("-p", "--platform"),
    type = "character",
    help = "What kind of scRNA-seq data? 10x or Smart-seq2"
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

process_gene_expression <- function(sce,
                                    platform) {
  # Given a SingleCellExperiment object, return a matrix of relative expression
  # values that can be passed to calculate_metaprogram_score()
  #
  # Input:
  #   sce: A SingleCellExperiment object
  #   platform: 10X or Smart-seq2? How the data is handled upstream is different
  #             between platforms, so how we calculate relative expression
  #             changes. For 10X, grab what is in logcounts(); for Smart-seq2,
  #             log2(TPM + 1) using the TPM in the counts slot.
  #
  # Returns: A mean-centered gene expression matrix -- "relative expression"
  #          per Hovestadt et al. and Tirosh et al. -- rownames are gene symbols

  # Extract the data to be used, depending on the platform
  if (tolower(platform) == "smart-seq2") {

    # The counts slot for the Smart-seq2 data contains TPM
    gene_exp <- SingleCellExperiment::counts(sce)
    # That we will log2 transform
    gene_exp <- log2(gene_exp + 1)

  } else if (tolower(platform) == "10x") {

    # Normalized, log transformed
    gene_exp <- SingleCellExperiment::logcounts(sce)

  } else {

    stop("Unsupported platform!")

  }

  # Convert to gene symbol
  gene_exp <- gene_exp |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    convert_gene_names(gene_column_before = "gene",
                       gene_column_after = "gene",
                       map_from = "ENSEMBL",
                       map_to = "SYMBOL") |>
    tibble::column_to_rownames("gene")

  # Mean center as in Hovestadt et al. and Tirosh et al. to get "relative
  # expression"
  gene_exp <- t(scale(t(gene_exp), center = TRUE, scale = FALSE))

  # Return relative expression matrix
  return(gene_exp)

}

calculate_metaprogram_score <- function(gene_exp,
                                        sample_name,
                                        program_name,
                                        program_genes,
                                        control_genes) {
  # Given a mean-centered gene expression matrix ("relative expression"), a set
  # of metaprogram genes, and a set of control genes, calculate a per-cell
  # metaprogram score, which is mean(program_genes) - mean(control_genes)
  #
  # Input:
  #   gene_exp: A mean-centered gene expression matrix
  #   sample_name: Sample identifier that will be used in the data frame that
  #                is returned
  #   program_name: The name of the metaprogram that will be used in the data
  #                 frame that is returned
  #   program_genes: A vector of gene symbols for genes in the metaprogram
  #   control_genes: A vector of control gene symbols for the metaprogram
  #
  # Returns: A data.frame with sample_accession, cell_index, metaprogram, and
  #          metaprogram_score

  # Subset to program genes
  program_genes <- intersect(rownames(gene_exp), program_genes)
  program_gene_exp <- gene_exp[program_genes, ]

  # We don't have to worry about control genes not being in the object, since
  # they are dataset-specific (derived from pseudobulk data)
  control_gene_exp <- gene_exp[control_genes, ]

  # The score is the average relative expression of metaprogram genes minus the
  # average relative expression of the control genes
  program_scores <- colMeans(program_gene_exp) - colMeans(control_gene_exp)

  # Return a data frame of the scores
  return(
    data.frame(
      sample_accession = sample_name,
      cell_index = 1:length(program_scores),
      metaprogram = program_name,
      metaprogram_score = program_scores
    )
  )

}

metaprogram_wrapper <- function(program) {
  # A convenient wrapper to pass to purrr::map()
  # Not intended to be used outside of this script

  # Grab the program gene symbols
  program_genes <- metaprogram_df |>
    dplyr::filter(metaprogram == program) |>
    dplyr::pull(gene)

  # Grab the relevant control set gene symbols
  control_genes <- control_genes_df |>
    dplyr::filter(metaprogram == program) |>
    dplyr::pull(gene)

  # Return a list of data frames with metaprogram scores for this program
  return(
    purrr::imap(gene_exp_list,
                \(gene_exp_mat, sample_name)
                calculate_metaprogram_score(
                  gene_exp = gene_exp_mat,
                  sample_name = sample_name,
                  program_name = program,
                  program_genes = program_genes,
                  control_genes = control_genes
                ))
  )
}

#### Set up samples ------------------------------------------------------------

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

#### Read in and prepare gene expression ---------------------------------------

# Read in SingleCellExperiment objects and get relative gene expression
gene_exp_list <- sce_files |>
  purrr::map(readr::read_rds) |>
  purrr::map(\(sce) process_gene_expression(sce, platform = opt$platform))

#### Read in program and control genes -----------------------------------------

metaprogram_df <- readr::read_tsv(opt$metaprogram_file)
control_genes_df <- readr::read_tsv(opt$controls_file)

# The metaprogram names will be a vector of unique values in the metaprogram
# column
metaprograms <- unique(metaprogram_df$metaprogram)

#### Calculate scores ----------------------------------------------------------

scores_df <- metaprograms |>
  # For each metaprogram, calculate scores for all samples
  purrr::map(\(program) metaprogram_wrapper(program)) |>
  # Bring list to one level
  purrr::flatten() |>
  # Bind rows together -- all data frames have the same format and contain
  # the metaprogram name already
  dplyr::bind_rows()

# Remove rownames
rownames(scores_df) <- NULL

# Write to output file
readr::write_tsv(scores_df, file = opt$output_file)
