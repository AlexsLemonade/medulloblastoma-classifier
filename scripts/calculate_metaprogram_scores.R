library(optparse)

#### Command line options ------------------------------------------------------

option_list <- list(
  make_option(
    opt_str = c("-i", "--sce_input_dir"),
    type = "character",
    help = "Directory that contains RDS of SingleCellExperiment objects"
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

calculate_metaprogram_score <- function(sce_filepath,
                                        sample_name,
                                        platform,
                                        program_name,
                                        program_genes,
                                        control_genes) {


  # Read in SingleCellExperiment object
  sce <- readr::read_rds(sce_filepath)

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
    purrr::imap(sce_files,
                \(sce_file, sample_name)
                calculate_metaprogram_score(
                  sce_filepath = sce_file,
                  sample_name = sample_name,
                  platform = opt$platform,
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

#### Read in program and control genes -----------------------------------------

metaprogram_df <- readr::read_tsv(opt$metaprogram_file)
control_genes_df <- readr::read_tsv(opt$controls_file)

# The metaprogram names will be a vector of unique values in the metaprogram
# column
metaprograms <- unique(metaprogram_df$metaprogram)

#### Calculate scores ----------------------------------------------------------

scores_list <- metaprograms |>
  purrr::map(\(program) metaprogram_wrapper(program))
