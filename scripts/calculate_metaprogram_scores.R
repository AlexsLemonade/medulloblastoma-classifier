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
    opt_str = c("-o", "--output_file"),
    type = "character",
    help = "File path for TSV file that contains metaprogram scores"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Create output directory if it doesn't exist
dir.create(dirname(opt$output_file), showWarnings = FALSE, recursive = TRUE)


#### Function ------------------------------------------------------------------

calculate_metaprogram_score <- function(sce,
                                        sample_name,
                                        program_name,
                                        program_genes,
                                        control_genes) {

  # Extract the data in the counts slot
  gene_exp <- SingleCellExperiment::counts(sce)

  # Subset to program genes
  program_genes <- intersect(rownames(gene_exp), program_genes)
  program_gene_exp <- gene_exp[program_genes, ]

  # We don't have to worry about control genes not being in the object, since
  # they are dataset-specific (derived from pseudobulk data)
  control_gene_exp <- gene_exp[control_genes, ]

  # The score is the average expression of metaprogram genes minus the
  # average expression of the control genes
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

#### Set up samples ------------------------------------------------------------

sce_files <- fs::dir_ls(path = fs::path(opt$sce_input_dir, "sce"),
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


