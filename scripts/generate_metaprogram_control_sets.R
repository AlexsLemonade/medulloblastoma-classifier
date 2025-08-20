# J. Taroni 2025
#
# For a set of metaprograms in a single-cell dataset, generate a control
# list of genes
#
# For each gene in the metaprogram, find the 100 genes outside of the
# metaprogram with the closest average pseudobulk expression across the whole
# cohort

library(optparse)

#### Command line options ------------------------------------------------------

option_list <- list(
  make_option(
    opt_str = c("-i", "--pseudobulk_input_file"),
    type = "character",
    help = "File path for TSV file that contains the pseudobulk expression data"
  ),
  make_option(
    opt_str = c("-m", "--metaprogram_file"),
    type = "character",
    help = "File path for TSV file that contains metaprograms",
    default = "processed_data/hovestadt-et-al-group-3-4-metaprogram-genes.tsv"
  ),
  make_option(
    opt_str = c("-p", "--platform"),
    type = "character",
    help = "Which type of scRNA-seq data is it? 10X or Smart-seq2"
  ),
  make_option(
    opt_str = c("-o", "--output_file"),
    type = "character",
    help = "File path for TSV file that will be output with the results"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Error handling for platform
stopifnot("Unsupported platform!" = (tolower(opt$platform)
                                     %in% c("10x", "smart-seq2")))

# Create output directory if it doesn't exist
dir.create(dirname(opt$output_file), showWarnings = FALSE, recursive = TRUE)

#### Functions -----------------------------------------------------------------

source(here::here("utils/convert_gene_names.R"))

pick_hundred_genes <- function(gene,
                               program_genes,
                               gene_means) {
  # Pick 100 control genes with the most similar average expression level to
  # the query gene
  #
  # Input:
  #   gene: A gene symbol of a gene to be used as a query
  #   program_genes: A vector of gene symbols in the program -- these will
  #                  be removed from consideration
  #   gene_means: A named vector of gene mean values where the names are
  #               gene symbols
  #
  # Output: A vector of 100 genes (not in the program) that are closest in their
  #         average expression level to the query gene

  # If you can't find the query gene, return NULL

  if (!(gene %in% names(gene_means))) {
    return(NULL)
  }

  # Subtract gene of interest average expression from all average expression
  gene_diff <- gene_means - gene_means[gene]

  # Remove genes that are in the program
  gene_diff <- gene_diff[setdiff(names(gene_means), program_genes)]

  # Closest 100 genes -- use the absolute value of the difference from the
  # query gene average expression level
  hundred_genes <- names(sort(abs(gene_diff), decreasing = FALSE))[1:100]

  # Return closest 100 genes
  return(hundred_genes)

}

#### Read in data --------------------------------------------------------------

pseudobulk_df <- readr::read_tsv(opt$pseudobulk_input_file)
metaprogram_df <- readr::read_tsv(opt$metaprogram_file)

#### Convert to gene symbols ---------------------------------------------------

pseudobulk_df <- convert_gene_names(genex_df = pseudobulk_df,
                                    gene_column_before = "gene",
                                    gene_column_after = "gene",
                                    map_from = "ENSEMBL",
                                    map_to = "SYMBOL") |>
  tibble::column_to_rownames("gene")

#### Set up metaprogram data ---------------------------------------------------

metaprogram_genes_list <- unique(metaprogram_df$metaprogram) |>
  purrr::map(\(program)
              metaprogram_df |>
                dplyr::filter(metaprogram == program) |>
                dplyr::pull(gene)
             ) |>
  purrr::set_names(unique(metaprogram_df$metaprogram))

#### Matrices and vectors ------------------------------------------------------

# If TPM, just log2 transform
if (tolower(opt$platform) == "smart-seq2") {
  pseudobulk_mat <- log2(as.matrix(pseudobulk_df) + 1)
} else if (tolower(opt$platform) == "10x") {
  # If 10X counts, convert counts to CPM first
  pseudobulk_mat <- edgeR::cpm(as.matrix(pseudobulk_df),
                               normalized.lib.sizes = FALSE)
  pseudobulk_mat <- log2(pseudobulk_mat + 1)
}

# Calculate gene_means
gene_means <- rowMeans(pseudobulk_mat)

#### Get control gene sets -----------------------------------------------------

metaprogram_controls_df <- metaprogram_genes_list |>
  purrr::map(\(program) program |>
               purrr::map(\(gene)
                          pick_hundred_genes(gene = gene,
                                             program_genes = program,
                                             gene_means = gene_means)
               ) |>
               purrr::flatten_chr() |>
               unique()
  ) |>
  purrr::imap(\(x, y) data.frame(gene = x, metaprogram = y)) |>
  dplyr::bind_rows()

#### Write to file -------------------------------------------------------------

readr::write_tsv(metaprogram_controls_df, file = opt$output_file)
