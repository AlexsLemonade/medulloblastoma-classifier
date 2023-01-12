# Read in, clean the pseudo-bulk expression files, and generate the pseudo-bulk
# matrix.
#
# Chante Bethell
# December 2022

suppressMessages(library(tidyverse))

data_dir <- here::here("data")
processed_data_dir <- here::here("processed_data")

# define input and output filepaths
input_metadata_filepath <- file.path(processed_data_dir, "pseudobulk_metadata.tsv")
output_matrix_filepath <- file.path(processed_data_dir, "pseudobulk_genex.tsv")

# read in pseudo-bulk metadata file
pseudobulk_metadata <- readr::read_tsv(input_metadata_filepath)

# grab the names of the individual expression files
expression_files <- file.path(data_dir, 
                              "GSE119926", 
                              str_c(pseudobulk_metadata$sample_accession, 
                                    "_", pseudobulk_metadata$title, ".txt.gz"))
names(expression_files) <- pseudobulk_metadata$title

# read in individual expression files
expression_df_list <- purrr::map(expression_files,
                                 function(x)
                                   readr::read_tsv(x, skip = 1, col_names = FALSE) %>%
                                   tibble::column_to_rownames(var = "X1"))

# revert log transformed-TPM values back to original TPM values using equation 
# 10*(2^x - 1) -- determined by working back from the equation log2(TPM/10+1)
# found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3905406
tpm_df_list <- lapply(expression_df_list, function(x) 10*(2^x - 1))

# average the TPM values across cells for each data frame
average_tpm_list <- lapply(tpm_df_list, rowMeans)

# combine the list of TPM values into a single matrix
pseudobulk_gene_names <- names(average_tpm_list[[1]])
pseudobulk_matrix <- dplyr::bind_cols(gene = pseudobulk_gene_names,
                                      average_tpm_list)

# save matrix object
readr::write_tsv(pseudobulk_matrix, output_matrix_filepath)
