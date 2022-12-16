# Read in and clean the pseudo-bulk expression files to prepare for the 
# pseudo-bulk matrix generation.
#
# Chante Bethell
# December 2022

suppressMessages(library(tidyverse))

data_dir <- here::here("data")
processed_data_dir <- here::here("processed_data")

# read in pseudo-bulk metadata file
pseudobulk_metadata <- readr::read_tsv(file.path(processed_data_dir,
                                                 "pseudobulk_metadata.tsv"))

# grab the names of the individual expression files
expression_files <-
  list.files(file.path(data_dir, "GSE119926"), "GSM", full.names = TRUE)
names(expression_files) <-
  pseudobulk_metadata$title[stringr::str_detect(expression_files, pseudobulk_metadata$title)]

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
average_tpm_list <- lapply(expression_df_list, rowMeans)
