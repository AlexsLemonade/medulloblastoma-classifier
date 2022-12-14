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

# read in individual expression files
expression_df_list <-
  purrr::map(expression_files[-c(4:5, 25, 28, 33)],
             function(x)
               readr::read_tsv(x, skip = 1, col_names = FALSE) %>%
               tibble::column_to_rownames(var = "X1"))
