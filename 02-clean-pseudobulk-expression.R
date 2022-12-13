# Read in and clean the pseudo-bulk expression files to prepare for the 
# pseudo-bulk matrix generation.
#
# Chante Bethell
# December 2022

processed_data_dir <- here::here("processed_data")

# read in pseudo-bulk metadata file
pseudobulk_metadata <- readr::read_tsv(file.path(processed_data_dir,
                                                 "pseudobulk_metadata.tsv"))

# grab the names of the individual expression files
expression_files <- list.files(file.path("data", "GSE119926"), "GSM")

# read in individual expression files
expression_df_list <- list()
for (file in expression_files[-c(4:5, 25, 28, 33)]) {
  file_df <-  data.frame(readr::read_tsv(file.path("data", "GSE119926", file))) %>%
    distinct()
  
  # set gene names as rownames  
  rownames(file_df) <- file_df[,1]
  file_df <- file_df[,-1]
  
  # add to list of expression data frames
  expression_df_list[[file]] <- file_df 
}
