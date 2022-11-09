# Gather gene expression data from each data source
#
# Steven Foltz
# November 2022

suppressMessages(library(tidyverse))

data_dir <- here::here("data")
processed_data_dir <- here::here("processed_data")

GSE124184_experiment_accessions_input_filepath <- file.path(data_dir,
                                                            "GSE124814_experiment_accessions.tsv")
combined_metadata_input_filepath <- file.path(processed_data_dir,
                                              "combined_metadata.tsv")

genex_df_output_filename <- file.path(processed_data,
                                      "genex_df.tsv")

################################################################################
# functions
################################################################################

get_genex_data <- function(genex_file_path,
                           mb_sample_accessions){
  
  genex_df_columns <- read_tsv(genex_file_path,
                               col_types = "c",
                               n_max = 0)
  
  select_these_samples_TF <- names(genex_df_columns)[-1] %in% mb_sample_accessions
  
  select_these_columns_types <- str_c(c("c", # for the gene column
                                        ifelse(select_these_samples_TF,
                                             "d",
                                             "-")),
                                      collapse = "")
  
  genex_df <- read_tsv(genex_file_path,
                       col_types = select_these_columns_types) %>%
    column_to_rownames(var = "Gene")
  
  return(genex_df)
  
}

get_common_genes <- function(genex_list){
  
  Reduce(intersect, lapply(genex_list, row.names))
  
}


################################################################################
# create a list containing each data set
################################################################################

# experiments part of GSE124814

experiment_accession_ids <- read_tsv(experiment_accessions_filepath,
                                     col_names = "experiment_accession",
                                     col_types = "c") %>%
  filter(experiment_accession != "E-MTAB-292")


genex_data_list <- purrr::map(experiment_accession_ids$experiment_accession,
                              function(x) get_genex_data(file.path(data_dir, x, x, str_c(x, ".tsv", sep = "")),
                                                         all_metadata$sample_accession))

names(genex_data_list) <- experiment_accession_ids$experiment_accession

### GSE164677

genex_data_list[["GSE164677"]] <- read_tsv("data/GSE164677/GSE164677_Asian_MB_RNA-seq.txt.gz",
                                           col_names = TRUE,
                                           show_col_types = FALSE,
                                           skip = 1) %>%
  left_join(ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                              keys = .$gene,
                              keytype = "SYMBOL",
                              columns = "GENEID"),
            by = c("gene" = "SYMBOL")) %>%
  filter(!duplicated(gene),
         !duplicated(GENEID),
         !is.na(GENEID)) %>%
  select(-gene) %>%
  column_to_rownames(var = "GENEID") %>%
  select(-ends_with("N")) # normals

### E-MTAB-292

### St. Jude
genex_data_list[["St. Jude"]] <- all_metadata %>%
  filter(study == "St. Jude") %>%
  pull(sample_accession) %>%
  purrr::map(function(x) read_tsv(file.path(data_dir,
                                            "stjudecloud",
                                            str_c(x,
                                                  ".RNA-Seq.feature-counts.txt")),
                                  col_names = c("SYMBOL", x), show_col_types = FALSE) %>%
               column_to_rownames("SYMBOL")) %>%
  bind_cols() %>%
  rownames_to_column("SYMBOL") %>%
  left_join(ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                              keys = .$SYMBOL,
                              keytype = "SYMBOL",
                              columns = "GENEID"),
            by = "SYMBOL") %>%
  filter(!duplicated(SYMBOL),
         !is.na(GENEID)) %>%
  select(-SYMBOL) %>%
  column_to_rownames(var = "GENEID")

# common genes
common_genes <- get_common_genes(genex_data_list)

################################################################################
# combine the list
################################################################################

genex_df <- lapply(genex_data_list,
                   function(x) x[common_genes,]) %>%
  bind_cols()

# set train/test
test_train_genex_list <- list()
test_train_metadata_list <- list()
set.seed(1)
for (i in 1:10) {
  
  test_train_metadata_list[[i]] <- list()
  test_train_genex_list[[i]] <- list()
  
  array_train <- all_metadata %>%
    filter(sample_accession %in% names(genex_df),
           platform == "Array") %>%
    select(study) %>%
    unique() %>%
    sample_n(size = 5) %>%
    pull(study)
    
  array_test <- all_metadata %>%
    filter(sample_accession %in% names(genex_df),
           platform == "Array") %>%
    select(study) %>%
    unique() %>%
    filter(!(study %in% array_train)) %>%
    pull(study)
  
  seq_train <- all_metadata %>%
    filter(sample_accession %in% names(genex_df),
           platform == "RNA-seq") %>%
    select(study) %>%
    unique() %>%
    sample_n(size = 1) %>%
    pull(study)
  
  seq_test <- all_metadata %>%
    filter(sample_accession %in% names(genex_df),
           platform == "RNA-seq") %>%
    select(study) %>%
    unique() %>%
    filter(!(study %in% seq_train)) %>%
    pull(study)
  
  train_samples_metadata <- all_metadata %>%
    filter(study %in% c(array_train, seq_train),
           sample_accession %in% names(genex_df))
  
  test_samples_metadata <- all_metadata %>%
    filter(study %in% c(array_test, seq_test),
           sample_accession %in% names(genex_df))
  
  test_train_metadata_list[[i]][["train"]] <- train_samples_metadata
  test_train_metadata_list[[i]][["test"]] <- test_samples_metadata
  
  test_train_genex_list[[i]][["train"]] <- genex_df %>%
    select(train_samples_metadata$sample_accession)
    
  test_train_genex_list[[i]][["test"]] <- genex_df %>%
    select(test_samples_metadata$sample_accession)
  
}


# genex and metadata by study

study_metadata_list <- all_metadata %>%
  filter(sample_accession %in% names(genex_df)) %>%
  arrange(study) %>%
  group_by(study) %>%
  group_split()

names(study_metadata_list) <- purrr::map(study_metadata_list, function(x) unique(x$study))


study_genex_list <- purrr::map(study_metadata_list,
                               function(x) genex_df %>%
                                 select(x$sample_accession))

names(study_genex_list) <- names(study_metadata_list)

