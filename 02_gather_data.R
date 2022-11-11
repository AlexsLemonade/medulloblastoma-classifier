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
                                              "bulk_metadata.tsv")
GSE164677_genex_input_filename <- file.path(data_dir,
                                            "GSE164677/GSE164677_Asian_MB_RNA-seq.txt.gz")
OpenPBTA_polya_genex_input_filename <- file.path(data_dir,
                                                 "OpenPBTA/pbta-gene-expression-rsem-tpm.polya.rds")
OpenPBTA_stranded_genex_input_filename <- file.path(data_dir,
                                                    "OpenPBTA/pbta-gene-expression-rsem-tpm.stranded.rds")

genex_df_output_filename <- file.path(processed_data_dir,
                                      "bulk_genex.tsv")

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
# Read in metadata
################################################################################

combined_metadata <- read_tsv(combined_metadata_input_filepath,
                              col_types = "c")

################################################################################
# create a list containing each data set
################################################################################

### experiments part of GSE124814

GSE124184_experiment_accession_ids <- read_tsv(GSE124184_experiment_accessions_input_filepath,
                                               col_names = "experiment_accession",
                                               col_types = "c") %>%
  filter(experiment_accession != "E-MTAB-292")


genex_data_list <- purrr::map(GSE124184_experiment_accession_ids$experiment_accession,
                              function(x) get_genex_data(file.path(data_dir, x, x, str_c(x, ".tsv", sep = "")),
                                                         combined_metadata$sample_accession))

names(genex_data_list) <- GSE124184_experiment_accession_ids$experiment_accession

### GSE164677

genex_data_list[["GSE164677"]] <- read_tsv(GSE164677_genex_input_filename,
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
  column_to_rownames(var = "GENEID")

### OpenPBTA

genex_data_list[["OpenPBTA"]] <- bind_cols(read_rds(OpenPBTA_polya_genex_input_filename),
                                           read_rds(OpenPBTA_stranded_genex_input_filename)[,-1]) %>%
  mutate(gene_id = stringr::str_split(gene_id, pattern = "\\.", simplify = TRUE)[,1]) %>%
  filter(!duplicated(gene_id)) %>%
  column_to_rownames(var = "gene_id") %>%
  select(combined_metadata %>%
           filter(study == "OpenPBTA") %>%
           pull(sample_accession))

### St. Jude

genex_data_list[["St. Jude"]] <- combined_metadata %>%
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

################################################################################
# combine the list
################################################################################

### common genes
common_genes <- get_common_genes(genex_data_list)

### column bind the studies together using common genes
lapply(genex_data_list,
       function(x) x[common_genes,]) %>%
  bind_cols() %>%
  rownames_to_column(var = "gene") %>%
  write_tsv(genex_df_output_filename)
