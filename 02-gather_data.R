# Gather gene expression data from each data source
#
# Steven Foltz
# November 2022

annotationhub_snapshot_date <- "2022-10-26" # reproducibility

suppressMessages(library(tidyverse))

data_dir <- here::here("data")
processed_data_dir <- here::here("processed_data")

GSE124184_experiment_accessions_input_filepath <- file.path(data_dir,
                                                            "GSE124814_experiment_accessions.tsv")
bulk_metadata_input_filepath <- file.path(processed_data_dir,
                                          "bulk_metadata.tsv")
GSE164677_genex_input_filename <- file.path(data_dir, "GSE164677",
                                            "GSE164677_Asian_MB_RNA-seq.txt.gz")
OpenPBTA_polya_genex_input_filename <- file.path(data_dir, "OpenPBTA",
                                                 "pbta-gene-expression-rsem-tpm.polya.rds")
OpenPBTA_stranded_genex_input_filename <- file.path(data_dir, "OpenPBTA",
                                                    "pbta-gene-expression-rsem-tpm.stranded.rds")

genex_df_output_filename <- file.path(processed_data_dir,
                                      "bulk_genex.tsv")

################################################################################
# set up gene name conversions
################################################################################

AnnotationHub::setAnnotationHubOption("ASK", FALSE) # download without asking
ah <- AnnotationHub::AnnotationHub()
AnnotationHub::snapshotDate(ah) <- annotationhub_snapshot_date
hs_orgdb <- AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]
map_ensembl_symbol_dedup_df <- AnnotationDbi::select(x = hs_orgdb,
                                                     keys = AnnotationDbi::keys(hs_orgdb, "ENSEMBL"),
                                                     columns = "SYMBOL",
                                                     keytype = "ENSEMBL") %>%
  dplyr::mutate(dup_ensembl = duplicated(ENSEMBL),
                dup_symbol = duplicated(SYMBOL)) %>%
  dplyr::filter(!dup_ensembl, !dup_symbol) %>%
  dplyr::select(ENSEMBL, SYMBOL)

################################################################################
# functions
################################################################################

get_genex_data <- function(genex_file_path,
                           mb_sample_accessions){
  
  # for a gene expression file located at genex_file_path,
  # read in columns associated with sample accessions in mb_sample_accessions
  
  genex_df_columns <- readr::read_tsv(genex_file_path,
                                      col_types = "c",
                                      n_max = 0)
  
  select_these_samples_TF <- names(genex_df_columns)[-1] %in% mb_sample_accessions
  
  select_these_columns_types <- stringr::str_c(c("c", # for the gene column
                                                 ifelse(select_these_samples_TF,
                                                        "d",
                                                        "-")),
                                               collapse = "")
  
  genex_df <- readr::read_tsv(genex_file_path,
                              col_types = select_these_columns_types) %>%
    tibble::column_to_rownames(var = "Gene")
  
  return(genex_df)
  
}

get_common_genes <- function(genex_list){
  
  # for a list of gene expression files with row names as genes,
  # find the common row names (genes)
  
  Reduce(intersect, lapply(genex_list, row.names))
  
}

################################################################################
# Read in metadata
################################################################################

bulk_metadata <- readr::read_tsv(bulk_metadata_input_filepath,
                                 col_types = "c")

################################################################################
# create a list containing each data set
################################################################################

### experiments part of GSE124814

GSE124184_experiment_accession_ids <- readr::read_tsv(GSE124184_experiment_accessions_input_filepath,
                                                      col_names = "experiment_accession",
                                                      col_types = "c") %>%
  dplyr::filter(experiment_accession != "E-MTAB-292")


genex_data_list <- purrr::map(GSE124184_experiment_accession_ids$experiment_accession,
                              function(x) get_genex_data(file.path(data_dir, x, x, stringr::str_c(x, ".tsv", sep = "")),
                                                         bulk_metadata$sample_accession))

names(genex_data_list) <- GSE124184_experiment_accession_ids$experiment_accession

### GSE164677

genex_data_list[["GSE164677"]] <- readr::read_tsv(GSE164677_genex_input_filename,
                                                  col_names = TRUE,
                                                  show_col_types = FALSE,
                                                  skip = 1) %>%
  dplyr::left_join(map_ensembl_symbol_dedup_df,
                   by = c("gene" = "SYMBOL")) %>%
  dplyr::filter(!duplicated(gene),
                !duplicated(ENSEMBL),
                !is.na(ENSEMBL)) %>%
  dplyr::select(-gene) %>%
  tibble::column_to_rownames(var = "ENSEMBL")

### OpenPBTA

genex_data_list[["OpenPBTA"]] <- dplyr::bind_cols(readr::read_rds(OpenPBTA_polya_genex_input_filename),
                                                  readr::read_rds(OpenPBTA_stranded_genex_input_filename)[,-1]) %>%
  dplyr::mutate(gene_id = stringr::str_split(gene_id, pattern = "\\.", simplify = TRUE)[,1]) %>%
  dplyr::filter(!duplicated(gene_id)) %>%
  tibble::column_to_rownames(var = "gene_id") %>%
  dplyr::select(bulk_metadata %>%
                  dplyr::filter(study == "OpenPBTA") %>%
                  dplyr::pull(sample_accession))

### St. Jude

genex_data_list[["St. Jude"]] <- bulk_metadata %>%
  dplyr::filter(study == "St. Jude") %>%
  dplyr::pull(sample_accession) %>%
  purrr::map(function(x) readr::read_tsv(file.path(data_dir, "stjudecloud",
                                                   stringr::str_c(x, ".RNA-Seq.feature-counts.txt")),
                                         col_names = c("SYMBOL", x), show_col_types = FALSE) %>%
               tibble::column_to_rownames("SYMBOL")) %>%
  dplyr::bind_cols() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::left_join(map_ensembl_symbol_dedup_df,
                   by = "SYMBOL") %>%
  dplyr::filter(!duplicated(SYMBOL),
                !duplicated(ENSEMBL),
                !is.na(ENSEMBL)) %>%
  dplyr::select(-SYMBOL) %>%
  tibble::column_to_rownames(var = "ENSEMBL")

################################################################################
# combine the list
################################################################################

### common genes
common_genes <- get_common_genes(genex_data_list)

### column bind the studies together using common genes
lapply(genex_data_list,
       function(x) x[common_genes,]) %>%
  dplyr::bind_cols() %>%
  tibble::rownames_to_column(var = "gene") %>%
  readr::write_tsv(genex_df_output_filename)
