# 1. Gather bulk gene expression data from each data source
# 2. Read in, clean the pseudo-bulk expression files, and generate the pseudo-bulk
# matrix.
#
# Chante Bethell, Steven Foltz
# November-December 2022

data_dir <- here::here("data")
processed_data_dir <- here::here("processed_data")
pseudobulk_data_dir <- here::here(processed_data_dir, "pseudobulk")
GSE119926_pseudobulk_dir <- here::here(pseudobulk_data_dir, "GSE119926")
GSE119926_pseudobulk_sce_dir <- here::here(GSE119926_pseudobulk_dir, "pseudobulk_sce")
utils_dir <- here::here("utils")

source(here::here(utils_dir, "convert_gene_names.R"))
source(here::here(utils_dir, "single-cell.R"))
source(here::here(utils_dir, "TPM_conversion.R"))

purrr::map(c(GSE119926_pseudobulk_sce_dir),
           function(dir) dir.create(dir,
                                    showWarnings = FALSE,
                                    recursive = TRUE))

################################################################################
# set input and output filepaths
################################################################################

# metadata inputs
bulk_metadata_input_filepath <- here::here(processed_data_dir,
                                           "bulk_metadata.tsv")
pseudobulk_metadata_input_filepath <- here::here(pseudobulk_data_dir,
                                                 "pseudobulk_metadata.tsv")

# bulk genex inputs
GSE124184_experiment_accessions_input_filepath <- here::here(data_dir,
                                                             "GSE124814_experiment_accessions.tsv")
GSE164677_genex_input_filepath <- here::here(data_dir, "GSE164677",
                                             "GSE164677_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
OpenPBTA_polya_genex_input_filepath <- here::here(data_dir, "OpenPBTA",
                                                  "pbta-gene-expression-rsem-tpm.polya.rds")
OpenPBTA_stranded_genex_input_filepath <- here::here(data_dir, "OpenPBTA",
                                                     "pbta-gene-expression-rsem-tpm.stranded.rds")

# gene map input
gene_map_input_filepath <- here::here(processed_data_dir,
                                      "gene_map.tsv")

# GENCODE annotations for St. Jude RNA-seq counts to TPM conversion
gencode_annotation_gtf_filepath <- here::here(data_dir, "GENCODE",
                                              "gencode.v31.annotation.gtf.gz")
GENCODE_gene_lengths_filepath <- here::here(data_dir, "GENCODE",
                                            "GENCODE_gene_lengths.tsv")

# genex outputs
bulk_genex_df_output_filepath <- here::here(processed_data_dir,
                                            "bulk_genex.tsv")
GSE119926_pseudobulk_genex_df_output_filepath <- here::here(processed_data_dir,
                                                            "GSE119926_pseudobulk_genex.tsv")

################################################################################
# functions
################################################################################

get_genex_data <- function(genex_filepath,
                           mb_sample_accessions){

  # for a gene expression file located at genex_filepath,
  # read in columns associated with sample accessions in the mb_sample_accessions vector

  genex_df_columns <- readr::read_tsv(genex_filepath,
                                      col_types = "c",
                                      n_max = 0)

  # we expect the structure of gene expression files read in this way to be:
  # genes (rows) x samples (columns)
  # gene names are kept in column 1 with column header "Gene"
  # sample names are found in columns 2-N
  if (names(genex_df_columns)[1] != "Gene") {

    stop("First column name of gene expression data should be 'Gene' in get_genex_data().")

  }

  select_these_samples_TF <- names(genex_df_columns)[-1] %in% mb_sample_accessions

  select_these_columns_types <- stringr::str_c(c("c", # for the gene column
                                                 ifelse(select_these_samples_TF,
                                                        "d",
                                                        "-")),
                                               collapse = "")

  genex_df <- readr::read_tsv(genex_filepath,
                              col_types = select_these_columns_types) |>
    dplyr::filter(!duplicated(Gene)) |>
    tibble::column_to_rownames(var = "Gene")

  return(genex_df)

}

################################################################################
# Read in metadata and gene map
################################################################################

bulk_metadata <- readr::read_tsv(bulk_metadata_input_filepath,
                                 col_types = "c")

################################################################################
# create a list containing each data set
################################################################################

### experiments part of GSE124814

GSE124184_experiment_accession_ids <- readr::read_tsv(GSE124184_experiment_accessions_input_filepath,
                                                      col_names = "experiment_accession",
                                                      col_types = "c") |>
  dplyr::filter(experiment_accession != "E-MTAB-292")


genex_data_list <- purrr::map(GSE124184_experiment_accession_ids$experiment_accession,
                              function(x) get_genex_data(here::here(data_dir, x, x, stringr::str_c(x, ".tsv", sep = "")),
                                                         bulk_metadata$sample_accession))

names(genex_data_list) <- GSE124184_experiment_accession_ids$experiment_accession

### GSE164677

genex_data_list[["GSE164677"]] <- readr::read_tsv(GSE164677_genex_input_filepath,
                                                  col_names = TRUE,
                                                  show_col_types = FALSE) |>
  dplyr::mutate(GeneID = as.character(GeneID)) |>
  dplyr::filter(!duplicated(GeneID)) |>
  convert_gene_names(gene_column_before = "GeneID",
                     gene_column_after = "gene",
                     map_from = "ENTREZID",
                     map_to = "ENSEMBL") |>
  tibble::column_to_rownames(var = "gene")

### OpenPBTA

genex_data_list[["OpenPBTA"]] <- dplyr::bind_cols(readr::read_rds(OpenPBTA_polya_genex_input_filepath),
                                                  readr::read_rds(OpenPBTA_stranded_genex_input_filepath)[,-1]) |>
  dplyr::mutate(gene_id = stringr::str_split(gene_id, pattern = "\\.", simplify = TRUE)[,1]) |>
  dplyr::filter(!duplicated(gene_id)) |>
  tibble::column_to_rownames(var = "gene_id") |>
  dplyr::select(bulk_metadata |>
                  dplyr::filter(study == "OpenPBTA") |>
                  dplyr::pull(sample_accession))

### St. Jude

GENCODE_gene_lengths_df <- get_GENCODE_gene_lengths(gtf_filepath = gencode_annotation_gtf_filepath,
                                                    GENCODE_gene_lengths_filepath = GENCODE_gene_lengths_filepath)

genex_data_list[["St. Jude"]] <- bulk_metadata |>
  dplyr::filter(study == "St. Jude") |>
  dplyr::pull(sample_accession) |>
  purrr::map(function(x) readr::read_tsv(here::here(data_dir, "stjudecloud",
                                                    stringr::str_c(x, ".RNA-Seq.feature-counts.txt")),
                                         col_names = c("SYMBOL", x), show_col_types = FALSE) |>
               tibble::column_to_rownames("SYMBOL")) |>
  dplyr::bind_cols() |>
  tibble::rownames_to_column("SYMBOL") |>
  dplyr::filter(!duplicated(SYMBOL)) |>
  convert_gene_names(gene_column_before = "SYMBOL",
                     gene_column_after = "gene",
                     map_from = "SYMBOL",
                     map_to = "ENSEMBL") |>
  tibble::column_to_rownames(var = "gene") |>
  convert_gene_counts_to_TPM(gene_lengths_df = GENCODE_gene_lengths_df)

################################################################################
# combine the list
################################################################################

### common genes
common_genes <- genex_data_list |>
  purrr::map(row.names) |>
  purrr::reduce(intersect)

### column bind the studies together using common genes
lapply(genex_data_list,
       function(x) x[common_genes,]) |>
  dplyr::bind_cols() |>
  tibble::rownames_to_column(var = "gene") |>
  readr::write_tsv(bulk_genex_df_output_filepath)

################################################################################
# Pseudobulk - generate pseudobulk data matrix from individual scRNA-seq data
################################################################################

# read in pseudo-bulk metadata file
pseudobulk_metadata <- readr::read_tsv(pseudobulk_metadata_input_filepath,
                                       show_col_types = FALSE)

# GSE119926, which has separate files for each sample

GSE119926_pseudobulk_metadata <- pseudobulk_metadata |>
  dplyr::filter(study == "GSE119926")

# grab the names of the individual expression files
GSE119926_sample_accession_ids <- GSE119926_pseudobulk_metadata$sample_accession
GSE119926_sample_titles <- GSE119926_pseudobulk_metadata$title

GSE119926_pseudobulk_df <- NULL

for (sample_iter in seq_along(GSE119926_sample_accession_ids)) {

  print(sample_iter)

  sample_acc <- GSE119926_sample_accession_ids[sample_iter]
  sample_title <- GSE119926_sample_titles[sample_iter]

  pseudobulk_expression_filepath <- here::here(data_dir,
                                               "GSE119926",
                                               stringr::str_c(sample_acc, "_", sample_title, ".txt.gz"))

  scrna_genex_df <- readr::read_tsv(pseudobulk_expression_filepath,
                                    skip = 1,
                                    col_names = FALSE,
                                    show_col_types = FALSE) |>
    dplyr::rename("gene" = "X1") |>
    dplyr::filter(!duplicated(gene)) |>
    convert_gene_names(gene_column_before = "gene",
                       gene_column_after = "gene",
                       map_from = "SYMBOL",
                       map_to = "ENSEMBL") |>
    tibble::column_to_rownames(var = "gene")

  # revert log transformed-TPM values back to original TPM values using equation
  # 10*(2^x - 1) -- determined by working back from the equation log2(TPM/10+1)
  # found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3905406
  tpm_genex_df <- 10*(2^scrna_genex_df - 1)

  # average the TPM values across cells
  average_tpm_vector <- rowMeans(tpm_genex_df)

  # combine the list of TPM values into a single matrix
  # if pseudobulk_df does not contain data yet (i.e. if sample_iter is 1),
  # then create a gene column and sample column for the first sample,
  # otherwise, add a column for each additional sample
  if (is.null(GSE119926_pseudobulk_df)) {

    GSE119926_pseudobulk_df <- tibble::tibble(gene = rownames(tpm_genex_df),
                                    "{sample_title}" := average_tpm_vector)

  } else {

    GSE119926_pseudobulk_df <- GSE119926_pseudobulk_df |>
      tibble::add_column("{sample_title}" := average_tpm_vector)

  }

  # define output file names for SCE objects
  sce_output_filepath = here::here(GSE119926_pseudobulk_sce_dir,
                                   stringr::str_c(sample_title, "_sce.rds"))

  # convert TPM matrix to SingleCellExperiment objects
  SingleCellExperiment::SingleCellExperiment(assays = list(counts = tpm_genex_df)) |>
    # calculate UMAP results
    add_sce_umap() |>
    # perform clustering
    perform_graph_clustering() |>
    # write to file
    readr::write_rds(file = sce_output_filepath)

}

# save pseudobulk df object
readr::write_tsv(GSE119926_pseudobulk_df,
                 GSE119926_pseudobulk_genex_df_output_filepath)
