# ESTIMATE tumor purity
#
# Steven Foltz
# November 2022

suppressMessages(library(magrittr))
suppressMessages(library(estimate))

processed_data_dir <- here::here("processed_data")
estimate_data_dir <- file.path(processed_data_dir, "estimate")

genex_input_filename <- file.path(processed_data_dir, "bulk_genex.tsv")
metadata_input_filename <- file.path(processed_data_dir, "bulk_metadata.tsv")

estimate_input_filename <- file.path(estimate_data_dir, "estimate_input.tsv")
estimate_input_gct_filename <- file.path(estimate_data_dir, "estimate_input.gct")

estimate_results_output_gct_filename <- file.path(estimate_data_dir, "estimate_results.gct")
estimate_results_output_filename <- file.path(estimate_data_dir, "estimate_results.tsv")

if (!file.exists(estimate_input_filename)) {
  
  # read in data
  genex_df <- readr::read_tsv(genex_input_filename) %>%
    tibble::column_to_rownames(var = "gene")
  
  metadata_df <- readr::read_tsv(metadata_input_filename) %>%
    dplyr::filter(sample_accession %in% names(genex_df))
  
  # convert gene names
  AnnotationHub::setAnnotationHubOption("ASK", FALSE) # download without asking
  ah <- AnnotationHub::AnnotationHub()
  AnnotationHub::snapshotDate(ah) <- "2022-10-26" # reproducibility
  hs_orgdb <- AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]
  map_ENSEMBL_SYMBOL_dedup_df <- AnnotationDbi::select(x = hs_orgdb,
                                                       keys = AnnotationDbi::keys(hs_orgdb, "ENSEMBL"),
                                                       columns = "SYMBOL",
                                                       keytype = "ENSEMBL") %>%
    dplyr::mutate(dup_ENSEMBL = duplicated(ENSEMBL),
                  dup_SYMBOL = duplicated(SYMBOL)) %>%
    dplyr::filter(!dup_ENSEMBL, !dup_SYMBOL) %>%
    dplyr::select(ENSEMBL, SYMBOL)
  
  genex_df_SYMBOL <- genex_df %>%
    tibble::rownames_to_column(var = "ENSEMBL") %>%
    dplyr::left_join(map_ENSEMBL_SYMBOL_dedup_df,
                     by = "ENSEMBL") %>%
    dplyr::filter(!duplicated(ENSEMBL),
                  !duplicated(SYMBOL),
                  !is.na(SYMBOL)) %>%
    dplyr::select(-ENSEMBL) %>%
    dplyr::select(SYMBOL, tidyselect::everything()) %>%
    readr::write_tsv(file = estimate_input_filename)  
  
}

# convert gene expression tsv to GCT file
if (!file.exists(estimate_input_gct_filename)) {
  
  estimate::filterCommonGenes(input.f = estimate_input_filename,
                              output.f = estimate_input_gct_filename,
                              id = "GeneSymbol") # "EntrezID" did not work
  
}

# run ESTIMATE
estimate::estimateScore(input.ds = estimate_input_gct_filename,
                        output.ds = estimate_results_output_gct_filename)

estimate_df <- readr::read_tsv(estimate_results_output_gct_filename,
                               skip = 2) %>% # skip two metadata lines
  dplyr::select(-Description) %>% # Description column is redundant
  tidyr::pivot_longer(cols = -c(NAME),
                      names_to = "sample_accession") %>%
  tidyr::pivot_wider(names_from = NAME,
                     values_from = value) %>%
  readr::write_tsv(file = estimate_results_output_filename)

# plot tumor purity violins
estimate_df %>%
  dplyr::left_join(metadata_df,
                   by = "sample_accession") %>%
  dplyr::filter(!is.na(subgroup)) %>%
  ggplot2::ggplot(ggplot2::aes(x = subgroup,
                               y = TumorPurity,
                               color = subgroup)) +
  ggplot2::geom_violin(show.legend = FALSE) +
  ggplot2::geom_jitter(shape = 16,
                       height = 0,
                       width = 0.1,
                       show.legend = FALSE) +
  ggplot2::facet_wrap(~ platform) +
  ggplot2::theme_bw()

# plot stromal vs. immune scores
estimate_df %>%
  dplyr::left_join(metadata_df,
                   by = "sample_accession") %>%
  dplyr::filter(!is.na(subgroup)) %>%
  ggplot2::ggplot(ggplot2::aes(x = StromalScore,
                               y = ImmuneScore,
                               color = subgroup)) +
  ggplot2::geom_point(shape = 16) +
  ggplot2::facet_wrap(~ platform) +
  ggplot2::theme_bw()

# plot estimate score vs. tumor purity
estimate_df %>%
  dplyr::left_join(metadata_df,
                   by = "sample_accession") %>%
  dplyr::filter(!is.na(subgroup)) %>%
  ggplot2::ggplot(ggplot2::aes(x = ESTIMATEScore,
                               y = TumorPurity,
                               color = subgroup)) +
  ggplot2::geom_point(shape = 16) +
  ggplot2::facet_wrap(~ platform) +
  ggplot2::theme_bw()
