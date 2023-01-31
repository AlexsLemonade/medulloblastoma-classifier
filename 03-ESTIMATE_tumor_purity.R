# ESTIMATE tumor purity of bulk expression data
#
# Steven Foltz
# November 2022 - January 2023

suppressMessages(library(magrittr))
suppressMessages(library(estimate))

# set up directories and file names

processed_data_dir <- here::here("processed_data")
plot_dir <- here::here("plots")
estimate_data_dir <- file.path(processed_data_dir, "estimate")
estimate_plot_dir <- file.path(plot_dir, "estimate")

# bulk gene expression and metadata file names
genex_input_filename <- file.path(processed_data_dir, "bulk_genex.tsv")
metadata_input_filename <- file.path(processed_data_dir, "bulk_metadata.tsv")

# gene map for converting gene names
gene_map_df_filename <- file.path(processed_data_dir, "gene_map.tsv")

# files used by ESTIMATE functions
estimate_input_filename <- file.path(estimate_data_dir, "estimate_input.tsv")
estimate_input_gct_filename <- file.path(estimate_data_dir, "estimate_input.gct")

# files produced by ESTIMATE functions (.gct) and after some reshaping (.tsv)
estimate_results_output_gct_filename <- file.path(estimate_data_dir, "estimate_results.gct")
estimate_results_output_filename <- file.path(estimate_data_dir, "estimate_results.tsv")

# output plot file names
tumor_purity_plot_filename <- file.path(estimate_plot_dir, "tumor_purity_by_subgroup.pdf")
stromal_immune_plot_filename <- file.path(estimate_plot_dir, "stromal_score_by_immune_score.pdf")
estimate_score_tumor_purity_plot_filename <- file.path(estimate_plot_dir, "estimate_score_by_tumor_purity.pdf")

# read in bulk genex and metadata
genex_df <- readr::read_tsv(genex_input_filename,
                            show_col_types = FALSE) %>%
  tibble::column_to_rownames(var = "gene")

metadata_df <- readr::read_tsv(metadata_input_filename,
                               show_col_types = FALSE)

# reformat genex data by converting ENSEMBL gene IDs to SYMBOL

# read in gene map
# includes ENSEMBL, ENTREZID, and SYMBOL columns for gene name conversion
gene_map_df <- readr::read_tsv(gene_map_df_filename,
                               show_col_types = FALSE)

# convert gene names
genex_df_SYMBOL <- genex_df %>%
  tibble::rownames_to_column(var = "ENSEMBL") %>%
  dplyr::left_join(gene_map_df %>% # do not include ENTREZID here
                     dplyr::select(ENSEMBL, SYMBOL),
                   by = "ENSEMBL") %>%
  dplyr::filter(!duplicated(ENSEMBL),
                !duplicated(SYMBOL),
                !is.na(SYMBOL)) %>%
  dplyr::select(SYMBOL, tidyselect::everything(), -ENSEMBL)

readr::write_tsv(x = genex_df_SYMBOL,
                 file = estimate_input_filename)

# convert gene expression tsv to GCT file
# filterCommonGenes reduces the input to only include genes used in ESTIMATE model 
estimate::filterCommonGenes(input.f = estimate_input_filename,
                            output.f = estimate_input_gct_filename,
                            id = "GeneSymbol") # "EntrezID" did not work

################################################################################
# run ESTIMATE
################################################################################

# ESTIMATE publication: https://www.nature.com/articles/ncomms3612
# ESTIMATE generates Stromal and Immune scores using ssGSEA, regardless of platform
# Stromal + Immune scores = ESTIMATE score
# The model converting ESTIMATE scores to tumor purity was trained on Affymetrix data
# That model is: tumor purity = cos(0.6049872018 + 0.0001467884 * ESTIMATE)
# To return tumor purity scores, we must set platform = "affymetrix" (default)
#   even though some data are RNA-seq, otherwise the TP conversion is not done

estimate::estimateScore(input.ds = estimate_input_gct_filename,
                        output.ds = estimate_results_output_gct_filename,
                        platform = "affymetrix")

# reshape ESTIMATE GCT file into something tidy
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

# There are four variables (Stromal, Immune, ESTIMATE, and tumor purity)
#   associated with each observation (one sample) but GCT is an untidy format
#   and we need to skip the first two metadata lines (version and data dimensions)

# We need to convert from this:
# #1.2
# n_features       n_samples
# NAME             Description    sample1 sample2
# StromalScore     StromalScore   w1      w2
# ImmuneScore      ImmuneScore    x1      x2
# ESTIMATEScore    ESTIMATEScore  y1      y2
# TumorPurity      TumorPurity    z1      z2

# Into this (intermediate step):
# NAME             sample_accession value
# StromalScore     sample1          w1
# ImmuneScore      sample1          x1
# ESTIMATEScore    sample1          y1
# TumorPurity      sample1          z1
# StromalScore     sample2          w2
# ImmuneScore      sample2          x2
# ESTIMATEScore    sample2          y2
# TumorPurity      sample2          z2

# Into this (tidy format):
# sample_accession StromalScore ImmuneScore ESTIMATEScore TumorPurity
# sample1          w1           x1          y1            z1
# sample2          w2           x2          y2            z2

estimate_df <- readr::read_tsv(estimate_results_output_gct_filename,
                               show_col_types = FALSE,
                               skip = 2) %>% # skip two metadata lines
  dplyr::select(-Description) %>% # Description column is redundant
  tidyr::pivot_longer(cols = -NAME, #create a column of sample accessions
                      names_to = "sample_accession") %>%
  tidyr::pivot_wider(names_from = NAME, # collect variables from each sample
                     values_from = value) # into a single row

readr::write_tsv(x = estimate_df,
                 file = estimate_results_output_filename)

################################################################################
# visualize results
################################################################################

estimate_metadata_plot_df <- estimate_df %>%
  dplyr::left_join(metadata_df,
                   by = "sample_accession") %>%
  dplyr::filter(!is.na(subgroup))

# set a consistent baseline theme
ggplot2::theme_set(ggplot2::theme_bw())

# plot tumor purity violins
tumor_purity_plot_object <- estimate_metadata_plot_df %>%
  ggplot2::ggplot(ggplot2::aes(x = subgroup,
                               y = TumorPurity,
                               color = subgroup)) +
  ggplot2::geom_violin(show.legend = FALSE) +
  ggplot2::geom_jitter(shape = 16,
                       height = 0,
                       width = 0.1,
                       show.legend = FALSE) +
  ggplot2::facet_wrap(ggplot2::vars(platform))

ggplot2::ggsave(filename = tumor_purity_plot_filename,
                plot = tumor_purity_plot_object,
                width = 7.5,
                height = 7.5)

# plot stromal vs. immune scores
stromal_immune_plot_object <- estimate_metadata_plot_df %>%
  ggplot2::ggplot(ggplot2::aes(x = StromalScore,
                               y = ImmuneScore,
                               color = subgroup)) +
  ggplot2::geom_point(shape = 16) +
  ggplot2::facet_wrap(ggplot2::vars(platform))

ggplot2::ggsave(filename = stromal_immune_plot_filename,
                plot = stromal_immune_plot_object,
                width = 7.5,
                height = 7.5)

# plot estimate score vs. tumor purity
estimate_score_tumor_purity_plot_object <- estimate_metadata_plot_df %>%
  ggplot2::ggplot(ggplot2::aes(x = ESTIMATEScore,
                               y = TumorPurity,
                               color = subgroup)) +
  ggplot2::geom_point(shape = 16) +
  ggplot2::facet_wrap(ggplot2::vars(platform))

ggplot2::ggsave(filename = estimate_score_tumor_purity_plot_filename,
                plot = estimate_score_tumor_purity_plot_object,
                width = 7.5,
                height = 7.5)
