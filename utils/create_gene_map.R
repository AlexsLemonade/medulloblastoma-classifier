# Create data frame to map gene names
#
# Steven Foltz
# November 2022

option_list <- list(
  optparse::make_option("--annotationhub_snapshot_date",
                        default = NA_character_,
                        help = "AnnotationHub snapshot date")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

library(magrittr)

# set up directories and output filepaths
processed_data_dir <- here::here("processed_data")

gene_map_filepath <- file.path(processed_data_dir,
                               "gene_map.tsv")

# set up AnnotationHub
AnnotationHub::setAnnotationHubOption("ASK", FALSE) # download without asking
ah <- AnnotationHub::AnnotationHub()
AnnotationHub::snapshotDate(ah) <- opt$annotationhub_snapshot_date # set date
hs_orgdb <- AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]] # humans

# starting with ENSEMBL, map to ENTREZID and SYMBOL
# remove duplicates
AnnotationDbi::select(x = hs_orgdb,
                      keys = AnnotationDbi::keys(hs_orgdb, "ENSEMBL"),
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "ENSEMBL") %>%
  dplyr::mutate(dup_ensembl = duplicated(ENSEMBL),
                dup_entrezid = duplicated(ENTREZID),
                dup_symbol = duplicated(SYMBOL)) %>%
  dplyr::filter(!dup_ensembl, !dup_entrezid, !dup_symbol) %>%
  dplyr::select(ENSEMBL, ENTREZID, SYMBOL) %>%
  readr::write_tsv(file = gene_map_filepath)
