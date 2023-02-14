# Add dimensionality reduction and clustering results to pseudobulk
# SingleCellExperiment objects
#
# Chante Bethell
# February 2023

suppressMessages(library(tidyverse))

processed_data_dir <- here::here("processed_data")
pseudobulk_sce_dir <- file.path(processed_data_dir, "pseudobulk_sce")
utils_dir <- here::here("utils")
pseudobulk_metadata_input_filepath <- file.path(processed_data_dir,
                                                "pseudobulk_metadata.tsv")
source(file.path(utils_dir, "single-cell.R"))

# read in pseudo-bulk metadata file
pseudobulk_metadata <- readr::read_tsv(pseudobulk_metadata_input_filepath,
                                       show_col_types = FALSE)

# read in individual SCE objects
pseudobulk_sce_list <- purrr::map(pseudobulk_metadata$title,
                                  function(x)
                                  readr::read_rds(file.path(pseudobulk_sce_dir, paste0(x, "_sce.rds"))))

names(pseudobulk_sce_list) <- pseudobulk_metadata$title

# define the names of the 3 funky samples, as well as one sample from each 
# MB subgroup for testing
sample_names_for_comparison <-
  c("BCH825",
    "Med2312FH",
    "SJ625",
    "BCH1205",
    "BCH1031",
    "MUV41",
    "BCH807")

# create plot list of UMAPs with cluster assignments for each defined sample
plot_list <- purrr::map(sample_names_for_comparison, function(x)
  scater::plotReducedDim(
    pseudobulk_sce_list[[x]],
    dimred = "UMAP",
    colour_by = "louvain_10",
    text_by = "louvain_10",
    text_size = 2,
    point_size = 0.4,
    point_alpha = 0.5,
    
  ) + theme_bw() +
    labs(caption = paste0("louvain clustering with a nearest neighbours value of 10 for sample ",
                          x)) +
    theme(text = element_text(size = 6),
          plot.caption = element_text(hjust = 0.5)) +
    guides(col = guide_legend(
      "Cluster assignment", override.aes = list(size = 0.5)
    )))

# use cowplot to combine plots and save as PDF
pdf(file.path(pseudobulk_sce_dir, "pseudobulk_umap_plots.pdf"))
cowplot::plot_grid(
  plotlist = plot_list,
  ncol = 2,
  byrow = FALSE,
  vjust = 0
) 
dev.off()
