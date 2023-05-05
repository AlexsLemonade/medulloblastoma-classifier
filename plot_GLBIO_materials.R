# GLBIO presentation and poster materials

library(tidyverse)
library(patchwork)

################################################################################
# color schemes
################################################################################

# palette.colors returns a named list, which must be unnamed to not match color names to category labels
okabe_ito_colors <- unname(grDevices::palette.colors(palette = "Okabe-Ito"))

# set subgroup colors
subgroup_colors <- c(G3 = okabe_ito_colors[2],
                     G4 = okabe_ito_colors[3],
                     SHH = okabe_ito_colors[4],
                     WNT = okabe_ito_colors[8])

# set platform colors
platform_colors = c(Array = okabe_ito_colors[1],
                    `RNA-seq` = okabe_ito_colors[9])

################################################################################
# data necessary for producing plots
################################################################################

# directories
processed_data_dir <- here::here("processed_data")

# input filepaths
bulk_metadata_filepath <- file.path(processed_data_dir,
                                    "bulk_metadata.tsv")
pseudobulk_metadata_filepath <- file.path(processed_data_dir,
                                          "pseudobulk_metadata.tsv")

# read in metadata files
bulk_metadata_df <- readr::read_tsv(bulk_metadata_filepath,
                                    show_col_types = FALSE)
pseudobulk_metadata_df <- readr::read_tsv(pseudobulk_metadata_filepath,
                                          show_col_types = FALSE)


#### output filepaths (go in GLBIO folder?)

# data set summaries
n_samples_per_subgroup_plot_filepath
bulk_mb_studies_plot_filepath
prop_samples_per_platform_plot_filepath

# baseline modeling

# pseudobulk and single cell

mb_subgroups <- c("G3", "G4", "SHH", "WNT")

################################################################################
# Data set summaries
################################################################################

# n samples per subgroup

subgroup_numbers_df <- bulk_metadata_df |>
  dplyr::bind_rows(pseudobulk_metadata_df |>
                     dplyr::mutate(subgroup = "MB scRNA",
                                   platform = "RNA-seq") |>
                     dplyr::select(sample_accession = title,
                                   subgroup,
                                   study,
                                   is_duplicate,
                                   platform)) |>
  dplyr::count(subgroup, platform) |>
  dplyr::mutate(mb_subgroup = dplyr::case_when(subgroup %in% mb_subgroups ~ "MB subgroups (bulk)",
                                               is.na(subgroup) ~ "MB subgroups (bulk)",
                                               TRUE ~ "Other data sets"))

n_samples_per_subgroup_plot_object <- subgroup_numbers_df |>
  ggplot2::ggplot(ggplot2::aes(x = subgroup,
                               y = n,
                               fill = platform)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(data = subgroup_numbers |>
                       dplyr::group_by(subgroup, mb_subgroup) |>
                       dplyr::summarize(total = sum(n),
                                        .groups = "drop"),
                     ggplot2::aes(x = subgroup,
                                  y = total,
                                  label = total,
                                  fill = NULL),
                     vjust = 0,
                     nudge_y = 5) +
  ggplot2::facet_grid(. ~ mb_subgroup,
                      scales = "free_x",
                      space = "free") +
  ggplot2::scale_fill_manual(values = platform_colors) +
  ggplot2::labs(x = "Subgroup",
                y = "Count",
                fill = "Platform") +
  ggplot2::theme_bw()

ggplot2::ggsave(filename = n_samples_per_subgroup_plot_filepath,
                plot = n_samples_per_subgroup_plot_object,
                height = 6,
                width = 7.5)

# collection of studies in this project

bulk_mb_studies_plot_object <- bulk_metadata_df |>
  dplyr::filter(subgroup %in% mb_subgroups) |>
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_infreq(study),
                               fill = subgroup)) +
  ggplot2::geom_bar() +
  ggplot2::facet_grid(. ~ platform,
                      scales = "free_x",
                      space = "free") +
  ggplot2::scale_fill_manual(values = subgroup_colors) +
  ggplot2::labs(x = "Study",
                y = "Count",
                fill = "Subgroup") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                     hjust = 1))

ggplot2::ggsave(filename = bulk_mb_studies_plot_filepath,
                plot = bulk_mb_studies_plot_object,
                height = 6,
                width = 7.5)

# Proportion of samples per platform

prop_samples_per_platform_plot_object <- bulk_metadata_df |>
  dplyr::filter(subgroup %in% mb_subgroups) |>
  ggplot2::ggplot(ggplot2::aes(x = subgroup,
                               fill = platform)) +
  ggplot2::geom_bar(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(x = "Subgroup",
                y = "Percentage",
                fill = "Platform") +
  ggplot2::scale_fill_manual(values = platform_colors) +
  ggplot2::theme_bw()

ggplot2::ggsave(filename = prop_samples_per_platform_plot_filepath,
                plot = prop_samples_per_platform_plot_object,
                height = 6,
                width = 7.5)

# Combined data set summary plot

(n_samples_per_subgroup_plot_object + prop_samples_per_platform_plot_object + patchwork::plot_layout(widths = c(3,1), guides = "collect")) / bulk_mb_studies_plot_object +
  patchwork::plot_annotation(tag_levels = "a")
