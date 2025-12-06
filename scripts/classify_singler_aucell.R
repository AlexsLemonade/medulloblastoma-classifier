# J. Taroni
# 2025
# Add SingleR and AUCell annotation information to SingleCellExperiment objects
# for the Hovestadt et al. Smart-seq2 data
#
# USAGE: Rscript classify_singler_aucell.R

#### Directories and files -----------------------------------------------------

# Directories
processed_data_dir <- here::here("processed_data")
single_cell_data_dir <- here::here(processed_data_dir, "single_cell")
smartseq_data_dir <- here::here(single_cell_data_dir, "GSE119926/sce")

# Reference files for annotation
singler_ref_file <- fs::path(
  single_cell_data_dir,
  "BlueprintEncode_SingleR_model.rds"
)
panglaodb_ref_file <- fs::path(
  single_cell_data_dir,
  "PanglaoDB_oligodendrocytes_ENSG.rds"
)

# Grab all the SingleCellExperiment files
sce_files <- fs::dir_ls(path = smartseq_data_dir,
                        glob = "*_sce.rds")

#### Functions -----------------------------------------------------------------

classify_cells <- function(sce) {
  # For a SingleCellExperiment object, classify cells with a SingleR model
  # using BlueprintEncode data as a reference and add AUCell scores for
  # PanglaoDB oligodendrocytes + assignments using automatic thresholding from
  # AUCell
  #
  # Returns the SingleCellExperiment with the cell type annotation data in
  # the colData

  # Classify cell types using BlueprintEncode reference
  singler_results <- SingleR::classifySingleR(
    test = sce,
    trained = singler_ref,
    fine.tune = TRUE
  )

  # Data frame version for joining to colData downstream
  singler_results_df <- singler_results |>
    as.data.frame() |>
    tibble::rownames_to_column("cell_id")

  # Run AUCell
  auc_results <- AUCell::AUCell_run(
    sce,
    geneSets = oligodendrocyte_gene_set
  )

  # Get AUCell scores as a data frame
  auc_df <- auc_results@assays@data$AUC |>
    # Transpose
    t() |>
    # Convert to data frame
    as.data.frame() |>
    # Make the cell identifiers a column
    tibble::rownames_to_column("cell_id")

  # Assign cells with automatically selected threshold
  auc_calls <- AUCell::AUCell_exploreThresholds(
    auc_results,
    plotHist = FALSE,
    assignCells = TRUE
  )

  # Add oligodendrocyte assignments using automated threshold
  auc_df <- auc_df |>
    dplyr::mutate(oligodendrocyte_assignment = dplyr::if_else(
      cell_id %in% auc_calls$panglaodb_oligodendrocyte$assignment,
      TRUE,
      FALSE
    ))

  # Put all the colData together
  coldata_df <- colData(sce) |>
    as.data.frame() |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::left_join(
      auc_df,
      by = "cell_id"
    ) |>
    dplyr::left_join(
      singler_results_df,
      by = "cell_id"
    )

  # Add back into SingleCellExperiment object
  colData(sce) <- DataFrame(
    coldata_df,
    row.names = coldata_df$cell_id
  )

  # Return SingleCellExperiment
  return(sce)

}

#### Read in reference files ---------------------------------------------------

# Trained SingleR model using BlueprintEncode data as a reference
singler_ref <- readr::read_rds(singler_ref_file)

# PanglaoDB oligodendrocytes marker genes in GSEABase::GeneSet format
oligodendrocyte_gene_set <- readr::read_rds(panglaodb_ref_file)

#### Classify ------------------------------------------------------------------

# Read in SingleCellExperiment objects into a list
sce_list <- sce_files |>
  purrr::map(readr::read_rds)

# Add classifications
sce_list <- sce_list |>
  purrr::map(classify_cells)

# Write to files
purrr::walk2(
  sce_list,
  sce_files,
  \(x, y) readr::write_rds(x, y)
)
