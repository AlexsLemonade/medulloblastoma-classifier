# J. Taroni
# 2025
#
# Extract marker genes from Aldinger et al. 2021. human cerebellar development
# data https://doi.org/10.1038/s41593-021-00872-y
# Downloaded from UCSC Cell Browser: https://cells.ucsc.edu/?ds=cbl-dev
#
# Usage: Rscript extract_cbl_dev_marker_genes.R

#### Set up --------------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)

#### Files and directories -----------------------------------------------------

# Aldinger et al. 2021. data as a Seurat object
srt_file <- here::here("data/CellBrowser/seurat.rds")

# Output file
markers_file <- here::here(
  "processed_data/single_cell/aldinger_cbl_dev_marker_genes.tsv"
)

#### Functions -----------------------------------------------------------------

# Gene conversion function to stay consistent within project
source(here::here("utils/convert_gene_names.R"))

extract_top_n_markers <- function(marker_dframe,
                                  cell_type_name,
                                  top_n = 100) {
  # Extract the top n markers (sorted by FDR) for marker Dframe and convert from
  # gene symbol to Ensembl gene identifier

  # Convert to data.frame and grab top n genes
  marker_df <- marker_dframe |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    dplyr::slice_min(order_by = FDR, n = top_n) |>
    # Select the columns that will be present in every marker gene df
    dplyr::select(
      gene,
      p.value,
      FDR,
      summary.logFC
    ) |>
    dplyr::mutate(cell_type = cell_type_name)

  # Convert to Ensembl gene identifiers
  marker_df <- convert_gene_names(marker_df,
                                  gene_column_before = "gene",
                                  gene_column_after = "gene",
                                  map_from = "SYMBOL",
                                  map_to = "ENSEMBL")
  return(marker_df)
}

#### Read in, convert, and normalize data --------------------------------------

# Read in Seurat object
srt <- readr::read_rds(srt_file)

# Convert to SingleCellExperiment, RNA only
sce <- as.SingleCellExperiment(srt, assay = "RNA")

# Remove large Seurat object
rm(srt)

# Normalize in a way that's consistent with the rest of the project
sce <- scater::logNormCounts(sce,
                             transform = "log")

# Remove reduced dims that were calculated using, presumably, the SCTransform
# data to make a bit smaller
reducedDim(sce) <- NULL

#### Find marker genes ---------------------------------------------------------

# Calculate marker genes for cell type vs. all other cell type
markers <- scran::findMarkers(
  sce,
  groups = sce$fig_cell_type,  # Use Aldinger et al. cell types
  pval.type = "all"  # Genes that differentiate from *all* other cell types
)

# Top 100 markers for every cell type
markers_df <- markers |>
  purrr::imap(
      \(df, ct) extract_top_n_markers(df, ct, top_n = 100)
    ) |>
  dplyr::bind_rows()

# Write to file
readr::write_tsv(markers_df, markers_file)
