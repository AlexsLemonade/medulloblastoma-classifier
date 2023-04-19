convert_gene_names <- function() {
  
  # Converts the gene column of a gene expression matrix to a different annotation style
  # e.g. from SYMBOL to ENSEMBL
  #
  # Many functions in this project require the gene column be converted to row name.
  # This function does NOT do that conversion. This function returns a "gene" column
  # along with the gene expression data columns.
  #
  # Inputs
  #  genex_df: gene expression matrix with a gene name column
  #  gene_column: name of the gene column (default: gene)
  #  gene_map_df: data frame mapping gene annotations across columns (one gene per row)
  #  map_from: the annotation style of the original gene column
  #  map_to: the annotation style of the new gene column
  #
  # Output
  #  Gene expression matrix with a "gene" column, now with new gene annotation style
  
  # check that gene_column is in genex_df
  
  # check that map_from and map_to are in gene_map_df
  
  genex_df |>
    dplyr::left_join(gene_map_df |> dplyr::select(map_from, map_to),
                     by = c(gene_column = map_from)) |>
    dplyr::filter(!duplicated(gene_column),
                  !duplicated(map_to),
                  !is.na(map_to)) |>
    dplyr::select(-gene_column,
                  gene = map_to,
                  tidyselect::everything()) |>
    return()
  
}
