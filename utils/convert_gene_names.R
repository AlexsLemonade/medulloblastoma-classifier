convert_gene_names <- function(genex_df,
                               gene_column_before,
                               gene_column_after,
                               gene_map_df,
                               map_from,
                               map_to) {
  
  # Converts the gene column of a gene expression matrix to a different annotation style
  # e.g. from SYMBOL to ENSEMBL
  #
  # Many functions in this project require the gene column be converted to row name.
  # This function does NOT do that conversion. This function returns a "gene" column
  # along with the gene expression data columns.
  #
  # Inputs
  #  genex_df: gene expression matrix with a gene name column
  #  gene_column_before: name of the gene column before conversion
  #  gene_column_after: name of the gene column after conversion
  #  gene_map_df: data frame mapping gene annotations across columns (one gene per row)
  #  map_from: the annotation style of the original gene column
  #  map_to: the annotation style of the new gene column
  #
  # Output
  #  Gene expression matrix with a "gene" column
  
  # check that gene_column_before is in genex_df
  
  if(!(gene_column_before %in% names(genex_df))) {
  
    stop("Gene column missing from gene expression matrix in convert_gene_names().")  
    
  }
  
  # check that map_from and map_to are in gene_map_df
  
  if(!all(c(map_from, map_to) %in% names(gene_map_df))) {
    
    stop("Annotation styles missing from gene map in convert_gene_names().")
    
  }
  
  genex_df |>
    dplyr::left_join(gene_map_df |> dplyr::select(map_from, map_to),
                     by = setNames(map_from, gene_column_before)) |>
    dplyr::filter(!duplicated(.data[[gene_column_before]]), # multi-mapped left to right
                  !duplicated(.data[[map_to]]), # multi-mapped right to left
                  !is.na(.data[[map_to]])) |> # no match left to right
    dplyr::select(-dplyr::all_of(gene_column_before)) |> # remove old gene column
    dplyr::relocate(!!gene_column_after := map_to) # rename and move new gene column to first position
  
}