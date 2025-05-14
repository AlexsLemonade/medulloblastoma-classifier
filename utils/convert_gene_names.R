set_up_AnnotationHub <- function(ah_date = "2022-10-30") {

  # Run this one time before doing anything with AnnotationHub in parallel
  #
  # Inputs
  #  ah_date: AnnotationHub snapshot date (default: 2022-10-30)

  # set up AnnotationHub
  AnnotationHub::setAnnotationHubOption("ASK", FALSE) # download without asking
  ah <- suppressMessages(AnnotationHub::AnnotationHub())
  AnnotationHub::snapshotDate(ah) <- ah_date # set AnnotationHub snapshot date
  hs_orgdb <- suppressMessages(AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]) # humans

}

convert_gene_names <- function(genex_df,
                               gene_column_before,
                               gene_column_after,
                               map_from,
                               map_to,
                               ah_date = "2022-10-30") {

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
  #  map_from: the annotation style of the original gene column (e.g., ENSEMBL, ENTREZID, SYMBOL)
  #  map_to: the annotation style of the new gene column (e.g., ENSEMBL, ENTREZID, SYMBOL)
  #  ah_date: AnnotationHub snapshot date (default: 2022-10-30)
  #
  # Output
  #  Gene expression matrix with a "gene" column

  # check that gene_column_before is in genex_df

  if(!(gene_column_before %in% names(genex_df))) {

    stop("Gene column missing from gene expression matrix in convert_gene_names().")

  }

  # set up AnnotationHub
  AnnotationHub::setAnnotationHubOption("ASK", FALSE) # download without asking
  ah <- suppressMessages(AnnotationHub::AnnotationHub())
  AnnotationHub::snapshotDate(ah) <- ah_date # set AnnotationHub snapshot date
  hs_orgdb <- suppressMessages(AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]) # humans

  gene_map_vector <- suppressMessages(AnnotationDbi::mapIds(x = hs_orgdb,
                                                            keys = genex_df[[gene_column_before]],
                                                            keytype = map_from,
                                                            column = map_to,
                                                            multiVals = "first"))

  # check that dimensions match up
  if (length(gene_map_vector) != nrow(genex_df)) {

    stop("In convert_gene_names(), length of gene map does not match number of genes in expression matrix.")

  }

  # check that the order of genes is the same in gene_map_vector and genex_df
  if (!all(names(gene_map_vector) == genex_df[[gene_column_before]])) {

    stop("In convert_gene_names(), order of gene map does not match order of genes in expression matrix.")

  }

  genex_df_mapped <- genex_df |>
    dplyr::select(-tidyselect::all_of(gene_column_before)) |>
    tibble::add_column("{gene_column_after}" := gene_map_vector,
                       .before = 1) |>
    dplyr::filter(!is.na(gene_map_vector), # if no map between key and value
                  !duplicated(gene_map_vector)) # if keys map to the same value

  return(genex_df_mapped)

}
