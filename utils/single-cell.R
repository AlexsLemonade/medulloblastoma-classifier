# Functions related to SingleCellExperiment operations
#
# Chante Bethell
# February 2023

convert_dataframe_list_to_sce <- function(df_list, results_dir) {
  # Purpose: Convert a list of data frame or matrix objects into
  # SingleCellExperiment objects and write to file
  
  # Args:
  #   df_list: the list of data frames or matrices to be converted into
  #            SingleCellExperiment objects; each object should be an expression
  #            matrix with genes as row names
  #   results_dir: the output directory for storing each individual 
  #                SingleCellExperiment RDS file
  
  # convert matrices to SingleCellExperiment objects
  sce_list <- lapply(df_list, function(x) SingleCellExperiment::SingleCellExperiment(x))
  
  # define the output file names
  output_filenames <- paste0(names(sce_list), "_sce.rds")
  
  # write each SingleCellExperiment object to file
  for (i in 1:length(output_filenames)) {
    if (is.null(sce_list[[i]]) |
        !("SingleCellExperiment" %in% class(sce_list[[i]]))) {
      stop(
        glue::glue(
          "SingleCellExperiment object conversion failed (result is NULL or is not a SingleCellExperiment). For 
          successful conversion, please ensure that the input for the following sample
          is an expression matrix with genes as rownames:
          {names(sce_list)[[i]]}"
        )
      )
    }
    write_rds(sce_list[[i]], file.path(results_dir, output_filenames[[i]]))
  }
  
  return(sce_list)
}
