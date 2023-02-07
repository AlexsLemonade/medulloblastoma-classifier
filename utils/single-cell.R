# Functions related to SingleCellExperiment operations
#
# Chante Bethell
# February 2023

convert_dataframe_list_to_sce <- function(df_list, results_dir) {
  # Purpose: Convert a list of data frame objects into SingleCellExperiment
  # objects and write to file
  
  # Args:
  #   df_list: the list of data frames to be converted into SingleCellExperiment
  #            objects
  #   results_dir: the output directory for storing each individual 
  #                SingleCellExperiment RDS file
  
  # convert matrices to SingleCellExperiment objects
  sce_list <- lapply(df_list, function(x) SingleCellExperiment::SingleCellExperiment(x))
  
  # define the output file names
  output_filenames <- paste0(names(sce_list), "_sce.rds")
  
  # write each SingleCellExperiment object to file
  for (i in 1:length(output_filenames)) {
    if (is.null(sce_list[[i]]) |
        !(class(sce_list[[i]]) == "SingleCellExperiment")) {
      stop(
        glue::glue(
          "There is at least one empty SingleCellExperiment object. For 
          successful conversion, please ensure that the input for the following 
          is an expression matrix with genes as rownames:
          {names(sce_list)[[i]]}"
        )
      )
    }
    write_rds(sce_list[[i]], file.path(results_dir, output_filenames[[i]]))
  }
  
  return(sce_list)
}
