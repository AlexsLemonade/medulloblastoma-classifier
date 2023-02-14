# Functions related to SingleCellExperiment operations
#
# Chante Bethell
# February 2023

convert_dataframe_list_to_sce <- function(df_list) {
  # Purpose: Convert a list of data frame or matrix objects into
  # SingleCellExperiment objects
  
  # Args:
  #   df_list: the list of data frames or matrices to be converted into
  #            SingleCellExperiment objects; each object should be an expression
  #            matrix with genes as row names
  
  # convert matrices to SingleCellExperiment objects
  sce_list <- lapply(df_list, function(x) SingleCellExperiment::SingleCellExperiment(assays=list(counts=x)))
  
  # write each SingleCellExperiment object to file
  for (i in 1:length(sce_list)) {
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
  }
  return(sce_list)
}

add_sce_umap <- function(sce) {
  # Purpose: Add UMAP results to SingleCellExperiment object
  
  # Args:
  #   sce: the SingleCellExperiment object with counts data
  
  # keep positive counts for `logNormCounts()`
  sce <- sce[, colSums(SingleCellExperiment::counts(sce)) > 0]
  
  # normalize and log transform
  normalized_sce <- scater::logNormCounts(sce)
  
  # model gene variance using `scran:modelGeneVar()`
  gene_variance <- scran::modelGeneVar(normalized_sce)
  
  # select the most variable genes
  subset_genes <- scran::getTopHVGs(gene_variance, n = 2000)
    
  # add PCA to normalized sce
  normalized_sce <- scater::runPCA(normalized_sce, subset_row = subset_genes)
  
  # calculate a UMAP matrix using the PCA results
  normalized_sce <- scater::runUMAP(normalized_sce, dimred = "PCA")
  
  return(normalized_sce)
}

perform_graph_clustering <- function(sce,
                                     nearest_neighbors = 10,
                                     cluster_type = "louvain",
                                     seed = 2021) {
  # Purpose: Perform the graph based clustering on a normalized SingleCellExperiment
  # object
  
  # Args:
  #   sce: normalized SingleCellExperiment object
  #   nearest_neighbors: number of nearest neighbors to use when
  #                      calculating/plotting the clustering results; default
  #                      is 10
  #   cluster_type: the type of graph-based clustering method that is being
  #     tested -- can be "walktrap" or "louvain"; the default is "louvain"
  #   seed: an integer to set the seed as for reproducibility
  
  # determine weighting type to use based on graph detection algorithm specified
  # if louvain is used, use jaccard
  # if walktrap is used, use rank
  weighting_type <-
    ifelse(cluster_type == "louvain", "jaccard", "rank")
  
  # set cluster name
  cluster_name <-
    sprintf("%s_%02d", cluster_type, nearest_neighbors)
  
  # set the seed for reproducible results
  set.seed(seed)
  
  # extract the principal components matrix
  pca_matrix <- SingleCellExperiment::reducedDim(sce, "PCA")
  
  # perform graph-based clustering
  clusters <- bluster::clusterRows(
    pca_matrix,
    bluster::NNGraphParam(
      k = nearest_neighbors,
      type = weighting_type,
      cluster.fun = cluster_type
    )
  )
  
  # store cluster results in the SCE object
  sce[[cluster_name]] <- factor(clusters)
  
  return(sce)
  
}
