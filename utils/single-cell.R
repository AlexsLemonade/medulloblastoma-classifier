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

add_sce_umap <- function(sce, seed = 2023) {
  # Purpose: Add UMAP results to SingleCellExperiment object

  # Args:
  #   sce: the SingleCellExperiment object with counts data
  #   seed: an integer to set the seed as for reproducibility

  # set the seed for reproducible results
  set.seed(seed)

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

test_single_cells <- function(sample_acc,
                              sce_filepath,
                              metadata_df,
                              gene_map_df,
                              labels,
                              classifier,
                              model_type,
                              study = "GSE119926",
                              platform = "scRNA-seq") {
  # Applies prediction model to gene expression matrix
  #
  # Inputs:
  #  sample_acc: sample accession used for filtering metadata out of metadata_df
  #  sce_filepath: file path to a single cell experiment object RDS file
  #  metadata_df: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #  classifier: classifier model
  #  model_type: prediction model used, must be one of 'ktsp' or 'rf'
  #  study: study ID for this sample (default: GSE119926)
  #  platform: gene expression platform for this sample (default: scRNA-seq)
  #
  # Outputs:
  #  Returns a test object

  if ( model_type %in% c("ktsp", "rf") ) {

    if ( model_type == "ktsp" & class(classifier) != "OnevsrestScheme_TSP" ) {

      stop("Type of model specified by classifier parameter in test_single_cells() does not match model_type, which must be 'ktsp' or 'rf'.")

    }

    if ( model_type == "rf" & class(classifier) != "rule_based_RandomForest" ) {

      stop("Type of model specified by classifier parameter in test_single_cells() does not match model_type, which must be 'ktsp' or 'rf'.")

    }

  } else {

    stop("model_type in test_single_cells() must be 'ktsp' or 'rf'.")

  }

  # read in single cell experiment object
  sce_object <- readr::read_rds(sce_filepath)

  # read in gene expression matrix
  genex_df_this_sample <- sce_object@assays@data@listData[[1]]

  # get number of cells to be tested
  n_cells <- ncol(genex_df_this_sample)

  # set column names as the same sample accession across all cells (columns)
  names(genex_df_this_sample) <- stringr::str_c(rep(sample_acc, n_cells),
                                                1:n_cells,
                                                sep = "_")

  # get subgroup of sample given
  sample_subgroup <- metadata_df |>
    dplyr::filter(sample_accession == sample_acc) |>
    dplyr::pull(subgroup)

  if (length(sample_subgroup) != 1) {

    stop("Wait -- sample accession is not unique in scRNA metadata")

  }

  # test_*() function requires a metadata file
  metadata_df_this_sample <- tibble::tibble(index = 1:n_cells,
                                            sample_accession = stringr::str_c(sample_acc,
                                                                              1:n_cells,
                                                                              sep = "_"),
                                            study = study,
                                            subgroup = sample_subgroup,
                                            platform = platform)

  # check_input_files() is sourced from utils/modeling.R. This function ensures
  #  the genex_df and metadata_df given to the test_*() function are properly
  #  formatted and consistent with each other
  check_input_files(genex_df = genex_df_this_sample,
                    metadata_df = metadata_df_this_sample)

  # predict the subgroup of each observation (individual cell) using given model
  if (model_type == "ktsp") {

    test_object <- test_ktsp(genex_df_test = genex_df_this_sample,
                             metadata_df_test = metadata_df_this_sample,
                             classifier = classifier,
                             labels = labels)

  } else if (model_type == "rf") {

    test_object <- test_rf(genex_df_test = genex_df_this_sample,
                           metadata_df_test = metadata_df_this_sample,
                           classifier = classifier)

  }

  return(test_object)

}
