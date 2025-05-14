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
                              genes_in_classifier,
                              prop_genes_detected = 0.25,
                              platform = "scRNA-seq") {
  # Applies prediction model to gene expression matrix
  #
  # Inputs:
  #  sample_acc: sample accession used for filtering metadata out of metadata_df
  #  sce_filepath: file path to a single cell experiment object RDS file
  #  metadata_df: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #  classifier: classifier model (model class must be OnevsrestScheme_TSP for kTSP or rule_based_RandomForest for random forest)
  #  genes_in_classifier: a vector of genes used in the rules-based classifier
  #  prop_genes_detected: what proportion of genes used in the classifier need to be detected to move forward with prediction? (default: 0.25)
  #  platform: gene expression platform for this sample (default: scRNA-seq)
  #
  # Outputs:
  #  Returns a test object

  if ( class(classifier) == "OnevsrestScheme_TSP" ) {

    model_type <- "ktsp"

  } else if ( class(classifier) == "rule_based_RandomForest" ) {

    model_type <- "rf"

  } else {

    stop("classifier model class must be OnevsrestScheme_TSP for kTSP or rule_based_RandomForest for random forest")

  }

  # read in single cell experiment object
  sce_object <- readr::read_rds(sce_filepath)

  # read in gene expression matrix
  genex_df_this_sample <- SingleCellExperiment::counts(sce_object)

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

  # get the study of a given sample
  sample_study <- metadata_df |>
    dplyr::filter(sample_accession == sample_acc) |>
    dplyr::pull(study)

  # If the sample_subgroup vector is longer than 1, more than one row has this
  # sample accession, meaning that sample accession is not a unique identifier
  # in this metadata (and this should be addressed by the user)
  
  if (length(sample_subgroup) != 1) {
    
    stop(glue::glue("Sample accession ", sample_acc, " is not unique in scRNA metadata."))
    
  }
  
  # test_*() function requires a metadata file
  metadata_df_this_sample <- tibble::tibble(index = 1:n_cells,
                                            sample_accession = stringr::str_c(sample_acc,
                                                                              1:n_cells,
                                                                              sep = "_"),
                                            study = sample_study,
                                            subgroup = sample_subgroup,
                                            platform = platform)
  
  # just genes in classifier
  classifier_genex <- genex_df_this_sample[which(rownames(genex_df_this_sample) %in% genes_in_classifier), ]
  # find the number of classifier-relevant genes detected
  num_genes_detected <- apply(classifier_genex, 2, function(x) sum(x > 0))
  # find cells to retain for prediction
  cells_to_retain <- which((num_genes_detected / length(genes_in_classifier)) > prop_genes_detected)
  # unclassified cell identifiers
  unclassified_cell_ids <- colnames(genex_df_this_sample)[-cells_to_retain]
  
  # for prediction, only retain the cells above the threshold
  genex_df_this_sample <- genex_df_this_sample[, cells_to_retain]
  metadata_df_this_sample <- metadata_df_this_sample |>
    dplyr::filter(sample_accession %in% colnames(genex_df_this_sample))

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
    
    if (length(unclassified_cell_ids) > 0) {
      
      # build out a prediction matrix for the unclassified cells to append
      # to test_object$model_output
      pred_mat_unclassified <- data.frame(
        cell_id = unclassified_cell_ids,
        G3 = NA_real_,
        G4 = NA_real_,
        SHH = NA_real_,
        WNT = NA_real_,
        max_score = "Unclassified",
        tie_flag = NA_character_
      ) |>
        tibble::column_to_rownames("cell_id")
      
      # tack on the unclassified cells
      test_object$model_output <- data.frame(test_object$model_output) |>
        dplyr::bind_rows(pred_mat_unclassified)

    }

  } else if (model_type == "rf") {

    test_object <- test_rf(genex_df_test = genex_df_this_sample,
                           metadata_df_test = metadata_df_this_sample,
                           classifier = classifier)
    
    if (length(unclassified_cell_ids) > 0) {
      # build out a prediction matrix for the unclassified cells to append to
      # test_object$model_output
      pred_mat_unclassified <- data.frame(
        cell_id = unclassified_cell_ids,
        G3 = NA_real_,
        G4 = NA_real_,
        SHH = NA_real_,
        WNT = NA_real_
      ) |>
        tibble::column_to_rownames("cell_id")
      
      # tack on the unclassified cells
      test_object$model_output <- data.frame(test_object$model_output) |>
        dplyr::bind_rows(pred_mat_unclassified)
    }
  }

  if (length(unclassified_cell_ids) > 0) {
    # make a data frame of the unclassified cells to rbind to the test object
    pred_unclassified_df <- tibble::tibble(
      sample_accession = unclassified_cell_ids,
      predicted_labels = "Unclassified"
    )
    
    # add unclassified cells & return
    test_object$predicted_labels_df <- test_object$predicted_labels_df |>
      dplyr::bind_rows(pred_unclassified_df)
  }
  
  return(test_object)

}
