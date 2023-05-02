suppressMessages(library(foreach))

run_many_models <- function(genex_df,
                            metadata_df,
                            labels,
                            model_types = c("ktsp", "rf", "mm2s", "lasso"),
                            initial_seed = 44,
                            n_repeats = 1,
                            n_cores = 1,
                            ktsp_featureNo = 1000,
                            ktsp_n_rules_min = 5,
                            ktsp_n_rules_max = 50,
                            ktsp_weighted = TRUE,
                            rf_num.trees = 500,
                            rf_genes_altogether = 50,
                            rf_genes_one_vs_rest = 50,
                            rf_gene_repetition = 1,
                            rf_rules_altogether = 50,
                            rf_rules_one_vs_rest = 50,
                            rf_weighted = TRUE,
                            mm2s_ah_date = "2022-10-30") {

  # Wrapper function to run many modeling jobs in parallel. New train/test sets
  # are created for each repeat. The same train/test data is used for each model.
  #
  # Inputs
  #  genex_df: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #  model_types: vector of model types (must be one or more of 'ktsp', 'rf', 'mm2s', or 'lasso') (default: all of them)
  #  initial_seed: seed used to set train/test seeds and modeling seeds (default: 44)
  #  n_repeats: how many times to repeat each modeling type (default: 1)
  #  n_cores: number of cores to use (default: 1)
  #
  #  kTSP parameters:
  #    ktsp_featureNo: number of most informative features to filter down to (kTSP only) (default: 1000)
  #    ktsp_n_rules_min: minimum number of rules allowed for kTSP modeling (default: 5)
  #    ktsp_n_rules_max: maximum number of rules allowed for kTSP modeling (default: 50)
  #    ktsp_weighted: logical, if TRUE (default) use one-vs-one and platform-wise comparisons to add weight to smaller subgroups and platforms
  #
  #  RF parameters:
  #    rf_num.trees: number of trees used in RF modeling (use more trees given more features) (default: 500)
  #    rf_genes_altogether: number of top genes used when comparing all classes together (default: 50)
  #    rf_genes_one_vs_rest: number of top genes used when comparing each class against rest of classes (default: 50)
  #    rf_gene_repetition: number of times a gene can be used throughout set of rules (default: 1)
  #    rf_rules_altogether: number of top rules used when comparing all classes together (default: 50)
  #    rf_rules_one_vs_rest: number of top rules used when comparing each class against rest of classes (default: 50)
  #    rf_weighted: logical, if TRUE (default) use one-vs-rest and platform-wise comparisons to add weight to smaller subgroups and platforms
  #
  #  MM2S parameters:
  #    mm2s_ah_date: AnnotationHub snapshot date (default: 2022-10-30)
  #
  # Output
  #  Model list with levels for repeat number and model type

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df,
                    metadata_df = metadata_df)

  # model types should be a list with elements limited to "ktsp", "rf", "mm2s", "lasso"
  if (!is.vector(model_types) |
      !all(model_types %in% c("ktsp", "rf", "mm2s", "lasso"))) {

    stop("model_types in run_models() should be a vector limited to 'ktsp', 'rf', 'mm2s', 'lasso'.")

  }

  # n_repeats should be a positive integer
  if (n_repeats < 1 | round(n_repeats) != n_repeats) {

    stop("n_repeats in run_many_models() should be a positive integer.")

  }

  # n_cores should not exceed parallel::detectCores() - 1
  n_cores <- min(n_cores, parallel::detectCores() - 1)

  # ktsp_n_rules_min is required for ktsp and is only used in ktsp
  if ("ktsp" %in% model_types) {

    if (is.na(ktsp_n_rules_min)) {
      stop("ktsp_n_rules_min in run_many_models() cannot be NA with ktsp model type.")
    }

    # ktsp should be a positive integer
    if (ktsp_n_rules_min < 1 | round(ktsp_n_rules_min) != ktsp_n_rules_min) {

      stop("ktsp_n_rules_min in run_models() should be a positive integer.")

    }

    # n_rules_max should be a positive integer and >= ktsp_n_rules_min
    # if n_rules_max is not given, set n_rules_max equal to ktsp_n_rules_min
    # this enables setting a specific number of rules for the model to use
    if (is.na(ktsp_n_rules_max)) {

      ktsp_n_rules_max <- ktsp_n_rules_min

    } else if (ktsp_n_rules_max < 1 | round(ktsp_n_rules_max) != ktsp_n_rules_max) {

      stop("ktsp_n_rules_max in run_many_models() should be NA or a positive integer.")

    }

  } else if (ktsp_n_rules_max < ktsp_n_rules_min) {

    stop("ktsp_n_rules_max cannot be less than ktsp_n_rules_min in run_many_models()")

  }

  # Set initial seed before creating test/train seeds and modeling seeds
  set.seed(initial_seed)

  # Seeds used for determining test/train split for each repeat
  # if n_repeats > 1000, then that becomes the max value for sampling seeds
  # alternatively we could just sample between 1:n_repeats every time
  train_test_seeds <- sample(1:max(1000, n_repeats), size = n_repeats)
  # Seeds used at start of each modeling step (same seed re-used for all model types within each repeat)
  modeling_seeds <- sample(1:max(1000, n_repeats), size = n_repeats)

  # parallel backend
  cl <- parallel::makeCluster(n_cores, outfile = "log") # use log file for troubleshooting
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,
                          c("convert_gene_names",
                            "get_train_test_samples",
                            "run_one_model",
                            "check_input_files",
                            "calculate_confusion_matrix",
                            "train_ktsp",
                            "test_ktsp",
                            "train_rf",
                            "test_rf",
                            "test_mm2s",
                            "train_lasso",
                            "test_lasso"))

  # these namespaces must be loaded for test_MM2S() and test_lasso() functions to work
  # done globally before splitting into parallel processes for consistency
  suppressMessages(library(MM2S))
  suppressMessages(library(glmnet))

  # run n_repeats in parallel
  model_list <- foreach(n = 1:n_repeats) %dopar% {

    # set up this repeat's train/test split
    train_test_samples_list <- get_train_test_samples(genex_df = genex_df,
                                                      metadata_df = metadata_df,
                                                      train_test_seed = train_test_seeds[n],
                                                      proportion_of_studies_train = 0.5)

    # split genex and metadata by train/test status
    genex_df_train <- genex_df |>
      dplyr::select(train_test_samples_list$train)
    genex_df_test <- genex_df |>
      dplyr::select(train_test_samples_list$test)
    metadata_df_train <- metadata_df |>
      dplyr::filter(sample_accession %in% train_test_samples_list$train)
    metadata_df_test <- metadata_df |>
      dplyr::filter(sample_accession %in% train_test_samples_list$test)

    # run model types one at a time
    repeat_list <- purrr::map(model_types,
                              function(x) run_one_model(type = x,
                                                        genex_df_train = genex_df_train,
                                                        genex_df_test = genex_df_test,
                                                        metadata_df_train = metadata_df_train,
                                                        metadata_df_test = metadata_df_test,
                                                        model_seed = modeling_seeds[n],
                                                        labels = labels,
                                                        ktsp_featureNo = ktsp_featureNo,
                                                        ktsp_n_rules_min = ktsp_n_rules_min,
                                                        ktsp_n_rules_max = ktsp_n_rules_max,
                                                        ktsp_weighted = ktsp_weighted,
                                                        rf_num.trees = rf_num.trees,
                                                        rf_genes_altogether = rf_genes_altogether,
                                                        rf_genes_one_vs_rest = rf_genes_one_vs_rest,
                                                        rf_gene_repetition = rf_gene_repetition,
                                                        rf_rules_altogether = rf_rules_altogether,
                                                        rf_rules_one_vs_rest = rf_rules_one_vs_rest,
                                                        rf_weighted = rf_weighted,
                                                        mm2s_ah_date = mm2s_ah_date)) |>
      purrr::set_names(model_types) # set names of each list element corresponding to model type

    # add metadata about this repeat (seeds used, train/test metadata)
    repeat_list[["train_test_seed"]] <- train_test_seeds[n]
    repeat_list[["modeling_seed"]] <- modeling_seeds[n]
    repeat_list[["train_metadata"]] <- metadata_df_train
    repeat_list[["test_metadata"]] <- metadata_df_test

    return(repeat_list)

  }

  # stop parallel backend
  parallel::stopCluster(cl)

  return(model_list)

}

get_train_test_samples <- function(genex_df,
                                   metadata_df,
                                   train_test_seed,
                                   proportion_of_studies_train = 0.5) {

  # Split data into training and test sets. Data gets split at project level.
  # Some proportion of projects become "training", remainder becomes "test".
  # Projects from different platforms (array, RNA-seq) are split separately.
  #
  # Inputs
  #  genex_df: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  train_test_seed: seed used for reproducibility given same input data
  #  proportion_of_studies_train: proportion of studies used as training data
  #
  # Outputs
  #  List of samples used for "train" and "test" sets
  #
  # Methodological choices
  #
  #  When determining the number of studies allocated to train and test sets,
  #  sometimes the number of studies is not a whole number. In that case, we use
  #  the ceiling() function to round up on the number of studies allocated to
  #  the training set, which necessarily rounds down the number of test studies.
  #  For example, with 3 RNA-seq studies and a 50% split for train and test,
  #  ceiling(3*0.5) = ceiling(1.5) = 2 studies are allocated to training, while
  #  3 - 2 = 1 study is allocated to testing. Note: ceiling() always rounds up,
  #  to the next whole number, not necessarily to the nearest whole number.

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df,
                    metadata_df = metadata_df)

  # check proportion_of_studies_train is numeric and 0 < p < 1
  if (!is.numeric(proportion_of_studies_train) |
      proportion_of_studies_train < 0 | proportion_of_studies_train > 1) {
    stop("proportion_of_studies_train must be numeric and in the range [0,1]")
  }

  set.seed(train_test_seed)

  # Vector of array studies
  array_studies <- metadata_df |>
    dplyr::filter(sample_accession %in% names(genex_df),
                  platform == "Array") |>
    dplyr::pull(study) |>
    unique()

  # Vector of RNA-seq studies
  rnaseq_studies <- metadata_df |>
    dplyr::filter(sample_accession %in% names(genex_df),
                  platform == "RNA-seq") |>
    dplyr::pull(study) |>
    unique()

  # Number of studies from each platform
  n_array_studies <- length(array_studies)
  n_rnaseq_studies <- length(rnaseq_studies)

  # Set number of training studies according to training proportion
  # ceiling() rounds the number of training studies up to the next whole number
  n_array_studies_train <- ceiling(n_array_studies*proportion_of_studies_train)
  n_rnaseq_studies_train <- ceiling(n_rnaseq_studies*proportion_of_studies_train)

  # Randomly sample studies for training set
  array_studies_train <- sample(array_studies, size = n_array_studies_train)
  rnaseq_studies_train <- sample(rnaseq_studies, size = n_rnaseq_studies_train)

  # Force other studies to be used as test set
  array_studies_test <- setdiff(array_studies, array_studies_train)
  rnaseq_studies_test <- setdiff(rnaseq_studies, rnaseq_studies_train)

  # Create output list
  train_test_samples_list <- list()

  # Save vector of samples used for training
  train_test_samples_list[["train"]] <- metadata_df |>
    dplyr::filter(study %in% c(array_studies_train,
                               rnaseq_studies_train)) |>
    dplyr::pull(sample_accession)

  # Save vector of samples used for testing
  train_test_samples_list[["test"]] <- metadata_df |>
    dplyr::filter(study %in% c(array_studies_test,
                               rnaseq_studies_test)) |>
    dplyr::pull(sample_accession)

  return(train_test_samples_list)

}

run_one_model <- function(type,
                          genex_df_train,
                          genex_df_test,
                          metadata_df_train,
                          metadata_df_test,
                          model_seed,
                          labels,
                          ktsp_featureNo = 1000,
                          ktsp_n_rules_min = 5,
                          ktsp_n_rules_max = 50,
                          ktsp_weighted = TRUE,
                          rf_num.trees = 500,
                          rf_genes_altogether = 50,
                          rf_genes_one_vs_rest = 50,
                          rf_gene_repetition = 1,
                          rf_rules_altogether = 50,
                          rf_rules_one_vs_rest = 50,
                          rf_weighted = TRUE,
                          mm2s_ah_date = "2022-10-30") {

  # Run a single model specified by the model type, input data, and parameters
  #
  # Inputs
  #  type: model type, one of "ktsp", "rf", "mm2s", or "lasso"
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  genex_df_test: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  metadata_df_test: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  model_seed: seed used for reproducibility in training step
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #
  #  kTSP parameters:
  #    ktsp_featureNo: number of most informative features to filter down to (kTSP only) (default: 1000)
  #    ktsp_n_rules_min: minimum number of rules allowed for kTSP modeling (default: 5)
  #    ktsp_n_rules_max: maximum number of rules allowed for kTSP modeling (default: 50)
  #    ktsp_weighted: logical, if TRUE use one-vs-one and platform-wise comparisons to add weight to smaller subgroups and platforms
  #
  #  RF parameters:
  #    rf_num.trees: number of trees used in RF modeling (use more trees given more features) (default: 500)
  #    rf_genes_altogether: number of top genes used when comparing all classes together (default: 50)
  #    rf_genes_one_vs_rest: number of top genes used when comparing each class against rest of classes (default: 50)
  #    rf_gene_repetition: number of times a gene can be used throughout set of rules (default: 1)
  #    rf_rules_altogether: number of top rules used when comparing all classes together (default: 50)
  #    rf_rules_one_vs_rest: number of top rules used when comparing each class against rest of classes (default: 50)
  #    rf_weighted: logical, if TRUE (default) use one-vs-rest and platform-wise comparisons to add weight to smaller subgroups and platforms
  #
  #  MM2S parameters:
  #    mm2s_ah_date: AnnotationHub snapshot date (default: 2022-10-30)


  # Outputs
  #  List of model elements, including the classifier, test results, and test confusion matrix

  if (type == "ktsp") {

    ktsp_classifier <- train_ktsp(genex_df_train = genex_df_train,
                                  metadata_df_train = metadata_df_train,
                                  model_seed = model_seed,
                                  ktsp_featureNo = ktsp_featureNo,
                                  ktsp_n_rules_min = ktsp_n_rules_min,
                                  ktsp_n_rules_max = ktsp_n_rules_max,
                                  ktsp_weighted = ktsp_weighted)

    ktsp_results <- test_ktsp(genex_df_test = genex_df_test,
                              metadata_df_test = metadata_df_test,
                              classifier = ktsp_classifier,
                              labels = labels)

    ktsp_cm <- calculate_confusion_matrix(predicted_labels = ktsp_results$predicted_labels_df$predicted_labels,
                                          true_labels = metadata_df_test$subgroup,
                                          labels = labels)

    model <- list(classifier = ktsp_classifier,
                  test_results = ktsp_results,
                  cm = ktsp_cm)

  } else if (type == "rf") {

    rf_classifier <- train_rf(genex_df_train = genex_df_train,
                              metadata_df_train = metadata_df_train,
                              model_seed = model_seed,
                              rf_num.trees = rf_num.trees,
                              rf_genes_altogether = rf_genes_altogether,
                              rf_genes_one_vs_rest = rf_genes_one_vs_rest,
                              rf_gene_repetition = rf_gene_repetition,
                              rf_rules_altogether = rf_rules_altogether,
                              rf_rules_one_vs_rest = rf_rules_one_vs_rest,
                              rf_weighted = rf_weighted)

    rf_results <- test_rf(genex_df_test = genex_df_test,
                          metadata_df_test = metadata_df_test,
                          classifier = rf_classifier)

    rf_cm <- calculate_confusion_matrix(predicted_labels = rf_results$predicted_labels_df$predicted_labels,
                                        true_labels = metadata_df_test$subgroup,
                                        labels = labels)

    model <- list(classifier = rf_classifier,
                  test_results = rf_results,
                  cm = rf_cm)

  } else if (type == "mm2s") {

    mm2s_results <- test_mm2s(genex_df_test = genex_df_test,
                              metadata_df_test = metadata_df_test,
                              model_seed = model_seed,
                              ah_date = mm2s_ah_date)

    mm2s_cm <- calculate_confusion_matrix(predicted_labels = mm2s_results$predicted_labels_df$predicted_labels,
                                          true_labels = metadata_df_test$subgroup,
                                          labels = c("G3", "G4", "NORMAL", "SHH", "WNT"))

    model <- list(test_results = mm2s_results,
                  cm = mm2s_cm)

  } else if (type == "lasso") {


    lasso_classifier <- train_lasso(genex_df_train = genex_df_train,
                                    metadata_df_train = metadata_df_train,
                                    model_seed = model_seed)

    lasso_results <- test_lasso(genex_df_test = genex_df_test,
                                metadata_df_test = metadata_df_test,
                                classifier = lasso_classifier)

    lasso_cm <- calculate_confusion_matrix(predicted_labels = lasso_results$predicted_labels_df$predicted_labels,
                                           true_labels = metadata_df_test$subgroup,
                                           labels = labels)

    model <- list(classifier = lasso_classifier,
                  test_results = lasso_results,
                  cm = lasso_cm)

  } else {

    stop("Type of model must be 'ktsp', 'rf', 'mm2s', or 'lasso'.")

  }

  return(model)

}

check_input_files <- function(genex_df,
                              metadata_df) {

  # Checks input files when performing modeling training and testing
  #
  # Inputs
  #  genex_df: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #
  # Output
  #  None (function calls 'stop' if anything is wrong)

  # Check that metadata_df has all necessary columns
  if (!all(c("sample_accession", "study", "subgroup", "platform") %in% names(metadata_df))) {
    stop("Metadata file must include sample_accession, study, subgroup, and platform columns")
  }

  # Check that the number of samples matches from gene expression and metadata
  if (ncol(genex_df) != nrow(metadata_df)) {
    stop("Different number of samples in genex and metadata")
  }

  # Check that order of samples matches from gene expression and metadata
  if (!all(names(genex_df) == metadata_df$sample_accession)) {
    stop("Sample order does not match")
  }

}

calculate_confusion_matrix <- function(predicted_labels,
                                       true_labels,
                                       labels) {

  # Calculates confusion matrix and statistics about model performance
  #
  # Inputs
  #  predicted_labels: vector of predicted best guess labels for each sample
  #  true_labels: vector of true labels given for each sample
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #
  # Outputs
  #  confusionMatrix object, including confusion matrix and model stats

  confusion_matrix <- caret::confusionMatrix(data = factor(predicted_labels,
                                                           levels = labels),
                                             reference = factor(true_labels,
                                                                levels = labels),
                                             mode = "everything")

  return(confusion_matrix)

}

train_ktsp <- function(genex_df_train,
                       metadata_df_train,
                       model_seed = 2988,
                       ktsp_featureNo = 1000,
                       ktsp_n_rules_min = 5,
                       ktsp_n_rules_max = 50,
                       ktsp_weighted = TRUE) {

  # Train a kTSP model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  model_seed: seed used for reproducibility in training step (default: 2988)
  #  ktsp_featureNo: number of most informative features to filter down to (kTSP only) (default: 1000)
  #  ktsp_n_rules_min: minimum number of rules allowed for kTSP modeling (default: 5)
  #  ktsp_n_rules_max: maximum number of rules allowed for kTSP modeling (default: 50)
  #  ktsp_weighted: logical, if TRUE use one-vs-one and platform-wise comparisons to add weight to smaller subgroups and platforms
  #
  # Outputs
  #  kTSP classifier object
  #
  # Methodological choices
  #
  #  Filtering with multiclassPairs::filter_genes_TSP()
  #    - "UpDown" (TRUE) considers an equal number of up and down regulated genes
  #    - ktsp_weighted controls the following parameters:
  #      - if ktsp_weighted = TRUE (default)
  #        - "filter" set to "one_vs_one", which gives more weight to smaller classes
  #        - "platform_wise" set to TRUE, which helps select genes relevant to all platforms
  #      - if ktsp_weighted = FALSE
  #        - "filter" set to "one_vs_rest"
  #        - "platform_wise" set to FALSE, which selects most relevant genes without considering relevance across platforms
  #
  #  Training kTSP with multiclassPairs::train_one_vs_rest_TSP()
  #    - "include_pivot" (FALSE) means only filtered features are used to make rules
  #    - ktsp_weighted controls the following parameters:
  #      - if ktsp_weighted = TRUE (default)
  #        - "one_vs_one_scores" set to TRUE, rule scores calculated as mean of one vs one comparisons (giving more weight to smaller classes)
  #        - "platform_wise_scores" set to TRUE, rule scores calculated as mean of within-platform scores (which gives more weight to smaller platforms)
  #      - if ktsp_weighted = FALSE
  #        - "one_vs_one_scores" set to FALSE, rule scores calculated from one vs rest comparison
  #        - "platform_wise_scores" set to FALSE, rule scores calculated without respect to platform
  #
  # More information on multiclassPairs R package
  #  https://cran.r-project.org/web/packages/multiclassPairs/index.html
  #  https://cran.r-project.org/web/packages/multiclassPairs/vignettes/Tutorial.html

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_train,
                    metadata_df = metadata_df_train)

  # create data object
  train_data_object <- multiclassPairs::ReadData(Data = genex_df_train,
                                                 Labels = metadata_df_train$subgroup,
                                                 Platform = metadata_df_train$platform,
                                                 verbose = TRUE)

  # reduce genes to most useful features
  filtered_genes <- multiclassPairs::filter_genes_TSP(data_object = train_data_object,
                                                      filter = ifelse(ktsp_weighted, "one_vs_one", "one_vs_rest"),
                                                      platform_wise = ktsp_weighted,
                                                      featureNo = ktsp_featureNo,
                                                      UpDown = TRUE,
                                                      verbose = TRUE)

  # train kTSP model
  classifier <- multiclassPairs::train_one_vs_rest_TSP(data_object = train_data_object,
                                                       filtered_genes = filtered_genes,
                                                       k_range = ktsp_n_rules_min:ktsp_n_rules_max,
                                                       include_pivot = FALSE,
                                                       one_vs_one_scores = ktsp_weighted,
                                                       platform_wise_scores = ktsp_weighted,
                                                       seed = model_seed,
                                                       verbose = TRUE)

  return(classifier)

}

test_ktsp <- function(genex_df_test,
                      metadata_df_test,
                      classifier,
                      labels) {

  # Test a kTSP model
  #
  # Inputs
  #  genex_df_test: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_test: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  classifier: kTSP classifier produced by train_ktsp()
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #
  # Outputs
  #  List containing "predicted_labels" and "model_output" elements
  #    "predicted_labels" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method
  #
  # Methodological choices
  #
  #  Predictions with multiclassPairs::predict_one_vs_rest_TSP()
  #    - "tolerate_missed_genes" (TRUE) allows test data to be missing some features
  #    - "weighted_votes" (TRUE) weights more informative (higher scoring) rules more than others
  #
  # More information on multiclassPairs R package
  #  https://cran.r-project.org/web/packages/multiclassPairs/index.html
  #  https://cran.r-project.org/web/packages/multiclassPairs/vignettes/Tutorial.html

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)

  # create data object
  test_data_object <- multiclassPairs::ReadData(Data = genex_df_test,
                                                Labels = metadata_df_test$subgroup,
                                                Platform = metadata_df_test$platform,
                                                verbose = TRUE)

  # use classifier to predict labels
  test_results <- multiclassPairs::predict_one_vs_rest_TSP(classifier = classifier,
                                                           Data = test_data_object,
                                                           tolerate_missed_genes = TRUE,
                                                           weighted_votes = TRUE,
                                                           classes = labels,
                                                           verbose = TRUE)

  # create df with sample names and predicted labels
  predicted_labels_df <- dplyr::tibble(sample_accession = metadata_df_test$sample_accession,
                                       predicted_labels = test_results$max_score) # best guess

  # create output list with predicted labels and the modeling object
  test_results_list <- list(predicted_labels_df = predicted_labels_df,
                            model_output = test_results)

  return(test_results_list)

}

train_rf <- function(genex_df_train,
                     metadata_df_train,
                     model_seed = 4032,
                     rf_num.trees = 500,
                     rf_genes_altogether = 50,
                     rf_genes_one_vs_rest = 50,
                     rf_gene_repetition = 1,
                     rf_rules_altogether = 50,
                     rf_rules_one_vs_rest = 50,
                     rf_weighted = TRUE) {

  # Train a Random Forest model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  model_seed: seed used to generate additional modeling seeds for reproducibility
  #  rf_num.trees: number of trees used in RF modeling (use more trees given more features) (default: 500)
  #  rf_genes_altogether: number of top genes used when comparing all classes together (default: 50)
  #  rf_genes_one_vs_rest: number of top genes used when comparing each class against rest of classes (default: 50)
  #  rf_gene_repetition: number of times a gene can be used throughout set of rules (default: 1)
  #  rf_rules_altogether: number of top rules used when comparing all classes together (default: 50)
  #  rf_rules_one_vs_rest: number of top rules used when comparing each class against rest of classes (default: 50)
  #  rf_weighted: logical, if TRUE (default) use one-vs-rest and platform-wise comparisons to add weight to smaller subgroups and platforms
  #
  # Outputs
  #  Random Forest classifier object
  #
  # Methodological choices
  #
  #  Sorting genes with multiclassPairs::sort_genes_RF()
  #    - "rank_data" (TRUE) within each sample because we assume data are not normalized
  #    - rf_weighted controls the following parameters:
  #      - if rf_weighted = TRUE (default)
  #        - "featureNo_one_vs_rest" set to number of genes (return all genes, sorted)
  #        - "platform_wise" set to TRUE, which favors genes that are important across all platforms
  #      - if rf_weighted = FALSE
  #        - "featureNo_one_vs_rest" set to 0 (do not do one vs rest comparison)
  #        - "platform_wise" set to FALSE, which does not account platform in sorting genes
  #
  #  Sorting rules with multiclassPairs::sort_rules_RF()
  #    - rf_weighted controls the following parameters:
  #      - if rf_weighted = TRUE (default)
  #        - "genes_one_vs_rest" set by rf_genes_one_vs_rest parameter
  #        - "run_one_vs_rest" set to TRUE, which sorts rules based on importance within each class
  #        - "platform_wise" set to TRUE, which favors rules that are important on every platform
  #      - if rf_weighted = FALSE
  #        - "genes_one_vs_rest" set to 0 (do not do one vs rest comparison)
  #        - "run_one_vs_rest" set to FALSE, which sorts rules based on overall importance
  #        - "platform_wise" set to FALSE, which does not account for platform when sorting rules
  #
  #  Training RF model with multiclassPairs::train_RF()
  #    - "run_boruta" (TRUE) use Boruta algorithm to remove unimportant rules
  #    - "probability" (TRUE) allows test data to get scores for each class
  #    - rf_weighted controls the following parameters:
  #      - if rf_weighted = TRUE (default)
  #        - "rules_one_vs_rest" set by rf_rules_one_vs_rest parameter
  #      - if rf_weighted = FALSE
  #        - "rules_one_vs_rest" set to 0 (do not use rules from one vs rest comparison)
  #
  # More information on multiclassPairs R package
  #  https://cran.r-project.org/web/packages/multiclassPairs/index.html
  #  https://cran.r-project.org/web/packages/multiclassPairs/vignettes/Tutorial.html

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_train,
                    metadata_df = metadata_df_train)

  # get additional seeds for modeling steps
  set.seed(model_seed)
  seeds <- sample(1:1000, size = 2, replace = FALSE)

  # create data object
  train_data_object <- multiclassPairs::ReadData(Data = genex_df_train,
                                                 Labels = metadata_df_train$subgroup,
                                                 Platform = metadata_df_train$platform,
                                                 verbose = TRUE)

  # identify top genes used downstream
  n_complete_genes <- sum(complete.cases(train_data_object$data$Data))
  genes <- multiclassPairs::sort_genes_RF(data_object = train_data_object,
                                          featureNo_one_vs_rest = ifelse(rf_weighted, n_complete_genes, 0),
                                          rank_data = TRUE,
                                          platform_wise = rf_weighted,
                                          num.trees = rf_num.trees,
                                          seed = seeds[1],
                                          verbose = TRUE)

  # identify top gene pairs used downstream
  rules <- multiclassPairs::sort_rules_RF(data_object = train_data_object,
                                          sorted_genes_RF = genes,
                                          genes_altogether = rf_genes_altogether,
                                          genes_one_vs_rest = ifelse(rf_weighted, rf_genes_one_vs_rest, 0),
                                          run_altogether = TRUE, # default
                                          run_one_vs_rest = rf_weighted,
                                          platform_wise = rf_weighted,
                                          num.trees = rf_num.trees,
                                          seed = seeds[2],
                                          verbose = TRUE)

  # train RF model
  classifier <- multiclassPairs::train_RF(data_object = train_data_object,
                                          sorted_rules_RF = rules,
                                          gene_repetition = rf_gene_repetition,
                                          rules_altogether = rf_rules_altogether,
                                          rules_one_vs_rest = ifelse(rf_weighted, rf_rules_one_vs_rest, 0),
                                          run_boruta = TRUE,
                                          probability = TRUE,
                                          num.trees = rf_num.trees,
                                          verbose = TRUE)

  return(classifier)

}

test_rf <- function(genex_df_test,
                    metadata_df_test,
                    classifier) {

  # Test a Random Forest model
  #
  # Inputs
  #  genex_df_test: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_test: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  classifier: RF classifier produced by train_rf()
  #
  # Outputs
  #  List containing "predicted_labels" and "model_output" elements
  #    "predicted_labels" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method
  #
  # Methodological choices
  #
  #  Predicting labels with multiclassPairs::predict_RF()
  #    - "impute" (TRUE) allows test data to be missing features by imputing with kNN
  #
  # More information on multiclassPairs R package
  #  https://cran.r-project.org/web/packages/multiclassPairs/index.html
  #  https://cran.r-project.org/web/packages/multiclassPairs/vignettes/Tutorial.html

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)

  # create data object
  test_data_object <- multiclassPairs::ReadData(Data = genex_df_test,
                                                Labels = metadata_df_test$subgroup,
                                                Platform = metadata_df_test$platform,
                                                verbose = TRUE)

  # predict labels with classifier
  test_results <- multiclassPairs::predict_RF(classifier = classifier,
                                              Data = test_data_object,
                                              impute = TRUE)

  # get the prediction matrix
  test_pred <- test_results$predictions

  # pick the column with maximum probability
  set.seed(1) # max.col randomly uses random selection in case of ties
  test_prediction_labels <- colnames(test_pred)[max.col(test_pred)]

  # create df with sample names and predicted labels
  predicted_labels_df <- dplyr::tibble(sample_accession = metadata_df_test$sample_accession,
                                       predicted_labels = test_prediction_labels) # best guess

  # create output list with predicted labels and the modeling object
  test_results_list <- list(predicted_labels_df = predicted_labels_df,
                            model_output = test_pred) # test_results is too much info

  return(test_results_list)

}

test_mm2s <- function(genex_df_test,
                      metadata_df_test,
                      model_seed = 4418,
                      ah_date = "2022-10-30") {

  # Test an MM2S model
  #
  # Inputs
  #  genex_df_test: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_test: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  model_seed: seed used for reproducibility in MM2S.human function
  #  ah_date: AnnotationHub snapshot date (default: 2022-10-30)
  #
  # Outputs
  #  List containing "predicted_labels" and "model_output" elements
  #    "predicted_labels" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)

  # load MM2S library (necessary for MM2S.human() function to work)
  library(MM2S)

  # convert genex_df gene names to from ENSEMBL to ENTREZID
  genex_df_test_ENTREZID <- genex_df_test |>
    tibble::rownames_to_column(var = "ENSEMBL") |>
    convert_gene_names(gene_column_before = "ENSEMBL",
                       gene_column_after = "gene",
                       map_from = "ENSEMBL",
                       map_to = "ENTREZID",
                       ah_date = ah_date) |>
    tibble::column_to_rownames(var = "gene")

  # predict labels with existing MM2S.human classifier
  mm2s_predictions <- MM2S::MM2S.human(InputMatrix = genex_df_test_ENTREZID,
                                       parallelize = 1,
                                       seed = model_seed)

  # give MM2S_Subtype data frame structure if only one test sample
  if(ncol(genex_df_test) == 1) {

    mm2s_predictions$MM2S_Subtype <- data.frame(as.list(mm2s_predictions$MM2S_Subtype))

  }

  # modify MM2S predictions to fit this project's medulloblastoma subgroup names
  test_results <- dplyr::bind_cols(mm2s_predictions$MM2S_Subtype,
                                   mm2s_predictions$Predictions) |>
    dplyr::mutate(MM2S_Prediction = dplyr::case_when(MM2S_Prediction == "Group3" ~ "G3",
                                                     MM2S_Prediction == "Group4" ~ "G4",
                                                     TRUE ~ MM2S_Prediction)) |>
    dplyr::select(SampleName,
                  MM2S_Prediction,
                  "G3" = Group3,
                  "G4" = Group4,
                  Normal,
                  SHH,
                  WNT)

  # create df with sample names and predicted labels
  predicted_labels_df <- dplyr::tibble(sample_accession = metadata_df_test$sample_accession,
                                       predicted_labels = test_results$MM2S_Prediction) # best guess

  # create output list with predicted labels and the modeling object
  test_results_list <- list(predicted_labels_df = predicted_labels_df,
                            model_output = test_results)

  return(test_results_list)

}

train_lasso <- function(genex_df_train,
                        metadata_df_train,
                        model_seed = 2064) {

  # Train a LASSO model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  model_seed: seed used for reproducibility in training step
  #
  # Outputs
  #  LASSO classifier object

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_train,
                    metadata_df = metadata_df_train)

  # do basic normalization: make each column sums to 1
  genex_df_train <- apply(genex_df_train, 2, function(x) x/sum(x))


  # set for reproducibility
  set.seed(model_seed)


  # train glmnet model with alpha = 1 for lasso
  classifier <- glmnet::cv.glmnet(x = t(genex_df_train),
                                  y = metadata_df_train$subgroup,
                                  family = "multinomial",
                                  type.measure = "class",
                                  alpha = 1) # lasso

  return(classifier)

}

test_lasso <- function(genex_df_test,
                       metadata_df_test,
                       classifier) {

  # Test a LASSO model
  #
  # Inputs
  #  genex_df_test: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_test: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  classifier: LASSO classifier produced by train_lasso()
  #
  # Outputs
  #  List containing "predicted_labels" and "model_output" elements
  #    "predicted_labels" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method

  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)

  # load glmnet library (necessary for predict() function to work on cv.glmnet object)
  library(glmnet)

  # do basic normalization: make each column sums to 1
  genex_df_test <- apply(genex_df_test, 2, function(x) x/sum(x))

  # get model subgroup labels
  lasso_subgroups <- row.names(classifier$glmnet.fit$a0)

  # predict using LASSO classifier
  test_results <- predict(classifier,
                          t(genex_df_test),
                          s = classifier$lambda.1se,
                          type = "response") |>
    as.data.frame() |>
    setNames(lasso_subgroups) |>
    tibble::rownames_to_column(var = "sample_accession") |>
    tibble::as_tibble()

  # add top scoring prediction column
  test_results <- test_results |>
    dplyr::mutate(prediction = lasso_subgroups[max.col(m = test_results[,-1])])

  # create df with sample names and predicted labels
  predicted_labels_df <- dplyr::tibble(sample_accession = metadata_df_test$sample_accession,
                                       predicted_labels = test_results$prediction) # best guess

  # create output list with predicted labels and the modeling object
  test_results_list <- list(predicted_labels_df = predicted_labels_df,
                            model_output = test_results)

  return(test_results_list)

}

return_model_metrics <- function(single_repeat,
                                 model_types,
                                 metadata_df,
                                 labels,
                                 platforms = NULL,
                                 studies = NULL) {

  # Returns the Kappa and Balanced Accuracy metrics from a set of models
  # associated with a single repeat (e.g., one repeat from run_many_models()).
  # Use this function inside purrr::map() with entire list output from run_many_models().
  #
  # Inputs
  #  single_repeat: single element of a models list corresponding to one repeat
  #    single_repeat is a list with one slot per model_type in addition to
  #    slots with metainfo about the repeat
  #  model_types: vector of model types to be assessed (e.g. c("ktsp", "rf)) (must be given)
  #  metadata_df: metadata data frame (must include sample_accession, study, subgroup, and platform columns)
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #  platforms: vector of platforms to be assessed separately (e.g., c("Array", "RNA-seq"))
  #    if platforms = NULL, results are presented in aggregate
  #  studies: vector of studies to be assessed separately (e.g., c("GSE85217", "OpenPBTA"))
  #    if studies = NULL, results are presented in aggregate
  #
  # Outputs
  #  Data frame of model metrics stratified by model types, platforms, and studies

  # check that all model_types are present in the repeat
  if (!all(model_types %in% names(single_repeat))) {

    stop("In return_model_metrics(), some model types missing from model list.")

  }

  # define the eventual output df
  metrics_df <- NULL

  # Loop over model_types
  for (model_type in model_types) {

    # Check that metadata exists for all test samples
    if (all(single_repeat[[model_type]]$test_results$predicted_labels_df$sample_accession %in%
            metadata_df$sample_accession)) {

      # merge test results with sample metadata
      df <- single_repeat[[model_type]]$test_results$predicted_labels_df |>
        dplyr::left_join(metadata_df,
                         by = "sample_accession")

    } else {

      stop("Sample accessions missing from metadata_df in return_model_metrics().")

    }

    # if no platforms or studies are specified, analyze everything together
    if (is.null(platforms) & is.null(studies)) {

      df <- df |>
        dplyr::mutate(platform_group = stringr::str_c(unique(platform), collapse = ","),
                      study_group = stringr::str_c(unique(study), collapse = ","))

    } else if (is.null(platforms) & !is.null(studies)) {

      df <- df |>
        dplyr::filter(study %in% studies) |>
        dplyr::group_by(study) |>
        dplyr::mutate(platform_group = stringr::str_c(unique(platform), collapse = ","),
                      study_group = study)

    } else if (!is.null(platforms) & is.null(studies)) {

      df <- df |>
        dplyr::filter(platform %in% platforms) |>
        dplyr::group_by(platform) |>
        dplyr::mutate(platform_group = platform,
                      study_group = stringr::str_c(unique(study), collapse = ","))

    } else {

      df <- df |>
        dplyr::filter(platform %in% platforms,
                      study %in% studies) |>
        dplyr::mutate(platform_group = platform,
                      study_group = study)

    }

    # create a list of data frames defined by platform and study group
    model_type_df <- split(df,
                           f = list(df$platform_group,
                                    df$study_group),
                           drop = TRUE) |>
      # map each data frame to calculate_model_metrics
      purrr::map(\(x) calculate_model_metrics(df = x,
                                              labels = labels)) |>
      # names of returned list are a combination of platform and study
      purrr::list_rbind(names_to = "platform_study") |>
      tidyr::separate(platform_study, into = c("platform", "study"), sep = "\\.") |>
      # add model type and official model status as metainfo
      tibble::add_column(model_type = model_type,
                         .before = "platform") |>
      tibble::add_column(official_model = single_repeat[["official_model"]],
                         .before = "platform")

    metrics_df <- dplyr::bind_rows(metrics_df, model_type_df)

  }

  return(metrics_df)

}

calculate_model_metrics <- function(df, labels) {

  # Calculate the Kappa and Balanced Accuracy metrics from a set of predictions
  #
  # Inputs
  #  df: data frame with "predicted_labels" and true "subgroup" labels
  #  labels: vector of possible sample labels (e.g., c("G3","G4","SHH","WNT"))
  #
  # Outputs
  #  Data frame of model metrics

  total_n <- nrow(df)

  # if at least one sample matches platform/study criteria
  if (total_n >= 1) {

    # confusion matrix
    cm <- calculate_confusion_matrix(predicted_labels = df$predicted_labels,
                                     true_labels = df$subgroup,
                                     labels = labels)

    # number of test samples from each subgroup
    n_subgroup_samples <- purrr::map_dbl(labels,
                                         \(x) sum(df$subgroup == x))

    # combine Kappa (one Overall value) with Balanced Accuracy (one value per subgroup)
    return_df <- tibble::tibble(total_samples = total_n,
                                subgroup = c("Overall", labels),
                                subgroup_samples = c(total_n, n_subgroup_samples),
                                metric = c("Kappa", rep("Balanced Accuracy", length(labels))),
                                value = c(cm$overall[["Kappa"]], cm$byClass[,"Balanced Accuracy"]))

    return(return_df)

  } else {

    return(NULL)

  }
}
