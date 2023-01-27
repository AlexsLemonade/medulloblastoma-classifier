check_input_files <- function(genex_df,
                              metadata_df) {
  
  # Checks input files when performing modeling training and testing
  #
  # Inputs
  #  genex_df: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #
  # Output
  #  None (function calls 'stop' if anything is wrong)
  
  # Check that metadata_df has all necessary columns
  if (!all(c("sample_accession", "subgroup", "platform") %in% names(metadata_df))) {
    stop("Metadata file must include sample_accession, subgroup, and platform columns")
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
                       n_rules_min = 5,
                       n_rules_max = 50) {
  
  # Train a kTSP model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #  model_seed: seed used for reproducibility in training step (default: 2988)
  #  ktsp_featureNo: number of most informative features to filter down to (kTSP only) (default: 1000)
  #  n_rules_min: minimum number of rules allowed for kTSP modeling (default: 5)
  #  n_rules_max: maximum number of rules allowed for kTSP modeling (default: 50)
  #
  # Outputs
  #  kTSP classifier object
  #
  # Methodological choices
  #
  #  Filtering with multiclassPairs::filter_genes_TSP()
  #    - "one_vs_one" filtering give more weight to smaller classes
  #    - "platform_wise" (TRUE) filtering helps select genes relevant to each platform
  #    - "UpDown" (TRUE) considers an equal number of up and down regulated genes
  #
  #  Training kTSP with multiclassPairs::train_one_vs_rest_TSP()
  #    - "include_pivot" (FALSE) means only filtered features are used to make rules
  #    - "one_vs_one_scores" (TRUE) gives more weight to small classes
  #    - "platform_wise_scores" (TRUE) gives more weight to small platforms
  
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
                                                      filter = "one_vs_one",
                                                      platform_wise = TRUE,
                                                      featureNo = ktsp_featureNo,
                                                      UpDown = TRUE,
                                                      verbose = TRUE)
  
  # train kTSP model
  classifier <- multiclassPairs::train_one_vs_rest_TSP(data_object = train_data_object,
                                                       filtered_genes = filtered_genes,
                                                       k_range = n_rules_min:n_rules_max,
                                                       include_pivot = FALSE,
                                                       one_vs_one_scores = TRUE,
                                                       platform_wise_scores = TRUE,
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
  #  metadata_df_test: metadata data frame (must include sample_accession, subgroup, and platform columns)
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
                     rf_rules_one_vs_rest = 50) {
  
  # Train a Random Forest model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #  model_seed: seed used to generate additional modeling seeds for reproducibility
  #  rf_num.trees: number of trees used in RF modeling (use more trees given more features)
  #  rf_genes_altogether: number of top genes used when comparing all classes together (default: 50)
  #  rf_genes_one_vs_rest: number of top genes used when comparing each class against rest of classes (default: 50)
  #  rf_gene_repetition: number of time a gene can be used throughout set of rules
  #  rf_rules_altogether: number of top rules used when comparing all classes together (default: 50)
  #  rf_rules_one_vs_rest: number of top rules used when comparing each class against rest of classes (default: 50)
  #
  # Outputs
  #  Random Forest classifier object
  #
  # Methodological choices
  #
  #  Sorting genes with multiclassPairs::sort_genes_RF()
  #    - "rank_data" (TRUE) within each sample because we assume data are not normalized
  #    - "platform_wise" (TRUE) favors genes that are important on every platform
  #
  #  Sorting rules with multiclassPairs::sort_rules_RF()
  #    - "platform_wise" (TRUE) favors rules that are important on every platform
  #
  #  Training RF model with multiclassPairs::train_RF()
  #    - "run_boruta" (TRUE) use Boruta algorithm to remove unimportant rules
  #    - "probability" (TRUE) allows test data to get scores for each class
  
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
  genes <- multiclassPairs::sort_genes_RF(data_object = train_data_object,
                                          rank_data = TRUE,
                                          platform_wise = TRUE,
                                          num.trees = rf_num.trees,
                                          seed = seeds[1],
                                          verbose = TRUE)
  
  # identify top gene pairs used downstream
  rules <- multiclassPairs::sort_rules_RF(data_object = train_data_object, 
                                          sorted_genes_RF = genes,
                                          genes_altogether = rf_genes_altogether,
                                          genes_one_vs_rest = rf_genes_one_vs_rest,
                                          platform_wise = TRUE,
                                          num.trees = rf_num.trees,
                                          seed = seeds[2],
                                          verbose = TRUE)
  
  # train RF model
  classifier <- multiclassPairs::train_RF(data_object = train_data_object,
                                          sorted_rules_RF = rules,
                                          gene_repetition = rf_gene_repetition,
                                          rules_altogether = rf_rules_altogether,
                                          rules_one_vs_rest = rf_rules_one_vs_rest,
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
  #  metadata_df_test: metadata data frame (must include sample_accession, subgroup, and platform columns)
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
                      gene_map_df) {
  
  # Test an MM2S model
  #
  # Inputs
  #  genex_df_test: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_test: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #  model_seed: seed used for reproducibility in MM2S.human function
  #  gene_map_df: gene map used to convert ENSEMBL IDs to ENTREZID to match MM2S model
  #
  # Outputs
  #  List containing "predicted_labels" and "model_output" elements
  #    "predicted_labels" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method
  
  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)
  
  # convert genex_df gene names to from ENSEMBL to ENTREZID
  genex_df_test_ENTREZID <- genex_df_test %>%
    tibble::rownames_to_column(var = "ENSEMBL") %>%
    dplyr::left_join(gene_map_df %>% dplyr::select(ENSEMBL, ENTREZID),
                     by = "ENSEMBL") %>%
    dplyr::filter(!duplicated(ENSEMBL),
                  !duplicated(ENTREZID),
                  !is.na(ENTREZID)) %>%
    dplyr::select(-ENSEMBL) %>%
    tibble::column_to_rownames(var = "ENTREZID")
  
  # predict labels with existing MM2S.human classifier
  mm2s_predictions <- MM2S::MM2S.human(InputMatrix = genex_df_test_ENTREZID,
                                       parallelize = 1,
                                       seed = model_seed)

  # modify MM2S predictions to fit this project's medulloblastoma subgroup names  
  test_results <- dplyr::bind_cols(mm2s_predictions$MM2S_Subtype,
                                   mm2s_predictions$Predictions) %>%
    dplyr::mutate(MM2S_Prediction = dplyr::case_when(MM2S_Prediction == "Group3" ~ "G3",
                                                     MM2S_Prediction == "Group4" ~ "G4",
                                                     TRUE ~ MM2S_Prediction)) %>%
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
  #  metadata_df_train: metadata data frame (must include sample_accession, subgroup, and platform columns)
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
  #  metadata_df_test: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #  classifier: LASSO classifier produced by train_lasso()
  #
  # Outputs
  #  List containing "predicted_labels" and "model_output" elements
  #    "predicted_labels" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method
  
  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)
  
  # do basic normalization: make each column sums to 1
  genex_df_test <- apply(genex_df_test, 2, function(x) x/sum(x))
  
  # predict using LASSO classifier
  test_results <- predict(classifier,
                          t(genex_df_test),
                          s = classifier$lambda.1se,
                          type = "response")[,,1] %>%
    as.data.frame() %>%
    dplyr::mutate(prediction = names(.)[max.col(.)]) %>%
    tibble::rownames_to_column(var = "sample_accession") %>%
    tibble::as_tibble()

  # create df with sample names and predicted labels
  predicted_labels_df <- dplyr::tibble(sample_accession = metadata_df_test$sample_accession,
                                       predicted_labels = test_results$prediction) # best guess
  
  # create output list with predicted labels and the modeling object
  test_results_list <- list(predicted_labels_df = predicted_labels_df,
                            model_output = test_results)
  
  return(test_results_list)
  
}
