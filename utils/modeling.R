check_input_files <- function(genex_df,
                              metatdata_df) {
  
  # Checks input files when performing modeling training and testin
  #
  # Inputs
  #  genex_df: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #
  # Output
  #  None (functions calls stop if anything is wrong)
  
  # Check that metadata_df has all necessary columns
  if (!all(c("sample_accesion", "subgroup", "platform") %in% names(metadata_df))) {
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
  
  # Create confusion matrix
  #
  # Inputs
  #  predicted_labels: vector of best guess labels for each sample
  #  true_labels: vector of true labels given for each sample
  #  labels: vector of possible label values that could exist
  #
  # Outputs
  #  Confusion matrix
  
  confusion_matrix <- caret::confusionMatrix(data = factor(predicted_labels, 
                                                           levels = labels),
                                             reference = factor(true_labels,
                                                                levels = labels),
                                             mode = "everything")
  
  return(confusion_matrix)
  
}

train_ktsp <- function(genex_df_train,
                       metadata_df_train,
                       model_seed,
                       n_rules_min,
                       n_rules_max) {
  
  # Train a kTSP model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #  model_seed: seed used for reproducibility in training step
  #  n_rules_min: minimum number of rules allowed for kTSP modeling
  #  n_rules_max: maximum number of rules allowed for kTSP modeling
  #
  # Outputs
  #  kTSP classifier object
  
  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_train,
                    metadata_df = metadata_df_train)
  
  train_data_object <- multiclassPairs::ReadData(Data = genex_df_train,
                                                 Labels = metadata_df_train$subgroup,
                                                 Platform = metadata_df_train$platform,
                                                 verbose = TRUE)
  
  filtered_genes <- multiclassPairs::filter_genes_TSP(data_object = train_data_object,
                                                      filter = "one_vs_one",
                                                      platform_wise = TRUE,
                                                      featureNo = 1000,
                                                      UpDown = TRUE,
                                                      verbose = TRUE)
  
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
  #  List containing "predicted_label" and "model_output" elements
  #    "predicted_label" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method
  
  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)
  
  test_data_object <- multiclassPairs::ReadData(Data = genex_df_test,
                                                Labels = metadata_df_test$subgroup,
                                                Platform = metadata_df_test$platform,
                                                verbose = TRUE)
  
  test_results <- multiclassPairs::predict_one_vs_rest_TSP(classifier = classifier,
                                                           Data = test_data_object,
                                                           tolerate_missed_genes = TRUE,
                                                           weighted_votes = TRUE,
                                                           classes = labels,
                                                           verbose = TRUE)

  predicted_label_df <- dplyr::tibble(sample_accesion = metadata_df_test$sample_accession,
                                      predicted_label = test_results$max_score) # best guess
    
  test_results_list <- list(predicted_label_df = predicted_label_df,
                            model_output = test_results)
  
  return(test_results_list)
  
}

train_rf <- function(genex_df_train,
                     metadata_df_train,
                     model_seed) {
  
  # Train a Random Forest model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (genes as row names and one column per sample)
  #  metadata_df_train: metadata data frame (must include sample_accession, subgroup, and platform columns)
  #  model_seed: seed used to generate additional modeling seeds for reproducibility
  #
  # Outputs
  #  Random Forest classifier object
  
  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_train,
                    metadata_df = metadata_df_train)
  
  # get additional seeds for modeling steps
  n_seeds <- 2
  set.seed(model_seed)
  seeds <- sample(1:1000, size = n_seeds, replace = FALSE)
  
  train_data_object <- multiclassPairs::ReadData(Data = genex_df_train,
                                                 Labels = metadata_df_train$subgroup,
                                                 Platform = metadata_df_train$platform,
                                                 verbose = TRUE)
  
  genes <- multiclassPairs::sort_genes_RF(data_object = train_data_object,
                                          rank_data = TRUE,
                                          platform_wise = TRUE,
                                          num.trees = 1000,
                                          seed = seeds[1],
                                          verbose = TRUE)
  
  rules <- multiclassPairs::sort_rules_RF(data_object = train_data_object, 
                                          sorted_genes_RF = genes,
                                          genes_altogether = 50,
                                          genes_one_vs_rest = 50,
                                          platform_wise = TRUE,
                                          num.trees = 1000,
                                          seed = seeds[2],
                                          verbose = TRUE)
  
  classifier <- multiclassPairs::train_RF(data_object = train_data_object,
                                          sorted_rules_RF = rules,
                                          gene_repetition = 1,
                                          rules_altogether = 50,
                                          rules_one_vs_rest = 50,
                                          run_boruta = TRUE,
                                          plot_boruta = FALSE,
                                          probability = TRUE,
                                          num.trees = 1000,
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
  #  classifier: kTSP classifier produced by train_ktsp()
  #
  # Outputs
  #  List containing "predicted_label" and "model_output" elements
  #    "predicted_label" contains a data frame with one row for each sample and its predicted label
  #    "model_output" is the prediction object returned by this method
  
  # ensure input files are properly formatted and sample orders match
  check_input_files(genex_df = genex_df_test,
                    metadata_df = metadata_df_test)
  
  test_data_object <- multiclassPairs::ReadData(Data = genex_df_test,
                                                Labels = metadata_df_test$subgroup,
                                                Platform = metadata_df_test$platform,
                                                verbose = TRUE)
  
  test_results <- multiclassPairs::predict_RF(classifier = classifier, 
                                              Data = test_data_object)
  
  # get the prediction matrix
  test_pred <- test_results$predictions
  
  # pick the column with maximum probability
  test_prediction_labels <- colnames(test_pred)[max.col(test_pred)]
  
  predicted_label_df <- dplyr::tibble(sample_accesion = metadata_df_test$sample_accession,
                                      predicted_label = test_prediction_labels) # best guess
  
  test_results_list <- list(predicted_label_df = predicted_label_df,
                            model_output = test_pred) # test_results is too much info
  
  return(test_results_list)
  
}

run_mm2s <- function(genex_df_test,
                     metadata_df_test,
                     model_seed,
                     gene_map_df) {
  
  # Run an MM2S model
  #
  # Inputs
  #  metadata_df_train: metadata data frame (train)
  #  metadata_df_test: metadata data frame (test)
  #  model_seed: seed re-used in each modeling step
  #  gene_map_df: gene map used to convert MM2S gene names
  #
  # Outputs
  #  List including classifier, test results, and test confusion matrix
  
  mb_subgroups <- c("G3", "G4", "NORMAL", "SHH", "WNT") # MM2S predicts "NORMAL" too
  
  set.seed(model_seed)
  
  genex_df_test_ENTREZID <- genex_df_test %>%
    tibble::rownames_to_column(var = "ENSEMBL") %>%
    dplyr::left_join(gene_map_df %>% dplyr::select(ENSEMBL, ENTREZID),
                     by = "ENSEMBL") %>%
    dplyr::filter(!duplicated(ENSEMBL),
                  !duplicated(ENTREZID),
                  !is.na(ENTREZID)) %>%
    dplyr::select(-ENSEMBL) %>%
    tibble::column_to_rownames(var = "ENTREZID")
  
  mm2s_predictions <- MM2S::MM2S.human(InputMatrix = genex_df_test_ENTREZID,
                                       parallelize = 1,
                                       seed = model_seed)
  
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
  
  # test accuracy
  test_cm <- caret::confusionMatrix(data = factor(test_results$MM2S_Prediction,
                                                  levels = mb_subgroups),
                                    reference = factor(metadata_df_test$subgroup,
                                                       levels = mb_subgroups),
                                    mode = "everything")
  
  list(classifier = NA,
       test_results = test_results,
       test_cm = test_cm)
  
}

run_lasso <- function(genex_df_train,
                      genex_df_test,
                      metadata_df_train,
                      metadata_df_test,
                      model_seed) {
  
  # Run a LASSO model
  #
  # Inputs
  #  genex_df_train: gene expression matrix (train)
  #  genex_df_test: gene expression matrix (test)
  #  metadata_df_train: metadata data frame (train)
  #  metadata_df_test: metadata data frame (test)
  #  model_seed: seed re-used in each modeling step
  #
  # Outputs
  #  List including classifier, test results, and test confusion matrix
  
  mb_subgroups <- c("G3", "G4", "SHH", "WNT")
  
  set.seed(model_seed)
  
  # do basic normalization: make each column sums to 1
  
  genex_df_train <- apply(genex_df_train, 2, function(x) x/sum(x))
  genex_df_test <- apply(genex_df_test, 2, function(x) x/sum(x))
  
  classifier <- glmnet::cv.glmnet(x = t(genex_df_train),
                                  y = metadata_df_train$subgroup,
                                  family = "multinomial",
                                  type.measure = "class",
                                  alpha = 1) # lasso
  
  test_results <- predict(classifier,
                          t(genex_df_test),
                          s = classifier$lambda.1se,
                          type = "response")[,,1] %>%
    as.data.frame() %>%
    dplyr::mutate(prediction = names(.)[max.col(.)]) %>%
    tibble::rownames_to_column(var = "sample_accession") %>%
    tibble::as_tibble()
  
  test_cm <- caret::confusionMatrix(factor(test_results$prediction,
                                           levels = mb_subgroups),
                                    factor(metadata_df_test$subgroup,
                                           levels = mb_subgroups),
                                    mode = "everything")
  
  list(classifier = classifier,
       test_results = test_results,
       test_cm = test_cm)
  
}
