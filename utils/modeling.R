suppressMessages(library(foreach))

run_many_models <- function(genex_df,
                            metadata_df,
                            model_types = c("ktsp", "rf", "mm2s", "lasso"),
                            initial_seed = 44,
                            n_repeats = 1,
                            n_cores = 1,
                            n_rules_min = 5,
                            n_rules_max = NA) {
  
  # model types should be a list with elements limited to "ktsp", "rf", "mm2s", "lasso"
  if (!is.vector(model_types) |
      !all(model_types %in% c("ktsp", "rf", "mm2s", "lasso"))) {
    
    stop("model_types in run_models() should be a vector limited to 'ktsp', 'rf', 'mm2s', 'lasso'.")
    
  }
  
  # n_repeats should be a positive integer <= 100
  if (n_repeats < 1 | n_repeats > 100 | round(n_repeats) != n_repeats) {
    
    stop("n_repeats in run_models() should be a positive integer <= 100.")
    
  }
  
  # n_cores should not exceed parallel::detectCores() - 1
  n_cores <- min(n_cores, parallel::detectCores() - 1)
  
  # n_rules_min should be a positive integer
  
  # n_rules_max should be a positive integer and >= n_rules_min
  # if n_rules_max is not given, set to n_rules_min
  if (is.na(n_rules_max)) {
    n_rules_max <- n_rules_min
  }
  
  set.seed(initial_seed)
  
  train_test_seeds <- sample(1:max(100, n_repeats), size = n_repeats)
  modeling_seeds <- sample(1:max(100, n_repeats), size = n_repeats)
  official_model_n <- sample(1:n_repeats, size = 1)
  
  # parallel backend
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,
                          c("get_train_test_samples",
                            "run_one_model",
                            "run_ktsp")) #,
                            #"run_rf",
                            #"run_mm2s",
                            #"run_lasso"))
  
  model_list <- foreach(n = 1:n_repeats) %dopar% {
    
    suppressMessages(library(magrittr))
    
    train_test_samples_list <- get_train_test_samples(genex_df,
                                                      metadata_df,
                                                      train_test_seed = train_test_seeds[n])
    
    genex_df_train <- genex_df %>%
      dplyr::select(train_test_samples_list$train)
    genex_df_test <- genex_df %>%
      dplyr::select(train_test_samples_list$test)
    metadata_df_train <- metadata_df %>%
      dplyr::filter(sample_accession %in% train_test_samples_list$train)
    metadata_df_test <- metadata_df %>%
      dplyr::filter(sample_accession %in% train_test_samples_list$test)
    
    repeat_list <- purrr::map(model_types,
                              function(x) run_one_model(x,
                                                        genex_df_train,
                                                        genex_df_test,
                                                        metadata_df_train,
                                                        metadata_df_test,
                                                        modeling_seeds[n],
                                                        n_rules_min,
                                                        n_rules_max))
    
    names(repeat_list) <- model_types
    
    repeat_list[["train_test_seed"]] <- train_test_seeds[n]
    repeat_list[["modeling_seed"]] <- modeling_seeds[n]
    
    return(repeat_list)
    
  }
  
  # stop parallel backend
  parallel::stopCluster(cl)
  
  # set official model
  model_list[[official_model_n]][["official_model"]] <- TRUE
  
  return(model_list)
  
}

get_train_test_samples <- function(genex_df,
                                   metadata_df,
                                   train_test_seed) {
  
  proportion_of_studies_train <- 0.5
  
  set.seed(train_test_seed)
  
  array_studies <- metadata_df %>%
    dplyr::filter(sample_accession %in% names(genex_df),
                  platform == "Array") %>%
    dplyr::pull(study) %>%
    unique()
  
  rnaseq_studies <- metadata_df %>%
    dplyr::filter(sample_accession %in% names(genex_df),
                  platform == "RNA-seq") %>%
    dplyr::pull(study) %>%
    unique()
  
  n_array_studies <- length(array_studies)
  n_rnaseq_studies <- length(rnaseq_studies)
  
  n_array_studies_train <- ceiling(n_array_studies*proportion_of_studies_train)
  n_rnaseq_studies_train <- ceiling(n_rnaseq_studies*proportion_of_studies_train)
  
  array_studies_train <- sample(array_studies, size = n_array_studies_train)
  rnaseq_studies_train <- sample(rnaseq_studies, size = n_rnaseq_studies_train)
  
  array_studies_test <- setdiff(array_studies, array_studies_train)
  rnaseq_studies_test <- setdiff(rnaseq_studies, rnaseq_studies_train)
  
  train_test_samples_list <- list()
  
  train_test_samples_list[["train"]] <- metadata_df %>%
    dplyr::filter(study %in% c(array_studies_train,
                               rnaseq_studies_train)) %>%
    dplyr::pull(sample_accession)
  
  train_test_samples_list[["test"]] <- metadata_df %>%
    dplyr::filter(study %in% c(array_studies_test,
                               rnaseq_studies_test)) %>%
    dplyr::pull(sample_accession)
  
  return(train_test_samples_list)
  
}

run_one_model <- function(type,
                          genex_df_train,
                          genex_df_test,
                          metadata_df_train,
                          metadata_df_test,
                          model_seed,
                          n_rules_min,
                          n_rules_max) {
  
  if (type == "ktsp") {
    
    model <- run_ktsp(genex_df_train,
                      genex_df_test,
                      metadata_df_train,
                      metadata_df_test,
                      model_seed,
                      n_rules_min,
                      n_rules_max)
    
  } else if (type == "rf") {
    
    model <- run_rf(genex_df_train,
                    genex_df_test,
                    metadata_df_train,
                    metadata_df_test,
                    model_seed,
                    n_rules_min,
                    n_rules_max)
    
  } else if (type == "mm2s") {
    
    model <- run_mm2s(genex_df_test,
                      metadata_df_test,
                      model_seed)
    
  } else if (type == "lasso") {
    
    model <- run_lasso(genex_df_train,
                       genex_df_test,
                       metadata_df_train,
                       metadata_df_test,
                       model_seed)
    
  } else {
    
    stop("Type of model must be ktsp, rf, mm2s, or lasso.")
    
  }
  
  return(model)
  
}


run_ktsp <- function(genex_df_train,
                     genex_df_test,
                     metadata_df_train,
                     metadata_df_test,
                     model_seed,
                     n_rules_min,
                     n_rules_max) {

  mb_subgroups <- c("G3", "G4", "SHH", "WNT")
  
  set.seed(model_seed)
  
  train_data_object <- multiclassPairs::ReadData(Data = genex_df_train,
                                                 Labels = metadata_df_train$subgroup,
                                                 Platform = metadata_df_train$platform,
                                                 verbose = TRUE)
  
  test_data_object <- multiclassPairs::ReadData(Data = genex_df_test,
                                                Labels = metadata_df_test$subgroup,
                                                Platform = metadata_df_test$platform,
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
  
  train_results <- multiclassPairs::predict_one_vs_rest_TSP(classifier = classifier,
                                                            Data = train_data_object,
                                                            tolerate_missed_genes = TRUE,
                                                            weighted_votes = TRUE,
                                                            classes = mb_subgroups,
                                                            verbose = TRUE)
  
  test_results <- multiclassPairs::predict_one_vs_rest_TSP(classifier = classifier,
                                                           Data = test_data_object,
                                                           tolerate_missed_genes = TRUE,
                                                           weighted_votes = TRUE,
                                                           classes = mb_subgroups,
                                                           verbose = TRUE)
  
  train_cm <- caret::confusionMatrix(data = factor(train_results$max_score, 
                                                   levels = mb_subgroups),
                                     reference = factor(train_data_object$data$Labels, 
                                                        levels = mb_subgroups),
                                     mode = "everything")
  
  test_cm <- caret::confusionMatrix(data = factor(test_results$max_score, 
                                                  levels = mb_subgroups),
                                    reference = factor(test_data_object$data$Labels, 
                                                       levels = mb_subgroups),
                                    mode = "everything")
  
  list(filtered_genes = filtered_genes,
       classifier = classifier,
       train_results = train_results,
       test_results = test_results,
       train_cm = train_cm,
       test_cm = test_cm)
  
}

run_rf <- function(genex_df_train,
                   genex_df_test,
                   metadata_df_train,
                   metadata_df_test,
                   model_seed,
                   n_rules_min,
                   n_rules_max) {

  mb_subgroups <- c("G3", "G4", "SHH", "WNT")
  
  set.seed(model_seed)
  
  train_data_object <- multiclassPairs::ReadData(Data = genex_df_train,
                                                 Labels = metadata_df_train$subgroup,
                                                 Platform = metadata_df_train$platform,
                                                 verbose = TRUE)
  
  test_data_object <- multiclassPairs::ReadData(Data = genex_df_test,
                                                Labels = metadata_df_test$subgroup,
                                                Platform = metadata_df_test$platform,
                                                verbose = TRUE)
  
  genes_RF <- multiclassPairs::sort_genes_RF(data_object = train_data_object,
                                             rank_data = TRUE,
                                             platform_wise = TRUE,
                                             num.trees = 500, # increase
                                             seed = model_seed,
                                             verbose = TRUE)
  
  rules_RF <- multiclassPairs::sort_rules_RF(data_object = train_data_object, 
                                             sorted_genes_RF = genes_RF,
                                             genes_altogether = 200, # more genes,
                                             genes_one_vs_rest = 200, # more rules
                                             platform_wise = TRUE, # pick rules that perform well across platforms
                                             num.trees = 500,# more rules, more tress are recommended 
                                             seed = model_seed,
                                             verbose = TRUE)
  
  RF_classifier <- multiclassPairs::train_RF(data_object = train_data_object,
                                             sorted_rules_RF = rules_RF,
                                             gene_repetition = 1, # * these parameters can be optimized
                                             rules_altogether = 10, # * using optimize_RF()
                                             rules_one_vs_rest = 10, # *
                                             run_boruta = TRUE, # * keep
                                             plot_boruta = FALSE, # * keep
                                             probability = TRUE, # * keep
                                             num.trees = 300, # *
                                             boruta_args = list(),
                                             verbose = TRUE)
  
  results <- multiclassPairs::predict_RF(classifier = RF_classifier, 
                                         Data = test_data_object)
  
  # get the prediction labels
  test_pred <- results$predictions
  
  # if the classifier trained using probability = FALSE
  if (is.factor(test_pred)) {
    x <- as.character(test_pred)
  }
  
  # if the classifier trained using probability = TRUE
  if (is.matrix(test_pred)) {
    x <- colnames(test_pred)[max.col(test_pred)]
  }
  
  # training accuracy
  caret::confusionMatrix(data = factor(x),
                         reference = factor(test_data_object$data$Labels, 
                                            levels = mb_subgroups),
                         mode = "everything")
  
}

#run_mm2s()
#run_lasso()