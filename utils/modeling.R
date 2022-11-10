suppressMessages(library(tidyverse))

create_and_run_models <- function(genex_df,
                                  metadata_df,
                                  model_types = list("ktsp", "rf", "mm2s", "lasso"),
                                  initial_seed = 44,
                                  n_repeats = 1,
                                  n_cores = parallel::detectCores() - 1,
                                  n_rules_min = 5,
                                  n_rules_max = n_rules_min) {
  
  set.seed(initial_seed)
  
  train_test_seeds <- sample(1:n_repeats, size = n_repeats)
  modeling_seeds <- sample(1:n_repeats, size = n_repeats)
  official_model_n <- sample(1:n_repeats, size = 1)
  
  # parallel backend
  cl <- parallel::makeCluster(n_cores)
  registerDoParallel(cl)
  
  model_list <- foreach(n = 1:n_repeats) %dopar% {
    
    
    
    test_train_samples_list <- get_test_train_samples(metadata_df,
                                                      test_train_seed = test_train_seeds[n])
    
    genex_df_test <- genex_df %>% select(test_train_samples_list$test)
    genex_df_train <- genex_df %>% select(test_train_samples_list$train)
    metadata_df_test <- metadata_df %>% select(test_train_samples_list$test)
    metadata_df_train <- metadata_df %>% select(test_train_samples_list$train)
    
    purrr::map(model_types,
               function(x) run_model(type = x,
                                     test_df = genex_df_test,
                                     train_df = genex_df_train,
                                     test_labels = metadata_df_test$subgroup,
                                     train_labels = metadata_df_train$subgroup,
                                     model_seed = model_seeds[n],
                                     n_rules_min,
                                     n_rules_max))
    
  }
  
  # stop parallel backend
  stopCluster(cl)
  
  # set official model
  model_list[[n]]$official_model <- TRUE
  
  return(model_list)
  
}


get_test_train_samples <- function(metadata_df,
                                   test_train_seed) {
  
}

run_model() <- function(type,
                        test_df,
                        train_df,
                        test_labels,
                        train_labels,
                        model_seed,
                        n_rules_min,
                        n_rules_max) {
  
  if (type == "ktsp") {
    
    model <- run_ktsp(test_df, train_df, test_labels, train_labels, model_seed, n_rules_min, n_rules_max)
    
  } else if (type == "rf") {
    
    model <- run_rf(test_df, train_df, test_labels, train_labels, model_seed, n_rules_min, n_rules_max)
    
  } else if (type == "mm2s") {
    
    model <- run_mm2s(train_df, train_labels, model_seed)
    
  } else if (type == "lasso") {
    
    model <- run_lasso(test_df, train_df, test_labels, train_labels, model_seed)
    
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

  set.seed(model_seed)
  
  train_data_object <- multiclassPairs::ReadData(Data = genex_df_train,
                                                 Labels = metadata_df_train$subgroup,
                                                 Platform = metadata_df_train$platform,
                                                 verbose = FALSE)
  
  test_data_object <- multiclassPairs::ReadData(Data = genex_df_test,
                                                Labels = metadata_df_test$subgroup,
                                                Platform = metadata_df_test$platform,
                                                verbose = FALSE)
  
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
                                                       verbose = FALSE)
  
  train_results <- multiclassPairs::predict_one_vs_rest_TSP(classifier = classifier,
                                                            Data = train_data_object,
                                                            tolerate_missed_genes = TRUE,
                                                            weighted_votes = TRUE,
                                                            classes = c("G3", "G4", "SHH", "WNT"),
                                                            verbose = TRUE)
  
  test_results <- multiclassPairs::predict_one_vs_rest_TSP(classifier = classifier,
                                                           Data = test_data_object,
                                                           tolerate_missed_genes = TRUE,
                                                           weighted_votes = TRUE,
                                                           classes = c("G3", "G4", "SHH", "WNT"),
                                                           verbose = TRUE)
  
  train_cm <- caret::confusionMatrix(data = factor(train_results$max_score, 
                                                   levels = unique(train_data_object$data$Labels)),
                                     reference = factor(train_data_object$data$Labels, 
                                                        levels = unique(train_data_object$data$Labels)),
                                     mode = "everything")
  
  test_cm <- caret::confusionMatrix(data = factor(test_results$max_score, 
                                                  levels = unique(test_data_object$data$Labels)),
                                    reference = factor(test_data_object$data$Labels, 
                                                       levels = unique(test_data_object$data$Labels)),
                                    mode = "everything")
  
  list(train_data_object = train_data_object,
       test_data_object = test_data_object,
       filtered_genes = filtered_genes,
       classifier = classifier,
       train_results = train_results,
       test_results = test_results,
       train_cm = train_cm,
       test_cm = test_cm)
  
}

run_rf()
run_mm2s()
run_lasso()