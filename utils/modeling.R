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
  if (n_repeats < 1 | round(n_repeats) != n_repeats) {
    
    stop("n_repeats in run_models() should be a positive integer.")
    
  }
  
  # n_cores should not exceed parallel::detectCores() - 1
  n_cores <- min(n_cores, parallel::detectCores() - 1)
  
  # n_rules_min is required for ktsp and is only used in ktsp
  if ("ktsp" %in% model_types) {
    
    if (is.na(n_rules_min)) {
      stop("n_rules_min in run_models() cannot be NA with ktsp model type.")
    }
    
    # ktsp should be a positive integer  
    if (n_rules_min < 1 | round(n_rules_min) != n_rules_min) {
      
      stop("n_rules_min in run_models() should be a positive integer.")
      
    }
    
    # n_rules_max should be a positive integer and >= n_rules_min
    # if n_rules_max is not given, set to n_rules_min
    if (is.na(n_rules_max)) {
      
      n_rules_max <- n_rules_min
      
    } else if (n_rules_max < 1 | round(n_rules_max) != n_rules_max) {
      
      stop("n_rules_max in run_models() should be NA or a positive integer.")
      
    }
    
  } else {
    
    n_rules_min <- NA
    n_rules_max <- NA
    
  }
  
  if ("mm2s" %in% model_types) {
    
    suppressMessages(library(MM2S))
  
    # set up gene name conversions
    
    AnnotationHub::setAnnotationHubOption("ASK", FALSE) # download without asking
    ah <- AnnotationHub::AnnotationHub()
    AnnotationHub::snapshotDate(ah) <- "2022-10-26" # reproducibility
    hs_orgdb <- AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]
    map_ENSEMBL_ENTREZID_dedup_df <- AnnotationDbi::select(x = hs_orgdb,
                                                           keys = AnnotationDbi::keys(hs_orgdb, "ENSEMBL"),
                                                           columns = "ENTREZID",
                                                           keytype = "ENSEMBL") %>%
      dplyr::mutate(dup_ENSEMBL = duplicated(ENSEMBL),
                    dup_ENTREZID = duplicated(ENTREZID)) %>%
      dplyr::filter(!dup_ENSEMBL, !dup_ENTREZID) %>%
      dplyr::select(ENSEMBL, ENTREZID)
    
  }
  
  set.seed(initial_seed)
  
  train_test_seeds <- sample(1:n_repeats, size = n_repeats)
  modeling_seeds <- sample(1:n_repeats, size = n_repeats)
  official_model_n <- sample(1:n_repeats, size = 1)
  
  # parallel backend
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,
                          c("get_train_test_samples",
                            "run_one_model",
                            "run_ktsp",
                            "run_rf",
                            "run_mm2s",
                            "run_lasso"))
  
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
                                   train_test_seed,
                                   proportion_of_studies_train = 0.5) {
  
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
                    model_seed)
    
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
  
  test_results <- multiclassPairs::predict_one_vs_rest_TSP(classifier = classifier,
                                                           Data = test_data_object,
                                                           tolerate_missed_genes = TRUE,
                                                           weighted_votes = TRUE,
                                                           classes = mb_subgroups,
                                                           verbose = TRUE)
  
  test_cm <- caret::confusionMatrix(data = factor(test_results$max_score, 
                                                  levels = mb_subgroups),
                                    reference = factor(metadata_df_test$subgroup,
                                                       levels = mb_subgroups),
                                    mode = "everything")
  
  list(classifier = classifier,
       test_results = test_results,
       test_cm = test_cm)
  
}

run_rf <- function(genex_df_train,
                   genex_df_test,
                   metadata_df_train,
                   metadata_df_test,
                   model_seed) {

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
  
  genes <- multiclassPairs::sort_genes_RF(data_object = train_data_object,
                                          rank_data = TRUE,
                                          platform_wise = TRUE,
                                          num.trees = 1000,
                                          seed = model_seed,
                                          verbose = TRUE)
  
  rules <- multiclassPairs::sort_rules_RF(data_object = train_data_object, 
                                          sorted_genes_RF = genes,
                                          genes_altogether = 1000,
                                          genes_one_vs_rest = 1000,
                                          platform_wise = TRUE,
                                          num.trees = 1000,
                                          seed = model_seed,
                                          verbose = TRUE)
  
  classifier <- multiclassPairs::train_RF(data_object = train_data_object,
                                          sorted_rules_RF = rules,
                                          gene_repetition = 1,
                                          rules_altogether = 1000,
                                          rules_one_vs_rest = 1000,
                                          run_boruta = TRUE,
                                          plot_boruta = FALSE,
                                          probability = TRUE,
                                          num.trees = 1000,
                                          verbose = TRUE)
  
  test_results <- multiclassPairs::predict_RF(classifier = classifier, 
                                              Data = test_data_object)
  
  # get the prediction matrix
  test_pred <- test_results$predictions
  
  # pick the column with maximum probability
  test_prediction_labels <- colnames(test_pred)[max.col(test_pred)]
  
  # test accuracy
  test_cm <- caret::confusionMatrix(data = factor(test_prediction_labels,
                                                  levels = mb_subgroups),
                                    reference = factor(metadata_df_test$subgroup, 
                                                       levels = mb_subgroups),
                                    mode = "everything")
  
  list(classifier = classifier,
       test_results = test_results,
       test_cm = test_cm)
  
}

run_mm2s <- function(genex_df_test,
                     metadata_df_test,
                     model_seed) {
  
  mb_subgroups <- c("G3", "G4", "NORMAL", "SHH", "WNT")
  
  genex_df_test_ENTREZID <- genex_df_test %>%
    tibble::rownames_to_column(var = "ENSEMBL") %>%
    dplyr::left_join(map_ENSEMBL_ENTREZID_dedup_df,
                     by = "ENSEMBL") %>%
    dplyr::filter(!duplicated(ENSEMBL),
                  !duplicated(ENTREZID),
                  !is.na(ENTREZID))
    dplyr::select(-ENSEMBL) %>%
    tibble::column_to_rownames(var = "ENTREZID")
  
  mm2s_predictions <- MM2S::MM2S.human(InputMatrix = genex_df_test,
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
  
  list(test_results = test_results,
       test_cm = test_cm)
  
}

run_lasso <- function(genex_df_train,
                      genex_df_test,
                      metadata_df_train,
                      metadata_df_test,
                      model_seed) {

  mb_subgroups <- c("G3", "G4", "SHH", "WNT")
  
  set.seed(model.seed)
  
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