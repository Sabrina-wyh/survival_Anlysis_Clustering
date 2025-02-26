####################### preparing for different data structure #################
set.seed(42)
prepare_for_survival_analysis_n_visits_backward <- function(df, n) {
  process_group <- function(group) {
    ad_indices <- which(group$event == 1)
    start_indices <- which(group$event == 0)
    
    if (length(ad_indices) > 0) {
      if (nrow(group) <= n) {
        return(group)
      } else {
        first_ad_idx <- ad_indices[1]
        first_ad_stop <- group$stop[first_ad_idx]
        retained <- tail(group, n)
        retained[nrow(retained), "event"] <- 1
        retained[nrow(retained), "stop"] <- first_ad_stop
        return(retained)
      }
    } else {
      retained <- tail(group, n)
      retained[nrow(retained), "stop"] <- group$stop[nrow(group)]
      return(retained)
    }
  }
  
  df_processed <- do.call(rbind, lapply(split(df, df$ID), process_group))
  return(df_processed)
}

prepare_for_survival_analysis_n_visits <- function(df, n) {
  process_group <- function(group) {
    ad_indices <- which(group$event == 1)
    start_indices <- which(group$event == 0)[1]
    
    if (length(ad_indices) > 0) {
      if (nrow(group) <= n) {
        return(group)
      } else {
        first_ad_idx <- ad_indices[1]
        first_ad_stop <- group$stop[first_ad_idx]
        retained <- group[1:min(start_indices + n - 1, nrow(group)), ]
        retained[nrow(retained), "event"] <- 1
        retained[nrow(retained), "stop"] <- first_ad_stop
        return(retained)
      }
    } else {
      retained <- head(group, n)
      retained[nrow(retained), "stop"] <- group$stop[nrow(group)]
      return(retained)
    }
  }
  
  df_processed <- do.call(rbind, lapply(split(df, df$ID), process_group))
  return(df_processed)
}

aggregate_for_standard_cox <- function(df_processed, static_features, dynamic_features, m) {
  aggregated_records <- list()
  grouped <- split(df_processed, df_processed$ID)
  
  for (ID in names(grouped)) {
    group <- grouped[[ID]]
    group <- group[order(group$start), ]
    n <- nrow(group)
    final_record <- tail(group, 1)
    final_record_stop <- final_record$stop
    final_record_event <- final_record$event
    
    if (m != 'all') {
      subset <- tail(group, m)
    } else {
      subset <- group
    }
    
    summary_stats <- list()
    for (feature in dynamic_features) {
      summary_stats[[paste(feature, 'min', sep = '_')]] <- min(subset[[feature]])
      summary_stats[[paste(feature, 'max', sep = '_')]] <- max(subset[[feature]])
      summary_stats[[paste(feature, 'median', sep = '_')]] <- median(subset[[feature]])
      summary_stats[[paste(feature, 'std', sep = '_')]] <- sd(subset[[feature]])
      
      last_record <- tail(subset, 1)
      first_record <- head(subset, 1)
      duration <- last_record$stop - first_record$start
      if (duration != 0) {
        slope <- (last_record[[feature]] - first_record[[feature]]) / duration
      } else {
        slope <- NA  # Handle division by zero if duration is zero
      }
      summary_stats[[paste(feature, 'slope', sep = '_')]] <- slope
    }
    
    static_data <- list()
    for (feature in static_features) {
      static_data[[feature]] <- subset[[1, feature]]
    }
    

    time <- final_record_stop - last_record$start
    aggregated_record <- list(
      ID = ID,
      time = time,
      event = final_record_event
    )
    aggregated_record <- c(aggregated_record, static_data, summary_stats)
    aggregated_records[[length(aggregated_records) + 1]] <- aggregated_record
  }
  df_aggregated <- do.call(rbind, lapply(aggregated_records, as.data.frame))
  return(df_aggregated)
}

missranger_impute <- function(df, columns_ignores) {
  if (!is.data.frame(df)) {
    stop("The input must be a dataframe.")
  }
  if (!all(columns_ignores %in% colnames(df))) {
    stop("Some ID columns specified are not present in the dataframe.")
  }
  
  
  ignore_cols <- df[, columns_ignores, drop = FALSE]
  data_to_impute <- df[, !colnames(df) %in% columns_ignores, drop = FALSE]
  set.seed(42)
  imputed_data <- missRanger(data_to_impute, pmm.k = 10)
  result <- cbind(ignore_cols, imputed_data)
  return(result)
}


######################### For Model Construction/Evaluation ####################


remove_high_corr_scales <- function(df_train, df_test, cols_consider=NULL, threshold = 0.9) { 
  ### first do the normalisation
  df_train[, cols_consider] <- as.data.frame(scale(df_train[, cols_consider])) 
  df_test[, cols_consider] <- as.data.frame(scale(df_test[, cols_consider])) 
  
  ### remove high correlated features
  correlation_matrix <- cor(df_train[, cols_consider], use = "pairwise.complete.obs")   
  high_corr <- findCorrelation(correlation_matrix, cutoff = threshold)   
  
  df_train <- df_train[, !(names(df_train) %in% high_corr)]  
  df_test <- df_test[, !(names(df_test) %in% high_corr)]   
  
  return(list(
    train = df_train,
    test = df_test
  ))
}

print_model_statistics <- function(all_c_index_list, auc_summary) {
  # Compute overall statistics for C-index
  all_c_indexes <- unlist(all_c_index_list)
  overall_c_median <- median(all_c_indexes, na.rm = TRUE)
  c_Q1 <- quantile(all_c_indexes, 0.25, na.rm = TRUE)
  c_Q3 <- quantile(all_c_indexes, 0.75, na.rm = TRUE)
  
  cat("Overall C-index Median:", overall_c_median, "\n")
  cat("C-index 25th Percentile (Q1):", c_Q1, "\n")
  cat("C-index 75th Percentile (Q3):", c_Q3, "\n\n")
  
  # Calculate overall statistics for AUC (across all iterations, folds, and time points)
  overall_auc_median <- median(auc_summary$AUC, na.rm = TRUE)
  overall_auc_q1 <- quantile(auc_summary$AUC, 0.25, na.rm = TRUE)
  overall_auc_q3 <- quantile(auc_summary$AUC, 0.75, na.rm = TRUE)
  cat("Overall AUC Median:", overall_auc_median, "\n")
  cat("AUC 25th Percentile (Q1):", overall_auc_q1, "\n")
  cat("AUC 75th Percentile (Q3):", overall_auc_q3, "\n\n")
}



do_internal_cv_time_varying_cox <- function(data, iterations = 10, test_ratio = 0.1, folds = 5, feature_selection = FALSE, shorten_visits = FALSE, n_visits = NULL) {
  all_predictors <- c()
  all_c_index_list <- list()
  all_auc_results <- list()
  num_selected_features <- c()
  
  # Set up parallel backend
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("remove_high_corr_scales", "prepare_for_survival_analysis_n_visits", 
                                "prepare_for_survival_analysis_n_visits_backward", 
                                "print_model_statistics", "timeROC"), envir = environment())
  on.exit(stopCluster(cl), add=TRUE)
  
  results <- foreach(iter = 1:iterations, .packages = c("survival", "MASS", "timeROC", "dplyr", "caret")) %dopar% {
    set.seed(iter)
    best_c_index <- 0
    best_features <- NULL
    unique_ids <- unique(data$ID)
    test_ids <- sample(unique_ids, size = floor(test_ratio * length(unique_ids)))
    train_ids <- setdiff(unique_ids, test_ids)
    train_data <- data[data$ID %in% train_ids, ]
    test_data <- data[data$ID %in% test_ids, ]
    fold_assignments <- sample(rep(1:folds, length.out = length(train_ids)))
    
    if (shorten_visits && !is.null(n_visits)) {
      test_data <- prepare_for_survival_analysis_n_visits_backward(test_data, n_visits)
    }
    
    for (fold in 1:folds) {
      fold_train_ids <- train_ids[fold_assignments != fold]
      fold_validation_ids <- train_ids[fold_assignments == fold]
      fold_train_data <- train_data[train_data$ID %in% fold_train_ids, ]
      fold_validation_data <- train_data[train_data$ID %in% fold_validation_ids, ]
      
      if (shorten_visits && !is.null(n_visits)) {
        fold_validation_data <- prepare_for_survival_analysis_n_visits_backward(fold_validation_data, n_visits)
      }
      
      predictors <- setdiff(names(fold_train_data), c("start", "stop", "event", "ID"))
      data_preprocess <- remove_high_corr_scales(fold_train_data, fold_validation_data, cols_consider = predictors, threshold = 0.9)
      fold_train_data <- data_preprocess$train
      fold_validation_data <- data_preprocess$test
      
      predictors <- setdiff(names(fold_train_data), c("start", "stop", "event", "ID"))
      
      if (feature_selection) {
        formula_no_frailty <- as.formula(paste("Surv(start, stop, event) ~", paste(predictors, collapse = " + ")))
        cox_model_no_frailty <- coxph(formula_no_frailty, data = fold_train_data, ties = "efron", control = coxph.control(iter.max = 100, eps = 1e-05))
        stepwise_model <- stepAIC(cox_model_no_frailty, direction = "both", trace = FALSE)
        predictors <- all.vars(formula(stepwise_model))[-1]
        predictors <- setdiff(predictors, c("start", "stop", "event", "ID"))
      }
      
      formula <- as.formula(paste("Surv(start, stop, event) ~", paste(predictors, collapse = " + "), "+ frailty(ID)"))
      cox_model <- coxph(formula, data = fold_train_data, id = fold_train_data$ID, ties = "efron", control = coxph.control(iter.max = 100, eps = 1e-05))
      
      risk_scores <- predict(cox_model, newdata = fold_validation_data, type = "risk")
      c_index <- concordance(cox_model, newdata = fold_validation_data)$concordance
      
      if (c_index > best_c_index) {
        best_c_index <- c_index
        best_features <- predictors
      }
    }
    
    data_preprocess_final <- remove_high_corr_scales(train_data, test_data, cols_consider = predictors, threshold = 0.9)
    train_data <- data_preprocess_final$train
    test_data <- data_preprocess_final$test
    
    final_train_formula <- as.formula(paste("Surv(start, stop, event) ~", paste(best_features, collapse = " + "), "+ frailty(ID)"))
    final_model <- coxph(final_train_formula, data = train_data, id = train_data$ID, ties = "efron", control = coxph.control(iter.max = 100, eps = 1e-05))
    
    final_risk_scores <- predict(final_model, newdata = test_data, type = "risk")
    final_c_index <- concordance(final_model, newdata = test_data)$concordance
    
    times <- na.omit(as.numeric(unique(test_data$stop)))
    times <- times[times > 0]
    
    auc_result <- timeROC(
      T = test_data$stop,
      delta = test_data$event,
      marker = final_risk_scores,
      cause = 1,
      times = times
    )
    
    auc_df <- data.frame(
      Iteration = iter,
      Time = auc_result$times,
      AUC = auc_result$AUC
    )
    
    list(
      C_Index = final_c_index,
      Best_Features = best_features,
      AUC_Results = auc_df
    )
  }
  
  
  all_c_index_list <- lapply(results, function(res) res$C_Index)
  all_auc_results <- do.call(rbind, lapply(results, function(res) res$AUC_Results))
  all_predictors <- unlist(lapply(results, function(res) res$Best_Features))
  
  auc_stats <- all_auc_results %>%
    group_by(Time) %>%
    summarise(
      Median_AUC = median(AUC, na.rm = TRUE),
      Q1_AUC = quantile(AUC, 0.25, na.rm = TRUE),
      Q3_AUC = quantile(AUC, 0.75, na.rm = TRUE)
    )
  
  print_model_statistics(all_c_index_list, all_auc_results)
  
  significant_features_freq <- data.frame(
    Feature = names(table(all_predictors)),
    Frequency = as.integer(table(all_predictors))
  ) %>% arrange(desc(Frequency))
  
  return(list(
    C_Index_List = unlist(all_c_index_list),
    AUC_list = unlist(all_auc_results$AUC),
    AUC_Stats = auc_stats,
    Significant_Features_Freq = significant_features_freq
  ))
}


do_external_time_varying_cox <- function(train_data, test_data, shorten_visits = FALSE, n_visits=NULL) {
  if (shorten_visits) {
    if (is.null(n_visits)) {
      stop("If shorten_visits is TRUE, n_visits must be specified.")
    }
    test_data <- prepare_for_survival_analysis_n_visits_backward(test_data, n_visits)
  }
  
  predictors <- setdiff(names(train_data), c("start", "stop", "event", "ID"))
  data_preprocess <- remove_high_corr_scales(train_data, test_data, cols_consider = predictors, threshold = 0.9)
  train_data <- data_preprocess$train
  test_data <- data_preprocess$test
  
  formula <- as.formula(paste("Surv(start, stop, event) ~", paste(predictors, collapse = " + "), "+ frailty(ID)"))
  cox_model <- coxph(formula, data = train_data, id = train_data$ID, ties = "efron",control = coxph.control(iter.max = 50, eps = 1e-04))
  print(summary(cox_model))
  
  
  # risk_scores <- predict(cox_model, newdata = test_data, type = "risk")
  #######
  risk_scores_training <- predict(cox_model, newdata = train_data, type = "risk")
  risk_scores_testing <- predict(cox_model, newdata = test_data, type = "risk")
  
  
  risk_score_final_train <- data.frame(
    ID = train_data$ID,
    time = train_data$start,
    risk_score = risk_scores_training
  )
  final_risk_score_train <- risk_score_final_train %>%
    arrange(ID, time) %>% 
    group_by(ID) %>%
    filter(row_number() == n()) %>%
    summarise(last_risk_score = risk_score)
  
  risk_score_final_test <- data.frame(
    ID = test_data$ID,
    time = test_data$start,
    risk_score = risk_scores_testing
  )
  final_risk_score_test <- risk_score_final_test %>%
    arrange(ID, time) %>% 
    group_by(ID) %>%
    filter(row_number() == n()) %>%
    summarise(last_risk_score = risk_score)
  
  ######
  c_index <- concordance(cox_model, newdata = test_data)$concordance
  times <- unique(test_data$stop)
  times <- as.numeric(times)
  times <- times[!is.na(times) & times > 0]  # Remove NA and non-positive times
  
  auc_result <- timeROC(
    T = test_data$stop,
    delta = test_data$event,
    marker = risk_scores_testing,
    cause = 1,
    times = times
  )
  
  auc_summary <- data.frame(
    Time = auc_result$times,
    AUC = auc_result$AUC
  )
  
  overall_auc_median <- median(auc_summary$AUC, na.rm = TRUE)
  overall_auc_q1 <- quantile(auc_summary$AUC, 0.25, na.rm = TRUE)
  overall_auc_q3 <- quantile(auc_summary$AUC, 0.75, na.rm = TRUE)
  cat("C-index:", c_index, "\n")
  cat("Overall AUC Median:", overall_auc_median, "\n")
  cat("AUC 25th Percentile (Q1):", overall_auc_q1, "\n")
  cat("AUC 75th Percentile (Q3):", overall_auc_q3, "\n")
  
  return(list(
    C_Index = c_index,
    AUC_list = unlist(auc_summary$AUC),
    AUC_Stats = auc_summary,
    model_summary = summary(cox_model),
    risk_train = final_risk_score_train,
    risk_test = final_risk_score_test
  ))
  
}


do_internal_cv_normal_cox <- function(data, static_ft = NULL, dynamic_ft = NULL, test_ratio = 0.1, folds = 5, iterations = 10,
                                      feature_selection = FALSE, shorten_visits = FALSE, n_visits = NULL) {
  
  
  all_predictors <- c()
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add=TRUE)
  clusterExport(cl, varlist = c("remove_high_corr_scales", "prepare_for_survival_analysis_n_visits", 
                                "prepare_for_survival_analysis_n_visits_backward", 
                                "aggregate_for_standard_cox",
                                "print_model_statistics", "timeROC"), envir = environment())
  
  results <- foreach(iter = 1:iterations, .packages = c("survival", "MASS", "timeROC", "dplyr", "caret")) %dopar% {
    set.seed(iter)
    unique_ids <- unique(data$ID)
    test_ids <- sample(unique_ids, size = floor(test_ratio * length(unique_ids)))
    train_ids <- setdiff(unique_ids, test_ids)
    train_data <- data[data$ID %in% train_ids, ]
    test_data <- data[data$ID %in% test_ids, ]
    
    fold_assignments <- sample(rep(1:folds, length.out = length(train_ids)))
    names(fold_assignments) <- train_ids
    best_c_index <- 0
    best_features <- NULL
    
    for (fold in 1:folds) {
      train_fold_ids <- train_ids[fold_assignments != fold]
      val_fold_ids <- train_ids[fold_assignments == fold]
      fold_train_data <- train_data[train_data$ID %in% train_fold_ids, ]
      fold_validation_data <- train_data[train_data$ID %in% val_fold_ids, ]
      
      if (!is.null(static_ft) || !is.null(dynamic_ft)) {
        fold_train_data <- aggregate_for_standard_cox(fold_train_data, static_ft, dynamic_ft, m = "all")
        fold_train_data_cols_valid <- names(fold_train_data[, colSums(is.na(fold_train_data)) == 0])
        predictors <- setdiff(fold_train_data_cols_valid, c("time", "event", "ID"))
        
        if (shorten_visits) {
          fold_validation_data <- aggregate_for_standard_cox(fold_validation_data, static_ft, dynamic_ft, m = n_visits)
        } else {
          fold_validation_data <- aggregate_for_standard_cox(fold_validation_data, static_ft, dynamic_ft, m = "all")
        }
        fold_validation_data_cols_valid <- names(fold_validation_data[, colSums(is.na(fold_validation_data)) == 0])
        predictors <- intersect(predictors, setdiff(fold_validation_data_cols_valid, c("time", "event", "ID")))
      } else {
        predictors <- setdiff(names(fold_train_data), c("time", "event", "ID"))
      }
      
      data_preprocess <- remove_high_corr_scales(fold_train_data, fold_validation_data, cols_consider = predictors, threshold = 0.9)
      fold_train_data <- data_preprocess$train
      fold_validation_data <- data_preprocess$test
      
      formula <- as.formula(paste("Surv(time, event) ~", paste(predictors, collapse = " + ")))
      cox_model <- coxph(formula, data = fold_train_data)
      
      if (feature_selection) {
        stepwise_model <- suppressWarnings(stepAIC(cox_model, direction = "both", trace = FALSE, steps = 5000, k = 10))
        predictors <- setdiff(all.vars(formula(stepwise_model))[-1], c("time", "event", "ID"))
      }
      
      risk_scores <- predict(cox_model, newdata = fold_validation_data, type = "risk")
      c_index <- concordance(cox_model, newdata = fold_validation_data)$concordance
      
      if (c_index > best_c_index) {
        best_c_index <- c_index
        best_features <- predictors
      }
    }
    
    if (!is.null(static_ft) || !is.null(dynamic_ft)) {
      train_data <- aggregate_for_standard_cox(train_data, static_ft, dynamic_ft, m = "all")
      if (shorten_visits) {
        test_data <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = n_visits)
      } else {
        test_data <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = "all")
      }
    }
    
    final_formula <- as.formula(paste("Surv(time, event) ~", paste(best_features, collapse = " + ")))
    final_model <- coxph(final_formula, data = train_data)
    
    final_risk_scores <- predict(final_model, newdata = test_data, type = "risk")
    final_c_index <- concordance(final_model, newdata = test_data)$concordance
    
    times <- unique(test_data$time)
    times <- as.numeric(times[!is.na(times) & times > 0])
    
    auc_result <- timeROC(
      T = test_data$time,
      delta = test_data$event,
      marker = final_risk_scores,
      cause = 1,
      times = times
    )
    
    auc_df <- data.frame(
      Iteration = iter,
      Time = auc_result$times,
      AUC = auc_result$AUC
    )
    
    list(
      C_Index = final_c_index,
      AUC_DF = auc_df,
      Predictors = best_features
    )
  }
  
  
  all_c_index_list <- sapply(results, function(x) x$C_Index)
  all_auc_results_df <- do.call(rbind, lapply(results, function(x) x$AUC_DF))
  all_predictors <- unlist(lapply(results, function(x) x$Predictors))
  
  auc_stats <- all_auc_results_df %>%
    group_by(Time) %>%
    summarise(
      Median_AUC = median(AUC, na.rm = TRUE),
      Q1_AUC = quantile(AUC, 0.25, na.rm = TRUE),
      Q3_AUC = quantile(AUC, 0.75, na.rm = TRUE)
    )
  
  print_model_statistics(all_c_index_list, all_auc_results_df)
  
  significant_features_freq <- data.frame(
    Feature = names(table(all_predictors)),
    Frequency = as.integer(table(all_predictors))
  ) %>% arrange(desc(Frequency))
  
  return(list(
    C_Index_List = unlist(all_c_index_list),
    AUC_list = unlist(all_auc_results_df$AUC),
    AUC_Stats = auc_stats,
    Significant_Features_Freq = significant_features_freq
  ))
}


do_external_normal_cox <- function(train_data, test_data, static_ft=NULL, dynamic_ft=NULL, shorten_visits = FALSE, n_visits=NULL) {
  
  
  if (!is.null(static_ft) || !is.null(dynamic_ft)) {
    train_data = aggregate_for_standard_cox(train_data, static_ft, dynamic_ft, m = "all")
    train_data <- train_data[, colSums(is.na(train_data)) == 0]
    predictors <- names(train_data)[-which(names(train_data) %in% c("time", "event", "ID"))]
    if (shorten_visits) {
      if (is.null(n_visits)) {
        stop("If shorten_visits is TRUE, n_visits must be specified.")
      }
      test_data <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = n_visits)
      test_data <- test_data[, colSums(is.na(test_data)) == 0]
      test_var = names(test_data)[-which(names(test_data) %in% c("time", "event", "ID"))]
      predictors <- intersect(predictors, test_var)
    }
    else{
      test_data <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = 'all')
      test_var = names(test_data)[-which(names(test_data) %in% c("time", "event", "ID"))]
      predictors <- intersect(predictors, test_var)
    }
  }
  else{
    predictors <- setdiff(names(train_data), c("time", "event", "ID"))
  }
  
  data_preprocess <- remove_high_corr_scales(train_data, test_data, cols_consider = predictors, threshold = 0.9)
  train_data <- data_preprocess$train
  test_data <- data_preprocess$test
  
  formula <- as.formula(paste("Surv(time, event) ~", paste(predictors, collapse = " + ")))
  cox_model <- suppressWarnings(coxph(formula, data = train_data, control = coxph.control(iter.max = 50, eps = 1e-04)))
  print(summary(cox_model))
  
  
  risk_scores <- predict(cox_model, newdata = test_data, type = "risk")
  c_index <- concordance(cox_model, newdata = test_data)$concordance
  times <- unique(test_data$time)
  times <- as.numeric(times)
  times <- times[!is.na(times) & times > 0]  # Remove NA and non-positive times
  
  auc_result <- timeROC(
    T = test_data$time,
    delta = test_data$event,
    marker = risk_scores,
    cause = 1,
    times = times
  )
  
  auc_summary <- data.frame(
    Time = auc_result$times,
    AUC = auc_result$AUC
  )
  
  overall_auc_median <- median(auc_summary$AUC, na.rm = TRUE)
  overall_auc_q1 <- quantile(auc_summary$AUC, 0.25, na.rm = TRUE)
  overall_auc_q3 <- quantile(auc_summary$AUC, 0.75, na.rm = TRUE)
  cat("C-index:", c_index, "\n")
  cat("Overall AUC Median:", overall_auc_median, "\n")
  cat("AUC 25th Percentile (Q1):", overall_auc_q1, "\n")
  cat("AUC 75th Percentile (Q3):", overall_auc_q3, "\n")
  
  return(list(
    C_Index = c_index,
    AUC_list = unlist(auc_summary$AUC),
    AUC_Stats = auc_summary,
    model_summary = summary(cox_model)
  ))
}

################################## For Clustering ##############################
analyze_risk_groups <- function(df, risk_group_df) {
  for (group in unique(risk_group_df$risk_group)) {
    cat(sprintf("\n--- Risk Group: %s ---\n", group))
    
    group_ids <- risk_group_df$ID[risk_group_df$risk_group == group]
    group_data <- df %>%
      filter(ID %in% group_ids) %>%
      arrange(ID, start)
    
    baseline <- group_data %>% group_by(ID) %>% slice(1) %>% ungroup()
    APOE4_positive_count <- sum(baseline$APOE4 == 1, na.rm = TRUE)
    APOE4_positive_percentage <- (APOE4_positive_count / nrow(baseline)) * 100
    
    
    
    gender_female_count <- sum(baseline$sex == 2, na.rm = TRUE)
    gender_female_percentage <- (gender_female_count / nrow(baseline)) * 100
    follow_up_intervals <- group_data %>% group_by(ID) %>%
      summarize(
        min_visit = min(start, na.rm = TRUE),
        max_visit = max(start, na.rm = TRUE)
      ) %>%
      ungroup()
    follow_up_range <- paste0(
      min(follow_up_intervals$max_visit, na.rm = TRUE),
      " - ",
      max(follow_up_intervals$max_visit, na.rm = TRUE),
      " months"
    )
    participants_with_event <- sum(group_data$event == 1, na.rm = TRUE)
    total_participants <- n_distinct(group_data$ID)
    AD_conversion_percentage <- (participants_with_event / total_participants) * 100
    if (participants_with_event > 0) {
      converters <- group_data %>% filter(event == 1)
      time_to_conversion <- converters %>%
        group_by(ID) %>%
        summarize(min_visit = min(start, na.rm = TRUE)) %>%
        ungroup() %>%
        pull(min_visit)
      time_to_conversion_mean <- mean(time_to_conversion, na.rm = TRUE)
      time_to_conversion_std <- sd(time_to_conversion, na.rm = TRUE)
    } else {
      time_to_conversion_mean <- NA
      time_to_conversion_std <- NA
    }
    
    # Print results
    total_participants <- length(unique(group_data$ID))
    cat(sprintf("Total participants: %d\n", total_participants))
    cat("APOE4 Positive (n, %):", APOE4_positive_count,
        "(", round(APOE4_positive_percentage, 2), "%)\n")
    cat("Gender (Female, n, %):", gender_female_count,
        "(", round(gender_female_percentage, 2), "%)\n")
    cat("Follow-up Interval:", follow_up_range, "\n")
    cat("AD Conversion (n, %):", participants_with_event,
        "(", round(AD_conversion_percentage, 2), "%)\n")
    if (!is.na(time_to_conversion_mean)) {
      cat("Time to AD Conversion (mean ± std):",
          round(time_to_conversion_mean, 2), "±",
          round(time_to_conversion_std, 2), "months\n")
    } else {
      cat("Time to AD Conversion (mean ± std): N/A\n")
    }
  }
}



##################################### plots ####################################

create_forest_plot <- function(summary_results, title = "Forest Plot for Cox Model", file_path_name=NULL) {
  coefficients <- summary_results$coefficients
  conf_intervals <- summary_results$conf.int
  
  # Exclude the frailty term if present
  if ("frailty(ID)" %in% rownames(coefficients)) {
    coefficients <- coefficients[rownames(coefficients) != "frailty(ID)", ]
    conf_intervals <- conf_intervals[rownames(conf_intervals) != "frailty(ID)", ]
  }
  
  
  forest_data <- data.frame(
    Predictor = rownames(coefficients),
    HR = conf_intervals[, "exp(coef)"],  # Hazard ratio
    Lower_CI = conf_intervals[, "lower .95"],  # Lower 95% CI
    Upper_CI = conf_intervals[, "upper .95"],  # Upper 95% CI
    # p_value = coefficients[, "Pr(>|z|)"]  # P-value
    p_value = coefficients[, "p"]  # P-value
  )
  
  forest_data <- forest_data %>%arrange(HR)
  forest_data$Predictor <- factor(forest_data$Predictor, levels = forest_data$Predictor)
  
  # print(forest_data)
  
  ggplot(forest_data, aes(y = Predictor, x = HR)) +
    geom_point(aes(color = p_value < 0.05), size = 3) +
    # geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI, color = p_value < 0.05), height = 0.2) +
    scale_color_manual(values = c("black", "red")) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = title,
      x = "Hazard Ratio (HR)",
      y = "Predictors"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    )
  ggsave(file_path_name, width = 10, height = 5, device = "pdf")
}

plot_auc_trajectory_group <- function(auc_stats_list, labels, validation_auc_stats = NULL, validation_label = "Validation", title_name = "", output_file = NULL) {
  combined_auc_stats <- do.call(rbind, lapply(1:length(auc_stats_list), function(i) {
    auc_stats_list[[i]] %>%
      filter(!is.na(Median_AUC) & !is.na(Q1_AUC) & !is.na(Q3_AUC)) %>%
      mutate(Model = labels[i])  # Add a label column for the model
  }))
  
  # If validation AUC stats are provided, format them and add to the combined data
  if (!is.null(validation_auc_stats)) {
    validation_auc_stats <- validation_auc_stats %>%
      filter(!is.na(AUC)) %>%
      rename(Median_AUC = AUC) %>%  # Rename AUC column for consistency
      mutate(
        Q1_AUC = NA,  # Validation doesn't have Q1_AUC
        Q3_AUC = NA,  # Validation doesn't have Q3_AUC
        Model = validation_label  # Add label for validation
      )
    
    combined_auc_stats <- bind_rows(combined_auc_stats, validation_auc_stats)
  }
  
  # Create the plot
  auc_plot <- ggplot(combined_auc_stats, aes(x = Time, y = Median_AUC, color = Model, group = Model)) +
    geom_line(linewidth = 1) +  # Line for Median AUC
    geom_point(data = combined_auc_stats %>% filter(!is.na(Q1_AUC)), size = 2) +  # Points only for models with IQR
    geom_errorbar(data = combined_auc_stats %>% filter(!is.na(Q1_AUC)), aes(ymin = Q1_AUC, ymax = Q3_AUC), width = 0.2) +  # Error bars for models with IQR
    geom_line(data = validation_auc_stats, linewidth = 1, linetype = "dotted") +  # Dotted line for validation data
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1) +  # Reference line at AUC = 0.5
    theme_minimal() +
    labs(
      title = paste0("Time-Varying AUC with IQR ", title_name),
      x = "Time",
      y = "AUC",
      color = "Model"  # Legend title
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
      # panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  
  # Save the plot as a PDF if an output file is provided
  if (!is.null(output_file)) {
    ggsave(output_file, plot = auc_plot, device = "pdf", width = 10, height = 5)
  }
  
  # Return the plot
  return(auc_plot)
}

plot_auc_trajectory_single <- function(auc_stats_list, labels, title_name = "", output_file = NULL) {
  # Combine AUC stats from all models into a single dataframe
  combined_auc_stats <- do.call(rbind, lapply(1:length(auc_stats_list), function(i) {
    auc_stats_list[[i]] %>%
      filter(!is.na(AUC)) %>%
      mutate(Model = labels[i])  # Add a label column for the model
  }))
  
  # Create the plot
  auc_plot <- ggplot(combined_auc_stats, aes(x = Time, y = AUC, color = Model, group = Model)) +
    geom_line(linewidth = 1) +  # Line for AUC
    geom_point(size = 2) +  # Points for AUC
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1) +  # Reference line at AUC = 0.5
    theme_minimal() +
    labs(
      title = paste0("Time-Varying AUC ", title_name),
      x = "Time",
      y = "AUC",
      color = "Model"  # Legend title
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    )
  
  # Save the plot as a PDF if an output file is provided
  if (!is.null(output_file)) {
    ggsave(output_file, plot = auc_plot, device = "pdf", width = 10, height = 5)
  }
  
  # Return the plot
  return(auc_plot)
}

plot_stacked_feature_importance <- function(dfs, labels, colors, file_path_name) {
  merged_df <- dfs[[1]] %>% rename(!!labels[1] := Frequency) %>% column_to_rownames(var = "Feature")
  for (i in 2:length(dfs)) {
    merged_df <- merge(
      merged_df,
      dfs[[i]] %>% rename(!!labels[i] := Frequency) %>% column_to_rownames(var = "Feature"),
      by = "row.names",
      all = TRUE
    )
    rownames(merged_df) <- merged_df$Row.names
    merged_df <- merged_df[,-1]
  }
  merged_df[is.na(merged_df)] <- 0
  merged_df$Total <- rowSums(merged_df)
  top_features_df <- merged_df %>% arrange(desc(Total)) %>% head(20) %>% select(-Total)
  top_features_long <- top_features_df %>% 
    rownames_to_column(var = "Feature") %>% 
    pivot_longer(-Feature, names_to = "Model", values_to = "Frequency")
  ggplot(top_features_long, aes(x = Frequency, y = fct_reorder(Feature, Frequency, .fun = sum), fill = Model)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    labs(
      title = "Top Feature Frequency for Time-Varying Cox Model",
      y = "Feature Name",  # Changed from x to y since coord_flip() switches them
      x = "Frequency (Count)",  # Changed from y to x since coord_flip() switches them
      fill = "Model"
      
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    ) +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = TRUE)) 
  
  ggsave(file_path_name, width = 10, height = 5, device = "pdf")
}

add_data <- function(final_df, model_name, feature_set_name, value_list) {
  temp_df <- data.frame(model = rep(model_name, length(value_list)),
                        feature_set = rep(feature_set_name, length(value_list)),
                        value = value_list,
                        stringsAsFactors = FALSE)
  final_df <- rbind(final_df, temp_df)
  return(final_df)
}

plot_boxplot <- function(df, y_label, save_path) {
  p <- ggplot(data = df, aes(x=model, y=value)) +
    geom_boxplot(aes(fill=feature_set), position=position_dodge(width=0.85)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis.text.x = element_blank(),  # Removes x-axis text labels
      # axis.ticks.x = element_blank(), # Removes x-axis ticks
      legend.text = element_text(size = 8),
      axis.title.x = element_blank()
    ) +
    labs(y = y_label) + # Change y-axis label dynamically
    scale_x_discrete(expand = expansion(mult = c(0.2, 0.2)))

  # Save the plot as a PDF
  ggsave(filename = save_path, plot = p, device = "pdf", width = 14, height = 4)

  return(p)
}



plot_centroid_trajectories <- function(df, risk_group, selected_features_trajectories, pdf_path) {
  common_time <- seq(min(df$start), max(df$start), length.out = 100)
  centroid_trajectories <- list()
  ci_trajectories <- list()
  
  for (feature in selected_features_trajectories) {
    centroid_trajectories[[feature]] <- lapply(unique(risk_group$risk_group), function(group) {
      group_ids <- risk_group$ID[risk_group$risk_group == group]
      group_data <- df %>%
        filter(ID %in% group_ids) %>%
        select(ID, start, all_of(feature))
      feature_matrix <- do.call(rbind, lapply(split(group_data, group_data$ID), function(df) {
        approx(x = df$start, y = df[[feature]], xout = common_time, rule = 2)$y
      }))
      
      list(mean = colMeans(feature_matrix, na.rm = TRUE),
           lower = apply(feature_matrix, 2, function(x) mean(x, na.rm = TRUE) - 1.96 * sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
           upper = apply(feature_matrix, 2, function(x) mean(x, na.rm = TRUE) + 1.96 * sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
      )
    })
  }
  
  plot_list <- list()
  for (feature in selected_features_trajectories) {
    plot_data <- data.frame(
      Time = rep(common_time, times = length(unique(risk_group$risk_group))),
      Mean = unlist(lapply(centroid_trajectories[[feature]], `[[`, "mean")),
      Lower = unlist(lapply(centroid_trajectories[[feature]], `[[`, "lower")),
      Upper = unlist(lapply(centroid_trajectories[[feature]], `[[`, "upper")),
      Risk_Group = rep(unique(risk_group$risk_group), each = length(common_time))
    )
    p <- ggplot(plot_data, aes(x = Time, y = Mean, color = Risk_Group, fill = Risk_Group)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
      labs(
        title = feature,
        x = "Time",
        y = feature
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    plot_list[[feature]] <- p
  }
  
  # Save plots as a PDF
  pdf(pdf_path, width = 10, height = 12)
  grid.arrange(grobs = plot_list, ncol = 2, nrow = 3)
  dev.off()
  
  return(centroid_trajectories)
}