library(parallel)
library(doParallel)
library(patchwork)

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
    
    # ── only compute dynamic stats if dynamic_features is non-empty ──────
    if (length(dynamic_features) > 0) {
      last_record  <- tail(subset, 1)
      first_record <- head(subset, 1)
      duration     <- last_record$stop - first_record$start
      
      for (feature in dynamic_features) {
        summary_stats[[paste(feature, 'min',    sep = '_')]] <- min(subset[[feature]])
        summary_stats[[paste(feature, 'max',    sep = '_')]] <- max(subset[[feature]])
        summary_stats[[paste(feature, 'median', sep = '_')]] <- median(subset[[feature]])
        summary_stats[[paste(feature, 'std',    sep = '_')]] <- sd(subset[[feature]])
        
        if (duration != 0) {
          slope <- (last_record[[feature]] - first_record[[feature]]) / duration
        } else {
          slope <- NA
        }
        summary_stats[[paste(feature, 'slope', sep = '_')]] <- slope
      }
      
      time <- final_record_stop - last_record$start
      
    } else {
      # no dynamic features — use final record start as time reference
      time <- final_record_stop - final_record$start
    }
    
    static_data <- list()
    for (feature in static_features) {
      static_data[[feature]] <- subset[[1, feature]]
    }
    
    aggregated_record <- list(
      ID    = ID,
      time  = time,
      event = final_record_event
    )
    aggregated_record <- c(aggregated_record, static_data, summary_stats)
    aggregated_records[[length(aggregated_records) + 1]] <- aggregated_record
  }
  
  df_aggregated <- do.call(rbind, lapply(aggregated_records, as.data.frame))
  return(df_aggregated)
}

# aggregate_for_standard_cox <- function(df_processed, static_features, dynamic_features, m) {
#   aggregated_records <- list()
#   grouped <- split(df_processed, df_processed$ID)
#   
#   for (ID in names(grouped)) {
#     group <- grouped[[ID]]
#     group <- group[order(group$start), ]
#     n <- nrow(group)
#     final_record <- tail(group, 1)
#     final_record_stop <- final_record$stop
#     final_record_event <- final_record$event
#     
#     if (m != 'all') {
#       subset <- tail(group, m)
#     } else {
#       subset <- group
#     }
#     
#     summary_stats <- list()
#     for (feature in dynamic_features) {
#       summary_stats[[paste(feature, 'min', sep = '_')]] <- min(subset[[feature]])
#       summary_stats[[paste(feature, 'max', sep = '_')]] <- max(subset[[feature]])
#       summary_stats[[paste(feature, 'median', sep = '_')]] <- median(subset[[feature]])
#       summary_stats[[paste(feature, 'std', sep = '_')]] <- sd(subset[[feature]])
#       
#       last_record <- tail(subset, 1)
#       first_record <- head(subset, 1)
#       duration <- last_record$stop - first_record$start
#       if (duration != 0) {
#         slope <- (last_record[[feature]] - first_record[[feature]]) / duration
#       } else {
#         slope <- NA  # Handle division by zero if duration is zero
#       }
#       summary_stats[[paste(feature, 'slope', sep = '_')]] <- slope
#     }
#     
#     static_data <- list()
#     for (feature in static_features) {
#       static_data[[feature]] <- subset[[1, feature]]
#     }
#     
# 
#     time <- final_record_stop - last_record$start
#     aggregated_record <- list(
#       ID = ID,
#       time = time,
#       event = final_record_event
#     )
#     aggregated_record <- c(aggregated_record, static_data, summary_stats)
#     aggregated_records[[length(aggregated_records) + 1]] <- aggregated_record
#   }
#   df_aggregated <- do.call(rbind, lapply(aggregated_records, as.data.frame))
#   return(df_aggregated)
# }

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

      # if (feature_selection) {
      #   formula_no_frailty <- as.formula(paste("Surv(start, stop, event) ~", paste(predictors, collapse = " + ")))
      #   cox_model_no_frailty <- coxph(formula_no_frailty, data = fold_train_data, ties = "efron", control = coxph.control(iter.max = 100, eps = 1e-05))
      #   stepwise_model <- stepAIC(cox_model_no_frailty, direction = "both", trace = FALSE)
      #   predictors <- all.vars(formula(stepwise_model))[-1]
      #   predictors <- setdiff(predictors, c("start", "stop", "event", "ID"))
      # }
      
      if (feature_selection) {
        formula_no_frailty <- as.formula(paste("Surv(start, stop, event) ~", paste(predictors, collapse = " + ")))
        cox_model_no_frailty <- coxph(formula_no_frailty, data = fold_train_data, ties = "efron", control = coxph.control(iter.max = 100, eps = 1e-05))
        stepwise_model <- stepAIC(cox_model_no_frailty, direction = "both", trace = FALSE)
        predictors <- all.vars(formula(stepwise_model))[-1]
        predictors <- setdiff(predictors, c("start", "stop", "event", "ID"))
        
        # VIF check: iteratively remove highest VIF feature until all VIF < 5
        repeat {
          if (length(predictors) <= 1) break
          
          vif_formula <- as.formula(paste("Surv(start, stop, event) ~", paste(predictors, collapse = " + ")))
          vif_model <- coxph(vif_formula, data = fold_train_data, ties = "efron", 
                             control = coxph.control(iter.max = 100, eps = 1e-05))
          vif_values <- tryCatch(vif(vif_model), error = function(e) NULL)
          
          if (is.null(vif_values)) break
          
          max_vif <- max(vif_values, na.rm = TRUE)
          if (max_vif < 5) break
          
          # Remove the predictor with the highest VIF
          worst_predictor <- names(which.max(vif_values))
          predictors <- setdiff(predictors, worst_predictor)
        }
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
  
  # return(list(
  #   C_Index = c_index,
  #   AUC_list = unlist(auc_summary$AUC),
  #   AUC_Stats = auc_summary,
  #   model_summary = summary(cox_model),
  #   risk_train = final_risk_score_train,
  #   risk_test = final_risk_score_test
  # ))
  return(list(
    C_Index = c_index,
    AUC_list = unlist(auc_summary$AUC),
    AUC_Stats = auc_summary,
    AUC_ROC = auc_result,        # raw timeROC object for ROC curve plotting
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
    
    # Wrap in tryCatch to identify exactly why a task fails
    tryCatch({
      unique_ids <- unique(data$ID)
      test_ids <- sample(unique_ids, size = floor(test_ratio * length(unique_ids)))
      train_ids <- setdiff(unique_ids, test_ids)
      train_data <- data[data$ID %in% train_ids, , drop = FALSE]
      test_data <- data[data$ID %in% test_ids, , drop = FALSE]
      
      fold_assignments <- sample(rep(1:folds, length.out = length(train_ids)))
      names(fold_assignments) <- train_ids
      best_c_index <- -1
      best_features <- NULL
      
      for (fold in 1:folds) {
        train_fold_ids <- train_ids[fold_assignments != fold]
        val_fold_ids <- train_ids[fold_assignments == fold]
        fold_train_data <- train_data[train_data$ID %in% train_fold_ids, , drop = FALSE]
        fold_validation_data <- train_data[train_data$ID %in% val_fold_ids, , drop = FALSE]
        
        if (!is.null(static_ft) || !is.null(dynamic_ft)) {
          fold_train_data <- aggregate_for_standard_cox(fold_train_data, static_ft, dynamic_ft, m = "all")
          
          # Check for columns that are ALL NA
          fold_train_data_cols_valid <- names(fold_train_data)[colSums(!is.na(fold_train_data)) > 0]
          predictors <- setdiff(fold_train_data_cols_valid, c("time", "event", "ID"))
          
          if (shorten_visits) {
            fold_validation_data <- aggregate_for_standard_cox(fold_validation_data, static_ft, dynamic_ft, m = n_visits)
          } else {
            fold_validation_data <- aggregate_for_standard_cox(fold_validation_data, static_ft, dynamic_ft, m = "all")
          }
          
          fold_validation_data_cols_valid <- names(fold_validation_data)[colSums(!is.na(fold_validation_data)) > 0]
          predictors <- intersect(predictors, setdiff(fold_validation_data_cols_valid, c("time", "event", "ID")))
        } else {
          predictors <- setdiff(names(fold_train_data), c("time", "event", "ID"))
        }
        
        # SAFETY: Check if we actually have predictors left
        if(length(predictors) == 0) stop("No valid predictors remaining after NA filtering.")
        
        # Preprocess and ensure output remains a data frame
        data_preprocess <- remove_high_corr_scales(fold_train_data, fold_validation_data, cols_consider = predictors, threshold = 0.9)
        fold_train_data <- as.data.frame(data_preprocess$train)
        fold_validation_data <- as.data.frame(data_preprocess$test)
        
        # Re-sync predictors after correlation removal
        predictors <- intersect(predictors, names(fold_train_data))
        
        formula_str <- paste("Surv(time, event) ~", paste(predictors, collapse = " + "))
        cox_model <- coxph(as.formula(formula_str), data = fold_train_data)
        
        # if (feature_selection) {
        #   stepwise_model <- suppressWarnings(stepAIC(cox_model, direction = "both", trace = FALSE, steps = 5000, k = 10))
        #   predictors <- setdiff(all.vars(formula(stepwise_model))[-1], c("time", "event", "ID"))
        # }
        
        if (feature_selection) {
          stepwise_model <- suppressWarnings(stepAIC(cox_model, direction = "both", trace = FALSE, steps = 5000, k = 10))
          predictors <- setdiff(all.vars(formula(stepwise_model))[-1], c("time", "event", "ID"))
          
          # VIF check: iteratively remove highest VIF feature until all VIF < 5
          repeat {
            if (length(predictors) <= 1) break
            
            vif_formula <- as.formula(paste("Surv(time, event) ~", paste(predictors, collapse = " + ")))
            vif_model <- coxph(vif_formula, data = fold_train_data)
            vif_values <- tryCatch(vif(vif_model), error = function(e) NULL)
            
            if (is.null(vif_values)) break
            
            max_vif <- max(vif_values, na.rm = TRUE)
            if (max_vif < 5) break
            
            # Remove the predictor with the highest VIF
            worst_predictor <- names(which.max(vif_values))
            predictors <- setdiff(predictors, worst_predictor)
          }
        }
        
        # Validation check
        c_index <- concordance(cox_model, newdata = fold_validation_data)$concordance
        
        if (c_index > best_c_index) {
          best_c_index <- c_index
          best_features <- predictors
        }
      }
      
      # Final Model fitting on full train set
      # train_data_agg <- aggregate_for_standard_cox(train_data, static_ft, dynamic_ft, m = "all")
      # if (shorten_visits) {
      #   test_data_agg <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = n_visits)
      # } else {
      #   test_data_agg <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = "all")
      # }
      # 
      
      if (!is.null(static_ft) || !is.null(dynamic_ft)) {
        train_data_agg <- aggregate_for_standard_cox(train_data, static_ft, dynamic_ft, m = "all")
        if (shorten_visits) {
          test_data_agg <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = n_visits)
        } else {
          test_data_agg <- aggregate_for_standard_cox(test_data, static_ft, dynamic_ft, m = "all")
        }
      } else {
        train_data_agg <- train_data
        test_data_agg  <- test_data
      }
      
      
      final_formula <- as.formula(paste("Surv(time, event) ~", paste(best_features, collapse = " + ")))
      final_model <- coxph(final_formula, data = train_data_agg)
      
      final_risk_scores <- predict(final_model, newdata = test_data_agg, type = "risk")
      final_c_index <- concordance(final_model, newdata = test_data_agg)$concordance
      
      # TimeROC evaluation
      times_to_eval <- sort(unique(test_data_agg$time[test_data_agg$event == 1]))
      
      auc_result <- timeROC(
        T = test_data_agg$time,
        delta = test_data_agg$event,
        marker = final_risk_scores,
        cause = 1,
        times = times_to_eval
      )
      
      list(
        C_Index = final_c_index,
        AUC_DF = data.frame(Iteration = iter, Time = auc_result$times, AUC = auc_result$AUC),
        Predictors = best_features,
        Status = "Success"
      )
      
    }, error = function(e) {
      # This will return the error message instead of crashing the whole loop
      list(Iteration = iter, Error = e$message, Status = "Failed")
    })
  }
  
  successful <- Filter(function(x) x$Status == "Success", results)
  if (length(successful) == 0) {
    print(results[[1]]$Error)
    stop("All iterations failed.")
  }
  
  all_c_index_list   <- sapply(successful, function(x) x$C_Index)
  all_auc_results_df <- do.call(rbind, lapply(successful, function(x) x$AUC_DF))
  all_predictors     <- unlist(lapply(successful, function(x) x$Predictors))
  
  # all_c_index_list <- sapply(results, function(x) x$C_Index)
  # all_auc_results_df <- do.call(rbind, lapply(results, function(x) x$AUC_DF))
  # all_predictors <- unlist(lapply(results, function(x) x$Predictors))
  
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
  
  # return(list(
  #   C_Index = c_index,
  #   AUC_list = unlist(auc_summary$AUC),
  #   AUC_Stats = auc_summary,
  #   model_summary = summary(cox_model)
  # ))
  
  return(list(
    C_Index = c_index,
    AUC_list = unlist(auc_summary$AUC),
    AUC_Stats = auc_summary,
    AUC_ROC = auc_result,        # raw timeROC object for ROC curve plotting
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

create_forest_plot <- function(summary_results, 
                               title = "Forest Plot for Cox Model", 
                               file_path_name = NULL,
                               feature_rename = NULL) {
  # feature_rename: named vector e.g. c("old_name" = "New Label")
  
  coefficients   <- summary_results$coefficients
  conf_intervals <- summary_results$conf.int
  
  # exclude frailty term
  if ("frailty(ID)" %in% rownames(coefficients)) {
    coefficients   <- coefficients[rownames(coefficients) != "frailty(ID)", ]
    conf_intervals <- conf_intervals[rownames(conf_intervals) != "frailty(ID)", ]
  }
  
  forest_data <- data.frame(
    Predictor = rownames(coefficients),
    HR        = conf_intervals[, "exp(coef)"],
    Lower_CI  = conf_intervals[, "lower .95"],
    Upper_CI  = conf_intervals[, "upper .95"],
    p_value   = coefficients[, "p"],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      # apply rename map if provided
      Predictor = ifelse(
        Predictor %in% names(feature_rename),
        feature_rename[Predictor],
        Predictor
      ),
      sig = case_when(
        p_value < 0.001 ~ "p < 0.001",
        p_value < 0.01  ~ "p < 0.01",
        p_value < 0.05  ~ "p < 0.05",
        TRUE            ~ "ns"
      ),
      CI_text = sprintf("%.3f (%.3f, %.3f)", HR, Lower_CI, Upper_CI),
      CI_text_star = ifelse(
        p_value < 0.05,
        paste0(CI_text, " *"),
        CI_text
      )
    ) %>%
    arrange(HR) %>%
    mutate(Predictor = factor(Predictor, levels = Predictor))
  
  # ── forest plot panel ────────────────────────────────────────────────────
  p_forest <- ggplot(forest_data,
                     aes(y = Predictor, x = HR,
                         xmin = Lower_CI, xmax = Upper_CI)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    geom_errorbarh(
      aes(colour = p_value < 0.05),
      height    = 0.25,
      linewidth = 0.8
    ) +
    geom_point(
      aes(colour = p_value < 0.05,
          shape  = p_value < 0.05),
      size = 3.5
    ) +
    scale_colour_manual(
      values = c("FALSE" = "black", "TRUE" = "black"),
      labels = c("FALSE" = "ns", "TRUE" = "p < 0.05"),
      name   = "Significance"
    ) +
    scale_shape_manual(
      values = c("FALSE" = 1, "TRUE" = 16),
      labels = c("FALSE" = "ns", "TRUE" = "p < 0.05"),
      name   = "Significance"
    ) +
    labs(
      x = "Hazard Ratio (HR)",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.title.x      = element_text(size = 14, face = "bold"),
      axis.text.x       = element_text(size = 13),
      axis.text.y       = element_text(size = 13),
      legend.position   = "bottom",
      legend.title      = element_text(size = 12, face = "bold"),
      legend.text       = element_text(size = 11)
    )
  
  # ── text panel (HR and CI on the right) ─────────────────────────────────
  p_text <- ggplot(forest_data,
                   aes(x = 1, y = Predictor, label = CI_text_star)) +
    geom_text(hjust = 0, size = 3.5, colour = "black") +
    xlim(1, 2) +
    labs(x = "HR (95% CI)", y = NULL) +
    theme_void(base_size = 12) +
    theme(
      axis.title.x     = element_text(size = 11, face = "bold",
                                      hjust = 0, vjust = 1),
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  
  # ── combine ──────────────────────────────────────────────────────────────
  p_final <- p_forest + p_text +
    plot_layout(widths = c(3, 1.5))
  
  if (!is.null(file_path_name)) {
    ggsave(file_path_name, plot = p_final,
           # width = 12, height = max(4, nrow(forest_data) * 0.5 + 2),
           width = 7, height = 5, 
           device = "pdf")
  }
  
  return(p_final)
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

get_sig_stars <- function(p_value) {
  if (p_value < 0.001) return("***")
  else if (p_value < 0.01) return("**")
  else if (p_value < 0.05) return("*")
  else return("ns")
}

plot_centroid_trajectories <- function(df, risk_group, selected_features_trajectories, 
                                       pdf_path,
                                       feature_rename = NULL,
                                       palette = c("#d4241f", "#3778ac")) {
  
  common_time <- seq(min(df$start), max(df$start), length.out = 100)
  centroid_trajectories <- list()
  
  # ── compute trajectories ─────────────────────────────────────────────────
  for (feature in selected_features_trajectories) {
    centroid_trajectories[[feature]] <- lapply(
      unique(risk_group$risk_group), function(group) {
        group_ids  <- risk_group$ID[risk_group$risk_group == group]
        group_data <- df %>%
          filter(ID %in% group_ids) %>%
          select(ID, start, all_of(feature))
        
        feature_matrix <- do.call(rbind, lapply(
          split(group_data, group_data$ID), function(d) {
            approx(x = d$start, y = d[[feature]],
                   xout = common_time, rule = 2)$y
          }
        ))
        
        list(
          mean  = colMeans(feature_matrix, na.rm = TRUE),
          lower = apply(feature_matrix, 2, function(x)
            mean(x, na.rm = TRUE) - 1.96 * sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
          upper = apply(feature_matrix, 2, function(x)
            mean(x, na.rm = TRUE) + 1.96 * sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
        )
      }
    )
  }
  
  # ── build plots ───────────────────────────────────────────────────────────
  risk_groups <- unique(risk_group$risk_group)
  plot_list   <- list()
  
  for (feature in selected_features_trajectories) {
    
    # apply rename if provided
    feature_label <- if (!is.null(feature_rename) && feature %in% names(feature_rename)) {
      feature_rename[[feature]]
    } else {
      feature
    }
    
    plot_data <- data.frame(
      Time       = rep(common_time, times = length(risk_groups)),
      Mean       = unlist(lapply(centroid_trajectories[[feature]], `[[`, "mean")),
      Lower      = unlist(lapply(centroid_trajectories[[feature]], `[[`, "lower")),
      Upper      = unlist(lapply(centroid_trajectories[[feature]], `[[`, "upper")),
      Risk_Group = rep(risk_groups, each = length(common_time))
    )
    
    p <- ggplot(plot_data,
                aes(x = Time, y = Mean,
                    colour = Risk_Group, fill = Risk_Group)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper),
                  alpha = 0.15, colour = NA) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = setNames(palette, risk_groups),
                          name   = "Risk Group") +
      scale_fill_manual(values   = setNames(palette, risk_groups),
                        name     = "Risk Group") +
      labs(
        title = feature_label,
        x     = "Time (months)",
        # y     = feature_label
        y     = " "
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title       = element_text(size = 15),
        axis.text        = element_text(size = 15),
        legend.position  = "bottom",
        legend.title     = element_text(size = 12, face = "bold"),
        legend.text      = element_text(size = 12),
        legend.key.size  = unit(0.4, "cm")
      )
    
    plot_list[[feature]] <- p
  }
  
  # # ── save as PDF ───────────────────────────────────────────────────────────
  # n_features <- length(selected_features_trajectories)
  # n_cols     <- 2
  # n_rows     <- ceiling(n_features / n_cols)
  # 
  # pdf(pdf_path, width = 10, height = n_rows * 4)
  # grid.arrange(grobs = plot_list, ncol = n_cols, nrow = n_rows)
  # dev.off()
  # ── save as PDF ───────────────────────────────────────────────────────────
  n_features <- length(selected_features_trajectories)
  n_cols     <- n_features    # all in one row
  n_rows     <- 1
  

  pdf(pdf_path, width = n_features * 4, height = 4)  # width scales with number of plots
  # pdf(pdf_path, width = 14, height = 4)  # width scales with number of plots
  grid.arrange(grobs = plot_list, ncol = n_cols, nrow = n_rows)
  dev.off()
  
  return(centroid_trajectories)
}

plot_boxplot_with_anno <- function(df, y_label, save_path,
                                   compare_group1, compare_group2,
                                   fill_var = "model",
                                   palette = NULL) {
  
  df[[fill_var]] <- factor(df[[fill_var]])
  df$feature_set <- factor(df$feature_set)
  
  # ── significance annotation (internal groups only) ───────────────────────
  # internal_levels <- levels(df$feature_set)[levels(df$feature_set) != "external (AIBL)"]
  internal_levels <- levels(df$feature_set)
  
  stat_anno <- do.call(rbind, lapply(internal_levels, function(fs) {
    grp1 <- df %>% filter(feature_set == fs, .data[[fill_var]] == compare_group1) %>% pull(value)
    grp2 <- df %>% filter(feature_set == fs, .data[[fill_var]] == compare_group2) %>% pull(value)
    
    if (length(grp1) < 3 || length(grp2) < 3) {
      return(data.frame(feature_set = fs, x = NA, xmin = NA, xmax = NA,
                        p_value = NA, stars = "ns", stringsAsFactors = FALSE))
    }
    
    p_val <- wilcox.test(grp1, grp2)$p.value
    stars <- get_sig_stars(p_val)
    x_pos <- as.numeric(factor(fs, levels = levels(df$feature_set)))
    
    n_groups  <- nlevels(df[[fill_var]])
    half_span <- 0.85 / 2 * ((n_groups - 1) / n_groups)
    
    data.frame(
      feature_set = fs,
      x           = x_pos,
      xmin        = x_pos - half_span,
      xmax        = x_pos + half_span,
      p_value     = p_val,
      stars       = stars,
      stringsAsFactors = FALSE
    )
  })) %>% filter(!is.na(x))
  
  # ── y positions for brackets ─────────────────────────────────────────────
  y_max_df <- df %>%
    group_by(feature_set) %>%
    summarise(y_max = max(value, na.rm = TRUE), .groups = "drop")
  
  y_range <- max(y_max_df$y_max, na.rm = TRUE) - min(df$value, na.rm = TRUE)
  y_step  <- y_range * 0.06
  
  stat_anno <- stat_anno %>%
    left_join(y_max_df, by = "feature_set") %>%
    mutate(
      y_bracket = y_max + y_step,
      y_text    = y_max + y_step * 1.8
    )
  
  # ── palette ───────────────────────────────────────────────────────────────
  # if (is.null(palette)) {
  #   if (fill_var == "model") {
  #     palette <- c("#d4241f", "#3778ac", "#4ea648", "#8e4b98")
  #   } else {
  #     palette <- c("#a05327", "#ed7b1c", "#e6df84", "#979797")
  #   }
  # }
  # if (is.null(palette)) {
  #   if (fill_var == "model") {
  #     palette <- c("#313772", "#2c4ca0", "#3c7fb1", "#75b5dc")
  #   } else {
  #     palette <- c("#376439", "#4d7e54", "#669877", "#81b095")
  #   }
  # }
  # 
  if (is.null(palette)) {
    if (fill_var == "model") {
      palette <- c("#96BDE4", "#67816D", "#AF7798", "#A8A1CA")
    } else {
      palette <- c("#7B92C7", "#F7C1CF", "#9bbf8a", "#CFAFD4")
    }
  }
  
  # ── whisker caps ─────────────────────────────────────────────────────────
  whisker_df <- df %>%
    group_by(feature_set, .data[[fill_var]]) %>%
    summarise(
      lower_whisker = boxplot.stats(value)$stats[1],
      upper_whisker = boxplot.stats(value)$stats[5],
      .groups = "drop"
    ) %>%
    mutate(
      fill_num  = as.numeric(.data[[fill_var]]),
      n_groups  = nlevels(df[[fill_var]]),
      fs_num    = as.numeric(factor(feature_set, levels = levels(df$feature_set))),
      dodge_x   = fs_num + (fill_num - (n_groups + 1) / 2) * (0.85 / n_groups),
      cap_width = 0.85 / n_groups * 0.4
    )
  
  # ── separate alpha for internal vs external ───────────────────────────────
  df <- df %>%
    mutate(is_external = ifelse(feature_set == "external (AIBL)", "external", "internal"))
  
  # ── plot ──────────────────────────────────────────────────────────────────
  p <- ggplot(data = df, aes(x = feature_set, y = value)) +
    
    # internal boxes — solid fill
    geom_boxplot(
      data     = df %>% filter(is_external == "internal"),
      aes(fill = .data[[fill_var]]),
      position      = position_dodge(width = 0.85),
      outlier.shape = NA,
      linewidth     = 0.8,
      width         = 0.7,
      fatten        = 1.5,
      alpha         = 1.0
    ) +
    
    # external boxes — same colour, lighter fill
    geom_boxplot(
      data     = df %>% filter(is_external == "external"),
      aes(fill = .data[[fill_var]]),
      position      = position_dodge(width = 0.85),
      outlier.shape = NA,
      linewidth     = 0.8,
      width         = 0.7,
      fatten        = 1.5,
      alpha         = 1
      # alpha         = 0.4          # lighter to distinguish from internal
    ) +
    
    # whisker cap lines lower
    geom_segment(
      data = whisker_df,
      aes(
        x    = dodge_x - cap_width,
        xend = dodge_x + cap_width,
        y    = lower_whisker,
        yend = lower_whisker
      ),
      linewidth = 0.9, inherit.aes = FALSE
    ) +
    
    # whisker cap lines upper
    geom_segment(
      data = whisker_df,
      aes(
        x    = dodge_x - cap_width,
        xend = dodge_x + cap_width,
        y    = upper_whisker,
        yend = upper_whisker
      ),
      linewidth = 0.9, inherit.aes = FALSE
    ) +
    
    # significance bracket horizontal
    geom_segment(
      data = stat_anno,
      aes(x = xmin, xend = xmax, y = y_bracket, yend = y_bracket),
      linewidth = 0.7, inherit.aes = FALSE
    ) +
    # left tick
    geom_segment(
      data = stat_anno,
      aes(x = xmin, xend = xmin,
          y = y_bracket - y_step * 0.3, yend = y_bracket),
      linewidth = 0.7, inherit.aes = FALSE
    ) +
    # right tick
    geom_segment(
      data = stat_anno,
      aes(x = xmax, xend = xmax,
          y = y_bracket - y_step * 0.3, yend = y_bracket),
      linewidth = 0.7, inherit.aes = FALSE
    ) +
    # stars
    geom_text(
      data = stat_anno,
      aes(x = x, y = y_text, label = stars),
      size = 4, fontface = "bold", inherit.aes = FALSE
    ) +
    
    # vertical dashed separator before external group
    # geom_vline(
    #   xintercept = length(internal_levels) + 0.5,
    #   linetype   = "dashed",
    #   colour     = "grey50",
    #   linewidth  = 0.5
    # ) +
    
    scale_fill_manual(values = palette) +
    scale_x_discrete(expand = expansion(mult = c(0.3, 0.15))) +
    coord_cartesian(ylim = c(
      min(df$value, na.rm = TRUE) - y_step,
      max(stat_anno$y_text, na.rm = TRUE) + y_step
    )) +
    
    theme_minimal(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 0.6),
      axis.ticks       = element_line(colour = "black", linewidth = 0.5),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = 12, face = "bold"),
      axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 14),
      axis.text.y      = element_text(size = 15),
      legend.title     = element_text(size = 15, face = "bold"),
      legend.text      = element_text(size = 14),
      legend.position  = "bottom",
      legend.key.size  = unit(0.4, "cm")
    ) +
    labs(
      y    = y_label,
      fill = ifelse(fill_var == "model", "Model", "Visit Count")
    )
  
  ggsave(filename = save_path, plot = p, device = "pdf", width = 14, height = 5)
  return(p)
}


