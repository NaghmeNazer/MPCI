library(data.table)
library(caret)
library(pROC)
library(ggplot2)
library(doParallel)
library(readr)
library(tidyr)
library(dplyr)
library(patchwork)

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# Clean SVM evaluation function
svm_nested_cv <- function(data, kernel_type = "linear", outer_folds = 5, inner_folds = 3, metric_name = "Unknown") {
  # Ensure proper factor levels
  data$Label <- factor(data$Label, levels = c("PRE", "POST"))
  
  # Create balanced folds
  outer_folds <- createMultiFolds(data$Label, k = outer_folds)
  
  results <- data.frame()
  all_predictions <- data.frame()
  
  for (i in seq_along(outer_folds)) {
    train_data <- data[outer_folds[[i]], ]
    test_data <- data[-outer_folds[[i]], ]
    
    ctrl <- trainControl(
      method = "cv",
      number = inner_folds,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      verboseIter = FALSE,
      savePredictions = "final"
    )
    
    # Train model with probability calibration
    svm_model <- train(
      Label ~ .,
      data = train_data,
      method = ifelse(kernel_type == "linear", "svmLinear", "svmRadial"),
      tuneGrid = if(kernel_type == "linear") {
        expand.grid(C = 10^seq(-3, 3, length.out = 7))
      } else {
        expand.grid(C = 10^seq(-3, 3, length.out = 5),
                    sigma = 10^seq(-3, 1, length.out = 5))
      },
      trControl = ctrl,
      metric = "ROC",
      preProcess = c("center", "scale"),
      probability = TRUE
    )
    
    # Get probabilities and predictions
    pred_prob <- predict(svm_model, test_data, type = "prob")[,"POST"]
    pred_class <- predict(svm_model, test_data)
    
    # Store predictions for ROC aggregation
    fold_predictions <- data.frame(
      Metric = metric_name,
      Fold = i,
      Sample = rownames(test_data),
      True_Label = test_data$Label,
      Predicted_Probability = pred_prob,
      Predicted_Class = pred_class,
      stringsAsFactors = FALSE
    )
    all_predictions <- rbind(all_predictions, fold_predictions)
    
    # ROC curve with automatic direction detection
    roc_obj <- roc(
      response = test_data$Label,
      predictor = pred_prob,
      levels = levels(test_data$Label),
      quiet = TRUE
    )
    
    # Confusion matrix
    conf_mat <- confusionMatrix(
      pred_class,
      test_data$Label,
      positive = "POST",
      mode = "everything"
    )
    
    # Extract AUC as numeric
    auc_value <- as.numeric(auc(roc_obj))
    
    results <- rbind(results, data.frame(
      Metric = metric_name,
      Fold = i,
      Kernel = kernel_type,
      AUC = auc_value,
      Sensitivity = as.numeric(conf_mat$byClass["Sensitivity"]),
      Specificity = as.numeric(conf_mat$byClass["Specificity"]),
      Accuracy = as.numeric(conf_mat$overall["Accuracy"]),
      PPV = as.numeric(conf_mat$byClass["Pos Pred Value"]),
      NPV = as.numeric(conf_mat$byClass["Neg Pred Value"]),
      F1 = as.numeric(conf_mat$byClass["F1"]),
      Prob_Mean = mean(pred_prob),
      Prob_SD = sd(pred_prob),
      Best_C = svm_model$bestTune$C,
      Best_sigma = ifelse(kernel_type == "rbf", svm_model$bestTune$sigma, NA),
      stringsAsFactors = FALSE
    ))
  }
  
  return(list(results = results, predictions = all_predictions))
}

prepare_data_for_classification <- function(data_matrix){
  # Check column names
  cat("First few column names:", names(data_matrix)[1:5], "\n")
  
  # Prepare data for classification
  features <- as.matrix(data_matrix[, -1])  # Remove first column (ID)
  rownames(features) <- data_matrix[[1]]  # Use first column as ID
  
  # Transpose so samples are rows, DMRs are columns
  features_t <- t(features)
  
  # Identify PRE and POST samples
  sample_names <- rownames(features_t)
  is_pre <- grepl("PRE", sample_names, ignore.case = TRUE)
  is_post <- grepl("POST", sample_names, ignore.case = TRUE)
  
  # Create labels
  labels <- rep(NA, length(sample_names))
  labels[is_pre] <- "PRE"
  labels[is_post] <- "POST"
  
  # Keep only PRE/POST samples
  keep_samples <- !is.na(labels)
  features_t <- features_t[keep_samples, ]
  features_t <- as.data.frame(features_t)
  labels <- factor(labels[keep_samples], levels = c("PRE", "POST"))
  features_t$Label <- labels
  
  cat("Samples prepared:", nrow(features_t), "(PRE:", sum(labels == "PRE"), 
      "POST:", sum(labels == "POST"), ")\n")
  return(features_t)
}

clean_high_na_data <- function(data, col_threshold = 0.8, row_threshold = 0.8) {
  # Make a copy of the data
  cleaned_data <- data
  
  # Step 1: Remove columns with >threshold% NA
  na_col_prop <- colMeans(is.na(cleaned_data))
  cols_to_keep <- na_col_prop <= col_threshold
  cleaned_data <- cleaned_data[, cols_to_keep, drop = FALSE]
  
  cat("Removed", sum(!cols_to_keep), "columns with >", 
      col_threshold*100, "% NA values\n")
  
  # Step 2: Remove rows with >threshold% NA
  na_row_prop <- rowMeans(is.na(cleaned_data))
  rows_to_keep <- na_row_prop <= row_threshold
  cleaned_data <- cleaned_data[rows_to_keep, , drop = FALSE]
  
  cat("Removed", sum(!rows_to_keep), "rows with >", 
      row_threshold*100, "% NA values\n")
  
  return(cleaned_data)
}

# Enhanced function to aggregate results with statistical tests
aggregate_results <- function(mpci_results, dmhl_results) {
  # Combine both results
  all_results <- rbind(mpci_results$results, dmhl_results$results)
  all_predictions <- rbind(mpci_results$predictions, dmhl_results$predictions)
  
  # Ensure all numeric columns are actually numeric
  numeric_cols <- c("AUC", "Sensitivity", "Specificity", "Accuracy", "PPV", "NPV", "F1")
  for (col in numeric_cols) {
    if (col %in% names(all_results)) {
      all_results[[col]] <- as.numeric(all_results[[col]])
    }
  }
  
  # Calculate summary statistics
  summary_stats <- all_results %>%
    group_by(Metric) %>%
    summarise(
      AUC_mean = mean(AUC, na.rm = TRUE),
      AUC_sd = sd(AUC, na.rm = TRUE),
      AUC_se = sd(AUC, na.rm = TRUE) / sqrt(n()),
      Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
      Sensitivity_sd = sd(Sensitivity, na.rm = TRUE),
      Sensitivity_se = sd(Sensitivity, na.rm = TRUE) / sqrt(n()),
      Specificity_mean = mean(Specificity, na.rm = TRUE),
      Specificity_sd = sd(Specificity, na.rm = TRUE),
      Specificity_se = sd(Specificity, na.rm = TRUE) / sqrt(n()),
      Accuracy_mean = mean(Accuracy, na.rm = TRUE),
      Accuracy_sd = sd(Accuracy, na.rm = TRUE),
      Accuracy_se = sd(Accuracy, na.rm = TRUE) / sqrt(n()),
      n_folds = n(),
      .groups = 'drop'
    )
  
  # Calculate Wilcoxon test p-values for each metric
  stat_tests <- list()
  metrics_to_test <- c("AUC", "Sensitivity", "Specificity", "Accuracy")
  
  for (metric_type in metrics_to_test) {
    mpci_vals <- all_results[all_results$Metric == "MPCI", metric_type]
    dmhl_vals <- all_results[all_results$Metric == "dMHL", metric_type]
    
    if (length(mpci_vals) >= 3 && length(dmhl_vals) >= 3) {
      wilcox_test <- wilcox.test(mpci_vals, dmhl_vals, paired = TRUE, exact = FALSE)
      stat_tests[[metric_type]] <- wilcox_test$p.value
    } else {
      stat_tests[[metric_type]] <- NA
    }
  }
  
  return(list(
    all_results = all_results,
    all_predictions = all_predictions,
    summary_stats = summary_stats,
    stat_tests = stat_tests
  ))
}

# Function to create ROC data properly
create_roc_data <- function(all_predictions) {
  roc_data <- data.frame()
  
  for (method in c("MPCI", "dMHL")) {
    method_predictions <- all_predictions[all_predictions$Metric == method, ]
    
    # Create ROC object
    roc_obj <- roc(method_predictions$True_Label, 
                   method_predictions$Predicted_Probability,
                   levels = c("PRE", "POST"), 
                   direction = "<", 
                   quiet = TRUE)
    
    # Extract coordinates from ROC object
    roc_coords <- coords(roc_obj, transpose = FALSE, ret = c("specificity", "sensitivity"))
    
    # Create data frame
    temp_df <- data.frame(
      Method = method,
      FPR = 1 - roc_coords$specificity,
      TPR = roc_coords$sensitivity,
      AUC = as.numeric(auc(roc_obj))
    )
    
    roc_data <- rbind(roc_data, temp_df)
  }
  
  return(roc_data)
}

# Updated function to visualize aggregated results with publication-ready formatting
visualize_comparison <- function(aggregated_results) {
  
  # Extract components
  all_results <- aggregated_results$all_results
  summary_stats <- aggregated_results$summary_stats
  all_predictions <- aggregated_results$all_predictions
  stat_tests <- aggregated_results$stat_tests
  
  # Ensure all metrics are numeric
  all_results$AUC <- as.numeric(all_results$AUC)
  all_results$Accuracy <- as.numeric(all_results$Accuracy)
  all_results$Sensitivity <- as.numeric(all_results$Sensitivity)
  all_results$Specificity <- as.numeric(all_results$Specificity)
  
  # 1. Enhanced bar plot of performance metrics with error bars
  # Prepare data for plotting
  plot_data <- summary_stats %>%
    select(Metric, AUC_mean, Sensitivity_mean, Specificity_mean, Accuracy_mean,
           AUC_se, Sensitivity_se, Specificity_se, Accuracy_se) %>%
    pivot_longer(cols = -Metric, 
                 names_to = c("Metric_Type", ".value"),
                 names_pattern = "(.*)_(mean|se)$") %>%
    filter(!is.na(mean)) %>%
    mutate(
      Metric_Type = factor(Metric_Type, 
                           levels = c("AUC", "Sensitivity", "Specificity", "Accuracy"),
                           labels = c("AUC", "Sensitivity", "Specificity", "Accuracy"))
    )
  
  # Create base plot
  p1 <- ggplot(plot_data, aes(x = Metric_Type, y = mean, fill = Metric)) +
    geom_bar(stat = "identity", 
             position = position_dodge(width = 0.8), 
             width = 0.7,
             color = "black",
             size = 0.3) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.8),
                  width = 0.25,
                  size = 0.5,
                  color = "black") +
    geom_text(aes(label = sprintf("%.2f", mean), 
                  y = mean + se + 0.05),
              position = position_dodge(width = 0.8),
              size = 6,
              fontface = "bold") +
    labs(x = "Performance Metric", 
         y = "Value",
         fill = "Method") +
    scale_fill_manual(values = c("dMHL" = "#ff7f0e","MPCI" = "#1f77b4"),
                      labels = c("dMHL", "MPCI")) +
    scale_y_continuous(limits = c(0, 1.15),
                       breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0, 0.05))) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22, face = "bold"),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 18),
      legend.position = "right",
      legend.key.size = unit(1.5, "cm"),
      plot.margin = margin(20, 20, 20, 20),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      panel.grid.minor.y = element_blank()
    )
  
  # Add significance annotations
  max_heights <- plot_data %>%
    group_by(Metric_Type) %>%
    summarise(max_height = max(mean + se, na.rm = TRUE))
  
  sig_y_pos <- max_heights$max_height + 0.08
  p_values <- unlist(stat_tests)[c("AUC", "Sensitivity", "Specificity", "Accuracy")]
  
  # Add significance stars
  for (i in 1:length(p_values)) {
    if (!is.na(p_values[i])) {
      sig_star <- ifelse(p_values[i] < 0.001, "***",
                         ifelse(p_values[i] < 0.01, "**",
                                ifelse(p_values[i] < 0.05, "*", "")))
      
      if (sig_star != "") {
        p1 <- p1 + annotate("text",
                            x = i,
                            y = sig_y_pos[i],
                            label = sig_star,
                            size = 10,
                            fontface = "bold",
                            vjust = 0.5)
      }
    }
  }
  
  # 2. Enhanced ROC curve comparison with publication-ready formatting
  roc_data <- create_roc_data(all_predictions)
  
  p2 <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Method)) +
    geom_line(size = 1.5) +
    geom_abline(slope = 1, intercept = 0, 
                linetype = "dashed", 
                color = "gray50",
                size = 0.8) +
    labs(x = "False Positive Rate (1 - Specificity)", 
         y = "True Positive Rate (Sensitivity)",
         color = "Method") +
    scale_color_manual(values = c("MPCI" = "#1f77b4", "dMHL" = "#ff7f0e"),
                       labels = c("MPCI", "dMHL")) +
    coord_equal() +
    theme_classic(base_size = 24) +
    theme(
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 24, face = "bold"),
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 20),
      legend.position = "right",
      legend.key.size = unit(1.5, "cm"),
      plot.margin = margin(20, 20, 20, 20),
      panel.grid.major = element_line(color = "grey90", size = 0.3)
    )
  
  # Add AUC annotations with publication formatting
  auc_summary <- roc_data %>%
    group_by(Method) %>%
    summarise(AUC = first(AUC))
  
  p2 <- p2 + 
    annotate("text", x = 0.65, y = 0.45, 
             label = "Pooled ROC AUC :",
             color = "darkgrey", 
             size = 6,
             fontface = "bold",
             hjust = 0) + 
    annotate("text", x = 0.65, y = 0.35, 
             label = paste("MPCI = ", sprintf("%.3f", auc_summary$AUC[auc_summary$Method == "MPCI"])),
             color = "#1f77b4", 
             size = 6,
             fontface = "bold",
             hjust = 0) +
    annotate("text", x = 0.65, y = 0.25, 
             label = paste("dMHL = ", sprintf("%.3f", auc_summary$AUC[auc_summary$Method == "dMHL"])),
             color = "#ff7f0e", 
             size = 6,
             fontface = "bold",
             hjust = 0)
  
  # Combine plots with patchwork
  combined_plot <- (p1 + p2) +
    plot_layout(ncol = 2, widths = c(1.2, 1)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 28, face = "bold"))
  
  # Print statistical summary
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("PERFORMANCE SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  cat("Wilcoxon Signed-Rank Tests (paired, non-parametric):\n")
  cat("===================================================\n")
  for (metric in names(stat_tests)) {
    p_val <- stat_tests[[metric]]
    if (!is.na(p_val)) {
      significance <- ifelse(p_val < 0.001, "***",
                             ifelse(p_val < 0.01, "**",
                                    ifelse(p_val < 0.05, "*", "ns")))
      cat(sprintf("%-12s: p = %.4f %s\n", metric, p_val, significance))
    }
  }
  
  cat("\nMean Performance ± Standard Error:\n")
  cat("=================================\n")
  for (metric in c("MPCI", "dMHL")) {
    stats <- summary_stats[summary_stats$Metric == metric, ]
    cat(sprintf("\n%s:\n", metric))
    cat(sprintf("  AUC:          %.3f ± %.3f\n", stats$AUC_mean, stats$AUC_se))
    cat(sprintf("  Sensitivity:  %.3f ± %.3f\n", stats$Sensitivity_mean, stats$Sensitivity_se))
    cat(sprintf("  Specificity:  %.3f ± %.3f\n", stats$Specificity_mean, stats$Specificity_se))
    cat(sprintf("  Accuracy:     %.3f ± %.3f\n", stats$Accuracy_mean, stats$Accuracy_se))
  }
  
  return(list(
    combined_plot = combined_plot,
    individual_plots = list(p1 = p1, p2 = p2),
    roc_data = roc_data,
    summary_stats = summary_stats,
    stat_tests = stat_tests
  ))
}

# Updated main analysis function
run_full_analysis <- function() {
  cat("Loading data...\n")
  
  # Load data
  MPCI_raw <- read_csv("~/GSE262275/high_cov_regions/400_most_covered/MPCI.csv")
  dMHL_raw <- read_csv("~/GSE262275/high_cov_regions/400_most_covered/dMHL.csv")
  
  cat("\n=== Preparing MPCI data ===\n")
  MPCI <- prepare_data_for_classification(MPCI_raw)
  MPCI <- clean_high_na_data(MPCI, 0.1, 0.8)
  MPCI <- na.omit(MPCI)
  cat("Final MPCI data: ", nrow(MPCI), "samples, ", ncol(MPCI)-1, "features\n")
  
  cat("\n=== Preparing dMHL data ===\n")
  dMHL <- prepare_data_for_classification(dMHL_raw)
  dMHL <- clean_high_na_data(dMHL, 0.1, 0.8)
  dMHL <- na.omit(dMHL)
  cat("Final dMHL data: ", nrow(dMHL), "samples, ", ncol(dMHL)-1, "features\n")
  
  # Check if we have enough data
  if (nrow(MPCI) < 10 || nrow(dMHL) < 10) {
    cat("\nERROR: Insufficient data after preprocessing!\n")
    return(NULL)
  }
  
  cat("\n=== Running SVM classification with linear kernel ===\n")
  set.seed(123)
  
  # Run SVM for both methods
  cat("Running MPCI classification...\n")
  MPCI_results <- svm_nested_cv(MPCI, "linear", 5, 3, "MPCI")
  
  cat("Running dMHL classification...\n")
  dMHL_results <- svm_nested_cv(dMHL, "linear", 5, 3, "dMHL")
  
  # Aggregate results
  cat("\n=== Aggregating results ===\n")
  aggregated_results <- aggregate_results(MPCI_results, dMHL_results)
  
  # Visualize comparison
  cat("\n=== Creating visualizations ===\n")
  visualizations <- visualize_comparison(aggregated_results)
  
  # Save results with publication quality
  cat("\n=== Saving results ===\n")
  
  # Save high-resolution plots
  ggsave("MPCI_vs_dMHL_performance.png", 
         visualizations$individual_plots$p1,
         width = 12, height = 8, dpi = 600, bg = "white")
  
  ggsave("MPCI_vs_dMHL_ROC.png", 
         visualizations$individual_plots$p2,
         width = 10, height = 8, dpi = 600, bg = "white")
  
  ggsave("MPCI_vs_dMHL_combined.png", 
         visualizations$combined_plot,
         width = 22, height = 10, dpi = 600, bg = "white")
  
  # Save data files
  write_csv(aggregated_results$all_results, "classification_results.csv")
  write_csv(aggregated_results$summary_stats, "performance_summary.csv")
  
  # Save statistical tests
  stat_tests_df <- data.frame(
    Metric = names(aggregated_results$stat_tests),
    p_value = unlist(aggregated_results$stat_tests)
  )
  write_csv(stat_tests_df, "statistical_tests.csv")
  
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("ANALYSIS COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("\nFiles saved:\n")
  cat("1. MPCI_vs_dMHL_performance.png (Performance bar plot)\n")
  cat("2. MPCI_vs_dMHL_ROC.png (ROC curves)\n")
  cat("3. MPCI_vs_dMHL_combined.png (Combined figure)\n")
  cat("4. classification_results.csv (Raw classification results)\n")
  cat("5. performance_summary.csv (Summary statistics)\n")
  cat("6. statistical_tests.csv (Wilcoxon test results)\n")
  
  return(list(
    MPCI_results = MPCI_results,
    dMHL_results = dMHL_results,
    aggregated_results = aggregated_results,
    visualizations = visualizations
  ))
}

# Run the analysis
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("MPCI vs dMHL CLASSIFICATION COMPARISON ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

results <- run_full_analysis()

# Stop parallel cluster
stopCluster(cl)

# Print final summary
if (!is.null(results)) {
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("FINAL SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  mpci_auc <- mean(results$MPCI_results$results$AUC)
  dmhl_auc <- mean(results$dMHL_results$results$AUC)
  
  cat(sprintf("\nMean AUC:\n"))
  cat(sprintf("  MPCI:  %.3f\n", mpci_auc))
  cat(sprintf("  dMHL:  %.3f\n", dmhl_auc))
  cat(sprintf("  Difference: %.3f\n", mpci_auc - dmhl_auc))
  
  # Calculate confidence intervals
  mpci_se <- sd(results$MPCI_results$results$AUC) / sqrt(nrow(results$MPCI_results$results))
  dmhl_se <- sd(results$dMHL_results$results$AUC) / sqrt(nrow(results$dMHL_results$results))
  
  cat(sprintf("\n95%% Confidence Intervals:\n"))
  cat(sprintf("  MPCI:  %.3f [%.3f, %.3f]\n", 
              mpci_auc, mpci_auc - 1.96*mpci_se, mpci_auc + 1.96*mpci_se))
  cat(sprintf("  dMHL:  %.3f [%.3f, %.3f]\n", 
              dmhl_auc, dmhl_auc - 1.96*dmhl_se, dmhl_auc + 1.96*dmhl_se))
}

results[["visualizations"]]
