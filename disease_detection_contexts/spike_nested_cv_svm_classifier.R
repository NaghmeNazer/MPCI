library(data.table)
library(caret)
library(pROC)
library(ggplot2)
library(doParallel)
library(readr)

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

# Clean SVM evaluation function
svm_nested_cv <- function(data, kernel_type = "linear", outer_folds = 5, inner_folds = 3) {
  # Ensure proper factor levels
  data$Label <- factor(data$Label, levels = c("Case", "Healthy"))
  
  # Create balanced folds
  outer_folds <- createMultiFolds(data$Label, k = outer_folds)
  
  results <- data.frame()
  
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
      probability = TRUE  # Ensure proper probability estimation
    )
    
    # Get probabilities and predictions
    pred_prob <- predict(svm_model, test_data, type = "prob")[,"Case"]
    pred_class <- predict(svm_model, test_data)
    
    # ROC curve with automatic direction detection
    roc_obj <- roc(
      response = test_data$Label,
      predictor = pred_prob,
      levels = levels(test_data$Label),
      quiet = TRUE
    )
    
    # Check probability distribution
    prob_hist <- hist(pred_prob, plot = FALSE)
    
    # Confusion matrix
    conf_mat <- confusionMatrix(
      pred_class,
      test_data$Label,
      positive = "Case",
      mode = "everything"
    )
    
    results <- rbind(results, data.frame(
      Fold = i,
      Kernel = kernel_type,
      AUC = auc(roc_obj),
      Sensitivity = conf_mat$byClass["Sensitivity"],
      Specificity = conf_mat$byClass["Specificity"],
      Accuracy = conf_mat$overall["Accuracy"],
      Prob_Mean = mean(pred_prob),
      Prob_SD = sd(pred_prob),
      Best_C = svm_model$bestTune$C,
      Best_sigma = ifelse(kernel_type == "rbf", svm_model$bestTune$sigma, NA),
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}
# Data processing function
load_and_process_data <- function(spike_ratio, spike_data.tissue, target_data.tissue, min_cov.thresh = 10, sample_num = 100) {
  # Convert percentage to decimal for file paths
  spike_ratio_decimal <- spike_ratio/100
  results_dir <- paste0("../simulate_spike_separate_healthy/12-24-2025_spike_generation_results_1percent/4/",
                        spike_ratio_decimal, "_", spike_data.tissue, "_in_", target_data.tissue)
  
  # Initialize data structures
  metrics <- list(
    dMHL = list(healthy = NULL, case = NULL),
    MPCI = list(healthy = NULL, case = NULL)
  )
  
  # Load data for all samples
  for (i in 1:sample_num) {
    file_path <- paste0(results_dir, "/",
                        spike_data.tissue, "_in_", target_data.tissue,
                        "_spikeratio", spike_ratio_decimal, "_S", i, ".csv")
    
    if (!file.exists(file_path)) stop(paste("File not found:", file_path))
    
    metric_df <- fread(file_path)
    metric_df$V1 <- NULL
    metric_df$dMHL_cfdna <- metric_df$MHL_cfdna - metric_df$uMHL_cfdna
    metric_df$dMHL_neuron <- metric_df$MHL_neuron - metric_df$uMHL_neuron
    metric_df$dMHL_simulated <- metric_df$MHL_simulated - metric_df$uMHL_simulated
    
    if (i == 1) {
      common_cols <- c("DMR_ID", "initial_coverage_cfdna", 
                       "coverage_neuron", "coverage_simulated")
      
      for (metric in names(metrics)) {
        target_col <- switch(metric,
                             "dMHL" = "dMHL_cfdna",
                             "MPCI" = "MPCI_cfdna")
        result_col <- switch(metric,
                             "dMHL" = "dMHL_simulated",
                             "MPCI" = "MPCI_simulated")
        
        metrics[[metric]]$healthy <- metric_df[, c(common_cols, target_col), with = FALSE]
        metrics[[metric]]$case <- metric_df[, c(common_cols, result_col), with = FALSE]
        
        names(metrics[[metric]]$healthy)[ncol(metrics[[metric]]$healthy)] <- paste0(metric, "_Healthy_S", i)
        names(metrics[[metric]]$case)[ncol(metrics[[metric]]$case)] <- paste0(metric, "_Case_S", i)
      }
    } else {
      for (metric in names(metrics)) {
        target_col <- switch(metric,
                             "dMHL" = "dMHL_cfdna",
                             "MPCI" = "MPCI_cfdna")
        result_col <- switch(metric,
                             "dMHL" = "dMHL_simulated",
                             "MPCI" = "MPCI_simulated")
        
        metrics[[metric]]$healthy[, paste0(metric, "_Healthy_S", i) := metric_df[[target_col]]]
        metrics[[metric]]$case[, paste0(metric, "_Case_S", i) := metric_df[[result_col]]]
      }
    }
  }
  
  # Process each metric
  result_dfs <- list()
  
  for (metric in names(metrics)) {
    # Select regions with sufficient coverage
    select_regions <- metrics[[metric]]$case$initial_coverage_cfdna > min_cov.thresh
    regions <- metrics[[metric]]$case$DMR_ID[select_regions]
    regions <- regions[!is.na(regions)]
    
    # Combine case and healthy data
    case_cols <- grep("_Case_", names(metrics[[metric]]$case), value = TRUE)
    healthy_cols <- grep("_Healthy_", names(metrics[[metric]]$healthy), value = TRUE)
    
    case_data <- metrics[[metric]]$case[select_regions, ..case_cols]
    healthy_data <- metrics[[metric]]$healthy[select_regions, ..healthy_cols]
    
    total_df <- cbind(case_data, healthy_data)
    total_df <- total_df[complete.cases(total_df), ]
    
    # Remove columns with zero variance
    col_vars <- sapply(total_df, function(x) var(x, na.rm = TRUE))
    total_df <- total_df[, which(col_vars > 0), with = FALSE]
    
    # Transpose and add labels
    total_mat <- t(as.matrix(total_df))
    total_df <- as.data.frame(total_mat)
    total_df$Label <- ifelse(grepl("Case", rownames(total_df)), "Case", "Healthy")
    total_df$Label <- factor(total_df$Label)
    
    # Ensure column names match the number of regions
    valid_regions <- regions[1:ncol(total_mat)]  # Only use regions that survived filtering
    colnames(total_df) <- c(valid_regions, "Label")
    
    result_dfs[[metric]] <- total_df
  }
  
  return(result_dfs)
}

# Main analysis function
run_full_analysis <- function() {
  spike_ratios <- c(1, 2, 3, 4, 5, 10)
  spike_data.tissue <- "Neuron"
  target_data.tissue <- "cfDNA"
  
  all_results <- data.table()
  
  for (ratio in spike_ratios) {
    cat("\nProcessing spike ratio:", ratio, "%\n")
    
    data_list <- load_and_process_data(ratio, spike_data.tissue, target_data.tissue)
    
    for (metric_name in names(data_list)) {
      cat("  Analyzing metric:", metric_name, "\n")
      
      set.seed(123)
      
      # Run both kernels
      linear_results <- svm_nested_cv(data_list[[metric_name]], "linear", 5, 3)
      
      # Combine results
      combined <- cbind(linear_results, Spike_Ratio = ratio, Metric = metric_name)
      
      all_results <- rbindlist(list(all_results, combined), fill = TRUE)
    }
  }
  
  # Save complete results
  fwrite(all_results, "svm_results_complete.csv")
  
  # Create and save summary statistics
  summary_stats <- all_results[, .(
    AUC = mean(AUC),
    Sensitivity = mean(Sensitivity),
    Specificity = mean(Specificity),
    Accuracy = mean(Accuracy),
    Best_C = median(Best_C),
    Best_sigma = median(Best_sigma, na.rm = TRUE)
  ), by = .(Spike_Ratio, Metric, Kernel)]
  
  fwrite(summary_stats, "svm_results_summary.csv")
  
  return(list(full_results = all_results, summary = summary_stats))
}

# Execute analysis
results <- run_full_analysis()
#results_backup <- results


# Visualization function for four separate line plots
create_four_line_plots <- function(all_results) {
  
  # Calculate summary statistics with standard error
  plot_summary <- all_results[, .(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_se = sd(AUC, na.rm = TRUE) / sqrt(.N),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Sensitivity_se = sd(Sensitivity, na.rm = TRUE) / sqrt(.N),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    Specificity_se = sd(Specificity, na.rm = TRUE) / sqrt(.N),
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Accuracy_se = sd(Accuracy, na.rm = TRUE) / sqrt(.N),
    N = .N
  ), by = .(Spike_Ratio, Metric, Kernel)]
  
  # Ensure Spike_Ratio is ordered correctly
  plot_summary$Spike_Ratio <- factor(plot_summary$Spike_Ratio, 
                                     levels = sort(unique(plot_summary$Spike_Ratio)))
  
  # Create color palette for metrics
  metric_colors <- c("dMHL" = "#E41A1C", "MPCI" = "#377EB8")
  
  # 1. AUC Plot
  auc_plot <- ggplot(plot_summary, aes(x = Spike_Ratio, group = Metric)) +
    geom_line(aes(y = AUC_mean, color = Metric), size = 1.2) +
    geom_point(aes(y = AUC_mean, color = Metric), size = 3) +
    geom_errorbar(aes(y = AUC_mean, ymin = AUC_mean - AUC_se, 
                      ymax = AUC_mean + AUC_se, color = Metric), 
                  width = 0.1, size = 0.8) +
    labs(
      title = "AUC Values Across Spike Ratios",
      x = "Spike Ratio (%)",
      y = "AUC Value",
      color = "Metric"
    ) +
    scale_color_manual(values = metric_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # 2. Accuracy Plot
  accuracy_plot <- ggplot(plot_summary, aes(x = Spike_Ratio, group = Metric)) +
    geom_line(aes(y = Accuracy_mean, color = Metric), size = 1.2) +
    geom_point(aes(y = Accuracy_mean, color = Metric), size = 3) +
    geom_errorbar(aes(y = Accuracy_mean, ymin = Accuracy_mean - Accuracy_se, 
                      ymax = Accuracy_mean + Accuracy_se, color = Metric), 
                  width = 0.1, size = 0.8) +
    labs(
      title = "Accuracy Values Across Spike Ratios",
      x = "Spike Ratio (%)",
      y = "Accuracy Value",
      color = "Metric"
    ) +
    scale_color_manual(values = metric_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # 3. Specificity Plot
  specificity_plot <- ggplot(plot_summary, aes(x = Spike_Ratio, group = Metric)) +
    geom_line(aes(y = Specificity_mean, color = Metric), size = 1.2) +
    geom_point(aes(y = Specificity_mean, color = Metric), size = 3) +
    geom_errorbar(aes(y = Specificity_mean, ymin = Specificity_mean - Specificity_se, 
                      ymax = Specificity_mean + Specificity_se, color = Metric), 
                  width = 0.1, size = 0.8) +
    labs(
      title = "Specificity Values Across Spike Ratios",
      x = "Spike Ratio (%)",
      y = "Specificity Value",
      color = "Metric"
    ) +
    scale_color_manual(values = metric_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # 4. Sensitivity Plot
  sensitivity_plot <- ggplot(plot_summary, aes(x = Spike_Ratio, group = Metric)) +
    geom_line(aes(y = Sensitivity_mean, color = Metric), size = 1.2) +
    geom_point(aes(y = Sensitivity_mean, color = Metric), size = 3) +
    geom_errorbar(aes(y = Sensitivity_mean, ymin = Sensitivity_mean - Sensitivity_se, 
                      ymax = Sensitivity_mean + Sensitivity_se, color = Metric), 
                  width = 0.1, size = 0.8) +
    labs(
      title = "Sensitivity Values Across Spike Ratios",
      x = "Spike Ratio (%)",
      y = "Sensitivity Value",
      color = "Metric"
    ) +
    scale_color_manual(values = metric_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Create a 2x2 grid of all four plots
  library(patchwork)
  
  # Arrange plots in a 2x2 grid
  combined_plots <- (auc_plot | accuracy_plot) / (specificity_plot | sensitivity_plot)
  
  # Add a main title to the combined grid
  combined_plots <- combined_plots + 
    plot_annotation(
      title = "SVM Performance Metrics Across Spike Ratios",
      subtitle = "Error bars represent standard error",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14)
      )
    )
  
  # Display the combined plot
  print(combined_plots)
  
  # Save individual plots
  ggsave("auc_plot.png", auc_plot, width = 8, height = 6, dpi = 300)
  ggsave("accuracy_plot.png", accuracy_plot, width = 8, height = 6, dpi = 300)
  ggsave("specificity_plot.png", specificity_plot, width = 8, height = 6, dpi = 300)
  ggsave("sensitivity_plot.png", sensitivity_plot, width = 8, height = 6, dpi = 300)
  
  # Save the combined grid
  ggsave("all_metrics_grid.png", combined_plots, width = 14, height = 10, dpi = 300)
  
  # Return all plots as a list
  return(list(
    auc_plot = auc_plot,
    accuracy_plot = accuracy_plot,
    specificity_plot = specificity_plot,
    sensitivity_plot = sensitivity_plot,
    combined_grid = combined_plots
  ))
}

create_four_line_plots_smart <- function(all_results) {
  
  # Calculate summary statistics with standard error
  plot_summary <- all_results[, .(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_se = sd(AUC, na.rm = TRUE) / sqrt(.N),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Sensitivity_se = sd(Sensitivity, na.rm = TRUE) / sqrt(.N),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    Specificity_se = sd(Specificity, na.rm = TRUE) / sqrt(.N),
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Accuracy_se = sd(Accuracy, na.rm = TRUE) / sqrt(.N),
    N = .N
  ), by = .(Spike_Ratio, Metric, Kernel)]
  
  # Ensure Spike_Ratio is ordered correctly
  plot_summary$Spike_Ratio <- factor(plot_summary$Spike_Ratio, 
                                     levels = sort(unique(plot_summary$Spike_Ratio)))
  
  # Create color palette for metrics
  metric_colors <- c("dMHL" = "#E41A1C", "MPCI" = "#377EB8")
  
  # Calculate MINIMUM values including error bars (mean - SE)
  auc_min <- min(plot_summary$AUC_mean - plot_summary$AUC_se, na.rm = TRUE)
  sens_min <- min(plot_summary$Sensitivity_mean - plot_summary$Sensitivity_se, na.rm = TRUE)
  spec_min <- min(plot_summary$Specificity_mean - plot_summary$Specificity_se, na.rm = TRUE)
  acc_min <- min(plot_summary$Accuracy_mean - plot_summary$Accuracy_se, na.rm = TRUE)
  
  # Find the overall minimum (considering error bars)
  overall_min <- min(c(auc_min, sens_min, spec_min, acc_min), na.rm = TRUE)
  
  # Calculate ceiling to nearest 0.1 (round DOWN to 0.1 increments)
  y_min <- floor(overall_min * 10) / 10
  
  # Ensure we don't go below 0
  y_min <- max(0, y_min)
  
  # Also ensure we leave some room - if y_min is very close to 0, start at 0
  if (y_min < 0.1) {
    y_min <- 0
  }
  
  # Create dynamic breaks based on the range
  create_y_scale <- function() {
    # Calculate appropriate breaks based on the range
    y_range <- 1 - y_min
    
    if (y_range <= 0.2) {
      # Very small range - use fine breaks
      breaks <- seq(y_min, 1, by = 0.02)
      minor_breaks <- seq(y_min, 1, by = 0.01)
    } else if (y_range <= 0.4) {
      # Medium range
      breaks <- seq(y_min, 1, by = 0.05)
      minor_breaks <- seq(y_min, 1, by = 0.01)
    } else {
      # Large range
      breaks <- seq(y_min, 1, by = 0.1)
      minor_breaks <- seq(y_min, 1, by = 0.05)
    }
    
    # Filter breaks to be within range
    breaks <- breaks[breaks >= y_min & breaks <= 1]
    
    list(
      limits = c(y_min, 1),
      breaks = breaks,
      minor_breaks = minor_breaks[minor_breaks >= y_min & minor_breaks <= 1]
    )
  }
  
  y_scale_info <- create_y_scale()
  
  # Create plotting function with dynamic y-axis
  create_plot <- function(data, y_mean, y_se, title, y_label) {
    ggplot(data, aes(x = Spike_Ratio, group = Metric)) +
      geom_line(aes(y = .data[[y_mean]], color = Metric), size = 1.2) +
      geom_point(aes(y = .data[[y_mean]], color = Metric), size = 3) +
      geom_errorbar(aes(y = .data[[y_mean]], 
                        ymin = .data[[y_mean]] - .data[[y_se]], 
                        ymax = .data[[y_mean]] + .data[[y_se]], 
                        color = Metric), 
                    width = 0.1, size = 0.8) +
      labs(
        title = title,
        x = "Spike Ratio (%)",
        y = y_label,
        color = "Metric"
      ) +
      scale_color_manual(values = metric_colors) +
      scale_y_continuous(
        limits = y_scale_info$limits,
        breaks = y_scale_info$breaks,
        minor_breaks = y_scale_info$minor_breaks
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "bottom",
        panel.grid.minor = element_line(color = "grey90"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = margin(20, 20, 20, 20)
      )
  }
  
  # Create all four plots
  auc_plot <- create_plot(plot_summary, "AUC_mean", "AUC_se", 
                          "AUC Values Across Spike Ratios", "AUC Value")
  
  accuracy_plot <- create_plot(plot_summary, "Accuracy_mean", "Accuracy_se",
                               "Accuracy Values Across Spike Ratios", "Accuracy Value")
  
  specificity_plot <- create_plot(plot_summary, "Specificity_mean", "Specificity_se",
                                  "Specificity Values Across Spike Ratios", "Specificity Value")
  
  sensitivity_plot <- create_plot(plot_summary, "Sensitivity_mean", "Sensitivity_se",
                                  "Sensitivity Values Across Spike Ratios", "Sensitivity Value")
  
  # Create combined plot
  library(patchwork)
  combined_plots <- (auc_plot | accuracy_plot) / (specificity_plot | sensitivity_plot)
  
  # Add annotation
  combined_plots <- combined_plots + 
    plot_annotation(
      title = "SVM Performance Metrics Across Spike Ratios",
      subtitle = sprintf("Error bars: standard error | Y-axis: %.1f to 1.0 (min: %.3f)", 
                         y_min, overall_min),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )
  
  # Save with informative filename
  suffix <- sprintf("_yaxis_%.1f-1.0", y_min)
  filenames <- list(
    auc = paste0("auc_plot", suffix, ".png"),
    accuracy = paste0("accuracy_plot", suffix, ".png"),
    specificity = paste0("specificity_plot", suffix, ".png"),
    sensitivity = paste0("sensitivity_plot", suffix, ".png"),
    combined = paste0("all_metrics_grid", suffix, ".png")
  )
  
  ggsave(filenames$auc, auc_plot, width = 8, height = 6, dpi = 300)
  ggsave(filenames$accuracy, accuracy_plot, width = 8, height = 6, dpi = 300)
  ggsave(filenames$specificity, specificity_plot, width = 8, height = 6, dpi = 300)
  ggsave(filenames$sensitivity, sensitivity_plot, width = 8, height = 6, dpi = 300)
  ggsave(filenames$combined, combined_plots, width = 14, height = 10, dpi = 300)
  
  # Return results
  return(list(
    plots = list(
      auc_plot = auc_plot,
      accuracy_plot = accuracy_plot,
      specificity_plot = specificity_plot,
      sensitivity_plot = sensitivity_plot,
      combined_grid = combined_plots
    ),
    y_axis = list(
      minimum_value = overall_min,
      y_min = y_min,
      y_max = 1,
      breaks = y_scale_info$breaks,
      filenames = filenames
    )
  ))
}

# Generate visualizations
if(exists("results")) {
  # Assuming results$full_results contains the data
  #plots <- create_four_line_plots(results$full_results)
  plots <- create_four_line_plots_smart(results$full_results)
  
  # Display summary message
  cat("Visualizations created and saved as PNG files:\n")
  cat("1. auc_plot.png - AUC values across spike ratios\n")
  cat("2. accuracy_plot.png - Accuracy values across spike ratios\n")
  cat("3. specificity_plot.png - Specificity values across spike ratios\n")
  cat("4. sensitivity_plot.png - Sensitivity values across spike ratios\n")
  cat("5. all_metrics_grid.png - Combined 2x2 grid of all plots\n")
}
