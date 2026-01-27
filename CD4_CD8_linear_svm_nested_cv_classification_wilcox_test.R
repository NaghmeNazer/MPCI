# Load required libraries
library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(grid)
library(parallel)
library(ggpubr)
library(broom)

# =============================================================================
# PARAMETERS CONFIGURATION
# =============================================================================

cat("=== PARAMETERS CONFIGURATION ===\n")

# Cross-validation parameters
PARAMS <- list(
  outer_folds = 5,
  inner_folds = 3,
  cv_repeats = 10,
  random_seed = 123,
  
  # SVM hyperparameter grids (linear kernel only)
  svm_linear_C = 10^seq(-3, 3, length.out = 7),
  
  # Visualization parameters
  plot_width = 15,
  plot_height = 10,
  plot_dpi = 300,
  point_size = 2,
  text_size = 10,
  
  # Color scheme
  colors_cell_types = c("#00AFBB", "#FC4E07"),
  colors_metrics = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                     "#8c564b", "#7f7f7f", "#bcbd22", "#38dccf")
)

# Set random seed for reproducibility
set.seed(PARAMS$random_seed)

# =============================================================================
# REGIONS TO REMOVE (by their ID strings)
# =============================================================================

# Directly specify the region IDs to remove
regions_to_remove <- c(
  "chr22-26188589-26188599",
  "chr22-26188578-26188593", 
  "chr22-26269670-26269683",
  "chr22-26191684-26191699",
  "chr22-26603546-26603558",
  "chr22-26269658-26269672",
  "chr22-26603554-26603571",
  "chr22-26191667-26191687"
)


# =============================================================================
# STEP 1: DATA LOADING AND CLEANING WITH REGION FILTERING
# =============================================================================

cat("\n=== STEP 1: DATA LOADING AND CLEANING ===\n")

# Read the generated CSV files
metric_files <- list(
  MPCI = "../random_regions/chr22/CD4_CD8/MPCI.csv",
  MHL = "../random_regions/chr22/CD4_CD8/MHL.csv",
  uMHL = "../random_regions/chr22/CD4_CD8/uMHL.csv",
  dMHL = "../random_regions/chr22/CD4_CD8/dMHL.csv",
  ME = "../random_regions/chr22/CD4_CD8/ME.csv",
  Epipoly = "../random_regions/chr22/CD4_CD8/Epipoly.csv",
  PDR = "../random_regions/chr22/CD4_CD8/PDR.csv",
  FDRP = "../random_regions/chr22/CD4_CD8/FDRP.csv",
  qFDRP = "../random_regions/chr22/CD4_CD8/qFDRP.csv"
)

# Function to safely read CSV files
read_metric_file <- function(file_path, metric_name) {
  if (file.exists(file_path)) {
    cat("Loading", metric_name, "...\n")
    return(read.csv(file_path))
  } else {
    cat("Warning: File not found -", file_path, "\n")
    return(NULL)
  }
}

# Load all metric data
metrics_raw <- map2(metric_files, names(metric_files), read_metric_file)
metrics_raw <- metrics_raw[!sapply(metrics_raw, is.null)]

# Function to simplify column names
simplify_colnames <- function(col_names) {
  simplified <- sapply(col_names, function(x) {
    if (grepl("chr", x)) {
      result <- sub(".*(chr.*?)\\.(GSM.*?)(\\.Z.*|$)", "\\1.\\2", x)
      result <- sub("\\.$", "", result)
      return(result)
    } else {
      return(x)
    }
  })
  return(simplified)
}

# Function to extract cell type from simplified column names
extract_cell_type <- function(col_names) {
  cell_types <- sapply(col_names, function(x) {
    if (grepl("\\.CD4", x) || grepl("CD4\\.", x)) {
      return("CD4")
    } else if (grepl("\\.CD8", x) || grepl("CD8\\.", x)) {
      return("CD8")
    } else {
      return(NA)
    }
  })
  return(cell_types)
}

# Function to preprocess data with region filtering
preprocess_metric_data <- function(metric_df) {
  if (is.null(metric_df) || ncol(metric_df) <= 1) return(NULL)
  
  # Simplify column names
  original_cols <- colnames(metric_df)
  simplified_cols <- c(original_cols[1], simplify_colnames(original_cols[-1]))
  colnames(metric_df) <- simplified_cols
  
  # Filter out specified regions
  cat("  Filtering specified regions...\n")
  original_count <- nrow(metric_df)
  metric_df_filtered <- metric_df[!metric_df$ID %in% regions_to_remove, ]
  removed_count <- original_count - nrow(metric_df_filtered)
  
  if (removed_count > 0) {
    cat("  Removed", removed_count, "specified regions\n")
    # Show which regions were actually removed
    removed_regions <- setdiff(metric_df$ID, metric_df_filtered$ID)
    cat("  Removed regions:", paste(removed_regions, collapse=", "), "\n")
  }
  
  # Remove regions with any NA values
  if (ncol(metric_df_filtered) > 1) {
    na_mask <- apply(metric_df_filtered[, -1, drop = FALSE], 1, function(x) !any(is.na(x)))
    metric_df_clean <- metric_df_filtered[na_mask, ]
  } else {
    metric_df_clean <- metric_df_filtered
  }
  
  return(metric_df_clean)
}

# Preprocess all metric dataframes
cat("Preprocessing data...\n")
metrics_clean <- map(metrics_raw, preprocess_metric_data)
metrics_clean <- metrics_clean[!sapply(metrics_clean, is.null)]

# Print preprocessing results
cat("\nOriginal vs Clean dimensions:\n")
walk2(names(metrics_clean), metrics_clean, function(name, df) {
  original_rows <- nrow(metrics_raw[[name]])
  clean_rows <- nrow(df)
  cat(name, ": ", original_rows, "->", clean_rows, "regions (", 
      round(clean_rows/original_rows * 100, 1), "% retained)\n")
})

# Prepare cell type labels from first available metric
if (length(metrics_clean) > 0) {
  first_metric <- metrics_clean[[1]]
  cell_names <- colnames(first_metric)[-1]
  cell_types <- extract_cell_type(cell_names)
  
  cat("\nCell type distribution:\n")
  print(table(cell_types))
} else {
  stop("No data available after preprocessing!")
}

# =============================================================================
# STEP 2: PCA ANALYSIS
# =============================================================================

cat("\n=== STEP 2: PCA ANALYSIS ===\n")

# Function to prepare data for analysis
prepare_data_for_analysis <- function(metric_df, cell_types) {
  data <- t(metric_df[, -1, drop = FALSE])
  colnames(data) <- metric_df$ID
  rownames(data) <- colnames(metric_df)[-1]
  data_clean <- data[complete.cases(data), , drop = FALSE]
  cell_types_clean <- cell_types[complete.cases(data)]
  return(list(data = data_clean, cell_types = cell_types_clean))
}

# Generate PCA plots
pca_plots <- list()
for (metric_name in names(metrics_clean)) {
  cat("Performing PCA for", metric_name, "\n")
  
  prepared_data <- prepare_data_for_analysis(metrics_clean[[metric_name]], cell_types)
  
  if (nrow(prepared_data$data) > 1 && length(unique(prepared_data$cell_types)) >= 2) {
    pca_result <- PCA(prepared_data$data, graph = FALSE)
    
    p <- fviz_pca_ind(pca_result,
                      geom.ind = "point",
                      col.ind = prepared_data$cell_types,
                      palette = PARAMS$colors_cell_types,
                      addEllipses = TRUE,
                      legend.title = "Cell Type",
                      title = paste(metric_name, "\n(", nrow(prepared_data$data), 
                                    "cells,", ncol(prepared_data$data), "regions)"),
                      pointsize = PARAMS$point_size) +
      theme(plot.title = element_text(size = PARAMS$text_size, hjust = 0.5),
            legend.position = "none")
    
    pca_plots[[metric_name]] <- p
  } else {
    cat("  Not enough data for PCA:", nrow(prepared_data$data), "cells\n")
  }
}

# Combine all plots into one grid if we have any
if (length(pca_plots) > 0) {
  n_plots <- length(pca_plots)
  ncol <- ifelse(n_plots >= 4, 3, 2)
  nrow <- ceiling(n_plots / ncol)
  
  combined_pca_plot <- grid.arrange(
    grobs = pca_plots,
    ncol = ncol,
    nrow = nrow,
    top = textGrob("PCA Analysis of Multiple Methylation Metrics", 
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )
  
  ggsave("../random_regions/12-13-2025_classifier_results/CD4_CD8_all_metrics_combined_pca_plots.png", combined_pca_plot, 
         width = PARAMS$plot_width, height = 5 * nrow, dpi = PARAMS$plot_dpi)
  cat("PCA plots saved as 'CD4_CD8_all_metrics_combined_pca_plots.png'\n")
} else {
  cat("No PCA plots generated - insufficient data\n")
}

# =============================================================================
# STEP 3: PREPARE DATA FOR CLASSIFICATION
# =============================================================================

cat("\n=== STEP 3: PREPARE DATA FOR CLASSIFICATION ===\n")

# Function to create data matrix for classification
create_data_matrix <- function(metric_df, cell_types, metric_name) {
  if (is.null(metric_df) || ncol(metric_df) <= 1) return(NULL)
  
  data_mat <- t(metric_df[, -1, drop = FALSE])
  colnames(data_mat) <- metric_df$ID
  data_df <- as.data.frame(data_mat)
  data_df$Sample <- rownames(data_mat)
  data_df$Label <- factor(cell_types)
  
  # Remove any rows with NA values
  data_clean <- data_df[complete.cases(data_df), ]
  
  cat(metric_name, ": ", nrow(data_clean), "samples,", 
      ncol(data_clean) - 2, "features\n")
  
  return(data_clean)
}

# Create data matrices for all metrics
data_matrices <- map2(metrics_clean, names(metrics_clean), 
                      ~create_data_matrix(.x, cell_types, .y))
data_matrices <- data_matrices[!sapply(data_matrices, is.null)]

# =============================================================================
# STEP 4: NESTED CV FUNCTION WITH LINEAR KERNEL ONLY
# =============================================================================

cat("\n=== STEP 4: RUNNING CLASSIFICATION ANALYSIS (LINEAR KERNEL) ===\n")

# SVM nested CV function with linear kernel
svm_nested_cv_linear <- function(data_matrix, metric_name) {
  # Remove Sample column and ensure Label is factor
  data <- data_matrix %>% select(-Sample) 
  data$Label <- factor(data$Label)
  
  # Check if we have at least 2 classes
  if (nlevels(data$Label) < 2) {
    warning(paste("Only one class present in", metric_name))
    return(data.frame())
  }
  
  results <- data.frame()
  
  for(rep in 1:PARAMS$cv_repeats) {
    set.seed(PARAMS$random_seed + rep)
    outer_folds <- createMultiFolds(data$Label, k = PARAMS$outer_folds, times = 1)
    
    for (i in seq_along(outer_folds)) {
      train_data <- data[outer_folds[[i]], ]
      test_data <- data[-outer_folds[[i]], ]
      
      # Check if both train and test have at least 2 classes
      if (nlevels(train_data$Label) < 2 || nlevels(test_data$Label) < 2) {
        next
      }
      
      # Set up control parameters
      ctrl <- trainControl(
        method = "cv",
        number = PARAMS$inner_folds,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        verboseIter = FALSE
      )
      
      # Set up grid for linear kernel only
      svm_grid <- expand.grid(C = PARAMS$svm_linear_C)
      
      # Train model with error handling
      svm_model <- tryCatch({
        train(
          Label ~ .,
          data = train_data,
          method = "svmLinear",
          tuneGrid = svm_grid,
          trControl = ctrl,
          metric = "ROC",
          preProcess = c("center", "scale")
        )
      }, error = function(e) {
        cat("  Error training model:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(svm_model)) next
      
      # Make predictions
      pred_prob <- tryCatch({
        predict(svm_model, test_data, type = "prob")[,2]
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(pred_prob)) next
      
      pred_class <- predict(svm_model, test_data)
      
      # Calculate metrics
      roc_obj <- tryCatch(roc(test_data$Label, pred_prob), error = function(e) NULL)
      conf_mat <- tryCatch(confusionMatrix(pred_class, test_data$Label), error = function(e) NULL)
      
      if (is.null(roc_obj) || is.null(conf_mat)) next
      
      results <- rbind(results, data.frame(
        Repeat = rep,
        Fold = i,
        Kernel = "linear",
        Metric = metric_name,
        AUC = auc(roc_obj),
        Accuracy = conf_mat$overall["Accuracy"],
        Sensitivity = conf_mat$byClass["Sensitivity"],
        Specificity = conf_mat$byClass["Specificity"]
      ))
    }
  }
  
  return(results)
}

# Run classification for all metrics with linear kernel only
classification_results <- data.frame()

for (metric_name in names(data_matrices)) {
  cat("Processing", metric_name, "with linear kernel...\n")
  
  linear_results <- svm_nested_cv_linear(data_matrices[[metric_name]], metric_name)
  
  if (nrow(linear_results) > 0) {
    classification_results <- rbind(classification_results, linear_results)
  }
}

# =============================================================================
# STEP 5: RESULTS ANALYSIS AND VISUALIZATION
# =============================================================================

cat("\n=== STEP 5: RESULTS ANALYSIS ===\n")

# Calculate summary statistics
summary_results <- classification_results %>%
  group_by(Metric) %>%
  summarise(
    N = n(),
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SD_AUC = sd(AUC, na.rm = TRUE),
    Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    SD_Accuracy = sd(Accuracy, na.rm = TRUE),
    Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
    SD_Sensitivity = sd(Sensitivity, na.rm = TRUE),
    Mean_Specificity = mean(Specificity, na.rm = TRUE),
    SD_Specificity = sd(Specificity, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nPerformance Summary (Linear Kernel Only):\n")
print(summary_results, n = Inf)

# =============================================================================
# STEP 6: NONPARAMETRIC STATISTICAL ANALYSIS BETWEEN METRICS
# =============================================================================

cat("\n=== STEP 6: NONPARAMETRIC STATISTICAL ANALYSIS ===\n")

# Function to perform nonparametric statistical comparisons
perform_nonparametric_comparisons <- function(classification_results) {
  test_results <- data.frame()
  
  metrics_to_test <- c("AUC", "Accuracy", "Sensitivity", "Specificity")
  
  # Get all metric pairs for comparison
  metric_pairs <- combn(unique(classification_results$Metric), 2, simplify = FALSE)
  
  for(metric in metrics_to_test) {
    for(pair in metric_pairs) {
      metric1 <- pair[1]
      metric2 <- pair[2]
      
      # Extract values for both metrics
      values1 <- classification_results %>% 
        filter(Metric == metric1) %>% 
        pull(!!sym(metric))
      
      values2 <- classification_results %>% 
        filter(Metric == metric2) %>% 
        pull(!!sym(metric))
      
      # Skip if not enough data
      if (length(values1) < 2 || length(values2) < 2) {
        cat("  Skipping", metric1, "vs", metric2, "for", metric, "- insufficient data\n")
        next
      }
      
      # Perform Wilcoxon rank-sum test (nonparametric alternative to t-test)
      wilcox_test <- wilcox.test(values1, values2, paired = FALSE, exact = FALSE)
      
      # Store results
      test_results <- rbind(test_results, data.frame(
        Metric_Type = metric,
        Metric1 = metric1,
        Metric2 = metric2,
        Mean_Metric1 = mean(values1),
        Mean_Metric2 = mean(values2),
        Median_Metric1 = median(values1),
        Median_Metric2 = median(values2),
        Difference = mean(values1) - mean(values2),
        W_Statistic = wilcox_test$statistic,
        P_Value = wilcox_test$p.value,
        Significance = ifelse(wilcox_test$p.value < 0.001, "***",
                              ifelse(wilcox_test$p.value < 0.01, "**",
                                     ifelse(wilcox_test$p.value < 0.05, "*", "NS")))
      ))
      
      cat("  ", metric1, "vs", metric2, "for", metric, ": p =", 
          format(wilcox_test$p.value, digits = 4), "\n")
    }
  }
  
  return(test_results)
}

# Perform nonparametric statistical tests
stat_test_results <- perform_nonparametric_comparisons(classification_results)

# Focus on MPCI vs dMHL comparisons
mpci_vs_dmhl <- stat_test_results %>%
  filter((Metric1 == "MPCI" & Metric2 == "dMHL") | (Metric1 == "dMHL" & Metric2 == "MPCI")) %>%
  mutate(Comparison = "MPCI_vs_dMHL")

cat("\nMPCI vs dMHL Nonparametric Statistical Comparisons:\n")
print(mpci_vs_dmhl)

# Save results
write.csv(classification_results, "../random_regions/12-13-2025_classifier_results/CD4_CD8_all_metrics_classification_results_linear.csv", row.names = FALSE)
write.csv(stat_test_results, "../random_regions/12-13-2025_classifier_results/CD4_CD8_all_metrics_statistical_comparisons_nonparametric.csv", row.names = FALSE)
write.csv(mpci_vs_dmhl, "../random_regions/12-13-2025_classifier_results/CD4_CD8_mpci_vs_dmhl_statistical_comparisons_nonparametric.csv", row.names = FALSE)

# =============================================================================
# STEP 7: COMPREHENSIVE VISUALIZATION
# =============================================================================

cat("\n=== STEP 7: CREATING VISUALIZATIONS ===\n")

# Set font sizes for all plots
large_text_theme <- theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 14, hjust = 0.5),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 11),
  strip.text = element_text(size = 12)
)

# Consistent color for all bars in AUC plot
consistent_color <- "#1f77b4"  # Professional blue color

# Calculate standard error function
std_error <- function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))

# Function to format p-values nicely
format_pvalue <- function(p_value) {
  if (p_value < 0.0001) {
    return("<0.0001")
  } else if (p_value < 0.001) {
    return(sprintf("%.4f", p_value))
  } else if (p_value < 0.01) {
    return(sprintf("%.4f", p_value))
  } else {
    return(sprintf("%.4f", p_value))
  }
}

# 1. Bar plot comparing AUCs with standard errors and sorted values
if (nrow(classification_results) > 0) {
  cat("Creating AUC comparison bar plot...\n")
  
  # Prepare AUC data with standard errors
  auc_summary <- classification_results %>%
    group_by(Metric) %>%
    summarise(
      Mean_AUC = mean(AUC, na.rm = TRUE),
      SE_AUC = std_error(AUC),
      N = n(),
      .groups = "drop"
    ) %>%
    arrange(Mean_AUC)  # Sort by AUC value
  
  auc_summary <- auc_summary[auc_summary$Metric != "dMHL",]
  
  # Convert Metric to factor with sorted levels for proper plotting
  auc_summary$Metric <- factor(auc_summary$Metric, levels = auc_summary$Metric)
  
  # Create the sorted AUC bar plot with standard errors
  auc_bar_plot <- ggplot(auc_summary, aes(x = Metric, y = Mean_AUC)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7, fill = consistent_color) +
    geom_errorbar(aes(ymin = Mean_AUC - SE_AUC, ymax = Mean_AUC + SE_AUC),
                  width = 0.2, color = "black", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Mean_AUC), y = Mean_AUC + SE_AUC + 0.02),
              size = 4, vjust = 0) +
    labs(title = "AUC Comparison Across Methylation Metrics",
         subtitle = "Linear SVM Kernel - Sorted by Performance (Mean ± SE)",
         x = "Methylation Metric", 
         y = "Area Under Curve (AUC)") +
    theme_bw() +
    large_text_theme +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank()) +
    ylim(0, 1.1)  # Ensure space for labels
  
  ggsave("../random_regions/12-13-2025_classifier_results/CD4_CD8_sorted_auc_comparison_se_linear.png", auc_bar_plot, 
         width = PARAMS$plot_width-7, height = PARAMS$plot_height-3, dpi = PARAMS$plot_dpi)
  cat("Saved: CD4_CD8_sorted_auc_comparison_se_linear.png\n")
  
  # 2. MPCI vs dMHL detailed comparison with statistical significance
  if (all(c("MPCI", "dMHL") %in% classification_results$Metric)) {
    cat("Creating MPCI vs dMHL detailed comparison...\n")
    
    mpci_dmhl_data <- classification_results %>%
      filter(Metric %in% c("MPCI", "dMHL")) %>%
      mutate(
        AUC = as.numeric(AUC),
        Accuracy = as.numeric(Accuracy),
        Sensitivity = as.numeric(Sensitivity),
        Specificity = as.numeric(Specificity)
      ) %>%
      select(Metric, AUC, Accuracy, Sensitivity, Specificity) %>%
      pivot_longer(cols = c(AUC, Accuracy, Sensitivity, Specificity), 
                   names_to = "Metric_Type", values_to = "Value")
    
    # Calculate summary statistics with standard errors
    mpci_dmhl_summary <- mpci_dmhl_data %>%
      group_by(Metric, Metric_Type) %>%
      summarise(
        Mean_Value = mean(Value, na.rm = TRUE),
        SE_Value = std_error(Value),
        N = n(),
        .groups = "drop"
      )
    
    # Get significance labels from statistical tests with formatted p-values
    significance_data <- mpci_vs_dmhl %>%
      select(Metric_Type, P_Value, Significance) %>%
      mutate(
        Metric_Type = case_when(
          Metric_Type == "AUC" ~ "AUC",
          Metric_Type == "Accuracy" ~ "Accuracy", 
          Metric_Type == "Sensitivity" ~ "Sensitivity",
          Metric_Type == "Specificity" ~ "Specificity"
        ),
        Formatted_P = sapply(P_Value, format_pvalue),
        Label = paste("p =", Formatted_P, Significance)
      )
    
    # Create the comparison plot
    mpci_dmhl_comparison <- ggplot(mpci_dmhl_summary, aes(x = Metric, y = Mean_Value, fill = Metric)) +
      geom_bar(stat = "identity", alpha = 0.8, width = 0.6, position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = Mean_Value - SE_Value, ymax = Mean_Value + SE_Value),
                    width = 0.2, position = position_dodge(0.8), color = "black", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.3f", Mean_Value), y = Mean_Value + SE_Value + 0.03),
                size = 4, position = position_dodge(0.8), vjust = 0) +
      # Add significance annotations with formatted p-values
      geom_text(data = significance_data, 
                aes(x = 1.5, y = 1.05,  # Fixed y-position for significance labels
                    label = Label),
                inherit.aes = FALSE, size = 4, fontface = "bold", check_overlap = TRUE) +
      facet_wrap(~Metric_Type, scales = "free_y", nrow = 1) +
      labs(title = "MPCI vs dMHL Performance Comparison (dMHL = MHL - uMHL)",
           subtitle = "Linear SVM Kernel - Mean ± Standard Error",
           x = "Methylation Metric", 
           y = "Performance Value") +
      scale_fill_manual(values = c("dMHL" = "#ff7f0e", "MPCI" = "#1f77b4")) +
      theme_bw() +
      large_text_theme +
      theme(legend.position = "bottom",
            strip.background = element_rect(fill = "lightgray"),
            panel.grid.major.x = element_blank()) +
      scale_y_continuous(limits = c(0, 1.1))  # Ensure consistent y-axis
    
    ggsave("../random_regions/12-13-2025_classifier_results/CD4_CD8_mpci_vs_dmhl_detailed_comparison_linear.png", mpci_dmhl_comparison, 
           width = PARAMS$plot_width-7, height = PARAMS$plot_height-3, dpi = PARAMS$plot_dpi)
    cat("Saved: CD4_CD8_mpci_vs_dmhl_detailed_comparison_linear.png\n")
  } else {
    cat("MPCI or dMHL not available for comparison plot\n")
  }
} else {
  cat("No classification results available for visualization\n")
}

