# GSE262275 Analysis Pipeline

This repository contains scripts for processing PAT files from the GSE262275 dataset, identifying high-coverage genomic regions, generating methylation metrics, and performing SVM-based classification to compare the performance of MPCI and dMHL in distinguishing between PRE and POST treatment conditions.

________________________________________
## Requirements

- R (version 4.0 or higher)
- Required R packages: `data.table`, `tidyverse`, `caret`, `e1071`, `pROC`, `ggplot2`, `doParallel`, `readr`, `tidyr`, `dplyr`, `patchwork`

________________________________________
## Repository Structure

### Region Identification Scripts

#### 1. **GSE262275_high_cov_regions_finding.R**
Identifies high-coverage regions from PAT files for downstream analysis.

**Functions:**

- **read_pat_data(path)**: Reads and processes PAT files
  - Input: Path to PAT file
  - Output: Data.table with chrom, start_cpg_index, methylation_pattern, count

- **find_high_coverage_regions(pat_data, min_reads = 10, min_cpgs = 3)**: Identifies consecutive CpGs with sufficient coverage using integer counting
  - Input: PAT data, minimum read threshold, minimum CpG sites
  - Output: Data.table of non-overlapping high-coverage regions with region scores

- **find_high_coverage_regions_fractional(pat_data, min_reads = 10, min_cpgs = 3)**: Identifies high-coverage regions using fractional coverage assignment
  - Input: PAT data, minimum read threshold, minimum CpG sites
  - Output: Data.table of non-overlapping high-coverage regions with fractional coverage calculation

________________________________________
### Metric Generation Scripts

#### 2. **GSE262275_metric_generation.R**
Generates methylation metrics for identified high-coverage regions across all samples.

**Core Binary Matrix Functions:**
- **subset_pat_in_window()**: Extracts reads overlapping specified genomic windows
- **truncated_to_binary()**: Converts methylation patterns to binary matrices
- **handle_excess_NAs_in_binary_DMR()**: Cleans binary matrices by removing all-NA rows/columns

**Methylation Metric Functions:**

- **give_sign_for_weights(x_1, x_2)**: Assigns weights based on average methylation status of two rows
  - Input: Two binary vectors (methylated = 1, unmethylated = 0)
  - Output: 1 (>50% methylated), -1 (<50% methylated), or random (±1 if exactly 50%)

- **signed_manhattan_sim(binary_dmr)**: Calculates signed Manhattan similarity
  - Input: Binary matrix (rows = reads, columns = CpG sites)
  - Output: Weighted average of pairwise similarities

- **MPCI(binary_dmr)**: Calculates Methylation Pattern Consistency Index
  - Input: Binary matrix
  - Output: Value from -1 to +1
    - Positive values: consistent methylation
    - Negative values: consistent unmethylation
    - Near zero: random methylation patterns

- **MHL(binary_dmr, max_hap_l = 9)**: Calculates Methylation Haplotype Load
  - Input: Binary matrix, maximum haplotype length
  - Output: Value from 0 to 1 indicating methylation block consistency

- **uMHL(binary_dmr, max_hap_l = 9)**: Calculates unmethylated MHL
  - Input: Binary matrix
  - Output: Value from 0 to 1 for unmethylated patterns

- **simple_uMHL(binary_dmr)**: Simplified uMHL calculation using complement matrix
  - Input: Binary matrix
  - Output: uMHL value

- **dMHL(MHL, uMHL)**: Calculates differential MHL (MHL - uMHL)
  - Input: MHL and uMHL values
  - Output: Difference metric

**Output Files:**
- `MPCI.csv`: Methylation Pattern Consistency Index values
- `MHL.csv`: Methylation Haplotype Load values
- `uMHL.csv`: Unmethylated haplotype load
- `dMHL.csv`: Differential MHL values
- `coverage.csv`: Read coverage per region

________________________________________
### Classification and Analysis Scripts

#### 3. **GSE262275_svm_nested_classifier_wilcoxon.R**
Performs SVM classification with nested cross-validation to compare MPCI and dMHL performance.

**Key Functions:**

- **svm_nested_cv(data, kernel_type = "linear", outer_folds = 5, inner_folds = 3, metric_name = "Unknown")**: 
  - Purpose: Runs nested cross-validation with SVM
  - Input: Prepared data matrix, kernel type, fold parameters
  - Output: List containing fold results and predictions

- **prepare_data_for_classification(data_matrix)**: 
  - Purpose: Transforms metric data into format suitable for classification
  - Input: Raw metric CSV with samples as columns
  - Output: Data frame with features and PRE/POST labels

- **clean_high_na_data(data, col_threshold = 0.1, row_threshold = 0.8)**:
  - Purpose: Removes features and samples with excessive missing values
  - Input: Data matrix, NA thresholds
  - Output: Cleaned data matrix

- **aggregate_results(mpci_results, dmhl_results)**:
  - Purpose: Combines results from both metrics and performs statistical tests
  - Input: Results from MPCI and dMHL SVM analyses
  - Output: List with combined results, predictions, summary statistics, and Wilcoxon test p-values

- **create_roc_data(all_predictions)**:
  - Purpose: Generates ROC curve coordinates from prediction probabilities
  - Input: Combined predictions from both methods
  - Output: Data frame with FPR, TPR, and AUC values

- **visualize_comparison(aggregated_results)**:
  - Purpose: Creates publication-ready visualizations
  - Output: List containing:
    - Combined plot with performance bars and ROC curves
    - Individual plot components
    - ROC data
    - Summary statistics
    - Statistical test results

- **run_full_analysis()**:
  - Purpose: Orchestrates complete analysis pipeline
  - Output: Comprehensive results including all metrics, visualizations, and statistical tests

________________________________________
## Output Files

### Visualization Outputs
- `MPCI_vs_dMHL_performance.png`: Bar plot comparing performance metrics (AUC, Sensitivity, Specificity, Accuracy) with error bars and significance stars
- `MPCI_vs_dMHL_ROC.png`: ROC curves comparing MPCI and dMHL classification performance
- `MPCI_vs_dMHL_combined.png`: Combined figure with both performance bar plot and ROC curves

### Data Outputs
- `classification_results.csv`: Raw classification results from nested CV folds
- `performance_summary.csv`: Summary statistics (mean, SD, SE) for all performance metrics
- `statistical_tests.csv`: Wilcoxon signed-rank test p-values comparing MPCI and dMHL

________________________________________
## Statistical Analysis

The pipeline implements statistical comparisons:

- **Paired Wilcoxon signed-rank tests** for non-parametric comparison of MPCI vs dMHL across all performance metrics
- **Significance annotation** with standard notation (* p < 0.05, ** p < 0.01, *** p < 0.001)

________________________________________
## Contact

naghme93@gmail.com
