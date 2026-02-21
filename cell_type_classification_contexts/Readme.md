## Overview

This repository contains scripts for calculating various methylation metrics from bisulfite sequencing PAT files, performing classification analyses using SVM, and generating visualizations for comparative evaluation of methylation metrics across different cell types and tissues.

________________________________________
## Requirements

- R (version 4.0 or higher)
- Required R packages: `tidyverse`, `caret`, `e1071`, `pROC`, `FactoMineR`, `factoextra`, `gridExtra`, `grid`, `parallel`, `ggpubr`, `broom`, `data.table`

________________________________________
## Repository Structure

### Metric Generation Scripts

#### 1. **random_regions_all_metric_generation_CD4_CD8.R**
Generates comprehensive methylation metrics for CD4+ and CD8+ T-cells across random genomic regions.

**Metrics calculated:**
- MPCI (Methylation Pattern Consistency Index)
- MHL (Methylation Haplotype Load)
- uMHL (unmethylated MHL)
- dMHL (differential MHL)
- Beta value (methylation frequency)
- ME (Methylation Entropy)
- Epipoly (Epipolymorphism)
- PDR (Pattern Discordance Rate)
- FDRP (Fragment Discordance Rate Proportion)
- qFDRP (quantitative FDRP)

#### 2. **Guo_chr22_all_metric_generation_CD4_CD8.R**
Similar to script #1 but specifically designed for Guo et al. MHB (Methylation Haplotype Block) regions on chromosome 22.

#### 3. **multi_tissue_random_regions_metric_generation.R**
Generates methylation metrics across multiple tissues using random genomic regions.

#### 4. **metric_generation_Guo_MHBs.R**
Focused on calculating MPCI, MHL, uMHL, dMHL, and Beta values for Guo et al. MHB regions.

#### 5. **Guo_MHBs_to_cpg_id_efficient_processing.R**
Utility script for converting genomic coordinates to CpG IDs for Guo et al. MHB regions using data.table for efficient processing.

________________________________________
### Classification and Analysis Scripts

#### 6. **CD4_CD8_all_metric_full_pipeline_pca_nested_cv_svm_linear_random_regions.R**
Complete analysis pipeline including PCA visualization and SVM classification with nested cross-validation.

**Key features:**
- PCA analysis for dimensionality reduction
- Linear SVM with nested cross-validation (10 repeats)
- Comprehensive performance metrics (AUC, Accuracy, Sensitivity, Specificity)
- Statistical comparisons between metrics
- Publication-ready visualizations

#### 7. **CD4_CD8_linear_svm_nested_cv_classification_wilcox_test.R**
Enhanced classification with region filtering and non-parametric statistical testing.

**Key features:**
- Removal of overlapping regions based on specific IDs
- Wilcoxon rank-sum tests for statistical comparisons
- Focus on MPCI vs dMHL comparisons
- Standard error visualization

#### 8. **Guo_MPCI_MHL_CD4_CD8_classifier_wilcoxon.R**
Focused analysis comparing only MPCI and MHL metrics for Guo et al. regions.

**Key features:**
- Two-metric comparative analysis
- Wilcoxon non-parametric testing
- Detailed performance visualization

________________________________________
## Core Functions

### Binary Matrix Processing
- **read_pat_data()**: Reads and processes PAT files from bisulfite sequencing
- **subset_pat_in_window()**: Extracts reads overlapping specified genomic windows
- **truncated_to_binary()**: Converts methylation patterns to binary matrices
- **handle_excess_NAs_in_binary_DMR()**: Cleans binary matrices by removing all-NA rows/columns

### Methylation Metrics

#### 1. **give_sign_for_weights(x_1, x_2)**
- **Purpose**: Assigns weights based on average methylation status of two rows
- **Input**: Two binary vectors (methylated = 1, unmethylated = 0)
- **Output**: 1 (>50% methylated), -1 (<50% methylated), or random (±1 if exactly 50%)

#### 2. **signed_manhattan_sim(binary_dmr)**
- **Purpose**: Calculates signed Manhattan similarity for a binary DMR matrix
- **Input**: Binary matrix (rows = reads, columns = CpG sites)
- **Output**: Weighted average of pairwise similarities between reads

#### 3. **MPCI(binary_dmr)**
- **Purpose**: Calculates Methylation Pattern Consistency Index
- **Input**: Binary matrix (rows = reads, columns = CpG sites)
- **Output**: Value from -1 to +1
  - Positive values: consistent methylation
  - Negative values: consistent unmethylation
  - Near zero: random methylation patterns

#### 4. **MHL(binary_dmr, max_hap_l = 9)**
- **Purpose**: Calculates Methylation Haplotype Load
- **Input**: Binary matrix and maximum haplotype length
- **Output**: Value from 0 to 1 indicating methylation block consistency

#### 5. **PDR(mat, min_cpgs_in_read = 4)**
- **Purpose**: Calculates Pattern Discordance Rate
- **Input**: Binary matrix, minimum CpGs per read
- **Output**: Proportion of discordant reads

#### 6. **epipoly_entropy(mat, k = 4)**
- **Purpose**: Calculates epipolymorphism and entropy
- **Input**: Binary matrix, window size k
- **Output**: List with epipoly and entropy values

#### 7. **FDRP_qFDRP(mat)**
- **Purpose**: Calculates Fragment Discordance metrics
- **Input**: Binary matrix
- **Output**: List with FDRP and qFDRP values

________________________________________
## Output Files

### Metric Files (.csv)
- Each metric saved as separate CSV file
- Format: First column = region IDs, subsequent columns = samples
- Values rounded to 3 decimal places

### Classification Results
- `*_classification_results.csv`: Raw classification results from nested CV
- `*_statistical_comparisons.csv`: Statistical test results between metrics
- `*_mpci_vs_*_comparisons.csv`: Focused comparisons

### Visualization Outputs
- `*_combined_pca_plots.png`: PCA visualizations for all metrics
- `*_sorted_auc_comparison_se.png`: AUC comparison bar plots with standard errors
- `*_detailed_comparison.png`: Detailed metric comparisons with significance

________________________________________
## Contact

naghme93@gmail.com
