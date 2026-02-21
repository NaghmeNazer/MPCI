# Methylation Pattern Analysis Suite

This suite provides a comprehensive collection of R scripts for analyzing DNA methylation patterns at the haplotype level. It includes functions for calculating various methylation metrics (MPCI, MHL, uMHL, dMHL, PDR, FDRP, Epipoly, Entropy), performing classification analyses with SVM, and evaluating detection sensitivity in cell type classification, liquid biopsy analysis, and simualted disease detection.


________________________________________
## Requirements

- R (version 4.0 or higher)
- Required R packages: `data.table`, `tidyverse`, `caret`, `e1071`, `pROC`, `FactoMineR`, `factoextra`, `ggplot2`, `doParallel`, `patchwork`, `gridExtra`, `ggpubr`, `broom`, `readr`, `tidyr`, `dplyr`

________________________________________
## Common Workflow

1. **Generate Binary Matrices** from PAT files using provided functions
2. **Calculate Methylation Metrics** (MPCI, MHL, uMHL, dMHL, etc.)
3. **Prepare Data** for classification by transposing and labeling
4. **Run Classification** with nested cross-validation SVM
5. **Visualize Results** with publication-ready plots
6. **Perform Statistical Tests** (Wilcoxon) for comparisons

________________________________________
## Repository Structure

## MPCI_calculation
- `MPCI_calculation.R` - Core functions for MPCI calculation

## cell_type_classification_contexts
- `random_regions_all_metric_generation_CD4_CD8.R` - Full metric generation for CD4/CD8 across random regions
- `Guo_chr22_all_metric_generation_CD4_CD8.R` - Metric generation for Guo MHB regions on chr22
- `multi_tissue_random_regions_metric_generation.R` - Metric generation across multiple tissues
- `metric_generation_Guo_MHBs.R` - Basic MPCI, MHL, uMHL, dMHL, Beta for Guo MHB regions
- `Guo_MHBs_to_cpg_id_efficient_processing.R` - Convert genomic coordinates to CpG IDs
- `CD4_CD8_all_metric_full_pipeline_pca_nested_cv_svm_linear_random_regions.R` - Full pipeline with PCA and nested CV SVM
- `CD4_CD8_linear_svm_nested_cv_classification_wilcox_test.R` - SVM classification with Wilcoxon tests and region filtering
- `Guo_MPCI_MHL_CD4_CD8_classifier_wilcoxon.R` - MPCI vs MHL comparison with Wilcoxon tests

## external_liquid_biopsy_data
- `GSE262275_high_cov_regions_finding.R` - Identify high-coverage regions from PAT files
- `GSE262275_metric_generation.R` - Generate methylation metrics for high-coverage regions
- `GSE262275_svm_nested_classifier_wilcoxon.R` - SVM classification comparing PRE vs POST treatment

## disease_detection_contexts
- `simulate_spike_with_healthy_separate_cfdna.R` - Simulate spiked-in methylation signals by mixing tissue reads into cfDNA background
- `spike_nested_cv_svm_classifier.R` - SVM classification on simulated data across multiple spike ratios

________________________________________
## Contact

naghme93@gmail.com

