# Spike Simulation and Classification Pipeline

This repository contains scripts for simulating realistic cfDNA mixtures by spiking methylation reads from tissue samples (e.g., Neuron) into healthy cfDNA background at controlled ratios. The pipeline then generates methylation metrics (MPCI, MHL, uMHL, dMHL) and performs SVM classification with nested cross-validation to assess the ability to detect the spiked-in signal across different concentrations.

________________________________________
## Requirements

- R (version 4.0 or higher)
- Required R packages: `data.table`, `caret`, `e1071`, `pROC`, `ggplot2`, `doParallel`, `readr`, `patchwork`

________________________________________
## Repository Structure

### Spike Simulation Script

#### 1. **simulate_spike_with_healthy_separate_cfdna.R**
Simulates spiked-in methylation signals by mixing reads from tissue samples into healthy cfDNA background.

**Core Spike Function:**

- **spike(target_region, spike_region, spike_ratio = 0.1, replace = TRUE)**: 
  - Purpose: Spikes reads from a source region into a target region at specified ratio
  - Input: 
    - `target_region`: Binary matrix of target region reads
    - `spike_region`: Binary matrix of source region reads
    - `spike_ratio`: Fraction of reads to replace/spike (default: 0.1)
    - `replace`: Whether to replace target reads with spike reads (TRUE) or add them (FALSE)
  - Output: Combined binary matrix with spiked-in reads

**Binary Matrix Processing Functions:**
- **read_pat_data(path)**: Reads and processes PAT files
- **subset_pat_in_window()**: Extracts reads overlapping specified genomic windows
- **truncated_to_binary()**: Converts methylation patterns to binary matrices
- **handle_excess_NAs_in_binary_DMR()**: Cleans binary matrices by removing all-NA rows/columns

**Methylation Metric Functions:**

- **give_sign_for_weights(x_1, x_2)**: Assigns weights based on average methylation status
  - Output: 1 (>50% methylated), -1 (<50% methylated), or random (±1 if exactly 50%)

- **signed_manhattan_sim(binary_dmr)**: Calculates signed Manhattan similarity
  - Output: Weighted average of pairwise similarities

- **MPCI(binary_dmr)**: Methylation Pattern Consistency Index
  - Output: Value from -1 to +1 indicating pattern consistency

- **MHL(binary_dmr, max_hap_l = 9)**: Methylation Haplotype Load
  - Output: Value from 0 to 1 indicating methylation block consistency

- **simple_uMHL(binary_dmr)**: Unmethylated haplotype load
  - Output: Value from 0 to 1 for unmethylated patterns

**Simulation Parameters:**
- `spike_data.tissue`: Source tissue for spike reads (e.g., "Neuron")
- `target_data.tissue`: Target background tissue (e.g., "cfDNA")
- `spike_ratios`: Vector of spike ratios to test (1%, 2%, 3%, 4%, 5%, 10%)
- `sample_num`: Number of simulation replicates (default: 100)
- `target_data.desired_coverage`: Target coverage depth (default: 100)
- `repeatition_number`: Replicate identifier for output organization

**Output Files:**
- CSV files per sample per spike ratio containing:
  - Region ID and coverage information
  - MPCI values for healthy, spike tissue, and simulated mixture
  - MHL values for all three conditions
  - uMHL values for all three conditions
  - dMHL calculated as MHL - uMHL

________________________________________
### Classification Script

#### 2. **spike_nested_cv_svm_classifier.R**
Performs SVM classification to evaluate detection of spiked-in signals across different concentrations.

**Key Functions:**

- **svm_nested_cv(data, kernel_type = "linear", outer_folds = 5, inner_folds = 3)**:
  - Purpose: Runs nested cross-validation with SVM
  - Input: Prepared data matrix, kernel type, fold parameters
  - Features:
    - Balanced fold creation using `createMultiFolds`
    - Probability calibration for ROC analysis
    - Automatic direction detection for ROC curves
    - Hyperparameter tuning for C (and sigma for RBF kernel)
  - Output: Data frame with per-fold results including:
    - AUC, Sensitivity, Specificity, Accuracy
    - Probability distributions
    - Optimal hyperparameters

- **load_and_process_data(spike_ratio, spike_data.tissue, target_data.tissue, min_cov.thresh = 10, sample_num = 100)**:
  - Purpose: Loads and processes simulation results for classification
  - Input: Spike ratio, tissue types, minimum coverage threshold
  - Processing steps:
    - Loads all simulation replicates for given spike ratio
    - Calculates dMHL from MHL and uMHL
    - Filters regions with sufficient coverage
    - Removes zero-variance features
    - Creates labeled dataset with Case (spiked) vs Healthy samples
  - Output: List containing data frames for MPCI and dMHL metrics

- **run_full_analysis()**:
  - Purpose: Orchestrates complete analysis across all spike ratios
  - Processes all spike ratios (1%, 2%, 3%, 4%, 5%, 10%)
  - Runs linear SVM for both MPCI and dMHL
  - Saves complete results and summary statistics
  - Output: List with full results and summary tables

**Visualization Functions:**

- **create_four_line_plots_smart(all_results)**:
  - Purpose: Creates publication-ready line plots with dynamic y-axis scaling
  - Features:
    - Four separate plots: AUC, Accuracy, Specificity, Sensitivity
    - Dynamic y-axis scaling based on minimum values (including error bars)
    - Standard error bars for all metrics
    - Color-coded by metric (dMHL: red, MPCI: blue)
    - 2x2 combined grid layout
  - Output: List containing all plots and y-axis parameters

**Output Files:**

*Data Files:*
- `svm_results_complete.csv`: Raw classification results for all folds
- `svm_results_summary.csv`: Summary statistics per spike ratio and metric

*Visualization Files:*
- `auc_plot_yaxis_X-X.png`: AUC values across spike ratios
- `accuracy_plot_yaxis_X-X.png`: Accuracy values across spike ratios
- `specificity_plot_yaxis_X-X.png`: Specificity values across spike ratios
- `sensitivity_plot_yaxis_X-X.png`: Sensitivity values across spike ratios
- `all_metrics_grid_yaxis_X-X.png`: Combined 2x2 grid of all plots

________________________________________
## Contact

naghme93@gmail.com
