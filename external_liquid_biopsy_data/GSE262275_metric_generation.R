library(data.table)
library(tidyverse)

## Functions to process pat files and transform them into binary matrix
read_pat_data <- function(path){
  pat.data <- fread(path)
  names(pat.data) <- c("chrom", "start_cpg_index", "methylation_pattern", "count")
  sort(table(pat.data$chrom))
  pat.data$end_cpg_index <- pat.data$start_cpg_index + nchar(pat.data$methylation_pattern)-1 #add end cpg_id to pat file
  return(pat.data)
}

subset_pat_in_window <- function(window.start.cpg.id, window.end.cpg.id, pat.data){
  ### select reads in window
  reads.in.window <- pat.data[(pat.data$start_cpg_index <= window.end.cpg.id & 
                                 pat.data$start_cpg_index >= window.start.cpg.id) |
                                (pat.data$end_cpg_index <= window.end.cpg.id &
                                   pat.data$end_cpg_index >= window.start.cpg.id),]
  if (nrow(reads.in.window) == 0) return(NULL)
  else {
    reads.in.window.binary <- truncated_to_binary(reads.in.window, window.start.cpg.id, window.end.cpg.id)
    return(reads.in.window.binary)
  }
}

truncate_read_df <- function(meth_pat, start_ind, end_ind, start_window, end_window){
  start_of <- start_ind - start_window
  end_of <- end_ind - end_window
  meth_pat <- strsplit(meth_pat, split = "")[[1]]
  l <- length(meth_pat)
  meth_pat <- meth_pat[(1-min(start_of,0)):(l-max(end_of,0))]
  start_ind <- max(start_ind, start_window)
  end_ind <- min(end_ind, end_window)
  trunc_read <- data.frame(matrix(NA, nrow = 1, ncol = end_window-start_window+1))
  names(trunc_read) <- start_window:end_window
  trunc_read[,as.character(start_ind:end_ind)] <- as.numeric(c("C" = "1", "T" = "0")[meth_pat])
  return(trunc_read)
}

truncated_to_binary <- function(reads, window_start_cpg, window_end_cpg){
  binary <- as.numeric(unlist(mapply(truncate_read_df,reads$methylation_pattern,reads$start_cpg_index,reads$end_cpg_index,window_start_cpg,window_end_cpg)))
  binary <- matrix(binary, nrow = nrow(reads), byrow = T)
  return(binary)
}

handle_excess_NAs_in_binary_DMR <- function(binary_matrix){
  # Remove rows with all NAs - handle single row case
  if (nrow(binary_matrix) == 1) {
    if (all(is.na(binary_matrix[1, ]))) {
      return(NULL)
    }
  } else {
    binary_matrix <- binary_matrix[!apply(binary_matrix, 1, function(x) all(is.na(x))), , drop = FALSE]
  }
  if (nrow(binary_matrix) == 0) {
    return(NULL)
  }
  
  # Remove leading columns with all NAs
  start_col <- 1
  while (start_col <= ncol(binary_matrix)) {
    if (all(is.na(binary_matrix[, start_col]))) {
      start_col <- start_col + 1
    } else {
      break
    }
  }
  
  # Remove trailing columns with all NAs
  end_col <- ncol(binary_matrix)
  while (end_col >= 1) {
    if (all(is.na(binary_matrix[, end_col]))) {
      end_col <- end_col - 1
    } else {
      break
    }
  }
  
  # Check if we have any valid columns left
  if (start_col > end_col) {
    return(NULL)
  }
  
  # Subset the matrix to remove leading/trailing all-NA columns
  binary_matrix <- binary_matrix[, start_col:end_col, drop = FALSE]
  return(binary_matrix)
}

# ---------- MPCI and MHL/uMHL/dMHL -------------
give_sign_for_weights <- function(x_1, x_2){
  thresh <- mean(c(unlist(x_1),unlist(x_2)), na.rm = T)
  if (is.nan(thresh)) return(NA)
  if (thresh > 0.5) return(1)
  if (thresh < 0.5) return(-1)
  if (thresh == 0.5) return(sample(c(1,-1), 1))
}

signed_manhattan_sim <- function(binary_dmr){
  pairwise.sim <- 1-(dist(binary_dmr, method = "manhattan")/ncol(binary_dmr))
  if(length(pairwise.sim) == 0){
    metric <- NA
  } else if(sum(is.na(pairwise.sim)) == length(pairwise.sim)) {
    metric <- NA
  } else{
    row_combinations <- combn(1:nrow(binary_dmr), 2)  # Generate all combinations of row indices
    weights <- sapply(1:ncol(row_combinations), function(i) {
      row_indices <- row_combinations[, i]
      row1 <- binary_dmr[row_indices[1], ]
      row2 <- binary_dmr[row_indices[2], ]
      give_sign_for_weights(row1, row2)
    })
    metric <- sum(pairwise.sim * weights, na.rm = T)/sum(abs(weights[!is.na(pairwise.sim)]))
  }
  return(metric)
}

MPCI <- function(binary_dmr){
  if (is.na(signed_manhattan_sim(binary_dmr)) | is.na(signed_manhattan_sim(t(binary_dmr)))){
    metric <- NA
  } else {
    metric <- mean(signed_manhattan_sim(binary_dmr), signed_manhattan_sim(t(binary_dmr)))
  }
  return(metric)
}

count_non_na <- function(row) {
  sum(!is.na(row))
}

# Function to find the length of the longest consecutive 1s in a vector
longest_consecutive_ones <- function(trunc_pat) {
  rle_result <- rle(trunc_pat)
  if (any(rle_result$values == 1, na.rm = T)) {
    run_lengths <- rle_result$lengths[rle_result$values == 1]
    max_run_length <- max(run_lengths, na.rm = T)
  } else {
    max_run_length <- 0
  }
  return(max_run_length)
}

# Function to find the length of the longest consecutive 0s in a vector
longest_consecutive_zeros <- function(trunc_pat) {
  rle_result <- rle(trunc_pat)
  if (any(rle_result$values == 0, na.rm = T)) {
    run_lengths <- rle_result$lengths[rle_result$values == 0]
    max_run_length <- max(run_lengths, na.rm = T)
  } else {
    max_run_length <- 0
  }
  return(max_run_length)
}

find_haplotypes_with_desired_length <- function(trunc_pat, l) {
  non_na_indices <- which(!is.na(trunc_pat))
  non_na_starts <- non_na_indices[non_na_indices <= length(trunc_pat) - l + 1]
  
  subsets <- lapply(non_na_starts, function(start_idx) {
    end_idx <- start_idx + l - 1
    if (any(is.na(trunc_pat[start_idx:end_idx]))) {
      return(NULL)
    } else {
      return(trunc_pat[start_idx:end_idx])
    }
  })
  
  subsets <- Filter(function(x) !is.null(x), subsets)
  return(subsets)
}

MHL <- function(binary_dmr, max_hap_l = 9){
  read_lengths <- apply(binary_dmr, 1, count_non_na)
  max_read_length <- max(read_lengths)
  longest_ful_meth_haps <- apply(binary_dmr, 1, longest_consecutive_ones)
  max_ful_meth_length <- max(longest_ful_meth_haps)
  numerator <- 0
  for (i in 1:min(max_ful_meth_length, max_hap_l)){
    haps <- unlist(apply(binary_dmr,1, find_haplotypes_with_desired_length, l = i), recursive = F)
    if (length(haps) != 0) full_meth_fraction <- sum(sapply(haps, function(subset) all(subset == 1)))/length(haps)
    if (length(haps) == 0) full_meth_fraction <- 0
    #print(paste0("i = ", i, " Full_meth_Fraction = ", round(full_meth_fraction,2), " i*full_meth_fraction = ", round(i*full_meth_fraction,2)))
    numerator <- numerator + i*full_meth_fraction
  }
  MHL <- numerator/sum(seq(1,min(max_read_length, max_hap_l)))
  return(MHL)
}

uMHL <- function(binary_dmr, max_hap_l = 9){
  read_lengths <- apply(binary_dmr, 1, count_non_na)
  max_read_length <- max(read_lengths)
  longest_ful_umeth_haps <- apply(binary_dmr, 1, longest_consecutive_zeros)
  max_ful_umeth_length <- max(longest_ful_umeth_haps)
  numerator <- 0
  for (i in 1:min(max_ful_umeth_length, max_hap_l)){
    haps <- unlist(apply(binary_dmr,1, find_haplotypes_with_desired_length, l = i), recursive = F)
    if (length(haps) != 0) full_umeth_fraction <- sum(sapply(haps, function(subset) all(subset == 0)))/length(haps)
    if (length(haps) == 0) full_umeth_fraction <- 0
    #print(paste0("i = ", i, " Full_meth_Fraction = ", round(full_meth_fraction,2), " i*full_meth_fraction = ", round(i*full_meth_fraction,2)))
    numerator <- numerator + i*full_umeth_fraction
  }
  uMHL <- numerator/sum(seq(1,min(max_read_length, max_hap_l)))
  return(uMHL)
}

simple_uMHL <- function(binary_dmr){
  contrast_dmr <- 1-binary_dmr
  uMHL <- MHL(contrast_dmr)
  return(uMHL)
}

dMHL <- function(MHL, uMHL){
  return(MHL-uMHL)
}

# --------- process regions over samples ----------
set.seed(123)
### read regions
#markers <- read.csv("../random_regions/400_random_regions.csv")
#marker_num <- 400
# markers <- markers %>%
#   arrange(desc(mean_coverage)) %>%  # arrange in descending order
#   slice(1:marker_num)                      # return rows 1 through 400

markers <- read.csv("../markers/GSE262275_high_cov_regions_fractional_minread30_mincpg_4.csv")
markers$ID <- paste0(markers$chrom,"-", markers$start_cpg,"-", markers$end_cpg)
markers <- markers[,c("chrom", "ID", "start_cpg", "end_cpg")]
names(markers) <- c("CHR", "ID", "START", "END")
#markers <- markers[sample.int(nrow(markers), marker_num),]
### read pat files
pat_files_path <- "/media/mehrmohammadi_hdd/hdd1/nazer/download/GSE262275/pat_files/POD_zero/"
pat_files <- paste0(pat_files_path, list.files(pat_files_path))

## Improved: Pre-allocate all metrics in a list for efficient storage
sample_names <- list.files(pat_files_path)
n_regions <- nrow(markers)
n_samples <- length(pat_files)

# Create a list to store all metrics
metrics_list <- list(
  MPCI = matrix(NA_real_, nrow = n_regions, ncol = n_samples),
  MHL = matrix(NA_real_, nrow = n_regions, ncol = n_samples),
  uMHL = matrix(NA_real_, nrow = n_regions, ncol = n_samples),
  dMHL = matrix(NA_real_, nrow = n_regions, ncol = n_samples),
  coverage = matrix(NA_real_, nrow = n_regions, ncol = n_samples)
)

for (j in seq_along(pat_files)) {
  fname <- pat_files[j]
  print(paste("Processing sample:", j, "/", n_samples, "-", basename(fname)))
  
  # Load pat data once per sample
  pat_data <- read_pat_data(path = fname)
  
  # Improved: Use vectorized operations where possible
  for (i in seq_len(n_regions)) {
    if (i %% 50 == 0) print(paste("  Region:", i, "/", n_regions))
    
    start_cpg_id <- as.numeric(markers[i, "START"])
    end_cpg_id <- as.numeric(markers[i, "END"])
    
    binary_DMR <- subset_pat_in_window(start_cpg_id, end_cpg_id, pat_data)
    if(!is.null(binary_DMR)){
      binary_DMR <- handle_excess_NAs_in_binary_DMR(binary_DMR)
    }
    if (is.null(binary_DMR)) {
      # All metrics remain NA (pre-allocated)
      next
    }

    # Improved: Calculate all metrics in batch to avoid repeated function calls
    metrics_list$MPCI[i, j] <- round(MPCI(binary_DMR), 3)
    metrics_list$MHL[i, j] <- round(MHL(binary_DMR), 3)
    metrics_list$uMHL[i, j] <- round(simple_uMHL(binary_DMR), 3)
    metrics_list$dMHL[i, j] <- round(metrics_list$MHL[i, j] - metrics_list$uMHL[i, j], 3)
    metrics_list$coverage[i, j] <- nrow(binary_DMR)
    # Clean up
    rm(binary_DMR, tmp_epi, tmp_fdrp)
  }
  
  # Force garbage collection after each sample
  gc()
}

## Improved: Convert to data frames and save all at once
output_dir <- "../GSE262275/high_cov_regions/all604_high_fractional_covered/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create data frames with IDs
result_dfs <- list()
for (metric_name in names(metrics_list)) {
  df <- as.data.frame(metrics_list[[metric_name]])
  colnames(df) <- sample_names
  df$ID <- markers$ID
  # Reorder columns to have ID first
  df <- df[, c("ID", sample_names)]
  result_dfs[[metric_name]] <- df
}

## Improved: Save all files in one operation
output_files <- paste0(output_dir, names(result_dfs), ".csv")
for (i in seq_along(result_dfs)) {
  write.csv(result_dfs[[i]], output_files[i], row.names = FALSE)
}

print("All metrics saved successfully!")
