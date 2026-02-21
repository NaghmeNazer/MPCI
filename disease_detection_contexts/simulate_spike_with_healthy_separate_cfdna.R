## code to spike reads in a region
library(data.table)

# function to spike reads from spike_region to the raget region
spike <- function(target_region, spike_region, spike_ratio = 0.1, replace = T){
  target_region.coverage <- nrow(target_region)
  spike_region.coverage <- nrow(spike_region)
  if (!replace) spike_ratio <- spike_ratio/(1-spike_ratio)
  spike_read_number <- ceiling(spike_ratio*target_region.coverage)
  if (spike_read_number <= spike_region.coverage){
    spike_reads.index <- sample(1:spike_region.coverage, size = spike_read_number)
    spike_reads <- spike_region[spike_reads.index,]
    if (replace){
      target_reads.index <- sample(1:target_region.coverage, size = spike_read_number)
      target_region <- target_region[-c(target_reads.index),]
    }
    result_region <- rbind(spike_reads, target_region)
    #return(result_region)
  }
  if (spike_read_number > spike_region.coverage) {
    print("The region from which you want to select spike reads has insufficient coverage, Spiking is performed with replacement. try to decrease the spike_ratio if needed.")
    spike_reads.index <- sample(1:spike_region.coverage, size = spike_read_number, replace = T)
    spike_reads <- spike_region[spike_reads.index,]
    if (replace){
      target_reads.index <- sample(1:target_region.coverage, size = spike_read_number)
      target_region <- target_region[-c(target_reads.index),]
    }
    result_region <- rbind(spike_reads, target_region)
    #return(result_region)
  }
  return(result_region)
}

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
      return(NA)
    }
  } else {
    binary_matrix <- binary_matrix[!apply(binary_matrix, 1, function(x) all(is.na(x))), , drop = FALSE]
  }
  if (nrow(binary_matrix) == 0) {
    return(NA)
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
    return(NA)
  }
  
  # Subset the matrix to remove leading/trailing all-NA columns
  binary_matrix <- binary_matrix[, start_col:end_col, drop = FALSE]
  return(binary_matrix)
}
# ---------- Metric Functions ----------
count_non_na <- function(row) sum(!is.na(row))

# returns TRUE if a vector has at least one 0 and one 1 (discordant)
is_discordant <- function(vec) {
  vec2 <- vec[!is.na(vec)]
  if (length(vec2) == 0) return(NA)
  return(any(vec2 == 0) && any(vec2 == 1))
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

simple_uMHL <- function(binary_dmr){
  contrast_dmr <- 1-binary_dmr
  uMHL <- MHL(contrast_dmr)
  return(uMHL)
}


##################
##
print("Reading Regions")
random.regions <- fread("../random_regions/400_random_regions.csv")

# Directly specify the overlapping region IDs to remove
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
random.regions <- random.regions[!(random.regions$ID %in% regions_to_remove),]
path <- "/media/mehrmohammadi_hdd/hdd2/nazer/home/data/loyfer_wgbs/pat file/chr22/"
all_files <- list.files(path)


metric_df <- data.frame(matrix(nrow=nrow(random.regions), ncol = 13))
names(metric_df) <- c("DMR_ID", 
                      "initial_coverage_cfdna", "coverage_neuron", "coverage_simulated",
                      "MPCI_cfdna", "MPCI_neuron", "MPCI_simulated",
                      "MHL_cfdna", "MHL_neuron", "MHL_simulated",
                      "uMHL_cfdna", "uMHL_neuron", "uMHL_simulated")

spike_data.tissue <- "Neuron"
target_data.tissue <- "cfDNA"
sample_num <- 100
target_data.desired_coverage <- 100
round_digits <- 2
spike_ratios <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
repeatition_number = 5
for (spike_ratio in spike_ratios){
  print(paste0("Spike Ratio = ", spike_ratio))
  results_dir <- paste0("../simulate_spike_separate_healthy/12-24-2025_spike_generation_results_1percent/",repeatition_number, "/",spike_ratio,"_",spike_data.tissue,"_in_",target_data.tissue)
  dir.create(results_dir, recursive = T)
  for (i in 1:sample_num){
    print(paste0("Calculating Metrics on ", nrow(random.regions)," DMRs in sample #", i))
    print(paste0("Read and Randomly Select Spike Tissue (", spike_data.tissue,") Pat File"))
    spike_data.fnames <- all_files[grepl(spike_data.tissue, all_files)]
    spike_data_i <- sample(1:length(spike_data.fnames),1)
    spike_data.pat <- read_pat_data(paste0(path,spike_data.fnames[spike_data_i]))
    
    print(paste0("Read and Randomly Select Target Tissue (", target_data.tissue,") Pat File"))
    target_data.fnames <- all_files[grepl(target_data.tissue, all_files)]
    target_data_i <- sample(1:length(target_data.fnames),1)
    target_data.pat <- read_pat_data(paste0(path,target_data.fnames[target_data_i]))
    healthy_data_i <- sample(1:length(target_data.fnames),1)
    healthy_data.pat <- read_pat_data(paste0(path,target_data.fnames[healthy_data_i]))
    
    gc()
    for (sig_dmr.i in 1:nrow(random.regions)){
      print(paste0("Spike Ratio = ", spike_ratio," <Sample #", i, "> Analysing DMR #", sig_dmr.i, " = ", random.regions[sig_dmr.i, "ID"]))
      spike_data.binary_DMR <- subset_pat_in_window(random.regions$START[sig_dmr.i], random.regions$END[sig_dmr.i], spike_data.pat)
      target_data.binary_DMR <- subset_pat_in_window(random.regions$START[sig_dmr.i], random.regions$END[sig_dmr.i], target_data.pat)
      healthy_data.binary_DMR <- subset_pat_in_window(random.regions$START[sig_dmr.i], random.regions$END[sig_dmr.i], healthy_data.pat)
      if (is.null(spike_data.binary_DMR) | is.null(target_data.binary_DMR)){
        result_data.binary_DMR <- NULL
      }
      if (!(is.null(spike_data.binary_DMR) | is.null(target_data.binary_DMR) | is.null(healthy_data.binary_DMR))){
        healthy_data.initial_coverage <- nrow(healthy_data.binary_DMR)
        healthy_data.read_index <- sample(1:healthy_data.initial_coverage, size = target_data.desired_coverage, replace = T)
        healthy_data.binary_DMR <- healthy_data.binary_DMR[healthy_data.read_index,]
        target_data.initial_coverage <- nrow(target_data.binary_DMR)
        target_data.read_index <- sample(1:target_data.initial_coverage, size = target_data.desired_coverage, replace = T)
        target_data.binary_DMR <- target_data.binary_DMR[target_data.read_index,]
        result_data.binary_DMR <- spike(target_data.binary_DMR, spike_data.binary_DMR, spike_ratio = spike_ratio)
      }
      if(is.null(result_data.binary_DMR) | is.null(healthy_data.binary_DMR)){
        metric_df[sig_dmr.i,] <- NA
      }
      else{
        metric_array <- c(random.regions$ID[sig_dmr.i],
                          healthy_data.initial_coverage, nrow(spike_data.binary_DMR), nrow(result_data.binary_DMR),
                          round(MPCI(healthy_data.binary_DMR), round_digits), 
                          round(MPCI(spike_data.binary_DMR), round_digits),
                          round(MPCI(result_data.binary_DMR), round_digits),
                          round(MHL(healthy_data.binary_DMR), round_digits),
                          round(MHL(spike_data.binary_DMR), round_digits),
                          round(MHL(result_data.binary_DMR), round_digits),
                          round(simple_uMHL(healthy_data.binary_DMR), round_digits),
                          round(simple_uMHL(spike_data.binary_DMR), round_digits),
                          round(simple_uMHL(result_data.binary_DMR), round_digits))
        metric_df[sig_dmr.i,1:13] <- metric_array
      }
    }
    write.csv(metric_df, paste0(results_dir,"/",spike_data.tissue,"_in_",target_data.tissue,"_spikeratio", spike_ratio,"_S",i,".csv"))
  }
}
