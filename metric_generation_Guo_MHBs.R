library(data.table)

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


## metric functions
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

beta_value <- function(binary_dmr){
  return(mean(binary_dmr, na.rm = T))
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

### read MHBs
MHBs <- read.csv("../Guo_MHBs/MHBs_with_cpg_id.csv")
MHBs$ID <- paste0(MHBs$chrom,"-",MHBs$start_cpg_id,"-",MHBs$end_cpg_id)
markers <- MHBs[,c("chrom", "ID", "start_cpg_id", "end_cpg_id")]
names(markers) <- c("CHR","ID","START","END")
markers <- markers[markers$CHR == "chr22",]

### read pat files
pat_files_path <- "/media/mehrmohammadi_hdd/hdd2/nazer/home/data/loyfer_wgbs/pat file/CD4_CD8/chr22/"
CD4_CD8_pat_files <- paste0(pat_files_path,list.files(pat_files_path))

## generate metrics
## define the metric dfs
MPCI_df <- data.frame(matrix(, nrow = nrow(markers), ncol = length(CD4_CD8_pat_files)+1))
names(MPCI_df) <- c("ID", gsub(".*Cap(\\d+).*", "\\1", CD4_CD8_pat_files))
MPCI_df$ID <- markers$ID
MHL_df <- data.frame(matrix(, nrow = nrow(markers), ncol = length(CD4_CD8_pat_files)+1))
names(MHL_df) <- c("ID", gsub(".*Cap(\\d+).*", "\\1", CD4_CD8_pat_files))
MHL_df$ID <- markers$ID
uMHL_df <- data.frame(matrix(, nrow = nrow(markers), ncol = length(CD4_CD8_pat_files)+1))
names(uMHL_df) <- c("ID", gsub(".*Cap(\\d+).*", "\\1", CD4_CD8_pat_files))
uMHL_df$ID <- markers$ID
dMHL_df <- data.frame(matrix(, nrow = nrow(markers), ncol = length(CD4_CD8_pat_files)+1))
names(dMHL_df) <- c("ID", gsub(".*Cap(\\d+).*", "\\1", CD4_CD8_pat_files))
dMHL_df$ID <- markers$ID
Beta_df <- data.frame(matrix(, nrow = nrow(markers), ncol = length(CD4_CD8_pat_files)+1))
names(Beta_df) <- c("ID", gsub(".*Cap(\\d+).*", "\\1", CD4_CD8_pat_files))
Beta_df$ID <- markers$ID

## calculate
for (j in 1:length(CD4_CD8_pat_files)){
  fname <- CD4_CD8_pat_files[j]
  print(fname)
  pat_data <- read_pat_data(path = fname)
  for (i in 1:nrow(markers)){
    print(paste0("Analyzing marker #", i," : ",markers[i, "ID"]))
    start_cpg_id <- as.numeric(markers[i,"START"])
    end_cpg_id <- as.numeric(markers[i,"END"])
    binary_DMR <- subset_pat_in_window(start_cpg_id, end_cpg_id, pat_data)
    if(is.null(binary_DMR)){
      print("NULL binary DMR")
      MPCI_df[i,j+1] <- NA
      MHL_df[i,j+1] <- NA
      uMHL_df[i,j+1] <- NA
      dMHL_df[i,j+1] <- NA
      Beta_df[i,j+1] <- NA
    }
    else{
      MPCI_df[i,j+1] <- round(MPCI(binary_DMR),3)
      MHL_df[i,j+1] <- round(MHL(binary_DMR),3)
      uMHL_df[i,j+1] <- round(simple_uMHL(binary_DMR),3)
      dMHL_df[i,j+1] <- round((MHL_df[i,j+1] - uMHL_df[i,j+1]),3)
      Beta_df[i,j+1] <- round(beta_value(binary_DMR),3)
    }
    gc()
  }
}

write.csv(MPCI_df, "../Guo_MHBs/chr22/CD4_CD8/MPCI.csv")
write.csv(Beta_df, "../Guo_MHBs/chr22/CD4_CD8/Beta.csv")
write.csv(MHL_df, "../Guo_MHBs/chr22/CD4_CD8/MHL.csv")
write.csv(uMHL_df, "../Guo_MHBs/chr22/CD4_CD8/uMHL.csv")
write.csv(dMHL_df, "../Guo_MHBs/chr22/CD4_CD8/dMHL.csv")

