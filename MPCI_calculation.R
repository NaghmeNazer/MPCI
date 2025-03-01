# Function to assign weights based on the average methylation status of two rows
give_sign_for_weights <- function(x_1, x_2){
  # Calculate the mean methylation status of the two rows
  thresh <- mean(c(unlist(x_1), unlist(x_2)), na.rm = T)
  
  # Handle cases where the mean is NaN (e.g., all values are NA)
  if (is.nan(thresh)) return(NA)
  
  # Assign weights based on the mean methylation status
  if (thresh > 0.5) return(1)       # More methylated CpGs
  if (thresh < 0.5) return(-1)      # More unmethylated CpGs
  if (thresh == 0.5) return(sample(c(1, -1), 1))  # Randomly assign +1 or -1 if equal
}

# Function to calculate the signed Manhattan similarity for a binary DMR matrix
signed_manhattan_sim <- function(binary_dmr){
  # Calculate pairwise Manhattan similarity between rows
  pairwise.sim <- 1 - (dist(binary_dmr, method = "manhattan") / ncol(binary_dmr))
  
  # Handle edge cases where pairwise similarity cannot be calculated
  if (length(pairwise.sim) == 0) {
    metric <- NA
  } else if (sum(is.na(pairwise.sim)) == length(pairwise.sim)) {
    metric <- NA
  } else {
    # Generate all pairwise combinations of row indices
    row_combinations <- combn(1:nrow(binary_dmr), 2)
    
    # Calculate weights for each pairwise combination
    weights <- sapply(1:ncol(row_combinations), function(i) {
      row_indices <- row_combinations[, i]
      row1 <- binary_dmr[row_indices[1], ]
      row2 <- binary_dmr[row_indices[2], ]
      give_sign_for_weights(row1, row2)  # Assign weights based on methylation status
    })
    
    # Calculate the weighted average of pairwise similarities
    metric <- sum(pairwise.sim * weights, na.rm = T) / sum(abs(weights[!is.na(pairwise.sim)]))
  }
  return(metric)
}

# Function to calculate the Methylation Pattern Consistency Index (MPCI)
MPCI <- function(binary_dmr){
  # Calculate signed Manhattan similarity for the original and transposed matrix
  sim_original <- signed_manhattan_sim(binary_dmr)
  sim_transposed <- signed_manhattan_sim(t(binary_dmr))
  
  # Handle cases where either similarity is NA
  if (is.na(sim_original) | is.na(sim_transposed)) {
    metric <- NA
  } else {
    # Calculate MPCI as the average of the two similarities
    metric <- mean(sim_original, sim_transposed)
  }
  return(metric)
}

# Test MPCI on synthetic binary DMR matrices
fully_methylated <- matrix(rep(1, 50), nrow = 10)  # Fully methylated DMR
MPCI(fully_methylated)

mixed_methylation_1 <- matrix(c(rep(0, 20), rep(1, 10), rep(0, 20)), nrow = 10, byrow = T)  # Mixed methylation
MPCI(mixed_methylation_1)

mostly_unmethylated <- matrix(rep(0, 50), nrow = 10)  # Mostly unmethylated DMR with some random methylation
mostly_unmethylated[6, 5] <- 1
mostly_unmethylated[10, 1] <- 1
mostly_unmethylated[2, 3] <- 1
mostly_unmethylated[9, 3] <- 1
mostly_unmethylated[9, 5] <- 1
mostly_unmethylated[5, 4] <- 1
mostly_unmethylated[4, 5] <- 1
mostly_unmethylated[4, 1] <- 1
mostly_unmethylated[6, 3] <- 1
mostly_unmethylated[1, 3] <- 1
MPCI(mostly_unmethylated)

fully_unmethylated <- matrix(rep(0, 50), nrow = 10)  # Fully unmethylated DMR
MPCI(fully_unmethylated)

mixed_methylation_2 <- matrix(c(rep(1, 20), rep(0, 10), rep(1, 20)), nrow = 10, byrow = T)  # Mixed methylation
MPCI(mixed_methylation_2)

mostly_methylated <- matrix(rep(1, 50), nrow = 10)  # Mostly methylated DMR with some random unmethylation
mostly_methylated[6, 5] <- 0
mostly_methylated[10, 1] <- 0
mostly_methylated[2, 3] <- 0
mostly_methylated[9, 3] <- 0
mostly_methylated[9, 5] <- 0
mostly_methylated[5, 4] <- 0
mostly_methylated[4, 5] <- 0
mostly_methylated[4, 1] <- 0
mostly_methylated[6, 3] <- 0
mostly_methylated[1, 3] <- 0
MPCI(mostly_methylated)

# Test MPCI on a real binary DMR matrix from a CSV file
real_binary_dmr <- read.csv("GSM5652222_Oligodendrocyte_chr6-9857061-9857066_binaryDMR.csv")
real_binary_dmr$X <- NULL  # Remove unnecessary column
MPCI(real_binary_dmr)
