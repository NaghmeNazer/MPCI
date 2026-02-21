library(data.table)

read_pat_data <- function(path){
  pat.data <- fread(path)
  names(pat.data) <- c("chrom", "start_cpg_index", "methylation_pattern", "count")
  sort(table(pat.data$chrom))
  pat.data$end_cpg_index <- pat.data$start_cpg_index + nchar(pat.data$methylation_pattern)-1 #add end cpg_id to pat file
  return(pat.data)
}

find_high_coverage_regions <- function(pat_data, min_reads = 10, min_cpgs = 3) {
  setDT(pat_data)
  
  # 1. Calculate pattern length
  pat_data[, pattern_length := nchar(methylation_pattern)]
  
  # 2. Create table of all CpG positions with coverage
  all_positions <- list()
  
  for (i in 1:nrow(pat_data)) {
    row <- pat_data[i, ]
    if (row$pattern_length == 1) {
      all_positions[[length(all_positions) + 1]] <- data.table(
        cpg_pos = row$start_cpg_index,
        count = row$count
      )
    } else {
      positions <- seq(row$start_cpg_index, row$end_cpg_index, by = 1)
      all_positions[[length(all_positions) + 1]] <- data.table(
        cpg_pos = positions,
        count = row$count
      )
    }
  }
  
  # 3. Aggregate coverage per CpG
  coverage_per_cpg <- rbindlist(all_positions)[, .(coverage = sum(count)), by = cpg_pos]
  setorder(coverage_per_cpg, cpg_pos)
  
  # 4. Find consecutive CpGs with sufficient coverage
  regions <- list()
  current_start <- NULL
  current_cpgs <- c()
  
  for (j in 1:nrow(coverage_per_cpg)) {
    if (coverage_per_cpg$coverage[j] >= min_reads) {
      if (is.null(current_start)) {
        current_start <- coverage_per_cpg$cpg_pos[j]
        current_cpgs <- coverage_per_cpg$cpg_pos[j]
      } else if (coverage_per_cpg$cpg_pos[j] == current_cpgs[length(current_cpgs)] + 1) {
        current_cpgs <- c(current_cpgs, coverage_per_cpg$cpg_pos[j])
      } else {
        # Save region if long enough
        if (length(current_cpgs) >= min_cpgs) {
          mean_cov <- mean(coverage_per_cpg[cpg_pos %in% current_cpgs, coverage])
          regions[[length(regions) + 1]] <- data.table(
            start_cpg = current_start,
            end_cpg = current_cpgs[length(current_cpgs)],
            cpgs_in_region = length(current_cpgs),
            mean_coverage = mean_cov,
            region_score = length(current_cpgs) * mean_cov
          )
        }
        # Start new region
        current_start <- coverage_per_cpg$cpg_pos[j]
        current_cpgs <- coverage_per_cpg$cpg_pos[j]
      }
    } else if (!is.null(current_start)) {
      # Save region if long enough
      if (length(current_cpgs) >= min_cpgs) {
        mean_cov <- mean(coverage_per_cpg[cpg_pos %in% current_cpgs, coverage])
        regions[[length(regions) + 1]] <- data.table(
          start_cpg = current_start,
          end_cpg = current_cpgs[length(current_cpgs)],
          cpgs_in_region = length(current_cpgs),
          mean_coverage = mean_cov,
          region_score = length(current_cpgs) * mean_cov
        )
      }
      current_start <- NULL
      current_cpgs <- c()
    }
  }
  
  # Check last region
  if (!is.null(current_start) && length(current_cpgs) >= min_cpgs) {
    mean_cov <- mean(coverage_per_cpg[cpg_pos %in% current_cpgs, coverage])
    regions[[length(regions) + 1]] <- data.table(
      start_cpg = current_start,
      end_cpg = current_cpgs[length(current_cpgs)],
      cpgs_in_region = length(current_cpgs),
      mean_coverage = mean_cov,
      region_score = length(current_cpgs) * mean_cov
    )
  }
  
  if (length(regions) == 0) {
    message("No regions found")
    return(data.table())
  }
  
  # 5. Select non-overlapping regions (greedy by score)
  all_regions <- rbindlist(regions)
  all_regions <- all_regions[order(-region_score, -cpgs_in_region)]
  
  non_overlapping <- list()
  covered_positions <- integer(0)
  
  for (k in 1:nrow(all_regions)) {
    region <- all_regions[k, ]
    region_positions <- seq(region$start_cpg, region$end_cpg, by = 1)
    
    if (!any(region_positions %in% covered_positions)) {
      non_overlapping[[length(non_overlapping) + 1]] <- region
      covered_positions <- c(covered_positions, region_positions)
    }
  }
  
  result <- rbindlist(non_overlapping)
  result[, region_id := 1:.N]
  result[, chrom := "chr22"]
  
  # Reorder columns
  setcolorder(result, c("chrom", "region_id", "start_cpg", "end_cpg", 
                        "cpgs_in_region", "mean_coverage", "region_score"))
  setorder(result, start_cpg)
  
  return(result)
}

find_high_coverage_regions_fractional <- function(pat_data,
                                                  min_reads = 10,
                                                  min_cpgs = 3) {
  setDT(pat_data)
  
  # 1. Calculate pattern length
  pat_data[, pattern_length := nchar(methylation_pattern)]
  
  # 2. Expand patterns into CpGs with FRACTIONAL coverage
  all_positions <- vector("list", nrow(pat_data))
  
  for (i in seq_len(nrow(pat_data))) {
    row <- pat_data[i]
    
    frac_count <- row$count / row$pattern_length
    
    positions <- seq(row$start_cpg_index, row$end_cpg_index, by = 1)
    
    all_positions[[i]] <- data.table(
      cpg_pos = positions,
      coverage = frac_count
    )
  }
  
  # 3. Aggregate fractional coverage per CpG
  coverage_per_cpg <- rbindlist(all_positions)[
    , .(coverage = sum(coverage)), by = cpg_pos
  ]
  
  setorder(coverage_per_cpg, cpg_pos)
  
  # 4. Find consecutive CpGs with sufficient coverage
  regions <- list()
  current_start <- NULL
  current_cpgs <- integer(0)
  
  for (j in seq_len(nrow(coverage_per_cpg))) {
    
    if (coverage_per_cpg$coverage[j] >= min_reads) {
      
      if (is.null(current_start)) {
        current_start <- coverage_per_cpg$cpg_pos[j]
        current_cpgs <- coverage_per_cpg$cpg_pos[j]
        
      } else if (coverage_per_cpg$cpg_pos[j] ==
                 current_cpgs[length(current_cpgs)] + 1) {
        current_cpgs <- c(current_cpgs, coverage_per_cpg$cpg_pos[j])
        
      } else {
        # finalize region
        if (length(current_cpgs) >= min_cpgs) {
          mean_cov <- mean(
            coverage_per_cpg[cpg_pos %in% current_cpgs, coverage]
          )
          
          regions[[length(regions) + 1]] <- data.table(
            start_cpg = current_start,
            end_cpg = current_cpgs[length(current_cpgs)],
            cpgs_in_region = length(current_cpgs),
            mean_coverage = mean_cov,
            region_score = length(current_cpgs) * mean_cov
          )
        }
        
        # start new region
        current_start <- coverage_per_cpg$cpg_pos[j]
        current_cpgs <- coverage_per_cpg$cpg_pos[j]
      }
      
    } else if (!is.null(current_start)) {
      
      # finalize region
      if (length(current_cpgs) >= min_cpgs) {
        mean_cov <- mean(
          coverage_per_cpg[cpg_pos %in% current_cpgs, coverage]
        )
        
        regions[[length(regions) + 1]] <- data.table(
          start_cpg = current_start,
          end_cpg = current_cpgs[length(current_cpgs)],
          cpgs_in_region = length(current_cpgs),
          mean_coverage = mean_cov,
          region_score = length(current_cpgs) * mean_cov
        )
      }
      
      current_start <- NULL
      current_cpgs <- integer(0)
    }
  }
  
  # check last region
  if (!is.null(current_start) && length(current_cpgs) >= min_cpgs) {
    mean_cov <- mean(
      coverage_per_cpg[cpg_pos %in% current_cpgs, coverage]
    )
    
    regions[[length(regions) + 1]] <- data.table(
      start_cpg = current_start,
      end_cpg = current_cpgs[length(current_cpgs)],
      cpgs_in_region = length(current_cpgs),
      mean_coverage = mean_cov,
      region_score = length(current_cpgs) * mean_cov
    )
  }
  
  if (length(regions) == 0) {
    message("No regions found")
    return(data.table())
  }
  
  # 5. Greedy selection of non-overlapping regions
  all_regions <- rbindlist(regions)
  all_regions <- all_regions[order(-region_score, -cpgs_in_region)]
  
  non_overlapping <- list()
  covered_positions <- integer(0)
  
  for (k in seq_len(nrow(all_regions))) {
    region <- all_regions[k]
    region_positions <- seq(region$start_cpg, region$end_cpg, by = 1)
    
    if (!any(region_positions %in% covered_positions)) {
      non_overlapping[[length(non_overlapping) + 1]] <- region
      covered_positions <- c(covered_positions, region_positions)
    }
  }
  
  result <- rbindlist(non_overlapping)
  result[, region_id := seq_len(.N)]
  result[, chrom := unique(pat_data$chrom)]
  
  setcolorder(result, c("chrom", "region_id", "start_cpg", "end_cpg",
                        "cpgs_in_region", "mean_coverage", "region_score"))
  setorder(result, start_cpg)
  
  return(result)
}

### read pat files
pat_files_path <- "/media/mehrmohammadi_hdd/hdd1/nazer/download/GSE262275/pat_files/POD_zero/"
pat_files <- paste0(pat_files_path, list.files(pat_files_path))
sample_names <- list.files(pat_files_path)
# Load one pat data 
j=1
fname <- pat_files[j]
pat_data <- read_pat_data(path = fname)
pat_data.chr22 <- pat_data[pat_data$chrom == "chr22",]
high_cov_regions <- find_high_coverage_regions(pat_data.chr22, min_reads = 50, min_cpgs = 5)
high_cov_regions_fractional <- find_high_coverage_regions_fractional(pat_data.chr22, min_reads = 30, min_cpgs = 4)
write.csv(high_cov_regions, "../markers/GSE262275_high_cov_regions_minread50_mincpg_5.csv")
write.csv(high_cov_regions_fractional, "../markers/GSE262275_high_cov_regions_fractional_minread30_mincpg_4.csv")
