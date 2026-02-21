library(data.table)
MHBs <- read.csv("../Guo_MHBs/MHL_MHBs.csv")
head(MHBs)
MHBs$start.hg19. <- MHBs$start.hg19. + 1
MHBs$end.hg19. <- MHBs$end.hg19. + 1

hg19.cpg.annot <- fread("/media/mehrmohammadi_hdd/hdd2/nazer/home/data/methylation_blocks/hg19CpG.bed.gz", data.table = F)
head(hg19.cpg.annot)
rownames(hg19.cpg.annot) <- hg19.cpg.annot$V3
names(hg19.cpg.annot) <- c("chrom", "pos", "cpg_id")

# Convert to data.tables
setDT(hg19.cpg.annot)
setDT(MHBs)

# Exact match for start positions
start_match <- hg19.cpg.annot[MHBs, .(chrom, start.hg19., end.hg19., cpg_id), 
                              on = .(chrom, pos = start.hg19.)]
setnames(start_match, "cpg_id", "start_cpg_id")

# Exact match for end positions
end_match <- hg19.cpg.annot[MHBs, .(chrom, start.hg19., end.hg19., cpg_id), 
                            on = .(chrom, pos = end.hg19.)]
setnames(end_match, "cpg_id", "end_cpg_id")

# Merge the results
MHBs_result <- start_match[end_match, on = .(chrom, start.hg19., end.hg19.)]

# Remove duplicate columns and keep the final result
MHBs_result <- MHBs_result[, .(chrom, start.hg19., end.hg19., start_cpg_id, end_cpg_id)]

# View the result
print(MHBs_result)
MHBs_result$start.hg19. <- MHBs_result$start.hg19. - 1
MHBs_result$end.hg19. <- MHBs_result$end.hg19. - 1
write.csv(MHBs_result, "../Guo_MHBs/MHBs_with_cpg_id.csv")
a<- read.csv("../Guo_MHBs/MHBs_with_cpg_id.csv")
