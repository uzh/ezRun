

require(ezRun)
require(stringr)

sce = readRDS("/export/local/scratch/gtan/p2497-SCCountQC/20171222.A-SiCSeq_SCs_P5_SCCountQC/sce.rds")
colData(sce)$nGenesDetected = colSums(assays(sce)$counts > 3)
cd = colData(sce)

platePos <- str_extract(cd$`PlatePosition [Characteristic]`,
                        "_[[:alpha:]]\\d+$")
plateRow <- str_extract(platePos, "[[:alpha:]]")
plateCol <- str_extract(platePos, "\\d+")

scores <- ezMatrix(NA, rows=LETTERS[1:ifelse(max(plateRow) > "I", 16, 8)],
                        cols=seq_len(ifelse(max(as.integer(plateCol)) > 12, 24, 12)))
mapToPlate = function(plateRow, plateCol, values){
  scoreMatrix <- ezMatrix(NA, rows=LETTERS[1:ifelse(max(plateRow) > "I", 16, 8)],
                     cols=seq_len(ifelse(max(as.integer(plateCol)) > 12, 24, 12)))
  for(j in seq_len(length(plateRow))) { scoreMatrix[plateRow[j], plateCol[j]] <- values[j]}
  return(scoreMatrix)
}


## for each score there should be three plots, e.g. for PCT_RIBOSOMAL_BASES
heatmapPlate(mapToPlate(plateRow, plateCol, cd$PCT_RIBOSOMAL_BASES), log10 = TRUE, center = TRUE)
plot(cd$PF_READS, cd$nGenesDetected, log="xy") ### with points colored by the same color as in the heatmap
plot(cd$PF_READS, cd$PCT_RIBOSOMAL_BASES, log="xy")


## scores that should be plotted the same way as the ribosomal bases
## (one still has to see whether one should take log and/or center)
# fraction aligned
plot(cd$PF_READS, cd$PF_READS_ALIGNED / cd$PF_READS, log="xy")
# mismatch rate
plot(cd$PF_READS, cd$PF_MISMATCH_RATE, log="x")
# pct adapter
plot(cd$PF_READS, cd$PCT_ADAPTER, log="x")
# cpt mrna bases
plot(cd$PF_READS, cd$PCT_MRNA_BASES, log="xy")
# median cv coverage
plot(cd$PF_READS, cd$MEDIAN_CV_COVERAGE, log="xy")
# median 5 prime bias
plot(cd$PF_READS, cd$MEDIAN_5PRIME_BIAS, log="x")
# median 3 prime biase
plot(cd$PF_READS, cd$MEDIAN_3PRIME_BIAS, log="x")

## percent duplication
plot(cd$PF_READS, cd$PERCENT_DUPLICATION, log="x")




######### stuff below can be ignored



plot(cd$TOTAL_READS, cd$PF_READS)
plot(cd$PF_READS, cd$PF_READS_ALIGNED+1, log="xy")
abline(a=0, b=1)
points(cd$PF_READS, cd$PF_HQ_ALIGNED_READS+1, col="red")

plot(cd$PF_READS, cd$PF_MISMATCH_RATE, log="x")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$PF_MISMATCH_RATE), log10 = FALSE)
plot(cd$PF_READS, cd$PF_HQ_ERROR_RATE, log="x")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$PF_HQ_ERROR_RATE), log10 = FALSE)
plot(cd$PF_READS, cd$PF_HQ_MEDIAN_MISMATCHES, log="x")
plot(cd$PF_READS, cd$BAD_CYCLES, log="x")
plot(cd$PF_READS, cd$STRAND_BALANCE, log="x")
## --> cells with unbalanced strand --> rRNA ???

plot(cd$PF_READS, cd$PCT_CHIMERAS, log="x")
plot(cd$PF_READS, cd$PCT_ADAPTER, log="x")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$PCT_ADAPTER), log10 = FALSE)
plot(cd$PF_READS, cd$RIBOSOMAL_BASES, log="xy")
plot(cd$CODING_BASES, cd$RIBOSOMAL_BASES / cd$CODING_BASES, log="x")

plot(cd$CODING_BASES, cd$UTR_BASES / cd$CODING_BASES, log="x")
plot(cd$UTR_BASES, cd$INTRONIC_BASES / cd$UTR_BASES, log="x")
plot(cd$UTR_BASES, cd$INTERGENIC_BASES/ cd$UTR_BASES, log="xy")

plot(cd$PF_READS, cd$IGNORED_READS, log="x")
plot(cd$PF_READS, cd$INCORRECT_STRAND_READS, log="x")
plot(cd$PF_READS, cd$NUM_R1_TRANSCRIPT_STRAND_READS, log="xy")
plot(cd$PF_READS, cd$NUM_R2_TRANSCRIPT_STRAND_READS, log="xy")
plot(cd$PF_READS, cd$NUM_UNEXPLAINED_READS, log="xy")
plot(cd$PF_READS, cd$PCT_RIBOSOMAL_BASES, log="xy")
plot(cd$PF_READS, cd$nGenesDetected, log="xy")

plot(cd$PF_READS, cd$PCT_CODING_BASES, log="xy")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$PCT_CODING_BASES), log10 = FALSE)
plot(cd$PCT_UTR_BASES, cd$PCT_MRNA_BASES, log="xy")
plot(cd$PF_READS, cd$PCT_UTR_BASES / cd$PCT_MRNA_BASES, log="xy")
plot(cd$PF_READS, cd$PCT_USABLE_BASES, log="xy")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$PCT_USABLE_BASES), log10 = FALSE)
plot(cd$PF_READS, cd$MEDIAN_CV_COVERAGE, log="xy")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$MEDIAN_CV_COVERAGE), log10 = FALSE, center=FALSE)
plot(cd$PF_READS, cd$MEDIAN_5PRIME_BIAS, log="x")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$MEDIAN_5PRIME_BIAS), log10 = FALSE, center=FALSE)
plot(cd$PF_READS, cd$MEDIAN_3PRIME_BIAS, log="x")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$MEDIAN_3PRIME_BIAS), log10 = FALSE, center=FALSE)
plot(cd$PF_READS, cd$MEDIAN_5PRIME_TO_3PRIME_BIAS, log="x")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$MEDIAN_5PRIME_TO_3PRIME_BIAS), log10 = FALSE, center=FALSE)
plot(cd$PF_READS, cd$SECONDARY_OR_SUPPLEMENTARY_RDS, log="x")
heatmapPlate(mapToPlate(plateRow, plateCol, cd$SECONDARY_OR_SUPPLEMENTARY_RDS / cd$PF_READS_ALIGNED), log10 = TRUE, center=TRUE)
plot(cd$PF_READS, cd$UNMAPPED_READS, log="xy")
abline(a=0, b=1)
heatmapPlate(mapToPlate(plateRow, plateCol, cd$UNMAPPED_READS / cd$PF_READS), log10 = TRUE, center=FALSE)

plot(cd$PF_READS, cd$UNPAIRED_READ_DUPLICATES, log="xy")
abline(a=0, b=1)
plot(cd$PF_READS, cd$READ_PAIR_OPTICAL_DUPLICATES, log="x")
abline(a=0, b=1)
plot(cd$PF_READS, cd$READ_PAIR_DUPLICATES, log="x")
abline(a=0, b=1)
plot(cd$PF_READS, cd$PERCENT_DUPLICATION, log="x")
abline(a=0, b=1)



# 
# c("Condition [Factor]", "Species", "FragmentSize [Characteristic]", 
#   "SampleConc [Characteristic]", "Tube [Characteristic]", "Sample Id [B-Fabric]", 
#   "PlatePosition [Characteristic]", "LibConc_100_800bp [Characteristic]", 
#   "LibConc_qPCR [Characteristic]", "Adapter1", "strandMode", "LibraryPrepKit", 
#   "EnrichmentMethod", "InputAmount [Characteristic]", "Read Count", 
#   "BAM [File]", "BAI [File]", "STARLog [File]", "PreprocessingLog [File]", 
#   "CountMatrix [File]", "Stats [File]", "CellCyclePhase [File]", 
#   "bias5", "bias3", "CATEGORY", "TOTAL_READS", "PF_READS", "PCT_PF_READS", 
#   "PF_NOISE_READS", "PF_READS_ALIGNED", "PCT_PF_READS_ALIGNED", 
#   "PF_ALIGNED_BASES", "PF_HQ_ALIGNED_READS", "PF_HQ_ALIGNED_BASES", 
#   "PF_HQ_ALIGNED_Q20_BASES", "PF_HQ_MEDIAN_MISMATCHES", "PF_MISMATCH_RATE", 
#   "PF_HQ_ERROR_RATE", "PF_INDEL_RATE", "MEAN_READ_LENGTH", "READS_ALIGNED_IN_PAIRS", 
#   "PCT_READS_ALIGNED_IN_PAIRS", "BAD_CYCLES", "STRAND_BALANCE", 
#   "PCT_CHIMERAS", "PCT_ADAPTER", "PF_BASES", "PF_ALIGNED_BASES", 
#   "RIBOSOMAL_BASES", "CODING_BASES", "UTR_BASES", "INTRONIC_BASES", 
#   "INTERGENIC_BASES", "IGNORED_READS", "CORRECT_STRAND_READS", 
#   "INCORRECT_STRAND_READS", "NUM_R1_TRANSCRIPT_STRAND_READS", "NUM_R2_TRANSCRIPT_STRAND_READS", 
#   "NUM_UNEXPLAINED_READS", "PCT_R1_TRANSCRIPT_STRAND_READS", "PCT_R2_TRANSCRIPT_STRAND_READS", 
#   "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", 
#   "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "PCT_USABLE_BASES", 
#   "PCT_CORRECT_STRAND_READS", "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS", 
#   "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS", "UNPAIRED_READS_EXAMINED", 
#   "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS", "UNMAPPED_READS", 
#   "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES", 
#   "PERCENT_DUPLICATION")

plot(cd$bias5, cd$bias3, xlim=c(0, 5), ylim=c(0,5))
cd$`PlatePosition [Characteristic]`
heatmapPlate(mapToPlate(plateRow, plateCol, cd$TOTAL_READS))
