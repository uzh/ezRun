library(ezRun)


library(rtracklayer)
library(GenomicRanges)
## assumes ENSEMBL style gtf files
gtf <- import(
  "/srv/GT/reference/Equus_caballus/Ensembl/EquCab3/Annotation/Release_111-2024-01-17/Genes/genes.gtf"
)
seqLengths <- readDNAStringSet(
  "/srv/GT/reference/Equus_caballus/Ensembl/EquCab3/Sequence/WholeGenomeFasta/genome.fa"
)
seqLengths <- setNames(width(seqLengths), names(seqLengths))
gtfMod <- extendGtfThreePrime(gtf, 2000, seqLengths)


length(gtf) == length(gtfMod)

isEqual <- poverlaps(gtf, gtfMod, type = "equal")
table(isEqual)

gtfD <- gtf[!isEqual]
gtfModD <- gtfMod[!isEqual]

table(gtfD$type)
table(gtfD$transcript_biotype)
table(gtfD$gene_biotype)

isInside <- poverlaps(gtfD, gtfModD, type = "within")
table(isInside)
isOutside <- poverlaps(gtfModD, gtfD, type = "within")
table(isOutside)


gtfDList <- split(gtfD, gtfD$gene_id)
gtfModDList <- split(gtfModD, gtfModD$gene_id)


ovlCount <- countOverlaps(gtfDList, gtfDList, type = "any")
ovlCount2 <- countOverlaps(gtfModDList, gtfModDList, type = "any")
table(ovlCount == ovlCount2)
table(ovlCount, ovlCount2)


delta <- psetdiff(gtfModD, gtfD)

downReg <- flank(gtfD, 2000, start = FALSE)
isInside <- poverlaps(delta, downReg, type = "within")
table(isInside)

upReg <- flank(gtfD, 2000, start = TRUE)
isInside <- poverlaps(delta, upReg, type = "within")
table(isInside)
