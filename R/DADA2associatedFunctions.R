
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Main DADA2 workflow (sample-wise)
##' @description Create  OTU table.
##' @param  a fastq or pairs of ffastq files .
##' @return Returns a DADA2 seqtab object.

DADA2CreateSeqTab <- function(sampleName,maxLen,file1PathInDataset,
                              minLen,file2PathInDataset){
  if(!missing(file2PathInDataset) & !missing(minLen)){
  fnFs <- file1PathInDataset
  fnRs <- file2PathInDataset
  filtFs <- paste("filt",sampleName,"R1.fastq.gz", sep = ".")
  filtRs <- paste("filt",sampleName,"R2.fastq.gz", sep = ".")
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(maxLen,minLen),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  seqtab <- makeSequenceTable(mergers)
  }else{
    fnFs <- file1PathInDataset
    filtFs <- paste("filt",sampleName,"R1.fastq.gz", sep = ".")
    out <- filterAndTrim(fnFs, filtFs, truncLen=minLen,
                         maxN=0, maxEE=1, truncQ=11, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)
    errF <- learnErrors(filtFs, multithread=TRUE)
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    seqtab <- makeSequenceTable(dadaFs)
  }
  rownames(seqtab) <- sampleName
  return(seqtab)
}


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title DADA2 seqtab merging, removing chimeras and assigning taxonomy
##' @description Merges   DADA2 seqTables.
##' @param  a list of  DADA2 seqTables.
##' @return Returns a merged and annotated DADA2 seqtab object, without chimeras and chimera stats.
##' 
DADA2mergeSeqTabs <- function(listOfSeqTabsFiles){
  listOfSeqTabs <- lapply(listOfSeqTabsFiles,readRDS)
  fullTableOfOTUs <- mergeSequenceTables(tables= listOfSeqTabs)
fullTableOfOTUsNoChim <- removeBimeraDenovo(fullTableOfOTUs, method="consensus", multithread=TRUE)
taxa <- assignTaxonomy(fullTableOfOTUsNoChim,SILVA_DB_DADA2, multithread=TRUE)
return(list(fullTableOfOTUsNoChimObj=fullTableOfOTUsNoChim,taxaObj=taxa))
}