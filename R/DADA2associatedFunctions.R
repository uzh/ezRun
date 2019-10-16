
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

DADA2CreateSeqTab <- function(sampleNames,minLen=0,file1PathInDataset,database){
    fnFs <- file1PathInDataset
    filtFs <- paste0(sampleNames,".filt.R1.fastq.gz")
    names(filtFs) <- sampleNames
    out <- filterAndTrim(fnFs, filtFs, truncLen=minLen,
                         maxN=0, maxEE=1, truncQ=11, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)
    errF <- learnErrors(filtFs,errorEstimationFunction = noqualErrfun, multithread=TRUE)
    derepFs <- derepFastq(filtFs, verbose=FALSE)
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    seqtab <- makeSequenceTable(dadaFs)
    fullTableOfOTUsNoChim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
    taxa <- assignTaxonomy(fullTableOfOTUsNoChim,database, multithread=TRUE)
    rownames(fullTableOfOTUsNoChim) <- sampleNames
  return(list(fullTableOfOTUsNoChimObj=fullTableOfOTUsNoChim,taxaObj=taxa))
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
DADA2mergeSeqTabs <- function(listOfSeqTabsFiles, database){
  listOfSeqTabs <- lapply(listOfSeqTabsFiles,readRDS)
  fullTableOfOTUs <- mergeSequenceTables(tables=listOfSeqTabs)
fullTableOfOTUsNoChim <- removeBimeraDenovo(fullTableOfOTUs, method="consensus", multithread=TRUE)
taxa <- assignTaxonomy(fullTableOfOTUsNoChim,database, multithread=TRUE)
return(list(fullTableOfOTUsNoChimObj=fullTableOfOTUsNoChim,taxaObj=taxa))
}