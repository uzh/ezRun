
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

DADA2CreateSeqTab <- function(sampleNames,maxLen=0,file1PathInDataset,database,seqTech,maxExpErr){
    fnFs <- file1PathInDataset
    filtFs <- paste0(sampleNames,".filt.R1.fastq.gz")
    names(filtFs) <- sampleNames
    out <- filterAndTrim(fnFs, filtFs, truncLen=maxLen,
                         maxN=0, maxEE=maxExpErr, truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)
    if (seqTech == "pacbio") {
    errF <-  learnErrors(filtFs, BAND_SIZE=32, multithread=TRUE, 
                         errorEstimationFunction=dada2:::PacBioErrfun)
    } else {
    errF <- learnErrors(filtFs, multithread=TRUE)
    }
    derepFs <- derepFastq(filtFs, verbose=FALSE)
    if (seqTech == "454" | seqTech == "ionTorrent"){
      dadaFs <- dada(derepFs, err=errF, multithread=TRUE,
                     HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
    }else if (seqTech == "pacbio") {
      dadaFs <- dada(derepFs, err=errF, multithread=TRUE, BAND_SIZE=32)
    }else {
      dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    }
    seqtab <- makeSequenceTable(dadaFs)
    fullTableOfOTUsNoChim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
    taxa <- assignTaxonomy(fullTableOfOTUsNoChim,database, multithread=FALSE)
    rownames(fullTableOfOTUsNoChim) <- sampleNames
  return(list(fullTableOfOTUsNoChimObj=fullTableOfOTUsNoChim,
              taxaObj=taxa,outFilt=out,dadaObj=dadaFs))
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