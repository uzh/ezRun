###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### Filter ATAC-seq bam file; process by sample;
#### 1. chrM removal
#### 2. low mapping quality filtering (< 10)
#### 3. shift the 5' start
#### 4. mark and remove duplicates: sambamba
atacBamProcess <- function(input=NA, output=NA, param=NA){
  require(GenomicAlignments)
  require(Rsamtools)
  require(ATACseqQC)
  require(rtracklayer)
  
  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")){
    output = input$copy()
    output$setColumn("BAM", paste0(getwd(), "/", input$getNames(), 
                                   "-shifted-nodup.bam"))
    output$setColumn("BAI", paste0(getwd(), "/", input$getNames(), 
                                   "-shifted-nodup.bam.bai"))
    output$dataRoot = NULL
  }
  
  bamFile = input$getFullPaths("BAM")
  localBamFile = .getBamLocally(bamFile)
  if(localBamFile != bamFile){
    on.exit(file.remove(c(localBamFile, paste0(localBamFile, ".bai"))),
            add=TRUE)
  }
  ## For now, we only have paired-end ATAC-seq
  stopifnot(param$paired)
  
  what <- c("qname", "flag", "mapq", "isize", "seq", "qual", "mrnm")
  tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
  flag <- scanBamFlag(isSecondaryAlignment = FALSE, 
                      isUnmappedQuery = FALSE, 
                      isNotPassingQualityControls = FALSE)
  
  ## mapq < 10: discarded
  scanBamParam <- ScanBamParam(flag=flag, tag=tags, what=what, mapqFilter=10L)
  
  reads <- readGAlignmentPairs(file=localBamFile, param=scanBamParam)
  
  ## Remove the reads on mitochondrial chr
  mitChrs <- c("M", "MT", "chrM")
  reads <- reads[!seqnames(reads) %in% mitChrs]
  
  ## shiftAlignments
  firstReads <- ATACseqQC:::shiftReads(first(reads),
                                       positive = 4L, negative = 5L)
  lastReads <- ATACseqQC:::shiftReads(last(reads),
                                      positive = 4L, negative = 5L)
  rm(reads)
  gc()
  shiftedReads <- GAlignmentPairs(first=firstReads, last=lastReads)
  
  tempShiftedBam <- tempfile(pattern="shifted_", tmpdir=".", fileext = ".bam")
  export(shiftedReads, tempShiftedBam)
  
  return(output)
}


dupBam <- function(inBam, outBam, operation=c("mark", "remove"),
                   cores=ezThreads()){
  operation <- match.arg(operation)
  setEnvironments("sambamba")
  
  if(operation == "mark"){
    cmd <- paste("sambamba markdup -t", cores, "-l 9 --tmpdir=.",
                 inBam, outBam)
  }else{
    cmd <- paste("sambamba markdup -r -t", cores, "-l 9 --tmpdir=.",
                 inBam, outBam)
  }
  ezSystem(cmd)
  invisible(outBam)
}

