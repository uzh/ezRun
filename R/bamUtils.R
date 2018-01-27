###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### Filter ATAC-seq bam file; process by sample;
#### 1. remove duplicates: sambamba
#### 2. chrM removal
#### 3. low mapping quality filtering (< 10)
#### 4. shift the 5' start
atacBamProcess <- function(input=NA, output=NA, param=NA){
  require(GenomicAlignments)
  require(Rsamtools)
  require(rtracklayer)
  
  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")){
    output = input$copy()
    output$setColumn("BAM", paste0(getwd(), "/", input$getNames(), 
                                   "_processed.bam"))
    output$setColumn("BAI", paste0(getwd(), "/", input$getNames(), 
                                   "_processed.bam.bai"))
    output$dataRoot = NULL
  }
  
  bamFile = input$getFullPaths("BAM")
  
  ## For now, we only have paired-end ATAC-seq
  stopifnot(param$paired)
  
  ## remove the duplicates in bam
  noDupBam <- tempfile(pattern="nodup_", tmpdir=".", fileext = ".bam")
  message("Remove duplicates...")
  dupBam(inBam=bamFile, outBam=noDupBam, operation="remove", 
         cores=param$cores)
  
  noDupNoMTNOLowBam <- basename(output$getColumn("BAM"))
  filterBam(inBam=noDupBam, outBam=noDupNoMTNOLowBam, 
            cores=param$cores, chrs=c("M", "MT", "chrM"), mapQ=10)
  file.remove(c(noDupBam, paste0(noDupBam, ".bai")))
  
  if(param$shiftATAC){
    require(ATACseqQC)
    what <- c("qname", "flag", "mapq", "isize", "seq", "qual", "mrnm")
    tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
    flag <- scanBamFlag(isSecondaryAlignment = FALSE, 
                        isUnmappedQuery = FALSE, 
                        isNotPassingQualityControls = FALSE)
    scanBamParam <- ScanBamParam(flag=flag, tag=tags, what=what)
    message("Reading nodup noMT bam file...")
    reads <- readGAlignmentPairs(file=noDupNoMTNOLowBam, param=scanBamParam)
    file.remove(c(noDupNoMTNOLowBam, paste0(noDupNoMTNOLowBam, ".bai")))
    ## shiftAlignments
    message("Shifting 5' start...")
    firstReads <- ATACseqQC:::shiftReads(first(reads),
                                         positive = 4L, negative = 5L)
    lastReads <- ATACseqQC:::shiftReads(last(reads),
                                        positive = 4L, negative = 5L)
    rm(reads)
    message("Exporting bam file...")
    export(GAlignmentPairs(first=firstReads, last=lastReads), 
           noDupNoMTNOLowBam)
    rm(firstReads)
    rm(lastReads)
    gc()
  }
  
  return(output)
}

### Make or remove duplicated in bam file
dupBam <- function(inBam, outBam, operation=c("mark", "remove"),
                   cores=ezThreads()){
  operation <- match.arg(operation)
  setEnvironments("sambamba")
  
  noDupBam <- tempfile(pattern="nodup_", tmpdir=".", fileext = ".bam")
  
  if(operation == "mark"){
    cmd <- paste("sambamba markdup -t", cores, "-l 9 --tmpdir=.",
                 inBam, outBam)
  }else{
    cmd <- paste("sambamba markdup -r -t", cores, "-l 9 --tmpdir=.",
                 inBam, outBam)
  }
  ezSystem(cmd)
  
  ## sort bam by coordinates
  #cmd <- paste("sambamba sort --tmpdir=. -l 9 -o", outBam, "-t", cores,
  #             noDupBam)
  #ezSystem(cmd)
  #file.remove(c(noDupBam, paste0(noDupBam, ".bai")))
  
  invisible(outBam)
}

### Filter bam by removing chrs
filterBam <- function(inBam, outBam, cores=ezThreads(), chrs=NULL, mapQ=NULL){
  require(Rsamtools)
  
  setEnvironments("samtools")
  
  cmd <- paste("samtools view -b", "-@", cores, "-o", outBam)
  
  if(!is.null(mapQ)){
    cmd <- paste(cmd, "-q", mapQ)
  }
  
  cmd <- paste(cmd, inBam)
  
  if(!is.null(chrs)){
    allChrs <- ezBamSeqNames(inBam)
    keepChrs <- setdiff(allChrs, chrs)
    cmd <- paste(cmd, paste(keepChrs, collapse=" "))
  }
  ezSystem(cmd)
  
  indexBam(outBam)
  
  invisible(outBam)
}

### Convert bam to bigwig file
bam2bw <- function(inBam, outBw, paired=FALSE,
                   #mode=c("RNA-seq", "DNA-seq"), ## TODO: add strandness for RNA-seq
                   method=c("deepTools", "Bioconductor")){
  method <- match.arg(method)
  
  if(method == "Bioconductor"){
    
  }else if(method == "deepTools"){
    
  }
  
  invisible(outBw)
}
