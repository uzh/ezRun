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
  filteroutBam(inBam=noDupBam, outBam=noDupNoMTNOLowBam, 
               cores=param$cores, chrs=c("M", "MT", "chrM", "chrMT"), mapQ=10)
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
                   program=c("sambamba", "picard"),
                   cores=ezThreads()){
  require(Rsamtools)
  
  operation <- match.arg(operation)
  program <- match.arg(program)
  
  noDupBam <- tempfile(pattern="nodup_", tmpdir=".", fileext = ".bam")
  
  headers <- scanBamHeader(inBam)
  if(length(headers[[1]]$targets) >= 1000){
    ## sambamba cannot work on highly fragmented genome
    ## 1000 is a random pick.
    program <- "picard"
  }
  
  message("Marking duplicates with ", program)
  
  
  if(program == "sambamba"){
    setEnvironments("sambamba")
    cmd <- paste("sambamba markdup -t", cores, "-l 9 --tmpdir=.",
                 ifelse(operation=="mark", "", "-r"), inBam, outBam)
    ezSystem(cmd)
  }else{
    setEnvironments("picard")
    tempBam <- tempfile(pattern="tempBam", tmpdir=".", fileext = ".bam")
    #ezSortIndexBam(inBam=inBam, bam=tempBam, removeBam=FALSE,
    #               method="picard")
    
    metricFn <- tempfile()
    tempMarkedBam <- tempfile(pattern="tempMarkedBam", tmpdir=".",
                              fileext = ".bam")
    cmd <- paste("java -Djava.io.tmpdir=. -jar", 
                 Sys.getenv("Picard_jar"), "MarkDuplicates",
                 paste0("I=", inBam),
                 paste0("O=", outBam),
                 paste0("M=", metricFn),
                 paste0("REMOVE_DUPLICATES=", 
                        ifelse(operation=="mark", "false", "true")),
                 "> /dev/null")
    ezSystem(cmd)
    #file.remove(c(tempBam, paste0(tempBam, ".bai")))
    
    #ezSortIndexBam(inBam=tempMarkedBam, bam=outBam, removeBam=TRUE)
    ezSystem(paste("samtools", "index", outBam))
  }
  
  invisible(outBam)
}

### Filter bam by removing chrs, low mapQ
filteroutBam <- function(inBam, outBam, cores=ezThreads(), chrs=NULL, mapQ=NULL){
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
  
  cmd <- paste("samtools index", outBam)
  ezSystem(cmd)
  
  invisible(outBam)
}

### Convert bam to bigwig file
bam2bw <- function(file,
                   destination=sub("\\.bam$", ".bw", file, ignore.case = TRUE), 
                   paired=FALSE,
                   #mode=c("RNA-seq", "DNA-seq"), ## TODO: add strandness for RNA-seq
                   method=c("deepTools", "Bioconductor"),
                   cores=ezThreads()){
  method <- match.arg(method)
  
  if(method == "Bioconductor"){
    ## Better compatability; single-base resolutio;
    ## Slower; higher RAM comsuption
    require(rtracklayer)
    require(GenomicAlignments)
    if (paired){
      aligns = readGAlignmentPairs(file)
    } else {
      aligns = readGAlignments(file)
    }
    cov = coverage(aligns)
    export.bw(cov, destination)
  }else if(method == "deepTools"){
    cmd <- paste("bamCoverage", "--bam", file, "-o", destination,
                 "--binSize 10 -of bigwig", "-p", cores
                 )
    ezSystem(cmd)
  }
  invisible(destination)
}

splitBamByRG <- function(inBam, mc.cores=ezThreads()){
  setEnvironments("samtools")
  require(Rsamtools)
  header <- scanBamHeader(inBam)
  readGroupNames <- header[[inBam]]$text[names(header[[inBam]]$text) == "@RG"]
  readGroupNames <- sub("ID:", "", sapply(readGroupNames, "[", 1))
  
  cmd <- paste("samtools split -@", mc.cores, "-f %!", inBam)
  ezSystem(cmd)
  
  bamFns <- paste0(readGroupNames, ".bam")
  file.rename(from=readGroupNames, to=bamFns)
  
  return(bamFns)
}
