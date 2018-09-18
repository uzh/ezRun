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
  # The following Rsamtools requires mapped Bam file.
  # require(Rsamtools)
  # header <- scanBamHeader(inBam)
  # readGroupNames <- header[[inBam]]$text[names(header[[inBam]]$text) == "@RG"]
  # readGroupNames <- sub("ID:", "", sapply(readGroupNames, "[", 1))
  cwd <- getwd()
  inBam <- normalizePath(inBam)
  
  # This is more general approach for uBam and mapped Bam.
  tempDIR <- file.path(cwd, paste("samtools-split", Sys.getpid(), sep="-"))
  setwdNew(tempDIR)
  
  cmd <- paste("samtools split -@", mc.cores, "-f %!", inBam)
  ezSystem(cmd)
  
  setwd(cwd)
  
  readGroupNames <- list.files(path=tempDIR)
  bamFns <- paste0(readGroupNames, ".bam")
  file.rename(from=file.path(tempDIR, readGroupNames), to=bamFns)
  names(bamFns) = readGroupNames
  
  # Clearning
  unlink(tempDIR, recursive=TRUE)
  
  return(bamFns)
}

mergeBamAlignments <- function(alignedBamFn, unmappedBamFn,
                               outputBamFn, fastaFn,
                               keepUnmapped=FALSE){
  setEnvironments("picard")
  ## Use . as tmp dir. Big bam generates big tmp files.
  cmd <- paste("java -Djava.io.tmpdir=. -jar",
               Sys.getenv("Picard_jar"), "MergeBamAlignment",
               paste0("ALIGNED=", alignedBamFn),
               paste0("UNMAPPED=", unmappedBamFn),
               paste0("O=", outputBamFn),
               paste0("R=", fastaFn)
  )
  ezSystem(cmd)
  if(!keepUnmapped){
    setEnvironments("samtools")
    tempBam <- tempfile(pattern="noUnmapped", tmpdir = ".",
                        fileext = ".bam")
    cmd <- paste("samtools view -F 4 -h -b", outputBamFn, ">",
                 tempBam)
    ezSystem(cmd)
    file.remove(outputBamFn)
    file.rename(from=tempBam, to=outputBamFn)
  }
  invisible(outputBamFn)
}

posSpecErrorBam <- function(bamGA, genomeFn){
  require(GenomicAlignments)
  require(stringr)
  require(GenomicRanges)
  require(Hmisc)
  require(BSgenome)
  require(Biostrings)
  require(Rsamtools)
  
  nMaxReads <- 100000
  
  what <- c("qname", "seq")
  stopifnot(all(what %in% colnames(mcols(bamGA))))
  
  hasGap <- grepl("I|D|N", cigar(bamGA))
  readLength <- qwidth(bamGA)
  genomeFnIndixe <- FaFile(genomeFn)
  genomeLengths <- seqlengths(genomeFnIndixe)
  
  isOutOfRange <- start(bamGA) + readLength - 1 > genomeLengths[as.character(seqnames(bamGA))] |
   start(bamGA) < readLength
  # this is very conservative; needed because there might be clipped bases in the beginning
  if (any(isOutOfRange)){
    ezWrite("#reads out of range: ", sum(isOutOfRange))
    idx = which(isOutOfRange)
    idx = idx[1:min(10, length(idx))]
    badAlignments = data.frame(pos=start(bamGA)[idx],
                               cigar=cigar(bamGA)[idx],
                               width=readLength[idx])
    print(badAlignments)
  }
  
  indexKeep <- which(!hasGap & !isOutOfRange)
  # indexKeep <- which(!hasGap)
  if (length(indexKeep) > nMaxReads){
    indexKeep <- sample(indexKeep, size=nMaxReads, replace=FALSE)
  }
  ezWrite("#alignments: ", length(bamGA),
          " #valid alignments: ", sum(!hasGap & !isOutOfRange),
          " #used:", length(indexKeep))
  if (length(indexKeep) == 0){
    return(NULL)
  }
  bamGA <- bamGA[indexKeep]
  bamGR <- granges(bamGA) ##used to record the actual mapped coordinates of seq
  
  ## adjust the start POS according to H and/or S
  tempCigar = str_extract(cigar(bamGA), "^(\\d+H)?(\\d+S)?\\d+M")
  ## get the number of H at the beginning
  noOfH <- explodeCigarOpLengths(tempCigar, ops="H")
  stopifnot(!any(lengths(noOfH) > 1))
  noOfH[lengths(noOfH) == 0] <- 0
  noOfH <- as.integer(noOfH)

  ## get the number of S at the beginning
  noOfS <- explodeCigarOpLengths(tempCigar, ops="S")
  stopifnot(!any(lengths(noOfS) > 1))
  noOfS[lengths(noOfS) == 0] <- 0
  noOfS <- as.integer(noOfS)
  nBeginClipped = noOfH + noOfS
  start(bamGR) = start(bamGR) - nBeginClipped
  
  ## add - to the begin and end of SEQ for revComplement
  Xbegin = makeNstr("-", noOfH)
  tempCigar = str_extract(cigar(bamGA), "\\d+[^HS](\\d+S)?(\\d+H)?$")
  noOfH <- explodeCigarOpLengths(tempCigar, ops="H")
  stopifnot(!any(lengths(noOfH) > 1))
  noOfH[lengths(noOfH) == 0] <- 0
  noOfH <- as.integer(noOfH)
  
  Xend = makeNstr("-", noOfH)
  noOfS <- explodeCigarOpLengths(tempCigar, ops="S")
  stopifnot(!any(lengths(noOfS) > 1))
  noOfS[lengths(noOfS) == 0] <- 0
  noOfS <- as.integer(noOfS)
  nEndClipped = noOfH + noOfS
  mcols(bamGR)$seq = DNAStringSet(paste0(Xbegin, mcols(bamGA)$seq, Xend))
  end(bamGR) <- end(bamGR) + nEndClipped
  
  indexNeg = as.logical(strand(bamGA) == "-")
  mcols(bamGR)$seq[indexNeg] <- reverseComplement(mcols(bamGR)$seq[indexNeg])
  
  seqChar = strsplit(as.character(mcols(bamGR)$seq), "")
  readLength <- lengths(seqChar)
  
  maxLength = quantile(readLength, 0.95)
  if (maxLength < max(readLength)){
    readLength[readLength > maxLength] = maxLength
    seqChar = mapply(function(x, l){x[1:l]}, seqChar, readLength)
    bamGR <- resize(bamGR, width=readLength, fix="start")
  }
  
  # referenceChar <- genome[bamGR]
  referenceChar <- getSeq(genomeFnIndixe, bamGR)
  stopifnot(all(width(referenceChar) == width(bamGR))) # check the position change is wright
  referenceChar = strsplit(as.character(referenceChar), "")
  
  # assuming we have unique read length and set it to the maximal read length here.
  nEndTrimmed = maxLength - readLength
  trimmedMatrix = mapply(function(readLength, nEndTrimmed){
                           rep(c(FALSE, TRUE), c(readLength, nEndTrimmed))}, 
                         readLength, nEndTrimmed, SIMPLIFY=FALSE)
  ## build a clippedMatrix to record the clipped character
  nNormal = readLength - nBeginClipped - nEndClipped
  clippedMatrixNoTrim = mapply(function(nBeginClipped, nNormal, nEndClipped){
    rep(c(TRUE, FALSE, TRUE),
        c(nBeginClipped, nNormal, nEndClipped))}, 
    nBeginClipped, nNormal, nEndClipped, SIMPLIFY=FALSE)
  
  matchMatrix = mapply("==", referenceChar, seqChar, SIMPLIFY=FALSE)
  
  clippedMatrixNoTrim[indexNeg] = lapply(clippedMatrixNoTrim[indexNeg], rev)
  clippedMatrix = mapply(function(clippedMatrixNoTrim, nEndTrimmed){
    c(clippedMatrixNoTrim, rep(FALSE, nEndTrimmed))}, 
    clippedMatrixNoTrim, nEndTrimmed, SIMPLIFY=FALSE)
  
  matchMatrix = sapply(matchMatrix, function(x){x[1:maxLength]})
  clippedMatrix = matrix(unlist(clippedMatrix), ncol=length(clippedMatrix))
  clippedRate = rowMeans(clippedMatrix, na.rm=TRUE)
  trimmedMatrix = matrix(unlist(trimmedMatrix), ncol=length(trimmedMatrix))
  trimmedRate = rowMeans(trimmedMatrix, na.rm=TRUE)
  matchMatrix[clippedMatrix] = NA
  errorRate = 1 - rowMeans(matchMatrix, na.rm=TRUE)
  names(errorRate) = 1:nrow(matchMatrix)
  names(clippedRate) = 1:nrow(matchMatrix)
  names(trimmedRate) = 1:nrow(matchMatrix)
  return(list(trimmedRate=trimmedRate, clippedRate=clippedRate, 
              errorRate=errorRate))
}

### Create the MD tag for BAM file and replace matches with "=" in seq
calmdBam <- function(bamFns, mdBamFns=sub("\\.bam$", "_md.bam", bamFns),
                     genomeFn, mc.cores=4L){
  setEnvironments("samtools")
  stopifnot(length(bamFns) == length(mdBamFns))
  
  cmdsCalmd <- paste("samtools calmd -b -e", bamFns, genomeFn, ">", mdBamFns)
  cmdsIndex <- paste("samtools index", mdBamFns)
  
  indicesExist <- file.exists(mdBamFns)
  if(any(indicesExist)){
    warning("Some md Bam files already exist!", 
            paste(mdBamFns[indicesExist], collapse=" "))
    cmdsCalmd <- cmdsCalmd[!indicesExist]
    cmdsIndex <- cmdsIndex[!indicesExist]
  }
  ezMclapply(cmdsCalmd, ezSystem, mc.preschedule = FALSE, mc.cores=mc.cores)
  ezMclapply(cmdsIndex, ezSystem, mc.preschedule = FALSE, mc.cores=mc.cores)
  return(mdBamFns)
}

### position specific error with md bam
posSpecErrorMDBam <- function(bamGA){
  require(Biostrings)
  hasIndel = grepl("I|D", cigar(bamGA))
  bamGA <- bamGA[!hasIndel]
  
  seqs <- mcols(bamGA)$seq
  negStrand <- strand(bamGA) == "-"
  seqs[negStrand] <- reverseComplement(seqs[negStrand])
  consensusCount <- consensusMatrix(seqs)
  errorRate <- 1 - consensusCount["-", ]/colSums(consensusCount)
  return(errorRate)
}