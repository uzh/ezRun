###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Loads the count dataset
##' @description Loads the count dataset with the given input.
##' @template input-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{dataRoot}{ the root directory of the files.}
##'   \item{expressionName}{ if specified, this will be used as the column name...}
##'   \item{knownExpressionNames}{ ...or otherwise known expression names that occur in the dataset will be used.}
##'   \item{ezRef@@refBuild}{ if specified, the sequence annotation will be extracted from \code{ezFeatureAnnotation()}.}
##'   \item{useTranscriptType}{ if specified, only the defined transcript type will be used.}
##'   \item{sigThresh}{ the threshold...}
##'   \item{useSigThresh}{ ...and whether it should be used.}
##'   \item{featureLevel}{ if equal to "gene" and the feature level of the dataset to "isoform", the rawdata will be passed to \code{aggregateCountsByGene()} before returning it.}
##' }
##' @template roxygen-template
##' @return Returns a list of raw data.
##' @seealso \code{\link{ezFeatureAnnotation}}
##' @seealso \code{\link{aggregateCountsByGene}}
##' @examples
##' param = ezParam()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' file = system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file)
##' cds = loadCountDataset(input, param)
loadCountDataset = function(input, param){

  files = input$getFullPaths(param, "Count")
  suffix = unique(toupper( sub(".*\\.", "", files)))
  if (length(suffix) > 1){
    return(list(error=paste("different  file suffixes not supported: <br>",
                            paste(files, collapse="<br>"))))
  }
  x = ezRead.table(files[1])
  if (ezIsSpecified(param$expressionName)){
    columnName = param$expressionName
  } else {
    columnName = intersect(param$knownExpressionNames, colnames(x))[1]
  }
  if (!columnName %in% colnames(x)){
    return(list(error=paste0("Specified column name not found in data!<br>columnName: '", columnName, "'\n",
                            "<br>Available column names:<br>\n",
                            paste0("'", colnames(x), "'", collapse="<br>"),
                            "<br>Set the option columnName to one of the names above!")))
  }
  dataFeatureLevel = unique(input$getColumn("featureLevel"))
  stopifnot(length(dataFeatureLevel) == 1)
  if (ezIsSpecified(param$ezRef@refBuild)){
    seqAnno = ezFeatureAnnotation(param, rownames(x), dataFeatureLevel)
  } else {
    seqAnno = x[ , intersect(c("type", "gene_name", "gene_id", "transcript_id", "Description", "GO BP", "GO MF", "GO CC"), colnames(x)), drop=FALSE]
  }
  signal = ezMatrix(0, rows=rownames(seqAnno), cols=names(files))
  columnNameStart = grep(paste(columnName, "[first"), colnames(x), fixed=TRUE, value=TRUE)
  if (length(columnNameStart) == 1){
    signalStart = signal
  } else {
    signalStart = NULL
  }
  columnNameEnd = grep(paste(columnName, "[last"), colnames(x), fixed=TRUE, value=TRUE)
  if (length(columnNameEnd) == 1){
    signalEnd = signal
  } else {
    signalEnd = NULL
  }
  for (i in 1:length(files)){
    message("loading file: ", files[i])
    x = ezRead.table(files[i], strip.white = FALSE)
    if(!setequal(rownames(x), rownames(seqAnno))){
      if (all( rownames(seqAnno) %in% rownames(x))){
        warning("inconsistent ID set")
      } else {
        stop("later arrays have IDs not present in the first array")
      }
    }
    y = x[rownames(seqAnno), columnName]
    y[is.na(y)] = 0
    signal[ , i] = y
    if (!is.null(signalStart)){
      y = x[rownames(seqAnno), columnNameStart]
      y[is.na(y)] = 0
      signalStart[ , i] = y
    }
    if (!is.null(signalEnd)){
      y = x[rownames(seqAnno), columnNameEnd]
      y[is.na(y)] = 0
      signalEnd[ , i] = y
    }
  }
  
  if (ezIsSpecified(param$useTranscriptType)){
    use = seqAnno$type == param$useTranscriptType
  } else {
    use = TRUE
  }
  
  if (param$useSigThresh){
    sigThresh = param$sigThresh
  } else {
    sigThresh = 0
  }
  
  rawData = list(counts=signal[use, ,drop=FALSE], countsStart=signalStart[use, ,drop=FALSE], 
                 countsEnd=signalEnd[use, ,drop=FALSE], isLog=FALSE,
                 presentFlag=signal[use, ,drop=FALSE] > sigThresh, 
                 seqAnno=seqAnno[use, , drop=FALSE], featureLevel=dataFeatureLevel,
                 type="Counts", countName=columnName, dataset=input$meta)
  if (dataFeatureLevel == "isoform" && param$featureLevel == "gene"){
    rawData = aggregateCountsByGene(param, rawData)
  }
  rawData$rpkm = getRpkm(rawData)  
  rawData$tpm = getTpm(rawData)
  return(rawData)
}

##' @title Writes the head of a file
##' @description Writes the head of a file into a newly created target file.
##' @param target a character specifying the path of the output.
##' @param x a character specifying the path to a file to read from.
##' @param n an integer specifying the amount of lines to return.
##' @template roxygen-template
##' @return Returns the name of the output file.
##' @examples 
##' description = system.file("DESCRIPTION", package="ezRun", mustWork=TRUE)
##' ezHead(x=description)
ezHead = function(target=paste0(x, "_head"), x, n=1000){
  ezSystem(paste("head -n", n, x, ">", target))
  return(target)
}

##' @title Reads NcPro results
##' @description Reads NcPro results and performs some checks on the data.
##' @param ncproDir a character representing the file path to the directory that contains the NcPro results.
##' @template roxygen-template
##' @return Returns the rawdata read from the NcPro results.
readNcProResult = function (ncproDir) {
  dataFiles = list.files(ncproDir, "_subfamcov.data$", full.names=TRUE)
  
  df = dataFiles[1]
  allData = NULL
  allTypes = character()
  for (df in dataFiles){
    type = sub("_all_samples_subfamcov.data", "", basename(df))
    x = ezRead.table(df)
    if (is.null(allData)){
      allData = x
    } else {
      stopifnot(colnames(allData) == colnames(x))
      #stopifnot(!rownames(x) %in% rownames(allData))
      allData = rbind(allData, x)
    }
    allTypes = c(allTypes, rep(type, nrow(x)))
  }
  seqAnno = data.frame(type=I(allTypes), row.names=rownames(allData))
  typeNames = c( "rfam_ACA_snoRNA"="ACA_snoRNA",
                 "rfam_CD_snoRNA"="CD_snoRNA",
                 "rmsk_L1"="L1",
                 "tRNA_tRNA"="tRNA",
                 "tRNA_tRNA_e_+30_+30"="tRNA_e",
                 "precursor_miRNA_miRNA"="pre-miRNA",
                 "mature_miRNA_miRNA_e_+2_+2"="miRNA")
  types = data.frame(row.names=rownames(seqAnno))
  stopifnot(allTypes %in% names(typeNames))
  for (nm in names(typeNames)){
    types[[typeNames[nm]]] = seqAnno$type == nm
  }
  rawData = list(counts=as.matrix(allData), presentFlag = allData > 0, seqAnno=seqAnno, featureLevel="isoform", type="Counts", countName="ncpro-counts", isLog=FALSE,
                 types=types)
}

##' @title Filters FastQ files by bam
##' @description Filters FastQ files by bam and writes them into new files.
##' @param fqFiles a character vector representing file paths to FastQ files.
##' @template bamFile-template
##' @param fqOutFiles an optional character vector representing file paths to the FastQ output files.
##' @param doGzip a logical indicating whether to archive the output files in a gzip archive.
##' @param keepUnmapped passed further to \code{scanBamFlag()}.
##' @param isProperPair passed further to \code{scanBamFlag()}.
##' @param readIds a character containing the read ID's from the \code{bamFile}.
##' @template roxygen-template
##' @seealso \code{\link[Rsamtools]{scanBam}}
##' @seealso \code{\link[Rsamtools]{ScanBamParam}}
##' @return Returns the names of the filtered FastQ files.
##' @examples 
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file)
##' param = list()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' fqFiles = input$getFullPaths(param, "Read1")
##' fqFiltered = filterFastqByBam(fqFiles, bamFile, doGzip = FALSE)
filterFastqByBam = function(fqFiles, bamFile, fqOutFiles=NULL, doGzip=TRUE, keepUnmapped=TRUE, isProperPair=NA){
  param = ScanBamParam(what=c("qname"),
                     flag=scanBamFlag(isUnmappedQuery=!keepUnmapped, isProperPair=isProperPair))
  bamReads = unique(scanBam(bamFile, param=param)[[1]]$qname)
  gc()
  fqFiltered = removeReadsFromFastq(fqFiles, bamReads, fqOutFiles=fqOutFiles, doGzip=doGzip)
  return(fqFiltered)
}

##' @describeIn filterFastqByBam Removes reads from the FastQ files according to the bam filter.
removeReadsFromFastq = function(fqFiles, readIds, fqOutFiles=NULL, doGzip=TRUE){
  requireNamespace("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  if (is.null(fqOutFiles)){
    fqOutFiles = sub(".gz$", "", basename(fqFiles))
    fqOutFiles = sub("R1.fastq", "clean_R1.fastq", fqOutFiles)
    fqOutFiles = sub("R2.fastq", "clean_R2.fastq", fqOutFiles)
  }
  for (i in 1:length(fqFiles)){
    fqs = FastqStreamer(fqFiles[i], 1e6)
    count = 0
    while (length(x <- yield(fqs))) {
      count = count + 1
      message(fqOutFiles[i], " ", count)
      ids = sub(" .*", "", as.character(id(x)))
      keep = !ids %in% readIds
      writeFastq(x[keep], file=fqOutFiles[i], mode="a", compress=FALSE)
    }
    if (doGzip){
      ezSystem(paste("gzip", fqOutFiles[i]))
      fqOutFiles[i] = paste0(fqOutFiles[i], ".gz")
    }
  }
  return(fqOutFiles)
}

##' @title Counts reads in FastQ files
##' @description Counts reads in FastQ files and returns the amount for each file.
##' @param fastqFiles a character vector representing file paths to FastQ files.
##' @template roxygen-template
##' @return Returns a named numeric vector
##' @examples 
##' file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file)
##' param = list()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' fqFiles = input$getFullPaths(param, "Read1")
##' result = countReadsInFastq(fqFiles)
countReadsInFastq = function(fastqFiles){
  nReads = c()
  for(fastq in fastqFiles){
    if(grepl(".gz$", fastq)){
      cmd = paste("gunzip -c", fastq, "| wc -l")
      nRead = ezSystem(cmd, intern=TRUE)
    }else{
      cmd = paste("wc -l", fastq)
      nRead = ezSystem(cmd , intern=TRUE)
    }
    nRead = trimWhiteSpace(nRead)
    nRead = ezSplit(nRead, " ")[1]
    stopifnot(as.integer(nRead) %% 4 == 0)
    nRead = as.integer(nRead) / 4
    nReads[fastq] = nRead
  }
  return(nReads)
}





.pairFastqReads = function(fqFile1, fqFile2, fqOut1, fqOut2, overwrite=FALSE,  doGzip=FALSE){
  requireNamespace("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  for (fn in c(fqOut1, fqOut2)){
    if (file.exists(fn)){
      if (overwrite){
        file.remove(fn)
      } else {
        stop("output file exists and overwrite is FALSE: ", fn)
      }
    }
  }
  fq1 = readFastq(fqFile1)
  ids1 = sub(" .*", "", as.character(id(fq1)))
  remove(fq1)
  gc()
  fq2 = readFastq(fqFile2)
  ids2 = sub(" .*", "", as.character(id(fq2)))
  remove(fq2)
  gc()
  ids = intersect(ids1, ids2)
  idx1 = match(ids, ids1)
  idx2 = match(ids, ids2)
  remove(ids, ids1, ids2)
  gc()
  writeFastq(readFastq(fqFile1)[idx1], fqOut1, compress=FALSE)
  gc()
  writeFastq(readFastq(fqFile2)[idx2], fqOut2, compress=FALSE)
  if (doGzip){
    ezSystem(paste("gzip", fqOut1, fqOut2))
  }
  return(NULL)
}

.getPairedReads = function(reads1){
  filePairs = data.frame(reads1=c("R1_001.fastq.gz", "R1.fastq.gz", "R1.fastq", "F3.csfasta", "F3.csfasta"),
                         reads2=c("R2_001.fastq.gz", "R2.fastq.gz", "R2.fastq", "F5-BC.csfasta", "F5-RNA.csfasta"))
  for (i in 1:nrow(filePairs)){
    #x = paste0(filePairs$reads1[i], "$")
    x = filePairs$reads1[i]
    if (ezGrepl(x, reads1)){
      reads2 = sub(x, filePairs$reads2[i], reads1)
      if(file.exists(reads2)){
        return(reads2)
      }
    }
  }
  return(NULL)
}

.getMeanQuality = function(phredQual, qualOffset=Biostrings:::offset(phredQual), batchSize=100000){
  ## too slow:
  #mg = lapply(phredQual, function(x){mean(as.integer(charToRaw(as.character(x))))-qualOffset})
  if (length(phredQual) <= batchSize){
    intQualConcatenated = as.integer(unlist(phredQual))
    #ezWriteElapsed(job, status="as integer")
    endPos = cumsum(width(phredQual))
    startPos = endPos - width(phredQual) + 1
    meanQual = unlist(mapply(function(s, e){mean(intQualConcatenated[s:e]) - qualOffset}, startPos, endPos, SIMPLIFY=FALSE))
    return(meanQual)
  } else {
    batchStart = seq(from=1, to=length(phredQual), by=batchSize)
    batchEnd = c(batchStart[-1]-1, length(phredQual))
    qo = qualOffset
    mqList = mapply(function(bs, be){.getMeanQuality(phredQual[bs:be], qualOffset=qo, batchSize=batchSize)},
                    batchStart, batchEnd, SIMPLIFY=FALSE)
    return(unlist(mqList))
  }
}

## will fail if sequences are too long
.phredQualToInt = function(phredQual){
  offset = Biostrings:::offset(phredQual)
  ## what we want to get is:
  #qList = lapply(as.character(phredQual), function(x){as.integer(charToRaw(x))-offset})
  intQualConcatenated = as.integer(unlist(phredQual)) - offset
  #ezWriteElapsed(job, status="as integer")
  endPos = cumsum(width(phredQual))
  startPos = endPos - width(phredQual) + 1
  qList = mapply(function(s, e){intQualConcatenated[s:e]}, startPos, endPos, SIMPLIFY=FALSE)
  #ezWriteElapsed(job, status="as list")
  ## option 2
  #ptm <- proc.time()
  #splitFactor = rep(1:length(qwidth), qwidth)
  #qList = split(intQualConcatenated, splitFactor)
  #proc.time() - ptm
  return(qList)  
}

.getQualities = function(shortReadQ){
  #job = ezJobStart(name="getQualities")
  qual = quality(quality(shortReadQ))
  offset = Biostrings:::offset(PhredQuality(qual))
  ## what we want to get is:
  #qList = lapply(as.character(qual), function(x){as.integer(charToRaw(x))-offset})
  intQualConcatenated = as.integer(unlist(qual)) - offset
  #ezWriteElapsed(job, status="as integer")
  endPos = cumsum(width(qual))
  startPos = endPos - width(qual) + 1
  qList = mapply(function(s, e){intQualConcatenated[s:e]}, startPos, endPos, SIMPLIFY=FALSE)
  #ezWriteElapsed(job, status="as list")
  ## option 2
  #ptm <- proc.time()
  #splitFactor = rep(1:length(qwidth), qwidth)
  #qList = split(intQualConcatenated, splitFactor)
  #proc.time() - ptm
  return(qList)
}

.getLastGoodBasePos = function(x, minTailQuality, qualityFilterWindow=4){
  x[1:5] = minTailQuality ## do not search for low qualities in the first 5 bases!!!!  
  lastGood1 = match(TRUE, caTools::runmean(x, qualityFilterWindow, align="left", alg="fast") < minTailQuality, nomatch=length(x)+1) - 1
  return(lastGood1)
}

.getSingleReadFileLocally = function(src, adapter=NA, param=NULL){
  
  job = ezJobStart(paste("get", src))  
  
  if (!is.na(adapter)){
    adapter = paste0(adapter, paste(rep("N", 100), collapse=""))
  }
  
  target=sub(".gz$", "", basename(src))
  stopifnot(target != src)
  if (file.exists(target)) file.remove(target)
  nYield = 5e5
  idx = NULL
  param$nReads = as.integer(param$nReads)
  
  fqs = FastqStreamer(src, nYield)
  if (ezIsSpecified(param$subsampleReads) && param$subsampleReads){
    idx = seq(from=1, to=nYield, by=param$subsampleReads)
  } else {
    if(param$nReads > 0){ ## nReads is set to -1L by default!
      nTotalReads = countReadsInFastq(src)
      ezWriteElapsed(job, status="reads counted")
      if (nTotalReads > param$nReads){
        nYield = min(nYield, nTotalReads)
        idx = round(seq(from=1, to=nYield, length=(round(nYield/nTotalReads * param$nReads))))
      }
    }
  }
  while(length(x <- yield(fqs))){
    if (!is.null(idx)){
      if (length(x) < nYield){
        ## when reading the last we want to truncate idx
        idx = idx[idx <= length(x)]
      }
      x = x[idx]
    }
    if (!is.null(param$minAvgQuality) && !is.na(param$minAvgQuality) && param$minAvgQuality > 0 ){
      minAvgQuality = param$minAvgQuality
      use = unlist(ezMclapply(.getQualities(x), function(qVec){mean(qVec) >= minAvgQuality}, mc.cores=round(ezThreads()/2)))
      if (any(is.na(use))){
        message("encountered na in minAvgQuality: ", sum(is.na(use)))
        use[is.na(use)] = FALSE  
      }
      x = x[use]
    }
    endPos = width(x)
    if (!is.null(param$minTailQuality) && !is.na(param$minTailQuality) && param$minTailQuality > 0 ){
      minTailQuality = param$minTailQuality
      lastGoodBase = unlist(ezMclapply(.getQualities(x), .getLastGoodBasePos, minTailQuality, 
                                       qualityFilterWindow=param$qualityFilterWindow, mc.cores=round(ezThreads()/2)))
      use = lastGoodBase >= param$minReadLength + param$trimLeft
      endPos[use] = lastGoodBase[use]
    }
    if (!is.null(param$trimAdapter) && param$trimAdapter){
      stopifnot(!is.null(adapter))
      stopifnot(adapter != "")
      endPos2 = end(trimLRPatterns(Rpattern=adapter, subject=narrow(sread(x), end=endPos), max.Rmismatch=0.15, with.Rindels=TRUE, ranges=TRUE))
      endPos2[endPos2 < param$minReadLength + param$trimLeft] = param$minReadLength  + param$trimLeft
      hasLongAdapt = endPos - endPos2 >= 5
      #ezWrite("trimming adapters: ", sum(hasLongAdapt), " / ", length(hasLongAdapt))
      endPos[hasLongAdapt] = endPos2[hasLongAdapt]
    }
    if (param$trimRight > 0) {
      endPos = endPos - param$trimRight ### always trim right bases ##pmin(endPos, width(x) - param$trimRight)
    }
    trimLeft = pmin(param$trimLeft+1, width(x))
    writeFastq(narrow(x, start=trimLeft, end=pmax(endPos, trimLeft)), file=target, mode="a", compress=F)
  }
  close(fqs)
  ezWriteElapsed(job)
  size = ezSystem(paste("du -h", target), echo=FALSE, intern=TRUE)
  size = sub("\t", " ", size)
  message("file size: ", size)
  return(target)
}

.ezGunzip = function(x, targetDir=".", verbose=FALSE){
  if (ezGrepl(".gz$", x)){
    start=proc.time()
    target = file.path(targetDir, sub(".gz$", "", basename(x)))
    ezSystem(paste("gunzip -c", x, ">", target))
    size = ezSystem(paste("du -h", target), echo=FALSE, intern=TRUE)
    size = sub("\t", " ", size)
    message("size: ", size, " time: ", signif((proc.time() - start)[ "elapsed"]/60, digits=4), " min")
    return(target)
  } else {
    return(x)
  }
}

.subsampleFromFastq = function(reads, outputDir, sampleno, size){
  requireNamespace("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  for(i in 1:length(reads)){
    reads1 = reads[i]
    reads2 = .getPairedReads(reads1)
    f1 = FastqFile(reads1)
    r1 = readFastq(f1)
    if(!is.null(reads2)){
      f2 = FastqFile(reads2)
      r2 = readFastq(f2)
    }
    outputDirNow = file.path(outputDir, sub(".fastq$", "", basename(reads1)))
    dir.create(outputDirNow, recursive=TRUE)
    for(j in 1:sampleno){
      output1 = file.path(outputDirNow, paste0("S", j, "_", sub(".fastq$", "", basename(reads1)), ".fastq"))
      if(is.null(reads2)){
        .samplingFastq(reads1=r1, output1=output1, size=size)
      }else{
        output2 = file.path(outputDirNow, paste0("S", j, "_", sub(".fastq$", "", basename(reads2)), ".fastq"))
        .samplingFastq(reads1=r1, reads2=r2, output1=output1, output2=output2, size=size)
      }
    }
    close(f1)
    if(!is.null(reads2)){
      close(f2)
    }
  }
}

.samplingFastq = function(reads1, reads2=NULL, output1, output2=NULL, size=NULL){
  nReads = length(reads1)
  nSampled = NULL
  if(size >=1 && size <= nReads){
    nSampled = size
  }else if(size >0 && size < 1){
    nSampled = size * nReads
  }else{
    stop("wrong size!")
  }
  if(is.null(reads2)){
    writeFastq(sample(reads1, nSampled), file=output1, mode="w", full=TRUE)
  }else{
    seed = sample(1:1e4,1)
    set.seed(seed)
    writeFastq(sample(reads1, nSampled), file=output1, mode="w", full=TRUE)
    set.seed(seed)
    writeFastq(sample(reads2, nSampled), file=output2, mode="w", full=TRUE)
  }
}

.getHomoPoloymerTails = function(readFile, maxTailWidth=50, nYield=1e6){
  
  requireNamespace("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  fqs = FastqStreamer(readFile, nYield)
  template = integer(maxTailWidth+1)
  names(template) = 0:maxTailWidth  
  lengthCounts = list()
  lengthCounts[["A"]] = list(left=template, right=template)
  lengthCounts[["T"]] = list(left=template, right=template)
  lengthCounts[["G"]] = list(left=template, right=template)
  lengthCounts[["C"]] = list(left=template, right=template)
  
  while(length(x <- yield(fqs)) > 0){
    for (base in (c("A", "T", "G", "C"))){
      tail = paste(rep(base, maxTailWidth), collapse="")
      rightSizes = width(x) - end(trimLRPatterns(Rpattern=tail, subject=sread(x), max.Rmismatch=0.1, ranges=TRUE))
      rst = table(rightSizes)
      lengthCounts[[base]]$right[names(rst)] = lengthCounts[[base]]$right[names(rst)] + rst
      leftSizes = start(trimLRPatterns(Lpattern=tail, subject=sread(x), max.Lmismatch=0.1, ranges=TRUE)) -1
      lst = table(leftSizes)
      lengthCounts[[base]]$left[names(lst)] = lengthCounts[[base]]$left[names(lst)] + lst
    }
  }
  close(fqs)
  return(lengthCounts)
}
