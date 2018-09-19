.subsampleFromFastq = function(reads, outputDir, sampleno, size){
  require("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
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


.getLastGoodBasePos = function(x, minTailQuality, qualityFilterWindow=4){
  x[1:5] = minTailQuality ## do not search for low qualities in the first 5 bases!!!!  
  lastGood1 = match(TRUE, caTools::runmean(x, qualityFilterWindow, align="left", alg="fast") < minTailQuality, nomatch=length(x)+1) - 1
  return(lastGood1)
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

.pairFastqReads = function(fqFile1, fqFile2, fqOut1, fqOut2, overwrite=FALSE,  doGzip=FALSE){
  require("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
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

.getHomoPoloymerTails = function(readFile, maxTailWidth=50, nYield=1e6){
  
  require("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
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

