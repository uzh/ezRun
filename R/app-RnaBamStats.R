###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodRnaBamStats = function(input=NA, output=NA, param=NA, 
                               htmlFile="00index.html"){
  
  setwdNew(basename(output$getColumn("Report")))
  param$featureLevel = "gene"
  param$projectId = sub("\\/.*", "", input$getColumn("BAM")[1]) ## project id is needed for the IGV link
    
  gff = ezLoadFeatures(param)
  if (!is.null(gff) && nrow(gff) == 0){
    writeErrorReport(htmlFile, param=param, error=paste("No features found in given feature file:<br>", 
                                                                   param$ezRef["refFeatureFile"]))
    return("Error")
  }
  result = computeBamStats(input, htmlFile, param, gff)
  return(result)
}

##' @template app-template
##' @templateVar method ezMethodRnaBamStats(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppRnaBamStats <-
  setRefClass("EzAppRnaBamStats",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodRnaBamStats
                  name <<- "EzAppRnaBamStats"
                  appDefaults <<- rbind(posErrorRates=ezFrame(Type="logical",	DefaultValue="TRUE",	Description="compute position specific error rates?"),
                                        dupRadar=ezFrame(Type="logical",	DefaultValue="TRUE",	Description="run dupradar"),
                                    fragSizeMax=ezFrame(Type="integer",  DefaultValue=500,	Description="maximum fragment size to plot in fragment size distribution"),
                                    writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"),
                                    ignoreDup=ezFrame(Type="logical", DefaultValue="NA", Description="should marked duplicates be ignored?"),
                                    skipCountQc=ezFrame(Type="logical", DefaultValue=FALSE, Description="should we skip the count QC as part of the report"))
                }
              )
  )

##' @title Computes the BAM statistics
##' @description Computes the BAM statistics.
##' @template input-template
##' @template htmlFile-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refBuild}{ a character containing the path of the reference build.}
##'   \item{seqNames}{ the sequence names to select.}
##'   \item{posErrorRates}{ a logical indicating whether to call \code{getPosErrorFromBam()}.}
##'   \item{saveImage}{ a logical indicating whether to save the image.}
##' }
##' @param gff an annotation data.frame in gtf or gff format.
##' @param resultList a list of results.
##' @template roxygen-template
##' @seealso \code{\link{plotBamStat}}
computeBamStats = function(input, htmlFile, param, gff, resultList=NULL){

  require("GenomicAlignments")
  require("S4Vectors")
  samples = input$getNames()
  files = input$getFullPaths("BAM")
  dataset = input$meta
  
  ## get the RNA_repeats if available
  refRepeatsFeat = file.path(GENOMES_ROOT, param$ezRef["refBuild"], 
                             "Repeats/RNA_repeats.gff")
  if (file.exists(refRepeatsFeat)){
    repeatsGff = ezReadGff(refRepeatsFeat)
  } else {
    repeatsGff = NULL
  }
  
  if (ezIsSpecified(param$seqNames)){
    gff = gff[gff$seqid %in% param$seqNames, ]
    if(!is.null(repeatsGff)){
      repeatsGff = repeatsGff[repeatsGff$seqid %in% param$seqNames, ]
    }
  }
  if (is.null(resultList)){
    resultList = list()
    for (sm in samples){
      message(sm)
      resultList[[sm]] = getStatsFromBam(param, files[sm], sm, gff=gff, 
                                         repeatsGff=repeatsGff, 
                                         nReads=dataset[sm, "Read Count"])
      if (isError(resultList[[sm]])){
        writeErrorReport(htmlFile, param=param, error=resultList[[sm]]$error)
        return()
      }
      print(gc())
    }
    seqLengths = resultList[[1]]$seqLengths
    
    if (is.null(param$posErrorRates) || param$posErrorRates == TRUE){
      errorRates = ezMclapply(files, getPosErrorFromBam, param, 
                              mc.preschedule=FALSE,
                              mc.cores = min(length(files), ezThreads()))
      for (sm in samples){
        resultList[[sm]][["ErrorRates"]] = errorRates[[sm]]
      }
      rm(errorRates)
      gc()
    }
    ## do the analysis from package rseqc
    junctionsResults = ezMclapply(files, getJunctionPlotsFromBam, param, 
                                  mc.preschedule=FALSE,
                                  mc.cores = min(length(files), ezThreads()))
    for(sm in samples){
      resultList[[sm]][["Junction"]] = junctionsResults[[sm]]
    }
    rm(junctionsResults)
    gc()
    
    ## do Assessment of duplication rates from package dupRadar
    if (is.null(param$dupRadar) || param$dupRadar == TRUE){
      dupRateResults <- ezMclapply(files, getDupRateFromBam, param, 
                                   mc.preschedule = FALSE,
                                   mc.cores = min(length(files), ezThreads()))
      for(sm in samples){
        resultList[[sm]][["dupRate"]] = dupRateResults[[sm]]
      }
      rm(dupRateResults)
      gc()
    }
  }
  
  if (param$saveImage){
    save(resultList, file="resultList.RData")
  }
  
  ## debug
  #save(dataset, param, file="dataParam.rda")
  
  #plotBamStat(resultList, dataset, param, htmlFile)
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "RNABamStats.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="RNABamStats.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  rm(resultList)
  gc()
  return("Success")
}

##' @describeIn computeBamStats Gets the error positions from the BAM file.
getPosErrorFromBam = function(bamFile, param){
  require("GenomicRanges", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  require("bitops", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  job = ezJobStart(paste("position error:", bamFile))
  ## heuristic to search for the chromosome with most reads
  seqLengths = ezBamSeqLengths(bamFile)
  seqLengths = seqLengths[seqLengths >= mean(seqLengths)]
  if (ezIsSpecified(param$seqNames)){
    seqLengths = seqLengths[intersect(param$seqNames, names(seqLengths))]
  }
  seqCounts = sapply(names(seqLengths), function(sn){
    ezScanBam(bamFile, seqname = sn, isDuplicate=!param$ignoreDup, countOnly=TRUE)$records[1]
  })
  chromSel = names(seqLengths)[which.max(seqCounts)]
  fai = fasta.index(param$ezRef["refFastaFile"])
  targetGenome = readDNAStringSet(fai[match(chromSel, sub(" .*", "", fai$desc)), ])[[1]]  ## the sub clips potentially present description terms from the read id
  
  result = list()
  ezWriteElapsed(job, "read reference genome done")
  #targetGenome = referenceGenome[[match(chromSel, sub(" .*", "", names(referenceGenome)))]]
  what = c("strand", "cigar", "seq", "rname", "pos")
  if (param$paired){
    reads = ezScanBam(bamFile, seqname=chromSel, isFirstMateRead=TRUE, isSecondMateRead=FALSE, isUnmappedQuery=FALSE, isDuplicate=!param$ignoreDup, what=what)
    result[[paste(chromSel, "Position Stats of First Read")]] = ezPosSpecErrorRate(reads, targetGenome)
    reads = ezScanBam(bamFile, seqname=chromSel, isFirstMateRead=FALSE, isSecondMateRead=TRUE, isUnmappedQuery=FALSE, isDuplicate=!param$ignoreDup, what=what)
    result[[paste(chromSel, "Position Stats of Second Read")]] = ezPosSpecErrorRate(reads, targetGenome)
  } else {
    reads = ezScanBam(bamFile, seqname=chromSel, isUnmappedQuery=FALSE, isDuplicate=!param$ignoreDup, what=what)
    result[[paste(chromSel, "Position Stats")]] = ezPosSpecErrorRate(reads, targetGenome)
  }
  ezWriteElapsed(job, "Position Error Rate done")
  rm(reads)
  rm(targetGenome)
  gc()
  return(result)
}

##' @describeIn computeBamStats Calculates the specific error rates for \code{getPosErrorFromBam()}.
ezPosSpecErrorRate = function(bam, ReferenceGenome, nMaxReads=100000){
  require("GenomicRanges", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  require("stringr", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  require("Hmisc", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  ## remove the reads containing the gaps, insertions, deletions
  hasGap = grepl("N|I|D", bam$cigar)
  readLength <- width(bam$seq)
  isOutOfRange = bam$pos + readLength - 1 > length(ReferenceGenome) | bam$pos < readLength ## this is very conservative; needed because there might be clipped bases in the beginning
  if (any(isOutOfRange)){
    ezWrite("#reads out of range: ", sum(isOutOfRange))
    ezWrite("last pos: ", length(ReferenceGenome))
    idx = which(isOutOfRange)
    idx = idx[1:min(10, length(idx))]
    badAlignments = data.frame(pos=bam$pos[idx], cigar=bam$cigar[idx], 
                               start=bam$pos[idx], width=readLength[idx])
    print(badAlignments)
  }
  indexKeep = which(!hasGap & !isOutOfRange)
  if (length(indexKeep) > nMaxReads){
    indexKeep = sample(indexKeep, size=nMaxReads, replace=FALSE)
  }
  ezWrite("#alignments: ", length(bam$cigar),
          " #valid alignments: ", sum(!hasGap & !isOutOfRange), 
          " #used:", length(indexKeep))
  if (length(indexKeep) == 0){
    return(NULL)
  }
  for (nm in setdiff(names(bam), "tag")){
    ## do the replacement in place in order to save memory
    bam[[nm]] = bam[[nm]][indexKeep]
  }
  ## treat the tag separately
  for (tagName in names(bam$tag)){
    bam$tag[[tagName]] = bam$tag[[tagName]][indexKeep]
  }
  
  ## adjust the start POS according to H and/or S
  tempCigar = str_extract(bam$cigar, "^(\\d+H)?(\\d+S)?\\d+M")
  ## get the number of H at the beginning
  clipCigar = str_extract(tempCigar, "^\\d+H")
  noOfH = as.integer(sub("H", "", clipCigar))
  noOfH[is.na(noOfH)] = 0
  ## get the number of S at the beginning
  clipCigar = str_extract(tempCigar, "\\d+S")
  noOfS = as.integer(sub("S", "", clipCigar))
  noOfS[is.na(noOfS)] = 0
  nBeginClipped = noOfH + noOfS
  bam$pos = bam$pos - nBeginClipped
  
  ## add X to the begin and end of SEQ
  Xbegin = makeNstr("X", noOfH)
  tempCigar = str_extract(bam$cigar, "(\\d+S|\\d+H)$")
  clipCigar = str_extract(tempCigar, "\\d+H$")
  noOfH = as.integer(sub("H", "", clipCigar))
  noOfH[is.na(noOfH)] = 0
  Xend = makeNstr("X", noOfH)
  clipCigar = str_extract(tempCigar, "\\d+S")
  noOfS = as.integer(sub("S", "", clipCigar))
  noOfS[is.na(noOfS)] = 0
  nEndClipped = noOfH + noOfS
  bam$seq = paste0(Xbegin, bam$seq, Xend)
  
  seqChar = strsplit(bam$seq,"")
  readLength <- lengths(seqChar)
  ## build the reference views object
  maxLength = quantile(readLength, 0.95)
  if (maxLength < max(readLength)){
    readLength[readLength > maxLength] = maxLength
    seqChar = mapply(function(x, l){x[1:l]}, seqChar, readLength)
  }
  ReferenceViews = Views(ReferenceGenome, start=bam$pos, width=readLength)
  referenceChar = strsplit(as.character(ReferenceViews), "")
  
  # assuming we have unique read length and set it to the maximal read length here.
  nEndTrimmed = maxLength - readLength
  trimmedMatrix = mapply(function(readLength, nEndTrimmed){rep(c(FALSE, TRUE), c(readLength, nEndTrimmed))}, readLength, nEndTrimmed, SIMPLIFY=FALSE)
  ## build a clippedMatrix to record the clipped character
  nNormal = readLength - nBeginClipped - nEndClipped
  clippedMatrix = mapply(function(nBeginClipped, nNormal, nEndClipped, nEndTrimmed){rep(c(TRUE, FALSE, TRUE, FALSE), c(nBeginClipped, nNormal, nEndClipped, nEndTrimmed))}, nBeginClipped, nNormal, nEndClipped, nEndTrimmed, SIMPLIFY=FALSE)
  
  
  matchMatrix = mapply("==", referenceChar, seqChar, SIMPLIFY=FALSE)
  indexNeg = which(bam$strand == "-")
  clippedMatrix[indexNeg] = lapply(clippedMatrix[indexNeg], rev)
  ###  trimmedMatrix[indexNeg] = lapply(trimmedMatrix[indexNeg], rev) ## trimmed matrix must not be reverted!!!
  matchMatrix[indexNeg] = lapply(matchMatrix[indexNeg], rev)
  lengthPadding = maxLength - readLength
  matchMatrix = sapply(matchMatrix, function(x){x[1:maxLength]})
  clippedMatrix = matrix(unlist(clippedMatrix), ncol=length(clippedMatrix))
  clippedRate = rowMeans(clippedMatrix, na.rm=TRUE)
  trimmedMatrix = matrix(unlist(trimmedMatrix), ncol=length(trimmedMatrix))
  trimmedRate = rowMeans(trimmedMatrix, na.rm=TRUE)
  ## To distinguish error rate and clipped rate, remove the clipped from the mismatch, make it as match
  matchMatrix[clippedMatrix] = NA
  errorRate = 1 - rowMeans(matchMatrix, na.rm=TRUE)
  names(errorRate) = 1:nrow(matchMatrix)
  names(clippedRate) = 1:nrow(matchMatrix)
  names(trimmedRate) = 1:nrow(matchMatrix)
  return(list(trimmedRate=trimmedRate, clippedRate=clippedRate, 
              errorRate=errorRate))
}

##' @describeIn computeBamStats Gets the result statistics from the BAM file.
getStatsFromBam = function(param, bamFile, sm, gff=NULL, repeatsGff=NULL, 
                           nReads=NA){
  require("bitops", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  job = ezJobStart(paste("count", bamFile))
  seqLengths = ezBamSeqLengths(bamFile)  
  if (ezIsSpecified(param$seqNames)){
    seqLengths = seqLengths[param$seqNames]
  }
  if (is.null(param$splitByChrom) || param$splitByChrom){
    result = getStatsFromBamParallel(seqLengths, param, bamFile, sm, 
                                     gff, repeatsGff, mc.cores=param$cores, 
                                     nReads=nReads)
  } else {
    result = getStatsFromBamSingleChrom(NULL, param, bamFile, sm, gff, 
                                        repeatsGff)
  }
  gc()
  ## TODO: this getBamMultiMatching should be moved to computeBamStats
  result$multiMatchInFileTable = getBamMultiMatching(param, bamFile, nReads)
  
  transcriptCov = result$transcriptCov
  transcriptCovRleList <- RleList(transcriptCov)
  # transcriptLengthCov = sapply(transcriptCov, function(x){sum(x>0)}) slow!
  transcriptLengthCov = sum(transcriptCovRleList > 0)
  # transcriptLengthTotal = sapply(transcriptCov, length) slow!
  transcriptLengthTotal <- elementNROWS(transcriptCovRleList)
  percentCovered = transcriptLengthCov / transcriptLengthTotal * 100
  percentCoveredHist = hist(percentCovered, 
                            breaks=c(0, seq(from=0, to=100, by=2)), plot=FALSE)
  result$TranscriptsCovered = percentCoveredHist
  
  ## Do the genebody_coverage
  #sampledTranscriptCov = sapply(transcriptCov, ## RleList will be slow. use list
  #                              function(x){as.integer(x[round(seq(1, length(x), length.out=101))])})
  sampledTranscriptCov <- ezMclapply(transcriptCov, ## RleList will be slow. use list
                                function(x){as.integer(x[round(seq(1, length(x), length.out=101))])},
                                mc.preschedule=TRUE, mc.cores=param$cores)
  sampledTranscriptCov <- do.call(cbind, sampledTranscriptCov)
  
  trUse = colSums(sampledTranscriptCov) > 0
  sampledTranscriptCov = sampledTranscriptCov[ , trUse, drop=FALSE]
  trLength = transcriptLengthTotal[trUse]
  lengthClasses = ezCut(trLength, breaks=c(399, 1000, 1200, 4000), 
                        labels=c("less than 400nt", "400 to 999nt", 
                                 "1000 to 1200nt", "1201 to 4000", "above 4000nt"))
  genebody_coverage = list()
  for (lc in levels(lengthClasses)){
    isInLc = lengthClasses == lc
    if (sum(isInLc) > 50){
      ltc = sampledTranscriptCov[ , isInLc, drop=FALSE]
      avgCov = colMeans(ltc)
      relativeCov = ezScaleColumns(ltc, 1/colSums(ltc)) ## normalize so that every transcripts adds the same weight
      avgCovQuant = unique(quantile(avgCov, c(0.25, 0.75)))
      if (length(avgCovQuant) == 2){
        covClasses = ezCut(avgCov, breaks=avgCovQuant,
                           labels=c("low expressed", "medium expressed", 
                                    "high expressed"))
        genebody_coverage[[lc]] = list()
        for (cc in levels(covClasses)){
          genebody_coverage[[lc]][[cc]] = rowMeans(relativeCov[ , covClasses == cc, drop=FALSE ])
        }
      }
    }
  }
  result$genebody_coverage = genebody_coverage
  result$transcriptCov = NULL
  gc()
  return(result)
}

##' @describeIn computeBamStats Gets parallel by chromosome statistics for \code{getStatsFromBam()} if the logical \code{param$splitByChrom} is true.
getStatsFromBamParallel = function(seqLengths, param, bamFile, sm, 
                                   gff=NULL, repeatsGff=NULL, 
                                   mc.cores=ezThreads(), nReads=NA){
  seqNames = names(sort(seqLengths, decreasing=TRUE))
  names(seqNames) = seqNames
  if (!is.na(nReads)){
    ## heuristic: reduce the number of cores so that we have at least 0.25GB RAM per chromosome per Million Reads in the total bam file
    #reduce the number of threads in case of
    maxCores = ceiling(param$ram / (nReads / 1e6) * 4)
    if (maxCores < mc.cores){
      message("too many reads --reducing the number of cores: ", maxCores)
      mc.cores = maxCores
    }
  }
  chromResults = ezMclapply(seqNames, getStatsFromBamSingleChrom, param, 
                            bamFile, sm, gff, repeatsGff,
                            mc.preschedule=FALSE, mc.cores=mc.cores)
  if (param$saveImage){
    save(chromResults, file=paste0(sm, "-chromResults.RData"))
  }
  gc()
  result = list()
  # merge the segmentCountHist
  idx = which(sapply(chromResults, function(x){!is.null(x$segmentCountHist)}))
  if (length(idx) > 0){
    for (i in idx){
      if (i == idx[1]){
        sch = chromResults[[i]]$segmentCountHist
        counts = sch$counts
      } else {
        counts = counts + chromResults[[i]]$segmentCountHist$counts
      }
    }
    sch$counts = counts
    sch$density = counts / sum(counts)
    result$segmentCountHist = sch
  }
  
  # merge the fragSizeHist
  idx = which(sapply(chromResults, function(x){!is.null(x$fragSizeHist)}))
  if (length(idx) > 0){
    for (i in idx){
      if (i == idx[1]){
        sch = chromResults[[i]]$fragSizeHist
        counts = sch$counts
      } else {
        counts = counts + chromResults[[i]]$fragSizeHist$counts
      }
    }
    sch$counts = counts
    sch$density = counts / sum(counts)
    result$fragSizeHist = sch
  }
  # merge the multiMatchTargetTypeCounts
  temp = data.frame(count=integer(0), width=integer(0))
  for (i in 1:length(chromResults)){
    newResult = chromResults[[i]]$multiMatchTargetTypeCounts
    ## extend the temp data frame if needed
    additionalRows = setdiff(rownames(newResult), rownames(temp))
    if (length(additionalRows) > 0){
      temp[additionalRows, ] = 0
    }
    temp[rownames(newResult), ] = temp[rownames(newResult), ] + newResult
  }
  if (any(is.na(temp))){
    message("na counts: ", sum(is.na(temp)))
    temp[is.na(temp)] = 0
  }
  
  #mytotal = colSums(temp[names(seqLengths), ])
  #result$multiMatchTargetTypeCounts = rbind(total=mytotal, temp)
  tempNamesOrdered = intersect(c(setdiff(rownames(temp), seqNames), seqNames), 
                               rownames(temp))
  result$multiMatchTargetTypeCounts = temp[tempNamesOrdered, ,drop=FALSE]
  rm(temp)
  gc()
  geneCountList = lapply(chromResults, function(x){x$geneCounts})
  geneCounts = do.call("c", geneCountList)
  targetNames = unlist(lapply(geneCountList, names))
  if (any(duplicated(targetNames))){
    geneCounts = tapply(geneCounts, targetNames, FUN=sum)
  } else {
    names(geneCounts) = targetNames ## undo the prepending of the list names
  }
  result$geneCounts = geneCounts
  result$seqLengths = seqLengths
  
  ## Merge the TranscriptsCovered results
  #transcriptCov <- list()
  transcriptCov <- vector("list", length(chromResults))
  i <- 1
  for(chrom in names(chromResults)){
    #transcriptCov = c(transcriptCov, chromResults[[chrom]]$transcriptCov)
    transcriptCov[[i]] <- chromResults[[chrom]]$transcriptCov
    i <- i + 1
  }
  #transcriptCov <- unlist(lapply(chromResults, "[[", "transcriptCov"))
  # list object
  result$transcriptCov = unlist(transcriptCov)
  
  return(result)
}

##' @describeIn computeBamStats Gets the statistics of a single chromosome for \code{getStatsFromBam()}.
getStatsFromBamSingleChrom = function(chrom, param, bamFile, sm, 
                                      gff=NULL, repeatsGff=NULL){
  require("bitops", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  message("Processing chr ", ifelse(is.null(chrom), "all", chrom))
  
  result = list()
  seqLengths = ezBamSeqLengths(bamFile)
  if (param$paired){
    reads = ezReadPairedAlignments(bamFile, seqname=chrom, 
                                   keepUnpaired=param$keepUnpaired,
                                   minMapQuality=param$minMapQuality,
                                   keepMultiHits=param$keepMultiHits)
  } else {
    reads = ezReadGappedAlignments(bamFile, seqname=chrom, 
                                   minMapQuality=param$minMapQuality,
                                   keepMultiHits=param$keepMultiHits)	
  }
  if (isError(reads)){
    return(reads)
  }
  #isMultiHit = my.duplicated(names(reads), mode = "all")
  if (length(reads) > 0){
    result$segmentCountHist = intHist(njunc(reads)+1, range=c(0.5, 10.5),
                                      plot=FALSE)
  }
  #ezWriteElapsed(job, "segment hist done")
  gc()
  
  if(param$paired && length(reads) > 0){
    pairedNames = ezScanBam(bamFile, seqname=chrom, 
                            isFirstMateRead=TRUE, isSecondMateRead=FALSE,
                            isProperPair=TRUE, isUnmappedQuery=FALSE,
                            isDuplicate=!param$ignoreDup, what="qname")$qname
    use = names(reads) %in% pairedNames
    result$fragSizeHist = intHist(width(reads)[use], 
                                  range=c(-0.5, param$fragSizeMax + 0.5), plot=FALSE)
    rm(pairedNames)
    gc()
  }
  
  result$multiMatchTargetTypeCounts = getTargetTypeCounts(param, gff, reads,
                                                          seqid=chrom, repeatsGff)
  result$geneCounts = getFeatureCounts(chrom=chrom, gff=gff[gff$type == "exon", ],
                                       reads=reads, param=param)
  #ezWriteElapsed(job, "multi match type counts")
  #if (any(isMultiHit)){
  #	result$uniqueMatchTargetTypeCounts = getTargetTypeCounts(param, gff, reads[!isMultiHit])
  #	ezWriteElapsed(job, "unique match type counts")
  #}
  gc()
  ## Do transcripts covered
  result$transcriptCov = getTranscriptCoverage(chrom, gff, reads, strandMode=param$strandMode)
  gc()
  rm(reads)
  #rm(isMultiHit)
  gc()
  return(result)
}

##' @describeIn computeBamStats Gets the counts of the target types for \code{getStatsFromBam()}.
getTargetTypeCounts = function(param, gff, rr, seqid=NULL, repeatsGff=NULL){
  require(data.table)
  if (class(rr) == "GRangesList"){
    #sn = unlist(sn, use.names=FALSE)[sn@partitioning@end]
    stop("GRangesList not supported")
  }
  seqNames = names(seqlengths(rr))
  if (!is.null(seqid)){
    stopifnot(length(seqid) == 1)
    seqNames = intersect(seqNames, seqid)
    readRefIsValid = as.character(seqnames(rr)) %in% seqNames
    if (!all(readRefIsValid)){
      rr = rr[readRefIsValid]
    }
  }
  #effWidth = sum(as.numeric(seqlengths(rr))) * ifelse(param$isStranded, 2, 1)
  effWidth = (seqlengths(rr) * ifelse(param$strandMode == "both", 1, 2))[seqNames]
  result = data.frame(count=0, width=effWidth, row.names=seqNames)
  #readCounts = table(as.character(seqnames(rr)))
  readCounts <- table(seqnames(rr)) # It also works with different order
  
  result[seqNames, "count"] = readCounts[seqNames]
  result$count[is.na(result$count)] = 0 ## if a chromosome has no reads the value would be na
  
  #hasAnyHit = rep(FALSE, length(rr))
  hasAnyHit <- logical(length(rr))
  repeatsRanges = NULL
  gffRanges = NULL
  
  ## the reads in the repeatsGff
  if(!is.null(repeatsGff)){
    repeatsGff = repeatsGff[repeatsGff$seqid %in% seqNames, ]
    if(nrow(repeatsGff) > 0){
      repeatsGff$strand = fixStrand(repeatsGff$strand, param$strandMode)
      repeatsRanges = gffToRanges(repeatsGff)
      classFam = ezGffAttributeField(repeatsGff$attributes, field="repClass")
      repFamily = ezGffAttributeField(repeatsGff$attributes, field="repFamily")
      use = repFamily != classFam
      classFam[use] = paste(classFam[use], repFamily[use], sep="--")
      for(type in unique(classFam)){
        use = classFam == type
        targetRanges = gffToRanges(repeatsGff[use, ])
        hitsTarget = overlapsAny(rr, targetRanges, minoverlap=10)
        result[type, ] = c(sum(hitsTarget), sum(width(reduce(targetRanges))))
        hasAnyHit = hasAnyHit | hitsTarget
      }
    }
  }
  
  if (!is.null(gff)){
    gff = gff[gff$seqid %in% seqNames, ]
    if (nrow(gff) > 0){
      gff$strand = fixStrand(gff$strand, param$strandMode)
      gffRanges = gffToRanges(gff)
      ensemblTypes = getEnsemblTypes(gff)
      if (!is.null(ensemblTypes)){
        # for (type in unique(ensemblTypes)){
        #   targetRanges = gffRanges[ensemblTypes == type]
        #   hitsTarget = overlapsAny(rr, targetRanges, minoverlap=10)
        #   result[type, ] = c(sum(hitsTarget), sum(width(reduce(targetRanges))))
        #   hasAnyHit = hasAnyHit | hitsTarget
        # }
        ## The following code is much faster than the loop above.
        ## The loop: 537.492 seconds;
        ## New implementation: 184 seconds
        hits <- findOverlaps(rr, gffRanges, minoverlap=10)
        hitsByType <- data.table(queryHits=queryHits(hits), 
                                 ensemblTypes=ensemblTypes[subjectHits(hits)])
        hitsByType <- unique(hitsByType)
        countsByType <- hitsByType[ , .N, by=ensemblTypes]
        if(nrow(countsByType) == 0L){
          ## When there is no hit for all reads
          countsByType <- data.table(ensemblTypes=unique(ensemblTypes),
                                     N=0)
        }
        ## some ensemblType has not hits and not included in findOverlaps.
        ## Add them with N=0
        missingTypes <- setdiff(unique(ensemblTypes), countsByType$ensemblTypes)
        if(length(missingTypes) != 0L){
          countsByType <- rbind(countsByType,
                                data.table(ensemblTypes=missingTypes, N=0))
        }
        
        widthByType <- sum(width(reduce(GenomicRanges::split(gffRanges, 
                                                             ensemblTypes))))
        
        hasAnyHit[hitsByType$queryHits] <- TRUE
        result[countsByType$ensemblTypes, ] <- cbind(countsByType$N,
                                                     widthByType[countsByType$ensemblTypes])
        
        isMsg = ensemblTypes == "protein_coding" & gff$type == "exon"
        msgRanges = gffGroupToRanges(gff[isMsg, ], gff$transcript_id[isMsg], 
                                     skipTransSpliced = TRUE)
        targetExonRanges = gffRanges[isMsg]
      } else {
        rootTypes = setdiff(gff$type, c("intron", "exon"))
        for (type in rootTypes){
          use = gff$type == type
          targetRanges = gffRanges[gff$type == type]
          hitsTarget = overlapsAny(rr, targetRanges, minoverlap=10)
          result[type, ] = c(sum(hitsTarget), sum(width(reduce(targetRanges))))
          hasAnyHit = hasAnyHit | hitsTarget
        }
        isExon = gff$type == "exon"
        msgRanges = gffGroupToRanges(gff[isExon, ], gff$transcript_id[isExon], 
                                     skipTransSpliced = TRUE)
        targetExonRanges = gffRanges[isExon]
      }
      ## Add seqlengths for msgRanges because of potential out-of-bound by flank
      ## TODO: now we use seqlengths from rr, which is not ideal.
      ## gff as data.frame has no seqlengths information.
      seqlengths(msgRanges) <- seqlengths(rr)[names(seqlengths(msgRanges))]
      
      ## check additionally for intron/exon/prom
      hitsTranscript = overlapsAny(rr, msgRanges, minoverlap=10)
      hasAnyHit = hasAnyHit | hitsTranscript  
      mRnaWidth = sum(width(reduce(msgRanges)))
      hitsTargetExons = overlapsAny(rr[hitsTranscript], 
                                    targetExonRanges, minoverlap=10)
      result["mRNA Exons", ] = c(sum(hitsTargetExons), 
                                 sum(width(reduce(targetExonRanges))))
      result["mRNA Introns", ] = c(sum(!hitsTargetExons), 
                                   mRnaWidth - result["mRNA Exons", "width"])
      ## suppressWarnings for out-of-bound ranges.
      promRanges = trim(suppressWarnings(flank(msgRanges, 2000)))
      hitsTargetProms = overlapsAny(rr, promRanges, minoverlap=10)
      result["mRNA Promoter 2kb", ] = c(sum(hitsTargetProms), 
                                        sum(width(reduce(promRanges))))
      hasAnyHit = hasAnyHit | hitsTargetProms
      downRanges = trim(suppressWarnings(flank(msgRanges, 2000, start=FALSE)))
      hitsTargetDown = overlapsAny(rr, downRanges, minoverlap=10)
      result["mRNA Downstream 2kb", ] = c(sum(hitsTargetDown), 
                                          sum(width(reduce(promRanges))))
      hasAnyHit = hasAnyHit | hitsTargetDown
      gffRanges = c(gffRanges, promRanges, downRanges)
    }
  }
  allRanges = GRanges()
  if (!is.null(gffRanges)){
    allRanges = c(allRanges, gffRanges)
  }
  if (!is.null(repeatsGff)){
    allRanges = c(allRanges, repeatsRanges)
  }
  if (length(allRanges) > 0){
    annotatedWidth = sum(width(reduce(allRanges)))
  } else {
    annotatedWidth = 0
  }
  #result["unannotated", ] = c(sum(!hasAnyHit), result[seqNames, "width"] - annotatedWidth)
  result["unannotated", ] = c(sum(!hasAnyHit), 
                              sum(result[seqNames, "width"]) - annotatedWidth)
  result["total", ] = colSums(result[seqNames, ], na.rm=TRUE)
  result = result[c("total", setdiff(rownames(result), "total")), ]
  return(result)
}

##' @describeIn computeBamStats Gets the junction results from the BAM file.
getJunctionPlotsFromBam = function(bamFile, param){
  pngFiles = list()
  ## do the junction annotation
  outputJunction = paste0("junction-", Sys.getpid())
  bed = getReferenceFeaturesBed(param)
  stopifnot(!is.null(bed))
  # junction_annotation.py is in RSeQC package available through Dev/Python2
  cmd = paste("junction_annotation.py", "--mapq=1", "-i", bamFile, 
              "-o", outputJunction, "-r", bed)
  res = ezSystem(cmd, stopOnFailure=FALSE)
  junctionFile = paste0(outputJunction, ".junction.xls")
  if (res == 0 && length(readLines(junctionFile)) > 1){
    juncsTable = read.table(junctionFile, header=TRUE)
    junctions = table(juncsTable$annotation) / length(juncsTable$annotation) *100
    foo = rep(juncsTable$annotation, juncsTable$read_count)
    events = table(foo) / length(foo) * 100
    pngFiles[["splice_events"]] = events
    pngFiles[["splice_junction"]] = junctions
    
    ## do the junction_saturation
    id = paste(juncsTable[["chrom"]], juncsTable[["intron_st.0.based."]], 
               juncsTable[["intron_end.1.based."]])
    juncReads = rep(id, juncsTable[["read_count"]])
    juncTypes = rep(juncsTable[["annotation"]], juncsTable[["read_count"]])
    juncTypeSet = c("annotated", "complete_novel", "partial_novel")
    nSim = 10
    quantiles = seq(0.05, 1, by=0.05)
    juncCounts = array(0, dim=c(length(juncTypeSet), length(quantiles), nSim), 
                       dimnames=list(types=juncTypeSet, 
                                     quantiles=as.character(quantiles), 
                                     sim=1:nSim))
    for(n in 1:nSim){
      idx = sample(1:length(juncReads), replace=FALSE)
      for(q in 1:length(quantiles)){
        idxUse = idx[1:round(quantiles[q] * length(idx))]
        cts = tapply(juncReads[idxUse], juncTypes[idxUse], 
                     function(x){length(unique(x))})
        juncCounts[names(cts), q, n] = as.vector(cts)
      }
    }
    ## TODO why can they be NA????
    juncCounts[is.na(juncCounts)] = 0
    juncCountMeans = apply(juncCounts, c(1,2), mean)
    juncCountMeans[is.na(juncCountMeans)] = 0
    junctionSaturations = list()
    junctionSaturations[["all junctions"]] = colSums(juncCountMeans[c("annotated", "complete_novel", "partial_novel"), ])
    junctionSaturations[["known junctions"]] = colSums(juncCountMeans[c("annotated"), , drop=FALSE])
    junctionSaturations[["novel junctions"]] = colSums(juncCountMeans[c("complete_novel", "partial_novel"), , drop=FALSE])
    pngFiles[["junctionSaturation"]] = junctionSaturations
  }
  
  ## do the cleaning
  file.remove(list.files(path=".", pattern=paste0(outputJunction, ".+")))
  return(pngFiles)
}


### classify reads according to where they match
### - assign reads to transcripts
### - split transcript reads into intronic/exonic
### - use remaining reads and assign to other root features
### - use remaining reads to assign to promoters
### - use remaining reads to downstream
### - use remaining reads to intergenic

## behaviour for strand-preserving library prep:
## - genome target size is 2 * genome size
## - match with ignoreStrand=false

## behaviour for non-strand-preserving library prep:
## - genome target size is 1 * genome size


# getTypePercentTable = function(resultList, name, minPercentage=1){
#   tbl = data.frame(row.names=rownames(resultList[[1]][[name]]))
#   for (sm in names(resultList)){
#     counts = resultList[[sm]][[name]]
#     tbl[sm] = round(counts[ rownames(tbl), "count"] / counts["total", "count"] * 100, digits=3)
#   }
#   useRow = apply(tbl, 1, max) > minPercentage & rownames(tbl) != "total"
#   return(tbl[useRow , ])
# }

getDupRateFromBam <- function(bamFile, param=NULL, gtfFn, 
                              stranded=c("both", "sense", "antisense"), 
                              paired=FALSE, #dupremover=c("picard", "bamutil"),
                              threads=1){
  if(!is.null(param)){
    gtfFn <- param$ezRef@refFeatureFile
    stranded <- param$strandMode
    paired <- param$paired
    threads <- ceiling(param$cores / 2)
  }
  require(dupRadar)
  #dupremover <- match.arg(dupremover)
  #setEnvironments(dupremover)
  
  ## Make the duplicates in bamFile
  inputBam <- paste(Sys.getpid(), basename(bamFile), sep="-")
  ### The bamFile may not be writable.
  file.symlink(from=bamFile, to=inputBam)
  
  ## intermediate files
  # picardMetricsFn <- gsub("\\.bam$", "_picard_metrics.txt", inputBam) # picard
  bamDuprmFn <- gsub("\\.bam$", "_duprm.bam", inputBam)
  # bamutilLogFn <- paste0(bamDuprmFn, ".log") # bamutil
  on.exit(file.remove(c(inputBam, bamDuprmFn)))#, picardMetricsFn, bamutilLogFn)))
  # 
  # if(dupremover == "bamutil"){
  #   dupremoverDir <- dirname(Sys.which("bam"))
  # }else if(dupremover == "picard"){
  #   dupremoverDir <- dirname(Sys.getenv("Picard_jar"))
  #   ## when modue load Tools/Picard.
  # }
  
  # bamDuprm <- markDuplicates(dupremover=dupremover,
  #                            bam=inputBam,
  #                            out=bamDuprmFn,
  #                            path=dupremoverDir,
  #                            rminput=FALSE)
  dupBam(inBam=inputBam, outBam=bamDuprmFn, operation="mark",
         cores=threads)
  ## Duplication rate analysis
  dm <- analyzeDuprates(bam=bamDuprmFn, gtf=gtfFn,
                        stranded=switch(stranded, "both"=0, "sense"=1, 
                                        "antisense"=2, 
                                        stop("unsupported strand mode: ", 
                                             stranded)), 
                        paired=paired, threads=threads)
  return(dm)
}

