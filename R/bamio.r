###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets sequence names
##' @description Gets sequence names from a bam file.
##' @template bamFile-template
##' @param sizeSorting a character. if equal to "decreasing" the names will be sorted decreasingly.
##' @template roxygen-template
##' @return Returns a character vector of sequence names.
##' @examples 
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' ezBamSeqLengths(bamFile)
##' ezBamSeqNames(bamFile)
ezBamSeqNames = function(bamFile, sizeSorting="decreasing"){
  seqLengths = ezBamSeqLengths(bamFile)
  if (sizeSorting == "decreasing"){
    seqLengths = sort(seqLengths, decreasing=TRUE)
  }
  return(names(seqLengths))
}

##' @describeIn ezBamSeqNames Gets the sequence lengths from a bam file.
ezBamSeqLengths = function(bamFile){
  require(Rsamtools)
  return(scanBamHeader(bamFile)[[1]]$targets)
}

ezSortIndexBam = function(inBam, bam, ram=2, removeBam=TRUE, cores=2,
                          method=c("samtools", "sambamba", "picard")){
  method <- match.arg(method)
  
  maxMem = paste0(as.integer(floor(ram * 0.7/cores*1000)), "M") 
  ## use only 70% --> 30% safety margin before crash
  setEnvironments(method)
  
  if(method == "sambamba"){
    cmd <- paste(method, "sort", "-l 9", "-m", maxMem, "-t", cores, inBam,
                 "-o", bam)
    ezSystem(cmd)
    cmd = paste(method, "index", "-t", cores, bam)
    ezSystem(cmd)
  }else if(method == "samtools"){
    cmd = paste(method, "sort", "-l 9", "-m", maxMem, "-@", cores, inBam, 
                "-o", bam)
    ezSystem(cmd)
    cmd = paste(method, "index", bam)
    ezSystem(cmd)
  }else{
    cmd <- paste("java -Djava.io.tmpdir=. -jar", 
                 Sys.getenv("Picard_jar"), "SortSam",
                 paste0("I=", inBam),
                 paste0("O=", bam),
                 "SORT_ORDER=coordinate")
    ezSystem(cmd)
    cmd <- paste("java -Djava.io.tmpdir=. -jar", 
                 Sys.getenv("Picard_jar"), "BuildBamIndex",
                 paste0("I=", bam),
                 paste0("OUTPUT=", bam, ".bai")
                 )
    ezSystem(cmd)
  }
  if (removeBam){
    file.remove(inBam)
  }
}

##' @title Scans a bam file
##' @description Scans a bam file with the option to select only a part of the output.
##' @template bamFile-template
##' @param seqname an optional character vector to keep only the specified sequence names in the returned reads.
##' @param start an optional integer vector to limit the start position of the range of returned reads.
##' @param end an optional integer vector to limit the start position of the range of returned reads.
##' @param strand a character specifying which strand to use (+, - or *).
##' @param tag passed further to \code{ScanBamParam()}.
##' @param what passed further to \code{ScanBamParam()}.
##' @param isFirstMateRead passed further to \code{scanBamFlag()}.
##' @param isSecondMateRead passed further to \code{scanBamFlag()}.
##' @param isUnmappedQuery passed further to \code{scanBamFlag()}.
##' @param isProperPair passed further to \code{scanBamFlag()}.
##' @param countOnly a logical specifying whether only counts should be returned or the whole scan of the bam file.
##' @template roxygen-template
##' @seealso \code{\link[Rsamtools]{scanBam}}
##' @seealso \code{\link[Rsamtools]{ScanBamParam}}
##' @return Returns a list representing a scanned bam file.
##' @examples 
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' bam = ezScanBam(bamFile)
ezScanBam = function(bamFile, seqname=NULL, start=NULL, end=NULL, strand="*",
                      tag=character(0), what=scanBamWhat(),
                      isFirstMateRead=NA, isSecondMateRead=NA, isUnmappedQuery=FALSE, isProperPair=NA,
                     isSecondaryAlignment=NA, isDuplicate = NA, countOnly=FALSE){
  ## initialize the parameters for scanBam
  require(Rsamtools)
  param = ScanBamParam(what=what, tag=tag)
  
  ### build and set the flag filter
  isMinusStrand = switch(strand, "-"=TRUE, "+"=FALSE, NA)
  bamFlag(param) = scanBamFlag(isMinusStrand=isMinusStrand, isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead,
                               isUnmappedQuery=isUnmappedQuery, isProperPair=isProperPair, isSecondaryAlignment=isSecondaryAlignment, isDuplicate=isDuplicate)
  
  ### limit the returned reads by chromosome and start/end if these are provided
  if (!is.null(seqname)){
    seqLengths = ezBamSeqLengths(bamFile)
    if (! seqname %in% names(seqLengths)){
      return(list(error=paste("Specified seqname: ", seqname, " not in chromosome list of ", bamFile, ":", paste(names(seqLengths), collapse=" "))))
    }
    if (is.null(start))  start = 1
    if (is.null(end))  end = seqLengths[seqname]
    bamWhich(param) = GRanges(seqnames=seqname, ranges=IRanges(start=start, end=end))
  }
  if (countOnly){
    return(countBam(bamFile, param=param))    
  } else {
    return(scanBam(bamFile, param=param)[[1]])
  }
}

##' @title Converts from bam to bigwig
##' @description Converts a bam file into a bigwig file.
##' @template bamFile-template
##' @param bigwigPrefix a character representing a prefix to add to the bigwig file names. 
##' @param param a list of parameter to extract the \code{strandMode} from.
##' @param paired a logical indicating whether the samples are paired.
##' @template roxygen-template
##' @seealso \code{\link[Rsamtools]{scanBamHeader}}
##' @examples 
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' ezBam2bigwig(bamFile, bigwigPrefix="bg", param=ezParam(), paired=TRUE)
ezBam2bigwig = function(bamFile, bigwigPrefix, param=NULL, paired=NULL){
  
  ## from bam to wig
  require(Rsamtools)
  sh = scanBamHeader(bamFile)
  seqLengths = sh[[1]]$targets
  chromSize = data.frame(seqLengths, row.names=names(seqLengths))
  chromSizeFile = paste0("chromSize-",Sys.getpid(), ".txt")
  write.table(chromSize, sep="\t", file=chromSizeFile, col.names=FALSE, quote=FALSE)
  outputWig = paste0("outputWig-", Sys.getpid())
  strand_rule = NULL
  if(paired){
    strand_rule = switch(param$strandMode,
                         both="",
                         sense="--strand='1++,1--,2+-,2-+'",
                         antisense="--strand='1+-,1-+,2++,2--'")
  } else {
    strand_rule = switch(param$strandMode,
                         both="",
                         sense="--strand='++,--'",
                         antisense="--strand='+-,-+'")    
  }
  # bam2wig.py is in RSeQC module
  cmd = paste("bam2wig.py", "-i", bamFile, "-s", chromSizeFile, "-o", outputWig, strand_rule)
  ezSystem(cmd)
  
  ## from wig to bigwig
  wigFiles = list.files(path=".", pattern=paste0(outputWig, ".*\\.wig$"))
  bigwigFiles = gsub(outputWig, bigwigPrefix, wigFiles)
  bigwigFiles = gsub(".wig$", ".bw", bigwigFiles)
  for(i in 1:length(wigFiles)){
    cmd = paste("wigToBigWig", wigFiles[i], chromSizeFile, bigwigFiles[i])
    try(ezSystem(cmd))
  }
  
  file.remove(wigFiles)
  file.remove(chromSizeFile)
}

##' @title Reads gapped alignments from bam
##' @description Reads gapped alignments from a bam file.
##' @template bamFile-template
##' @param seqname an optional character vector to keep only the specified sequence names in the returned reads.
##' @param start an optional integer vector to limit the start position of the range of returned reads.
##' @param end an optional integer vector to limit the start position of the range of returned reads.
##' @param strand a character specifying which strand to use (+, - or *).
##' @param tag passed further to \code{ScanBamParam()}.
##' @param what passed further to \code{ScanBamParam()}.
##' @param use.names passed further to \code{readGAlignments()}.
##' @param isFirstMateRead passed further to \code{scanBamFlag()}.
##' @param isSecondMateRead passed further to \code{scanBamFlag()}.
##' @param isUnmappedQuery passed further to \code{scanBamFlag()}.
##' @param isProperPair passed further to \code{scanBamFlag()}.
##' @param minMapQuality an integer specifying the minimal mapping quality to use.
##' @param keepMultiHits a logical specifying whether to keep multiple hits.
##' @template roxygen-template
##' @seealso \code{\link[Rsamtools]{scanBam}}
##' @seealso \code{\link[Rsamtools]{ScanBamParam}}
##' @seealso \code{\link[GenomicAlignments]{readGAlignments}}
##' @return Returns an object of the class GAlignments.
##' @examples 
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' ga = ezReadGappedAlignments(bamFile)
ezReadGappedAlignments = function(bamFile, seqname=NULL, start=NULL, end=NULL, strand="*",
                                  tag=character(0), what=character(0), use.names=TRUE,
                                  isFirstMateRead=NA, isSecondMateRead=NA, 
                                  isUnmappedQuery=FALSE, isProperPair=NA,
                                  isSecondaryAlignment = NA,
                                  minMapQuality=0, 
                                  tagFilter=list(),
                                  keepMultiHits=TRUE){
  ## initialize the parameters for scanBam
  require(Rsamtools)
  require(GenomicAlignments)
  if (minMapQuality > 0){
    what = union(what, "mapq")
  }
  if (!keepMultiHits){
    tag = union(tag, c("NH", "IH"))
  }
  param = ScanBamParam(what=what, tag=tag, tagFilter = tagFilter)
  
  ### build and set the flag filter
  isMinusStrand = switch(strand, "-"=TRUE, "+"=FALSE, NA)
  bamFlag(param) = scanBamFlag(isMinusStrand=isMinusStrand, 
                               isFirstMateRead=isFirstMateRead, 
                               isSecondMateRead=isSecondMateRead,
                               isUnmappedQuery=isUnmappedQuery, 
                               isProperPair=isProperPair, 
                               isSecondaryAlignment=isSecondaryAlignment)
  
  ### limit the returned reads by chromosome and start/end if these are provided
  if (!is.null(seqname)){
    seqLengths = ezBamSeqLengths(bamFile[1])
    if (! seqname %in% names(seqLengths)){
      return(list(error=paste("Specified seqname: ", seqname, " not in chromosome list of ", bamFile, ":", paste(names(seqLengths), collapse=" "))))
    }
    if (is.null(start))  start = 1
    if (is.null(end))  end = seqLengths[seqname]
    bamWhich(param) = GRanges(seqnames=seqname, ranges=IRanges(start=start, end=end))
  }
  
  if (length(bamFile) > 1){
    ga = Reduce(c, lapply(bamFile, readGAlignments, use.names=use.names, param=param))
  } else {
    ga = readGAlignments(bamFile, use.names=use.names, param=param)
  }
  if (minMapQuality > 0){
    mc = mcols(ga)
    use = mc$mapq >= minMapQuality
    use[is.na(use)] = TRUE   ## accept also mapping quality of 255
    ga = ga[use]
  }
  if (!keepMultiHits){
    mc = mcols(ga)
    numberOfHits = NULL
    if (all(!is.na(mc$NH))){
      numberOfHits = mc$NH
    } else {
      if (all(!is.na(mc$IH))){
        numberOfHits = mc$IH
      }
    }
    if (is.null(numberOfHits)){
      stop(paste("multi-hit information not available in bam file:", bamFile))
    }
    ga = ga[numberOfHits == 1]
  }
  return(ga)
}

##' @title Reads paired alignments from bam
##' @description Reads paired alignments from a bam file.
##' @template bamFile-template
##' @param seqname an optional character vector to keep only the specified sequence names.
##' @param start passed further to \code{ezReadGappedAlignments()}.
##' @param end passed further to \code{ezReadGappedAlignments()}.
##' @param strand a character specifying which strand to use (+, - or *).
##' @param tag passed further to \code{ezReadGappedAlignments()}.
##' @param keepUnpaired a character indicating which aligned but unpaired reads to keep. Possible options: "none", "first", "second", "both".
##' @param fillGap a character to fill into the gaps produced by \code{ezMergeLeftRightAlignments()}.
##' @param minMapQuality passed further to \code{ezReadGappedAlignments()}.
##' @param keepMultiHits passed further to \code{ezReadGappedAlignments()}.
##' @param gaLeft the left gapped alignments.
##' @param gaRight the right gapped alignments.
##' @template roxygen-template
##' @seealso \code{\link{ezReadGappedAlignments}}
##' @return Returns an object of the class GAlignments.
##' @examples 
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' ezReadPairedAlignments(bamFile)
ezReadPairedAlignments = function(bamFile, seqname=NULL, start=NULL, end=NULL, strand="*",
                                   tag=c("NH"), keepUnpaired="both", fillGap="N", minMapQuality=0, keepMultiHits=TRUE){
  
  require(Rsamtools)
  require(GenomicAlignments)
  .loadPairedSingleChrom = function(chrom){
    return(ezReadPairedAlignments(bamFile, seqname=chrom, start=start, end=end, strand=strand,
                                  tag=tag, keepUnpaired=keepUnpaired, minMapQuality=minMapQuality, keepMultiHits=keepMultiHits))
  }
  
  if (is.null(seqname)){
    ## if we load the entire genome we run it multi-threaded and combine then the different
    ## GappedAlignment objects
    job = ezJobStart(paste("readPaired", bamFile))
    seqNames = ezBamSeqNames(bamFile)
    gaList = ezMclapply(seqNames, .loadPairedSingleChrom,  mc.preschedule=FALSE)
    names(gaList) = seqNames
    ezWriteElapsed(job, "chromosomes loaded")
    ga = GAlignments()
    for (sn in seqNames){
      ga = c(ga, gaList[[sn]])
      gaList[[sn]] = NULL
      gc()
    }
    ezWriteElapsed(job, "chromosomes merged")
    return(ga)
  }
  
  seqLengths = ezBamSeqLengths(bamFile)  
  gaAll = GAlignments(seqnames=Rle(factor( levels=names(seqLengths))), seqlengths=seqLengths)
  if (strand %in% c("+", "*")){
    gaLeft = ezReadGappedAlignments(bamFile, seqname=seqname, start=start, end=end, strand="+", 
                                    isFirstMateRead=TRUE, isSecondMateRead=FALSE, isProperPair=TRUE, what=c("mpos"), tag=tag,
                                    minMapQuality=minMapQuality, keepMultiHits=keepMultiHits)
    #ezWriteElapsed(job, "loaded paired gaLeft for positive strand")
    gaRight = ezReadGappedAlignments(bamFile, seqname=seqname, start=start, end=end, strand="-", 
                                     isFirstMateRead=FALSE, isSecondMateRead=TRUE, isProperPair=TRUE, what=c("mpos"), tag=tag,
                                     minMapQuality=minMapQuality, keepMultiHits=keepMultiHits)
    #ezWriteElapsed(job, "loaded paired gaRight for positive strand")
    ## strand of gaRight will be wrong but we will never use it!!
    gaPositive <- .ezMergeLeftRightAlignments(gaLeft, gaRight, fillGap=fillGap)
    rm(gaLeft, gaRight)
    gc()
    gaAll = c(gaAll, gaPositive)
    rm(gaPositive)
    gc()
  }
  if (strand %in% c("-", "*")){
    gaLeft = ezReadGappedAlignments(bamFile, seqname=seqname, start=start, end=end, strand="+", 
                                    isFirstMateRead=FALSE, isSecondMateRead=TRUE, isProperPair=TRUE, what=c("mpos"), tag=tag,
                                    minMapQuality=minMapQuality, keepMultiHits=keepMultiHits)
    #ezWriteElapsed(job, "loaded paired gaLeft for negative strand")
    gaRight = ezReadGappedAlignments(bamFile, seqname=seqname, start=start, end=end, strand="-", 
                                     isFirstMateRead=TRUE, isSecondMateRead=FALSE, isProperPair=TRUE, what=c("mpos"), tag=tag,
                                     minMapQuality=minMapQuality, keepMultiHits=keepMultiHits)
    #ezWriteElapsed(job, "loaded paired gaRight for negative strand")
    gaNegative = .ezMergeLeftRightAlignments(gaLeft, gaRight)
    if (!is.null(gaNegative)){
      strand(gaNegative) = "-"
    }
    rm(gaLeft, gaRight)
    gc()
    gaAll = c(gaAll, gaNegative)
    rm(gaNegative)
    gc()
  }
  ## by default we load also read that were aligned but not paired with their mate
  stopifnot(keepUnpaired %in% c("none", "first", "second", "both"))
  if (keepUnpaired != "none"){
    if (keepUnpaired == "both"){
      isFirstMateRead = NA
      isSecondMateRead = NA
    }
    if (keepUnpaired == "first"){
      isFirstMateRead = TRUE
      isSecondMateRead = FALSE
    }
    if (keepUnpaired == "second"){
      isFirstMateRead = FALSE
      isSecondMateRead = TRUE
    }
    gaSingle = ezReadGappedAlignments(bamFile, seqname=seqname, start=start, end=end, strand=strand,
                                      isProperPair=FALSE,
                                      isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead,
                                      what=c("flag"), tag=tag,
                                      minMapQuality=minMapQuality, keepMultiHits=keepMultiHits)
    isAlsoPaired = names(gaSingle) %in% names(gaAll)
    if (any(isAlsoPaired)){
      gaSingle = gaSingle[!isAlsoPaired]
    }
    doFlipStrand = bamFlagTest(values(gaSingle)$flag, "isSecondMateRead")
    if (any(doFlipStrand)){
      myStrand = strand(gaSingle)
      setToPos = myStrand == "-" & doFlipStrand
      setToNeg = myStrand == "+" & doFlipStrand
      if (any(setToPos)){
        myStrand[setToPos] = "+"
      }
      if (any(setToNeg)){
        myStrand[setToNeg] = "-"
      }
      strand(gaSingle) = myStrand
    }
    values(gaSingle)$flag = NULL
    gaAll = c(gaAll, gaSingle)
    rm(gaSingle)
    gc()
  }
  return(gaAll)
}

##' @describeIn ezReadPairedAlignments Merges left and right alignments and fills gaps with \code{fillGap}.
.ezMergeLeftRightAlignments <- function(gaLeft, gaRight, fillGap="N"){
  ## include the multi-paired-hits
  idx <- match(paste(names(gaLeft), values(gaLeft)[["mpos"]]), 
               paste(names(gaRight), start(gaRight)))
  isNaIdx = is.na(idx)
  gaAll = gaLeft[isNaIdx]
  values(gaAll)$mpos = NULL
  if (any(isNaIdx)){
    message(paste("found invalid pairs:", sum(isNaIdx), "/", length(isNaIdx)))
    gaLeft = gaLeft[!isNaIdx]
    gaRight = gaRight[idx[!isNaIdx]]
  } else {
    gaRight <- gaRight[idx]
  }
  values(gaRight)$mpos = NULL
  values(gaLeft)$mpos = NULL
  
  
  dist = start(gaRight) - end(gaLeft) -1    
  ## some useful index
  # gaRight is totally to the right of gaLeft
  hasGap = dist >= 0
  # has overhang
  hasOverhang = end(gaRight) > end(gaLeft)
  ## check if there are insertions or deletions
  hasIndel = ezGrepl("D|I", cigar(gaRight))
  
  ## situation one: gap or zero distance between gaLeft and gaRight
  if (any(hasGap)){
    gapString = sub("^0N$", "", paste0(format(dist[hasGap], scientific=FALSE, trim=TRUE), fillGap))
    mergedCigar = paste0(cigar(gaLeft)[hasGap], gapString, cigar(gaRight)[hasGap])
    gaGap = GAlignments(seqnames=seqnames(gaLeft)[hasGap], pos=start(gaLeft)[hasGap], cigar=mergedCigar, strand=strand(gaLeft)[hasGap],
                        names=names(gaLeft)[hasGap])
    values(gaGap) = values(gaLeft)[hasGap, ,drop=FALSE]
    gaAll = c(gaGap, gaAll)
  }
  
  ## situation two: overlap and overhang no indel problems
  useNarrow = !hasGap & hasOverhang & !hasIndel
  if (any(useNarrow)){
    gaRightNarrowed = narrow(gaRight[useNarrow], -dist[useNarrow] + 1)
    distNarrowed = start(gaRightNarrowed) - end(gaLeft[useNarrow]) -1
    gapString = sub("^0N$", "", paste0(distNarrowed, "N"))
    mergedCigar = paste0(cigar(gaLeft)[useNarrow], gapString, cigar(gaRightNarrowed))
    gaNew = GAlignments(seqnames=seqnames(gaLeft)[useNarrow], pos=start(gaLeft)[useNarrow],
                        cigar=mergedCigar, strand=strand(gaLeft)[useNarrow],
                        names=names(gaLeft)[useNarrow])
    values(gaNew) = values(gaLeft)[useNarrow, ,drop=FALSE]
    gaAll = c(gaNew, gaAll)
  }
  useNarrowIndel = !hasGap & hasOverhang & hasIndel
  if (any(useNarrowIndel)){
    ## TODO: do the narrowing of the right in the presence of indel in gaRight
    ## Current: uses only the left as a workaround!!!!
    gaAll = c(gaLeft[useNarrowIndel], gaAll)
  }
  
  ## situation three: overlap but no overhang --> use left reads
  useLeft = !hasGap & !hasOverhang
  if (any(useLeft)){
    gaAll = c(gaLeft[useLeft], gaAll)
  }
  return(gaAll)
}


##' @title Gets bam multi matching
##' @description Gets bam multi matching and returns the result as an integer vector.
##' @param param a list of parameters:
##' \itemize{
##'   \item{paired}{ a logical indicating whether the samples are paired.}
##' }
##' @template bamFile-template
##' @param nReads an integer specifying the number of reads. These additional reads will be counted in "0"
##' @template roxygen-template
##' @return Returns a integer vector specifying the amount of matching reads.
##' @examples 
##' param = ezParam()
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' getBamMultiMatching(param, bamFile, nReads=10000)
getBamMultiMatching = function(param, bamFile, nReads=NULL){
  require(data.table)
  require(Rsamtools)
  
  # This data.table based implementation is the fastest and most elegant!
  # test file: /srv/gstore/projects/p2438/STAR_18564_2017-06-12--13-46-30/26EV_d3_A.bam
  ## Shell based is ugly and slow. Writing temps to disks. 1534.172 seconds
  ## table(table()) method is more tidy but slow. 1521.440 seconds + reading qname
  ## plyr::count is slow: 1770.432 seconds + reading qname
  ## rle(sort()) is slow: 1781.552 seconds + reading qname
  ## data.table: 4.036 seconds + reading qname
  
  #  job = ezJobStart(paste("bam multimatch", basename(bamFile)))
  paramBam = ScanBamParam(what="qname")
  bamFlag(paramBam) = scanBamFlag(isUnmappedQuery=FALSE)
  if(isTRUE(param$paired)){
  #     bamFlag(param) = switch(param$paired,
  #                       paired=scanBamFlag(isFirstMateRead=TRUE, isProperPair=TRUE, isUnmappedQuery=FALSE),
  #                       first=scanBamFlag(isFirstMateRead=TRUE, isUnmappedQuery=FALSE),
  #                       second=scanBamFlag(isSecondMateRead=TRUE, isUnmappedQuery=FALSE))
    bamFlag(paramBam) = scanBamFlag(isFirstMateRead=TRUE, isProperPair=TRUE, 
                                    isUnmappedQuery=FALSE)
  }
  bamReads = scanBam(bamFile, param=paramBam)[[1]]$qname
  bamReadsDT <- as.data.table(bamReads)
  bamReadsDTCount <- bamReadsDT[, .N, by=bamReads]
  colnames(bamReadsDTCount)[2] <- "qnameN"
  bamReadsDTCount2 <- bamReadsDTCount[, .N, by=qnameN]
  bamReadsDTCount2 <- bamReadsDTCount2[order(qnameN)]
  result <- setNames(bamReadsDTCount2$N,
                     bamReadsDTCount2$qnameN)
  #result = table(table(bamReads))
  #if (param$paired){
  #  flagOption = "-f 3 -F 132"
#     switch(param$pairedMode,
#                         paired="-f 3 -F 132",
#                         first="-F 132",
#                         second="-F 68")    
  #} else {
  #  flagOption = "-F 4"
  #}
  #countFile = paste0(Sys.getpid(), "-BamMultiMatch.txt")
  #cmd = paste("samtools", "view", flagOption, bamFile, "| cut -f1 | sort | uniq -c | cut -c 1-7 | sort | uniq -c >", countFile)
  ## direct usage of regular expression, slow
  #cmd = paste("samtools", "view", flagOption, bamFile, "| cut -f1 | sort | uniq -c | grep -o "[[:blank:]]*[[:digit:]]\+[ ]" | sort | uniq -c >", countFile)
  ## trim the leading spaces first and then cut the first field, faster
  ## set the temp directory to be the current one, because /tmp may be too small
  #cmd = paste("samtools", "view", flagOption, bamFile, "| cut -f1 | sort --temporary-directory=. | uniq -c | sed -e \"s/^[ \t]*//\" | cut -f1 -d\" \" |sort | uniq -c >", countFile)
  
  #ezSystem(cmd)
  # temp = read.table(countFile)
  # tempOrdered = temp[order(temp[,2]), ]
  # result = tempOrdered[ ,1]
  # names(result) = tempOrdered[, 2]
  #file.remove(countFile)
  
  if (!is.null(nReads)){
    nReads = as.integer(nReads)
    if (!is.na(nReads)){
      nReadsUnmapped = nReads - sum(result)
      stopifnot(nReadsUnmapped >= 0)
      result = c("0"=nReadsUnmapped, result)
    }
  }
  return(result)
}

getBamLocally = function(src, toSam=FALSE){
  if (normalizePath(dirname(src)) %in% c(".", normalizePath(getwd()))){
    return(src)
  }
  target = basename(src)
  stopifnot(target != src)
  if (toSam){
    target = sub(".bam$", ".sam", target)
    if (file.exists(target)) file.remove(target)
    cmd = paste("samtools", "view -h", "-o", target, src)
    ezSystem(cmd)
  } else {
    file.copy(from=src, to=target)
    file.copy(from=paste0(src, ".bai"),
              to=paste0(target, ".bai"))
  }
  return(target)
}
