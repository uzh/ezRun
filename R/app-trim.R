###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## order of the trimming
## - adapter
## - quality window (SLIDING_WINDOW)
## - avg qual
## - minlen
## - fixed
##' @title Trims input reads
##' @description Trims input reads. There are several options to influence trimming with parameters.
##' @template input-template
##' @param output an object of the class EzDataset or NA. If it is NA, it will be copied from the input.
##' @param param a list of parameters:
##' \itemize{
##'   \item{paired}{ a logical specifying whether the samples have paired ends.}
##'   \item{subsampleReads}{ an integer specifying how many subsamples there are. This will call \code{ezMethodSubsampleReads()} if > 1.}
##'   \item{trimAdapter}{ a logical specifying whether to use a trim adapter.}
##'   \item{minTailQuality}{ an integer specifying the minimal tail quality to accept. Only used if > 0.}
##'   \item{minAvgQuality}{ an integer specifying the minimal average quality to accept. Only used if > 0.}
##'   \item{minReadLength}{ an integer specifying the minimal read length to accept.}
##'   \item{dataRoot}{ a character specifying the path of the data root to get the full column paths from.}
##' }
##' @template roxygen-template
##' @return Returns the output after trimming as an object of the class EzDataset.
ezMethodTrim = function(input=NA, output=NA, param=NA){
  
  ## if output is not defined, set it!
  ## TODO: gives warning message if outpu is S4 class
  if (is.na(output)){
    output = input$copy()
    output$setColumn("Read1", paste0(input$getNames(), "-trimmed-R1.fastq"))
    if (param$paired){
      output$setColumn("Read2", paste0(input$getNames(), "-trimmed-R2.fastq"))
    } else {
      if ("Read2" %in% input$colNames){
        output$setColumn("Read2", NULL)
      }
    }
  }
  
  ## if there are multiple samples loop through them
  if (input$getLength() > 1){
    for (nm in input$getNames()){
      ezMethodTrim(input$copy()$subset(nm), output$copy()$subset(nm), param)
      ## NOTE: potential risk, temp files might be overwritten
    }
    return(output)
  }
  
  ## now we deal only with one sample!
  param$trimSeedMismatches = 1
  param$trimPalindromClipThresh = 20
  param$trimSimpleClipThresh = 7
  param$trimMinAdaptLength = 5
  param$trimKeepBothReads = "true"
  param$trimQualWindowWidth = 4
  adaptFile = "/usr/local/ngseq/src/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa"
  
  if (param$subsampleReads > 1 || param$nReads > 0){
    input = ezMethodSubsampleReads(input=input, param=param)
  }
  
  if (param$trimAdapter){
    trimAdaptOpt =  paste("ILLUMINACLIP", adaptFile, param$trimSeedMismatches, param$trimPalindromClipThresh, 
                          param$trimSimpleClipThresh, param$trimMinAdaptLength, param$trimKeepBothReads, sep=":")
  } else {
    trimAdaptOpt = ""
  }
  if (param$minTailQuality > 0){
    tailQualOpt = paste("SLIDINGWINDOW", param$trimQualWindowWidth, param$minTailQuality, sep=":")
  } else {
    tailQualOpt = ""
  }
  if (param$minAvgQuality > 0){
    minAvgQualOpt = paste("AVGQUAL", param$minAvgQuality, sep=":")
  } else {
    minAvgQualOpt = ""
  }
  r1TmpFile = "trimmed-R1.fastq"
  r2TmpFile = "trimmed-R2.fastq"
  if (any(c(trimAdaptOpt, tailQualOpt, minAvgQualOpt) != "")){
    if (param$paired){
      method = "PE"
      readOpts = paste(
        input$getFullPaths(param, "Read1"),
        input$getFullPaths(param, "Read2"),
        r1TmpFile,
        "unpaired-R1.fastq",
        r2TmpFile,
        "unpaired-R2.fastq")
    } else {
      method = "SE"
      readOpts = paste(
        input$getFullPaths(param, "Read1"), r1TmpFile)
    }
    cmd = paste(TRIMMOMATIC, method,
                "-threads", min(ezThreads(), 8),
                #"-trimlog", paste0(input$getNames(), "-trimmomatic.log"),
                readOpts, trimAdaptOpt, tailQualOpt, minAvgQualOpt,
                #               paste("HEADCROP", param$trimLeft, sep=":"),
                #               paste("CROP", param$trimRight, sep=":"),
                paste("MINLEN", param$minReadLength, sep=":"),
                "> trimmomatic.out 2> trimmomatic.err")
    ezSystem(cmd)
  } else {
    ezSystem(paste("gunzip -c", input$getFullPaths(param, "Read1"), ">", r1TmpFile))
    if (param$paired){
      ezSystem(paste("gunzip -c", input$getFullPaths(param, "Read2"), ">", r2TmpFile))      
    }
  }
  
  ## the flexbar call is done separately because we do want to make sure that fixed trimming is done on top of adapter trimming
  ## this is needed for STAR to be able to pair the reads properly
  if (param$trimLeft > 0 || param$trimRight > 0){
    if (param$paired){
      pairedOpt = paste("-p", r2TmpFile)
    } else {
      pairedOpt = ""
    }
    cmd = paste(FLEXBAR,
                "--threads", min(ezThreads(), 8),
                "-r", r1TmpFile,
                pairedOpt,
                "--format", "i1.8",
                "-u", 20, ##### max uncalled bases
                "--pre-trim-left", param$trimLeft,
                "--pre-trim-right", param$trimRight,
                "--target", "flexbar",
                "> flexbar.out 2> flexbar.err")
    ezSystem(cmd)
    if (param$paired) {
      file.remove(r1TmpFile)
      file.remove(r2TmpFile)
      r1TmpFile = "flexbar_1.fastq"
      r2TmpFile = "flexbar_2.fastq"
    } else {
      file.remove(r1TmpFile)
      r1TmpFile = "flexbar.fastq"
    }
  }
  

  if (param$paired){
    ezSystem(paste("mv", r1TmpFile, basename(output$getColumn("Read1"))))
    ezSystem(paste("mv", r2TmpFile, basename(output$getColumn("Read2"))))
  } else {
    ezSystem(paste("mv", r1TmpFile, basename(output$getColumn("Read1"))))
  }
  return(output)
}

##' @describeIn ezMethodTrim Gets the subsample files, calls \code{ezSubsampleFastq()} on them and returns the output, which is an object of the class EzDataset.
ezMethodSubsampleReads = function(input=NA, output=NA, param=NA){
  ## TODO: gives warning message if outpu is S4 class
  if (is.na(output)){
    output = input$copy()
    subsampleFiles = sub(".fastq.*", "-subsample.fastq", basename(input$getColumn("Read1")))
    output$setColumn(name="Read1", values = file.path(getwd(), subsampleFiles))
    if (param$paired){
      subsampleFiles = sub(".fastq.*", "-subsample.fastq", basename(input$getColumn("Read2")))
      output$setColumn(name="Read2", values = file.path(getwd(), subsampleFiles))      
    }
  }
  if (param$nReads > 0){
    totalReads = input$getColumn("Read Count")
    subsampleFactor = shrinkToRange(totalReads / param$nReads, c(1, Inf))
    if (param$subsampleReads > 1){
      message("subsampleReads setting will be overwritten by nReads parameter")
    }
  } else {
    subsampleFactor = param$subsampleReads
  }
  newReadCounts = ezSubsampleFastq(input$getFullPaths(param, "Read1"), output$getColumn("Read1"), subsampleFactor = subsampleFactor)
  output$setColumn("Read Count", newReadCounts)
  if (param$paired){
    ezSubsampleFastq(input$getFullPaths(param, "Read2"), output$getColumn("Read2"), subsampleFactor = subsampleFactor)
  }
  return(output)
}

##' @describeIn ezMethodTrim Performs the fastq for the subsamples using the package ShortRead.
##' @examples 
##'  inputFile = system.file(package = "ezRun", "extdata/yeast_10k/wt_1_R1.fastq.gz")
##'  subsampledFile = "sub_R1.fastq"
##'  ezSubsampleFastq(inputFile, subsampledFile, subsampleFactor = 5)
ezSubsampleFastq = function(full, sub, subsampleFactor=NA, nYield=1e5, overwrite=FALSE){
  stopifnot(full != sub)
  stopifnot(length(full) == length(sub))
  if (any(file.exists(sub))){
    filesToRemove = sub[file.exists(sub)]
    if (!overwrite){
      stop("files do exist: ", filesToRemove)
    }
    warning("removing first: ", filesToRemove)
    file.remove(filesToRemove)
  }
  requireNamespace("ShortRead")
  nms = names(full)
  if (is.null(nms)){
    nms = full
  }
  nReadsVector = integer()
  for (i in seq_along(full)){ 
    nReads = 0
    fqs = FastqStreamer(full[i], n = nYield) 
    idx = seq(from=1, to=nYield, by=subsampleFactor)
    while(length(x <- yield(fqs))){
      if (length(x) >= nYield){
        writeFastq(x[idx], file=sub[i], mode="a", full=F, compress=F)
        nReads = nReads + length(idx)
      } else {
        writeFastq(x[idx[idx<length(x)]], file=sub[i], mode="a", full=F, compress=F)
        nReads = nReads + sum(idx<length(x))
      }
    }
    close(fqs)
    nReadsVector[nms[i]] = nReads
  }
  return(nReadsVector)
}
