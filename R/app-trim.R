###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Trims input reads
##' @description Trims input reads. The trimming happens in the following order:
##' \itemize{
##'   \item{adapter}
##'   \item{quality window}
##'   \item{avg quality}
##'   \item{minimum length}
##'   \item{fixed trimming}
##' }
##' This means if you specify a fixed trimming this comes on top of the adaptive trimming
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

  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")){
    output = input$copy()
    output$setColumn("Read1", paste0(getwd(), "/", input$getNames(), "-trimmed-R1.fastq"))
    if (param$paired){
      output$setColumn("Read2", paste0(getwd(), "/", input$getNames(), "-trimmed-R2.fastq"))
    } else {
      if ("Read2" %in% input$colNames){
        output$setColumn("Read2", NULL)
      }
    }
    output$dataRoot = NULL
  }


  ## if there are multiple samples loop through them
  if (input$getLength() > 1){
    for (nm in input$getNames()){
      ezSystem(paste0('touch ', nm, '_preprocessing.log'))
      ezMethodTrim(input$subset(nm), output$subset(nm), param)
      ## NOTE: potential risk, temp files might be overwritten
    }
    return(output)
  }

  ## now we deal only with one sample!


  ## make a local copy of the dataset and check the md5sum
  if (param$paired){
    reads = c("Read1", "Read2")
  } else{
    reads = "Read1"
  }
  for (rds in reads){
    readFileIn = input$getFullPaths(rds)
    ezSystem(paste("cp -n", readFileIn, "."))
    input$setColumn(rds, basename(readFileIn))
    if (Sys.info()["user"] == "trxcopy") { ## only run the check for the user trxcopy!!!
      md5Local = ezSystem(paste("md5sum", basename(readFileIn)), intern = TRUE)
      md5Local = sub(" .*", "", md5Local)
      md5File = file.path(dirname(readFileIn), "md5.txt")
      md5Remote = NA
      if (file.exists(md5File)){
        md5Set = ezRead.table(md5File)
        md5Remote = md5Set[basename(readFileIn), 1]
      }
      if (is.na(md5Remote)){
        md5Remote = ezSystem(paste("ssh fgcz-s-022 md5sum", readFileIn), intern = TRUE)
        md5Remote = sub(" .*", "", md5Remote)
      }
      stopifnot(md5Local == md5Remote)
    }
  }
  input$dataRoot = NULL


  param$trimSeedMismatches = 1
  param$trimPalindromClipThresh = 20
  param$trimSimpleClipThresh = 7
  param$trimMinAdaptLength = 5
  param$trimKeepBothReads = "true"
  param$trimQualWindowWidth = 4

  if (param$subsampleReads > 1 || param$nReads > 0){
    input = ezMethodSubsampleReads(input=input, param=param)
  }

  if (param$trimAdapter){
    if (!is.null(input$meta$Adapter1) && !is.na(input$meta$Adapter1) && input$meta$Adapter1 != ""){
      adapter1 = DNAStringSet(input$meta$Adapter1)
      names(adapter1) = "GivenAdapter1"
    } else {
      adapter1 = DNAStringSet()
    }
    if (param$paired && !is.null(input$meta$Adapter2) && !is.na(input$meta$Adapter2) && input$meta$Adapter2 != ""){
      adapter2 = DNAStringSet(input$meta$Adapter2)
      names(adapter2) = "GivenAdapter2"
    } else {
      adapter2 = DNAStringSet()
    }
    # take only adapter from dataset and ignore the ones from TRIMMOMATIC_ADAPTERS
    if (!is.null(param$onlyAdapterFromDataset) && param$onlyAdapterFromDataset){
      adapters = c(adapter1, adapter2)
    } else {
      adapters = c(adapter1, adapter2, readDNAStringSet(TRIMMOMATIC_ADAPTERS))
    }
    adaptFile = "adapters.fa"
    writeXStringSet(adapters, adaptFile)
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
  if (any(c(trimAdaptOpt, tailQualOpt, minAvgQualOpt) != "") || param$minReadLength > 0){
    if (param$paired){
      method = "PE"
      readOpts = paste(
        input$getFullPaths("Read1"),
        input$getFullPaths("Read2"),
        r1TmpFile,
        "unpaired-R1.fastq",
        r2TmpFile,
        "unpaired-R2.fastq")
    } else {
      method = "SE"
      readOpts = paste(
        input$getFullPaths("Read1"), r1TmpFile)
    }
    cmd = paste(TRIMMOMATIC, method,
                "-threads", min(ezThreads(), 8), "-phred33", ## hardcode phred33 quality encoding
                #"-trimlog", paste0(input$getNames(), "-trimmomatic.log"),
                readOpts, trimAdaptOpt, tailQualOpt, minAvgQualOpt,
                #               paste("HEADCROP", param$trimLeft, sep=":"),
                #               paste("CROP", param$trimRight, sep=":"),
                paste("MINLEN", param$minReadLength, sep=":"),
                ">> trimmomatic.out 2>> trimmomatic.err")
    ezSystem(cmd)
    ## TRIMOMMATIC may throw exception but still return status 0
    exceptionCount = length(grep("Exception", readLines("trimmomatic.err")))
    stopifnot(exceptionCount == 0)
    cmd = paste0('cat trimmomatic.err >>',input$getNames(),'_preprocessing.log')
    ezSystem(cmd)
  } else {
    ezSystem(paste("gunzip -c", input$getFullPaths("Read1"), ">", r1TmpFile))
    if (param$paired){
      ezSystem(paste("gunzip -c", input$getFullPaths("Read2"), ">", r2TmpFile))
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
    if(param$minReadLength > 0){
      minReadLengthOpt = paste("--min-read-length", param$minReadLength)
    }
    else {
      minReadLengthOpt = ""
    }
    cmd = paste(FLEXBAR,
                "--threads", min(ezThreads(), 8),
                "-r", r1TmpFile,
                pairedOpt,
                "--format", "i1.8",
                "-u", 20, ##### max uncalled bases
                "--pre-trim-left", param$trimLeft,
                "--pre-trim-right", param$trimRight,
                minReadLengthOpt,
                "--target", "flexbar",
                "> flexbar.out 2> flexbar.err")
    ezSystem(cmd)
    cmd = paste0('cat flexbar.out >>',input$getNames(),'_preprocessing.log')
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

  ## filter by max read length
  if (!is.null(param$maxReadLength) && !is.na(as.integer(param$maxReadLength))){
    newFile = "lengthTrimmed_R1.fastq"
    maxLength = as.integer(param$maxReadLength)
    fqs = FastqStreamer(r1TmpFile, n=1e6)
    while(length(x <- yield(fqs))){
      writeFastq(x[width(x) <= maxLength], file = newFile, mode="a", compress=FALSE)
    }
    close(fqs)
    file.remove(r1TmpFile)
    r1TmpFile = newFile
    if (param$paired){
      newFile = "lengthTrimmed_R2.fastq"
      maxLength = as.integer(param$maxReadLength)
      fqs = FastqStreamer(r2TmpFile, n=1e6)
      while(length(x <- yield(fqs))){
        writeFastq(x[width(x) <= maxLength], file = newFile, mode="a")
      }
      close(fqs)
      file.remove(r2TmpFile)
      r2TmpFile = newFile
    }
  }

  ezSystem(paste("mv", r1TmpFile, basename(output$getColumn("Read1"))))
  if (param$paired){
    ezSystem(paste("mv", r2TmpFile, basename(output$getColumn("Read2"))))
  }
  return(output)
}

##' @title Subsample reads in a fastq dataset file
##' @description The subsampled reads are equally distributed across the original files
##' @param param a list of parameters where the following entries are used
##' \itemize{
##'   \item{dataRoot} the prefix of the file paths
##'   \item{nReads} the number of reads to keep; if given will be used to compute subsampleFactor; it is not guaranteed that the number of reads kept is exact
##'   \item{subsampleFactor} the factor by which subsampling has been done. if \code{nReads} is specified subsampleFactor will not be used
##'   \item{paired} whether these are paired-end reads
##' }
##' @examples
##' inputDatasetFile = system.file(package = "ezRun", "extdata/yeast_10k/dataset.tsv")
##' param = ezParam(list(dataRoot=system.file(package = "ezRun"), subsampleFactor=5))
##' input = EzDataset(file=inputDatasetFile, dataRoot=param$dataRoot)
##' xSubsampled = ezMethodSubsampleReads(input=input, param=param)
##' # NOTE: the subsampled files will not be gzip compressed!
ezMethodSubsampleReads = function(input=NA, output=NA, param=NA){
  if (!is(output, "EzDataset")){
    output = input$copy()
    subsampleFiles = sub(".fastq.*", "-subsample.fastq", basename(input$getColumn("Read1")))
    output$setColumn(name="Read1", values = file.path(getwd(), subsampleFiles))
    if (param$paired){
      subsampleFiles = sub(".fastq.*", "-subsample.fastq", basename(input$getColumn("Read2")))
      output$setColumn(name="Read2", values = file.path(getwd(), subsampleFiles))
    }
    output$dataRoot = NULL
  }
  if (param$nReads > 0){
    totalReads = as.numeric(input$getColumn("Read Count"))
    subsampleFactor = shrinkToRange(totalReads / param$nReads, c(1, Inf))
    if (param$subsampleReads > 1){
      message("subsampleReads setting will be overwritten by nReads parameter")
    }
  } else {
    subsampleFactor = param$subsampleReads
  }
  newReadCounts = ezSubsampleFastq(input$getFullPaths("Read1"), output$getColumn("Read1"), subsampleFactor = subsampleFactor)
  output$setColumn("Read Count", newReadCounts)
  if (param$paired){
    ezSubsampleFastq(input$getFullPaths("Read2"), output$getColumn("Read2"), subsampleFactor = subsampleFactor)
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
  require("ShortRead")
  nms = names(full)
  if (is.null(nms)){
    nms = full
  }
  nReadsVector = integer()
  for (i in seq_along(full)){
    nReads = 0
    fqs = FastqStreamer(full[i], n = nYield)
    idx = seq(from=1, to=nYield, by=subsampleFactor)
    tmpFile = sub("\\.fastq.*", "_temp.fastq", sub[i])
    while(length(x <- yield(fqs))){
      if (length(x) >= nYield){
        writeFastq(x[idx], file=tmpFile, mode="a", full=F, compress=F)
        nReads = nReads + length(idx)
      } else {
        writeFastq(x[idx[idx<length(x)]], file=tmpFile, mode="a", full=F, compress=F)
        nReads = nReads + sum(idx<length(x))
      }
    }
    close(fqs)
    nReadsVector[nms[i]] = nReads
    if (grepl(".gz$", sub[i])){
      ezSystem(paste("pigz -p 2 --best -c", tmpFile, ">", sub[i]))
      file.remove(tmpFile)
    } else {
      file.rename(tmpFile, sub[i])
    }
  }
  return(nReadsVector)
}
