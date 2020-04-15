###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Trims input reads using fastp (Chen et al. 2018)
##' @description Trims input reads. The trimming happens in the following order:
##' \itemize{
##'   \item{global trimming fron and tail}
##'   \item{quality window}
##'   \item{avg quality}
##'   \item{adapter}
##'   \item{trim polyG (by default, only for NextSeq and NovaSeq)}
##'   \item{fixed trimming}
##'   \item{length filtering}
##' }
##' This means if you specify a fixed trimming this comes on top of the adapter trimming

ezMethodFastpTrim = function(input=NA, output=NA, param=NA){
  # Input/Output Preparation
  ## stop early
  if (any(grepl("bam$", input$getFullPaths("Read1")))){
    stop("cannot process unmapped bam as input")
  }
  ## if output is not an EzDataset, set it! (when ezMethodFastpTrim is used inside another app)
  if (!is(output, "EzDataset")){
    output = input$copy()
    output$setColumn("Read1", paste0(getwd(), "/", input$getNames(), "-trimmed_R1.fastq.gz"))
    if (param$paired){
      output$setColumn("Read2", paste0(getwd(), "/", input$getNames(), "-trimmed_R2.fastq.gz"))
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
      ezMethodFastpTrim(input$subset(nm), output$subset(nm), param)
      ## NOTE: potential risk, temp files might be overwritten
    }
    return(output)
  }
  
  ## now we deal only with one sample!
  
  ## make a local copy of the dataset and check the md5sum
  if (!is.null(param$copyReadsLocally) && param$copyReadsLocally){
    input = copyReadsLocally(input, param)
  }
  
  if (param$subsampleReads > 1 || param$nReads > 0){
    input = ezMethodSubsampleReads(input=input, param=param)
    inputRead1 <- input$getFullPaths("Read1")
    on.exit(file.remove(inputRead1), add=TRUE)
    if(param$paired){
      inputRead2 <- input$getFullPaths("Read2")
      on.exit(file.remove(inputRead2), add=TRUE)
    }
  }
  
  # Prepare shell command
  ## binary
  fastpBin = 'fastp'
  ## input/output file names:
  r1TmpFile = "trimmed_R1.fastq.gz"
  if(param$paired){
    r2TmpFile = "trimmed_R2.fastq.gz"
    readsInOut = paste('--in1', input$getFullPaths("Read1"),
                       '--in2', input$getFullPaths("Read2"),
                       '--out1', r1TmpFile,
                       '--out2', r2TmpFile)
  }else{
    readsInOut = paste('--in1', input$getFullPaths("Read1"),
                       '--out1', r1TmpFile)
  }
  ## adapter trimming options
  if (param[['trimAdapter']]){
    # read1
    if (!is.null(input$meta$Adapter1) && !is.na(input$meta$Adapter1) && input$meta$Adapter1 != ""){
      adapter1 = DNAStringSet(input$meta$Adapter1)
      names(adapter1) = "GivenAdapter1"
    } else {
      adapter1 = DNAStringSet()
    }
    # read2 (if paired)
    if (param$paired && !is.null(input$meta$Adapter2) && !is.na(input$meta$Adapter2) && input$meta$Adapter2 != ""){
      adapter2 = DNAStringSet(input$meta$Adapter2)
      names(adapter2) = "GivenAdapter2"
    } else {
      adapter2 = DNAStringSet()
    }
    adaptFile = "adapters.fa"
    adapters = c(adapter1, adapter2)
    if (!is.null(param$onlyAdapterFromDataset) && param$onlyAdapterFromDataset){
      # take only adapter from dataset and ignore the ones from TRIMMOMATIC_ADAPTERS
      writeXStringSet(adapters, adaptFile)
    } else {
      file.copy(from=TRIMMOMATIC_ADAPTERS,
                to=adaptFile)
      writeXStringSet(adapters, adaptFile, append=TRUE)
    }

    trimAdapt = paste('--adapter_fasta', adaptFile)
    
  } else {
    ezSystem('touch adapters.fa')
    adaptFile = "adapters.fa"
    trimAdapt = "--disable_adapter_trimming"
  }

  ## paste command
  cmd = paste(fastpBin,
              readsInOut,
              # general options
              paste('--thread',param[['cores']]),
              # global trimming
              paste('--trim_front1',param[['trim_front1']]),
              paste('--trim_tail1',param[['trim_tail1']]),
              # quality-based trimming per read
              if(param[['cut_front']]) paste("--cut_front", "--cut_front_window_size", param[['cut_front_window_size']], "--cut_front_mean_quality", param[['cut_front_mean_quality']]), # like Trimmomatic's LEADING
              if(param[['cut_right']]) paste("--cut_right", "--cut_right_window_size", param[['cut_right_window_size']], "--cut_right_mean_quality", param[['cut_right_mean_quality']]), # like Trimmomatic's SLIDINGWINDOW
              if(param[['cut_tail']])  paste("--cut_tail", "--cut_tail_window_size", param[['cut_tail_window_size']], "--cut_tail_mean_quality", param[['cut_tail_mean_quality']]), # like Trimmomatic's TRAILING
              paste("--average_qual", param[['average_qual']]),
              # adapter trimming
              trimAdapt,
              # read length trimming
              paste('--max_len1',param[['max_len1']]),
              paste('--max_len2',param[['max_len2']]),
              # polyX
              paste("--disable_trim_poly_g","--trim_poly_x",
                    "--poly_x_min_len",param[['poly_x_min_len']]),
              # read length filtering
              paste('--length_required',param[['length_required']]),
              # compression output
              paste('--compression 4')
  )
  
  ## run
  if(ezIsSpecified(param[['cmdOptionsFastp']])){
    cmd = paste(cmd, param[['cmdOptionsFastp']])
  }
  ezSystem(paste0(cmd,' 2> fastp.err'))
  
  ## remove reports
  ezSystem("rm fastp.json fastp.html")
  
  ## rename adapters.fa (standalone) or not (within another app)
  if("Adapters" %in% output$colNames){
    renamedAdaptFile = paste0(input$getNames(),"_",adaptFile)
    ezSystem(paste("mv",adaptFile,renamedAdaptFile))
  }else{
    on.exit(file.remove(adaptFile), add=TRUE)
  }
  
  ## rename log
  ezSystem(paste0('mv fastp.err ',input$getNames(),'_preprocessing.log'))

  ## rename trimmed output
  ezSystem(paste("mv", r1TmpFile, basename(output$getColumn("Read1"))))
  if (param$paired){
    ezSystem(paste("mv", r2TmpFile, basename(output$getColumn("Read2"))))
  }
  
  return(output)
}

##' @title EzAppFastp app
##' @description fast read pre-processing.
##' @author Miquel Anglada Girotto
EzAppFastp =
  setRefClass( "EzAppFastp",
               contains = "EzApp",
               methods = list(
                 initialize = function(){
                   "Initializes the application using its specific defaults."
                   runMethod <<- ezMethodFastpTrim
                   name <<- "EzAppFastp"
                 }
               )
               
  )

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
##'   \item{minTrailingQuality}{ an integer specifying the minimal trailing quality to accept. Only used if > 0.}
##'   \item{minAvgQuality}{ an integer specifying the minimal average quality to accept. Only used if > 0.}
##'   \item{minReadLength}{ an integer specifying the minimal read length to accept.}
##'   \item{dataRoot}{ a character specifying the path of the data root to get the full column paths from.}
##' }
##' @template roxygen-template
##' @return Returns the output after trimming as an object of the class EzDataset.
ezMethodTrim = function(input=NA, output=NA, param=NA){
  
  if (any(grepl("bam$", input$getFullPaths("Read1")))){
    stop("can not process unmapped bam as input")
  }
  
  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")){
    output = input$copy()
    output$setColumn("Read1", paste0(getwd(), "/", input$getNames(), "-trimmed_R1.fastq"))
    if (param$paired){
      output$setColumn("Read2", paste0(getwd(), "/", input$getNames(), "-trimmed_R2.fastq"))
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
  if (!is.null(param$copyReadsLocally) && param$copyReadsLocally){
    input = copyReadsLocally(input, param)
  }
  
  param$trimSeedMismatches = 1
  param$trimPalindromClipThresh = 20
  param$trimSimpleClipThresh = 7
  param$trimMinAdaptLength = 5
  param$trimKeepBothReads = "true"
  param$trimQualWindowWidth = 4
  
  if (param$subsampleReads > 1 || param$nReads > 0){
    input = ezMethodSubsampleReads(input=input, param=param)
    inputRead1 <- input$getFullPaths("Read1")
    on.exit(file.remove(inputRead1), add=TRUE)
    if(param$paired){
      inputRead2 <- input$getFullPaths("Read2")
      on.exit(file.remove(inputRead2), add=TRUE)
    }
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
    adaptFile = "adapters.fa"
    adapters = c(adapter1, adapter2)
    if (!is.null(param$onlyAdapterFromDataset) && param$onlyAdapterFromDataset){
      writeXStringSet(adapters, adaptFile)
    } else {
      file.copy(from=TRIMMOMATIC_ADAPTERS,
                to=adaptFile)
      writeXStringSet(adapters, adaptFile, append=TRUE)
    }
    on.exit(file.remove(adaptFile), add=TRUE)
    
    trimAdaptOpt =  paste("ILLUMINACLIP", adaptFile, param$trimSeedMismatches, param$trimPalindromClipThresh,
                          param$trimSimpleClipThresh, param$trimMinAdaptLength, param$trimKeepBothReads, sep=":")
  } else {
    trimAdaptOpt = ""
  }
  
  if (param$minTailQuality > 0){
    tailQualOpt = paste("SLIDINGWINDOW", param$trimQualWindowWidth,
                        param$minTailQuality, sep=":")
  } else {
    tailQualOpt = ""
  }
  
  if(param$minLeadingQuality >0){
    leadingQualOpt <- paste("LEADING", param$minLeadingQuality, sep=":")
  }else{
    leadingQualOpt <- ""
  }

    if(param$minTrailingQuality >0){
    trailingQualOpt <- paste("TRAILING", param$minTrailingQuality, sep=":")
  }else{
    trailingQualOpt <- ""
  }
  
  if (param$minAvgQuality > 0){
    minAvgQualOpt = paste("AVGQUAL", param$minAvgQuality, sep=":")
  } else {
    minAvgQualOpt = ""
  }
  
  r1TmpFile = "trimmed_R1.fastq"
  r2TmpFile = "trimmed_R2.fastq"
  if (any(c(trimAdaptOpt, tailQualOpt, minAvgQualOpt) != "") || param$minReadLength > 0){
    if (param$paired){
      method = "PE"
      readOpts = paste(
        input$getFullPaths("Read1"),
        input$getFullPaths("Read2"),
        r1TmpFile,
        "unpaired_R1.fastq",
        r2TmpFile,
        "unpaired_R2.fastq")
        on.exit(file.remove(c("unpaired_R1.fastq", "unpaired_R2.fastq")),
                add=TRUE)
    } else {
      method = "SE"
      readOpts = paste(
        input$getFullPaths("Read1"), r1TmpFile)
    }
    cmd = paste(prepareTrimmomatic(), method,
                ## hardcode phred33 quality encoding
                "-threads", min(param$cores, 8), "-phred33",
                readOpts, trimAdaptOpt, tailQualOpt, leadingQualOpt, trailingQualOpt, minAvgQualOpt,
                paste("MINLEN", param$minReadLength, sep=":"),
                "> trimmomatic.out 2> trimmomatic.err")
    on.exit(file.remove(c("trimmomatic.out", "trimmomatic.err")), add=TRUE)
    
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
    cmd = paste("flexbar",
                "--threads", min(param$cores, 8),
                "-r", r1TmpFile,
                pairedOpt,
                "-u", 20, ##### max uncalled bases
                "--pre-trim-left", param$trimLeft,
                "--pre-trim-right", param$trimRight,
                minReadLengthOpt,
                "--target", "flexbar",
                "> flexbar.out 2> flexbar.err")
    ezSystem(cmd)
    on.exit(file.remove(c("flexbar.out", "flexbar.err", "flexbar.log")),
            add=TRUE)
    
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


## copies the read files locally and computes md5 sums
copyReadsLocally = function(input, param){
  if (param$paired){
    reads = c("Read1", "Read2")
  } else{
    reads = "Read1"
  }
  for (rds in reads){
    readFileIn = input$getFullPaths(rds)
    
    file.copy(from=readFileIn, to=".")
    
    input$setColumn(rds, basename(readFileIn))
    # if (Sys.info()["user"] == "trxcopy") { ## only run the check for the user trxcopy!!!
    #   md5Local = ezSystem(paste("md5sum", basename(readFileIn)), intern = TRUE)
    #   md5Local = sub(" .*", "", md5Local)
    #   md5File = file.path(dirname(readFileIn), "md5.txt")
    #   md5Remote = NA
    #   if (file.exists(md5File)){
    #     md5Set = ezRead.table(md5File)
    #     md5Remote = md5Set[basename(readFileIn), 1]
    #   }
    #   if (is.na(md5Remote)){
    #     md5Remote = ezSystem(paste("ssh fgcz-s-022 md5sum", readFileIn), intern = TRUE)
    #     md5Remote = sub(" .*", "", md5Remote)
    #   }
    #   stopifnot(md5Local == md5Remote)
    # }
  }
  input$dataRoot = NULL
  return(input)
}
### -----------------------------------------------------------------
### ezMethodSubsampleReads: subsample fastq
###
ezMethodSubsampleReads = function(input=NA, output=NA, param=NA){
  nYield=1e5
  require(ShortRead)
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
  totalReads = as.numeric(input$getColumn("Read Count"))
  if (param$nReads > 0){
    #subsampleFactor = shrinkToRange(totalReads / param$nReads, c(1, Inf))
    if (param$subsampleReads > 1){
      message("subsampleReads setting will be overwritten by nReads parameter")
    }
    nReads <- pmin(param$nReads, totalReads)
  } else {
    nReads <- as.integer(1/param$subsampleReads * totalReads)
  }
  for(i in 1:length(input$getFullPaths("Read1"))){
    idxMaster = sample.int(size=nReads[i], n=totalReads[i])
    fullFile = input$getFullPaths("Read1")[i]
    subFile = output$getColumn("Read1")[i]
    idx = idxMaster
    fqs <- FastqStreamer(fullFile, n=nYield)
    while(length(x <- yield(fqs))){
      use = idx <= length(x)
      writeFastq(x[idx[use]], file=sub(".gz$", "", subFile), mode="a", full=F, compress=F)
      idx = idx[!use] - length(x)
    }
    close(fqs)
    if (grepl(".gz$", subFile)){
      ezSystem(paste("pigz -p 2 ", sub(".gz$", "", subFile)))
    }
    if (param$paired){
      fullFile = input$getFullPaths("Read2")[i]
      subFile = output$getColumn("Read2")[i]
      idx = idxMaster
      fqs <- FastqStreamer(fullFile, n=nYield)
      while(length(x <- yield(fqs))){
        use = idx <= length(x)
        writeFastq(x[idx[use]], file=sub(".gz$", "", subFile), mode="a", full=F, compress=F)
        idx = idx[!use] - length(x)
      }
      close(fqs)
      if (grepl(".gz$", subFile)){
        ezSystem(paste("pigz -p 2 ", sub(".gz$", "", subFile)))
      }
    }
  }
  return(output)
}

##' @describeIn ezMethodTrim Performs the fastq for the subsamples using the package ShortRead.
##' @examples
##'  inputFile = system.file(package = "ezRun", "extdata/yeast_10k/wt_1_R1.fastq.gz")
##'  subsampledFile = "sub_R1.fastq"
##'  ezSubsampleFastq(inputFile, subsampledFile, subsampleFactor = 5)
# ezSubsampleFastq = function(full, sub, subsampleFactor=NA, nYield=1e5, overwrite=FALSE){
#   stopifnot(full != sub)
#   stopifnot(length(full) == length(sub))
#   if (any(file.exists(sub))){
#     filesToRemove = sub[file.exists(sub)]
#     if (!overwrite){
#       stop("files do exist: ", filesToRemove)
#     }
#     warning("removing first: ", filesToRemove)
#     file.remove(filesToRemove)
#   }
#   require("ShortRead")
#   nms = names(full)
#   if (is.null(nms)){
#     nms = full
#   }
#   nReadsVector = integer()
#   for (i in seq_along(full)){
#     nReads = 0
#     fqs = FastqStreamer(full[i], n = nYield)
#     idx = seq(from=1, to=nYield, by=subsampleFactor)
#     tmpFile = sub("\\.fastq.*", "_temp.fastq", sub[i])
#     while(length(x <- yield(fqs))){
#       if (length(x) >= nYield){
#         writeFastq(x[idx], file=tmpFile, mode="a", full=F, compress=F)
#         nReads = nReads + length(idx)
#       } else {
#         writeFastq(x[idx[idx<length(x)]], file=tmpFile, mode="a", full=F, compress=F)
#         nReads = nReads + sum(idx<length(x))
#       }
#     }
#     close(fqs)
#     nReadsVector[nms[i]] = nReads
#     if (grepl(".gz$", sub[i])){
#       ezSystem(paste("pigz -p 2 --best -c", tmpFile, ">", sub[i]))
#       file.remove(tmpFile)
#     } else {
#       file.rename(tmpFile, sub[i])
#     }
#   }
#   return(nReadsVector)
# }
