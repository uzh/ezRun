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

ezMethodFastpTrim <- function(input = NA, output = NA, param = NA) {
  require(withr)
  if (any(str_detect(input$getFullPaths("Read1"), "bam$"))) {
    stop("cannot process unmapped bam as input")
  }

  if (!ezIsSpecified(param$gzipTrimmed)) {
    param$gzipTrimmed <- TRUE
  }

  ## if output is not an EzDataset, set it! (when ezMethodFastpTrim is used inside another app)
  if (!is(output, "EzDataset")) {
    output <- input$copy()
    output$setColumn("Read1", file.path(getwd(), str_c(input$getNames(), "-trimmed_R1.fastq",
                                                       if_else(param$gzipTrimmed, ".gz", ""))))
    if (param$paired) {
      output$setColumn("Read2", file.path(getwd(), str_c(input$getNames(), "-trimmed_R2.fastq",
                                                         if_else(param$gzipTrimmed, ".gz", ""))))
    } else {
      if ("Read2" %in% input$colNames) {
        output$setColumn("Read2", NULL)
      }
    }
    output$dataRoot <- NULL
  }
  ## if there are multiple samples loop through them
  if (input$getLength() > 1) {
    for (nm in input$getNames()) {
      ezSystem(str_c("touch ", nm, "_preprocessing.log"))
      ezMethodFastpTrim(input$subset(nm), output$subset(nm), param)
      ## NOTE: potential risk, temp files might be overwritten
    }
    return(output)
  }
  ## now we deal only with one sample!

  ## make a local copy of the dataset and check the md5sum
  if (!is.null(param$copyReadsLocally) && param$copyReadsLocally) {
    input <- copyReadsLocally(input, param)
  }

  if (param$subsampleReads > 1 || param$nReads > 0) {
    if (param$subsampleReads > 1 && param$nReads > 0){
      stop("can not have subsampleReads and nReads")
    }
    totalReads <- as.numeric(input$getColumn("Read Count"))
    if (param$nReads > 0) {
      nReads <- min(param$nReads, totalReads)
    }
    if (param$subsampleReads > 1){
      nReads <- as.integer(1 / param$subsampleReads * totalReads)
    }
    input <- ezMethodSubsampleFastq(input, param=param, n=nReads) 
  }

  ## fastp
  r1TmpFile <- "trimmed_R1.fastq.gz"
  if (param$paired) {
    r2TmpFile <- "trimmed_R2.fastq.gz"
    readsInOut <- str_c(
      "--in1", input$getFullPaths("Read1"),
      "--in2", input$getFullPaths("Read2"),
      "--out1", r1TmpFile,
      "--out2", r2TmpFile,
      sep = " "
    )
  } else {
    readsInOut <- str_c(
      "--in1", input$getFullPaths("Read1"),
      "--out1", r1TmpFile,
      sep = " "
    )
  }
  ## adapter trimming options
  if (param[["trimAdapter"]]) {
    if (!is.null(input$meta$Adapter1) && !is.na(input$meta$Adapter1) &&
        input$meta$Adapter1 != "") {
      adapter1 <- DNAStringSet(input$meta$Adapter1)
      names(adapter1) <- "GivenAdapter1"
    } else {
      adapter1 <- DNAStringSet()
    }
    # read2 (if paired)
    if (param$paired && !is.null(input$meta$Adapter2) && !is.na(input$meta$Adapter2) &&
        input$meta$Adapter2 != "") {
      adapter2 <- DNAStringSet(input$meta$Adapter2)
      names(adapter2) <- "GivenAdapter2"
    } else {
      adapter2 <- DNAStringSet()
    }
    adaptFile <- "adapters.fa"
    adapters <- c(adapter1, adapter2)
    if (!is.null(param$onlyAdapterFromDataset) && param$onlyAdapterFromDataset) {
      # take only adapter from dataset and ignore the ones from TRIMMOMATIC_ADAPTERS
      writeXStringSet(adapters, adaptFile)
    } else {
      file.copy(
        from = TRIMMOMATIC_ADAPTERS,
        to = adaptFile
      )
      writeXStringSet(adapters, adaptFile, append = TRUE)
    }

    trimAdapt <- str_c("--adapter_fasta", adaptFile, sep = " ")
  } else {
    ezSystem("touch adapters.fa")
    adaptFile <- "adapters.fa"
    trimAdapt <- "--disable_adapter_trimming"
  }

  ## paste command
  cmd <- str_c(
    "fastp",
    readsInOut,
    # general options
    "--thread", param$cores,
    # global trimming
    "--trim_front1", param$trim_front1,
    "--trim_tail1", param$trim_tail1,
    # quality-based trimming per read
    if_else(param$cut_front,
            str_c("--cut_front", "--cut_front_window_size", param$cut_front_window_size,
                  "--cut_front_mean_quality", param$cut_front_mean_quality, sep = " "),
            ""), # like Trimmomatic's LEADING
    if_else(param$cut_right,
            str_c("--cut_right", "--cut_right_window_size", param$cut_right_window_size,
                  "--cut_right_mean_quality", param$cut_right_mean_quality, sep = " "),
            ""), # like Trimmomatic's SLIDINGWINDOW
    if_else(param$cut_tail,
            str_c("--cut_tail", "--cut_tail_window_size", param$cut_tail_window_size,
                  "--cut_tail_mean_quality", param$cut_tail_mean_quality, sep = " "),
            ""), # like Trimmomatic's TRAILING
    "--average_qual", param$average_qual,
    # adapter trimming
    trimAdapt,
    # read length trimming
    "--max_len1", param$max_len1,
    "--max_len2", param$max_len2,
    # polyX
    "--trim_poly_x", "--poly_x_min_len", param$poly_x_min_len,
    # read length filtering
    "--length_required", param$length_required,
    # compression output
    "--compression 4",
    sep = " "
  )

  ## run
  if (ezIsSpecified(param$cmdOptionsFastp)) {
    cmd <- str_c(cmd, param$cmdOptionsFastp, sep = " ")
  }
  if ("PreprocessingLog" %in% output$colNames) {
    logFile <- basename(output$getColumn("PreprocessingLog"))
  } else {
    logFile <- str_c(output$getNames(), "_preprocessing.log")
  }
  ezSystem(str_c(cmd, "2>", logFile, sep = " "))
  ezSystem(str_c("cat", "fastp.json", ">>", logFile, sep = " "))

  ## rename trimmed output
  if (param$gzipTrimmed) {
    file.rename(from = r1TmpFile, to = basename(output$getColumn("Read1")))
    if (param$paired) {
      file.rename(from = r2TmpFile, to = basename(output$getColumn("Read2")))
    }
  } else {
    ezSystem(str_c("zcat", r1TmpFile, ">", basename(output$getColumn("Read1")), sep = " "))
    file.remove(r1TmpFile)
    if (param$paired) {
      ezSystem(str_c("zcat", r2TmpFile, ">", basename(output$getColumn("Read2")), sep = " "))
      file.remove(r2TmpFile)
    }
  }

  return(output)
}

##' @title EzAppFastp app
##' @description fast read pre-processing.
##' @author Miquel Anglada Girotto
EzAppFastp <-
  setRefClass("EzAppFastp",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodFastpTrim
        name <<- "EzAppFastp"
        appDefaults <<- rbind(gzipTrimmed = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "whether to return gzipped reads"))
      }
    )
  )

