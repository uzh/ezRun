# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodPbmm2 <- function(input = NA, output = NA, param = NA) {
  ref <- getPbmm2Reference(param)
  bamFile <- output$getColumn("BAM")
  sampleName <- sub(".bam", "", basename(bamFile))
  Input <- input$getFullPaths("Read1")
  defOpt <- paste("--sort -m 2000M -J", param$cores)
  inputFileType = param$ReadOpt
  readGroupOpt <- paste0(
    "--rg-id ", sampleName, " --rg SM:", sampleName,
    " --rg LB:RGLB_", sampleName,
    " --rg PL:pacbio", " --rg PU:RGPU_", sampleName
  )
  cmd <- paste(
    "/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbmm2 align", param$cmdOptions, defOpt, "--preset", inputFileType, "--rg", readGroupOpt, paste0(ref, ".", inputFileType, ".mmi"), Input, bamFile, "2>", paste0(sampleName, "_pbmm2.log"))
  ezSystem(cmd)


  cmd <- paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbindex", bamFile)
  ezSystem(cmd)
  cmd <- paste("samtools index", bamFile)
  ezSystem(cmd)

  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Pbmm2
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refBuildDir}{ a character specifying the directory of the reference build.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getPbmm2Reference <- function(param) {
  refBase <- ifelse(param$ezRef["refIndex"] == "",
    file.path(param$ezRef["refBuildDir"], "Sequence/Pbmm2Index/genome"),
    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    setwd(dirname(refBase))

    fastaFile <- param$ezRef["refFastaFile"]
    ezSystem(paste("ln -s", fastaFile, "."))
    cmd <- paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbmm2 index --preset HIFI", basename(fastaFile), paste0(basename(refBase), "HIFI.mmi"))
    ezSystem(cmd)
    cmd <- paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbmm2 index --preset SUBREAD", basename(fastaFile), paste0(basename(refBase),"SUBREAD.mmi"))
    ezSystem(cmd)
    cmd <- paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbmm2 index --preset ISOSEQ", basename(fastaFile), paste0(basename(refBase),"ISOSEQ.mmi"))
    ezSystem(cmd)
    cmd <- paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbmm2 index --preset UNROLLED", basename(fastaFile), paste0(basename(refBase),"UNROLLED.mmi"))
    ezSystem(cmd)
    # ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles <- list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) < 3) {
    ## we assume the index is built and complete
    stop(paste("index not available: ", refBase))
  }
  return(refBase)
}

##' @template app-template
##' @templateVar method ezMethodPbmm2(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getPbmm2Reference}}
EzAppPbmm2 <-
  setRefClass("EzAppPbmm2",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodPbmm2
        name <<- "EzAppPbmm2"
        appDefaults <<- rbind(
        ReadOpt = ezFrame(Type="character",  DefaultValue="HIFI",  Description="input read types: SUBREAD, CCS, HIFI, ISOSEQ, UNROLLED. Default is HIFI")
	)
      }
    )
  )

