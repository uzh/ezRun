###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodSalmon = function(input=NA, output=NA, param=NA){
  sampleName = input$getNames()
  ref = getSalmonReference(param)
  refIdx = paste0(ref, ".idx")

  iftrue <- function(p, yes, no = "") {
    if (!is.null(p) && p) { yes } else { no }
  }
  condCharAdd <- function(name, val) {
    if (ezIsSpecified(val)) { paste(name, val) } else { "" }
  }
  condNumAdd <- function(name, val) {
    if (!is.null(val) && val != 0) { paste(name, val) } else { "" }
  }

  libType = paste0(
    iftrue(param$paired, "I", ""),
    switch(param$strandMode, both = "U", "S"),
    switch(param$strandMode, sense = "F", antisense = "R", "")
  )
  if(param$strandMode == 'NA'){
    libType = 'A'
  }
  trimmedInput = ezMethodTrim(input = input, param = param)

  opt = paste(
      "-i", refIdx,
      "-o", param$outputDir,
      "-p", ezThreads(),
      "-l", libType,
      iftrue(param$paired,
        paste("-1", trimmedInput$getColumn("Read1"), "-2", trimmedInput$getColumn("Read2")),
        paste("-r", trimmedInput$getColumn("Read1"))
      ),
      condNumAdd("--fldMean", param$fldMean),
      condNumAdd("--fldSD", param$fldSD),
      iftrue(param$useVBOpt, "--useVBOpt", ""),
      condNumAdd("--numBootstraps", param$numBootstraps),
      condNumAdd("--numGibbsSamples", param$numGibbsSamples),
      iftrue(param$seqBias, "--seqBias"),
      iftrue(param$gcBias, "--gcBias"),
      iftrue(param$posBias, "--posBias"),
      if (ezIsSpecified(param$specialParams)) { param$specialParams } else { "" }
  )

  pathQuants = file.path(param$outputDir, "quant.sf")

  cmd = paste("salmon", "quant", opt, "2> salmon.stderr > salmon.stdout")
  ezSystem(cmd)

  ezSystem(paste("mv", pathQuants, basename(output$getColumn("Count"))))

  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodSalmon(input=NA, output=NA, param=NA)
##' @description Use this reference class to run Salmon
##' @seealso \code{\link{getSalmonReference}}
##' @seealso \code{\link{ezMethodTrim}}
##' @author Roman Briskine
EzAppSalmon <-
  setRefClass("EzAppSalmon",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSalmon
                  name <<- "EzAppSalmon"
                  appDefaults <<- rbind(
                    fldMean = ezFrame(
                      Type = "integer",
                      DefaultValue = 0,
                      Description = 'expected mean fragment length for single-end reads'
                    ),
                    fldSD = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0,
                      Description = 'expected standard deviation of the fragment length distribution for single-end reads'
                    ),
                    useVBOpt = ezFrame(
                      Type = "logical",
                      DefaultValue = F,
                      Description = 'use the variational Bayesian EM algorithm rather than the "standard" EM algorithm to optimize abundance estimates'
                    ),
                    numBootstraps = ezFrame(
                      Type = "integer",
                      DefaultValue = 0,
                      Description = 'number of bootstrapped samples; by default, Salmon does not use bootstrapping'
                    ),
                    numGibbsSamples = ezFrame(
                      Type = "integer",
                      DefaultValue = 0,
                      Description = 'an alternative to bootstrapping for estimating the variance in abundance estimates'
                    ),
                    seqBias = ezFrame(
                      Type = "logical",
                      DefaultValue = F,
                      Description = 'enable sequence bias correction'
                    ),
                    gcBias = ezFrame(
                      Type = "logical",
                      DefaultValue = F,
                      Description = 'enable GC bias correction'
                    ),
                    posBias = ezFrame(
                      Type = "logical",
                      DefaultValue = F,
                      Description = 'enable position bias correction'
                    ),
                    specialParams = ezFrame(
                      Type = "character",
                      DefaultValue = '',
                      Description = 'additional command line parameters to pass to Salmon'
                    ),
                    outputDir = ezFrame(
                      Type = "character",
                      DefaultValue = ".",
                      Description = "Output directory"
                    )
                  )
                }
              )
  )

getSalmonReference = function(param){

  if (ezIsSpecified(param$transcriptFasta)){
    refBase = file.path(getwd(), "salmonIndex/transcripts") 
    #paste0(file_path_sans_ext(param$trinityFasta), "_salmonIndex/transcripts")
  } else {
    if(ezIsSpecified(param$transcriptTypes)){
      salmonBase <- paste(sort(param$transcriptTypes), collapse="-")
      ## This is a combination of transcript types to use.
    }else{
      salmonBase <- ""
    }
    refBase = ifelse(param$ezRef["refIndex"] == "",
                     sub(".gtf$",
                         paste0("_", salmonBase, "_salmonIndex/transcripts"),
                         #"_salmonIndex/transcripts", 
                         param$ezRef["refFeatureFile"]),
                     param$ezRef["refIndex"])
  }
  lockFile = file.path(dirname(refBase), "lock")
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait
    Sys.sleep( 60)
    i = i + 1
  }
  if (file.exists(lockFile)){
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles = list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) > 0 ){
    ## we assume the index is built and complete
    return(refBase)
  }
  ## we have to build the reference
  wd = getwd()
  dir.create(dirname(refBase))
  setwd(dirname(refBase))
  ezWrite(Sys.info(), con=lockFile)
  on.exit(file.remove(lockFile))
  on.exit(setwd(wd), add=TRUE)

  job = ezJobStart("salmon index")
  if (ezIsSpecified(param$transcriptFasta)){
    pathTranscripts = param$transcriptFasta
  } else{
    pathTranscripts <- paste0(refBase, ".fa")
    transcripts = getTranscriptSequences(param)
    if(ezIsSpecified(param$transcriptTypes)){
      ## Prepare the subset of transcipts sequence
      seqAnno = ezFeatureAnnotation(param$ezRef@refAnnotationFile,
                                    dataFeatureType="transcript")
      transcriptsUse = rownames(seqAnno)[seqAnno$type %in% param$transcriptTypes]
      transcripts <- transcripts[unique(transcriptsUse)]
    }
    writeXStringSet(transcripts, pathTranscripts)
  }
  cmdTemplate = "salmon index -i %s.idx -t %s"
  cmd = sprintf(cmdTemplate, refBase, pathTranscripts)
  ezSystem(cmd)
  ezWriteElapsed(job, "done")
  return(refBase)
}
