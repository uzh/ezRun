###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodKallisto = function(input=NA, output=NA, param=NA){
  sampleName = input$getNames()
  ref = getKallistoReference(param)
  refIdx = paste0(ref, ".idx")
  if(!param$paired){
    if(param$"fragment-length" == 0){
      param$"fragment-length" = 180
    }
    if(param$sd == 0){
      param$sd = 50
    }
  }
  
  iftrue <- function(p, yes, no = "") {
    if (!is.null(p) && p) { yes } else { no }
  }
  condCharAdd <- function(name, val) {
    if (ezIsSpecified(val)) { paste(name, val) } else { "" }
  }
  condNumAdd <- function(name, val) {
    if (!is.null(val) && val != 0) { paste(name, val) } else { "" }
  }

  opt = paste(
      "-i", refIdx,
      "-o", param$outputDir,
      "-t", ezThreads(),
      iftrue(param$bias, "--bias"),
      condCharAdd("--bootstrap-samples", param$"bootstrap-samples"),
      condCharAdd("--seed", param$seed),
      iftrue(param$paired, "", "--single"),
      strandOpt = switch(param$strandMode,
                         "sense"="--fr-stranded",
                         "antisense"="--rf-stranded",
                         "both"=""),
      condNumAdd("--fragment-length", param$"fragment-length"),
      condNumAdd("--sd", param$sd),
      iftrue(param$pseudobam, "--pseudobam")
  )

  trimmedInput = ezMethodTrim(input = input, param = param)

  pathFastqFiles = paste(
      trimmedInput$getColumn("Read1"),
      iftrue(param$paired, trimmedInput$getColumn("Read2"))
  )

  pathAbundance.h5 = file.path(param$outputDir, "abundance.h5")
  pathAbundance.tsv = file.path(param$outputDir, "abundance.tsv")
  pathRunInfo = file.path(param$outputDir, "run_info.json")
  pathPseudobam = file.path(param$outputDir, "pseudoalignments.bam")

  cmd = paste(
      "kallisto",
      "quant",
      opt,
      pathFastqFiles,
      "2> kallisto.stderr",
      " > kallisto.stdout"
  )
  ezSystem(cmd)

  ezSystem(paste("mv", pathAbundance.tsv, basename(output$getColumn("Count"))))
  if (!is.null(param$"bootstrap-samples") && param$"bootstrap-samples" > 0) {
    ezSystem(paste("mv", pathAbundance.h5, basename(output$getColumn("bootstrappedCount"))))
  }
  ezSystem(paste("mv", pathRunInfo, basename(output$getColumn("runInfo"))))
  if (!is.null(param$pseudobam) && param$pseudobam){
    fnBam = basename(output$getColumn("BAM"))
    ezSortIndexBam(pathPseudobam, fnBam, ram=param$ram, cores=min(param$cores, 8))
  }

  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodKallisto(input=NA, output=NA, param=NA)
##' @description Use this reference class to run kallisto
##' @seealso \code{\link{getKallistoReference}}
##' @seealso \code{\link{ezMethodTrim}}
##' @author Roman Briskine
EzAppKallisto <-
  setRefClass("EzAppKallisto",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodKallisto
                  name <<- "EzAppKallisto"
                  appDefaults <<- rbind(
                    "bootstrap-samples" = ezFrame(
                      Type = "integer",
                      DefaultValue = 10,
                      Description = "Number of bootstrap samples"
                    ),
                    seed = ezFrame(
                      Type = "numeric",
                      DefaultValue = 42,
                      Description = "seed for the bootstrap sampling"
                    ),
                    "fragment-length" = ezFrame(
                      Type = "integer",
                      DefaultValue = 0,
                      Description = "estimated average fragment length (required for single-end reads but should be set to 0 for paired-end reads)"
                    ),
                    sd = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0,
                      Description = 'estimated fragment length standard deviation (required for single-end reads but should be set to 0 for paired-end reads)'
                    ),
                    bias = ezFrame(
                      Type = "logical",
                      DefaultValue = T,
                      Description = 'perform sequence based bias correction'
                    ),
                    pseudobam = ezFrame(
                      Type = "logical",
                      DefaultValue = F,
                      Description = 'generate a bam file with pseudoalignments'
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


getKallistoReference = function(param){

  if (ezIsSpecified(param$transcriptFasta)){
    refBase = file.path(getwd(), "kallistoIndex/transcripts") 
    #paste0(file_path_sans_ext(param$trinityFasta), "_kallistoIndex/transcripts")
  } else {
    if(ezIsSpecified(param$transcriptTypes)){
      kallistoBase <- paste(sort(param$transcriptTypes), collapse="-")
      ## This is a combination of transcript types to use.
    }else{
      kallistoBase <- ""
    }
    refBase = ifelse(param$ezRef["refIndex"] == "",
                     sub(".gtf$", 
                         paste0("_", kallistoBase, "_kallistoIndex/transcripts"),
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

  job = ezJobStart("kallisto index")
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
      transcriptsUse <- intersect(transcriptsUse, names(transcripts))
      transcripts <- transcripts[unique(transcriptsUse)]
    }
    writeXStringSet(transcripts, pathTranscripts)
  }
  cmdTemplate = "kallisto index -i %s.idx %s"
  cmd = sprintf(cmdTemplate, refBase, pathTranscripts)
  ezSystem(cmd)
  ezWriteElapsed(job, "done")
  return(refBase)
}
