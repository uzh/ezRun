###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodGatkRna = function(input=NA, output=NA, param=NA){
  meta = input$meta
  pathSampleInfo = "00_samples.txt"
  pathRef = param$ezRef["refFastaFile"]
  vcfName = basename(output$getColumn("VCF"))
  sampleNames = rownames(meta)
  sampleInfo = data.frame(
    Sample = sampleNames,
    Library = sampleNames,
    Platform = rep(param$rgpl, length(sampleNames)),
    RGID = sampleNames,
    FileBAM = input$getFullPaths("BAM")
  )
  write.table(sampleInfo, pathSampleInfo, quote = F, sep = "\t", row.names = F)
  
  pathFilter = system.file("extdata/GATK/RnaFilters.txt", package = "ezRun", mustWork = T)
  pathGtFilter = system.file("extdata/GATK/RnaGtFilters.txt", package = "ezRun", mustWork = T)
  pathScript = system.file("extdata/GATK/VariantCallingRna.scala", package = "ezRun", mustWork = T)
  noSoftClip = ! is.null(param$dontUseSoftClippedBases) && as.logical(param$dontUseSoftClippedBases)
  # All jobs should be submitted to the same node, so that they could access the data
  cmd = paste("java -jar $Queue_jar",
              "-R", pathRef,
              "-libs", pathSampleInfo,
              "-t", param$threads,
              "-ct", param$cthreads,
              "-sj", param$sj,
              "-maxMem", param$ram,
              "-maxConcurrentRun", param$maxConcurrentRun,
              ifelse(noSoftClip, "--dontUseSoftClippedBases", ""),
              "-filters", pathFilter,
              "-gtfilters", pathGtFilter,
              "-pathOut", vcfName,
              "-qsub",
              "-jobQueue", "GT",
              "-jobNative", '"', "-l", paste0("h=", Sys.info()["nodename"]), '"',
              "-jobParaEnv", "smp",
              "-runDir", ".",
              "--script", pathScript,
              "-run")
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodGatkRna(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppGatkRna <-
  setRefClass("EzAppGatkRna",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGatkRna
                  name <<- "EzAppGatkRna"
                }
              )
  )
