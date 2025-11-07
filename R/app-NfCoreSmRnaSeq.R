###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNfCoreSmRnaSeq <- function(input = NA, output = NA, param = NA) {
  refbuild = param$refBuild
  outFolder = output$getColumn("Result") |> basename()
  
  nfSampleFile <- file.path('datasetSmRnaSeq.csv')
  nfSampleInfo = getSmRnaSeqSampleSheet(input)
  write_csv(nfSampleInfo, nfSampleFile)
  
  
  setNFCacheDir()
  cmd = paste(
    "nextflow run nf-core/smrnaseq",
    ## i/o
    "--input", nfSampleFile,
    "--outdir", outFolder,
    ## genome 
    "--genome", param[['referenceGenome']],
    ## parameters
    "--mirtrace_species", param[['mirtraceSpecies']],
    ## configuration
    "-work-dir work",
    "-profile apptainer",
    "-r 2.4.0"
  )
  ezSystem(cmd)
  
  return("Success")
}


EzAppNfCoreSmRnaSeq <- setRefClass(
  "EzAppNfCoreSmnRnaSeq",
  contains = "EzApp",
  methods = list(
    initialize = function()
    {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodNfCoreSmRnaSeq
      name <<- "EzAppNfCoreSmRnaSeq"
    }
  )
)

##' @description get an nf-core/smrnaseq-formatted csv file
getSmRnaSeqSampleSheet <- function(input){
  nfSampleInfo <- ezFrame(
    sample = names(input$getFullPaths("Read1")),
    fastq_1 = input$getFullPaths("Read1")
  )
  
  return(nfSampleInfo)
}

