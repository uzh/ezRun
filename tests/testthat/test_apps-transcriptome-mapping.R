context("Mapping apps with example data")

## this tests do take long therefore we only run them if the environment variable RUN_LONG_TEST is set to TRUE
# Sys.setenv(RUN_LONG_TEST=TRUE)

cwd = getwd()

skipLong = function() {
  if (Sys.getenv("RUN_LONG_TEST") == "TRUE") {
    return()
  } else {
    skip("not running lengthy tests")
  }
}

yeastCommonMapParam = function() {
  param = list()
  param[['cores']] = '8'
  param[['ram']] = '20'
  param[['scratch']] = '100'
  param[['node']] = ''
  param[['process_mode']] = 'SAMPLE'
  param[[
    'refBuild'
  ]] = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07'
  param[['paired']] = 'true'
  param[['cmdOptions']] = ''
  param[['strandMode']] = 'sense'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['trimAdapter']] = 'false'
  param[['trimLeft']] = '0'
  param[['trimRight']] = '0'
  param[['minTailQuality']] = '0'
  param[['specialOptions']] = ''
  param[['mail']] = ''
  param[['dataRoot']] = system.file(package = "ezRun", mustWork = TRUE)
  param[['resultDir']] = 'p1001/Map_Result'
  return(param)
}


test_that("Map_Bowtie2Transcriptome", {
  skipLong()
  ezSystem("rm -fr /scratch/test_bowtie2transcriptomeMapping/*")
  setwdNew("/scratch/test_bowtie2transcriptomeMapping")
  param = yeastCommonMapParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  output = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_Bowtie2Transcriptome/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  param[['cmdOptions']] = '--no-unal'
  myApp = EzAppBowtie2Transcriptome$new()
  myApp$run(
    input = input$copy()$subset(1),
    output = output$copy()$subset(1),
    param = param
  )
  setwd(cwd)
})


test_that("Profiles_TranscriptCoverage", {
  skipLong()
  ezSystem("rm -fr /scratch/test_TranscriptCoverage/*")
  setwdNew("/scratch/test_TranscriptCoverage")
  param = yeastCommonMapParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_Bowtie2Transcriptome/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  meta = input$meta
  meta$"Count [File]" = paste0(input$getNames(), ".txt")
  meta$"Profiles [File]" = paste0(input$getNames(), "-profiles.RData")
  output = EzDataset$new(meta = meta, dataRoot = param$dataRoot)
  param[['minReadLength']] = 35
  param[['maxReadLength']] = 37
  param[["paired"]] = FALSE
  param[["getCoverageByReadLength"]] = TRUE

  myApp = EzAppTranscriptCoverage$new()
  myApp$run(
    input = input$copy()$subset(1),
    output = output$copy()$subset(1),
    param = param
  )
  setwd(cwd)
})
