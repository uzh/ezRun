context("small RNA apps with example data")

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

yeastCommonSmallRnaParam = function() {
  param = list()
  param[['cores']] = '8'
  param[['ram']] = '16'
  param[['scratch']] = '100'
  param[['node']] = ''
  param[['process_mode']] = 'DATASET'
  param[[
    'refBuild'
  ]] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  param[['mail']] = ''
  param[['dataRoot']] = '/srv/gstore/projects'
  param[['resultDir']] = 'p1001/smRNA_Result'
  return(param)
}

test_that("smRNA Ncpro", {
  skipLong()
  ezSystem("rm -fr /scratch/test_ncpro/*")
  setwdNew("/scratch/test_ncpro")
  param = yeastCommonSmallRnaParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/smRNA_250k/datasetWithAdapter.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  output = list()
  output[['Name']] = 'ncPRO_Result'
  output[['Species']] = ''
  output[[
    'refBuild'
  ]] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  output[['refBuild']] = 'Mus_musculus/UCSC/mm10/Annotation/Version-2012-05-23'
  output[[
    'Report [File]'
  ]] = 'p1001/Count_ncPRO_Report_5750_2015-12-18--12-33-40/ncPRO_Result'
  output[[
    'Html [Link]'
  ]] = 'p1001/Count_ncPRO_Report_5750_2015-12-18--12-33-40/ncPRO_Result/ncpro/report.html'
  output[[
    'TrimCounts [Link]'
  ]] = 'p1001/Count_ncPRO_Report_5750_2015-12-18--12-33-40/ncPRO_Result/trimCounts-barplot.png'
  param[['name']] = 'ncPRO_Result'
  param[['mail']] = ''
  param[['refBuild']] = 'Mus_musculus/UCSC/mm10/Annotation/Version-2012-05-23'
  myApp = EzAppNcpro$new()
  myApp$run(input = input, output = output, param = param)
  setwd(cwd)
})
