context("ChIP apps with example data")

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

mouseCommonChipParam = function() {
  param = list()
  param[['cores']] = '1'
  param[['ram']] = '16'
  param[['scratch']] = '100'
  param[['node']] = ''
  param[['process_mode']] = 'SAMPLE'
  param[['refBuild']] = 'Mus_musculus/UCSC/mm10/Annotation/Version-2012-05-23'
  param[['paired']] = 'false'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['mail']] = ''
  param[['dataRoot']] = "/srv/gstore/projects/" ##system.file(package="ezRun", mustWork = TRUE)
  param[['specialOptions']] = ""
  return(param)
}

test_that("ChIP_Macs2", {
  skipLong()
  ezSystem("rm -fr /scratch/test_macs2/*")
  setwdNew("/scratch/test_macs2")
  param = mouseCommonChipParam()
  output = EzDataset$new(
    file = "/scratch/ExampleData_Macs2/PeakCalling_MACS2_7312_TEST_2015-09-01--17-43-41/dataset.tsv",
    dataRoot = param$dataRoot
  )
  input = EzDataset$new(
    file = "/scratch/ExampleData_Macs2/PeakCalling_MACS2_7312_TEST_2015-09-01--17-43-41/input_dataset.tsv",
    dataRoot = param$dataRoot
  )
  param[['useControl']] = "false"
  param[['cmdOptions']] = "--nomodel --extsize 147 -g hs --bw 200"
  param[['dataRoot']] = "/srv/gstore/projects"
  myApp = EzAppMacs2$new()
  myApp$run(input = input$subset(1), output = output$subset(1), param = param)
  setwd(cwd)
})
