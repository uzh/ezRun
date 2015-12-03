context("Runs all function examples and cleans them up afterwards.")

test_that("Function examples", {
  cwd = getwd()
  setwdNew("run_examples")
  devtools::run_examples()
  setwd(cwd)
  ezSystem(paste("rm -fr", system.file("./run_examples", package="ezRun", mustWork=TRUE)))
  remove1 = system.file("DESCRIPTION_head", package="ezRun", mustWork=TRUE)
  remove2 = system.file("extdata/genes.bed", package="ezRun", mustWork=TRUE)
  remove3 = system.file("extdata/genesWithPrespliced.gtf", package="ezRun", mustWork=TRUE)
  cmd = paste("rm -fr", remove1, remove2, remove3)
  ezSystem(cmd)
})
