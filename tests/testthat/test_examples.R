context("Runs all function examples and cleans them up afterwards.")

cwd = getwd()
rm0 = system.file("tests/testthat/run_examples", package="ezRun")
rm1 = system.file("DESCRIPTION_head", package="ezRun")
rm2 = system.file("extdata/genesWithPrespliced.gtf", package="ezRun")
cmd = paste("rm -fr", rm0, rm1, rm2)

test_that("Function examples", {
  setwdNew("./run_examples/")
  devtools::run_examples()
  setwd(cwd)
  ezSystem(cmd)
})
