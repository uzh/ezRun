#options(error=recover)
require(ezRun)

ds = EzDataset$new(file=system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE))

ds2 = ds$copy()

ds$getNames()
ds$file
ds$getColumn("Read1")
ds$meta
ds$columnHasTag("File")
#ds[,!ds$columnHasTag("File")]
ds$meta[,!ds$columnHasTag("File")]
EzDataset$methods()

demoApp = EzApp$new()
demoApp$runMethod

demoMethod = function(input, output, param){
  print(input)
}


DemoApp <-
  setRefClass("demoApp",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  runMethod <<- demoMethod
                  name <<- "demoApp"
                }
              ))

demoApp = DemoApp$new()
demoApp$run(input=ds, output=ds, param=list(process_mode="DATASET"))



