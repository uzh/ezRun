
.libPaths("/srv/GT/analysis/course_sushi/lib")

library(ezRun)


######### some dummy tests
param = ezParam(userParam = list("process_mode"="DATASET"))

tagList = ezTagListFromNames(c("foo", "bar", "tagged [thetag]", "double tagged [tag1,tag2]"))
grepl("tag2", tagList)

### test the dataset class
ds = EzDataset(file="~/R/ezRun/inst/extdata/ventricles_10k/dataset.tsv")
ds2 = EzDataset(meta=ds$meta)
EzDataset()

options(error=recover)
ds2[, 3]


demoAppMethod = function(input=NULL, output=NULL, param=NULL){
  print(input)
  return("Success")
}



ezApp = EzApp(runMethod=demoAppMethod, name="demoAppName")
ezApp$runMethod(5)

input = list(a=5, b=10)
output = list("Result [File]"="foo.txt")

## directly calling the method works
ezApp$runMethod(input, output, param)

## with the generic we do add some pre and post-processing
#options(error=recover)
ezRunApp(ezApp, input, output, param)

## now an app that fails
failAppMethod = function(input=NULL, output=NULL, param=NULL){
  log("a")
  return("Success")
}
failApp = EzApp(runMethod=failAppMethod, name="failApp")
#options(error=NULL)
ezRunApp(failApp, input, output, param)


## method needs to be defined first
# fooApp = EzApp(runMethod=fooAppMethod, name="demoApp")
# getParam()



## test the real ngs count qc app
# input = list(a=5, b=10)
# output = list("Result [File]"="foo.txt")
# param = ezParam()
# ezRunApp(ezAppCountQC, input, output, param)
# 
# 
# tryCatch(1, error=function(e){traceback(); deparse(sys.calls()[[sys.nframe()-1]])})
# tryCatch(log("a"), error=function(e){traceback(); deparse(sys.calls()[[sys.nframe()-1]])})
# tryCatch(log("a"), error=function(e){print(e);traceback(); sys.calls()[[sys.nframe()-1]]})
# dump.frames




####### test the ref class ----------



EzRef(param=list(refBuild="foo"))


