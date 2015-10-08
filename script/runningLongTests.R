
## settings for running the tests as Hubert
.libPaths("/srv/GT/analysis/course_sushi/lib")
require(ezRun)
EZ_PARAM_DEFAULTS["adminMail", "DefaultValue"] = "Hubert.Rehrauer@fgcz.ethz.ch"
Sys.setenv(RUN_LONG_TEST=TRUE)
devtools::test()


## settings for running the tests as Peter
.libPaths("/home/pdschmid/R-libs")
setwd("/home/pdschmid/ezRun")
require(ezRun)
EZ_PARAM_DEFAULTS["adminMail", "DefaultValue"] = "peter.schmid@ieu.uzh.ch"
Sys.setenv(RUN_LONG_TEST=TRUE)
devtools::test()
