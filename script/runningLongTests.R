
## settings for running the tests as Hubert
.libPaths("/srv/GT/analysis/course_sushi/lib")
setwd("~/R/git/ezRun")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
require(ezRun)
EZ_PARAM_DEFAULTS["adminMail", "DefaultValue"] = "Hubert.Rehrauer@fgcz.ethz.ch"
Sys.setenv(RUN_LONG_TEST=TRUE)
devtools::test()


## settings for running the tests as Peter
.libPaths("/home/pdschmid/R-libs")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
setwd("/home/pdschmid/Git/ezRun")
require(ezRun)
EZ_PARAM_DEFAULTS["adminMail", "DefaultValue"] = "peter.schmid@ieu.uzh.ch"
Sys.setenv(RUN_LONG_TEST=TRUE)
devtools::test()
