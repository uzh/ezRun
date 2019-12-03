
## settings for running the tests as Hubert
.libPaths("/srv/localdata/hubert/R-libs")
setwd("~/R/git/ezRun")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
require(ezRun)
EZ_PARAM_DEFAULTS["adminMail", "DefaultValue"] = "Hubert.Rehrauer@fgcz.ethz.ch"
Sys.setenv(RUN_LONG_TEST=TRUE)
devtools::test()


