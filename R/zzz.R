###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


EZ_PARAM_DEFAULTS <<- ezRead.table(system.file("extdata/EZ_PARAM_DEFAULTS.txt", package="ezRun", mustWork = TRUE), comment.char="#")  


.onLoad = function(libname, pkgname){
  ## NICE_TO_HAVE: put global variables in a list ezGlobals  
  ## users should have their own global parameters
  if (!exists("EZ_GLOBAL_VARIABLES")){
    EZ_GLOBAL_VARIABLES <<- system.file("extdata/EZ_GLOBAL_VARIABLES.txt", package=pkgname, mustWork = TRUE)
  }
  if (file.exists(EZ_GLOBAL_VARIABLES)){
    ezWrite("loading EZ_GLOBAL_VARIABLES  from: ", EZ_GLOBAL_VARIABLES)
    source(EZ_GLOBAL_VARIABLES, local = .GlobalEnv)
  } else {
    ezWrite("EZ_GLOBAL_VARIABLES file defined but does not exist: ", EZ_GLOBAL_VARIABLES)    
  }
  # users should not overwrite the existing default parameter list; this could break apps in the package; they can extend the parameter list
  #EZ_PARAM_DEFAULTS <<- ezRead.table("~/R/ezRun/inst/extdata/EZ_PARAM_DEFAULTS.txt",  comment.char="#")  
  EZ_PARAM_DEFAULTS <<- ezRead.table(system.file("extdata/EZ_PARAM_DEFAULTS.txt", package=pkgname, mustWork = TRUE), comment.char="#")  
}



## DEFINE AN .onLoad Method here
## could load stuff from the subdirectory extdata
## see example from BSgenome:
# .onLoad <- function(libname, pkgname)
# {
#   if (pkgname != .pkgname)
#     stop("package name (", pkgname, ") is not ",
#          "the expected name (", .pkgname, ")")
#   extdata_dirpath <- system.file("extdata", package=pkgname,
#                                  lib.loc=libname, mustWork=TRUE)
#   
#   
