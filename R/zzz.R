###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

.onLoad = function(libname, pkgname){

  EZ_PARAM_DEFAULTS <<- ezRead.table(system.file("extdata/EZ_PARAM_DEFAULTS.txt", 
                                                 package=pkgname, mustWork = TRUE), comment.char="#")  
}

if (!exists("EZ_GLOBAL_VARIABLES")){
  EZ_GLOBAL_VARIABLES <<- system.file("extdata/EZ_GLOBAL_VARIABLES.txt",
                                      package="ezRun", mustWork = TRUE)
}
if (file.exists(EZ_GLOBAL_VARIABLES)){
  ezWrite("loading EZ_GLOBAL_VARIABLES  from: ", EZ_GLOBAL_VARIABLES)
  source(EZ_GLOBAL_VARIABLES, local = TRUE)
} else {
  ezWrite("EZ_GLOBAL_VARIABLES file defined but does not exist: ", EZ_GLOBAL_VARIABLES)    
}
