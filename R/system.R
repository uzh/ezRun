###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Invokes a system command
##' @description Invokes the system command specified by \code{cmd}.
##' @param cmd a system command input as a character.
##' @param echo a logical defining whether \code{cmd} should be written by \code{ezWrite}.
##' @param intern a logical which indicates whether to capture the output of the command as an R character vector.
##' @param stopOnFailure defaults to the opposite of intern and specifies whether the command should stop.
##' @return Returns input as character vector if \code{echo} is set to TRUE. Returns an error code if \code{intern} is set to FALSE.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun system()
##' @seealso \code{\link[base]{system}}
##' @examples
##' try(ezSystem("who"))
# bash -c 'set -o pipefail; ls -a'
ezSystem = function(cmd, echo=TRUE, intern=FALSE, stopOnFailure=!intern, ...){
  if (echo){
    ezWrite(paste("EXECUTED CMD:", cmd))
  }
  if (grepl("|", cmd, fixed=TRUE) | grepl(";", cmd, fixed=TRUE)){
    if(grepl("'", cmd)){
      stop(paste("single quotes not supported in command if there is a | or a ; : ", cmd))
    }
    res = system(paste("bash -c 'set -e; set -o pipefail; ", cmd, "'"), intern=intern, ...)
  } else {
    res = system(cmd, intern=intern, ...)    
  }
  #res = system(paste("set -o pipefail; ", cmd), intern=intern, ...)
  if (stopOnFailure){
    if (res != 0){
      stop(paste(cmd, "\n", "failed"))
    }
  }
  return(res)
}

##' @title Determines the number of CPU cores to be used
##' @description The function will try to find the amount of CPU cores on the current host or how many should be used. First, it will try \code{getOption("cores")}, then \code{Sys.getenv("NSLOTS")} and lastly \code{detectCores(logical=TRUE)}.
##' @return Returns an integer specifying the amount of CPU cores R is using.
##' @template roxygen-template
##' @seealso \code{\link[parallel]{detectCores}}
##' @examples
##' ezThreads()
ezThreads = function(){
  coreOpt = getOption("cores")
  if (!is.null(coreOpt)){
    return(as.integer(coreOpt))
  }
  nslots = as.integer(Sys.getenv("NSLOTS"))
  if (!is.na(nslots)){
    return(nslots)
  }
  nslots = parallel::detectCores(logical=FALSE)
  return(nslots)
}

