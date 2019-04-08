###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch



##' @title Build parameter object
##' @description Parameters can be specified in the global defaults, app-specific defaults and user-given.
##' If a parameter is specified in multiple places, the user parameters override the app defaults which again
##' override the global defaults
##' @param userParam a list of parameters defined by the user
##' @param globalDefaults a data.frame containing the global defaults for parameters. Includes for each parameter a name, 
##' a default value, the type, and a description. The global defaults are read from a file. 
##' See the file EZ_PARAM_DEFAULTS.txt in the package directory.
##' @param appDefaults a data.frame containing application specific defaults.
##' Must have the same columns as the globalDefaults data.frame.
##' @template roxygen-template
##' @return Returns the merged list of parameters
##' @examples
##' library(ezRun)
##' head(EZ_PARAM_DEFAULTS)
##' globalDefaultTypedParams = ezParam()
##' modifiedParams = ezParam(list(ram=40))
##' modifiedParams$ram
##' globalDefaultTypedParams$ram
##' parseOptions("a=5 b='foo with space' c=foo")
## current implementation: every param needs a default value
ezParam = function(userParam=list(), globalDefaults=getGlobalDefaults(),
                   appDefaults=ezFrame()){

  specialParam = parseOptions(userParam$specialOptions)
  userParam[names(specialParam)] = specialParam
  defaults = rbind(globalDefaults[setdiff(rownames(globalDefaults), 
                                          rownames(appDefaults)), ], 
                   appDefaults)
  unknownParams = setdiff(names(userParam), rownames(defaults))
  sapply(unknownParams, function(x){message("unknown param: ", x)})
  for (nm in rownames(defaults)){
    if(!is.null(userParam[[nm]])){
      value <- userParam[[nm]]
    } else {
      value <- defaults[nm, "DefaultValue"]
    }
    userParam[[nm]] <- switch(defaults[nm, "Type"],
                              integer=as.integer(value),                         
                              numeric=as.numeric(value),
                              character=as.character(value),
                              charVector=if (length(value) > 1) value else unlist(strsplit(value, ",", fixed=TRUE)),
                              charList=if (length(value) > 1) value else parseListOptions(value),
                              logical=as.logical(value),
                              stop("unsupported type: ", defaults[nm, "Type"]))
  }
  
  ## avoid special characters in any option
  lapply(userParam, function(optString){
    if (class(optString) == "character" && any(grepl("[;\\{}$%#!*]", optString))){
      stop("special characters not allowed in option string: ", optString)
    }
  })
  
  # we build the ezRef object from hints in the general parameters
  if (is.null(userParam$ezRef)){
    userParam$ezRef = EzRef(userParam)
  }
  
  return(userParam)
}

getGlobalDefaults = function(){
  if (exists("EZ_PARAM_DEFAULTS")){
    return(EZ_PARAM_DEFAULTS)
  } else {
    ezRead.table(system.file("extdata/EZ_PARAM_DEFAULTS.txt", package="ezRun", mustWork = TRUE), comment.char="#")  
  }
}

# 'DC-like=Lgals3,Napsa B;cells=Cd79a,Ly6d;' to a named list
parseListOptions <- function(optString){
  params <- strsplit(optString, ";")[[1]]
  paramsList <- list()
  param <- params[1]
  for(param in params){
    param=strsplit(param, "=")[[1]]
    paramsList[[param[1]]] <- trimws(strsplit(param[2], ",")[[1]])
  }
  return(paramsList)
}

##' @describeIn ezParam Used to parse additional options specified in \code{userParam$specialOptions}.
##' Converts an option specification in the from "key1=value1 key2=value2" into a named list
parseOptions = function(optString){
  param = list()
  if (is.null(optString) || optString == ""){
    return(param)
  }
  if (length(grep("[[:cntrl:]]", optString))> 0){
    stop("control characters not allowed in option string")
  }
  if (any(grepl("[;\\{}$%#!*]", optString))){
    stop("special characters not allowed in option string")
  }
  params = unlist(strsplit(optString, " ", fixed=TRUE))
  quoteIdx = grep("\'", params)
  for (i in 1:length(params)){
    if (ezGrepl("=\'", params[i])){
      j = min(quoteIdx[quoteIdx > i])
      message(i, "- ", j, " ", params[i])
      params[i] = gsub("\'", "", paste(params[i:j], collapse=" "))
      params[(i+1):j] = "" 
    }
  }
  params = params[params != ""]
  for (p in params){
    pv = unlist(strsplit(p, "=", fixed=TRUE))
    if (length(pv) == 2){
      name = pv[1]
      if (is.na(suppressWarnings(as.numeric(pv[2])))){
        value = pv[2]
      } else {
        value = as.numeric(pv[2])
      }
      cat("found parameter in name: ", name, " = ", value, "\n")
      param[[name]] = value
    } else {
      stop("invalid syntax in option string: ", p)
    }
  }
  return(param)
}

##' @title Check if a value is specified
##' @description A value is specified if it is not \code{NULL}, it is not an empty vector, empty list or empty string.
##' If the value is a vector, the first element must be different from the empty string.
##' @param x usually a parameter to check.
##' @template roxygen-template
##' @return Returns FALSE or TRUE.
##' @examples
##' ezIsSpecified(5)
##' ezIsSpecified(c("this","is"))
##' ezIsSpecified(c("","this isn't"))
ezIsSpecified = function(x){
  !is.null(x) && length(x) > 0 && x[1] != "" && !is.na(x) && x[1] != "NA"
}

##' @title Wrapper for data.frame suitable for data processing rather than stastistical modelling
##' @description The original \code{data.frame()} is suitable for statistical modelling with factors etc.
##' Our version is designed to be more suitable for data processing, i.e. we
##' * do keep column names the way they are; we do not require that they represent valid variable names
##' * by default we do not convert strings to factors
##' There exist, similarities with data_frame 
##' in Hadleys package. The difference is that we encourage to use rownames while he discourages this.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun data.frame()
##' @return Returns a data.frame.
##' @examples
##' ezFrame(first=1:3, second=5, "with space"="text", row.names=letters[1:3])
ezFrame = function(..., row.names = NULL, check.rows=TRUE,
                   check.names=FALSE,
                   stringsAsFactors=FALSE){
  x = data.frame(..., check.rows=check.rows, check.names=check.names, stringsAsFactors=stringsAsFactors)
  if (!is.null(row.names)){
    rownames(x) = row.names
  }
  return(x)
}
