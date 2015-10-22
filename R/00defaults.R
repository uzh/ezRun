###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets and processes all parameters
##' @description Gets and processes all parameters. Removes any duplicates and recognizes the correct datatype.
##' @param userParam a list of parameters defined by the user.
##' @param globalDefaults global package defaults. These should normally not be changed.
##' @param appDefaults a data.frame containing application specific defaults.
##' @template roxygen-template
##' @return Returns the complete list of parameters used by the package.
##' @examples
##' ezParam()
## current implementation: every param needs a default value
ezParam = function(userParam=list(), globalDefaults=EZ_PARAM_DEFAULTS, appDefaults=ezFrame()){

  specialParam = parseOptions(userParam$specialOptions)
  userParam[names(specialParam)] = specialParam
  defaults = rbind(globalDefaults[setdiff(rownames(globalDefaults), rownames(appDefaults)), ], appDefaults)
  unknownParams = setdiff(names(userParam), rownames(defaults))
  sapply(unknownParams, function(x){message("unknown param: ", x)})
  for (nm in rownames(defaults)){
    value = ifelse(!is.null(userParam[[nm]]), userParam[[nm]], defaults[nm, "DefaultValue"])
    userParam[[nm]] = switch(defaults[nm, "Type"],
                         integer=as.integer(value),                         
                         numeric=as.numeric(value),
                         character=as.character(value),
                         charVector=unlist(strsplit(value, ",", fixed=TRUE)),
                         logical=as.logical(value),
                         stop("unsupported type: ", defaults[nm, "Type"]))
  }
  if (is.null(userParam$ezRef)){
    userParam$ezRef = EzRef(userParam)
  }
  return(userParam)
}

##' @title Is x specified?
##' @description Checks whether \code{x} is an existing parameter and whether its first entry not an empty character.
##' @param x usually a parameter to check.
##' @template roxygen-template
##' @return Returns FALSE or TRUE.
##' @examples
##' ezIsSpecified(c("this","is"))
##' ezIsSpecified(c("","this isn't"))
ezIsSpecified = function(x){
  !is.null(x) && length(x) > 0 && x[1] != ""
}

##' @title Modified default of data.frame
##' @description Modified version of \code{data.frame()} with a different default.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun data.frame
##' @return Returns a data.frame.
##' @examples
##' ezFrame(a=1:10,b=5)
ezFrame = function(..., row.names = NULL, check.rows=TRUE,
                   check.names=FALSE,
                   stringsAsFactors=FALSE){
  data.frame(..., row.names=row.names, check.rows=check.rows, check.names=check.names, stringsAsFactors=stringsAsFactors)
}

##' @describeIn ezParam Used to parse additional options specified in \code{userParam$specialOptions}.
parseOptions = function(optString){
  param = list()
  if (is.null(optString) || optString == ""){
    return(param)
  }
  if (length(grep("[[:cntrl:]]", optString))> 0){
    stop("control characters not allowed in option string")
  }
  if (length(grep("[;\\(){}$%:#!?*]", optString))> 0){
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
