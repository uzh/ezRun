###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Mothur Summary Table
##' @description Create summary table from mothur summary files.
##' @param  summary, mothur summary files.
##' @return Returns a data.frame.

convertDiamondAnnotationToAbund <- function(annFile){
  annDF <- read.delim(annFile, header = F, stringsAsFactors = F)
  foundIds <- annDF$V2
  funcCatTemp <- sapply(foundIds, function(x) unlist(strsplit(x,"\\.1_"))[2], USE.NAMES = F)
  funcCatTemp <- sapply(funcCatTemp, function(x) unlist(strsplit(x,"_\\["))[1], USE.NAMES = F)
  funcCatTemp <- gsub("MULTISPECIES:_","",funcCatTemp)
  funcCat <- gsub(",_partial","",funcCatTemp)
  orgTemp <- sapply(foundIds, function(x) unlist(strsplit(x,"_\\["))[2], USE.NAMES = F)
  org <- gsub("\\]","",orgTemp)
  funcAbundDF <- data.frame(table(funcCat))
  orgAbundDF <- data.frame(table(org))
  return(list(funcAbundDF=funcAbundDF,orgAbundDF=orgAbundDF))
}

