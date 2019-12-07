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

convertDiamondAnnotationToAbund <- function(annFile,feature){
  annDF <- read.delim(annFile, header = F, stringsAsFactors = F)
  foundIds <- annDF$V2
  if (feature == "function"){
  funcCatTemp <- sapply(foundIds, function(x) unlist(strsplit(x,"\\.1_"))[2], USE.NAMES = F)
  funcCatTemp <- sapply(funcCatTemp, function(x) unlist(strsplit(x,"_\\["))[1], USE.NAMES = F)
  funcCatTemp <- gsub("MULTISPECIES:_","",funcCatTemp)
  funcCat <- gsub(",_partial","",funcCatTemp)
  funcAbundDF <- data.frame(table(funcCat))
  rownames(funcAbundDF) <- funcAbundDF$funcCat
  funcAbundDF <- subset(funcAbundDF,select=-c(funcCat))
  names(funcAbundDF) <- gsub(".RefSeq.annotated.txt","",basename(annFile))
  return(data.frame(funcAbundDF, stringsAsFactors = F))
  } else if (feature == "organism"){
  orgTemp <- sapply(foundIds, function(x) unlist(strsplit(x,"_\\["))[2], USE.NAMES = F)
  orgTemp2 <- gsub("\\]","",orgTemp)
  org <- gsub("\\[","",orgTemp2)
  orgAbundDF <- data.frame(table(org))
  rownames(orgAbundDF) <- orgAbundDF$org
  orgAbundDF <- subset(orgAbundDF,select=-c(org))
  names(orgAbundDF) <- gsub(".RefSeq.annotated.txt","",basename(annFile))
  return(data.frame(orgAbundDF, stringsAsFactors = F))
  }
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title listOfAbundMerge
##' @description Create summary table from list of individual abundance.
##' @param  a list of data frame, one per sample.
##' @return Returns a data.frame.

listOfAbundMerge <- function(listOfAbund,names){
mergedDF <- Reduce(function(x, y) {
  z <- merge(x, y, all=TRUE,by=0)
  rownames(z) <- z$Row.names
  z <- subset(z, select=-c(Row.names))
  return(z)},
  listOfAbund)
names(mergedDF) <- names
mergedDF[is.na(mergedDF)] = 0
mergedDF <- mergedDF[apply(mergedDF,1,sd)>0,]
return(mergedDF)
}
