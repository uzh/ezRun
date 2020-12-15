###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets a strand value
##' @description Gets a strand value of plus or minus from an integer.
##' @param x an integer to select the strand value: 1 equals plus and 2 equals minus.
##' @template roxygen-template
##' @return Returns a named character "+" or "-".
##' @examples 
##' strandValue(1)
strandValue = function(x){
  c(plus="+", minus="-")[x]
}

##' @describeIn strandValue Does the same, but returns the character as "plus" or "minus".
strandName = function(x){
  c("+"="plus", "-"="minus")[x]
}

##' @title Fixes the strand information
##' @description Fixes the strand information by either leaving it as it is, flipping it or setting both.
##' @param strandValues a set of strand values + - and *.
##' @param strandMode a character defining what should be done with the values:
##' \itemize{
##'   \item{sense}{ returns \code{strandValues} without modifying them.}
##'   \item{antisense}{ returns the flipped \code{strandValues}.}
##'   \item{both}{ sets all the \code{strandValues} to "*".}
##' }
##' @template roxygen-template
##' @return Returns a character vector or Rle object containing the strand information.
##' @examples
##' strandValues = c("+","+","-","*")
##' fixStrand(strandValues, "sense")
##' flipStrand(strandValues)
fixStrand = function(strandValues, strandMode="sense"){
  return(switch(strandMode,
         sense=return(strandValues),
         antisense=flipStrand(strandValues),
         both=Rle(factor("*", levels=c("+", "-", "*")), length(strandValues))))
}

##' @describeIn fixStrand This function performs the flipping of the strand values.
flipStrand = function(strandValues){
  strandMap = c("+"="-", "-"="+", "*"="*")
  if (class(strandValues) == "Rle"){
    return(Rle(factor(as.character(strandMap[runValue(strandValues)]), levels=c("+", "-", "*")), runLength(strandValues)))
  } else {
    return(strandMap[match(strandValues, names(strandMap))])
  }
}

## strand option for the tuxedo suite
##' @title Gets the tuxedo library type
##' @description Gets the tuxedo library type from the strand mode.
##' @param strandMode a character defining the strand mode and what to return.
##' \itemize{
##'   \item{sense}{ returns "fr-secondstrand".}
##'   \item{antisense}{ returns "fr-firststrand".}
##'   \item{both}{ returns "fr-unstranded".}
##' }
##' @template roxygen-template
##' @return Returns a character.
##' @examples 
##' getTuxedoLibraryType("sense")
getTuxedoLibraryType = function(strandMode){
  return(switch(strandMode,
         "both"="fr-unstranded",
         "sense"="fr-secondstrand",
         "antisense"="fr-firststrand",
         stop(paste("bad strandMode: ", strandMode))))
}

##' @title Is \code{x} a valid cigar?
##' @description Tests whether \code{x} is a valid cigar by checking the format of the character.
##' @param x a character vector to check.
##' @template roxygen-template
##' @return Returns a logical vector answering which elements are valid cigars.
##' @seealso \code{\link[GenomicAlignments]{cigarWidthAlongReferenceSpace}}
##' @examples
##' isValidCigar("3M5G")
isValidCigar = function(x){
  sapply(x, function(y){!inherits(try(cigarWidthAlongReferenceSpace(y), silent=TRUE), "try-error")})
}

##' @title Shifts zeros
##' @description Shifts zeros by changing values lower than \code{minSignal} to a random numeric in the range of (1/4 to 3/4) * \code{minSignal}. This is done to prevent zero signals to stack at 0.
##' @param counts a set of numeric or integer values.
##' @param minSignal a numeric or integer specifying the minimal signal amount.
##' @template roxygen-template
##' @return Returns the modified count values.
##' @examples 
##' shiftZeros(1:10, 5)
shiftZeros = function(counts, minSignal){
  isLow = counts < minSignal
  isLow[is.na(isLow)] = TRUE
  counts[isLow] = runif(sum(isLow), min=0.25 * minSignal, max=0.75 * minSignal)
  return(counts)
}
