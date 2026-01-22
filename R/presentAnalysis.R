###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Which values to keep?
##' @description Returns TRUE for values that are considered as measured correctly and that are over a threshold, if one is used.
##' @param x a matrix with row and column names.
##' @param presentFlag a binary matrix with the same size as \code{x} which indicates if a values is considered as measured correctly.
##' @param param contains a threshold and whether it should be used.
##' \itemize{
##'  \item{sigThresh}{ the threshold...}
##'  \item{useSigThresh}{ ...and whether it should be used.}
##' }
##' @param isLog a logical indicating whether the \code{log2()} of the threshold should be used.
##' @return Returns a logical matrix.
##' @template roxygen-template
##' @examples
##' m1 = matrix(1:20,5)
##' rownames(m1) = letters[1:5]
##' colnames(m1) = letters[6:9]
##' ezPresentFlags(m1, param=list(useSigThresh=TRUE, sigThresh=10))
ezPresentFlags = function(x, presentFlag = NULL, param = NULL, isLog = FALSE) {
  isPresent = ezMatrix(TRUE, rows = rownames(x), cols = colnames(x))
  if (param$useSigThresh) {
    if (isLog) {
      thresh = log2(param$sigThresh)
    } else {
      thresh = param$sigThresh
    }
    isPresent[x < thresh] = FALSE
    isPresent[is.na(x)] = FALSE
  }
  if (!is.null(presentFlag)) {
    isPresent[!presentFlag] = FALSE
  }
  isPresent
}
