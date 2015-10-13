###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Deseq2
##' @seealso \code{\link{EzAppDeseq2}}
ezMethodDeseq2 = function(input=NA, output=NA, param=NA){
  param$testMethod = "deseq2"
  param$normMethod = ""
  param$runGfold = FALSE ## no need to compute moderated ratios; deseq2 does this already
  if (!is.null(param$markOutliers) && param$markOutliers){
    stop("DESeq2 does not support marking outliers because marked outliers would still be used in dispersion estimates")
  }
  ngsTwoGroupAnalysis(input, output, param=param)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDeseq2()
##' @seealso \code{\link{ezMethodDeseq2}}
EzAppDeseq2 <-
  setRefClass("EzAppDeseq2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  runMethod <<- ezMethodDeseq2
                  name <<- "EzAppDeseq2"
                }
              )
  )


## NOTEP: gets called from edgerTwoGroups. the Deseq2 method calls ngsTwoGroupAnalysis() from edgerTwoGroups,
## maybe it makes sense to combine these apps in one R file.
##' @describeIn twoGroupCountComparison Runs the Deseq2 test method.
runDeseq2 = function(x, sampleGroup, refGroup, grouping, batch=NULL, isPresent=NULL){
  library(DESeq2, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  if (ezIsSpecified(batch)){
    if (!is.numeric(batch)){
      batch = as.factor(batch)
    } else {
      message("using numeric batch factor")
    }
    colData = data.frame(grouping=as.factor(grouping), batch=batch, row.names=colnames(x))
    dds = DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping + batch)
  } else {
    colData = data.frame(grouping=as.factor(grouping), row.names=colnames(x))
    dds = DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping)
  }
  dds = estimateSizeFactors(dds, controlGenes=isPresent)
  dds = DESeq(dds, quiet=FALSE)
  res = results(dds, contrast=c("grouping", sampleGroup, refGroup), cooksCutoff=FALSE)
  res=as.list(res)
  res$sf = 1/colData(dds)$sizeFactor
  return(res)
}
