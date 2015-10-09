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
## maybe it makes sense to combine these apps (and perhaps also edgerMulti) in one R file.
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
  
## NOTEP: seems unused
runDeseq = function(x, sampleGroup, refGroup, grouping){
  
  library(DESeq, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  cds = newCountDataSet( x, grouping)
  #sf = ezLogmeanScalingFactor(x, presentFlag=x > 0)
  #sizeFactors(cds) = 1/sf
  cds = estimateSizeFactors(cds)
  sf = 1/sizeFactors(cds)
  if (min(table(grouping)) > 1){
    cds = estimateDispersions( cds , method="pooled",sharingMode="maximum")
  } else {
    cds = estimateDispersions( cds , method="blind",sharingMode="fit-only")
  }
  
  res = nbinomTest(cds, refGroup, sampleGroup)
  res = as.list(res)
  #bmvA <- getBaseMeansAndVariances( counts(cds)[,colsA], sizeFactors(cds)[colsA] )
  #bmvB <- getBaseMeansAndVariances( counts(cds)[,colsB], sizeFactors(cds)[colsB] )
  
  #res$log2Expr = log2(res$baseMean)
  ## do not return groupMeans now!
  #res$groupMeans = cbind(res$baseMeanB, res$baseMeanA)
  #colnames(res$groupMeans) = c(sampleGroup, refGroup)
  #rownames(res$groupMeans) = res$id
  #res$resVar = cbind(res$resVarB, res$resVarA)
  #colnames(res$resVar) = c(sampleGroup, refGroup)
  #rownames(res$resVar) = res$id
  res$sf = sf
  return(res)
}

