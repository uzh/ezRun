###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Deseq2
##' @seealso \code{\link{EzAppDeseq2}}
ezMethodDeseq2 = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  param$testMethod = "deseq2"
  param$normMethod = ""
  param$runGfold = FALSE ## no need to compute moderated ratios; deseq2 does this already
  if (!is.null(param$markOutliers) && param$markOutliers){
    stop("DESeq2 does not support marking outliers because marked outliers would still be used in dispersion estimates")
  }
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(basename(output$getColumn("Report")))
  stopifnot(param$sampleGroup != param$refGroup)
  dataset = input$meta
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset, param), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  }
  input$meta = dataset
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  
  modifiedParams = modifyParameters(dataset, param)
  param$grouping = modifiedParams$grouping
  param$batch = modifiedParams$batch
  
  result = twoGroupCountComparison(rawData, param)
  if (isError(result)){
    writeErrorReport(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  result$featureLevel = rawData$featureLevel
  result$countName = rawData$countName
  
  writeNgsTwoGroupReport(dataset, result, htmlFile, param=param, rawData=rawData)
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
                  "Initializes the application using its specific defaults."
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
