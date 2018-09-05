###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Gets the signal
##' @description Gets the signal from raw data and modifies it if necessary due to logarithm.
##' @template rawData-template
##' @template roxygen-template
##' @return Returns the signal.
##' @examples
##' param = ezParam()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' file = system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file, dataRoot=param$dataRoot)
##' rawData = loadCountDataset(input, param)
##' sig = getSignal(rawData)
## TODOEXAMPLE: get example with proper return and improve description.
getSignal = function(rawData){
  if (metadata(rawData)$isLog){
    return(2^assays(rawData)$signal)
  } else {
    return(assays(rawData)$signal)
  }
}

##' @describeIn getSignal Does the same but returns the log2 instead.
getLog2Signal = function(rawData){
  if (rawData$isLog){
    return(rawData$signal)
  } else {
    return(log2(rawData$signal))
  }
}

getRpkm = function(rawData){
  require(Matrix)
  require(SummarizedExperiment)
  
  if (!is.null(assays(rawData)$rpkm)){
    return(assays(rawData)$rpkm)
  }
  if (is.null(rowData(rawData)$width)){
    if(is.null(rowData(rawData)$featWidth)){
      warning("The `width` is not available in annotation.")
      return(NULL)
    }else{
      rowData(rawData)$width <- rowData(rawData)$featWidth
    }
  }
  libSize = Matrix::colSums(assays(rawData)$counts)
  # rpkm <- edgeR::rpkm(assays(rawData)$counts, gene.length=rowData(rawData)$width,
  #                     normalized.lib.sizes=FALSE, log=FALSE, prior.count=0)
  # Didn't test whether edgeR::rpkm works on sparse matrxi from single cell
  rpkm <- sweep(assays(rawData)$counts * 1e9, MARGIN=1,
                STATS=rowData(rawData)$width, FUN="/")
  rpkm <- sweep(rpkm, MARGIN=2, STATS=libSize, FUN="/")
  return(rpkm)
}

getTpm = function(rawData) {
  
  if (!is.null(assays(rawData)$tpm)){
    return(assays(rawData)$tpm)
  }
  if (is.null(rowData(rawData)$width)){
    if(is.null(rowData(rawData)$featWidth)){
      warning("The `width` is not available in annotation.")
      return(NULL)
    }else{
      rowData(rawData)$width <- rowData(rawData)$featWidth
    }
  }
  tpm = assays(rawData)$counts
  for (i in 1:ncol(tpm)){
    rpk = (assays(rawData)$counts[,i] * 1e3) /(rowData(rawData)$width)
    scalingFactor = sum(rpk)/1e6
    tpm[, i] = rpk / scalingFactor
  }
  return(tpm)
}


##' @title Aggregates counts by gene
##' @description Aggregates counts by gene.
##' @param param a list of parameters:
##' \itemize{
##'   \item{geneColumnSet}{ a character naming the columns.}
##'   \item{sigThresh}{ the threshold...}
##'   \item{useSigThresh}{ ...and whether it should be used.}
##' }
##' @template rawData-template
##' @template roxygen-template
##' @seealso \code{\link{getGeneMapping}}
##' @return Returns the aggregates raw data.
aggregateCountsByGene = function(param, rawData){
  
  genes = getGeneMapping(param, rawData$seqAnno)
  
  if (is.null(genes)){
    return(list(error=paste("gene summaries requested but not gene column available. did you specify the build?<br>column names tried:<br>",
                            paste(param$geneColumnSet, collapse="<br>"))))
  }
  seqAnnoNew = data.frame(row.names=na.omit(unique(genes)))
  #   commonCols = intersect(c("gene_id", "type", "seqid", "strand", "GO BP", "GO MF", "GO CC"), 
  #                          colnames(rawData$seqAnno))
  for (nm in colnames(rawData$seqAnno)){
    seqAnnoNew[[nm]] = tapply(rawData$seqAnno[[nm]], genes, 
                              ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE, na.rm=TRUE)[rownames(seqAnnoNew)]
  }
  ## special merging for special columns
  if (!is.null(rawData$seqAnno$start)){
    gStart = tapply(as.integer(sub(";.*", "", rawData$seqAnno$start)), genes, min)
    seqAnnoNew$start = gStart[rownames(seqAnnoNew)]
  }
  if (!is.null(rawData$seqAnno$end)){
    gEnd = tapply(as.integer(sub(".*;", "", rawData$seqAnno$end)), genes, max)
    seqAnnoNew$end = gEnd[rownames(seqAnnoNew)]
  }
  if (!is.null(rawData$seqAnno$width)){
    seqAnnoNew$width = tapply(rawData$seqAnno$width, genes, mean)[rownames(seqAnnoNew)]
  }
  if (!is.null(rawData$seqAnno$gc)){
    seqAnnoNew$gc = tapply(as.numeric(rawData$seqAnno$gc), genes, mean)[rownames(seqAnnoNew)]
  }
  rawData$seqAnno = seqAnnoNew
  if (rawData$isLog){
    stop("not supported")
  }
  for (nm in c("counts", "countsStart", "countsEnd", "rpkm")){
    if (!is.null(rawData[[nm]])){
      rawData[[nm]] = as.matrix(averageRows(rawData[[nm]], genes, func=sum))[rownames(seqAnnoNew), ]
    }
  }
  if (!is.null(rawData[["presentFlag"]])){
    if (param$useSigThresh){
      rawData[["presentFlag"]] = rawData[["counts"]] > param$sigThresh
    } else {
      rawData[["presentFlag"]] = rawData[["counts"]] > 0     
    }
  }
  
  rawData$signal = NULL
  rawData$featureLevel = "gene"
  return(rawData)
}

aggregateCountsByGeneSE <- function(rawData){
  require(SummarizedExperiment)
  param <- metadata(rawData)$param
  seqAnno <- data.frame(rowData(rawData), row.names=rownames(rawData),
                        check.names = FALSE, stringsAsFactors=FALSE)
  genes = getGeneMapping(param, seqAnno)
  
  if (is.null(genes)){
    return(list(error=paste("gene summaries requested but not gene column available. did you specify the build?<br>column names tried:<br>",
                            paste(param$geneColumnSet, collapse="<br>"))))
  }
  
  seqAnnoNew <- aggregateFeatAnno(seqAnno)
  
  if (metadata(rawData)$isLog){
    stop("Counts in logarithm are not supported!")
  }
  
  newRawCounts <- SimpleList()
  for (nm in setdiff(names(assays(rawData)), "presentFlag")){
      newRawCounts[[nm]] = as.matrix(averageRows(assays(rawData)[[nm]], genes, 
                                       func=sum))[rownames(seqAnnoNew), ]
  }
  if (param$useSigThresh){
    newRawCounts[["presentFlag"]] = newRawCounts[["counts"]] > param$sigThresh
  } else {
    newRawCounts[["presentFlag"]] = newRawCounts[["counts"]] > 0
  }
  
  ## Get rid of signal matrix in case it exists
  ## signal will not be valid after the aggregation
  newRawCounts[["signal"]] <- NULL
    
  newRawData <- SummarizedExperiment(
    assays=newRawCounts,
    rowData=seqAnnoNew, colData=colData(rawData),
    metadata=list(isLog=FALSE, featureLevel="gene",
                  type="Counts", countName=metadata(rawData)$countName,
                  param=param)
  )
  return(newRawData)
}
