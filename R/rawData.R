###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

getSignal = function(rawData){
  if (metadata(rawData)$isLog){
    return(2^assays(rawData)$signal)
  } else {
    return(assays(rawData)$signal)
  }
}

getRpkm <- function(rawData){
  require(Matrix)
  require(SummarizedExperiment)
  require(SingleCellExperiment)
  # if(is(rawData, "SingleCellExperiment")){
  #   require(scater)
  #   rpkm <- calculateFPKM(rawData, effective_length=rowData(rawData)$featWidth)
  # }else{
    # for bulk, but also valid for single cell
  # scater implementation fails when no counts from one cell
    libSize <- Matrix::colSums(assays(rawData)$counts)
    rpkm <- sweep(assays(rawData)$counts * 1e9, MARGIN=1,
                  STATS=rowData(rawData)$featWidth, FUN="/")
    rpkm <- sweep(rpkm, MARGIN=2, STATS=libSize, FUN="/")
  # }
  return(rpkm)
}

getTpm <- function(rawData){
  require(SummarizedExperiment)
  require(SingleCellExperiment)
  # if(is(rawData, "SingleCellExperiment")){
  #   # scater's implementation seems to be faster.   
  #   if(!ezIsSpecified(metadata(rawData)$param$scProtocol)){
  #     stop("scProtocol must be specified in param.")
  #   }
  if(ezIsSpecified(metadata(rawData)$param$scProtocol) && toupper(metadata(rawData)$param$scProtocol) == "10X"){
    require(scater)
    tpm <- calculateTPM(rawData, effective_length=NULL)
    # }else if(metadata(rawData)$param$scProtocol == "smart-Seq2"){
    #   tpm <- calculateTPM(rawData, effective_length=rowData(rawData)$featWidth)
    # }
  }else{
    # for bulk, but also valid for single cell
    # scater implementation fails when no counts from one cell
    tpm <- sweep(assays(rawData)$counts * 1e3, MARGIN=1,
                 STATS=rowData(rawData)$featWidth, FUN="/")
    tpm <- sweep(tpm * 1e6, MARGIN=2, STATS=Matrix::colSums(tpm),
                 FUN="/")
  }
  return(tpm)
}

getCpm <- function(rawData){
  require(SummarizedExperiment)
  require(SingleCellExperiment)
  # if(is(rawData, "SingleCellExperiment")){
  #   # scater's implementation seems to be faster.
  #   require(scater)
  #   ans <- calculateCPM(rawData, exprs_values = "counts",
  #                       use_size_factors = TRUE)
  # }else{
    # for bulk, but also valid for single cell
    ans <- sweep(assays(rawData)$counts*1e6, MARGIN=2, 
                 STATS=Matrix::colSums(assays(rawData)$counts),
                 FUN="/")
  # }
  return(ans)
}

aggregateCountsByGene <- function(rawData){
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

## a feature will typically be a gene, isoform, or a microarray probe
## all matrices/data.frames in rawData must have rownames!
##' @title Selects features from raw data
##' @description Selects features from raw data by looking at "signal" (if specified) or "counts" from \code{rawData}.
##' @template rawData-template
##' @param keep a character vector specifying which signals or counts to keep.
##' @template roxygen-template
##' @return Returns the selected raw data.
##' @examples
##' param = ezParam()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' file = system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file, dataRoot=param$dataRoot)
##' rawData = loadCountDataset(input, param)
##' keep = "YFR014C"
##' rd2 = selectFeatures(rawData, keep)
selectFeatures = function(rawData, keep){
    
    if (!is.null(rawData$signal)){
        stopifnot(keep %in% rownames(rawData$signal))
        idx = match(keep, rownames(rawData$signal))
    } else {
        stopifnot(keep %in% rownames(rawData$counts))
        idx = match(keep, rownames(rawData$counts))
    }
    for (nm in names(rawData)){
        if (is.matrix(rawData[[nm]]) || is.data.frame(rawData[[nm]])){
            if (all(keep %in% rownames(rawData[[nm]]))){
                rawData[[nm]] = rawData[[nm]][idx, , drop=FALSE]      
            }
        }
    }
    return(rawData)
}

