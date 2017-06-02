###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Merges raw data
##' @description Merges raw data by combining two lists of raw data into one.
##' @param rawData1 a list of raw data. This is the main dataset where the information from the second gets added to.
##' @param rawData2 the second list of raw data. Both are obtained from \code{loadCountDataset()}.
##' @template roxygen-template
##' @return Returns a list containing the merged raw data.
##' @examples
##' param = ezParam()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' file = system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file, dataRoot=param$dataRoot)
##' rawData1 = loadCountDataset(input$subset(1), param)
##' rawData2 = loadCountDataset(input$subset(2), param)
##' rdm = mergeRawData(rawData1, rawData2)
mergeRawData = function(rawData1, rawData2){
  for (nm in names(rawData1)){
    if (is.matrix(rawData1[[nm]]) || is.data.frame(rawData1[[nm]])){
      if (setequal(rownames(rawData1[[nm]]), rownames(rawData2[[nm]]))){
        rawData1[[nm]] = cbind(rawData1[[nm]], rawData2[[nm]][rownames(rawData1[[nm]]), ,drop=FALSE])
      } else {
        if (setequal(colnames(rawData1[[nm]]), colnames(rawData2[[nm]]))){
          rawData1[[nm]] = rbind(rawData1[[nm]], rawData2[[nm]][ , colnames(rawData1[[nm]]), drop=FALSE])
        }
      }
    }
  }
  return(rawData1)
}

## a sample can figure as a  colum name (e.g. signal) and as a rowname (e.g. dataset)
##' @title Select samples from raw data
##' @description Select samples from raw data by specifying either the signal names or the sample names.
##' @template rawData-template
##' @param samples a character vector specifying which samples to select.
##' @template roxygen-template
##' @return Returns a list of the selected samples.
##' @examples
##' param = ezParam()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' file = system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file, dataRoot=param$dataRoot)
##' rawData = loadCountDataset(input, param)
##' rd2 = selectSamples(rawData, c("wt_1","wt_2"))
selectSamples = function(rawData, samples){
  if (!is.null(rawData$signal)){
    stopifnot(samples %in% colnames(rawData$signal))
  } else {
    stopifnot(samples %in% colnames(rawData$counts))
  }
  for (nm in names(rawData)){
    if (is.matrix(rawData[[nm]]) || is.data.frame(rawData[[nm]])){
      if (all(samples %in% colnames(rawData[[nm]]))){
        rawData[[nm]] = rawData[[nm]][, samples, drop=FALSE]
      }
      if (all(samples %in% rownames(rawData[[nm]]))){
        rawData[[nm]] = rawData[[nm]][samples, , drop=FALSE]
      }
    }
  }
  return(rawData)
}


renameSamples = function(rawData, samplesNew){
  samplesOld = rownames(rawData$dataset)
  stopifnot(length(samplesOld) == length(samplesNew))
  for (nm in names(rawData)){
    if (is.matrix(rawData[[nm]]) || is.data.frame(rawData[[nm]])){
      if (all(samplesOld %in% colnames(rawData[[nm]]))){
        colnames(rawData[[nm]]) = samplesNew
      }
      if (all(samplesOld %in% rownames(rawData[[nm]]))){
        rownames(rawData[[nm]]) = samplesNew
      }
    }
  }
  return(rawData)
}


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
  if (rawData$isLog){
    return(2^rawData$signal)
  } else {
    return(rawData$signal)
  }
}

getSignalSE = function(rawData){
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


##' @title Gets the rpkm measurement
##' @description Gets the rpkm measurement.
##' @template rawData-template
##' @template roxygen-template
##' @return Returns the rpkm measurement.
##' @examples
##' param = ezParam()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' file = system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file, dataRoot=param$dataRoot)
##' rawData = loadCountDataset(input, param)
##' rpkm = getRpkm(rawData)
## TODOEXAMPLE: get example with proper return
getRpkm = function(rawData){
  #edgeR::rpkm.default
  #edgeR::cpm.default
  if (!is.null(rawData$rpkm)){
    return(rawData$rpkm)
  }
  if (is.null(rawData$seqAnno$width)){
    return(NULL)
  }
  libSize = colSums(rawData$counts)
  rpkm = rawData$counts
  for (i in 1:ncol(rpkm)){
    rpkm[, i] = (rawData$counts[,i] * 1e9) /(rawData$seqAnno$width * libSize[i])
    #rpkm[, i] = rawData$counts[,i] /(libSize[i]* 1e6) /(rawData$seqAnno$width / 1e3)
  }
  return(rpkm)
}

getRpkmSE = function(rawData){
  #edgeR::rpkm.default
  #edgeR::cpm.default
  if (!is.null(assays(rawData)$rpkm)){
    return(assays(rawData)$rpkm)
  }
  if (is.null(rowData(rawData)$width)){
    warning("The `width` is not available in annotation.")
    return(NULL)
  }
  libSize = colSums(assays(rawData)$counts)
  rpkm = assays(rawData)$counts
  for (i in 1:ncol(rpkm)){
    rpkm[, i] = (assays(rawData)$counts[,i] * 1e9) /(rowData(rawData)$width * libSize[i])
    #rpkm[, i] = rawData$counts[,i] /(libSize[i]* 1e6) /(rawData$seqAnno$width / 1e3)
  }
  return(rpkm)
}


##' @title Gets the tpm measurement
##' @description Gets the transcripts per million measurement.
##' @template rawData-template
##' @template roxygen-template
##' @return Returns the tpm measurement.
getTpm = function(rawData) {
  if (!is.null(rawData$tpm)){
    return(rawData$tpm)
  }
  if (is.null(rawData$seqAnno$width)){
    return(NULL)
  }
  tpm = rawData$counts
  for (i in 1:ncol(tpm)){
    rpk = (rawData$counts[,i] * 1e3) /(rawData$seqAnno$width)
    scalingFactor = sum(rpk)/1e6
    tpm[, i] = rpk / scalingFactor
  }
  return(tpm)
}

getTpmSE = function(rawData) {
  if (!is.null(assays(rawData)$tpm)){
    return(assays(rawData)$tpm)
  }
  if (is.null(rowData(rawData)$width)){
    warning("The `width` is not available in annotation.")
    return(NULL)
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

aggregateCountsByGeneSE <- function(param, rawData){
  require(SummarizedExperiment)
  genes = getGeneMapping(param, rowData(rawData))
  
  if (is.null(genes)){
    return(list(error=paste("gene summaries requested but not gene column available. did you specify the build?<br>column names tried:<br>",
                            paste(param$geneColumnSet, collapse="<br>"))))
  }
  
  ## `GO BP` will be replaced by `GO.BP` after as.data.frame
  seqAnnoNew <- setNames(as.data.frame(rowData(rawData)),
                         colnames(rowData(rawData)))
  seqAnnoNew <- aggregateFeatAnno(seqAnnoNew)
  
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
                  type="Counts", countName=metadata(rawData)$countName)
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
