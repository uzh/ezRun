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

getRpkm = function(rawData){
  require(Matrix)
  require(SummarizedExperiment)
  
  if (!is.null(assays(rawData)$rpkm)){
    return(assays(rawData)$rpkm)
  }
  libSize = Matrix::colSums(assays(rawData)$counts)
  # Didn't test whether edgeR::rpkm works on sparse matrxi from single cell
  rpkm <- sweep(assays(rawData)$counts * 1e9, MARGIN=1,
                STATS=rowData(rawData)$featWidth, FUN="/")
  rpkm <- sweep(rpkm, MARGIN=2, STATS=libSize, FUN="/")
  return(rpkm)
}

getTpm = function(rawData) {
  
  if (!is.null(assays(rawData)$tpm)){
    return(assays(rawData)$tpm)
  }
  tpm <- sweep(assays(rawData)$counts * 1e3, MARGIN=1,
               STATS=rowData(rawData)$featWidth, FUN="/")
  tpm <- sweep(tpm * 1e6, MARGIN=2, STATS=Matrix::colSums(tpm),
               FUN="/")
  return(tpm)
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
