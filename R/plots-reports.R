###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Adds QC scatter plots
##' @description Adds QC scatter plots to an html file.
##' @template doc-template
##' @templateVar object plots
##' @param param a list of parameters.
##' @param design a data.frame containing the factorial design.
##' @param conds a named character vector containing the conditions of the factorial design.
##' @template rawData-template
##' @param signalCond a set of values containing signals averaged by the conditions.
##' @param isPresentCond Either NULL or a data.frame containing coloring information.
##' @template types-template
##' @template roxygen-template
countQcScatterPlots = function(param, design, conds, rawData, signalCond, 
                               isPresentCond, types=NULL){
  
  samples = rownames(design)
  nConds = length(unique(conds))
  signal = getSignal(rawData) + param$backgroundExpression
  #signal[signal <= 0] = NA
  isPresent = ezPresentFlags(signal, presentFlag=assays(rawData)$presentFlag, 
                             param=param, isLog=metadata(rawData)$isLog)
  signalRange = range(signal, na.rm=TRUE)
  if (!is.null(rowData(rawData)$gc)){
    gcTypes = data.frame("GC < 0.4"=as.numeric(rowData(rawData)$gc) < 0.4,
                         "GC > 0.6"=as.numeric(rowData(rawData)$gc) > 0.6,
                         check.names=FALSE)
  } else {
    gcTypes = NULL
  }
  if (!is.null(rowData(rawData)$featWidth)){
    widthTypes = data.frame("width < 500nt"=as.numeric(rowData(rawData)$featWidth) < 500, 
                            "width > 5000nt"=as.numeric(rowData(rawData)$featWidth) > 5000,
                            check.names=FALSE)
    widthColors = brewPalette(12, alpha=1)[c(7,6)]
  } else {
    widthTypes = NULL
  }
  
  widePlots <- list()
  # if (nConds > 1 & nConds <=  param$allPairsMaxCondNumber){
  #   widePlots[["allPairs"]] <- list()
  #   widePlots[["allPairs"]][["std"]] <- 
  #     ezAllPairScatter(signalCond, isPresent=isPresentCond, types=types)
  #   if (!is.null(gcTypes)){
  #     widePlots[["allPairs"]][["gc"]] <- 
  #       ezAllPairScatter(signalCond, main="color by GC", isPresent=isPresentCond, types=gcTypes)
  #   }
  #   if (!is.null(widthTypes)){
  #     widePlots[["allPairs"]][["width"]] <-
  #       ezAllPairScatter(signalCond, main="color by width", isPresent=isPresentCond, types=widthTypes, colors = widthColors)
  #   }
  # }
  narrowPlots <- list()
  nPlots <- 0
  for (factorName in head(colnames(design), 4)){ ## take the first 4 factors
    for (factorLevel in unique(design[, factorName])){
      idx = which(factorLevel == design[, factorName])
      if (length(idx) == 1){
        next
      }
      if (length(idx) == 2){
        nPlots = 1
      } else {
        nPlots = length(idx)
      }
      plotName = paste(factorName, factorLevel)
      narrowPlots[[plotName]] <- list()
      idx = idx[order(samples[idx])] ## order alphabetically
      narrowPlots[[plotName]][["std"]] <-
        ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=types, lim=signalRange, xlab=paste("Avg of", factorLevel), ylab=NULL, mode="ggplot2")
      if (!is.null(gcTypes)){
        narrowPlots[[plotName]][["gc"]] <-
          ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=gcTypes, lim=signalRange, xlab=paste("Avg of", factorLevel), ylab=NULL, mode="ggplot2")
      }
      if (!is.null(widthTypes)){
        narrowPlots[[plotName]][["width"]] <-
          ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=widthTypes, lim=signalRange, xlab=paste("Avg of", factorLevel), ylab=NULL, mode="ggplot2")
      }
    }
  }
  nPlots <- max(c(1, nPlots))
  return(list(widePlots=widePlots, narrowPlots=narrowPlots, nPlots=nPlots))
}


### -----------------------------------------------------------------
### Prepare the group means and types for producing the scatter plots
###
makeTestScatterData <- function(param, se, types=NULL){
  seqAnno <- data.frame(rowData(se), row.names=rownames(se),
                        check.names = FALSE, stringsAsFactors=FALSE)
  
  logSignal <- log2(shiftZeros(assays(se)$xNorm, param$minSignal))
  groupMeans <- cbind(rowMeans(logSignal[ , param$grouping == param$sampleGroup, 
                                          drop=FALSE]),
                      rowMeans(logSignal[ , param$grouping == param$refGroup, 
                                          drop=FALSE])
  )
  colnames(groupMeans) = c(param$sampleGroup, param$refGroup)
  
  if (is.null(types)){
    types = data.frame(row.names=rownames(se))
    if ("IsControl" %in% colnames(seqAnno)){
      types$Controls = seqAnno[rownames(se), "IsControl"]
    }
  }
  if (!is.null(param$pValueHighlightThresh)){
    significants = rowData(se)$pValue <= param$pValueHighlightThresh & rowData(se)$usedInTest
    types$Significants = significants
    if (!is.null(param$log2RatioHighlightThresh)){
      if (!is.null(rowData(se)$log2Ratio)){
        types$Significants = types$Significants & 
          abs(rowData(se)$log2Ratio) >= param$log2RatioHighlightThresh
      } else {
        types$Significants = types$Significants & rowData(se)$log2Effect >= param$log2RatioHighlightThresh
      }
    }
  }
  return(list(logSignal=logSignal, groupMeans=groupMeans, types=types))
}
