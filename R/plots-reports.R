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
addQcScatterPlots = function(doc, param, design, conds, rawData, signalCond, isPresentCond, types=NULL){
  samples = rownames(design)
  nConds = length(unique(conds))
  signal = getSignal(rawData)
  signal[signal <= 0] = NA
  isPresent = ezPresentFlags(signal, presentFlag=assays(rawData)$presentFlag, 
                             param=param, isLog=metadata(rawData)$isLog)
  signalRange = range(signal, na.rm=TRUE)
  qcScatterTitles = list()
  qcScatterTitles[["Scatter Plots by Conditions"]] = "Scatter Plots by Conditions"
  addTitle(doc, qcScatterTitles[[length(qcScatterTitles)]], 2, id=qcScatterTitles[[length(qcScatterTitles)]])
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
  } else {
    widthTypes = NULL
  }
  if (nConds > 1 & nConds <=  param$allPairsMaxCondNumber){
    plotCmd = expression({
      ezAllPairScatter(signalCond, isPresent=isPresentCond, types=types)
    })
    imgLinks = character()
    imgLinks["def"] = ezImageFileLink(plotCmd, file="allPairs-scatter.png",
                                      width=min(max(ncol(signalCond) * 200, 400), 2000),
                                      height=min(max(ncol(signalCond) * 200, 400), 2000)) # dynamic png with possibly many plots
    if (!is.null(gcTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by GC", isPresent=isPresentCond, types=gcTypes)
      })
      imgLinks["gc"] = ezImageFileLink(plotCmd, file="allPairs-scatter-byGc.png",
                                       width=min(max(ncol(signalCond) * 200, 400), 2000),
                                       height=min(max(ncol(signalCond) * 200, 400), 2000)) # dynamic png with possibly many plots
    }
    if (!is.null(widthTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by width", isPresent=isPresentCond, types=widthTypes)
      })
      imgLinks["width"] = ezImageFileLink(plotCmd, file="allPairs-scatter-byWidth.png",
                                          width=min(max(ncol(signalCond) * 200, 400), 2000),
                                          height=min(max(ncol(signalCond) * 200, 400), 2000)) # dynamic png with possibly many plots
    }
    addFlexTable(doc, ezGrid(t(imgLinks)))
  }
  for (i in 1:min(4, ncol(design))){
    for (cond in unique(design[,i])){
      idx = which(cond == design[,i])
      if (length(idx) > 1){
        idx = idx[order(samples[idx])] ## order alphabetically
        condName = paste(colnames(design)[i], cond)
        nPlots = length(idx)
        if (nPlots == 2) nPlots = 1
        qcScatterTitles[[paste(condName, cond, i)]] = condName
        addTitle(doc, qcScatterTitles[[length(qcScatterTitles)]], 3, id=qcScatterTitles[[length(qcScatterTitles)]])
        pngName = ezValidFilename(paste0(condName, "-scatter.png"))
        plotCmd = expression({
          ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=types, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
        })
        imgLinks = character()
        imgLinks["def"] = ezImageFileLink(plotCmd, file=pngName,
                                  width=min(nPlots, 6) * 300,
                                  height=ceiling(nPlots/6) * 330) # dynamic png with possibly many plots
        if (!is.null(gcTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByGcScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=gcTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          })
          imgLinks["gc"] = ezImageFileLink(plotCmd, file=pngName,
                                   width=min(nPlots, 6) * 300,
                                   height=ceiling(nPlots/6) * 330) # dynamic png with possibly many plots
        }
        if (!is.null(widthTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByWidthScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=widthTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          })
          imgLinks["width"] = ezImageFileLink(plotCmd, file=pngName,
                                      width=min(nPlots, 6) * 300,
                                      height=ceiling(nPlots/6) * 330) # dynamic png with possibly many plots
        }
        addFlexTable(doc, ezGrid(imgLinks))
      }
    }
  }
  return(qcScatterTitles)
}

countQcScatterPlots = function(param, design, conds, rawData, signalCond, 
                               isPresentCond, types=NULL){
  samples = rownames(design)
  nConds = length(unique(conds))
  signal = getSignal(rawData)
  signal[signal <= 0] = NA
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
  
  plotsToInclude <- list()
  
  if (nConds > 1 & nConds <=  param$allPairsMaxCondNumber){
    plotCmd = expression({
      ezAllPairScatter(signalCond, isPresent=isPresentCond, types=types)
    })
    imgLinks = character()
    imgLinks["def"] = ezImageFileLink(plotCmd, file="allPairs-scatter.png",
                                      width=min(max(ncol(signalCond) * 200, 400), 2000),
                                      height=min(max(ncol(signalCond) * 200, 400), 2000)) # dynamic png with possibly many plots
    plotsToInclude[["allPairs"]] <- c()
    plotsToInclude[["allPairs"]] <- c(plotsToInclude[["allPairs"]], "allPairs-scatter.png")
    
    if (!is.null(gcTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by GC", 
                         isPresent=isPresentCond, types=gcTypes)
      })
      imgLinks["gc"] = ezImageFileLink(plotCmd, file="allPairs-scatter-byGc.png",
                                       width=min(max(ncol(signalCond) * 200, 400), 2000),
                                       height=min(max(ncol(signalCond) * 200, 400), 2000)) # dynamic png with possibly many plots
      plotsToInclude[["allPairs"]] <- c(plotsToInclude[["allPairs"]], 
                                        "allPairs-scatter-byGc.png")
    }
    if (!is.null(widthTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by width", 
                         isPresent=isPresentCond, types=widthTypes, colors = widthColors)
      })
      imgLinks["width"] = ezImageFileLink(plotCmd, file="allPairs-scatter-byWidth.png",
                                          width=min(max(ncol(signalCond) * 200, 400), 2000),
                                          height=min(max(ncol(signalCond) * 200, 400), 2000)) # dynamic png with possibly many plots
      plotsToInclude[["allPairs"]] <- c(plotsToInclude[["allPairs"]], 
                                        "allPairs-scatter-byWidth.png")
    }
  }
  for (i in 1:min(4, ncol(design))){
    for (cond in unique(design[,i])){
      idx = which(cond == design[,i])
      if (length(idx) > 1){
        idx = idx[order(samples[idx])] ## order alphabetically
        condName = paste(colnames(design)[i], cond)
        plotsToInclude[[condName]] <- c()
        nPlots = length(idx)
        if (nPlots == 2) nPlots = 1
        pngName = ezValidFilename(paste0(condName, "-scatter.png"))
        plotCmd = expression({
          ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=types, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
        })
        imgLinks = character()
        imgLinks["def"] = ezImageFileLink(plotCmd, file=pngName,
                                          width=min(nPlots, 6) * 300,
                                          height=ceiling(nPlots/6) * 330) # dynamic png with possibly many plots
        plotsToInclude[[condName]] <- c(plotsToInclude[[condName]], pngName)
        
        if (!is.null(gcTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByGcScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=gcTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          })
          imgLinks["gc"] = ezImageFileLink(plotCmd, file=pngName,
                                           width=min(nPlots, 6) * 300,
                                           height=ceiling(nPlots/6) * 330) # dynamic png with possibly many plots
          plotsToInclude[[condName]] <- c(plotsToInclude[[condName]], pngName)
        }
        if (!is.null(widthTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByWidthScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=widthTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL, colors=widthColors)
          })
          imgLinks["width"] = ezImageFileLink(plotCmd, file=pngName,
                                              width=min(nPlots, 6) * 300,
                                              height=ceiling(nPlots/6) * 330) # dynamic png with possibly many plots
          plotsToInclude[[condName]] <- c(plotsToInclude[[condName]], pngName)
        }
      }
    }
  }
  return(plotsToInclude)
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
