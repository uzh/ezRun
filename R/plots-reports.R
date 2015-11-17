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
##' @param colors a character vector containing rgb codes. The default is a scale from blue to red.
##' @param types a character vector containing the types.
##' @template roxygen-template
addQcScatterPlots = function(doc, param, design, conds, rawData, signalCond, isPresentCond, types=NULL){
  samples = rownames(design)
  nConds = length(unique(conds))
  signal = getSignal(rawData)
  signal[signal <= 0] = NA
  isPresent = ezPresentFlags(signal, presentFlag=rawData$presentFlag, param=param, isLog=rawData$isLog)
  signalRange = range(signal, na.rm=TRUE)
  qcScatterTitles = list()
  qcScatterTitles[["Scatter Plots by Conditions"]] = "Scatter Plots by Conditions"
  addTitleWithAnchor(doc, qcScatterTitles[[length(qcScatterTitles)]], 2)
  if (!is.null(rawData$seqAnno$gc)){
    gcTypes = data.frame("GC < 0.4"=as.numeric(rawData$seqAnno$gc) < 0.4,
                         "GC > 0.6"=as.numeric(rawData$seqAnno$gc) > 0.6,
                         check.names=FALSE)
  } else {
    gcTypes = NULL
  }
  if (!is.null(rawData$seqAnno$width)){
    widthTypes = data.frame("width < 500nt"=as.numeric(rawData$seqAnno$width) < 500, 
                            "width > 5000nt"=as.numeric(rawData$seqAnno$width) > 5000,
                            check.names=FALSE)
  } else {
    widthTypes = NULL
  }
  if (nConds > 1 & nConds <=  param$allPairsMaxCondNumber){
    plotCmd = expression({
      ezAllPairScatter(signalCond, isPresent=isPresentCond, types=types)
    })
    defLink = ezImageFileLink(plotCmd, file="allPairs-scatter.png",
                              width=min(max(ncol(signalCond) * 200, 480), 2000),
                              height=min(max(ncol(signalCond) * 200, 480), 2000))
    if (!is.null(gcTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by GC", isPresent=isPresentCond, types=gcTypes)
      })
      gcLink = ezImageFileLink(plotCmd, file="allPairs-scatter-byGc.png",
                               width=min(max(ncol(signalCond) * 200, 480), 2000),
                               height=min(max(ncol(signalCond) * 200, 480), 2000))
    }
    if (!is.null(widthTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by width", isPresent=isPresentCond, types=widthTypes)
      })
      widthLink = ezImageFileLink(plotCmd, file="allPairs-scatter-byWidth.png",
                                  width=min(max(ncol(signalCond) * 200, 480), 2000),
                                  height=min(max(ncol(signalCond) * 200, 480), 2000))
    }
    doc = addFlexTable(doc, ezGrid(cbind(defLink, ifelse(!is.null(gcTypes), gcLink, NULL) , ifelse(!is.null(widthTypes), widthLink, NULL))))
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
        addTitleWithAnchor(doc, qcScatterTitles[[length(qcScatterTitles)]], 3)
        pngName = ezValidFilename(paste0(condName, "-scatter.png"))
        plotCmd = expression({
          ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=types, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
        })
        defLink = ezImageFileLink(plotCmd, file=pngName,
                                  width=min(nPlots, 6) * 480,
                                  height=ceiling(nPlots/6) * 480)
        if (!is.null(gcTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByGcScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=gcTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          })
          gcLink = ezImageFileLink(plotCmd, file=pngName,
                                   width=min(nPlots, 6) * 480,
                                   height=ceiling(nPlots/6) * 480)
        }
        if (!is.null(widthTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByWidthScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=widthTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          })
          widthLink = ezImageFileLink(plotCmd, file=pngName,
                                      width=min(nPlots, 6) * 480,
                                      height=ceiling(nPlots/6) * 480)
        }
        doc = addFlexTable(doc, ezGrid(cbind(defLink, ifelse(!is.null(gcTypes), gcLink, NULL) , ifelse(!is.null(widthTypes), widthLink, NULL))))
      }
    }
  }
  return(qcScatterTitles)
}

##' @title Adds test scatter plots
##' @description Adds test scatter plots to an html file.
##' @template doc-template
##' @templateVar object plots
##' @param param a list of parameters.
##' @param x a vector containing the signal.
##' @template result-template
##' @param seqAnno the sequence annotation. This is used if types is NULL.
##' @param types a character vector containing the types.
##' @template roxygen-template
addTestScatterPlots = function(doc, param, x, result, seqAnno, types=NULL){
  if (is.null(types)){
    types = data.frame(row.names=rownames(x))
    if ("IsControl" %in% colnames(seqAnno)){
      types$Controls = seqAnno[ rownames(x), "IsControl"]
    }
  }
  msg = "Highlighting significants with: "
  if (!is.null(param$pValueHighlightThresh)){
    significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
    types$Significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
    msg = paste(msg, "p <= ", param$pValueHighlightThresh)
    if (!is.null(param$log2RatioHighlightThresh)){
      msg = paste(msg, "and log ratio >= ", param$log2RatioHighlightThresh)
      if (!is.null(result$log2Ratio)){
        types$Significants = types$Significants & abs(result$log2Ratio) >= param$log2RatioHighlightThresh
      } else {
        types$Significants = types$Significants & result$log2Effect >= param$log2RatioHighlightThresh
      }
    }
  }
  
  testScatterTitles = list()
  testScatterTitles[["Scatter Plots"]] = "Scatter Plots"
  addTitleWithAnchor(doc, testScatterTitles[[length(testScatterTitles)]], 2)
  doc = addParagraph(doc, msg)
  testScatterTitles[["Between-group Comparison"]] = "Between-group Comparison"
  addTitleWithAnchor(doc, testScatterTitles[[length(testScatterTitles)]], 3)
  
  links = character()
  if (ncol(result$groupMeans) == 2 & !is.null(param$sampleGroup) & !is.null(param$refGroup)){
    sampleValues = 2^result$groupMeans[ , param$sampleGroup]
    refValues = 2^result$groupMeans[ , param$refGroup]
    plotCmd = expression({
      ezScatter(x=refValues, y=sampleValues, isPresent=result$usedInTest, types=types, xlab=param$refGroup, ylab=param$sampleGroup)
    })
    links["scatter"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                       width=min(ncol(as.matrix(sampleValues)), 6) * 480,
                                       height=ceiling(ncol(as.matrix(sampleValues))/6) * 480)
    plotCmd = expression({
      ezVolcano(log2Ratio=result$log2Ratio, pValue=result$pValue, isPresent=result$usedInTest, types=types, main=param$comparison)
    })
    links["volcano"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-volcano.png"))
    plotCmd = expression({
      ezVolcano(log2Ratio=result$log2Ratio, pValue=result$fdr, isPresent=result$usedInTest, types=types, main=param$comparison, yType="FDR")
    })
    links["volcanoFdr"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-FDR-volcano.png"))
  } else {
    plotCmd = expression({
      ezAllPairScatter(x=2^result$groupMeans, isPresent=result$usedInTest, types=types)
    })
    links["allPair"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                       width=min(max(ncol(result$groupMeans) * 200, 480), 2000),
                                       height=min(max(ncol(result$groupMeans) * 200, 480), 2000))
  }
  
  plotCmd = expression({
    myBreaks = seq(0, 1, by=0.002)
    histUsed = hist(result$pValue[result$usedInTest], breaks=myBreaks, plot=FALSE)
    histAbs = hist(result$pValue[!result$usedInTest], breaks=myBreaks, plot=FALSE)
    xx = rbind(used=histUsed$counts, absent=histAbs$counts)
    xx = shrinkToRange(xx, c(0, max(xx["used", ])))
    barplot(xx, space=0, border=NA, col=c("blue", "darkorange"), 
            xlab="p-value", ylab="counts", ylim=c(0, max(xx["used", ])),
            main="p-value histogram")
    abline(h=sum(result$usedInTest)/ncol(xx))
    at = c(0.01, 0.1, 0.25, 0.5, 0.75, 1)
    axis(1, at=at*ncol(xx), labels = at)
    legend("top", c("used", "absent"), col=c("blue", "darkorange"), pch=20, cex=1)
  })
  links["pValueHist"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-pValueHist.png"), height=400, width=800)
  doc = addFlexTable(doc, ezGrid(rbind(links)))
  
  if (is.null(x)){
    return()
  }
  
  theRange = 2^(range(x, na.rm=TRUE))
  advancedTitles = list()
  advancedTitles[["Advanced Plots"]] = "Advanced Plots"
  advancedDoc = openBsdocReport(title=advancedTitles[[length(advancedTitles)]])
  if (!ezIsSpecified(param$batch)){ ## TODO: we no longer use pairing, we now use batch which is more general; however these plots only work if batch is a real pairing
    for (group in unique(c(param$refGroup, colnames(result$groupMeans)))){
      idx = which(group == param$grouping)
      if (length(idx) > 1){
        advancedTitles[[paste("Intra-group Comparison:", group)]] = paste("Intra-group Comparison:", group)
        addTitleWithAnchor(advancedDoc, advancedTitles[[length(advancedTitles)]], 3)
        pngName = paste0(group, "-scatter.png")
        xlab = paste("Avg of", group)
        refValues = result$groupMeans[ , group]
        plotCmd = expression({
          ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], isPresent=result$isPresent[, idx, drop=FALSE], types=types, lim=theRange, xlab=xlab)
        })
        advancedDoc = addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                                    width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                                    height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480))
        if (ncol(result$groupMeans) == 2){
          otherGroup = setdiff(colnames(result$groupMeans), group)
          pngName = paste0(group, "-over-", otherGroup, "-scatter.png")
          xlab = paste("Avg of", otherGroup)
          refValues = result$groupMeans[ , otherGroup]
          plotCmd = expression({
            ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], isPresent=result$isPresent[, idx, drop=FALSE], types=types, lim=theRange, xlab=xlab)
          })
          advancedDoc = addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                                      width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                                      height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480))
        }
      }
    }
  } else {
    advancedTitles[["Pairs ... over ..."]] = paste("Pairs:", param$sampleGroup, "over", param$refGroup)
    addTitleWithAnchor(advancedDoc, advancedTitles[[length(advancedTitles)]], 3)
    use = param$grouping %in% c(param$sampleGroup, param$refGroup)
    if (all(table(param$batch[use], param$grouping[use]) == 1)){
      groups = paste(param$grouping, param$batch, sep="--")
      sampleGroups = sort(unique(groups[param$grouping == param$sampleGroup]))
      refGroups = sort(unique(groups[param$grouping == param$refGroup]))
      avgValues = averageColumns(x[ ,use], groups[use], mean)
      avgPresent= averageColumns(x[ ,use], groups[use], function(x){mean(x) > 0.5})
      sampleValues = avgValues[ , sampleGroups, drop=FALSE]
      refValues = avgValues[ , refGroups, drop=FALSE]
      samplePresent = avgPresent[ ,sampleGroups, drop=FALSE]
      refPresent = avgPresent[ , refGroups, drop=FALSE]
      pngName = paste0(param$sampleGroup, "-over-", param$refGroup, "-pairs.png")
      plotCmd = expression({
        ezScatter(x=2^refValues, y=2^sampleValues, isPresent=samplePresent | refPresent, types=types, lim=theRange, xlab=colnames(refValues))
      })
      advancedDoc = addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                                  width=min(ncol(as.matrix(sampleValues)), 6) * 480,
                                                                  height=ceiling(ncol(as.matrix(sampleValues))/6) * 480))
    }
  }
  closeBsdocReport(advancedDoc, "advancedPlots.html", advancedTitles)
  doc = addParagraph(doc, pot("advancedPlots.html", hyperlink = "advancedPlots.html"))
  return(testScatterTitles)
}
