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
                         isPresent=isPresentCond, types=widthTypes)
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
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=widthTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
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


##' @title Adds test scatter plots
##' @description Adds test scatter plots to an html file.
##' @template doc-template
##' @templateVar object plots
##' @param param a list of parameters.
##' @param x a vector containing the signal.
##' @template result-template
##' @param seqAnno the sequence annotation. This is used if types is NULL.
##' @param resultFile a character representing the link to the result file.
##' @template types-template
##' @template roxygen-template
addTestScatterPlots = function(doc, param, x, result, seqAnno, resultFile, types=NULL){
  if (is.null(types)){
    types = data.frame(row.names=rownames(x))
    if ("IsControl" %in% colnames(seqAnno)){
      types$Controls = seqAnno[ rownames(x), "IsControl"]
    }
  }
  if (!is.null(param$pValueHighlightThresh)){
    significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
    types$Significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
    if (!is.null(param$log2RatioHighlightThresh)){
      if (!is.null(result$log2Ratio)){
        types$Significants = types$Significants & abs(result$log2Ratio) >= param$log2RatioHighlightThresh
      } else {
        types$Significants = types$Significants & result$log2Effect >= param$log2RatioHighlightThresh
      }
    }
  }
  
  testScatterTitles = list()
  testScatterTitles[["Inspection of significant genes"]] = "Inspection of significant genes"
  addTitle(doc, testScatterTitles[[length(testScatterTitles)]], 2, id=testScatterTitles[[length(testScatterTitles)]])
  paragraphs = character()
  if (!is.null(param$pValueHighlightThresh)){
    paragraphs["P-value threshold:"] = paste("p <=", param$pValueHighlightThresh)
  }
  if (!is.null(param$log2RatioHighlightThresh)){
    paragraphs["Log ratio threshold:"] = paste("log ratio >=", param$log2RatioHighlightThresh)
  }
  paragraphs["Number of significant genes:"] = length(which(types$Significants))
  paragraphs["Subsequent plots highlight significant genes in "] = as.html(pot("red", format = textProperties(color="red")))
  addFlexTable(doc, ezGrid(paragraphs, add.rownames=TRUE))
  testScatterTitles[["Between-group Comparison"]] = "Between-group Comparison"
  addTitle(doc, testScatterTitles[[length(testScatterTitles)]], 3, id=testScatterTitles[[length(testScatterTitles)]])
  
  links = ezMatrix("", rows=1:2, cols=1:2)
  advancedLinks = character()
  if (ncol(result$groupMeans) == 2 & !is.null(param$sampleGroup) & !is.null(param$refGroup)){
    sampleValues = 2^result$groupMeans[ , param$sampleGroup]
    refValues = 2^result$groupMeans[ , param$refGroup]
    
    sortedGenes = names(sort(result$pValue[types$Significants]))
    if (length(sortedGenes) > param$maxInteractivePoints){
      sortedGenes = sortedGenes[1:param$maxInteractivePoints]
    }
    useForInteractivePoints = names(result$pValue) %in% sortedGenes
    
    if (ezIsSpecified(seqAnno)){
      if (is.null(seqAnno$gene_name)) {
        if (!is.null(seqAnno$gene_id)){
          seqAnno$gene_name = seqAnno$gene_id
        } else {
          seqAnno$gene_name = rownames(seqAnno)
        }
      } 
      clickActions = paste0("window.open('http://www.ihop-net.org/UniPub/iHOP/?search=%22", seqAnno$gene_name[useForInteractivePoints],
                            "%22&field=synonym&ncbi_tax_id=0');")
      popupLabels = seqAnno$gene_name[useForInteractivePoints]
    } else {
      clickActions = paste0("window.open('http://www.ihop-net.org/UniPub/iHOP/?search=%22", rownames(result$groupMeans)[useForInteractivePoints],
                            "%22&field=synonym&ncbi_tax_id=0');")
      popupLabels = rownames(result$groupMeans)[useForInteractivePoints]
    }
    
    interactiveDoc = openBsdocReport("Interactive plots")
#     if (!is.null(types)){
#       plotCmd = expression({
#         ezLegend(legend=c("too few reads", "enough reads but not significant", colnames(types)), fill=c("grey", "black", rainbow(ncol(types))))
#       })
#       legendLink = ezImageFileLink(plotCmd, file="scatterAndVolcanoLegend.png", height=length(ncol(types))*15 + 55, width=300, addPdfLink=FALSE)
#       addParagraph(interactiveDoc, legendLink)
#     }
    .interactiveSmoothScatter = function(){
      ezSmoothScatter(x=refValues, y=sampleValues, isPresent=result$usedInTest, types=types,
                      xlab=param$refGroup, ylab=param$sampleGroup, legendPos=NULL, nbin=32)
      ## TODO replace this with working code using the new package
#       add.plot.interactivity(fun=points, col="red", pch=16,
#                              x=log2(refValues[useForInteractivePoints]),
#                              y=log2(sampleValues[useForInteractivePoints]),
#                              popup.labels=popupLabels,
#                              click.actions=clickActions)
    }
    addPlot(interactiveDoc, .interactiveSmoothScatter, par.properties=parLeft())
#     .interactiveVolcanoPlot = function(){
#       ezSmoothScatter(x=result$log2Ratio, y=result$pValue, isPresent=result$usedInTest, types=types, main=param$comparison, legendPos=NULL)
#       add.plot.interactivity(fun=points, col="red", pch=16,
#                              x=result$log2Ratio[useForInteractivePoints],
#                              y=result$pValue[useForInteractivePoints],
#                              popup.labels=popupLabels,
#                              click.actions=clickActions)
#       values = ezVolcano(log2Ratio=result$log2Ratio, pValue=result$pValue, isPresent=result$usedInTest, types=types, main=param$comparison, legendPos=NULL)
#       add.plot.interactivity(fun=points, col="red", pch=16,
#                              x=values$x[useForInteractivePoints],
#                              y=values$y[useForInteractivePoints],
#                              popup.labels=popupLabels,
#                              click.actions=clickActions)
#     }
#     addPlot(interactiveDoc, .interactiveVolcanoPlot, fontname="", par.properties=parLeft()) ## the plots are plotted next to each other. I don't know why, but I like it.
    addParagraph(interactiveDoc, "Significant genes are plotted in red and clicking on them will open an iHop search of the gene.")
    closeBsdocReport(interactiveDoc, "interactivePlots.html")
    links[2, 1] = as.html(pot("Comparison plot of clickable significant genes", hyperlink = "interactivePlots.html"))
    
    plotCmd = expression({
      ezXYScatter(xVec=refValues, yVec=sampleValues, isPresent=result$usedInTest, types=types,
                         xlab=param$refGroup, ylab=param$sampleGroup, main="Comparison of average expression")
    })
    links[1, 1] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                       width=min(ncol(as.matrix(sampleValues)), 6) * 400,
                                       height=ceiling(ncol(as.matrix(sampleValues))/6) * 400) # dynamic png with possibly many plots
    plotCmd = expression({
      ezVolcano(log2Ratio=result$log2Ratio, pValue=result$pValue, isPresent=result$usedInTest, types=types, main=param$comparison)
    })
    links[1, 2] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-volcano.png"))
    plotCmd = expression({
      ezVolcano(log2Ratio=result$log2Ratio, pValue=result$fdr, isPresent=result$usedInTest, types=types, main=param$comparison, yType="FDR")
    })
    advancedLinks["volcanoFdr"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-FDR-volcano.png"))
  } else {
    plotCmd = expression({
      ezAllPairScatter(x=2^result$groupMeans, isPresent=result$usedInTest, types=types)
    })
    links[1, 1] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                       width=min(max(ncol(result$groupMeans) * 200, 400), 2000),
                                       height=min(max(ncol(result$groupMeans) * 200, 400), 2000)) # dynamic png with possibly many plots
  }
  
  tableLink = sub(".txt", "-viewTopSignificantGenes.html", resultFile)
  links[2, 2] = as.html(ezLink(tableLink, label="Interactive table of significant genes", target = "_blank"))
  
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
    legend("top", c("used", "not expressed"), col=c("blue", "darkorange"), pch=20, cex=1)
  })
  advancedLinks["pValueHist"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-pValueHist.png"), width=800)
  
  if (is.null(x)){
    return()
  }
  
  theRange = 2^(range(x, na.rm=TRUE))
  advancedTitles = list()
  advancedTitles[["Advanced Plots"]] = "Advanced Plots"
  advancedDoc = openBsdocReport(title=advancedTitles[[length(advancedTitles)]])
  addFlexTable(advancedDoc, ezGrid(rbind(advancedLinks)))
  if (!ezIsSpecified(param$grouping2)){ ## TODO: we no longer use pairing, we now use batch which is more general; however these plots only work if batch is a real pairing
    for (group in unique(c(param$refGroup, colnames(result$groupMeans)))){
      idx = which(group == param$grouping)
      if (length(idx) > 1){
        advancedTitles[[paste("Intra-group Comparison:", group)]] = paste("Intra-group Comparison:", group)
        addTitle(advancedDoc, advancedTitles[[length(advancedTitles)]], 3, id=advancedTitles[[length(advancedTitles)]])
        pngName = paste0(group, "-scatter.png")
        xlab = paste("Avg of", group)
        refValues = result$groupMeans[ , group]
        plotCmd = expression({
          ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], isPresent=result$isPresent[, idx, drop=FALSE], types=types, lim=theRange, xlab=xlab)
        })
        addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                                    width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                                    height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480)) # dynamic png with possibly many plots
        if (ncol(result$groupMeans) == 2){
          otherGroup = setdiff(colnames(result$groupMeans), group)
          pngName = paste0(group, "-over-", otherGroup, "-scatter.png")
          xlab = paste("Avg of", otherGroup)
          refValues = result$groupMeans[ , otherGroup]
          plotCmd = expression({
            ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], isPresent=result$isPresent[, idx, drop=FALSE], types=types, lim=theRange, xlab=xlab)
          })
          addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                                      width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                                      height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480)) # dynamic png with possibly many plots
        }
      }
    }
  } else {
    advancedTitles[["Pairs ... over ..."]] = paste("Pairs:", param$sampleGroup, "over", param$refGroup)
    addTitle(advancedDoc, advancedTitles[[length(advancedTitles)]], 3, id=advancedTitles[[length(advancedTitles)]])
    use = param$grouping %in% c(param$sampleGroup, param$refGroup)
    if (all(table(param$grouping2[use], param$grouping[use]) == 1)){
      groups = paste(param$grouping, param$grouping2, sep="--")
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
      addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                                  width=min(ncol(as.matrix(sampleValues)), 6) * 400,
                                                                  height=ceiling(ncol(as.matrix(sampleValues))/6) * 400)) # dynamic png with possibly many plots
    }
  }
  closeBsdocReport(advancedDoc, "advancedPlots.html", advancedTitles)
#   testScatterTitles[["Significant Genes"]] = "Significant Genes"
#   addTitle(doc, testScatterTitles[[length(testScatterTitles)]], 3, id=testScatterTitles[[length(testScatterTitles)]])
  addFlexTable(doc, ezGrid(links))
  return(testScatterTitles)
}

addTestScatterPlotsSE = function(doc, param, x, se, resultFile, types=NULL){
  seqAnno <- data.frame(rowData(se), row.names=rownames(se),
                        check.names = FALSE, stringsAsFactors=FALSE)
  if (is.null(types)){
    types = data.frame(row.names=rownames(x))
    if ("IsControl" %in% colnames(seqAnno)){
      types$Controls = seqAnno[rownames(x), "IsControl"]
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
  
  testScatterTitles = list()
  testScatterTitles[["Inspection of significant genes"]] = "Inspection of significant genes"
  addTitle(doc, testScatterTitles[[length(testScatterTitles)]], 2, 
           id=testScatterTitles[[length(testScatterTitles)]])
  paragraphs = character()
  if (!is.null(param$pValueHighlightThresh)){
    paragraphs["P-value threshold:"] = paste("p <=", param$pValueHighlightThresh)
  }
  if (!is.null(param$log2RatioHighlightThresh)){
    paragraphs["Log ratio threshold:"] = paste("log ratio >=", param$log2RatioHighlightThresh)
  }
  paragraphs["Number of significant genes:"] = length(which(types$Significants))
  paragraphs["Subsequent plots highlight significant genes in "] = as.html(pot("red", format = textProperties(color="red")))
  addFlexTable(doc, ezGrid(paragraphs, add.rownames=TRUE))
  testScatterTitles[["Between-group Comparison"]] = "Between-group Comparison"
  addTitle(doc, testScatterTitles[[length(testScatterTitles)]], 3, id=testScatterTitles[[length(testScatterTitles)]])
  
  links = ezMatrix("", rows=1:2, cols=1:2)
  advancedLinks = character()
  
  ## We compute groupMeans for plotting
  logSignal = log2(shiftZeros(assays(se)$xNorm, param$minSignal))
  groupMeans = cbind(rowMeans(logSignal[ , param$grouping == param$sampleGroup, drop=FALSE]),
                     rowMeans(logSignal[ , param$grouping == param$refGroup, drop=FALSE]))
  colnames(groupMeans) = c(param$sampleGroup, param$refGroup)
  
  if (ncol(groupMeans) == 2 & !is.null(param$sampleGroup) & !is.null(param$refGroup)){
    sampleValues = 2^groupMeans[ , param$sampleGroup]
    refValues = 2^groupMeans[ , param$refGroup]
    
    seSignificant <- se[types$Significants, ]
    sortedGenes <- rownames(seSignificant[order(rowData(seSignificant)$pValue), ])
    
    if (length(sortedGenes) > param$maxInteractivePoints){
      sortedGenes = sortedGenes[1:param$maxInteractivePoints]
    }
    useForInteractivePoints = rownames(se) %in% sortedGenes
    
    if (ezIsSpecified(seqAnno)){
      if (is.null(seqAnno$gene_name)) {
        if (!is.null(seqAnno$gene_id)){
          seqAnno$gene_name = seqAnno$gene_id
        } else {
          seqAnno$gene_name = rownames(seqAnno)
        }
      } 
      clickActions = paste0("window.open('http://www.ihop-net.org/UniPub/iHOP/?search=%22", seqAnno$gene_name[useForInteractivePoints],
                            "%22&field=synonym&ncbi_tax_id=0');")
      popupLabels = seqAnno$gene_name[useForInteractivePoints]
    } else {
      clickActions = paste0("window.open('http://www.ihop-net.org/UniPub/iHOP/?search=%22", rownames(groupMeans)[useForInteractivePoints],
                            "%22&field=synonym&ncbi_tax_id=0');")
      popupLabels = rownames(groupMeans)[useForInteractivePoints]
    }
    
    interactiveDoc = openBsdocReport("Interactive plots")
    #     if (!is.null(types)){
    #       plotCmd = expression({
    #         ezLegend(legend=c("too few reads", "enough reads but not significant", colnames(types)), fill=c("grey", "black", rainbow(ncol(types))))
    #       })
    #       legendLink = ezImageFileLink(plotCmd, file="scatterAndVolcanoLegend.png", height=length(ncol(types))*15 + 55, width=300, addPdfLink=FALSE)
    #       addParagraph(interactiveDoc, legendLink)
    #     }
    .interactiveSmoothScatter = function(){
      ezSmoothScatter(x=refValues, y=sampleValues, isPresent=rowData(se)$usedInTest, types=types,
                      xlab=param$refGroup, ylab=param$sampleGroup, legendPos=NULL, nbin=32)
      ## TODO replace this with working code using the new package
      #       add.plot.interactivity(fun=points, col="red", pch=16,
      #                              x=log2(refValues[useForInteractivePoints]),
      #                              y=log2(sampleValues[useForInteractivePoints]),
      #                              popup.labels=popupLabels,
      #                              click.actions=clickActions)
    }
    addPlot(interactiveDoc, .interactiveSmoothScatter, par.properties=parLeft())
    #     .interactiveVolcanoPlot = function(){
    #       ezSmoothScatter(x=result$log2Ratio, y=result$pValue, isPresent=result$usedInTest, types=types, main=param$comparison, legendPos=NULL)
    #       add.plot.interactivity(fun=points, col="red", pch=16,
    #                              x=result$log2Ratio[useForInteractivePoints],
    #                              y=result$pValue[useForInteractivePoints],
    #                              popup.labels=popupLabels,
    #                              click.actions=clickActions)
    #       values = ezVolcano(log2Ratio=result$log2Ratio, pValue=result$pValue, isPresent=result$usedInTest, types=types, main=param$comparison, legendPos=NULL)
    #       add.plot.interactivity(fun=points, col="red", pch=16,
    #                              x=values$x[useForInteractivePoints],
    #                              y=values$y[useForInteractivePoints],
    #                              popup.labels=popupLabels,
    #                              click.actions=clickActions)
    #     }
    #     addPlot(interactiveDoc, .interactiveVolcanoPlot, fontname="", par.properties=parLeft()) ## the plots are plotted next to each other. I don't know why, but I like it.
    addParagraph(interactiveDoc, "Significant genes are plotted in red and clicking on them will open an iHop search of the gene.")
    closeBsdocReport(interactiveDoc, "interactivePlots.html")
    links[2, 1] = as.html(pot("Comparison plot of clickable significant genes", hyperlink = "interactivePlots.html"))
    
    plotCmd = expression({
      ezXYScatter(xVec=refValues, yVec=sampleValues, isPresent=rowData(se)$usedInTest, types=types,
                  xlab=param$refGroup, ylab=param$sampleGroup, main="Comparison of average expression")
    })
    links[1, 1] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                  width=min(ncol(as.matrix(sampleValues)), 6) * 400,
                                  height=ceiling(ncol(as.matrix(sampleValues))/6) * 400) # dynamic png with possibly many plots
    plotCmd = expression({
      ezVolcano(log2Ratio=rowData(se)$log2Ratio, pValue=rowData(se)$pValue, 
                isPresent=rowData(se)$usedInTest, types=types, 
                main=param$comparison)
    })
    links[1, 2] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-volcano.png"))
    plotCmd = expression({
      ezVolcano(log2Ratio=rowData(se)$log2Ratio, pValue=rowData(se)$fdr, 
                isPresent=rowData(se)$usedInTest, types=types, 
                main=param$comparison, yType="FDR")
    })
    advancedLinks["volcanoFdr"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-FDR-volcano.png"))
  } else {
    plotCmd = expression({
      ezAllPairScatter(x=2^groupMeans, isPresent=rowData(se)$usedInTest, 
                       types=types)
    })
    links[1, 1] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                  width=min(max(ncol(groupMeans) * 200, 400), 2000),
                                  height=min(max(ncol(groupMeans) * 200, 400), 2000)) # dynamic png with possibly many plots
  }
  
  tableLink = sub(".txt", "-viewTopSignificantGenes.html", resultFile)
  links[2, 2] = as.html(ezLink(tableLink, label="Interactive table of significant genes", target = "_blank"))
  
  plotCmd = expression({
    myBreaks = seq(0, 1, by=0.002)
    histUsed = hist(rowData(se)$pValue[rowData(se)$usedInTest], breaks=myBreaks, plot=FALSE)
    histAbs = hist(rowData(se)$pValue[!rowData(se)$usedInTest], breaks=myBreaks, plot=FALSE)
    xx = rbind(used=histUsed$counts, absent=histAbs$counts)
    xx = shrinkToRange(xx, c(0, max(xx["used", ])))
    barplot(xx, space=0, border=NA, col=c("blue", "darkorange"), 
            xlab="p-value", ylab="counts", ylim=c(0, max(xx["used", ])),
            main="p-value histogram")
    abline(h=sum(rowData(se)$usedInTest)/ncol(xx))
    at = c(0.01, 0.1, 0.25, 0.5, 0.75, 1)
    axis(1, at=at*ncol(xx), labels = at)
    legend("top", c("used", "not expressed"), col=c("blue", "darkorange"), pch=20, cex=1)
  })
  advancedLinks["pValueHist"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-pValueHist.png"), width=800)
  
  if (is.null(x)){
    return()
  }
  
  theRange = 2^(range(x, na.rm=TRUE))
  advancedTitles = list()
  advancedTitles[["Advanced Plots"]] = "Advanced Plots"
  advancedDoc = openBsdocReport(title=advancedTitles[[length(advancedTitles)]])
  addFlexTable(advancedDoc, ezGrid(rbind(advancedLinks)))
  if (!ezIsSpecified(param$grouping2)){ ## TODO: we no longer use pairing, we now use batch which is more general; however these plots only work if batch is a real pairing
    for (group in unique(c(param$refGroup, colnames(groupMeans)))){
      idx = which(group == param$grouping)
      if (length(idx) > 1){
        advancedTitles[[paste("Intra-group Comparison:", group)]] = paste("Intra-group Comparison:", group)
        addTitle(advancedDoc, advancedTitles[[length(advancedTitles)]], 3, id=advancedTitles[[length(advancedTitles)]])
        pngName = paste0(group, "-scatter.png")
        xlab = paste("Avg of", group)
        refValues = groupMeans[ , group]
        plotCmd = expression({
          ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], 
                    isPresent=assays(se)$isPresent[, idx, drop=FALSE], 
                    types=types, lim=theRange, xlab=xlab)
        })
        addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                  width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                  height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480)) # dynamic png with possibly many plots
        if (ncol(groupMeans) == 2){
          otherGroup = setdiff(colnames(groupMeans), group)
          pngName = paste0(group, "-over-", otherGroup, "-scatter.png")
          xlab = paste("Avg of", otherGroup)
          refValues = groupMeans[ , otherGroup]
          plotCmd = expression({
            ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], 
                      isPresent=assays(se)$isPresent[, idx, drop=FALSE], 
                      types=types, lim=theRange, xlab=xlab)
          })
          addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                    width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                    height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480)) # dynamic png with possibly many plots
        }
      }
    }
  } else {
    advancedTitles[["Pairs ... over ..."]] = paste("Pairs:", param$sampleGroup, "over", param$refGroup)
    addTitle(advancedDoc, advancedTitles[[length(advancedTitles)]], 3, id=advancedTitles[[length(advancedTitles)]])
    use = param$grouping %in% c(param$sampleGroup, param$refGroup)
    if (all(table(param$grouping2[use], param$grouping[use]) == 1)){
      groups = paste(param$grouping, param$grouping2, sep="--")
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
      addParagraph(advancedDoc, ezImageFileLink(plotCmd, file=pngName,
                                                width=min(ncol(as.matrix(sampleValues)), 6) * 400,
                                                height=ceiling(ncol(as.matrix(sampleValues))/6) * 400)) # dynamic png with possibly many plots
    }
  }
  closeBsdocReport(advancedDoc, "advancedPlots.html", advancedTitles)
  #   testScatterTitles[["Significant Genes"]] = "Significant Genes"
  #   addTitle(doc, testScatterTitles[[length(testScatterTitles)]], 3, id=testScatterTitles[[length(testScatterTitles)]])
  addFlexTable(doc, ezGrid(links))
  return(testScatterTitles)
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
