###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Count QC
##' @template htmlFile-template
##' @seealso \code{\link{EzAppCountQC}}
ezMethodCountQC = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){

  if (is.null(param$runGO)){
    param$runGO = TRUE
  }
  
  cwd = getwd()
  on.exit(setwd(cwd))
  dataset = input$meta
  setwdNew(basename(output$getColumn("Report")))
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  }
  input$meta = dataset
      
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport(htmlFile, param=param, dataset=dataset, error=rawData$error)
    return("Error")
  }
  rawData$signal = ezNorm(rawData$counts, presentFlag=rawData$presentFlag, method=param$normMethod)
  runNgsCountQC(dataset, htmlFile, param, rawData=rawData)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodCountQC()
##' @seealso \code{\link{ezMethodCountQC}}
EzAppCountQC <-
  setRefClass("EzAppCountQC",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCountQC
                  name <<- "EzAppCountQC"
                }
              )
  )

##' @title Runs the NGS count QC
##' @description Runs the NGS count quality control, does various plots and creates a report.
##' @template dataset-template
##' @template htmlFile-template
##' @param param a list of parameters, possibly passed to other functions as well:
##' \itemize{
##'   \item{name}{ a character representing the name of the app.}
##'   \item{minSignal}{ a numeric or integer specifying the minimal signal amount.}
##'   \item{normMethod}{ the normalization method.}
##'   \item{sigThresh}{ the threshold...}
##'   \item{useSigThresh}{ ...and whether it should be used.}
##'   \item{doZip}{ a logical indicating whether to archive the result file.}
##'   \item{bgExpression}{ a numeric specifying the expression baseline value that is added before heatmap plots.}
##'   \item{topGeneSize}{ an integer specifying the number of high variance genes to consider in gene clustering.}
##'   \item{highVarThreshold}{ a numeric specifying the threshold for minimum standard deviation of log2 signal across samples.}
##'   \item{maxGenesForClustering}{ an integer specifying the maximum amount of genes for clustering. If the amount is higher, the least significant genes get removed first.}
##'   \item{minGenesForClustering}{ an integer specifying the minimum amount of genes for clustering. If the amount is lower, no clustering will be done.}
##'   \item{logColorRange}{ a logarithmic color range.}
##'   \item{writeScatterPlots}{ a logical indicating whether to write scatter plots.}
##' }
##' @template rawData-template
##' @param writeDataFiles a logical indicating whether to write the data files into seperate tables.
##' @param types a character vector containing the types.
##' @template roxygen-template
runNgsCountQC = function(dataset, htmlFile="00index.html", param=param, rawData=NULL,
                         writeDataFiles=TRUE, types=NULL){
  
  #library(affy, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  seqAnno = rawData$seqAnno
  if (is.null(types) && !is.null(seqAnno$type)){
    types = data.frame(row.names=rownames(seqAnno))
    for (nm in setdiff(na.omit(unique(seqAnno$type)), "")){
      types[[nm]] = seqAnno$type == nm
    }
  }
  
  design = ezDesignFromDataset(dataset, param)
  samples = rownames(design)
  nSamples = length(samples)
  conds = ezConditionsFromDesign(design, maxFactors = 2)
  nConds = length(unique(conds))
  sampleColors = getSampleColors(conds)
  
  titles = list()
  titles[["Analysis"]] = paste("Analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]], dataset=dataset)
  
  if (nSamples < 2){
    titles[["Note"]] = "Note: Statistics and Plots are not available for single sample experiments"
    addTitleWithAnchor(doc, titles[[length(titles)]], 2)
    doc = addParagraph(doc, "Run the report again with multiple samples selected.")
    return("Success")
  }
  
  signal = shiftZeros(getSignal(rawData), param$minSignal)
  presentFlag = rawData$presentFlag
  signalRange = range(signal, na.rm=TRUE)
  log2Signal = log2(signal)
  isPresent = ezPresentFlags(signal, presentFlag=presentFlag, param=param, isLog=rawData$isLog)
  signalCond = 2^averageColumns(log2Signal, by=conds)
  isPresentCond = averageColumns(isPresent, by=conds) >= 0.5
  isPresentStudy = apply(isPresentCond, 1, mean) >= 0.5
  
  titles[["Parameters"]] = "Parameters"
  addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  settings = character()
  settings["Normalization method:"] = param$normMethod
  if (param$useSigThresh){
    settings["Log2 signal threshold:"] = signif(log2(param$sigThresh), digits=4)
    settings["Linear signal threshold:"] = signif(param$sigThresh, digits=4)
  }
  settings["Feature level:"] = rawData$featureLevel
  settings["Number of features:"] = nrow(signal)
  settings["Data Column Used:"] = rawData$countName
  doc = addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
  
  if (writeDataFiles){
    if (!is.null(rawData$presentFlag)){
      combined = interleaveMatricesByColumn(rawData$signal, rawData$presentFlag)			
    } else {
      combined = rawData$signal
    }
    if (!is.null(seqAnno)){
      combined = cbind(seqAnno[rownames(combined), ,drop=FALSE], combined)
    }
    countFile = paste0(ezValidFilename(param$name), "-raw-count.txt")
    ezWrite.table(rawData$counts, file=countFile, head="Feature ID", digits=4)
    signalFile = paste0(ezValidFilename(param$name), "-normalized-signal.txt")
    ezWrite.table(combined, file=signalFile, head="Feature ID", digits=4)
    rpkmFile = paste0(ezValidFilename(param$name), "-rpkm.txt")
    ezWrite.table(getRpkm(rawData), file=rpkmFile, head="Feature ID", digits=4)
    if (param$doZip){
      dataFiles =c(zipFile(countFile), zipFile(signalFile), zipFile(rpkmFile))
      titles[["Data Files"]] = "Data Files"
      addTitleWithAnchor(doc, titles[[length(titles)]], 2)
      addTxtLinksToReport(doc, dataFiles, mime="application/zip")
    } else {
      titles[["Data Files"]] = "Data Files"
      addTitleWithAnchor(doc, titles[[length(titles)]], 2)
      addTxtLinksToReport(doc, c(countFile, signalFile, rpkmFile), mime="application/txt")
    }
  }
  
  titles[["Count Statistics"]] = "Count Statistics"
  addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  totalCounts = signif(apply(rawData$counts, 2, sum) / 1e6, digits=3)
  presentCounts = apply(isPresent, 2, sum)
  names(totalCounts) = samples
  names(presentCounts) = samples
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    barplot(totalCounts, las=2, ylab="Counts [Mio]", main="Total reads")
  })
  totalLink = ezImageFileLink(plotCmd, file="totalCounts.png", height=600,
                              width=min(600 + (length(samples)-10)* 30, 2000))
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    bplot = barplot(presentCounts, las=2, ylab="Counts", main="Genomic Features with Reads above threshold")
    percentages = paste(signif(100*presentCounts/nrow(isPresent), digits=2), "%")
    text(x=bplot, y=3, labels=percentages)
  })
  presentLink = ezImageFileLink(plotCmd, file="presentCounts.png", height=600,
                                width=min(600 + (length(samples)-10)* 30, 2000))
  
  doc = addFlexTable(doc, ezGrid(cbind(totalLink, presentLink)))
  
  rawData$signal = signal
  
  #################################
  ## correlation plot
  
  isValid = isPresentStudy
  if (!is.null(seqAnno$IsControl)){
    isValid = isValid & !seqAnno$IsControl
  }
  
  pngLinks = ezMatrix("", rows=1, cols=1:2)
  pngAdvancedLinks = pngLinks
  x = log2(2^log2Signal[isValid, ] + param$bgExpression)
  xNormed = sweep(x, 1 , rowMeans(x));
  xSd = apply(x, 1, sd, na.rm=TRUE)
  ord = order(xSd, decreasing=TRUE)
  topGenes = ord[1:min(param$topGeneSize, length(ord))]
  
  pngName = "sampleCorrelation-AllPresent.png"
  plotCmd = expression({
    ezCorrelationPlot(cor(x, use="complete.obs"), cond=conds, condLabels=conds, 
                      main=paste0("all present genes (", sum(isValid), ")"), sampleColors=sampleColors)
  })
  height = nrow(cor(x, use="complete.obs")) * 20
  if (height < 500) height= 500
  if (height > 2000) height = 2000
  try({
    pngLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=round(height*1.25), height=height)
  })
  
  pngName = "sampleCorrelation-AllPresentNormalized.png"
  plotCmd = expression({
    ezCorrelationPlot(cor(xNormed, use="complete.obs"), cond=conds, condLabels=conds, 
                      main=paste0("all present genes (", sum(isValid), ") gene-wise normalized"))
  })
  height = nrow(cor(xNormed, use="complete.obs")) * 20
  if (height < 500) height= 500
  if (height > 2000) height = 2000
  try({
    pngAdvancedLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=round(height*1.25), height=height)
  })
  
  pngName = "sampleCorrelation-TopGenes.png"
  plotCmd = expression({
    ezCorrelationPlot(cor(x[topGenes,], use="complete.obs"), cond=conds, condLabels=conds, 
                      main=paste0("top genes (", length(topGenes), ")"))
  })
  height = nrow(cor(x[topGenes,], use="complete.obs")) * 20
  if (height < 500) height= 500
  if (height > 2000) height = 2000
  try({
    pngLinks[1, 2] = ezImageFileLink(plotCmd, file=pngName, width=round(height*1.25), height=height)
  })
  
  pngName = "sampleCorrelation-TopGenesNormalized.png"
  plotCmd = expression({
    ezCorrelationPlot(cor(xNormed[topGenes,], use="complete.obs"), cond=conds, condLabels=conds, 
                      main=paste0("top genes (", length(topGenes), ") gene-wise normalized"))
  })
  height = nrow(cor(xNormed[topGenes,], use="complete.obs")) * 20
  if (height < 500) height= 500
  if (height > 2000) height = 2000
  try({
    pngAdvancedLinks[1, 2] = ezImageFileLink(plotCmd, file=pngName, width=round(height*1.25), height=height)
  })
  
  titles[["Sample correlation"]] = "Sample correlation"
  addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  doc = addFlexTable(doc, ezGrid(pngLinks))
  
  advancedTitles = list()
  advancedTitles[["Advanced Plots"]] = "Advanced Plots"
  advancedDoc = openBsdocReport(title=advancedTitles[[length(advancedTitles)]])
  advancedTitles[["Sample correlation"]] = "Sample correlation"
  addTitleWithAnchor(advancedDoc, advancedTitles[[length(advancedTitles)]], 2)
  advancedDoc = addFlexTable(advancedDoc, ezGrid(pngAdvancedLinks))
  
  
  ################################
  ## cluster plots
  
  if(nSamples > 3){
    pngLinks = ezMatrix("", rows=1, cols=1:2)
    pngAdvancedLinks = pngLinks
    pngName = "sampleClustering-AllPresent.png"
    plotCmd = expression({
      d = as.dist(1-cor(x, use="complete.obs"));
      hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
      hcd = colorClusterLabels(hcd, sampleColors)
      mai = par("mai")
      mai[1] = 3
      par(mai=mai)
      plot(hcd, main="all present genes", xlab="")
    })
    pngLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    
    pngName = "sampleClustering-AllPresentNormalized.png"
    plotCmd = expression({
      d = as.dist(1-cor(xNormed, use="complete.obs"));
      hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
      hcd = colorClusterLabels(hcd, sampleColors)
      mai = par("mai")
      mai[1] = 3
      par(mai=mai)
      plot(hcd, main="all present genes; gene-wise normalized", xlab="")
    })
    pngAdvancedLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    
    pngName = "sampleClustering-TopGenes.png"
    plotCmd = expression({
      d = as.dist(1-cor(x[topGenes, ], use="complete.obs"));
      hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
      hcd = colorClusterLabels(hcd, sampleColors)
      mai = par("mai")
      mai[1] = 3
      par(mai=mai)
      plot(hcd, main=paste("top", length(topGenes), " genes"), xlab="")
    })
    pngLinks[1, 2] = ezImageFileLink(plotCmd, file=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    
    pngName = "sampleClustering-TopGenesNormalized.png"
    plotCmd = expression({
      d = as.dist(1-cor(xNormed[topGenes, ], use="complete.obs"));
      hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
      hcd = colorClusterLabels(hcd, sampleColors)
      mai = par("mai")
      mai[1] = 3
      par(mai=mai)
      plot(hcd, main=paste("top", length(topGenes), "genes; gene-wise normalized"), xlab="")
    })
    pngAdvancedLinks[1, 2] = ezImageFileLink(plotCmd, file=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    
    titles[["Sample Clustering"]] = "Sample Clustering"
    addTitleWithAnchor(doc, titles[[length(titles)]], 2)
    doc = addFlexTable(doc, ezGrid(pngLinks))
    
    advancedTitles[["Sample Clustering"]] = "Sample Clustering"
    addTitleWithAnchor(advancedDoc, advancedTitles[[length(advancedTitles)]], 2)
    advancedDoc = addFlexTable(advancedDoc, ezGrid(pngAdvancedLinks))
    
    ## gene clustering
    use = xSd > param$highVarThreshold & apply(!is.na(x), 1, all)
    sdThresh = param$highVarThreshold
    if (sum(use, na.rm=TRUE) > param$maxGenesForClustering){
      use[use] = rank(-1 * xSd[use], ties.method="max") <= param$maxGenesForClustering
      sdThresh = signif(min(xSd[use]), digits=3)
    }
    
    if (sum(use, na.rm=TRUE) > param$minGenesForClustering){
      clusterPng = "cluster-heatmap.png"
      clusterColors = c("red", "yellow", "orange", "green", "blue", "cyan")
      
      clusterResult = clusterResults(xNormed[use, ], nClusters=6, clusterColors=clusterColors)
      plotCmd = expression({
        clusterHeatmap(xNormed[use, ], param, clusterResult, file=clusterPng, margins=c(18, 9),
                       colColors=sampleColors, lim=c(-param$logColorRange, param$logColorRange),
                       doClusterColumns=TRUE)
      })
      clusterLink = ezImageFileLink(plotCmd, file=clusterPng, width=max(800, 400 + 10 * ncol(xNormed[use, ])), height=1000)
      
      if (doGo(param, seqAnno)){
        clusterResult = goClusterResults(xNormed[use, ], param, clusterResult, seqAnno=seqAnno,
                                         universeProbeIds=rownames(seqAnno))
      }
      
      titles[["Clustering of High Variance Features"]] = "Clustering of High Variance Features"
      addTitleWithAnchor(doc, titles[[length(titles)]], 2)
      doc = addParagraph(doc, paste("Threshold for std. dev. of log2 signal across samples:", sdThresh))
      doc = addParagraph(doc, paste("Number of features with high std. dev.:", sum(use)))
      
      if (!is.null(clusterResult$GO)){
        goLink = as.html(ezGrid(c("Background color corresponds to the color of the feature cluster in the heatmap plot.",
                                  as.html(goClusterTable(param, clusterResult)))))
      } else {
        goLink = as.html(pot("No information available"))
      }
      tbl = ezGrid(t(c("Cluster Plot"=clusterLink, "GO categories of feature clusters"=goLink)), header.columns = TRUE)
      doc = addFlexTable(doc, tbl)
    }
    
    ##########################################
    ## mds plot
    
    titles[["MDS-Plot"]] = "MDS-Plot"
    addTitleWithAnchor(doc, titles[[length(titles)]], 2)
    pngLinks=vector(length=2)
    
    pngName = "mdsPlot_PresentGenes.png"
    plotCmd = expression({
      myMdsPlot(signal=x, sampleColors=sampleColors, main=sub('.png','',pngName))
    })
    presentLink = ezImageFileLink(plotCmd, file=pngName, width=480, height=480)
    
    pngName = "mdsPlot_TopGenes.png"
    plotCmd = expression({
      myMdsPlot(signal=x[topGenes,], sampleColors=sampleColors, main=sub('.png','',pngName))
    })
    topLink = ezImageFileLink(plotCmd, file=pngName, width=480, height=480)
    
    doc = addFlexTable(doc, ezGrid(cbind(presentLink, topLink)))
    
    if (param$writeScatterPlots){
      qcScatterTitles = addQcScatterPlots(doc, param, design, conds,
                                            rawData, signalCond, isPresentCond, types=types)
      titles = append(titles, qcScatterTitles)
    }
    
    if (!is.null(rawData$countsStart) & !is.null(rawData$countsEnd)){
      titles[["3' Bias analysis"]] = "3' Bias analysis"
      addTitleWithAnchor(doc, titles[[length(titles)]], 2)
      
      valStart = shiftZeros(rawData$countsStart, param$minSignal)
      colnames(valStart) = paste(colnames(valStart), "[transcript start]")
      valEnd = shiftZeros(rawData$countsEnd, param$minSignal)
      colnames(valEnd) = paste(colnames(valEnd), "[transcript end]")
      
      pngName = "start-end-countScatter.png"
      plotCmd = expression({
        ezScatter(x=valStart, y=valEnd, types=types)
      })
      pngLink = ezImageFileLink(plotCmd, file=pngName,
                                width=min(ncol(as.matrix(valEnd)), 6) * 480,
                                height=ceiling(ncol(as.matrix(valEnd))/6) * 480)
      doc = addParagraph(doc, pngLink)
      
      pngName = "start-end-countSmoothScatter.png"
      plotCmd = expression({
        ezSmoothScatter(x=valStart, y=valEnd)
      })
      pngLink = ezImageFileLink(plotCmd, file=pngName,
                                width=min(ncol(as.matrix(valEnd)), 6) * 480,
                                height=ceiling(ncol(as.matrix(valEnd))/6) * 480)
      doc = addParagraph(doc, pngLink)
    }
    
    ##########################################
    ## count density plots
    
    pngName = "signalDens.png"
    plotCmd = expression({
      countDensPlot(param, signal, sampleColors, main="all transcripts", bw=0.7)
    })
    pngLink = ezImageFileLink(plotCmd, file=pngName, width=640, height=640)
    
    titles[["Expression densities"]] = "Expression densities"
    addTitleWithAnchor(doc, titles[[length(titles)]], 2)
    doc = addParagraph(doc, "Zero or negative counts are not represented by the area!")
    doc = addParagraph(doc, pngLink)
  }
  closeBsdocReport(advancedDoc, "advancedPlots.html", advancedTitles)
  doc = addParagraph(doc, pot("advancedPlots.html", hyperlink = "advancedPlots.html"))
  closeBsdocReport(doc, htmlFile, titles)
}
