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
##'   \item{backgroundExpression}{ a numeric specifying the expression baseline value that is added before heatmap plots.}
##'   \item{topGeneSize}{ an integer specifying the number of high variance genes to consider in gene clustering.}
##'   \item{highVarThreshold}{ a numeric specifying the threshold for minimum standard deviation of log2 signal across samples.}
##'   \item{maxGenesForClustering}{ an integer specifying the maximum amount of genes for clustering. If the amount is higher, the least significant genes get removed first.}
##'   \item{minGenesForClustering}{ an integer specifying the minimum amount of genes for clustering. If the amount is lower, no clustering will be done.}
##'   \item{logColorRange}{ a logarithmic color range.}
##'   \item{writeScatterPlots}{ a logical indicating whether to write scatter plots.}
##' }
##' @template rawData-template
##' @param writeDataFiles a logical indicating whether to write the data files into separate tables.
##' @template types-template
##' @template roxygen-template
runNgsCountQC = function(htmlFile="00index.html",
                         rawData=SummarizedExperiment::SummarizedExperiment(),
                         writeDataFiles=TRUE, types=NULL, output=NULL){
  param <- metadata(rawData)$param
  if (is.null(assays(rawData)$signal)){
    assays(rawData)$signal = ezNorm(assays(rawData)$counts,
                                    presentFlag=assays(rawData)$presentFlag,
                                    method=param$normMethod)
  }
  seqAnno <- data.frame(rowData(rawData), row.names=rownames(rawData),
                        check.names = FALSE, stringsAsFactors=FALSE)
  dataset <- data.frame(colData(rawData), 
                        row.names=colnames(rawData), check.names = FALSE,
                        stringsAsFactors=FALSE)
  metadata(rawData)$analysis <- "Count_QC"
  
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
  titles[["Analysis"]] = paste("Analysis:", ifelse(grepl("bam", param$name, ignore.case = TRUE), sub("$", "_Count_QC", param$name), param$name))
  doc = openBsdocReport(title=titles[[length(titles)]])
  
  if (nSamples < 2){
    titles[["Note"]] = "Note: Statistics and Plots are not available for single sample experiments"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    addParagraph(doc, "Run the report again with multiple samples selected.")
    return("Success")
  }
  
  signal = shiftZeros(getSignal(rawData), param$minSignal)
  presentFlag = assays(rawData)$presentFlag
  signalRange = range(signal, na.rm=TRUE)
  log2Signal = log2(signal)
  isPresent = ezPresentFlags(assays(rawData)$counts, presentFlag=presentFlag, 
                             param=param, isLog=metadata(rawData)$isLog)
  signalCond = 2^averageColumns(log2Signal, by=conds)
  isPresentCond = averageColumns(isPresent, by=conds) >= 0.5
  isPresentStudy = apply(isPresentCond, 1, mean) >= 0.5
  
  addDataset(doc, dataset, param)
  
  settings = character()
  settings["Normalization method:"] = param$normMethod
  if (param$useSigThresh){
    settings["Log2 signal threshold:"] = signif(log2(param$sigThresh), digits=4)
    settings["Linear signal threshold:"] = signif(param$sigThresh, digits=4)
  }
  settings["Feature level:"] = metadata(rawData)$featureLevel
  settings["Number of features:"] = nrow(signal)
  settings["Data Column Used:"] = metadata(rawData)$countName
  addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
  
  if (!is.null(output)){
    liveReportLink = output$getColumn("Live Report")
    result = EzResult(se=rawData)
    result$saveToFile(basename(output$getColumn("Live Report")))
    addParagraph(doc, ezLink(liveReportLink,
                             "Live Report and Visualizations",
                             target = "_blank"))
  }  
  
  if (writeDataFiles){
    if (!is.null(assays(rawData)$presentFlag)){
      combined = interleaveMatricesByColumn(assays(rawData)$signal, 
                                            assays(rawData)$presentFlag)
    } else {
      combined = assays(rawData)$signal
    }
    if (!is.null(seqAnno)){
      combined = cbind(seqAnno[rownames(combined), ,drop=FALSE], combined)
    }
    if (!is.null(combined$featWidth)){
      combined$featWidth = as.integer(combined$featWidth)
    }
    countFile = paste0(ezValidFilename(param$name), "-raw-count.txt")
    ezWrite.table(assays(rawData)$counts, file=countFile, 
                  head="Feature ID", digits=NA)
    signalFile = paste0(ezValidFilename(param$name), "-normalized-signal.txt")
    ezWrite.table(combined, file=signalFile, head="Feature ID", digits=4)
    
    selectSignals = grepl("Signal", colnames(combined))
    combined$"Mean signal" = rowMeans(combined[, selectSignals])
    combined$"Maximum signal" = apply(combined[, selectSignals], 1, max)
    topGenesPerSample = apply(combined[, selectSignals], 2, function(col){
      col = sort(col, decreasing = TRUE)
      if (length(col) > 100) col = col[1:100]
      return(names(col))
    })
    
    topGenes = unique(as.character(topGenesPerSample))
    
    combined = combined[order(combined$"Maximum signal", decreasing = TRUE), , drop=FALSE]
    useInInteractiveTable = c("seqid", "gene_name", "Maximum signal", "Mean signal", "description", "featWidth", "gc")
    useInInteractiveTable = intersect(useInInteractiveTable, colnames(combined))
    tableLink = sub(".txt", "-viewHighExpressionGenes.html", signalFile)
    combinedTopGenes = combined[which(rownames(combined) %in% topGenes),] ## select top genes
    interactiveTable = head(combinedTopGenes[, useInInteractiveTable, drop=FALSE], param$maxTableRows) ## restrict number of table rows if necessary
    nRows = ifelse(length(topGenes)>=param$maxTableRows, param$maxTableRows, length(topGenes))
    ezInteractiveTable(interactiveTable, tableLink=tableLink, digits=3, colNames=c("ID", colnames(interactiveTable)),
                       title=paste("Showing the top", nRows, "genes with the highest expression"))
    
    rpkmFile = paste0(ezValidFilename(param$name), "-rpkm.txt")
    ezWrite.table(getRpkm(rawData), file=rpkmFile, head="Feature ID", digits=4) 
    
    tpmFile = paste0(ezValidFilename(param$name), "-tpm.txt")
    ezWrite.table(getTpm(rawData), file=tpmFile, head="Feature ID", digits=4)
    
    dataFiles = c(countFile, signalFile, rpkmFile, tpmFile)
    titles[["Data Files"]] = "Data Files"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    addTxtLinksToReport(doc, dataFiles, param$doZip)
    addParagraph(doc, ezLink(tableLink, target = "_blank"))
  }
  
  titles[["Count Statistics"]] = "Count Statistics"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  totalCounts = signif(colSums(assays(rawData)$counts) / 1e6, digits=3)
  presentCounts = colSums(isPresent)
  names(totalCounts) = samples
  names(presentCounts) = samples
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    barplot(totalCounts, las=2, ylab="Counts [Mio]", main="Total reads",
            names.arg=ezSplitLongLabels(names(totalCounts)))
  })
  totalLink = ezImageFileLink(plotCmd, file="totalCounts.png", width=min(600 + (nSamples-10)* 30, 2000)) # nSamples dependent width
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    bplot = barplot(presentCounts, las=2, ylab="Counts", main="Genomic Features with Reads above threshold",
                    names.arg=ezSplitLongLabels(names(presentCounts)))
    percentages = paste(signif(100*presentCounts/nrow(isPresent), digits=2), "%")
    text(x=bplot, y=0, pos=3, offset=2, labels=percentages)
  })
  presentLink = ezImageFileLink(plotCmd, file="presentCounts.png", width=min(600 + (nSamples-10)* 30, 2000)) # nSamples dependent width
  
  addFlexTable(doc, ezGrid(cbind(totalLink, presentLink)))
  
  assays(rawData)$signal = signal
  
  #################################
  ## correlation plot
  
  isValid = isPresentStudy
  if (!is.null(seqAnno$IsControl)){
    isValid = isValid & !seqAnno$IsControl
  }
  
  if (sum(isValid) < 10){
    doc = addParagraph(doc, "Not enough valid features for further plots")
    closeBsdocReport(doc, htmlFile, titles)
    return("SUCCESS")
  }
  
  pngLinks = ezMatrix("", rows=1, cols=1:2)
  pngAdvancedLinks = pngLinks
  x = log2(2^log2Signal[isValid, ] + param$backgroundExpression)
  xNormed = sweep(x, 1 , rowMeans(x));
  xSd = apply(x, 1, sd, na.rm=TRUE)
  ord = order(xSd, decreasing=TRUE)
  topGenes = ord[1:min(param$topGeneSize, length(ord))]
  
  pngName = "sampleCorrelation-AllPresent.png"
  plotCmd = expression({
    zValues = cor(x, use="complete.obs")
    ezCorrelationPlot(zValues, cond=conds, condLabels=conds,
                      main=paste0("all present genes (", sum(isValid), ")"), colors=sampleColors)
  })
  height = nrow(cor(x, use="complete.obs")) * 20
  if (height < 500) height = 500
  if (height > 2000) height = 2000
  try({
    pngLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=round(height*1.25), height=height)
  })
  
  pngName = "sampleCorrelation-AllPresentNormalized.png"
  plotCmd = expression({
    zValues = cor(xNormed, use="complete.obs")
    ezCorrelationPlot(zValues, cond=conds, condLabels=conds,
                      main=paste0("all present genes (", sum(isValid), ") gene-wise normalized"))
  })
  height = nrow(cor(xNormed, use="complete.obs")) * 20
  if (height < 500) height = 500
  if (height > 2000) height = 2000
  try({
    pngAdvancedLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=round(height*1.25), height=height)
  })
  
  pngName = "sampleCorrelation-TopGenes.png"
  plotCmd = expression({
    zValues = cor(x[topGenes,], use="complete.obs")
    ezCorrelationPlot(zValues, cond=conds, condLabels=conds,
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
    zValues = cor(xNormed[topGenes,], use="complete.obs")
    ezCorrelationPlot(zValues, cond=conds, condLabels=conds,
                      main=paste0("top genes (", length(topGenes), ") gene-wise normalized"))
  })
  height = nrow(cor(xNormed[topGenes,], use="complete.obs")) * 20
  if (height < 500) height= 500
  if (height > 2000) height = 2000
  try({
    pngAdvancedLinks[1, 2] = ezImageFileLink(plotCmd, file=pngName, width=round(height*1.25), height=height)
  })
  
  titles[["Sample correlation"]] = "Sample correlation"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addFlexTable(doc, ezGrid(pngLinks))
  
  advancedTitles = list()
  advancedTitles[["Advanced Plots"]] = "Advanced Plots"
  advancedDoc = openBsdocReport(title=advancedTitles[[length(advancedTitles)]])
  advancedTitles[["Sample correlation"]] = "Sample correlation"
  addTitle(advancedDoc, advancedTitles[[length(advancedTitles)]], 2, id=advancedTitles[[length(advancedTitles)]])
  addFlexTable(advancedDoc, ezGrid(pngAdvancedLinks))
  
  ################################
  ## cluster plots
  
  if (nSamples > 3){
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
    pngLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=min(600 + (nSamples-10)*20, 2000)) # nSamples dependent width
    
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
    pngAdvancedLinks[1, 1] = ezImageFileLink(plotCmd, file=pngName, width=min(600 + (nSamples-10)*20, 2000)) # nSamples dependent width
    
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
    pngLinks[1, 2] = ezImageFileLink(plotCmd, file=pngName, width=min(600 + (nSamples-10)*20, 2000)) # nSamples dependent width
    
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
    pngAdvancedLinks[1, 2] = ezImageFileLink(plotCmd, file=pngName, width=min(600 + (nSamples-10)*20, 2000)) # nSamples dependent width
    
    titles[["Sample Clustering"]] = "Sample Clustering"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    addFlexTable(doc, ezGrid(pngLinks))
    
    advancedTitles[["Sample Clustering"]] = "Sample Clustering"
    addTitle(advancedDoc, advancedTitles[[length(advancedTitles)]], 2, id=advancedTitles[[length(advancedTitles)]])
    addFlexTable(advancedDoc, ezGrid(pngAdvancedLinks))
    
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
      clusterLink = ezImageFileLink(plotCmd, file=clusterPng, width=max(800, 400 + 10 * ncol(xNormed[use, ])), height=1000) # HEATMAP
      
      if (doGo(param, seqAnno)){
        clusterResult = goClusterResults(xNormed[use, ], param, clusterResult, seqAnno=seqAnno,
                                         universeProbeIds=rownames(seqAnno))
      }
      
      titles[["Clustering of High Variance Features"]] = "Clustering of High Variance Features"
      addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
      addParagraph(doc, paste("Threshold for std. dev. of log2 signal across samples:", sdThresh))
      addParagraph(doc, paste("Number of features with high std. dev.:", sum(use)))
      
      jsFile = system.file("extdata/enrichr.js", package="ezRun", mustWork=TRUE)
      addJavascript(doc, jsFile)
      if (!is.null(clusterResult$GO)){
        goTables = goClusterTable(param, clusterResult, seqAnno)
        # addFlexTable(doc, ezFlexTable(goTables$linkTable, add.rownames=TRUE))
        if (doEnrichr(param)){
          goAndEnrichr = cbind(goTables$linkTable, goTables$enrichrTable)
        } else {
          goAndEnrichr = goTables$linkTable
        }
        goAndEnrichrFt = ezFlexTable(goAndEnrichr, border = 2, header.columns = TRUE, add.rownames=TRUE)
        bgColors = rep(gsub("FF$", "", clusterResult$clusterColors))
        goAndEnrichrFt = setFlexTableBackgroundColors(goAndEnrichrFt, j=1, colors=bgColors)
        goAndEnrichrTableLink = as.html(ezGrid(rbind("Background color corresponds to the row colors in the heatmap plot.",
                                                     as.html(goAndEnrichrFt))))
        goLink = as.html(ezGrid(rbind("Background color corresponds to the row colors in the heatmap plot.",
                                      as.html(goTables$ft))))
      } else {
        goAndEnrichrTableLink = as.html(pot("No information available"))
        goLink =  as.html(pot("No information available"))
      }
      goClusterTableDoc = openBsdocReport("GO Cluster tables")
      tbl = ezGrid(cbind("Cluster Plot"=clusterLink, "GO categories of feature clusters"=goLink), header.columns = TRUE)
      addFlexTable(goClusterTableDoc, tbl)
      closeBsdocReport(goClusterTableDoc, "goClusterTable.html")
      addParagraph(doc, pot("GO cluster tables", hyperlink="goClusterTable.html"))
      tbl = ezGrid(cbind("Cluster Plot"=clusterLink, "GO categories of feature clusters"=goAndEnrichrTableLink), header.columns = TRUE)
      addFlexTable(doc, tbl)
    }
    
    ##########################################
    ## mds plot
    
    titles[["MDS-Plot"]] = "MDS-Plot"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    pngLinks=vector(length=2)
    
    pngName = "mdsPlot_PresentGenes.png"
    plotCmd = expression({
      ezMdsPlot(signal=x, sampleColors=sampleColors, main=sub('.png','',pngName))
    })
    presentLink = ezImageFileLink(plotCmd, file=pngName)
    
    pngName = "mdsPlot_TopGenes.png"
    plotCmd = expression({
      ezMdsPlot(signal=x[topGenes,], sampleColors=sampleColors, main=sub('.png','',pngName))
    })
    topLink = ezImageFileLink(plotCmd, file=pngName)
    
    addFlexTable(doc, ezGrid(cbind(presentLink, topLink)))
    
    if (param$writeScatterPlots){
      qcScatterTitles = addQcScatterPlots(doc, param, design, conds, rawData,
                                          signalCond, isPresentCond, types=types)
      titles = append(titles, qcScatterTitles)
    }
    
    ##########################################
    ## count density plots
    
    pngName = "signalDens.png"
    plotCmd = expression({
      #countDensPlot(signal, sampleColors, main="all transcripts", bw=0.7)
      p = countDensGGPlot(cts=data.frame(signal,stringsAsFactors = F),
                          colors=sampleColors, alpha=0.4)
      print(p)
    })
    
    pngLink = ezImageFileLink(plotCmd, file=pngName, width=700, height=550)
    
    titles[["Expression densities"]] = "Expression densities"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    addParagraph(doc, "Zero or negative counts are not represented by the area!")
    addParagraph(doc, pngLink)
  }
  closeBsdocReport(advancedDoc, "advancedPlots.html", advancedTitles)
  addParagraph(doc, pot("advanced Plots", hyperlink = "advancedPlots.html"))
  closeBsdocReport(doc, htmlFile, titles)
}
