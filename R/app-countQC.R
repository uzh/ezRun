###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Count QC
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
    writeErrorHtml(htmlFile, param=param, dataset=dataset, error=rawData$error)
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

runNgsCountQC = function(dataset, htmlFile="00index.html", param=param, rawData=NULL,
 writeDataFiles=TRUE, types=NULL, annoFile=NULL){


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

  html = openHtmlReport(htmlFile, param=param, title=paste("Analysis:", param$name),
                        dataset=dataset)  
  on.exit(closeHTML(html))
  
  
  if (nSamples < 2){
    ezWrite("<h2>Note: Statistics and Plots are not available for single sample experiments</h2>", con=html)
    ezWrite("<p>Run the report again with multiple samples selected.</p>", con=html)
    return("Success");
  }

  signal = shiftZeros(getSignal(rawData), param$minSignal)
  presentFlag = rawData$presentFlag
  signalRange = range(signal, na.rm=TRUE)
  log2Signal = log2(signal)
  isPresent = ezPresentFlags(signal, presentFlag=presentFlag, param=param, isLog=rawData$isLog)
  signalCond = 2^averageColumns(log2Signal, by=conds)
  isPresentCond = averageColumns(isPresent, by=conds) >= 0.5
  isPresentStudy = apply(isPresentCond, 1, mean) >= 0.5
  
  
  ezWrite("<h2>Parameters</h2>", con=html)
	ezWrite("<table border='0'>", con=html)
  ezWrite("<tr><td>Normalization method:</td><td>", param$normMethod, "</td></tr>", con=html)
  if (param$useSigThresh){
     ezWrite("<tr><td>Log2 signal threshold:</td><td>", signif(log2(param$sigThresh), digits=4), "</td></tr>", con=html)
     ezWrite("<tr><td>Linear signal threshold:</td><td>", signif(param$sigThresh, digits=4), "</td></tr>", con=html)
  }
  ezWrite("<tr><td>Feature level:</td><td>", rawData$featureLevel, "</td></tr>", con=html)
  ezWrite("<tr><td>Number of features:</td><td>", nrow(signal), "</td></tr>", con=html)
  ezWrite("<tr><td>Data Column Used:</td><td>", rawData$countName, "</td></tr>", con=html)
  ezWrite("</table>", con=html)

  if (writeDataFiles){
    if (!is.null(rawData$presentFlag)){
      combined = interleaveMatricesByColumn(rawData$signal, rawData$presentFlag)			
    } else {
      combined = rawData$signal
    }
    if (!is.null(seqAnno)){
      combined = cbind(seqAnno[rownames(combined), ,drop=FALSE], combined)
    }
    countFile = paste(ezValidFilename(param$name), "-raw-count.txt", sep="")
    ezWrite.table(rawData$counts, file=countFile, head="Feature ID", digits=4)
    signalFile = paste(ezValidFilename(param$name), "-normalized-signal.txt", sep="")
    ezWrite.table(combined, file=signalFile, head="Feature ID", digits=4)
    rpkmFile = paste(ezValidFilename(param$name), "-rpkm.txt", sep="")
    ezWrite.table(getRpkm(rawData), file=rpkmFile, head="Feature ID", digits=4)
    if (param$doZip){
      dataFiles =c(zipFile(countFile), zipFile(signalFile), zipFile(rpkmFile))
      ezWrite("<h2>Data Files</h2>", con=html)
      writeTxtLinksToHtml(dataFiles, con=html, mime="application/zip")
    } else {
      ezWrite("<h2>Data Files</h2>", con=html)
      writeTxtLinksToHtml(c(countFile, signalFile, rpkmFile), con=html, mime="application/txt")
    }
  }
  
  ezWrite("<h2>Count Statistics</h2>", con=html)
	totalCounts = signif(apply(rawData$counts, 2, sum) / 1e6, digits=3)
  presentCounts = apply(isPresent, 2, sum)
  presentFrame = data.frame("Total Read Counts [Mio]"=totalCounts,
														"Genomic Features <br>with Reads above threshold"=presentCounts,
														"Genomic Features <br>with Reads above threshold [%]"=signif(100*presentCounts/nrow(isPresent), digits=3),
														check.names=FALSE)
  rownames(presentFrame) = samples
  writeTableToHtml(presentFrame, con=html)

  rawData$signal = signal



  #################################
  ## correlation plot

  isValid = isPresentStudy
  if (!is.null(seqAnno$IsControl)){
   isValid = isValid & !seqAnno$IsControl
  }

  pngNames = ezMatrix("", rows=1:2, cols=1:2)
  x = log2(2^log2Signal[isValid, ] + param$bgExpression)
  xNormed = sweep(x, 1 , rowMeans(x));
  xSd = apply(x, 1, sd, na.rm=TRUE)
  ord = order(xSd, decreasing=TRUE)
  topGenes = ord[1:min(param$topGeneSize, length(ord))]

  pngName = "sampleCorrelation-AllPresent.png"
  try({ezCorrelationPlot(png=pngName, cor(x, use="complete.obs"), cond=conds, condLabels=conds, 
    main=paste("all present genes (", sum(isValid), ")", sep=""), sampleColors=sampleColors)})
  pngNames[1, 1] = pngName

  pngName = "sampleCorrelation-AllPresentNormalized.png"
  try({ezCorrelationPlot(png=pngName, cor(xNormed, use="complete.obs"), cond=conds, condLabels=conds, 
                         main=paste("all present genes (", sum(isValid), ") gene-wise normalized", sep=""))})
  pngNames[1,2] = pngName

  pngName = "sampleCorrelation-TopGenes.png"
  try({ezCorrelationPlot(png=pngName, cor(x[topGenes,], use="complete.obs"), cond=conds, condLabels=conds, 
                         main=paste("top genes (", length(topGenes), ")", sep=""))})
  pngNames[2,1] = pngName

  pngName = "sampleCorrelation-TopGenesNormalized.png"
  try({ezCorrelationPlot(png=pngName, cor(xNormed[topGenes,], use="complete.obs"), cond=conds, condLabels=conds, 
                         main=paste("top genes (", length(topGenes), ") gene-wise normalized", sep=""))})
  pngNames[2,2] = pngName

  ezWrite("<h2>Sample correlation</h2>", con=html)
  writeImageTableToHtml(pngNames, con=html)
  flush(html)


  ################################
  ## cluster plots

  if(nSamples > 3){

    pngNames = ezMatrix("", rows=1:2, cols=1:2)

    pngName = "sampleClustering-AllPresent.png"
    pngNames[ 1, 1] = pngName
    png(filename=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    d = as.dist(1-cor(x, use="complete.obs"));
    hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
    hcd = colorClusterLabels(hcd, sampleColors)
    mai = par("mai")
    mai[1] = 3
    par(mai=mai)
    plot(hcd, main="all present genes", xlab="")
    dev.off()

    pngName = "sampleClustering-AllPresentNormalized.png"
    pngNames[ 1,2] = pngName
    png(filename=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    d = as.dist(1-cor(xNormed, use="complete.obs"));
    hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
    hcd = colorClusterLabels(hcd, sampleColors)
    mai = par("mai")
    mai[1] = 3
    par(mai=mai)
    plot(hcd, main="all present genes; gene-wise normalized", xlab="")
    dev.off()

    pngName = "sampleClustering-TopGenes.png"
    pngNames[ 2, 1] = pngName
    png(filename=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    d = as.dist(1-cor(x[topGenes, ], use="complete.obs"));
    hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
    hcd = colorClusterLabels(hcd, sampleColors)
    mai = par("mai")
    mai[1] = 3
    par(mai=mai)
    plot(hcd, main=paste("top", length(topGenes), " genes"), xlab="")
    dev.off()

    pngName = "sampleClustering-TopGenesNormalized.png"
    pngNames[ 2, 2] = pngName
    png(filename=pngName, width=800 + max(0, 10 * (nSamples-20)), height=500)
    d = as.dist(1-cor(xNormed[topGenes, ], use="complete.obs"));
    hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
    hcd = colorClusterLabels(hcd, sampleColors)
    mai = par("mai")
    mai[1] = 3
    par(mai=mai)
    plot(hcd, main=paste("top", length(topGenes), "genes; gene-wise normalized"), xlab="")
    dev.off()
    
    ezWrite("<h2>Sample Clustering</h2>", con=html)
    writeImageTableToHtml(pngNames, con=html)
    
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
		  clusterResult = clusterHeatmap(param, xNormed[use, ], file=clusterPng, nClusters=6, 
           lim=c(-param$logColorRange, param$logColorRange), margins=c(18, 9),
		       colColors=sampleColors, clusterColors=clusterColors, doClusterColumns=TRUE,
					 doGO=doGo(param, seqAnno), seqAnno=seqAnno,
		       universeProbeIds=rownames(seqAnno))
			ezWrite("<h2>Clustering of High Variance Features</h2>", con=html)
			ezWrite("<p>Threshold for std. dev. of log2 signal across samples: ", sdThresh, "<br>", con=html)
			ezWrite("Number of features with high std. dev.: ", sum(use), "</p>", con=html)
			ezWrite("<table border=0><tr><th>Cluster Plot</th><th>GO categories of feature clusters</th></tr>", con=html)
			ezWrite("<tr valign=top><td>", con=html)
		  writeImageRowToHtml(clusterPng, con=html)
		  ezWrite("</td><td>", con=html)
			if (!is.null(clusterResult$GO)){
   			ezWrite("Background color corresponds to the color of the feature cluster in the heatmap plot.<br>", con=html)
		    writeGOClusterResult(html, param, clusterResult)
			} else {
				ezWrite("No information available", con=html)
    	}
			ezWrite("</td></tr></table>", con=html)
		}
		##########################################
		###### mds plot
		ezWrite("<h2>MDS-Plot</h2>", con=html)
		pngNames=vector(length=2)
    pngNames[1] = "mdsPlot_PresentGenes.png"
		myMdsPlot(signal=x, sampleColors=sampleColors, pngNames[1])
    
    pngNames[2] = "mdsPlot_TopGenes.png"
		myMdsPlot(signal=x[topGenes,], sampleColors=sampleColors, pngNames[2])
		
    writeImageRowToHtml(pngNames, con=html)
    flush(html)
    
    if (param$writeScatterPlots){
      writeQcScatterPlots(html, param, design, conds,
                          rawData, signalCond, isPresentCond, types=types)
    }
		
		if (!is.null(rawData$countsStart) & !is.null(rawData$countsEnd)){
		  ezWrite("<h2>3' Bias analysis</h2>", con=html)
		  pngName = "start-end-countScatter.png"
		  valStart = shiftZeros(rawData$countsStart, param$minSignal)
		  colnames(valStart) = paste(colnames(valStart), "[transcript start]")
		  valEnd = shiftZeros(rawData$countsEnd, param$minSignal)
		  colnames(valEnd) = paste(colnames(valEnd), "[transcript end]")
		  ezScatter(x=valStart, y=valEnd, file=pngName, types=types)
		  writeImageRowToHtml(pngName, con=html)    
		  pngName = "start-end-countSmoothScatter.png"
		  ezSmoothScatter(x=valStart, y=valEnd, file=pngName)
		  writeImageRowToHtml(pngName, con=html)    
		}
		
		##########################################
		###### count density plots
		pngNames = character()
		pngName = "signalDens.png"
		pngNames["all"] = pngName
		countDensPlot(param, signal, pngName, sampleColors, main="all transcripts", bw=0.7)
		ezWrite("<h2>Expression densities</h2>", con=html)
		ezWrite("<p>Zero or negative counts are not represented by the area!</p>", con=html)
		writeImageRowToHtml(pngNames, con=html)
		flush(html)
		
  
  }
  ezSessionInfo()
  writeTxtLinksToHtml('sessionInfo.txt',con=html)
  flush(html)
}
