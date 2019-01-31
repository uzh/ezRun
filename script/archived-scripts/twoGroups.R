##' @title Writes the report for the two group analysis
##' @description Writes the report for the two group analysis.
##' @template dataset-template
##' @template result-template
##' @template output-template
##' @template htmlFile-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{logColorRange}{ an integer or numeric specifying the log color range.}
##'   \item{minSignal}{ a numeric or integer specifying the minimal signal amount.}
##'   \item{grouping}{ a character specifying the grouping.}
##'   \item{sampleGroup}{ a character specifying the group to sample.}
##'   \item{refGroup}{ a character specifying the reference group.}
##'   \item{pValueHighlightThresh}{ a numeric specifying the threshold for highlighting p-values.}
##'   \item{log2RatioHighlightThresh}{ a numeric specifying the threshold for highlighting log2 ratios.}
##'   \item{maxGenesForClustering}{ an integer specifying the maximum amount of genes for clustering. If the amount is higher, the least significant genes get removed first.}
##'   \item{minGenesForClustering}{ an integer specifying the minimum amount of genes for clustering. If the amount is lower, no clustering will be done.}
##'   \item{doZip}{ a logical indicating whether to archive the result file.}
##'   \item{goseqMethod}{ a character specifying the method for the GO analysis. Accepted values: Wallenius, Sampling or Hypergeometric.}
##'   \item{maxNumberGroupsDisplayed}{ an integer specifying the maximum amount of rows for each result file of the GO analysis.}
##' }
##' @template rawData-template
##' @template types-template
##' @template roxygen-template
writeNgsTwoGroupReport = function(deResult, output, 
                                  htmlFile="00index.html", types=NULL) {
  param = deResult$param
  se <- deResult$se
  
  dataset <- setNames(as.data.frame(colData(se)),
                      colnames(colData(se)))
  
  # TODO: this is temporary fix to make it compatible with GO and GAGE analysis
  ## In the future, we should always use SummarizedExperiement rawDtaa.
  rawData <- deResult$rawData
  
  ## TODO: this is temporary fix to make it compatible with GO and GAGE analysis
  result = list()
  result$analysis = metadata(se)$analysis
  result$method = metadata(se)$method
  result$isPresentProbe = rowData(se)$isPresentProbe
  result$isPresent = assays(se)$isPresent
  result$log2Ratio = setNames(rowData(se)$log2Ratio, 
                              rownames(assays(se)$counts))
  result$fitGlm = metadata(se)$fitGlm
  result$sf = colData(se)$sf
  result$gfold = rowData(se)$gfold
  result$nativeResult = metadata(se)$nativeResult
  result$pValue = setNames(rowData(se)$pValue,
                           rownames(assays(se)$counts))
  result$fdr = rowData(se)$fdr
  result$usedInTest = rowData(se)$usedInTest
  result$xNorm = assays(se)$xNorm
  result$featureLevel <- metadata(se)$featureLeve
  result$countName <- metadata(se)$countName
  result$summary <- metadata(se)$summary
  deResult$result <- result  ## For deResult output
  
  seqAnno <- data.frame(rowData(se), row.names=rownames(se),
                        check.names = FALSE, stringsAsFactors=FALSE)
  
  titles = list()
  titles[["Analysis"]] = paste("Analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  
  addDataset(doc, dataset, param)
  addCountResultSummarySE(doc, param, se)
  
  titles[["Result Summary"]] = "Result Summary"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  settings = character()
  settings["Number of features:"] = length(rowData(se)$pValue)
  if (!is.null(rowData(se)$isPresentProbe)){
    settings["Number of features with counts above threshold:"] = 
      sum(rowData(se)$isPresentProbe)
  }
  addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
  
  titles[["Number of significants by p-value and fold-change"]] = "Number of significants by p-value and fold-change"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addSignificantCountsSE(doc, se)
  
  resultFile = addResultFileSE(doc, param, se)
  ezWrite.table(colData(se)$sf, file="scalingfactors.txt", head="Name", digits=4)
  
  liveReportLink = output$getColumn("Live Report")
  #resultObjFile = paste0("result--", param$comparison, "--", ezRandomString(length=12), "--EzResult.RData")
  deResult$saveToFile(basename(output$getColumn("Live Report")))
  addParagraph(doc, ezLink(liveReportLink,
                           "Live Report and Visualizations",
                           target = "_blank"))
  
  ## for scatter plots we show the highly variable low counts
  logSignal = log2(shiftZeros(assays(se)$xNorm, param$minSignal))
  result$groupMeans = cbind(rowMeans(logSignal[ , param$grouping == param$sampleGroup, drop=FALSE]),
                            rowMeans(logSignal[ , param$grouping == param$refGroup, drop=FALSE]))
  colnames(result$groupMeans) = c(param$sampleGroup, param$refGroup)
  
  if (param$writeScatterPlots){
    testScatterTitles = addTestScatterPlotsSE(doc, param, logSignal, se,
                                              resultFile$resultFile, types)
    titles = append(titles, testScatterTitles)
  }
  
  use = rowData(se)$pValue < param$pValueHighlightThresh & 
    abs(rowData(se)$log2Ratio) > param$log2RatioHighlightThresh & rowData(se)$usedInTest
  use[is.na(use)] = FALSE
  if (sum(use) > param$maxGenesForClustering){
    use[use] = rank(rowData(se)$pValue[use], ties.method="max") <= param$maxGenesForClustering
  }
  ## for clustering we use a moderated logSignal
  logSignal = log2(assays(se)$xNorm + param$backgroundExpression)
  
  if (sum(use) > param$minGenesForClustering){
    xCentered = logSignal[use , ]
    if (!is.null(param$useRefGroupAsBaseline) && param$useRefGroupAsBaseline){
      xCentered = xCentered - rowMeans(xCentered[ , param$grouping == param$refGroup])
    } else {
      xCentered = xCentered - rowMeans(xCentered)
    }
    xCentered = xCentered[, order(param$grouping)]
    sampleColors = getSampleColors(param$grouping)[order(param$grouping)]
    clusterPng = "cluster-heatmap.png"
    clusterColors = c("red", "yellow", "orange", "green", "blue", "cyan")
    clusterResult = clusterResults(xCentered, nClusters=6, clusterColors=clusterColors)
    plotCmd = expression({
      clusterHeatmap(xCentered, param, clusterResult, file=clusterPng,
                     colColors=sampleColors, lim=c(-param$logColorRange, param$logColorRange))
    })
    clusterLink = ezImageFileLink(plotCmd, file=clusterPng, width=max(800, 400 + 10 * ncol(xCentered)), height=1000) # HEATMAP
    
    if (doGo(param, seqAnno)){
      clusterResult = goClusterResults(xCentered, param, clusterResult, seqAnno=seqAnno,
                                       universeProbeIds=rownames(seqAnno)[rowData(se)$isPresentProbe])
    }
    
    ## append the result file with the cluster colors
    resultLoaded = ezRead.table(resultFile$resultFile)
    resultLoaded$Cluster = clusterResult$clusterColors[clusterResult$clusterNumbers[rownames(resultLoaded)]]
    ezWrite.table(resultLoaded, file=resultFile$resultFile)
    if (param$doZip){
      zipFile(resultFile$resultFile)
    }
    titles[["Clustering of Significant Features"]] = "Clustering of Significant Features"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    paragraphs = character()
    paragraphs["Significance threshold:"] = param$pValueHighlightThresh
    if (param$log2RatioHighlightThresh > 0){
      paragraphs["log2 Ratio threshold:"] = param$log2RatioHighlightThresh
    }
    paragraphs["Number of significant features:"] = sum(use)
    addFlexTable(doc, ezGrid(paragraphs, add.rownames=TRUE))
    
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
      bgColors = gsub("FF$", "", clusterResult$clusterColors)
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
  ## only do GO if we have enough genes
  if (doGo(param, seqAnno)){
    if (doEnrichr(param)){
      isSig = rowData(se)$pValue < param$pValThreshGO & rowData(se)$usedInTest
      isUp = rowData(se)$log2Ratio > param$log2RatioThreshGO & isSig
      isDown = rowData(se)$log2Ratio < -param$log2RatioThreshGO & isSig
      regulatedGenes = list()
      regulatedGenes$upGenes = na.omit(unique(seqAnno[isUp, "gene_name"]))
      regulatedGenes$downGenes = na.omit(unique(seqAnno[isDown, "gene_name"]))
      regulatedGenes$bothGenes = union(regulatedGenes$upGenes, regulatedGenes$downGenes)
      
      enrichrLinks = ezMatrix("", rows=names(regulatedGenes), cols=c("External", "Precomputed"))
      maxResultsPerLibrary = 5
      for (row in rownames(enrichrLinks)){
        genesToUse = regulatedGenes[[row]]
        jsCall = paste0('enrich({list: "', paste(genesToUse, collapse="\\n"), '", popup: true});')
        enrichrLinks[row, "External"] = as.html(pot(paste0("<a href='javascript:void(0)' onClick='", jsCall, "'>Analyse at Enrichr website</a>")))
        resMerged <- NA
        if (!is.null(genesToUse) && length(genesToUse) > 3 && param$doPrecomputeEnrichr) {
          resList = runEnrichr(genesToUse, maxResult = maxResultsPerLibrary)
          resList = lapply(names(resList), function(nm){return(cbind("Gene-set library"=nm, resList[[nm]][, c(2:5, 7:10)]))}) ## add the name as a first column
          if (length(resList) > 0) {
            resMerged = do.call("rbind", resList)
            resMerged <- resMerged[order(-resMerged[,5]), ]
            resMerged[, c(3,6:8)] <- apply(resMerged[, c(3,6:8)], 2, sprintf, fmt = "%0.2e")
            resMerged[, c(4,5)] <- apply(resMerged[, c(4,5)], 2, sprintf, fmt = "%0.3f")
            enrichrTablePath <- paste0("enrichrTable_", row, ".html")
            ezInteractiveTable(resMerged, tableLink=enrichrTablePath, title=paste("Enrichr report for ", row))
            enrichrLinks[row, "Precomputed"] = as.html(ezLink(link = enrichrTablePath, label = "Report", target = "_blank"))
          } else {
            enrichrLinks[row, "Precomputed"] = "No significant results"
          }
        } else {
          enrichrLinks[row, "Precomputed"] = "Not run (too few genes)"
        }
      }
      titles[["Enrichr"]] = "Enrichr"
      addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
      addFlexTable(doc, ezFlexTable(enrichrLinks, valign="middle", add.rownames = TRUE, header.columns = T))
    }
    
    goResult = twoGroupsGO(param, result, seqAnno, normalizedAvgSignal=rowMeans(result$groupMeans), method=param$goseqMethod)
    titles[["GO Enrichment Analysis"]] = "GO Enrichment Analysis"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    revigoTitle = addGoUpDownResult(doc, param, goResult)
    titles = append(titles, revigoTitle)
    
  }
  
  ## Run Gage
  if(param[['GAGEanalysis']] ) {
    gageRes = runGageAnalysis(result, param=param, output=output, rawData=rawData)
    titles[["GAGE Enrichment Analysis"]] = "GAGE Enrichment Analysis"
    addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
    addGageTables(doc, param, gageRes)
  }
  addParagraph(doc, pot("advanced Plots", hyperlink = "advancedPlots.html"))
  closeBsdocReport(doc, htmlFile, titles)
}


twoGrouplsRoast = function(param, testResult, seqAnno){
  job = ezJobStart("twoGroupsRoast")
  require("GOstats", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  require("annotate", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  if (param$featureLevel != "gene"){
    stop("only feature level: gene is supported")
  }
  ontologies = c("BP", "MF", "CC")
  #goResults = list()
  #for (onto in ontologies){
  goResults = ezMclapply(ontologies, function(onto){
    gene2goList = goStringsToList(seqAnno[[paste("GO", onto)]], listNames=rownames(seqAnno))[presentGenes]
    if (param$includeGoParentAnnotation){
    }
  })
}
