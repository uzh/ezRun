###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## TODO: REFAC into RC inheriting from EzApp
edgerMultiGroupApp = function(inputDatasetFile=NA, output=NA, param=NA){
  param = fillWithDefaults(param)
  on.exit({
    if(!is.null(param$mail) && grepl('@',param$mail)){
      text=paste0('http://fgcz-gstore.uzh.ch/projects/',output[['Report [File]']],'/00index.html')
      subject=paste0("EdgeR ", sub("/.*", "", output[['Report [File]']]),' done.')
      ezMail(text=text,subject=subject,to=param$mail)
    }
  })
  #waitUntilFileExists(inputDatasetFile, maxWaitSeconds=300)
  checkFreeDiskSpace(param)
  cwd = getwd()
  on.exit(setwd(cwd), add=TRUE)
  setwdNew(basename(output$Report))
  param$name = basename(output$Report)
  ngsMultiGroupAnalysis(inputDatasetFile,  param=param)
  return("Success")
}


ngsMultiGroupAnalysis = function(datasetFile, htmlFile="00index.html", param=NULL){
  
  param = fillWithDefaults(param)
  logMessage("ngsMultiGroupAnalysis", param, "Starting") ## TODO: on.exit() and logStart() will run in EzAPP that this app should eventually inherit from.
  on.exit({traceback(); logMessage("ngsMultiGroupAnalysis", param, "Finished")}) 
  
  if (is.null(param$runGO)){
    param$runGO = FALSE
  }
  
  dataset = ezRead.table(datasetFile)
  if (param$useFactorsAsSampleName){
    sampleName = apply(dataset[ , grepl("Factor", colnames(dataset)), drop=FALSE], 1, paste, collapse="_")
    rownames(dataset) = paste(sampleName, rownames(dataset))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""') == TRUE, ]
  }
  
  rawData = loadCountDataset(param, dataset) ## TODO: REFAC to loadCountDataset(input, param, dataset) or this call won't work
  param$comparison = paste("glm fit for", param$grouping)
  if (ezIsSpecified(param$batch)){
    param$comparison = paste(param$comparison, "with second factor", param$batch)
  }
  
  if (!is.null(dataset[[param$grouping]])){
    param$grouping = dataset[[param$grouping]]    
  } else {
    if (!is.null(dataset[[paste(param$grouping, "[Factor]")]])){
      param$grouping = dataset[[paste(param$grouping, "[Factor]")]]    
    } else {
      stop("column not found: ", param$grouping)
    }
  }
  if (!is.null(dataset[[param$batch]])){
    param$batch = dataset[[param$batch]]    
  } else {
    if (!is.null(dataset[[paste(param$batch, "[Factor]")]])){
      param$batch = dataset[[paste(param$batch, "[Factor]")]]    
    } else {
      stop("column not found: ", param$batch)
    }
  }
  
  result = runNgsMultiGroupAnalysis(dataset, rawData=rawData, param=param)
  return(result)
}

runNgsMultiGroupAnalysis = function(dataset, htmlFile="00index.html", param=param, rawData=NULL, types=NULL){
  
  if (isError(rawData)){
    writeErrorHtml(htmlFile, param=param, error=rawData$error)
    return("Error")
  }  
  x = rawData$counts
  result = multiGroupCountComparison(x, rawData$presentFlag, param)
  if (isError(result)){
    writeErrorHtml(htmlFile, param=param, error=rawData$error)
    return("Error")
  }  
  result$featureLevel = rawData$featureLevel
  result$countName = rawData$countName
  writeNgsMultiGroupReport(dataset, result, htmlFile, param=param, rawData=rawData, types=types)
}


writeNgsMultiGroupReport = function(dataset, result, htmlFile, param=NA, rawData=NA, types=NULL) {
  keggOrganism = NA
  seqAnno = rawData$seqAnno
  html = openHtmlReport(htmlFile, param=param, title=paste("Analysis:", param$name),
                        dataset=dataset)  
  on.exit(closeHTML(html))
  
  writeCountResultSummary(html, param, result)                ## REFAC addCountResultSummary
  writeResultCounts(html, param, result)                      ## REFAC addSignificantCounts
  resultFile = writeResultFile(html, param, result, rawData)  ## REFAC addResultFile
  
  cr = c(-param$logColorRange, param$logColorRange)
  logSignal = log2(shiftZeros(result$xNorm, param$minSignal))
  result$groupMeans = averageColumns(logSignal, by=param$grouping)
  
  if (param$writeScatterPlots){
    writeTestScatterPlots(html, param, logSignal, result, seqAnno, types=types, colorRange=cr) ## REFAC addTestScatterPlots
  }
  
  use = result$pValue < param$pValueHighlightThresh & abs(result$log2Ratio) > param$log2RatioHighlightThresh & result$usedInTest
  use[is.na(use)] = FALSE
  if (sum(use) > param$maxGenesForClustering){
    use[use] = rank(result$pValue[use], ties.method="max") <= param$maxGenesForClustering
  }
  
  if (sum(use) > param$minGenesForClustering){
    xCentered = (logSignal - rowMeans(logSignal))[use, order(param$grouping)]
    sampleColors = getSampleColors(param$grouping)[order(param$grouping)]
    clusterPng = "cluster-heatmap.png"
    clusterColors = c("red", "yellow", "orange", "green", "blue", "cyan")
    doGO = doGo(param, seqAnno)
    clusterResult = clusterHeatmap(param, xCentered, file=clusterPng, nClusters=6, ## REFAC (see clusterHeatmap in twoGroups)
                                   lim=c(-param$logColorRange, param$logColorRange),
                                   colColors=sampleColors, clusterColors=clusterColors,
                                   doGO=doGO, seqAnno=seqAnno,
                                   universeProbeIds=rownames(seqAnno)[result$isPresentProbe],
                                   keggOrganism=keggOrganism)
    
    ## append the result file with the cluster colors
    resultLoaded = ezRead.table(resultFile$resultFile)
    resultLoaded$Cluster = clusterResult$clusterColors[clusterResult$clusterNumbers[rownames(resultLoaded)]]
    ezWrite.table(resultLoaded, file=resultFile$resultFile)
    if (param$doZip){
      zipFile(resultFile$resultFile)
    }
    ezWrite("<h2>Clustering of Significant Features</h2>", con=html)
    ezWrite("<p>Significance threshold: ", param$pValueHighlightThresh, "<br>", con=html)
    if (param$log2RatioHighlightThresh > 0){
      ezWrite("log2 Ratio threshold: ", param$log2RatioHighlightThresh, "<br>", con=html)      
    }
    ezWrite("Number of significant features: ", sum(use), "</p>", con=html)
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
    if (!is.null(clusterResult$Kegg)){
      writeKeggClusterResult(html, param, clusterResult, keggOrganism)
    }
    ezWrite("</td></tr></table>", con=html)
  }
  ## only do GO if we have enough genes
  if (doGo(param, seqAnno)){
    goResult = twoGroupsGO(param, result, seqAnno, normalizedAvgSignal=rowMeans(result$groupMeans), method=param$goseqMethod)
    writeGOTables(html, param, goResult) ## REFAC addGoUpDownResult
  } 
  ezSessionInfo()
  writeTxtLinksToHtml('sessionInfo.txt',con=html)
}



multiGroupCountComparison = function(x, presentFlag=NULL, param){
  
  result = list()
  result$analysis ="NGS multi group analysis"
  if (!ezIsSpecified(param$testMethod)){
    param$testMethod = "glm"    
  }
  
  result$method = param$testMethod
  isRef = param$grouping == param$refGroup
  res = switch(param$testMethod,
               glm = runEdgerGlmMultiGroup(round(x), param$refGroup, param$grouping, param$normMethod, batch=param$batch),
               stop("unsupported testMethod: ", param$testMethod)
  )
  result$log2Ratio = res$log2FoldChange  
  pValue = res$pval
  pValue[is.na(pValue)] = 1
  
  ## compute which probes are present
  groups = unique(param$grouping)
  isPresent = ezPresentFlags(x, presentFlag=presentFlag, param=param, isLog=TRUE)
  useProbe = rep(FALSE, nrow(x))
  for (group in groups){
    #groupMeans[ , group] = apply(x[ , group == param$grouping], 1, mean, na.rm=TRUE)
    useProbe[ apply(isPresent[ , group == param$grouping, drop=FALSE], 1, mean) > 0.5] = TRUE
  }
  isPresentProbe = useProbe
  ## further quality filtering of useProbe, e.g. variance outlier, see two group analysis
  fdr = rep(NA, length(pValue))
  fdr[useProbe] = p.adjust(pValue[useProbe], method="fdr")
  
  names(pValue) = rownames(x)
  names(fdr) = rownames(x)
  
  result$pValue = pValue  
  result$fdr=fdr
  result$isPresent = isPresent
  result$isPresentProbe = isPresentProbe
  result$usedInTest = useProbe
  result$xNorm = ezScaleColumns(x, res$sf)
  return(result)
}



runDeseq2MultiGroup = function(x, sampleGroup, refGroup, grouping, batch=NULL){
  if (is.null(batch)){
    colData = data.frame(grouping=as.factor(grouping), row.names=colnames(x))
    dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping)
  } else {
    colData = data.frame(grouping=as.factor(grouping), batch=as.factor(batch), row.names=colnames(x))
    dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping + batch)
  }
  dds = DESeq2::DESeq(dds, quiet=FALSE)
  res = results(dds, contrast=c("grouping", sampleGroup, refGroup), cooksCutoff=FALSE)
  res=as.list(res)
  res$sf = 1/colData(dds)$sizeFactor
  return(res)
}




runEdgerGlmMultiGroup = function(x, refGroup, grouping, normMethod, batch=NULL){
  requireNamespace("edgeR", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  ## get the scaling factors for the entire data set
  cds = DGEList(counts=x, group=grouping)
  cds = calcNormFactors(cds, method=normMethod)
  sf = 1/(cds$samples$norm.factors * cds$samples$lib.size)
  sf = sf / ezGeomean(sf)
  #sf = ezLogmeanScalingFactor(x, presentFlag=x>0)
  
  groupFactor = factor(grouping, levels=c(refGroup, setdiff(grouping, refGroup)))
  if (is.null(batch)){
    design = model.matrix( ~ groupFactor)
    #colnames(design) = c("Intercept", "Grouping")
  } else {
    design = model.matrix( ~ groupFactor + factor(batch))
    #colnames(design) = c("Intercept", "Grouping", paste("Batch", 1:(ncol(design)-2), sep="_"))
  }
  
  ## dispersion estimation
  cds = estimateGLMTrendedDisp(cds, design)
  cds = estimateGLMTagwiseDisp(cds, design)
  
  ## Testing for DE genes
  fitGlm = glmFit(cds, design)
  lrt = glmLRT(fitGlm, coef=2:length(levels(groupFactor)))
  res = list()
  res$id = rownames(lrt$table)
  logFC = lrt$table[ , grep("^logFC", colnames(lrt$table))]
  logFC$maxChange = apply(logFC, 1, max) - apply(logFC, 1, min)
  res$log2FoldChange = apply(abs(logFC), 1, max)
  res$pval = lrt$table$PValue
  #res$log2Expr = lrt$table$logCPM
  res$sf = sf
  ## do not return groupMeans now.
  #res$groupMeans = cbind(apply(cds$count[ , grouping == sampleGroup, drop=FALSE], 1, mean), apply(cds$count[ , grouping == refGroup, drop=FALSE], 1, mean))
  #colnames(res$groupMeans) = c(sampleGroup, refGroup)
  #rownames(res$groupMeans) = rownames(cds$count)
  return(res)
}

