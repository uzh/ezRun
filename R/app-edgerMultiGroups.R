###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Edger Multi
##' @template htmlFile-template
##' @seealso \code{\link{EzAppEdgerMulti}}
ezMethodEdgerMulti = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(basename(output$getColumn("Report")))
  stopifnot(param$sampleGroup != param$refGroup)
  
  if (is.null(param$runGO)){
    param$runGO = FALSE
  }
  
  input = cleanupMultiGroupsInput(input, param)
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$batch) && length(param$batch) == 1){
    param$batch = input$meta[[param$batch]]
  }
  param$comparison = paste("glm fit for", param$grouping)
  if (ezIsSpecified(param$batch)){
    param$comparison = paste(param$comparison, "with second factor", param$batch)
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  
  result = multiGroupCountComparison(rawData, param)
  if (isError(result)){
    writeErrorReport(htmlFile, param=param, error=result$error)
    return("Error")
  }
  result$featureLevel = rawData$featureLevel
  result$countName = rawData$countName
  
  writeNgsMultiGroupReport(input$meta, result, htmlFile, param=param, rawData=rawData)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodEdgerMulti()
##' @seealso \code{\link{ezMethodEdgerMulti}}
EzAppEdgerMulti <-
  setRefClass("EzAppEdgerMulti",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodEdgerMulti
                  name <<- "EzAppEdgerMulti"
                  appDefaults <<- rbind(testMethod=ezFrame(Type="character",  DefaultValue="glm",  Description="which test method in edgeR to use: glm or exactTest"),
                                        normMethod=ezFrame(Type="character", DefaultValue="TMM", Description="edgeR's norm method: TMM, upperquartile, RLE, or none"))
                }
              )
  )


# functions for multiGroups.R
cleanupMultiGroupsInput = function(input, param){
  dataset = input$meta
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset, param), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]  ## CHECK: should FALSE be there?
  }
  inputMod = EzDataset(meta=dataset)
  if (!is.null(param$markOutliers) && param$markOutliers){
    stopifnot(!is.null(dataset$Outlier))
    grouping = inputMod$getColumn(param$grouping)
    isOut = dataset$Outlier %in% c("", "NO", '""', "FALSE") == FALSE
    grouping[isOut] = paste(grouping[isOut], "OUTLIER", sep="_")
    inputMod$setColumn(grouping)
  }
  return(inputMod)
}

multiGroupCountComparison = function(rawData , param){
  x = rawData$counts
  presentFlag = rawData$presentFlag
  
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
  if (!ezIsSpecified(batch)){
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
  logFC = lrt$table$logFC
  logFC$maxChange = apply(logFC, 1, max) - apply(logFC, 1, min) ## doesn't work and I'm not sure what it's supposed to do? perhaps the input data doesn't make sense.
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

writeNgsMultiGroupReport = function(dataset, result, htmlFile, param=NA, rawData=NA, types=NULL) {
  seqAnno = rawData$seqAnno
  
  titles = list()
  titles[["Analysis"]] = paste("Analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  
  titles[["Parameters"]] = "Parameters"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addDataset(doc, dataset, param)
  
  titles[["Result Summary"]] = "Result Summary"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addCountResultSummary(doc, param, result)
  
  titles[["Significant Counts"]] = "Significant Counts"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addSignificantCounts(doc, result)
  
  resultFile = addResultFile(doc, param, result, rawData)
  
  logSignal = log2(shiftZeros(result$xNorm, param$minSignal))
  result$groupMeans = averageColumns(logSignal, by=param$grouping)
  
  if (param$writeScatterPlots){
    testScatterTitles = addTestScatterPlots(doc, param, logSignal, result, seqAnno, types)
    titles = append(titles, testScatterTitles)
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
    clusterResult = clusterResults(xCentered, nClusters=6, clusterColors=clusterColors)
    plotCmd = expression({
      clusterHeatmap(xCentered, param, clusterResult, file=clusterPng,
                     colColors=sampleColors, lim=c(-param$logColorRange, param$logColorRange))
    })
    clusterLink = ezImageFileLink(plotCmd, file=clusterPng, width=max(800, 400 + 10 * ncol(xCentered)), height=1000) # HEATMAP
    
    if (doGo(param, seqAnno)){
      clusterResult = goClusterResults(xCentered, param, clusterResult, seqAnno=seqAnno,
                                       universeProbeIds=rownames(seqAnno)[result$isPresentProbe])
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
    addParagraph(doc, paste("Significance threshold:", param$pValueHighlightThresh))
    if (param$log2RatioHighlightThresh > 0){
      addParagraph(doc, paste("log2 Ratio threshold:", param$log2RatioHighlightThresh))
    }
    addParagraph(doc, paste("Number of significant features:", sum(use)))
    
    if (!is.null(clusterResult$GO)){
      goTables = goClusterTable(param, clusterResult)
      addFlexTable(doc, ezFlexTable(goTables$linkTable, add.rownames=TRUE))
      goLink = as.html(ezGrid(cbind("Background color corresponds to the row colors in the heatmap plot.",
                                    as.html(goTables$ft))))
      #goLink[[2]] = addGOClusterResult(doc, param, clusterResult)
    } else {
      goLink = as.html(pot("No information available"))
    }
    tbl = ezGrid(cbind("Cluster Plot"=clusterLink, "GO categories of feature clusters"=goLink), header.columns = TRUE)
    addFlexTable(doc, tbl)
  }
  ## only do GO if we have enough genes
  if (doGo(param, seqAnno)){
    goResult = twoGroupsGO(param, result, seqAnno, normalizedAvgSignal=rowMeans(result$groupMeans), method=param$goseqMethod) ## should probably be multiGroupsGO()
    titles[["GO Enrichment Analysis"]] = "GO Enrichment Analysis"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    addGoUpDownResult(doc, param, goResult)
  } 
  closeBsdocReport(doc, htmlFile, titles)
}

