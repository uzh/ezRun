###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodEdgerMulti = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  setwdNew(basename(output$getColumn("Report")))
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupMultiGroupsInput(input, param)
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1){
    param$grouping2 = input$getColumn(param$grouping2)
  }
  param$comparison = paste("glm fit for", param$grouping)
  if (ezIsSpecified(param$grouping2)){
    param$comparison = paste(param$comparison, "with second factor", param$grouping2)
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
##' @templateVar method ezMethodEdgerMulti(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
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
                                        normMethod=ezFrame(Type="character", DefaultValue="TMM", Description="edgeR's norm method: TMM, upperquartile, RLE, or none"),
                                        runGO=ezFrame(Type="logical", DefaultValue=FALSE, Description="whether to run the GO analysis"))
                }
              )
  )

##' @title Cleans up input from the edgeR multi groups app
##' @description Cleans up input from  the edgeR multi groups app.
##' @inheritParams cleanupTwoGroupsInput
##' @template roxygen-template
##' @return Returns an object of the class EzDataset that is the modified \code{input}.
cleanupMultiGroupsInput = function(input, param){
  dataset = input$meta
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset, param), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]  ## CHECK: should FALSE be there?
  }
  inputMod = EzDataset(meta=dataset, dataRoot=param$dataRoot)
  if (!is.null(param$markOutliers) && param$markOutliers){
    stopifnot(!is.null(dataset$Outlier))
    grouping = inputMod$getColumn(param$grouping)
    isOut = dataset$Outlier %in% c("", "NO", '""', "FALSE") == FALSE
    grouping[isOut] = paste(grouping[isOut], "OUTLIER", sep="_")
    inputMod$setColumn(grouping)
  }
  return(inputMod)
}

##' @title Compares the counts of many groups
##' @description Compares the counts of many groups with the option to choose from several methods to test them.
##' @template rawData-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{testMethod}{ defines the method to run: deseq2, exactTest, glm, sam or limma. Defaults to glm.}
##'   \item{grouping}{ a character specifying the grouping.}
##'   \item{grouping2}{ a character vector specifying a secondary factor.}
##'   \item{refGroup}{ a character specifying the reference group.}
##'   \item{normMethod}{ a character specifying the normalization method for the edger and glm test methods.}
##' }
##' @template roxygen-template
##' @return Returns a list containing the results of the comparison.
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
               glm = runEdgerGlmMultiGroup(round(x), param$refGroup, param$grouping, param$normMethod, grouping2=param$grouping2),
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

##' @describeIn multiGroupCountComparison Runs the Deseq2 test method for many groups.
runDeseq2MultiGroup = function(x, sampleGroup, refGroup, grouping, grouping2=NULL){
  if (is.null(grouping2)){
    colData = data.frame(grouping=as.factor(grouping), row.names=colnames(x))
    dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping)
  } else {
    if (ezTagListFromNames(grouping2) == "Factor") {
    colData = data.frame(grouping=as.factor(grouping), grouping2=as.factor(grouping2), row.names=colnames(x))
    dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping + grouping2)
    } else if (ezTagListFromNames(grouping2) == "Numeric") {
      colData = data.frame(grouping=as.factor(grouping), grouping2=as.numeric(grouping2), row.names=colnames(x))
      dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping + grouping2)
    } else {
      stop("Column header of grouping2 must have the tag [Factor] or [Numeric]")
    }
  }
  dds = DESeq2::DESeq(dds, quiet=FALSE)
  res = DESeq2::results(dds, contrast=c("grouping", sampleGroup, refGroup), cooksCutoff=FALSE)
  res=as.list(res)
  res$sf = 1/colData(dds)$sizeFactor
  return(res)
}

##' @describeIn multiGroupCountComparison Runs the EdgeR GLM test method for many groups.
runEdgerGlmMultiGroup = function(x, refGroup, grouping, normMethod, grouping2=NULL){
  require("edgeR", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  ## get the scaling factors for the entire data set
  cds = DGEList(counts=x, group=grouping)
  cds = calcNormFactors(cds, method=normMethod)
  sf = 1/(cds$samples$norm.factors * cds$samples$lib.size)
  sf = sf / ezGeomean(sf)
  #sf = ezLogmeanScalingFactor(x, presentFlag=x>0)
  
  groupFactor = factor(grouping, levels=c(refGroup, setdiff(grouping, refGroup)))
  if (!ezIsSpecified(grouping2)){
    design = model.matrix( ~ groupFactor)
    #colnames(design) = c("Intercept", "Grouping")
  } else {
    if (ezTagListFromNames(grouping2) == "Factor") {
      design = model.matrix( ~ groupFactor + factor(grouping2))
    } else if (ezTagListFromNames(grouping2) == "Numeric") {
      design = model.matrix( ~ groupFactor + as.numeric(grouping2))
    } else {
      stop("Column header of grouping2 must have the tag [Factor] or [Numeric]")
    }
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

##' @title Writes the report for the multi group analysis
##' @description Writes the report for the multi group analysis.
##' @template dataset-template
##' @template result-template
##' @template htmlFile-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{name}{ a character specifying the name of the run.}
##'   \item{logColorRange}{ an integer or numeric specifying the log color range.}
##'   \item{minSignal}{ a numeric or integer specifying the minimal signal amount.}
##'   \item{grouping}{ a character specifying the grouping.}
##'   \item{pValueHighlightThresh}{ a numeric specifying the threshold for highlighting p-values.}
##'   \item{log2RatioHighlightThresh}{ a numeric specifying the threshold for highlighting log2 ratios.}
##'   \item{maxGenesForClustering}{ an integer specifying the maximum amount of genes for clustering. If the amount is higher, the least significant genes get removed first.}
##'   \item{minGenesForClustering}{ an integer specifying the minimum amount of genes for clustering. If the amount is lower, no clustering will be done.}
##'   \item{doZip}{ a logical indicating whether to archive the result file.}
##'   \item{goseqMethod}{ a character specifying the method for the GO analysis. Accepted values: Wallenius, Sampling or Hypergeometric.}
##' }
##' @template rawData-template
##' @template types-template
##' @template roxygen-template
writeNgsMultiGroupReport = function(dataset, result, htmlFile, param=NA, rawData=NA, types=NULL) {
  seqAnno = rawData$seqAnno
  
  titles = list()
  titles[["Analysis"]] = paste("Analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  
  addDataset(doc, dataset, param)
  addCountResultSummary(doc, param, result)
  
  titles[["Result Summary"]] = "Result Summary"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  settings = character()
  settings["Number of features:"] = length(result$pValue)
  if (!is.null(result$isPresentProbe)){
    settings["Number of features with counts above threshold:"] = sum(result$isPresentProbe)
  }
  addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
  
  titles[["Number of significants by p-value and fold-change"]] = "Number of significants by p-value and fold-change"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addSignificantCounts(doc, result)
  
  resultFile = addResultFile(doc, param, result, rawData)
  
  logSignal = log2(shiftZeros(result$xNorm, param$minSignal))
  result$groupMeans = averageColumns(logSignal, by=param$grouping)
  
  if (param$writeScatterPlots){
    testScatterTitles = addTestScatterPlots(doc, param, logSignal, result, seqAnno, resultFile$resultFile, types)
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
    
    jsFile = system.file("extdata/enrichr.js", package="ezRun", mustWork=TRUE)
    addJavascript(doc, jsFile)
    if (!is.null(clusterResult$GO)){
      goTables = goClusterTable(param, clusterResult, seqAnno)
      # addFlexTable(doc, ezFlexTable(goTables$linkTable, add.rownames=TRUE))
      if (any(c(grepl("Homo_", getOrganism(param$ezRef)), grepl("Mus_", getOrganism(param$ezRef))))){
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
  ## only do GO if we have enough genes
  if (doGo(param, seqAnno)){
    goResult = twoGroupsGO(param, result, seqAnno, normalizedAvgSignal=rowMeans(result$groupMeans), method=param$goseqMethod) ## TODO -> multiGroupsGO()
    titles[["GO Enrichment Analysis"]] = "GO Enrichment Analysis"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    revigoTitle = addGoUpDownResult(doc, param, goResult)
    titles = append(titles, revigoTitle)
    
    ## enrichrLink for mice and humans
    if (any(c(grepl("Homo_", getOrganism(param$ezRef)), grepl("Mus_", getOrganism(param$ezRef))))){
      
      isSig = result$pValue < param$pValThreshGO & result$usedInTest
      isUp = result$log2Ratio > param$log2RatioThreshGO & isSig
      isDown = result$log2Ratio < -param$log2RatioThreshGO & isSig
      regulatedGenes = list()
      regulatedGenes$upGenes = na.omit(unique(seqAnno[isUp, "gene_name"]))
      regulatedGenes$downGenes = na.omit(unique(seqAnno[isDown, "gene_name"]))
      regulatedGenes$bothGenes = union(regulatedGenes$upGenes, regulatedGenes$downGenes)
      
      enrichrLinks = ezMatrix("", rows=c('bothGenes', 'downGenes', 'upGenes'), cols=1)
      for (row in rownames(enrichrLinks)){
        genesToUse = seqAnno$gene_name[which(seqAnno$gene_name %in% regulatedGenes[[row]])]
        genesList = paste(genesToUse, collapse="\\n")
        jsCall = paste0('enrich({list: "', genesList, '", popup: true});')
        enrichrLinks[row, 1] = as.html(pot(paste0("<a href='javascript:void(0)' onClick='", jsCall, "'>Enrichr</a>")))
      }
      titles[["Enrichr"]] = "Enrichr"
      addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
      addFlexTable(doc, ezFlexTable(enrichrLinks, valign="middle", add.rownames = TRUE))
    }
  }
  addParagraph(doc, pot("advanced Plots", hyperlink = "advancedPlots.html"))
  closeBsdocReport(doc, htmlFile, titles)
}
