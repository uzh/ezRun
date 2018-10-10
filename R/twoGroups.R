###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch




##' @title Cleans up input from two group apps
##' @description Cleans up input from two group apps.
##' @template input-template
##' @param param a list of parameters to use or clean up.
##' @template roxygen-template
##' @return Returns an object of the class EzDataset that is the modified \code{input}.
cleanupTwoGroupsInput = function(input, param){
  dataset = input$meta
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset, param), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  }
  inputMod = EzDataset(meta=dataset, dataRoot=param$dataRoot)
  if (!is.null(param$markOutliers) && param$markOutliers){
    stopifnot(!is.null(dataset$Outlier))
    grouping = inputMod$getColumn(param$grouping)
    isOut = dataset$Outlier %in% c("", "NO", '""', "FALSE") == FALSE
    grouping[isOut] = paste(grouping[isOut], "OUTLIER", sep="_")
    inputMod$setColumn(grouping)
  }
  if (!is.null(param$removeOtherGroups) && param$removeOtherGroups){
    grouping = inputMod$getColumn(param$grouping)
    keep = grouping %in% c(param$sampleGroup, param$refGroup)
    inputMod = inputMod$subset(keep)
  }
  return(inputMod)
}

## TODO: add runLimma function.
##' @title Compares the counts of two groups
##' @description Compares the counts of two groups with the option to choose from several methods to test them.
##' @template rawData-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{testMethod}{ defines the method to run: deseq2, exactTest, glm, sam or limma. Defaults to glm.}
##'   \item{grouping}{ a character specifying the grouping.}
##'   \item{grouping2}{ a character vector specifying a secondary grouping.}
##'   \item{sampleGroup}{ a character specifying the group to sample.}
##'   \item{refGroup}{ a character specifying the reference group.}
##'   \item{normMethod}{ a character specifying the normalization method for the edger and glm test methods.}
##'   \item{runGfold}{ a logical indicating whether to run Gfold.}
##' }
##' @template roxygen-template
##' @return Returns a list containing the results of the comparison.
twoGroupCountComparison = function(rawData){
  require(SummarizedExperiment)
  x = assays(rawData)$counts
  param <- metadata(rawData)$param
  presentFlag = assays(rawData)$presentFlag
  job = ezJobStart("twoGroupCountComparison")
  metadata(rawData)$analysis <- "NGS two group analysis"
  if (is.null(param$testMethod)){
    param$testMethod = "glm"
    metadata(rawData)$param <- param
  }
  if(param$testMethod != "glm"){
    ## deTest option is only for glm method.
    param$deTest <- NULL
    metadata(rawData)$param <- param
  }
  
  metadata(rawData)$method <- param$testMethod
  
  if (ezIsSpecified(param$grouping2)){
    if (param$testMethod %in% c("glm", "sam","deseq2")){
      metadata(rawData)$method <- paste(metadata(rawData)$method, 
                                        "using secondary factor")
    } else {
      return(list(error=paste("Second factor only supported for the test methods glm, sam and deseq2")))
    }
  }
  
  isSample = param$grouping == param$sampleGroup
  isRef = param$grouping == param$refGroup
  ## compute which probes are present
  isPresent = ezPresentFlags(x, presentFlag=presentFlag, param=param,
                             isLog=FALSE)
  useProbe = rep(FALSE, nrow(x))
  useProbe[rowMeans(isPresent[, isRef, drop=FALSE]) >= 0.5] = TRUE
  useProbe[rowMeans(isPresent[, isSample, drop=FALSE]) >= 0.5] = TRUE
  rowData(rawData)$isPresentProbe <- useProbe
  assays(rawData)$isPresent <- isPresent
  
  res = switch(param$testMethod,
               deseq2 = runDeseq2(round(x), param$sampleGroup, param$refGroup, 
                                  param$grouping, grouping2=param$grouping2, 
                                  isPresent=useProbe,
                                  cooksCutoff=ezIsSpecified(param$cooksCutoff) && param$cooksCutoff),
               exactTest = runEdger(round(x), param$sampleGroup, param$refGroup,
                                    param$grouping, param$normMethod,
                                    priorCount=param$backgroundExpression),
               glm = runGlm(round(x), param$sampleGroup, param$refGroup, 
                            param$grouping, param$normMethod, 
                            grouping2=param$grouping2,
                            priorCount=param$backgroundExpression,
                            deTest=param$deTest),
               limma = runLimma(x, param$sampleGroup, param$refGroup, 
                                param$grouping, grouping2=param$grouping2),
               stop("unsupported testMethod: ", param$testMethod)
  )
  rowData(rawData)$log2Ratio <- res$log2FoldChange
  metadata(rawData)$fitGlm = res$fitGlm
  colData(rawData)$sf <- res$sf
  pValue = res$pval
  pValue[is.na(pValue)] = 1
  
  if (!is.null(param$runGfold) && param$runGfold && 
      !is.null(rowData(rawData)$featWidth) && !is.null(rowData(rawData)$gene_name)){
    rowData(rawData)$gfold <- runGfold(rawData, colData(rawData)$sf, 
                                       isSample, isRef)
  }
  metadata(rawData)$nativeResult <- res
  useProbe[is.na(useProbe)] = FALSE
  fdr = rep(NA, length(pValue))
  fdr[useProbe] = p.adjust(pValue[useProbe], method="fdr")
  rowData(rawData)$pValue <- pValue
  rowData(rawData)$fdr <- fdr
  ## usedInTest was used before to pre-filter. Not used for now.
  rowData(rawData)$usedInTest = useProbe
  assays(rawData)$xNorm = ezScaleColumns(x, colData(rawData)$sf)
  
  ezWriteElapsed(job, status="done")
  metadata(rawData)$summary = c("Name"=param$name,
                                "Reference Build"=param$refBuild,
                                "Feature Level"=metadata(rawData)$featureLevel,
                                "Normalization"=param$normMethod)
  
  deResult = EzResult(se=rawData)
  return(deResult)
}

##' @describeIn twoGroupCountComparison Runs the Gfold test method.
runGfold = function(rawData, scalingFactors, isSample, isRef){
  message("running gfold ")
  # prepare input data for gfold
  .writeGfoldInput = function(sampleName){
    gene_name = rowData(rawData)$gene_name
    if (is.null(gene_name)) gene_name = "NA"
    gfoldData = data.frame(gene_name=gene_name, 
                           count=assays(rawData)$counts[, sampleName], 
                           rowData(rawData)$featWidth,
                           assays(rawData)$rpkm[, sampleName], 
                           row.names=rownames(assays(rawData)$counts), 
                           check.names=FALSE, stringsAsFactors=FALSE)
    gfoldFile = paste0(sampleName, ".read_cnt")
    ezWrite.table(gfoldData, file=gfoldFile, col.names = FALSE)
    return(gfoldFile)
  }
  sampleFiles = sapply(colnames(assays(rawData)$counts)[isSample], 
                       .writeGfoldInput)
  refFiles = sapply(colnames(assays(rawData)$counts)[isRef], 
                    .writeGfoldInput)
  
  # run gfold
  cmd = paste ("gfold", "diff -v 0",
               "-norm", paste(signif(1/c(scalingFactors[isRef],
                                         scalingFactors[isSample]), digits=4),
                              collapse=","),
               "-s1", paste(sub(".read_cnt", "", refFiles), collapse=","), 
               "-s2", paste(sub(".read_cnt", "", sampleFiles), collapse=","),
               "-suf .read_cnt -o out.diff")
  ezSystem(cmd)
  gfoldRes = ezRead.table("out.diff", header=FALSE, comment.char="#")
  names(gfoldRes)=c('geneName','gfold','e_FDR', 'log2fdc', '1stRPKM', '2ndRPKM')
  gfold = gfoldRes$gfold
  names(gfold) = rownames(gfoldRes)
  # remove gfold input / output files
  file.remove(c("out.diff", "out.diff.ext", sampleFiles, refFiles))
  return(gfold)
}

##' @describeIn twoGroupCountComparison Runs the Deseq2 test method.
runDeseq2 = function(x, sampleGroup, refGroup, grouping, grouping2=NULL, isPresent=NULL, cooksCutoff=FALSE){
  ## get size factors -- grouping2 not needed
  colData = data.frame(grouping=as.factor(grouping), row.names=colnames(x))
  dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping)
  dds = DESeq2::estimateSizeFactors(dds, controlGenes=isPresent)
  sf = 1/dds@colData$sizeFactor

  
  ## remove the samples that do not participate in the comparison
  isSample = grouping == sampleGroup
  isRef = grouping == refGroup
  grouping = grouping[isSample|isRef]
  x = x[ ,isSample|isRef]
  if (ezIsSpecified(grouping2)){
    grouping2 = grouping2[isSample|isRef]
  }    
  
  ## run the analysis
  if (ezIsSpecified(grouping2)){
    colData = data.frame(grouping=as.factor(grouping), grouping2=grouping2, row.names=colnames(x))
    dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping + grouping2)
  } else {
    colData = data.frame(grouping=as.factor(grouping), row.names=colnames(x))
    dds = DESeq2::DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping)
  }
  dds = DESeq2::estimateSizeFactors(dds, controlGenes=isPresent)
  dds = DESeq2::DESeq(dds, quiet=FALSE, minReplicatesForReplace=Inf)
  res = DESeq2::results(dds, contrast=c("grouping", sampleGroup, refGroup), cooksCutoff=cooksCutoff)
  res = as.list(res)
  res$sf = sf
  return(res)
}

##' @describeIn twoGroupCountComparison Runs the EdgeR test method.
runEdger = function(x, sampleGroup, refGroup, grouping, normMethod, priorCount=0.125){
  require("edgeR", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  cds = DGEList(counts=x, group=grouping)
  cds = calcNormFactors(cds, method=normMethod)
  sf = 1/(cds$samples$norm.factors * cds$samples$lib.size)
  sf = sf / ezGeomean(sf)
  #sf = ezLogmeanScalingFactor(x, presentFlag=x>0)
  if (sum(grouping == refGroup) >=2 & sum(grouping == sampleGroup) >=2){
    cds <- estimateDisp(cds)
  } else {
    cds$common.dispersion = 0.1
  }
  et <- exactTest(cds, pair=c(refGroup, sampleGroup), prior.count = priorCount)
  
  res = list()
  res$id = rownames(et$table)
  res$log2FoldChange = et$table$logFC
  res$pval = et$table$PValue
  #res$log2Expr = et$table$logCPM
  res$sf = sf
  ## the following should be double checked
  ## do not return the groupMeans now
  #   res$groupMeans = log2(cbind(apply(cds$counts[ , grouping == sampleGroup, drop=FALSE], 1, mean), 
  #                          apply(cds$counts[ , grouping == refGroup, drop=FALSE], 1, mean)))
  #   colnames(res$groupMeans) = c(sampleGroup, refGroup)
  #   rownames(res$groupMeans) = res$id
  
  return(res)
}

##' @describeIn twoGroupCountComparison Runs the Glm test method.
runGlm = function(x, sampleGroup, refGroup, grouping, normMethod, grouping2=NULL,
                  priorCount=0.125, deTest=c("QL", "LR")){
  require("edgeR", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  ## differential expression test by quasi-likelihood (QL) F-test or 
  ## likelihood ratio test.
  ## QL as default.
  deTest <- match.arg(deTest)
  
  robust = ezIsSpecified(param$robust) && param$robust

  ## get the scaling factors for the entire data set
  cds = DGEList(counts=x, group=grouping)
  cds = calcNormFactors(cds, method=normMethod)
  sf = 1/(cds$samples$norm.factors * cds$samples$lib.size)
  sf = sf / ezGeomean(sf)
  #sf = ezLogmeanScalingFactor(x, presentFlag=x>0)
  
  ## run analysis and especially dispersion estimates only on subset of the data
  isSample = grouping == sampleGroup
  isRef = grouping == refGroup
  grouping = grouping[isSample|isRef]
  x2 = x[ ,isSample|isRef]
  cds = DGEList(counts=x2, group=grouping)
  cds = calcNormFactors(cds, method=normMethod)
  groupFactor = factor(grouping, levels = c(refGroup, sampleGroup))
  if (ezIsSpecified(grouping2)){
    design = model.matrix( ~ groupFactor + grouping2[isSample|isRef])
    colnames(design) = c("Intercept", "Grouping", paste("Grouping2", 1:(ncol(design)-2), sep="_"))
  } else {
    design = model.matrix( ~ groupFactor)
    colnames(design) = c("Intercept", "Grouping")
  }
  
  ## dispersion estimation
  if (sum(grouping == refGroup) >=2 & sum(grouping == sampleGroup) >=2){
    if (robust){
      cds = estimateGLMRobustDisp(cds, design)
    } else {
      cds <- estimateDisp(cds, design)
    }
  } else {
    cds$common.dispersion = 0.1
  }
  
  ## Testing for DE genes
  if(deTest == "QL"){
    ## quasi-likelihood (QL) F-test
    message("Using quasi-likelihood (QL) F-test!")
    fitGlm = glmQLFit(cds, design, prior.count=priorCount, robust = robust)
    lrt.2vs1 = glmQLFTest(fitGlm, coef=2)
  }else{
    ## likelihood ratio test
    message("Using likelihood ratio test!")
    fitGlm = glmFit(cds, design, prior.count=priorCount, robust=robust)
    lrt.2vs1 = glmLRT(fitGlm, coef=2)
  }
  res = list()
  res$id = rownames(lrt.2vs1$table)
  res$log2FoldChange = lrt.2vs1$table$logFC
  res$pval = lrt.2vs1$table$PValue
  res$fitGlm = fitGlm
  #res$log2Expr = lrt.2vs1$table$logCPM
  res$sf = sf
  ## do not return groupMeans now.
  #res$groupMeans = cbind(apply(cds$count[ , grouping == sampleGroup, drop=FALSE], 1, mean), apply(cds$count[ , grouping == refGroup, drop=FALSE], 1, mean))
  #colnames(res$groupMeans) = c(sampleGroup, refGroup)
  #rownames(res$groupMeans) = rownames(cds$count)
  return(res)
}

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
