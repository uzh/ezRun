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

## TODO: add runLimma function.
##' @title Compares the counts of two groups
##' @description Compares the counts of two groups with the option to choose from several methods to test them.
##' @template rawData-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{testMethod}{ defines the method to run: deseq2, exactTest, glm, sam or limma. Defaults to glm.}
##'   \item{batch}{ a character vector specifying the batch groups.}
##'   \item{grouping}{ a character specifying the grouping.}
##'   \item{sampleGroup}{ a character specifying the group to sample.}
##'   \item{refGroup}{ a character specifying the reference group.}
##'   \item{normMethod}{ a character specifying the normalization method for the edger and glm test methods.}
##'   \item{runGfold}{ a logical indicating whether to run Gfold.}
##' }
##' @template roxygen-template
##' @return Returns a list containing the results of the comparison.
twoGroupCountComparison = function(rawData, param){
  x = rawData$counts
  presentFlag = rawData$presentFlag
  job = ezJobStart("twoGroupCountComparison")
  result = list()
  result$analysis ="NGS two group analysis"
  if (is.null(param$testMethod)){
    param$testMethod = "glm"    
  }
  
  result$method = param$testMethod
  if (ezIsSpecified(param$batch)){
    if (param$testMethod %in% c("glm", "sam","deseq2")){
      result$method = paste(result$method, "batch")
    } else {
      return(list(error=paste("Second factor only supported for the test methods glm, sam and deseq2")))
    }
  }
  
  isSample = param$grouping == param$sampleGroup
  isRef = param$grouping == param$refGroup
  ## compute which probes are present
  isPresent = ezPresentFlags(x, presentFlag=presentFlag, param=param, isLog=FALSE)
  useProbe = rep(FALSE, nrow(x))
  useProbe[apply(isPresent[, isRef, drop=FALSE], 1, mean) >= 0.5] = TRUE
  useProbe[apply(isPresent[, isSample, drop=FALSE], 1, mean) >= 0.5] = TRUE
  result$isPresentProbe = useProbe
  result$isPresent = isPresent
  res = switch(param$testMethod,
               deseq2 = runDeseq2(round(x), param$sampleGroup, param$refGroup, param$grouping, batch=param$batch, isPresent=useProbe),
               exactTest = runEdger(round(x), param$sampleGroup, param$refGroup, param$grouping, param$normMethod),
               glm = runGlm(round(x), param$sampleGroup, param$refGroup, param$grouping, param$normMethod, batch=param$batch),
               limma = runLimma(x, param$sampleGroup, param$refGroup, param$grouping, param$batch),
               stop("unsupported testMethod: ", param$testMethod)
  )
  result$log2Ratio = res$log2FoldChange  
  result$fitGlm = res$fitGlm
  result$sf = res$sf
  pValue = res$pval
  pValue[is.na(pValue)] = 1
  
  if (!is.null(param$runGfold) && param$runGfold && !is.null(rawData$seqAnno$width) && !is.null(rawData$seqAnno$gene_name)){
    result$gfold = runGfold(rawData, result$sf, isSample, isRef)
  }
  
  useProbe[is.na(useProbe)] = FALSE
  fdr = rep(NA, length(pValue))
  fdr[useProbe] = p.adjust(pValue[useProbe], method="fdr")
  names(pValue) = rownames(x)
  names(fdr) = rownames(x)
  result$pValue = pValue  
  result$fdr=fdr
  result$usedInTest = useProbe
  result$xNorm = ezScaleColumns(x, result$sf)
  ezWriteElapsed(job, status="done")
  return(result)
}

##' @describeIn twoGroupCountComparison Runs the Gfold test method.
runGfold = function(rawData, scalingFactors, isSample, isRef){
  message("running gfold ")
  # prepare input data for gfold
  .writeGfoldInput = function(sampleName){
    gene_name = rawData$seqAnno$gene_name
    if (is.null(gene_name)) gene_name = "NA"
    gfoldData = data.frame(gene_name=gene_name, count=rawData$counts[, sampleName], rawData$seqAnno$width, 
                           rawData$rpkm[, sampleName], row.names=rownames(rawData$seqAnno), check.names=FALSE, stringsAsFactors=FALSE)
    gfoldFile = paste0(sampleName, ".read_cnt")
    ezWrite.table(gfoldData, file=gfoldFile, col.names = FALSE)
    return(gfoldFile)
  }
  sampleFiles = sapply(colnames(rawData$counts)[isSample], .writeGfoldInput)
  refFiles = sapply(colnames(rawData$counts)[isRef], .writeGfoldInput)
  
  # run gfold
  cmd = paste (file.path(GFOLD_DIR, "gfold"), "diff -v 0",
               "-norm", paste(signif(1/c(scalingFactors[isRef], scalingFactors[isSample]), digits=4), collapse=","),
               "-s1", paste(sub(".read_cnt", "", refFiles), collapse=","), 
               "-s2", paste(sub(".read_cnt", "", sampleFiles), collapse=","),
               "-suf .read_cnt -o out.diff")
  ezSystem(cmd)
  gfoldRes = ezRead.table("out.diff", header=FALSE, comment.char="#")
  names(gfoldRes)=c('geneName','gfold','e_FDR', 'log2fdc', '1stRPKM', '2ndRPKM')
  gfold = gfoldRes$gfold
  names(gfold) = rownames(gfoldRes)
  # remove gfold input / output files
  cmd = paste("rm out.diff out.diff.ext", paste(c(sampleFiles, refFiles), collapse=" "))
  ezSystem(cmd)
  return(gfold)
}

##' @describeIn twoGroupCountComparison Runs the Deseq2 test method.
runDeseq2 = function(x, sampleGroup, refGroup, grouping, batch=NULL, isPresent=NULL){
  require("DESeq2", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  if (ezIsSpecified(batch)){
    if (!is.numeric(batch)){
      batch = as.factor(batch)
    } else {
      message("using numeric batch factor")
    }
    colData = data.frame(grouping=as.factor(grouping), batch=batch, row.names=colnames(x))
    dds = DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping + batch)
  } else {
    colData = data.frame(grouping=as.factor(grouping), row.names=colnames(x))
    dds = DESeqDataSetFromMatrix(countData=x, colData=colData, design= ~ grouping)
  }
  dds = estimateSizeFactors(dds, controlGenes=isPresent)
  dds = DESeq(dds, quiet=FALSE)
  res = results(dds, contrast=c("grouping", sampleGroup, refGroup), cooksCutoff=FALSE)
  res=as.list(res)
  res$sf = 1/colData(dds)$sizeFactor
  return(res)
}

##' @describeIn twoGroupCountComparison Runs the Edger test method.
runEdger = function(x, sampleGroup, refGroup, grouping, normMethod){
  library(edgeR, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  cds = DGEList(counts=x, group=grouping)
  cds = calcNormFactors(cds, method=normMethod)
  sf = 1/(cds$samples$norm.factors * cds$samples$lib.size)
  sf = sf / ezGeomean(sf)
  #sf = ezLogmeanScalingFactor(x, presentFlag=x>0)
  if (sum(grouping == refGroup) >=2 & sum(grouping == sampleGroup) >=2){
    cds = estimateCommonDisp(cds)
    cds = estimateTagwiseDisp(cds)
  } else {
    cds$common.dispersion = 0.1
  }
  et <- exactTest(cds, pair=c(refGroup, sampleGroup))
  
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
runGlm = function(x, sampleGroup, refGroup, grouping, normMethod, batch=NULL){
  library(edgeR, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
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
  if (ezIsSpecified(batch)){
    design = model.matrix( ~ groupFactor + factor(batch[isSample|isRef]))
    colnames(design) = c("Intercept", "Grouping", paste("Batch", 1:(ncol(design)-2), sep="_"))
  } else {
    design = model.matrix( ~ groupFactor)
    colnames(design) = c("Intercept", "Grouping")
  }
  
  ## dispersion estimation
  if (sum(grouping == refGroup) >=2 & sum(grouping == sampleGroup) >=2){
    cds = estimateGLMTrendedDisp(cds, design)
    cds = estimateGLMTagwiseDisp(cds, design)
  } else {
    cds$common.dispersion = 0.1
  }
  
  ## Testing for DE genes
  fitGlm = glmFit(cds, design)
  lrt.2vs1 = glmLRT(fitGlm, coef=2)
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
writeNgsTwoGroupReport = function(dataset, result, htmlFile="00index.html", param=NA, rawData=NA, types=NULL) {
  seqAnno = rawData$seqAnno
  
  titles = list()
  titles[["Analysis"]] = paste("Analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]], dataset=dataset)
  
  titles[["Result Summary"]] = "Result Summary"
  addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  addCountResultSummary(doc, param, result)
  
  titles[["Significant Counts"]] = "Significant Counts"
  addTitleWithAnchor(doc, titles[[length(titles)]], 3)
  addSignificantCounts(doc, result)
  
  resultFile = addResultFile(doc, param, result, rawData)
  ezWrite.table(result$sf, file="scalingfactors.txt", head="Name", digits=4)
  
  logSignal = log2(shiftZeros(result$xNorm, param$minSignal))
  result$groupMeans = cbind(rowMeans(logSignal[ , param$grouping == param$sampleGroup, drop=FALSE]),
                            rowMeans(logSignal[ , param$grouping == param$refGroup, drop=FALSE]))
  colnames(result$groupMeans) = c(param$sampleGroup, param$refGroup)
  
  if (param$writeScatterPlots){
    testScatterTitles = addTestScatterPlots(doc, param, logSignal, result, seqAnno, types) ## colorRange was also not used in the old function
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
    clusterLink = ezImageFileLink(plotCmd, file=clusterPng, width=max(800, 400 + 10 * ncol(xCentered)), height=1000)
    
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
    addTitleWithAnchor(doc, titles[[length(titles)]], 2)
    doc = addParagraph(doc, paste("Significance threshold:", param$pValueHighlightThresh))
    if (param$log2RatioHighlightThresh > 0){
      doc = addParagraph(doc, paste("log2 Ratio threshold:", param$log2RatioHighlightThresh))
    }
    doc = addParagraph(doc, paste("Number of significant features:", sum(use)))
    
    if (!is.null(clusterResult$GO)){
      goLink = as.html(ezGrid(cbind("Background color corresponds to the row colors in the heatmap plot.",
                                as.html(goClusterTable(param, clusterResult)))))
      #goLink[[2]] = addGOClusterResult(doc, param, clusterResult)
    } else {
      goLink = as.html(pot("No information available"))
    }
    tbl = ezGrid(cbind("Cluster Plot"=clusterLink, "GO categories of feature clusters"=goLink), header.columns = TRUE)
    doc = addFlexTable(doc, tbl)
  }
  ## only do GO if we have enough genes
  if (doGo(param, seqAnno)){
    goResult = twoGroupsGO(param, result, seqAnno, normalizedAvgSignal=rowMeans(result$groupMeans), method=param$goseqMethod)
    titles[["GO Enrichment Analysis"]] = "GO Enrichment Analysis"
    addTitleWithAnchor(doc, titles[[length(titles)]], 2)
    doc = addGoUpDownResult(doc, param, goResult)
  }
  
  ## Run Gage
  if(param[['GAGEanalysis']] ) {
    gageRes <- runGageAnalysis(result, param=param, output=output, rawData=rawData)
    titles[["GAGE Enrichment Analysis"]] = "GAGE Enrichment Analysis"
    addTitleWithAnchor(doc, titles[[length(titles)]], 3)
    addGageTables(doc, param, gageRes)
  }
  closeBsdocReport(doc, htmlFile, titles)
}
