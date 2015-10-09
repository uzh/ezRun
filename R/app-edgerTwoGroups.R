###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Edger
##' @seealso \code{\link{EzAppEdger}}
ezMethodEdger = function(input=NA, output=NA, param=NA){
  ngsTwoGroupAnalysis(input, output, param=param)
}

##' @template app-template
##' @templateVar method ezMethodEdger()
##' @seealso \code{\link{ezMethodEdger}}
EzAppEdger <-
  setRefClass("EzAppEdger",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  runMethod <<- ezMethodEdger
                  name <<- "EzAppEdger"
                  appDefaults <<- rbind(testMethod=ezFrame(Type="character",  DefaultValue="glm",  Description="which test method in edgeR to use: glm or exactTest"),
                                        normMethod=ezFrame(Type="character", DefaultValue="TMM", Description="edgeR's norm method: TMM, upperquartile, RLE, or none"))
                }
              )
  )

# DESeq2App = function(input=NA, output=NA, param=NA){
#   param$testMethod = "deseq2"
#   param$normMethod = ""
#   if (!is.null(param$markOutliers) && param$markOutliers){
#     stop("DESeq2 does not support marking outliers because marked outliers would still be used in dispersion estimates")
#   }
#   ngsTwoGroupAnalysis(input, output, param=param)
# }


## NOTEP: all functions below go only into here, plus runDeseq2 (from deseq2TwoGroups). missing: runSam, runLimma
##' @title 1
##' @description 1
##' @param input an object of the class EzDataset.
##' @param output an object of the class EzDataset.
##' @param param a list of parameters:
##' \itemize{
##'   \item{sampleGroup}{}
##'   \item{refGroup}{}
##'   \item{useFactorsAsSampleName}{}
##'   \item{removeOutliers}{}
##'   \item{grouping}{}
##'   \item{batch}{}
##'   \item{markOutliers}{}
##' }
##' @param htmlFile
##' @template roxygen-template
##' @return Returns
##' @seealso \code{\link{twoGroupCountComparison}}
##' @seealso \code{\link{writeNgsTwoGroupReport}}
##' @examples
##' 1
ngsTwoGroupAnalysis = function(input=NA, output=NA, param=NULL, htmlFile="00index.html"){
  
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(basename(output$getColumn("Report")))
  stopifnot(param$sampleGroup != param$refGroup)
  dataset = input$meta
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset, param), 1, paste, collapse="_")) # refactorhelp0
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  } 
  input$meta = dataset
  
  rawData = loadCountDataset(input, param)
  
  if (!is.null(dataset[[param$grouping]])){
    param$grouping = dataset[[param$grouping]]    
  } else {
    if (!is.null(dataset[[paste(param$grouping, "[Factor]")]])){
      param$grouping = dataset[[paste(param$grouping, "[Factor]")]]    
    } else {
      stop("column not found: ", param$grouping)
    }
  }
  if (ezIsSpecified((param$batch)) && length(param$batch) == 1){
    if (is.null(dataset[[param$batch]])){
      stop("column not found: ", param$batch)
    }
    param$batch = dataset[[param$batch]]
  }

  if (!is.null(param$markOutliers) && param$markOutliers){
    stopifnot(!is.null(dataset$Outlier))
    isOut = dataset$Outlier %in% c("", "NO", '""', "FALSE") == FALSE
    param$grouping[isOut] = paste(param$grouping[isOut], "OUTLIER", sep="_")
  }  
  
  result = runNgsTwoGroupAnalysis(dataset, rawData=rawData, param=param)
  return(result)
}

##' @describeIn ngsTwoGroupAnalysis Checks for errors and calls \code{twoGroupCountComparison()}, then \code{writeNgsTwoGroupReport()}.
runNgsTwoGroupAnalysis = function(dataset, htmlFile="00index.html", param=param, rawData=NULL, types=NULL){
  
#   checkResult = checkTwoGroupAnalysisConfig(param)
  #if (isError(checkResult)){
  #	writeErrorHtml(htmlFile, param=param, error=checkResult$error)
  #	return("Error")
  #}
  
  if (isError(rawData)){
    writeErrorHtml(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  
  x = rawData$counts ## TODOP: x unused
  result = twoGroupCountComparison(rawData, param)
  if (isError(result)){
    writeErrorHtml(htmlFile, param=param, error=rawData$error)
    return("Error")
  }  
  result$featureLevel = rawData$featureLevel
  result$countName = rawData$countName

  writeNgsTwoGroupReport(dataset, result, htmlFile, param=param, rawData=rawData, types=types)
}


##' @title 1
##' @description 1
##' @param dataset
##' @param result
##' @param htmlFile
##' @param param a list of parameters:
##' \itemize{
##'   \item{logColorRange}{}
##'   \item{minSignal}{}
##'   \item{grouping}{}
##'   \item{sampleGroup}{}
##'   \item{refGroup}{}
##'   \item{pValueHighlightThresh}{}
##'   \item{log2RatioHighlightThresh}{}
##'   \item{maxGenesForClustering}{}
##'   \item{minGenesForClustering}{}
##'   \item{doZip}{}
##'   \item{goseqMethod}{}
##'   \item{maxNumberGroupsDisplayed}{}
##' }
##' @param rawData
##' @param types
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
writeNgsTwoGroupReport = function(dataset, result, htmlFile, param=NA, rawData=NA, types=NULL) {
  keggOrganism = NA
  seqAnno = rawData$seqAnno
  html = openHtmlReport(htmlFile, param=param, title=paste("Analysis:", param$name),
                        dataset=dataset)  
  on.exit(closeHTML(html))
  
  writeCountResultSummary(html, param, result, rawData$type)
  writeResultCounts(html, param, result)
  resultFile = writeResultFile(html, param, result, rawData)
  ezWrite.table(result$sf, file="scalingfactors.txt", head="Name", digits=4)
  
  cr = c(-param$logColorRange, param$logColorRange)
  logSignal = log2(shiftZeros(result$xNorm, param$minSignal))
  result$groupMeans = cbind(rowMeans(logSignal[ , param$grouping == param$sampleGroup, drop=FALSE]),
                            rowMeans(logSignal[ , param$grouping == param$refGroup, drop=FALSE]))
  colnames(result$groupMeans) = c(param$sampleGroup, param$refGroup)
  writeTestScatterPlots(html, param, logSignal, result, seqAnno, types=types, colorRange=cr)
  
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
    clusterResult = clusterHeatmap(param, xCentered, file=clusterPng, nClusters=6, 
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
    writeGOTables(html, param, goResult)
    goFiles = list.files('.',pattern='enrich.*txt')
    keepCols = c('GO.ID','Pvalue')
    revigoLinks = matrix(nrow=3,ncol=3)
    colnames(revigoLinks) = c('BP','CC','MF')
    rownames(revigoLinks) = c('Both','Down','Up')
    for(j in 1:length(goFiles)){
      goResult = read.table(goFiles[j], sep='\t', stringsAsFactors = F, quote='',
                            comment.char = '', header = T)[,keepCols]
      goResult = goResult[which(goResult$Pvalue < param$pValThreshFisher),]
      if(nrow(goResult) > param$maxNumberGroupsDisplayed) {
        goResult = goResult[1:param$maxNumberGroupsDisplayed,]
      }
      revigoLinks[j] = paste('http://revigo.irb.hr/?inputGoList=',
                             paste(goResult[,'GO.ID'],goResult[,'Pvalue'],
                                   collapse='%0D%0A'),sep='')
    }
    ezWrite(paste("<h3>ReViGO </h3>",sep=''), con=html)
    revigoResult = capture.output(writeTableToHtml(revigoLinks))          
    revigoResult = gsub("<td valign='middle' bgcolor='#ffffff'>","<td valign='middle' bgcolor='#ffffff'><a target='_blank' href='",revigoResult)
    revigoResult = gsub("</td>","' type='text/plain'>Link2ReViGo</a></td>",revigoResult)
    ezWrite(revigoResult, con=html)
  }
  
  ## Run Gage
  if(param[['GAGEanalysis']] ) {
    gageRes <- runGageAnalysis(result, param=param, output=output, rawData=rawData)
    writeGageTables(html, param, gageRes)
  }
  
  ezSessionInfo()
  writeTxtLinksToHtml('sessionInfo.txt',con=html)
  return("Success")
}


##' @title 1
##' @description 1
##' @param rawData
##' @param param a list of parameters:
##' \itemize{
##'   \item{testMethod}{}
##'   \item{batch}{}
##'   \item{grouping}{}
##'   \item{sampleGroup}{}
##'   \item{refGroup}{}
##'   \item{normMethod}{}
##'   \item{runGfold}{}
##' }
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
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
               sam = runSam(round(x), param$sampleGroup, param$refGroup, param$grouping, param$batch), 
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
    gfoldFile = paste(sampleName, ".read_cnt", sep="")
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
