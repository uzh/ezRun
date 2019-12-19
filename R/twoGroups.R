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
  inputMod = EzDataset(meta=dataset, dataRoot=param$dataRoot)
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
    if (param$testMethod %in% c("glm", "limma","deseq2")){
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
  useProbe = logical(nrow(x))
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
                            deTest=param$deTest,
                            robust=ezIsSpecified(param$robust) && param$robust),
               limma = runLimma(x, param$sampleGroup, param$refGroup, 
                                param$grouping, grouping2=param$grouping2,
                                priorCount=param$priorCount,
                                modelMethod=param$modelMethod),
               stop("unsupported testMethod: ", param$testMethod)
  )
  rowData(rawData)$log2Ratio <- res$log2FoldChange
  metadata(rawData)$fitGlm = res$fitGlm
  colData(rawData)$sf <- res$sf
  pValue = res$pval
  pValue[is.na(pValue)] = 1
  
  # if (!is.null(param$runGfold) && param$runGfold && 
  #     !is.null(rowData(rawData)$featWidth) && !is.null(rowData(rawData)$gene_name)){
  #   rowData(rawData)$gfold <- runGfold(rawData, colData(rawData)$sf, 
  #                                      isSample, isRef)
  # }
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
runDeseq2 = function(x, sampleGroup, refGroup, grouping, grouping2=NULL, 
                     isPresent=NULL, cooksCutoff=FALSE){
  require(DESeq2)
  ## get size factors -- grouping2 not needed
  colData = data.frame(grouping=as.factor(grouping), 
                       row.names=colnames(x))
  dds = DESeqDataSetFromMatrix(countData=x, colData=colData, 
                               design= ~ grouping)
  dds = estimateSizeFactors(dds, controlGenes=isPresent)
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
    colData = data.frame(grouping=as.factor(grouping), grouping2=grouping2, 
                         row.names=colnames(x))
    dds = DESeqDataSetFromMatrix(countData=x, colData=colData, 
                                 design= ~ grouping + grouping2)
  } else {
    colData = data.frame(grouping=as.factor(grouping), row.names=colnames(x))
    dds = DESeqDataSetFromMatrix(countData=x, colData=colData,
                                 design= ~ grouping)
  }
  dds = estimateSizeFactors(dds, controlGenes=isPresent)
  dds = DESeq(dds, quiet=FALSE, minReplicatesForReplace=Inf)
  res = results(dds, contrast=c("grouping", sampleGroup, refGroup),
                cooksCutoff=cooksCutoff)
  res = as.list(res)
  res$sf = sf
  return(res)
}

##' @describeIn twoGroupCountComparison Runs the EdgeR test method.
runEdger = function(x, sampleGroup, refGroup, grouping, normMethod, 
                    priorCount=0.125){
  require("edgeR")
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
                  priorCount=0.125, deTest=c("QL", "LR"), robust=FALSE){
  require("edgeR")
  
  ## differential expression test by quasi-likelihood (QL) F-test or 
  ## likelihood ratio test.
  ## QL as default.
  deTest <- match.arg(deTest)

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
    colnames(design) = c("Intercept", "Grouping", 
                         paste("Grouping2", 1:(ncol(design)-2), sep="_"))
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

runLimma = function(x, sampleGroup, refGroup, grouping, grouping2=NULL,
                    priorCount=3, modelMethod=c("limma-trend", "voom")){
  require(limma)
  require(edgeR)
  
  modelMethod <- match.arg(modelMethod)

  cds <- DGEList(counts=x)
  cds <- calcNormFactors(cds)
  
  sf = 1/(cds$samples$norm.factors * cds$samples$lib.size)
  sf = sf / ezGeomean(sf)
  
  ## run analysis and especially dispersion estimates only on subset of the data
  isSample = grouping == sampleGroup
  isRef = grouping == refGroup
  grouping = grouping[isSample|isRef]
  if(ezIsSpecified(grouping2)){
    grouping2 <- grouping2[isSample|isRef]
  }
  
  x2 = x[ ,isSample|isRef]
  cds = DGEList(counts=x2, group=grouping)
  cds = calcNormFactors(cds)
  
  groupFactor = factor(grouping, levels = c(refGroup, sampleGroup))
  design = model.matrix( ~ groupFactor)
  
  if(modelMethod == "limma-trend"){
    logCPM <- cpm(cds, log=TRUE, prior.count=priorCount)
    if(ezIsSpecified(grouping2)){
      corfit <- duplicateCorrelation(logCPM, design, block=grouping2)
      fit <- lmFit(logCPM, design, block=grouping2, correlation=corfit$consensus)
    }else{
      fit <- lmFit(logCPM, design)
    }
    fit <- eBayes(fit, trend=TRUE)
    topDT <- topTable(fit, num=Inf, coef=ncol(design))
  }else if(modelMethod == "voom"){
    v <- voom(cds, design, plot=FALSE)
    if(ezIsSpecified(grouping2)){
      corfit <- duplicateCorrelation(v, design, block=grouping2)
      fit <- lmFit(v, design, block=grouping2, correlation=corfit$consensus)
    }else{
      fit <- lmFit(v, design)
    }
    fit <- eBayes(fit)
    topDT <- topTable(fit, num=Inf, coef=ncol(design))
  }
  res = list()
  res$id = rownames(topDT)
  res$log2FoldChange = topDT$logFC
  res$pval = topDT$P.Value
  res$sf = sf
  return(res)
}
