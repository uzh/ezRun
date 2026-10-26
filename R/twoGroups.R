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
cleanupTwoGroupsInput <- function(input, param) {
  dataset <- input$meta
  if (param$useFactorsAsSampleName) {
    dataset$Name <- rownames(dataset)
    rownames(dataset) <- addReplicate(apply(
      ezDesignFromDataset(dataset, param),
      1,
      paste,
      collapse = "_"
    ))
  }
  inputMod <- EzDataset(meta = dataset, dataRoot = param$dataRoot)
  if (!is.null(param$removeOtherGroups) && param$removeOtherGroups) {
    grouping <- inputMod$getColumn(param$grouping)
    keep <- grouping %in% c(param$sampleGroup, param$refGroup)
    inputMod <- inputMod$subset(keep)
    ezLogMethods("Samples outside the sample group and the reference group were removed before the analysis.")
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
twoGroupCountComparison <- function(rawData) {
  require(SummarizedExperiment)
  x <- assays(rawData)$counts
  param <- metadata(rawData)$param
  presentFlag <- assays(rawData)$presentFlag
  job <- ezJobStart("twoGroupCountComparison")
  metadata(rawData)$analysis <- "NGS two group analysis"
  if (is.null(param$testMethod)) {
    param$testMethod <- "glm"
    metadata(rawData)$param <- param
  }
  if (param$testMethod != "glm") {
    ## deTest option is only for glm method.
    param$deTest <- NULL
    metadata(rawData)$param <- param
  }

  metadata(rawData)$method <- param$testMethod

  if (ezIsSpecified(param$grouping2)) {
    if (param$testMethod %in% c("glm", "limma", "deseq2")) {
      metadata(rawData)$method <- paste(
        metadata(rawData)$method,
        "using two factors"
      )
    } else {
      return(list(
        error = paste(
          "Second factor only supported for the test methods glm, sam and deseq2"
        )
      ))
    }
  }

  isSample <- param$grouping == param$sampleGroup
  isRef <- param$grouping == param$refGroup
  ## compute which probes are present
  isPresent <- ezPresentFlags(
    x,
    presentFlag = presentFlag,
    param = param,
    isLog = FALSE
  )
  useProbe <- logical(nrow(x))
  useProbe[rowMeans(isPresent[, isRef, drop = FALSE]) >= 0.5] <- TRUE
  useProbe[rowMeans(isPresent[, isSample, drop = FALSE]) >= 0.5] <- TRUE
  ezLogMethods(
    "Genes were considered expressed when their count exceeded ",
    if (isTRUE(param$useSigThresh)) param$sigThresh else "0",
    " in at least half of the samples of either group."
  )
  rowData(rawData)$isPresentProbe <- useProbe
  assays(rawData)$isPresent <- isPresent

  res <- switch(
    param$testMethod,
    deseq2 = runDeseq2(
      round(x),
      param$sampleGroup,
      param$refGroup,
      param$grouping,
      grouping2 = param$grouping2,
      isPresent = useProbe,
      cooksCutoff = ezIsSpecified(param$cooksCutoff) && param$cooksCutoff,
      useLfcShrink = ezIsSpecified(param$useLfcShrink) && param$useLfcShrink
    ),
    exactTest = runEdger(
      round(x),
      param$sampleGroup,
      param$refGroup,
      param$grouping,
      param$normMethod,
      priorCount = param$backgroundExpression
    ),
    glm = runGlm(
      round(x),
      sampleGroup = param$sampleGroup,
      refGroup = param$refGroup,
      grouping = param$grouping,
      sampleGroupBaseline = param$sampleGroupBaseline,
      refGroupBaseline = param$refGroupBaseline,
      normMethod = param$normMethod,
      grouping2 = param$grouping2,
      priorCount = param$backgroundExpression,
      deTest = param$deTest,
      robust = ezIsSpecified(param$robust) && param$robust
    ),
    limma = runLimma(
      x,
      param$sampleGroup,
      param$refGroup,
      param$grouping,
      grouping2 = param$grouping2,
      priorCount = param$priorCount,
      modelMethod = param$modelMethod
    ),
    stop("unsupported testMethod: ", param$testMethod)
  )
  pValue <- res$pval
  pValue[is.na(pValue)] <- 1
  useProbe[is.na(useProbe)] <- FALSE
  fdr <- rep(NA, length(pValue))
  fdr[useProbe] <- p.adjust(pValue[useProbe], method = "fdr")
  ezLogMethods(
    "P-values were adjusted with the Benjamini-Hochberg method over the expressed genes."
  )

  rowData(rawData)$log2Ratio <- res$log2FoldChange
  colData(rawData)$sf <- res$sf
  rowData(rawData)$pValue <- pValue
  rowData(rawData)$fdr <- fdr
  rowData(rawData)$usedInTest <- useProbe
  metadata(rawData)$nativeResult <- res
  assays(rawData)$xNorm <- ezScaleColumns(x, colData(rawData)$sf)

  ezWriteElapsed(job, status = "done")
  metadata(rawData)$summary <- c(
    "Name" = param$name,
    "Reference Build" = param$refBuild,
    "Feature Level" = metadata(rawData)$featureLevel,
    "Normalization" = param$normMethod
  )

  seqAnno <- data.frame(
    rowData(rawData),
    row.names = rownames(rawData),
    check.names = FALSE
  )

  if (doGo(param, seqAnno)) {
    metadata(rawData)$enrichInput <- compileEnrichmentInput(param, rawData)
    metadata(rawData)$enrichResult <- ezEnricher(
      metadata(rawData)$enrichInput,
      param
    )
    metadata(rawData)$enrichResultGSEA <- ezGSEA(
      metadata(rawData)$enrichInput,
      param
    )
  }
  return(rawData)
}


##' @describeIn twoGroupCountComparison Runs the Deseq2 test method.
runDeseq2 <- function(
  x,
  sampleGroup,
  refGroup,
  grouping,
  grouping2 = NULL,
  isPresent = NULL,
  cooksCutoff = FALSE,
  useLfcShrink = FALSE
) {
  require(DESeq2)
  ## get size factors -- grouping2 not needed
  colData <- data.frame(
    grouping = as.factor(grouping),
    row.names = colnames(x)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = x,
    colData = colData,
    design = ~grouping
  )
  dds <- estimateSizeFactors(dds, controlGenes = isPresent)
  sf <- 1 / dds@colData$sizeFactor
  ezLogMethods(
    "Normalized counts were computed with DESeq2 size factors estimated on all samples",
    if (is.null(isPresent)) "" else ", using the expressed genes as control genes",
    "."
  )

  ## remove the samples that do not participate in the comparison
  isSample <- grouping == sampleGroup
  isRef <- grouping == refGroup
  grouping <- grouping[isSample | isRef]
  x <- x[, isSample | isRef]
  if (ezIsSpecified(grouping2)) {
    grouping2 <- grouping2[isSample | isRef]
  }
  ezLogMethods("Only the samples of the sample group and the reference group were used in the DESeq2 fit.")

  ## run the analysis
  if (ezIsSpecified(grouping2)) {
    colData <- data.frame(
      grouping = as.factor(grouping),
      grouping2 = grouping2,
      row.names = colnames(x)
    )
    colData$grouping <- factor(
      x = colData$grouping,
      levels = c(refGroup, sampleGroup)
    )
    dds <- DESeqDataSetFromMatrix(
      countData = x,
      colData = colData,
      design = ~ grouping + grouping2
    )
    ezLogMethods(
      "The DESeq2 model included the second factor as an additive ",
      if (is.numeric(grouping2)) "numeric covariate" else "categorical covariate",
      "."
    )
  } else {
    colData <- data.frame(
      grouping = as.factor(grouping),
      row.names = colnames(x)
    )
    colData$grouping <- factor(
      x = colData$grouping,
      levels = c(refGroup, sampleGroup)
    )
    dds <- DESeqDataSetFromMatrix(
      countData = x,
      colData = colData,
      design = ~grouping
    )
    ezLogMethods("The DESeq2 model had the comparison factor as its only term.")
  }
  dds <- estimateSizeFactors(dds, controlGenes = isPresent)
  dds <- DESeq(dds, quiet = FALSE, minReplicatesForReplace = Inf)
  ezLogMethods(
    "Differential expression was tested with DESeq2 (Wald test) on these samples, with size factors re-estimated on them and without outlier replacement."
  )
  res <- results(
    dds,
    contrast = c("grouping", sampleGroup, refGroup),
    cooksCutoff = cooksCutoff
  )
  ezLogMethods(
    "Log2 fold changes are for the sample group over the reference group; p-values were ",
    if (isTRUE(cooksCutoff)) "filtered" else "not filtered",
    " by Cook's distance."
  )
  if (useLfcShrink){
    message("use lfcShrink")
    resRaw <- res
    res <- lfcShrink(dds, res=res, contrast = c("grouping", sampleGroup, refGroup), 
                             type="ashr")
    ezLogMethods("Log2 fold changes were shrunk with ashr (DESeq2 lfcShrink).")
    res$log2FoldChangeNoShrink <- resRaw$log2FoldChange
    res$lfcSENoShrink <- resRaw$lfcSE
  }
  res$gene_id <- rownames(res)
  res <- as.list(res)
  res$sf <- sf
  res$dds = dds
  return(res)
}

##' @describeIn twoGroupCountComparison Runs the EdgeR test method.
runEdger <- function(
  x,
  sampleGroup,
  refGroup,
  grouping,
  normMethod,
  priorCount = 0.125
) {
  require("edgeR")
  cds <- DGEList(counts = x, group = grouping)
  cds <- calcNormFactors(cds, method = normMethod)
  sf <- 1 / (cds$samples$norm.factors * cds$samples$lib.size)
  sf <- sf / ezGeomean(sf)
  ezLogMethods(
    "edgeR normalization factors were computed with calcNormFactors (method ", normMethod,
    ") on all samples of the dataset and, with the library sizes, gave the normalized counts."
  )
  # sf = ezLogmeanScalingFactor(x, presentFlag=x>0)
  if (sum(grouping == refGroup) >= 2 & sum(grouping == sampleGroup) >= 2) {
    cds <- estimateDisp(cds)
    ezLogMethods("Negative-binomial dispersions (common, trended and tagwise) were estimated with edgeR estimateDisp over all groups of the dataset.")
  } else {
    cds$common.dispersion <- 0.1
    ezLogMethods("Because the sample group or the reference group had fewer than two samples, the negative-binomial dispersion was fixed at 0.1.")
  }
  et <- exactTest(
    cds,
    pair = c(refGroup, sampleGroup),
    prior.count = priorCount
  )
  ezLogMethods(
    "The sample group was compared with the reference group with the edgeR exact test using the ",
    if (!is.null(cds$tagwise.dispersion)) "tagwise dispersions" else "fixed common dispersion",
    ", and log2 fold changes were computed with a prior count of ", priorCount, "."
  )

  res <- list()
  res$id <- rownames(et$table)
  res$log2FoldChange <- et$table$logFC
  res$pval <- et$table$PValue
  # res$log2Expr = et$table$logCPM
  res$sf <- sf
  ## the following should be double checked
  ## do not return the groupMeans now
  #   res$groupMeans = log2(cbind(apply(cds$counts[ , grouping == sampleGroup, drop=FALSE], 1, mean),
  #                          apply(cds$counts[ , grouping == refGroup, drop=FALSE], 1, mean)))
  #   colnames(res$groupMeans) = c(sampleGroup, refGroup)
  #   rownames(res$groupMeans) = res$id

  return(res)
}

##' @describeIn twoGroupCountComparison Runs the Glm test method.
runGlm <- function(
  x,
  sampleGroup,
  refGroup,
  grouping,
  sampleGroupBaseline = NULL,
  refGroupBaseline = NULL,
  normMethod,
  grouping2 = NULL,
  priorCount = 0.125,
  deTest = c("QL", "LR"),
  robust = FALSE
) {
  require(edgeR)
  deTest <- match.arg(deTest)
  if (!ezIsSpecified(sampleGroupBaseline)) {
    sampleGroupBaseline <- NULL
  }
  if (!ezIsSpecified(refGroupBaseline)) {
    refGroupBaseline <- NULL
  }
  ## get the scaling factors for the entire data set
  cds <- DGEList(counts = x, group = grouping)
  cds <- calcNormFactors(cds, method = normMethod)
  sf <- 1 / (cds$samples$norm.factors * cds$samples$lib.size)
  sf <- sf / ezGeomean(sf)
  ezLogMethods(
    "The normalized counts were computed from edgeR calcNormFactors (method ", normMethod,
    ") normalization factors and library sizes over all samples of the dataset."
  )

  ## run analysis and especially dispersion estimates only on subset of the data
  indices2Keep <- grouping %in%
    c(refGroup, refGroupBaseline, sampleGroup, sampleGroupBaseline)
  grouping <- grouping[indices2Keep]
  x2 <- x[, indices2Keep]

  cds <- DGEList(counts = x2, group = grouping)
  cds <- calcNormFactors(cds, method = normMethod)
  ezLogMethods(
    "For the edgeR GLM, only the samples of the ",
    if (!is.null(refGroupBaseline) || !is.null(sampleGroupBaseline)) {
      "sample group, the reference group and the baseline groups"
    } else {
      "sample group and the reference group"
    },
    " were kept, and normalization factors were recomputed on them with calcNormFactors (method ", normMethod, ")."
  )
  groupFactor <- factor(
    grouping,
    levels = c(refGroup, sampleGroup, refGroupBaseline, sampleGroupBaseline)
  )
  if (ezIsSpecified(grouping2)) {
    grouping2 <- grouping2[indices2Keep]
    design <- model.matrix(~ 0 + groupFactor + grouping2)
    ezLogMethods(
      "The edgeR GLM used a no-intercept design with one coefficient per group plus the second factor (grouping2) as an additive term without interaction, ",
      "entered as a ",
      if (is.numeric(grouping2)) "numeric covariate" else "categorical factor",
      "."
    )
  } else {
    design <- model.matrix(~ 0 + groupFactor)
    ezLogMethods("The edgeR GLM used a no-intercept design with one coefficient per group.")
  }
  colnames(design) <- str_replace(colnames(design), "^groupFactor", "")

  if (sum(grouping == refGroup) >= 2 & sum(grouping == sampleGroup) >= 2) {
    if (robust) {
      cds <- estimateGLMRobustDisp(cds, design)
      if (deTest == "LR") {
        ezLogMethods("Negative-binomial dispersions were estimated robustly with edgeR estimateGLMRobustDisp on the design matrix, and the GLM fit used its tagwise dispersions and observation weights.")
      }
    } else {
      cds <- estimateDisp(cds, design)
      if (deTest == "LR") {
        ezLogMethods("Negative-binomial dispersions (common, trended and tagwise) were estimated with edgeR estimateDisp on the design matrix, and the GLM fit used the tagwise dispersions.")
      }
    }
  } else {
    cds$common.dispersion <- 0.1
    if (deTest == "LR") {
      ezLogMethods("Because the sample group or the reference group had fewer than two samples, the negative-binomial dispersion was fixed at 0.1.")
    }
  }

  if (deTest == "QL") {
    ezLog("Using quasi-likelihood (QL) F-test!")
    fitGlm <- glmQLFit(cds, design, prior.count = priorCount, robust = robust)
    ezLogMethods(
      "Features were tested with the edgeR quasi-likelihood F-test",
      if (isTRUE(robust)) " with robust moderation of the dispersions",
      "; log2 fold changes were computed with a prior count of ", priorCount, "."
    )
    if (is.null(refGroupBaseline) || is.null(sampleGroupBaseline)) {
      contrastsIndices <- c(-1, 1, rep(0, ncol(fitGlm) - 2))
    } else {
      contrastsIndices <- c(-1, 1, 1, -1, rep(0, ncol(fitGlm) - 4))
    }
    lrt.2vs1 <- glmQLFTest(fitGlm, contrast = contrastsIndices)
  } else {
    ezLog("Using likelihood ratio test!")
    fitGlm <- glmFit(cds, design, prior.count = priorCount, robust = robust)
    ezLogMethods(
      "Features were tested with the edgeR likelihood-ratio test (glmFit, glmLRT); log2 fold changes were computed with a prior count of ",
      priorCount, "."
    )
    if (is.null(refGroupBaseline) || is.null(sampleGroupBaseline)) {
      contrastsIndices <- c(-1, 1, rep(0, ncol(fitGlm) - 2))
    } else {
      contrastsIndices <- c(-1, 1, 1, -1, rep(0, ncol(fitGlm) - 4))
    }
    lrt.2vs1 <- glmLRT(fitGlm, contrast = contrastsIndices)
  }
  ezLogMethods(
    "Log2 fold changes and p-values are for the contrast ",
    if (is.null(refGroupBaseline) || is.null(sampleGroupBaseline)) {
      "sample group minus reference group"
    } else {
      "(sample group minus its baseline group) minus (reference group minus its baseline group)"
    },
    "."
  )
  res <- list()
  res$id <- rownames(lrt.2vs1$table)
  res$log2FoldChange <- lrt.2vs1$table$logFC
  res$pval <- lrt.2vs1$table$PValue
  res$fitGlm <- fitGlm
  res$sf <- sf
  return(res)
}

runLimma <- function(
  x,
  sampleGroup,
  refGroup,
  grouping,
  grouping2 = NULL,
  priorCount = 3,
  modelMethod = c("limma-trend", "voom")
) {
  require(limma)
  require(edgeR)

  modelMethod <- match.arg(modelMethod)

  cds <- DGEList(counts = x)
  cds <- calcNormFactors(cds)

  sf <- 1 / (cds$samples$norm.factors * cds$samples$lib.size)
  sf <- sf / ezGeomean(sf)

  ## run analysis and especially dispersion estimates only on subset of the data
  isSample <- grouping == sampleGroup
  isRef <- grouping == refGroup
  grouping <- grouping[isSample | isRef]
  if (ezIsSpecified(grouping2)) {
    grouping2 <- grouping2[isSample | isRef]
  }

  x2 <- x[, isSample | isRef]
  cds <- DGEList(counts = x2, group = grouping)
  cds <- calcNormFactors(cds)

  groupFactor <- factor(grouping, levels = c(refGroup, sampleGroup))
  design <- model.matrix(~groupFactor)

  if (modelMethod == "limma-trend") {
    logCPM <- cpm(cds, log = TRUE, prior.count = priorCount)
    if (ezIsSpecified(grouping2)) {
      corfit <- duplicateCorrelation(logCPM, design, block = grouping2)
      fit <- lmFit(
        logCPM,
        design,
        block = grouping2,
        correlation = corfit$consensus
      )
    } else {
      fit <- lmFit(logCPM, design)
    }
    fit <- eBayes(fit, trend = TRUE)
    topDT <- topTable(fit, num = Inf, coef = ncol(design))
  } else if (modelMethod == "voom") {
    v <- voom(cds, design, plot = FALSE)
    if (ezIsSpecified(grouping2)) {
      corfit <- duplicateCorrelation(v, design, block = grouping2)
      fit <- lmFit(v, design, block = grouping2, correlation = corfit$consensus)
    } else {
      fit <- lmFit(v, design)
    }
    fit <- eBayes(fit)
    topDT <- topTable(fit, num = Inf, coef = ncol(design))
  }
  res <- list()
  res$id <- rownames(topDT)
  res$log2FoldChange <- topDT$logFC
  res$pval <- topDT$P.Value
  res$sf <- sf
  return(res)
}
