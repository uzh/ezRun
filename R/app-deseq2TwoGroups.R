###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDeseq2 = function(input=NA, output=NA, param=NA){
  if (ezIsSpecified(param$samples)){
    input = input$subset(param$samples)
  }
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupTwoGroupsInput(input, param)
  param$groupingName <- param$grouping
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1){
    param$grouping2Name <- param$grouping2
    param$grouping2 = input$getColumn(param$grouping2)
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport("00index.html", param=param, error=rawData$error)
    return("Error")
  }
  
  deResult = twoGroupCountComparison(rawData)
  if (isError(deResult)){
    writeErrorReport("00index.html", param=param, error=deResult$error)
    return("Error")
  }
  dds = metadata(deResult)$nativeResult$dds
  dataset <- data.frame(colData(deResult), check.names = FALSE)
  dataset <- dataset[rownames(dataset) %in% rownames(dds@colData), ]
  seqAnno <- data.frame(rowData(deResult),
                        row.names = rownames(deResult),
                        check.names = FALSE)
  
  glimmaResults <- list()
  glimmaResults[["res"]] <- results(
    object = dds, 
    contrast = c("grouping", param$sampleGroup, param$refGroup))
  glimmaResults[["status"]] <- as.numeric(
    glimmaResults[["res"]]$pvalue < param$pValueHighlightThresh & 
    abs(glimmaResults[["res"]]$log2FoldChange) > param$log2RatioHighlightThresh)
  glimmaResults[["counts"]] <- as.data.frame(
    assay(varianceStabilizingTransformation(dds)))
  if (!is.null(param$grouping2Name)) {
    glimmaResults[["counts"]] <- limma::removeBatchEffect(
      glimmaResults[["counts"]], dataset[[param$grouping2Name]]
    )
  }
  glimmaResults[["samples"]] <- colnames(dds)
  glimmaResults[["anno"]] <- seqAnno[, c("gene_id", "gene_name", "description")]
  glimmaResults[["anno"]]$GeneID <- glimmaResults[["anno"]]$gene_id
  glimmaResults[["groups"]] <- colData(dds)$grouping
  glimmaResults[["resDF"]] <- as.data.frame(glimmaResults[["res"]])
  glimmaResults[["resDF"]] <- glimmaResults[["resDF"]][
    rownames(glimmaResults[["resDF"]]) %in% 
      seqAnno$gene_id[seqAnno$usedInTest == TRUE], ]
  glimmaResults[["countsXY"]] <- glimmaResults[["counts"]][
    rownames(glimmaResults[["counts"]]) %in% 
      rownames(glimmaResults[["resDF"]]), ]
  glimmaResults[["statusXY"]] <- as.numeric(
    glimmaResults[["resDF"]]$pvalue < 
      param$pValueHighlightThresh & 
      abs(glimmaResults[["resDF"]]$log2FoldChange) > 
      param$log2RatioHighlightThresh)
  glimmaResults[["annoXY"]] <- glimmaResults[["anno"]][
    glimmaResults[["anno"]]$gene_id %in% 
      rownames(glimmaResults[["resDF"]]), ]
  
  Glimma::glMDPlot(
    x = glimmaResults[["res"]],
    counts = glimmaResults[["counts"]], 
    anno = glimmaResults[["anno"]], 
    groups = glimmaResults[["groups"]], 
    samples = glimmaResults[["samples"]],
    status = glimmaResults[["status"]], 
    main = param$comparison, 
    html = paste0(param$comparison, "_MA"), 
    side.main = "gene_name",
    launch = FALSE)
  
  Glimma::glXYPlot(
    x = glimmaResults[["resDF"]]$log2FoldChange,
    y = -log10(glimmaResults[["resDF"]]$pvalue), 
    xlab = "Log2 Fold Change", 
    ylab = "-log10 p-value",
    counts = glimmaResults[["countsXY"]],
    groups = glimmaResults[["groups"]],
    samples = glimmaResults[["samples"]], 
    status = glimmaResults[["statusXY"]],
    anno = glimmaResults[["annoXY"]],
    main = param$comparison, 
    html = paste0(param$comparison, "_Volcano"), 
    launch = FALSE)
  
  Glimma::glMDSPlot(
    x = glimmaResults[["counts"]], 
    top = 2000, 
    labels = glimmaResults[["samples"]],
    groups = glimmaResults[["groups"]],
    main = param$comparison,
    html = paste0(param$comparison, "_MDS"), 
    launch = FALSE)
  
  system(paste0("zip -r ", param$comparison, "_glimma-plots glimma-plots"))
  
  makeRmdReport(output=output, param=param, deResult=deResult, rmdFile="twoGroups.Rmd", reportTitle = param$comparison)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDeseq2(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppDeseq2 <-
  setRefClass("EzAppDeseq2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDeseq2
                  name <<- "EzAppDeseq2"
                  appDefaults <<- rbind(testMethod=ezFrame(Type="character",  DefaultValue="deseq2",  Description="which test method in DESeq to use: deseq2"),
                                        normMethod=ezFrame(Type="character", DefaultValue="DESeq2_MedianRatio", Description="Deseq2's default norm method; this is actually not read"),
                                        useRefGroupAsBaseline=ezFrame(Type="logical", DefaultValue=FALSE, Description="should the log-ratios be centered at the reference samples"),
                                        onlyCompGroupsHeatmap=ezFrame(Type="logical", DefaultValue=FALSE, Description="Only show the samples from comparison groups in heatmap")
                                        )
                }
              )
  )
