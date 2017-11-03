###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodScater = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  require(scater)
  require(grid)
  if (is.null(param$minExpressedCells)) {
    param$minExpressedCells <- 4
  }
  
  ## input refers to a meta input
  dsFile = input$getFullPaths("CountDataset")
  cmFile = input$getFullPaths("CountMatrix")
  phases = NULL
  if ("CellCyclePhase" %in% input$colNames && any(nchar(input$getColumn("CellCyclePhase")) > 0)) {
	  phaseFile = input$getFullPaths("CellCyclePhase")
	  phases = read.delim(phaseFile, header = T, stringsAsFactors = F)
  }
  

  ## we replace the input with the actual input
  input = EzDataset(file=dsFile, dataRoot=param$dataRoot)
  meta = input$meta
  setwdNew(basename(output$getColumn("Report")))
  if (param$useFactorsAsSampleName){
    meta$Name = rownames(meta)
    rownames(meta) = addReplicate(apply(ezDesignFromDataset(meta), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(meta$Outlier)){
    meta = meta[toupper(meta$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  }
  if (!is.null(phases)) {
    meta$Phase = phases$Phase
  }
  
  input$meta = meta

  titles = list()
  titles[["scater"]] = paste("scater analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  addDataset(doc, meta, param)

  # Prepare the data for SCE set
  phenoData = new("AnnotatedDataFrame", data = meta)
  rownames(phenoData) <- row.names(meta)
  countData = ezRead.table(cmFile)
  featureNames = row.names(countData)
  featureData <- new("AnnotatedDataFrame", data = data.frame(Feature = featureNames))
  rownames(featureData) <- featureNames

  sceset <- newSCESet(countData = as.matrix(countData), phenoData = phenoData, featureData = featureData)
  # When countData is provided, scater ignores is_exprsData parameter value in the constructor. This
  # is probably because it wants to rely on lowerDetectionLimit, which is not properly implemented.
  # Setting lowerDetectionLimit does not generate is_exprs attrbute regardless whether 
  # lowerDetectionLimit is set in the constructor or afterwards. Therefore, is_exprs should be set
  # after sceset object is created.
  is_exprs(sceset) = countData > param$sigThresh
  sceset <- calculateQCMetrics(sceset)
  # Remove features that are not expressed
  expressedMask <- rowSums(is_exprs(sceset)) > param$minExpressedCells
  sceset <- sceset[expressedMask, ]
  
  
  if (!is.null(output)){
    liveReportLink = output$getColumn("Live Report")
    summary = c("Name"=param$name,
                "Reference Build"=param$refBuild,
                "Feature Level"=param$featureLevel) ## the feature level should be in the count dataset but is not!!
    result = EzResult(param=param, sceset=sceset, result=list(summary=summary, analysis="scater"))
    result$saveToFile(basename(output$getColumn("Live Report")))
    addParagraph(doc, ezLink(liveReportLink,
                             "Live Report and Visualizations",
                             target = "_blank"))
  }  
  

  plotCmd <- expression(grid.draw(plot(sceset, exprs_values = "counts")))
  cumPropLink <- ezImageFileLink(plotCmd, file = "libCumProp.png",
                                 mouseOverText = "Cumulative proportion of counts for features with highest expression")
  plotCmd <- expression(grid.draw(plotQC(sceset, type = "highest-expression")))
  topExprLink <- ezImageFileLink(plotCmd, file = "topExpressed.png",
                                 mouseOverText = "Proportion of counts in the top 50 most-expressed features")
  plotCmd <- expression(grid.draw(plotQC(sceset, type = "exprs-freq-vs-mean")))
  exprFreqLink <- ezImageFileLink(plotCmd, file = "exprsFreq.png",
                                  mouseOverText = "Expression frequency vs mean expression level")
  plotCmd <- expression(grid.draw(plotQC(sceset, type = "find-pcs")))
  pcsLink <- ezImageFileLink(plotCmd, file = "pcs.png",
                             mouseOverText = "Top principal components")
  addFlexTable(doc, ezGrid(matrix(cbind(cumPropLink, topExprLink, exprFreqLink, pcsLink), ncol = 1)))

  closeBsdocReport(doc=doc, file=htmlFile, titles)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodScater(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
EzAppScater <-
  setRefClass("EzAppScater",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScater
                  name <<- "EzAppScater"
                  appDefaults <<- rbind(runGO=ezFrame(Type="logical", DefaultValue=TRUE, Description="whether to run the GO analysis"))
                }
              )
  )

EzAppScaterMtx <-
  setRefClass("EzAppScaterMtx",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScaterMtx
                  name <<- "EzAppScaterMtx"
                  appDefaults <<- rbind(runGO=ezFrame(Type="logical",
                                                      DefaultValue=TRUE,
                                                      Description="whether to run the GO analysis"))
                }
              )
  )

ezMethodScaterMtx <- function(input=NA, output=NA, param=NA,
                              htmlFile="00index.html"){
  require(scater)
  require(grid)
  if (is.null(param$minExpressedCells)) {
    param$minExpressedCells <- 4
  }
  
  ## input refers to a meta input
  dsFile = input$getFullPaths("CountDataset")
  cmFile = input$getFullPaths("CountMatrix")
  phases = NULL
  if ("CellCyclePhase" %in% input$colNames && any(nchar(input$getColumn("CellCyclePhase")) > 0)) {
    phaseFile = input$getFullPaths("CellCyclePhase")
    phases = read.delim(phaseFile, header = T, stringsAsFactors = F)
  }
  
  
  ## we replace the input with the actual input
  input = EzDataset(file=dsFile, dataRoot=param$dataRoot)
  meta = input$meta
  setwdNew(basename(output$getColumn("Report")))
  if (param$useFactorsAsSampleName){
    meta$Name = rownames(meta)
    rownames(meta) = addReplicate(apply(ezDesignFromDataset(meta), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(meta$Outlier)){
    meta = meta[toupper(meta$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  }
  if (!is.null(phases)) {
    meta$Phase = phases$Phase
  }
  
  input$meta = meta
  
  titles = list()
  titles[["scater"]] = paste("scater analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  addDataset(doc, meta, param)
  
  # Prepare the data for SCE set
  phenoData = new("AnnotatedDataFrame", data = meta)
  rownames(phenoData) <- row.names(meta)
  countData = ezRead.table(cmFile)
  featureNames = row.names(countData)
  featureData <- new("AnnotatedDataFrame", data = data.frame(Feature = featureNames))
  rownames(featureData) <- featureNames
  
  sceset <- newSCESet(countData = as.matrix(countData), phenoData = phenoData, featureData = featureData)
  # When countData is provided, scater ignores is_exprsData parameter value in the constructor. This
  # is probably because it wants to rely on lowerDetectionLimit, which is not properly implemented.
  # Setting lowerDetectionLimit does not generate is_exprs attrbute regardless whether 
  # lowerDetectionLimit is set in the constructor or afterwards. Therefore, is_exprs should be set
  # after sceset object is created.
  is_exprs(sceset) = countData > param$sigThresh
  sceset <- calculateQCMetrics(sceset)
  # Remove features that are not expressed
  expressedMask <- rowSums(is_exprs(sceset)) > param$minExpressedCells
  sceset <- sceset[expressedMask, ]
  
  
  if (!is.null(output)){
    liveReportLink = output$getColumn("Live Report")
    summary = c("Name"=param$name,
                "Reference Build"=param$refBuild,
                "Feature Level"=param$featureLevel) ## the feature level should be in the count dataset but is not!!
    result = EzResult(param=param, sceset=sceset, result=list(summary=summary, analysis="scater"))
    result$saveToFile(basename(output$getColumn("Live Report")))
    addParagraph(doc, ezLink(liveReportLink,
                             "Live Report and Visualizations",
                             target = "_blank"))
  }  
  
  
  plotCmd <- expression(grid.draw(plot(sceset, exprs_values = "counts")))
  cumPropLink <- ezImageFileLink(plotCmd, file = "libCumProp.png",
                                 mouseOverText = "Cumulative proportion of counts for features with highest expression")
  plotCmd <- expression(grid.draw(plotQC(sceset, type = "highest-expression")))
  topExprLink <- ezImageFileLink(plotCmd, file = "topExpressed.png",
                                 mouseOverText = "Proportion of counts in the top 50 most-expressed features")
  plotCmd <- expression(grid.draw(plotQC(sceset, type = "exprs-freq-vs-mean")))
  exprFreqLink <- ezImageFileLink(plotCmd, file = "exprsFreq.png",
                                  mouseOverText = "Expression frequency vs mean expression level")
  plotCmd <- expression(grid.draw(plotQC(sceset, type = "find-pcs")))
  pcsLink <- ezImageFileLink(plotCmd, file = "pcs.png",
                             mouseOverText = "Top principal components")
  addFlexTable(doc, ezGrid(matrix(cbind(cumPropLink, topExprLink, exprFreqLink, pcsLink), ncol = 1)))
  
  closeBsdocReport(doc=doc, file=htmlFile, titles)
  return("Success")
}