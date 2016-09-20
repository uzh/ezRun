###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodScater = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  if (is.null(param$minExpressedCells)) {
    param$minExpressedCells <- 4
  }

  countMatrix = EzDataset(file=input$getFullPaths("CountMatrix"), dataRoot=param$dataRoot)
  input = EzDataset(file=input$getFullPaths("CountDataset"), dataRoot=param$dataRoot)
  dataset = input$meta

  setwdNew(basename(output$getColumn("Report")))
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  }
  input$meta = dataset

  titles = list()
  titles[["scater"]] = paste("scater analysis:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  addDataset(doc, dataset, param)

  # Prepare the data for SCE set
  phenoData = new("AnnotatedDataFrame", data = dataset$meta)
  rownames(phenoData) <- row.names(dataset$meta)
  countData = countMatrix$meta
  featureNames = row.names(countData)
  featureData <- new("AnnotatedDataFrame", data = data.frame(Feature = featureNames))
  rownames(featureData) <- featureNames

  sceset <- newSCESet(countData = countData, phenoData = phenoData, featureData = featureData)
  sceset <- calculateQCMetrics(sceset)
  # Remove features that are not expressed
  expressedMask <- rowSums(is_exprs(sceset)) > param$minExpressedCells
  sceset <- sceset[expressedMask, ]

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

