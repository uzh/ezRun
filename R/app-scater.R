###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodScater = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  
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
  
  ## load the single cell data
  
  ## create the scater graphics here
  
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

