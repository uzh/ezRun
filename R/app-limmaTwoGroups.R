###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## for now only gene set testing from limma

EzAppLimma <-
  setRefClass("EzAppLimma",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodLimma
                  name <<- "EzAppLimma"
                  appDefaults <<- rbind(testMethod=ezFrame(Type="character",  DefaultValue="limma",  Description="which test method in limma to use: limma"),
                                        modelMethod=ezFrame(Type="character", DefaultValue="limma-trend", Description="which mean-variance relationship model method in limma to use: limma-trend or voom"),
                                        useRefGroupAsBaseline=ezFrame(Type="logical", DefaultValue=FALSE, Description="should the log-ratios be centered at the reference samples"),
                                        onlyCompGroupsHeatmap=ezFrame(Type="logical", DefaultValue=FALSE, Description="Only show the samples from comparison groups in heatmap"),
                                        priorCount=ezFrame(Type="numeric", DefaultValue=10, Description="prior count to be added to shrink the log-fold-changes")
                  )
                }
              )
  )

ezMethodLimma = function(input=NA, output=NA, param=NA,
                         htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupTwoGroupsInput(input, param)
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1){
    param$grouping2 = input$getColumn(param$grouping2)
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  
  deResult = twoGroupCountComparison(rawData)
  if (isError(deResult)){
    writeErrorReport(htmlFile, param=param, error=deResult$error)
    return("Error")
  }
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "twoGroups.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="twoGroups.Rmd", envir=new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  return("Success")
}
