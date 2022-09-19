###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodPostSamsa2Analysis = function(input=NA, output=NA, param=NA, 
                                      htmlFile="00index.html"){
  
  library(DESeq2)
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(data.table)
  library(genefilter)
  library(ggplot2)
  library(gridExtra)
  library(knitr)
  library(optparse)
  library(plyr)
  library(reshape2)
  library(scales)
  library(vegan)
  require(rmarkdown)
  
  dataset = input$meta
  sampleNames = input$getNames()
  isGroupThere = param$grouping != ""
  group = param$grouping
  sampleGroup = param$sampleGroup
  refGroup = param$refGroup
  ORGannotationFiles <- input$getFullPaths("annotationORGFileRefSeq")
  FUNCannotationFiles <- input$getFullPaths("annotationFUNCFileRefSeq")
  
  ezSystem(paste("mkdir org_results"))
  ezSystem(paste("mkdir func_results"))
  ezSystem(paste("cp", ORGannotationFiles, "org_results/"))
  ezSystem(paste("cp", FUNCannotationFiles, "func_results/"))
  
  rename_input1 <- paste0("rename s/",refGroup,"/control_",refGroup,"/g */TP1-*")
  ezSystem(rename_input1)
  rename_input2 <- paste0("rename s/",sampleGroup,"/experimental_",sampleGroup,"/g */TP5-*")
  ezSystem(rename_input2)
  
  samsa2RscriptsToBeExecCmd1 <- paste("bash /usr/local/ngseq/src/samsa2/R_scripts/run_all_Rscripts.sh org_results/")
  ezSystem(samsa2RscriptsToBeExecCmd1)
  samsa2RscriptsToBeExecCmd2 <- paste("bash /usr/local/ngseq/src/samsa2/R_scripts/run_all_Rscripts.sh func_results/")
  ezSystem(samsa2RscriptsToBeExecCmd2)
  
  setwdNew(basename(output$getColumn("Report")))
  markdownFile <- "PostSamsa2.Rmd"
  ## Copy the style files and templates
  makeRmdReport(param=param, output=output, rmdFile = markdownFile, reportTitle = "Samsa2 Report")
}
##' @template app-template
##' @templateVar method ezMethodPostSamsa2Analysis()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppPostSamsa2Analysis<-
  setRefClass("EzAppPostSamsa2Analysis",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodPostSamsa2Analysis
                  name <<- "EzAppPostSamsa2Analysis"
                }
              )
  )