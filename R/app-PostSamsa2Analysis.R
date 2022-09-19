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
  
  for (i in ORGannotationFiles) {
    ezSystem(paste("cp", i, "org_results/"))
  }
  
  for (x in FUNCannotationFiles) {
    ezSystem(paste("cp", x, "func_results/"))
  }
  
  rename_input1 <- paste0("rename s/",refGroup,"/control_",refGroup,"/g */TP1-*")
  ezSystem(rename_input1)
  rename_input2 <- paste0("rename s/",sampleGroup,"/experimental_",sampleGroup,"/g */TP5-*")
  ezSystem(rename_input2)
  
  samsa2RscriptsToBeExecCmd1 <- paste("bash /usr/local/ngseq/src/samsa2/R_scripts/run_all_Rscripts.sh org_results/")
  ezSystem(samsa2RscriptsToBeExecCmd1)
  samsa2RscriptsToBeExecCmd2 <- paste("bash /usr/local/ngseq/src/samsa2/R_scripts/run_all_Rscripts.sh func_results/")
  ezSystem(samsa2RscriptsToBeExecCmd2)
  
  org_Shannon_Simpson <- readLines("org_results/Shannon_Simpson_diversity.txt")
  func_Shannon_Simpson <- readLines("func_results/Shannon_Simpson_diversity.txt")
  org_div_graph <- readRDS("org_results/diversity_graph.rds")
  func_div_graph <- readRDS("func_results/diversity_graph.rds")
  org_PCA <- readRDS("org_results/PCA_plot.tab.rds")
  func_PCA <- readRDS("func_results/PCA_plot.tab.rds")
  org_heatmap <- readRDS("org_results/DESeq_heatmap.rds")
  func_heatmap <- readRDS("func_results/DESeq_heatmap.rds")
  org_combined <- readRDS("org_results/combined_graph.rds")
  func_combined <- readRDS("func_results/combined_graph.rds")
  org_volcano <- readRDS("org_results/volcano_plot.rds")
  func_volcano <- readRDS("func_results/volcano_plot.rds")
  org_deseq_results <- read.csv("org_results/DESeq_results.tab", sep = "\t")
  func_deseq_results <- read.csv("func_results/DESeq_results.tab", sep = "\t")
  
  setwdNew(basename(output$getColumn("Report")))
  markdownFile <- "PostSamsa2.Rmd"
  ## Copy the style files and templates
  makeRmdReport(param=param, org_Shannon_Simpson = org_Shannon_Simpson,
                func_Shannon_Simpson = func_Shannon_Simpson, org_div_graph = org_div_graph,
                func_div_graph = func_div_graph, org_PCA = org_PCA,
                func_PCA = func_PCA, org_heatmap = org_heatmap,
                func_heatmap = func_heatmap, org_combined = org_combined,
                func_combined = func_combined, org_volcano = org_volcano,
                func_volcano = func_volcano, org_deseq_results = org_deseq_results,
                func_deseq_results = func_deseq_results, output=output, 
                rmdFile = markdownFile, reportTitle = "Samsa2 Report")
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