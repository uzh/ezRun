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
  
  ezSystem(paste("mkdir orgresults"))
  ezSystem(paste("mkdir funcresults"))
  
  for (i in ORGannotationFiles) {
    ezSystem(paste("cp", i, "orgresults/"))
  }
  
  for (x in FUNCannotationFiles) {
    ezSystem(paste("cp", x, "funcresults/"))
  }
  
  setwdNew(getwd())
  write.csv(dataset, "check.metadata.tsv")
  dataset <- read.csv("check.metadata.tsv")
  for (directory in c("orgresults/", "funcresults/")) {
    print(list.files(directory))
    factor_list <- ""
    for (i in as.list(dataset[dataset[ , grepl( group , names( dataset ) ) ] == refGroup,]$annotationFileRefSeq..File.)) { 
      factor_list <- c(factor_list, unlist(strsplit(i, split='/', fixed=TRUE))[3]) 
      }
    factor_list <- factor_list[-1]
    print(as.list(dataset[dataset[ , grepl( group , names( dataset ) ) ] == refGroup,]$annotationFileRefSeq..File.))
    print(factor_list)
    print(refGroup)
    print(group)
    factor_list2 <- ""
    for (i in factor_list) { 
      factor_list2 <- c(factor_list2, unlist(strsplit(i, split='.', fixed=TRUE))[1]) 
    }
    factor_list2 <- factor_list2[-1]
    print(factor_list2)
    list_of_files_torename <- ""
    for (i in factor_list2) { 
      list_of_files_torename <- c(list_of_files_torename, list.files(path = directory, pattern = i)) 
      }
    list_of_files_torename <- list_of_files_torename[-1]
    print(list_of_files_torename)
    add_prefix <- function(x, path = directory, sep = "_"){ paste(paste0(directory,"control"), x, sep = sep) }
    file.rename(paste0(directory,list_of_files_torename), add_prefix(list_of_files_torename))
    
  }
  
  for (directory in c("orgresults/", "funcresults/")) {
    factor_list <- ""
    for (i in as.list(dataset[dataset[ , grepl( group , names( dataset ) ) ] == sampleGroup,]$annotationFileRefSeq..File.)){ 
      factor_list <- c(factor_list, unlist(strsplit(i, split='/', fixed=TRUE))[3]) 
      }
    factor_list <- factor_list[-1]
    factor_list2 <- ""
    for (i in factor_list) { 
      factor_list2 <- c(factor_list2, unlist(strsplit(i, split='.', fixed=TRUE))[1]) 
      }
    factor_list2 <- factor_list2[-1]
    list_of_files_torename <- ""
    for (i in factor_list2) { 
      list_of_files_torename <- c(list_of_files_torename, list.files(path = directory, pattern = i)) 
      }
    list_of_files_torename <- list_of_files_torename[-1]
    add_prefix <- function(x, path = directory, sep = "_"){ 
      paste(paste0(directory,"experimental"), x, sep = sep) 
      }
    file.rename(paste0(directory,list_of_files_torename), add_prefix(list_of_files_torename))
    
  }
  ezSystem(paste("rename s/_/-/g */*.txt"))
  samsa2RscriptsToBeExecCmd1 <- paste("bash /usr/local/ngseq/src/samsa2/R_scripts/run_all_Rscripts.sh orgresults/")
  ezSystem(samsa2RscriptsToBeExecCmd1)
  samsa2RscriptsToBeExecCmd2 <- paste("bash /usr/local/ngseq/src/samsa2/R_scripts/run_all_Rscripts.sh funcresults/")
  ezSystem(samsa2RscriptsToBeExecCmd2)
  
  org_Shannon_Simpson <- readLines("orgresults/Shannon_Simpson_diversity.txt")
  func_Shannon_Simpson <- readLines("funcresults/Shannon_Simpson_diversity.txt")
  org_div_graph <- readRDS("orgresults/diversity_graph.rds")
  func_div_graph <- readRDS("funcresults/diversity_graph.rds")
  org_PCA <- readRDS("orgresults/PCA_plot.tab.rds")
  func_PCA <- readRDS("funcresults/PCA_plot.tab.rds")
  org_heatmap <- readRDS("orgresults/DESeq_heatmap.rds")
  func_heatmap <- readRDS("funcresults/DESeq_heatmap.rds")
  org_combined <- readRDS("orgresults/combined_graph.rds")
  func_combined <- readRDS("funcresults/combined_graph.rds")
  org_volcano <- readRDS("orgresults/volcano_plot.rds")
  func_volcano <- readRDS("funcresults/volcano_plot.rds")
  org_deseq_results <- read.csv("orgresults/DESeq_results.tab", sep = "\t")
  func_deseq_results <- read.csv("funcresults/DESeq_results.tab", sep = "\t")
  
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