###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetagenomeAtlas = function(input=NA, output=NA, param=NA, 
                                   htmlFile="00index.html"){
  library(xml2)
  library(ape)
  library(tidyverse)
  
  opt = param$cmdOptions
  sampleName = input$getNames()
  cmdPrepare1 = paste("mkdir sub_folder")
  ezSystem(cmdPrepare1)
  for (i in input$getFullPaths("Read1")){
    ezSystem(paste("cp", i, "sub_folder/"))
  } 
  for (i in input$getFullPaths("Read2")){
    ezSystem(paste("cp", i, "sub_folder/"))
  } 
  cmdInitiate = paste("atlas init", "--db-dir /srv/GT/databases/metagenome_atlas_db/", "sub_folder/")
  ezSystem(cmdInitiate)
  tx  <- readLines("config.yaml")
  tx2  <- gsub(pattern = "final_binner: DAS Tool", replace = "final_binner: maxbin", x = tx)
  tx2  <- gsub(pattern = "- gtdb_tree", replace = "#- gtdb_tree", x = tx2)
  tx2  <- gsub(pattern = "- gtdb_taxonomy", replace = "#- gtdb_taxonomy", x = tx2)
  tx2  <- gsub(pattern = "- gtdb_tree", replace = "#- gtdb_tree", x = tx2)
  tx2  <- gsub(pattern = "- gtdb_taxonomy", replace = "#- gtdb_taxonomy", x = tx2)
  tx2  <- gsub(pattern = "#  - checkm_taxonomy", replace = "- checkm_taxonomy", x = tx2)
  tx2  <- gsub(pattern = "#  - checkm_tree", replace = "- checkm_tree", x = tx2)
  writeLines(tx2, con="config.yaml")
  cmd = paste("atlas run qc --profile cluster --latency-wait 400000")
  ezSystem(cmd)
  
  QCHTML <- "reports/QC_report.html"
  #ASSEMBLYHTML <- "reports/assembly_report.html"
  #BINNINGHTML <- "reports/bin_report_maxbin.html"
  #genomic_bins_table <- read_tsv("reports/genomic_bins_maxbin.tsv")
  #taxonomy_table <- read_tsv("genomes/checkm/taxonomy.tsv")
  #newick <- read_file("genomes/tree/checkm.nwk")
  
  setwdNew(basename(output$getColumn("Report")))
  markdownFile <- "MetagenomeAtlasSummary.Rmd"
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", markdownFile, 
                            "fgcz_header.html", "banner.png"))
  rmarkdown::render(input=markdownFile, envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}
##' @template app-template
##' @templateVar method ezMethodMetagenomeAtlas()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMetagenomeAtlas <-
  setRefClass("EzAppMetagenomeAtlas",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMetagenomeAtlas
                  name <<- "EzAppMetagenomeAtlas"
                }
              )
  )