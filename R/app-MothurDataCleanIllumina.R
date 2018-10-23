###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurDataCleanIllumina = function(input=NA, output=NA, param=NA, 
                                   htmlFile="00index.html"){
  
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  library(scales)
  dataset = input$meta
  ### read fastq files and prepare inputs for Mothur
  prepareFilesLocallyForMothur(dataset,param)
  projNum <- dirname(param$resultDir)
  
  ### cp silva reference locally
  cpSilvaRefCmd <- "cp /srv/GT/databases/silva/silva.bacteria.forMothur.fasta silva.bacteria.fasta"
  ezSystem(cpSilvaRefCmd)
  
  ### update batch file Illumina with parameters and run mothur: step 1, identify region
  updateBatchCmdIllumina <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen_Illumina, "/g",
                                   " -e s/\"MAX_LEN\"/", param$maxLen_Illumina, "/g",
                                   " -e s/\"Mothur\"/\"Illumina\"/g ",
                                   MOTHUR_DATA_CLEAN_BATCH_TEMPLATE_STEP1, " >", MOTHUR_DATA_CLEAN_BATCH_ILLUMINA_STEP1)
  ezSystem(updateBatchCmdIllumina)
  cmdMothurIllumina = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_ILLUMINA_STEP1)
  ezSystem(cmdMothurIllumina)
  ### extract region
  summaryOfAlign <- read.delim('Illumina.unique.good.good.summary', stringsAsFactors = FALSE, header = T)
  regionStartCoord <- quantile(summaryOfAlign$start, probs = seq(0, 1, 0.025))["5%"]
  regionEndCoord <- quantile(summaryOfAlign$end, probs = seq(0, 1, 0.025))["95%"]
  
  ### update batch file Illumina with parameters and run mothur: step 2, precluster and remove non-bacterial reads
  updateBatchCmdIllumina <- paste0("sed -e s/\"CUTOFF\"/", param$cutOff, "/g",
                                   " -e s/\"DIFFS\"/", param$diffs_Illumina, "/g",
                                   " -e s/\"START_COORD\"/", regionStartCoord, "/g",
                                   " -e s/\"END_COORD\"/", regionEndCoord, "/g",
                                   " -e s/\"Mothur\"/\"Illumina\"/g ",
                                   MOTHUR_DATA_CLEAN_BATCH_TEMPLATE_STEP2,
                                   " >", MOTHUR_DATA_CLEAN_BATCH_ILLUMINA_STEP2)
  ezSystem(updateBatchCmdIllumina)
  cmdMothurIllumina = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_ILLUMINA_STEP2)
  ezSystem(cmdMothurIllumina)
  
  
  ## Define input for rmd file
  rawIllumina <- ezRead.table("Illumina.summary")
  deduppedIllumina <- ezRead.table("Illumina.unique.summary")
  lengthAndHomopFilteredIllumina <- ezRead.table("Illumina.unique.good.summary")
  mappedFilteredIllumina <- ezRead.table("Illumina.unique.good.good.good.summary")
  chimeraIllumina <- read.delim("Illumina.unique.good.good.good.filter.unique.precluster.denovo.vsearch.chimeras", header = FALSE)
  preClusteredAndChimeraFilteredIllumina <- ezRead.table("Illumina.unique.good.good.good.filter.unique.precluster.pick.summary")
  preClusteredAndChimeraCountIllumina <- ezRead.table("Illumina.unique.good.good.good.filter.unique.precluster.denovo.vsearch.pick.count_table")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurDataCleanIllumina.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="MothurDataCleanIllumina.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}

##' @template app-template
##' @templateVar method ezMethodMothurDataCleanIllumina()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurDataCleanIllumina <-
  setRefClass("EzAppMothurDataCleanIllumina",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurDataCleanIllumina
                  name <<- "EzAppMothurDataCleanIllumina"
                  appDefaults <<- rbind(cutOff = ezFrame(Type="integer",  DefaultValue="80",Description="Cut-off for taxonomy assignment"),
                                        diffs_Illumina = ezFrame(Type="integer",  DefaultValue="2",Description="Differences allowed in the pre.cluster step.Should be 1 every 100 bases"),
                                        minLen_Illumina = ezFrame(Type="integer",  DefaultValue="290",Description="Min length Illumina"),     
                                        maxLen_Illumina = ezFrame(Type="integer",  DefaultValue="330",Description="Max length Illumina")
                  )
                }
              )
  )
