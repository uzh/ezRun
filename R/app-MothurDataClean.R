###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurDataClean = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
#  mothurExe = "/usr/local/ngseq/src/mothur-1.39.5/mothur"
#  mothurBatchFile = "/home/grusso/Rcodes/giancarlo/genericScipts/mothurSingleEndCleanApp.batch"
  setwdNew(basename(output$getColumn("Static Report")))
  dataset = input$meta
  ### read fastq files and prepare inputs for Mothur
  datasetToMothur(dataset)

### update batch file pacbio with parameters and run mothur
updateBatchCmdPacbio <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen_PacBio, "/g",
                               " -e s/\"MAX_LEN\"/", param$maxLen_PacBio, "/g",
                               " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                               " -e s/\"DIFFS\"/", param$diffs_PacBio, "/g",
                               " -e s/\"Mothur\"/\"PacBio\/g",
                               MOTHUR_DATA_CLEAN_BATCH, " >", MOTHUR_DATA_CLEAN_BATCH_PACBIO)
ezSystem(updateBatchCmdPacbio)
cmdMothurPacBio = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_PACBIO)
ezSystem(cmdMothurPacBio)
### update batch file Illumina with parameters and run mothur
updateBatchCmdPacbio <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen_Illumina, "/g",
                               " -e s/\"MAX_LEN\"/", param$maxLen_Illumina, "/g",
                               " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                               " -e s/\"DIFFS\"/", param$diffs_Illumina, "/g",
                               " -e s/\"Mothur\"/\"Illumina\/g",
                               MOTHUR_DATA_CLEAN_BATCH, " >", MOTHUR_DATA_CLEAN_BATCH_ILLUMINA)
ezSystem(updateBatchCmdIllumina)
cmdMothurPacBio = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_ILLUMINA)
ezSystem(cmdMothurIllumina)

## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurDataClean.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
}
 
##' @template app-template
##' @templateVar method ezMethodMothurDataClean()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurDataClean <-
  setRefClass("EzAppMothurDataClean",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurDataClean
                  name <<- "EzAppMothurDataClean"
                  appDefaults <<- rbind(cutOff = ezFrame(Type="numeric",  DefaultValue="80",Description="Cut-off for taxonomy assignment"),
                                        diffs_Illumina = ezFrame(Type="integer",  DefaultValue="2",Description="Differences allowed in the pre.cluster step. 
                                                                 Should be 1 every 100 bases"),
                                        diffs_PacBio = ezFrame(Type="integer",  DefaultValue="15",Description="Differences allowed in the pre.cluster step. 
                                                               Should be 1 every 100 bases"),
                                        minLen_Illumina = ezFrame(Type="integer",  DefaultValue="290",Description="Min length Illumina"),     
                                        maxLen_Illumina = ezFrame(Type="integer",  DefaultValue="290",Description="Max length Illumina"),
                                        minLen_PacBio = ezFrame(Type="integer",  DefaultValue="290",Description="Min length PacBio"),
                                        maxLen_PacBio = ezFrame(Type="integer",  DefaultValue="290",Description="Max length PacBio"),
                  )
                }
              )
  )
