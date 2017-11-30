###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurDataClean = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
 # cwd <- getwd()
#  setwdNew(basename(output$getColumn("Report")))
#  on.exit(setwd(cwd))
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  dataset = input$meta
  ### read fastq files and prepare inputs for Mothur
  datasetToMothur(dataset,param)

### update batch file pacbio with parameters and run mothur
updateBatchCmdPacbio <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen_PacBio, "/g",
                               " -e s/\"MAX_LEN\"/", param$maxLen_PacBio, "/g",
                               " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                               " -e s/\"DIFFS\"/", param$diffs_PacBio, "/g",
                               " -e s/\"Mothur\"/\"PacBio\"/g ",
                               " -e s/\"pcr.seqs\"/\"#pcr.seqs\"/g ",
                               " -e s/\"rename.file\"/\"#rename.file\"/g ",
                               " -e s/\"start=1968, end=11540\"/\"start=1046, end=43116\"/g ",
                               " -e s/\"v4\"/\"bacteria\"/g ",
                               MOTHUR_DATA_CLEAN_BATCH_TEMPLATE, " > ", MOTHUR_DATA_CLEAN_BATCH_PACBIO)
ezSystem(updateBatchCmdPacbio)
cmdMothurPacBio = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_PACBIO)
ezSystem(cmdMothurPacBio)
### update batch file Illumina with parameters and run mothur
updateBatchCmdIllumina <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen_Illumina, "/g",
                               " -e s/\"MAX_LEN\"/", param$maxLen_Illumina, "/g",
                               " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                               " -e s/\"DIFFS\"/", param$diffs_Illumina, "/g",
                               " -e s/\"Mothur\"/\"Illumina\"/g ",
                               MOTHUR_DATA_CLEAN_BATCH_TEMPLATE, " >", MOTHUR_DATA_CLEAN_BATCH_ILLUMINA)
ezSystem(updateBatchCmdIllumina)
cmdMothurIllumina = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_ILLUMINA)
ezSystem(cmdMothurIllumina)

## Define input for rmd file
rawPacbio <- ezRead.table("PacBio.summary")
lengthDeduppedPacbio <- ezRead.table("PacBio.good.unique.summary")

mappedAndHomopFilteredPacbio <- ezRead.table("PacBio.good.unique.good.summary")
preClusteredAndChimeraFilteredPacbio <- ezRead.table("PacBio.good.unique.good.filter.unique.precluster.pick.summary")

rawIllumina <- ezRead.table("Illumina.summary")
lengthDeduppedIllumina <- ezRead.table("Illumina.good.unique.summary")
mappedAndHomopFilteredIllumina <- ezRead.table("Illumina.good.unique.good.summary")
preClusteredAndChimeraFilteredIllumina <- ezRead.table("Illumina.good.unique.good.filter.unique.precluster.pick.summary")

## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurDataClean.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="MothurDataClean.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
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
                  appDefaults <<- rbind(cutOff = ezFrame(Type="integer",  DefaultValue="80",Description="Cut-off for taxonomy assignment"),
                                        diffs_Illumina = ezFrame(Type="integer",  DefaultValue="2",Description="Differences allowed in the pre.cluster step.Should be 1 every 100 bases"),
                                        diffs_PacBio = ezFrame(Type="integer",  DefaultValue="15",Description="Differences allowed in the pre.cluster step.Should be 1 every 100 bases"),
                                        minLen_Illumina = ezFrame(Type="integer",  DefaultValue="290",Description="Min length Illumina"),     
                                        maxLen_Illumina = ezFrame(Type="integer",  DefaultValue="330",Description="Max length Illumina"),
                                        minLen_PacBio = ezFrame(Type="integer",  DefaultValue="1400",Description="Min length PacBio"),
                                        maxLen_PacBio = ezFrame(Type="integer",  DefaultValue="1700",Description="Max length PacBio")
                  )
                }
              )
  )
