###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurErrorEstimateAndCluster = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  dataset = input$meta
  ### copy files locally
  copyRefCmd <- paste("cp", param$referenceFasta,"./", sep = " ")
  ezSystem(copyRefCmd)
  copyCountTablePacbioCmd <- paste("cp", input$getFullPaths("CountTablePacBio"),"./", sep = " ")
  ezSystem(copyCountTablePacbioCmd)
  copyClusteredFastaFilePacbioCmd <- paste("cp", input$getFullPaths("PreClusteredFastaFilePacBio"),"./", sep = " ")
  ezSystem(copyCClusteredFastaFilePacbioCmd)
  copyCountTableIllCmd <- paste("cp", input$getFullPaths("CountTableIllumina"),"./", sep = " ")
  ezSystem(copyCountTableIllCmd)
  copyClusteredFastaFileIllCmd <- paste("cp", input$getFullPaths("PreClusteredFastaFileIllumina"),"./", sep = " ")
  ezSystem(copyClusteredFastaFileIllCmd)
  ### update batch file pacbio with parameters and run mothur
  updateBatchCmdPacbio <- paste0("sed -e s/\"REFERENCE\"/", basename(param$referenceFasta), "/g",
                                 " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                                 " -e s/\"Mothur\"/\"PacBio\"/g ",
                                 " -e s/\"GROUPS\"/", param$group, "/g",
                                 " -e s/\"INPUT_COUNT\"/", basename(dataset$CountTablePacbio), "/g",
                                 " -e s/\"INPUT_FASTA\"/", basename(dataset$PreClusteredFastaFilePacbio), "/g ",                                 
                                 MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_TEMPLATE, " > ", MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_PACBIO)
  ezSystem(updateBatchCmdPacbio)
  cmdMothurPacBio = paste(MOTHUR_EXE,MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_PACBIO)
  ezSystem(cmdMothurPacBio)
  ### update batch file Illumina with parameters and run mothur
  updateBatchCmdPacbio <- paste0("sed -e s/\"REFERENCE\"/", basename(param$referenceFasta), "/g",
                                 " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                                 " -e s/\"Mothur\"/\"Illumina\"/g ",
                                 " -e s/\"GROUPS\"/", param$group, "/g",
                                 " -e s/\"INPUT_COUNT\"/", basename(dataset$CountTableIllumina), "/g",
                                 " -e s/\"INPUT_FASTA\"/", basename(dataset$PreClusteredFastaFileIllumina), "/g ",  
                                 MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_TEMPLATE, " >", MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_ILLUMINA)
  ezSystem(updateBatchCmdIllumina)
  cmdMothurIllumina = paste(MOTHUR_EXE,MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_ILLUMINA)
  ezSystem(cmdMothurIllumina)
  
  ## Define input for rmd file
  errorCountFileNamePacbio <- "PacBio.good.unique.good.filter.unique.precluster.pick.pick.error.count"
  stepFilePacbio <- "PacBio.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.steps"
  sharedFilePacbio <- "PacBio.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
  
  errorCountFileNameIllumina <- "Illumina.good.unique.good.filter.unique.precluster.pick.pick.error.count"
  stepFileIllumina <- "Illumina.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.steps"
  sharedFileIllumina <- "Illumina.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurErrorEstimateAndCluster.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="MothurErrorEstimateAndCluster.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}
 
##' @template app-template
##' @templateVar method ezMethodMothurErrorEstimateAndCluster()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurErrorEstimateAndCluster <-
  setRefClass("EzAppMothurErrorEstimateAndCluster",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurErrorEstimateAndCluster
                  name <<- "EzAppMothurErrorEstimateAndCluster"
                  appDefaults <<- rbind(cutOff = ezFrame(Type="numeric",  DefaultValue="0.03",Description="Cut-off for OTU clustering."),
                                        group = ezFrame(Type="character",  DefaultValue="Mock",Description="Mock group."),
                                        referenceFasta = ezFrame(Type="character",  DefaultValue="",Description="Mock reference seqeuences.")
                  )
                }
              )
  )
