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
  mothurExe = "/usr/local/ngseq/src/mothur-1.39.5/mothur"
  mothurBatchFile = "/home/grusso/Rcodes/giancarlo/genericScipts/mothurSingleEndErrorEstimateAndClusterApp.batch"
  mothurInput = "datasetNoHeaderForMothur.tsv"
  setwdNew(basename(output$getColumn("Static Report")))
  dataset = input$meta
  
  ### update batch file pacbio with parameters and run mothur
  updateBatchCmdPacbio <- paste0("sed -e s/\"REFERENCE\"/", param$mockReference, "/g",
                                 " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                                 " -e s/\"Mothur\"/\"PacBio\"/g ",
                                 " -e s/\"GROUPS\"/", param$group, "/g",
                                 " -e s/\"INPUT_COUNT\"/", input$CountTablePacBio, "/g",
                                 " -e s/\"INPUT_FASTA\"/", input$PreClusteredFastaFilePacBio, "/g",                                 
                                 MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_TEMPLATE, " > ", MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_PACBIO)
  ezSystem(updateBatchCmdPacbio)
  cmdMothurPacBio = paste(MOTHUR_EXE,MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_PACBIO)
  ezSystem(cmdMothurPacBio)
  ### update batch file Illumina with parameters and run mothur
  updateBatchCmdPacbio <- paste0("sed -e s/\"REFERENCE\"/", param$mockReference, "/g",
                                 " -e s/\"CUTOFF\"/", param$cutOff, "/g",
                                 " -e s/\"Mothur\"/\"Illumina\"/g ",
                                 " -e s/\"GROUPS\"/", param$group, "/g",
                                 " -e s/\"INPUT_COUNT\"/", input$CountTableIllumina, "/g",
                                 " -e s/\"INPUT_FASTA\"/", input$PreClusteredFastaFileIllumina, "/g",  
                                 MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_TEMPLATE, " >", MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_ILLUMINA)
  ezSystem(updateBatchCmdIllumina)
  cmdMothurIllumina = paste(MOTHUR_EXE,MOTHUR_ERROR_ESTIMATE_AND_CLUSTER_BATCH_ILLUMINA)
  ezSystem(cmdMothurIllumina)
  

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "FastQC.Rmd", "FastQC_overview.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  
  plots = c("Summary of sequences in each sample.png",
            "Per sequence quality scores"="per_sequence_quality.png",
            "Per tile sequence quality"="per_tile_quality.png",
            "Per base sequence content"="per_base_sequence_content.png",
            "Per sequence GC content"="per_sequence_gc_content.png",
            "Per base N content"="per_base_n_content.png",
            "Sequence Length Distribution"="sequence_length_distribution.png",
            "Sequence Duplication Levels"="duplication_levels.png",
            "Adapter Content"="adapter_content.png",
            "Kmer Content"="kmer_profiles.png")
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
                                        reference = ezFrame(Type="character",  DefaultValue="",Description="Mock reference seqeuences")
                  )
                }
              )
  )
