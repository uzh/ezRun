###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurDataClean = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  require(rmarkdown)
  mothurExe = "/usr/local/ngseq/src/mothur-1.39.5/mothur"
  mothurBatchFile = "/home/grusso/Rcodes/giancarlo/genericScipts/mothurMiSeqSOPtest.batch"
  mothurInput = "datasetNoHeaderForMothur.tsv"
  
  setwdNew(basename(output$getColumn("Static Report")))
  ## remove header to creat input files for mothur and copy it 
  cmdNoHead = paste("grep -v Name", "dataset.tsv", ">", mothurInput)
  ezSystem(cmdNoHead)
  cmdMothur = paste(mothurExe,mothurBatchFile)
  ezSystem(cmdMothur)
  
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
                  appDefaults <<- rbind(cutOff = ezFrame(Type="numeric",  DefaultValue="0,03",Description="Cut-off for OTU clustering.")
                  )
                }
              )
  )
