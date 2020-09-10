###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMegahit = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  ### de novo metagenome assemby with Megahit
  
  library(plyr)
  library(dplyr)
  
  sampleName = input$getNames()
  file1PathInDatset <- input$getFullPaths("Read1")
  fastqName1 <- paste0(sampleName,".R1.fastq")
  cpCmd1 <- paste0("gunzip -c ", file1PathInDatset, "  > ", fastqName1)
  ezSystem(cpCmd1)  
  if(param$paired){
    file2PathInDatset <- input$getFullPaths("Read2")
    fastqName2 <- paste0(sampleName,".R2.fastq")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", fastqName2)
    ezSystem(cpCmd2)
    inputStringAss <- paste("\"-1", fastqName1, "-2", fastqName2,"\"")
  } else {
    inputStringAss <- paste("-s",fastqName1)
  }
  if (param$noMercy){
    noMercyString <- "--no-mercy"
  } else  {
    noMercyString <- ""
  }
  if (param$kmin1pass){
    kmin1passString <- "--kmin-1pass"
  } else  {
    kmin1passString <- ""
  }
  

  ##update template
  updateTemplateScriptCmd <- paste0("sed -e s/\"SAMPLE_NAME\"/", sampleName, "/g",
                                    " -e s/\"KMER_MIN\"/", param$kmerMin, "/g ",
                                    " -e s/\"KMER_MAX\"/", param$kmerMax, "/g ",
                                    " -e s/\"KMER_STEP\"/", param$kmerStep, "/g ",
                                    " -e s/\"MIN_COUNT\"/", param$minCount, "/g ",
                                    " -e s/\"NO_MERCY\"/", noMercyString, "/g ",
                                    " -e s/\"KMIN_1_PASS\"/", kmin1passString, "/g ",
                                    " -e s/\"INPUT_FILE_STRING_ASS\"/", inputStringAss, "/g ",
                                    " -e s/\"NUM_CPU\"/", param$cores, "/g ",
                                    " -e s/\"PAIRED\"/", param$paired, "/g ",
                                    file.path(METAGENOMICS_ROOT,MEGAHIT_TEMPLATE_SCRIPT), " >",
                                    MEGAHIT_TEMPLATE_SCRIPT)

  ezSystem(updateTemplateScriptCmd)
  ## run script
  megahitCmd <- paste("bash",MEGAHIT_TEMPLATE_SCRIPT)
  ezSystem(megahitCmd)
  ### megahit completes
  ## place output files
  #1) contigs
  oldContigFile <- "megahitResults/final.contigs.fa"
  
  newContigFile <- basename(output$getColumn("contigFile"))
  ezSystem(paste("cp",oldContigFile,newContigFile))
}
##' @template app-template
##' @templateVar method ezMethodMegahit()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMegahit <-
  setRefClass("EzAppMegahit",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMegahit
                  name <<- "EzAppMegahit"
                  appDefaults <<- rbind(kmerMin = ezFrame(Type="numeric",
                                                                DefaultValue="31",
                                                                Description="Minimum k-mer for the assembly."),
                                        kmerMax = ezFrame(Type="numeric",
                                                          DefaultValue="101",
                                                          Description="Maximum k-mer for the assembly."),
                                        kmerStep = ezFrame(Type="numeric",
                                                          DefaultValue="10",
                                                          Description="Step k-mer for the assembly."),
                                        minCount = ezFrame(Type="numeric",
                                                          DefaultValue="2",
                                                          Description="(kmin+1)-mer with multiplicity lower than this will be discarded."),
                                        kmin1pass = ezFrame(Type="boolean",
                                                          DefaultValue="false",
                                                          Description="Enabled, makes  memory more efficient for ultra low-depth datasets"),
                                        noMercy = ezFrame(Type="boolean",
                                                                  DefaultValue="false",
                                                                  Description="Wheter or not to add mercy k-mers.")
                                        
                                        
                  )
                }
              )
  )

