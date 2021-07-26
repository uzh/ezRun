###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### ATAC-seq ENCODE pipeline
### https://github.com/kundajelab/atac_dnase_pipelines

EzAppAtacENCODE <-
  setRefClass("EzAppAtacENCODE",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodAtacENCODE
                  name <<- "EzAppAtacENCODE"
                }
              )
  )

ezMethodAtacENCODE <- function(input=NA, output=NA, param=NA, 
                               htmlFile="00index.html"){
  
  ## Environments for this pipeline
  setEnvironments("conda")
  Sys.setenv("_JAVA_OPTIONS"="-Xms256M -Xmx728M -XX:ParallelGCThreads=1")
  
  setwdNew(basename(output$getColumn("Report")))
  dataset = input$meta
  samples = rownames(dataset)
  
  ## This pipeline only supports human and mouse
  species <- unique(input$getColumn("Species"))
  stopifnot(length(species)==1L)
  
  if(grepl("(Homo|human)", species, ignore.case = TRUE)){
    refBuild <- "hg38"
  }else if(grepl("(Mus|mouse)", species, ignore.case = TRUE)){
    refBuild <- "mm10"
  }else{
    stop("Only human and mouse are supported by this pipeline!")
  }
  
  ## We only work on paired-end data for now; although the pipeline supporst SE
  fastqR1Fns <- input$getFullPaths("Read1")
  if(param$paired){
  fastqR2Fns <- input$getFullPaths("Read2")
  cmd <- paste(paste0("-fastq", 1:length(fastqR1Fns), "_1"), fastqR1Fns,
               paste0("-fastq", 1:length(fastqR1Fns), "_2"), fastqR2Fns)
  } else {
    cmd <- paste(paste0("-fastq", 1:length(fastqR1Fns), "_1"), fastqR1Fns)
}
  cmd <- paste0(cmd, collapse=" ")
  cmd <- paste("bds", ATACENCODE, "-species", refBuild, 
               "-auto_detect_adapter -nth", param$cores, cmd)
  
  Sys.setenv("LD_LIBRARY_PATH"="/usr/local/ngseq/lib")
  ## This is needed due to library clash in /usr/local/ngseq/lib and /usr/lib/x86_64-linux-gnu
  
  ## This ATAC ENCODE pipeline can fail randomly. 
  ## We force it to rerun until it succeed.
  attempt <- 1L
  status <- try(ezSystem(cmd))
  while(class(status) == "try-error"){
    message("Attemp: ", attempt)
    attempt <- attempt + 1L
    status <- try(ezSystem(cmd))
    if(attempt >= 10L){
      stop("The pipeline still fails after 10 attempts.")
    }
  }
  return("Success")
}
