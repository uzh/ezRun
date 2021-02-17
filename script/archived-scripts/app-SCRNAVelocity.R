###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCRNAVelocity <-
  setRefClass("EzAppSCRNAVelocity",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCRNAVelocity
                  name <<- "EzAppSCRNAVelocity"
                  appDefaults <<- rbind(markersToCheck=ezFrame(Type="charVector", DefaultValue="", Description="The markers to check"),
                                        scProtocol=ezFrame(Type="character", DefaultValue="", Description="Which single cell protocol?")
                                        )
                }
              )
  )

ezMethodSCRNAVelocity <- function(input=NA, output=NA, param=NA, 
                                  htmlFile="00index.html"){
  require(velocyto.R)
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))  
  on.exit(setwd(cwd), add=TRUE)
  
  reportCwd <- getwd()
  
  param$name <- paste(param$name, input$getNames(), sep=": ")
  
  # pythonPATH <- Sys.getenv("PYTHONPATH")
  # Sys.setenv("PYTHONPATH"="")
  # on.exit(Sys.setenv("PYTHONPATH"=pythonPATH), add=TRUE)
  
  if(toupper(param$scProtocol) == "SMART-SEQ2"){
    bamFn <- list.files(path=input$getFullPaths("ResultDir"),
                        pattern="\\.bam$", full.names = TRUE)
    stopifnot(length(bamFn) == 1L)
    
    bams <- splitBamByRG(bamFn, mc.cores=min(param$cores, 8L))
    # run velocyto
    cmd <- paste("velocyto run_smartseq2 -v", paste(bams, collapse=" "),
                 param$ezRef['refFeatureFile'])
    ezSystem(cmd)
    file.remove(bams)
  }else if(toupper(param$scProtocol) == "10X"){
    cellRangerDir <- input$getFullPaths("ResultDir")
    file.copy(from=cellRangerDir, to=".",
              recursive = TRUE)
    # run velocyto
    cmd <- paste("velocyto run10x --samtools-threads 8 --samtools-memory 512 -v",
                 basename(cellRangerDir), param$ezRef['refFeatureFile'])
    ezSystem(cmd)
    
    file.copy(from=file.path(basename(cellRangerDir), "velocyto"),
              to=".", recursive = TRUE)
    unlink(basename(cellRangerDir), recursive = TRUE)
  }else{
    stop("Unsupported single cell protocol.")
  }
  
  # gene-relative model
  loomFn <- list.files("velocyto", "\\.loom$", full.names = TRUE)
  ldat <- read.loom.matrices(loomFn)
  ldat <- lapply(ldat,function(x) {
    colnames(x) <-  gsub("\\.bam","",gsub(".*:","",colnames(x)))
    x
  })
  
  scResults <- readRDS(file.path(input$getFullPaths("Report"),
                                 basename(input$getColumn("Live Report"))))
  tSNE_data <- scResults$tSNE_data
  cell.dist <- 1-armaCor(t(scResults$scData@dr$pca@cell.embeddings))
  
  ## Remove plate/run name
  tSNE_data$cells <- sub(".*___", "", tSNE_data$cells)
  rownames(cell.dist) <- colnames(cell.dist) <- sub(".*___", "", 
                                                    colnames(cell.dist))
  
  cell.dist <- as.dist(cell.dist)
  if(toupper(param$scProtocol) == "SMART-SEQ2"){
    # Use Pearson linear correlation distance on all genes (log scale) to 
    # find k closest cells for the SMART-seq2 datasets
    cell.dist <- NULL
  }
  
  ## save object for report
  ans <- list(ldat=ldat, tSNE_data=tSNE_data, 
              cell.dist=cell.dist, param=param)
  saveRDS(ans, file="ans.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCRNAVelocity.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCRNAVelocity.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}
