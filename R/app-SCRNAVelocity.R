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
  
  if(param$scProtocol == "smart-Seq2"){
    bams <- splitBamByRG(input$getFullPaths("BAM"),
                         mc.cores=min(param$cores, 8L))
    # run velocyto
    cmd <- paste("velocyto run_smartseq2", paste(bams, collapse=" "),
                 param$ezRef['refFeatureFile'])
    ezSystem(cmd)
    file.remove(bams)
  }else if(param$scProtocol == "10X"){
    stop("Not implemented yet!")
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
  
  scResults <- readRDS(file.path(input$getFullPath("Report"),
                                 basename(input$getColumn("Live Report"))))
  cell.colorsPalette <- setNames(gg_color_hue(length(levels(scResults$tSNE_data$cluster))),
                                 levels(scResults$tSNE_data$cluster))
  cell.colors <- cell.colorsPalette[scResults$tSNE_data$cluster]
  names(cell.colors) <- sub(".*___", "", scResults$tSNE_data$cells)
  
  emb <- data.frame(row.names=sub(".*___", "", scResults$tSNE_data$cells),
                    scResults$tSNE_data[ ,c("X", "Y")])
  emb <- as.matrix(emb)
  
  ## save object for report
  ans <- list(ldat=ldat, emb=emb, param=param)
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
