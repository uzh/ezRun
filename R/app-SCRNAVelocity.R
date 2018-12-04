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
                #   appDefaults <<- rbind(minReadsPerCell=ezFrame(Type="numeric", DefaultValue=1500, Description="Filter cells with less reads counted on genes"),
                #                         minReadsPerGene=ezFrame(Type="numeric", DefaultValue=3, Description="Minimal number of reads per gene to be expressed"),
                #                         minGenesPerCell=ezFrame(Type="numeric", DefaultValue=500, Description="Filter cells with less genes expressed"))
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
  
  param$scProtocol <- ifelse("STARLog" %in% input$colNames, "smart-Seq2", "10x")
  
  if(param$scProtocol == "smart-Seq2"){
    bams <- splitBamByRG(input$getFullPath("BAM"),
                         mc.cores=min(param$cores, 8L))
    # run velocyto
    cmd <- paste("velocyto run_smartseq2", paste(bams, collapse=" "),
                 param$ezRef['refFeatureFile'])
    ezSystem(cmd)
    file.remove(bams)
  }else{
    stop("Not implemented yet!")
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
  
}
