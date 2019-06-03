###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


EzAppSCScran <-
  setRefClass("EzAppSCScran",
              contains = "EzApp",
              methods = list(
                initialize = function(){
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCScran
                  name <<- "EzAppSCScran"
                  appDefaults <<- rbind(scProtocol=ezFrame(Type="character", DefaultValue="10X", 
                                                           Description="Which single cell protocol? 10X or Smart-seq2."),
                                        snnK=ezFrame(Type="integer", DefaultValue=10, 
                                                     Description="An integer scalar specifying the number of nearest neighbors to consider during graph construction. Larger value, more fined clusters. Used in 10X."),
                                        visMethod=ezFrame(Type="character", DefaultValue="TSNE", 
                                                          Description="Which visualisation method? TSNE, UMAP or DiffusionMap."),
                                        knownMarkers=ezFrame(Type="charList", DefaultValue="", 
                                                             Description="Known markers to plot."),
                                        runPseudoTime=ezFrame(Type="logical", 
                                                              DefaultValue=FALSE,
                                                              Description="Run PseudoTime for single cell data?"))
                }
                )
              )

ezMethodSCScran <- function(input=NA, output=NA, param=NA, 
                            htmlFile="00index.html"){
  require(SummarizedExperiment)
  # SAMPLE mode
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  sce <- loadSCCountDataset(input, param)
  metadata(sce)$output <- output
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse=", "),
                                    sep=": ")
  # debug
  # saveRDS(sce,  "sce.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCScran.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCScran.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)

  return("Success")
}
