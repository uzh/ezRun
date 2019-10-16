###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCScranIntegration <-
  setRefClass("EzAppSCScranIntegration",
              contains = "EzApp",
              methods = list(
                initialize = function(){
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCScranIntegration
                  name <<- "EzAppSCScranIntegration"
                  appDefaults <<- rbind(batchCorrection=ezFrame(Type="character", 
                                                                DefaultValue="MNN",
                                                                Description="Which batch correction method to use? None or MNN"),
                                        runPseudoTime=ezFrame(Type="logical", 
                                                              DefaultValue=FALSE,
                                                              Description="Run PseudoTime for single cell data?"))
                }
              )
  )

ezMethodSCScranIntegration <- function(input=NA, output=NA, param=NA, 
                                       htmlFile="00index.html"){
  # Dataset mode
  
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  param$name <- paste(param$name, paste(input$getNames(), collapse=", "),
                      sep=": ")
  
  sceURLs <- input$getColumn("Static Report")
  
  saveRDS(param, file = "param.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCScranIntegration.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCScranIntegration.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}
