###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCSeurat <-
  setRefClass("EzAppSCSeurat",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCSeurat
                  name <<- "EzAppSCSeurat"
                  appDefaults <<- rbind(scProtocol=ezFrame(Type="character", DefaultValue="10x", Description="Which single cell protocol?"),
                                        minCellsPerGene=ezFrame(Type="numeric", 
                                                                DefaultValue=3, 
                                                                Description="Minimum number of cells per gene for creating Seurat object"),
                                        minGenesPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=1000,
                                                                Description="Minimal number of genes per cell for Seurat filtering"),
                                        maxGenesPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=7000, 
                                                                Description="Maximal number of genes per cell for Seurat filtering"),
                                        maxMitoPercent=ezFrame(Type="numeric", 
                                                                DefaultValue=25, 
                                                                Description="Maximal fraction of mitochondrial reads per cell for Seurat filtering"),
                                        pcs=ezFrame(Type="numeric",
                                                    DefaultValue=50,
                                                    Description="The maximal dimensions to use for reduction"),
                                        vars.to.regress=ezFrame(Type="charVector", 
                                                                DefaultValue="", # nFeature_RNA,percent.mt,CellCycle
                                                                Description="Variables to regress out"),
                                        resolution=ezFrame(Type="numeric",
                                                           DefaultValue=0.5,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        markersToShow=ezFrame(Type="numeric", 
                                                              DefaultValue=10, 
                                                              Description="The markers to show in the heatmap of cluster marker genes"),
                                        knownMarkers=ezFrame(Type="charList", 
                                                               DefaultValue="", 
                                                               Description="The markers to check"),
                                        runPseudoTime=ezFrame(Type="logical", 
                                                              DefaultValue=FALSE,
                                                              Description="Run PseudoTime for single cell data?"))
                }
              )
  )

ezMethodSCSeurat = function(input=NA, output=NA, param=NA, 
                            htmlFile="00index.html"){
  require(Seurat)
  require(scDblFinder)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  sce <- loadSCCountDataset(input, param)
  metadata(sce)$input <- input
  metadata(sce)$output <- output
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse=", "),
                                    sep=": ")
  
  # Doublet detection
  pdf(NULL) ## scDblFinder plots a figure automatically.
  sce <- scDblFinder(sce)
  dev.off()
  
  # saveRDS(sce, file = "sce.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCSeurat.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCSeurat.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}



