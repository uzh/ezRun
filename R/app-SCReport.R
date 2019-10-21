###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCReport <-
  setRefClass("EzAppSCReport",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCReport
                  name <<- "EzAppSCReport"
                  appDefaults <<- rbind(scProtocol=ezFrame(Type="character", DefaultValue="10X", Description="Which single cell protocol?"),
                                        minCellsPerGene=ezFrame(Type="numeric", 
                                                                DefaultValue=5, 
                                                                Description="Minimum number of cells per gene for creating Seurat object"),
                                        minGenesPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=500, 
                                                                Description="Minimal number of genes per cell for Seurat filtering"),
                                        maxGenesPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=3000, 
                                                                Description="Maximal number of genes per cell for Seurat filtering"),
                                        maxMitoFraction=ezFrame(Type="numeric", 
                                                                DefaultValue=0.25, 
                                                                Description="Maximal fraction of mitochondrial reads per cell for Seurat filtering"),
                                        minReadsPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=5e4, 
                                                                Description="Minimal reads per cell of smart-Seq2 for Seurat filtering"),
                                        pcs=ezFrame(Type="numeric", 
                                                    DefaultValue=20,
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes=ezFrame(Type="charVector", 
                                                        DefaultValue="", 
                                                        Description="The genes used in supvervised clustering"),
                                        x.low.cutoff=ezFrame(Type="numeric", 
                                                             DefaultValue=0.125,
                                                             Description="Bottom cutoff on x-axis for identifying variable genes"),
                                        x.high.cutoff=ezFrame(Type="numeric", 
                                                              DefaultValue=3,
                                                              Description="Top cutoff on x-axis for identifying variable genes"),
                                        y.cutoff=ezFrame(Type="numeric", 
                                                         DefaultValue=0.5,
                                                         Description="Bottom cutoff on y-axis for identifying variable genes"),
                                        vars.to.regress=ezFrame(Type="charVector", 
                                                                DefaultValue="nUMI,perc_mito", 
                                                                Description="Variables to regress out"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        markersToShow=ezFrame(Type="numeric", 
                                                              DefaultValue=10, 
                                                              Description="The markers to show in the heatmap of cluster marker genes"),
                                        markersToCheck=ezFrame(Type="charList", 
                                                               DefaultValue="", 
                                                               Description="The markers to check"),
                                        runPseudoTime=ezFrame(Type="logical", 
                                                              DefaultValue=FALSE,
                                                              Description="Run PseudoTime for single cell data?"),
                                        all2allMarkers=ezFrame(Type="logical", 
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"))
                }
              )
  )

ezMethodSCReport = function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  ## subset the selected sample names
  # samples <- param$samples
  # input <- input$subset(samples)
  
  sce <- loadSCCountDataset(input, param)
  metadata(sce)$output <- output
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse=", "),
                                    sep=": ")
  sce <- seuratPreProcess(sce)
  
  saveRDS(sce, file="sce.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCReport.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCReport.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}
