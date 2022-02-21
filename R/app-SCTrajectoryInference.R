###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCTrajectoryInference <-
  setRefClass("EzAppSCTrajectoryInference",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCTrajectoryInference
                  name <<- "EzAppSCTrajectoryInference"
                  appDefaults <<- rbind(start_id=ezFrame(Type="charVector", 
                                                         DefaultValue="", 
                                                         Description="Start cluster(s)"),
                                        end_id=ezFrame(Type="character", 
                                                         DefaultValue='', 
                                                         Description="End cluster(s)"),
                                        start_n=ezFrame(Type="numeric", 
                                                         DefaultValue="1", 
                                                         Description="The number of start states"),
                                        end_n=ezFrame(Type="numeric", 
                                                        DefaultValue="1", 
                                                        Description="The number of end states"),
                                        TI_method=ezFrame(Type="charVector", 
                                                          DefaultValue="none", 
                                                          Description="Trajectory inference method(s)"),
                                        show_genes=ezFrame(Type="character", 
                                                               DefaultValue='none', 
                                                               Description="Genes to show along the trajectory"),
                                        root_expression=ezFrame(Type="character", 
                                                         DefaultValue='none', 
                                                         Description="Genes that are highly expressed at the start of the trajectory"),
                                        diff_Branch=ezFrame(Type="character", 
                                                            DefaultValue='none', 
                                                            Description="Method and branch name to extract dysregulated genes from. (For example: Slingshot,3)"),
                                        diff_Branch_Point=ezFrame(Type="character", 
                                                                  DefaultValue='none', 
                                                                  Description="Method and branching point name to extract dysregulated genes from. (For example: Slingshot,3)"))
                } 
                  )
              )
ezMethodSCTrajectoryInference <- function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  library(HDF5Array)
  library(Seurat)
  objectURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(objectURLs)), "sce_h5")
  filePath_course <- file.path(paste0("/srv/GT/analysis/course_sushi/public/projects/", dirname(objectURLs)), "/sce_h5")
  
  if(file.exists(filePath)) {
    object <- loadHDF5SummarizedExperiment(filePath)
    #if it is an rds object it has been likely generated from old reports, so we need to update the seurat version before using the clustering functions below.                                             )
  } else if (file.exists(filePath_course)) {
    object <- loadHDF5SummarizedExperiment(filePath_course)
    } else {
    filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(objectURLs)), "scData.rds")
    scData <- readRDS(filePath)
  }
  
  if (is(object, "SingleCellExperiment")) {
    counts = counts(object)
    logcounts = logcounts(object)
    cells_meta = colData(object)
  } else {
    counts = GetAssayData(object, "counts")
    logcounts = GetAssayData(object, "data")
    cells_meta = object@meta.data
  }
  
  
  library(dyno)
  #library(dyncli)
 
  #Prepare the data
  dyno_dataset <- wrap_expression(expression = t(as.matrix(logcounts)), counts = t(as.matrix(counts)))
  
  #Selecting the best 2 methods predicted by dyno in case no method is specified by the user
  if(param$TI_method=="none") {
     guidelines <- dynguidelines::guidelines(dyno_dataset)
     TI_method <- guidelines$methods_selected[1:2]
  } else {
    TI_method <- param$TI_method
  }
  #Add priors
  if(param$start_id == "none") 
    param$start_id <- levels(object$ident)[1]
  start_cells <- rownames(cells_meta[object$ident %in% param$start_id,])
  end_cells <- rownames(cells_meta[object$ident %in% param$end_id,])
  dyno_dataset <- dyno_dataset %>% add_prior_information(start_id = start_cells, 
                                                         end_id = end_cells,
                                                         start_n = param$start_n,
                                                         end_n = param$end_n)
  
  #Running the best 2 methods selected by dyno or the ones specified by the user
  if(!identical(end_cells, character(0))) 
     priors <- c("start_id", "end_id", "start_n", "end_n")
  else
    priors <- c("start_id", "start_n", "end_n")
  model <- infer_trajectories(dyno_dataset, TI_method, 
                            give_priors = priors,
                            seed=38, verbose = TRUE)
  
  #Model the trajectory using dimensionality reduction
  for (i in 1:length(TI_method)) 
     model[i,]$model[[1]] <- model[i,]$model[[1]]  %>% add_dimred(dyndimred::dimred_mds, expression_source = dyno_dataset$expression)
  
  #save the dyno dataset and the trajectories
  saveRDS(dyno_dataset, "dyno_dataset.rds")
  saveRDS(model, "model.rds")
  
## Copy the style files and templates
styleFiles <- file.path(system.file("templates", package="ezRun"),
                        c("fgcz.css", "SCTrajectoryInference.Rmd",
                          "fgcz_header.html", "banner.png"))
file.copy(from=styleFiles, to=".", overwrite=TRUE)
rmarkdown::render(input="SCTrajectoryInference.Rmd", envir = new.env(),
                  output_dir=".", output_file=htmlFile, quiet=TRUE)

return("Success")
  
}
  