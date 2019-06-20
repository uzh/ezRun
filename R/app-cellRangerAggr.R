###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppCellRangerAggr <-
  setRefClass("EzAppCellRangerAggr",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCellRangerAggr
                  name <<- "EzAppCellRangerAggr"
                  appDefaults <<- rbind(normalize=ezFrame(Type="character",
                                                          DefaultValue="mapped",
                                                          Description="String specifying how to normalize depth across the input libraries. mapped or none."))
                }
              )
  )

ezMethodCellRangerAggr = function(input=NA, output=NA, param=NA){
  ## dataset mode
  require(tibble)
  require(dplyr)
  require(readr)
  
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  # if(!any(input$columnHasTag("Factor"))){
  #   meta <- input$meta
  #   meta$`Condition [Factor]` <- input$getNames()
  #   input <- EzDataset$new(meta=meta, dataRoot=input$dataRoot)
  # }
  
  aggr_input <- tibble(library_id=input$getNames(),
                       molecule_h5=file.path(input$getFullPaths("ResultDir"),
                                             "molecule_info.h5")
                       )
  if(any(input$columnHasTag("Factor"))){
    aggr_input2 <- as_tibble(input$meta[ ,input$columnHasTag("Factor"), drop=FALSE],
                             rownames="library_id")
    colnames(aggr_input2) <- sub(" \\[.*", "", colnames(aggr_input2))
    aggr_input <- left_join(aggr_input, aggr_input2)
  }
  
  aggr_input_fn <- tempfile(pattern="aggr_input", fileext = ".csv")
  write_csv(aggr_input, path=aggr_input_fn)
  
  cellRangerFolder = paste0(param$name, "-cellRanger")
  
  cmd <- paste(CELLRANGER, "aggr", paste0("--id=", cellRangerFolder),
               paste0("--csv=", aggr_input_fn),
               paste0("--normalize=", param$normalize))
  ezSystem(cmd)
  ezSystem(paste("mv ", file.path(cellRangerFolder, "outs"),  param$name))
  
  return("Success")
}
