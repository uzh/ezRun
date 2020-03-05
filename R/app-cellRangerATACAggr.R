###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppCellRangerATACAggr <-
  setRefClass("EzAppCellRangerATACAggr",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCellRangerATACAggr
                  name <<- "EzAppCellRangerATACAggr"
                  appDefaults <<- rbind(normalize=ezFrame(Type="character",
                                                          DefaultValue="depth",
                                                          Description="String specifying how to normalize the input libraries. Valid values: depth (default), signal, or none."))
                }
              )
  )

ezMethodCellRangerATACAggr = function(input=NA, output=NA, param=NA){
  ## dataset mode
  require(tibble)
  require(dplyr)
  require(readr)
  
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  aggr_input <- tibble(library_id=input$getNames(),
                       fragments=file.path(dirname(input$getFullPaths("Report")),
                                           "fragments.tsv.gz"),
                       cells=file.path(dirname(input$getFullPaths("Report")),
                                       "singlecell.csv")

                       )
  if(any(input$columnHasTag("Factor"))){
    aggr_input2 <- as_tibble(input$meta[ ,input$columnHasTag("Factor"), drop=FALSE],
                             rownames="library_id")
    colnames(aggr_input2) <- sub(" \\[.*", "", colnames(aggr_input2))
    aggr_input <- left_join(aggr_input, aggr_input2)
  }
  
  aggr_input_fn <- tempfile(pattern="aggr_input", tmpdir=getwd(),
                            fileext = ".csv")
  write_csv(aggr_input, path=aggr_input_fn)
  
  cellRangerFolder = paste0(param$name, "-cellRanger")
  
  refDir <- getCellRangerATACReference(param)
  message("Using the reference: ", refDir)
  
  cmd <- paste(CELLRANGERATAC, "aggr", paste0("--id=", cellRangerFolder),
               paste0("--csv=", aggr_input_fn),
               paste0("--normalize=", param$normalize),
               paste0("--reference=", refDir))
  
  ezSystem(cmd)
  
  file.rename(file.path(cellRangerFolder, "outs"),  param$name)
  unlink(cellRangerFolder, recursive=TRUE)
  
  return("Success")
}
