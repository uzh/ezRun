###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpaceRangerAggr <-
  setRefClass("EzAppSpaceRangerAggr",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSpaceRangerAggr
        name <<- "EzAppSpaceRangerAggr"
        appDefaults <<- rbind(normalize = ezFrame(
          Type = "character",
          DefaultValue = "mapped",
          Description = "String specifying how to normalize depth across the input libraries. mapped or none."
        ))
      }
    )
  )

ezMethodSpaceRangerAggr <- function(input = NA, output = NA, param = NA) {
  ## dataset mode

  aggr_input <- tibble(
    library_id = input$getNames(),
    molecule_h5 = file.path(
      dirname(input$getFullPaths("Report")), "molecule_info.h5"
    ),
    cloupe_file = file.path(
        dirname(input$getFullPaths("Report")), "cloupe.cloupe"
    ),
    spatial_folder = file.path(
        dirname(input$getFullPaths("Report")), "spatial"
    )
  )
  if (any(input$columnHasTag("Factor"))) {
    aggr_input2 <- as_tibble(input$meta[, input$columnHasTag("Factor"), drop = FALSE],
      rownames = "library_id"
    )
    colnames(aggr_input2) <- str_replace(colnames(aggr_input2), " \\[.*", "")
    for (i in 1:ncol(aggr_input2)){
        if(any(aggr_input2[,i] == ''))
            aggr_input2[[colnames(aggr_input2)[i]]] = NULL
    }
    
    aggr_input <- left_join(aggr_input, aggr_input2)
  }

  aggr_input_fn <- tempfile(
    pattern = "aggr_input", tmpdir = getwd(), fileext = ".csv"
  )
  write_csv(aggr_input, file = aggr_input_fn)

  spaceRangerFolder <- str_c(param$name, "-spaceRanger")
  cmd <- str_c(
    "spaceranger aggr", str_c("--id=", spaceRangerFolder),
    str_c("--csv=", aggr_input_fn),
    str_c("--normalize=", param$normalize), sep=" "
  )
  ezSystem(cmd)

  file.rename(file.path(spaceRangerFolder, "outs"), param$name)
  unlink(spaceRangerFolder, recursive = TRUE)

  return("Success")
}
