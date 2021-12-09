ezMethodMageckTest = function(input=NA, output=NA, param=NA){
  require(Herper)

  # Loading the variables
  dataset <- input$meta
  dir.create(param$outputDir, showWarnings=FALSE)
  outputPrefix <- file.path(param$outputDir, output$getNames())
  
  mergedCountFileName <- paste0(output$getNames(), ".merged.count.tsv")
  mergedCountFileLoc <- file.path(param$outputDir, mergedCountFileName)
  sampleNames <- rownames(dataset)
  
  # Combining the count files into a single count file
  countFilePaths <- file.path(param$dataRoot, dataset$`Count [File]`)
  counts <- lapply(countFilePaths, data.table::fread)
  mergeCounts <- counts %>% reduce(inner_join, by="sgRNA", keep=FALSE)
  mergeCounts <- mergeCounts %>%
    rename("Gene"="Gene.x") %>%
    select(!starts_with("Gene."))
  colnames(mergeCounts)[3:(3+length(sampleNames)-1)] <- sampleNames
  
  ezWrite.table(mergeCounts, file=mergedCountFileLoc, row.names=FALSE)
  
  # We give the design of the experiment as indices corresponding to the
  # columns (skipping the first 2 positions) which are sample vs ref groups
  rId <- paste(which(dataset$`Condition [Factor]` == param$refGroup) - 1, collapse=",")
  sId <- paste(which(dataset$`Condition [Factor]` == param$sampleGroup) - 1, collapse=",")
  
  local_CondaEnv("mageckenv", pathToMiniConda = "/usr/local/ngseq/miniconda3")

  opt <- c(
    "test",
    "-k",
    mergedCountFileLoc,
    "-t",
    sId,
    "-c",
    rId,
    "-n",
    outputPrefix
  )
  
  system2("mageck", args=opt)
  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMageckTest(input=NA, output=NA, param=NA)
##' @description Use this reference class to run Mageck Test
##' @author Falko Noé
EzAppMageckTest <-
  setRefClass("EzAppMageckTest",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMageckTest
                  name <<- "EzAppMageckTest"
                  appDefaults <<- rbind(
                    outputDir = ezFrame(
                      Type = "character",
                      DefaultValue = ".",
                      Description = "Output directory"
                    )
                  )
                }
              )
  )
