ezMethodMageckTest = function(input=NA, output=NA, param=NA){
  require(Herper)
  
  inputTsv <- data.table::fread(input)
  mergedCountFileName <- paste0(output$Name, ".merged.count.tsv")
  mergedCountFileLoc <- file.path(param$outputDir, mergedCountFileName)
  print(mergedCountFileLoc)
  sampleNames <- inputTsv$Name
  
  countFilePaths <- file.path(param$dataRoot, inputTsv$`Count [File]`)
  print(countFilePaths)
  counts <- lapply(countFilePaths, data.table::fread)
  mergeCounts <- counts %>% reduce(inner_join, by="sgRNA", keep=FALSE)
  mergeCounts <- mergeCounts %>%
    rename("Gene"="Gene.x") %>%
    select(!starts_with("Gene."))
  colnames(mergeCounts)[3:(3+length(sampleNames)-1)] <- sampleNames
  
  data.table::fwrite(mergeCounts, mergedCountFileLoc, sep="\t", na="NA")
  
  rId <- paste(which(inputTsv$`Condition [Factor]` == param$refGroup) - 1, collapse=",")
  sId <- paste(which(inputTsv$`Condition [Factor]` == param$sampleGroup) - 1, collapse=",")
  
  local_CondaEnv("mageckenv", pathToMiniConda = "/usr/local/ngseq/miniconda3")

  opt <- paste(
    "mageck",
    "test",
    "-k",
    mergedCountFileLoc,
    "-t",
    sId,
    "-c",
    rId,
    "-n",
    "foo"
  )
  
  ezSystem(opt)
  
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

source("./R/app-mageck-vars.R")
EzAppMageckTest$new()$run(input=input, output=output, param=param)
#ezMethodMageckTest(input=input, output=output, param=param)