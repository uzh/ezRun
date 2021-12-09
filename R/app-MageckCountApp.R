###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMageckCount <- function(input, output, param){
    require(Herper)
    local_CondaEnv("mageckenv", pathToMiniConda = "/usr/local/ngseq/miniconda3")
    sampleName <- input$getNames()
    inputFile  <- input$getFullPaths("Read1")
    param[['dictPath']] <- list.files(param[['libName']], pattern = 'MAGeCK.csv$', full.names = TRUE)
    ##TODO: handle empty dictPath
    system2("mageck", args = c("count", "-l", param[['dictPath']], "--fastq", inputFile, "-n", sampleName))
}

##' @template app-template
##' @templateVar method ezMethodMageckCount(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppMageckCount <-
    setRefClass("EzAppMageckCount",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodMageckCount
                        name <<- "EzAppMageckCount"
                        appDefaults <<- rbind(
                            libName = ezFrame(
                                Type = "character",
                                DefaultValue = "",
                                Description = "sgRNA Library Name"
                            )
                        )
                    }
                )
    )