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
    param <- getMageckReference(param)
    if(identical(param[['ctrlFile']], character(0))){
        system2("mageck", args = c("count", "-l", param[['dictFile']], "--fastq", inputFile, "-n", sampleName))
    } else {
        system2("mageck", args = c("count", "-l", param[['dictFile']], "--control-sgrna", param[['ctrlFile']], "--fastq", inputFile, "-n", sampleName))
    }
}

getMageckReference <- function(param){
    param[['dictFile']] <- list.files(param[['libName']], pattern = 'MAGeCK.csv$', full.names = TRUE)
    param[['ctrlFile']] <- list.files(param[['libName']], pattern = 'MAGeCK_Ctrl.csv$', full.names = TRUE)
    
    lockFile <- file.path(param[['libName']], "lock")
    if(identical(param[['dictFile']], character(0)) & !file.exists(lockFile)){
        ###Create LockFile
        ezWrite(Sys.info(), con = lockFile)
        
        ###FIND csv file
        basicFile <- list.files(param[['libName']], pattern = '.csv$', full.names = TRUE)
        if(length(basicFile) == 1L){
            myRef <- ezRead.table(basicFile, row.names = NULL, sep = ',', header = FALSE)
            colnames(myRef) <- c('TranscriptName', 'Sequence', 'GeneSymbol', 'isControl')
            
            myRef[['ID']] <- paste(myRef[['TranscriptName']],myRef[['Sequence']], sep = '_')
            myRefList <- split(myRef, f = myRef$isControl)
            refFile <- sub('.csv', '_MAGeCK.csv', basicFile)
            ezWrite.table(myRef[,c('ID', 'Sequence', 'GeneSymbol')], refFile, col.names = FALSE, row.names = FALSE, sep = ',')
            param[['dictFile']] <- refFile
            if(length(myRefList) == 2L){
                ctrlFile <- sub('.csv', '_MAGeCK_Ctrl.csv', basicFile)
                ezWrite.table(myRefList[['TRUE']][,c('ID')], ctrlFile, col.names = FALSE, row.names = FALSE)
                param[['ctrlFile']] <- ctrlFile
            } else {
                param[['ctrlFile']] <- list.files(param[['libName']], pattern = 'MAGeCK_Ctrl.csv$', full.names = TRUE)
            }
        } else {
            file.remove(lockFile)
            stop('no or multiple basic reference file(s) available')
        }
        file.remove(lockFile)
    }
    return(param)
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