###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetaquast = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  library(Biostrings)
  outFileName = param$Name
    ## copy everything locally 
    ## refs
   copyLoopOverFiles <- function(x){ 
    sapply(x,function(x) ezSystem(paste("cp",x,"./")))
   }
   isThereRef <- param$isThereRef
   refFile <-  param$fileWithListOfRefs
   if (isThereRef) {
    listOfRef <- readDNAStringSet(refFile)
    listOfRefNames <- vector()
    for (fileName in names(listOfRef)) {
      writeXStringSet(listOfRef[fileName],paste0(fileName,".fasta"))
      listOfRefNames[fileName] <- paste0(fileName,".fasta")
    }
    listOfRefToParse <- paste(listOfRefNames, collapse = ",")
   }
    ## draft
    draftList = input$getFullPaths("contigFile")
    localDraft <- vector()
    for (k in 1:length(draftList)){
      localDraft[k] <- basename(draftList[k])
      ezSystem(paste("cp",draftList[k],localDraft[k]))
    }
    localDraftCollapsed <- paste(localDraft,collapse = " ")
    sampleNameList <- localDraftCollapsed
    if (isThereRef) {
      cmd = paste("metaquast.py", "-R", listOfRefToParse, "-o", outFileName, 
                  '-t', ezThreads(), sampleNameList, "1> ", "metaQuast.log")
    }else{
      cmd = paste("metaquast.py", "-o", outFileName, 
                  '-t', ezThreads(), sampleNameList, "1> ", "metaQuast.log")
    }
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMetaquast()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMetaquast <-
  setRefClass("EzAppMetaquast",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMetaquast
                  name <<- "EzAppMetaquast"
                  appDefaults <<- rbind(fileWithListOfRefs = ezFrame(Type="character",DefaultValue="", Description="full path to a 
                                                            comma-sep list of refs or a directory with multiple dasta files")
                  )
                                    

                }
              )
  )
