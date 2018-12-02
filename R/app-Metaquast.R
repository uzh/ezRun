###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetaquast = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  sampleName = input$getNames()
    ## copy everything locally 
    ## refs
   copyLoopOverFiles <- function(x){ 
    lapply(x,function(x) ezSystem(paste("cp",x,"./")))
   }
   refListFile <-  param$fileWithListOfRefs
   refList = as.list(read.delim(refListFile))
    lapply(refList,copyLoopOverFiles)
    localrefListFile <- basename(refListFile)
    ezSystem(paste("cp",refListFile,localrefListFile))
    
    ## draft
    draft = input$getFullPaths("contigFile")
    localDraft <- basename(draft)
    ezSystem(paste("cp",draft,basename(draft)))
    
      cmd = paste("quast", "--references-list", localrefListFile, "-o", sampleName, 
                  '-t', ezThreads(), localDraft, "1> ", paste0(sampleName,"_quast.log"))
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
