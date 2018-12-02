###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetaquast = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  sampleName = input$getNames()
  draft = input$getFullPaths("Draft")
    refList = param$fileWithListOfRefs
      cmd = paste("quast", "--references-list", ref, "-o", sampleName, 
                  '-t', ezThreads(), draft, "1> ", paste0(sampleName,"_quast.log"))
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
