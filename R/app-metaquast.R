###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetaQuast = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  draft = input$getFullPaths("Draft")
  cmd = paste("metaquast.py", "--blast-db /srv/GT/databases/silva/release_138_1/SILVA_138.1_LSURef_NR99_tax_silva.fasta -o", sampleName, '-t', ezThreads(), opt, draft, "1> ", paste0(sampleName,"_metaquast.log"))
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMetaQuast <-
  setRefClass("EzAppMetaQuast",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMetaQuast
                  name <<- "EzAppMetaQuast"
                }
              )
  )
