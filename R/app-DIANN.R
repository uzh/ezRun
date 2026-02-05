###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDIANN = function(input=NA, output=NA, param=NA,
                           htmlFile="00index.html"){
  setwdNew(basename(output$getColumn("Result")))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodCountQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppDIANN <-
  setRefClass("EzAppDIANN",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCountQC
                  name <<- "EzAppCountQC"
                  appDefaults <<- rbind(runGO=ezFrame(Type="logical", DefaultValue=TRUE, Description="whether to run the GO analysis"),
                                        nSampleClusters=ezFrame(Type="numeric", DefaultValue=6, Description="Number of SampleClusters, default value 6"),
                                        selectByFtest=ezFrame(Type="logical", DefaultValue=FALSE, Description="select topGenes by Test instead of SD"),
                                        topGeneSize=ezFrame(Type="numeric", DefaultValue=100, Description="number of genes to consider in gene clustering, mds etc"))
                }
              )
  )
