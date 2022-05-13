###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch
ezMethodSubSampleReads <- function(input, output, param){
    input <- ezMethodSubsampleFastq(input = input, param = param, n = param$nReads)
    system('rename s/-subsample// *.gz')
    return('success')
} 



##' @template app-template
##' @templateVar method EzAppSubSampleReads(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppSubSampleReads <-
    setRefClass("EzAppSubSampleReads",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodSubSampleReads
                        name <<- "EzAppSubSampleReads"
                        appDefaults <<- rbind(nReads = ezFrame(Type="numeric",  DefaultValue=FALSE, Description="max number of reads after subsampling"),
                                              paired = ezFrame(Type="logical",  DefaultValue=TRUE, Description="pairedEnd reads"))
                    }
                )
    )
