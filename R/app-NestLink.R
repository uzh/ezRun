###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch
ezMethodNestLink <- function(input=NA, output=NA, param=NA){
    require(NestLink)
    knownNB_data <- read.table(param[['knownNBPath']], sep='\t', header = TRUE, 
                               row.names = 1, stringsAsFactors = FALSE)
    knownNB <- Biostrings::translate(DNAStringSet(knownNB_data$Sequence))
    names(knownNB) <- rownames(knownNB_data)
    param[['knownNB']] <- sapply(knownNB, toString)
    
    file <- input$getFullPaths("Read1")
    sampleName <- input$getNames()
    setwdNew(sampleName)
    NB2FC <- runNGSAnalysis(file = file, param)
    write.table(NB2FC, paste0(sampleName, '_NB2FC.txt'), sep = '\t')
    nanobodyFlycodeLinking.as.fasta(NB2FC, paste0(sampleName, '_NB2FC.fasta'), name = sampleName)
}

##' @template app-template
##' @templateVar method ezMethodNestLink(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodTrim}}
EzAppNestLink <-
    setRefClass("EzAppNestLink",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodNestLink
                        name <<- "EzAppNestLink"
                    }
                )
    )