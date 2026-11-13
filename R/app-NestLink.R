###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch
ezMethodNestLink <- function(input = NA, output = NA, param = NA) {
  require(NestLink)

  if (grepl('>', readLines(param[['knownNBPath']], n = 1))) {
    knownNB_data <- readDNAStringSet(param[['knownNBPath']])
    knownNB <- translate(knownNB_data, if.fuzzy.codon = 'solve')
    param[['knownNB']] <- sapply(knownNB, toString)
  } else {
    knownNB_data <- read.table(
      param[['knownNBPath']],
      sep = '\t',
      header = TRUE,
      row.names = 1,
      stringsAsFactors = FALSE
    )
    knownNB <- translate(
      DNAStringSet(knownNB_data$Sequence),
      if.fuzzy.codon = 'solve'
    )
    names(knownNB) <- rownames(knownNB_data)
    param[['knownNB']] <- sapply(knownNB, toString)
  }

  file <- input$getFullPaths("Read1")
  sampleName <- input$getNames()
  setwdNew(sampleName)
  NB2FC <- runNGSAnalysis(file = file, param)
  write.table(NB2FC, paste0(sampleName, '_NB2FC.txt'), sep = '\t')
  nanobodyFlycodeLinking.as.fasta(
    NB2FC,
    paste0(sampleName, '_NB2FC.fasta'),
    name = sampleName
  )
}

##' @template app-template
##' @templateVar method ezMethodNestLink(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppNestLink <-
  setRefClass(
    "EzAppNestLink",
    contains = "EzApp",
    methods = list(
      ## NestLink, unconditional: runNGSAnalysis() of the NestLink package always runs
      ## (ezMethodNestLink, :33); the package asks to cite this paper. The reads come merged
      ## from FlashApp, so fastp and FLASH are cited there, not here.
      citation = function(param = list()) {
        c(
          "Egloff, P. et al. Engineered peptide barcodes for in-depth analyses of binding protein libraries. Nature Methods 16(5), 421-428 (2019). https://doi.org/10.1038/s41592-019-0389-8"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodNestLink
        name <<- "EzAppNestLink"
      }
    )
  )
