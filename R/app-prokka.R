###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodProkka = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  draft = input$getFullPaths("Draft")
  cmd = paste("prokka --outdir", sampleName, "--prefix", sampleName, "--locustag", sampleName, "--compliant --centre FGCZ", "--kingdom", param$kingdom, "--cpus", ezThreads(), opt, draft, "1>", paste0(sampleName,"_prokka.log"))
  ezSystem(cmd)
  wddir <- "."
  outfile <- paste0("sampleName", ".faa")
  outfile <- file.path(wddir, sampleName, outfile)
  ezSystem(paste("cp", outfile, basename(output$getColumn("Proteins"))))
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodProkka
##' @description Use this reference class to run
##' @seealso \code{\link{getPbmm2Reference}}
EzAppProkka <-
  setRefClass("EzAppProkka",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodProkka
        name <<- "EzAppProkka"
        appDefaults <<- rbind(
        kingdom = ezFrame(Type="character",  DefaultValue="Bacteria",  Description="annotation mode and genetic code. Default is Bacteria")
	)
      }
    )
  )

