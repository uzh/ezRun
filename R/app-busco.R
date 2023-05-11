###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

#busco -i /srv/gstore/projects/p25396/o30319_Spades_2023-02-06--14-53-26/DAMREC0027.fasta -o busco_test -m geno -l gammaproteobacteria_odb10 -c 4 
#generate_plot.py -wd busco_test/
#/srv/GT/analysis/qiwei/p25396/o30319/o30319_Prodigal_2023-02-07--12-29-31

ezMethodBusco = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  draft = input$getFullPaths("Draft")
  cmd = paste("busco -i ", draft, "-o", sampleName, "-l", param$lineage, "-m", param$mode, "-c", ezThreads(), opt, "1>", paste0(sampleName,"_busco.log"))
  ezSystem(cmd)
  cmd = paste ("generate_plot.py -wd", sampleName)
  ezSystem(cmd)
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodProkka
##' @description Use this reference class to run
##' @seealso \code{\link{getPbmm2Reference}}
EzAppBusco <-
  setRefClass("EzAppBusco",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBusco
        name <<- "EzAppBusco"
        appDefaults <<- rbind(
        kingdom = ezFrame(Type="character",  DefaultValue="Bacteria",  Description="annotation mode and genetic code. Default is Bacteria")
	)
      }
    )
  )

