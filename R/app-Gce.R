###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodGce = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  setwdNew(sampleName)
  cmd = paste("echo", input$getFullPaths("Read1"), ">", paste0(sampleName, ".lib"))
  ezSystem(cmd)
  cmd = paste("/usr/local/ngseq/src/kmerfreq/kmerfreq -k", param$kSize, "-t", ezThreads(), "-f 1 -w 1 -c 1 -q 1 -p", sampleName, paste0(sampleName, ".lib"))
  ezSystem(cmd)
  cmd = paste("perl /usr/local/ngseq/src/GCE/gce-alternative/estimate_genome_character_real.pl", paste0(sampleName, "kmer.freq.stat"))
  ezSystem(cmd)
  cmd = paste("convert", paste0(sampleName, "kmer.freq.stat_real_theory_distribution.svg"), paste0(sampleName, "kmer.freq.stat_real_theory_distribution.png"))
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodGce()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppGce <-
  setRefClass("EzAppGce",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGce
                  name <<- "EzAppGce"
                }
              )
)
              
              
