###################################################################
### git-init-24 (ruby script in /srv/GT/analysis/masaomi/2022/FGCZ/ check it out) create a test folder on github to test scripts
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodHifiasm = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  Hifi_Path = input$getFullPaths("Read1")
  Hifi_File = basename(input$getColumn("Read1"))
  inputFileType = param$inputType
  cmd = paste("hifiasm", "-o", sampleName,
                paste0("-t", ezThreads(),
                paste0("-l", param$purgingLevel), 
                paste0("-s", param$similarityThreshold),
                paste0("--n-hap", param$haplotypeNumber)), opt, 
                input$getFullPaths("Read1"), 
                "1> ", paste0(sampleName,"_hifiasm.log")) #(filename.join(basename(output$getColumn("column"), "vcfstats"), (input$getColumn()) input from gstore (check VcfStats) 
    ezSystem(cmd)
    OutputFile <- paste0(sampleName,".fasta")
    return("Success")
  }
}

##' @template app-template
##' @templateVar method ezMethodHifiasm()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppHifiasm <-
  setRefClass("EzAppHifiasm",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodHifiasm
                  name <<- "EzAppHifiasm"
                  appDefaults <<- rbind(purgingLevel = ezFrame(Type="integer",  DefaultValue="3",  Description="purge level"), 
                                        similarityThreshold = ezFrame(Type="integer",  DefaultValue="0.55",  Description="similarity threshold"), 
                                        haplotypeNumber = ezFrame(Type="integer",  DefaultValue="2",  Description="haplotype number"))
                }
              )
  )