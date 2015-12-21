###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template app-template
##' @templateVar method ezMethodTeqc
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppTeqc <-
  setRefClass("EzAppTeqc",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodTeqc
                  name <<- "EzAppTeqc"
                  appDefaults <<- rbind(designFile=ezFrame(Type="character",  DefaultValue="",  Description="file describing the regions selected by the enrichment kit"),
                                    covUniformityPlot=ezFrame(Type="logical", DefaultValue="TRUE", Description="generate plots for coverage uniformity?"),
                                    covTargetLengthPlot=ezFrame(Type="logical", DefaultValue="TRUE", Description="generate plots for coverage vs. target length"),
                                    duplicatesPlot=ezFrame(Type="logical", DefaultValue="TRUE", Description="generate plots for duplicates"))
                }
              )
  )

ezMethodTeqc = function(input=NA, output=NA, param=NA){
  setwdNew(basename(output$getColumn("Report")))
  if(basename(param$designFile) == param$designFile){
    param$designFile = list.files(file.path(TEQC_DESIGN_DIR, param$designFile), 
                                  pattern='Covered\\.bed$', full.names = T)[1]
  }
  samples = input$getNames()
  jobList = input$getFullPaths(param, "BAM")
  #Create one Report per Sample:
  destDirs = ezMclapply(jobList, runTEQC, param, mc.cores=ezThreads())
  
  #Create MultiSampleReport:
  destDirs = unlist(destDirs)
  TEQC::multiTEQCreport(singleReportDirs=destDirs,
                    samplenames=samples,
                    projectName=param$name,
                    targetsName=basename(dirname(param$designFile)),
                    referenceName="Human Genome hg19",
                    destDir="multiTEQCreport",
                    k = c(1,5,10,20,30,50),
                    figureFormat = c("png"))
  
  capture.output(print(sessionInfo()), file = "sessionInfo.txt")
  htmlFile="00index.html"
  
  titles = list()
  titles[["TEQC-Report"]] = paste("TEQC-Report:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  
  addDataset(doc, input$meta, param)
  
  titles[["MultiSample-Report"]] = "MultiSample-Report"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addTxtLinksToReport(doc, "multiTEQCreport/index.html")
  titles[["Individual Reports"]] = "Individual Reports"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addTxtLinksToReport(doc, paste0(destDirs, '/index.html'))
  closeBsdocReport(doc, htmlFile, titles)
  return("Success")
}

##' @title Runs the target enrichment quality control
##' @description Performs a target enrichment quality control and creates reports of the outcome.
##' @param param a list of parameters:
##' \itemize{
##'  \item{designFile}{ a file describing the regions selected by the enrichment kit.}
##'  \item{covUniformityPlot}{ a logical indicating whether to generate plots for coverage uniformity.}
##'  \item{covTargetLengthPlot}{ a logical indicating whether to generate plots for coverage vs. target length.}
##'  \item{duplicatesPlot}{ a logical indicating whether to generate plots for duplicates.}
##'  \item{paired}{ a logical indicating whether the samples are paired.}
##' }
##' @param file a character representing the path to the file containing the reads.
##' @template roxygen-template
runTEQC = function(file, param){
  sampleName = gsub('\\.bam', '', basename(file))
  destDir = paste0("report_", sampleName)
  targetsfile = param$designFile
  TEQC::TEQCreport(sampleName=sampleName,
                   CovUniformityPlot = param$covUniformityPlot, CovTargetLengthPlot = param$covTargetLengthPlot, duplicatesPlot=param$duplicatesPlot,#CovGCPlot = T,
                   k = c(1,5,10,20,30,50),
                   targetsName=basename(dirname(targetsfile)),
                   referenceName='hg19',
                   pairedend=param$paired,
                   destDir=destDir,
                   reads=TEQC::get.reads(file, filetype="bam"),
                   targets=TEQC::get.targets(targetsfile, 
                                       skip=grep("^track", readLines(targetsfile, n=200))),
                   genome='hg19',figureFormat = c("png"))
  return(destDir)
}
