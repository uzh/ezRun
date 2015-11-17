###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Teqc
##' @seealso \code{\link{EzAppTeqc}}
ezMethodTeqc = function(input=NA, output=NA, param=NA){
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(basename(output$getColumn("Report")))
  teqc(input, param)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodTeqc()
##' @seealso \code{\link{ezMethodTeqc}}
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

##' @title Target enrichment quality control
##' @description Performs a target enrichment quality control and creates reports of the outcome.
##' @template input-template
##' @param param a list of parameters, some of which are passed to \code{runTEQC()}:
##' \itemize{
##'  \item{designFile}{ a file describing the regions selected by the enrichment kit.}
##'  \item{name}{ a character representing the name of the app.}
##'  \item{covUniformityPlot}{ a logical indicating whether to generate plots for coverage uniformity.}
##'  \item{covTargetLengthPlot}{ a logical indicating whether to generate plots for coverage vs. target length.}
##'  \item{duplicatesPlot}{ a logical indicating whether to generate plots for duplicates.}
##'  \item{paired}{ a logical indicating whether the samples are paired.}
##' }
##' @template roxygen-template
teqc = function(input, param=NULL){
  require(TEQC)
  if(basename(param$designFile) == param$designFile){
    path = file.path("/srv/GT/databases/targetEnrichment_designs", param$designFile)
    param$designFile = list.files(path, pattern='Covered\\.bed$', full.names = T)[1]
  }
  samples = input$getNames()
  jobList = input$getFullPaths(param, "BAM")
  #Create one Report per Sample:
  ezMclapply(jobList, runTEQC, param, mc.cores=ezThreads())
  
  #Create MultiSampleReport:
  reportDirs = unlist(jobList)
  reportDirs = paste0("report_",gsub('\\.bam', '', basename(reportDirs)))
  multiTEQCreport(singleReportDirs=reportDirs,
                    samplenames=samples,
                    projectName=param$name,
                    targetsName=basename(dirname(param$designFile)),
                    referenceName="Human Genome hg19",
                    destDir="multiTEQCreport",
                    k = c(1,5,10,20,30,50),
                    figureFormat = c("png"))
  
  capture.output(print(sessionInfo()), file = "sessionInfo.txt")
  htmlFile="00index.html"
#   titles = list()
#   titles[["TEQC-Report"]] = paste("TEQC-Report:", param$name)
#   doc = openBsdocReport(title=titles[[length(titles)]], dataset=input$meta)
  html = openHtmlReport(htmlFile, param=param, title=paste("TEQC-Report:", param$name),
                        dataset=input$meta)
#   titles[["MultiSample-Report"]] = "MultiSample-Report"
#   addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  ezWrite("<h2>MultiSample-Report</h2>",con=html)
#   addTxtLinksToReport(doc, "multiTEQCreport/index.html")
  writeTxtLinksToHtml('multiTEQCreport/index.html', con=html)
#   titles[["Individual Reports"]] = "Individual Reports"
#   addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  ezWrite("<h2>Individual Reports</h2>",con=html)
#   addTxtLinksToReport(doc, paste0(reportDirs, '/index.html'))
  writeTxtLinksToHtml(paste0(reportDirs, '/index.html'), con=html)
#   >>>remove this title<<<
  ezWrite("<h2>Misc</h2>",con=html)
#   closeBsdocReport(doc, htmlFile, titles)
  writeTxtLinksToHtml('sessionInfo.txt',con=html)
  flush(html)
  closeHTML(html)
  return("Success")
}

##' @describeIn teqc Runs the target enrichment QC.
runTEQC = function(file, param){
  readsfile = file
  targetsfile = param$designFile
  TEQCreport(sampleName=gsub('\\.bam','',basename(file)),
               CovUniformityPlot = param$covUniformityPlot, CovTargetLengthPlot = param$covTargetLengthPlot, duplicatesPlot=param$duplicatesPlot,#CovGCPlot = T,
               k = c(1,5,10,20,30,50),
               targetsName=basename(dirname(targetsfile)),
               referenceName='hg19',
               pairedend=param$paired,
               destDir=paste0("report_",gsub('\\.bam', '', basename(file))),
               reads=get.reads(readsfile,filetype="bam"),
               targets=get.targets(targetsfile, 
                                   skip=grep("^track", readLines(targetsfile, n=200))),
               genome='hg19',figureFormat = c("png"))
  return("Success")
}

#prepareEnvironment = function(inputDatasetFile=NA, output=NA, param=NA){
#  param = fillWithDefaults(param)
#  options(cores=param$cores)
  #waitUntilFileExists(inputDatasetFile, maxWaitSeconds=300)
#  checkFreeDiskSpace(param)
#  setwdNew(basename(output$Report))
#  param$Name = basename(output$Report)
#  return(param)
#}
