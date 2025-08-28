###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDelly <- function(input = NA, output = NA, param = NA, htmlFile="00index.html") {
  require("VariantAnnotation")
  genomeSeq = param$ezRef["refFastaFile"]
  sampleName = input$getNames()
  setwdNew(sampleName)

  bamFile <- input$getFullPaths("BAM")
  ovf = paste0(sampleName, ".vcf")
  ovzf = paste0(sampleName, ".vcf.gz")
  DellyLogFile = paste0(sampleName, "_delly.log")
  DellyStatsFile = paste0(sampleName, ".stats.txt")
  SurvivorLogFile = paste0(sampleName, ".SURVIVOR.log")
  cmd = "delly"
  if(param$ReadOpt=="sr"){
	    cmd=paste(cmd, "call")
  }else if (param$ReadOpt=="pb"){
	    cmd=paste(cmd, "lr -y pb")
  }else if (param$ReadOpt=="ont"){
	    cmd=paste(cmd, "lr -y ont")
  }
  cmd=paste(cmd, "-g", param$ezRef["refFastaFile"], bamFile, ">", ovf, "2>", DellyLogFile)
  ezSystem(cmd)
  ezSystem(paste("bgzip -c", ovf, ">", ovzf)) 
  indexTabix(basename(ovzf),format = "vcf")

  SurvivorCmd=paste("SURVIVOR stats", ovf, "20 -1 3",  DellyStatsFile, "1>", SurvivorLogFile) 
  ezSystem(SurvivorCmd) 

  ##html file  
  #setwd(start_path)
  htmlFile = output$getColumn("OutReport")
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "Pbsv.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  params = list(sample=sampleName)
  rmarkdown::render(input="Pbsv.Rmd", envir = new.env(), output_dir=".", output_file=htmlFile, quiet=TRUE)

  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodPbsv(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getPbmm2Reference}}
EzAppPbsv <-
  setRefClass("EzAppPbsv",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodPbsv
        name <<- "EzAppPbsv"
        appDefaults <<- rbind(
        ReadOpt = ezFrame(Type="character",  DefaultValue="HIFI",  Description="input read types: SUBREAD, CCS, HIFI. Default is HIFI")
	)
      }
    )
  )

