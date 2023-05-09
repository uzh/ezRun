###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodPbsv <- function(input = NA, output = NA, param = NA) {
  require("VariantAnnotation")
  require(tidyr)
  require(ggplot2)
  genomeSeq = param$ezRef["refFastaFile"]
  sampleName = input$getNames()
  bamFile <- input$getFullPaths("BAM")
  osf = paste0(sampleName, ".svsig.gz")
  ovf = paste0(sampleName, ".vcf")
  ovzf = paste0(sampleName, ".vcf.gz")
  PbsvLogFile = paste0(sampleName, "_pbsv.log")
  PbsvStatsFile = paste0(sampleName, ".stats.txt")
  PbsvStatsPlot = paste0(sampleName, ".stats.pdf")
  cmd = paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbsv discover", param$afOptions)
  if(param$ReadOpt=="HIFI"){
	    cmd=paste(cmd, "--hifi")
  }else if (param$ReadOpt=="CCS"){
	    cmd=paste(cmd, "--ccs")
  }
  if (!is.null(param$region) && param$region != ""){
	    cmd=paste(cmd, "-r", param$region)
  }
  cmd=paste(cmd, bamFile, osf)
  ezSystem(cmd)

  
  PbsvCallCmd = paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbsv call --log-file", PbsvLogFile, "--types", param$types, "--min-sv-length", param$minL, param$callOptions, param$filterOptions, param$cmdOptions)
  
  if(param$ReadOpt=="HIFI"){
	    PbsvCallCmd=paste(PbsvCallCmd, "--hifi")
  }else if (param$ReadOpt=="CCS"){
	    PbsvCallCmd=paste(PbsvCallCmd, "--ccs")
  }
  if (!is.null(param$region) && param$region != ""){
	    cmd=paste(cmd, "-r", param$region)
  }
 

  PbsvCallCmd=paste(PbsvCallCmd, genomeSeq, osf, ovf)
  ezSystem(PbsvCallCmd)
  ezSystem(paste("bgzip -c", ovf, ">", ovzf)) 
  indexTabix(basename(ovzf),format = "vcf")
  ezSystem(paste("SURVIVOR stats", ovf, "20 -1 3",  PbsvStatsFile)) 

  plotPbsv(PbsvStatsFile, sampleName)

  return("Success")
}

plotPbsv <- function(PbsvStatsFile, sampleName){
  	tbl<-read.table(PbsvStatsFile, header=TRUE)
  	data_long<- gather(tbl, svType, count, Del:TRA, factor_key=TRUE)
  	data_long$Len<-factor(data_long$Len, levels = c("0-50bp", "50-100bp", "100-1000bp", "1000-10000bp", "10000+bp"))
  	p1<-ggplot(data=data_long, aes(x=Len, y=count, color=svType, fill=svType)) 
  	p1<-p1 + geom_bar(stat="identity") 
  	p1<-p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  	p1<-p1 + facet_grid(. ~ svType) + scale_y_continuous(trans='log10') + ggtitle(sampleName)
  	pdf(paste0(sampleName, '.stats.pdf'))
    		print(p1)
  	dev.off()
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

