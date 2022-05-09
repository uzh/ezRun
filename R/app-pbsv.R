###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodPbsv = function(input=NA, output=NA, param=NA){

  require("VariantAnnotation")
  
  vcfOutputFile = output$getColumn("VCF")
  PbsvLogFile = output$getColumn("PbsvLog")
  
  bamFiles = input$getFullPaths("BAM")
  bamDataset = input$meta
  genomeSeq = param$ezRef["refFastaFile"]
  nBamsInParallel = min(4, param$cores)
  bamFilesDiscover = ezMclapply(names(bamFiles), function(sampleName){
    setwdNew(paste(sampleName, "discover", sep="-"))
    bf = bamFiles[sampleName]
    osf = file.path(getwd(), sub(".bam", ".svsig.gz", basename(bf)))
    cmd = paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbsv discover", param$afOptions)
    if(param$ReadOpt=="CCS/HIFI"){
	    cmd=paste(cmd, "--hifi")
    }
    if (!is.null(param$region) && param$region != ""){
	    cmd=paste(cmd, "-r", param$region)
    }
    cmd=paste(cmd, bf, osf)

    ezSystem(cmd)
    
    setwd("..")
    
    return(osf)
  }, mc.cores=nBamsInParallel, mc.preschedule=FALSE)
  
  PbsvCallCmd = paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pbsv call --log-file", PbsvLoFile, param$types, param$minL, param$afOptions, param$callOptions, param$filterOptions, param$cmdOptions)
  
  if(param$ReadOpt=="CCS/HIFI"){
	    PbsvCallCmd=paste(PbsvCallCmd, "--hifi")
  }

  if (!is.null(param$region) && param$region != ""){
	    PbsvCallCmd=paste(PbsvCallCmd, "-r", param$region)
  }

  PbsvCallCmd=paste(PbsvCallCmd, genomeSeq, paste(bamFilesClean, collapse=" "), "Pbsv_Variants.vcf")
  ezSystem(PbsvCallCmd)
  ezSystem("bgzip Pbsv_Variants.vcf") 
  indexTabix(basename(vcfOutputFile),format = "vcf")
  ezSystem("SURVIVOR stats Pbsv_Variants.vcf.gz -1 -1 -1 Pbsv_Variants.stats.txt") 
  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMpileup(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppPbsv <-
  setRefClass("EzAppPbsv",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodPbsv
                  name <<- "EzAppPbsv"
                  appDefaults <<- rbind(region=ezFrame(Type="character",  DefaultValue="",  Description="The region of the genome. You can give either a chromosome name or a region on a chromosome like chr1 or chr1:1000-2000"),
                                    ReadOpt=ezFrame(Type="character",  DefaultValue="",	Description="input read types: SUBREAD, CCS/HIFI. Default is CCS/HIFI"),
                                    types=ezFrame(Type="character", DefaultValue="", Description="SV types to call: DEL,INS,INV,DUP,BND"),
                                    minL=ezFrame(Type="integer", DefaultValue="", Description="Ignore SV with length shorter than the given length (bp)"),
                                    afOptions=ezFrame(Type="character", DefaultValue="", Description="The options to filter alignments"),
                                    callOptions=ezFrame(Type="character", DefaultValue="", Description="The options to call"),
                                    filterOptions=ezFrame(Type="character", DefaultValue="", Description="The options to filter"))
                }
              )
  )
