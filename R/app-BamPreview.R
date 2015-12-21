###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodBamPreview = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){

  if (ezIsSpecified(param$samples)){
    input$subset(param$samples)
  }
  
  bamMeta = input$meta[ , !input$columnHasTag("File")]
  bamMeta[["BAM [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), ".bam")
  bamMeta[["BAI [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), ".bam.bai")
  bamMeta[["Read Count"]] = ceiling(bamMeta[["Read Count"]] / param$subsampleReads)
  bamOutput = EzDataset(meta=bamMeta)
  bamParam = param
  bamParam$mail = ""
  switch(param$mapMethod,
         STAR={
           mappingApp = EzAppSTAR$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif")
         },
         bowtie={
           mappingApp = EzAppBowtie$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "")
         },
         bowtie2={
           mappingApp = EzAppBowtie2$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--no-unal")
         },
         tophat={
           mappingApp = EzAppTophat$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--mate-inner-dist 100 --mate-std-dev 150")
         },
         "bwa-mem"={
           mappingApp = EzAppBWA$new()
           bamParam$algorithm = "mem"
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "")
         },
         stop("unsupported mapMethod: ", param$mapMethod)
  )
  for (i in 1:nrow(input$meta)){
    setwdNew(input$getNames()[i])
    mappingApp$run(input=input$copy()$subset(i), output=bamOutput$copy()$subset(i), param=bamParam) ## TODO: Read Count of the bamOutput should be set by the mapping app.
    setwd("..")
  }
  param$dataRoot = ""
  result = ezMethodRnaBamStats(bamOutput, output, param)
  return(result)
}

##' @template app-template
##' @templateVar method ezMethodBamPreview
##' @templateVar htmlArg , htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppBamPreview <-
  setRefClass("EzAppBamPreview",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBamPreview
                  name <<- "EzAppBamPreview"
                  appDefaults <<- rbind(mapMethod=ezFrame(Type="character",	DefaultValue="STAR",	Description="the mapper to use"),
                                        mapOptions=ezFrame(Type="character", DefaultValue="", Description="options passed to the mapper"),
                                        fragSizeMax=ezFrame(Type="integer",  DefaultValue=500,	Description="maximum fragment size to plot in fragment size distribution"),
                                        writeIgvSessionLink=ezFrame(Type="logical",  DefaultValue=FALSE,	Description="whether to write IGV session links."))
                }
              )
  )
