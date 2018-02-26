###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCCountQC <-
  setRefClass("EzAppSCCountQC",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCCountQC
                  name <<- "EzAppSCCountQC"
                  #appDefaults <<- rbind(runGO=ezFrame(Type="logical", DefaultValue=TRUE, Description="whether to run the GO analysis"))
                }
              )
  )

ezMethodSCCountQC = function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  if(length(input$getNames()) > 1L)
    stop("Currently we support one pooled bam file!")
  
  setwdNew(basename(output$getColumn("Report")))
  rawData = loadSCCountDataset(input, param)
  
  ## CollectAlignmentSummaryMetrics
  alnMetrics <- CollectAlignmentSummaryMetrics(inBam=input$getFullPaths("BAM"),
                                               fastaFn=param$ezRef['refFastaFile'],
                                               metricLevel="SAMPLE")
  
  ## CollectRnaSeqMetrics
  rnaSeqMetrics <- CollectRnaSeqMetrics(inBam=input$getFullPaths("BAM"),
                                        gtfFn=param$ezRef['refFeatureFile'],
                                        featAnnoFn=param$ezRef['refAnnotationFile'],
                                        strandMode=param$strandMode,
                                        metricLevel="SAMPLE")
  
  ## DuplicationMetrics
  
}
