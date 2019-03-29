###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @template app-template
##' @templateVar method ezMethodSingleCellCounts(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppSCCounts <-
  setRefClass("EzAppSCCounts",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCCounts
                  name <<- "EzAppSCCounts"
                  appDefaults <<- rbind(mapMethod=ezFrame(Type="character",	
                                                          DefaultValue="STAR",
                                                          Description="the mapper to use"),
                                        mapOptions=ezFrame(Type="character",
                                                           DefaultValue="",
                                                           Description="options passed to the mapper"),
                                        controlSeqs=ezFrame(Type="charVector",
                                                            DefaultValue="",
                                                            Description="control sequences to add"),
                                        writeIgvSessionLink=ezFrame(Type="logical",
                                                                    DefaultValue=FALSE,
                                                                    Description="whether to write IGV session links.")
                  )
                }
              )
  )


ezMethodSCCounts = function(input=NA, output=NA, param=NA,
                            htmlFile="00index.html"){
  require(readr)
  require(tibble)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("ResultDir")))
  on.exit(setwd(cwd), add=TRUE)
  
  metainput = EzDataset(file=input$getFullPaths("CellDataset"),
                        dataRoot=param$dataRoot)
  
  cellMeta = metainput$meta[ , !metainput$columnHasTag("File")]
  cellMeta[['CountMatrix [File]']] = output$getColumn("CountMatrix")
  
  output$setColumn("BAM [File]", sub("-counts\\.mtx$", ".bam",
                                     output$getColumn("CountMatrix")))
  cellMeta[["BAM [File]"]] = output$getColumn("BAM")
  
  output$setColumn("BAI [File]", sub("-counts\\.mtx$", ".bai",
                                     output$getColumn("CountMatrix")))
  cellMeta[["BAI [File]"]] = output$getColumn("BAI")
  
  output$setColumn("STARLog [File]", sub("-counts\\.mtx$", "_STAR.log",
                                     output$getColumn("CountMatrix")))
  cellMeta[["STARLog [File]"]] = output$getColumn("STARLog")
  
  output$setColumn("PreprocessingLog [File]", sub("-counts\\.mtx$", "_preprocessing.log",
                                                  output$getColumn("CountMatrix")))
  cellMeta[["PreprocessingLog [File]"]] = output$getColumn("PreprocessingLog")
  
  output$setColumn("Stats [File]", sub("-counts\\.mtx$", "_stats.txt",
                                       output$getColumn("CountMatrix")))
  cellMeta[['Stats [File]']] = output$getColumn("Stats")
  
  output$setColumn("CellCyclePhase [File]", sub("-counts\\.mtx$", "-CellCyclePhase.txt",
                                                output$getColumn("CountMatrix")))
  cellMeta[['CellCyclePhase [File]']] = output$getColumn("CellCyclePhase")
  
  write_tsv(as_tibble(cellMeta, rownames="Name"),
            path=sub("-counts\\.mtx$", "-dataset.tsv",
                     basename(cellMeta[['CountMatrix [File]']])))
  
  bamParam = param
  #bamParam$mail = "" Let's send several emails during SCCounts. It's a long computation.

  switch(param$mapMethod,
         STAR={
           mappingApp = EzAppSingleCellSTAR$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All")
           if (!grepl("outSAMattributes", bamParam$cmdOptions)){
             bamParam$cmdOptions = paste(bamParam$cmdOptions, "--outSAMattributes All")
           }
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

  mappingApp$run(input=input, output=output, param=bamParam)
  
  ## Prepare the input for featurecounts
  featurecountsInput = output$copy()
  featurecountsInput$dataRoot <- ""
  featurecountsInput$setColumn("BAM",
                               basename(output$getColumn("BAM")))
  EzAppSingleCellFeatureCounts$new()$run(input=featurecountsInput,
                                         output=output, param=bamParam)
  return("SUCCESS")
}
