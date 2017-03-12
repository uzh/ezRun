###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodSingleCellCounts = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){

  input = EzDataset(file=input$getFullPaths("ReadDataset"), dataRoot=param$dataRoot)
  if (ezIsSpecified(param$samples)){
    input = input$subset(param$samples)
  }
  countDir = basename(output$getColumn('CountFolder'))
  dir.create(countDir)
  bamMeta = input$meta[ , !input$columnHasTag("File")]
  bamMeta[["BAM [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), ".bam")
  bamMeta[["BAI [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), ".bam.bai")
  bamMeta[["Read Count"]] = ceiling(bamMeta[["Read Count"]] / param$subsampleReads)
  bamMeta[["STARLog [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), "_STAR.log")
  bamOutput = EzDataset(meta=bamMeta, param$dataRoot)
  bamParam = param
  bamParam$mail = ""
  countMeta = input$meta[ , !input$columnHasTag("File")]
  countMeta[['Count [File]']] = paste0(rownames(countMeta), '.txt')
  countMeta[['Stats [File]']] = paste0(rownames(countMeta), '-stats.txt')
  ezWrite.table(countMeta, file=basename(output$getColumn('CountDataset')), head='Name')
  countDs = EzDataset(meta=countMeta)
  
  switch(param$mapMethod,
         STAR={
           refDir = getSTARReference(param)
           mappingApp = EzAppSTAR$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All")
           if (!grepl("genomeLoad", bamParam$cmdOptions)){
             bamParam$cmdOptions = paste("--genomeLoad LoadAndKeep" , bamParam$cmdOptions)
           }
           if (!grepl("outSAMattributes", bamParam$cmdOptions)){
             bamParam$cmdOptions = paste(bamParam$cmdOptions, "--outSAMattributes All")
           }

           on.exit({
             if (ezIsSpecified(param$ezRef["refIndex"])){
               refDir = param$ezRef["refIndex"]
             } else {
               if (!ezIsSpecified(param$ezRef["refFeatureFile"])){
                 stop("refFeatureFile not defined")
               }
               refDir = sub(".gtf$", "_STARIndex", param$ezRef["refFeatureFile"])
             }
             ezSystem(paste(STAR, '--genomeDir', refDir, '--genomeLoad Remove'))
           })
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
  
  
  for (nm in input$getNames()){
    setwdNew(nm)
    readI = input$subset(nm)
    bamI = bamOutput$subset(nm)
    countI = countDs$subset(nm)
    mappingApp$run(input=readI, output=bamI, param=bamParam)
    starLogfile = paste0("../", countDir, '/', nm, '-STAR.log')
    file.rename('Log.final.out', to = starLogfile)
    EzAppFeatureCounts$new()$run(input=bamI, output=countI, param=bamParam)
    file.rename(basename(countI$getColumn("Count")), paste0("../", countDir, "/", basename(countI$getColumn("Count"))))
    file.rename(basename(countI$getColumn("Stats")), paste0("../", countDir, "/", basename(countI$getColumn("Stats"))))
    ezSystem(paste("mv", basename(bamI$getColumn("BAM")), file.path("..", countDir) ))
    ezSystem(paste("mv", basename(bamI$getColumn("BAI")), file.path("..", countDir) ))
    setwd('..')
    unlink(nm, recursive=TRUE, force=TRUE)
    
  }
  
  countFiles = paste0(countDir, "/", basename(countDs$getColumn("Count")))
  templateFile = file.path(countDir, countDs$subset( input$getNames()[1] )$getColumn("Count"))
  x = ezRead.table(templateFile)
  counts = ezMatrix(0, rows=rownames(x), cols=countDs$getNames())
  for (nm in input$getNames()) {
    countFile = file.path(countDir, countDs$subset(nm)$getColumn("Count"))
    x = ezRead.table(countFile)
    stopifnot(setequal(rownames(x), rownames(counts)))
    counts[rownames(x), nm] = x[ , "matchCounts"]
  }
  ezWrite.table(counts, head=paste0(param$featureLevel, "_id"), file=basename(output$getColumn('CountMatrix')))
  return("SUCCESS")
}

##' @template app-template
##' @templateVar method ezMethodSingleCellCounts(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppSingleCellCounts <-
  setRefClass("EzAppSingleCellCounts",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSingleCellCounts
                  name <<- "EzAppSingleCellCounts"
                  appDefaults <<- rbind(mapMethod=ezFrame(Type="character",	DefaultValue="STAR",	Description="the mapper to use"),
                                        mapOptions=ezFrame(Type="character", DefaultValue="", Description="options passed to the mapper"),
                                        writeIgvSessionLink=ezFrame(Type="logical", DefaultValue=FALSE, Description="whether to write IGV session links.")
                                        )
                }
              )
  )
