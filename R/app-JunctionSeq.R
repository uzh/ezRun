###QoRTs/JunctionSeqApp
ezMethodRunJunctionSeqPipeline <- function(input, output, param, htmlFile="00index.html"){
  setwdNew(param$comparison)
  runQoRTs(input, output, param)
  runJunctionSeq(input, output, param)
  ###TODO: multigroup support, order of levels in dataset (resort dataset by groupname?), usage of flattend gff?? only known SJ, result annotation -> simple report,
  ##add parameters for junctionSeq call
  return('success')
}

runQoRTs <- function(input, output, param){
  samples = input$getNames()
  dataset = input$meta
  javaCall = paste0('java -Xmx', param$ram, 'G', ' -jar /usr/local/ngseq/src/QoRTs/QoRTs.jar')
  #1. Get BAM-Files:
  bamFiles = input$getFullPaths("BAM")
  bamFileList = as.list(bamFiles)
  ezWrite.table(data.frame(sample.ID = rownames(dataset), dataset, stringsAsFactors = F, check.names = F), 'decoderFile.txt', row.names = FALSE)
  gtfFile = param$ezRef@refFeatureFile
  
  #Make flattend gtfFile
  gffOutput = 'JunctionSeq.gff'
  cmd = paste(javaCall, 'makeFlatGff', gtfFile, gffOutput)
  ezSystem(cmd)
  
  #Get Junctions per Sample
  for (i in c(1:length(bamFiles))){
    setwdNew(samples[i])
    if(param$paired){
      cmd = paste(javaCall, 'QC --stranded --maxReadLength', param$maxReadLength, '--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon',
                bamFiles[i], gtfFile, '.')
    } else {
      cmd = paste(javaCall, 'QC --stranded --singleEnded --maxReadLength', param$maxReadLength, '--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon',
                  bamFiles[i], gtfFile, '.')
    }
    ezSystem(cmd)
    setwd('..')
  }
  
  #Merge Junctions and create new GFF which includes Novel Junctions
  cmd = paste(javaCall, 'mergeNovelSplices  --minCount', param$minCount, '--stranded .', 'decoderFile.txt', gtfFile, '.')
  ezSystem(cmd)
  
  return('success')
}

#Run junctionSeq
runJunctionSeq <- function(input, output, param, htmlFile="00index.html"){
  require("JunctionSeq")
  samples = input$getNames()
  dataset = input$meta
  gffOutputIncludingNovelJunctions = 'withNovel.forJunctionSeq.gff.gz'
  # including novel junctions
  countFiles <- file.path(samples,"QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz");
  param[['grouping']] = sub('.\\[.*', '', param[['grouping']])
  colnames(dataset) = sub('.\\[.*', '', colnames(dataset))
  jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                                 sample.names = samples,
                                 condition = factor(dataset[[param$grouping]]),
                                 flat.gff.file = gffOutputIncludingNovelJunctions,
                                 nCores = param$cores,
                                 analysis.type = "junctionsAndExons")
  
  writeSizeFactors(jscs, file = "sizeFactors.txt")
  
  writeCompleteResults(jscs,
                       outfile.prefix="./JuncSeq_Result-File",
                       save.jscs = TRUE, 
                       FDR.threshold = as.numeric(param$fdr));
  
  buildAllPlots(jscs=jscs,
                outfile.prefix = "./finalOutput",
                use.plotting.device = "CairoPNG",
                FDR.threshold = as.numeric(param$fdr))
  return('sucess')
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodRunQoRTs(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run JunctionSeq on bamFiles for alternative Splicing analysis
EzAppJunctionSeq <-
  setRefClass(Class = "EzAppJunctionSeq",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodRunJunctionSeqPipeline
                  name <<- "EzAppJunctionSeq"
                  appDefaults <<- rbind(maxReadLength = ezFrame(Type="numeric", DefaultValue=151, Description="define max read length"),
                                        minCount = ezFrame(Type="numeric", DefaultValue=6, Description="min Junction Count"))
                }
              ))
