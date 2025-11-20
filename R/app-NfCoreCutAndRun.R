###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNfCoreCutAndRun <- function(input = NA, output = NA, param = NA) {
  refbuild = param$refBuild
  outFolder = output$getColumn("CutAndRun_Result") |> basename()

  fullGenomeSize <- param$ezRef@refFastaFile %>% Rsamtools::FaFile() %>% GenomeInfoDb::seqlengths() %>% sum()
  effectiveGenomeSize <- (fullGenomeSize * 0.8 ) %>% round()

  nfSampleFile <- file.path('dataset.csv')
  nfSampleInfo = getCutAndRunSampleSheet(input, param)
  write_csv(nfSampleInfo, nfSampleFile)


  blackListFile <- getBlackListFile(input, param)
  prepNFCoreEnv()
  configFile <- writeNextflowLimits(param)
  cmd = paste(
    "nextflow run nf-core/cutandrun",
     ## i/o
    "--input", nfSampleFile,
    "--outdir", outFolder,
    ## genome files
    "--fasta", param$ezRef@refFastaFile,
    "--gtf", param$ezRef@refFeatureFile,
    "--gene_bed", str_replace(param$ezRef@refAnnotationFile,
                              basename(param$ezRef@refAnnotationFile),
                              'genes.bed'),
    "--blacklist", blackListFile,
    # if (param[['blacklist']] == "")  "" else paste0("--blacklist ", param[['blacklist']]),
    ## parameters
    "--dt_calc_all_matrix false",
    "--macs_gsize", effectiveGenomeSize,
    "--peakcaller", param[['peakCaller']],
    "--spikein_genome", param[['spikeinGenome']],
    "--normalisation_mode", param[['normalization']],
    if (param[['peakStyle']] == 'broad')  "" else "--macs2_narrow_peak",
    ## configuration
    "-work-dir work",
    "-profile apptainer",
    "-r", param$pipelineVersion,
    "-c", configFile,
    param$cmdOptions
  )
  ezSystem(cmd)
  ezSystem(paste('mv', configFile, outFolder))
  getFastaFromBedFiles(outFolder, refFile = param$ezRef["refFastaFile"])
  getAnnotatedPeaks(gtfFile = param$ezRef@refFeatureFile, outFolder)
  jsonFile = writeCutAndRunIgvSession(param, outFolder, jsonFileName = paste0(outFolder, "/igv_session.json"), bigwigRelPath = "/04_reporting/igv/",
                      baseUrl = file.path(PROJECT_BASE_URL, output$getColumn("CutAndRun_Result")))
  writeNfCoreIgvHtml(param, jsonFile, title = "NfCoreCutAndRun MultiSample Coverage Tracks", htmlTemplate = "templates/igvNfCoreTemplate.html", htmlFileName = paste0(outFolder, "/igv_session.html"))
  makeRmdReportWrapper(outFolder, rmdFile="NfCoreCutAndRun.Rmd", reportTitle="NfCoreCutAndRun")

  if(ezIsSpecified(param$keepBams)){
    keepBams <- param$keepBams
  } else {
    keepBams <- TRUE
  }
  cleanupOutFolder(outFolder, keepBams)

  return("Success")
}

EzAppNfCoreCutAndRun <- setRefClass(
  "EzAppNfCoreCutAndRun",
  contains = "EzApp",
  methods = list(
    initialize = function()
    {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodNfCoreCutAndRun
      name <<- "EzAppNfCoreCutAndRun"
      ## minimum nf-core parameters
      appDefaults <<- rbind(
        peakCaller = ezFrame(Type="character", DefaultValue="macs2", Description="Select the peak caller for the pipeline"),
        spikeinGenome = ezFrame(Type="character", DefaultValue="macs2", Description="Select the reference for the spike-in genome"),
        normalization = ezFrame(Type="character", DefaultValue="macs2", Description="Select the target read normalization mode"),
        peakStyle  = ezFrame(Type="character", DefaultValue="broad", Description="Run MACS2 in broadPeak mode, otherwise in narrowPeak mode"),
        keepBams = ezFrame(Type="logical", DefaultValue = FALSE, Description= "Should bam files be stored"),
        pipelineVersion = ezFrame(Type="character", DefaultValue = '3.2.2', Description= "specify pipeline version")
      )
    }
  )
)

##' @description fetch nfcore blacklist files
getBlackListFile <- function(input, param){
  buildName = param$ezRef@refBuildName
  basePath = '/srv/GT/databases/nf-core/cutandrun/'
  blackListPath <- case_when(
    grepl('GRCm39', buildName, ignore.case = TRUE) ~ paste0(basePath, 'GRCm39-blacklist.bed'),
    grepl('GRCm38', buildName, ignore.case = TRUE) ~ paste0(basePath, 'GRCm38-blacklist.bed'),
    grepl('GRCh38', buildName, ignore.case = TRUE) ~ paste0(basePath, 'GRCh38-blacklist.bed'),
    grepl('GRCh37', buildName, ignore.case = TRUE) ~ paste0(basePath, 'GRCh37-blacklist.bed'),
    TRUE ~ ""
  )
  return(blackListPath)
}

##' @description get an nf-core/cutandrun-formatted csv file
getCutAndRunSampleSheet <- function(input, param){
  ## TODO: does not yet support comma-separated file paths
  nfSampleInfo <- ezFrame(
    group = input$getColumn(param$grouping),
    replicate = ezReplicateNumber(input$getColumn(param$grouping)),
    fastq_1 = input$getFullPaths("Read1"),
    fastq_2 = input$getFullPaths("Read2"),
    control = input$getColumn(param$controlColumn)
  )

  return(nfSampleInfo)
}

##' @description clean up NfCoreCutAndRun_result directory
cleanupOutFolder <- function(outFolder, keepBams=TRUE){
  if(!keepBams){
    bamPath <- paste0(outFolder,"/02_alignment/")
    bamsToDelete <- dir(path=bamPath, pattern="*.bam(.bai)?$", recursive=TRUE)
    file.remove(file.path(bamPath, bamsToDelete))
    cat("Deleted bam and bam.bai files form bwa directory.\n")
  }
  genomePath <- paste0(outFolder,"/04_reporting/igv/")
  filesToDelete <- dir(path=genomePath, pattern="genome.fa(.fai)?$")
  file.remove(file.path(genomePath, filesToDelete))
  cat("Deleted genome.fa and genome.fai files from 04_reporting/igv/ directory.\n")
}

##' @description write IGV session in json format
writeCutAndRunIgvSession <- function(param, outFolder, jsonFileName, bigwigRelPath, baseUrl){
  refBuildName = param$ezRef@refBuildName
  refUrlBase = file.path(REF_HOST, param$ezRef@refBuild)
  fastaUrl = sub("Annotation.*", "Sequence/WholeGenomeFasta/genome.fa", refUrlBase)
  faiUrl = paste0(fastaUrl, ".fai")

  bigwigPath=file.path(outFolder, bigwigRelPath)
  bigwigFiles <- dir(path=bigwigPath, pattern="*.bigWig$")
  trackNames <- bigwigFiles |> str_replace("\\..*", "")
  tracks <- list()
  tracks[[1]] <- list(type=	"sequence")
  for (i in 1:length(bigwigFiles)){
    tracks[[i+1]] <- list(id = trackNames[[i]],
                          url = paste0(baseUrl,file.path(bigwigRelPath, bigwigFiles[[i]])),
                          format =	"bigWig",
                          name	= trackNames[[i]])

  }
  jsonLines <- list( version =	"3.5.3",
                     showSampleNames = FALSE,
                     reference = list(id = refBuildName , fastaUrl = fastaUrl, indexURL = faiUrl),
                     tracks = tracks
  )
  jsonFile <- rjson::toJSON(jsonLines, indent=5, method="C")
  write(jsonFile, jsonFileName)
  return(jsonFile)
}

##' @description write html wrapper for IGV session
writeNfCoreIgvHtml = function(param, jsonFile, title, htmlTemplate, htmlFileName){
  htmlLines = readLines(system.file(htmlTemplate, package="ezRun", mustWork = TRUE))
  htmlLines = gsub("TITLE", title, htmlLines)
  htmlLines = gsub("IGV_JSON_CONTENT", jsonFile, htmlLines)
  writeLines(htmlLines, htmlFileName)
}


##' @description generate fasta files from BED files
getFastaFromBedFiles <- function(outFolder, refFile){
  bedFilePath <- paste0(outFolder,"/04_reporting/igv/")
  bedFileNames <- dir(path=bedFilePath, pattern=".bed$", recursive=TRUE, full.names = TRUE)
  for (name in bedFileNames){
    peakSeqFile = paste0(name, "_peaks.fa")
    cmd = paste("bedtools", " getfasta -fi", refFile, "-bed", name, " -name -fo ",peakSeqFile)
    ezSystem(cmd)
  }
}

##' @description annotate peaks in BED files
getAnnotatedPeaks <- function(gtfFile, outFolder){
  require(ChIPpeakAnno)
  require(GenomicRanges)
  require(rtracklayer)
  gtf <- rtracklayer::import(gtfFile)
  if(grepl('gtf$',gtfFile)){
    names_gtf = make.unique(gtf$'gene_id')
  } else {
    names_gtf = make.unique(gtf$'ID')
  }
  names(gtf) = names_gtf

  bedFilePath <- paste0(outFolder,"/04_reporting/igv/")
  bedFileNames <- dir(path=bedFilePath, pattern=".bed$", recursive=TRUE)
  for (name in bedFileNames){
    peakBedFile <- file.path(bedFilePath, name)
    peakAnnFile = paste0(peakBedFile, ".xlsx")
    myPeaks = ezRead.table(peakBedFile, row.names = NULL, header = F)
    peaksRD = makeGRangesFromDataFrame(myPeaks, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field="V1")
    annotatedPeaks <- annotatePeakInBatch(peaksRD, AnnotationData = gtf, output='nearestStart', multiple=FALSE, FeatureLocForDistance='TSS')
    annotatedPeaks = merge(myPeaks, annotatedPeaks,by.x = c('V6','V4','V5','V2','V3','V1'), by.y = c('V6','V4','V5','start','end','seqnames'), all.x = TRUE)
    writexl::write_xlsx(annotatedPeaks, peakAnnFile)
  }
}

##' @description write HTML report
makeRmdReportWrapper <- function(outFolder, rmdFile, reportTitle){
  plotsPath <- paste0(outFolder,"/04_reporting/deeptools_heatmaps/")
  filesToPlot <- dir(path=plotsPath, pattern=".pdf$", recursive=TRUE)
  cd = getwd()
  setwdNew(paste0(outFolder,"/04_reporting/"))
  makeRmdReport(filesToPlot=file.path("./deeptools_heatmaps", filesToPlot), rmdFile=rmdFile,
                reportTitle=reportTitle, selfContained = TRUE)
  setwd(cd)
}