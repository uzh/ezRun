###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNfCoreCutAndRun <- function(input = NA, output = NA, param = NA) {
  refbuild = param$refBuild
  outFolder = output$getColumn("Result") |> basename()

  fullGenomeSize <- param$ezRef@refFastaFile %>% Rsamtools::FaFile() %>% GenomeInfoDb::seqlengths() %>% sum()
  effectiveGenomeSize <- (fullGenomeSize * 0.8 ) %>% round()

  nfSampleFile <- file.path('dataset.csv')
  nfSampleInfo = getCutAndRunSampleSheet(input, param)
  write_csv(nfSampleInfo, nfSampleFile)


  blackListFile <- getBlackListFile(input, param)
  setNFCacheDir()
  cmd = paste(
    "nextflow run nf-core/cutandrun",
     ## i/o
    "--input", nfSampleInfo,
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
    "-r 3.2.2"
  )
  ezSystem(cmd)

  writeAtacIgvSession(param, outFolder, jsonFileName = paste0("/scratch/mderrico/TestNfCoreCutNRun", "/igv_session.json"), bigwigRelPath = "/04_reporting/igv/",
                      baseUrl = file.path(PROJECT_BASE_URL, output$getColumn("Result")))
  
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
        keepBams = ezFrame(Type="logical", DefaultValue = FALSE, Description= "Should bam files be stored")
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
}