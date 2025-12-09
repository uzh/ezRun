###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNfCoreAtacSeq <- function(input = NA, output = NA, param = NA) {
  refbuild = param$refBuild
  outFolder = output$getColumn("ATAC_Result") |> basename()
  
  fullGenomeSize <- param$ezRef@refFastaFile %>% Rsamtools::FaFile() %>% GenomeInfoDb::seqlengths() %>% sum()
  effectiveGenomeSize <- (fullGenomeSize * 0.8 ) %>% round()

  nfSampleFile <- file.path('dataset.csv')
  nfSampleInfo = getAtacSampleSheet(input, param)
  write_csv(nfSampleInfo, nfSampleFile)
  prepNFCoreEnv()
  configFile <- writeNextflowLimits(param)
  cmd = paste(
    "nextflow run nf-core/atacseq",
     ## i/o
    "--input", nfSampleFile,
    "--outdir", outFolder,
    ## genome files
    "--fasta", param$ezRef@refFastaFile,
    "--gtf", param$ezRef@refFeatureFile,
    "--gene_bed", str_replace(param$ezRef@refAnnotationFile,
                              basename(param$ezRef@refAnnotationFile),
                              'genes.bed'),
    ## parameters
    "--macs_gsize", effectiveGenomeSize,
    if (param[['peakStyle']] == 'broad')  "" else "--narrow_peak",
    if (param[['varStabilizationMethod']] != 'vst') "--deseq2_vst false"  else "",
    ## configuration
    "-work-dir nfatacseq_work",
    "-profile apptainer",
    "-r", param$pipelineVersion,
    "-c", configFile,
    param$cmdOptions #,
    # "-resume"  ## for testing
  )
  ezSystem(cmd)
  ezSystem(paste('mv', configFile, outFolder))
  ## multiple fastq files per library have been merged by the processing (if any)
  ## now we work with the library names and reduce the dataset
  nfSampleInfo$libName <- paste0(nfSampleInfo$sample, "_REP", nfSampleInfo$replicate)
  nfSampleInfo <- nfSampleInfo[!duplicated(nfSampleInfo$sid), ]
  
  sampleCountFiles <- writePerSampleCountPeaksFiles(nfSampleInfo, countDir=paste0(outFolder, "/bwa/merged_library/macs2/", param$peakStyle, "_peak/consensus/"))
  
  jsonFile <- writeAtacIgvSession(param, outFolder, jsonFileName = paste0(outFolder, "/igv_session.json"), bigwigRelPath = "/bwa/merged_library/bigwig/",
                      baseUrl = file.path(PROJECT_BASE_URL, output$getColumn("ATAC_Result")))
  writeNfCoreIgvHtml(param, jsonFile, title = "NfCoreAtacSeq MultiSample Coverage Tracks", htmlTemplate = "templates/igvNfCoreTemplate.html", htmlFileName = paste0(outFolder, "/igv_session.html"))
  

  if(param[['runTwoGroupAnalysis']]){
    library(DESeq2)
    nfCoreOutDir <- paste0(param$name, '_results', '/bwa/merged_replicate/macs2/', param$peakStyle, '_peak/consensus')
    peakAnno <- vroom::vroom(paste0(nfCoreOutDir, '/consensus_peaks.mRp.clN.annotatePeaks.txt'), delim="\t", col_types = cols()) %>%
      rename(c("PeakID"=1))
    
    grouping <- input$getColumn(param$grouping)
    dds <- getDdsFromConcensusPeaks(output, param, grouping)
    outDir <- file.path(basename(output$getColumn('Result')), 'diffpeak_analysis')
    cd = getwd()
    setwdNew(outDir)
    makeRmdReport(
      output = output, param = param, peakAnno=peakAnno, dds=dds, selfContained = TRUE,
      rmdFile = "DiffPeak.Rmd", htmlFile = "DifferentialPeakAnalysisReport.html",
      reportTitle = 'Differential Peak Analysis', use.qs2 = TRUE
    )
    setwd(cd)
  }

  dirsToRemove <- c("genome", "trimgalore", "fastqc", "igv")
  if(ezIsSpecified(param$keepBams)){
    keepBams <- param$keepBams
  } else {
    keepBams <- TRUE
  }
  cleanupAtacOutFolder(outFolder, dirsToRemove, keepBams)

  return("Success")
}

EzAppNfCoreAtacSeq <- setRefClass(
  "EzAppNfCoreAtacSeq",
  contains = "EzApp",
  methods = list(
    initialize = function()
    {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodNfCoreAtacSeq
      name <<- "EzAppNfCoreAtacSeq"
      ## minimum nf-core parameters
      appDefaults <<- rbind(
        runTwoGroupAnalysis = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "Run two group analysis"),
        peakStyle  = ezFrame(Type="character", DefaultValue="broad", Description="Run MACS2 in broadPeak mode, otherwise in narrowPeak mode"),
        varStabilizationMethod = ezFrame(Type="character", DefaultValue="vst", Description="Use rlog transformation or vst (DESeq2)"),
        keepBams = ezFrame(Type="logical", DefaultValue = FALSE, Description= "Should bam files be stored"),
        pipelineVersion = ezFrame(Type="character", DefaultValue = '2.1.2', Description= "specify pipeline version")
      )
    }
  )
)

##' @description get an nf-core/atacseq-formatted csv file
getAtacSampleSheet <- function(input, param){
  groups <- input$getColumn(param$grouping)
  if(any(groups == "") || any(is.na(groups)))
    stop("No conditions detected. Please add them in the dataset before calling NfCoreAtacSeqApp.")


  listFastq1 <- input$getFullPathsList("Read1")
  listFastq2 <- input$getFullPathsList("Read2")
  
  nfSampleInfo <- ezFrame(
    sample = rep(input$getColumn(param$grouping), lengths(listFastq1)),
    fastq_1 = unlist(listFastq1),
    fastq_2 = unlist(listFastq2),
    replicate = rep(ezReplicateNumber(input$getColumn(param$grouping)), lengths(listFastq1)),
    sid = rep(input$getNames(), lengths(listFastq1))
  )
  
  # input$meta |> 
  #   arrange(`Condition [Factor]`, `Read1 [File]`, `Read2 [File]`) |>
  #   rownames_to_column(var = 'SampleID [Factor]') |>
  #   group_by(`Condition [Factor]`) |>
  #   mutate(`Replicate [Factor]` = row_number()) |>
  #   ungroup() |>
  #   select('Condition [Factor]', 'Read1 [File]', 'Read2 [File]', 'Replicate [Factor]', 'SampleID [Factor]') |>
  #   ## the first 4 columns of the header must be: sample,fastq_1,fastq_2,replicate
  #   rename(sample    = 'Condition [Factor]', 
  #          fastq_1   = 'Read1 [File]', 
  #          fastq_2   = 'Read2 [File]', 
  #          replicate = 'Replicate [Factor]',
  #          sid       = 'SampleID [Factor]') |>
  #   mutate(fastq_1 = replace(fastq_1, sid %in% names(input$getFullPaths('Read1')), input$getFullPaths('Read1')[sid]),
  #          fastq_2 = replace(fastq_2, sid %in% names(input$getFullPaths('Read2')), input$getFullPaths('Read2')[sid])) |>
  #   write_csv(csvPath)

  return(nfSampleInfo)
}

writePerSampleCountPeaksFiles <- function(nfSampleInfo, countDir="."){
  libColumnNames <- paste0(nfSampleInfo$libName, ".mLb.clN.sorted.bam")
  sampleNames <- nfSampleInfo$sid
  sampleCountFiles <- paste0(countDir, "/", sampleNames, ".txt")
  annoColumnNames <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  x <- data.table::fread(file.path(countDir, "consensus_peaks.mLb.clN.featureCounts.txt"), data.table=FALSE)
  for (i in 1:nrow(nfSampleInfo)){
    xSel <- x[ , c(annoColumnNames, libColumnNames[i])] |> dplyr::rename("matchCounts" := !!libColumnNames[i])
    ezWrite.table(xSel, file=sampleCountFiles[i], row.names = FALSE)
  }
  return(sampleCountFiles)
}



writeHtmlWrapper <- function(htmlFile, igvAppLink){
  library(htmltools)
  
  page <- htmlTemplate(
    text_ = "
  <!DOCTYPE html>
  <html>
    <head><title>{{title}}</title></head>
    <body>
      <h1>{{header}}</h1>
      <a href={{igvAppLink}}>{{igvAppLink}}</a>
    </body>
  </html>
  ",
    title = "IGV Starter",
    header = "IGV Starter Link",
    igvAppLink = igvAppLink
  )
  
  # Save to file
  save_html(page, htmlFile)
  
}


getDdsFromConcensusPeaks <- function(output, param, grouping){
  nfCoreOutDir <- paste0(param$name, '_results', '/bwa/merged_replicate/macs2/', 
                         param$peakStyle, '_peak/consensus')
  
  dds <- readRDS(paste0(nfCoreOutDir, '/deseq2/consensus_peaks.mRp.clN.rds'))
  
  featureCounts <- vroom::vroom(paste0(nfCoreOutDir, '/consensus_peaks.mRp.clN.featureCounts.txt'), 
                                delim="\t", comment="#", col_types = cols())
  
  rowData(dds) <- featureCounts[, c("Chr", "Start", "End", "Strand", "Length")]
  ## samples and grouping must be in the same order
  samples <- colData(dds)$sample %>% str_remove(., '_REP\\d+')
  grouping <- grouping[match(samples, grouping)]
  dds$Condition <- grouping %>% as.factor()
  design(dds) <- ~ Condition

  return(dds)  
}

##' @description clean up NfCoreAtacSeq_result directory
cleanupAtacOutFolder <- function(outFolder, dirsToRemove, keepBams=TRUE){
  if(!keepBams){
    bamPath <- paste0(outFolder,"/bwa/")
    bamsToDelete <- dir(path=bamPath, pattern="*.bam(.bai)?$", recursive=TRUE)
    file.remove(file.path(bamPath, bamsToDelete))
    cat("Deleted .bam and .bam.bai files from the bwa directory.\n")
  }
  absolutePaths <- paste(outFolder, dirsToRemove, sep="/")
  unlink(absolutePaths, recursive=TRUE)
  cat(paste0("Deleted subdirectory: ",dirsToRemove, "\n"))
}

##' @description write IGV session in json format
writeAtacIgvSession <- function(param, outFolder, jsonFileName, bigwigRelPath, baseUrl){
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
  tracks[[i+2]] <- list(
      id = "genes",
      url = file.path(REF_HOST, param$ezRef@refBuild,'Genes/transcripts.only.gtf'),
      format =	"gtf",
      type = "annotation",
      name = "genes")
  
  tracks[[i+3]] <- list(
      id = "exons",
      url = file.path(REF_HOST, param$ezRef@refBuild,'Genes/genes.bed'),
      format =	"bed",
      type = "annotation",
      name = "exons")
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
