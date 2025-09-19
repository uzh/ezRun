###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNfCoreAtacSeq <- function(input = NA, output = NA, param = NA) {
  sampleDataset = getAtacSampleSheet(input, param)
  refbuild = param$refBuild
  outFolder = paste0(param$name, '_results')
  
  fullGenomeSize <- param$ezRef@refFastaFile %>% Rsamtools::FaFile() %>% GenomeInfoDb::seqlengths() %>% sum()
  effectiveGenomeSize <- (fullGenomeSize * 0.8 ) %>% round()

  cmd = paste(
    "nextflow run nf-core/atacseq",
     ## i/o
    "--input", sampleDataset,
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
    "-r 2.1.2" #,
    # "-resume"  ## for testing
  )

  ezSystem(cmd)

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

  dirsToRemove <- c("genome", "trimgalore", "fastqc")
  if(ezIsSpecified(param$keepBams)){
    keepBams <- param$keepBams
  } else {
    keepBams <- TRUE
  }
  cleanupOutFolder(outFolder, dirsToRemove, keepBams)

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
        varStabilizationMethod = ezFrame(Type="character", DefaultValue="vst", Description="Use rlog transformation or vst (DESeq2)")
      )
    }
  )
)

##' @description get an nf-core/atacseq-formatted csv file
getAtacSampleSheet <- function(input, param){
  groups <- input$getColumn(param$grouping)
  if(any(groups == "") || any(is.na(groups)))
    stop("No conditions detected. Please add them in the dataset before calling NfCoreAtacSeqApp.")

  # oDir <- '.' ## param[['resultDir']]
  #if(!dir.exists(oDir)) dir.create(path = oDir)

  csvPath <- file.path('dataset.csv')

  ## TODO: does not yet support comma-separated file paths
  nfSampleInfo <- ezFrame(
    sample = input$getColumn(param$grouping),
    fastq_1 = input$getFullPaths("Read1"),
    fastq_2 = input$getFullPaths("Read2"),
    replicate = ezReplicateNumber(input$getColumn(param$grouping)),
    sid = input$getNames()
  )
  write_csv(nfSampleInfo, csvPath)
  
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

  return(csvPath)
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

cleanupOutFolder <- function(outFolder, dirsToRemove, keepBams=TRUE){
  if(!keepBams){
    bamPath <- paste0(outFolder,"/bwa/merged_library/")
    bamsToDelete <- dir(path=bamPath, pattern="*.bam*")
    file.remove(file.path(bamPath, bamsToDelete))
    cat("Deleted bam and bam.bai files form bwa directory.\n")
  }
  absolutePaths <- paste(outFolder, dirsToRemove, sep="/")
  unlink(absolutePaths, recursive=TRUE)
  cat(paste0("Deleted subdirectory: ",dirsToRemove, "\n"))
}
