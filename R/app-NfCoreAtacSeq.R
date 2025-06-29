###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNfCoreAtacSeq <- function(input = NA, output = NA, param = NA) {
  sampleDataset = getSampleSheet(input, param)
  refbuild = param$refBuild
  outFolder = paste0(param$name, '_results')
  
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
    "--macs_gsize", getGenomeSize(param),
    if (param[['peakStyle']] == 'broad')  "" else "--narrow_peak",
    if (param[['varStabilizationMethod']] != 'vst') "--deseq2_vst false"  else "",
    ## configuration
    "-work-dir nfatacseq_work",
    "-profile apptainer",
    "-r 2.1.2"
  )

  ezSystem(cmd)

  cwd <- getwd()
  if(param[['runTwoGroupAnalysis']]){
    require(DESeq2)
    getData(output, param)

    makeRmdReport(
      output = output, param = param, selfContained = TRUE,
      rmdFile = "DiffPeak.Rmd", htmlFile = "DifferentialPeakAnalysisReport.html",
      reportTitle = 'Differential Peak Analysis'
    )
    #file.remove(list.files(pattern="^(dds|peak)\\.qs2$"))
    on.exit(setwd(cwd))
  }
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
getSampleSheet <- function(input, param){
  groups <- input$getColumn(param$grouping)
  if(any(groups == "") || any(is.na(groups)))
    stop("No conditions detected. Please add them in the dataset before calling NfCoreAtacSeqApp.")
  
  if(any(str_detect(groups, "[^a-zA-Z0-9]"))){
    separator <- table(unlist(str_extract_all(groups, "[^a-zA-Z0-9]"))) |> which.max() |> names()
    ngroups <- ncol(str_split(groups, separator, simplify = T))
    if(ngroups >2)
      stop('Values in the Condition column cannot be splitted in two groups for pairwise comparison. Please use a proper separator.')
  }
  
  oDir <- '.' ## param[['resultDir']]
  #if(!dir.exists(oDir)) dir.create(path = oDir)

  csvPath <- file.path(oDir, 'dataset.csv')
  
  input$meta |> 
    arrange(`Condition [Factor]`, `Read1 [File]`, `Read2 [File]`) |>
    rownames_to_column(var = 'SampleID [Factor]') |>
    group_by(`Condition [Factor]`) |>
    mutate(`Replicate [Factor]` = row_number()) |>
    ungroup() |>
    select('Condition [Factor]', 'Read1 [File]', 'Read2 [File]', 'Replicate [Factor]', 'SampleID [Factor]') |>
    ## the first 4 columns of the header must be: sample,fastq_1,fastq_2,replicate
    rename(sample    = 'Condition [Factor]', 
           fastq_1   = 'Read1 [File]', 
           fastq_2   = 'Read2 [File]', 
           replicate = 'Replicate [Factor]',
           sid       = 'SampleID [Factor]') |>
    mutate(fastq_1 = replace(fastq_1, sid %in% names(input$getFullPaths('Read1')), input$getFullPaths('Read1')[sid]),
           fastq_2 = replace(fastq_2, sid %in% names(input$getFullPaths('Read2')), input$getFullPaths('Read2')[sid])) |>
    write_csv(csvPath)
    return(csvPath)
}

##' @description estimate genome size
getGenomeSize <- function(param){
  fastaFile <- param$ezRef@refFastaFile
  gsize <- sum(as.numeric(fasta.seqlengths(fastaFile)))
  gsize <- round(gsize * 0.8)
  return(gsize)
}

getData <- function(output, param){
  nfCoreOutDir <- paste0(param$name, '_results', '/bwa/merged_replicate/macs2/', param$peakStyle, '_peak/consensus')
  
  dds <- readRDS(paste0(nfCoreOutDir, '/deseq2/consensus_peaks.mRp.clN.rds'))
  
  featureCounts <- vroom::vroom(paste0(nfCoreOutDir, '/consensus_peaks.mRp.clN.featureCounts.txt'), 
                                delim="\t", skip = 1, col_types = cols())
  
  peakAnno <- vroom::vroom(paste0(nfCoreOutDir, '/consensus_peaks.mRp.clN.annotatePeaks.txt'), delim="\t", col_types = cols()) |> 
    rename(c("PeakID"=1))
  
  rowData(dds) <- featureCounts[,1:5]
  dsgn <- data.frame(Condition = factor(paste(colData(dds)$Group1, colData(dds)$Group2, sep = "_")))
  rownames(dsgn) <- rownames(colData(dds))
  dds$Condition <- dsgn$Condition
  design(dds) <- ~ Condition
  
  outDir <- file.path(basename(output$getColumn('Result')),'diffpeak_analysis')
  if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
  setwdNew(outDir)
  
  qs2::qs_save(dds, file = 'dds.qs2')
  qs2::qs_save(peakAnno, file = 'peakAnno.qs2')
  # qs2::qs_save(list(dds = dds, peakAnno = peakAnno),  file = paste0(peakDir, '/results.qs2'))
  # return(list(dds=dds, peakAnno=peakAnno)) ## wo storing results, datasets need to be reloaded in DiffPeak.Rmd
}
