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
    "-r 2.1.2",
    ## testing config
    # "-bg",
    # "-resume",
    # "-c ~/ezRun/test-nfcore/nf.config"
  )
  
  ezSystem(cmd)
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
        peakStyle  = ezFrame(Type="character", DefaultValue="broad", Description="Run MACS2 in broadPeak mode, otherwise in narrowPeak mode"),
        varStabilizationMethod = ezFrame(Type="character", DefaultValue="vst", Description="Use rlog transformation or vst (DESeq2)")
      )
    }
  )
)

##' @description get an nf-core/atacseq-formatted csv file
getSampleSheet <- function(input, param){
  if(any(input$getColumn(param$grouping) == "") || any(is.na(input$getColumn(param$grouping))))
    stop("No conditions detected. Please add them in the dataset before calling NfCoreAtacSeqApp.")
  
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
