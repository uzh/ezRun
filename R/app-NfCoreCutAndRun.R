###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNfCoreCutAndRun <- function(input = NA, output = NA, param = NA) {
  sampleDataset = getSampleSheet(input, param)
  refbuild = param$refBuild
  outFolder = paste0(param$name, '_results')
  
  cmd = paste(
    "nextflow run nf-core/cutandrun",
     ## i/o
    "--input", sampleDataset,
    "--outdir", outFolder,
    ## genome files
    "--fasta", param$ezRef@refFastaFile,
    "--gtf", param$ezRef@refFeatureFile,
    "--gene_bed", str_replace(param$ezRef@refAnnotationFile,
                              basename(param$ezRef@refAnnotationFile),
                              'genes.bed'),
    if (param[['blacklist']] == "")  "" else paste0("--blacklist ", param[['blacklist']]),
    ## parameters
    "--dt_calc_all_matrix false",
    "--macs_gsize", getGenomeSize(param),
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
        blacklist = ezFrame(Type="character", DefaultValue="", Description="Path to genome blacklist"),
        peakCaller = ezFrame(Type="character", DefaultValue="macs2", Description="Select the peak caller for the pipeline"),
        spikeinGenome = ezFrame(Type="character", DefaultValue="macs2", Description="Select the reference for the spike-in genome"),
        normalization = ezFrame(Type="character", DefaultValue="macs2", Description="Select the target read normalization mode"),
        peakStyle  = ezFrame(Type="character", DefaultValue="broad", Description="Run MACS2 in broadPeak mode, otherwise in narrowPeak mode")
      )
    }
  )
)

##' @description get an nf-core/cutandrun-formatted csv file
getSampleSheet <- function(input, param){
  oDir <- '.'
  csvPath <- file.path(oDir, 'dataset.csv')
  
  input$meta |> 
    arrange(`Condition [Factor]`, `Read1 [File]`, `Read2 [File]`) |>
    rownames_to_column(var = 'SampleID [Factor]') |>
    group_by(`Condition [Factor]`) |>
    mutate(`Replicate [Factor]` = row_number()) |>
    ungroup() |>
    select('Condition [Factor]', 'Replicate [Factor]', 'Read1 [File]', 'Read2 [File]', 'Control [Factor]', 'SampleID [Factor]') |>
    ## header must have the following columns and in this order: group, replicate, fastq_1, fastq_2, control
    rename(group     = 'Condition [Factor]', 
           replicate = 'Replicate [Factor]',
           fastq_1   = 'Read1 [File]',
           fastq_2   = 'Read2 [File]', 
           control   = 'Control [Factor]',
           sid       = 'SampleID [Factor]') |>
    mutate(fastq_1 = replace(fastq_1, sid %in% names(input$getFullPaths('Read1')), input$getFullPaths('Read1')[sid]),
           fastq_2 = replace(fastq_2, sid %in% names(input$getFullPaths('Read2')), input$getFullPaths('Read2')[sid])) |>
    select(-sid) |> ## max 5 columns allowed ...
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
