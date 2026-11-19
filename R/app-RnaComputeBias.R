###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppRnaComputeBias <-
  setRefClass(
    "EzAppRnaComputeBias",
    contains = "EzApp",
    methods = list(
      ## kallisto and fastp unconditional: ezComputeBias (R/qc-functions.R) runs EzAppKallisto, which trims
      ## with fastp, on every sample, and the GC and length bias statistics come from those counts. STAR
      ## is not cited: its genome alignments feed only the gene-body coverage panel, which coverage() of a
      ## genome BAM computes per chromosome, not per transcript. Picard, RSeQC and dupRadar (inside the
      ## STAR step) and samtools are not cited either.
      citation = function(param = list()) {
        c(
          "Bray, N.L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525-527 (2016). https://doi.org/10.1038/nbt.3519",
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodRnaComputeBias
        name <<- "EzAppRnaComputeBias"
        appDefaults <<- rbind(
          minReadsPerSample = ezFrame(
            Type = "numeric",
            DefaultValue = 30000,
            Description = ""
          ),
          maxReadsPerSample = ezFrame(
            Type = "numeric",
            DefaultValue = 1e6,
            Description = ""
          ),
          minReadsPerGene = ezFrame(
            Type = "numeric",
            DefaultValue = 3,
            Description = ""
          ),
          minPresentFraction = ezFrame(
            Type = "numeric",
            DefaultValue = 0.2,
            Description = ""
          )
        )
      }
    )
  )


ezMethodRnaComputeBias <- function(input, output, param) {
  system('export HDF5_DISABLE_VERSION_CHECK=1 kallisto')
  ezWrite.table(input$meta, 'input_dataset.tsv')
  dsFile <- 'input_dataset.tsv'
  qcSummaryDir = getwd()
  dsName = 'report-RNAseq'
  dir.create(dsName)
  plateId <- unique(input$getColumn("PlateName"))
  orderId <- unique(input$getColumn("Order Id"))
  tag <- paste0("Plate_", plateId, '_o', orderId)
  ezComputeBias(
    dsFile = dsFile,
    dsName = dsName,
    param = param,
    qcSummaryDir = qcSummaryDir,
    minReadsPerSample = param$minReadsPerSample,
    maxReadsPerSample = param$maxReadsPerSample,
    minReadsPerGene = param$minReadsPerGene,
    minPresentFraction = param$minPresentFraction,
    tag = tag
  )
  return('success')
}
