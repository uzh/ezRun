# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodONTwfSc <- function(input = NA, output = NA, param = NA) {
  sampledataset = input$getFullPaths("Read1")
  samplename = input$getNames()
  refbuild = param$refBuild
  kitname = param$kitname
  kitversion = param$kitVersion
  expCells = param$expCells

  cmd = paste(
    "nextflow run /srv/GT/software/epi2me-labs/wf-single-cell/",
    "-w",
    paste0(samplename, "/workspace"),
    "--fastq",
    sampledataset,
    "--kit_name",
    kitname,
    "--kit_version",
    kitversion,
    "--expected_cells",
    expCells,
    "--ref_genome_dir",
    paste0("/srv/GT/reference/ont-wf-single-cell/", refbuild),
    "--out_dir",
    samplename,
    "-profile singularity",
    "--plot_umaps"
  )
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodONTwfSc()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppONTwfSc <-
  setRefClass(
    "EzAppONTwfSc",
    contains = "EzApp",
    methods = list(
      ## All unconditional: ezMethodONTwfSc always runs the whole EPI2ME wf-single-cell Nextflow
      ## workflow. Its tools that define the result: vsearch (adapter scan and read orientation),
      ## minimap2 (genome and transcriptome alignment), StringTie2 long-read mode (-L, per-sample
      ## transcriptome), gffcompare/gffread (one paper), UMI-tools directional clustering (UMI
      ## deduplication). seqkit, samtools, parasail/rapidfuzz libraries and the UMAP plots are
      ## helpers or report content: not cited.
      citation = function(param = list()) {
        c(
          "wf-single-cell. Oxford Nanopore Technologies. https://github.com/epi2me-labs/wf-single-cell",
          "Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35(4), 316-319 (2017). https://doi.org/10.1038/nbt.3820",
          "Rognes, T., Flouri, T., Nichols, B., Quince, C. & Mahé, F. VSEARCH: a versatile open source tool for metagenomics. PeerJ 4, e2584 (2016). https://doi.org/10.7717/peerj.2584",
          "Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34(18), 3094-3100 (2018). https://doi.org/10.1093/bioinformatics/bty191",
          "Kovaka, S., Zimin, A.V., Pertea, G.M., Razaghi, R., Salzberg, S.L. & Pertea, M. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biology 20, 278 (2019). https://doi.org/10.1186/s13059-019-1910-1",
          "Pertea, G. & Pertea, M. GFF Utilities: GffRead and GffCompare. F1000Research 9, 304 (2020). https://doi.org/10.12688/f1000research.23297.1",
          "Smith, T., Heger, A. & Sudbery, I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research 27(3), 491-499 (2017). https://doi.org/10.1101/gr.209601.116"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodONTwfSc
        name <<- "EzAppONTwfSc"
        appDefaults <<- rbind(
          expected_cells = ezFrame(
            Type = "integer",
            DefaultValue = "500",
            Description = "number of expected cells"
          ),
          kitname = ezFrame(
            Type = "character",
            DefaultValue = "3prime",
            Description = "10x kit name"
          )
        )
      }
    )
  )
