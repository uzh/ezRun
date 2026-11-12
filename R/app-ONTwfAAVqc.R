# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodONTwfAAVqc <- function(input = NA, output = NA, param = NA) {
  sampledataset = input$getFullPaths("Read1")
  samplename = input$getNames()
  cores = param$cores
  refHelper = param$refHelper
  refRepCap = param$refRepCap
  refTrans = param$refTrans
  itr1Start = param$itr1Start
  itr1End = param$itr1End
  itr2Start = param$itr2Start
  itr2End = param$itr2End
  refHost = file.path(
    dirname(dirname(file.path('/srv/GT/reference/', param$refBuild))),
    'Sequence/WholeGenomeFasta/genome.fa'
  )
  opt = param$cmdOptions
  cmd = paste(
    "nextflow run /srv/GT/software/epi2me-labs/wf-aav-qc/",
    "-w",
    paste0(samplename, "_workspace"),
    "--threads",
    cores,
    "--fastq",
    sampledataset,
    "--itr1_start",
    itr1Start,
    "--itr1_end",
    itr1End,
    "--itr2_start",
    itr2Start,
    "--itr2_end",
    itr2End,
    "--ref_helper",
    refHelper,
    "--ref_host",
    refHost,
    "--ref_rep_cap",
    refRepCap,
    "--ref_transgene_plasmid",
    refTrans,
    opt,
    "--out_dir",
    samplename,
    "-profile singularity"
  )
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodONTwfAAVqc()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppONTwfAAVqc <-
  setRefClass(
    "EzAppONTwfAAVqc",
    contains = "EzApp",
    methods = list(
      ## All unconditional: ezMethodONTwfAAVqc always runs the whole EPI2ME wf-aav-qc Nextflow
      ## workflow. Its tools that define the result: minimap2 (mapping to the combined
      ## reference), samtools depth (ITR-ITR coverage) and bcftools consensus (one SAMtools/
      ## BCFtools paper), medaka (transgene plasmid variants and consensus; no paper).
      ## fastcat ingress and seqkit bam (alignment summary tables) are helpers: not cited.
      citation = function(param = list()) {
        c(
          "wf-aav-qc. Oxford Nanopore Technologies. https://github.com/epi2me-labs/wf-aav-qc",
          "Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35(4), 316-319 (2017). https://doi.org/10.1038/nbt.3820",
          "Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34(18), 3094-3100 (2018). https://doi.org/10.1093/bioinformatics/bty191",
          "medaka. Oxford Nanopore Technologies. https://github.com/nanoporetech/medaka",
          "Danecek, P. et al. Twelve years of SAMtools and BCFtools. GigaScience 10(2), giab008 (2021). https://doi.org/10.1093/gigascience/giab008"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodONTwfAAVqc
        name <<- "EzAppONTwfAAVqc"
        appDefaults <<- rbind(
          itr1Start = ezFrame(
            Type = "integer",
            DefaultValue = "11",
            Description = "The start position of ITR1."
          ),
          itr1End = ezFrame(
            Type = "integer",
            DefaultValue = "156",
            Description = "The end position of ITR1."
          ),
          itr2Start = ezFrame(
            Type = "integer",
            DefaultValue = "2156",
            Description = "The start position of ITR2."
          ),
          itr2End = ezFrame(
            Type = "integer",
            DefaultValue = "2286",
            Description = "The end position of ITR2."
          )
        )
      }
    )
  )
