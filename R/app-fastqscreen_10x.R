###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFastqScreen_10x <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  dataset <- input$meta
  sampleDirs <- input$getFullPaths("RawDataDir")
  stopifnot(all(grepl("\\.tar$", sampleDirs)))

  taredfiles <- lapply(
    sampleDirs,
    untar,
    list = TRUE,
    tar = system("which tar", intern = TRUE)
  )
  if (any(str_detect(unlist(taredfiles), "_R3_"))) {
    ## ATAC has R3 for real data
    taredfiles_R2 <- sapply(
      taredfiles,
      function(x) {
        grep("_R3_", x, value = TRUE) %>% head(1)
      }
    )
  } else if (any(str_detect(unlist(taredfiles), "_R2_"))) {
    ## RNA has R2 for real data
    taredfiles_R2 <- sapply(
      taredfiles,
      function(x) {
        grep("_R2_", x, value = TRUE) %>% head(1)
      }
    )
  }
  for (i in 1:length(sampleDirs)) {
    untar(
      sampleDirs[i],
      files = taredfiles_R2[i],
      tar = system("which tar", intern = TRUE)
    )
  }
  taredfiles_R2 <- normalizePath(taredfiles_R2)
  dataset$`Read1` <- taredfiles_R2
  input <- EzDataset(meta = dataset, dataRoot = NULL)
  ezMethodFastqScreen(input = input, output = output, param = param)

  return("Success")
}

EzAppFastqScreen_10x <-
  setRefClass(
    "EzAppFastqScreen_10x",
    contains = "EzApp",
    methods = list(
      ## Thin wrapper: unpacks the 10x tar bundle then delegates 100% to
      ## ezMethodFastqScreen -- identical citations to the non-10x FastqScreenApp.
      citation = function() {
        c(
          "Wingett, S.W. & Andrews, S. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research 7, 1338 (2018). https://doi.org/10.12688/f1000research.15931.2",
          "Langmead, B. & Salzberg, S.L. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357-359 (2012). https://doi.org/10.1038/nmeth.1923",
          "Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0",
          "Morgan, M., Anders, S., Lawrence, M., Aboyoun, P., Pagès, H. & Gentleman, R. ShortRead: a bioconductor package for input, quality assessment and exploration of high-throughput sequence data. Bioinformatics 25(19), 2607-2608 (2009). https://doi.org/10.1093/bioinformatics/btp450",
          "Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient manipulation of biological strings. R package version 2.80.1. https://doi.org/10.18129/B9.bioc.Biostrings",
          "Bembom, O. & Ivanek, R. seqLogo: Sequence logos for DNA sequence alignments. R package version 1.78.0. https://doi.org/10.18129/B9.bioc.seqLogo",
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodFastqScreen_10x
        name <<- "EzAppFastqScreen_10x"
        appDefaults <<- rbind(
          nTopSpecies = ezFrame(
            Type = "integer",
            DefaultValue = 10,
            Description = "number of species to show in the plots"
          ),
          confFile = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "the configuration file for fastq screen"
          ),
          virusCheck = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "check for viruses in unmapped data"
          ),
          minAlignmentScore = ezFrame(
            Type = "integer",
            DefaultValue = "-20",
            Description = "the min alignment score for bowtie2"
          ),
          trimAdapter = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "whether to search for the adapters and trim them"
          ),
          readFileToUse = ezFrame(
            Type = "character",
            DefaultValue = "Read1",
            Description = "which read file(s) to use for the analysis: Read1, Read2, or both"
          )
        )
      }
    )
  )
