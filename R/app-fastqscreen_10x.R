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
