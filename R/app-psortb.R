###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodPsortb = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  require(data.table)
  require(readr)
  require(dplyr)
  require(reshape2)
  require(tidyverse)

  opt = param$cmdOptions
  org = param$org
  sampleName = input$getNames()
  proteins = input$getFullPaths("Proteins")
  ezSystem(paste("mkdir", sampleName))
  cmd = paste(
    "/usr/local/ngseq/src/psortb/psortb_app -i",
    proteins,
    "-r",
    sampleName,
    org,
    "--output terse",
    "-s /usr/local/ngseq/src/psortb/psortb.sif",
    opt,
    "1>",
    paste0(sampleName, "_psortb.log")
  )
  ezSystem(cmd)
  wddir <- "."
  outfile <- list.files(
    paste0(wddir, "/", sampleName),
    pattern = "_psortb_.*\\.txt"
  )
  outfile <- file.path(wddir, sampleName, outfile)
  ezSystem(paste("cp", outfile, basename(output$getColumn("PsortbOut"))))
  df = read_delim(outfile[1])
  list_of_OM_Ids = df[df$Localization == "OuterMembrane", c(1), drop = FALSE]
  write.table(
    sapply(str_split(list_of_OM_Ids$SeqID, " "), .subset, 1),
    paste0(wddir, "/", sampleName, "psortb.om.txt"),
    row.names = FALSE,
    quote = FALSE
  )
  cmd = paste(
    "/usr/local/ngseq/packages/Tools/seqtk/1.4/bin/bin subseq",
    proteins,
    paste0(wddir, "/", sampleName, "psortb.om.txt"),
    ">",
    basename(output$getColumn("Proteins"))
  )
  ezSystem(cmd)
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodPsortb
##' @description Use this reference class to run
EzAppPsortb <-
  setRefClass(
    "EzAppPsortb",
    contains = "EzApp",
    methods = list(
      ## PSORTb 3.0: always runs (psortb_app call in ezMethodPsortb). seqtk only subsets the
      ## protein FASTA and is not cited. Kept identical to the EzAppPsortb copy in app-signalp.R.
      citation = function(param = list()) {
        c(
          "Yu, N.Y. et al. PSORTb 3.0: improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes. Bioinformatics 26(13), 1608-1615 (2010). https://doi.org/10.1093/bioinformatics/btq249"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodPsortb
        name <<- "EzAppPsortb"
        appDefaults <<- rbind(
          org = ezFrame(
            Type = "character",
            DefaultValue = "--negative",
            Description = "type of organism: gram negative/ gram positive bacteria or archea"
          )
        )
      }
    )
  )
