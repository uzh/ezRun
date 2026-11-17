###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDelly <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  require("VariantAnnotation")
  genomeSeq = param$ezRef["refFastaFile"]
  sampleName = input$getNames()
  setwdNew(sampleName)

  bamFile <- input$getFullPaths("BAM")
  ovf = paste0(sampleName, ".vcf")
  ovzf = paste0(sampleName, ".vcf.gz")
  DellyLogFile = paste0(sampleName, "_delly.log")
  DellyStatsFile = paste0(sampleName, ".stats.txt")
  SurvivorLogFile = paste0(sampleName, ".SURVIVOR.log")
  cmd = "delly"
  if (param$ReadOpt == "sr") {
    cmd = paste(cmd, "call")
  } else if (param$ReadOpt == "pb") {
    cmd = paste(cmd, "lr -y pb")
  } else if (param$ReadOpt == "ont") {
    cmd = paste(cmd, "lr -y ont")
  } else {
    stop("ReadOpt must be sr, pb or ont, not: ", param$ReadOpt)
  }
  ## delly -t takes one type or a comma-separated list (e.g. DEL,INV)
  svTypes <- c("DEL", "INS", "DUP", "INV", "BND", "ALL")
  if (!all(strsplit(param$types, ",", fixed = TRUE)[[1]] %in% svTypes)) {
    stop(
      "types must be one of ",
      paste(svTypes, collapse = ", "),
      ", not: ",
      param$types
    )
  }
  cmd = paste(
    cmd,
    "-t",
    param$types,
    "-g",
    param$ezRef["refFastaFile"],
    param$cmdOptions,
    bamFile,
    ">",
    ovf,
    "2>",
    DellyLogFile
  )
  ezSystem(cmd)
  ezSystem(paste("bgzip -c", ovf, ">", ovzf))
  indexTabix(basename(ovzf), format = "vcf")

  SurvivorCmd = paste(
    "SURVIVOR stats",
    ovf,
    "-1 -1 -1",
    DellyStatsFile,
    "1>",
    SurvivorLogFile
  )
  ezSystem(SurvivorCmd)

  ##html file
  #setwd(start_path)
  htmlFile = output$getColumn("OutReport")
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c("fgcz.css", "Pbsv.Rmd", "fgcz_header.html", "banner.png")
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)
  params = list(sample = sampleName)
  rmarkdown::render(
    input = "Pbsv.Rmd",
    envir = new.env(),
    output_dir = ".",
    output_file = htmlFile,
    quiet = TRUE
  )

  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodPbsv(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getPbmm2Reference}}
EzAppDelly <-
  setRefClass(
    "EzAppDelly",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDelly
        name <<- "EzAppDelly"
        appDefaults <<- rbind(
          ReadOpt = ezFrame(
            Type = "character",
            DefaultValue = "sr",
            Description = "input read types: sr, pb, ont. Default is sr"
          ),
          types = ezFrame(
            Type = "character",
            DefaultValue = "ALL",
            Description = "SV type to call (delly -t): DEL, INS, DUP, INV, BND or ALL"
          )
        )
      }
    )
  )
