###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodSubreads = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  sampleName = input$getNames()
  SMRT_Path = input$getFullPaths("Reads")
  SMRT_File = basename(input$getColumn("Reads"))
  ezSystem(paste("cp -r", SMRT_Path, "."))
  ezSystem(paste("mkdir", "smrt_input"))
  ezSystem(paste("tar -zxf", SMRT_File, "--strip-components=4 -C smrt_input"))
  readFile = file.path(getwd(), "smrt_input", "Analysis_Results", "*.bax.h5")
  ezSystem(paste("mkdir", sampleName))
  start_path = getwd()
  workdir = file.path(getwd(), sampleName)
  setwd(workdir)
  ezSystem(paste("ls", readFile, ">input.fofn"))
  require(xml2)
  inputXmlFile <- system.file(
    package = "ezRun",
    "extdata/subreads_settings.xml"
  )
  setting <- read_xml(inputXmlFile)
  if (param$minSubReadLength != 50) {
    new_node <- read_xml(paste0(
      "<new><value>",
      param$minSubReadLength,
      "</value></new>"
    ))
    xml_replace(
      xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[2])[
        1
      ],
      xml_children(new_node)
    )
  }
  if (param$readScore != 75) {
    new_node <- read_xml(paste0(
      "<new><value>",
      param$readScore,
      "</value></new>"
    ))
    xml_replace(
      xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[3])[
        1
      ],
      xml_children(new_node)
    )
  }
  if (param$minLength != 50) {
    new_node <- read_xml(paste0(
      "<new><value>",
      param$minLength,
      "</value></new>"
    ))
    xml_replace(
      xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[4])[
        1
      ],
      xml_children(new_node)
    )
  }
  write_xml(setting, "settings.xml")
  ezSystem(SMRT_CMD)
  ezSystem(paste("cp", "data/filtered_subreads.fastq", start_path))
  setwd(start_path)
  cmd = "gzip filtered_subreads.fastq"
  ezSystem(cmd)
  ezSystem(paste(
    "mv",
    "filtered_subreads.fastq.gz",
    basename(output$getColumn("Read1"))
  ))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppSubreads <-
  setRefClass(
    "EzAppSubreads",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSubreads
        name <<- "EzAppSubreads"
        appDefaults <<- rbind(
          minSubReadLength = ezFrame(
            Type = "integer",
            DefaultValue = "50",
            Description = "Minimum Subread Length. Subreads shorter than this value (in bps) are filtered out and excluded from pre-assembly."
          ),
          readScore = ezFrame(
            Type = "integer",
            DefaultValue = "75",
            Description = "Minimum Polymerase Read Quality. Polymerase reads with lower quality than this value are filtered out and excluded from pre-assembly."
          ),
          minLength = ezFrame(
            Type = "integer",
            DefaultValue = "50",
            Description = "Minimum Polymerase Read Length. Polymerase reads shorter than this value (in bps) are filtered out and excluded from pre-assembly."
          )
        )
      }
    )
  )
