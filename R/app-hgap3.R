###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodHGAP = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  sampleName = input$getNames()
  SMRT_Path = input$getFullPaths("Reads")
  SMRT_File = basename(input$getColumn("Reads"))
  ezSystem(paste("cp -r", SMRT_Path, "."))
  ezSystem(paste("mkdir", "smrt_input"))
  ezSystem(paste("tar -zxf", SMRT_File, "--strip-components=4 -C smrt_input"))
  readFile = file.path(getwd(), "smrt_input", "Analysis_Results", "*.subreads.fastq") 
  ezSystem(paste("ls", readFile, ">input.fofn"))
  library(xml2)
  inputXmlFile <-system.file(package = "ezRun", "extdata/hgap3_settings.xml")
  setting<-read_xml(inputXmlFile)
  new_size<-read_xml(paste0("<new><value>", param$genomeSize, "</value></new>"))
  xml_replace( xml_children(xml_children(xml_children(xml_children(setting)[5])[2])[2])[1], xml_children(new_size))
  new_coverage<-read_xml(paste0("<new><value>", param$xCoverage, "</value></new>"))
  xml_replace( xml_children(xml_children(xml_children(xml_children(setting)[5])[2])[5])[1], xml_children(new_coverage))
  write_xml(setting, "settings.xml")
  ezSystem("export SMRT=/misc/ngseq8/opt/smrtanalysis.2.3.0/install/smrtanalysis_2.3.0.140936")
  ezSystem("source $SMRT/etc/setup.sh")
  ezSystem("fofnToSmrtpipeInput.py input.fofn > input.xml")
  ezSystem("smrtpipe.py --params=settings.xml xml:input.xml")
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppHGAP <-
  setRefClass("EzAppHGAP",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodHGAP
                  name <<- "EzAppHGAP"
                  appDefaults <<- rbind(genomeSize = ezFrame(Type="character",  DefaultValue="5000000",  Description="The approximate genome size, in base pairs."),
                                        xCoverage = ezFrame(Type="integer",  DefaultValue="25",  Description="Fold coverage to target for when picking the minimum fragment length for assembly; typically 15 to 25."))
                }
              )
)
              
              
