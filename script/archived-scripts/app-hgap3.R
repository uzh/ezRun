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
  readFile = file.path(getwd(), "smrt_input", "Analysis_Results", "*.bax.h5") 
  ezSystem(paste("mkdir", sampleName))
  start_path = getwd()
  workdir = file.path(getwd(), sampleName)
  setwd(workdir)
  ezSystem(paste("ls", readFile, ">input.fofn"))
  require(xml2)
  inputXmlFile <-system.file(package = "ezRun", "extdata/hgap3_settings.xml")
  setting<-read_xml(inputXmlFile)
  if (param$minSubReadLength != 500){
	new_node<-read_xml(paste0("<new><value>", param$minSubReadLength, "</value></new>"))
	xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[2])[1], xml_children(new_node))
  }
  if(param$readScore != 0.8){
	new_node<-read_xml(paste0("<new><value>", param$readScore, "</value></new>"))
	xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[3])[1], xml_children(new_node))
  }
  if(param$minLength != 100){
        new_node<-read_xml(paste0("<new><value>", param$minLength, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[4])[1], xml_children(new_node))
  }
  if (param$genomeSize != 5000000){
  	new_node<-read_xml(paste0("<new><value>", param$genomeSize, "</value></new>"))
  	xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[2])[2])[1], xml_children(new_node))
  }
  if (param$xCoverage != 25){
  	new_coverage<-read_xml(paste0("<new><value>", param$xCoverage, "</value></new>"))
  	xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[2])[5])[1], xml_children(new_node))
  }
  if(param$targetChunks != 6){
        new_node<-read_xml(paste0("<new><value>", param$targetChunks, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[1])[4])[1], xml_children(new_node))
  }
  if(param$splitBestn != 10){
        new_node<-read_xml(paste0("<new><value>", param$splitBestn, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[1])[5])[1], xml_children(new_node))
  }
  if(param$totalBestn != 24){
        new_node<-read_xml(paste0("<new><value>", param$totalBestn, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[1])[6])[1], xml_children(new_node))
  }
  if(param$minCorCov != 6){
        new_node<-read_xml(paste0("<new><value>", param$minCorCov, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[1])[7])[1], xml_children(new_node))
  }
  if(param$ovlErrorRate != 0.06){
        new_node<-read_xml(paste0("<new><value>", param$ovlErrorRate, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[2])[6])[1], xml_children(new_node))
  }
  if(param$ovlMinLen != 40){
        new_node<-read_xml(paste0("<new><value>", param$ovlMinLen, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[2])[7])[1], xml_children(new_node))
  }
  if(param$merSize != 14){
        new_node<-read_xml(paste0("<new><value>", param$merSize, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[2])[8])[1], xml_children(new_node))
  }
  if(param$maxDivergence != 30){
        new_node<-read_xml(paste0("<new><value>", param$maxDivergence, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[7])[1])[3])[1], xml_children(new_node))
  }
  if(param$minAnchorSize != 12){
        new_node<-read_xml(paste0("<new><value>", param$minAnchorSize, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[7])[1])[4])[1], xml_children(new_node))
  }
  write_xml(setting, "settings.xml")
  ezSystem(SMRT_CMD)
  ezSystem(paste("cp", "data/polished_assembly.fasta.gz", start_path))
  setwd(start_path)
  cmd = "gunzip -d polished_assembly.fasta.gz"
  ezSystem(cmd)
  ezSystem(paste("mv", "polished_assembly.fasta", basename(output$getColumn("Draft")))) 
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
                  appDefaults <<- rbind(minSubReadLength = ezFrame(Type="integer",  DefaultValue="500",  Description="Minimum Subread Length. Subreads shorter than this value (in bps) are filtered out and excluded from pre-assembly."),
					readScore  = ezFrame(Type="numeric", DefaultValue="0.8",  Description="Minimum Polymerase Read Quality. Polymerase reads with lower quality than this value are filtered out and excluded from pre-assembly."),
     					minLength = ezFrame(Type="integer",  DefaultValue="100",  Description="Minimum Polymerase Read Length. Polymerase reads shorter than this value (in bps) are filtered out and excluded from pre-assembly."),
					genomeSize = ezFrame(Type="character",  DefaultValue="5000000",  Description="The approximate genome size, in base pairs."),
                                        xCoverage = ezFrame(Type="integer",  DefaultValue="25",  Description="Fold coverage to target for when picking the minimum fragment length for assembly; typically 15 to 25."),
					targetChunks = ezFrame(Type="integer",  DefaultValue="6", Description="The number of pieces to split the data files into while running pre-assembly."),
					splitBestn = ezFrame(Type="integer",  DefaultValue="10", Description="The number of alignments to consider for each read for a particular chunk."),
					totalBestn = ezFrame(Type="integer",  DefaultValue="24", Description="The number of potential alignments to consider across all chunks for a particular read."),
					minCorCov = ezFrame(Type="integer",  DefaultValue="6", Description="The minimum coverage to maintain correction for a read.  If the coverage falls below that threshold, the read will be broken at that juntion."),
					ovlErrorRate = ezFrame(Type="numeric", DefaultValue="0.06",  Description="Trimming and assembly overlaps above this error limit are not computed."),
					ovlMinLen = ezFrame(Type="integer",  DefaultValue="40",  Description="Overlaps shorter than this length (in bps) are not computed."),
					merSize = ezFrame(Type="integer",  DefaultValue="14",  Description="The length of the seeds (in base pairs) used by the seed-and-extend algorithm."),
					maxDivergence = ezFrame(Type="integer",  DefaultValue="30",  Description="The maximum allowed divergence (in %) of a read from the reference sequence."),
					minAnchorSize = ezFrame(Type="integer",  DefaultValue="12",  Description="The minimum size of the read (in bps) that must match against the reference.")
)
                }
              )
)
              
              
