###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodResequencing = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
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
  inputXmlFile <-system.file(package = "ezRun", "extdata/resequencing_settings.xml")
  setting<-read_xml(inputXmlFile)
  if (param$minSubReadLength != 50){
	new_node<-read_xml(paste0("<new><value>", param$minSubReadLength, "</value></new>"))
	xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[2])[1], xml_children(new_node))
  }
  if(param$readScore != 75){
	new_node<-read_xml(paste0("<new><value>", param$readScore, "</value></new>"))
	xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[3])[1], xml_children(new_node))
  }
  if(param$minLength != 50){
        new_node<-read_xml(paste0("<new><value>", param$minLength, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[3])[1])[4])[1], xml_children(new_node))
  }
  if(param$maxDivergence != 30){
        new_node<-read_xml(paste0("<new><value>", param$maxDivergence, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[1])[3])[1], xml_children(new_node))
  }
  if(param$minAnchorSize != 12){
        new_node<-read_xml(paste0("<new><value>", param$minAnchorSize, "</value></new>"))
        xml_replace(xml_children(xml_children(xml_children(xml_children(setting)[5])[1])[4])[1], xml_children(new_node))
  }
  new_node<-read_xml(paste0("<new><value>common/userdata.d/references/", param$refFile, "</value></new>"))
  xml_replace(xml_children(xml_children(xml_children(setting)[1])[5]), xml_children(new_node))
  write_xml(setting, "settings.xml")
  ezSystem(SMRT_CMD)
  ezSystem(paste("cp", "data/aligned_reads.bam", start_path))
  ezSystem(paste("cp", "data/aligned_reads.bam.bai", start_path))
  ezSystem(paste("cp", "data/variants.vcf", start_path))
  ezSystem(paste("cp", "data/consensus.fasta.gz", start_path))
  setwd(start_path)
  cmd = "gunzip -d consensus.fasta.gz"
  ezSystem(cmd)
  cmd = "bgzip -c variants.vcf > variants.vcf.gz"
  ezSystem(cmd)
  cmd = "tabix -p vcf variants.vcf.gz"
  ezSystem(cmd)
  ezSystem(paste("mv", "aligned_reads.bam", basename(output$getColumn("BAM"))))
  ezSystem(paste("mv", "aligned_reads.bam.bai", basename(output$getColumn("BAI"))))
  ezSystem(paste("mv", "variants.vcf.gz", basename(output$getColumn("VCF"))))
  ezSystem(paste("mv", "variants.vcf.gz.tbi", basename(output$getColumn("VCFINDEX"))))
  ezSystem(paste("mv", "consensus.fasta", basename(output$getColumn("Consensus")))) 
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppResequencing <-
  setRefClass("EzAppResequencing",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodResequencing
                  name <<- "EzAppResequencing"
                  appDefaults <<- rbind(minSubReadLength = ezFrame(Type="integer",  DefaultValue="50",  Description="Minimum Subread Length. Subreads shorter than this value (in bps) are filtered out and excluded from pre-assembly."),
					readScore  = ezFrame(Type="numeric", DefaultValue="75",  Description="Minimum Polymerase Read Quality. Polymerase reads with lower quality than this value are filtered out and excluded from pre-assembly."),
     					minLength = ezFrame(Type="integer",  DefaultValue="50",  Description="Minimum Polymerase Read Length. Polymerase reads shorter than this value (in bps) are filtered out and excluded from pre-assembly."),
					maxDivergence = ezFrame(Type="integer",  DefaultValue="30",  Description="The maximum allowed divergence (in %) of a read from the reference sequence."),
					minAnchorSize = ezFrame(Type="integer",  DefaultValue="12",  Description="The minimum size of the read (in bps) that must match against the reference.")
)
                }
              )
)
              
              
