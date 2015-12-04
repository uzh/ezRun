###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title The S4 class representing a collection of information about the reference build
##' @description Objects of this class are a collection of references derived from a list of parameters.
##' @slot refBuild a character specifying the file path to the reference build.
##' @slot refBuildName a character specifying the name of the reference build.
##' @slot refBuildDir a character specifying the directory of the reference build.
##' @slot refIndex a character specifying the location of the index that is used in the alignment.
##' @slot refFeatureFile a character specifying the file path to the annotation feature file (.gtf).
##' @slot refAnnotationFile a character specifying the file path to the annotation file (.txt).
##' @slot refFastaFile a character specifying the file path to the fasta file.
##' @slot refChromDir a character specifying the file path to the directory of the chromosome information.
##' @slot refChromSizesFile a character specifying the file path to the file containing the chromosome sizes.
##' @slot refAnnotationVersion a character specifying the annotation version.
##' @template roxygen-template
##' @seealso \code{\link{cleanGenomeFiles}}
##' @examples
##' refBuild = "Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18"
##' gtf = system.file("extdata/genes.gtf", package="ezRun", mustWork = TRUE)
##' fp = system.file("extdata/genome.fa", package="ezRun", mustWork = TRUE)
##' genomesRoot = "./refExample"
##' myRef = EzRef(param=ezParam(list(refBuild=refBuild)), genomesRoot=genomesRoot)
##' buildRefDir(myRef, fp, gtf)
##' buildIgvGenome(myRef)
EzRef = setClass("EzRef",
                 slots = c(refBuild="character",
                           refBuildName="character",
                           refBuildDir="character",
                           refIndex="character",
                           refFeatureFile="character",
                           refAnnotationFile="character",
                           refFastaFile="character",
                           refChromDir="character",
                           refChromSizesFile="character",
                           refAnnotationVersion="character"))

##' @describeIn EzRef Initializes the slots of EzRef. It will also try to specify some fields and if necessary get full file paths.
setMethod("initialize", "EzRef", function(.Object, param=list(), genomesRoot=GENOMES_ROOT){
  #   if (!ezIsSpecified(param$refBuild)){
  #     return(.Object)
  #   }
  .Object@refBuild = param$refBuild
  refFields = strsplit(.Object@refBuild, "/", fixed=TRUE)[[1]]
  if (ezIsSpecified(param$refBuildName)){
    .Object@refBuildName = param$refBuildName
  } else {
    .Object@refBuildName = refFields[3]
  }
  if (ezIsSpecified(param$refBuildDir)){
    .Object@refBuildDir = param$refBuildDir
  } else {
    .Object@refBuildDir = file.path(genomesRoot, paste(refFields[1:3], collapse="/"))
  }
  if (length(refFields) == 5 && grepl("^Version", refFields[5])){
    .Object@refAnnotationVersion = refFields[5]
  } else {
    NULL
  }
  if (ezIsSpecified(param$refIndex)){
    .Object@refIndex = param$refIndex
  } else {
    .Object@refIndex = ""
  }
  if (ezIsAbsolutePath(param$refFeatureFile)){
    .Object@refFeatureFile = param$refFeatureFile
  } else {
    .Object@refFeatureFile =  file.path(.Object@refBuildDir, "Annotation", .Object@refAnnotationVersion, "Genes", param$refFeatureFile)
  }
  if (ezIsAbsolutePath(param$refAnnotationFile)){
    .Object@refAnnotationFile = param$refAnnotationFile
  } else {
    .Object@refAnnotationFile =  sub(".gtf$", "_annotation.txt", .Object@refFeatureFile)
  }
  if (ezIsAbsolutePath(param$refFastaFile)){
    .Object@refFastaFile = param$refFastaFile
  } else {
    .Object@refFastaFile =  file.path(.Object@refBuildDir, param$refFastaFile)
  }
  if (ezIsAbsolutePath(param$refChromDir)){
    .Object@refChromDir = param$refChromDir
  } else {
    .Object@refChromDir =  file.path(.Object@refBuildDir, "Sequence/Chromosomes")
  }
  if (ezIsAbsolutePath(param$refChromSizesFile)){
    .Object@refChromSizesFile = param$refChromSizesFile
  } else {
    .Object@refChromSizesFile =  file.path(.Object@refChromDir, "chromsizes.txt")
  }
  return(.Object)
})

##' @describeIn EzRef Allows selection with square brackets [ ].
setMethod("[", "EzRef", function(x, i){
  slot(x, i)
})

setGeneric("getOrganism", function(.Object){
  standardGeneric("getOrganism")
})
##' @describeIn EzRef Gets the organism name from the reference build.
setMethod("getOrganism", "EzRef", function(.Object){
  strsplit(.Object@refBuild,"/")[[1]][1]
})

setGeneric("buildRefDir", function(.Object, genomeFile, genesFile, genomesRoot = "."){
  standardGeneric("buildRefDir")
})
##' @describeIn EzRef Builds the reference directory and copies the annotation and fasta file into the right folders.
setMethod("buildRefDir", "EzRef", function(.Object, genomeFile, genesFile, genomesRoot = "."){
  cd = getwd()
  on.exit(setwd(cd))
  setwdNew(genomesRoot)
  
  gtfPath = dirname(.Object@refFeatureFile)
  fastaPath = dirname(.Object@refFastaFile)
  dir.create(gtfPath, recursive=T)
  dir.create(fastaPath, recursive=T)
  dir.create(.Object@refChromDir)
  if (!is.null(.Object@refAnnotationVersion)){
    ezSystem(paste("cd", file.path(.Object@refBuildDir, "Annotation"), "; ", "ln -s",
                   file.path(.Object@refAnnotationVersion, "*"), "."))
  }
  
  genomeInfoList = cleanGenomeFiles(genomeFile, genesFile)
  writeXStringSet(genomeInfoList$genomeSeq, .Object@refFastaFile)
  ezWriteGff(genomeInfoList$gtf, .Object@refFeatureFile)
  cmd = paste(SAMTOOLS, "faidx", .Object@refFastaFile)
  ezSystem(cmd)
  cmd = paste("java -jar", PICARD_JAR, "CreateSequenceDictionary",
              paste0("R=", .Object@refFastaFile), paste0("O=", sub(".fa$", ".dict", .Object@refFastaFile)))
  ezSystem(cmd)
})

## should be called after buildRefDir created the folder structure with genes.gtf and genome.fa
## I did not include the following commented out functionality from Masa's ruby script:
# ENTRY_FILE = '/srv/GT/reference/igv_genomes.txt'
# if ADD_LIST
#   zip_file = "http://fgcz-gstore.uzh.ch#{zip_file.gsub('/srv/GT','')}"
#   open(ENTRY_FILE, 'a') do |file|
#     file.print name, "\t", zip_file, "\t", id_base, "\n" 
#   end
# end
setGeneric("buildIgvGenome", function(.Object, param=NULL){
  standardGeneric("buildIgvGenome")
})
##' @describeIn EzRef Builds the IGV genome.
setMethod("buildIgvGenome", "EzRef", function(.Object, param=NULL){
  
  ## create transcript.only.gtf
  gtfFile = .Object@refFeatureFile
  genomeFile = .Object@refFastaFile
  stopifnot(file.exists(gtfFile))
  stopifnot(file.exists(genomeFile))
  tryCatch({
    gtf = ezLoadFeatures(param=param, featureFile=gtfFile, types="exon")
    transcriptGtf = groupGff(gtf, grouping=gtf$transcript_id, type="transcript")
    transcriptGtf$attributes = ezBuildAttributeField(transcriptGtf[ , c("transcript_id", "gene_id", "gene_name")])
    trxFile = file.path(.Object@refBuildDir, "transcripts.only.gtf")
    ezWriteGff(transcriptGtf, trxFile)
    cmd = paste0("sort -k1,1 -k4,4n -k5,5n  ", trxFile, " -o ", trxFile)
    ezSystem(cmd)
  }, error=function(e){
    message("Could not load features. Copy the annotation file instead.")
    trxFile = file.path(.Object@refBuildDir, "transcripts.only.gtf")
    cmd = paste("cp", gtfFile, trxFile)
    ezSystem(cmd)
  })
  
  ## sort and index genes.gtf
  sortedGtfFile = file.path(dirname(gtfFile), "genes.sorted.gtf")
  cmd = paste("/usr/local/ngseq/bin/igvtools sort", gtfFile, sortedGtfFile)
  ezSystem(cmd)
  cmd = paste("/usr/local/ngseq/bin/igvtools index", sortedGtfFile)
  ezSystem(cmd)
  
  ## make chrom_alias.tab
  chromFile = file.path(.Object@refBuildDir, "chrom_alias.tab")
  cmd = paste("grep '>'", genomeFile, ">", chromFile)
  ezSystem(cmd)
  
  ## make property.txt
  propertyFile = file.path(.Object@refBuildDir, "property.txt")
  id = .Object@refBuildName
  name = paste(getOrganism(.Object), id, sep="_")
  properties = c("fasta=true", "fastaDirectory=false", "ordered=true")
  properties = c(properties, paste0("id=", id), paste0("name=", name), paste0("sequenceLocation=", genomeFile))
  properties = c(properties, "geneFile=transcripts.only.gtf", "chrAliasFile=chrom_alias.tab")
  writeLines(properties, con=propertyFile)
  
  ## make zip file and clean up
  zipFile = paste0("igv_", id, ".genome")
  filesToZip = c(trxFile, chromFile, propertyFile)
  cd = getwd()
  setwd(.Object@refBuildDir)
  zipFile(basename(filesToZip), zipFile)
  setwd(cd)
  cmd = paste("rm -fr", paste(filesToZip, collapse=" "))
  ezSystem(cmd)
})
