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
##' @template roxygen-template
##' @examples
##' refBuild = "Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18"
##' gtf = "genes.gtf"
##' fasta = "genome.fa"
##' genomesRoot = "~/refExample"
##' myRef = EzRef(param=ezParam(list(refBuild=refBuild)), genomesRoot=genomesRoot)
##' ezBuildRefdir(myRef, gtf, fasta, genomesRoot)
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

setGeneric("ezBuildRefdir", function(.Object, refBuild, gtf, fasta, genomesRoot = "."){
  standardGeneric("ezBuildRefdir")
})

##' @describeIn EzRef Builds the reference directory and copies the annotation and fast file into the right folders.
setMethod("ezBuildRefdir", "EzRef", function(.Object, refBuild, gtf, fasta, genomesRoot = "."){
  cd = getwd()
  setwdNew(genomesRoot)
  gtfPath = dirname(.Object@refFeatureFile)
  fastaPath = dirname(.Object@refFastaFile)
  dir.create(gtfPath, recursive=T)
  dir.create(fastaPath, recursive=T)
  dir.create(.Object@refChromDir)
  #gtf = list.files(pattern=".gtf$")
  #fasta = list.files(pattern=".fa$")
  file.copy(gtf, gtfPath)
  file.copy(fasta, fastaPath)
  setwd(cd)
})
