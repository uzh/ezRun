
# Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18


ezBuildRefDir = function(refBuild){
  ## make the directory structure
}



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
##' EzRef(param=list(refBuild="foo"))
EzRef = setClass("EzRef",
                 slots = c(refBuild="character",
                           refBuildName="character",
                           refBuildDir="character",
                           refIndex="character",
                           refFeatureFile="character",
                           refAnnotationFile="character",
                           refFastaFile="character",
                           refChromDir="character",
                           refChromSizesFile="character"))

##' @describeIn EzRef Initializes the slots of EzRef. It will fill also try to specify some fields and if necessary get full file paths.
setMethod("initialize", "EzRef", function(.Object, param=list()){
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
    .Object@refBuildDir = file.path(GENOMES_ROOT, paste(refFields[1:3], collapse="/"))
  }  
  if (ezIsSpecified(param$refIndex)){
    .Object@refIndex = param$refIndex
  } else {
    .Object@refIndex = ""
  }
  if (ezIsAbsolutePath(param$refFeatureFile)){
    .Object@refFeatureFile = param$refFeatureFile
  } else {
    .Object@refFeatureFile =  file.path(.Object@refBuildDir, "Annotation/Genes", param$refFeatureFile)
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
