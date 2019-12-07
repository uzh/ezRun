###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets a bed file
##' @description Gets a bed file from a gtf annotation file.
##' @param param a list of parameters to extract \code{ezRef@@refFeatureFile} from.
##' @template roxygen-template
##' @return Returns the path to the created bed file.
##' @examples
##' \dontrun{
##' param = ezParam()
##' param$ezRef@@refFeatureFile = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' rfbed = getReferenceFeaturesBed(param)
##' }
getReferenceFeaturesBed = function(param){
  bedFile = sub(".gtf$", ".bed", param$ezRef["refFeatureFile"])
  ## bedFile exists
  if (file.exists(bedFile)){
    return(bedFile)
  }
  lockFile = sub(".bed", ".bed.lock", bedFile)
  ## I build the bed file
  if (!file.exists(lockFile)){
    ezWrite(Sys.info(), con=lockFile)
    require(GenomicFeatures)
    message("generating bed file from gtf")
    txdb = makeTxDbFromGFF(param$ezRef["refFeatureFile"], dataSource="FGCZ", organism=NA, taxonomyId=NA, chrominfo = NULL)
    require(rtracklayer)
    export(txdb, bedFile)
    ezSystem(paste("chmod", "g+w", bedFile))
    file.remove(lockFile)
    return(bedFile)
  }
  ## I wait until it is build
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait
    Sys.sleep( 60)
    i = i + 1
  }
  if (file.exists(bedFile)){
    return(bedFile)
  } else {
    stop("bed file unavailable: ", bedFile)
  }
}

# 
# ##' @title Gets the file containing the chromosome sizes
# ##' @description Gets the file containing the chromosome sizes either directly or by the fasta files in the chromosome directory.
# ##' @param param a list of parameters:
# ##' \itemize{
# ##'   \item{ezRef@@refChromSizesFile}{ a character representing the path to the file to get.}
# ##'   \item{ezRef@@refChromDir}{ a character representing the path to the directory to get the information to calculate chromosome sizes. Used if the file does not yet exist.}
# ##' }
# ##' @template roxygen-template
# ##' @return Returns the path to the file containing the chromosome sizes.
# ##' @examples
# ##' param = ezParam()
# ##' param$ezRef@@refChromSizesFile = "example.txt"
# ##' param$ezRef@@refChromDir = "./script/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes"
# ##' getRefChromSizesFile(param)
# getRefChromSizesFile = function(param){
#   if (!file.exists(param$ezRef@refChromSizesFile)){
#     faFiles = list.files(path = param$ezRef@refChromDir, pattern = ".fa$", full.names = TRUE)
#     chromSizes = character()
#     for (ff in faFiles){
#       message(ff)
#       dss = readDNAStringSet(ff)
#       chromSizes[names(dss)] = width(dss)
#     }
#     ezWrite.table(chromSizes, file=param$ezRef@refChromSizesFile, col.names = FALSE)    
#   }
#   return(param$ezRef@refChromSizesFile)
# }

##' @title Cleans genome files
##' @description Removes from the seqence files all descriptions in the header line.
##' @param genomeFile a character specifying the path to a fasta file.
##' @param genesFile a character specifying the path to a gtf file.
##' @template roxygen-template
##' @return Returns a list containing a fasta and a gtf object.
##' @examples
##' gtf = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' fasta = system.file("extdata/genome.fa", package="ezRun", mustWork=TRUE)
##' cg = cleanGenomeFiles(fasta, gtf)
cleanGenomeFiles = function(genomeFile, genesFile){
  require(Biostrings)
  require(rtracklayer)
  genome = readDNAStringSet(genomeFile)
  names(genome) = sub(" .*", "", names(genome))
  gtf <- import(genesFile)
  use <- S4Vectors::`%in%`(seqnames(gtf), names(genome))
  stopifnot(any(use))
  gtf <- gtf[use]
  return(list(genomeSeq=genome, gtf=gtf))
}
