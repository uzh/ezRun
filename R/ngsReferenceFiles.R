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
getReferenceFeaturesBed <- function(param) {
  bedFile <- str_replace(param$ezRef["refFeatureFile"], "\\.gtf$", ".bed")
  if (file.exists(bedFile)) {
    return(bedFile)
  }
  lockFile <- str_replace(bedFile, "\\.bed$", ".bed.lock")
  ## I build the bed file
  if (!file.exists(lockFile)) {
    require(rtracklayer)
    require(withr)
    write_lines(Sys.info(), file = lockFile)
    defer(file.remove(lockFile))
    message("generating bed file from gtf")
    ## the gtf file must contain at least CDS and exons, otherwise there is an error using itemrgb ... blah
    stopifnot(c("CDS", "exon") %in% import(param$ezRef["refFeatureFile"])$type)
    txdb <- GenomicFeatures::makeTxDbFromGFF(param$ezRef["refFeatureFile"], dataSource="FGCZ", organism=NA, taxonomyId=NA, chrominfo = NULL)
    export.bed(txdb, bedFile)
    return(bedFile)
  }
  ## I wait until it is build
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(bedFile)) {
    return(bedFile)
  } else {
    stop("bed file unavailable: ", bedFile)
  }
}

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
cleanGenomeFiles <- function(genomeFile, genesFile) {
  require(rtracklayer)
  genome <- readDNAStringSet(genomeFile)
  names(genome) <- sub(" .*", "", names(genome))
  gtf <- import(genesFile)
  use <- S4Vectors::`%in%`(seqnames(gtf), names(genome))
  stopifnot(any(use))
  gtf <- gtf[use]
  return(list(genomeSeq = genome, gtf = gtf))
}
