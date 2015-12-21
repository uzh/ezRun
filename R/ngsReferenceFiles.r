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
##' getReferenceFeaturesBed(param)
##' }
getReferenceFeaturesBed = function(param){
  bedFile = sub(".gtf$", ".bed", param$ezRef["refFeatureFile"])
  if (!file.exists(bedFile)){
    ezSystem(paste(GTF2BED, "--do-not-sort", "<", param$ezRef["refFeatureFile"], ">", bedFile))
    ezSystem(paste("chmod", "g+w", bedFile))
  }
  return(bedFile)
}

##' @title Gets the file containing the chromosome sizes
##' @description Gets the file containing the chromosome sizes either directly or by the fasta files in the chromosome directory.
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refChromSizesFile}{ a character representing the path to the file to get.}
##'   \item{ezRef@@refChromDir}{ a character representing the path to the directory to get the information to calculate chromosome sizes. Used if the file does not yet exist.}
##' }
##' @template roxygen-template
##' @return Returns the path to the file containing the chromosome sizes.
##' @examples
##' param = ezParam()
##' param$ezRef@@refChromSizesFile = "example.txt"
##' param$ezRef@@refChromDir = "./script/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes"
##' getRefChromSizesFile(param)
getRefChromSizesFile = function(param){
  if (!file.exists(param$ezRef@refChromSizesFile)){
    faFiles = list.files(path = param$ezRef@refChromDir, pattern = ".fa$", full.names = TRUE)
    chromSizes = character()
    for (ff in faFiles){
      message(ff)
      dss = readDNAStringSet(ff)
      chromSizes[names(dss)] = width(dss)
    }
    ezWrite.table(chromSizes, file=param$ezRef@refChromSizesFile, col.names = FALSE)    
  }
  return(param$ezRef@refChromSizesFile)
}

##' @title Cleans genome files
##' @description Removes from the seqence files all descriptionso in the header line.
##' RemovesCleans the fasta and gtf files before they get written into the folder structure.
##' @param genomeFile a character specifying the path to a fasta file.
##' @param genesFile a character specifying the path to a gtf file.
##' @param patchPattern a character specifying the pattern of patches to remove from the genome.
##' @template roxygen-template
##' @return Returns a list containing a fasta and a gtf object.
##' @examples
##' gtf = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' fasta = system.file("extdata/genome.fa", package="ezRun", mustWork=TRUE)
##' cleanGenomeFiles(fasta, gtf)
cleanGenomeFiles = function(genomeFile, genesFile, patchPattern="PATCH"){
  
  genome = readDNAStringSet(genomeFile)
  names(genome) = sub(" .*", "", names(genome))
  genome = genome[!grepl(patchPattern, names(genome))]
  gtf = ezReadGff(genesFile)
  use = gtf$seqid %in% names(genome)
  stopifnot(any(use))
  gtf = gtf[use,]
  return(list(genomeSeq=genome, gtf=gtf))
}
