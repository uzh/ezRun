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
##' param = ezParam()
##' param$ezRef@@refFeatureFile = "./inst/extdata/genes.gtf"
##' getReferenceFeaturesBed(param)
getReferenceFeaturesBed = function(param){
  bedFile = sub(".gtf$", ".bed", param$ezRef["refFeatureFile"])
  if (!file.exists(bedFile)){
    ezSystem(paste(GTF2BED, param$ezRef["refFeatureFile"], ">", bedFile))
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
