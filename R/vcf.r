###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets chromosome sizes from a VCF file
##' @description Gets chromosome sizes from a VCF file by accessing its \code{contig} column.
##' @param vcfFile the VCF file to get the chromosome sizes from.
##' @template roxygen-template
##' @return Returns a names vector of the chromosome sizes.
## TODOEXAMPLE: get vcf file with $contig information.
ezChromSizesFromVcf = function(vcfFile){
  require(VariantAnnotation)
  vh = scanVcfHeader(vcfFile)
  contigs = header(vh)$contig
  chromSizes = as.integer(contigs[ , "length"])
  names(chromSizes) = rownames(contigs)
  return(chromSizes)
}


# load the vcf filter and write back
##' @title Filters a VCF file and writes it back
##' @description Filters a VCF file and writes it back to the specified connection.
##' @param vcfFile a character representing the file path of the input VCF file
##' @param vcfFiltFile a character representing the file path of the output VCF file.
##' @param discardMultiAllelic a logical. Removes multi allelic variants if set to TRUE.
##' @param bamDataset a data.frame to help renaming and reordering the fields in genotype.
##' @param param a list of parameters to extract the value from \code{vcfFilt.minAltCount}.
##' @param vcf an object of the class VCF.
##' @template roxygen-template
##' @return Returns a filtered VCF file.
## TODOEXAMPLE: get working .vcf file
ezFilterVcf = function(vcfFile, vcfFiltFile, discardMultiAllelic=TRUE, 
                       bamDataset=bamDataset, param=NULL){
  require(VariantAnnotation)
  vcf = readVcf(vcfFile, genome="genomeDummy")
  genotype = geno(vcf)
  if (discardMultiAllelic){
    isMultiAllelic = apply(genotype$AD, 1, 
                           function(x){any(sapply(x, length) > 2)})
    table(isMultiAllelic)
    vcf = vcf[!isMultiAllelic, ]
    genotype = geno(vcf)
  }
  altCount = apply(genotype$AD, 2, function(x){
    sapply(x, function(y){if (length(y) == 0) return(NA); return(max(y))})
    })
  hasHighAltCount = apply(altCount > param$vcfFilt.minAltCount, 1, any)
  vcf = vcf[which(hasHighAltCount), ]
  genotype = geno(vcf)
  
  ## reorder and rename the fields in genotype
  for (nm in names(genotype)){
    #colnames(genotype[[nm]]) = sub("^RGSM_", "", colnames(genotype[[nm]]))
    stopifnot(setequal(colnames(genotype[[nm]]), rownames(bamDataset)))
    genotype[[nm]] = genotype[[nm]][ , rownames(bamDataset)] ## establish the original order
  }
  geno(vcf) = genotype
  colData(vcf) = DataFrame(Samples=1:ncol(genotype$AD), 
                           row.names=colnames(genotype$AD))
  ezWriteVcf(vcf, vcfFiltFile)
}

##' @describeIn ezFilterVcf Writes a new VCF file.
ezWriteVcf = function(vcf, vcfFiltFile){
  #vcfTemp = sub("-haplo.vcf", "-haplo-filt.vcf", vcfFiles[dsName])
  stopifnot(grepl(".gz$", vcfFiltFile))
  nm = sub(".vcf", "", sub(".gz", "", basename(vcfFiltFile)))
  vcfTemp = paste0(nm , "-", Sys.getpid(), ".vcf")
  writeVcf(vcf, filename=vcfTemp, index=TRUE)
  stopifnot(file.rename(paste0(vcfTemp, ".bgz"), vcfFiltFile))
  stopifnot(file.rename(paste0(vcfTemp, ".bgz.tbi"), paste0(vcfFiltFile, ".tbi")))
  return(invisible(vcfFiltFile))
}
