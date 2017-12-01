
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq OTU object
##' @description Create Phyloseq OTU object mothur OTU files.
##' @param  otuFileName, mothur shared clustered OTU  files.
##' @return Returns a Phyloseq OTU object.

phyloSeqOTU <- function(otuFileName){
otuFile <- read.table(otuFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
rownames(otuFile) <- otuFile$Group
colToDrop <- c("label","Group","numOtus")
otuFile1 <- as.matrix(otuFile[,!names(otuFile)%in%colToDrop])
otuObject <- otu_table(otuFile1, taxa_are_rows = FALSE)
return(otuObject)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq Taxa object
##' @description Create Phyloseq taxa object from mothur taxonomy files.
##' @param  taxaFileName, mothur taxonomy file.
##' @return Returns a Phyloseq Taxa object.

phyloSeqTaxa <- function(taxaFileName){
taxaFile <- read.table(taxaFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
tempList <- lapply(taxaFile$Taxonomy,function(y) unlist(strsplit(y,";")))
taxaMatrix <- as.matrix(ldply(tempList))
rownames(taxaMatrix) <- taxaFile$OTU
colnames(taxaMatrix) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")[1:(ncol(taxaMatrix))]
taxaObject <- tax_table(taxaMatrix)
return(taxaObject)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq Sample object
##' @description Create Phyloseq sample object from dataset.
##' @param  taxaFileName, mothur taxonomy file.
##' @return Returns a Phyloseq Taxa object.
phyloSeqSample <- function(sampleFileName){
sampleFile <- ezRead.table(sampleFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
colToKeep <- grep("Factor",colnames(sampleFile))
colnames(sampleFile) <- sub('\\s\\[Factor\\]',"",colnames(sampleFile))
sampleObject <- sample_data(sampleFile)
return(sampleObject)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq preprocess
##' @description Preprocesses a phyloseq object.
##' @param  phyloseqObj, a phyloseq object.
##' @return Returns a  filtered Phyloseq  object.
phyloSeqPreprocess <- function(phyloseqObj){
  ## Standardize abundances to the median sequencing depth
  total = median(sample_sums(phyloseqObj))
  standf = function(x, t=total) round(t * (x / sum(x)))
  gps = transform_sample_counts(phyloseqObj, standf)
  ## transformed to relative abundance and filtered by abundance
  GPr  = transform_sample_counts(gps, function(x) x / sum(x) )
  GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
  filteredPhyloseqObj = filter_taxa(GPfr, function(x) sum(x > 3*1e-5) > (0.2*length(x)), TRUE)
  filteredPhyloseqObj = subset_taxa(filteredPhyloseqObj, Domain=="Bacteria")
  return(filteredPhyloseqObj)
}

