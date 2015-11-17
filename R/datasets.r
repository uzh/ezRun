###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## code for handling dataset data.frames
##' @title Makes a minimal single end dataset
##' @description Makes a minimal single end dataset by combining the arguments to a data.frame.
##' @param fqDir a character specifying the path to the fastq files.
##' @param species a character specifying the species name.
##' @param adapter1 a character representing the adapter sequence. There is also an adapter2 argument for the paired end dataset.
##' @param strandMode a character specifying the strand mode for the dataset.
##' @template roxygen-template
##' @return Returns a data.frame containing the provided information.
##' @examples
##' fqDir = "inst/extdata/yeast_10k"
##' species = "Example"
##' ds = makeMinimalSingleEndReadDataset(fqDir, species)
##' ds2 = makeMinimalPairedEndReadDataset(fqDir, species)
makeMinimalSingleEndReadDataset = function(fqDir, species="", adapter1="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                                       strandMode="both"){
  gzFiles = list.files(fqDir, ".gz$", full.names=TRUE)
  samples = sub(".*-", "", sub(".fastq.gz", "", gzFiles))
  ds = data.frame("Read1 [File]"=gzFiles, row.names=samples, stringsAsFactors=FALSE, check.names=FALSE)
  ds$Species = species
  ds$Adapter1 = adapter1
  ds$strandMode = strandMode
  ds$"Read Count"=""
  #ezWrite.table(ds, file=file.path(fqDir, "dataset.tsv"), head="Name")
  return(ds)
}

##' @describeIn makeMinimalSingleEndReadDataset Does the same for paired end reads.
makeMinimalPairedEndReadDataset = function(fqDir, species="", adapter1="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                                       strandMode="both"){
  gzFiles = list.files(fqDir, "R1.fastq.gz$", full.names=TRUE)
  samples = sub(".*-", "", sub("R1.fastq.gz", "", gzFiles))
  ds = data.frame("Read1 [File]"=gzFiles, row.names=samples, stringsAsFactors=FALSE, check.names=FALSE)
  r2Files = sub("R1.fastq", "R2.fastq", gzFiles)
  stopifnot(file.exists(r2Files))
  ds$"Read2 [File]" = r2Files
  ds$Adapter1 = adapter1
  ds$Adapter2 = adapter2
  ds$strandMode = strandMode
  ds$Species = species
  ds$"Read Count"=""
  #ezWrite.table(ds, file=file.path(fqDir, "dataset.tsv"), head="Name")
  return(ds)
}

##' @title Gets the design from the dataset
##' @description Gets the design from the dataset or parameters, if specified. 
##' @template dataset-template
##' @param param a list of parameters to specify the factors directly using \code{param$factors}.
##' @template roxygen-template
##' @return Returns the factorial design of the dataset.
##' @examples
##' ds = EzDataset$new(file=system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE))
##' design = ezDesignFromDataset(ds$meta)
##' cond = ezConditionsFromDesign(design)
##' addReplicate(apply(design, 1, paste, collapse="_"))
##' addReplicate(cond)
ezDesignFromDataset = function(dataset, param=NULL){
  if (!ezIsSpecified(param$factors)){
    design = dataset[ , grepl("Factor", colnames(dataset)), drop=FALSE]
    if (ncol(design) == 0){
      design$Condition = rownames(design)
    }
  }
  else {
    factorNames = unlist(strsplit(param$factors, split = ","))
    design = dataset[ , paste(factorNames, "[Factor]"), drop=FALSE]
  }
  factorLevelCount = apply(design, 2, function(x){length(unique(x))})
  if (any(factorLevelCount > 1)){
    for (nm in colnames(design)){
      if (length(unique(design[[nm]])) == 1){
        design[[nm]] = NULL
      }
    }
  }
  return(design)
}

##' @describeIn ezDesignFromDataset Gets the conditions from the design. Use \code{maxFactors} to limit the amount of factors.
ezConditionsFromDesign = function(design, maxFactors=2){
  apply(design[ , 1:min(ncol(design), maxFactors), drop=FALSE], 1, function(x){paste(x, collapse="_")})
}

##' @describeIn ezDesignFromDataset A wrapper to get the conditions directly from the dataset.
ezConditionsFromDataset = function(dataset, param=NULL, maxFactors=2){
  ezConditionsFromDesign(ezDesignFromDataset(dataset, param), maxFactors=maxFactors)
}

##' @describeIn ezDesignFromDataset add a replicate id to make values unique
addReplicate = function(x, sep="_", repLabels=1:length(x)){
  repId = ezReplicateNumber(x)
  paste(x, repLabels[repId], sep=sep)
}
