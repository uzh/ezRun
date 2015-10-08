###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## code for handling dataset data.frames

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


## TODOP: It should be possible to refactor the next four functions to some extent. search in all files for refactorhelp0 to find relevant places
##' @title Gets the design from the dataset
##' @description Gets the design from the dataset or parameters, if specified. 
##' @param dataset a data.frame to get the design from, usually the field meta of an EzDataset object.
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

## this one is used the most, maybe delete the two functions above and refactor where necessary?
##' @describeIn ezDesignFromDataset A wrapper to get the conditions directly from the dataset.
ezConditionsFromDataset = function(dataset, param=NULL, maxFactors=2){
  ezConditionsFromDesign(ezDesignFromDataset(dataset, param), maxFactors=maxFactors)
}

##' @describeIn ezDesignFromDataset Adds the names into the input character vector.
addReplicate = function(x, sep="_", repLabels=1:length(x)){
  
  idx = order(x)
  xc = table(x)
  idReps = lapply(xc, function(x){repLabels[1:x]})
  idBase = rep(names(xc), xc)
  idNew = paste(idBase, unlist(idReps), sep=sep)
  x[idx] = idNew
  return(x)
}
