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
##' @param adapter1 a character representing the adapter sequence.
##' @param adapter2 a character representing the second adapter sequence in the case of a paired end dataset.
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
##' file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
##' ds = EzDataset$new(file=file)
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

##' @title Combine the reads from two datasets in a single dataset
##' @description Takes the union of the samples in both input datasets and generates a new dataset.
##' @param ds1 a data.frame from the meta field of an EzDataset.
##' @param ds2 a data.frame from the meta field of an EzDataset.
##' @param dataRoot a character specifying the file root of the datasets.
##' @param newDsDir a character specifying the directory to save the new dataset in.
##' @template roxygen-template
##' @details 
##' If a sample is present in both datasets, the read files are concatenated and a new file is written.
##' If a sample is present in only one dataset it is simply copied
##' The Read Count column must be present and is updated if two files are combined.
##' A new dataset is written.
##' @examples 
##' ds1 = ezRead.table(system.file("extdata/yeast_10k/dataset.tsv", package = "ezRun", mustWork = TRUE))
##' ds2 = ds1
##' dataRoot = normalizePath("./inst")
##' newDsDir = "./scratch"
##' ezCombineReadDatasets(ds1, ds2, dataRoot, newDsDir)
ezCombineReadDatasets = function(ds1, ds2, dataRoot="/srv/gstore/projects", newDsDir=NULL){
  
  # these are used a lot
  rowDs1 = rownames(ds1)
  rowDs2 = rownames(ds2)
  colDs1 = colnames(ds1)
  colDs2 = colnames(ds2)
  cols = intersect(colDs1, colDs2)
  rowdiff1 = setdiff(rowDs1, rowDs2)
  rowdiff2 = setdiff(rowDs2, rowDs1)
  
  # create new dataset starting with the rows in ds1 and adding those only found in ds2 and set the right rownames
  dsNew = rbind(ds1[, cols], ds2[rowdiff2, cols])
  rownames(dsNew) = c(rowDs1, rowdiff2)
  rowDsNew = rownames(dsNew)
  
  # if colnames are not equal, we're not finished with this part:
  # this code should be able to merge columns correctly while filling in NA's no matter how column names occur in ds1 and ds2
  if (!setequal(colDs1, colDs2)){
    # add ds1 columns
    cdiff1 = setdiff(colDs1, colDs2)
    for (i in 1:length(cdiff1)){
      dsNew = cbind(dsNew, c(ds1[, cdiff1[i]], rep(NA, length(rowdiff2))))
    }
    # add ds2 columns while making sure to add everything in the right place
    cdiff2 = setdiff(colDs2, colDs1)
    for (i in 1:length(cdiff2)){
      for (j in 1:nrow(dsNew)){
        colToAdd[j] = ifelse(rowDsNew[j] %in% rowDs2,
                             ds2[rowDs2==rowDsNew[j], cdiff2[i]],
                             NA)
      }
      dsNew = cbind(dsNew, colToAdd)
    }
    # set the colnames for dsNew
    colnames(dsNew) = c(cols, cdiff1, cdiff2)
  }
  
  # adjust the read count and set the directory
  dsNew$"Read Count" = NA
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(newDsDir)
  
  # loop through rows of dsNew to apply the merging
  for (nm in rowDsNew){
    read2IsNull.Ds1 = is.null(ds1[nm, "Read2 [File]"])
    read2IsNull.Ds2 = is.null(ds2[nm, "Read2 [File]"])
    if (nm %in% rowDs1 && !(nm %in% rowDs2)){
      # nm is in ds1, but not in ds2
      dsNew[nm, "Read Count"] = ds1[nm, "Read Count"]
      fileRead1 = file.path(dataRoot, ds1[nm, "Read1 [File]"])
      ezSystem(paste("cp", fileRead1, "."))
      dsNew[nm, "Read1 [File]"] = file.path(newDsDir, basename(fileRead1))
      if (!read2IsNull.Ds1){
        fileRead2 = file.path(dataRoot, ds1[nm, "Read2 [File]"])
        ezSystem(paste("cp", fileRead2, "."))
        dsNew[nm, "Read2 [File]"] = file.path(newDsDir, basename(fileRead2))
      }
    } else if (nm %in% rowDs2 && !(nm %in% rowDs1)){
      # nm is in ds2, but not in ds1
      dsNew[nm, "Read Count"] = ds2[nm, "Read Count"]
      fileRead1 = file.path(dataRoot, ds2[nm, "Read1 [File]"])
      ezSystem(paste("cp", fileRead1, "."))
      dsNew[nm, "Read1 [File]"] = file.path(newDsDir, basename(fileRead1))
      if (!read2IsNull.Ds2){
        fileRead2 = file.path(dataRoot, ds2[nm, "Read2 [File]"])
        ezSystem(paste("cp", fileRead2, "."))
        dsNew[nm, "Read2 [File]"] = file.path(newDsDir, basename(fileRead2))
      }
    } else {
      # nm is in ds1 and ds2, thus they need to be merged. there should be no other case.
      dsNew[nm, "Read Count"] = ds1[nm, "Read Count"] + ds2[nm, "Read Count"]
      fileRead1.1 = file.path(dataRoot, ds1[nm, "Read1 [File]"])
      fileRead1.2 = file.path(dataRoot, ds2[nm, "Read1 [File]"])
      fileMerged = paste0("combined-", nm, "_R1.fastq.gz")
      cmd = paste("gunzip -c", fileRead1.1, fileRead1.2, "|", "pigz -p4 --best >", fileMerged)
      cat(cmd, "\n")
      ezSystem(cmd)
      dsNew[nm, "Read1 [File]"] = file.path(newDsDir, fileMerged)
      if (!read2IsNull.Ds1 && !read2IsNull.Ds2){
        # Read2 exists in both datasets and needs to be merged as well
        fileRead2.1 = file.path(dataRoot, ds1[nm, "Read2 [File]"])
        fileRead2.2 = file.path(dataRoot, ds2[nm, "Read2 [File]"])
        fileMerged = paste0("combined-", nm, "_R2.fastq.gz")
        cmd = paste("gunzip -c", fileRead2.1, fileRead2.2, "|", "pigz -p4 --best >", fileMerged)
        cat(cmd, "\n")
        ezSystem(cmd)
        dsNew[nm, "Read2 [File]"] = file.path(newDsDir, fileMerged)
      } else if (!read2IsNull.Ds1 && read2IsNull.Ds2){
        # Read2 only exists in ds1
        fileRead2 = file.path(dataRoot, ds1[nm, "Read2 [File]"])
        ezSystem(paste("cp", fileRead2, "."))
        dsNew[nm, "Read2 [File]"] = file.path(newDsDir, basename(fileRead2))
      } else if (!read2IsNull.Ds2 && read2IsNull.Ds1){
        # Read2 only exists in ds2
        fileRead2 = file.path(dataRoot, ds2[nm, "Read2 [File]"])
        ezSystem(paste("cp", fileRead2, "."))
        dsNew[nm, "Read2 [File]"] = file.path(newDsDir, basename(fileRead2))
      } 
    }
  }
  return(dsNew)
}
