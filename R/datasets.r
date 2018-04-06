###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## code for handling dataset data.frames
##' @title Makes a minimal single end dataset
##' @description Makes a minimal single end dataset by combining the arguments to a data.frame.
##' You have to switch directory to the root of the data store first.
##' @param fqDir a character specifying the path to the directory holding the fastq files. The paths stored in the dataset will be relative to the dataRoot
##' @param species a character specifying the species name.
##' @param adapter1 a character representing the adapter sequence.
##' @param adapter2 a character representing the second adapter sequence in the case of a paired end dataset.
##' @param strandMode a character specifying the strand mode for the dataset.
##' @param dataRoot a string specifying the root directory of the data
##' @template roxygen-template
##' @return Returns a data.frame containing the provided information.
##' @examples
##' fqDir = system.file("extdata/yeast_10k", package="ezRun", mustWork=TRUE)
##' species = "Example"
##' ds = makeMinimalSingleEndReadDataset(fqDir, species)
##' ds2 = makeMinimalPairedEndReadDataset(fqDir, species)
makeMinimalSingleEndReadDataset = function(fqDir, species="", adapter1="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                                       strandMode="both", readCount=NULL, dataRoot=DEFAULT_DATA_ROOT){
  if (!grepl("/$", dataRoot)){
    dataRoot = paste0(dataRoot, "/")
  }
  gzFiles = list.files(fqDir, ".gz$", full.names=TRUE)
  samples = sub(".*-", "", sub(".fastq.gz", "", basename(gzFiles)))
  stopifnot(!duplicated(samples))
  ds = data.frame("Read1 [File]"=sub(dataRoot, "", gzFiles),
                  row.names=samples, stringsAsFactors=FALSE, check.names=FALSE)
  ds$Species = species
  ds$Adapter1 = adapter1
  ds$strandMode = strandMode
  if (is.null(readCount)){
    ds$"Read Count"=countReadsInFastq(gzFiles)
  } else {
    ds$"Read Count" = readCount
  }
  #ezWrite.table(ds, file=file.path(fqDir, "dataset.tsv"), head="Name")
  return(ds)
}

##' @describeIn makeMinimalSingleEndReadDataset Does the same for paired end reads.
makeMinimalPairedEndReadDataset = function(fqDir, species="", adapter1="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                                       strandMode="both", readCount=NULL, readTypeSuffix=c("-R1.fastq.gz", "-R2.fastq.gz"),
                                       dataRoot=DEFAULT_DATA_ROOT){
  if (!grepl("/$", dataRoot)){
    dataRoot = paste0(dataRoot, "/")
  }
  gzFiles = list.files(fqDir, readTypeSuffix[1], full.names=TRUE)
  samples = sub(readTypeSuffix[1], "", basename(gzFiles))
  ds = data.frame("Read1 [File]"=sub(dataRoot, "", gzFiles), row.names=samples, stringsAsFactors=FALSE, check.names=FALSE)
  r2Files = sub(readTypeSuffix[1], readTypeSuffix[2], gzFiles)
  stopifnot(file.exists(r2Files))
  ds$"Read2 [File]" = sub(dataRoot, "", r2Files)
  ds$Adapter1 = adapter1
  ds$Adapter2 = adapter2
  ds$strandMode = strandMode
  ds$Species = species
  if (is.null(readCount)){
    ds$"Read Count"=countReadsInFastq(gzFiles)
  } else {
    ds$"Read Count" = readCount
  }
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
##' ds = EzDataset$new(file=file, dataRoot=NULL)
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
  design = design[ ,factorLevelCount > 1, drop = FALSE]
  if (ncol(design) == 0){
    design$Condition = rownames(design)
  }
  colnames(design) = sub(" \\[.*", "", colnames(design))
  return(design)
}


ezColorsFromDesign = function(design){
  ## see http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
  colorFrame = design[NULL] ## make a data frame with rows but without columns
  cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") ## without black
  for (nm in names(design)){
    conds = design[[nm]]
    n = length(unique(conds))
    if (n > length(cbbPalette)){
      condColors = rainbow(n)
    } else {
      condColors = cbbPalette[1:n]
    }
    names(condColors) = unique(conds)
    colorFrame[[nm]] = condColors[conds]
  }
  return(colorFrame)
}

ezColorMapFromDesign = function(design){
  ## see http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
  colorMap = list()
  cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") ## without black
  for (nm in names(design)){
    conds = design[[nm]]
    n = length(unique(conds))
    if (n > length(cbbPalette)){
      condColors = rainbow(n)
    } else {
      condColors = cbbPalette[1:n]
    }
    names(condColors) = unique(conds)
    colorMap[[nm]] = condColors
  }
  return(colorMap)
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
##' dataRoot = system.file("./inst", package = "ezRun", mustWork = TRUE)
##' newDsDir = "./scratch"
##' ezCombineReadDatasets(ds1, ds2, dataRoot, newDsDir)
ezCombineReadDatasets = function(ds1, ds2, dataRoot="/srv/gstore/projects", newDsDir=NULL){
  rowDs1 = rownames(ds1)
  rowDs2 = rownames(ds2)
  commonCols = intersect(colnames(ds1), colnames(ds2))
  ds1[ , setdiff(colnames(ds2), commonCols)] = NA
  ds2[ , setdiff(colnames(ds1), commonCols)] = NA
  rowdiff1 = setdiff(rowDs1, rowDs2)
  rowdiff2 = setdiff(rowDs2, rowDs1)
  
  # create new dataset starting with the rows in ds1 and adding those only found in ds2 and set the right rownames
  dsNew = rbind(ds1, ds2[rowdiff2, colnames(ds1)])
  rownames(dsNew) = c(rowDs1, rowdiff2)
  rowDsNew = rownames(dsNew)
  
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
    }
    if (nm %in% rowDs2 && !(nm %in% rowDs1)){
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
    }
    if (nm %in% rowDs2 && nm %in% rowDs1){
      # nm is in ds1 and ds2, thus they need to be merged. there should be no other case.
      dsNew[nm, "Read Count"] = ds1[nm, "Read Count"] + ds2[nm, "Read Count"]
      fileRead1.1 = file.path(dataRoot, ds1[nm, "Read1 [File]"])
      fileRead1.2 = file.path(dataRoot, ds2[nm, "Read1 [File]"])
      fileMerged = paste0("combined-", nm, "_R1.fastq.gz")
      cmd = paste("gunzip -c", fileRead1.1, fileRead1.2, "|", "pigz -p4 --best >", fileMerged)
      ezSystem(cmd)
      dsNew[nm, "Read1 [File]"] = file.path(newDsDir, fileMerged)
      if (!read2IsNull.Ds1 && !read2IsNull.Ds2){
        # Read2 exists in both datasets and needs to be merged as well
        fileRead2.1 = file.path(dataRoot, ds1[nm, "Read2 [File]"])
        fileRead2.2 = file.path(dataRoot, ds2[nm, "Read2 [File]"])
        fileMerged = paste0("combined-", nm, "_R2.fastq.gz")
        cmd = paste("gunzip -c", fileRead2.1, fileRead2.2, "|", "pigz -p4 --best >", fileMerged)
        ezSystem(cmd)
        dsNew[nm, "Read2 [File]"] = file.path(newDsDir, fileMerged)
      }
    }
  }
  return(dsNew)
}
