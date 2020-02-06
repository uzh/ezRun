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
                                       strandMode="both", readCount=NULL, readTypeSuffix=c("_R1.fastq.gz", "_R2.fastq.gz"),
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

### -----------------------------------------------------------------
### Combine fastqs from multiple sequencing runs
###
ezCombineReadDatasets = function(..., dataRoot="/srv/gstore/projects",
                                 newDsDir=NULL){
  require(dplyr)
  stopifnot(!is.null(newDsDir))
  dir.create(newDsDir, recursive = TRUE)
  
  ds <- bind_rows(...)
  dsNew <- NULL
  isPaired <- "Read2 [File]" %in% colnames(ds)
  uniqNames <- names(which(table(ds$Name) == 1L))
  dupNames <- unique(ds$Name[duplicated(ds$Name)])
  
  ## uniqNames: copy to the destination
  if(length(uniqNames) > 0L){
    message("Copying files without merging.")
    dsUniq <- filter(ds, Name %in% uniqNames)
    dsUniqNew <- mutate(dsUniq,
                        `Read1 [File]`=file.path(newDsDir,
                                                 basename(`Read1 [File]`)))
    file.copy(from=file.path(dataRoot, dsUniq$`Read1 [File]`),
              to=dsUniqNew$`Read1 [File]`)
    if(isPaired){
      dsUniqNew <- mutate(dsUniqNew, 
                          `Read2 [File]`=file.path(newDsDir, 
                                                   basename(`Read2 [File]`)))
      file.copy(from=file.path(dataRoot, dsUniq$`Read2 [File]`),
                to=dsUniqNew$`Read2 [File]`)
    }
    dsNew <- dsUniqNew
  }
  
  ## dupNames: merge to the destination
  if(length(dupNames) > 0L){
    message("Merging files.")
    dsDupList <- list()
    dupName <- dupNames[1]
    for(dupName in dupNames){
      dsDup <- filter(ds, Name == dupName)
      dsDupNew <- head(dsDup, 1) %>%
        mutate(`Read1 [File]`=file.path(newDsDir, paste0("combined-", dupName, "_R1.fastq.gz")))
      dsDupNew$`Read Count` <- sum(dsDup$`Read Count`)
      cmd <- paste("gunzip -c", 
                   paste(file.path(dataRoot, dsDup$`Read1 [File]`), 
                         collapse=" "),
                   "| pigz -p8 --best >",
                   dsDupNew$`Read1 [File]`)
      ezSystem(cmd)
      if(isPaired){
        dsDupNew <- mutate(dsDupNew,
                           `Read2 [File]`=file.path(newDsDir, paste0("combined-", dupName, "_R2.fastq.gz")))
        cmd <- paste("gunzip -c", 
                     paste(file.path(dataRoot, dsDup$`Read2 [File]`), 
                           collapse=" "),
                     "| pigz -p8 --best >",
                     dsDupNew$`Read2 [File]`)
        ezSystem(cmd)
      }
      dsDupList[[dupName]] <- dsDupNew
    }
    dsDupNew <- bind_rows(dsDupList)
    if(!is.null(dsNew)){
      dsNew <- bind_rows(dsNew, dsDupNew)
    }else{
      dsNew <- dsDupNew
    }
  }
  return(dsNew)
}
