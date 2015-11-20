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
##' If a sample is present in both datasets, the read files are concatenated and a new file is written.
##' If a sample is present in only one dataset it is simply copied
##' The Read Count column must be present and is updated if two files are combined.
##' A new dataset is written.
## newDsDir = "/srv/GT/analysis/p1314/HiSeq-20151116-Combined"
## newDsDir = "/scratch/test_merging_datasets"
.ezCombineReadDatasets = function(ds1, ds2, newDsDir){
  ds1 = ezRead.table("/srv/GT/analysis/hubert/ds1.tsv")
  ds2 = ezRead.table("/srv/GT/analysis/hubert/ds2.tsv")
  colnames(ds2) == colnames(ds1)
  intersect(rownames(ds1), rownames(ds2))
  #vennFromSets(list(ds1=rownames(ds1), ds2=rownames(ds2)))
  
  cols = intersect(colnames(ds1), colnames(ds2))
  ds = rbind(ds1[ , cols],
             ds2[setdiff(rownames(ds2), rownames(ds1)), cols])
  ds$"Read Count" = NA
  
  setwdNew(newDsDir)
  ## below implements only the special case where the ds1 samples are a subset of the ds2 samples
  for (nm in rownames(ds)){
    file2 = file.path("/srv/gstore/projects", ds2[nm, "Read1 [File]"])
    if (nm %in% rownames(ds1)){
      file1 = file.path("/srv/gstore/projects", ds1[nm, "Read1 [File]"])
      fileMerged = paste0("20151116.combined-", nm, "_R1.fastq.gz")
      cmd = paste("zcat", file1, file2, "|", "pigz -p4 --best >", fileMerged)
      cat(cmd, "\n")
      ezSystem(cmd)
      ds[nm, "Read Count"] = ds1[nm, "Read Count"] + ds2[nm, "Read Count"]
      ds[nm, "Read1 [File]"] = file.path(newDsDir, fileMerged)
    } else {
      ezSystem(paste("cp", file2, "."))
      ds[nm, "Read1 [File]"] = file.path(newDsDir, basename(file2))
    }
  }
  
  ezWrite.table(ds, file=file.path(".", "dataset.tsv"), head="Name")
  return(ds)
}


ezCombineReadDatasets = function(df1, df2, newDsDir=NULL){
  
  # these are used a lot
  rowDf1 = rownames(df1)
  rowDf2 = rownames(df2)
  colDf1 = colnames(df1)
  colDf2 = colnames(df2)
  cols = intersect(colDf1, colDf2)
  rowdiff1 = setdiff(rowDf1, rowDf2)
  rowdiff2 = setdiff(rowDf2, rowDf1)
  
  # create new df starting with the rows in df1 and adding those only found in df2 and set the right rownames
  dfNew = rbind(df1[, cols], df2[rowdiff2, cols])
  rownames(dfNew) = c(rowDf1, rowdiff2)
  rowDfNew = rownames(dfNew)
  
  # if colnames are not equal, we're not finished with this part:
  # this code should be able to merge columns correctly while filling in NA's no matter how column names occur in df1 and df2
  if (!setequal(colDf1, colDf2)){
    
    # add df1 columns
    cdiff1 = setdiff(colDf1, colDf2)
    for (i in 1:length(cdiff1)){
      dfNew = cbind(dfNew, c(df1[, cdiff1[i]], rep(NA, length(rowdiff2))))
    }
    
    # add df2 columns while making sure to add everything in the right place
    cdiff2 = setdiff(colDf2, colDf1)
    for (i in 1:length(cdiff2)){
      for (j in 1:nrow(dfNew)){
        colToAdd[j] = ifelse(rowDfNew[j] %in% rowDf2,
                             df2[rowDf2==rowDfNew[j], cdiff2[i]],
                             NA)
      }
      dfNew = cbind(dfNew, colToAdd)
    }
    
    # set the colnames for dfNew
    colnames(dfNew) = c(cols, cdiff1, cdiff2)
  }
  
  # adjust the read count and set the directory
  dfNew$"Read Count" = NA
  cwd = getwd()
  setwdNew(newDsDir)
  
  # loop through rows of dfNew to apply the merging
  for (nm in rowDfNew){
    if (nm %in% rowDf1 && !(nm %in% rowDf2)){
      
      # nm is in df1, but not in df2
      dfNew[nm, "Read Count"] = df1[nm, "Read Count"]
      file = file.path("/srv/gstore/projects", df1[nm, "Read1 [File]"])
      ezSystem(paste("cp", file, "."))
      ds[nm, "Read1 [File]"] = file.path(newDsDir, basename(file))
    } else if (nm %in% rowDf2 && !(nm %in% rowDf1)){
      
      # nm is in df2, but not in df1
      dfNew[nm, "Read Count"] = df2[nm, "Read Count"]
      file = file.path("/srv/gstore/projects", df2[nm, "Read1 [File]"])
      ezSystem(paste("cp", file, "."))
      ds[nm, "Read1 [File]"] = file.path(newDsDir, basename(file))
    } else {
      
      # nm is in df1 and df2, thus they need to be merged. there should be no other case.
      dfNew[nm, "Read Count"] = df1[nm, "Read Count"] + df2[nm, "Read Count"]
      file1 = file.path("/srv/gstore/projects", df1[nm, "Read1 [File]"])
      file2 = file.path("/srv/gstore/projects", df2[nm, "Read1 [File]"])
      fileMerged = paste0("combined-", nm, "_R1.fastq.gz")
      cmd = paste("gunzip -c", file1, file2, "|", "pigz -p4 --best >", fileMerged)
      cat(cmd, "\n")
      ezSystem(cmd)
      ds[nm, "Read1 [File]"] = file.path(newDsDir, fileMerged)
    }
  }
  setwd(cwd)
  return(dfNew)
}

