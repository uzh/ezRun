###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


# ezSupportedParam = function(appParam=NULL){
#   commonParam = c("resultDir", "processMode", "mail", "adminMail")
#   union(commonParam, appParam)
# }



##' @title The R5 class representing a dataset
##' @description Use this to create an object of the class EzDataset that contains the necessary fields for the input and output datasets.
##' @field file a character representing the file path of the dataset's contents.
##' @field meta a data.frame containing the information about the samples.
##' @field colNames a character vector of the column names in \code{meta}.
##' @field tags a list with possible tags for \code{colNames}.
##' @field isModified whether the dataset has been modified.
##' @section Functions:
##' \itemize{
##'   \item{\code{ezTagListFromNames(names): }}{Gets the tags from \code{names} that are in the format [tag].}
##' }
##' @template roxygen-template
##' @examples 
##' file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
##' dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' ds = EzDataset$new(file=file, dataRoot=dataRoot)
##' ds$file
##' ds$meta
##' ds$getColumn("Read1")
##' ds$getFullPaths("Read1")
##' ds2 = ds$copy()
##' ds2$setColumn("Read1","replacement")
##' ds$columnHasTag("File")
##' ds$getNames()
##' ds$meta$"Genotype [Factor]"[1] = "a"
##' ds2 = EzDataset$new(meta = ds$meta, dataRoot=dataRoot)
EzDataset <-
  setRefClass("EzDataset",
              fields = c("file", "meta", "colNames", "tags", "isModified", "dataRoot"),
              methods = list(
                initialize = function(fileNew=character(0), metaNew=list(), dataRoot=NULL)
                {
                  if (length(metaNew) > 0){
                    if (is.data.frame(metaNew)){
                      meta <<- metaNew ## we accept a data frame for multiple samples
                    } else {
                      stopifnot(length(metaNew$Name) == 1) ## we accept a list for single samples; a sample must always have a Name
                      meta <<- data.frame(metaNew, stringsAsFactors=FALSE, 
                                          check.names=FALSE, row.names=metaNew$Name)
                      meta$Name <<- NULL
                    }
                  } else {
                    if (length(fileNew) == 1){
                      if (length(metaNew) == 0){
                        .waitUntilFileExists(fileNew, maxWaitSeconds=60, interval=1)
                        stopifnot(file.exists(fileNew))
                        file <<- fileNew
                        meta <<- ezRead.table(fileNew)
                      } else {
                        file <<- fileNew
                      }
                    }
                  }
                  tags <<- ezTagListFromNames(base::names(meta))
                  ## reorder the meta-information such that factors come first!
                  if (class(meta) != "uninitializedField"){
                    meta <<- meta[ , order(.self$columnHasTag("Factor"), decreasing = TRUE), drop=FALSE]
                  }
                  tags <<- ezTagListFromNames(base::names(meta))
                  colNames <<- sub(" \\[.*", "", base::names(meta))
                  for (i in which(.self$columnHasTag("Factor"))){
                    if (is.logical(meta[[i]])){
                      toReplace = is.na(meta[[i]])
                      meta[toReplace, i] <<- ''
                    }
                    meta[ ,i] <<- as.character(meta[ ,i])
                    hasBadCharacter = !hasFilesafeCharacters(meta[ ,i])
                    if (any(hasBadCharacter)){
                      stop("Invalid character in: ", colnames(meta)[i], " - ", paste("'", meta[hasBadCharacter ,i], "'", sep="", collapse=" "))
                    }
                  }
                  dataRoot <<- dataRoot
                  isModified <<- FALSE
                },
                getColumn = function(names)
                {
                  "Gets the column(s) selected with \\code{names}."
                  idx = match(names, colNames)
                  if (any(is.na(idx))){
                    stop("Column not found in dataset: ", paste(names[is.na(idx)], collapse=" "),
                         "\nAvailable columns: ", paste(colNames, collapse=" "))
                  }
                  x = meta[[idx]]
                  names(x) = rownames(meta)
                  return(x)
                },
                setColumn = function(name, values)
                {
                  "Sets the column selected with \\code{name} to \\code{values}. If \\code{values} is \\code{NULL} the column gets removed"
                  idx = match(name, colNames)
                  if (any(is.na(idx))){
                    stop("Column not found in dataset: ", paste(name[is.na(idx)], collapse=" "),
                       "\nAvailable columns: ", paste(colNames, collapse=" "))
                    meta[[name]] <<- NA
                    colNames <<- sub(" \\[.*", "", base::names(meta))
                    idx = match(sub(" \\[.*", "", name), colNames)
                  }
                  meta[ ,idx] <<- values
                  if (is.null(values)){
                    ## if the values are NULL the column gets remove
                    colNames <<- colNames[-idx]
                    tags <<- tags[-idx]
                  }
                  isModified <<- TRUE
                },
                columnHasTag = function(tag)
                {
                  "Checks each column whether its \\code{tags} matches \\code{tag}."
                  return(grepl(tag, tags))
                },
                subset = function(samples)
                {
                  "Subsets the meta field keeping \\code{samples} and generates a new EzDataset"
                  # meta <<- meta[samples, , drop=FALSE]
                  # isModified <<- TRUE
                  # return(.self)
                  return(EzDataset(meta=meta[samples, , drop=FALSE], dataRoot=dataRoot))
                },
                getNames = function()
                {
                  "Gets the row names."
                  return(rownames(meta))
                },
                getLength = function()
                {
                  "Gets the number of samples."
                  return(length(rownames(meta)))
                },
                getFullPaths = function(name)
                {
                  "Gets the files in the nameed column prepended with the \\code{dataRoot}."
                  ### ok = ezSystem(paste("cd", dataRoot, "; pwd")) ### workaround to make sure the drive where the data sits is mounted by the automounter
                  files = .self$getColumn(name)
                  if (is.null(dataRoot) || dataRoot == "" ){
                    fullPaths = files
                  } else {
                    fullPaths = file.path(dataRoot, files)
                    names(fullPaths) = names(files)
                  }
                  isInvalid = file.access(fullPaths) != 0
                  if (any(isInvalid)){
                    stop("Files are not readable using root:\n", paste(dataRoot, collapse="\n"), "\nfiles:\n", paste(files[isInvalid], collapse="\n"))
                  }
                  return(fullPaths)
                },
                readType = function(){
                  if("Read1" %in% colNames){
                    isFastq <- all(grepl("\\.(fastq|fq)(\\.gz){0,1}$", 
                                         .self$getColumn("Read1")))
                    isBam <- all(grepl("bam$", .self$getColumn("Read1"),
                                       ignore.case = TRUE))
                    stopifnot(isFastq || isBam)
                    if(isTRUE(isFastq)){
                      return("fastq")
                    }else if(isTRUE(isBam)){
                      return("bam")
                    }
                  }else{
                    return(NA)
                  }
                }
              )
  )
# require(ezRun)
# options(error=recover)

# @describeIn doesn't work in RC classes. It is described "manually" in EzDataset.
ezTagListFromNames = function(names){
  lapply(names, function(nm){
    if (grepl("\\[", nm)){
      unlist(strsplit(sub(".*\\[(.*)\\]", "\\1", nm), ","))
    } else {
      return(NULL)
    }
  })
}

##' @title The R5 class representing a runnable app
##' @description This reference class is the basis of all other apps that inherit from it. It sets the framework to run different apps.
##' @field runMethod the function that will be executed in the \code{run} method.
##' @field name the name of the app.
##' @field appDefaults the defaults to run the application with.
##' @param input a list, file path or an object of the class EzDataset containing the input.
##' @param output a list, file path or an object of the class EzDataset containing the output information.
##' @param param a list of parameters to customize the application run.
##' @section Applications inheriting from \code{EzApp}:
##' \itemize{
##'   \item{\code{EzAppBamPreview: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppBismark: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppBowtie: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppBowtie2: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppBWA: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppChipStats: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppCountOverlaps: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppCountQC: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppDeseq2: }}
##'   {Use this reference class to run a differential expression analysis with the application deseq2 on two groups.}
##'   \item{\code{EzAppEdger: }}
##'   {Use this reference class to run a differential expression analysis with the application edgeR on two groups.}
##'   \item{\code{EzAppEdgerMulti: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppFastqc: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppFastqScreen: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppFeatureCounts: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppFlash: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppGatkRnaHaplotyper: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppMacs2: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppMEME: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppMpileup: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppNcpro: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppRSEM: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppSTAR: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppTeqc: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppTophat: }}
##'   {Use this reference class to run }
##'   \item{\code{EzAppTrinity: }}
##'   {Use this reference class to run }
##' }
##' @template roxygen-template
##' @seealso \code{\link{EzDataset}}
##' @seealso \code{\link{waitForFreeDiskSpace}}
##' @seealso \code{\link{EzAppBamPreview}}
##' @seealso \code{\link{EzAppBismark}}
##' @seealso \code{\link{EzAppBowtie}}
##' @seealso \code{\link{EzAppBowtie2}}
##' @seealso \code{\link{EzAppBWA}}
##' @seealso \code{\link{EzAppChipStats}}
##' @seealso \code{\link{EzAppCountOverlaps}}
##' @seealso \code{\link{EzAppCountQC}}
##' @seealso \code{\link{EzAppDeseq2}}
##' @seealso \code{\link{EzAppEdger}}
##' @seealso \code{\link{EzAppEdgerMulti}}
##' @seealso \code{\link{EzAppFastqc}}
##' @seealso \code{\link{EzAppFastqScreen}}
##' @seealso \code{\link{EzAppFeatureCounts}}
##' @seealso \code{\link{EzAppFlash}}
##' @seealso \code{\link{EzAppGatkRnaHaplotyper}}
##' @seealso \code{\link{EzAppMacs2}}
##' @seealso \code{\link{EzAppMEME}}
##' @seealso \code{\link{EzAppMpileup}}
##' @seealso \code{\link{EzAppNcpro}}
##' @seealso \code{\link{EzAppRSEM}}
##' @seealso \code{\link{EzAppSTAR}}
##' @seealso \code{\link{EzAppTeqc}}
##' @seealso \code{\link{EzAppTophat}}
##' @seealso \code{\link{EzAppTrinity}}
##' @examples
##' require("ezRun")
##' file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
##' ds = EzDataset$new(file=file, dataRoot=NULL)
##' NULLApp = EzApp$new(runMethod=function(input, output, param){},name="NULLApp")
##' NULLApp$run(input=ds, output=ds, param=list(process_mode="DATASET"))
EzApp <- 
  setRefClass("EzApp",
              fields = list(runMethod="function",
                            name="character",
                            appDefaults="data.frame",
                            stackTrace="character"),
              methods = list(
                run = function (input, output, param) 
                {
                  "Runs the app with the provided \\code{input}, \\code{output} and \\code{param}."
                  if (is.list(input)){
                    input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
                  } else {
                    if (is.character(input)){
                      input = EzDataset$new(file=input, dataRoot=param$dataRoot)
                    }
                  }
                  if (is.list(output)){
                    output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
                  } else {
                    if (is.character(output)){
                      output = EzDataset$new(file=output, dataRoot=param$dataRoot)
                    }
                  }
                  on.exit(.self$appExitAction(param, output, appName=name))
                  withCallingHandlers({
                    if (param$process_mode == "SAMPLE"){
                      if (input$getLength() > 1){
                        stop("process mode is SAMPLE but the input contains more than one sample.")
                      }
                    }
                    options(cores=param$cores)
                    param$appName = name
                    logMessage(name, param, "Starting")
                    param = ezParam(param, appDefaults=appDefaults)
                    cleanForFreeDiskSpace(param)
                    #waitForFreeDiskSpace(param)
                    jobDir = getwd()
                    result = runMethod(input=input$copy(), output=output$copy(), param=param)
                    setwd(jobDir)
                    return(result)
                  }, error=function(e){
                    dump.frames(format(Sys.time(), format="dump_%Y%m%d%H%M%S"), to.file=TRUE);
                    stackTrace <<- limitedLabels(sys.calls(), maxwidth = 200);
                  }
                  )
                },
                appExitAction = function(param, output, appName="unknown")
                {
                  "Executes actions on exit of an application. This includes links to the output and possibly sending an e-mail."
                  text=.self$outputLinks(output, param)
                  resultName = switch(param$process_mode,
                                      SAMPLE=output$getNames()[1],
                                      DATASET=param$name)
                  subject=paste(appName, resultName, 'done.', sep=' ')
                  .self$exitMail(text, subject, param)
                  logMessage(appName, param, "Finished")
                },
                outputLinks = function(output, param)
                {
                  "Returns URLs, that are tagged as Links, specified in the output list together with relevant metadata."
                  use = grepl("Link", output$tags)
                  relUrls = c(param$resultDir, unlist(output$meta[use])) ## always show the link to the resultdir and to all Links if available.
                  return(paste(PROJECT_BASE_URL, relUrls, sep="/"))
                },
                exitMail = function(text, subject, param)
                {
                  "Sends a report e-mail to the specified e-mail address. If not valid, an e-mail will be sent to the administrator if there was an error."
                  if (ezIsSpecified(stackTrace)){
                    if (ezValidMail(param$mail)){
                      recipient = param$mail
                    } else {
                      recipient = param$adminMail
                    }
                    message("error exists: ", recipient)
                    if (ezValidMail(recipient)){
                      ezMail(subject = paste("Error: ", subject),
                             text=c(text, " ", geterrmessage(), " ", stackTrace[1:(length(stackTrace)-2)]), 
                             to=recipient)
                      message("mail sent to: ", recipient)
                    } else {
                      message(c(text, " ", geterrmessage(), " ", stackTrace[1:(length(stackTrace)-2)]))
                    }
                  } else {
                    if (ezValidMail(param$mail)){
                      ezMail(subject=subject, text=text, to=param$mail)      
                    }
                  }
                  return()
                }
              )
  )

##' @title Checks if there is enough free disk space
##' @description Checks if there is enough free disk space. If there is not enough disk space, an e-mail will be sent and the job will be put on hold for up to two hours.
##' @param param a list of parameters:
##' \itemize{
##'   \item{\code{scratch}}{ the required disk space in gigabytes.}
##'   \item{\code{mail}}{ the e-mail address of the recipient.}
##' }
##' @param dirPath a character specifying the path of the directory to check the disk space in.
##' @template roxygen-template
##' @examples
##' param = list()
##' param[['mail']] = ''
##' param[['scratch']] = '100'
##' waitForFreeDiskSpace(param)
waitForFreeDiskSpace = function(param){
  if (is.null(param$scratch)){
    return()
  }
  freeSpace = getGigabyteFree(".")
  if (freeSpace < param$scratch){
    if (ezValidMail(param$mail)){
      recipient = param$mail
    } else{
      recipient = param$adminMail
    }
    ezMail(to=recipient,
           subject=paste("Alert: not enough disk space ", Sys.info()["nodename"], "-", getwd()),
           text="Please free up space! Job is on hold for 2 hours and will be terminated afterwards if the issue persists.")
    cat('Wait for free disk space') 
    i = 0
    while(getGigabyteFree(".") < param$scratch & i < 60){
      Sys.sleep( 120)
      i = i + 1
    }
    if (getGigabyteFree(".") < param$scratch) stop("actual free disk space is less than required")
  }
  return()
}

### Check scratch for enough space and clean it
### When the working director is other than scratch, no cleaning.
cleanForFreeDiskSpace <- function(param){
  if (is.null(param$scratch) || !grepl("^(/scratch|/export/local/scratch)", getwd())){
    message("Scratch is not specificed or the current working directory is not under /scratch. No cleaning.")
    return(TRUE)
  }
  
  freeSpace = getGigabyteFree(".")
  i = 0
  while(getGigabyteFree(".") < param$scratch & i < 200){
    if(getGigabyteTotal(".") > 1024){
      ## For big nodes with more than 1TB scratch, only clean for trxcopy
      message("Clean for trxcopy!")
      cleanOldestDir(dirPath="/scratch", user="trxcopy")
    }else{
      message("Clean for all users!")
      cleanOldestDir(dirPath="/scratch", user=NULL)
    }
    Sys.sleep(5)
    i = i + 1
  }
  if (getGigabyteFree(".") < param$scratch){
    if (ezValidMail(param$mail)){
      recipient = param$mail
    } else{
      recipient = param$adminMail
    }
    ezMail(to=recipient,
           subject=paste("Alert: not enough disk space ", Sys.info()["nodename"], "-", getwd()),
           text="Please free up space manually!")
    stop("actual free disk space is less than required")
  }

  return(TRUE)
}

##' @describeIn waitForFreeDiskSpace Gets the number of free gigabytes.
##' @examples 
##' getGigabyteFree(".")
##' getGigabyteFree("/")
getGigabyteFree = function(dirPath){
  as.numeric(strsplit(ezSystem(paste("df", dirPath), intern=TRUE, echo=FALSE), " +")[[2]][4]) / 1e6
}
getGigabyteTotal = function(dirPath){
  as.numeric(strsplit(ezSystem(paste("df", dirPath), intern=TRUE, echo=FALSE), " +")[[2]][2]) / 1e6
}

### Clean the oldest, not used dir
cleanOldestDir <- function(dirPath, user=NULL){
  allDirs <- list.dirs(path=dirPath, recursive=FALSE)
  
  ## Don't clean symlinks
  allDirs <- allDirs[Sys.readlink(allDirs) == ""]
  
  ## Don't clean smrt* , pacbio stuff
  allDirs <- grep("(smrt|pacbio)", allDirs, invert = TRUE, value=TRUE)

  ## Don't clean rstudio folders
  allDirs <- grep("rstudio$", allDirs, invert = TRUE, value=TRUE)
  
  ## Don't clean **.GT folder; created by grid engine
  allDirs <- grep("GT$", allDirs, invert = TRUE, value=TRUE)
  
  allInfo <- file.info(allDirs)
  if(!is.null(user)){
    allInfo <- allInfo[allInfo$uname %in% user, ]
  }
  
  ## Check being used or not
  isUsed <- suppressWarnings(lapply(paste("lsof", rownames(allInfo)), 
                                    system, intern=TRUE))
  isUsed <- lengths(isUsed) != 0L
  if(!all(isUsed)){
    allInfo <- allInfo[!isUsed, ]
    ## order by ctime
    allInfo <- allInfo[order(allInfo$ctime), ]
    message("Deleting ", rownames(allInfo)[1])
    unlink(rownames(allInfo)[1], recursive=TRUE, force=TRUE)
  }
}
