###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodBdRhapsodySA <- function(input = NA, output = NA, param = NA) {
  sampleName <- input$getNames()
  
  #1. Get the reference
  if (ezIsSpecified(param$refBuild)){
    bdRef <- getBdWtaReference(param)
  } else {
    bdRef <- NULL
  }
  bdRhapsodyFolder <- str_sub(sampleName, 1, 45) %>% str_c("-BD-Rhapsody")
  
  #2. Generate the yml file
  bdYmlFile <- makeBdYmlFile(input, param, bdRef)
  
  #3. Build command
  cmd <- paste(
    "rhapsody",
    "pipeline",
    "--tmpdir-prefix ./bd_sa_tmp/",
    "--outdir", sampleName,
    bdYmlFile
  )
  
  #3. Execute the command
  ezSystem(cmd)
  
  #4. Post-process results
  #4.1 Unzip general results
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  setwd(sampleName)
  cellsMexOutput <- basename(output$getColumn("CountMatrix"))
  ezSystem(sprintf("unzip %s.zip -d %s", cellsMexOutput, cellsMexOutput))
  cellsMexUnfilteredOutput <- basename(output$getColumn("UnfilteredCountMatrix"))
  ezSystem(sprintf("unzip %s.zip -d %s", cellsMexUnfilteredOutput, cellsMexUnfilteredOutput))
  #4.2 Unzip tag-specific count matrices
  if (ezIsSpecified(param$sampleTagsVersion)) {
    postProcessTagResults(param, output, sampleName)
  }
  return("Success")
}

getJobConfig <- function(param) {
  require(yaml)
  jobParams <- list(
    "ResourceRequirement"=list(
      ramMin=param$ram, 
      tmpdirMin=param$scratch, 
      outdirMin=param$scratch
    )
  )
  jobYamlFile <- "job.yml"
  write_yaml(jobParams, jobYamlFile)
  return(jobYamlFile)
}

makeBdYmlFile <- function(input, param, bdRef) {
  require(yaml)
  
  # define read file locations and general parameters
  bdParams <- list(
    "cwl:tools"="Rhapsody",
    "Reads"=list(
      list("class"="File",
           "location"=input$getFullPaths("Read1")),
      list("class"="File",
           "location"=input$getFullPaths("Read2"))
    ),
    "Enable_Refined_Cell_Call"=param$enableRefinedPutativeCellCalling,
    "Generate_Bam"=param$generateBamOutput,
    "Exclude_Intronic_Reads"=param$excludeIntronicReads,
    "Putative_Cell_Call"=param$putativeCellCalling,
    "Run_Name"=input$getNames(),
    "Maximum_Threads"=param$cores
  )
  if (ezIsSpecified(param$exactCellCount)) {
    bdParams$Exact_Cell_Count <- as.integer(param$exactCellCount)
  }
  if (ezIsSpecified(param$expectedCellCount)) {
    bdParams$Expected_Cell_Count <- as.integer(param$expectedCellCount)
  }
  if (ezIsSpecified(param$sampleTagsVersion)) {
    bdParams$Sample_Tags_Version <- param$sampleTagsVersion
    if (ezIsSpecified(param$tagNames)) {
      bdParams$Tag_Names <- param$tagNames
    }
  }
  if (ezIsSpecified(param$vdjSpeciesVersion)) {
    bdParams$VDJ_Version <- param$vdjSpeciesVersion
  }
  
  # lastly the references
  if (ezIsSpecified(param$refBuild)) {
    bdParams$Reference_Archive <- list(
      "class"="File",
      "location"=bdRef
    )
  }
  if (ezIsSpecified(param$abSeqReference)) {
    bdParams$AbSeq_Reference <- list(
      list(
        "class"="File",
        "location"=file.path(param$dataRoot, param$abSeqReference)
      )
    )
  }
  if (ezIsSpecified(param$targetedReference)) {
    bdParams$Targeted_Reference <- list(
      list(
      "class"="File", 
      "location"=file.path(param$dataRoot, param$targetedReference)
      )
    )
  }
  
  bdYaml <- "bd_sa.yaml"
  write_yaml(bdParams, bdYaml, handlers=list(
    logical = function(x) {
      result <- ifelse(x, "true", "false")
      class(result) <- "verbatim"
      return(result)
    }
  ))
  return(bdYaml)
}


makeBdRhapRefYmlFile <- function(fasta_path, gtf_path, archive_prefix, param) {
  require(yaml)
  
  # define file locations and general parameters
  bdParams <- list(
    "cwl:tools"="Rhapsody",
    "Genome_fasta"=list(
      list("class"="File",
           "location"=fasta_path)
    ),
    "Gtf"=list(
      list("class"="File",
           "location"=gtf_path)
    ),
    "Archive_prefix"=archive_prefix,
    "Maximum_Threads"=param$cores,
    "Filtering_off"=FALSE,
    "Extra_STAR_params"="",
    "ResourceRequirement"=list(
      ramMin=param$ram, 
      tmpdirMin=param$scratch, 
      outdirMin=param$scratch
    )
  )
  
  bdYaml <- "bd_sa_make_ref.yaml"
  write_yaml(bdParams, bdYaml, handlers=list(
    logical = function(x) {
      result <- ifelse(x, "true", "false")
      class(result) <- "verbatim"
      return(result)
    }
  ))
  return(bdYaml)
}

getFastqDirs <- function(input, column, sampleName) {
  fastqDirs <- strsplit(input$getColumn(column), ",")[[sampleName]]
  fastqDirs <- file.path(input$dataRoot, fastqDirs)
  return(fastqDirs)
}

getBdWtaReference <- function(param) {
  require(rtracklayer)
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  
  if (ezIsSpecified(param$controlSeqs)) {
    refDir <- file.path(getwd(), "BD-WTA_customised_Ref")
  } else {
    if (ezIsSpecified(param$transcriptTypes)) {
      bdWtaBase <- paste(sort(param$transcriptTypes), collapse = "-")
      ## This is a combination of transcript types to use.
    } else {
      bdWtaBase <- ""
    }
    refDir <- sub(
      "\\.gtf$", paste0("_BD-WTA_SC_", bdWtaBase, "_Index"),
      param$ezRef["refFeatureFile"]
    )
  }
  refArchive <- paste0(refDir, "_Rhap_reference.tar.gz")  # the actual output of the workflow
  
  lockFile <- paste0(refDir, ".lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste(
      "reference building still in progress after",
      INDEX_BUILD_TIMEOUT, "min"
    ))
  }
  ## there is no lock file
  if (file.exists(refArchive)) {
    ## we assume the index is built and complete
    return(refArchive)
  }
  
  ## we have to build the reference
  ezWrite(Sys.info(), con = lockFile)
  on.exit(file.remove(lockFile), add = TRUE)
  
  job <- ezJobStart("BD WTA build")
  
  if (ezIsSpecified(param$controlSeqs)) {
    ## make reference genome
    genomeLocalFn <- tempfile(
      pattern = "genome", tmpdir = getwd(),
      fileext = ".fa"
    )
    file.copy(from = param$ezRef@refFastaFile, to = genomeLocalFn)
    writeXStringSet(getControlSeqs(param$controlSeqs),
                    filepath = genomeLocalFn,
                    append = TRUE
    )
    on.exit(file.remove(genomeLocalFn), add = TRUE)
  } else {
    genomeLocalFn <- param$ezRef@refFastaFile
  }
  
  ## make gtf
  gtfFile <- tempfile(
    pattern = "genes", tmpdir = getwd(),
    fileext = ".gtf"
  )
  if (ezIsSpecified(param$transcriptTypes)) {
    export.gff2(gtfByTxTypes(param, param$transcriptTypes),
                con = gtfFile
    )
  } else {
    file.copy(from = param$ezRef@refFeatureFile, to = gtfFile)
  }
  if (ezIsSpecified(param$controlSeqs)) {
    extraGR <- makeExtraControlSeqGR(param)
    gtfExtraFn <- tempfile(
      pattern = "extraSeqs", tmpdir = getwd(),
      fileext = ".gtf"
    )
    on.exit(file.remove(gtfExtraFn), add = TRUE)
    export.gff2(extraGR, con = gtfExtraFn)
    ezSystem(paste("cat", gtfExtraFn, ">>", gtfFile))
  }
  
  # Download the cwl for building the reference
  cmd <- paste(
    "rhapsody", 
    "makeRhapReference",
    "--tmpdir-prefix ./tmp_ref/",
    makeBdRhapRefYmlFile(genomeLocalFn, gtfFile, basename(refDir), param)
  )
  ezSystem(cmd)
  
  # Move the result over and remove files
  ezSystem(sprintf("mv %s %s", basename(refArchive), dirname(refDir)))
  file.remove(gtfFile)
  
  return(refArchive)
}

postProcessTagResults <- function(param, output, sampleName) {
  datasetByTagFn <- "datasetByTag.tsv"
  # Get the tag numbers and tag names
  if (ezIsSpecified(param$tagNames)) {
    tagsSplit <- str_split(param$tagNames, pattern="-", simplify=TRUE)
    tagNums <- tagsSplit[,1] %>% as.integer()
    tagNames <- tagsSplit[,2]
  } else {
    tagNums <- 1:12  # The maximum number of possible tags
    tagNames <- sprintf("Sample%02d", tagNums)
  }
  # Recursively unzip the by-tag count matrix results
  sampleNameMut <- str_replace_all(sampleName, "_", "-")
  sampleTagZips <- c()
  for (i in seq_along(tagNums)){
    myZipFile <- Sys.glob(sprintf("%s_SampleTag%02d*.zip", sampleNameMut, tagNums[i]))
    sampleTagZips[i] <- switch(as.character(length(myZipFile)),"0"='', "1"=myZipFile)
  }
  tagIsFound <- sampleTagZips != ""
  mtxFolders <- sapply(sampleTagZips[tagIsFound], function(tagZip) {
    sampleTagFolder <- tools::file_path_sans_ext(tagZip)
    ezSystem(sprintf("unzip %s -d %s", tagZip, sampleTagFolder))
    mtxZip <- Sys.glob(file.path(sampleTagFolder, "*.zip"))
    mtxFolder <- basename(tools::file_path_sans_ext(mtxZip))
    ezSystem(sprintf("unzip %s -d %s", mtxZip, mtxFolder))
    unlink(sampleTagFolder, recursive=TRUE)
    return(mtxFolder)
  })
  # Create a dataset file split by sample
  bySampleOutput <- tibble(
    `Name`=tagNames[tagIsFound], 
    `Condition [Factor]`=output$getColumn("Condition"),
    `Species`=output$getColumn("Species"),
    `refBuild`=output$getColumn("refBuild"),
    `SCDataOrigin`=output$getColumn("SCDataOrigin"),
    `CountMatrix [Link]`=file.path(output$getColumn("ResultDir"), mtxFolders),
    `UnfilteredCountMatrix [File]`=output$getColumn("UnfilteredCountMatrix"),
    `ResultDir [File]`=output$getColumn("ResultDir")
  )
  ezWrite.table(bySampleOutput, file=datasetByTagFn, row.names=FALSE)
  return(datasetByTagFn)
}

##' @author Noé, Falko
##' @template app-template
##' @templateVar method ezMethodBdRhapsodySA(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppBDRhapsodySA <-
  setRefClass("EzAppBDRhapsodySA",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBdRhapsodySA
                  name <<- "EzAppBDRhapsodySA"
                  appDefaults <<- rbind(
                    controlSeqs = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "control sequences to add"
                    ),
                    enableRefinedPutativeCellCalling = ezFrame(
                      Type = "logical",
                      DefaultValue = FALSE,
                      Description = "enable refined putative cell calling"
                    ),
                    exactCellCount = ezFrame(
                      Type = "integer",
                      DefaultValue = 10000,
                      Description = "exact number of cells"
                    ),
                    expectedCellCount = ezFrame(
                      Type = "integer",
                      DefaultValue = 10000,
                      Description = "Expected number of cells"
                    ),
                    generateBamOutput = ezFrame(
                      Type = "logical",
                      DefaultValue = FALSE,
                      Description = "Should we output the bam file"
                    ),
                    excludeIntronicReads = ezFrame(
                      Type = "logical",
                      DefaultValue = FALSE,
                      Description = "Should we exclude intronic reads"
                    ),
                    putativeCellCalling = ezFrame(
                       Type = "character",
                       DefaultValue = "mRNA",
                       Description = "Specify the data to be used for putative cell calling"
                    ),
                    tagNames = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "Names of the tags"
                    )
                  )
                }
              )
  )
