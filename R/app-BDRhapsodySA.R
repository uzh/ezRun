###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodBdRhapsodySA <- function(input = NA, output = NA, param = NA) {
  sampleName <- input$getNames()
  
  #1. Get the reference
  bdRef <- getBdWtaReference(param)
  bdRhapsodyFolder <- str_sub(sampleName, 1, 45) %>% str_c("-BD-Rhapsody")
  
  #2. Generate the yml file
  bdYmlFile <- makeBdYmlFile(input, param, bdRef)
  
  #3.1 Copy over the necessary workflow
  pipelineCwl <- sprintf("rhapsody_pipeline_%s.cwl", param$version)
  ezSystem(sprintf(
    "cp /srv/GT/software/BD_Rhapsody/cwl/v%s/%s .", 
    param$version, pipelineCwl)
  )
  
  #3. Build command
  cmd <- paste(
    "cwltool",
    "--tmpdir-prefix ./bd_sa_tmp/",
    "--outdir", sampleName,
    "--singularity",
    pipelineCwl,
    bdYmlFile
  )
  
  #3. Execute the command
  ezSystem(cmd)
  
  #4. Process result
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  setwd(sampleName)
  cellsMexOutput <- basename(output$getColumn("CountMatrix"))
  ezSystem(sprintf("unzip %s.zip -d %s", cellsMexOutput, cellsMexOutput))

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
    bdParams$Exact_Cell_Count <- param$exactCellCount
  }
  if (ezIsSpecified(param$expectedCellCount)) {
    bdParams$Expected_Cell_Count <- param$expectedCellCount
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
      "class"="File", 
      "location"=file.path(param$dataRoot, param$targetedReference)
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
    extraGR <- makeExtraControlSeqGR(param$controlSeqs)
    gtfExtraFn <- tempfile(
      pattern = "extraSeqs", tmpdir = getwd(),
      fileext = ".gtf"
    )
    on.exit(file.remove(gtfExtraFn), add = TRUE)
    export.gff2(extraGR, con = gtfExtraFn)
    ezSystem(paste("cat", gtfExtraFn, ">>", gtfFile))
  }
  
  # Download the cwl for building the reference
  referenceCwlFile <- sprintf("make_rhap_reference_%s.cwl", param$version)
  ezSystem(sprintf("cp /srv/GT/software/BD_Rhapsody/cwl/v%s/Extra_Utilities/%s .", param$version,referenceCwlFile))
  cmd <- paste(
    "cwltool",
    "--tmpdir-prefix ./tmp_ref/",
    "--singularity",
    referenceCwlFile,
    getJobConfig(param),
    "--Genome_fasta", genomeLocalFn,
    "--Gtf", gtfFile,
    "--Archive_prefix", basename(refDir),
    "--Maximum_threads", param$cores
  )
  ezSystem(cmd)
  
  # Move the result over and remove files
  ezSystem(sprintf("mv %s %s", basename(refArchive), dirname(refDir)))
  file.remove(gtfFile)
  
  return(refArchive)
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
                      Type = "numeric",
                      DefaultValue = 10000,
                      Description = "exact number of cells"
                    ),
                    expectedCellCount = ezFrame(
                      Type = "numeric",
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
