###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCellRanger = function(input=NA, output=NA, param=NA){
  sampleName = input$getNames()
  sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
  sampleDirs <- file.path(input$dataRoot, sampleDirs)
  if(all(grepl("\\.tar$", sampleDirs))){
    # This is new .tar folder
    lapply(sampleDirs, untar)
    sampleDirs <- sub("\\.tar$", "", basename(sampleDirs))
  }
  sampleDir <- paste(sampleDirs, collapse=",")
  cellRangerFolder = paste0(sampleName, "-cellRanger")
  
  if(param$TenXLibrary == "GEX"){
    refDir <- getCellRangerGEXReference(param)
    cmd <- paste(CELLRANGER, "count", paste0("--id=", cellRangerFolder),
                 paste0("--transcriptome=", refDir),
                 paste0("--fastqs=", sampleDir),
                 paste0("--sample=", sampleName),
                 paste0("--localmem=", param$ram),
                 paste0("--localcores=", param$cores),
                 paste0("--chemistry=", param$chemistry))
  }else if(param$TenXLibrary == "VDJ"){
    refDir <- getCellRangerVDJReference(param)
    cmd <- paste(CELLRANGER, "vdj", paste0("--id=", cellRangerFolder),
                 paste0("--reference=", refDir),
                 paste0("--fastqs=", sampleDir),
                 paste0("--sample=", sampleName),
                 paste0("--localmem=", param$ram),
                 paste0("--localcores=", param$cores))
  }
  
  if(ezIsSpecified(param$cmdOptions)){
    cmd = paste(cmd, param$cmdOptions)
  }
  ezSystem(cmd)
  
  unlink(basename(sampleDirs), recursive=TRUE)
  file.rename(file.path(cellRangerFolder, "outs"),  sampleName)
  unlink(cellRangerFolder, recursive=TRUE)
  
  if(ezIsSpecified(param$controlSeqs)){
    unlink(refDir, recursive = TRUE)
  }
  
  if(param$TenXLibrary == "GEX"){
    require(DropletUtils)
    require(Matrix)
    require(readr)
    countMatrixFn <- list.files(path=file.path(sampleName, 'filtered_feature_bc_matrix'),
                                pattern="\\.mtx(\\.gz)*$", recursive=TRUE,
                                full.names=TRUE)
    sce <- read10xCounts(dirname(countMatrixFn), col.names=TRUE)
    
    cellPhase <- getCellCycle(sce, param$refBuild)
    write_tsv(cellPhase,
              path=file.path(dirname(countMatrixFn), "CellCyclePhase.txt"))
  }
  return("Success")
}

getCellRangerGEXReference <- function(param){
  require(Biostrings)
  require(rtracklayer)
  cwd <- getwd()
  on.exit(setwd(cwd), add=TRUE)
  
  if(ezIsSpecified(param$controlSeqs)){
    refDir = file.path(getwd(), "10X_customised_Ref")
  }else{
    if(ezIsSpecified(param$transcriptTypes)){
      cellRangerBase <- paste(sort(param$transcriptTypes), collapse="-")
      ## This is a combination of transcript types to use.
    }else{
      cellRangerBase <- ""
    }
    if(param$scMode == "SN"){
      refDir = sub("\\.gtf$", paste0("_10XGEX_SN_", cellRangerBase, "_Index"),
                   param$ezRef["refFeatureFile"])
    }else if(param$scMode == "SC"){
      refDir = sub("\\.gtf$", paste0("_10XGEX_SC_", cellRangerBase, "_Index"),
                   param$ezRef["refFeatureFile"])
    }
  }
  
  lockFile = paste0(refDir, ".lock")
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait
    Sys.sleep(60)
    i = i + 1
  }
  if (file.exists(lockFile)){
    stop(paste("reference building still in progress after", 
               INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  if (file.exists(refDir)){
    ## we assume the index is built and complete
    return(refDir)
  }
  
  ## we have to build the reference
  setwd(dirname(refDir))
  ezWrite(Sys.info(), con=lockFile)
  on.exit(file.remove(lockFile), add=TRUE)
  
  job = ezJobStart("10X CellRanger build")
  
  if(ezIsSpecified(param$controlSeqs)){
    ## make reference genome
    genomeLocalFn <- tempfile(pattern="genome", tmpdir=getwd(),
                              fileext = ".fa")
    file.copy(from=param$ezRef@refFastaFile, to=genomeLocalFn)
    writeXStringSet(getControlSeqs(param$controlSeqs), filepath=genomeLocalFn,
                    append=TRUE)
    on.exit(file.remove(genomeLocalFn, add=TRUE))
  }else{
    genomeLocalFn <- param$ezRef@refFastaFile
  }
  
  ## make gtf
  gtfFile <- tempfile(pattern="genes", tmpdir=getwd(),
                      fileext = ".gtf")
  if(ezIsSpecified(param$transcriptTypes)){
    export.gff2(gtfByTxTypes(param, param$transcriptTypes),
                con=gtfFile)
  }else{
    file.copy(from=param$ezRef@refFeatureFile, to=gtfFile)
  }
  if(ezIsSpecified(param$controlSeqs)){
    extraGR <- makeExtraControlSeqGR(param$controlSeqs)
    gtfExtraFn <- tempfile(pattern="extraSeqs", tmpdir=getwd(),
                           fileext = ".gtf")
    on.exit(file.remove(gtfExtraFn), add=TRUE)
    export.gff2(extraGR, con=gtfExtraFn)
    ezSystem(paste("cat", gtfExtraFn, ">>", gtfFile))
  }
  if(param$scMode == "SN"){
    gtf <- import(gtfFile)
    gtf <- gtf[gtf$type == "transcript"]
    gtf$type <- "exon"
    export.gff2(gtf, gtfFile)
  }
  
  cmd <- paste(CELLRANGER, "mkref",
               paste0("--genome=", basename(refDir)),
               paste0("--fasta=", genomeLocalFn),
               paste0("--genes=", gtfFile),
               paste0("--nthreads=", param$cores))
  ezSystem(cmd)
  
  file.remove(gtfFile)
  
  return(refDir)
}

getCellRangerVDJReference <- function(param){
  require(Biostrings)
  require(rtracklayer)
  cwd <- getwd()
  on.exit(setwd(cwd), add=TRUE)
  
  refDir = sub("\\.gtf$", "_10XVDJ_Index",
               param$ezRef["refFeatureFile"])
  
  lockFile = paste0(refDir, ".lock")
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait
    Sys.sleep(60)
    i = i + 1
  }
  if (file.exists(lockFile)){
    stop(paste("reference building still in progress after", 
               INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  if (file.exists(refDir)){
    ## we assume the index is built and complete
    return(refDir)
  }
  
  ## we have to build the reference
  setwd(dirname(refDir))
  ezWrite(Sys.info(), con=lockFile)
  on.exit(file.remove(lockFile), add=TRUE)
  
  job = ezJobStart("10X CellRanger build")
  
  cmd <- paste(CELLRANGER, "mkvdjref",
               paste0("--genome=", basename(refDir)),
               paste0("--fasta=", param$ezRef@refFastaFile),
               paste0("--genes=", param$ezRef@refFeatureFile))
  ezSystem(cmd)
  
  return(refDir)
}

## not used any more
getCellRangerReference <- function(param){
  if(ezIsSpecified(param$controlSeqs)){
    if(param$TenXLibrary == "VDJ"){
      stop("VDJ library with extra control sequences is not implemented yet!")
    }
    if(param$scMode == "SN"){
      stop("Single-nuclei with extra control sequences is not implemented yet!")
    }
    require(Biostrings)
    require(rtracklayer)
    ## make reference genome
    genomeLocalFn <- tempfile(pattern="genome", tmpdir=getwd(),
                              fileext = ".fa")
    file.copy(from=param$ezRef@refFastaFile, to=genomeLocalFn)
    writeXStringSet(getControlSeqs(param$controlSeqs), filepath=genomeLocalFn,
                    append=TRUE)
    on.exit(file.remove(genomeLocalFn, add=TRUE))
    
    ## make gtf
    gtfFile <- tempfile(pattern="genes", tmpdir=getwd(),
                        fileext = ".gtf")
    on.exit(file.remove(gtfFile, add=TRUE))
    
    file.copy(from=param$ezRef@refFeatureFile, to=gtfFile)
    extraGR <- makeExtraControlSeqGR(param$controlSeqs)
    gtfExtraFn <- tempfile(pattern="extraSeqs", tmpdir=getwd(),
                           fileext = ".gtf")
    on.exit(file.remove(gtfExtraFn), add=TRUE)
    export.gff2(extraGR, con=gtfExtraFn)
    ezSystem(paste("cat", gtfExtraFn, ">>", gtfFile))
    
    ## build the index
    refDir <- file.path(getwd(), "10X_customised_Ref")
    cmd <- paste(CELLRANGER, "mkref",
                 paste0("--genome=", basename(refDir)),
                 paste0("--fasta=", genomeLocalFn),
                 paste0("--genes=", gtfFile),
                 paste0("--nthreads=", param$cores))
    ezSystem(cmd)
  }else{
    ## TODO: automate the reference building
    refDir <- dirname(param$ezRef["refFeatureFile"])
    if(param$TenXLibrary == "VDJ"){
      refDirs <- list.files(path=refDir, pattern="^10X_Ref.*_VDJ_", 
                            full.names = TRUE)
    }else if(param$TenXLibrary == "GEX"){
      if(param$scMode == "SC"){
        refDirs <- list.files(path=refDir, pattern="^10X_Ref.*_GEX_", 
                              full.names = TRUE)
      }else if(param$scMode == "SN"){
        refDirs <- list.files(path=refDir, pattern="^10X_Ref.*_premRNA_", 
                              full.names = TRUE)
      }
    }else{
      stop("Unsupported 10X library: ", param$TenXLibrary)
    }
    
    if(length(refDirs) == 0){
      stop("No 10X_Ref folder found in", refDir)
    }
    if(length(refDirs) > 1){
      warning("Multiple 10X_Ref folders in ", refDir)
    }
    refDir <- refDirs[1]
  }
  return(refDir)
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodCellRanger(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppCellRanger <-
  setRefClass("EzAppCellRanger",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCellRanger
                  name <<- "EzAppCellRanger"
                  appDefaults <<- rbind(TenXLibrary=ezFrame(Type="charVector",
                                                            DefaultValue="GEX",
                                                            Description="Which 10X library? GEX or VDJ."),
                                        scMode=ezFrame(Type="character",
                                                       DefaultValue="SC",
                                                       Description="Single cell or single nuclei?"),
                                        chemistry=ezFrame(Type="character",
                                                          DefaultValue="auto",
                                                          Description="Assay configuration."),
                                        controlSeqs=ezFrame(Type="charVector",
                                                            DefaultValue="",
                                                            Description="control sequences to add"))
                }
              )
  )
