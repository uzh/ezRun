###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##ToD: --localmem, --localcores, use opt-Parameters
ezMethodCellRanger = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  sampleName = input$getNames()
  sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
  sampleDirs <- file.path(input$dataRoot, sampleDirs)
  
  if(length(sampleDirs) == 1L){
    ## Not merged sample dataset
    sampleDir <- sampleDirs[1]
  }else{
    ## Merged sample dataset
    ## soft link to a temp dir
    sampleDir <- file.path(getwd(), paste0("CellRangerDataDir-", Sys.getpid()))
    dir.create(sampleDir)
    on.exit(unlink(sampleDir, recursive=TRUE), add=TRUE)
    for(eachSampleDir in sampleDirs){
      res <- file.symlink(from=list.files(path=eachSampleDir,
                                          pattern=sampleName, full.names=TRUE),
                          to=sampleDir)
      stopifnot(all(res))
    }
  }
  
  refDir <- getCellRangerReference(param)
  
  cmd = paste(CELLRANGER,"count", paste0("--id=", sampleName),
              paste0("--transcriptome=", refDir),
              paste0("--fastqs=", sampleDir),
              paste0("--sample=", sampleName),
              paste0("--localmem=",param$ram),
              paste0("--localcores=",param$cores))
  if(opt!=''){
    cmd = paste(cmd, opt)
  }
  ezSystem(cmd)
  
  if(ezIsSpecified(param$controlSeqs)){
    unlink(refDir, recursive = TRUE)
  }
  
  return("Success")
}


getCellRangerReference <- function(param){
  if(ezIsSpecified(param$controlSeqs)){
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
    refDirs <- list.files(path=refDir, pattern="^10X_Ref", full.names = TRUE)
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
                  appDefaults <<- rbind(controlSeqs=ezFrame(Type="charVector",
                                                            DefaultValue="",
                                                            Description="control sequences to add"))
                }
              )
  )
