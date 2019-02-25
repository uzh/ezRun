###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##ToD: --localmem, --localcores, use opt-Parameters
ezMethodCellRanger = function(input=NA, output=NA, param=NA){
  sampleName = input$getNames()
  sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
  sampleDirs <- file.path(input$dataRoot, sampleDirs)
  sampleDir <- paste(sampleDirs, collapse=",")
  refDir <- getCellRangerReference(param)
  if(param$TenXLibrary == "GEX"){
    cmd <- paste(CELLRANGER, "count", paste0("--id=", sampleName),
                 paste0("--transcriptome=", refDir),
                 paste0("--fastqs=", sampleDir),
                 paste0("--sample=", sampleName),
                 paste0("--localmem=", param$ram),
                 paste0("--localcores=", param$cores),
                 paste0("--chemistry=", param$chemistry))
  }else if(param$TenXLibrary == "VDJ"){
    cmd <- paste(CELLRANGER, "vdj", paste0("--id=", sampleName),
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
  
  if(ezIsSpecified(param$controlSeqs)){
    unlink(refDir, recursive = TRUE)
  }
  
  if(param$TenXLibrary == "GEX" &&
     param$refBuild %in% c("Homo_sapiens/Ensembl", "Mus_musculus/Ensembl")){
    require(DropletUtils)
    countMatrixFn <- list.files(path=sampleName,
                                pattern="\\.mtx(\\.gz)*$", recursive=TRUE,
                                full.names=TRUE)
    sce <- read10xCounts(dirname(countMatrixFn), col.names=TRUE)
    cellPhase <- getCellCycle(sce, param$refBuild)
    write_tsv(cellPhase,
              path=file.path(dirname(countMatrixFn), "CellCyclePhase.txt"))
  }
  
  return("Success")
}


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

## Make 10X reference
### VDJ 
#### /usr/local/ngseq/opt/cellranger-3.0.1/cellranger mkvdjref --genome=10X_Ref_Mouse_VDJ_GRCm38.p5_20192019_Release_91 --fasta=../../../Sequence/WholeGenomeFasta/genome.fa --genes=genes.gtf


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
