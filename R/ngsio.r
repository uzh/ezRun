###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

loadCountDataset <- function(input, param){
  require(tools)
  require(SummarizedExperiment)
  require(readr)
  require(dplyr)
  require(readr)
  files = input$getFullPaths("Count")
  
  dataFeatureLevel = unique(input$getColumn("featureLevel"))
  stopifnot(length(dataFeatureLevel) == 1)
  
  x1 <- read_tsv(files[1], guess_max=1e6)
  
  if (ezIsSpecified(param$expressionName)){
    columnName = param$expressionName
  } else {
    columnName = intersect(param$knownExpressionNames, colnames(x1))[1]
  }
  if (!columnName %in% colnames(x1)){
    return(list(error=paste0("Specified column name not found in data!<br>columnName: '", columnName, "'\n",
                             "<br>Available column names:<br>\n",
                             paste0("'", colnames(x1), "'", collapse="<br>"),
                             "<br>Set the option columnName to one of the names above!")))
  }
  identifier <- 1 #colnames(x1)[1]
  
  x <- mapply(function(x, y){
    message("loading file: ", x)
    tempTibble <- read_tsv(x, progress=FALSE, guess_max=1e6)
    #stopifnot(setequal(pull(tempTibble[1]), pull(x1[1])))
    tempTibble %>%
      dplyr::select(identifier, columnName) %>%
      dplyr::rename("id":= 1, !! y := columnName)
  }, files, names(files), SIMPLIFY=FALSE)
  x <- Reduce(function(x,y){left_join(x, y, by="id")}, x)
  
  if(dataFeatureLevel == "isoform" && param$featureLevel == "gene"){
    ## aggregate from isoform to gene level
    seqAnnoDFData <- ezFeatureAnnotation(param, pull(x["id"]), 
                                         dataFeatureLevel)
    stopifnot(identical(seqAnnoDFData$transcript_id, x[[identifier]]))
    
    x$gene_id <- seqAnnoDFData$gene_id
    x <- select(x, -id) %>% group_by(gene_id) %>% 
      summarise_all(funs(sum))
    ## TODO: consider using rowsum()
  }
  signal <- as.matrix(x[ ,-1])
  rownames(signal) <- pull(x[ ,1])

  if(ezIsSpecified(param$ezRef@refBuild)){
    seqAnnoDFFeature <- ezFeatureAnnotation(param, rownames(signal),
                                            param$featureLevel)
  }else{
    seqAnnoDFFeature <- data.frame(row.names = rownames(signal),
                                   gene_id=rownames(signal),
                                   transcript_id=rownames(signal),
                                   gene_name=rownames(signal),
                                   type="protein_coding",
                                   strand="*",
                                   seqid=1,
                                   biotypes="protein_coding",
                                   description=rownames(signal),
                                   start=1, end=100, gc=NA, featWidth=NA,
                                   "GO BP"="", "GO CC"="", "GO MF"="")
    if(any(colnames(x1) %in% c('GeneID', 'Chr', 'Start', 'End', 'Strand'))){
      x1 <- rename(x1, gene_id=GeneID, seqid=Chr, start=Start, end=End,
                 strand=Strand)
    }
    colsInData <- intersect(colnames(seqAnnoDFFeature), colnames(x1))
    if(length(colsInData) > 0){
      seqAnnoDFFeature[ ,colsInData] <- x1[ ,colsInData]
    }
    if(!all(seqAnnoDFFeature$strand %in% c("+", "-", "*"))){
      wrongStrands <- ! seqAnnoDFFeature$strand %in% c("+", "-", "*")
      seqAnnoDFFeature$strand[wrongStrands] <- "*"
    }
  }
  stopifnot(identical(rownames(seqAnnoDFFeature), rownames(signal)))
  
  if (ezIsSpecified(param$correctBias) && param$correctBias){
    ## output will be floating point, but we don't round; input might already be floating point
    signal = ezCorrectBias(signal, gc=seqAnnoDFFeature$gc,
                           width=seqAnnoDFFeature$featWidth)$correctedCounts
  }
  
  seqAnno <- makeGRangesFromDataFrame(seqAnnoDFFeature, keep.extra.columns=TRUE)
  
  if (param$useSigThresh){
    sigThresh = param$sigThresh
  } else {
    sigThresh = 0
  }
  
  ## assays: counts, presentFlag, RPKM, TPM, (signal)
  ## rowData, colData
  ## meta: isLog, featureLevel, type, countName, param
  rawData <- SummarizedExperiment(
    assays=SimpleList(counts=signal, presentFlag=signal > sigThresh),
    rowRanges=seqAnno, colData=input$meta,
    metadata=list(isLog=FALSE, featureLevel=param$featureLevel,
                  type="Counts", countName=columnName,
                  param=param)
    )
  
  if (ezIsSpecified(param$transcriptTypes)){
    use = seqAnno$type %in% param$transcriptTypes
  } else {
    use = TRUE
  }
  rawData <- rawData[use, ]

  assays(rawData)$rpkm = getRpkm(rawData)
  assays(rawData)$tpm = getTpm(rawData)
  return(rawData)
}

loadSCCountDataset <- function(input, param){
  require(SingleCellExperiment)
  require(SummarizedExperiment)
  require(DropletUtils)
  require(tools)
  require(Matrix)

  if(length(input$getNames()) > 1L){
    ## Recursively load and combine sce
    sce <- SingleCellExperiment::cbind(loadSCCountDataset(input$subset(input$getNames()[1]), param),
                                       loadSCCountDataset(input$subset(input$getNames()[-1]), param))
    return(sce)
  }

  if(any(grepl("___", input$getNames()))){
    stop("___ is not allowed in sample name.")
  }
  
  dataFeatureLevel <- unique(input$getColumn("featureLevel"))
  stopifnot(length(dataFeatureLevel) == 1)
  
  if(toupper(param$scProtocol) == "SMART-SEQ2"){
    countMatrixFn <- input$getFullPaths("CountMatrix")
    if(file_ext(countMatrixFn) == "mtx"){
      countMatrix <- readSCMM(countMatrixFn)
    }else if(file_ext(countMatrixFn) == "txt"){
      countMatrix <- Matrix(as.matrix(ezRead.table(countMatrixFn)))
    }
    cellDataSet <- ezRead.table(sub("-counts\\.mtx$", "-dataset.tsv",
                                    countMatrixFn))
    
    cellCycleFn <- sub("-counts\\.mtx$", "-CellCyclePhase.txt",
                       countMatrixFn)
    if(file.exists(cellCycleFn)){
      ## TODO: remove this test as CellCyclePhase file should always exist
      cellCycle <- ezRead.table(cellCycleFn)
      cellDataSet$CellCycle <- cellCycle[rownames(cellDataSet), "Phase"]
    }
    
    ## TODO: this is a temporary solution to fix the discrepency of sample names
    if(!setequal(colnames(countMatrix), rownames(cellDataSet))){
      ## fix the colnames in countMatrix
      ## It's only possible if sample names are part of fastq file names
      matches <- lapply(rownames(cellDataSet), grep, colnames(countMatrix), 
                        fixed=TRUE, value=FALSE)
      stopifnot(all(lengths(matches) <= 1L)) ## stop when mutiple matches happen
      cellDataSet = cellDataSet[sapply(matches, length) == 1, ]
      matches <- sapply(rownames(cellDataSet), grep, colnames(countMatrix), 
                        fixed=TRUE, value=FALSE)
      colnames(countMatrix)[unlist(matches)] = names(matches)
      ## TODO this is a temporary solution to support cells with zero reads
    }
    ## Reorder the countMatrix columns
    ## This should be unnecessary if we retain the order of RG when creating unmapped bam
    countMatrix <- countMatrix[ , rownames(cellDataSet)]
    
    seqAnnoDF <- ezFeatureAnnotation(param, rownames(countMatrix),
                                     dataFeatureLevel)
    seqAnno <- makeGRangesFromDataFrame(seqAnnoDF, keep.extra.columns=TRUE)
    
    if (ezIsSpecified(param$correctBias) && param$correctBias){
      countMatrix <- ezCorrectBias(countMatrix, gc = seqAnno$gc, 
                                   width=seqAnno$featWidth)$correctedCounts
    }
    
    sce <- SingleCellExperiment(assays=list(counts=countMatrix),
                                rowRanges=seqAnno, colData=cellDataSet,
                                metadata=list(isLog=FALSE,
                                              featureLevel=dataFeatureLevel,
                                              type="Counts", param=param))
  }else if(toupper(param$scProtocol) == "10X"){
    countMatrixFn <- list.files(path=input$getFullPaths("CountMatrix"),
                                pattern="\\.mtx(\\.gz)*$", recursive=TRUE, 
                                full.names=TRUE)
    sce <- read10xCounts(dirname(countMatrixFn), col.names=TRUE)
    cellCycleFn <- list.files(path=input$getFullPaths("CountMatrix"),
                              pattern="CellCyclePhase\\.txt$", recursive=TRUE,
                              full.names=TRUE)
    if(length(cellCycleFn) == 1){ 
      ## TODO: remove this test as CellCyclePhase file should always exist
      cellCycle <- ezRead.table(cellCycleFn)
      colData(sce)$CellCycle <- cellCycle[colnames(sce), "Phase"]
      colData(sce)$CellCycleG1 <- cellCycle[colnames(sce), "G1"]
      colData(sce)$CellCycleS <- cellCycle[colnames(sce), "S"]
      colData(sce)$CellCycleG2M <- cellCycle[colnames(sce), "G2M"]
    }
    seqAnnoDF <- ezFeatureAnnotation(param, rownames(sce),
                                     dataFeatureLevel)
    seqAnno <- makeGRangesFromDataFrame(seqAnnoDF, keep.extra.columns=TRUE)
    rowRanges(sce) <- seqAnno
    metadata(sce)$param <- param
  }else{
    stop("Unsupported single cell protocol!")
  }
  
  if (ezIsSpecified(param$transcriptTypes)){
    use = seqAnno$type %in% param$transcriptTypes
    sce <- sce[use, ]
  }
  
  if (dataFeatureLevel == "isoform" && param$featureLevel == "gene"){
    sce <- aggregateCountsByGene(sce)
  }
  ## unique cell names when merging two samples
  colnames(sce) <- paste(input$getNames(), colnames(sce), sep="___")
  
  colData(sce)$Batch <- input$getNames()
  
  return(sce)
}

##' @title Writes the head of a file
##' @description Writes the head of a file into a newly created target file.
##' @param target a character specifying the path of the output.
##' @param x a character specifying the path to a file to read from.
##' @param n an integer specifying the amount of lines to return.
##' @template roxygen-template
##' @return Returns the name of the output file.
##' @examples 
##' description = system.file("DESCRIPTION", package="ezRun", mustWork=TRUE)
##' ezHead(x=description)
ezHead = function(target=paste0(x, "_head"), x, n=1000){
  writeLines(readLines(x, n=n), con=target)
  return(target)
}

##' @title Reads NcPro results
##' @description Reads NcPro results and performs some checks on the data.
##' @param ncproDir a character representing the file path to the directory that contains the NcPro results.
##' @template roxygen-template
##' @return Returns the rawdata read from the NcPro results.
readNcProResult = function (ncproDir) {
  dataFiles = list.files(ncproDir, "_subfamcov.data$", full.names=TRUE)
  dataFiles = grep("precursor", dataFiles, invert = TRUE, value = TRUE)
  
  df = dataFiles[1]
  allData = NULL
  allTypes = character()
  for (df in dataFiles){
    type = sub("_all_samples_subfamcov.data", "", basename(df))
    x = ezRead.table(df)
    if (is.null(allData)){
      allData = x
    } else {
      stopifnot(colnames(allData) == colnames(x))
      #stopifnot(!rownames(x) %in% rownames(allData))
      allData = rbind(allData, x)
    }
    allTypes = c(allTypes, rep(type, nrow(x)))
  }
  seqAnno = data.frame(type=I(allTypes), row.names=rownames(allData))
  typeNames = c( "rfam_ACA_snoRNA"="ACA_snoRNA",
                 "rfam_CD_snoRNA"="CD_snoRNA",
                 "rmsk_L1"="L1",
                 "tRNA_tRNA"="tRNA",
                 "tRNA_tRNA_e_+30_+30"="tRNA_e",
                 "precursor_miRNA_miRNA"="pre-miRNA",
                 "mature_miRNA_miRNA_e_+2_+2"="miRNA")
  types = data.frame(row.names=rownames(seqAnno))
  stopifnot(allTypes %in% names(typeNames))
  for (nm in names(typeNames)){
    types[[typeNames[nm]]] = seqAnno$type == nm
  }
  rawData = list(counts=as.matrix(allData), presentFlag = allData > 0, seqAnno=seqAnno, featureLevel="isoform", type="Counts", countName="ncpro-counts", isLog=FALSE,
                 types=types)
}

##' @title Filters FastQ files by bam
##' @description Filters FastQ files by bam and writes them into new files.
##' @param fqFiles a character vector representing file paths to FastQ files.
##' @template bamFile-template
##' @param fqOutFiles an optional character vector representing file paths to the FastQ output files.
##' @param doGzip a logical indicating whether to archive the output files in a gzip archive.
##' @param keepUnmapped passed further to \code{scanBamFlag()}.
##' @param isProperPair passed further to \code{scanBamFlag()}.
##' @param readIds a character containing the read ID's from the \code{bamFile}.
##' @template roxygen-template
##' @seealso \code{\link[Rsamtools]{scanBam}}
##' @seealso \code{\link[Rsamtools]{ScanBamParam}}
##' @return Returns the names of the filtered FastQ files.
##' @examples 
##' bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
##' file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
##' param = list()
##' param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
##' input = EzDataset$new(file=file, dataRoot=param$dataRoot)
##' fqFiles = input$getFullPaths("Read1")
##' fqFiltered = filterFastqByBam(fqFiles, bamFile, doGzip = FALSE)
filterFastqByBam = function(fqFiles, bamFile, fqOutFiles=NULL, doGzip=TRUE, keepUnmapped=TRUE, isProperPair=NA){
  param = ScanBamParam(what=c("qname"),
                     flag=scanBamFlag(isUnmappedQuery=!keepUnmapped, isProperPair=isProperPair))
  bamReads = unique(scanBam(bamFile, param=param)[[1]]$qname)
  gc()
  fqFiltered = removeReadsFromFastq(fqFiles, bamReads, fqOutFiles=fqOutFiles, doGzip=doGzip)
  return(fqFiltered)
}

##' @describeIn filterFastqByBam Removes reads from the FastQ files according to the bam filter.
removeReadsFromFastq = function(fqFiles, readIds, fqOutFiles=NULL, doGzip=TRUE){
  require("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  if (is.null(fqOutFiles)){
    fqOutFiles = sub(".gz$", "", basename(fqFiles))
    fqOutFiles = sub("R1.fastq", "clean_R1.fastq", fqOutFiles)
    fqOutFiles = sub("R2.fastq", "clean_R2.fastq", fqOutFiles)
  }
  for (i in 1:length(fqFiles)){
    fqs = FastqStreamer(fqFiles[i], 1e6)
    count = 0
    while (length(x <- yield(fqs))) {
      count = count + 1
      message(fqOutFiles[i], " ", count)
      ids = sub(" .*", "", as.character(id(x)))
      keep = !ids %in% readIds
      writeFastq(x[keep], file=fqOutFiles[i], mode="a", compress=FALSE)
    }
    if (doGzip){
      ezSystem(paste("gzip", fqOutFiles[i]))
      fqOutFiles[i] = paste0(fqOutFiles[i], ".gz")
    }
  }
  return(fqOutFiles)
}

writeSCMM <- function(x, file){
  require(Matrix)
  require(readr)
  
  stopifnot(grepl("\\.mtx$", file))
  
  write_lines(rownames(x), path=sub("\\.mtx$", ".rowNames", file))
  write_lines(colnames(x), path=sub("\\.mtx$", ".colNames", file))

  writeMM(Matrix(x), file=file)
  
  invisible(file)
}

readSCMM <- function(file){
  require(Matrix)
  require(readr)
  ans <- readMM(file)
  ans <- as(ans, "dgCMatrix")
  colnames(ans) <- read_lines(sub("\\.mtx$", ".colNames", file))
  rownames(ans) <- read_lines(sub("\\.mtx$", ".rowNames", file))
  return(ans)
}

saveExternalFiles = function(sce, ...) {
  tr_cnts <- expm1(logcounts(sce))
  geneMeans <- rowsum(t(as.matrix(tr_cnts)), group=colData(sce)[,"ident"])
  geneMeans <- sweep(geneMeans, 1, STATS=table(colData(sce)[,"ident"])[rownames(geneMeans)], FUN="/")
  geneMeans <- log1p(t(geneMeans))
  colnames(geneMeans) <- paste("cluster", colnames(geneMeans), sep="_")
  geneMeanPerClusterFn = "gene_means_per_cluster.txt"
  ezWrite.table(geneMeans, geneMeanPerClusterFn)
  
  geneMeans <- Matrix::rowMeans(tr_cnts)
  geneMeans <- log1p(geneMeans)
  geneMeansFn = "gene_means.txt"
  ezWrite.table(geneMeans, geneMeansFn)
  
  tSNE_data <- as_tibble(reducedDims(sce)$TSNE,
                         rownames="cells")
  tSNE_data <- dplyr::rename(tSNE_data, X=`tSNE_1`, Y=`tSNE_2`)
  tSNE_data$cluster <- colData(sce)[,"ident"]
  tSNEFn = "tSNE_data.tsv"
  write_tsv(tSNE_data, path=tSNEFn)
  
  add_results = list(...)
  for(i in 1:length(add_results[[1]])) {
    if(!is.null(add_results[[1]][[i]])) {
      file_name = names(add_results[[1]][i])
      write_tsv(add_results[[1]][[i]], path=paste0(file_name, ".tsv"))
    }
  }
}