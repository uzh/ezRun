###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

#' RunMethod for reference class EzAppDEXSeqAnalysis
#' 
#' @description 
#' Differential exon usage is assessed using the function \code{DEXSeq}
#' from the DEXSeq package.
#' 
#' @param input either EzDataSet reference object or path to input dataset file
#' @param output 
#' @param param
#' 
ezMethodDEXSeqAnalysis <- function(input=NA, output=NA, param=NA){
  param[['BPPARAM']] = BiocParallel::MulticoreParam(workers=param$cores)
  ### # get count files based on the name of the bamfiles  
  sCountfileExt <- 'count'
  if (ezIsSpecified(param$countfile_ext))
    sCountfileExt <- param$countfile_ext
  countFiles <- gsub("bam$", replacement = sCountfileExt, basename(input$getColumn("BAM")))
  ### # if count files do not exist, generate them
  if(!all(file.exists(countFiles)))
    DEXSeqCounting(input = input, output = output, param = param)
    #EzAppDEXSeqCounting$new()$run(input = input, output = output, param = param)
  
  ### # check whether conditions are specified
  colnames(input$meta) = gsub(' \\[.*','',colnames(input$meta))
  if (param$grouping %in% colnames(input$meta))
    condition <- input$meta[[param$grouping]]
  #if (ezIsSpecified(param$grouping))
  #  condition <- param$grouping
  ### # if conditions are not specified, then we have to stop here
  if (is.null(condition))
    stop(" * No conditions were specified in ezMethodDEXSeqAnalysis")
  ### # row indices of samples and reference
  vIdxSample <- which(condition == param$sampleGroup)
  vIdxRef <- which(condition == param$refGroup)
  ### # define order, from the vignette, it seams that 
  ### #  first come the sample rows then the reference rows
  vCompOrder <- c(vIdxSample, vIdxRef)
  if(length(countFiles)>length(vCompOrder)){
    countFiles = countFiles[vCompOrder]
  }
  ### # sample table from rownames and conditions
  sampleTable <- data.frame(
    row.names = rownames(input$meta)[vCompOrder],
    condition = condition[vCompOrder]
  )
  ### # check the reference
  sRefFeatGff <- gsub("gtf$", "gff", basename(param[['ezRef']]@refFeatureFile))
  if(ezIsSpecified(param$gff_file))
    sRefFeatGff <- param$gff_file
  stopifnot(file.exists(sRefFeatGff))
  
  ### # check whether special design was specified, o/w use minimal default design
  if (ezIsSpecified(param$design)) {
    design <- param$design
  } else {
    design <- ~ sample + exon + condition:exon
    param$design <- design
  }
    
  ### # create the initial DEXSeqDataSet object
  dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData    = sampleTable,
    design        = design,
    flattenedfile = sRefFeatGff )
  
  ### # estimate size factors and dispersion
  dxd <- DEXSeq::estimateSizeFactors( dxd )
  dxd <- DEXSeq::estimateDispersions( dxd, BPPARAM = param[['BPPARAM']] )

  ### # testing for differential usage
  dxd  <- DEXSeq::testForDEU( dxd, BPPARAM = param[['BPPARAM']] )
  
  ### # fold changes
  dxd <- DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar = tolower(param$grouping), BPPARAM = param[['BPPARAM']])

  ### # generate a report
  writeDEXSeqReport(dataset = input$meta, dexResult = list(param = param, dxd=dxd), output = output, sResultDir = basename(output$meta[['Report [File]']]))
  return("Success")  
}


##' @template app-template
##' @templateVar method ezMethodDEXSeqCounting
##' @templateVar htmlArg, htmlFile="00index.html" )
##' @description Use this reference class to run analysis on differential exon usage
EzAppDEXSeqAnalysis <- 
  setRefClass(Class = "EzAppDEXSeqAnalysis",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDEXSeqAnalysis
                  name <<- "EzAppDEXSeqAnalysis"
                  appDefaults <<- rbind(disp_plot      = ezFrame(Type="character", DefaultValue="dispersion_estimate_plot", Description="which test method in DESeq to use: deseq2"),
                                        ma_plot        = ezFrame(Type="character", DefaultValue="ma_plot",    Description="no need to compute moderated ratios; deseq2 does this already"),
                                        countfile_ext  = ezFrame(Type="character", DefaultValue="count",      Description="extension of count files"),
                                        countfile_path = ezFrame(Type="character", DefaultValue=".",          Description="path where count files should be stored"),
                                        gff_file       = ezFrame(Type="character", DefaultValue="genes.gff",  Description="name of the gff annotation file"),
                                        fdr            = ezFrame(Type="numeric",   DefaultValue=0.1,          Description="false discovery rate below which genes are reported"))
                }
              ))


#' Addition experimental conditions to input files
#'
addDEXSeqCondition = function(psInput, pvCondition){
  ezObjInput <- EzDataset$new(file=psInput)
  ezObjInput$meta$Condition <- pvCondition
  write.table(ezObjInput$meta, file = ezObjInput$file, quote = FALSE, sep = "\t")
}


#' @title Writing a report for a DEXSeq analysis
#' 
#' @description 
writeDEXSeqReport <- function(dataset, dexResult, output, htmlFile="00index.html", types=NULL, sResultDir = "html") {
  ### # retrieve parameters 
  param <- dexResult$param
  ### # extract name appearing in the report
  name <- param$name
  ### # result dataframe and generate html results to be included later
  dxd <- dexResult$dxd
  dxr <- DEXSeq::DEXSeq(dxd)

  ### # write tsv file from results
  if (ezIsSpecified(param$ResultFile)) {
    sResultFile <- param$ResultFile
  } else {
    sResultFile <- "DexSeqResult.tsv"
  }
  write.table(dxr, file = sResultFile, quote = FALSE, sep = "\t")
  
  ### # put the results into a different subdirectory
  sCurWd <- getwd()
  setwdNew(sResultDir)

  ### # write that generic report for a given FDR, using 0.1 as the default
#   nFdr <- 0.1
#   if (ezIsSpecified(param$fdr))
    nFdr <- param$fdr
    DEXSeq::DEXSeqHTML(dxr, FDR = nFdr, BPPARAM = param[['BPPARAM']])

  ### # put a title to the report using name in output
  titles <- list()
  titles[["Analysis"]]  <- paste("Analysis:", name)
  ### # create a report instance
  doc <- openBsdocReport(title=titles[[length(titles)]])
  ### # adding the dataset meta information
  addDataset(doc, dataset, param)
  
  ### # result summary
  titles[["Result Summary"]] = "Result Summary"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  settings = character()
  settings["Grouping:"] = param$grouping
  settings["Sample group:"] = param$sampleGroup
  settings["Reference group:"] = param$refGroup
  settings["Design:"] = paste(as.character(param$design), collapse = " ")
  settings["FDR:"] = as.character(nFdr)
  settings["Number of result features:"] = nrow(dxr)
  addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
  
  ### # experimental design
  titles[["Experimental Design"]] = "Experimental Design"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  ### # put together experimental condition
  sampleData <- dxr@sampleData
  fitExpToVar <- tolower(param$grouping)
  numcond <- length(unique(sampleData[[fitExpToVar]]))
  cond <- as.data.frame(sampleData[, !colnames(sampleData) %in% "sizeFactor"])
  addFlexTable(doc, ezFlexTable(cond, add.rownames=FALSE, header.columns = TRUE))
  
  ### # Dispersion plot
  titles[["Dispersion-Plot"]] = "Dispersion Plot"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  if (ezIsSpecified(param$disp_plot)) {
    sDispPlotFile <- ifelse(identical(tools::file_ext(param$disp_plot), "png"), param$disp_plot, paste(param$disp_plot, "png", sep = "."))
    addParagraph(doc, 
                 ezImageFileLink(plotCmd = expression(DESeq2::plotDispEsts( dxd )), 
                                 file=sDispPlotFile, 
                                 name="Dispersion Plot",
                                 mouseOverText = "Dispersion Plot"))
  }
  
  ### # MA-Plot
  titles[["MA-Plot"]] = "MA Plot"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  if (ezIsSpecified(param$ma_plot)) {
    sMaPlotPdfFile <- ifelse(identical(tools::file_ext(param$ma_plot), "png"), param$ma_plot, paste(param$ma_plot, "png", sep = "."))
    addParagraph(doc, 
                 ezImageFileLink(plotCmd = expression(DEXSeq::plotMA( dxr )), 
                                 file=sMaPlotPdfFile, 
                                 name="MA Plot",
                                 mouseOverText = "MA Plot"))
  }
  
  ### # Put simply a link to the already existing report
  titles[["DEXSeq differential exon usage test"]] = "DEXSeq differential exon usage test"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addParagraph(doc, pot("Test results for differential exon usage", hyperlink = "DEXSeqReport/testForDEU.html"))
                                                                                                
  ### # closing the report leads to writing it to htmlFile
  closeBsdocReport(doc, htmlFile, titles)
  setwd(sCurWd)
}

DEXSeqCounting <- function(input = input, output = output, param = param){
  ### # check whether GFF formatted annotation is available
  sGtfFile <- param$ezRef@refFeatureFile
  ### # gff will be placed in actual working directory, hence no 
  ### #  soft links will be required
  sGffFile <- gsub("gtf$", "gff", basename(sGtfFile))
  if(ezIsSpecified(param$gff_file))
    sGffFile <- param$gff_file
  if (!file.exists(sGffFile))
    convertGtfToGff(psGtfFile = sGtfFile, psGffFile = sGffFile)
  
  ### # do the counting, get the bam files from input
  bamFiles = as.list(input$getFullPaths(param, "BAM"))
  
  ### # determine extension for count files
  sCountfileExt <- 'count'
  if (ezIsSpecified(param$countfile_ext))
    sCountfileExt <- param$countfile_ext
  ### # call counting routine
  vCountFiles <- ezMclapply(bamFiles, runCountSingleBam, sGffFile, sCountfileExt,mc.cores = param[['Cores']])
  
  return("Success")
}

#' Convert annotation file from GTF format to GFF
#' 
#' @description 
#' \code{convertGtfToGff} converts an annotation file from 
#' the GTF format into the GFF format which is required 
#' by the package \code{DEXSeq}. Input file name and the 
#' name of the file to be generated are both given as 
#' function parameters. The conversion is done by a python 
#' script that is given by the content of \code{DEXSEQ_PREPARE} 
#' which is either taken as a global variable or from the 
#' result of function \code{lGetPyScriptPaths}
#' 
#' @param psGtfFile   name of the GTF annotation file
#' @param psGffFile   name of the GFF file to be generated
convertGtfToGff <- function(psGtfFile, psGffFile) {
  cat(" * Converting GTF to GFF ...\n")
  ### # check whether the path exists
  if (!exists("DEXSEQ_PREPARE")) {
    DEXSEQ_PREPARE <- lGetPyScriptPaths()$DEXSEQ_PREPARE
  }
  sPyConvCmd <- paste(DEXSEQ_PREPARE, psGtfFile, psGffFile)
  ezSystem(sPyConvCmd)
  cat("  ==> created: ", psGffFile, "\n")
  invisible(TRUE)
}

#' Run counts for a single BAM file
#' 
runCountSingleBam <- function(psBamFile, psGffFile, psCountfileExt){
  sSamCmd <- paste(SAMTOOLS, "view -h", psBamFile)
  ### # run counting on sam file
  sCountBaseFn <- gsub("bam$", psCountfileExt, basename(psBamFile))
  if (!exists("DEXSEQ_COUNT")){
    DEXSEQ_COUNT <- lGetPyScriptPaths()$DEXSEQ_COUNT
  }
  sPyCountCmd <- paste(sSamCmd, "|", DEXSEQ_COUNT, psGffFile, "-", sCountBaseFn)
  ezSystem(sPyCountCmd)
  sCountDir <- getwd()
  return(file.path(sCountDir, sCountBaseFn))
}


#' Write names of countfiles back into the input file
#' 
writeCountFilesToMeta <- function(pvCountFiles, input) {
  ### # add column with counts to the meta information
  input$meta$Count <- pvCountFiles
  ### # write the extended meta information back to the file
  write.table(input$meta, file = input$file, quote = FALSE, sep = "\t")
}


#' Get list with required python script paths
#' 
lGetPyScriptPaths <- function(){
  if (!exists("PYTHON_CMD")){
    PYTHON_CMD='PYTHONPATH="/usr/local/ngseq/lib/python/:/usr/local/ngseq/lib/python2.7:/usr/local/ngseq/lib/python2.7/dist-packages" /usr/local/ngseq/bin/python'
  }
  return(list(DEXSEQ_PREPARE = paste(PYTHON_CMD, file.path(system.file(package = "DEXSeq", "python_scripts"), "dexseq_prepare_annotation.py")),
              DEXSEQ_COUNT = paste(PYTHON_CMD, file.path(system.file(package = "DEXSeq", "python_scripts"), "dexseq_count.py"))))
  
}
