###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

#' RunMethod for reference class EzAppDEXSeqAnalysis
#' 
#' @description 
#' Differential exon usage is assessed using the same steps that are done in 
#' the wrapper function \code{DEXSeq} from the DEXSeq package. The wrapper 
#' is not used here, because dispersion plots and MA-plots are integrated 
#' in the report and those cannot be produced from the resulting object of 
#' the wrapper function. The following steps are done in the analysis
#' \itemize{
#'   \item{DEXSeqDataSetFromHTSeq: }{instantiates a DEXSeqDataSet object}
#'   \item{estimateSizeFactors: }{normalizes to account for different coverage}
#'   \item{estimateDispersions: }{assesses variability within and between experimental groups}
#'   \item{testForDEU: }{getting fdr statistics for deu}
#'   \item{estimateExonFoldChanges: }{fold changes are estimated}
#' }
#' Once these steps are done the results are written to a report.
#' 
#' @details 
#' Before running the DEU-analysis, it is verified that count data 
#' is available. It is expected that the current working directory 
#' contains files with the same name as the BAM-files specified in the 
#' input, but with a file-extension that can be specified and which 
#' defaults to 'count'. In case, the count files are not available 
#' they are generated using the function \code{DEXSeqCounting}.
#' 
#' @param input EzDataSet reference object specifying input
#' @param output EzDataSet reference object specifying output 
#' @param param EzParam reference object specifying additional parameters
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
  
  sampleTable <- data.frame(
    row.names = rownames(input$meta),
    condition = condition
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
  ### # define the reference group in the condition levels
  dxd$condition <- relevel(dxd$condition, param$refGroup)
  
  ### # estimate size factors and dispersion
  dxd <- DEXSeq::estimateSizeFactors( dxd )
  dxd <- DEXSeq::estimateDispersions( dxd, BPPARAM = param[['BPPARAM']] )

  ### # testing for differential usage
  dxd  <- DEXSeq::testForDEU( dxd, BPPARAM = param[['BPPARAM']] )
  
  ### # fold changes
  dxd <- DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar = tolower(param$grouping), BPPARAM = param[['BPPARAM']])

  ### # generate a report
  writeDEXSeqReport(dataset = input$meta, dexResult = list(param = param, dxd=dxd), sResultDir = basename(output$meta[['Report [File]']]))
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
                                        fdr            = ezFrame(Type="numeric",   DefaultValue=0.1,          Description="false discovery rate below which genes are reported"),
                                        dexseq_report_path = ezFrame(Type="character", DefaultValue="DEXSeqReport",  Description="path DEXSeqHTML report is written to"),
                                        dexseq_report_file = ezFrame(Type="character", DefaultValue="testForDEU.html",  Description="file name for DEXSeqHTML report")   )
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
writeDEXSeqReport <- function(dataset, dexResult, htmlFile="00index.html", sResultDir = "html") {
  ### # retrieve parameters 
  param <- dexResult$param
  ### # extract name appearing in the report
  name <- param$name
  ### # extract DEXSeqResults object
  dxd <- dexResult$dxd
  dxr <- DEXSeq::DEXSeqResults(dxd)

  ### # put the results into a different subdirectory
  sCurWd <- getwd()
  setwdNew(sResultDir)
  
  ### # write tsv file from results
  if (ezIsSpecified(param$ResultFile)) {
    sResultFile <- param$ResultFile
  } else {
    sResultFile <- "DexSeqResult.tsv"
  }
  write.table(dxr, file = sResultFile, quote = FALSE, sep = "\t")
  
  ### # if parameters about biomart are specified, use them
  if (ezIsSpecified(param$bio_mart) & 
      ezIsSpecified(param$mart_dataset) &
      ezIsSpecified(param$mart_filter) & 
      ezIsSpecified(param$mart_attributes) ) {
    ensembl_mart <- biomaRt::useMart(biomart = param$bio_mart, dataset=param$mart_dataset)  
    DEXSeq::DEXSeqHTML(dxr, path = param$dexseq_report_path, file = param$dexseq_report_file, FDR = param$fdr,
                       mart = ensembl_mart, filter = param$mart_filter, attributes = param$mart_attributes)
  } else {
    ### # write that generic report for a given FDR, using 0.1 as the default
    DEXSeq::DEXSeqHTML(dxr, path = param$dexseq_report_path, file = param$dexseq_report_file, FDR = param$fdr)
  }
  
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
  settings["FDR:"] = as.character(param$fdr)
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

  ### # table with annotations
  genetable <- getGeneTable(pdxr = dxr, param = param)
  titles[["DEXSeq differential exon usage test"]] = "DEXSeq differential exon usage test"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addFlexTable(doc, ezFlexTable(genetable, add.rownames=FALSE, header.columns = TRUE))
  
  
  ### # Put simply a link to the already existing report
  titles[["Report generated by DEXSeq"]] = "Report generated by DEXSeq"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addParagraph(doc, pot("Test results for differential exon usage", hyperlink = "DEXSeqReport/testForDEU.html"))
  
  
                                                                                                
  ### # closing the report leads to writing it to htmlFile
  closeBsdocReport(doc, htmlFile, titles)
  setwd(sCurWd)
}


#' Generate gene table from DEXSeqResults object
#' 
#' @param pdxr       DEXSeqResults object
#' @param param      EzParam object
getGeneTable <- function(pdxr, param){
  ### # check that argument is a DEXSeqResults object
  stopifnot(is(pdxr, "DEXSeqResults"))
  
  ### # put together a result data consisting of genomic data and the modelling
  genomicData <- as.data.frame(pdxr$genomicData)
  results <- data.frame(pdxr[, c("groupID", 
                                 "featureID", 
                                 "exonBaseMean", 
                                 "dispersion", 
                                 "pvalue", 
                                 "padj")], 
                        stringsAsFactors = FALSE)
  results <- cbind(results, genomicData)
  results[, c("dispersion", "pvalue", "padj")] <- round(results[, c("dispersion", "pvalue", "padj")], 3)
  dexseqR <- elementMetadata(pdxr)$type == "DEXSeq results"
  if (sum(dexseqR, na.rm = TRUE) > 0) {
    results <- cbind(results, round(as.data.frame(pdxr[, which(dexseqR)]), 3))
  }
  rownames(results) <- NULL

  ### # from the results, generate the genetable which seams to be the basis for the
  ### #  table on the results page
  gns <- as.character(unique(results$groupID[which(results$padj < param$fdr)]))
  results <- results[as.character(results$groupID) %in% gns,]
  splitCols <- split(seq_len(nrow(results)), results$groupID)
  genetable <- lapply(splitCols, function(x) {
    data.frame(chr = unique(results$seqnames[x]), start = min(results$start[x]),
               end = max(results$end[x]), total_exons = length(x),
               exon_changes = sum(results$padj[x] < param$fdr, na.rm = TRUE))
  })
  
  ### # seams to convert the list "genetable" to a data.frame
  genetable <- do.call(rbind, genetable)
  genetable <- cbind(geneID = rownames(genetable), genetable)
  
  
  ### # reading gene annotations from annotation file
  ### #  extract name of annotation file from param
  sGnAnFn <- param[['ezRef']]@refAnnotationFile
  dfGnsAnnot <- ezRead.table(file = sGnAnFn, row.names = NULL)

  ### # extract gene_names and descriptions for all genes in the whole genetable
  gene_name <- sapply(as.character(genetable$geneID), 
                      function(x) {
                        ### # some entries in geneID can contain multiple
                        ### #  gene_ids pasted together with "+"
                        sSingleId <- unlist(strsplit(x, split = "+", fixed = TRUE))
                        nNrIds <- length(sSingleId)
                        if (nNrIds > 1) {
                          sResultGeneName <- paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == sSingleId[1],"gene_name"]), sep = "", collapse = " | ")
                          for(nIdx in 2:nNrIds){
                            sResultGeneName <- paste(sResultGeneName, paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == sSingleId[nIdx],"gene_name"]), sep = "", collapse = " | "),
                                                     sep = "+")
                          }
                        } else {
                          sResultGeneName <- paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == x,"gene_name"]), sep = "", collapse = " | ")
                        }
                        return(sResultGeneName)
                      }, 
                      USE.NAMES = FALSE)
  gene_description <- sapply(as.character(genetable$geneID), 
                             function(x) {
                               ### # some entries in geneID can contain multiple
                               ### #  gene_ids pasted together with "+"
                               sSingleId <- unlist(strsplit(x, split = "+", fixed = TRUE))
                               nNrIds <- length(sSingleId)
                               if (nNrIds > 1) {
                                 sResultGeneDesc <- paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == sSingleId[1],"description"]), sep = "", collapse = " | ")
                                 for(nIdx in 2:nNrIds){
                                   sResultGeneDesc <- paste(sResultGeneDesc, paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == sSingleId[nIdx],"description"]), sep = "", collapse = " | "),
                                                            sep = "+")
                                 }
                               } else {
                                 sResultGeneDesc <- paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == x,"description"]), sep = "", collapse = " | ")
                               }
                               return(sResultGeneDesc)
                             }, 
                             USE.NAMES = FALSE)
  ### # add extracted columns and return genetable
  genetable <- cbind(genetable, gene_name, gene_description)
  
  ### # add links to result files
  genetable$geneID <- as.character(genetable$geneID)
  genetable$geneID <- getGeneIdExprLinks(pvGeneIds = genetable$geneID, psdexseq_report_path = param$dexseq_report_path)
  
  return(genetable)
}


#' Get Links from GeneIds to Expression result files
#' 
#' @param pvGeneIds              vector of gene ids
#' @param psdexseq_report_path   path to where DEXSeqHTML report is saved
getGeneIdExprLinks <- function(pvGeneIds, psdexseq_report_path){
  ### # elements in pvGeneIds can have multiple GeneIds separated with a "+"
  ### # the following local function will take a single element in pvGeneIds 
  ### # split it up, if required and generate the link to the result file
  getGeneLink <- function(psGeneId) {
    vSplitIds <- unlist(strsplit(psGeneId, split = "+", fixed = TRUE))
    nrIds <- length(vSplitIds)
    if (nrIds > 1) {
      ### # in case there are multiple GeneIds in psGeneId, the first determines the link, hence we save it
      sLinkToFirstGene <- paste0(psdexseq_report_path, "/files/", vSplitIds[1], "expression.html")
      sResultLink <- as.html(pot(vSplitIds[1], hyperlink = sLinkToFirstGene))
      for (nIdx in 2:nrIds){
        sResultLink <- paste(sResultLink, 
                             as.html(pot(vSplitIds[nIdx], hyperlink = sLinkToFirstGene)),
                             sep = "+")
      }
    } else {
      sResultLink <- as.html(pot(psGeneId, hyperlink = paste0(psdexseq_report_path, "/files/", psGeneId, "expression.html")))
    }
    return(sResultLink)
  }
  return(sapply(pvGeneIds, getGeneLink, USE.NAMES = FALSE))
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
  bamFiles = as.list(input$getFullPaths("BAM"))
  
  ### # determine extension for count files
  sCountfileExt <- 'count'
  if (ezIsSpecified(param$countfile_ext))
    sCountfileExt <- param$countfile_ext
  ### # call counting routine
  vCountFiles <- ezMclapply(bamFiles, runCountSingleBam, sGffFile, sCountfileExt,mc.cores = param[['cores']])
  
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
