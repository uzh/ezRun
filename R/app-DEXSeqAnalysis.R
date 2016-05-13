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
  ### # get count files based on the name of the bamfiles  
  sCountfileExt <- 'count'
  if (ezIsSpecified(param$countfile_ext))
    sCountfileExt <- param$countfile_ext
  countFiles <- gsub("bam$", replacement = sCountfileExt, basename(input$getColumn("BAM")))
  ### # if count files do not exist, generate them
  if(!all(file.exists(countFiles)))
    EzAppDEXSeqCounting$new()$run(input = input, output = output, param = param)
  
  ### # check whether conditions are specified
  if (param$grouping %in% colnames(input$meta))
    condition <- input$meta$Condition
  if (ezIsSpecified(param$condition))
    condition <- param$condition
  ### # if conditions are not specified, then we have to stop here
  if (is.null(condition))
    stop(" * No conditions were specified in ezMethodDEXSeqAnalysis")
  ### # row indices of samples and reference
  vIdxSample <- which(condition == param$sampleGroup)
  vIdxRef <- which(condition == param$refGroup)
  ### # define order, from the vignette, it seams that 
  ### #  first come the sample rows then the reference rows
  vCompOrder <- c(vIdxSample, vIdxRef)
  
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
  dxd <- DEXSeq::estimateDispersions( dxd )

  ### # testing for differential usage
  dxd  <- DEXSeq::testForDEU( dxd )
  
  ### # fold changes
  dxd <- DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar=tolower(param$grouping))
  
  ### # generate a report
  writeDEXSeqReport(dataset = input$meta, dexResult = list(param = param, dxd=dxd), output=output)
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
  ### # extract DEXSeqResults object
  dxd <- dexResult$dxd
  dxr <- DEXSeq::DEXSeqResults(dxd)

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
  #nFdr <- param$fdr
  DEXSeq::DEXSeqHTML(dxr, FDR = param$fdr)

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
  dfGnsAnnot <- read.table(file = sGnAnFn, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL)

  ### # extract gene_names and descriptions for all genes in the whole genetable
  gene_name <- sapply(genetable$geneID, 
                      function(x) 
                        paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == x,"gene_name"]), sep = "", collapse = " | "), 
                      USE.NAMES = FALSE)
  gene_description <- sapply(genetable$geneID, 
                             function(x) 
                               paste(unique(dfGnsAnnot[dfGnsAnnot[, "gene_id"] == x,"description"]), sep = "", collapse = " | "), 
                             USE.NAMES = FALSE)
  ### # add extracted columns and return genetable
  genetable <- cbind(genetable, gene_name, gene_description)
  return(genetable)
}

