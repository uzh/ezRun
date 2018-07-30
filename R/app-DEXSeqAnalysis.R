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
  require(DEXSeq)
  require(ReporteRs)

  BPPARAM = BiocParallel::MulticoreParam(workers=max(param$cores -1, 1))
  register(BPPARAM) ## register it because exon fold changes does use the default
  
  ### # check whether conditions are specified
  condition = make.names(input$getColumn(param$grouping))
  param$sampleGroup = make.names(param$sampleGroup)
  param$refGroup = make.names(param$refGroup)
  # colnames(input$meta) = gsub(' \\[.*','',colnames(input$meta))
  # if (param$grouping %in% colnames(input$meta))
  #   condition <- input$meta[[param$grouping]]

  ### # if conditions are not specified, then we have to stop here
  if (is.null(condition))
    stop(" * No conditions were specified in ezMethodDEXSeqAnalysis")

  ###Count only relevant bam-files
  samplesUse = condition %in% c(param$sampleGroup,param$refGroup)
  input = input$subset(samplesUse)
  condition = factor(condition[samplesUse], levels=c(param$refGroup, param$sampleGroup))


  sampleTable <- data.frame(
    row.names = rownames(input$meta),
    condition = condition
  )

  ### # get count files based on the name of the bamfiles
  sCountfileExt <- 'count'
  if (ezIsSpecified(param$countfile_ext))
    sCountfileExt <- param$countfile_ext
  countFiles <- gsub("bam$", replacement = sCountfileExt, basename(input$getColumn("BAM")))
  ### # if count files do not exist, generate them
  if(!all(file.exists(countFiles)))
    DEXSeqCounting(input = input, output = output, param = param)

  ### # check the reference
  if(ezIsSpecified(param$gff_file)){
    sRefFeatGff <- param$gff_file
  } else {
    sRefFeatGff <- gsub("gtf$", "gff", basename(param[['ezRef']]@refFeatureFile))
  }
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

  countData = counts(dxd)[,1:length(countFiles)]
  countDataPerGene = averageRows(countData, by=gsub(':.*','',rownames(countData)), func = sum)
  # rownames(countData) = gsub(':.*','',rownames(countData))
  # countDataPerGene = matrix(0,length(unique(rownames(countData))),ncol(countData))
  # colnames(countDataPerGene) = colnames(countData)
  # rownames(countDataPerGene) = unique(rownames(countData))
  # 
  # for (j in 1:ncol(countDataPerGene)){
  #    countDataPerGene[,j] = tapply(countData[,j],INDEX = rownames(countData),sum)
  # }

  presentGenes = rownames(countDataPerGene)[rowSums(countDataPerGene > param[['minGeneExprCount']]) > 1]
  #presentGenes = rownames(countDataPerGene)[which(rowMax(countDataPerGene)>param[['minGeneExprCount']] & apply(countDataPerGene,1,aboveMinExprSamples,minExpr=param[['minGeneExprCount']])>1)]
  useRow = gsub(':.*','',rownames(countData)) %in% presentGenes
  # filteredCountData = countData[useRow , ] + param$countOffset
  # transcripts = rowData(dxd)$transcripts[useRow]
  # featureRanges = rowRanges(dxd)[useRow]
  # 
  # dxd <- DEXSeq::DEXSeqDataSet(
  #   filteredCountData,
  #   sampleData    = sampleTable,
  #   design        = design,
  #   featureID     = gsub('.*:','',rownames(filteredCountData)),
  #   groupID       = gsub(':.*','',rownames(filteredCountData)),
  #   transcripts   = transcripts,
  #   featureRanges = featureRanges
  # )


  ### # estimate size factors and dispersion
  dxd <- DEXSeq::estimateSizeFactors( dxd )
  dxd <- DEXSeq::estimateDispersions( dxd, BPPARAM = BPPARAM, quiet = T )

  ### # testing for differential usage
  dxd  <- DEXSeq::testForDEU( dxd, BPPARAM = BPPARAM )

  ### # fold changes
  dxd <- DEXSeq::estimateExonFoldChanges( dxd, BPPARAM = BPPARAM, denominator = param$refGroup, independentFiltering = TRUE)

  ### # generate a report
  writeDEXSeqReport(dataset = input$meta, dexResult = list(param = param, dxd=dxd, useRow=useRow), 
                    sResultDir = basename(output$meta[['Report [File]']]))
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
                  appDefaults <<- rbind(disp_plot      = ezFrame(Type="character", DefaultValue="dispersion_estimate_plot", Description="which test method in DEXSeq to use: deseq2"),
                                        ma_plot        = ezFrame(Type="character", DefaultValue="ma_plot",    Description="no need to compute moderated ratios; deseq2 does this already"),
                                        countfile_ext  = ezFrame(Type="character", DefaultValue="count",      Description="extension of count files"),
                                        countfile_path = ezFrame(Type="character", DefaultValue=".",          Description="path where count files should be stored"),
                                        gff_file       = ezFrame(Type="character", DefaultValue="genes.gff",  Description="name of the gff annotation file"),
                                        paired         = ezFrame(Type="character",   DefaultValue='false',          Description="different counting for paired end data"),
                                        strandMode     = ezFrame(Type="character",   DefaultValue='both',          Description="read orientiation"),
                                        fdr            = ezFrame(Type="numeric",   DefaultValue=0.05,          Description="false discovery rate below which genes are reported"),
                                        minGeneExprCount   = ezFrame(Type="numeric",   DefaultValue=20,          Description="minimal Mean GeneCount for candidate selection"),
                                        minExonExprCount   = ezFrame(Type="numeric",   DefaultValue=10,          Description="minimal Mean ExonCount for candidate selection"),
                                        minExonLog2Ratio = ezFrame(Type="numeric",   DefaultValue=0.3,          Description="minimal log2Ratio for diff. exon for candidate selection"),
                                        bgExpression = ezFrame(Type="numeric",   DefaultValue=10,          Description="add pseudo counts to shrink logRatios"),
                                        dexseq_report_path = ezFrame(Type="character", DefaultValue="DEXSeqReport",  Description="path DEXSeqHTML report is written to"),
                                        dexseq_report_file = ezFrame(Type="character", DefaultValue="testForDEU.html",  Description="file name for DEXSeqHTML report")   )
                }
              ))


#' Addition experimental conditions to input files
#'
# addDEXSeqCondition = function(psInput, pvCondition){
#   x = ezRead.table(psInput)
#   x$Condition = pvCondition
#   write.table(x, file = psInput, quote = FALSE, sep = "\t")
# }


#' @title Writing a report for a DEXSeq analysis
#'
writeDEXSeqReport <- function(dataset, dexResult, htmlFile="00index.html", sResultDir = "html", BPPARAM=bpparam()) {
  ### # retrieve parameters
  param <- dexResult$param
  ### # extract name appearing in the report
  name <- param$name
  ### # extract DEXSeqResults object
  dxd <- dexResult$dxd
  dxr <- DEXSeq::DEXSeqResults(dxd,independentFiltering = TRUE)

  ### # put the results into a different subdirectory
  sCurWd <- getwd()
  setwdNew(sResultDir)

  ###Save dexResult as RData-Object for Shiny
  resultObj = list(dxd = dxd, param=param)
  resultObjFile = paste0("result--", param$comparison, "--", ezRandomString(length=12), "--EzDEXSeqResult.rds")
  saveRDS(resultObj,file=resultObjFile)

  ### # write tsv file from results
  sResultFile <- paste0("result--", param$comparison, "--", "DEXSeqResult.txt")
  dxr$transcripts = sapply(dxr$transcripts, function(trList){paste(trList, collapse="; ")})
  write.table(dxr, file = sResultFile, quote = FALSE, sep = "\t")

  ### # if parameters about biomart are specified, use them
  if (ezIsSpecified(param$bio_mart) &
      ezIsSpecified(param$mart_dataset) &
      ezIsSpecified(param$mart_filter) &
      ezIsSpecified(param$mart_attributes) ) {
    ensembl_mart <- biomaRt::useMart(biomart = param$bio_mart, dataset=param$mart_dataset)
    DEXSeq::DEXSeqHTML(dxr, path = param$dexseq_report_path, file = param$dexseq_report_file, FDR = param$fdr,
                       mart = ensembl_mart, filter = param$mart_filter, attributes = param$mart_attributes,
                       BPPARAM = BPPARAM)
  } else {
    ### # write that generic report for a given FDR, using 0.1 as the default
    DEXSeq::DEXSeqHTML(dxr, path = param$dexseq_report_path, file = param$dexseq_report_file, FDR = param$fdr,
                       BPPARAM = BPPARAM)
  }
  geneTableInfo = getGeneTable(pdxr = dxr, param = param)
  candidateReportFile = geneTableInfo[['candidateReportFile']]

  ### # put a title to the report using name in output
  titles <- list()
  titles[["DEXSeq Analysis"]]  <- paste("Analysis:", name)
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
  settings["Number of result features: "] = nrow(dxr)
  settings["Number of significant exons: "] = length(which(dxr$padj < param$fdr))
  settings["Number of candidate genes: "] = geneTableInfo[['nCandidates']]
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
  titles[["DEXSeq differential exon usage results"]] = "DEXSeq differential exon usage results"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addParagraph(doc, pot("Interactive Candidate Table", hyperlink = candidateReportFile))
  addParagraph(doc, pot("Txt-Result-File Exon-Level", hyperlink = sResultFile))

  ### # Put simply a link to the already existing report
  titles[["Report generated by DEXSeq"]] = "Report generated by DEXSeq"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  addParagraph(doc, pot("Test results for differential exon usage", hyperlink = "DEXSeqReport/testForDEU.html"))



  ### # closing the report leads to writing it to htmlFile
  titles[["Misc"]] = "Misc"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
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
  log2column = grep('log2fold',colnames(pdxr), value=TRUE)
  results <- data.frame(pdxr[, c("groupID",
                                 "featureID",
                                 "exonBaseMean",
                                 "dispersion",
                                 log2column,
                                 "pvalue",
                                 "padj")],
                        stringsAsFactors = FALSE)
  results <- cbind(results, genomicData)
  results[, c("dispersion",log2column, "pvalue", "padj")] <- round(results[, c("dispersion",log2column, "pvalue", "padj")], 3)
  dexseqR <- elementMetadata(pdxr)$type == "DEXSeq results"
  if (sum(dexseqR, na.rm = TRUE) > 0) {
    results <- cbind(results, round(as.data.frame(pdxr[, which(dexseqR)]), 3))
  }
  rownames(results) <- NULL

  ### # from the results, generate the genetable which seems to be the basis for the
  ### #  table on the results page
  #gns <- as.character(unique(results$groupID[which(results$padj < param$fdr)]))
  exns <- results[which(results$padj < param$fdr &
                          abs(results[[log2column]]) > param$minExonLog2Ratio
                        & results[['exonBaseMean']] > param$minExonExprCount), c('groupID','featureID',log2column,'padj')]
  gns  <- names(perGeneQValue(pdxr)) ##[perGeneQValue(pdxr) < param$fdr] -- report all genes
  ##gns  <- intersect(gns,unique(exns$groupID)) --- completely unclear why there should be a non-intersecting gene id
  #names(perGeneQValue(pdxr)< param$fdr)
  #results <- results[as.character(results$groupID) %in% gns,]
  dfGnsAnnot <- ezFeatureAnnotation(param, ids=NULL, "gene")
  genetable = ezFrame(gene_id=unique(results$groupID))
  rownames(genetable) = genetable$gene_id
  genetable$gene_name <- sapply(genetable$gene_id,
                                function(x) {
                                  geneIds <- unlist(strsplit(x, split = "+", fixed = TRUE))
                                  geneNames = ezCollapse(dfGnsAnnot[geneIds, "gene_name"], uniqueOnly=TRUE, empty.rm=TRUE, sep="+")
                                  return(geneNames)
                                },
                                USE.NAMES=TRUE)[genetable$gene_id]
  genetable$gene_description <- sapply(genetable$gene_id,
                                       function(x) {
                                         geneIds <- unlist(strsplit(x, split = "+", fixed = TRUE))
                                         geneNames = ezCollapse(dfGnsAnnot[geneIds, "description"], uniqueOnly=TRUE, empty.rm=TRUE, sep="+")
                                         return(geneNames)
                                       },
                                       USE.NAMES = TRUE)[genetable$gene_id]
  genetable$seqnames=tapply(as.character(results$seqnames), results$groupID, unique)[genetable$gene_id]
  genetable$end = tapply(results$end, results$groupID, max)[genetable$gene_id]
  genetable$exon_changes = tapply(results$padj < param$fdr & 
                                    results$exonBaseMean>param$minExonExprCount & 
                                    abs(results[[log2column]]) > param$minExonLog2Ratio,
                                  results$groupID, sum, na.rm = TRUE)[genetable$gene_id]
  genetable$meanRawCount = round(tapply(results$exonBaseMean, results$groupID, sum)[genetable$gene_id],3)
  genetable = genetable[genetable$exon_changes > 0, ]
  ####TODO: get ExonLog2FC only for significant exons with realistic log2FC calculation (pseudo count +10), only exon longer than >50nt
  ####Filter based on ExonLog2FC  
  genetable$max_ExonLog2FC = round(tapply(results[[log2column]], results$groupID, 
                                          function(x){
                                            idx = which.max(abs(x))
                                            if (length(idx) == 0){
                                              return(0)
                                            } else {
                                              return(x[idx])
                                            }
                                          })[rownames(genetable)], 3)
  ### # add extracted columns and return genetable
  genetable$fdr = round(perGeneQValue(pdxr),5)[genetable$gene_id]
  genetable$gene_id <- getGeneIdExprLinks(pvGeneIds = genetable$gene_id, psdexseq_report_path = param$dexseq_report_path)
  genetable = genetable[order(genetable$fdr),]
  ezWrite.table(genetable, file=paste0("result--", param$comparison, "--DEXSeqGeneResult.txt"), head="gene_id")
  require(DT)
  geneTableSel = genetable[!is.na(genetable$fdr) & genetable$fdr < param$fdr & abs(genetable$max_ExonLog2FC) > param$minExonLog2Ratio, ]
  x = datatable(geneTableSel, escape = F,rownames = FALSE, filter = 'bottom',extensions = c('ColReorder','Buttons'),
                caption = paste('Candidates DEXSeq:',param$comparison, sep=''),
                options = list(
                  initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': '#0000A0', 'color': '#fff'});",
                    "}"),
                  dom = c('Bfrtip'),
                  buttons = c('colvis','copy', 'csv', 'excel', 'pdf', 'print')))
  candidateReportFile = paste0('candidates_',param$comparison,'.html')
  saveWidget(x,candidateReportFile)


  return(list(candidateReportFile=candidateReportFile,nCandidates=nrow(genetable)))
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
  ramPerJob = round((param[['ram']]*1000)/param[['cores']] * 0.25) ## keep a reserve of 25% RAM ... for samtools being greedhy and the Dexseq python script
  vCountFiles <- ezMclapply(bamFiles, runCountSingleBam, sGffFile, sCountfileExt, param$strandMode, param$paired, ramPerJob, param$bgExpression, mc.cores = param[['cores']])

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
#' which is taken from the result of function
#' \code{lGetPyScriptPaths}
#'
#' @param psGtfFile   name of the GTF annotation file
#' @param psGffFile   name of the GFF file to be generated
convertGtfToGff <- function(psGtfFile, psGffFile) {
  cat(" * Converting GTF to GFF ...\n")
  ezSystem(paste('cp', psGtfFile, '.'))
  ezSystem(paste('grep protein_coding', basename(psGtfFile), '>protein_coding_genes.gtf'))
  sPyConvCmd <- paste(
    "python",
    file.path(system.file(package = "DEXSeq", "python_scripts"), "dexseq_prepare_annotation.py"),
    'protein_coding_genes.gtf',
    psGffFile)
  ezSystem(sPyConvCmd)
  cat("  ==> created: ", psGffFile, "\n")
  invisible(TRUE)
}

#' Run counts for a single BAM file
#'
runCountSingleBam <- function(psBamFile, psGffFile, psCountfileExt, strandMode, Paired, ramPerJob, bgExpression){
  if(Paired){
    # samCmd = paste("samtools", "sort -n" , psBamFile, "-m", paste0(ramPerJob,"M"), "-O SAM")
    stopifnot(psBamFile != basename(psBamFile))
    ezSystem(paste("samtools", "sort -n" ,psBamFile, "-m", paste0(ramPerJob,"M"), "-o", basename(psBamFile)))
    psBamFile = basename(psBamFile)
  }
  samCmd = paste("samtools", "view -h", psBamFile)

  ### # run counting on sam file
  sCountBaseFn <- gsub("bam$", psCountfileExt, basename(psBamFile))
  cmd <- paste("python", file.path(system.file(package = "DEXSeq", "python_scripts"), "dexseq_count.py"))
  if(Paired){
    cmd = paste(cmd, '--paired yes')
  }

  if(strandMode=='antisense'){
    cmd = paste(cmd,'--stranded reverse')
  } else if(strandMode=='both'){
    cmd = paste(cmd,'--stranded no')
  }

  sPyCountCmd <- paste(samCmd, "|", cmd, psGffFile, "-", sCountBaseFn, "2>", paste0(sCountBaseFn, ".err"))
  ezSystem(sPyCountCmd)
  myCounts = ezRead.table(sCountBaseFn, row.names = NULL, header = F)
  myCounts[,2] = myCounts[,2] + bgExpression
  ezWrite.table(myCounts, sCountBaseFn, col.names = F, row.names = F)
  return(file.path(getwd(), sCountBaseFn))
}


#' Write names of countfiles back into the input file
#'
writeCountFilesToMeta <- function(pvCountFiles, input) {
  ### # add column with counts to the meta information
  input$meta$Count <- pvCountFiles
  ### # write the extended meta information back to the file
  write.table(input$meta, file = input$file, quote = FALSE, sep = "\t")
}


aboveMinExprSamples = function(x,minExpr=5){
  nSamples = sum(x>minExpr)
  return(nSamples)
}
