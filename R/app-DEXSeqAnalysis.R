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
  ### # if count files are not available, generate them
  if(!is.element("Count", colnames(input$meta)) || !all(file.exists(input$meta$Count))) {
    EzAppDEXSeqCounting$new()$run(input = input, output = output, param = param)
    input <- EzDataset$new(file = input$file)    
  }
    
  ### # get count files
  countFiles <- input$getFullPaths(param, "Count")
  ### # check whether conditions are specified
  if (param$grouping %in% colnames(input$meta))
    condition <- input$meta$Condition
  if (!is.null(param$condition))
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
  sRefFeatGff <- gsub("gtf$", "gff", param[['ezRef']]@refFeatureFile)
  stopifnot(file.exists(sRefFeatGff))
  
  ### # check whether special design was specified, o/w use minimal default design
  if (!is.null(param$design)) {
    design <- param$design
  } else {
    design <- ~ sample + exon + condition:exon
  }
    
  ### # create the initial DEXSeqDataSet object
  dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData    = sampleTable,
    design        = design,
    flattenedfile = sRefFeatGff )
  
  ### # plot disersion estimates, if indicated by the parameter
  if (!is.null(param$disp_plot)) {
    sDispPlotPdfFile <- ifelse(identical(tools::file_ext(param$disp_plot), "pdf"), param$disp_plot, paste(param$disp_plot, "pdf", sep = "."))
    pdf(file = sDispPlotPdfFile)
    dxd %>% 
      DEXSeq::estimateSizeFactors %>% 
      DEXSeq::estimateDispersions %>% 
      DEXSeq::plotDispEsts
    dev.off()
  }
  ### # extract the result from the differential analysis
  dxr <- DEXSeq::DEXSeq(dxd)
  
  ### # plotting MA, if indicated
  if (!is.null(param$ma_plot)){
    sMaPlotPdfFile <- ifelse(identical(tools::file_ext(param$ma_plot), "pdf"), param$ma_plot, paste(param$ma_plot, "pdf", sep = "."))
    pdf(file = sMaPlotPdfFile)
    DEXSeq::plotMA( dxr )
    dev.off()
  }
  
  ### # write tsv file from results
  if (!is.null(param$ResultFile)) {
    sResultFile <- param$ResultFile
  } else {
    sResultFile <- "DexSeqResult.tsv"
  }
  write.table(dxr, file = sResultFile, quote = FALSE, sep = "\t")

  ### # generate a predefined html-report
  if (!is.null(param$output_format) && identical(tolower(param$output_format), 'html'))
    DEXSeq::DEXSeqHTML(dxr)
  
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
                }
              ))


#' Addition experimental conditions to input files
#'
addDEXSeqCondition = function(psInput, pvCondition){
  ezObjInput <- EzDataset$new(file=psInput)
  ezObjInput$meta$Condition <- pvCondition
  write.table(ezObjInput$meta, file = ezObjInput$file, quote = FALSE, sep = "\t")
}

