###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

#' RunMethod for reference class EzAppDEXSeqAnalysis
#' 
ezMethodDEXSeqAnalysis <- function(input=NA, output=NA, param=NA){
  ### # if count files are not available, generate them
  if(!is.element("Count", colnames(input$meta))) {
    EzAppDEXSeqCounting$new()$run(input = input, output = output, param = param)
    input <- EzDataset$new(file = input$file)    
  }
    
  ### # get count files
  countFiles <- input$getFullPaths(param, "Count")
  ### # check whether conditions are specified
  if ("Condition" %in% colnames(input$meta))
    condition <- input$meta$Condition
  if (!is.null(param$condition))
    condition <- param$condition
  ### # if conditions are not specified, then we have to stop here
  if (is.null(condition))
    stop(" * No conditions were specified in ezMethodDEXSeqAnalysis")
  ### # sample table from rownames and conditions
  sampleTable <- data.frame(
    row.names = rownames(input$meta),
    condition = condition
  )
  ### # check the reference
  sRefFeatGff <- gsub("gtf$", "gff", param[['ezRef']]@refFeatureFile)
  stopifnot(file.exists(sRefFeatGff))
  
  ### # check whether special design was specified, o/w use default
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
  
  ### # extract the result from the differential analysis
  dxr <- DEXSeq::DEXSeq(dxd)
  
  ### # generate a report
  DEXSeq::DEXSeqHTML(dxr)
  
  return("Success")  
}


##' @template app-template
##' @templateVar method ezMethodDEXSeqCounting
##' @templateVar htmlArg )
##' @description Use this reference class to run 
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
