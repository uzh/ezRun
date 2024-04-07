###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch
ezMethodMageckCountQC = function(input=NA, output=NA, param=NA,
                           htmlFile="00index.html"){
    setwdNew(basename(output$getColumn("Report")))
    param <- ezParam(list(refBuild='Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30'))
    param$expressionName <- 'sample1'
    rawData <- loadMageckCountDataset(input, param)
    metadata(rawData)$output <- output
    makeRmdReport(rawData=rawData, rmdFile="CrisprCountQC.Rmd", 
                  reportTitle="ExploreCounts", selfContained = TRUE)
    return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMageckCountQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppExploreMageckCounts <-
    setRefClass("EzAppExploreMageckCounts",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodMageckCountQC
                        name <<- "EzAppExploreMageckCounts"
                        appDefaults <<- rbind(runGO=ezFrame(Type="logical", DefaultValue=FALSE, Description="whether to run the GO analysis"),
                                              nSampleClusters=ezFrame(Type="numeric", DefaultValue=6, Description="Number of SampleClusters, max value 6"))
                    }
                )
    )

loadMageckCountDataset <- function(input, param){
    require(tools)
    require(SummarizedExperiment)
    files <- input$getFullPaths("Count")
    
    dataFeatureLevel <- 'miRNA' #unique(input$getColumn("featureLevel"))
    stopifnot(length(dataFeatureLevel) == 1)
    
    x1 <- read_tsv(files[1], guess_max=1e6,  col_types = cols())
    ## col type messages are suppressed by col_types = cols()
    ## Alternatively: The message is triggered by readr:::show_cols_spec. 
    ##To suppress the message, put this at the top of your script: options(readr.num_columns = 0)
    
    if (ezIsSpecified(param$expressionName)){
        columnName <- param$expressionName
    } else {
        columnName <- intersect(param$knownExpressionNames, colnames(x1))[1]
    }
    if (!columnName %in% colnames(x1)){
        return(list(error=paste0("Specified column name not found in data!<br>columnName: '", columnName, "'\n",
                                 "<br>Available column names:<br>\n",
                                 paste0("'", colnames(x1), "'", collapse="<br>"),
                                 "<br>Set the option columnName to one of the names above!")))
    }
    identifier <- 1
    
    x <- mapply(function(x, y){
        message("loading file: ", x)
        tempTibble <- read_tsv(x, progress=FALSE, guess_max=1e6, col_types = cols())
        tempTibble %>%
            dplyr::select(identifier, columnName) %>%
            dplyr::rename("id":= 1, !! y := columnName)
    }, files, names(files), SIMPLIFY=FALSE)
    x <- Reduce(function(x,y){full_join(x, y, by="id")}, x)
    
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
    counts <- as.matrix(x[ ,-identifier])
    rownames(counts) <- x[[identifier]]
    counts[is.na(counts)] <- 0
    
    if(ezIsSpecified(param$ezRef@refBuild)){
        seqAnnoDFFeature <- ezFeatureAnnotation(param, rownames(counts),
                                                param$featureLevel)
    }else{
        seqAnnoDFFeature <- .makeSeqAnnoFromCounts(x, identifier)
    }
    stopifnot(identical(rownames(seqAnnoDFFeature), rownames(counts)))
    
    if (ezIsSpecified(param$correctBias) && param$correctBias){
        ## output will be floating point, but we don't round; input might already be floating point
        counts <- ezCorrectBias(counts, gc=seqAnnoDFFeature$gc,
                                width=seqAnnoDFFeature$featWidth)$correctedCounts
    }
    
    if (ezIsSpecified(param$runRUV) && as.logical(param$runRUV)){
        library(RUVSeq)
        differences <- makeGroups(param$grouping)
        ruvCorr = RUVs(counts, cIdx=rownames(counts), k=as.integer(param$kRUVFactors), scIdx=differences, epsilon=10)
        counts <- ruvCorr$normalizedCounts
    }
    
    
    seqAnno <- makeGRangesFromDataFrame(seqAnnoDFFeature, keep.extra.columns=TRUE)
    
    if (param$useSigThresh){
        sigThresh = param$sigThresh
    } else {
        sigThresh = 0
    }
    
    ## assays: counts, presentFlag
    ## rowData, colData
    ## meta: isLog, featureLevel, type, countName, param
    rawData <- SummarizedExperiment(
        assays=SimpleList(counts=counts, presentFlag=counts > sigThresh),
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
    return(rawData)
}

#library(ezRun)
#dataset <- ezRead.table('/srv/gstore/projects/p33200/MageckCount_Baggen_2021_2024-04-04--15-31-22/dataset.tsv')
#input <- EzDataset$new(meta=dataset, dataRoot = DEFAULT_DATA_ROOT)
#param <- ezParam(list(refBuild='Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30'))
#param$expressionName <- 'sample1'

#rawData <- loadMageckCountDataset(input, param)

#metadata(rawData)$output <- output
#setwdNew('CountQC')
#makeRmdReport(rawData=rawData, rmdFile="CountQC.Rmd", 
#              reportTitle="CountQC", selfContained = TRUE)