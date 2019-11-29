
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Create and save SGSea database 
##' @description Create and save SGSea database 
##' @param  a gtf file
##' @return a makeTxDbFromGFF database object
 
makeSGSeqtxdbObject <- function(refBuil,refFeat,txdbDatabase) {
txdbGtfOrigin <- file.path("/srv/GT/reference",refBuil,"Genes",refFeat)
organism <- unlist(strsplit(refBuil,"/"))[1]
organismForDB <- gsub("_"," ",organism)
dataSource <- unlist(strsplit(refBuil,"/"))[2]
build <- unlist(strsplit(refBuil,"/"))[3]
genomeIndexFile <- file.path("/srv/GT/reference",organism,dataSource,build,
                             "Sequence","WholeGenomeFasta","genome.fa.fai")
genomeIndex <- read.delim(genomeIndexFile, stringsAsFactors = F, header = F)
chrInfo <- data.frame(chrom=genomeIndex$V1,length=genomeIndex$V2)
txDbObject <- makeTxDbFromGFF(txdbGtfOrigin,format = "gtf", organism = organismForDB,
                              dataSource = dataSource, chrominfo = chrInfo)
saveDb(txDbObject, file=txdbDatabase)
return(txDbObject)
}


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Generate flattened file for SGSeq-based differential expression 
##' @description Generate flattened file for SGSeq-based differential expression 
##' @param  a list of count files generated with SGseq functions analyzeFeatures, analyzeVariants and counts
##' @return a text file

generateFlattenedFileForSGeq <- function(x){
  listOfCountFiles <- list.files(path = './', pattern = "inputForDexSeq.")
  allSampleNames <- ldply(strsplit(listOfCountFiles,"\\."))[,2]
  RdataToLoad <- paste("sgSplice",allSampleNames[1],"RData",sep = ".")
  load(RdataToLoad)
  mainTable <- mcols(sgfc_ucscV)
  ### functions to create the appropriate row format from sgSeq variant file
  ## gene row
  buildGeneRowsForFlattenedFile <- function(test){
    chrName <- strsplit(test$from,":")[[1]][2]
    script <- "flattenedFileForSgSeq"
    featType <- "aggregate_gene"
    start <- min(sapply(strsplit(test$from,":"), function(x) x[3]))
    end <-  max(sapply(strsplit(test$to,":"), function(x) x[3]))
    effStart <- min(start,end)
    effEnd <- max(start,end)
    dot1 <- "."
    strand <- strsplit(test$from,":")[[1]][4]
    dot2 <- "."
    descr <- paste0("gene_id ", "\"",unique(test$eventID),"\"")
    tempEventDF <- data.frame(chrName,script,featType,effStart,stringsAsFactors = FALSE)
    return(tempEventDF)
  }
  ## exon rows
  buildExonRowsForFlattenedFile <- function(test){
    chrName <- sapply(strsplit(test$from,":"), function(x) x[2])
    script <- rep("flattenedFileForSgSeq", nrow(test))
    featType <- rep("exonic_part", nrow(test))
    effStart <- sapply(strsplit(test$from,":"), function(x) x[3])
    effEnd <-  sapply(strsplit(test$to,":"), function(x) x[3])
    dot1 <- rep(".", nrow(test))
    strand <- sapply(strsplit(test$from,":"),function(x) x[4])
    dot2 <- rep(".", nrow(test))
    geneToAnnotate <- unique(unlist(test$geneName))
    if (length(geneToAnnotate) > 0){
      geneNameToMap <- AnnotationDbi::mget(x=geneToAnnotate,envir=org.Mm.egSYMBOL,ifnotfound = NA)[[1]]
    }else{
      geneNameToMap <- "noGeneName"
    }
    descr <- paste0("transcripts ", "\"", geneNameToMap, "\"; ",
                    "exonic_part_number ","\"", test$variantID ,"\"; ",
                    "gene_id ", "\"",unique(test$eventID),"\"")
    tempExonDF <- data.frame(cbind(chrName,script,featType,effStart,
                                   effEnd,dot1,strand,dot2,descr),
                             stringsAsFactors = FALSE)
    return(tempExonDF)
  }
  ## gene-wide mergin
  listOfEventIds <- levels(as.factor(mainTable$eventID))
  listOfTempDF <- lapply(listOfEventIds, function(x) {
    eventIdVar <- mainTable[mainTable$eventID == x,]
    geneRow <- buildGeneRowsForFlattenedFile(eventIdVar)
    exonRows <- buildExonRowsForFlattenedFile(eventIdVar)
    geneRow$effStart <- as.numeric(geneRow$effStart)
    geneRow$effEnd <- as.numeric(geneRow$effEnd)
    exonRows$effStart <- as.numeric(exonRows$effStart)
    exonRows$effEnd <- as.numeric(exonRows$effEnd)
    rbind(geneRow,exonRows)
  })
  ## global merging
  finalFlattened <- do.call("rbind",listOfTempDF)
  indexToSwap <- which(finalFlattened$effStart > finalFlattened$effEnd)
  finalFlattened$temp1 <- finalFlattened$effStart
  finalFlattened$temp2 <- finalFlattened$effEnd
  finalFlattened[indexToSwap,]$temp1 <- finalFlattened[indexToSwap,]$effEnd
  finalFlattened[indexToSwap,]$temp2 <- finalFlattened[indexToSwap,]$effStart
  finalFlattened$effStart <- finalFlattened$temp1
  finalFlattened$effEnd <- finalFlattened$temp2
  finalFlattened <- subset(finalFlattened, select=chrName:descr)
  write.table(finalFlattened,'dexSeqFlattenedFileForSgSeq.txt', 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}