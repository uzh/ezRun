###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodInlineBarcodeDmx <- function(input = NA, output = NA, param = NA) {
    require(ShortRead)
    setwdNew(param[['name']])
    sampleName <- input$getNames()
    inputFile  <- input$getFullPaths("Read1")
    indexInfo <- ezRead.table(input$meta[['indexFile']], row.names = NULL, sep = ',')
    index <- indexInfo$Barcode
    index <- setNames(index, indexInfo$Name)
    dataChunks <- 5*10^6
    
    strm <- FastqStreamer(inputFile, n = dataChunks)
    processedReads <- 0
    resultFiles <- paste(names(index), index, 'R1.fastq.gz', sep = '_')
    k = 1
    repeat {
        currentReads <- yield(strm)
        reads <- sread(currentReads)
        maxReadLength <- max(width(reads))
        qual <- quality(currentReads)
        ids <- id(currentReads)
        reads <- DNAStringSet(substr(reads, param$barcodePos, maxReadLength))
        qual <- BStringSet(substr(quality(qual),param$barcodePos, maxReadLength))
        
        for (i in 1:length(index)){
            file <- paste(names(index)[i], index[i], 'R1.fastq.gz', sep = '_')
            
            leftPattern = index[i]
            vp <- vmatchPattern(leftPattern, reads, max.mismatch = 0)
            leftEnd <- vapply(endIndex(vp),
                              .dummyFunction, c('endIndex' = 0))
            pos <- which(leftEnd == nchar(leftPattern))
            filteredReads <- reads[pos]
            filteredQual <- qual[pos]
            filteredID <- ids[pos]
            vp <- vmatchPattern(param$rightPattern, filteredReads, max.mismatch = param$maxMismatch)
            rightStart <- vapply(startIndex(vp),
                                 .dummyFunction, c('startIndex' = 0))
            rightStart[is.na(rightStart)] <- maxReadLength
            toSmall <- which(rightStart < (nchar(leftPattern)))
            rightStart[toSmall] <- maxReadLength
            filteredReads <- DNAStringSet(substr(filteredReads, nchar(leftPattern) + param$leftLinkerSize, rightStart))
            filteredQual <- BStringSet(substr(filteredQual, nchar(leftPattern) + param$leftLinkerSize, rightStart))
            filteredReads <- filteredReads[width(filteredReads) >= param$minReadLength & width(filteredReads) <= param$maxReadLength]
            filteredID <- filteredID[width(filteredQual) >= param$minReadLength & width(filteredQual) <= param$maxReadLength]
            filteredQual <- filteredQual[width(filteredQual) >= param$minReadLength & width(filteredQual) <= param$maxReadLength]
            
            if(k == 1){
                writeFastq(ShortReadQ(sread = filteredReads, quality = filteredQual, id = filteredID), file, mode="w", full=FALSE, compress=TRUE)
            } else {
                writeFastq(ShortReadQ(sread = filteredReads, quality = filteredQual, id = filteredID), file, mode="a", full=FALSE, compress=TRUE)
            }
        }
        processedReads = processedReads + dataChunks
        ezLog(paste0(processedReads/10^6, 'M reads processed'))
        if (length(currentReads) == 0)
            break
        
        k <- k + 1
        gc()
    }
    gc()
    
    ###Generate Dataset:
    dataset <- data.frame(Name = names(index), Read1 = resultFiles, Species = '', ReadCount = 0)
    resultFiles <- file.path(param[['resultDir']], param[['name']], resultFiles)
    
    for (i in 1:nrow(dataset)){
        ezLog(paste('Count Reads for ',dataset$Name[i], '...'))
        myCount <- system(paste('zcat', dataset$Read1[i], '|wc -l'), intern = TRUE)
        ezLog(paste('Count Reads for ',dataset$Name[i], ' done'))
        dataset$ReadCount[i] <- as.numeric(myCount)/4
    }
    dataset$Read1 <- resultFiles
    
    colnames(dataset)[2] <- 'Read1 [File]' 
    colnames(dataset)[4] <- 'Read Count'
    
    ezWrite.table(dataset, paste0(sampleName, '_dataset.tsv'), row.names = FALSE)
    return('sucess')
    }

.dummyFunction <- function(x) {
    if (is.null(x))
        NA
    else
        x[1]
}


##' @template app-template
##' @templateVar method ezMethodInlineBarcodeDmx(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppInlineBarcodeDmx <-
    setRefClass("EzAppInlineBarcodeDmx",
                contains = "EzApp",
                methods = list(
                    initialize = function() {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodInlineBarcodeDmx
                        name <<- "EzAppInlineBarcodeDmx"
                        appDefaults <<- rbind(barcodePos = ezFrame(Type = "numeric", DefaultValue = 8, Description = "StartPosition of Inline Barcode in Read1"),
                                              leftLinkerSize = ezFrame(Type = "numeric", DefaultValue = 5, Description = "Size of 5-Prime Linker next to Inline Barcode"),
                                              rightPattern = ezFrame(Type = "character", DefaultValue = "GTGTCAGTCACTTCCAG", Description = "3-Prime linker sequence for trimming"),
                                              maxMismatch = ezFrame(Type = "numeric", DefaultValue = 1, Description = "Accepted mismatches in 3-Prime linker sequence"),
                                              minReadLength = ezFrame(Type = "numeric", DefaultValue = 20, Description = "minimal read length after dmx, filtering and trimming"),
                                              maxReadLength = ezFrame(Type = "numeric", DefaultValue = 50, Description = "maximal read length after dmx, filtering and trimming"))
                    }
                )
    )