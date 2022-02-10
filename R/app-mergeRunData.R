###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch
ezMethodMergeRunData <- function(input=NA, output=NA, param=NA){
    setwdNew(param[['Name']])
    project = dirname(param[['resultDir']])
    matchCol = param[['matchingColumn']]
    
    dataset1 = input$meta
    dataset1$Name = rownames(dataset1)
    rownames(dataset1) <- NULL
    dataset1 <- dataset1[,c('Name', colnames(dataset1)[1:(length(colnames(dataset1))-1)])]
    
    inputDir2 = file.path(param[['dataRoot']], project, param[['DataSetName2']])
    datasetFile2 = list.files(inputDir2, pattern='^dataset.tsv$', full.names = TRUE)
    input2 <- EzDataset(datasetFile2)
    dataset2 <-input2$meta
    dataset2$Name = rownames(dataset2)
    rownames(dataset2) <- NULL
    dataset2 <- dataset2[,c('Name', colnames(dataset2)[1:(length(colnames(dataset2))-1)])]
    
    #Find Factor columns and replace NA with '':
    #factorColPositions <- grep('\\[Factor\\]$', colnames(dataset2))
    #if(length(factorColPositions) > 0){
    #    for (i in 1:length(factorColPositions)){
    #        if(all(is.na(dataset2[,factorColPositions[i]]))){
    #            dataset2[,factorColPositions[i]] = ''
    #        }
    #    }
    #}
    #Remove almost empty dataFiles:
    dataset1 = dataset1[dataset1[['Read Count']] >= param$minReadCount,]
    dataset2 = dataset2[dataset2[['Read Count']] >= param$minReadCount,]
    
    commonCols = intersect(colnames(dataset1), colnames(dataset2))
    
    dataset1 = dataset1[,commonCols]
    dataset2 = dataset2[,commonCols]
    dataset1[[matchCol]] = gsub('/', '_', dataset1[[matchCol]])
    dataset2[[matchCol]] = gsub('/', '_', dataset2[[matchCol]])

    #SampleData to Pool:
    intersectNames = intersect(dataset1[[matchCol]], dataset2[[matchCol]])
    uniqSet1 = setdiff(dataset1[[matchCol]], dataset2[[matchCol]])
    uniqSet2 = setdiff(dataset2[[matchCol]], dataset1[[matchCol]])
    outputRunName = paste(gsub('-','',Sys.Date()),'.X-',sep='')
    
    cat("#Files to merge:", length(intersectNames),"\n")
    cat("#Files only in ", param[['FirstDataSet']], ':', length(uniqSet1), "\n")
    cat("#Files only in ", basename(inputDir2), ':', length(uniqSet2), "\n")
    
    for (i in 1:length(intersectNames)){
        file1 = file.path(param[['dataRoot']], dataset1[dataset1[[matchCol]] == intersectNames[i],'Read1 [File]'])
        file2 = file.path(param[['dataRoot']], dataset2[dataset2[[matchCol]] == intersectNames[i],'Read1 [File]'])
        
        mergedFile = paste0(outputRunName, intersectNames[i], '_R1.fastq.gz')
        cmd = paste('cat', file1, file2, '>', mergedFile)
        ezSystem(cmd)
        
        if(param[['paired']]){
            file1 = file.path(param[['dataRoot']], dataset1[dataset1[[matchCol]] == intersectNames[i], 'Read2 [File]'])
            file2 = file.path(param[['dataRoot']], dataset2[dataset2[[matchCol]] == intersectNames[i], 'Read2 [File]'])
            
            mergedFile = paste0(outputRunName, intersectNames[i], '_R2.fastq.gz')
            cmd = paste('cat', file1, file2, '>', mergedFile)
            ezSystem(cmd)
        }
    }
    
    ###Create new dataset
    dataset = rbind(dataset1, dataset2)
    dataset = dataset[order(dataset[[matchCol]]), ]
    dataset = dataset[dataset[[matchCol]] %in% intersectNames, ]
    colsToRemove <- c('Read1 [File]', 'Read Count', 'SampleConc [Characteristic]', 'InputAmount [Characteristic]', 'Sample Id [B-Fabric]', 'RIN [Characteristic]', 'PlatePosition [Characteristic]',
      'LibConc_100_800bp [Characteristic]', 'LibraryPrepKit', 'PlateName [Characteristic]')
    
    if('Order Id [B-Fabric]' %in% colnames(dataset)){
        if(length(unique(dataset[['Order Id [B-Fabric]']])) > 1){
            colsToRemove <- c(colsToRemove, 'Order Id [B-Fabric]', 'Tube [Characteristic]', 'EnrichmentMethod')
        }
    }
    
    if('Species' %in% colnames(dataset)){
        if(length(unique(dataset[['Species']])) > 1){
            colsToRemove <- c(colsToRemove, 'Species')
        }
    }
    
    uniqueDataset = unique(dataset[ ,-which(colnames(dataset) %in% colsToRemove)])
    if(param[['paired']]){
        uniqueDataset = unique(uniqueDataset[, -which(colnames(uniqueDataset) %in% c('Read2 [File]'))])
    }
    uniqueDataset[['Read Count']] = tapply(dataset[['Read Count']], dataset[[matchCol]], sum)
    
    ###Add Read1 Column
    uniqueDataset[['Read1 [File]']] = file.path(param[['resultDir']], param[['Name']], paste0(outputRunName, uniqueDataset[[matchCol]], '_R1.fastq.gz'))
    if(param[['paired']]){
        uniqueDataset[['Read2 [File]']] = file.path(param[['resultDir']], param[['Name']], paste0(outputRunName, uniqueDataset[[matchCol]], '_R2.fastq.gz'))
    }
    
    datasetKeep = rbind(dataset1[which(dataset1[[matchCol]] %in% uniqSet1),], dataset2[which(dataset2[[matchCol]] %in% uniqSet2),])
    if(nrow(datasetKeep) > 0){
        ##add unique rows to dataset and transfer run specific samples to gstore
        datasetKeep = datasetKeep[,intersect(colnames(uniqueDataset),colnames(datasetKeep))]
        uniqueDataset = rbind(uniqueDataset, datasetKeep)
        cmd <- paste('rsync', file.path(param[['dataRoot']], datasetKeep[['Read1 [File]']]), '.')
        for (k in 1:length(cmd)){
                ezSystem(cmd[k])
        }
        if(param[['paired']]){
            cmd <- paste('rsync', file.path(param[['dataRoot']], datasetKeep[['Read2 [File]']]), '.')
            for (k in 1:length(cmd)){
                ezSystem(cmd[k])
            }
        }
    }
    ezWrite.table(uniqueDataset, 'dataset.tsv', row.names = FALSE)
    ##compute md5 sums
    ezSystem('md5sum *.gz > md5.txt')
    return('success')
}

##' @template app-template
##' @templateVar method ezMethodMergeRunData(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppMergeRunData <-
    setRefClass("EzAppMergeRunData",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodMergeRunData
                        name <<- "EzAppMergeRunData"
                        appDefaults <<- rbind(
                            minReadCount = ezFrame(
                                Type = "integer",
                                DefaultValue = 10000,
                                Description = "minimal Read Count to consider a sample for merging in a run"
                            )
                        )
                    }
                )
    )
