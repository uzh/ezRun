ezMethodGetEnaData <- function(input=NA, output=NA, param=NA){
    require(XML)
    containerId <- sub('p', '', strsplit(output$getColumn('ENA Result'),'/')[[1]][1])
    cmd = paste0("curl -o fastqLinks.txt -X GET ","\'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",param[['projectID']],"&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,fastq_md5&format=tsv&download=true\'")
    ezSystem(cmd)
    fastqInfo = ezRead.table('fastqLinks.txt', row.names = NULL)
    fastqInfo[['Name']] = ''
    fastqInfo[['ReadCount']] = 0
    colnames(fastqInfo)[grep('scientific_name', colnames(fastqInfo))] = 'Species'
    
    
    ##Filter for included or excluded samples:
    if(param$excludedSamples != ''){
        excludedSamples <- as.vector(gsub(' ', '', limma::strsplit2(param$excludedSamples, split = ',')))
        toExclude <- which(fastqInfo[['sample_accession']] %in% excludedSamples)
        if(length(toExclude) > 0){
            cat('Samples excluded:')
            print(fastqInfo[toExclude,'sample_accession'])
            fastqInfo <- fastqInfo[-c(toExclude),]
        }
    }
    
    if(param$includedSamples != ''){
        includedSamples <- as.vector(gsub(' ', '', limma::strsplit2(param$includedSamples, split = ',')))
        toInclude <- which(fastqInfo[['sample_accession']] %in% includedSamples)
        if(length(toInclude) > 0){
            cat('Samples included:')
            print(fastqInfo[toInclude,'sample_accession'])
            fastqInfo <- fastqInfo[c(toInclude),]
        }
    }
    
    ###Filter fastqInfo for empty fastq_ftp entries
    toRemove <- which(fastqInfo[['fastq_ftp']] == '')
    
    if(length(toRemove) > 0){
        cat('Samples removed because of missing raw data:')
        print(fastqInfo[toRemove,])
        fastqInfo <- fastqInfo[-c(toRemove),]
    }
    
    ##Check sequencing mode and keep only supported mode in mixed cases
    pairedDataPos <- grep(';',fastqInfo$fastq_ftp)
    if(length(pairedDataPos) > 0 & length(pairedDataPos) != nrow(fastqInfo)){
        if(param$supportedMode == 'paired'){
            fastqInfo <- fastqInfo[pairedDataPos,]
        } else if(param$supportedMode == 'single'){
            fastqInfo <- fastqInfo[-c(pairedDataPos),]
        }
    }
    
    
    setwdNew(rownames(input$meta))
    
    fastqInfo$Name <- fastqInfo$sample_accession ## use the accession as a fallback
    for (i in 1:nrow(fastqInfo)){
        #download ERR xml File
        cmd = paste0("curl -o ", fastqInfo$run_accession[i],".xml ", "-X GET \'https://www.ebi.ac.uk/ena/browser/api/xml/",fastqInfo$run_accession[i],"?download=true\'")
        ezSystem(cmd)
        #Extract read number from xml, extract sampleID from xml
        runInfo <- xmlParse(paste0(fastqInfo$run_accession[i], '.xml'))
        
        runAttr = xmlToList(runInfo)$RUN$RUN_ATTRIBUTES
        for (k in 1:length(runAttr)){
            if(runAttr[[k]]$TAG == 'ENA-SPOT-COUNT'){
                fastqInfo[['ReadCount']][i] <- runAttr[[k]]$VALUE}
            }
        sampleID <- xmlToList(runInfo)$RUN$RUN_LINKS[[2]]$XREF_LINK$ID
        
        if(is.list(sampleID)|any(grepl('^ERP',sampleID))){
            sampleID <- fastqInfo$sample_accession[i]
        }
        cmd = paste0("curl -o ", sampleID,".xml ", "-X GET \'https://www.ebi.ac.uk/ena/browser/api/xml/",sampleID,"?download=true\'")
        ezSystem(cmd)
        xml <- xmlParse(paste0(sampleID, '.xml'))
        sampleInfo <- xmlToList(xml)
        
        if (!is.null(sampleInfo$SAMPLE$TITLE)){
          cleanName <- gsub("[\\# \\(\\):,;]", "_", sampleInfo$SAMPLE$TITLE)
          fastqInfo[['Name']][i] <- paste(cleanName, sampleID, sep = '_')
          sampleAttributes <- xmlToDataFrame(xmlRoot(xml)[[1]][[5]])
          sampleAttributes <- sampleAttributes[grep('ENA',sampleAttributes$TAG, invert = TRUE),]
          
          if(!is.null(nrow(sampleAttributes)) && nrow(sampleAttributes) > 0){
            for (j in 1:nrow(sampleAttributes)){
              attrName <- as.character(sampleAttributes[j, 1])
              if (is.null(fastqInfo[[attrName]])){
                fastqInfo[[attrName]] = '-'
              }
              fastqInfo[i, attrName] = as.character(sampleAttributes[j,2])
            }
          }
        } else {
            tryCatch({
                experimentID <- xmlToList(runInfo)[[1]]$EXPERIMENT_REF$.attrs
                if(length(experimentID) > 1){
                    experimentID <- experimentID['accession']
                }
                cmd = paste0("curl -o ", experimentID,".xml ", "-X GET \'https://www.ebi.ac.uk/ena/browser/api/xml/",experimentID,"?download=true\'")
                ezSystem(cmd)
                xml <- xmlParse(paste0(experimentID, '.xml'))
                experimentInfo <- xmlToList(xml)
                cleanName <-  gsub("[\\# \\(\\):,;]", "_", experimentInfo$EXPERIMENT$TITLE)
                fastqInfo[['Name']][i] <- paste(cleanName, sampleID, sep = '_')},
                warning=function(w) {
                    message('sample title not available in experiment accession file')
                    print(w)
                    return(NA)
                })
            }
        
        if(grepl(';',fastqInfo$fastq_ftp[i])) {
            paired = TRUE
            fastqList <- strsplit(fastqInfo$fastq_ftp[i],';')
            cmd = paste0('wget --quiet -t 0 ', 'ftp://', fastqList[[1]][1],'; wget --quiet -t 0 ', 'ftp://', fastqList[[1]][2])
        } else {
            paired = FALSE
            cmd = paste0('wget --quiet -t 0 ','ftp://', fastqInfo$fastq_ftp[i])
        }
        ezSystem(cmd)
    }
    myPath = output$meta[['ENA Result [File]']]
    dataset <- createDataset(fastqInfo, myPath, paired = paired)
    
    if (sum(dataset[['Read Count']])==0) {
        dataset[['Read Count']] <- countReadsInFastq(basename(fastqInfo$fastq_ftp))
    }
    
    
    if(length(dataset$Name) == length(unique(dataset$Name))){
        ezWrite.table(dataset, 'dataset.tsv', row.names = FALSE)
    } else {
        message('Data pooling needed')
        sampleNames = unique(dataset$Name)
        for (j in 1:length(sampleNames)){
            files_R1 <- basename(dataset[['Read1 [File]']][dataset[['Name']] == sampleNames[j]])
            sampleNames[j] = gsub('\\/', '_', gsub(' ', '_', sampleNames[j]))
            sampleNames[j] = gsub(';', '_', sampleNames[j])
            pooledR1_File <- paste0(sampleNames[j],'_R1.fastq.gz')
            ezSystem(paste('touch', pooledR1_File))
            for (k in 1:length(files_R1)){
                ezSystem(paste0('cat ', files_R1[k],'>>', pooledR1_File))
                ezSystem(paste('rm', files_R1[k]))
                dataset[['Read1 [File]']] = sub(files_R1[k], pooledR1_File, dataset[['Read1 [File]']])
            }
            if(paired){
                files_R2 <- basename(dataset[['Read2 [File]']][dataset[['Name']] == sampleNames[j]])
                sampleNames[j] = gsub('\\/', '_', gsub(' ', '_', sampleNames[j]))
                sampleNames[j] = gsub(';', '_', sampleNames[j])
                pooledR2_File <- paste0(sampleNames[j],'_R2.fastq.gz')
                ezSystem(paste('touch', pooledR2_File))
                for (k in 1:length(files_R2)){
                    ezSystem(paste0('cat ', files_R2[k],'>>', pooledR2_File))
                    ezSystem(paste('rm', files_R2[k]))
                    dataset[['Read2 [File]']] = sub(files_R2[k], pooledR2_File, dataset[['Read2 [File]']])
                }
            }
            dataset[['Read Count']][dataset[['Name']] == sampleNames[j]] = rep(sum(as.numeric(dataset[['Read Count']][dataset[['Name']] == sampleNames[j]])), length(files_R1))
        }
        dataset = dataset[,-grep('md5sum', colnames(dataset))]
        dataset = unique(dataset)
        ##TODO: Recompute md5sums after pooling & handle PairedEnd data
        ezWrite.table(dataset, 'dataset.tsv', row.names = FALSE)
    }
     
    if(param$tarOutput){
    for (i in 1:nrow(dataset)){
        mySample <- dataset[i,'Name']
        dir.create(mySample)
        R1File <-  basename(dataset[['Read1 [File]']][i])
        R2File <- basename(dataset[['Read2 [File]']][i])
        cmd <- paste('mv', file.path(R1File), file.path(R2File), mySample)
        system(cmd)
        ##Rename files
        system(paste('mv',  file.path(mySample, basename(R1File)), file.path(mySample, paste0(mySample, '_S1_L001_R1_001.fastq.gz'))))
        system(paste('mv',  file.path(mySample, basename(R2File)), file.path(mySample, paste0(mySample, '_S1_L001_R2_001.fastq.gz'))))
        ##tar folder
        cmd <- paste('tar vcf', paste0(mySample, '.tar'), mySample)
        system(cmd)
        ##rm folder
        system(paste('rm -Rf', mySample))
    }
    
    ##Update dataset
    dataset[['RawDataDir [File]']] <- file.path(dirname(dataset[['Read1 [File]']])[1], paste0(dataset$Name, '.tar'))
    dataset[['Read1 [File]']] <- NULL
    dataset[['Read2 [File]']] <- NULL
    ezWrite.table(dataset, 'dataset.tsv', row.names = FALSE)
    }
    
    ##Register dataset/resources in BF as short term storage
    datasetName = paste0('ENA_App_', output$getColumn('projectID'))
    myCmd <- paste('register_sushi_dataset_into_bfabric', containerId, 'dataset.tsv', datasetName, '-b ~/.bfabricpy.yml --skip-file-check -a 372')
    tryCatch({
    ezSystem(myCmd)}, warning=function(w) {
        message('please in B-Fabric if registration was succesful')
        print(w)
        return(NA)
    })
}

createDataset <- function(fastqInfo, myPath, paired = FALSE){
    if(ncol(fastqInfo) > 10){
        xtraCols <- fastqInfo[,11:ncol(fastqInfo)]
        colnames(fastqInfo[,11:ncol(fastqInfo)])
    } else {
        xtraCols <- data.frame(sample_accession = fastqInfo$sample_accession)
        
    }
    if (!paired){
        dataset = data.frame(Name = gsub(' ', '_', fastqInfo$Name), Read1 = file.path(myPath, basename(fastqInfo$fastq_ftp)), 
                             md5sum = fastqInfo$fastq_md5, Species = fastqInfo$Species, 
                             ReadCount = fastqInfo$ReadCount, xtraCols, stringsAsFactors = FALSE)
        colnames(dataset) = c('Name', 'Read1 [File]', 'md5sum', 'Species', 'Read Count', colnames(xtraCols))
    } else {
        dataset = data.frame(Name = gsub(' ', '_', fastqInfo$Name), Read1 = file.path(myPath, sapply(strsplit(fastqInfo$fastq_ftp, ';'), basename)[1,]),
                             Read2 = file.path(myPath, sapply(strsplit(fastqInfo$fastq_ftp, ';'), basename)[2,]),
                             md5sum = fastqInfo$fastq_md5, Species = fastqInfo$Species, 
                             ReadCount = fastqInfo$ReadCount, xtraCols, stringsAsFactors = FALSE)
        colnames(dataset) = c('Name', 'Read1 [File]', 'Read2 [File]', 'md5sum', 'Species', 'Read Count', colnames(xtraCols))
    }
    return(dataset)
}

##' @template app-template
##' @templateVar method ezMethodGetEnaData(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppENA <-
    setRefClass("EzAppENA",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodGetEnaData
                        name <<- "EzAppENA"
                        appDefaults <<- rbind(
                            excludedSamples = ezFrame(
                                Type = "character",
                                DefaultValue = "",
                                Description = 'sample_accession SAMN-ids for samples to exclude'),
                            includedSamples = ezFrame(
                                Type = "character",
                                DefaultValue = "",
                                Description = 'sample_accession SAMN-ids for samples to include')
                            )
                    }
                )
    )
