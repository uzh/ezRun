ezMethodGetEnaData <- function(input=NA, output=NA, param=NA){
    #param[['projectID']] = 'PRJNA163241'
    require(XML)
    cmd = paste0("curl -o fastqLinks.txt -X GET ","\'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",param[['projectID']],"&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,fastq_md5&format=tsv&download=true\'")
    ezSystem(cmd)
    fastqInfo = ezRead.table('fastqLinks.txt', row.names = NULL)
    fastqInfo[['Name']] = ''
    fastqInfo[['ReadCount']] = 0
    colnames(fastqInfo)[grep('scientific_name', colnames(fastqInfo))] = 'Species'
    setwdNew(rownames(input$meta))
    
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
        
        if(is.list(sampleID)){
            sampleID <- fastqInfo$sample_accession[i]
        }
        cmd = paste0("curl -o ", sampleID,".xml ", "-X GET \'https://www.ebi.ac.uk/ena/browser/api/xml/",sampleID,"?download=true\'")
        ezSystem(cmd)
        xml <- xmlParse(paste0(sampleID, '.xml'))
        sampleInfo <- xmlToList(xml)
        fastqInfo[['Name']][i] <- sampleInfo$SAMPLE$TITLE
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
            #ezSystem(paste('pigz --best', pooledR1_File))
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
                #ezSystem(paste('pigz --best', pooledR2_File))
            }
            #files_R2 <- dataset[['Read2 [File]']][dataset[['Name']] == sampleNames[j]]
            dataset[['Read Count']][dataset[['Name']] == sampleNames[j]] = rep(sum(as.numeric(dataset[['Read Count']][dataset[['Name']] == sampleNames[j]])), length(files_R1))
        }
        dataset = dataset[,-grep('md5sum', colnames(dataset))]
        dataset = unique(dataset)
        ##TODO: Recompute md5sums after pooling & handle PairedEnd data
        ezWrite.table(dataset, 'dataset.tsv', row.names = FALSE)
    }
}

createDataset <- function(fastqInfo, myPath, paired = FALSE){
    if (!paired){
        dataset = data.frame(Name = gsub(' ', '_', fastqInfo$Name), Read1 = file.path(myPath, basename(fastqInfo$fastq_ftp)), 
                             md5sum = fastqInfo$fastq_md5, Species = fastqInfo$Species, 
                             ReadCount = fastqInfo$ReadCount, fastqInfo[,11:ncol(fastqInfo)], stringsAsFactors = FALSE)
        colnames(dataset) = c('Name', 'Read1 [File]', 'md5sum', 'Species', 'Read Count', colnames(fastqInfo)[11:ncol(fastqInfo)])
        #dataset = dataset[,-c(6)]
    } else {
        dataset = data.frame(Name = gsub(' ', '_', fastqInfo$Name), Read1 = file.path(myPath, sapply(strsplit(fastqInfo$fastq_ftp, ';'), basename)[1,]),
                             Read2 = file.path(myPath, sapply(strsplit(fastqInfo$fastq_ftp, ';'), basename)[2,]),
                             md5sum = fastqInfo$fastq_md5, Species = fastqInfo$Species, 
                             ReadCount = fastqInfo$ReadCount, fastqInfo[,11:ncol(fastqInfo)], stringsAsFactors = FALSE)
        colnames(dataset) = c('Name', 'Read1 [File]', 'Read2 [File]', 'md5sum', 'Species', 'Read Count', colnames(fastqInfo)[11:ncol(fastqInfo)])
        #dataset = dataset[,-c(7)]
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
                    }
                )
    )
