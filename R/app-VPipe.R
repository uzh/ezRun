###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### V-Pipe
### https://github.com/cbg-ethz/V-pipe

EzAppVPipe <-
    setRefClass("EzAppVPipe",
                contains = "EzApp",
                methods = list(
                    initialize = function()
                    {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodVPipe
                        name <<- "EzAppVPipe"
                        appDefaults <<- rbind(
                            readLength = ezFrame(
                                Type = "integer",
                                DefaultValue = 150,
                                Description = "Read Length"
                            ),
                            samplePrefix = ezFrame(
                                Type = "character",
                                DefaultValue = "",
                                Description = "prefix to remove from sample name"
                            )
                        )
                    }
                )
    )

ezMethodVPipe <- function(input=NA, output=NA, param=NA, 
                               htmlFile="00index.html"){

    dataset <- input$meta
    samples <- rownames(dataset)
    dataDir <- dirname(input$getFullPaths("Read1")[1])
    orderId <- paste0('o', unique(dataset[['Order Id [B-Fabric]']]))
    
    samples <- sub(param$samplePrefix, '', samples)
    samples <- limma::strsplit2(samples, '_')
    if(ncol(samples) == 1){
        samples <- data.frame(samples)
        samples[['InternalNumber']] = paste(samples[,1], '_S', 1:nrow(samples), sep = '')
    }
    samples <- data.frame(ID = samples[,1], InternalNumber = samples[,2], RL = param$readLength, stringsAsFactors = FALSE)
    rownames(samples) <- rownames(dataset)
    ezWrite.table(samples, 'samples.tsv', row.names = FALSE, col.names = FALSE)
    
    
    setwdNew('samples')
    oDir <- "outputDir"
    cDir <- "consensus_majority_dels"
    aDir <- "consensus_ambig_dels"
    vDir <- "variants_lofreq"
    
    dir.create(oDir)
    dir.create(file.path(oDir, cDir))
    dir.create(file.path(oDir, aDir))
    dir.create(file.path(oDir, vDir))
    
    for (i in rownames(dataset)){
        dir.create(samples[i, 'ID'])
        dir.create(file.path(samples[i, 'ID'],samples[i, 'InternalNumber']))
        dest <- file.path(samples[i, 'ID'],samples[i, 'InternalNumber'], 'raw_data')
        dir.create(dest)
        system(paste('ln -s', input$getFullPaths("Read1")[i], file.path(dest, basename(input$getFullPaths("Read1")[i]))))
        system(paste('ln -s', input$getFullPaths("Read2")[i], file.path(dest, basename(input$getFullPaths("Read2")[i]))))
    }
    setwd('..')
    system('ln -s /srv/GT/analysis/p24680/references references')
    system('cp  /srv/GT/analysis/p24680/vpipe.config .')
    system('cp /srv/GT/analysis/p24680/vpipe .')
    dir.create('/scratch/tmp',showWarnings = FALSE)
    cmd <- paste('export TMPDIR=/scratch/tmp; ./vpipe --use-conda -p --keep-going --rerun-incomplete --cores', param$cores, ' 2>&1 | tee', paste0(orderId, '.log'))
    system(cmd)
    gc()
    
    setwd('samples')
    for (j in 1:nrow(samples)){
        cmd <- paste('mv', file.path(samples[j, 'ID'],samples[j,'InternalNumber'], 'references/ref_ambig_dels.fasta'), 
                       paste0(file.path(oDir, aDir),'/',samples[j, 'ID'],'_',samples[j,'InternalNumber'],'.ambig.fasta'))
        ezSystem(cmd)
        cmd <- paste('mv', file.path(samples[j, 'ID'],samples[j,'InternalNumber'], 'references/ref_majority_dels.fasta'), 
                     paste0(file.path(oDir, cDir),'/',samples[j, 'ID'],'_',samples[j,'InternalNumber'],'.fasta'))
        ezSystem(cmd)
        cmd <- paste('mv', file.path(samples[j, 'ID'],samples[j,'InternalNumber'], 'variants/SNVs/snvs.vcf'), 
                     paste0(file.path(oDir, vDir),'/',samples[j, 'ID'],'_',samples[j,'InternalNumber'],'.vcf'))
        ezSystem(cmd)
    }
    setwd(file.path(oDir, cDir))
    
    consensusFile <- paste0(orderId, '.consensus.fasta')
    system(paste0('cat *.fasta >', '../', consensusFile))
    setwd('..')
    createPangolinScript(consensusFile)
    system('bash runPangolin.sh')
    file.remove('runPangolin.sh')
    gc()
    
    lineageReport <- ezRead.table('lineage_report.csv', sep = ',', row.names = NULL)
    x = DT::datatable(lineageReport, escape = FALSE, rownames = FALSE, filter = 'bottom',
                  caption = paste('',sep=''),extensions = c('Buttons'),
                  options = list(
                                 initComplete = DT::JS(
                                     "function(settings, json) {",
                                     "$(this.api().table().header()).css({'background-color': '#0000A0', 'color': '#fff'});",
                                     "}"),
                                 dom = c('Bfrtip'),buttons = c('colvis','copy', 'csv', 'excel', 'pdf', 'print'), pageLength=200, autoWidth=TRUE))
    DT::saveWidget(x, 'lineage_report.html')
    
    dir.create('bedtools_coverage')
    dir.create('samtools_depth')
    for (j in 1:nrow(samples)){
        cmd <- paste('bedtools coverage -a /srv/gstore/projects/p24680/gffs/artic.expected_amplicons.bed -b', 
                     file.path('..', samples[j, 'ID'],samples[j,'InternalNumber'], 'variants/SNVs/REF_aln_indelqual.bam'), '-mean -bed >', file.path('bedtools_coverage', paste0(samples[j, 'ID'],'_',samples[j,'InternalNumber'], '.bed')))
        ezSystem(cmd)             
        cmd <- paste('samtools depth -a', file.path('..', samples[j, 'ID'],samples[j,'InternalNumber'], 'variants/SNVs/REF_aln_indelqual.bam >'), file.path('samtools_depth', paste0(samples[j, 'ID'],'_',samples[j,'InternalNumber'], '.tsv')))
        ezSystem(cmd)
    }
    setwd('../..')
    dir.create(param[['name']])
    system(paste('mv samples/*', param[['name']]))
    return('success')
}

createPangolinScript <- function(fileName){
    cat('#!/bin/bash \n',  file = 'runPangolin.sh')
    cat('. /usr/local/ngseq/miniconda3/etc/profile.d/conda.sh \n', file = 'runPangolin.sh', append = TRUE)
    cat('conda activate pangolin \n', file = 'runPangolin.sh', append = TRUE)
    cat(paste('pangolin', fileName), file = 'runPangolin.sh', append = TRUE)
}