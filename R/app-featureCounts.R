###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodFeatureCounts = function(input=NA, output=NA, param=NA){
  require(GenomicRanges)
  require(Rsubread)
  require(rtracklayer)
  
  bamFile = input$getFullPaths("BAM")
  localBamFile = getBamLocally(bamFile)
  if(localBamFile != bamFile){
    on.exit(file.remove(c(localBamFile, paste0(localBamFile, ".bai"))),
            add=TRUE)
  }
  
  outputFile = basename(output$getColumn("Count"))
  statFile = basename(output$getColumn("Stats"))
  
  if(!is.null(param$aroundTSSCounting) && param$aroundTSSCounting){
    gtf <- rtracklayer::import(param$ezRef@refFeatureFile)
    gtf <- gtf[gtf$type == "gene"]
    if (ezIsSpecified(param$transcriptTypes)){
      seqAnno = ezFeatureAnnotation(param$ezRef@refAnnotationFile,
                                    dataFeatureType="gene")
      if( !is.null(seqAnno$type) && any(seqAnno$type != "")){
        genesUse <- rownames(seqAnno)[seqAnno$type %in% param$transcriptTypes]
        gtf <- gtf[gtf$gene_id %in% genesUse]
        
        seqlengthsRef <- fasta.seqlengths(param$ezRef@refFastaFile)
        names(seqlengthsRef) <- sub('[[:blank:]].*$','',names(seqlengthsRef))
        seqlengths(gtf) <- seqlengthsRef[names(seqlengths(gtf))]
        gtf <- trim(suppressWarnings(promoters(gtf, upstream=param$upstreamFlanking,
                                               downstream=param$downstreamFlanking)))
      }
    }
    
    gtfFile = paste(Sys.getpid(), "genes.gtf", sep="-")
    rtracklayer::export(gtf, gtfFile)
    on.exit(file.remove(gtfFile), add=TRUE)

    countResult = featureCounts(localBamFile, annot.inbuilt=NULL,
                                          annot.ext=gtfFile, isGTFAnnotationFile=TRUE,
                                          GTF.featureType='gene',
                                          GTF.attrType= 'gene_id',
                                          useMetaFeatures=param$useMetaFeatures,
                                          allowMultiOverlap=param$allowMultiOverlap,
                                          isPairedEnd=param$paired, 
                                          requireBothEndsMapped=FALSE,
                                          checkFragLength=FALSE,minFragLength=50,
                                          maxFragLength=600,
                                          nthreads=param$cores, 
                                          strandSpecific=switch(param$strandMode, "both"=0, "sense"=1, "antisense"=2, stop("unsupported strand mode: ", param$strandMode)),
                                          minMQS=param$minMapQuality,
                                          readExtension5=0,readExtension3=0,
                                          read2pos=NULL,
                                          minOverlap=param$minFeatureOverlap,
                                          ignoreDup=param$ignoreDup,
                                          splitOnly=FALSE,
                                          countMultiMappingReads=param$keepMultiHits,
                                          fraction=param$keepMultiHits & !param$countPrimaryAlignmentsOnly,
                                          primaryOnly=param$countPrimaryAlignmentsOnly,
                                          countChimericFragments=TRUE,
                                          chrAliases=NULL,reportReads=NULL)
  }else{
    ## Count exons by gene
    gtfFile = param$ezRef@refFeatureFile
    if (ezIsSpecified(param$transcriptTypes)){
      seqAnno = ezFeatureAnnotation(param$ezRef@refAnnotationFile,
                                    dataFeatureType="transcript")
      if ( !is.null(seqAnno$type) && any(seqAnno$type != "")){
        transcriptsUse = rownames(seqAnno)[seqAnno$type %in% param$transcriptTypes]
        require(data.table)
        require(stringr)
        gtf <- ezReadGff(param$ezRef@refFeatureFile)
        transcripts <- ezGffAttributeField(gtf$attributes,
                                           field="transcript_id", 
                                           attrsep="; *", valuesep=" ")
        gtf = gtf[transcripts %in% transcriptsUse, ]
        ## write.table is much faster than write_tsv and export.
        ### write_tsv: 38.270s
        ### export: 347.998s
        ### write.table: 8.787s
        gtfFile = paste(Sys.getpid(), "genes.gtf", sep="-")
        on.exit(file.remove(gtfFile), add=TRUE)
        write.table(gtf, gtfFile, quote=FALSE, sep="\t", 
                    row.names=FALSE, col.names=FALSE)
      }
    }
    countResult = featureCounts(localBamFile, annot.inbuilt=NULL,
                              annot.ext=gtfFile, isGTFAnnotationFile=TRUE,
                              GTF.featureType=param$gtfFeatureType,
                              GTF.attrType=switch(param$featureLevel,
                                                  "gene"="gene_id",
                                                  "transcript"="transcript_id",
                                                  "isoform"="transcript_id",
                                                  stop("unsupported feature level: ", param$featureLevel)),
                              useMetaFeatures=param$useMetaFeatures,
                              allowMultiOverlap=param$allowMultiOverlap,
                              isPairedEnd=param$paired, 
                              requireBothEndsMapped=FALSE,
                              checkFragLength=FALSE,minFragLength=50,
                              maxFragLength=600,
                              nthreads=param$cores, 
                              strandSpecific=switch(param$strandMode, "both"=0, "sense"=1, "antisense"=2, stop("unsupported strand mode: ", param$strandMode)),
                              minMQS=param$minMapQuality,
                              readExtension5=0,readExtension3=0,read2pos=NULL,
                              minOverlap=param$minFeatureOverlap,
                              ignoreDup=param$ignoreDup,
                              splitOnly=FALSE,
                              countMultiMappingReads=param$keepMultiHits,
                              fraction=param$keepMultiHits & !param$countPrimaryAlignmentsOnly,
                              primaryOnly=param$countPrimaryAlignmentsOnly,
                              countChimericFragments=TRUE,chrAliases=NULL,
                              reportReads=NULL)
  }
  
  colnames(countResult$counts) = "matchCounts"
  ezWrite.table(countResult$counts, file=outputFile)
  colnames(countResult$stat) = c("Status", "Count")
  ezWrite.table(countResult$stat, file=statFile, row.names=FALSE)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFeatureCounts(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppFeatureCounts <-
  setRefClass("EzAppFeatureCounts",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFeatureCounts
                  name <<- "EzAppFeatureCounts"
                  appDefaults <<- rbind(gtfFeatureType = ezFrame(Type="character",  DefaultValue="exon",  Description="which gtf feature types to use; with Ensembl GTF files; use 'transcript' to count also intronic reads"),
                                        allowMultiOverlap = ezFrame(Type="logical",  DefaultValue="TRUE",  Description="make sure that every read counts only as one"),
                                        countPrimaryAlignmentsOnly = ezFrame(Type="logical",  DefaultValue="TRUE",  Description="count only the primary alignment"),
                                        minFeatureOverlap=ezFrame(Type="integer", DefaultValue="10", Description="the number of bases overlap are need to generate a count"),
                                        useMetaFeatures=ezFrame(Type="logical", DefaultValue="TRUE", Description="should counts be summarized to meta-features"),
                                        ignoreDup=ezFrame(Type="logical", DefaultValue="FALSE", Description="ignore reads marked as duplicates"),
                                        aroundTSSCounting=ezFrame(Type="logical",  DefaultValue="FALSE",  Description="count reads only around TSS of genes"),
                                        downstreamFlanking=ezFrame(Type="integer", DefaultValue="250", Description="the number of bases downstream of TSS"),
                                        upstreamFlanking=ezFrame(Type="integer", DefaultValue="250", Description="the number of bases upstream of TSS"))
                }
              )
  )


### EzAppSingleCellFeatureCounts
EzAppSingleCellFeatureCounts <-
  setRefClass("EzAppSingleCellFeatureCounts",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSingleCellFeatureCounts
                  name <<- "EzAppSingleCellFeatureCounts"
                  appDefaults <<- 
                    rbind(gtfFeatureType = ezFrame(Type="character",
                                                   DefaultValue="exon",
                                                   Description="which gtf feature types to use; with Ensembl GTF files; use 'transcript' to count also intronic reads"),
                          allowMultiOverlap = ezFrame(Type="logical",
                                                      DefaultValue="TRUE",
                                                      Description="make sure that every read counts only as one"),
                          countPrimaryAlignmentsOnly = ezFrame(Type="logical",
                                                               DefaultValue="TRUE",
                                                               Description="count only the primary alignment"),
                          minFeatureOverlap=ezFrame(Type="integer",
                                                    DefaultValue="10",
                                                    Description="the number of bases overlap are need to generate a count"),
                          useMetaFeatures=ezFrame(Type="logical",
                                                  DefaultValue="TRUE",
                                                  Description="should counts be summarized to meta-features"),
                          ignoreDup=ezFrame(Type="logical",
                                            DefaultValue="FALSE",
                                            Description="ignore reads marked as duplicates"),
                          aroundTSSCounting=ezFrame(Type="logical",
                                                    DefaultValue="FALSE",
                                                    Description="count reads only around TSS of genes"),
                          downstreamFlanking=ezFrame(Type="integer",
                                                     DefaultValue="250",
                                                     Description="the number of bases downstream of TSS"),
                          upstreamFlanking=ezFrame(Type="integer",
                                                   DefaultValue="250",
                                                   Description="the number of bases upstream of TSS"),
                          controlSeqs=ezFrame(Type="charVector",
                                              DefaultValue="",
                                              Description="control sequences to add")
                          )
                }
              )
  )

ezMethodSingleCellFeatureCounts <- function(input=NA, output=NA, param=NA){
  require(GenomicRanges)
  require(Rsubread)
  require(rtracklayer)
  
  bamFile = input$getFullPaths("BAM")
  localBamFile = getBamLocally(bamFile)
  if(localBamFile != bamFile){
    on.exit(file.remove(c(localBamFile, paste0(localBamFile, ".bai"))),
            add=TRUE)
  }
  
  outputFile = basename(output$getColumn("CountMatrix"))
  statFile = basename(output$getColumn("Stats"))
  
  require(Rsamtools)
  bamHeaders <- scanBamHeader(localBamFile)
  hasRG <- "@RG" %in% names(bamHeaders[[1]]$text)
  
  if(!is.null(param$aroundTSSCounting) && param$aroundTSSCounting){
    gtf <- rtracklayer::import(param$ezRef@refFeatureFile)
    seqlengthsRef <- fasta.seqlengths(param$ezRef@refFastaFile)
    seqlengths(gtf) <- seqlengthsRef[names(seqlengths(gtf))]
    gtf <- gtf[gtf$type == "gene"]
    if (ezIsSpecified(param$transcriptTypes)){
      seqAnno = ezFeatureAnnotation(param$ezRef@refAnnotationFile,
                                    dataFeatureType="gene")
      genesUse <- rownames(seqAnno)[seqAnno$type %in% param$transcriptTypes]
      gtf <- gtf[gtf$gene_id %in% genesUse]
    }
    
    gtfFile = paste(Sys.getpid(), "genes.gtf", sep="-")
    gtf <- trim(suppressWarnings(promoters(gtf, upstream=param$upstreamFlanking,
                                           downstream=param$downstreamFlanking)))
    rtracklayer::export(gtf, gtfFile)
    on.exit(file.remove(gtfFile), add=TRUE)
    
    countResult = featureCounts(localBamFile, annot.inbuilt=NULL,
                                          annot.ext=gtfFile, isGTFAnnotationFile=TRUE,
                                          GTF.featureType='gene',
                                          GTF.attrType= 'gene_id',
                                          useMetaFeatures=param$useMetaFeatures,
                                          allowMultiOverlap=param$allowMultiOverlap,
                                          isPairedEnd=param$paired, 
                                          requireBothEndsMapped=FALSE,
                                          checkFragLength=FALSE,minFragLength=50,
                                          maxFragLength=600,
                                          nthreads=param$cores, 
                                          strandSpecific=switch(param$strandMode, "both"=0, "sense"=1, "antisense"=2, stop("unsupported strand mode: ", param$strandMode)),
                                          minMQS=param$minMapQuality,
                                          readExtension5=0,readExtension3=0,
                                          read2pos=NULL,
                                          minOverlap=param$minFeatureOverlap,
                                          ignoreDup=param$ignoreDup,
                                          splitOnly=FALSE,
                                          countMultiMappingReads=param$keepMultiHits,
                                          fraction=param$keepMultiHits & !param$countPrimaryAlignmentsOnly,
                                          primaryOnly=param$countPrimaryAlignmentsOnly,
                                          countChimericFragments=TRUE,
                                          chrAliases=NULL,reportReads=NULL,
                                          byReadGroup=ifelse(hasRG, TRUE, FALSE))
  }else{
    gtfFile = paste(Sys.getpid(), "genes.gtf", sep="-")
    on.exit(file.remove(gtfFile), add=TRUE)
    if (ezIsSpecified(param$transcriptTypes)){
      seqAnno = ezFeatureAnnotation(param$ezRef@refAnnotationFile,
                                    dataFeatureType="transcript")
      transcriptsUse = rownames(seqAnno)[seqAnno$type %in% param$transcriptTypes]
      require(data.table)
      require(stringr)
      gtf <- ezReadGff(param$ezRef@refFeatureFile)
      transcripts <- ezGffAttributeField(gtf$attributes,
                                         field="transcript_id", 
                                         attrsep="; *", valuesep=" ")
      gtf = gtf[transcripts %in% transcriptsUse, ]
      write.table(gtf, gtfFile, quote=FALSE, sep="\t", 
                  row.names=FALSE, col.names=FALSE)
    } else {
      file.copy(from=param$ezRef@refFeatureFile, to=gtfFile)
    }
    
    if(ezIsSpecified(param$controlSeqs)){
      ## control sequences
      genomesRoot <- strsplit(GENOMES_ROOT, ":")[[1]]
      controlSeqsFn <- file.path(genomesRoot, "controlSeqs.fa")
      controlSeqsFn <- head(controlSeqsFn[file.exists(controlSeqsFn)], 1)
      controlSeqs <- readDNAStringSet(controlSeqsFn)
      names(controlSeqs) <- sub(" .*$", "", names(controlSeqs))
      controlSeqs <- controlSeqs[param$controlSeqs]
      txids <- rep(paste0("Transcript_",names(controlSeqs)), each=4)
      txids[seq(1, length(txids), by=4)] <- NA
      transcript_biotype <- rep("protein_coding", length(txids))
      transcript_biotype[seq(1, length(transcript_biotype), by=4)] <- NA
      mcols=DataFrame(source="NCBI", 
                      type=rep(c("gene", "transcript",
                                 "exon", "CDS"), length(controlSeqs)),
                      score=NA, phase=rep(c(NA, NA, NA, 0), length(controlSeqs)),
                      gene_id=rep(paste0("Gene_",names(controlSeqs)),
                                  each=4),
                      gene_version=NA,
                      gene_name=rep(names(controlSeqs), each=4),
                      gene_source="NCBI",
                      gene_biotype="protein_coding",
                      havana_gene=NA, havana_gene_version=NA,
                      transcript_id=txids, transcript_version=NA,
                      transcript_name=txids, transcript_source=NA,
                      transcript_biotype=transcript_biotype,
                      havana_transcript=NA,
                      havana_transcript_version=NA,
                      tag=NA, transcript_support_level=NA,
                      exon_number=NA, exon_id=NA, exon_version=NA,
                      ccds_id=NA, protein_id=NA, protein_version=NA
      )
      extraGR <- GRanges(seqnames=rep(names(controlSeqs), each=4),
                         ranges=IRanges(start=1,
                                        end=rep(width(controlSeqs), each=4)),
                         strand="+")
      mcols(extraGR) <- mcols
      gtfExtraFn <- tempfile(pattern="extraSeqs", tmpdir=getwd(),
                             fileext = ".gtf")
      on.exit(file.remove(gtfExtraFn), add=TRUE)
      export.gff2(extraGR, con=gtfExtraFn)
      ezSystem(paste("cat", gtfExtraFn, ">>", gtfFile))
    }
    
    countResult = featureCounts(localBamFile, annot.inbuilt=NULL,
                                annot.ext=gtfFile, isGTFAnnotationFile=TRUE,
                                GTF.featureType=param$gtfFeatureType,
                                GTF.attrType=switch(param$featureLevel,
                                                    "gene"="gene_id",
                                                    "transcript"="transcript_id",
                                                    "isoform"="transcript_id",
                                                    stop("unsupported feature level: ", param$featureLevel)),
                                useMetaFeatures=param$useMetaFeatures,
                                allowMultiOverlap=param$allowMultiOverlap,
                                isPairedEnd=param$paired, 
                                requireBothEndsMapped=FALSE,
                                checkFragLength=FALSE,minFragLength=50,
                                maxFragLength=600,
                                nthreads=param$cores,
                                strandSpecific=switch(param$strandMode, "both"=0, "sense"=1, "antisense"=2, stop("unsupported strand mode: ", param$strandMode)),
                                minMQS=param$minMapQuality,
                                readExtension5=0,readExtension3=0,read2pos=NULL,
                                minOverlap=param$minFeatureOverlap,
                                ignoreDup=param$ignoreDup,
                                splitOnly=FALSE,
                                countMultiMappingReads=param$keepMultiHits,
                                fraction=param$keepMultiHits & !param$countPrimaryAlignmentsOnly,
                                primaryOnly=param$countPrimaryAlignmentsOnly,
                                countChimericFragments=TRUE,chrAliases=NULL,
                                reportReads=NULL,
                                byReadGroup=ifelse(hasRG, TRUE, FALSE))
  }
  
  ## The count matrix from featurecounts has colnames messed up
  ## recover them here
  countsFixed <- countResult$counts
  colnames(countsFixed) <- sub(paste0(make.names(normalizePath(localBamFile)), "."), "",
                               colnames(countsFixed), fixed=TRUE)
  #colnames(countsFixed) <- sub(paste0(".*", make.names(basename(localBamFile)), "\\."), "",
  #                             colnames(countsFixed))
  tagsRG <- sub("ID:", "",
                sapply(bamHeaders[[1]]$text[names(bamHeaders[[1]]$text) == "@RG"], "[", 1))
  ## RG starts with numbers will have X after make.names
  ## But featureCounts doesn't have this X.
  fixNameMapping <- setNames(tagsRG, make.names(tagsRG))
  indexStartNumber <- grep("^\\d", fixNameMapping)
  names(fixNameMapping)[indexStartNumber] <- sub("^X", "", 
                                                 names(fixNameMapping)[indexStartNumber])
  colnames(countsFixed) <- fixNameMapping[colnames(countsFixed)]
  
  ## wirteMM doesn't hold the colnames and rownames in mtx
  # ezWrite.table(countsFixed, head=paste0(param$featureLevel, "_id"),
  #               file=outputFile)
  writeSCMM(countsFixed, file=outputFile)
  ezWrite.table(countResult$stat, file=statFile, row.names=FALSE)
  
  # Determine cell cycle phases. The training data is only available for Hsap and Mmus Ensembl
  trainData = NULL
  if (startsWith(param$refBuild, "Homo_sapiens/Ensembl")) {
    trainData = readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                    package = "scran", mustWork=TRUE))
  } else if (startsWith(param$refBuild, "Mus_musculus/Ensembl")) {
    trainData = readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                    package = "scran", mustWork=TRUE))
  }
  if (!is.null(trainData)) {
    cellCycleData = scran::cyclone(countResult$counts, trainData)
    cellPhase = data.frame(Name = colnames(countResult$counts),
                           Phase = cellCycleData$phases)
    write.table(cellPhase, file = basename(output$getColumn('CellCyclePhase')),
                quote = F, sep = "\t", row.names = F)
  }
  
  return("Success")
}
