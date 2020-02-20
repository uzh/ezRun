###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezFeatureAnnotation = function(param, ids=NULL,
                               dataFeatureType=c("gene", "transcript", 
                                                 "isoform")){
  require(data.table)
  require(rtracklayer)
  require(tibble)
  require(readr)
  dataFeatureType <- match.arg(dataFeatureType)
  if(is.list(param)){
    featAnnoFn <- param$ezRef["refAnnotationFile"]
  }else{
    featAnnoFn <- param
  }
  
  if(dataFeatureType == "gene"){
    refAnnoGeneFn <- sub("(_byTranscript)*\\.txt$", "_byGene.txt", featAnnoFn)
    
    if(file.exists(refAnnoGeneFn)){
      message("Using gene level annotation: ", refAnnoGeneFn)
      seqAnno <- fread(refAnnoGeneFn, data.table=FALSE)
      ## replace the NA character columns with ""; 
      ## numeric columns shall not have NA.
      seqAnno[is.na(seqAnno)] <- ""
      rownames(seqAnno) <- seqAnno$gene_id
    }else{
      message("Using isoform level annotation and aggregating.")
      ## For compatibility of old annotation without _byGene.txt
      seqAnnoTx <- fread(featAnnoFn)
      ## replace the NA character columns with ""; 
      ## numeric columns shall not have NA.
      seqAnnoTx[is.na(seqAnnoTx)] <- ""
      ## historical reason: replace Identifier with transcript_id
      colnames(seqAnnoTx)[colnames(seqAnnoTx)=="Identifier"] <- "transcript_id"
      seqAnno <- aggregateFeatAnno(seqAnnoTx)
    }
  }else if(dataFeatureType %in% c("transcript", "isoform")){
    seqAnno <- fread(featAnnoFn, data.table=FALSE)
    ## replace the NA character columns with ""; 
    ## numeric columns shall not have NA.
    seqAnno[is.na(seqAnno)] <- ""
    ## historical reason: replace Identifier with transcript_id
    colnames(seqAnno)[colnames(seqAnno)=="Identifier"] <- "transcript_id"
    rownames(seqAnno) <- seqAnno$transcript_id
  }else{
    stop("Only support dataFeatureType in 'transcript', 'isoform', 'gene'")
  }
  minimalCols <- c("gene_id", "transcript_id", "gene_name", "type", "strand",
                   "seqid", "biotypes", "description", "start", "end",
                   "gc", "featWidth", "GO BP", "GO MF", "GO CC")
  if(!"featWidth" %in% colnames(seqAnno) && "width" %in% colnames(seqAnno)){
    # For back compatibility
    seqAnno$featWidth <- seqAnno$width
    seqAnno$width <- NULL
  }
  if(!"description" %in% colnames(seqAnno) || all(seqAnno$description == "")){
    message("Assigning description with gene_id.")
    seqAnno$description <- seqAnno$gene_id
  }
  if(!"type" %in% colnames(seqAnno) || all(seqAnno$type == "")){
    message("Assigning type with protein coding.")
    seqAnno$type <- "protein_coding"
  }
  if(!"biotypes" %in% colnames(seqAnno) || all(seqAnno$biotypes == "")){
    message("Assigning biotypes with protein coding.")
    seqAnno$biotypes <- "protein_coding"
  }
  if(!all(minimalCols %in% colnames(seqAnno))){
    stop(paste(minimalCols[!minimalCols %in% colnames(seqAnno)], collapse="; "),
         " must exist in annotation file!")
  }
  if(!is.null(ids)){
    if(all(ids %in% rownames(seqAnno))){
      seqAnno <- seqAnno[ids, , drop=FALSE]
    }else{
      extraIds <- setdiff(ids, rownames(seqAnno))
      extraIds <- sub("^(Gene|Transcript)_", "", extraIds)
      # there is control sequences added.
      extraGR <- makeExtraControlSeqGR(extraIds)
      extraSeqs <- getControlSeqs(extraIds)
      gtfExtraFn <- tempfile(fileext = ".gtf")
      export.gff2(extraGR, con=gtfExtraFn)
      controlSeqsTempFn <- tempfile(fileext = ".fa")
      writeXStringSet(getControlSeqs(extraIds),
                      filepath=controlSeqsTempFn)
      extraMart <- tibble("Transcript stable ID"=unique(extraGR$transcript_id),
                          "Gene description"=NA, "GO term accession"=NA,
                          "GO domain"=NA)
      biomartFile <- tempfile(fileext = ".mart")
      write_tsv(extraMart, path=biomartFile)
      
      featAnnoExtra <- makeFeatAnnoEnsembl(featureFile=gtfExtraFn,
                                           genomeFile=controlSeqsTempFn,
                                           biomartFile=biomartFile)
      featAnnoExtra <- featAnnoExtra[[dataFeatureType]]
      seqAnno <-  rbind(seqAnno, featAnnoExtra)[ids, , drop=FALSE]
    }
  }
  return(seqAnno)
}

### -----------------------------------------------------------------
### make the feature annotation file <name>_annotation.txt
### for Ensembl gtf.
makeFeatAnnoEnsembl <- function(featureFile,
                                genomeFile,
                                biomartFile=NULL,
                                organism=NULL,
                                host=NULL,
                                mart='ENSEMBL_MART_ENSEMBL'){
  require(rtracklayer)
  require(data.table)
  
  featAnnoFile <- sub(".gtf", "_annotation_byTranscript.txt", featureFile)
  featAnnoGeneFile <- sub(".gtf", "_annotation_byGene.txt", featureFile)
  
  feature <- import(featureFile)
  transcripts <- feature[feature$type=="transcript"]
  if(length(transcripts) == 0L){
    ## Incomplete gtf with only exons.
    ## Try to reconstruct the transcripts.
    exons <- feature[feature$type == "exon"]
    exonsByTx <- GenomicRanges::split(exons, exons$transcript_id)
    transcripts <- unlist(range(exonsByTx))
    transcripts$transcript_id <- names(transcripts)
    names(transcripts) <- NULL
    transcripts$gene_id <- exons$gene_id[match(transcripts$transcript_id, 
                                               exons$transcript_id)]
    transcripts$gene_name <- exons$gene_name[match(transcripts$transcript_id, 
                                                   exons$transcript_id)]
    transcripts$gene_biotype <- exons$gene_biotype[match(transcripts$transcript_id, 
                                                         exons$transcript_id)]
  }
  
  transcripts <- transcripts[!duplicated(transcripts$transcript_id)]
  ## This is to deal with the cases of duplicates transcripts from GENCODE annotation
  ## Example: ENST00000399012 can be on chrX and chrY.
  ## Ensembl only keeps the ones on chrX.
  
  ## Calculate gc and featWidth
  gw <- getTranscriptGcAndWidth(genomeFn=genomeFile,
                                featureFn=featureFile)
  featAnno <- data.table(transcript_id=transcripts$transcript_id,
                         gene_id=transcripts$gene_id,
                         gene_name=transcripts$gene_name,
                         type=transcripts$gene_biotype,
                         strand=as.character(strand(transcripts)),
                         seqid=as.character(seqnames(transcripts)),
                         start=start(transcripts),
                         end=end(transcripts),
                         biotypes=transcripts$gene_biotype,
                         gc=gw$gc[transcripts$transcript_id],
                         featWidth=gw$featWidth[transcripts$transcript_id]
                        )
  ## The numeric columns should not have NAs.
  stopifnot(!any(is.na(featAnno[ ,.(start, end, gc, featWidth)])))
  
  ## Group the biotype into more general groups
  stopifnot(all(featAnno$biotypes %in% listBiotypes("all")))
  isProteinCoding <- featAnno$biotypes %in% listBiotypes("protein_coding")
  isLNC <- featAnno$biotypes %in% listBiotypes("long_noncoding")
  isSHNC <- featAnno$biotypes %in% listBiotypes("short_noncoding")
  isrRNA <- featAnno$biotypes %in% listBiotypes("rRNA")
  istRNA <- featAnno$biotypes %in% listBiotypes("tRNA")
  isMtrRNA <- featAnno$biotypes %in% listBiotypes("Mt_rRNA")
  isMttRNA <- featAnno$biotypes %in% listBiotypes("Mt_tRNA")
  isPseudo <- featAnno$biotypes %in% listBiotypes("pseudogene")
  featAnno$type[isPseudo] <- "pseudogene"
  featAnno$type[isLNC] <- "long_noncoding"
  featAnno$type[isSHNC] <- "short_noncoding"
  featAnno$type[isProteinCoding] <- "protein_coding"
  ### rRNA and tRNA have to be after noncoding
  ### since they are subset of noncoding
  featAnno$type[isrRNA] <- "rRNA"
  featAnno$type[istRNA] <- "tRNA"
  featAnno$type[isMtrRNA] <- "Mt_rRNA"
  featAnno$type[isMttRNA] <- "Mt_tRNA"
  
  ## additional information from Ensembl or downloaded biomart file
  attributes <- c("ensembl_transcript_id", "description", 
                  "go_id", "namespace_1003")
                  
  names(attributes) <- c("Transcript stable ID", "Gene description",
                         "GO term accession", "GO domain")
  ## Older web-page biomart has different names
  attributesOld <- setNames(attributes,
                            c("Ensembl Transcript ID", "Description",
                              "GO Term Accession", "GO domain"))
  if(!is.null(biomartFile)){
    message("Using local biomart file!")
    ### Use the downloaded biomartFile when availble
    stopifnot(file.exists(biomartFile))
    require(readr)
    
    # fread cannot handle compressed file
    mapping <- as.data.table(read_tsv(biomartFile)) 
    if(all(names(attributes) %in% colnames(mapping))){
      mapping <- mapping[ ,names(attributes), with=FALSE]
      # To make it consistent with biomaRt
      colnames(mapping) <- attributes[colnames(mapping)] 
    }else if(all(names(attributesOld) %in% colnames(mapping))){
      mapping <- mapping[ ,names(attributesOld), with=FALSE]
      # To make it consistent with biomaRt
      colnames(mapping) <- attributesOld[colnames(mapping)]
    }else{
      stop("Make sure ", paste(names(attributes), collapse="; "), 
           "are downloaded from web biomart!")
    }
  }else if(!is.null(organism)){
    ### Query Biomart from R
    message("Query via biomaRt package!")
    require(biomaRt)
    if(is.null(host)){
      ensembl <- useMart(mart)
    }else{
      ensembl <- useMart(mart, host=host)
    }
    ensembl <- useDataset(organism, mart=ensembl)
    mapping1 <-
      getBM(attributes=setdiff(attributes, c("go_id", "namespace_1003")),
            filters=c("ensembl_transcript_id"),
            values=featAnno$transcript_id, mart=ensembl)
    mapping1 <- as.data.table(mapping1)
    mapping2 <-
      getBM(attributes=c("ensembl_transcript_id", "go_id", "namespace_1003"),
            filters=c("ensembl_transcript_id"),
            values=featAnno$transcript_id, mart=ensembl)
    mapping2 <- as.data.table(mapping2)
    mapping <- merge(mapping1, mapping2, all=TRUE)
  }else{
    message("Not using any additional annotation!")
    mapping <- data.table(ensembl_transcript_id=featAnno$transcript_id,
                          description="",
                          go_id="",
                          namespace_1003="")
  }
  
  if(!all(featAnno$transcript_id %in% mapping$ensembl_transcript_id)){
    warning("Some transcript ids don't exist in biomart file!") #Normal for GENCODE
  }
  
  ### description
  txid2description <- mapping[!duplicated(ensembl_transcript_id), 
                              .(transcript_id=ensembl_transcript_id,
                                description=description)]
  
  featAnno <- merge(featAnno, txid2description, all.x=TRUE, all.y=FALSE)
  
  ### GO
  GOMapping <- c("biological_process"="GO BP",
                 "molecular_function"="GO MF",
                 "cellular_component"="GO CC",
                 "ensembl_transcript_id"="transcript_id")
  go <- mapping[!(is.na(go_id) | go_id == "") & namespace_1003 %in% c("biological_process", "molecular_function", "cellular_component"),
                ## They can be NA, "", or weird character (EFO)
                .(go_id=ezCollapse(go_id, na.rm=TRUE, empty.rm=TRUE, 
                                   uniqueOnly=TRUE)),
                by=.(ensembl_transcript_id, namespace_1003)]
  if(nrow(go)==0L){
    ## If go is an empty data.table
    go <- data.table(ensembl_transcript_id=featAnno$transcript_id,
                     biological_process="",
                     molecular_function="",
                     cellular_component="")
  }else{
    go <- dcast(go, ensembl_transcript_id~namespace_1003, value.var ="go_id")
  }
  colnames(go) <- GOMapping[colnames(go)]
  featAnno <- merge(featAnno, go, all.x=TRUE, all.y=FALSE)
  
  featAnno[is.na(featAnno)] <- ""  ## replace NA in GO with ""
  
  ## output annotation file on transcript level
  ezWrite.table(featAnno, file=featAnnoFile, row.names=FALSE)
  cwd <- getwd()
  ### For compatibility, create _annotation.txt symlink
  setwd(dirname(featAnnoFile))
  setwd(cwd)
  
  ## make annotation at gene level
  featAnnoGene <- aggregateFeatAnno(featAnno)
  ezWrite.table(featAnnoGene, file=featAnnoGeneFile, row.names=FALSE)
  
  invisible(list("transcript"=featAnno, "gene"=featAnnoGene))
}

##' @describeIn ezFeatureAnnotation Aggregates the Go annotation.
aggregateGoAnnotation = function(seqAnno, genes,
                                 goColumns=c("GO BP", "GO CC", "GO MF")){
  if (setequal(genes, rownames(seqAnno))){
    return(seqAnno)
  }
  geneAnno = data.frame(row.names=na.omit(unique(genes)))
  mergeGo = function(x){
    ezCollapse(strsplit(x, "; "), na.rm=TRUE, empty.rm=TRUE, uniqueOnly=TRUE)
  }
  for (nm in intersect(goColumns, colnames(seqAnno))){
    geneAnno[nm] = ""
    x = tapply(seqAnno[[nm]], genes, mergeGo)
    geneAnno[names(x), nm] = x
  }
  return(geneAnno)
}

### -----------------------------------------------------------------
### aggregate feature annotation at isoform level into gene level
### 
aggregateFeatAnno <- function(featAnno){
  ## featAnno is the content of *_annotation.txt
  ## it's expected to contain the columns: transcript_id, gene_id, gene_name, 
  ## type, strand, seqid, start, end, biotypes, description, gc, featWidth, GO BP,
  ## GO MF, GO CC
  features <- c("gene_id", "transcript_id", "gene_name", "type", "strand", 
                "seqid", "start", "end", "biotypes", "description", "gc", 
                "featWidth", "GO BP", "GO MF", "GO CC", 
                ## below is for compatibility with old _annotation.txt
                "gene_source", "transcript_name", "hgnc_symbol", "orignal type",
                "uniprot", "Short_description", "Curator_summary", "GO",
                "Blast.Score.Eval", "EnsemblGeneID_HS", "Gene.Symbol",
                "Gene.Description")
  goColumns=c("GO BP", "GO MF", "GO CC")
  if(!all(colnames(featAnno) %in% features)){
    stop("`featAnno` can only have the columns: ", ezCollapse(features))
  }
  
  features <- intersect(features, colnames(featAnno))
  
  require(data.table)
  featAnno <- as.data.table(featAnno)
  ## Aggregate the character columns
  featAnnoGene <- featAnno[ ,
                            lapply(.SD, ezCollapse, empty.rm=TRUE, 
                                   uniqueOnly=TRUE, na.rm=TRUE),
                            by=.(gene_id), 
                            .SDcols = setdiff(features, c("gene_id", "start", 
                                                          "end", "gc", "featWidth",
                                                          goColumns))]
  ## Aggregate the numeric columns
  if (all(c("start", "end", "gc", "featWidth") %in% colnames(featAnno))){
    featAnnoGeneNumeric <- featAnno[ , .(start=min(start),
                                         end=max(end),
                                         gc=signif(mean(gc), digits = 4),
                                         featWidth=signif(mean(featWidth), 
                                                          digits = 4)),
                                     by=.(gene_id)
                                     ]
    featAnnoGene <- merge(featAnnoGene, featAnnoGeneNumeric)
  }

  ## Aggregate the GO columns which reuqire more processing
  mergeGo = function(x){
    ezCollapse(strsplit(x, "; "), na.rm=TRUE, empty.rm=TRUE, uniqueOnly=TRUE)
  }
  if(all(goColumns %in% colnames(featAnno))){
    featAnnoGeneGO <- featAnno[ , lapply(.SD, mergeGo),
                                by=.(gene_id),
                                .SDcols=goColumns]
    featAnnoGene <- merge(featAnnoGene, featAnnoGeneGO)
  }else{
    ## Some annotation has no GO terms.
    featAnnoGene$"GO BP" = ""
    featAnnoGene$"GO MF" = ""
    featAnnoGene$"GO CC" = ""
  }
  ## TODO: in the future, maybe we want to return featAnnoGene as data.table
  featAnnoGene <- as.data.frame(featAnnoGene)
  rownames(featAnnoGene) <- featAnnoGene$gene_id
  featAnnoGene$strand[!featAnnoGene$strand %in% c("+", "-")] <- "*"
  return(featAnnoGene)
}

##' @title Gets the isoform-to-gene mapping
##' @description The annotation in the reference folders is at the isoform level but for gene-level analyses we need to map the annotation to the gene level.
##' Especially all enrichment analyses should be run at the gene-level
##' @param param the global parameter object. We need the entries \code{geneColumn} or \code{geneColumnSet} to figure out 
##' which columns of the annotation data frame 
##' should be used for the mapping
##' @param seqAnnoDF a data.frame containing the sequence annotation.
##' @template roxygen-template
##' @return Returns a named character vector containing the gene mapping.
##' @examples
##' param = ezParam()
##' annoFile = system.file("extdata/genes_annotation.txt", package="ezRun", mustWork=TRUE)
##' param$ezRef["refAnnotationFile"] = annoFile
##' param$ezRef["refFeatureFile"] = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' param$ezRef["refFastaFile"] = system.file("extdata/genome.fa", package="ezRun", mustWork=TRUE)
##' seqAnno = writeAnnotationFromGtf(param)
##' seqAnno = ezFeatureAnnotation(param, rownames(seqAnno), dataFeatureType="gene")
##' gm = getGeneMapping(param, seqAnno)
##' hasGeneMapping(param, seqAnno)
getGeneMapping = function(param, seqAnnoDF){
  
  if (is.null(seqAnnoDF)){
    return(NULL)
  }
  
  if (ezIsSpecified(param[["geneColumn"]])){
    geneCol = param[["geneColumn"]]
  } else {
    geneCol = intersect(param$geneColumnSet, colnames(seqAnnoDF))[1]  	
  }
  if (is.null(geneCol) || length(geneCol) == 0 || is.na(geneCol)){
    message("getGeneMapping -- no geneCol provided: ", geneCol)
    x = NULL
  } else {
    message("getGeneMapping -- using: ", geneCol)
    x = seqAnnoDF[[geneCol]]
    x[ x == ""] = NA
    names(x) = rownames(seqAnnoDF)
  }
  return(x)
}

##' @describeIn getGeneMapping Returns a logical indicating if the set parameters have gene mapping or not.
hasGeneMapping = function(param, seqAnnoDF){
  any(param$geneColumnSet %in% colnames(seqAnnoDF) | param$featureLevel == "gene")
}

### Get the black listed genomic regions for functional genomics analysis
getBlacklistedRegions <- function(refBuild=c("hg38", "hg19", "mm10", "mm9",
                                             "ce10", "dm3"),
                                  format=c("UCSC", "Ensembl")){
  require(rtracklayer)
  require(GenomeInfoDb)
  refBuild <- match.arg(refBuild)
  format <- match.arg(format)
  
  url <- switch(refBuild,
                hg38="https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz",
                hg19="https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz",
                mm10="https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz",
                mm9="http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed.gz",
                ce10="http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz",
                dm3="http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz")
  bedGR <- import.bed(url)
  
  if(format == "Ensembl"){
    seqlevelsStyle(bedGR) <- "Ensembl"
  }
  return(bedGR)
}
 
### -----------------------------------------------------------------
### Functions for control sequences
###
makeExtraControlSeqGR <- function(ids=NULL){
  controlSeqs <- getControlSeqs(ids)
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
  return(extraGR)
}