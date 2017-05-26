###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets the sequence annotation
##' @description Gets the sequence annotation from annotation (gtf or gff) and sequence (fasta) information.
##' Sequence annotation files are on the isoform level. If the analysis is to be done at the gene level, the annotation for
##' the isoforms is aggregated.
##' @param param the parameters to load the annotation and sequence files from and possibly write to.
##' @param ids a character vector containing the gene ID's to return.
##' @param dataFeatureType either \code{"isoform"} or \code{"gene"}.
##' @template roxygen-template
##' @return Returns a data.frame containing information about the genes in an easily readable way.
##' @examples
##' param = ezParam()
##' param$ezRef["refFeatureFile"] = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' param$ezRef["refFastaFile"] = system.file("extdata/genome.fa", package="ezRun", mustWork=TRUE)
##' annoFile = system.file("extdata/genes_annotation.txt", package="ezRun", mustWork=TRUE)
##' param$ezRef["refAnnotationFile"] = annoFile
##' seqAnno = writeAnnotationFromGtf(param)
##' seqAnno2 = ezFeatureAnnotation(param, rownames(seqAnno), dataFeatureType="gene")
ezFeatureAnnotation = function(param, ids, dataFeatureType){
  
  stopifnot(dataFeatureType %in% c("transcript", "isoform", "gene")) ## not yet supported
  if (!file.exists(param$ezRef["refAnnotationFile"])){
    seqAnno = writeAnnotationFromGtf(param)
  } else {
    seqAnno = ezRead.table(param$ezRef["refAnnotationFile"], colClasses="character")
    if (!is.null(seqAnno$gc)){
      seqAnno$gc = signif(as.numeric(seqAnno$gc), digits=4)
    }
    if (!is.null(seqAnno$width)){
      seqAnno$width = as.numeric(seqAnno$width)
    }
  }
  if (dataFeatureType == "gene"){
    ## TODO: use the aggregated _annotation_byGene.txt file
    if (!all(ids %in% rownames(seqAnno))){
      stopifnot(ids %in% seqAnno$gene_id)
      trNames = tapply(rownames(seqAnno), seqAnno$gene_id, ezCollapse, uniqueOnly=TRUE)
      if (!is.null(seqAnno$gene_name)){
        gene_name = tapply(seqAnno$gene_name, seqAnno$gene_id, ezCollapse, uniqueOnly=TRUE)
      } else {
        gene_name = NULL
      }
      seqAnnoGene = data.frame(row.names=names(trNames), stringsAsFactors=FALSE)
      if (!is.null(seqAnno$gene_name)){
        seqAnnoGene$gene_name = tapply(seqAnno$gene_name, seqAnno$gene_id, ezCollapse, uniqueOnly=TRUE)
      }
      seqAnnoGene$transcript_id=as.character(trNames)
      for (nm in intersect(colnames(seqAnno), c("type", "Description", "description", "hgnc_symbol", "Entrez Gene ID", "gene_name", "gene_symbol"))){
        seqAnnoGene[[nm]] = tapply(seqAnno[[nm]], seqAnno$gene_id,
                                   ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE, na.rm=TRUE)[rownames(seqAnnoGene)]
      }  
      goAnno = aggregateGoAnnotation(seqAnno, seqAnno$gene_id)
      if (!is.null(seqAnno$width)){
        seqAnnoGene$width = signif(tapply(seqAnno$width, seqAnno$gene_id, mean)[rownames(seqAnnoGene)], digits=3)
      }
      if (!is.null(seqAnno$gc)){
        seqAnnoGene$gc = signif(tapply(seqAnno$gc, seqAnno$gene_id, mean)[rownames(seqAnnoGene)], digits=3)
      }
      seqAnno = cbind(seqAnnoGene, goAnno[rownames(seqAnnoGene), ], stringsAsFactors=FALSE)
    }
  }
  seqAnno = seqAnno[ids, ,drop=FALSE]
  return(seqAnno)
}

##' @describeIn ezFeatureAnnotation Gets the annotation from a .gtf file and transforms it into a tab-separated tabular .txt file.
writeAnnotationFromGtf = function(param, featureFile=param$ezRef@refFeatureFile, featAnnoFile=param$ezRef@refAnnotationFile){
  gtf = ezLoadFeatures(param, featureFile=featureFile)
  gtf = gtf[gtf$type == "exon", ]
  seqAnno = transcriptAnnoFromGtf(gtf, attributes=c("gene_name"))
  trProps = getTranscriptGcAndWidth(param)
  seqAnno$gc = trProps$gc[rownames(seqAnno)]
  ezWrite.table(seqAnno, file=featAnnoFile)
  return(seqAnno)
}

### -----------------------------------------------------------------
### make the feature annotation file <name>_annotation.txt
### for Ensembl gtf.
makeFeatAnnoEnsembl <- function(featureFile, 
                                featAnnoFile=sub(".gtf", "_annotation.txt", 
                                                 featureFile),
                                biomartFile=NULL,
                                organism="hsapiens_gene_ensembl",
                                host=NULL){
  require(plyr)
  
  require(rtracklayer)
  feature <- import(featureFile)
  transcripts <- feature[feature$type=="transcript"]
  featAnno <- ezFrame(transcript_id=transcripts$transcript_id,
                      gene_id=transcripts$gene_id,
                      gene_name=transcripts$gene_name,
                      type=transcripts$gene_biotype,
                      strand=strand(transcripts),
                      seqid=seqnames(transcripts),
                      start=start(transcripts),
                      end=end(transcripts),
                      biotypes=transcripts$gene_biotype,
                      row.names=transcripts$transcript_id
                      )
  
  ## Group the biotype into more general groups
  stopifnot(all(featAnno$type %in% listBiotypes("all")))
  isProteinCoding <- featAnno$type %in% listBiotypes("protein_coding")
  isLNC <- featAnno$type %in% listBiotypes("long_noncoding")
  isSHNC <- featAnno$type %in% listBiotypes("short_noncoding")
  isrRNA <- featAnno$type %in% listBiotypes("rRNA")
  istRNA <- featAnno$type %in% listBiotypes("tRNA")
  isPseudo <- featAnno$type %in% listBiotypes("pseudogene")
  featAnno$type[isPseudo] <- "pseudogene"
  featAnno$type[isLNC] <- "long_noncoding"
  featAnno$type[isSHNC] <- "short_noncoding"
  featAnno$type[isProteinCoding] <- "protein_coding"
  ### rRNA and tRNA have to be after noncoding
  ### since they are subset of noncoding
  featAnno$type[isrRNA] <- "rRNA"
  featAnno$type[istRNA] <- "tRNA"
  
  ## additional information from Ensembl or downloaded biomart file
  attributes <- c("ensembl_transcript_id", "description", 
                  "go_id", "namespace_1003",
                  "percentage_gene_gc_content", "transcript_length")
  names(attributes) <- c("Transcript stable ID", "Gene description",
                         "GO term accession", "GO domain",
                         "% GC content", 
                         "Transcript length (including UTRs and CDS)")
  if(!is.null(biomartFile)){
    require(readr)
    mapping <- read_tsv(biomartFile)
    if(!setequal(colnames(mapping), names(attributes))){
      stop("Make sure ", paste(names(attributes), collapse="; "), 
           "are downloaded from web biomart!")
    }
    colnames(mapping) <- attributes[colnames(mapping)]
  }else{
    require(biomaRt)
    if(is.null(host)){
      ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
    }else{
      ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host=host)
    }
    ensembl <- useDataset(organism, mart=ensembl)
    mapping <-
      getBM(attributes=attributes,
            filters=c("ensembl_transcript_id"),
            values=rownames(featAnno), mart=ensembl)
  }
  
  ### description, gc, width
  txid2description <- mapping[c("ensembl_transcript_id", "description", 
                                "percentage_gene_gc_content", 
                                "transcript_length")]
  colnames(txid2description) <- c("transcript_id", "description",
                                  "gc", "width")
  txid2description$gc <- txid2description$gc / 100
  featAnno <- join(featAnno, txid2description, by="transcript_id", 
                   match="first")
  
  ### GO
  GOMapping <- c("biological_process"="GO BP",
                 "molecular_function"="GO MF",
                 "cellular_component"="GO CC")
  for(i in 1:length(GOMapping)){
    isBP <- mapping$namespace_1003 == names(GOMapping)[i]
    txid2BP <- lapply(split(mapping$go_id[isBP], 
                            mapping$ensembl_transcript_id[isBP]), unique)
    txid2BP <- sapply(txid2BP, paste, collapse="; ")
    #### Some transcript_id may not have GO IDs. we use match.
    featAnno[[GOMapping[i]]] <- txid2BP[match(featAnno$transcript_id, 
                                              names(txid2BP))]
  }
  featAnno[is.na(featAnno)] <- ""
  featAnno$transcript_id <- NULL
  ezWrite.table(featAnno, file=featAnnoFile)
  invisible(featAnno)
}

##' @describeIn ezFeatureAnnotation Aggregates the Go annotation.
aggregateGoAnnotation = function(seqAnno, genes, goColumns=c("GO BP", "GO CC", "GO MF")){
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
