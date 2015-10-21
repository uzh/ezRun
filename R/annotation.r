###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets the sequence annotation
##' @description Gets the sequence annotation from annotation (gtf or gff) and sequence (fasta) information.
##' @param param the parameters to load the annotation and sequence files from and possibly write to.
##' @param ids a character vector containing the gene ID's to return.
##' @param dataFeatureType either \code{"isoform"} or \code{"gene"}.
##' @template roxygen-template
##' @return Returns a data.frame containing information about the genes in an easily readable way.
##' @examples
##' param = ezParam()
##' param$ezRef@@refFeatureFile = "./inst/extdata/genes.gtf"
##' param$ezRef@@refAnnotationFile = "./inst/extdata/genes_annotation_example.txt"
##' param$ezRef@@refFastaFile = "./script/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa"
##' x = writeAnnotationFromGtf(param)
##' ezFeatureAnnotation(param, rownames(x), "gene")
ezFeatureAnnotation = function(param, ids, dataFeatureType){
  
  stopifnot(dataFeatureType %in% c("isoform", "gene")) ## not yet supported
  if (!file.exists(param$ezRef["refAnnotationFile"])){
    seqAnno = writeAnnotationFromGtf(param)
  } else {
    seqAnno = ezRead.table(param$ezRef["refAnnotationFile"], colClasses="character")
    if (!is.null(seqAnno$gc)){
      seqAnno$gc = signif(as.numeric(seqAnno$gc), digits=4)
    }
    if (!is.null(seqAnno$wdith)){
      seqAnno$width = as.numeric(seqAnno$width)
    }
  }
  if (dataFeatureType == "gene"){
    if (!all(ids %in% rownames(seqAnno))){
      stopifnot(ids %in% seqAnno$gene_id)
      trNames = tapply(rownames(seqAnno), seqAnno$gene_id, collapse, uniqueOnly=TRUE)
      if (!is.null(seqAnno$gene_name)){
        gene_name = tapply(seqAnno$gene_name, seqAnno$gene_id, collapse, uniqueOnly=TRUE)
      } else {
        gene_name = NULL
      }
      seqAnnoGene = data.frame(row.names=names(trNames), stringsAsFactors=FALSE)
      if (!is.null(seqAnno$gene_name)){
        seqAnnoGene$gene_name = tapply(seqAnno$gene_name, seqAnno$gene_id, collapse, uniqueOnly=TRUE)
      }
      seqAnnoGene$transcript_id=as.character(trNames)
      for (nm in intersect(colnames(seqAnno), c("type", "Description", "description", "hgnc_symbol", "Entrez Gene ID", "gene_name", "gene_symbol"))){
        seqAnnoGene[[nm]] = tapply(seqAnno[[nm]], seqAnno$gene_id,
                                   collapse, empty.rm=TRUE, uniqueOnly=TRUE, na.rm=TRUE)[rownames(seqAnnoGene)]
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

##' @describeIn ezFeatureAnnotation Gets the annotation from a .gtf file and transforms it into a more readable .txt file.
writeAnnotationFromGtf = function(param, featureFile=param$ezRef["refFeatureFile"], featAnnoFile=param$ezRef["refAnnotationFile"]){
  gtf = ezLoadFeatures(param, featureFile=featureFile)
  gtf = gtf[gtf$type == "exon", ]
  seqAnno = transcriptAnnoFromGtf(gtf, attributes=c("gene_name"))
  trProps = getTranscriptGcAndWidth(param)
  seqAnno$gc = trProps$gc[rownames(seqAnno)]
  ezWrite.table(seqAnno, file=featAnnoFile)
  return(seqAnno)
}

##' @describeIn ezFeatureAnnotation Aggregates the Go annotation.
aggregateGoAnnotation = function(seqAnno, genes, goColumns=c("GO BP", "GO CC", "GO MF")){
  if (setequal(genes, rownames(seqAnno))){
    return(seqAnno)
  }
  geneAnno = data.frame(row.names=na.omit(unique(genes)))
  mergeGo = function(x){
    collapse(strsplit(x, "; "), na.rm=TRUE, empty.rm=TRUE, uniqueOnly=TRUE)
  }
  for (nm in intersect(goColumns, colnames(seqAnno))){
    geneAnno[nm] = ""
    x = tapply(seqAnno[[nm]], genes, mergeGo)
    geneAnno[names(x), nm] = x
  }
  return(geneAnno)
}

##' @title Gets the gene mapping
##' @description Gets the gene mapping.
##' @param param the parameters to specify which names to use. This is taken from \code{geneColumn} or \code{geneColumnSet} if the first one is not specified.
##' @param seqAnnoDF a data.frame containing the sequence annotation.
##' @template roxygen-template
##' @return Returns a named character vector containing the gene mapping.
##' @examples
##' param = ezParam()
##' param$ezRef@@refFeatureFile = "./inst/extdata/genes.gtf"
##' param$ezRef@@refAnnotationFile = "./inst/extdata/genes_annotation_example.txt"
##' param$ezRef@@refFastaFile = "./script/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa"
##' seqAnno = writeAnnotationFromGtf(param)
##' getGeneMapping(param,seqAnno)
##' hasGeneMapping(param,seqAnno)
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
