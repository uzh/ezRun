###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezLoadFeatures = function(param=NULL, featureFile=param$ezRef["refFeatureFile"], 
                          types=NULL){
  
  if (!ezIsSpecified(featureFile)){
    ezWrite("no features specified in param")
    return(NULL)
  }
  
  gff = ezReadGff(featureFile)
  
  gff = gff[!ezGrepl(c("CDS", "codon", "UTR", "protein"), gff$type, 
                     ignore.case=TRUE), ]
  if (getSuffix(featureFile) == "gtf" ){
    gff$transcript_id = ezGffAttributeField(gff$attributes, 
                                             field="transcript_id", 
                                             attrsep="; *", valuesep=" ")
    stopifnot(!is.na(gff$transcript_id[gff$type == "exon"]))
    gff$gene_id = 	ezGffAttributeField(gff$attributes, field="gene_id", attrsep=";", valuesep=" ")
    use = is.na(gff$gene_id) & gff$type == "exon"
    gff$gene_id[use] = gff$transcript_id[use]
    stopifnot(!is.na(gff$gene_id[gff$type == "exon"]))
    gff$gene_name =   ezGffAttributeField(gff$attributes, field="gene_name", 
                                          attrsep="; *", valuesep=" ")
    gff$Parent = gff$transcript_id
    gff$ID = gff$transcript_id
    gff$Name = gff$gene_id
    isChild = gff$type %in% c("exon", "intron")
    # 		gff$ID[isChild] = NA
    # 		gff$Name[isChild] = NA
    gff$Parent[!isChild] = NA
    gff$tss_id = ezGffAttributeField(gff$attributes, field="tss_id", 
                                     attrsep="; *", valuesep=" ")
  }
  if (getSuffix(featureFile) == "gff"){
    gff$Parent = ezGffAttributeField(gff$attributes, field="Parent")
    gff$ID = ezGffAttributeField(gff$attributes, field="ID")
    gff$Name = ezGffAttributeField(gff$attributes, field="Name")
  }
  if (ezIsSpecified(param$addPromoters) && param$addPromoters){
    gff = addPromotersToGff(gff, param$addPromoters)
  }
  if (!is.null(types)){
    gff = gff[gff$type %in% types, ]
  }
  return(gff)
}

ezGffAttributeField = function (x, field, attrsep = ";", valuesep="=") {
  require(stringr)
  ## gene_id "ENSG00000278267"; gene_version "1"; gene_name "MIR6859-1"; gene_source "mirbase"; gene_biotype "miRNA"; transcript_id "ENST00000619216"; transcript_version "1"; transcript_name "MIR6859-1-201"; transcript_source "mirbase"; transcript_biotype "miRNA"; tag "basic"; transcript_support_level "NA";
  #x <- str_extract(x, paste(field, paste0("\"{0,1}[^\"]+\"{0,1}", attrsep),
  #                          sep=valuesep)
  #                 )
  x <- str_extract(paste0(x, attrsep), paste(field, paste0('"{0,1}[^"',attrsep, ']+"{0,1}', attrsep),
                            sep=valuesep)
                   )
  
  x <- sub(paste(field, "\"{0,1}", sep=valuesep), "", x)
  x <- sub(paste0("\"{0,1}", attrsep, "$"), "", x)
  # This implementation is 2 times faster than the old implementation below.
  # Old: 30.419 seconds; New: 6.353 seconds
  
#  x[ !ezGrepl(paste0(field, valuesep), x)] = NA
#  x = sub(paste0(".*", field, valuesep), "", x)
#  x = sub(paste0(attrsep, ".*"), "", x)
#  x = sub("^\"", "", sub("\"$", "", x))
  return(x)
}

## TODO: we should use promoters() to replace the following code.
## The boundaries of chromosoems should always be checked.
addPromotersToGff = function(gff, promWidth){
  isTranscript = gff$type %in% c("transcript", "mRNA", "gene")
  if (any(isTranscript)){
    geneGff = gff[isTranscript, ]
  } else {
    geneGff = groupGff(gff)
  }
  promGff = geneGff
  promGff$type = "promoter"
  isPos = geneGff$strand == "+"
  promGff$start[isPos] = pmax(geneGff$start[isPos] - promWidth, 1)
  promGff$end[isPos] = pmax(geneGff$start[isPos] - 1, 1)
  promGff$start[!isPos] = geneGff$end[!isPos] +1
  promGff$end[!isPos] = geneGff$end[!isPos] + promWidth
  for (nm in intersect(colnames(geneGff), c("ID", "Name"))){
    promGff[[nm]] = paste0(geneGff[[nm]], "_", promWidth, "B")
  }
  result = rbind(gff, promGff)
  return(result)
}

##' @title Builds an annotation attributes field
##' @description Builds a compatible attributes field for either the gtf or gff annotation format.
##' @param x a subset of an annotation data.frame to combine the attributes from.
##' @param format the annotation format to build \code{x} with. Either \code{"gtf"} or \code{"gff"}.
##' @template roxygen-template
##' @return Returns a gtf or gff attribute field.
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' attrNew = ezBuildAttributeField(gtf[ , c("gene_id", "gene_name", "transcript_id", "tss_id")],"gtf")
##' gtf$attributes2 = attrNew
ezBuildAttributeField = function(x, format="gtf"){
  if (format == "gtf") {
    for (nm in colnames(x)){
      x[[nm]] = paste0(nm, " \"", x[[nm]], "\";")
    }
    return(apply(x, 1, paste, collapse=" "))
  }
  else if (format == "gff") {
    for (nm in colnames(x)){
      x[[nm]] = paste0(nm, "=\"", x[[nm]], "\"")
    }
    return(apply(x, 1, paste, collapse=";"))
  }
  else stop("Format not supported. Use either 'gtf' or 'gff'.")
}

#### consider also
## library(genomeIntervals)
##  gff = readGff3(gffFile)
##' @title Reads an annotation table from a file.
##' @description Reads an annotation table from a gtf or gff file into a data.frame.
##' @param gffFile the gtf or gff file to read the information from.
##' @param nrows an integer specifying the maximum number of rows to read in.
##' @param x an annotation matrix or data.frame.
##' @param file a character representing the connection to write the table to.
##' @template roxygen-template
##' @return Returns a data.frame containing the annotation information.
##' @examples
##' gtf = ezReadGff(system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' ezWriteGff(gtf,"newgtf")
ezReadGff = function(gffFile, nrows=-1){
  require(data.table)
  gff <- fread(gffFile, sep="\t", header=FALSE, data.table=FALSE,
               quote="", col.names=c("seqid", "source", "type", "start", "end",
                                     "score", "strand", "phase", "attributes"),
               colClasses=c("character", "character", "character", "integer", 
                            "integer", "character", "character", "character",
                            "character"),
               nrows=nrows)
  ## fread is faster than read_tsv, read.table and import.
  ## some benchmarks for human gtf on Mac i7-3615QM
  ### read_tsv: 27.662s
  ### ezRead.table: 127.461s
  ### import: 33.419s
  
  # fread is much faster than read.table.
  # human genes.gtf, read.table: 97.124 seconds. fread: 5.976 seconds
  # gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
  #                 header=FALSE, comment.char="#", nrows = nrows,
  #                 colClasses=c("character", "character", "character", "integer", "integer",
  #                              "character", "character", "character", "character"))
  # colnames(gff) = c("seqid", "source", "type", "start", "end",
  #                  "score", "strand", "phase", "attributes")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

##' @describeIn ezReadGff Writes a gtf or gff annotation table passed by with \code{x} to the connection \code{file}, keeping only the original columns.
ezWriteGff = function(x, file){
  
  gff = x[ ,c("seqid", "source", "type", "start", "end",
              "score", "strand", "phase", "attributes")]
  ezWrite.table(gff, file=file, row.names=FALSE, col.names=FALSE)
}

## usage: gff = ezReadGff("gffFile.gff")
## gff$Parent = ezGffAttributeField(gff$attributes, "Parent")
## gffTRimmed = gffTrimTranscripts(gff[gff$type == "exon", ])
##' @title Trims transcripts in annotation data.frames
##' @description Trims transcripts in annotation data.frames of gtf or gff format.
##' @param gff an annotation data.frame in gtf or gff format.
##' @param gffSingle the annotation information for a single transcript.
##' @param maxLength an integer specifying the maximum length of transcripts after trimming.
##' @param start a logical indicating the start direction of each transcript to be trimmed.
##' @template roxygen-template
##' @return Returns
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' gtfTrimmed = gffTrimTranscripts(gtf[gtf$type == "exon", ])
gffTrimTranscripts = function(gff, maxLength=100, start=TRUE){
  
  gff$width = gff$end -  gff$start + 1
  gffList = split(gff, gff$Parent)
  gffListTrimmed = ezMclapply(gffList, function(x){gffTrimSingleTranscript(x, maxLength=maxLength, start=start)}, mc.cores=min(ezThreads(), 8))
  result = do.call("rbind", gffListTrimmed)
  return(result)
}

##' @describeIn gffTrimTranscripts Trims a single transcript according to \code{maxLength} and in reversed order depending on the strand (+ / -).
gffTrimSingleTranscript = function(gffSingle, maxLength=100, start=TRUE){
  
  reverseOrder = xor(gffSingle$strand[1] == "+", start)
  stopifnot(length(gffSingle$start) > 0 && (all(diff(gffSingle$start)> 0) || all(diff(gffSingle$start) < 0)))
  ## and place the order according the "start"
  gffSingle = gffSingle[order(gffSingle$start, decreasing=reverseOrder), , drop=FALSE]
  cs = cumsum(gffSingle$width)
  trimThis = which(cs > maxLength)[1]
  if (is.na(trimThis)){
    gffTrim = gffSingle
  } else {
    gffTrim = gffSingle[1:trimThis, ]
    if (trimThis > 1){
      maxLength = maxLength - cs[trimThis-1]
    }
    if(maxLength == 0){
      gffTrim = gffTrim[1:(trimThis-1),]
    }else{
      if(reverseOrder){
        gffTrim[trimThis, "start"] = gffTrim[trimThis, "end"] - as.integer(maxLength - 1)
      } else {
        gffTrim[trimThis, "end"] = gffTrim[trimThis, "start"] + as.integer(maxLength - 1)
      }
    }
  }
  gffTrim$width = NULL
  if(reverseOrder){
    gffTrim = gffTrim[order(gffTrim$start, decreasing=FALSE), , drop=FALSE]
  }
  return(gffTrim)
  
}

##' @title Gets GRanges from annotation
##' @description Gets GRanges from annotation of the format gtf or gff.
##' @param gff the annotation data.frame to get the GRanges from.
##' @template roxygen-template
##' @return Returns a GRanges object
##' @seealso \code{\link[GenomicRanges]{GRanges}}
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' gRanges = gffToRanges(gtf)
##' gRanges
gffToRanges = function(gff){
  require("GenomicRanges", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  gff$strand[gff$strand == "."] = "*"
  if (is.null(gff$ID)){
    gff$ID = ezGffAttributeField(gff$attributes, field="ID")
  } 
  if (any(is.na(gff$ID))){
    id = paste("gffid", 1:nrow(gff), sep="_")
  }
  return(GRanges(seqnames=gff$seqid,
                 ranges=IRanges(start=gff$start, end=gff$end, names=gff$ID),
                 strand=gff$strand))
}

## get a gtf file with the exon coordinates
# gtf = ezLoadFeatures(featureFile="/srv/GT/reference/Homo_sapiens/Ensembl/GRCh37/Annotation/Version-2014-03-29/Genes/genes.gtf", types="exon")
# transcriptGtf = groupGff(gtf, grouping=gtf$transcript_id, type="transcript")
# transcriptGtf$attributes = ezBuildAttributeField(transcriptGtf[ , c("transcript_id", "gene_id", "gene_name")],"gtf")
# ezWriteGff(transcriptGtf, "transcript.gtf")
##' @title Groups annotation
##' @description Groups annotation according to the \code{grouping} factor. Duplicated values are dropped and the function checks for transspliced elements.
##' @param gtf the annotation data.frame to group.
##' @param grouping a vector of the annotation data.frame to group it by.
##' @param type a character representing the new \code{type} name.
##' @param skipTransSpliced a logical indicating whether to skip transspliced elements if they exist.
##' @template roxygen-template
##' @return Returns a grouped annotation data.frame.
##' @seealso \code{\link[GenomicRanges]{GRanges}}
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' groupedGtf = groupGff(gtf)
##' gffGroupToRanges(gtf, grouping = gtf$Parent)
groupGff = function(gtf, grouping=gtf$Parent, type="transcript", skipTransSpliced=FALSE){
  
  isTransSpliced = tapply(paste(gtf$seqid, gtf$strand), grouping, function(x){length(unique(x)) > 1})
  if (any(isTransSpliced)){
    tsNames = names(isTransSpliced)[isTransSpliced]
    if (skipTransSpliced){
      use = !grouping %in% tsNames
      gtf = gtf[ use, ]
      grouping = grouping[use]
    } else {
      stop("transspliced elements exist in GTF: ", paste(tsNames, collapse=" "))
    }
  }
  stopifnot(tapply(paste(gtf$seqid, gtf$strand), grouping, function(x){length(unique(x))}) == 1)
  start = tapply(gtf$start, grouping, function(x){min(x)})
  end = tapply(gtf$end, grouping, function(x){max(x)})
  isDup = duplicated(grouping)
  gtfGrouped = gtf[!isDup, ]
  gtfGrouped$ID = grouping[!isDup]
  gtfGrouped$start = start[gtfGrouped$ID]
  gtfGrouped$end = end[gtfGrouped$ID]
  gtfGrouped$type = type
  if (!is.null(gtfGrouped$Parent)){
    gtfGrouped$Parent = NA
  }
  return(gtfGrouped)
}

##' @describeIn groupGff Groups annotation and returns a GRanges object from the grouped annotation.
gffGroupToRanges = function(gtf, grouping, skipTransSpliced=FALSE){
  gtfGrouped = groupGff(gtf, grouping=grouping, skipTransSpliced=skipTransSpliced)
  return(gffToRanges(gtfGrouped))
}

##' @title Gets the transcript annotation
##' @description Gets the transcript annotation from a gtf annotation.
##' @param gtf an annotation data.frame of the gtf format.
##' @param id a column of the gtf data.frame that specifies a factor to group transcripts by.
##' @template types-template
##' @param attributes additional annotation attributes to pass in.
##' @template roxygen-template
##' @return Returns a transposed transcript annotation data.frame.
##' @examples
##' param = ezParam()
##' gtfFile = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' gtf = ezLoadFeatures(param, gtfFile)
##' ta = transcriptAnnoFromGtf(gtf)
transcriptAnnoFromGtf = function(gtf, id=gtf$transcript_id, types="exon", attributes=NULL){
  stopifnot(length(id) == nrow(gtf))
  use = gtf$type %in% types
  gtf = gtf[use, ]
  id = id[use]
  seqAnno = data.frame(row.names=unique(id))
  seqAnno$gene_id = tapply(gtf$gene_id, id, ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)] 
  seqAnno$seqid = tapply(gtf$seqid, id, ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)]
  seqAnno$strand = tapply(gtf$strand, id, ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)]
  seqAnno$start = tapply(gtf$start, id, ezCollapse)[rownames(seqAnno)]        
  seqAnno$end = tapply(gtf$end, id, ezCollapse)[rownames(seqAnno)]        
  seqAnno$width = tapply(gtf$end - gtf$start + 1, id, sum)[rownames(seqAnno)]
  if (!is.null(attributes)){
    for (attr in attributes){
      aValues = ezGffAttributeField(x=gtf$attributes, field=attr, attrsep="; *", valuesep=" ")
      seqAnno[[attr]] = tapply(aValues, id, ezCollapse, na.rm=TRUE, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)]
    }
  }
  return(seqAnno)
}

##' @title Adds transcripts from exons to the annotation
##' @description Adds transcripts from exons to the annotation using a gtf or gff data.frame.
##' @param gff the annotation data.frame to get transcripts from.
##' @template roxygen-template
##' @return Returns an annotation data.frame with transcripts added to it.
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' gtf2 = addTranscriptsToGffExons(gtf)
addTranscriptsToGffExons = function(gff){
  
  gffExon = gff[gff$type == "exon", ]
  if (is.null(gffExon$transcript_id)){
    gffExon$transcript_id = ezGffAttributeField(gffExon$attributes, field="transcript_id", attrsep="; *", valuesep=" ")
  }
  gffExon$id = paste(gffExon$transcript_id, gffExon$seqid, gffExon$strand)
  gffTranscript = gffExon[!duplicated(gffExon$id), ]
  trStart = by(gffExon$start, gffExon$id, min)
  trEnd = by(gffExon$end, gffExon$id, max)
  gffTranscript$start = trStart[gffTranscript$id]
  gffTranscript$end = trEnd[gffTranscript$id]
  gffTranscript$type = "transcript"
  gffTranscript$transcript_id = NULL
  gffTranscript$id = NULL
  gffExon$transcript_id = NULL
  gffExon$id = NULL
  gffAll = rbind(gffTranscript, gffExon, gff[gff$type != "exon", ])
  return(gffAll)
}

##' @title Gets the exon numbers
##' @description Gets the exon numbers of an annotation data.frame of the gtf or gff format.
##' @param gtf the annotation data.frame to get transcripts from.
##' @template roxygen-template
##' @return Returns a vector containing the exon numbers.
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' en = getExonNumber(gtf)
getExonNumber = function(gtf){
  ## check if the same transcript id occurs at multiple loci
  locusId = paste(gtf$transcript_id, gtf$seqid, gtf$strand)
  stopifnot(tapply(locusId, gtf$transcript_id, function(x){length(unique(x))}) == 1) ## every transcript must occur only in one locus!
  exonNumbers = rep(NA, nrow(gtf))
  use = gtf$type == "exon"
  exonNumbers[use] = unsplit(tapply(gtf$start[use], gtf$transcript_id[use], order), gtf$transcript_id[use])
  return(exonNumbers)
}

##' @title Gets the transcripts database from the reference
##' @description Gets the transcripts database from the reference.
##' @param reference a reference to read the information from.
##' @param dataSource a character specifying the source of the data.
##' @template roxygen-template
##' @return Returns an object of the class TxDb
##' @seealso \code{\link[GenomicFeatures]{makeTxDbFromGFF}}
##' @examples 
##' refBuild = "Mus_musculus/UCSC/mm10/Annotation/Version-2012-05-23"
##' param = ezParam(list(refBuild=refBuild))
##' tdb = ezTranscriptDbFromRef(param$ezRef)
ezTranscriptDbFromRef = function(reference, dataSource="FGCZ"){
  require("GenomicFeatures", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  genomeFastaIndexFile = paste0(reference@refFastaFile, ".fai")
  fai = ezRead.table(genomeFastaIndexFile, header=FALSE)
  chrominfo = data.frame(chrom=rownames(fai), length=fai$V2, is_circular=FALSE)
  makeTxDbFromGFF(reference@refFeatureFile, dataSource=dataSource, chrominfo=chrominfo)
}

##' @title Gets gene names from annotation
##' @description Gets gene names from an annotation data.frame of gtf or gff format.
##' @param gtf the annotation data.frame to get the gene data from.
##' @param tr2Gene a character vector containing the transcript ID's obtained from the gene ID's.
##' @template roxygen-template
##' @return Returns a character vector with the gene names.
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' tr2Gene = getTranscript2Gene(gtf)
##' ipg = countIsoformsPerGene(tr2Gene)
getTranscript2Gene = function(gtf){
  stopifnot(!is.null(gtf$transcript_id))
  stopifnot(!is.null(gtf$gene_id))
  use = !duplicated(gtf$transcript_id)
  tr2Gene = gtf$gene_id[use]
  names(tr2Gene) = gtf$transcript_id[use]
  return(tr2Gene)
}

##' @describeIn getTranscript2Gene Returns a contingency table of the transcript names counting isoforms per gene.
countIsoformsPerGene = function(tr2Gene){
  return(table(tr2Gene))
}

##' @title Writes pre-spliced annotation
##' @description Writes pre-spliced annotation directly from the parameters into a new .gtf file.
##' @param param contains the feature file to load the annotation from.
##' @param featureFile the file to load the features in the \code{ezLoadFeatures()} call from.
##' @template roxygen-template
##' @seealso \code{\link{ezLoadFeatures}}
##' @examples
##' param = ezParam()
##' writePresplicedGtf(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
writePresplicedGtf <- function (param, featureFile=param$ezRef["refFeatureFile"]) {
  gtfFile = featureFile
  gtf = ezLoadFeatures(param=param, featureFile=featureFile)
  gtf = gtf[gtf$type == "exon", ]
  preSplicedId = paste(gtf$gene_id, gtf$tss_id, sep="_")
  gtfTr = groupGff(gtf, grouping=preSplicedId, type="exon")
  gtfTr$transcript_id = paste0(gtfTr$preSplicedId, "_pre")
  gtfBoth = rbind(gtf, gtfTr)
  attrNew = ezBuildAttributeField(gtfBoth[ , c("gene_id", "gene_name", "transcript_id", "tss_id")],"gtf")
  gtfBoth$attributes = attrNew
  gtfBothFile = sub(".gtf", "WithPrespliced.gtf", gtfFile)
  stopifnot(gtfBothFile != gtfFile)
  ezWriteGff(gtfBoth, file=gtfBothFile)
}


getTranscriptSequences = function(param=NULL, genomeFn=NULL, featureFn=NULL){
  require(GenomicFeatures)
  require(Rsamtools)
  if(!is.null(param)){
    genomeFn <- param$ezRef["refFastaFile"]
    featureFn <- param$ezRef["refFeatureFile"]
  }
  txdb = makeTxDbFromGFF(featureFn,
                         dataSource="FGCZ", taxonomyId = "10090")
  # organism=organism, chrominfo=NULL)
  exonRgList = exonsBy(txdb, by="tx", use.names=TRUE)
  genomeFasta = FaFile(genomeFn)
  trSeqs = extractTranscriptSeqs(genomeFasta, exonRgList)
  return(trSeqs)
}

##' @title Gets UTR sequences
##' @description Gets sequences from untranslated regions.
##' @param gtf the annotation data.frame to get the gene data from.
##' @param chromSeqs the chromosome sequences to get the Views from.
##' @template roxygen-template
##' @return Returns a list with two objects of the class DNAStringSet from the packages Biostrings that represent the 3 and 5 primer.
##' @seealso \code{\link[Biostrings]{DNAStringSet}}
##' @seealso \code{\link[IRanges]{Views}}
##' @seealso \code{\link[Biostrings]{reverseComplement}}
# param = ezParam()
# gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
# param$ezRef@@refFeatureFile = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
# fp = "/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa"
# param$ezRef@@refFastaFile = fp
# trSeqs = getTranscriptSequences(param)
# ezUtrSequences(gtf)
ezUtrSequences = function(gtf, chromSeqs){
  
  if (is.null(gtf$transcript_id)){
    gtf$transcript_id = ezGffAttributeField(gtf$attributes, 
                                            field = "transcript_id", attrsep = "; *", valuesep = " ")
  }
  gtf = gtf[order(gtf$start), ]
  utr = gtf[gtf$type == "UTR", ] ## UTRs may contain introns!!!
  cds = gtf[gtf$type == "CDS", ]
  codSeqid = tapply(cds$seqid, cds$transcript_id, unique)
  codStrand = tapply(cds$strand, cds$transcript_id, unique)
  codLeft = tapply(cds$start, cds$transcript_id, min)
  codRight = tapply(cds$end, cds$transcript_id, max)
  codMid = (codLeft + codRight)/2
  isLeft = utr$start < codMid[utr$transcript_id]
  isThreePrime = xor( isLeft, utr$strand == "+")
  
  threePrimeUtr = DNAStringSet()
  fivePrimeUtr = DNAStringSet()
  for (nm in names(chromSeqs)){
    message(nm)
    utr3 = utr[utr$seqid == nm & isThreePrime, ,drop=FALSE]
    utr5 = utr[utr$seqid == nm & !isThreePrime, ,drop=FALSE]
    if (nrow(utr3) > 0){
      seqs = as.character(Views(chromSeqs[[nm]], 
                                start=utr3$start,
                                end=utr3$end, names=utr3$transcript_id))
      newUtrChar = tapply(seqs, names(seqs), paste, collapse="")
      newUtr = DNAStringSet(newUtrChar)
      names(newUtr) = names(newUtrChar)
      threePrimeUtr = c(threePrimeUtr, newUtr)
    }
    if (nrow(utr5) > 0){
      seqs = as.character(Views(chromSeqs[[nm]], 
                                start=utr5$start,
                                end=utr5$end, names=utr5$transcript_id))
      newUtrChar = tapply(seqs, names(seqs), paste, collapse="")
      newUtr = DNAStringSet(newUtrChar)
      names(newUtr) = names(newUtrChar)
      fivePrimeUtr = c(fivePrimeUtr, newUtr)
    }
  }
  revTranscripts = unique(utr$transcript_id[utr$strand == "-"])
  xn = intersect(revTranscripts, names(threePrimeUtr))
  threePrimeUtr[xn] = reverseComplement(threePrimeUtr[xn])
  xn = intersect(revTranscripts, names(fivePrimeUtr))
  fivePrimeUtr[xn] = reverseComplement(fivePrimeUtr[xn])
  return(list(threePrime=threePrimeUtr, fivePrime=fivePrimeUtr))
}

##' @title Gets GC proportions and gene widths
##' @description Gets GC proportions and gene widths from annotation (gtf or gff) and sequence (fasta) information.
##' @param param the parameters to load the annotation and sequence files from.
##' @template roxygen-template
##' @return Returns a list containing named vectors for the GC proportion and the gene width. 
##' @seealso \code{\link{ezLoadFeatures}}
##' @seealso \code{\link[Biostrings]{readDNAStringSet}}
##' @seealso \code{\link[IRanges]{Views}}
##' @seealso \code{\link[Biostrings]{letterFrequency}}
##' @examples
##' param = ezParam()
##' param$ezRef@@refFeatureFile = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' fp = "/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa"
##' param$ezRef@@refFastaFile = fp
##' gw = getTranscriptGcAndWidth(param)
getTranscriptGcAndWidth = function(param=NULL, genomeFn=NULL, featureFn=NULL){
  require(data.table)
  require(Biostrings)
  require(GenomicRanges)
  require(stringr)
  require(Rsamtools)
  
  if(!is.null(param)){
    genomeFn <- param$ezRef["refFastaFile"]
    featureFn <- param$ezRef["refFeatureFile"]
  }
  genomeFasta = genomeFn
  if (!file.exists(genomeFasta)){
    return(NULL)
  }
  #genomeSeq = readDNAStringSet(genomeFasta)
  #names(genomeSeq) = sub(" .*", "", names(genomeSeq))
  gff <- fread(featureFn, header=FALSE, sep="\t")
  gff <- gff[V3=="exon", .(chr=V1, start=V4, end=V5, strand=V7, attribute=V9)]
  exonsByTx <- GRanges(seqnames=gff$chr,
                       ranges=IRanges(start=gff$start,
                                      end=gff$end),
                       strand=gff$strand)
  transcripts <- ezGffAttributeField(gff$attribute,
                                     field="transcript_id", 
                                     attrsep="; *", valuesep=" ")
  #exonsByTx <- genomeSeq[exonsByTx] ## This can results in error of XString
  #                                     object of length 2^31 or more
  exonsByTx <- getSeq(FaFile(genomeFasta), exonsByTx) 
  # FaFile can deal with spaces witin header name.
  # ">chr1 test1" can be subset by chr1 GRanges.
  gcCount <- letterFrequency(exonsByTx, letters="GC", as.prob = FALSE)[ ,"G|C"]
  txWidth <- width(exonsByTx)
  data <- data.table(tx=transcripts, txWidth=txWidth, gcCount=gcCount)
  data <- data[ , lapply(.SD, sum), by=.(tx)]
  ans <- list(gc=signif(setNames(data[ , gcCount/txWidth], data[ ,tx]),
                        digits=4),
              width=setNames(data[ , txWidth], data[ ,tx]))
  return(ans)
}

##' @title Gets the ensembl types
##' @description Gets the ensembl types of an annotation data.frame either directly from the source, the gene_biotype or the gene_type.
##' @param gff an annotation data.frame in gtf or gff format.
##' @template roxygen-template
##' @return Returns a character vector containing the types of the elements or NULL if they are not provided in \code{gff}.
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' et = getEnsemblTypes(gtf)
getEnsemblTypes = function(gff){
  if ("protein_coding" %in% gff$source){
    return(gff$source)
  }
  types = ezGffAttributeField(x=gff$attributes, field="gene_biotype", attrsep="; *", valuesep=" ")
  if ("protein_coding" %in% types){
    return(types)
  }
  types = ezGffAttributeField(x=gff$attributes, field="gene_type", attrsep="; *", valuesep=" ")
  if ("protein_coding" %in% types){
    return(types)
  }
  return(NULL)
}

.checkGtfForExons = function(){
  ## TODO 
  gtfFile = "/srv/GT/reference/Homo_sapiens/Ensembl/GRCh37/Annotation/Version-2014-03-28/Genes/genes.gtf"
  gtf = ezReadGff(gtfFile)
  gtf$Parent = ezGffAttributeField(x=gtf$attributes, field="transcript_id", attrsep="; *", valuesep=" ")
  cdsGtf = gtf[gtf$type == "CDS", ]
  exonGtf = gtf[gtf$type == "exon", ]
  cdsRanges = gffToRanges(cdsGtf)
  exonRanges = gffToRanges(exonGtf)
  cds2exon = findOverlaps(cdsRanges, exonRanges, type="within")
  cdsIdx = queryHits(cds2exon)
  exonIdx = subjectHits(cds2exon)
  table(1:length(cdsRanges) %in% cdsIdx)
  cdsParent = cdsGtf$Parent[cdsIdx]
  exonParents = exonGtf$Parent[exonIdx]
  isMatch = cdsParent == exonParents
  hasExon = tapply(isMatch, cdsParent, any)
  table(hasExon, useNA="always")
}
