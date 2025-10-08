###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title The S4 class representing a collection of information about the reference build
##' @description Objects of this class are a collection of references derived from a list of parameters.
##' @slot refBuild a character specifying the file path to the reference build.
##' @slot refBuildName a character specifying the name of the reference build.
##' @slot refBuildDir a character specifying the directory of the reference build.
##' @slot refIndex a character specifying the location of the index that is used in the alignment.
##' @slot refFeatureFile a character specifying the file path to the annotation feature file (.gtf).
##' @slot refAnnotationFile a character specifying the file path to the annotation file (.txt). Subsetting species are automatically removed
##' @slot refFastaFile a character specifying the file path to the fasta file.
##' @slot refChromSizesFile a character specifying the file path to the file containing the chromosome sizes.
##' @slot refAnnotationVersion a character specifying the annotation version.
##' @template roxygen-template
##' @seealso \code{\link{cleanGenomeFiles}}
##' @examples
##' refBuild = "Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18"
##' genomesRoot = "./refExample"
##' param = ezParam(list(refBuild=refBuild, genomesRoot=genomesRoot))
##' gtf = system.file("extdata/genes.gtf", package="ezRun", mustWork = TRUE)
##' fp = system.file("extdata/genome.fa", package="ezRun", mustWork = TRUE)
##' buildRefDir(param$ezRef, fp, gtf)
##' buildIgvGenome(param$ezRef)
##' seqAnno = writeAnnotationFromGtf(param=param)
## featureFile=param$ezRef["refFeatureFile"], featAnnoFile=param$ezRef["refAnnotationFile"])
setClass("EzRef", slots = c(refBuild="character",
                            refBuildName="character",
                            refBuildDir="character",
                            refIndex="character",
                            refFeatureFile="character",
                            refAnnotationFile="character",
                            refFastaFile="character",
                            refChromSizesFile="character",
                            refAnnotationVersion="character",
                            refVariantsDir="character")
)

##' @describeIn EzRef Initializes the slots of EzRef. It will also try to specify some fields and if necessary get full file paths.
EzRef <- function(param){
  if(ezIsSpecified(param$genomesRoot)){
    genomesRoot <- param$genomesRoot
  }else{
    genomesRoot <- strsplit(GENOMES_ROOT, ":")[[1]]
  }
  refBuild <- param$refBuild
  refFields <- str_split(refBuild, fixed("/"))[[1]]
  if (ezIsSpecified(param$refBuildName)){
    refBuildName <- param$refBuildName
  } else {
    refBuildName <- refFields[3]
  }
  if (ezIsSpecified(param$refBuildDir)){
    refBuildDir <- param$refBuildDir
  } else {
    for (gr in genomesRoot){
      stopifnot(file.exists(gr))
      rbd <- file.path(gr, paste(refFields, collapse="/"))
      if (file.exists(rbd)){
        break
      }
    }
    refBuildDir <- file.path(gr, str_c(refFields[1:3], collapse="/"))
  }
  refVariantsDir <- file.path(refBuildDir, "Variants")
  if (length(refFields) == 5 && str_detect(refFields[5], "^Version|^Release")){
    refAnnotationVersion <- refFields[5]
  } else {
    refAnnotationVersion <- ""
  }
  if(ezIsSpecified(param$refIndex)){
    refIndex <- param$refIndex
  }else{
    refIndex <- ""
  }
  if(ezIsAbsolutePath(param$refFeatureFile)){
    refFeatureFile <- param$refFeatureFile
  }else{
    if(ezIsSpecified(refAnnotationVersion)){
      refFeatureFile <- file.path(refBuildDir, "Annotation",
                                  refAnnotationVersion, "Genes",
                                  param$refFeatureFile)
    }else{
      refFeatureFile <-  file.path(refBuildDir, "Annotation", "Genes",
                                   param$refFeatureFile)
    }
  }
  if(ezIsAbsolutePath(param$refAnnotationFile)){
    refAnnotationFile <- param$refAnnotationFile
  }else{
    annoBaseName <- str_replace(basename(refFeatureFile), "\\.gtf$", "")
    annoBaseName <- str_replace(annoBaseName, "_.*", "")
    annoBaseName <- str_c(annoBaseName, "_annotation_byTranscript.txt")
    refAnnotationFile <- file.path(dirname(refFeatureFile), annoBaseName)
  }
  if(ezIsAbsolutePath(param$refFastaFile)){
    refFastaFile <- param$refFastaFile
  }else{
    refFastaFile <- file.path(refBuildDir, param$refFastaFile)
  }
  refChromSizesFile <- str_replace(refFastaFile, "\\.fa$", "-chromsizes.txt")
  
  new("EzRef", refBuild=refBuild, refBuildName=refBuildName,
      refBuildDir=refBuildDir, refVariantsDir=refVariantsDir,
      refAnnotationVersion=refAnnotationVersion, refIndex=refIndex,
      refAnnotationFile=refAnnotationFile, refFeatureFile=refFeatureFile,
      refFastaFile=refFastaFile, refChromSizesFile=refChromSizesFile
  )
}

##' @describeIn EzRef Access of slots by name with square brackets [ ]
setMethod("[", "EzRef", function(x, i){
  slot(x, i)
})

##' @describeIn EzRef Assignment to slots by name with square brackets [ ]
setMethod("[<-", "EzRef", function(x, i, value){
  slot(x, i) = value
  x
})

getOrganism <- function(x){
  str_split(x@refBuild, "/")[[1]][1]
}

buildRefDir <- function(x, genomeFile, genesFile=NULL, keepOriginalIDs = FALSE){
  # x is EzRef object
  require(rtracklayer)
  require(Rsamtools)
  
  gtfPath <- dirname(x@refFeatureFile)
  fastaPath <- dirname(x@refFastaFile)
  dir.create(gtfPath, recursive=TRUE)
  dir.create(fastaPath, recursive=TRUE)
  if(!is.null(x@refAnnotationVersion)){
    unlink(file.path(x@refBuildDir, "Annotation", "Genes"), recursive = TRUE)
    file.symlink(file.path(x@refAnnotationVersion, "Genes"),
                 file.path(x@refBuildDir, "Annotation", "Genes"))
  }
  
  ## fasta
  genome <- readBStringSet(genomeFile)
  ### remove everything after chr id
  names(genome) <- str_replace(names(genome), " .*$", "")
  writeXStringSet(genome, x@refFastaFile)
  indexFa(x@refFastaFile)
  
  ## create the chromsizes file
  fai <- read_tsv(str_c(x@refFastaFile, ".fai"), col_names = FALSE)
  write_tsv(fai %>% dplyr::select(1:2), file = x@refChromSizesFile, col_names = FALSE)
  
  dictFile <- str_replace(x@refFastaFile, "\\.fa$", ".dict")
  if (file.exists(dictFile)) {
    file.remove(dictFile)
  }
  cmd <- paste("java -Xms1g -Xmx10g -Djava.io.tmpdir=. -jar /misc/ngseq12/packages/Tools/Picard/3.2.0/picard.jar", "CreateSequenceDictionary",
               "-R", x@refFastaFile, "-O", dictFile)
  ezSystem(cmd)
  
  
  if (!is.null(genesFile)){  
    ## 2 GTF files:
    ### features.gtf
    gtf <- import(genesFile)
    
    #### some controls over gtf
    if(is.null(gtf$gene_biotype)){
      if(is.null(gtf$gene_type)){
        message("gene_biotype is not available in gtf. Assigning protein_coding.")
        gtf$gene_biotype <- "protein_coding"
      }else{
        ## In GENCODE gtf, there is gene_type, instead of gene_biotype.
        gtf$gene_biotype <- gtf$gene_type
        gtf$gene_type <- NULL
      }
    }
    #### GENCODE.gtf: remove the version number from gene_id, transcript_id but keep additional information if available, e.g. ENST00000429181.6_PAR_Y -> ENST00000429181_PAR_Y
    if(!keepOriginalIDs){
      gtf$gene_id <- str_replace(gtf$gene_id, "\\.\\d+", "")
      gtf$transcript_id <- str_replace(gtf$transcript_id, "\\.\\d+", "")
    }
    
    if(is.null(gtf$gene_name)){
      message("gene_name is not available in gtf. Assigning gene_id.")
      gtf$gene_name <- gtf$gene_id
    }
    
    export(gtf, con=file.path(gtfPath, "features.gtf"))
    ### genes.gtf
    export(gtf[gtf$gene_biotype %in% listBiotypes("genes")],
           con=file.path(gtfPath, "genes.gtf"))
    ### transcripts.only.gtf
    export(gtf[gtf$type %in% "transcript"],
           con=file.path(gtfPath, "transcripts.only.gtf"))
  }    
}

##' @describeIn listBiotypes returns the Ensembl gene_biotypes according to more general groups.
listBiotypes <- function(select=c("genes", "protein_coding", "long_noncoding",
                                  "short_noncoding", "rRNA", "tRNA", "Mt_tRNA", 
                                  "Mt_rRNA", "pseudogene", "all")) {
  require(yaml)
  
  ## Based on http://www.ensembl.org/Help/Faq?id=468
  ## http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
  ## http://www.ensembl.org/Help/Glossary
  select <- match.arg(select)
  
  ## Load biotype mappings from YAML configuration file
  configFile <- system.file("extdata", "biotypes.yml", package="ezRun")
  if (!file.exists(configFile)) {
    stop("Biotypes configuration file not found at: ", configFile)
  }
  biotypesConfig <- read_yaml(configFile)
  
  ## Extract biotype categories from config
  proteinCoding <- biotypesConfig$protein_coding
  pseudogene <- biotypesConfig$pseudogene
  lnc <- biotypesConfig$long_noncoding
  shnc <- biotypesConfig$short_noncoding
  
  unionBiotypes <- unique(c(proteinCoding, pseudogene, lnc, shnc))
  ## genes.gtf
  genes <- unique(c(setdiff(unionBiotypes, pseudogene),
                    grep("transcribed", pseudogene, value=TRUE)
  )
  )
  
  types <- switch(select,
                  "genes"=genes,
                  "protein_coding"=proteinCoding,
                  "long_noncoding"=lnc,
                  "short_noncoding"=shnc,
                  "pseudogene"=pseudogene,
                  "rRNA"=c("rRNA"),
                  "tRNA"=c("tRNA"),
                  "Mt_tRNA"=c("Mt_tRNA"),
                  "Mt_rRNA"=c("Mt_rRNA"),
                  "all"=unionBiotypes
  )
  return(types)
}

buildIgvGenome <- function(x){
  # x is a EzRef object
  
  ## create transcript.only.gtf
  gtfFile <- x@refFeatureFile
  genomeFile <- x@refFastaFile
  genomeFileURL <- file.path(REF_HOST, str_replace(x@refBuild, "/Annotation.*$", ""),
                             str_replace(genomeFile, "^.*/Sequence", "Sequence"))
  stopifnot(file.exists(gtfFile))
  stopifnot(file.exists(genomeFile))
  trxFile <- file.path(x@refBuildDir, "transcripts.only.gtf")
  tryCatch({
    gtf <- ezLoadFeatures(featureFile=gtfFile, types="exon")
    transcriptGtf <- groupGff(gtf, grouping=gtf$transcript_id, type="transcript")
    transcriptGtf$attributes <- ezBuildAttributeField(transcriptGtf[ , c("transcript_id", "gene_id", "gene_name")])
    transcriptGtf <- transcriptGtf[order(transcriptGtf$seqid, transcriptGtf$start, transcriptGtf$end), ]
    ezWriteGff(transcriptGtf, trxFile)
  }, error=function(e){
    message("Could not load features. Copy the annotation file instead.")
    file.copy(gtfFile, trxFile)
  })
  
  ## sort and index genes.gtf
  sortedGtfFile <- file.path(dirname(gtfFile), "genes.sorted.gtf")
  cmd <- paste("igvtools", "sort", gtfFile, sortedGtfFile)
  ezSystem(cmd)
  cmd <- paste("igvtools", "index", sortedGtfFile)
  ezSystem(cmd)
  
  ## make chrom_alias.tab
  chromFile <- file.path(x@refBuildDir, "chrom_alias.tab")
  write_lines(fasta.index(genomeFile)$desc, chromFile)
  
  ## make property.txt
  propertyFile <- file.path(x@refBuildDir, "property.txt")
  id <- x@refBuildName
  name <- str_c(getOrganism(x), id, sep="_")
  properties <- c("fasta=true", "fastaDirectory=false", "ordered=true")
  properties <- c(properties, str_c("id=", id), str_c("name=", name),
                  str_c("sequenceLocation=", genomeFileURL))
  properties <- c(properties, "geneFile=transcripts.only.gtf",
                  "chrAliasFile=chrom_alias.tab")
  write_lines(properties, propertyFile)
  
  ## make zip file and clean up
  zipFile <- paste0("igv_", id, ".genome")
  filesToZip <- c(trxFile, chromFile, propertyFile)
  cd <- getwd()
  setwd(x@refBuildDir)
  zip::zip(zipFile, basename(filesToZip))
  setwd(cd)
  unlink(filesToZip)
}

### -----------------------------------------------------------------
### Fetch the control sequences from https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa
###
getControlSeqs <- function(ids=NULL){
  genomesRoot <- strsplit(GENOMES_ROOT, ":")[[1]]
  controlSeqsFn <- file.path(genomesRoot, "controlSeqs.fa")
  controlSeqsFn <- head(controlSeqsFn[file.exists(controlSeqsFn)], 1)
  controlSeqs <- readDNAStringSet(controlSeqsFn)
  names(controlSeqs) <- sub(" .*$", "", names(controlSeqs))
  if(!is.null(ids)){
    controlSeqs <- controlSeqs[ids]
  }
  return(controlSeqs)
}

makeBiomartFile <- function(gtf){
  require(httr)
  require(jsonlite)
  require(xml2)
  require(GO.db)
  
  server <- "http://rest.ensembl.org"
  
  txids <- unique(na.omit(gtf$transcript_id))
  tb0 <- tibble("Transcript stable ID"=txids)
  
  ## txid, GO ID
  goList <- list()
  for(txid in txids){
    ext <- paste0("/xrefs/id/", txid, "?external_db=GO;all_levels=0;object_type=transcript")
    r <- GET(paste0(server, ext), content_type("application/json"))
    ans_go <- fromJSON(toJSON(httr::content(r))) ## GO.db masks content
    goList[[txid]] <- unique(na.omit(unlist(ans_go$primary_id)))
  }
  tb1 <- stack(goList)
  tb1 <- as_tibble(tb1) %>% mutate(ind=as.character(ind)) %>%
    rename("Transcript stable ID"=ind,
           "GO term accession"=values)
  
  
  ## txid, gene name
  tb2 <- mcols(gtf)[ ,c("gene_name", "transcript_id")]
  tb2 <- as_tibble(tb2) %>% filter(!is.na(transcript_id)) %>%
    rename("Transcript stable ID"=transcript_id, "Gene description"=gene_name)
  
  ## GO ID, GO domain
  # uniqueGOIDs <- unique(tb1$`GO term accession`)
  # goDomain <- setNames(character(length(uniqueGOIDs)), uniqueGOIDs)
  # for(goid in uniqueGOIDs){
  #   ext <- paste0("/ontology/id/", goid, "?")
  #   r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  #   ans_go <- fromJSON(toJSON(httr::content(r)))
  #   goDomain[goid] <- ans_go$namespace
  # }
  # tb3 <- tibble(`GO term accession`=names(goDomain),
  #               `GO domain`=goDomain)
  # Information from Ensembl REST is not consistent. Some GO IDs are missing.
  tb3 <- AnnotationDbi::select(GO.db, columns=c("GOID", "ONTOLOGY"),
                               keys=unique(tb1$`GO term accession`),
                               keytype="GOID")
  tb3 <- as_tibble(tb3) %>% mutate(ONTOLOGY=c("BP"="biological_process",
                                              "MF"="molecular_function",
                                              "CC"="cellular_component")[ONTOLOGY])
  tb3 <- rename(tb3, "GO term accession"=GOID, "GO domain"=ONTOLOGY)
  
  tb <- left_join(left_join(left_join(tb0, tb1), tb2), tb3)
  tb <- dplyr:: select(tb, `Transcript stable ID`, `Gene description`,
                       `GO term accession`, `GO domain`)
  tb <- mutate(tb, `GO term accession`=ifelse(is.na(`GO domain`), NA,
                                              `GO term accession`))
  return(tb)
}

