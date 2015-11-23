###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


.dexseqFlattenGtf = function(gtf, param){
  
  gtfFileModified = sub(".gtf$", "", param$ezRef["refFeatureFile"])
  gtfFileModified = paste0(gtfFileModified, "-dexseq-", param$featureLevel, ".gtf")
  gtfFileFlattened = sub(".gtf", "-flattened.gtf", gtfFileModified)
  lockFile = paste0(gtfFileFlattened, ".lock")
  i = 0
  while(file.exists(lockFile) && i < 30){
    ### somebody else builds and we wait at most 30min
    Sys.sleep( 60)
    i = i + 1
  }
  if (file.exists(lockFile)){
    stop("reference building still in progress after 30 min")    
  }
  if (file.exists(gtfFileFlattened)){
    return(gtfFileFlattened)
  }
  ezWrite(Sys.info(), con=lockFile)
  on.exit(file.remove(lockFile), add=TRUE)
  gtf = ezLoadFeatures(param)
  stopifnot(param$featureLevel %in% c("gene", "tss"))
  if (param$featureLevel == "tss"){
    gtf$attributes = sub("gene_id", "orig_gene_id", gtf$attributes)
    gtf$attributes = sub("tss_id", "gene_id", gtf$attributes)
  }
  ## fix the gene ids
  chromStrand = paste0(gtf$seqid, strandName(gtf$strand))
  geneToChromStrand = tapply(chromStrand, gtf$gene_id, unique)
  ## append to the loci a the chromStrand
  isDup = sapply(geneToChromStrand, length) > 1
  dupGenes = names(geneToChromStrand)[isDup]
  use = gtf$gene_id %in% dupGenes
  gtf$gene_id[use] = paste(gtf$gene_id[use], chromStrand[use], sep="_")
  gtf$attributes = ezBuildAttributeField(gtf[ , c("transcript_id", "gene_id", "gene_name", "tss_id")],"gtf")
  ezWriteGff(gtf, gtfFileModified)
  cmd = paste(DEXSEQ_PREPARE, gtfFileModified, gtfFileFlattened)
  ezSystem(cmd)
  return(gtfFileFlattened)
}

.countDexseq = function(input=NA, output=NA, param=NA){
  
  param = fillWithDefaults(param) ## TODO: function doesn't exist
  checkFreeDiskSpace(param)
  options(cores=param$cores)
  if (param$paired == FALSE){
    bamFile = input$"BAM"
    paired = "no"
  } else {
    ## currently only support for full paired!!!!
    maxMem = "1000000000" ## sort with 10GB RAM
    sampleName = paste0("mySample",Sys.getpid())
    cmd = paste(SAMTOOLS, "sort", "-n", "-m", maxMem, input, sampleName)
    ezSystem(cmd)
    bamFile = paste0(sampleName, ".bam")
    paired = "yes"
  }
  
  gtfFile = param$ezRef["refFeatureFile"]
  gtfFileFlattened = dexseqFlattenGtf(gtfFile, param)  
  strandOpt = switch(param$strandMode, sense="yes", antisense="reverse", both="no", 
                     stop("unsupported strand mode: ", param$strandMode))
  cmd = paste(SAMTOOLS, "view", bamFile, "|", DEXSEQ_COUNT, "-p", paired, "-s", strandOpt, gtfFileFlattened, "-", 
              output$"Count [File]")
  ezSystem(cmd)
  return("Success")
}


