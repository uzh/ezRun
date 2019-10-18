###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDADA2Step1Sample = function(input=NA, output=NA, param=NA, 
                                     htmlFile="00index.html"){
  
  require(dada2)
  require(purrr)
  require(phyloseq)
  dataset = input$meta
  sampleNames = input$getNames() 
  databaseParam <- param$database
  kingdomParam <- param$Kingdom
  if (databaseParam == "silva") {
    if (kingdomParam == "Bacteria"){
    database <- SILVA_BACTERIA_DADA2
    }else if (kingdomParam == "Archea"){
    database <- SILVA_ARACHAEA_DADA2
    }else if (kingdomParam == "Eukaryota"){
    database <- SILVA_EUKARYOTA_DADA2
    }else if (kingdomParam == "All"){
      database <- SILVA_ALL_DADA2    
    }
  } else if (databaseParam == "RDP" & kingdomParam == "Bacteria" ) {
    database <- RDP_DB
  } else if (databaseParam == "greenGenes" & kingdomParam == "Bacteria") {
    database <- GREENGENES_DB
  } else {
    stop("Currently from RDP and greenGenes we bacterial databases.")
  }
  minLen <- param$minLen
  isPaired <- param$paired
  concat <- param$concatenateReads
  ### read fastq files and prepare inputs for DADA2
  ### if reads are paired, they are first joined. The DADA2 inbulit joining works only 
  ### if there is only one V-region (almost never the case)
  file1PathInDataset <- input$getFullPaths("Read1")
  if(isPaired){
    file2PathInDataset <- input$getFullPaths("Read2")
    fastqJoin="/usr/local/ngseq/src/ea-utils.1.1.2-686/fastq-join"

    fastqJoinFun <- function(x,y,z){
    joinedFileName <- paste0(z, ".temp.")
    fastqJoinCmd <- paste(fastqJoin, x,y, "-o", joinedFileName)
    ezSystem(fastqJoinCmd)
    joinedFileName <- paste0(joinedFileName,"join")
    joinedFile <- file.path(getwd(),joinedFileName)
    return(joinedFile)
    }
    listOfJoinedFiles <- mapply(fastqJoinFun,file1PathInDataset,file2PathInDataset,
                                sampleNames)
    DADA2mainSeqTabObj <- DADA2CreateSeqTab(sampleNames = sampleNames,
                                          minLen = minLen,
                                          file1PathInDataset = listOfJoinedFiles,
                                          database)
  }else{
    DADA2mainSeqTabObj <- DADA2CreateSeqTab(sampleNames= sampleNames,
                                            minLen = minLen,
                                            file1PathInDataset = file1PathInDataset,
                                            database)
  }
  
  ## rename output files
  # taxonomy
  newOTUsToTaxFileName <- basename(output$getColumn("OTUsToTaxonomyFile"))
  taxaOTUs <- data.frame(DADA2mainSeqTabObj$taxaObj, stringsAsFactors = F)
  taxaOTUs$OTU <- paste0("OTU",seq(1:nrow(taxaOTUs)))
  rownames(taxaOTUs) <- NULL
  write.table(taxaOTUs,newOTUsToTaxFileName,
              row.names = F, col.names = T, quote = F,sep = "\t")
  # OTU count
  newOTUsToCountFileName <- basename(output$getColumn("OTUsCountTable"))
  countOTUs <- data.frame(DADA2mainSeqTabObj$fullTableOfOTUsNoChimObj, stringsAsFactors = F)
  colnames(countOTUs) <- paste0("OTU",seq(1:ncol(countOTUs)))
  countOTUs$sample <- rownames(countOTUs)
  rownames(countOTUs) <- NULL
  write.table(countOTUs,newOTUsToCountFileName,
              row.names = F, col.names = T, quote = F,sep = "\t")
  ## design Matrix 
  if (param$group){
    factorCols <- grep("Factor",colnames(dataset))
    designMatrix <- data.frame(dataset[,factorCols])
    colnames(designMatrix) <- gsub(" \\[Factor\\]","",colnames(dataset)[factorCols])
    rownames(designMatrix) <- rownames(dataset)
    designMatrixFile <-  basename(output$getColumn("sampleDescriptionFile"))
    write.table(designMatrix,designMatrixFile,row.names = F, col.names = T, quote = F,sep = "\t")
  }
  
  ## create phyloseqObject
  phyloseqObjectRdata <-  basename(output$getColumn("RObjectPhyloseq"))
  if (param$group){
  phyloseqObject <- phyloseq(otu_table(DADA2mainSeqTabObj$fullTableOfOTUsNoChimObj, taxa_are_rows=FALSE), 
                 sample_data(designMatrix), 
                 tax_table(DADA2mainSeqTabObj$taxaObj))
  }else{
  phyloseqObject <- phyloseq(otu_table(DADA2mainSeqTabObj$fullTableOfOTUsNoChimObj, taxa_are_rows=FALSE), 
                               tax_table(DADA2mainSeqTabObj$taxaObj))
  }
  saveRDS(phyloseqObject,phyloseqObjectRdata)
  
}

##' @template app-template
##' @templateVar method ezMethodDADA2Step1Sample()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppDADA2Step1Sample <-
  setRefClass("EzAppDADA2Step1Sample",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDADA2Step1Sample
                  name <<- "EzAppDADA2Step1Sample"
                  appDefaults <<- rbind(minLen= ezFrame(Type="integer",  DefaultValue="150",Description="Min length")
                  )
                }
              )
  )












