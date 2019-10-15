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
  dataset = input$meta
  sampleNames = input$getNames() 
  databaseParam <- param$database
  if (databaseParam == "silva") {
    database <- SILVA_DB_DADA2
  } else if (databaseParam == "RDP") {
    database <- RDP_DB_DADA2
  }  else if (databaseParam == "greenGenes") {
    database <- GREENGENES_DB_DADA2
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

    fastqJoin <- function(x,y,z){
    joinedFileName <- paste0(z, ".temp.")
    fastqJoinCmd <- paste(fastqJoin, x,y, "-o", joinedFileName)
    ezSystem(fastqJoinCmd)
    joinedFileName <- paste0(joinedFileName,"join")
    joinedFile <- file.path(getwd(),joinedFileName)
    return(joinedFile)
    }
    listOfJoinedFiles <- mapply(fastqJoin,file1PathInDataset,file2PathInDataset,
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
  if (param$Group){
    designMatrixFile <-  basename(output$getColumn("sampleDescriptionFile"))
    write.table(designMatrix,designMatrixFile,row.names = F, col.names = T, quote = F,sep = "\t")
  }
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












