###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDADA2Step2Dataset = function(input=NA, output=NA, param=NA, 
                                      htmlFile="00index.html"){
  
  require(dada2)
  dataset = input$meta
  sampleName = "merged"
  databaseParam <- param$database
  if (databaseParam == "silva") {
    database <- SILVA_DB_DADA2
  } else if (databaseParam == "RDP") {
    database <- RDP_DB_DADA2
  }  else if (databaseParam == "greenGenes") {
    database <- GREENGENES_DB_DADA2
  }
  ### create a merged object, remove chimeras and assign taxa 
  RfilesWithSeqTabs <- input$getFullPaths("RObjectWithSeqTab")
  mergeDADA2Obj <- DADA2mergeSeqTabs(RfilesWithSeqTabs, database)
  
  ### select groups from dataset
  dataset$Name <- rownames(dataset)
  colsToKeep <- c("Name",grep("Factor",colnames(dataset), value = T))
  designMatrix <- dataset[,colsToKeep]
  colnames(designMatrix) <- gsub(" \\[Factor\\]","",colnames(designMatrix))

  ##  output files

  ### Files needed for Phyloseq
  # taxonomy
  newOTUsToTaxFileName <- basename(output$getColumn("OTUsToTaxonomyFile"))
  taxaOTUs <- data.frame(mergeDADA2Obj$taxaObj, stringsAsFactors = F)
  taxaOTUs$OTU <- paste0("OTU",seq(1:nrow(taxaOTUs)))
  rownames(taxaOTUs) <- NULL
  write.table(taxaOTUs,newOTUsToTaxFileName,
              row.names = F, col.names = T, quote = F,sep = "\t")
  # OTU count
  newOTUsToCountFileName <- basename(output$getColumn("OTUsCountTable"))
  countOTUs <- data.frame(mergeDADA2Obj$fullTableOfOTUsNoChimObj, stringsAsFactors = F)
  colnames(countOTUs) <- paste0("OTU",seq(1:ncol(countOTUs)))
  countOTUs$sample <- rownames(countOTUs)
  rownames(countOTUs) <- NULL
  write.table(countOTUs,newOTUsToCountFileName,
              row.names = F, col.names = T, quote = F,sep = "\t")
  ## design Matrix 
  if (param$group){
    designMatrixFile <-  basename(output$getColumn("sampleDescriptionFile"))
    write.table(designMatrix,designMatrixFile,row.names = F, col.names = T, quote = F,sep = "\t")
}
}
##' @template app-template
##' @templateVar method ezMethodDADA2Step2Dataset()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppDADA2Step2Dataset <-
  setRefClass("EzAppDADA2Step2Dataset",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDADA2Step2Dataset
                  name <<- "EzAppDADA2Step2Dataset"
                  appDefaults <<- rbind(database = ezFrame(Type="character",  
                                                                 DefaultValue="silva",
                                                           Description="16S database")
                                      
                  )
                }
              )
  )
