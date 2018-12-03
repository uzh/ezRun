###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodUparse = function(input=NA, output=NA, param=NA, 
                                                 htmlFile="00index.html"){
### run UPARSE and reformat output for phylotseq 
  
library(plyr)
library(dplyr)

  sampleName = input$getNames()
  sampleNameString = paste("\"", paste(sampleName, collapse = " "),"\"")
  file1PathInDatset <- input$getFullPaths("Read1")
  for (k in 1:length(sampleName)){
    cpCmd1 <- paste0("gunzip -c ", file1PathInDatset[k], "  > ", sampleName[k],"_R1",".fastq")
  ezSystem(cpCmd1)  
  if(param$paired){
    file2PathInDatset <- input$getFullPaths("Read2")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset[k], "  > ", sampleName[k],"_R2",".fastq")
    ezSystem(cpCmd2)
  }
  }
uparseTemplScript <- USEARCH_TEMPLATE_SCRIPT
uparseToBeExec <- "uparse.sh"
##update template
updateTemplateScriptCmd <- paste0("sed -e s/\"SAMPLE_LIST\"/", sampleNameString, "/g",
                               " -e s/\"MAX_EE\"/", param$fastqErrorMax, "/g ",
                               uparseTemplScript, " >",
                               uparseToBeExec)
ezSystem(updateTemplateScriptCmd)
## run script
uparseCmd <- paste("bash",uparseToBeExec)
ezSystem(uparseCmd)

### 1 OTU abundance table 
inputTabFile <- "all.OTU.tab.txt"
otuTabRaw <- ezRead.table(inputTabFile, header = T, stringsAsFactors = F, 
                        check.names = F)
phyloseqAppFormatOtuAbundTable <- data.frame(t(otuTabRaw), stringsAsFactors = F, 
                        check.names = F, check.rows = F)
nOTUS <- ncol(phyloseqAppFormatOtuAbundTable)
phyloseqAppFormatOtuAbundTable$label <- "0.03"
phyloseqAppFormatOtuAbundTable$numOtus <- nOTUS
phyloseqAppFormatOtuAbundTable$Group <- rownames(phyloseqAppFormatOtuAbundTable)

### 2 OTU-taxonomy table 
inputTaxFile <- paste0("all.OTU.taxonomy.txt")
otuTaxonomyRaw <- read.delim(inputTaxFile, header = F, stringsAsFactors = F)
otuTaxonomyRaw <- otuTaxonomyRaw[,1:2]
names(otuTaxonomyRaw) <- c("OTU","Taxonomy")
bactRowsToKeep <- grep("d:Bacteria",otuTaxonomyRaw$Taxonomy)
otuTaxonomyRawOnlyBact <- otuTaxonomyRaw[bactRowsToKeep,]
taxRanks <- c("d:","p:","c:","o:","f:","g:","s:")
createAppropriateRows <- function(inputTaxTable){
  x <- inputTaxTable[2]
  x1 <- strsplit(x,",")[[1]]
  taxClassif <- vector()
  taxClassif[1] <- "Bacteria"
  for (k in 2:7){
    tax <- taxRanks[k]
    isDomThere <- grep(tax,x1, value = T)
    if (length(isDomThere) >0){
      taxClassif[k] <- unlist(strsplit(gsub(tax,"",isDomThere),"\\("))[1]
    } else {
      k1=k-1
      alreadyUnclassified <- grep(taxClassif[k1],"unclassified")
      if (length(alreadyUnclassified) > 0 ) {
        taxClassif[k] <- taxClassif[k1]
      } else{
      taxClassif[k] <- paste(taxClassif[k1],"unclassified",sep = "_")
    }
    }
  }
 return(c(inputTaxTable[1],paste(taxClassif, collapse = ";")))
}

### compact everything 
formattedRows <- apply(otuTaxonomyRawOnlyBact,1,createAppropriateRows)
formattedDF <- data.frame(t(formattedRows), stringsAsFactors = F)
colnames(formattedDF) <- c("OTU","Taxonomy")
OTUsToKeep <- colnames(phyloseqAppFormatOtuAbundTable)[
  colnames(phyloseqAppFormatOtuAbundTable)%in%formattedDF$OTU]
finalOtuAbundTable <- phyloseqAppFormatOtuAbundTable[,c("Group", "label","numOtus",OTUsToKeep)]
finalOtuTaxon <- formattedDF[formattedDF$OTU%in%OTUsToKeep,]
OTUsTotAbund <- apply(phyloseqAppFormatOtuAbundTable[,OTUsToKeep],2,sum)
finalOtuTaxon$Size <- OTUsTotAbund

## sort final files 
newOtuAbundFile <- basename(output$getColumn("OTUsCountTable"))
newTaxAbundFile <- basename(output$getColumn("OTUsToTaxonomyFile"))
write.table(finalOtuAbundTable,newOtuAbundFile,row.names = F, col.names = T, quote = F,sep = "\t")
write.table(finalOtuTaxon,newTaxAbundFile,row.names = F, col.names = T, quote = F,sep = "\t")

## eventual design Matrix 
if (param$Group){
designMatrix <- data.frame(Name = sampleName,Group=input$getColumn("Group"), 
                           check.names = F)
designMatrixFile <-  basename(output$getColumn("sampleDescriptionFile"))
write.table(designMatrix,designMatrixFile,row.names = F, col.names = T, quote = F,sep = "\t")
}

}

##' @template app-template
##' @templateVar method ezMethodUparse()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppUparse <-
  setRefClass("EzAppUparse",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodUparse
                  name <<- "EzAppUparse"
                  appDefaults <<- rbind(fastqErrorMax = ezFrame(Type="numeric",
                                                                 DefaultValue="1",
                                                                 Description="Max EE in fastx_truncate")

                  )
                }
              )
  )

