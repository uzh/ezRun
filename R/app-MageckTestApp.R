ezMethodMageckTest = function(input=NA, output=NA, param=NA){
  require(Herper)
  require(stringr)
  require(MAGeCKFlute)
  require(clusterProfiler) 
  require(ggplot2)
  require(limma)
    
  # Loading the variables
  dataset <- input$meta
  dir.create(param$comparison, showWarnings=FALSE)
  outputPrefix <- file.path(param$comparison, output$getNames())
  
  mergedCountFileName <- paste0(output$getNames(), ".merged.count.tsv")
  mergedCountFileLoc <- file.path(param$comparison, mergedCountFileName)
  sampleNames <- rownames(dataset)
  
  # Combining the count files into a single count file
  countFilePaths <- file.path(param$dataRoot, dataset$`Count [File]`)
  counts <- lapply(countFilePaths, data.table::fread)
  mergeCounts <- counts %>% reduce(inner_join, by="sgRNA")
  mergeCounts[['Gene']] = NULL
  mergeCounts <- mergeCounts %>% rename("Gene"="Gene.x") %>% select(!starts_with("Gene."))
  
  # Rename the sample columns to the actual sample names
  sampleColumns <- !(colnames(mergeCounts) %in% c("sgRNA", "Gene"))
  colnames(mergeCounts)[sampleColumns] <- sampleNames
  mergeCounts$Gene <- gsub(' ', '_', mergeCounts$Gene) #bug in Mageck regarding Spaces in GeneNames
  
  ezWrite.table(mergeCounts, file=mergedCountFileLoc, row.names=FALSE)
  
  # We give the design of the experiment as indices corresponding to the
  # columns (skipping the first 2 positions) which are sample vs ref groups
  fullColumnName <- paste(param$grouping, "[Factor]")
  #which(input$getColumn(param$grouping) == param$refGroup)-1
  rId <- paste(which(dataset[[fullColumnName]] == param$refGroup) - 1, collapse=",")
  sId <- paste(which(dataset[[fullColumnName]] == param$sampleGroup) - 1, collapse=",")
  
  opt <- c(
          "test",
          "-k",
          mergedCountFileLoc,
          "-t",
          sId,
          "-c",
          rId,
          "-n",
          outputPrefix,
          "--pdf-report",
          as.vector(str_split(param$cmdOptions, "\ +", simplify = TRUE)))
  
  # Load the conda environment
  local_CondaEnv("gi_mageck", pathToMiniConda = "/usr/local/ngseq/miniforge3")
  
  ctrlFile <- list.files(param$libName, pattern = 'MAGeCK_Ctrl.csv$', full.names = TRUE)
  if(length(ctrlFile) == 1L){
    opt <- c(opt, "--control-sgrna", ctrlFile)
  }
  # Execute the command
  system2("mageck", args=opt)
  
  # add official gene symbol to gene_summary file for human samples
  if(param$species %in% c('hsa','mmu')){
    dat <- ezRead.table(file.path(param$comparison,paste0(param$comparison, '.gene_summary.txt')), row.names=NULL)
    dat[['GeneSymbol_Addgene']] = dat[['id']]
        for (j in 1:nrow(dat)){
            gene <- c()
            if(param$species == 'hsa'){
                gene <- alias2Symbol(dat$id[j], species = "Hs")
            } else if(param$species == 'mmu') {
                gene <- alias2Symbol(dat$id[j], species = "Mm")
            }
            if(length(gene)==1L)
            dat[['id']][j] <- gene
        }
    dat <- dat[!duplicated(dat$id),]
    ezWrite.table(dat, file.path(param$comparison,paste0(param$comparison, '.gene_summary.txt')), row.names = FALSE)
  }
  
  # We convert the raw outputs to xlsx files
  lapply(c(".sgrna_summary", ".gene_summary"), function(fileComp) {
    dat <- ezRead.table(paste0(outputPrefix, fileComp, ".txt"), row.names=NULL)
    writexl::write_xlsx(dat, paste0(outputPrefix, fileComp, ".xlsx"))
  })
  
  #run MAGECK FLUTE
  setwd(param$comparison)
  file1 =  paste0(param$comparison, '.gene_summary.txt')
  file2 =  paste0(param$comparison, '.sgrna_summary.txt')
  
  out <- tryCatch(FluteRRA(file1, file2, proj="output", organism=param$species, outdir = "./", omitEssential = FALSE), error = function(e) return(message('Error in running MAGECK FLUTE')))
  
  saveRDS(param, 'param.rds')
  makeRmdReport(param=param, output=output,
                rmdFile = "MageckTest.Rmd", reportTitle = paste0(param$name))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMageckTest(input=NA, output=NA, param=NA)
##' @description Use this reference class to run Mageck Test
##' @author Falko Noé
EzAppMageckTest <-
  setRefClass("EzAppMageckTest",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMageckTest
                  name <<- "EzAppMageckTest"
                  appDefaults <<- rbind(
                    outputDir = ezFrame(
                      Type = "character",
                      DefaultValue = ".",
                      Description = "Output directory"
                    )
                  )
                }
              )
  )
