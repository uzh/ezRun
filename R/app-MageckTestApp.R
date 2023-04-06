ezMethodMageckTest = function(input=NA, output=NA, param=NA){
  require(Herper)
  require(stringr)
  require(MAGeCKFlute)
  require(clusterProfiler) 
  require(ggplot2)
    
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
  local_CondaEnv("mageckenv", pathToMiniConda = "/usr/local/ngseq/miniconda3")
  
  ctrlFile <- list.files(param$libName, pattern = 'MAGeCK_Ctrl.csv$', full.names = TRUE)
  if(length(ctrlFile) == 1L){
    opt <- c(opt, "--control-sgrna", ctrlFile)
  }
  # Execute the command
  system2("mageck", args=opt)
  
  # We convert the raw outputs to xlsx files
  lapply(c(".sgrna_summary", ".gene_summary"), function(fileComp) {
    dat <- ezRead.table(paste0(outputPrefix, fileComp, ".txt"), row.names=NULL)
    writexl::write_xlsx(dat, paste0(outputPrefix, fileComp, ".xlsx"))
  })
  
  #run MAGECK FLUTE
  setwd(param$comparison)
  file1 =  paste0(param$comparison, '.gene_summary.txt')
  file2 =  paste0(param$comparison, '.sgrna_summary.txt')
  FluteRRA(file1, file2, proj="output", organism=param$species, outdir = "./", omitEssential = FALSE)
  
  
#  rmdFile <- paste0(param$comparison, '.report.Rmd')
#  htmlFile <- sub('.Rmd', '.html', rmdFile)
  
 # rmarkdown::render(
  #    input = rmdFile, envir = new.env(),
  #    output_dir = ".", output_file = htmlFile, quiet = TRUE
  #)
  
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
