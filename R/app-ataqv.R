###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppAtaqv <-
  setRefClass("EzAppAtaqv",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodAtaqv
                  name <<- "EzAppAtaqv"
                }
              )
  )

ezMethodAtaqv = function(input=NA, output=NA, param=NA, 
                         htmlFile="index.html"){
  require(GenomeInfoDb)
  require(rtracklayer)
  require(taxize)
  
  bamFiles <- input$getFullPaths("BAM")
  
  ## Get the autosomes
  #species <- unique(input$getColumn("Species"))
  #stopifnot(length(species)==1L)
  ## in the format of Homo_sapiens
  speciesLatin <- strsplit(param$refBuild, split="/")[[1]][1]
  #speciesLatin <- sub(" ", "_", sub(" \\(.*$", "", species))
  speciesNormal <- strsplit(sci2comm(scinames=speciesLatin, db='ncbi')[[1]],
                            split=" ")[[1]]
  speciesNormal <- tail(speciesNormal, 1)
  #sub("\\)$", "", sub(".*\\(", "", species))
  
  if(speciesLatin %in% names(genomeStyles())){
    # Use the autosomes from GenomeInfoDb
    # Guess chr format
    if(grepl("(UCSC|GENCODE)", param$refBuild, ignore.case = TRUE)){
      chrFormat <- "UCSC"
      mitChr <- "chrM"
    }else if(grepl("Ensembl", param$refBuild, ignore.case = TRUE)){
      chrFormat <- "Ensembl"
      mitChr <- "MT"
    }else if(grepl("NCBI", param$refBuild, ignore.case = TRUE)){
      chrFormt <- "NCBI"
      mitChr <- "MT"
    }else{
      stop("Unknown refBuild: ", param$refBuild)
    }
    message("Using autosomes from GenomeInfoDb package.")
    autosomes <- extractSeqlevelsByGroup(speciesLatin, style=chrFormat, 
                                         group="auto")
  }else{
    message("Use all the chromosomes from the genome!")
    autosomes <- ezBamSeqNames(bamFiles[1])
    mitChr <- "chrM"
  }
  automsomesFn <- tempfile(pattern="autosomes-", fileext=".txt")
  writeLines(autosomes, con=automsomesFn)
  
  ## TSS from gtf
  tx <- ezFeatureAnnotation(param, dataFeatureType = "transcript")
  tss <- promoters(makeGRangesFromDataFrame(tx), upstream=0, downstream = 1)
  tssFn <- tempfile(pattern="tss-", fileext = ".txt")
  export.bed(tss, con=tssFn)
  
  ## ataqv
  jsonFns <- paste0(names(bamFiles), ".ataqv.json")
  cmd <- paste("ataqv", "--autosomal-reference-file", automsomesFn,
               "--mitochondrial-reference-name", mitChr,
               "--metrics-file", jsonFns,
               "--tss-file", tssFn)
  
  if(any(grepl("BED", input$colNames))){
    cmd <- paste(cmd, "--peak-file", input$getFullPaths("BED"))
  }
  cmd <- paste(cmd, speciesNormal, bamFiles)
  ezMclapply(cmd, ezSystem, mc.cores=param$cores)
  
  ## mkarv
  cmd <- paste("mkarv", param$name, paste(jsonFns, collapse=" "))
  ezSystem(cmd)
  
  file.remove(jsonFns)
  
  return("Success")
}
