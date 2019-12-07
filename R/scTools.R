###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

getCellCycle <- function(counts, refBuild){
  require(scran)
  require(tibble)
  require(dplyr)
  # The training data is only available for Hsap and Mmus Ensembl
  if(startsWith(refBuild, "Homo_sapiens")){
    species <- "human"
    hasTrainData <- TRUE
  }else if(startsWith(refBuild, "Mus_musculus")){
    species <- "mouse"
    hasTrainData <- TRUE
  }else{
    hasTrainData <- FALSE
  }
  if(isTRUE(hasTrainData)){
    trainData <- readRDS(system.file("exdata",
                                     paste0(species, "_cycle_markers.rds"), 
                                     package = "scran", mustWork=TRUE))
    cellCycleData <- cyclone(counts, trainData)
    cellPhase <- tibble(Name = colnames(counts),
                        Phase = cellCycleData$phases)
    cellPhase <- bind_cols(cellPhase, cellCycleData$scores)
  }else{
    cellPhase <- tibble(Name = colnames(counts),
                        Phase = NA)
  }
  return(cellPhase)
}

getPerplexity <- function(n){
  ifelse(n > 200, 30, 10)
}
