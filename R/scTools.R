###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

getCellCycle <- function(counts, refBuild=c("Homo_sapiens/Ensembl",
                                            "Mus_musculus/Ensembl")){
  require(scran)
  # The training data is only available for Hsap and Mmus Ensembl
  refBuild <- match.arg(refBuild)
  
  if (startsWith(param$refBuild, "Homo_sapiens/Ensembl")) {
    trainData = readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                    package = "scran", mustWork=TRUE))
  } else if (startsWith(param$refBuild, "Mus_musculus/Ensembl")) {
    trainData = readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                    package = "scran", mustWork=TRUE))
  }
  cellCycleData <- cyclone(counts, trainData)
  cellPhase <- tibble(Name = colnames(counts),
                      Phase = cellCycleData$phases)
  return(cellPhase)
}